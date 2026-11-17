import hashlib
import json
import pathlib
import re
import time
import xml.etree.ElementTree as ET

import requests

from . import clients as K
from . import config as C

VEP_URL = "https://rest.ensembl.org/vep/human/region"
ACMG = re.compile(r"\b(PVS1|PS[1-4]|PM[1-6]|PP[1-5]|BA1|BS[1-4]|BP[1-7])"
                  r"(?:[_ ](Very ?Strong|Strong|Moderate|Supporting|Stand[- ]?Alone))?\b")
FUNCTIONAL_METHODS = {"in vitro", "in vivo", "research"}
VEP_ERRORS: list[dict] = []


def genpept_features(gene: str, np_acc: str) -> list[dict]:
    url = C.EUTILS + "efetch.fcgi"
    p = K._ncbi_params(db="protein", id=np_acc, rettype="gp", retmode="xml")
    xml = K._get(url, p).text
    time.sleep(C.NCBI_DELAY)
    K._log("NCBI protein", "E-utilities efetch (GenPept XML)", np_acc, K._public_url(url, p), 1)
    root = ET.fromstring(xml)
    rows = []
    for f in root.iter("GBFeature"):
        key = f.findtext("GBFeature_key")
        if key in ("source",):
            continue
        q = {}
        for qq in f.iter("GBQualifier"):
            q.setdefault(qq.findtext("GBQualifier_name"), []).append(qq.findtext("GBQualifier_value") or "")
        ivs = [(i.findtext("GBInterval_from") or i.findtext("GBInterval_point"),
                i.findtext("GBInterval_to") or i.findtext("GBInterval_point")) for i in f.iter("GBInterval")]
        rows.append(dict(
            gen=gene, NP=np_acc, tipo_feature=key, ubicacion=f.findtext("GBFeature_location"),
            inicio=ivs[0][0] if ivs else None, fin=ivs[-1][1] if ivs else None,
            nombre=";".join(q.get("region_name", []) + q.get("site_type", []) + q.get("product", [])),
            nota=" | ".join(q.get("note", [])), db_xref=";".join(q.get("db_xref", [])),
            fuente="NCBI Protein GenPept", enlace=C.LINKS["protein"].format(np_acc)))
    return rows


def clinvar_submissions(gene: str, variation_ids: list[str], batch=100) -> list[dict]:
    url = C.EUTILS + "efetch.fcgi"
    rows = []
    for i in range(0, len(variation_ids), batch):
        chunk = variation_ids[i:i + batch]
        p = K._ncbi_params(db="clinvar", id=",".join(chunk), rettype="vcv", is_variationid="true")
        r = K._post(url, p)
        time.sleep(C.NCBI_DELAY)
        root = ET.fromstring(r.text)
        for va in root.iter("VariationArchive"):
            vid, vcv = va.get("VariationID"), va.get("Accession")
            for ca in va.iter("ClinicalAssertion"):
                acc = ca.find("ClinVarAccession")
                cl = ca.find("Classification")
                comment = " ".join((c.text or "") for c in ca.iter("Comment"))
                methods = sorted({(m.text or "").strip() for m in ca.iter("MethodType")})
                crit = sorted({m.group(0).replace(" ", "_") for m in ACMG.finditer(comment)})
                base = sorted({c.split("_")[0] for c in crit})
                pmids = sorted({c.text for c in ca.iter("ID") if c.get("Source") == "PubMed" and c.text})
                rows.append(dict(
                    gen=gene, variation_id=vid, VCV=vcv,
                    SCV=acc.get("Accession") if acc is not None else None,
                    remitente=acc.get("SubmitterName") if acc is not None else None,
                    clasificacion_SCV=cl.findtext("GermlineClassification") if cl is not None else None,
                    revision_SCV=cl.findtext("ReviewStatus") if cl is not None else None,
                    fecha_evaluacion=cl.get("DateLastEvaluated") if cl is not None else None,
                    metodos=";".join(methods),
                    metodo_funcional=bool(set(methods) & FUNCTIONAL_METHODS),
                    criterios_ACMG=";".join(crit),
                    cita_PS3_BS3=bool({"PS3", "BS3"} & set(base)),
                    cita_PP3_BP4=bool({"PP3", "BP4"} & set(base)),
                    PubMed=";".join(pmids),
                    comentario=comment.strip()[:2000],
                ))
        K._log("ClinVar", "E-utilities efetch (VCV XML)", f"{gene}: {len(chunk)} VariationID",
               f"{url}?db=clinvar&rettype=vcv&is_variationid=true&id=<{len(chunk)} ids>", len(chunk))
    return rows

def spdi_to_vcf(chrom, spdi):
    """'NC_000004.12:987861:A:C' -> '4 987862 . A C . . .' (solo sustituciones de igual longitud)."""
    try:
        _, pos, ref, alt = spdi.split(":")
    except (AttributeError, ValueError):
        return None
    if not ref or not alt or len(ref) != len(alt):
        return None
    return f"{chrom} {int(pos) + 1} {chrom}_{int(pos) + 1}_{ref}_{alt} {ref} {alt} . . ."


def _vep_post(lines, retries=4):
    body = dict(variants=lines, AlphaMissense=1, REVEL=1, CADD=1, SpliceAI=1, mane=1,
                canonical=1, numbers=1, protein=1, hgvs=1)
    for k in range(retries + 2):
        try:
            r = requests.post(VEP_URL, json=body, timeout=(20, 120),
                              headers={"Content-Type": "application/json", "Accept": "application/json"})
        except requests.exceptions.RequestException as e:   # cortes de red / proxy transitorios
            print(f"  VEP reintento {k + 1}: {type(e).__name__}", flush=True)
            time.sleep(15 * (k + 1)); continue
        if r.status_code == 200:
            return r.json()
        if r.status_code in (429, 500, 502, 503, 504):
            time.sleep(float(r.headers.get("Retry-After", 5 * (k + 1)))); continue
        if r.status_code == 400:      # alguna variante rechazada: dividir el lote para aislarla
            if len(lines) == 1:
                VEP_ERRORS.append(dict(linea=lines[0], error=r.text[:300]))
                return []
            mid = len(lines) // 2
            return _vep_post(lines[:mid]) + _vep_post(lines[mid:])
        r.raise_for_status()
    raise RuntimeError("VEP no respondió tras varios reintentos")


def vep_annotate(variants: list[dict], batch=200, cache_dir=None) -> list[dict]:
    """variants: dicts con gen, chr, spdi, MANE_ENST. Devuelve una fila por variante (consecuencia MANE)."""
    cache = pathlib.Path(cache_dir) if cache_dir else None
    if cache: cache.mkdir(parents=True, exist_ok=True)
    todo = {}
    for v in variants:
        line = spdi_to_vcf(v["chr"], v.get("spdi"))
        if line:
            k = line.split()[2]
            if k in todo: todo[k][1].append(v)
            else: todo[k] = (line, [v])
    keys = list(todo)
    out = []
    for i in range(0, len(keys), batch):
        chunk = keys[i:i + batch]
        h = hashlib.md5("|".join(chunk).encode()).hexdigest()[:12]
        cf = cache / f"vep_{h}.json" if cache else None
        if cf and cf.exists():
            res = json.loads(cf.read_text())
        else:
            lines = [todo[k][0] for k in chunk]
            res = []
            for j in range(0, len(lines), 50):     # sub-lotes cortos: respuestas rápidas y reintentos baratos
                res += _vep_post(lines[j:j + 50])
            if cf: cf.write_text(json.dumps(res))
            K._log("Ensembl VEP", "REST POST /vep/human/region", f"{len(chunk)} variantes (VCF GRCh38)",
                   VEP_URL + " (AlphaMissense, REVEL, CADD, SpliceAI, mane, hgvs)", len(res))
        for rec in res:
            key = rec.get("id")
            if key not in todo:
                continue
            for v in todo[key][1]:
                enst = (v.get("MANE_ENST") or "").split(".")[0]
                tcs = rec.get("transcript_consequences", [])
                # Solo el ENST MANE del gen del panel (nunca el MANE de un gen vecino)
                tc = next((t for t in tcs if t.get("transcript_id") == enst), {})
                sp = tc.get("spliceai") or {}
                am = tc.get("alphamissense") or {}
                am_tx = tc.get("transcript_id") if am else None
                if not am and tc:   # AlphaMissense se calculó sobre la isoforma UniProt (p. ej. HGSNAT 663 aa)
                    alt = next((t for t in tcs if t.get("gene_symbol") == v["gen"] and t.get("alphamissense")
                                and t.get("amino_acids") == tc.get("amino_acids")), None)
                    if alt:
                        am, am_tx = alt["alphamissense"], f"{alt['transcript_id']} (p.{alt.get('protein_start')})"
                out.append(dict(
                    gen=v["gen"], variation_id=v.get("variation_id"), id_VEP=key,
                    transcrito_VEP=tc.get("transcript_id"), MANE_RefSeq=tc.get("mane_select"),
                    consecuencia=";".join(tc.get("consequence_terms", [])),
                    hgvsc=tc.get("hgvsc"), hgvsp=tc.get("hgvsp"),
                    pos_proteina=tc.get("protein_start"), aminoacidos=tc.get("amino_acids"),
                    SIFT=tc.get("sift_score"), SIFT_pred=tc.get("sift_prediction"),
                    PolyPhen2=tc.get("polyphen_score"), PolyPhen2_pred=tc.get("polyphen_prediction"),
                    CADD_phred=tc.get("cadd_phred"), REVEL=tc.get("revel"),
                    SpliceAI_max=max([sp.get(k, 0) or 0 for k in ("DS_AG", "DS_AL", "DS_DG", "DS_DL")]) if sp else None,
                    AlphaMissense=am.get("am_pathogenicity"), AlphaMissense_clase=am.get("am_class"),
                    AlphaMissense_transcrito=am_tx,
                ))
        print(f"  VEP {min(i + batch, len(keys))}/{len(keys)}", flush=True)
    return out


def vep_release():
    r = requests.get("https://rest.ensembl.org/info/software", headers={"Content-Type": "application/json"}, timeout=30)
    K.SOURCE_VERSIONS["Ensembl VEP release"] = r.json().get("release")
    return K.SOURCE_VERSIONS["Ensembl VEP release"]


def alphamissense_uniprot(gene: str, acc: str) -> dict:
    """{'M1A': (0.2651, 'LBen'), ...} en numeración UniProt. Fuente alternativa cuando VEP no puntúa el MANE."""
    api = f"https://alphafold.ebi.ac.uk/api/prediction/{acc}"
    meta = K._get(api).json()[0]
    K.SOURCE_VERSIONS["AlphaFold DB model version"] = meta.get("latestVersion")
    url = f"https://alphafold.ebi.ac.uk/files/AF-{acc}-F1-aa-substitutions.csv"
    txt = K._get(url).text
    out = {}
    for line in txt.splitlines()[1:]:
        parts = line.strip().split(",")
        if len(parts) >= 2 and parts[1]:
            out[parts[0]] = (float(parts[1]), parts[2] if len(parts) > 2 else None)
    K._log("AlphaFold DB", "files/AF-{acc}-F1-aa-substitutions.csv", f"{gene} ({acc})", url, len(out))
    return out
