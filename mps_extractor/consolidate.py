import re

import pandas as pd

AA3 = dict(Ala="A", Arg="R", Asn="N", Asp="D", Cys="C", Gln="Q", Glu="E", Gly="G", His="H", Ile="I",
           Leu="L", Lys="K", Met="M", Phe="F", Pro="P", Ser="S", Thr="T", Trp="W", Tyr="Y", Val="V",
           Sec="U", Ter="*")
TITLE = re.compile(r"^(NM_\d+\.\d+)\(([^)]+)\):c\.\S+ \(p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})\)")
STARS = {"practice guideline": 4, "reviewed by expert panel": 3,
         "criteria provided, multiple submitters, no conflicts": 2,
         "criteria provided, conflicting classifications": 1, "criteria provided, single submitter": 1}


def _p_from_title(title, mane_nm):
    m = TITLE.match(str(title))
    if not m or m.group(1) != mane_nm:
        return None
    ref, pos, alt = AA3.get(m.group(3)), int(m.group(4)), AA3.get(m.group(5))
    return (ref, pos, alt) if ref and alt else None


def _p_from_vep(row):
    aa, pos = row.get("aminoacidos"), row.get("pos_proteina")
    if not isinstance(aa, str) or "/" not in aa or pd.isna(pos):
        return None
    ref, alt = aa.split("/")
    return (ref, int(pos), alt) if len(ref) == 1 and len(alt) == 1 else None


def build(ids, conc, seqs, clinvar, scv, vep, up_feats, am_up=None):
    am_up = am_up or {}
    ids = ids.set_index("gen"); off = conc.set_index("gen")["desplazamiento_UniProt_menos_NP"].to_dict()
    vep_by = vep.set_index("variation_id").to_dict("index") if len(vep) else {}
    scv_g = scv.groupby("variation_id") if len(scv) else None
    recs = {}

    def get(gen, ref, pos, alt):
        k = (gen, ref, pos, alt)
        if k not in recs:
            o = int(off.get(gen, 0)); pu = pos + o; s = seqs.get(gen, "")
            recs[k] = dict(gen=gen, variante=f"{ref}{pos}{alt}",
                           HGVS_p=f"{ids.loc[gen, 'MANE_NP']}:p.{ref}{pos}{alt}",
                           pos_NP=pos, pos_UniProt=pu,
                           residuo_ref_coincide=(0 < pu <= len(s) and s[pu - 1] == ref),
                           VCV=[], clasificacion_ClinVar=[], estrellas_max=None, n_SCV=0,
                           SCV_metodo_funcional=False, SCV_cita_PS3_BS3=False, SCV_cita_PP3_BP4=False,
                           criterios_ACMG=set(), UniProt_variante=[], UniProt_descripcion=[],
                           UniProt_experimental=False, UniProt_PubMed=set(), mutagenesis=[],
                           fuentes=set())
        return recs[k]

    long_rows = []

    def ev(r, tipo, nivel, fuente, medida, valor, ident=None, enlace=None, detalle=None):
        long_rows.append(dict(gen=r["gen"], variante=r["variante"], HGVS_p=r["HGVS_p"], pos_UniProt=r["pos_UniProt"],
                              tipo_evidencia=tipo, nivel=nivel, fuente=fuente, medida=medida, valor=valor,
                              identificador=ident, enlace=enlace, detalle=detalle))

    sin_mapeo = []
    for _, c in clinvar.iterrows():
        g, vid = c["gen"], str(c["variation_id"])
        v = vep_by.get(int(vid)) or vep_by.get(vid) or {}
        p = _p_from_title(c["nombre_HGVS"], ids.loc[g, "MANE_NM"]) or (_p_from_vep(v) if v else None)
        if not p:
            m = re.match(r"^\S+?\(([^)]+)\)", str(c["nombre_HGVS"]))
            tg = m.group(1) if m else None
            csq = v.get("consecuencia") if v else None
            if tg and tg != g:
                motivo = f"missense en otro gen ({tg}); en {g}: {csq or 'sin anotación VEP'}"
            else:
                motivo = f"sin sustitución simple sobre el NP_ MANE (VEP: {csq or 'sin anotación'})"
            sin_mapeo.append(dict(gen=g, VCV=c["VCV"], nombre_HGVS=c["nombre_HGVS"], gen_del_titulo=tg,
                                  consecuencia_en_gen_MPS=csq, motivo=motivo))
            continue
        r = get(g, *p); r["fuentes"].add("ClinVar")
        r["VCV"].append(c["VCV"]); r["clasificacion_ClinVar"].append(str(c["clasificacion"]))
        st = STARS.get(c["estado_revision"], 0)
        r["estrellas_max"] = max(r["estrellas_max"] or 0, st)
        ev(r, "clínica", "transcrito", "ClinVar (agregada)", "clasificación", c["clasificacion"], c["VCV"], c["enlace"],
           f"{c['estado_revision']} | {c['condiciones']}")
        if scv_g is not None and vid in scv_g.groups:
            sub = scv_g.get_group(vid)
            r["n_SCV"] += len(sub)
            r["SCV_metodo_funcional"] |= bool(sub.metodo_funcional.any())
            r["SCV_cita_PS3_BS3"] |= bool(sub.cita_PS3_BS3.any())
            r["SCV_cita_PP3_BP4"] |= bool(sub.cita_PP3_BP4.any())
            for _, s in sub.iterrows():
                if isinstance(s.criterios_ACMG, str) and s.criterios_ACMG:
                    r["criterios_ACMG"].update(s.criterios_ACMG.split(";"))
                if s.metodo_funcional or s.cita_PS3_BS3:
                    ev(r, "experimental", "transcrito", "ClinVar (submission)",
                       "ensayo funcional citado" if s.cita_PS3_BS3 else f"método {s.metodos}",
                       s.clasificacion_SCV, s.SCV, f"https://www.ncbi.nlm.nih.gov/clinvar/{s.SCV}/" if s.SCV else None,
                       f"{s.remitente} | {s.criterios_ACMG} | PubMed {s.PubMed}")
        if v:
            r["fuentes"].add("VEP")
            for col, nivel in [("SIFT", "transcrito"), ("PolyPhen2", "transcrito"), ("CADD_phred", "transcrito"),
                               ("REVEL", "transcrito"), ("SpliceAI_max", "transcrito"), ("AlphaMissense", "proteína")]:
                if v.get(col) is not None and not pd.isna(v.get(col)):
                    r[col] = v[col]
                    ev(r, "in silico", nivel, f"Ensembl VEP ({col})", col, v[col], v.get("hgvsc"), None,
                       v.get("AlphaMissense_clase") if col == "AlphaMissense" else None)
            r["AlphaMissense_clase"] = v.get("AlphaMissense_clase")
            r["VEP_hgvsp"] = v.get("hgvsp")
            r["VEP_concuerda_p"] = (_p_from_vep(v) == p)

    # UniProt (variantes naturales y mutagénesis)
    uf = up_feats[up_feats.tipo_feature.isin(["Natural variant", "Mutagenesis"])]
    for _, f in uf.iterrows():
        g = f["gen"]; o = int(off.get(g, 0))
        if not isinstance(f["original"], str) or len(f["original"]) != 1 or pd.isna(f["inicio"]):
            continue
        for alt in str(f["alternativa"]).split(";"):
            if len(alt) != 1:
                continue
            r = get(g, f["original"], int(f["inicio"]) - o, alt); r["fuentes"].add("UniProt")
            exp = bool(f["evidencia_experimental"])
            if f["tipo_feature"] == "Natural variant":
                r["UniProt_variante"].append(str(f["id_feature"])); r["UniProt_descripcion"].append(str(f["descripcion"]))
                r["UniProt_experimental"] |= exp
                medida = "variante natural"
            else:
                r["mutagenesis"].append(str(f["descripcion"])); medida = "mutagénesis"
            if isinstance(f["PubMed"], str):
                r["UniProt_PubMed"].update(f["PubMed"].split(";"))
            ev(r, "experimental" if exp else "in silico", "proteína", "UniProtKB", medida, f["descripcion"],
               f["id_feature"], f"https://www.uniprot.org/uniprotkb/{f['UniProt']}/entry", f"{f['ECO']} | PubMed {f['PubMed']}")

    # AlphaMissense desde AlphaFold DB (numeración UniProt) 
    for r in recs.values():
        if r.get("AlphaMissense") is not None and not pd.isna(r.get("AlphaMissense")):
            r["AlphaMissense_fuente"] = "Ensembl VEP"
            continue
        ref, alt = r["variante"][0], r["variante"][-1]
        hit = am_up.get(r["gen"], {}).get(f"{ref}{r['pos_UniProt']}{alt}")
        if hit:
            r["AlphaMissense"], r["AlphaMissense_clase"] = hit
            r["AlphaMissense_fuente"] = "AlphaFold DB (isoforma UniProt)"
            ev(r, "in silico", "proteína", "AlphaFold DB (AlphaMissense)", "AlphaMissense", hit[0],
               f"{ref}{r['pos_UniProt']}{alt} (numeración UniProt)", None, hit[1])

    wide = []
    for r in recs.values():
        r = dict(r)
        for k in ("VCV", "clasificacion_ClinVar", "UniProt_variante", "UniProt_descripcion", "mutagenesis"):
            r[k] = " | ".join(dict.fromkeys(r[k]))
        for k in ("criterios_ACMG", "UniProt_PubMed", "fuentes"):
            r[k] = ";".join(sorted(x for x in r[k] if x))
        wide.append(r)
    cols = ["gen", "variante", "HGVS_p", "pos_NP", "pos_UniProt", "residuo_ref_coincide", "fuentes", "VCV",
            "clasificacion_ClinVar", "estrellas_max", "n_SCV", "SCV_metodo_funcional", "SCV_cita_PS3_BS3",
            "SCV_cita_PP3_BP4", "criterios_ACMG", "UniProt_variante", "UniProt_descripcion", "UniProt_experimental",
            "UniProt_PubMed", "mutagenesis", "SIFT", "PolyPhen2", "CADD_phred", "REVEL", "SpliceAI_max",
            "AlphaMissense", "AlphaMissense_clase", "AlphaMissense_fuente", "VEP_hgvsp", "VEP_concuerda_p"]
    W = pd.DataFrame(wide).reindex(columns=cols).sort_values(["gen", "pos_NP", "variante"])
    return W, pd.DataFrame(long_rows), pd.DataFrame(sin_mapeo)
