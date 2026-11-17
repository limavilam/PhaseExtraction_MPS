from . import clients as K
from . import config as C

L = C.LINKS
REFSEQ_PREFIX = {
    "NC": "cromosoma (genómico)", "NG": "RefSeqGene (región genómica)",
    "NM": "ARNm curado", "NR": "ARN no codificante curado",
    "XM": "ARNm predicho (modelo)", "XR": "ARN no codificante predicho",
    "NW": "contig/parche genómico (GRCh38 patch)", "NT": "contig genómico",
}


def identifiers(gene: str) -> dict:
    """GeneID, NC, NG, MANE Select (NM/NP/ENST/ENSP), UniProt, OMIM y enlaces."""
    K.esearch("gene", f"{gene}[sym] AND {C.TAXID}[taxid]", base="NCBI Gene")  
    reps = K.datasets(f"gene/symbol/{gene}/taxon/{C.TAXID}", f"símbolo {gene}")["reports"]
    match = [r["gene"] for r in reps if r["gene"]["symbol"] == gene]
    if len(match) != 1:
        raise ValueError(f"{gene}: {len(match)} coincidencias exactas entre {[r['gene']['symbol'] for r in reps]}")
    g = match[0]
    gid = g["gene_id"]
    pr = K.datasets(f"gene/id/{gid}/product_report",
                    f"GeneID {gid} → select_category=MANE_SELECT")["reports"][0]["product"]
    mane = next((t for t in pr.get("transcripts", []) if t.get("select_category") == "MANE_SELECT"), {})
    ann = g["annotations"][0]
    K.SOURCE_VERSIONS["NCBI RefSeq annotation"] = f"{ann.get('annotation_name')} ({ann.get('annotation_release_date')}), {ann.get('assembly_name')}"
    loc = ann["genomic_locations"][0]
    nc = loc["genomic_accession_version"]
    ng = [x["gene_range"]["accession_version"] for x in g.get("reference_standards", [])]
    nm = mane.get("accession_version"); npp = mane.get("protein", {}).get("accession_version")
    enst = mane.get("ensembl_transcript"); ensp = mane.get("protein", {}).get("ensembl_protein")
    ensg = (g.get("ensembl_gene_ids") or [None])[0]
    up = (g.get("swiss_prot_accessions") or [None])[0]
    omim = (g.get("omim_ids") or [None])[0]
    rng = loc["genomic_range"]
    return dict(
        gen=gene, MPS=C.PANEL.get(gene), descripcion=g.get("description"),
        GeneID=gid, GeneID_enlace=L["gene"].format(gid),
        cromosoma=g.get("chromosomes", [None])[0], NC=nc, NC_enlace=L["nuccore"].format(nc),
        inicio_GRCh38=rng["begin"], fin_GRCh38=rng["end"], hebra=rng.get("orientation"),
        NG=";".join(ng), NG_enlace=";".join(L["nuccore"].format(x) for x in ng),
        MANE_NM=nm, MANE_NM_enlace=L["nuccore"].format(nm) if nm else None,
        MANE_NP=npp, MANE_NP_enlace=L["protein"].format(npp) if npp else None,
        ENSG=ensg, ENSG_enlace=L["ensembl_gene"].format(ensg) if ensg else None,
        MANE_ENST=enst, MANE_ENST_enlace=L["ensembl_tx"].format(enst.split(".")[0]) if enst else None,
        MANE_ENSP=ensp,
        UniProt=up, UniProt_enlace=L["uniprot"].format(up) if up else None,
        AlphaFoldDB_enlace=L["alphafold"].format(up) if up else None,
        OMIM=omim, OMIM_enlace=L["omim"].format(omim) if omim else None,
        n_transcritos=g.get("transcript_count"), n_proteinas=g.get("protein_count"),
        fuente_API=C.DATASETS + f"gene/id/{gid}/product_report",
    )


def refseq_records(gene: str, chrom: str | None = None) -> list[dict]:
    """Todos los registros RefSeq del gen en Nucleotide, con prefijo y decisión de inclusión."""
    n_all, _ = K.esearch("nuccore", f'{gene}[Gene Name] AND "Homo sapiens"[Organism]',
                         base="NCBI Nucleotide")
    _, ids = K.esearch("nuccore", f'{gene}[Gene Name] AND "Homo sapiens"[Organism] AND refseq[filter]',
                       base="NCBI Nucleotide (RefSeq)")
    summ = K.esummary("nuccore", ids)
    rows = []
    for u, s in summ.items():
        acc = s["accessionversion"]; pre = acc[:2]
        tipo = REFSEQ_PREFIX.get(pre, "otro")
        if pre == "NC":
            tipo += " — GRCh38" if acc.startswith("NC_0000") else " — T2T-CHM13v2.0"
        title = s.get("title") or ""
        # La búsqueda [Gene Name] también devuelve genes que comparten alias (p. ej. ARSB -> ANKRD44, SLURP1).
        if pre in ("NC", "NW", "NT"):
            pertenece = chrom is None or f"chromosome {chrom}," in title or f"chromosome {chrom} " in title
        else:
            pertenece = f"({gene})" in title or f" {gene} " in title
        rows.append(dict(gen=gene, accesion=acc, prefijo=pre + "_", tipo=tipo,
                         pertenece_al_gen=pertenece,
                         motivo_exclusion=None if pertenece else "otro gen (alias compartido)", titulo=title,
                         longitud=s.get("slen"), enlace=L["nuccore"].format(acc),
                         total_nuccore_gen=n_all))
    return rows

def _comments(e, ctype):
    return [c for c in e.get("comments", []) if c["commentType"] == ctype]


def _eco(evs):
    codes = sorted({v["evidenceCode"] for v in evs})
    pmids = sorted({v["id"] for v in evs if v.get("source") == "PubMed"})
    tipo = sorted({C.ECO.get(c, c) for c in codes})
    experimental = any(c in ("ECO:0000269", "ECO:0007744") for c in codes)
    return ";".join(codes), "; ".join(tipo), experimental, ";".join(pmids)


def uniprot_protein(gene: str):
    """Devuelve (resumen, features, estructuras_PDB, secuencia) de la entrada Swiss-Prot."""
    hit = K.uniprot_search(f"gene_exact:{gene} AND organism_id:{C.TAXID} AND reviewed:true",
                           fields="accession", size=1)[0]
    acc = hit["primaryAccession"]
    e = K.uniprot_entry(acc)
    seq = e["sequence"]["value"]
    desc = e["proteinDescription"]
    rec = desc.get("recommendedName", {})
    txt = lambda cs: " | ".join(t["value"] for c in cs for t in c.get("texts", []))
    summary = dict(
        gen=gene, UniProt=acc, enlace=L["uniprot"].format(acc),
        entrada=e.get("uniProtkbId"), version_entrada=e.get("entryAudit", {}).get("entryVersion"),
        ultima_actualizacion=e.get("entryAudit", {}).get("lastAnnotationUpdateDate"),
        nombre=rec.get("fullName", {}).get("value"),
        EC=";".join(x["value"] for x in rec.get("ecNumbers", [])),
        longitud_aa=e["sequence"]["length"], masa_Da=e["sequence"]["molWeight"],
        funcion=txt(_comments(e, "FUNCTION")),
        actividad_catalitica=" | ".join(c["reaction"]["name"] for c in _comments(e, "CATALYTIC ACTIVITY")),
        cofactor=" | ".join(f["name"] for c in _comments(e, "COFACTOR") for f in c.get("cofactors", [])),
        localizacion=" | ".join(l["location"]["value"] for c in _comments(e, "SUBCELLULAR LOCATION")
                                for l in c.get("subcellularLocations", [])),
        enfermedades=" | ".join(c["disease"]["diseaseId"] for c in _comments(e, "DISEASE") if "disease" in c),
        isoformas=" | ".join(i.get("name", {}).get("value", "") for c in _comments(e, "ALTERNATIVE PRODUCTS")
                             for i in c.get("isoforms", [])),
        anotacion_score=e.get("annotationScore"),
        fuente_API=C.UNIPROT + f"{acc}.json",
    )
    feats = []
    for f in e.get("features", []):
        codes, tipo, exp, pmids = _eco(f.get("evidences", []))
        alt = f.get("alternativeSequence", {})
        xr = ";".join(f"{x['database']}:{x['id']}" for x in f.get("featureCrossReferences", []))
        feats.append(dict(
            gen=gene, UniProt=acc, tipo_feature=f["type"], id_feature=f.get("featureId"),
            inicio=f["location"]["start"].get("value"), fin=f["location"]["end"].get("value"),
            original=alt.get("originalSequence"), alternativa=";".join(alt.get("alternativeSequences", [])),
            descripcion=f.get("description"), refs_cruzadas=xr,
            ECO=codes, tipo_evidencia=tipo, evidencia_experimental=exp, PubMed=pmids,
        ))
    pdbs = []
    for x in e.get("uniProtKBCrossReferences", []):
        if x["database"] == "PDB":
            pr = {p["key"]: p["value"] for p in x.get("properties", [])}
            pdbs.append(dict(gen=gene, UniProt=acc, PDB=x["id"], metodo=pr.get("Method"),
                             resolucion=pr.get("Resolution"), cadenas=pr.get("Chains"),
                             enlace=L["pdb"].format(x["id"])))
    return summary, feats, pdbs, seq


def seq_concordance(gene, np_acc, up_acc, up_seq):
    np_seq = K.efetch_fasta("protein", np_acc)
    same = np_seq == up_seq
    offset, nota = 0, ""
    if not same and len(up_seq) > len(np_seq) and up_seq.endswith(np_seq):
        offset = len(up_seq) - len(np_seq)
        nota = f"UniProt tiene {offset} aa N-terminales extra: posición UniProt = posición NP + {offset}"
    elif not same and len(np_seq) > len(up_seq) and np_seq.endswith(up_seq):
        offset = -(len(np_seq) - len(up_seq))
        nota = f"NP tiene {-offset} aa N-terminales extra: posición UniProt = posición NP - {-offset}"
    elif not same and len(np_seq) == len(up_seq):
        diffs = [f"{a}{i+1}{b}" for i, (a, b) in enumerate(zip(up_seq, np_seq)) if a != b]
        nota = "sustituciones (residuo UniProt, posición, residuo NP): " + ";".join(diffs[:20])
    elif not same:
        nota = "longitud distinta sin desplazamiento simple: requiere alineamiento"
    return dict(gen=gene, MANE_NP=np_acc, UniProt=up_acc, long_NP=len(np_seq), long_UniProt=len(up_seq),
                identicas=same, desplazamiento_UniProt_menos_NP=offset, nota=nota)


def variant_length(spdi):
    """Longitud y coordenadas de una variante a partir de su SPDI (secuencia:posición desde 0:eliminado:insertado).
    longitud = número de bases eliminadas; inicio = posición SPDI + 1 (numeración desde 1); fin = inicio + longitud - 1.
    """
    try:
        _, pos, ref, alt = str(spdi).split(":")
        pos = int(pos)
    except ValueError:
        return dict(longitud_variante_pb=None, inicio_variante_GRCh38=None, fin_variante_GRCh38=None, tipo_variante=None)
    L = len(ref)
    tipo = ("SNV (1 nucleótido)" if L == 1 and len(alt) == 1 else
            "MNV (varios nucleótidos contiguos)" if L == len(alt) and L > 1 else "indel")
    return dict(longitud_variante_pb=L, inicio_variante_GRCh38=pos + 1,
                fin_variante_GRCh38=pos + max(L, 1), tipo_variante=tipo)


def clinvar_missense(gene: str) -> list[dict]:
    term = f'{gene}[gene] AND "missense variant"[molecular consequence]'
    _, ids = K.esearch("clinvar", term, base="ClinVar")
    rows = []
    for u, s in K.esummary("clinvar", ids).items():
        syms = [g["symbol"] for g in s.get("genes", [])]
        if gene not in syms:
            continue
        gc = s.get("germline_classification", {}) or {}
        vs = (s.get("variation_set") or [{}])[0]
        loc = next((l for l in vs.get("variation_loc", []) if l.get("assembly_name") == "GRCh38"), {})
        rows.append(dict(
            gen=gene, VCV=s.get("accession"), variation_id=u, enlace=L["clinvar"].format(u),
            nombre_HGVS=s.get("title"), cambio_cDNA=vs.get("cdna_change"),
            cambio_proteina=s.get("protein_change"),
            consecuencias=";".join(s.get("molecular_consequence_list", [])),
            clasificacion=gc.get("description"), estado_revision=gc.get("review_status"),
            ultima_evaluacion=gc.get("last_evaluated"),
            condiciones=" | ".join(t["trait_name"] for t in gc.get("trait_set", [])),
            genes_en_registro=";".join(syms), chr=loc.get("chr"), pos_GRCh38=loc.get("start"),
            rsID=";".join(x["db_id"] for x in vs.get("variation_xrefs", []) if x.get("db_source") == "dbSNP"),
            spdi=vs.get("canonical_spdi"),
            tipo_variante_ClinVar=vs.get("variant_type"),
            **variant_length(vs.get("canonical_spdi")),
        ))
    return rows
