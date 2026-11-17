import datetime as dt
import json
import pathlib
import shutil
import tempfile
import time

import pandas as pd

from . import clients as K
from . import config as C
from . import extract as X
from . import evidence as V
from . import consolidate as CO
from . import phenotype as PH

DICCIONARIO = {
    "Identificadores": ("Identificadores canónicos por gen y enlace web verificable de cada uno",
        "NCBI E-utilities (gene) + NCBI Datasets v2 (gene, product_report)", {
        "GeneID": "Identificador NCBI Gene", "NC": "Secuencia cromosómica RefSeq GRCh38",
        "NG": "RefSeqGene: región genómica de referencia para nomenclatura HGVS g./c. intrónica",
        "MANE_NM": "Transcrito MANE Select (RefSeq), referencia para HGVS c.",
        "MANE_NP": "Proteína del MANE Select (RefSeq), referencia para HGVS p.",
        "MANE_ENST": "Transcrito MANE Select equivalente en Ensembl", "MANE_ENSP": "Proteína Ensembl equivalente",
        "UniProt": "Accesión UniProtKB/Swiss-Prot revisada", "OMIM": "Número MIM del gen",
        "*_enlace": "Enlace web a la ficha del identificador en su base de origen",
        "fuente_API": "Solicitud API exacta de la que se obtuvo el registro"}),
    "Registros_RefSeq": ("Todos los registros RefSeq del gen en Nucleotide con su tipo según prefijo",
        "NCBI E-utilities esearch/esummary db=nuccore", {
        "prefijo": "NC_/NG_/NM_/NR_/XM_/XR_", "tipo": "Significado del prefijo",
        "seleccionado": "True si es el registro usado como referencia (MANE NM o NG)",
        "pertenece_al_gen": "False si el registro corresponde a otro gen que comparte alias con el símbolo buscado",
        "total_nuccore_gen": "Registros totales del gen en Nucleotide (RefSeq + GenBank)"}),
    "Proteina_UniProt": ("Resumen funcional de cada enzima", "UniProt REST /uniprotkb/{acc}.json", {
        "EC": "Número de la Enzyme Commission", "funcion": "Comentario FUNCTION",
        "actividad_catalitica": "Reacción(es) catalizada(s)", "localizacion": "Localización subcelular",
        "enfermedades": "Enfermedades asociadas (comentario DISEASE)",
        "anotacion_score": "Puntuación de anotación UniProt (1-5)"}),
    "Features_UniProt": ("Anotaciones por posición (sitios, dominios, PTM, variantes, mutagénesis)",
        "UniProt REST /uniprotkb/{acc}.json, campo features", {
        "tipo_feature": "Tipo de anotación UniProt (Active site, Binding site, Natural variant, Mutagenesis, ...)",
        "inicio/fin": "Posición en la secuencia canónica UniProt",
        "original/alternativa": "Residuo de referencia y sustituto (solo variantes y mutagénesis)",
        "ECO": "Códigos Evidence & Conclusion Ontology", "tipo_evidencia": "Traducción del código ECO",
        "evidencia_experimental": "True si algún código es ECO:0000269 o ECO:0007744",
        "PubMed": "Artículos que respaldan la anotación"}),
    "Variantes_UniProt": ("Subconjunto de Features: variantes naturales y mutagénesis (evidencia a nivel de proteína)",
        "Derivado de Features_UniProt", {}),
    "Estructuras_PDB": ("Estructuras experimentales vinculadas a cada entrada UniProt",
        "UniProt uniProtKBCrossReferences (PDB)", {"metodo": "X-ray / EM / NMR", "resolucion": "Å",
        "cadenas": "Cadenas y rango de residuos cubiertos"}),
    "Concordancia_secuencia": ("Comparación de la proteína MANE (NP_) con la secuencia canónica UniProt",
        "NCBI efetch db=protein + UniProt", {"identicas": "True si ambas secuencias son idénticas",
        "desplazamiento_UniProt_menos_NP": "Constante a sumar a la posición HGVS p. (NP) para obtener la posición UniProt/AlphaFold",
        "nota": "Descripción de la discrepancia"}),
    "ClinVar_missense": ("Registros ClinVar con consecuencia missense en el gen (evidencia clínica a nivel de transcrito)",
        "NCBI E-utilities esearch/esummary db=clinvar", {
        "VCV": "Accesión de variación ClinVar", "nombre_HGVS": "Nomenclatura HGVS sobre el transcrito",
        "cambio_cDNA": "HGVS c.", "cambio_proteina": "Cambio de aminoácido",
        "clasificacion": "Clasificación germinal agregada", "estado_revision": "Nivel de revisión (estrellas)",
        "condiciones": "Fenotipos asociados", "pos_GRCh38": "Posición genómica GRCh38",
        "spdi": "Variante en formato SPDI: secuencia:posición (desde 0):bases eliminadas:bases insertadas",
        "tipo_variante_ClinVar": "Tipo de variante según ClinVar (single nucleotide variant, Indel…)",
        "longitud_variante_pb": "Nucleótidos que cambian = número de bases eliminadas en el SPDI",
        "inicio_variante_GRCh38": "Primera base de la variante (numeración desde 1) = posición SPDI + 1",
        "fin_variante_GRCh38": "Última base de la variante = inicio + longitud − 1",
        "tipo_variante": "SNV (1 nucleótido), MNV (varios contiguos sustituidos por el mismo número) o indel"}),
    "GenPept_NP": ("Anotaciones del registro GenPept de la proteína MANE (regiones, sitios, péptidos)",
        "NCBI E-utilities efetch db=protein rettype=gp retmode=xml", {
        "tipo_feature": "Clave GenPept (Region, Site, sig_peptide, mat_peptide, Protein, CDS)",
        "nombre": "region_name / site_type / product", "db_xref": "Referencias (p. ej. CDD:...)",
        "inicio/fin": "Posición en el NP_ MANE (numeración HGVS p.)"}),
    "ClinVar_submissions": ("Cada submission (SCV) de los registros ClinVar missense",
        "NCBI E-utilities efetch db=clinvar rettype=vcv", {
        "SCV": "Accesión de la submission", "remitente": "Laboratorio o panel que la envió",
        "metodos": "MethodType de las observaciones (clinical testing, literature only, in vitro, in vivo, research, curation)",
        "metodo_funcional": "True si algún método es in vitro, in vivo o research",
        "criterios_ACMG": "Criterios ACMG/AMP citados en el comentario público (extraídos por expresión regular)",
        "cita_PS3_BS3": "True si cita evidencia funcional (PS3/BS3)",
        "cita_PP3_BP4": "True si cita predicción in silico (PP3/BP4): riesgo de circularidad",
        "comentario": "Comentario público (primeros 2000 caracteres)"}),
    "VEP_MANE": ("Anotación Ensembl VEP de cada registro ClinVar sobre el transcrito MANE",
        "Ensembl VEP REST POST /vep/human/region", {
        "SIFT": "0-1; ≤0.05 deletéreo", "PolyPhen2": "0-1; alto = dañino", "CADD_phred": "Escala phred",
        "REVEL": "0-1", "SpliceAI_max": "Máximo de DS_AG, DS_AL, DS_DG, DS_DL (0-1)",
        "AlphaMissense": "Patogenicidad 0-1 (nivel proteína)", "AlphaMissense_clase": "likely_benign / ambiguous / likely_pathogenic",
        "AlphaMissense_transcrito": "Transcrito del que proviene AlphaMissense (otro distinto al MANE solo si el MANE no tiene puntuación, p. ej. HGSNAT)"}),
    "Variantes_consolidadas": ("Una fila por sustitución proteica (numeración NP_ MANE) que une ClinVar, submissions, VEP y UniProt",
        "Derivado (consolidate.py)", {
        "variante": "Sustitución en código de una letra sobre el NP_ MANE",
        "pos_UniProt": "Posición en UniProt/AlphaFold = pos_NP + desplazamiento",
        "residuo_ref_coincide": "True si el residuo de referencia coincide con UniProt en pos_UniProt",
        "estrellas_max": "Máximo nivel de revisión ClinVar (0-4)",
        "SCV_cita_PP3_BP4": "Alguna submission usó predictores in silico como evidencia",
        "AlphaMissense_fuente": "Ensembl VEP (transcrito MANE) o AlphaFold DB (isoforma UniProt, por pos_UniProt)",
        "VEP_concuerda_p": "True si la sustitución de VEP coincide con la del título ClinVar"}),
    "Efecto_fenotipico": ("Una fila por variante: longitud, efecto molecular, evidencia clínica, fenotipo asociado y efecto potencial predicho",
        "Derivado de ClinVar_missense, VEP_MANE y Variantes_consolidadas (phenotype.py)", {
        "longitud_variante_pb / tipo_variante": "Máxima longitud y tipos entre los registros ClinVar de la variante",
        "impacto_VEP": "Categoría de impacto de Ensembl de la consecuencia más grave: HIGH, MODERATE (missense), LOW, MODIFIER",
        "categoria_clinica": "P/LP, B/LB, VUS, conflicto, mixta o sin clasificación",
        "fenotipos_ClinVar": "Enfermedades declaradas en ClinVar (sin 'not provided', 'not specified', 'Inborn genetic diseases', 'See cases')",
        "grupo_fenotipo_ClinVar": "MPS; gangliosidosis GM1; Charcot-Marie-Tooth 2V; retinitis pigmentosa 73 / distrofia retiniana; trastorno relacionado con el gen; otro; no especificado",
        "fenotipo_UniProt / gravedad_UniProt / efecto_actividad_UniProt": "Enfermedad, forma clínica (severe, intermediate, mild) y efecto sobre la actividad según UniProt",
        "predictores_daninos": "Predictores que indican daño / con valor. Umbrales: SIFT ≤0,05; PolyPhen-2 ≥0,447; CADD ≥20; REVEL ≥0,644; AlphaMissense >0,564",
        "consenso_in_silico": "dañino (≥3 y todos menos uno), benigno (≥3 con valor y ≤1 indica daño), discordante",
        "potencial_efecto_fenotipico": "Si ClinVar es P/LP o B/LB manda la evidencia clínica; si no, el consenso in silico como efecto potencial. No es una clasificación clínica"}),
    "Evidencia_por_variante": ("Tabla larga: una fila por elemento de evidencia, etiquetada por tipo y nivel",
        "Derivado (consolidate.py)", {
        "tipo_evidencia": "in silico / experimental / clínica", "nivel": "transcrito / proteína",
        "fuente": "Base y medida de origen", "valor": "Valor o descripción", "identificador": "VCV, SCV, HGVS c. o id UniProt"}),
    "Sin_mapeo": ("Registros ClinVar que no se pudieron mapear a una sustitución sobre el NP_ MANE", "Derivado", {}),
    "Registro_consultas": ("Log automático de cada solicitud API: fecha, base, consulta exacta, URL y conteo",
        "Generado por mps_extractor.clients", {}),
}


def run(genes=None, outdir="salida", clinvar=True, save_raw=True, vep=True, submissions=True):
    genes = genes or list(C.PANEL)
    out = pathlib.Path(outdir); out.mkdir(exist_ok=True)
    K.QUERY_LOG.clear(); K.SOURCE_VERSIONS.clear(); K.capture_versions()
    ids, refseq, prot, feats, pdbs, conc, cv = [], [], [], [], [], [], []
    gp, scv, seqs = [], [], {}
    stamp = dt.date.today().isoformat()
    gcache = out / "cache_genes"; gcache.mkdir(exist_ok=True)
    for g in genes:
        cf = gcache / f"{g}_{stamp}_{int(clinvar)}{int(submissions)}.json"
        if cf.exists():   
            d = json.loads(cf.read_text())
            ids.append(d["i"]); refseq += d["rs"]; prot.append(d["s"]); feats += d["f"]; pdbs += d["p"]
            seqs[g] = d["seq"]; conc.append(d["conc"]); gp += d["gp"]; cv += d["cv"]; scv += d["scv"]
            K.QUERY_LOG.extend(d.get("log", []))
            print(f"{g}: desde caché del día", flush=True); continue
        n_cv, n_scv, n_log = len(cv), len(scv), len(K.QUERY_LOG)
        i = X.identifiers(g); ids.append(i)
        rs = X.refseq_records(g, chrom=i["cromosoma"])
        for r in rs:
            r["seleccionado"] = r["accesion"] in (i["MANE_NM"], *i["NG"].split(";"))
        refseq += rs
        s, f, p, seq = X.uniprot_protein(g); prot.append(s); feats += f; pdbs += p; seqs[g] = seq
        conc.append(X.seq_concordance(g, i["MANE_NP"], s["UniProt"], seq))
        gp_g = V.genpept_features(g, i["MANE_NP"]); gp += gp_g
        if clinvar:
            rows_g = X.clinvar_missense(g)
            for r in rows_g:
                r["MANE_ENST"] = i["MANE_ENST"]
            cv += rows_g
            if submissions:
                scv += V.clinvar_submissions(g, [r["variation_id"] for r in rows_g])
        cf.write_text(json.dumps(dict(i=i, rs=rs, s=s, f=f, p=p, seq=seq, conc=conc[-1], gp=gp_g,
                                      cv=cv[n_cv:], scv=scv[n_scv:], log=K.QUERY_LOG[n_log:]), default=str))
        print(f"{g}: {len(f)} features UniProt, {len(p)} PDB, ClinVar {len(cv)}, SCV {len(scv)}", flush=True)
    vep_rows = []; V.VEP_ERRORS.clear()
    if clinvar and vep:
        V.vep_release()
        vep_rows = V.vep_annotate(cv, cache_dir=out / "cache_vep")
    feats_df = pd.DataFrame(feats)
    cv_df = pd.DataFrame(cv).drop(columns=["MANE_ENST"], errors="ignore")
    conc_df = pd.DataFrame(conc)
    am_up = {}
    if vep:
        for r in prot:
            try:
                am_up[r["gen"]] = V.alphamissense_uniprot(r["gen"], r["UniProt"])
            except Exception as e:
                print(f"  AlphaMissense AFDB {r['gen']}: {e}", flush=True)
    W, Lg, SM = CO.build(pd.DataFrame(ids), conc_df, seqs, cv_df, pd.DataFrame(scv), pd.DataFrame(vep_rows), feats_df, am_up)
    cv_df["mapeada_a_NP_MANE"] = ~cv_df["VCV"].isin(SM["VCV"]) if len(SM) else True
    sheets = {
        "Identificadores": pd.DataFrame(ids),
        "Registros_RefSeq": pd.DataFrame(refseq),
        "Proteina_UniProt": pd.DataFrame(prot),
        "Features_UniProt": feats_df,
        "Variantes_UniProt": feats_df[feats_df.tipo_feature.isin(["Natural variant", "Mutagenesis"])],
        "Estructuras_PDB": pd.DataFrame(pdbs),
        "Concordancia_secuencia": conc_df,
        "GenPept_NP": pd.DataFrame(gp),
        "ClinVar_missense": cv_df,
        "ClinVar_submissions": pd.DataFrame(scv),
        "VEP_MANE": pd.DataFrame(vep_rows),
        "Variantes_consolidadas": W,
        "Efecto_fenotipico": PH.build(W, cv_df, pd.DataFrame(vep_rows)),
        "Evidencia_por_variante": Lg,
        "Sin_mapeo": pd.concat([SM, pd.DataFrame(V.VEP_ERRORS).assign(motivo="rechazada por VEP")], ignore_index=True),
        "Registro_consultas": pd.DataFrame(K.QUERY_LOG),
    }
    stamp = dt.date.today().isoformat()
    if save_raw:  
        for k, df in sheets.items():
            df.to_csv(out / f"{k}.csv", index=False)
        (out / "versiones_fuentes.json").write_text(json.dumps(K.SOURCE_VERSIONS, indent=2))
    xlsx = out / f"extraccion_MPS_{stamp}.xlsx"
    xlsx = write_excel(sheets, xlsx, genes, stamp)
    return sheets, xlsx


def write_excel(sheets, path, genes, stamp):
    readme = pd.DataFrame([
        ("Proyecto", "Extracción de evidencia de proteína y transcrito para enzimas MPS"),
        ("Fecha de extracción", stamp),
        ("Genes", ", ".join(genes)),
        ("Organismo", f"Homo sapiens (taxid {C.TAXID})"),
        *[(k, v) for k, v in K.SOURCE_VERSIONS.items()],
        ("Solicitudes API registradas", len(K.QUERY_LOG)),
        ("Criterio de referencia", "Transcrito MANE Select (NM_/ENST) y su proteína (NP_/ENSP) ↔ UniProtKB/Swiss-Prot"),
        ("Exclusiones", "Registros XM_/XR_ (predichos) y NR_ (no codificantes); GenBank no RefSeq solo se cuenta"),
    ], columns=["campo", "valor"])
    dic_rows = []
    for hoja, (desc, fuente, cols) in DICCIONARIO.items():
        n = len(sheets[hoja]) if hoja in sheets else 0
        dic_rows.append(dict(hoja=hoja, columna="(hoja)", descripcion=desc, fuente=fuente, filas=n))
        for c, d in cols.items():
            dic_rows.append(dict(hoja=hoja, columna=c, descripcion=d, fuente=fuente, filas=None))

    path = pathlib.Path(path)
    tmpdir = pathlib.Path(tempfile.mkdtemp(prefix="mps_xlsx_"))
    tmp = tmpdir / path.name
    with pd.ExcelWriter(tmp, engine="openpyxl") as w:
        readme.to_excel(w, sheet_name="LEEME", index=False)
        pd.DataFrame(dic_rows).to_excel(w, sheet_name="Diccionario", index=False)
        for k, df in sheets.items():
            df.to_excel(w, sheet_name=k[:31], index=False)
        for ws in w.book.worksheets:         
            ws.freeze_panes = "A2"
            for col in ws.columns:
                width = max(len(str(c.value or "")) for c in list(col)[:200])
                ws.column_dimensions[col[0].column_letter].width = min(max(10, width + 2), 60)
            for row in ws.iter_rows(min_row=2):
                for c in row:
                    if isinstance(c.value, str) and c.value.startswith("https://") and ";" not in c.value:
                        c.hyperlink = c.value; c.style = "Hyperlink"
    for intento in range(1, 6):   
        try:
            shutil.copyfile(tmp, path)
            break
        except OSError as e:
            print(f"  copia a {path} falló (intento {intento}): {e}", flush=True)
            time.sleep(5 * intento)
    else:
        print(f"  No se pudo copiar a OneDrive; el Excel quedó en: {tmp}", flush=True)
        return tmp
    shutil.rmtree(tmpdir, ignore_errors=True)
    return path
