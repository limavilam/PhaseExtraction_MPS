"""Descriptores estructurales sobre los modelos de AlphaFold DB (estructura silvestre, numeración UniProt).

Se describe la POSICIÓN silvestre (no se construyen modelos mutantes) y el CAMBIO fisicoquímico de la
sustitución. La posición de cada variante en el modelo es `pos_UniProt` (Variantes_consolidadas).

"""
import os
import pathlib
import re
import shutil

import numpy as np
import pandas as pd
import requests
from scipy.spatial import cKDTree

AF_API = "https://alphafold.ebi.ac.uk/api/prediction/{}"

# Sitios catalíticos no anotados en UniProt ni GenPept, tomados de la literatura (numeración UniProt)
SITIOS_MANUALES = {
    "NAGLU": [(316, "Glu catalítico (ácido/base)", "PDB 4XWH; doi:10.3390/ijms27104404"),
              (446, "Glu catalítico (nucleófilo/base)", "PDB 4XWH; doi:10.3390/ijms27104404")],
}

KD = dict(A=1.8, R=-4.5, N=-3.5, D=-3.5, C=2.5, Q=-3.5, E=-3.5, G=-0.4, H=-3.2, I=4.5, L=3.8, K=-3.9,
          M=1.9, F=2.8, P=-1.6, S=-0.8, T=-0.7, W=-0.9, Y=-1.3, V=4.2)            # Kyte y Doolittle 1982
VOL = dict(A=88.6, R=173.4, N=114.1, D=111.1, C=108.5, Q=143.8, E=138.4, G=60.1, H=153.2, I=166.7,
           L=166.7, K=168.6, M=162.9, F=189.9, P=112.7, S=89.0, T=116.1, W=227.8, Y=193.6, V=140.0)  # Zamyatnin 1972, Å3
CHG = dict(D=-1, E=-1, K=1, R=1)                                                     # pH 7 (His = 0)
GR = dict(S=(1.42, 9.2, 32), R=(0.65, 10.5, 124), L=(0, 4.9, 111), P=(0.39, 8.0, 32.5), T=(0.71, 8.6, 61),
          A=(0, 8.1, 31), V=(0, 5.9, 84), G=(0.74, 9.0, 3), I=(0, 5.2, 111), F=(0, 5.2, 132), Y=(0.20, 6.2, 136),
          C=(2.75, 5.5, 55), H=(0.58, 10.4, 96), Q=(0.89, 10.5, 85), N=(1.33, 11.6, 56), K=(0.33, 11.3, 119),
          D=(1.38, 13.0, 54), E=(0.92, 12.3, 83), M=(0, 5.7, 105), W=(0.13, 5.4, 170))   # Grantham 1974
SS3 = {"H": "hélice", "G": "hélice", "I": "hélice", "P": "hélice", "E": "lámina", "B": "lámina"}


def _dssp_env():
    exe = shutil.which("mkdssp")
    if exe is None:
        raise RuntimeError("No se encontró 'mkdssp' en el PATH. Instálalo con: "
                           "conda install -c conda-forge dssp")
    if "LIBCIFPP_DATA_DIR" not in os.environ:
        d = pathlib.Path(exe).resolve().parent.parent / "share" / "libcifpp"
        if d.exists():
            os.environ["LIBCIFPP_DATA_DIR"] = str(d)
    return exe


def fetch_model(acc, outdir):
    """Descarga el modelo AlphaFold DB (PDB) de una accesión UniProt; devuelve (ruta, metadatos)."""
    outdir = pathlib.Path(outdir); outdir.mkdir(parents=True, exist_ok=True)
    m = requests.get(AF_API.format(acc), timeout=60).json()[0]
    path = outdir / m["pdbUrl"].split("/")[-1]
    if not path.exists():
        path.write_bytes(requests.get(m["pdbUrl"], timeout=300).content)
    meta = dict(UniProt=acc, modelo=m["entryId"], version=m.get("latestVersion"),
                pLDDT_global=m.get("globalMetricValue"), inicio=m.get("uniprotStart"), fin=m.get("uniprotEnd"),
                fecha_modelo=m.get("modelCreatedDate"), url=m["pdbUrl"], archivo=path.name)
    return path, meta


def _ranges(df, tipo):
    x = df[df.tipo_feature == tipo]
    return [(int(a), int(b), str(d)) for a, b, d in zip(x.inicio, x.fin, x.descripcion)
            if not (pd.isna(a) or pd.isna(b))]


def _cdd_active(gp_gene, offset):
    """Residuos de 'Site active' del registro GenPept (CDD, numeración NP) -> numeración UniProt."""
    out = set()
    for loc in gp_gene[(gp_gene.tipo_feature == "Site") & (gp_gene.nombre == "active")].ubicacion.dropna():
        for a, b in re.findall(r"(\d+)(?:\.\.(\d+))?", loc):
            a = int(a); b = int(b) if b else a
            if b - a < 3:                       # puntos o intervalos cortos (no el rango global)
                out.update(range(a + offset, b + offset + 1))
    return sorted(out)


def position_descriptors(gene, acc, pdb_path, feats, genpept, offset):
    """Una fila por residuo del modelo silvestre."""
    from Bio.PDB import PDBParser
    from Bio.PDB.DSSP import DSSP
    from Bio.PDB.Polypeptide import three_to_index, index_to_one
    exe = _dssp_env()
    model = PDBParser(QUIET=True).get_structure(gene, str(pdb_path))[0]
    chain = next(iter(model))
    residues = [r for r in chain if r.id[0] == " "]
    pos = np.array([r.id[1] for r in residues])
    aa = [index_to_one(three_to_index(r.get_resname())) for r in residues]
    ca = np.array([r["CA"].coord for r in residues])
    plddt = np.array([r["CA"].get_bfactor() for r in residues])
    # átomos pesados con índice de residuo
    at_xyz, at_res = [], []
    for i, r in enumerate(residues):
        for a in r:
            if a.element != "H":
                at_xyz.append(a.coord); at_res.append(i)
    at_xyz = np.array(at_xyz); at_res = np.array(at_res)

    dssp = DSSP(model, str(pdb_path), dssp=exe, acc_array="Wilke")
    ss = {k[1][1]: dssp[k][2] for k in dssp.keys()}
    rsa = {k[1][1]: dssp[k][3] for k in dssp.keys()}

    tree_ca = cKDTree(ca)
    n8 = np.array([len(x) - 1 for x in tree_ca.query_ball_point(ca, 8.0)])
    n12 = np.array([len(x) - 1 for x in tree_ca.query_ball_point(ca, 12.0)])
    vec = tree_ca.query_ball_point(ca, 10.0)
    plddt_vec = np.array([plddt[v].mean() for v in vec])

    f = feats[feats.gen == gene]
    idx = {p: i for i, p in enumerate(pos)}

    def min_dist(site_pos):
        site_pos = [p for p in site_pos if p in idx]
        if not site_pos:
            return np.full(len(residues), np.nan)
        mask = np.isin(at_res, [idx[p] for p in site_pos])
        d, _ = cKDTree(at_xyz[mask]).query(at_xyz)
        out = np.full(len(residues), np.inf)
        np.minimum.at(out, at_res, d)
        return out

    cat = [int(p) for p in f[f.tipo_feature == "Active site"].inicio.dropna()]
    cat_src = "UniProt"
    if not cat and gene in SITIOS_MANUALES:
        cat = [p for p, _, _ in SITIOS_MANUALES[gene]]; cat_src = "literatura (" + SITIOS_MANUALES[gene][0][2] + ")"
    bind = sorted({int(p) for a, b, _ in _ranges(f, "Binding site") for p in range(a, b + 1)})
    cdd = _cdd_active(genpept[genpept.gen == gene], offset)
    glyc = [int(p) for p in f[f.tipo_feature == "Glycosylation"].inicio.dropna()]
    ss_pairs = [(int(a), int(b)) for a, b, _ in _ranges(f, "Disulfide bond")]
    ss_cys = sorted({p for ab in ss_pairs for p in ab})
    sig = _ranges(f, "Signal"); prop = _ranges(f, "Propeptide"); chains = _ranges(f, "Chain")
    tm = _ranges(f, "Transmembrane"); topo = _ranges(f, "Topological domain")

    def inr(p, rr):
        return any(a <= p <= b for a, b, _ in rr)

    def topo_of(p):
        if inr(p, tm):
            return "transmembrana"
        for a, b, d in topo:
            if a <= p <= b:
                return d.lower()
        return None

    d_cat, d_bind, d_cdd, d_gly, d_ss = (min_dist(cat), min_dist(bind), min_dist(cdd), min_dist(glyc), min_dist(ss_cys))
    sequon = set(glyc) | {g + 2 for g in glyc}
    rows = []
    for i, p in enumerate(pos):
        s = ss.get(p, "-")
        rows.append(dict(
            gen=gene, UniProt=acc, pos_UniProt=int(p), aa_silvestre=aa[i],
            pLDDT=round(float(plddt[i]), 2), pLDDT_vecindad_10A=round(float(plddt_vec[i]), 2),
            confianza_modelo=("muy alta" if plddt[i] > 90 else "alta" if plddt[i] > 70 else "baja" if plddt[i] > 50 else "muy baja"),
            RSA=round(float(rsa.get(p, np.nan)), 3) if rsa.get(p) not in (None, "NA") else np.nan,
            enterrado=bool(rsa.get(p, 1) < 0.2) if rsa.get(p) not in (None, "NA") else None,
            estructura_secundaria_DSSP=s, estructura_secundaria=SS3.get(s, "lazo/giro"),
            contactos_CA_8A=int(n8[i]), contactos_CA_12A=int(n12[i]),
            dist_sitio_catalitico=round(float(d_cat[i]), 2), es_sitio_catalitico=p in cat,
            dist_sitio_union=round(float(d_bind[i]), 2), es_sitio_union=p in bind,
            dist_bolsillo_CDD=round(float(d_cdd[i]), 2),
            dist_glicosilacion=round(float(d_gly[i]), 2), es_sitio_glicosilacion=p in glyc,
            en_secuon_glicosilacion=p in sequon,
            dist_cys_disulfuro=round(float(d_ss[i]), 2), es_cys_disulfuro=p in ss_cys,
            en_peptido_senal=inr(p, sig), en_propeptido=inr(p, prop), en_cadena_madura=inr(p, chains),
            topologia=topo_of(p), fuente_sitio_catalitico=cat_src,
        ))
    sites = ([dict(gen=gene, tipo="catalítico", pos_UniProt=p, fuente=cat_src) for p in cat] +
             [dict(gen=gene, tipo="unión", pos_UniProt=p, fuente="UniProt") for p in bind] +
             [dict(gen=gene, tipo="bolsillo activo CDD", pos_UniProt=p, fuente="NCBI GenPept (CDD)") for p in cdd] +
             [dict(gen=gene, tipo="glicosilación", pos_UniProt=p, fuente="UniProt") for p in glyc] +
             [dict(gen=gene, tipo="cisteína en disulfuro", pos_UniProt=p, fuente="UniProt") for p in ss_cys])
    return rows, sites


def substitution_descriptors(ref, alt):
    from Bio.Align import substitution_matrices
    if not hasattr(substitution_descriptors, "B"):
        substitution_descriptors.B = substitution_matrices.load("BLOSUM62")
    B = substitution_descriptors.B
    if ref not in KD or alt not in KD:
        return {}
    (c1, p1, v1), (c2, p2, v2) = GR[ref], GR[alt]
    return dict(
        delta_volumen=round(VOL[alt] - VOL[ref], 1), delta_hidrofobicidad=round(KD[alt] - KD[ref], 2),
        delta_carga=CHG.get(alt, 0) - CHG.get(ref, 0), BLOSUM62=int(B[ref][alt]),
        Grantham=round(50.723 * np.sqrt(1.833 * (c1 - c2) ** 2 + 0.1018 * (p1 - p2) ** 2 + 0.000399 * (v1 - v2) ** 2), 1),
        introduce_prolina=alt == "P", pierde_prolina=ref == "P", introduce_glicina=alt == "G",
        pierde_glicina=ref == "G", introduce_cisteina=alt == "C", pierde_cisteina=ref == "C",
    )


def etiqueta(cl):
    cl = str(cl) if isinstance(cl, str) else ""
    if not cl:
        return "sin ClinVar"
    if "Conflicting" in cl or ("athogenic" in cl and "enign" in cl):
        return "conflicto"
    if "athogenic" in cl and "Uncertain" not in cl:
        return "P/LP"
    if "enign" in cl and "Uncertain" not in cl:
        return "B/LB"
    if "Uncertain" in cl and "athogenic" not in cl and "enign" not in cl:
        return "VUS"
    return "mixta"


def run(salida="salida", out_xlsx=None):
    salida = pathlib.Path(salida)
    W = pd.read_csv(salida / "Variantes_consolidadas.csv")
    F = pd.read_csv(salida / "Features_UniProt.csv"); G = pd.read_csv(salida / "GenPept_NP.csv")
    conc = pd.read_csv(salida / "Concordancia_secuencia.csv").set_index("gen")
    prot = pd.read_csv(salida / "Proteina_UniProt.csv")
    pos_rows, site_rows, metas = [], [], []
    for _, r in prot.iterrows():
        g, acc = r.gen, r.UniProt
        path, meta = fetch_model(acc, salida / "modelos_af")
        meta["gen"] = g; metas.append(meta)
        rows, sites = position_descriptors(g, acc, path, F, G, int(conc.loc[g, "desplazamiento_UniProt_menos_NP"]))
        seq_model = "".join(x["aa_silvestre"] for x in rows)
        meta["longitud_modelo"] = len(rows); meta["longitud_UniProt"] = int(r.longitud_aa)
        pos_rows += rows; site_rows += sites
        print(f"{g}: {len(rows)} residuos, pLDDT global {meta['pLDDT_global']}", flush=True)
    P = pd.DataFrame(pos_rows)
    keep = ["gen", "variante", "HGVS_p", "pos_NP", "pos_UniProt", "residuo_ref_coincide", "fuentes", "VCV",
            "clasificacion_ClinVar", "estrellas_max", "SCV_cita_PP3_BP4", "SCV_metodo_funcional", "SCV_cita_PS3_BS3",
            "UniProt_experimental", "AlphaMissense", "REVEL", "CADD_phred", "SIFT", "PolyPhen2", "SpliceAI_max"]
    V = W[keep].copy()
    V.insert(9, "etiqueta_ClinVar", V.clasificacion_ClinVar.map(etiqueta))
    V["aa_ref"] = V.variante.str[0]; V["aa_alt"] = V.variante.str[-1]
    no_missense = V[~V.aa_alt.isin(list(KD))]           # p. ej. codones de parada etiquetados como missense
    if len(no_missense):
        print(f"  excluidas {len(no_missense)} sustituciones no missense: {', '.join(no_missense.gen + ' ' + no_missense.variante)}")
    V = V[V.aa_alt.isin(list(KD))].reset_index(drop=True)
    S = pd.DataFrame([substitution_descriptors(a, b) for a, b in zip(V.aa_ref, V.aa_alt)], index=V.index)
    D = V.merge(P.drop(columns=["UniProt"]), on=["gen", "pos_UniProt"], how="left")
    D = pd.concat([D, S.set_index(D.index)], axis=1)
    D["aa_modelo_coincide"] = D.aa_silvestre == D.aa_ref
    M = pd.DataFrame(metas)
    P.to_csv(salida / "Descriptores_posicion.csv", index=False)
    D.to_csv(salida / "Descriptores_variante.csv", index=False)
    pd.DataFrame(site_rows).to_csv(salida / "Sitios_funcionales.csv", index=False)
    M.to_csv(salida / "Modelos_AlphaFold.csv", index=False)
    return P, D, pd.DataFrame(site_rows), M
