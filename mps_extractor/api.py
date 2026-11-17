import pathlib
import tempfile

from fastapi import FastAPI, HTTPException, Query
from fastapi.responses import FileResponse

from . import clients as K
from . import config as C
from . import extract as X
from . import pipeline as P
from . import evidence as V
from . import consolidate as CO
from . import phenotype as PH
import pandas as pd

app = FastAPI(title="MPS extractor",
              description="Extracción de identificadores, evidencia de proteína (UniProt) "
                          "y de transcrito (NCBI/ClinVar) para las enzimas de las mucopolisacaridosis.",
              version="0.1.0")


def _check(gene):
    gene = gene.upper()
    if gene not in C.PANEL:
        raise HTTPException(404, f"{gene} no está en el panel: {', '.join(C.PANEL)}")
    return gene


@app.get("/genes")
def genes():
    return C.PANEL


@app.get("/genes/{gene}/identificadores")
def identificadores(gene: str):
    return X.identifiers(_check(gene))


@app.get("/genes/{gene}/refseq")
def refseq(gene: str):
    return X.refseq_records(_check(gene))


@app.get("/genes/{gene}/proteina")
def proteina(gene: str, tipo: list[str] | None = Query(None, description="Filtrar features, p. ej. 'Natural variant'")):
    s, f, p, seq = X.uniprot_protein(_check(gene))
    if tipo:
        f = [x for x in f if x["tipo_feature"] in tipo]
    return dict(resumen=s, features=f, pdb=p, secuencia=seq)


@app.get("/genes/{gene}/clinvar")
def clinvar(gene: str):
    return X.clinvar_missense(_check(gene))


@app.get("/genes/{gene}/genpept")
def genpept(gene: str):
    g = _check(gene)
    return V.genpept_features(g, X.identifiers(g)["MANE_NP"])


@app.get("/genes/{gene}/clinvar/submissions")
def clinvar_scv(gene: str):
    g = _check(gene)
    return V.clinvar_submissions(g, [r["variation_id"] for r in X.clinvar_missense(g)])


@app.get("/genes/{gene}/evidencia")
def evidencia(gene: str, formato: str = Query("ancha", pattern="^(ancha|larga|efecto)$")):
    g = _check(gene)
    i = X.identifiers(g); s, f, _, seq = X.uniprot_protein(g)
    conc = pd.DataFrame([X.seq_concordance(g, i["MANE_NP"], s["UniProt"], seq)])
    cv = X.clinvar_missense(g)
    for r in cv: r["MANE_ENST"] = i["MANE_ENST"]
    scv = V.clinvar_submissions(g, [r["variation_id"] for r in cv])
    vep = V.vep_annotate(cv)
    W, L, _ = CO.build(pd.DataFrame([i]), conc, {g: seq}, pd.DataFrame(cv), pd.DataFrame(scv),
                       pd.DataFrame(vep), pd.DataFrame(f))
    if formato == "efecto":
        df = PH.build(W, pd.DataFrame(cv), pd.DataFrame(vep))
    else:
        df = W if formato == "ancha" else L
    return df.astype(object).where(df.notna(), None).to_dict("records")


@app.get("/registro")
def registro():
    return K.QUERY_LOG


@app.get("/exportar")
def exportar(genes: list[str] = Query(None), clinvar: bool = True, vep: bool = True):
    """Corre la extracción completa y devuelve el Excel documentado."""
    gl = [_check(g) for g in genes] if genes else None
    tmp = pathlib.Path(tempfile.mkdtemp())
    _, xlsx = P.run(gl, outdir=tmp, clinvar=clinvar, save_raw=False, vep=vep)
    return FileResponse(xlsx, filename=xlsx.name,
                        media_type="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
