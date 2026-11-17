import re

import numpy as np
import pandas as pd

IMPACTO = {"stop_gained": 4, "start_lost": 4, "stop_lost": 4, "splice_donor_variant": 4, "splice_acceptor_variant": 4,
           "frameshift_variant": 4, "missense_variant": 3, "inframe_insertion": 3, "inframe_deletion": 3,
           "protein_altering_variant": 3, "splice_region_variant": 2, "synonymous_variant": 2,
           "stop_retained_variant": 2, "start_retained_variant": 2}
IMPACTO_NOMBRE = {4: "HIGH", 3: "MODERATE", 2: "LOW", 1: "MODIFIER"}

UNIPROT_FENOTIPO = {
    "MPS1H": "MPS I Hurler (grave)", "MPS1H/S": "MPS I Hurler-Scheie (intermedia)", "MPS1S": "MPS I Scheie (atenuada)",
    "MPS2": "MPS II Hunter", "MPS3A": "MPS IIIA Sanfilippo A", "MPS3B": "MPS IIIB Sanfilippo B",
    "MPS3C": "MPS IIIC Sanfilippo C", "MPS3D": "MPS IIID Sanfilippo D", "MPS4A": "MPS IVA Morquio A",
    "MPS4B": "MPS IVB Morquio B", "MPS6": "MPS VI Maroteaux-Lamy", "MPS7": "MPS VII Sly", "MPS9": "MPS IX",
    "GM1G1": "Gangliosidosis GM1 tipo 1 (infantil)", "GM1G2": "Gangliosidosis GM1 tipo 2 (infantil tardía/juvenil)",
    "GM1G3": "Gangliosidosis GM1 tipo 3 (adulta)", "GM1": "Gangliosidosis GM1",
    "RP73": "Retinitis pigmentosa 73", "CMT2V": "Charcot-Marie-Tooth axonal tipo 2V",
}
FENOTIPO_NO_INFORMATIVO = {"not provided", "not specified", "inborn genetic diseases", "see cases", ""}

UMBRALES = {
    "SIFT": ("SIFT", lambda x: x <= 0.05),
    "PolyPhen2": ("PolyPhen2", lambda x: x >= 0.447),
    "CADD": ("CADD_phred", lambda x: x >= 20),
    "REVEL": ("REVEL", lambda x: x >= 0.644),
    "AlphaMissense": ("AlphaMissense", lambda x: x > 0.564),
}


def categoria_clinica(cl):
    cl = cl if isinstance(cl, str) else ""
    if not cl:
        return "sin clasificación"
    if "Conflicting" in cl or ("athogenic" in cl and "enign" in cl):
        return "conflicto"
    if "athogenic" in cl and "Uncertain" not in cl:
        return "P/LP"
    if "enign" in cl and "Uncertain" not in cl:
        return "B/LB"
    if "Uncertain" in cl and "athogenic" not in cl and "enign" not in cl:
        return "VUS"
    return "mixta"


def grupo_fenotipo(f):
    f = f.lower()
    if any(k in f for k in ["mucopolysacchar", "hurler", "scheie", "hunter", "sanfilippo", "morquio", "maroteaux",
                            "sly syndrome", "hyaluronoglucosaminidase"]):
        return "MPS"
    if "gangliosid" in f:
        return "Gangliosidosis GM1"
    if "charcot" in f:
        return "Charcot-Marie-Tooth 2V"
    if "retinitis" in f or "retinal" in f:
        return "Retinitis pigmentosa 73 / distrofia retiniana"
    if "related disorder" in f:
        return "Trastorno relacionado con el gen (sin especificar)"
    return "Otro (signo clínico o enfermedad lisosomal sin especificar)"


def consenso(row):
    votos = [(k, f(row[c])) for k, (c, f) in UMBRALES.items() if c in row and pd.notna(row[c])]
    n_dano, n = sum(v for _, v in votos), len(votos)
    if n == 0:
        c = "sin predicciones"
    elif n_dano >= 3 and n_dano >= n - 1:
        c = "dañino"
    elif n >= 3 and n_dano <= 1:
        c = "benigno"
    else:
        c = "discordante"
    return (f"{n_dano}/{n}" if n else None), (", ".join(k for k, v in votos if v) or None), c


def build(W, clinvar, vep):
    """W = Variantes_consolidadas; clinvar = ClinVar_missense (con columnas de longitud); vep = VEP_MANE."""
    cv = clinvar.set_index("VCV")
    vep_idx = vep.set_index(vep.variation_id.astype(str)) if len(vep) else pd.DataFrame()
    rows = []
    for _, w in W.iterrows():
        vcvs = [x.strip() for x in str(w.VCV).split("|")] if isinstance(w.VCV, str) else []
        recs = [cv.loc[v] for v in vcvs if v in cv.index]
        longs = [r.longitud_variante_pb for r in recs if pd.notna(r.longitud_variante_pb)]
        tipos = sorted({r.tipo_variante for r in recs if isinstance(r.tipo_variante, str)})
        csq = set()
        for r in recs:
            vid = str(r.variation_id)
            if vid in vep_idx.index and isinstance(vep_idx.loc[vid, "consecuencia"], str):
                csq.update(vep_idx.loc[vid, "consecuencia"].split(";"))
        imp = max([IMPACTO.get(t, 1) for t in csq], default=None)
        fen = sorted({f.strip() for r in recs if isinstance(r.condiciones, str) for f in r.condiciones.split("|")
                      if f.strip().lower() not in FENOTIPO_NO_INFORMATIVO})
        desc = w.UniProt_descripcion if isinstance(w.UniProt_descripcion, str) else ""
        codigos = [c for c in re.findall(r"\bin ([A-Z0-9/]+)", desc) if c in UNIPROT_FENOTIPO]
        gravedad = sorted(set(re.findall(r"(severe|intermediate|mild|attenuated)", desc.lower())))
        actividad = sorted(set(re.findall(r"(loss of activity|no effect|(?:reduced|decreased|decrease in)[^;]*activity)",
                                          desc.lower())))
        n_dano, cuales, cons = consenso(w)
        cat = categoria_clinica(w.clasificacion_ClinVar)
        if cat == "P/LP":
            res = "Asociada a enfermedad (evidencia clínica)"
        elif cat == "B/LB":
            res = "Sin efecto patogénico (evidencia clínica)"
        elif cons == "dañino":
            res = "Efecto potencialmente dañino (predicción, sin clasificación clínica concluyente)"
        elif cons == "benigno":
            res = "Efecto potencialmente benigno (predicción, sin clasificación clínica concluyente)"
        else:
            res = "Indeterminado (predicciones discordantes o ausentes)"
        rows.append(dict(
            gen=w.gen, variante=w.variante, HGVS_p=w.HGVS_p, pos_NP=w.pos_NP, pos_UniProt=w.pos_UniProt,
            VCV=w.VCV if isinstance(w.VCV, str) else None,
            longitud_variante_pb=max(longs) if longs else None, tipo_variante="; ".join(tipos) or None,
            consecuencia_VEP=";".join(sorted(csq)) or None, impacto_VEP=IMPACTO_NOMBRE.get(imp),
            clasificacion_ClinVar=w.clasificacion_ClinVar, categoria_clinica=cat, estrellas_max=w.estrellas_max,
            fenotipos_ClinVar=" | ".join(fen) or None,
            grupo_fenotipo_ClinVar="; ".join(sorted({grupo_fenotipo(f) for f in fen})) or "no especificado",
            fenotipo_UniProt="; ".join(UNIPROT_FENOTIPO[c] for c in dict.fromkeys(codigos)) or None,
            gravedad_UniProt="; ".join(gravedad) or None, efecto_actividad_UniProt="; ".join(actividad) or None,
            predictores_daninos=n_dano, predictores_que_indican_dano=cuales, consenso_in_silico=cons,
            potencial_efecto_fenotipico=res))
    return pd.DataFrame(rows)
