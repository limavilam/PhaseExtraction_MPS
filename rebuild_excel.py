import json
import pathlib
import sys

import pandas as pd

from mps_extractor import clients as K
from mps_extractor.pipeline import write_excel

out = pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else "salida")
orden = ["Identificadores", "Registros_RefSeq", "Proteina_UniProt", "Features_UniProt", "Variantes_UniProt",
         "Estructuras_PDB", "Concordancia_secuencia", "GenPept_NP", "ClinVar_missense", "ClinVar_submissions",
         "VEP_MANE", "Variantes_consolidadas", "Efecto_fenotipico", "Evidencia_por_variante", "Sin_mapeo",
         "Registro_consultas"]   # mismo orden que pipeline.run
nombres = [n for n in orden if (out / f"{n}.csv").exists()]
sheets = {n: pd.read_csv(out / f"{n}.csv", low_memory=False) for n in nombres}
K.SOURCE_VERSIONS.update(json.loads((out / "versiones_fuentes.json").read_text()))
K.QUERY_LOG.extend(sheets.get("Registro_consultas", pd.DataFrame()).to_dict("records"))
genes = sheets["Identificadores"]["gen"].tolist()
stamp = K.SOURCE_VERSIONS.get("ClinVar consultado", "")
xlsx = write_excel(sheets, out / f"extraccion_MPS_{stamp}.xlsx", genes, stamp)
print("Excel:", xlsx)
for k, df in sheets.items():
    print(f"  {k}: {len(df)} filas")
