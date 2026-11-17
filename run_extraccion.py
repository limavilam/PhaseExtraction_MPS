""" python run_extraccion.py               (11 genes, con ClinVar)
       python run_extraccion.py IDUA GNS      (solo esos genes)
"""
import sys
from mps_extractor.pipeline import run

if __name__ == "__main__":
    sheets, xlsx = run(sys.argv[1:] or None, outdir="salida")
    print("Excel:", xlsx)
    for k, df in sheets.items():
        print(f"  {k}: {len(df)} filas")
