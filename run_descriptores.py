import pathlib, runpy, sys
from mps_extractor.structure import run

if __name__ == "__main__":
    salida = sys.argv[1] if len(sys.argv) > 1 else "salida"
    P, D, S, M = run(salida)
    print(f"Descriptores_posicion: {len(P)} residuos | Descriptores_variante: {len(D)} variantes | Sitios: {len(S)}")
    sys.argv = ["armar_excel_descriptores.py", salida]
    runpy.run_path(str(pathlib.Path(__file__).with_name("armar_excel_descriptores.py")), run_name="__main__")
