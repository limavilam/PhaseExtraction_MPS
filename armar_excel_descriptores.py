import sys, pathlib
SALIDA = pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else "salida")
import pandas as pd
import openpyxl
from openpyxl.styles import PatternFill, Font

src_sheets = {s: pd.read_csv(SALIDA / f"{s}.csv")
              for s in ["Descriptores_variante", "Descriptores_posicion", "Sitios_funcionales", "Modelos_AlphaFold"]}

xl = dict(src_sheets)

import numpy as np
from scipy.stats import mannwhitneyu

D = xl["Descriptores_variante"]
L = D[D.etiqueta_ClinVar.isin(["P/LP", "B/LB"]) & D.aa_modelo_coincide].copy()
L["y"] = (L.etiqueta_ClinVar == "P/LP").astype(int)
L2 = L[~L.SCV_cita_PP3_BP4.astype(bool)]

cols = ["RSA", "contactos_CA_12A", "contactos_CA_8A", "dist_sitio_catalitico", "dist_sitio_union",
        "dist_bolsillo_CDD", "dist_cys_disulfuro", "dist_glicosilacion", "pLDDT", "pLDDT_vecindad_10A",
        "Grantham", "BLOSUM62", "delta_volumen", "delta_hidrofobicidad", "delta_carga",
        "AlphaMissense", "REVEL", "CADD_phred"]

def auc(df, c):
    x = df[[c, "y"]].dropna()
    if x.y.nunique() < 2:
        return np.nan, len(x)
    u = mannwhitneyu(x[c][x.y == 1], x[c][x.y == 0]).statistic
    return u / (x.y.sum() * (len(x) - x.y.sum())), len(x)

res = []
for c in cols:
    a, n = auc(L, c)
    a2, n2 = auc(L2, c)
    x = L[[c, "y"]].dropna()
    res.append(dict(
        descriptor=c,
        tipo=("predictor genómico (comparador)" if c in ["AlphaMissense", "REVEL", "CADD_phred"]
              else "sustitución" if c in ["Grantham", "BLOSUM62", "delta_volumen", "delta_hidrofobicidad", "delta_carga"]
              else "posición (estructura AlphaFold)"),
        n=n, n_PLP=int(x.y.sum()), n_BLB=int((1 - x.y).sum()),
        mediana_PLP=round(x[c][x.y == 1].median(), 3), mediana_BLB=round(x[c][x.y == 0].median(), 3),
        AUC=round(a, 3), AUC_orientada=round(max(a, 1 - a), 3),
        direccion=("mayor en P/LP" if a >= 0.5 else "menor en P/LP"),
        AUC_orientada_sin_PP3_BP4=round(max(a2, 1 - a2), 3), n_sin_PP3_BP4=n2
    ))
R = pd.DataFrame(res).sort_values(["tipo", "AUC_orientada"], ascending=[True, False])
xl["Resumen_univariado"] = R

DL = "Descargado"
DT = "Descargado (predicción de terceros)"
LIT = "Literatura"
DER = "Derivado (organización)"
CP = "Cálculo propio"
AN = "Análisis propio (exploratorio)"

UP = "UniProtKB REST, release 2026_03"
AF = "AlphaFold DB, modelos v6"
CV = "ClinVar (E-utilities), consulta 29-10-2026"
VEP = "Ensembl VEP REST, release 116"

M = {
    "gen": ("Símbolo del gen del panel", DER, "Panel de 11 enzimas MPS definido en el proyecto"),
    "UniProt": ("Accesión UniProtKB/Swiss-Prot", DL, UP),
    "variante": ("Sustitución en código de una letra sobre el NP_ MANE", DER, "Reescrita desde el HGVS de ClinVar, VEP o UniProt"),
    "HGVS_p": ("Nomenclatura HGVS p. sobre el NP_ MANE", DER, "Reescrita desde ClinVar/VEP/UniProt"),
    "pos_NP": ("Posición en la proteína MANE (NP_)", DER, "Leída del HGVS de ClinVar o de VEP"),
    "pos_UniProt": ("Posición en UniProt y en el modelo AlphaFold (= pos_NP + desplazamiento)", DER, "Desplazamiento de la comparación NP_ vs UniProt (HGSNAT +28)"),
    "residuo_ref_coincide": ("El residuo de referencia coincide con UniProt en pos_UniProt", CP, "Control de calidad: comparación de secuencias"),
    "fuentes": ("Bases en las que aparece la variante", DER, "Registro de procedencia"),
    "VCV": ("Accesión(es) de variación ClinVar", DL, CV),
    "clasificacion_ClinVar": ("Clasificación germinal agregada", DL, CV),
    "etiqueta_ClinVar": ("Clasificación agrupada: P/LP, B/LB, VUS, conflicto, mixta, sin ClinVar", DER, "Agrupación de clasificacion_ClinVar"),
    "estrellas_max": ("Nivel de revisión máximo (0-4 estrellas)", DER, "Conversión del estado de revisión de ClinVar a estrellas (convención ClinVar)"),
    "SCV_cita_PP3_BP4": ("Alguna submission cita PP3/BP4 (predicción in silico) como criterio", CP, "Búsqueda por expresión regular en el comentario público de cada submission (ClinVar VCV XML)"),
    "SCV_metodo_funcional": ("Alguna submission tiene método in vitro, in vivo o research", DER, "Campo MethodType descargado de ClinVar"),
    "SCV_cita_PS3_BS3": ("Alguna submission cita PS3/BS3 (ensayo funcional) como criterio", CP, "Búsqueda por expresión regular en el comentario público"),
    "UniProt_experimental": ("La variante natural de UniProt tiene código ECO experimental", DER, "Códigos ECO descargados de UniProt (ECO:0000269 / ECO:0007744)"),
    "AlphaMissense": ("Patogenicidad AlphaMissense (0-1)", DT, VEP + "; AlphaFold DB (archivo de sustituciones) para HGSNAT y variantes solo en UniProt"),
    "REVEL": ("Puntuación REVEL (0-1)", DT, VEP),
    "CADD_phred": ("Puntuación CADD en escala phred", DT, VEP),
    "SIFT": ("Puntuación SIFT (0-1; ≤0,05 deletéreo)", DT, VEP),
    "PolyPhen2": ("Puntuación PolyPhen-2 (0-1)", DT, VEP),
    "SpliceAI_max": ("Máximo de las cuatro puntuaciones delta de SpliceAI", DT, VEP + "; el máximo de las 4 se toma en este proyecto"),
    "aa_ref": ("Aminoácido de referencia de la variante", DER, "Leído de la variante"),
    "aa_alt": ("Aminoácido alterno de la variante", DER, "Leído de la variante"),
    "aa_silvestre": ("Aminoácido del modelo en esa posición", DL, AF),
    "pLDDT": ("Confianza local del modelo (0-100)", DL, AF + " (columna B-factor del Cα)"),
    "pLDDT_vecindad_10A": ("Media del pLDDT de los residuos con Cα a ≤10 Å", CP, "Coordenadas y pLDDT de AlphaFold DB"),
    "confianza_modelo": ("Categoría de confianza: muy alta >90, alta >70, baja >50, muy baja", DER, "Umbrales de AlphaFold DB aplicados al pLDDT"),
    "RSA": ("Accesibilidad relativa al solvente", CP, "DSSP 4.6.1 sobre el modelo AlphaFold, escala de Tien et al. 2013"),
    "enterrado": ("RSA < 0,2", CP, "Umbral aplicado a RSA"),
    "estructura_secundaria_DSSP": ("Código DSSP de 8 estados", CP, "DSSP 4.6.1 sobre el modelo AlphaFold"),
    "estructura_secundaria": ("Hélice / lámina / lazo-giro", CP, "Agrupación del código DSSP"),
    "contactos_CA_8A": ("Residuos con Cα a ≤8 Å", CP, "Coordenadas AlphaFold (SciPy cKDTree)"),
    "contactos_CA_12A": ("Residuos con Cα a ≤12 Å", CP, "Coordenadas AlphaFold (SciPy cKDTree)"),
    "dist_sitio_catalitico": ("Distancia mínima entre átomos pesados al residuo catalítico más cercano (Å)", CP, "Coordenadas AlphaFold + sitio activo de UniProt (NAGLU: literatura)"),
    "es_sitio_catalitico": ("La posición es un residuo catalítico", DL, UP + "; NAGLU: literatura (PDB 4XWH)"),
    "dist_sitio_union": ("Distancia mínima a un sitio de unión (Å)", CP, "Coordenadas AlphaFold + sitios de unión de UniProt"),
    "es_sitio_union": ("La posición es un sitio de unión", DL, UP),
    "dist_bolsillo_CDD": ("Distancia mínima a un residuo del sitio activo del dominio conservado (Å)", CP, "Coordenadas AlphaFold + sitios 'active' del registro GenPept (CDD), solo sulfatasas"),
    "dist_glicosilacion": ("Distancia mínima a la Asn glicosilada más cercana (Å)", CP, "Coordenadas AlphaFold + glicosilaciones de UniProt"),
    "es_sitio_glicosilacion": ("La posición es una Asn glicosilada", DL, UP),
    "en_secuon_glicosilacion": ("La posición es N o N+2 de un secuón N-X-S/T glicosilado", CP, "Regla aplicada a las glicosilaciones de UniProt"),
    "dist_cys_disulfuro": ("Distancia mínima a una cisteína de puente disulfuro (Å)", CP, "Coordenadas AlphaFold + disulfuros de UniProt"),
    "es_cys_disulfuro": ("La posición es una cisteína de puente disulfuro", DL, UP),
    "en_peptido_senal": ("La posición está en el péptido señal", DL, UP),
    "en_propeptido": ("La posición está en un propéptido", DL, UP),
    "en_cadena_madura": ("La posición está en la cadena madura", DL, UP),
    "topologia": ("Lumenal / citoplasmático / transmembrana (solo HGSNAT)", DL, UP),
    "fuente_sitio_catalitico": ("Procedencia del sitio catalítico usado", DER, "Registro de procedencia"),
    "delta_volumen": ("Volumen alterno − silvestre (Å³)", CP, "Tabla de volúmenes de Zamyatnin 1972"),
    "delta_hidrofobicidad": ("Hidrofobicidad alterna − silvestre", CP, "Escala de Kyte y Doolittle 1982"),
    "delta_carga": ("Carga alterna − silvestre a pH 7 (D, E = −1; K, R = +1)", CP, "Regla de carga"),
    "BLOSUM62": ("Puntuación BLOSUM62 de la sustitución", CP, "Matriz BLOSUM62 (Henikoff y Henikoff 1992, Biopython)"),
    "Grantham": ("Distancia de Grantham de la sustitución", CP, "Fórmula y parámetros de Grantham 1974"),
    "introduce_prolina": ("El aminoácido alterno es prolina", CP, "Regla"),
    "pierde_prolina": ("El silvestre es prolina", CP, "Regla"),
    "introduce_glicina": ("El alterno es glicina", CP, "Regla"),
    "pierde_glicina": ("El silvestre es glicina", CP, "Regla"),
    "introduce_cisteina": ("El alterno es cisteína", CP, "Regla"),
    "pierde_cisteina": ("El silvestre es cisteína", CP, "Regla"),
    "aa_modelo_coincide": ("El aminoácido del modelo coincide con el de referencia", CP, "Control de calidad"),
}

SM = {
    "tipo": ("Tipo de sitio: catalítico, unión, bolsillo CDD, glicosilación, disulfuro", DER, "Clasificación de las anotaciones"),
    "pos_UniProt": ("Posición del sitio en numeración UniProt", DL, "UniProt; GenPept (CDD, convertido a numeración UniProt); NAGLU: literatura"),
    "fuente": ("Fuente del sitio", DER, "Registro de procedencia"),
}

MM = {
    "modelo": ("Identificador del modelo", DL, AF),
    "version": ("Versión del modelo", DL, AF),
    "pLDDT_global": ("pLDDT medio del modelo", DL, AF),
    "inicio": ("Primer residuo UniProt cubierto", DL, AF),
    "fin": ("Último residuo UniProt cubierto", DL, AF),
    "fecha_modelo": ("Fecha de creación del modelo", DL, AF),
    "url": ("URL de descarga del modelo", DL, AF),
    "archivo": ("Archivo local del modelo", DER, "Nombre del archivo descargado"),
    "longitud_modelo": ("Residuos contados en el modelo", CP, "Conteo sobre el archivo descargado"),
    "longitud_UniProt": ("Longitud de la secuencia UniProt", DL, UP),
}

RM = {
    "descriptor": ("Descriptor comparado", AN, "—"),
    "tipo": ("Grupo del descriptor", AN, "—"),
    "n": ("Variantes con valor", AN, "Variantes P/LP y B/LB"),
    "n_PLP": ("Variantes P/LP", AN, "—"),
    "n_BLB": ("Variantes B/LB", AN, "—"),
    "mediana_PLP": ("Mediana en P/LP", AN, "—"),
    "mediana_BLB": ("Mediana en B/LB", AN, "—"),
    "AUC": ("Área bajo la curva (U de Mann-Whitney normalizada)", AN, "SciPy; todas las enzimas juntas, sin particiones por gen"),
    "AUC_orientada": ("max(AUC, 1 − AUC)", AN, "—"),
    "direccion": ("Sentido del efecto", AN, "—"),
    "AUC_orientada_sin_PP3_BP4": ("AUC orientada excluyendo etiquetas que citaron PP3/BP4", AN, "—"),
    "n_sin_PP3_BP4": ("Variantes en ese subconjunto", AN, "—"),
}

rows = []
for sh, mp in [("Descriptores_variante", M), ("Descriptores_posicion", M),
               ("Sitios_funcionales", {**M, **SM}), ("Modelos_AlphaFold", {**M, **MM}),
               ("Resumen_univariado", RM)]:
    for c in xl[sh].columns:
        d = mp.get(c) or M.get(c)
        rows.append(dict(hoja=sh, columna=c, descripcion=d[0], origen=d[1], fuente_o_metodo=d[2]))
DIC = pd.DataFrame(rows)

defs = pd.DataFrame([
    ("Origen: Descargado", "Valor tal como lo publica la fuente (UniProt, ClinVar, AlphaFold DB); solo se copió"),
    ("Origen: Descargado (predicción de terceros)", "Puntuación calculada por otra herramienta (SIFT, PolyPhen-2, CADD, REVEL, SpliceAI, AlphaMissense) y entregada ya hecha por la API"),
    ("Origen: Derivado (organización)", "Dato descargado reescrito, agrupado o renumerado sin cálculo nuevo (p. ej. etiqueta agrupada, posición UniProt)"),
    ("Origen: Cálculo propio", "Calculado en este proyecto con software o fórmula sobre datos descargados (DSSP, distancias, contactos, escalas fisicoquímicas)"),
    ("Origen: Análisis propio (exploratorio)", "Estadística exploratoria de este proyecto; no es un dato ni un resultado de evaluación"),
    ("Excepción de literatura", "Sitio catalítico de NAGLU (Glu316, Glu446) transcrito de la estructura cristalina PDB 4XWH; señalado en fuente_o_metodo"),
], columns=["campo", "valor"])

leeme_base = pd.DataFrame([
    ("Proyecto", "Descriptores estructurales de variantes missense en enzimas MPS"),
    ("Estructuras", "AlphaFold DB, modelos v6 (monómero, apo), numeración UniProt"),
    ("Posición de la variante", "pos_UniProt = posición HGVS p. sobre el NP_ MANE + desplazamiento (HGSNAT +28)"),
    ("Principio", "Se describe la posición silvestre y el cambio fisicoquímico; no se construyen modelos mutantes"),
    ("Software", "DSSP 4.6.1 (mkdssp), Biopython 1.88, SciPy 1.18"),
    ("Sitio catalítico NAGLU", "Glu316 y Glu446, de la estructura cristalina PDB 4XWH (no anotados en UniProt)"),
    ("Excluidas", "4 sustituciones a codón de parada etiquetadas como missense en ClinVar (GLB1 Y154*, Q175*, Q184*; GUSB W111*)"),
    ("Limitaciones", "Modelos monoméricos (GALNS y GLB1 son dímeros, GUSB tetrámero: las interfaces no se describen); sin glicanos (RSA sobreestimada cerca de sitios glicosilados); sin modificación a formilglicina ni Ca2+ en sulfatasas; regiones con pLDDT < 70 (p. ej. péptido señal) tienen geometría poco fiable"),
], columns=["campo", "valor"])

LE = pd.concat([leeme_base, defs], ignore_index=True)

fills = {DL: "DDEBF7", DT: "E2EFDA", DER: "F2F2F2", CP: "FCE4D6", AN: "FFF2CC"}

out = SALIDA / "descriptores_estructurales_MPS.xlsx"
# Escribir primero en disco local (OneDrive puede bloquear/expirar escrituras directas) y luego copiar
import tempfile, shutil
tmp = pathlib.Path(tempfile.mkdtemp()) / out.name
with pd.ExcelWriter(tmp, engine="openpyxl") as w:
    LE.to_excel(w, sheet_name="LEEME", index=False)
    DIC.to_excel(w, sheet_name="Diccionario", index=False)
    for s in ["Descriptores_variante", "Descriptores_posicion", "Sitios_funcionales", "Modelos_AlphaFold", "Resumen_univariado"]:
        xl[s].to_excel(w, sheet_name=s, index=False)
    for ws in w.book.worksheets:
        ws.freeze_panes = "A2"
        for c in ws[1]:
            c.font = Font(bold=True)
        for colc in ws.columns:
            ws.column_dimensions[colc[0].column_letter].width = min(
                max(10, max(len(str(c.value or "")) for c in list(colc)[:120]) + 2), 70)
    ws = w.book["Diccionario"]
    for row in ws.iter_rows(min_row=2):
        f = fills.get(row[3].value)
        if f:
            row[3].fill = PatternFill("solid", fgColor=f)
    ws = w.book["LEEME"]
    for row in ws.iter_rows(min_row=2):
        for k, f in fills.items():
            if row[0].value == f"Origen: {k}":
                row[0].fill = PatternFill("solid", fgColor=f)
try:
    shutil.copyfile(tmp, out)
    print(f"Escrito {out}")
except OSError as e:
    print(f"No se pudo copiar a {out} ({e}). ¿Está abierto en Excel o OneDrive sincronizando?\nCopia local en: {tmp}")
