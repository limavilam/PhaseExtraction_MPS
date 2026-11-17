# mps_extractor

Extracción reproducible de identificadores, evidencia a nivel de proteína (UniProtKB) y a nivel de
transcrito (NCBI RefSeq / ClinVar) para las 11 enzimas de las mucopolisacaridosis
(IDUA, IDS, SGSH, NAGLU, HGSNAT, GNS, GALNS, GLB1, ARSB, GUSB, HYAL1).

## Estructura

```
mps_extractor/
├── mps_extractor/
│   ├── config.py     # panel, endpoints, plantillas de enlaces, códigos ECO
│   ├── clients.py    # clientes NCBI E-utilities, NCBI Datasets v2, UniProt REST + log automático
│   ├── extract.py    # identificadores, registros RefSeq, proteína UniProt, ClinVar, concordancia
│   ├── evidence.py   # GenPept del NP_, submissions ClinVar (SCV), Ensembl VEP (+AlphaMissense)
│   ├── consolidate.py# tabla por variante (ancha) y tabla de evidencia por tipo/nivel (larga)
│   ├── phenotype.py  # hoja Efecto_fenotipico: longitud, impacto VEP, fenotipo (ClinVar/UniProt) y consenso in silico
│   ├── pipeline.py   # orquestación + Excel documentado (LEEME, Diccionario, hojas de datos, Registro)
│   └── api.py        # servicio FastAPI
└── run_extraccion.py # ejecución por línea de comandos
```

## Instalación y uso

```bash
pip install requests pandas openpyxl fastapi uvicorn
export NCBI_API_KEY=...      # opcional: sube el límite de 3 a 10 solicitudes/s
export NCBI_EMAIL=...        # recomendado por NCBI

python run_extraccion.py               # 11 genes -> salida/extraccion_MPS_<fecha>.xlsx + CSV (~15-30 min con VEP)
python run_extraccion.py IDUA GNS      # subconjunto

uvicorn mps_extractor.api:app --reload # servicio; documentación en http://127.0.0.1:8000/docs
```

## Endpoints

| Método | Ruta | Devuelve |
|---|---|---|
| GET | `/genes` | Panel y tipo de MPS |
| GET | `/genes/{gen}/identificadores` | GeneID, NC, NG, MANE NM/NP/ENST/ENSP, UniProt, OMIM + enlaces |
| GET | `/genes/{gen}/refseq` | Registros RefSeq de Nucleotide con tipo según prefijo |
| GET | `/genes/{gen}/proteina?tipo=Natural variant` | Resumen UniProt, features con ECO/PubMed, PDB, secuencia |
| GET | `/genes/{gen}/clinvar` | Registros ClinVar missense |
| GET | `/genes/{gen}/genpept` | Anotaciones GenPept del NP_ MANE |
| GET | `/genes/{gen}/clinvar/submissions` | Submissions (SCV): método, criterios ACMG citados, PubMed |
| GET | `/genes/{gen}/evidencia?formato=ancha\|larga\|efecto` | Tabla consolidada por variante o por elemento de evidencia (incluye VEP) |
| GET | `/registro` | Log de solicitudes API |
| GET | `/exportar?genes=IDUA&genes=GNS` | Excel documentado |

## Longitud y efecto fenotípico

- ClinVar_missense incluye `longitud_variante_pb`, `inicio_variante_GRCh38`, `fin_variante_GRCh38` y `tipo_variante`, calculados del SPDI (longitud = bases eliminadas; inicio = posición SPDI + 1), y `tipo_variante_ClinVar`.
- La hoja Efecto_fenotipico resume por variante la evidencia clínica, la forma clínica de UniProt y un consenso de cinco predictores; `potencial_efecto_fenotipico` es una regla de resumen, no una clasificación clínica.

## Notas

- Las respuestas de VEP se guardan en `salida/cache_vep/` (por contenido del lote): si la conexión se corta, al volver a ejecutar solo se piden los lotes faltantes.
- Numeración: `pos_UniProt = pos_NP + desplazamiento` (HGSNAT: +28). Ver hoja Concordancia_secuencia.
- `cita_PP3_BP4` marca submissions de ClinVar que usaron predictores in silico como evidencia (riesgo de circularidad al entrenar con esas etiquetas).

## Pendiente

- gnomAD v4.1 (frecuencias, homocigotos) como clase benigna de apoyo.
