import os

# Panel de 11 enzimas asociadas a mucopolisacaridosis (gen -> tipo de MPS)
PANEL = {
    "IDUA": "MPS I", "IDS": "MPS II", "SGSH": "MPS IIIA", "NAGLU": "MPS IIIB",
    "HGSNAT": "MPS IIIC", "GNS": "MPS IIID", "GALNS": "MPS IVA", "GLB1": "MPS IVB",
    "ARSB": "MPS VI", "GUSB": "MPS VII", "HYAL1": "MPS IX",
}
TAXID = 9606

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
DATASETS = "https://api.ncbi.nlm.nih.gov/datasets/v2/"
UNIPROT = "https://rest.uniprot.org/uniprotkb/"

NCBI_API_KEY = os.environ.get("NCBI_API_KEY")
NCBI_EMAIL = os.environ.get("NCBI_EMAIL")      
NCBI_TOOL = "mps_extractor"
NCBI_DELAY = 0.11 if NCBI_API_KEY else 0.35

LINKS = {
    "gene":    "https://www.ncbi.nlm.nih.gov/gene/{}",
    "nuccore": "https://www.ncbi.nlm.nih.gov/nuccore/{}",
    "protein": "https://www.ncbi.nlm.nih.gov/protein/{}",
    "ensembl_tx": "https://www.ensembl.org/Homo_sapiens/Transcript/Summary?t={}",
    "ensembl_gene": "https://www.ensembl.org/Homo_sapiens/Gene/Summary?g={}",
    "uniprot": "https://www.uniprot.org/uniprotkb/{}/entry",
    "omim":    "https://omim.org/entry/{}",
    "alphafold": "https://alphafold.ebi.ac.uk/entry/{}",
    "pdb":     "https://www.rcsb.org/structure/{}",
    "clinvar": "https://www.ncbi.nlm.nih.gov/clinvar/variation/{}/",
    "pubmed":  "https://pubmed.ncbi.nlm.nih.gov/{}/",
}

ECO = {
    "ECO:0000269": "experimental (artículo)",
    "ECO:0007744": "combinada computacional+experimental (manual; p. ej. PDB)",
    "ECO:0000305": "curador (inferencia)",
    "ECO:0000250": "por similitud de secuencia",
    "ECO:0000255": "predicción de secuencia (in silico)",
    "ECO:0000256": "anotación automática (in silico)",
    "ECO:0000303": "autor sin respaldo experimental",
    "ECO:0000312": "importado de otra base",
    "ECO:0007829": "combinada computacional+experimental (automática)",
}
