import datetime as dt
import time
from urllib.parse import urlencode

import requests

from . import config as C

QUERY_LOG: list[dict] = []
SOURCE_VERSIONS: dict = {}
_session = requests.Session()


def _log(base, endpoint, query, url, n):
    QUERY_LOG.append(dict(
        fecha_hora=dt.datetime.now().isoformat(timespec="seconds"),
        base=base, endpoint=endpoint, consulta=query, url=url, registros=n))


def _request(method, url, retries=5, **kw):
    for i in range(retries):
        try:
            r = _session.request(method, url, timeout=kw.pop("timeout", 120) if i == 0 else 300, **kw)
        except requests.exceptions.RequestException as e:
            print(f"  reintento {i + 1} ({type(e).__name__}) {url}", flush=True)
            time.sleep(5 * (i + 1)); continue
        if r.status_code in (429, 500, 502, 503, 504):
            print(f"  reintento {i + 1} (HTTP {r.status_code}) {url}", flush=True)
            time.sleep(5 * (i + 1)); continue
        r.raise_for_status()
        return r
    raise RuntimeError(f"Sin respuesta válida tras {retries} intentos: {url}")


def _get(url, params=None, retries=5):
    return _request("GET", url, retries, params=params)


def _post(url, data, retries=5):
    return _request("POST", url, retries, data=data)


def _ncbi_params(**kw):
    p = dict(kw, tool=C.NCBI_TOOL)
    if C.NCBI_API_KEY: p["api_key"] = C.NCBI_API_KEY
    if C.NCBI_EMAIL: p["email"] = C.NCBI_EMAIL
    return p


def _public_url(url, params):
    safe = {k: v for k, v in params.items() if k not in ("api_key", "email")}
    return f"{url}?{urlencode(safe)}"


def esearch(db, term, retmax=10000, base=None):
    url = C.EUTILS + "esearch.fcgi"
    p = _ncbi_params(db=db, term=term, retmode="json", retmax=retmax)
    res = _get(url, p).json()["esearchresult"]
    time.sleep(C.NCBI_DELAY)
    n = int(res["count"])
    _log(base or f"NCBI {db}", "E-utilities esearch", term, _public_url(url, p), n)
    return n, res["idlist"]


def esummary(db, ids, batch=200):
    out = {}
    url = C.EUTILS + "esummary.fcgi"
    for i in range(0, len(ids), batch):
        chunk = ids[i:i + batch]
        p = _ncbi_params(db=db, id=",".join(chunk), retmode="json")
        r = _post(url, p)
        res = r.json()["result"]
        out.update({u: res[u] for u in res.get("uids", [])})
        time.sleep(C.NCBI_DELAY)
    return out


def efetch_fasta(db, acc):
    url = C.EUTILS + "efetch.fcgi"
    p = _ncbi_params(db=db, id=acc, rettype="fasta", retmode="text")
    txt = _get(url, p).text
    time.sleep(C.NCBI_DELAY)
    _log(f"NCBI {db}", "E-utilities efetch (FASTA)", acc, _public_url(url, p), 1)
    return "".join(l.strip() for l in txt.splitlines() if not l.startswith(">"))

def datasets(path, query_desc):
    url = C.DATASETS + path
    js = _get(url).json()
    time.sleep(C.NCBI_DELAY)
    _log("NCBI Datasets v2", path.split("/")[0] + "/…/" + path.split("/")[-1], query_desc, url,
         len(js.get("reports", [])))
    return js


# ---------------- UniProt REST ----------------
def uniprot_search(query, fields=None, size=10):
    url = C.UNIPROT + "search"
    p = dict(query=query, format="json", size=size)
    if fields: p["fields"] = fields
    r = _get(url, p)
    js = r.json()
    _log("UniProtKB", "REST /uniprotkb/search", query, f"{url}?{urlencode(p)}",
         int(r.headers.get("x-total-results", len(js["results"]))))
    return js["results"]


def uniprot_entry(acc):
    url = C.UNIPROT + f"{acc}.json"
    r = _get(url); js = r.json()
    SOURCE_VERSIONS["UniProtKB release"] = r.headers.get("X-UniProt-Release")
    SOURCE_VERSIONS["UniProtKB release date"] = r.headers.get("X-UniProt-Release-Date")
    _log("UniProtKB", "REST /uniprotkb/{accession}.json", acc, url, 1)
    return js


def capture_versions():
    """Versiones de RefSeq/UniProt en cada ejecución (también cuando los genes vienen de caché)."""
    r = _get(C.UNIPROT + "P35475.json")
    SOURCE_VERSIONS["UniProtKB release"] = r.headers.get("X-UniProt-Release")
    SOURCE_VERSIONS["UniProtKB release date"] = r.headers.get("X-UniProt-Release-Date")
    ann = _get(C.DATASETS + "gene/id/3425").json()["reports"][0]["gene"]["annotations"][0]
    SOURCE_VERSIONS["NCBI RefSeq annotation"] = (f"{ann.get('annotation_name')} ({ann.get('annotation_release_date')}), "
                                                 f"{ann.get('assembly_name')}")
    SOURCE_VERSIONS["ClinVar consultado"] = dt.date.today().isoformat()
