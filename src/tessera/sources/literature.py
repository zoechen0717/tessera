"""Deterministic literature routes: LitVar2 gene variants + Ensembl variant_recoder.

Verified 2026-10-05:
- LitVar2 GET https://www.ncbi.nlm.nih.gov/research/litvar2-api/variant/search/gene/<SYMBOL>
  returns one Python-literal dict per line: {'_id': 'litvar@rs123##' | 'litvar@#<geneid>#<hgvs>',
  'pmids_count': n, 'rsid'?}. GET /variant/get/<id>/publications returns {pmids, pmcids}.
- Ensembl POST /variant_recoder/homo_sapiens {"ids": [...], "fields": "spdi"} returns,
  per resolvable input, {ALT: {input, spdi: [...]}}; unresolvable inputs are omitted
  (their warnings may be attached to a neighbouring record), so results are matched
  by `input`, never by position. rsIDs expand to every ALT; protein changes expand to
  every compatible SNV; REF disagreement is rejected.

Resolution rule (DEC-21): exactly one compatible genomic allele → candidate;
several → ambiguous; none → unmapped. If the literature text states the deleted
sequence (e.g. `delAG`) its length must match the mapped deletion, else ambiguous.
A literature-derived candidate is a concrete allele to consider; evidence from the
same paper still links only at protein_equivalent_only unless the paper reports
the genomic/cDNA allele (R-09).
"""

from __future__ import annotations

import ast
import re

from tessera.sources.http import Fetcher, SourceError

LITVAR = "https://www.ncbi.nlm.nih.gov/research/litvar2-api"
ENSEMBL = "https://rest.ensembl.org"
RECODER_BATCH = 50

def litvar_gene_variants(f: Fetcher, symbol: str) -> tuple[list[dict], str]:
    r = f.get("litvar2", f"search/gene/{symbol}", "live", f"{LITVAR}/variant/search/gene/{symbol}",
              media="text/plain", ttl_days=7)
    rows = []
    for line in r.text().splitlines():
        line = line.strip()
        if line:
            rows.append(ast.literal_eval(line))  # literals only; never eval
    return rows, r.snapshot_id


def litvar_name(row: dict) -> str | None:
    """'litvar@#9739#p.R913C' → 'p.R913C'; rsID rows and the gene aggregate → None."""
    vid = row["_id"]
    if row.get("rsid"):
        return None
    name = vid.split("#")[-1]
    return name or None


def recoder_input(name: str, transcript: str, protein: str) -> str | None:
    if name.startswith("p."):
        if "fs" in name:
            return None  # frameshift protein notation does not determine a DNA change
        return f"{protein}:{name}"
    if name.startswith("c."):
        return f"{transcript}:{name}"
    return None


def recoder(f: Fetcher, ids: list[str], rel: str, failed: set[str] | None = None) -> dict[str, list[str]]:
    """{input: [spdi, ...]} for every input the recoder resolved. Inputs of a batch whose
    request failed are added to `failed` (if given) instead of aborting the run."""
    out: dict[str, list[str]] = {}
    for i in range(0, len(ids), RECODER_BATCH):
        chunk = ids[i : i + RECODER_BATCH]
        try:
            r = f.post_json("ensembl", f"variant_recoder/batch{i // RECODER_BATCH}", rel,
                            f"{ENSEMBL}/variant_recoder/homo_sapiens", {"ids": chunk, "fields": "spdi"}, ttl_days=90)
        except SourceError:
            if failed is None:
                raise
            failed.update(chunk)
            continue
        for rec in r.json():
            if not isinstance(rec, dict):
                continue  # the recoder returns null for some unparseable inputs
            for alt, v in rec.items():
                if isinstance(v, dict) and "input" in v:
                    out.setdefault(v["input"], []).extend(v.get("spdi") or [])
    return {k: sorted(set(v)) for k, v in out.items()}


_DEL_SEQ = re.compile(r"del([ACGT]+)$")


def stated_deletion_ok(name: str, spdi: str) -> bool:
    m = _DEL_SEQ.search(name)
    if not m:
        return True
    _, _, dele, ins = spdi.split(":")
    return len(dele) - len(ins) == len(m.group(1))


# ------------------------------------------------------------- document pool
# Verified 2026-10-05:
# - Europe PMC GET /webservices/rest/search?query=…&format=json&resultType=core&pageSize≤1000
#   &cursorMark=* — core records include abstractText, pubTypeList, isOpenAccess, inPMC,
#   hasSuppl, pmcid, doi.
# - PubTator3 GET /research/pubtator3-api/search/?text=@GENE_X AND @DISEASE_Y&page=n
#   (10 results per page; `count`, `total_pages`).
# - LitVar2 GET /variant/get/<id>/publications → {pmids, pmcids}.

import json as _json
import urllib.parse as _up

EUROPEPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest"
PUBTATOR = "https://www.ncbi.nlm.nih.gov/research/pubtator3-api"


def europepmc_search(f: Fetcher, query: str, max_records: int = 5000) -> tuple[list[dict], list[str]]:
    out, snaps, cursor = [], [], "*"
    while len(out) < max_records:
        r = f.get("europepmc", f"search/{query}/{cursor}", "live", f"{EUROPEPMC}/search",
                  {"query": query, "format": "json", "resultType": "core", "pageSize": 1000, "cursorMark": cursor})
        d = r.json()
        snaps.append(r.snapshot_id)
        out.extend(d["resultList"]["result"])
        nxt = d.get("nextCursorMark")
        if not d["resultList"]["result"] or not nxt or nxt == cursor:
            break
        cursor = nxt
    return out, snaps


def europepmc_by_pmids(f: Fetcher, pmids: list[str], batch: int = 50) -> tuple[list[dict], list[str]]:
    out, snaps = [], []
    for i in range(0, len(pmids), batch):
        chunk = pmids[i : i + batch]
        q = "(" + " OR ".join(f"EXT_ID:{p}" for p in chunk) + ") AND SRC:MED"
        r = f.get("europepmc", f"pmids/{i // batch}", "live", f"{EUROPEPMC}/search",
                  {"query": q, "format": "json", "resultType": "core", "pageSize": 1000}, ttl_days=30)
        out.extend(r.json()["resultList"]["result"])
        snaps.append(r.snapshot_id)
    return out, snaps


def pubtator_search(f: Fetcher, text: str, max_pages: int = 100) -> tuple[list[str], list[str], int]:
    pmids, snaps, page, total = [], [], 1, 0
    while page <= max_pages:
        r = f.get("pubtator3", f"search/{text}/{page}", "live", f"{PUBTATOR}/search/", {"text": text, "page": page})
        d = r.json()
        snaps.append(r.snapshot_id)
        total = d.get("count", 0)
        pmids.extend(str(x["pmid"]) for x in d.get("results", []) if x.get("pmid"))
        if page >= int(d.get("total_pages") or 0):
            break
        page += 1
    return pmids, snaps, total


def litvar_publications(f: Fetcher, litvar_id: str) -> list[str]:
    r = f.get("litvar2", f"publications/{litvar_id}", "live",
              f"{LITVAR}/variant/get/{_up.quote(litvar_id, safe='')}/publications", ttl_days=30)
    return [str(p) for p in _json.loads(r.content).get("pmids", [])]


def europepmc_fulltext(f: Fetcher, pmcid: str) -> str | None:
    """Open-access full text XML (Europe PMC); None when not available."""
    try:
        r = f.get("europepmc", f"fullTextXML/{pmcid}", "live", f"{EUROPEPMC}/{pmcid}/fullTextXML",
                  media="text/plain", ttl_days=float("inf"))
    except SourceError:
        return None
    return r.text()


# ------------------------------------------------------- notation cleanup
import html as _html

_HYPHENS = dict.fromkeys(map(ord, "‐‑‒–—−"), "-")
_TOKEN = re.compile(r"(?:(NM_\d+(?:\.\d+)?)\s*:\s*)?(?:[A-Z0-9]+\s*:\s*)?([cp])\.\s*\(?([^\s,;()]+(?:\s+[ACGT]+\b)?)")


def clean_notations(text: str) -> list[tuple[str, str | None]]:
    """Split literature allele text into (notation, refseq transcript or None).

    Handles HTML entities, typographic hyphens/asterisks, spaces inside notations
    ("c.3930_3940del CCCTGCGCCAG", "c.1067C > T"), zygosity suffixes ("/+"),
    transcript prefixes ("NM_014712.2 :c.…") and several alleles in one string.
    When a c. and a p. description of the same allele are given, only c. is kept.
    Bare genomic descriptions without an assembly are dropped (cannot be placed).
    """
    t = _html.unescape(text).translate(_HYPHENS).replace("∗", "*")
    t = re.sub(r"\s*>\s*", ">", t)
    t = re.sub(r"/(?=\s*[cp]\.)", " ", t)  # "c.1067C>T/p.Ser356Phe" → two descriptions
    t = re.sub(r"(del|dup|ins)\s+([ACGT]+)\b", r"\1\2", t)
    t = re.sub(r"\s*\*\s*", "*", t)
    out = []
    for m in _TOKEN.finditer(t):
        tx, kind, body = m.group(1), m.group(2), m.group(3).replace(" ", "")
        body = re.sub(r"/\+$|/\+?wt$|/WT$", "", body).rstrip(".")
        out.append((f"{kind}.{body}", tx))
    out = [(expand_single_position_deletion(n), tx) for n, tx in out]
    cs = [o for o in out if o[0].startswith("c.")]
    return cs or [o for o in out if o[0].startswith("p.")]


_SINGLE_DEL = re.compile(r"^c\.(\d+)([-+]\d+)?del([ACGT]{2,})$")


def expand_single_position_deletion(note: str) -> str:
    """Common shorthand 'c.4582-2delAG' (one position, two deleted bases) →
    'c.4582-2_4582-1delAG'. Only expanded when the span stays on one side of the
    exon/intron boundary; otherwise the notation is returned unchanged (and will be
    judged ambiguous by the stated-deletion check)."""
    m = _SINGLE_DEL.match(note)
    if not m:
        return note
    base, off, seq = int(m.group(1)), m.group(2), m.group(3)
    n = len(seq) - 1
    if off is None:
        return f"c.{base}_{base + n}del{seq}"
    o = int(off)
    end = o + n
    if (o < 0 and end >= 0) or (o > 0 and end <= 0):
        return note  # would cross into the exon: do not guess
    return f"c.{base}{o:+d}_{base}{end:+d}del{seq}"
