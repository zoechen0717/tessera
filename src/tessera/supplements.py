"""Layer 3a: deterministic supplementary-material scan (DEC-22). No LLM.

For each deep-read paper with open supplements: download the Europe PMC
supplementary zip (verified 2026-10-05: GET /webservices/rest/<PMCID>/supplementaryFiles;
the stream can arrive truncated, so archive integrity is checked and the download
retried), split every member into addressable units (xlsx row, csv row, docx
paragraph/table row, pdf page), and keep variant-like strings from units that
also name the gene. Each hit carries a locator (file, sheet, row / page).

The scan proposes allele MENTIONS; identity is decided downstream by the same
variant_recoder + reference rules as every other route. A c. notation is
assumed to be on the MANE transcript unless the unit names another transcript;
that assumption is recorded on the mention.
"""

from __future__ import annotations

import csv
import io
import re
import zipfile
from dataclasses import dataclass

from tessera.sources.http import Fetcher, SourceError

EUROPEPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest"

C_NOTATION = re.compile(r"\bc\.[-*]?\d+(?:[-+]\d+)?(?:_[-*]?\d+(?:[-+]\d+)?)?(?:[ACGT]>[ACGT]|del[ACGT]*|dup[ACGT]*|ins[ACGT]+|delins[ACGT]+)")
P_NOTATION = re.compile(r"\bp\.\(?(?:[A-Z][a-z]{2}|[A-Z])\d+(?:[A-Z][a-z]{2}|[A-Z]|Ter|\*|=)(?!fs)\)?(?=[\s,;)]|$)")
TRANSCRIPT = re.compile(r"\b(NM_\d+\.\d+|ENST\d+\.\d+)\b")
# Gene-symbol-like token: upper-case alphanumeric containing a digit (SCN8A, SETD1A, KMT2D),
# optionally with one internal space as some tables write "SETD 1A". ACMG codes excluded.
GENE_LIKE = re.compile(r"\b[A-Z][A-Z0-9]{1,8} ?[0-9][A-Z0-9]{0,4}\b")
ACMG = re.compile(r"^(PVS|PS|PM|PP|BA|BS|BP)\d$")
OWNER_WINDOW = 150


@dataclass(frozen=True)
class Unit:
    file: str
    where: str  # e.g. "sheet=Sheet1;row=12" or "page=3" or "para=40"
    text: str


@dataclass(frozen=True)
class SupplementHit:
    pmid: str
    file: str
    where: str
    notation: str
    transcript: str | None
    context: str


def fetch_zip(f: Fetcher, pmcid: str, attempts: int = 3) -> tuple[bytes | None, str]:
    last = "not_attempted"
    for i in range(attempts):
        try:
            # first try may come from cache; a bad cached archive forces a network refetch
            r = f.get("europepmc", f"supplementaryFiles/{pmcid}", "live",
                      f"{EUROPEPMC}/{pmcid}/supplementaryFiles", media="application/zip",
                      ttl_days=float("inf") if i == 0 else None)
        except SourceError as e:
            last = e.code
            continue
        try:
            z = zipfile.ZipFile(io.BytesIO(r.content))
            if z.testzip() is None:
                return r.content, r.snapshot_id
            last = "corrupt_member"
        except zipfile.BadZipFile:
            last = "truncated_archive"
    return None, last


def units(zbytes: bytes) -> list[Unit]:
    out: list[Unit] = []
    z = zipfile.ZipFile(io.BytesIO(zbytes))
    for name in z.namelist():
        low = name.lower()
        data = z.read(name)
        try:
            if low.endswith((".xlsx", ".xlsm")):
                import openpyxl
                wb = openpyxl.load_workbook(io.BytesIO(data), read_only=True, data_only=True)
                for ws in wb.worksheets:
                    for i, row in enumerate(ws.iter_rows(values_only=True), start=1):
                        cells = [str(c) for c in row if c is not None]
                        if cells:
                            out.append(Unit(name, f"sheet={ws.title};row={i}", " | ".join(cells)))
            elif low.endswith((".csv", ".tsv", ".txt")):
                text = data.decode("utf-8", errors="replace")
                delim = "\t" if low.endswith(".tsv") or text.count("\t") > text.count(",") else ","
                for i, row in enumerate(csv.reader(io.StringIO(text), delimiter=delim), start=1):
                    if row:
                        out.append(Unit(name, f"row={i}", " | ".join(row)))
            elif low.endswith(".docx"):
                doc = zipfile.ZipFile(io.BytesIO(data)).read("word/document.xml").decode("utf-8", errors="replace")
                for i, para in enumerate(re.findall(r"<w:(?:p|tr)[ >].*?</w:(?:p|tr)>", doc, re.S), start=1):
                    t = re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", para)).strip()
                    if t:
                        out.append(Unit(name, f"block={i}", t))
            elif low.endswith(".pdf"):
                from pypdf import PdfReader
                for i, page in enumerate(PdfReader(io.BytesIO(data)).pages, start=1):
                    t = re.sub(r"\s+", " ", page.extract_text() or "")
                    if t:
                        out.append(Unit(name, f"page={i}", t))
        except Exception as e:  # unreadable member: recorded as a unit-level failure
            out.append(Unit(name, "unreadable", f"__unreadable__ {type(e).__name__}"))
    return out


def _gene_in(text: str, aliases: list[str]) -> bool:
    squashed = re.sub(r"\s+", "", text).upper()
    return any(a.upper() in squashed for a in aliases)


def owner_is_target(text: str, at: int, aliases: list[str]) -> bool:
    """The nearest gene-like token before position `at` (within OWNER_WINDOW chars)
    must be the target gene. Multi-gene tables and PDF pages otherwise attach other
    genes' variants to the target -- and a protein change of another gene can pass
    the reference check by coincidence (same residue number, same amino acid)."""
    targets = {re.sub(r"\s+", "", a).upper() for a in aliases}
    last = None
    for m in GENE_LIKE.finditer(text, max(0, at - OWNER_WINDOW), at):
        tok = re.sub(r"\s+", "", m.group(0)).upper()
        if not ACMG.match(tok):
            last = tok
    return last in targets


def scan(pmid: str, us: list[Unit], aliases: list[str], window: int = 300) -> list[SupplementHit]:
    """Variant-like strings attributed to the gene: the unit must name the gene and
    the nearest preceding gene symbol of each notation must be the gene."""
    hits: list[SupplementHit] = []
    for u in us:
        if u.text.startswith("__unreadable__") or not _gene_in(u.text, aliases):
            continue
        if u.where.startswith(("page=", "block=")) and len(u.text) > 2 * window:
            spans = []
            squashed_positions = [m.start() for a in aliases for m in re.finditer(re.escape(a), u.text, re.I)]
            for p in squashed_positions:
                spans.append(u.text[max(0, p - window) : p + window])
            texts = spans or [u.text]
        else:
            texts = [u.text]
        tx = TRANSCRIPT.search(u.text)
        for t in texts:
            for m in list(C_NOTATION.finditer(t)) + list(P_NOTATION.finditer(t)):
                if not owner_is_target(t, m.start(), aliases):
                    continue
                hits.append(SupplementHit(pmid, u.file, u.where, m.group(0).replace("(", "").replace(")", "").strip(),
                                          tx.group(1) if tx else None, t[:200]))
    seen, uniq = set(), []
    for h in hits:
        k = (h.file, h.where, h.notation)
        if k not in seen:
            seen.add(k)
            uniq.append(h)
    return uniq
