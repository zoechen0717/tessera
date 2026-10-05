"""ClinVar (E-utilities), UniProt, ClinGen dosage.

Verified 2026-10-05:
- ClinVar esearch/esummary (db=clinvar, retmode=json): records carry
  accession, germline_classification {description, review_status, trait_set},
  variation_set[] with canonical_spdi on GRCh38.
- UniProt REST https://rest.uniprot.org/uniprotkb/<acc>.json: features + sequence.
- ClinGen https://ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh38.tsv:
  Haploinsufficiency Score per gene, with a disease ID column.
"""

from __future__ import annotations

from datetime import date

from tessera.sources.http import Fetcher, SourceError

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
UNIPROT = "https://rest.uniprot.org/uniprotkb"
CLINGEN_TSV = "https://ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh38.tsv"


def clinvar_records(f: Fetcher, symbol: str, batch: int = 200) -> list[tuple[dict, str]]:
    version = f"eutils-live-{date.today().isoformat()}"
    s = f.get("clinvar", f"esearch/{symbol}", version, f"{EUTILS}/esearch.fcgi",
              {"db": "clinvar", "term": f"{symbol}[gene]", "retmax": 10000, "retmode": "json"})
    res = s.json()["esearchresult"]
    ids = res["idlist"]
    if int(res["count"]) != len(ids):
        raise SourceError("clinvar", "truncated", f"count {res['count']} but {len(ids)} ids returned")
    out = []
    for i in range(0, len(ids), batch):
        chunk = ids[i : i + batch]
        r = f.get("clinvar", f"esummary/{symbol}/{i // batch}", version, f"{EUTILS}/esummary.fcgi",
                  {"db": "clinvar", "id": ",".join(chunk), "retmode": "json"})
        result = r.json()["result"]
        for uid in result.get("uids", []):
            out.append((result[uid], r.snapshot_id))
    return out


def clinvar_genes(rec: dict) -> set[str]:
    return {g.get("symbol") for g in rec.get("genes", [])}


def parse_spdi(spdi: str) -> tuple[str, int, str, str]:
    """`ACC:pos0:del:ins` → (accession, 0-based position, deleted, inserted)."""
    acc, pos0, dele, ins = spdi.split(":")
    return acc, int(pos0), dele, ins


def spdi_to_vcf(pos0: int, dele: str, ins: str, anchor_base) -> tuple[int, str, str]:
    """SPDI (0-based, interbase deletion) → VCF (1-based, anchored when an allele is empty)."""
    if dele and ins:
        return pos0 + 1, dele, ins
    anchor = anchor_base(pos0)  # 1-based position pos0 is the base before the event
    return pos0, anchor + dele, anchor + ins


def clinvar_classification(rec: dict) -> dict:
    g = rec.get("germline_classification") or {}
    return {
        "accession": rec.get("accession"),
        "description": g.get("description") or None,
        "review_status": g.get("review_status") or None,
        "traits": sorted({t.get("trait_name") for t in g.get("trait_set", []) if t.get("trait_name")}),
    }


def uniprot_entry(f: Fetcher, accession: str) -> tuple[dict, str]:
    r = f.get("uniprot", accession, "live", f"{UNIPROT}/{accession}.json")
    d = r.json()
    return d, r.snapshot_id


# Declared rule for `in_functional_domain` (provisional): UniProt Domain and
# Motif features, and Region features annotated as an interaction.
def functional_regions(entry: dict) -> list[tuple[int, int, str]]:
    out = []
    for ft in entry.get("features", []):
        kind, desc = ft.get("type"), ft.get("description") or ""
        if kind in ("Domain", "Motif") or (kind == "Region" and desc.startswith("Interaction with")):
            loc = ft["location"]
            out.append((int(loc["start"]["value"]), int(loc["end"]["value"]), f"{kind}:{desc}"))
    return sorted(out)


def clingen_row(f: Fetcher, symbol: str) -> tuple[dict | None, str]:
    r = f.get("clingen", "gene_curation_list_GRCh38", "live", CLINGEN_TSV, media="text/plain")
    header = None
    for line in r.text().splitlines():
        if line.startswith("#Gene Symbol"):
            header = line[1:].split("\t")
            continue
        if line.startswith("#") or not header:
            continue
        cols = line.split("\t")
        if cols[0] == symbol:
            return dict(zip(header, cols)), r.snapshot_id
    return None, r.snapshot_id
