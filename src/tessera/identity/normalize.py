"""Genomic allele normalization (ARCHITECTURE §5, steps 3, 7, 8).

M0 accepts concrete VCF-style alleles only (DEC-11). Protein shorthand and bare
rsIDs are recorded and left unresolved: neither identifies one genomic allele.
Normalization is the standard trim-and-left-align: REF is validated against the
pinned reference before and after.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field

from tessera.identity.reference import ReferenceBundle
from tessera.schemas.identity import (
    IDENTITY_VERSION,
    AlleleType,
    IdentityStatus,
    compute_variant_id,
)

NORMALIZER = "tessera.left_align"
NORMALIZER_VERSION = "1"

_DNA = re.compile(r"^[ACGT]+$")
_PROTEIN = re.compile(r"^p\.\(?[A-Z][a-z]{2}\d+")
_RSID = re.compile(r"^rs\d+")


@dataclass(frozen=True)
class Capabilities:
    allowed_types: frozenset[AlleleType] = frozenset(
        {AlleleType.SNV, AlleleType.DELETION, AlleleType.INSERTION, AlleleType.DELINS}
    )
    max_indel_length: int = 50


@dataclass(frozen=True)
class NormalizedAllele:
    chrom: str
    pos_1based: int
    ref: str
    alt: str
    allele_type: AlleleType
    reference_accession: str
    variant_id: str


@dataclass
class NormalizationResult:
    status: IdentityStatus
    allele: NormalizedAllele | None = None
    errors: list[str] = field(default_factory=list)
    justification: str = ""


def classify_text(raw_text: str) -> IdentityStatus | None:
    """Status for mentions that carry no genomic coordinates."""
    text = raw_text.strip()
    if _PROTEIN.match(text):
        return IdentityStatus.UNSUPPORTED
    if _RSID.match(text):
        return IdentityStatus.INSUFFICIENT
    return None


def allele_type(ref: str, alt: str) -> AlleleType:
    if len(ref) == len(alt):
        return AlleleType.SNV if len(ref) == 1 else AlleleType.MNV
    if len(alt) == 1 and ref[0] == alt:
        return AlleleType.DELETION
    if len(ref) == 1 and alt[0] == ref:
        return AlleleType.INSERTION
    return AlleleType.DELINS


def _left_align(ref_bundle: ReferenceBundle, chrom: str, pos: int, ref: str, alt: str) -> tuple[int, str, str]:
    """Trim shared suffix (extending left when an allele empties), then shared prefix."""
    while True:
        changed = False
        if ref and alt and ref[-1] == alt[-1]:
            ref, alt = ref[:-1], alt[:-1]
            changed = True
        if not ref or not alt:
            if pos <= 1:
                raise ValueError("cannot extend allele left of contig start")
            pos -= 1
            base = ref_bundle.fetch(chrom, pos, 1)
            ref, alt = base + ref, base + alt
            changed = True
        if not changed:
            break
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt = ref[1:], alt[1:]
        pos += 1
    return pos, ref, alt


def normalize(
    ref_bundle: ReferenceBundle,
    chrom: str,
    pos_1based: int,
    ref: str,
    alt: str,
    *,
    assembly: str | None = None,
    capabilities: Capabilities = Capabilities(),
) -> NormalizationResult:
    ref, alt = ref.upper(), alt.upper()
    if assembly is not None and assembly != ref_bundle.assembly:
        return NormalizationResult(
            IdentityStatus.UNSUPPORTED,
            errors=[f"reported assembly {assembly} != bundle {ref_bundle.assembly}; "
                    "variant-aware assembly mapping is not an M0 capability"],
        )
    if not ref_bundle.has_contig(chrom):
        return NormalizationResult(IdentityStatus.UNMAPPED, errors=[f"unknown contig {chrom}"])
    if not _DNA.match(ref) or not _DNA.match(alt):
        return NormalizationResult(
            IdentityStatus.UNSUPPORTED, errors=["symbolic or non-ACGT allele"]
        )
    if ref == alt:
        return NormalizationResult(IdentityStatus.INSUFFICIENT, errors=["ref equals alt"])
    try:
        observed = ref_bundle.fetch(chrom, pos_1based, len(ref))
    except Exception as exc:  # out of range
        return NormalizationResult(IdentityStatus.REF_MISMATCH, errors=[str(exc)])
    if observed != ref:
        return NormalizationResult(
            IdentityStatus.REF_MISMATCH,
            errors=[f"REF {ref} != reference {observed} at {chrom}:{pos_1based}"],
        )

    pos, nref, nalt = _left_align(ref_bundle, chrom, pos_1based, ref, alt)
    if ref_bundle.fetch(chrom, pos, len(nref)) != nref:  # defensive re-check (step 8)
        return NormalizationResult(IdentityStatus.REF_MISMATCH, errors=["post-normalization REF check failed"])

    atype = allele_type(nref, nalt)
    if atype not in capabilities.allowed_types:
        return NormalizationResult(
            IdentityStatus.UNSUPPORTED, errors=[f"allele type {atype.value} not in capability set"]
        )
    if abs(len(nref) - len(nalt)) > capabilities.max_indel_length:
        return NormalizationResult(IdentityStatus.UNSUPPORTED, errors=["indel exceeds capability length"])

    accession = ref_bundle.accession(chrom)
    vid = compute_variant_id(
        identity_version=IDENTITY_VERSION,
        assembly=ref_bundle.assembly,
        reference_accession=accession,
        pos_1based=pos,
        ref=nref,
        alt=nalt,
    )
    moved = (pos, nref, nalt) != (pos_1based, ref, alt)
    return NormalizationResult(
        IdentityStatus.RESOLVED,
        allele=NormalizedAllele(chrom, pos, nref, nalt, atype, accession, vid),
        justification="REF validated; " + ("trimmed/left-aligned" if moved else "already normalized"),
    )
