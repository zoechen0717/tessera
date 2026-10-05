"""Pinned reference sequence bundle.

A bundle is a FASTA plus a manifest naming the assembly, the checksum and the
versioned accession of every contig. The checksum is verified on load: a
silently different FASTA would make every REF check meaningless.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import yaml

from tessera.schemas.base import sha256_hex


class ReferenceError_(Exception):
    pass


@dataclass(frozen=True)
class ReferenceBundle:
    bundle_id: str
    assembly: str
    is_synthetic: bool
    fasta_sha256: str
    contig_accessions: dict[str, str]
    sequences: dict[str, str]
    offsets: dict[str, int]

    @classmethod
    def load(cls, manifest_path: Path) -> "ReferenceBundle":
        manifest = yaml.safe_load(manifest_path.read_text())
        fasta_path = manifest_path.parent / manifest["fasta"]
        raw = fasta_path.read_bytes()
        digest = sha256_hex(raw)
        if digest != manifest["fasta_sha256"]:
            raise ReferenceError_(
                f"FASTA checksum mismatch for {fasta_path.name}: manifest "
                f"{manifest['fasta_sha256'][:12]}…, file {digest[:12]}…"
            )
        sequences = _parse_fasta(raw.decode("ascii"))
        # A contig is either `name: accession` (whole sequence) or
        # `name: {accession, start}` (a slice whose first base is `start`).
        accessions, offsets = {}, {}
        for name, spec in manifest["contigs"].items():
            if isinstance(spec, dict):
                accessions[str(name)], offsets[str(name)] = spec["accession"], int(spec["start"]) - 1
            else:
                accessions[str(name)], offsets[str(name)] = spec, 0
        missing = set(accessions) - set(sequences)
        if missing:
            raise ReferenceError_(f"contigs declared but absent from FASTA: {sorted(missing)}")
        return cls(
            bundle_id=manifest["bundle_id"],
            assembly=manifest["assembly"],
            is_synthetic=bool(manifest["is_synthetic"]),
            fasta_sha256=digest,
            contig_accessions=accessions,
            sequences=sequences,
            offsets=offsets,
        )

    def has_contig(self, chrom: str) -> bool:
        return chrom in self.contig_accessions

    def accession(self, chrom: str) -> str:
        return self.contig_accessions[chrom]

    def span(self, chrom: str) -> tuple[int, int]:
        """1-based inclusive coordinates covered for this contig."""
        off = self.offsets[chrom]
        return off + 1, off + len(self.sequences[chrom])

    def contains(self, chrom: str, pos_1based: int, length: int = 1) -> bool:
        lo, hi = self.span(chrom)
        return lo <= pos_1based and pos_1based + length - 1 <= hi

    def fetch(self, chrom: str, pos_1based: int, length: int) -> str:
        """Reference bases [pos, pos+length) in 1-based genomic coordinates."""
        if not self.contains(chrom, pos_1based, length):
            raise ReferenceError_(f"{chrom}:{pos_1based}+{length} outside the reference bundle")
        i = pos_1based - 1 - self.offsets[chrom]
        return self.sequences[chrom][i : i + length]


def _parse_fasta(text: str) -> dict[str, str]:
    out: dict[str, list[str]] = {}
    name: str | None = None
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            name = line[1:].split()[0]
            out[name] = []
        elif name is None:
            raise ReferenceError_("FASTA sequence before header")
        else:
            out[name].append(line.upper())
    return {k: "".join(v) for k, v in out.items()}
