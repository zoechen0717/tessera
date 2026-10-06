"""Mentions → identity decisions → canonical registry.

Equivalent representations merge on canonical identity; distinct alleles from
one multiallelic record stay distinct. Every mention that does not resolve is
kept as an exclusion with its reason (R-02, INV-01).
"""

from __future__ import annotations

from dataclasses import dataclass, field

from tessera.identity.normalize import (
    NORMALIZER,
    NORMALIZER_VERSION,
    Capabilities,
    classify_text,
    normalize,
)
from tessera.identity.reference import ReferenceBundle
from tessera.schemas.base import content_digest
from tessera.schemas.identity import (
    Alias,
    CanonicalVariant,
    Exclusion,
    IdentityDecision,
    IdentityStatus,
    VariantMention,
)


def split_multiallelic(mention: VariantMention) -> list[VariantMention]:
    alts = (mention.reported_alt or "").split(",")
    if len(alts) < 2:
        return [mention]
    return [
        mention.model_copy(
            update={
                "mention_id": f"{mention.mention_id}.{i}",
                "reported_alt": alt,
                "parent_mention_id": mention.mention_id,
            }
        )
        for i, alt in enumerate(alts, start=1)
    ]


@dataclass
class Registry:
    variants: dict[str, CanonicalVariant] = field(default_factory=dict)
    decisions: list[IdentityDecision] = field(default_factory=list)
    exclusions: list[Exclusion] = field(default_factory=list)
    mentions: list[VariantMention] = field(default_factory=list)
    roles: dict[str, set[str]] = field(default_factory=dict)
    mention_to_variant: dict[str, str] = field(default_factory=dict)
    mention_to_decision: dict[str, str] = field(default_factory=dict)


def build_registry(
    mentions: list[VariantMention],
    bundle: ReferenceBundle,
    gene_id: str,
    capabilities: Capabilities = Capabilities(),
) -> Registry:
    reg = Registry()
    resources = {"reference_bundle": bundle.bundle_id, "normalizer": f"{NORMALIZER}:{NORMALIZER_VERSION}"}

    for parent in mentions:
        for m in split_multiallelic(parent):
            reg.mentions.append(m)
            decision_id = "idd_" + content_digest({"mention": m.mention_id, "resources": resources})[:24]
            reg.mention_to_decision[m.mention_id] = decision_id

            has_coords = m.reported_chrom and m.reported_pos and m.reported_ref and m.reported_alt
            if not has_coords:
                status = m.resolution_status or classify_text(m.raw_text) or IdentityStatus.INSUFFICIENT
                why = m.resolution_note or {
                    IdentityStatus.UNSUPPORTED: "protein-level description; compatible with "
                    "multiple DNA edits and not reverse-translated",
                    IdentityStatus.INSUFFICIENT: "no concrete genomic allele (rsIDs and free "
                    "text are not unique identifiers)",
                }.get(status, "no coordinates")
                _exclude(reg, m, decision_id, status, why, resources)
                continue

            result = normalize(
                bundle,
                m.reported_chrom,
                m.reported_pos,
                m.reported_ref,
                m.reported_alt,
                assembly=m.reported_assembly,
                capabilities=capabilities,
            )
            if result.status is not IdentityStatus.RESOLVED or result.allele is None:
                _exclude(reg, m, decision_id, result.status, "; ".join(result.errors), resources)
                continue

            a = result.allele
            reg.decisions.append(
                IdentityDecision(
                    decision_id=decision_id,
                    mention_id=m.mention_id,
                    status=IdentityStatus.RESOLVED,
                    variant_ids=[a.variant_id],
                    mapping_method=m.resolution_method or "vcf_left_align",
                    resource_versions=resources,
                    justification=(f"{m.resolution_note}; " if m.resolution_note else "") + result.justification,
                )
            )
            alias = Alias(
                raw_value=f"{m.reported_chrom}:{m.reported_pos}:{m.reported_ref}>{m.reported_alt}",
                alias_type="vcf",
                mention_id=m.mention_id,
                decision_id=decision_id,
            )
            existing = reg.variants.get(a.variant_id)
            aliases = (list(existing.aliases) if existing else []) + [alias]
            reg.variants[a.variant_id] = CanonicalVariant(
                variant_id=a.variant_id,
                assembly=bundle.assembly,
                reference_accession=a.reference_accession,
                reference_bundle_id=bundle.bundle_id,
                chrom=a.chrom,
                pos_1based=a.pos_1based,
                ref=a.ref,
                alt=a.alt,
                allele_type=a.allele_type,
                normalization_tool=NORMALIZER,
                normalization_version=NORMALIZER_VERSION,
                gene_ids=[gene_id],
                aliases=aliases,
            )
            reg.roles.setdefault(a.variant_id, set()).update(m.roles)
            reg.mention_to_variant[m.mention_id] = a.variant_id
    return reg


def _exclude(reg: Registry, m: VariantMention, decision_id: str, status: IdentityStatus, why: str, resources) -> None:
    reg.decisions.append(
        IdentityDecision(
            decision_id=decision_id,
            mention_id=m.mention_id,
            status=status,
            mapping_method=m.resolution_method or ("vcf_left_align" if m.reported_chrom else "text_classification"),
            resource_versions=resources,
            justification=why,
            validation_errors=[why],
        )
    )
    reg.exclusions.append(Exclusion(mention_id=m.mention_id, decision_id=decision_id, status=status, reason=why))
