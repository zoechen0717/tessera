"""Closed vocabularies.

Every field in the evidence model whose value space is bounded uses one of these
rather than a free string. This is what makes the distinctions in the spec
mechanically enforceable instead of conventional -- an LLM cannot invent a new
evidence level, and a scorer cannot silently treat two states as equivalent.
"""

from __future__ import annotations

from enum import StrEnum


class EvidenceLevel(StrEnum):
    """Spec section 12.2. Strictly ordered: lower ordinal = more direct.

    The single most important distinction in the system is VARIANT_DIRECT vs
    GENE_LEVEL. Collapsing them is the failure mode that turns "SETD1A loss of
    function affects neuronal biology" into "variant Q191H has direct functional
    validation".
    """

    VARIANT_DIRECT = "variant_direct"
    VARIANT_ASSOCIATION = "variant_association"
    GENE_LEVEL = "gene_level"
    LOCUS_LEVEL = "locus_level"
    PREDICTION = "prediction"

    @property
    def rank(self) -> int:
        return _EVIDENCE_LEVEL_RANK[self]

    def is_variant_specific(self) -> bool:
        return self in (EvidenceLevel.VARIANT_DIRECT, EvidenceLevel.VARIANT_ASSOCIATION)


_EVIDENCE_LEVEL_RANK = {
    EvidenceLevel.VARIANT_DIRECT: 1,
    EvidenceLevel.VARIANT_ASSOCIATION: 2,
    EvidenceLevel.GENE_LEVEL: 3,
    EvidenceLevel.LOCUS_LEVEL: 4,
    EvidenceLevel.PREDICTION: 5,
}


class DiseaseMatchType(StrEnum):
    """Spec section 12.1.1.

    A gene may have abundant, high-quality clinical records whose condition is a
    DIFFERENT disease from the query. Presenting those as evidence for the query
    disease is the most likely way for this system to be confidently wrong.
    """

    EXACT = "exact"
    NARROWER = "narrower"
    BROADER = "broader"
    RELATED = "related"
    UNRELATED = "unrelated"


class EvidenceDirection(StrEnum):
    """Spec section 12.3. Conflict is preserved, not averaged away."""

    SUPPORTS_PRIORITY = "supports_priority"
    CONTRADICTS_PRIORITY = "contradicts_priority"
    NEUTRAL = "neutral"
    UNCERTAIN = "uncertain"


class Missingness(StrEnum):
    """Spec section 42. These five states are NOT interchangeable.

    NOT_SEARCHED scored as 0.0 is the silent failure that fills a panel with
    variants the pipeline never looked at: the novelty multiplier rewards low
    saturation, so an unresearched variant would be promoted as a discovery
    candidate. See ranking.components and section 42.1.
    """

    PRESENT = "present"
    SEARCHED_NOT_FOUND = "searched_not_found"
    NOT_SEARCHED = "not_searched"
    NOT_AVAILABLE = "not_available"
    NOT_APPLICABLE = "not_applicable"
    FAILED_RETRIEVAL = "failed_retrieval"


class ConsequenceClass(StrEnum):
    """Spec section 19.3. Selects which M branch is used, and the rank stratum.

    Cross-class score comparison is not meaningful: M_coding and M_splice are
    computed from different predictors on different scales.
    """

    MISSENSE = "missense"
    SPLICE = "splice"
    PTV = "ptv"
    SYNONYMOUS = "synonymous"
    NONCODING = "noncoding"
    INFRAME_INDEL = "inframe_indel"
    OTHER = "other"


class Mechanism(StrEnum):
    """Spec section 9.8. Gates which variant classes are worth screening at all."""

    HAPLOINSUFFICIENCY = "haploinsufficiency"
    GAIN_OF_FUNCTION = "gain_of_function"
    DOMINANT_NEGATIVE = "dominant_negative"
    MIXED = "mixed"
    UNKNOWN = "unknown"


class MechanismConfidence(StrEnum):
    HIGH = "high"
    MODERATE = "moderate"
    LOW = "low"


class VariantStatus(StrEnum):
    """Spec section 17.1.2. Editability gate.

    EDITABILITY_UNKNOWN is distinct from NOT_EDITABLE_* and must not be scored
    as zero -- the host line simply has not been genotyped at that position.
    """

    EDITABLE = "editable"
    NOT_EDITABLE_NO_PAM = "not_editable_no_pam"
    NOT_EDITABLE_IN_HOST_LINE = "not_editable_in_host_line"
    NOT_EDITABLE_LENGTH = "not_editable_length"
    EDITABILITY_UNKNOWN = "editability_unknown"

    def blocks_panel_inclusion(self) -> bool:
        return self.startswith("not_editable")


class VerificationStatus(StrEnum):
    """Spec sections 14, 50."""

    PENDING = "pending"
    REJECTED_DETERMINISTIC = "rejected_deterministic"
    AGENT_VERIFIED = "agent_verified"
    HUMAN_VERIFIED = "human_verified"
    HUMAN_REJECTED = "human_rejected"

    def usable_in_scoring(self) -> bool:
        return self in (
            VerificationStatus.AGENT_VERIFIED,
            VerificationStatus.HUMAN_VERIFIED,
        )


class LabelLeakage(StrEnum):
    """Spec sections 9.5.2, 47.3, 48.10.

    Benchmarking a ClinVar-supervised predictor against a ClinVar-derived answer
    key is circular. Metrics are reported per leakage class, never pooled.

    ALLELE_FREQUENCY covers predictors trained on population-frequency proxy labels
    (AVI is trained to separate filtering allele frequency above vs below 0.1%).
    Such a predictor is relatively fair against a ClinVar key, but it is NOT
    independent of any score component that also uses allele frequency -- which in
    this system means the rarity term in D.
    """

    NONE = "none"
    CLINVAR = "clinvar"
    DMS = "dms"
    HGMD = "hgmd"
    ALLELE_FREQUENCY = "allele_frequency"
    UNKNOWN = "unknown"


class RankingMode(StrEnum):
    DISCOVERY = "discovery"
    VALIDATION = "validation"
    BALANCED = "balanced"
    CONTROL = "control"
    PANEL = "panel"


class ScoreComponent(StrEnum):
    """The four weighted components. Novelty is deliberately absent: it is a
    multiplier, not a component (spec 19.5)."""

    DISEASE = "disease"
    FUNCTIONAL = "functional"
    MECHANISM = "mechanism"
    EDITABILITY = "editability"


class AnchorType(StrEnum):
    """Spec section 19.6. Every mapping is absolute; none is candidate-set
    dependent, because that would make one variant's score a function of which
    other variants happen to be in the run."""

    IDENTITY = "identity"
    CATEGORICAL = "categorical"
    PIECEWISE = "piecewise"
    LINEAR_CLAMP = "linear_clamp"
    BOOLEAN = "boolean"


class RunMode(StrEnum):
    FROZEN = "frozen"
    REFRESH = "refresh"
    COMPARE = "compare"


class ControlType(StrEnum):
    """Spec section 8.5."""

    KNOWN_PATHOGENIC = "known_pathogenic"
    KNOWN_BENIGN = "known_benign"
    COMMON_VARIANT = "common_variant"
    SYNONYMOUS = "synonymous"
    CONSEQUENCE_MATCHED = "consequence_matched"
    EDITABILITY_MATCHED = "editability_matched"
    OUTSIDE_DOMAIN = "outside_domain"
