from __future__ import annotations

from pathlib import Path

import pytest

from tessera.config.scoring_config import ScoringConfig
from tessera.enums import (
    ConsequenceClass,
    DiseaseMatchType,
    EvidenceLevel,
    Mechanism,
    MechanismConfidence,
    VerificationStatus,
)
from tessera.schemas.annotations import VariantAnnotations
from tessera.schemas.common import Quantity
from tessera.schemas.evidence import Evidence, ObservationKey
from tessera.schemas.gene import GeneMechanism
from tessera.schemas.variant import CanonicalVariant

CONFIG_DIR = Path(__file__).resolve().parents[1] / "config"


@pytest.fixture(scope="session")
def cfg() -> ScoringConfig:
    return ScoringConfig.load(CONFIG_DIR)


@pytest.fixture
def mechanism_hi() -> GeneMechanism:
    return GeneMechanism(
        gene_symbol="SETD1A",
        mechanism=Mechanism.HAPLOINSUFFICIENCY,
        confidence=MechanismConfidence.HIGH,
    )


def make_variant(
    vrs_id: str = "ga4gh:VA.1",
    consequence_class: ConsequenceClass = ConsequenceClass.MISSENSE,
    **kw,
) -> CanonicalVariant:
    base = dict(
        vrs_id=vrs_id,
        chrom="16",
        pos=30_976_000,
        ref="G",
        alt="A",
        gene_symbol="SETD1A",
        hgvs_p="p.Gln191His",
        consequence_class=consequence_class,
    )
    base.update(kw)
    return CanonicalVariant(**base)


def make_annotations(
    vrs_id: str = "ga4gh:VA.1",
    *,
    alphamissense: float | None = 0.85,
    esm1b: float | None = -9.0,
    domain: str | None = "annotated_domain",
    rmc: float | None = 0.45,
    avi: float | None = 14.0,
    af: float | None = 3e-6,
    clinvar_class: str | None = "likely_pathogenic",
    stars: float | None = 2,
    pip: float | None = None,
    efficiency: float | None = 0.62,
    designs: float | None = 4,
    spliceai: float | None = 0.05,
) -> VariantAnnotations:
    a = VariantAnnotations(vrs_id=vrs_id)
    if alphamissense is not None:
        a.protein_effect.alphamissense_score = Quantity.present(alphamissense)
    if esm1b is not None:
        a.protein_effect.esm1b_llr = Quantity.present(esm1b)
    if rmc is not None:
        a.protein_effect.regional_missense_constraint = Quantity.present(rmc)
    if domain is not None:
        a.protein_context.domain_context = Quantity.category(domain)
    if avi is not None:
        # avi is given as AVI_PHRED (0 to ~80; 20 = top 1% genome-wide).
        a.alphagenome.avi_score = Quantity.present(avi, source_release="avi-2025.1")
        a.alphagenome.avi_raw = Quantity.present(avi / 20.0)
        a.alphagenome.avi_regulatory_shap = Quantity.present(avi / 30.0)
        a.alphagenome.avi_table_release = "avi-2025.1"
        a.alphagenome.avi_variant_in_table = True
    if af is not None:
        a.frequency.global_af = Quantity.present(af)
    if clinvar_class is not None:
        a.clinvar.germline_classification = clinvar_class
    if stars is not None:
        a.clinvar.review_stars = Quantity.present(stars)
    if pip is not None:
        a.gwas.posterior_inclusion_probability = Quantity.present(pip)
    if spliceai is not None:
        a.splice.spliceai_max_delta = Quantity.present(spliceai)
        a.splice.distance_to_nearest_junction = Quantity.present(120)
    if efficiency is not None:
        a.editability.design_tool = "PRIDICT2.0"
        a.editability.design_tool_version = "2.0"
        a.editability.predicted_efficiency = Quantity.present(efficiency)
        a.editability.pegrna_count = Quantity.present(designs or 1)
        a.editability.bystander_risk = Quantity.present(0.1)
        a.editability.off_target_risk = Quantity.present(0.05)
        a.editability.indel_length = Quantity.present(1)
        a.editability.pam_available = True
        a.editability.host_line = "WTC11"
        a.editability.host_line_genotype_verified = True
        a.editability.host_line_ref_allele_confirmed = True
    a.consequence.position_in_transcript = Quantity.present(0.35)
    a.consequence.nmd_predicted = Quantity.category(True)
    a.protein_context.conservation_phylop = Quantity.present(6.2)
    return a


def make_evidence(
    vrs_id: str | None = "ga4gh:VA.1",
    *,
    evidence_id: str = "ev_1",
    level: EvidenceLevel = EvidenceLevel.VARIANT_DIRECT,
    directness: str = "variant_direct_relevant_human_cell",
    model_system: str = "exact_cell_type",
    match: DiseaseMatchType = DiseaseMatchType.EXACT,
    proband: str | None = "P1",
    cohort: str = "cohort_a",
    cohort_size: int = 4000,
    source_type: str = "publication",
    source_id: str = "PMID:1",
    verified: bool = True,
) -> Evidence:
    return Evidence(
        evidence_id=evidence_id,
        vrs_id=vrs_id,
        gene_symbol="SETD1A",
        evidence_type="functional_experiment",
        evidence_level=level,
        directness_key=directness,
        model_system_key=model_system,
        disease_match_type=match,
        observation=ObservationKey(
            first_report_source_id=source_id,
            cohort_id=cohort,
            proband_id=proband,
            proband_resolution="proband_level" if proband else "cohort_only",
        ),
        cohort_size=cohort_size,
        source_type=source_type,
        source_id=source_id,
        verification_status=(
            VerificationStatus.AGENT_VERIFIED if verified else VerificationStatus.PENDING
        ),
    )
