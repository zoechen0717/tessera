"""Deterministic annotations (spec 9). Every numeric field is a Quantity.

This is the reason the module is verbose: a bare `float | None` cannot say
whether an absent AlphaMissense score means "not a coding variant", "table
lookup failed", or "we never looked". Those three must score differently.
"""

from __future__ import annotations

from datetime import datetime

from pydantic import Field

from tessera.enums import LabelLeakage, Missingness, VariantStatus
from tessera.schemas.common import Quantity, SourceSnapshot, Strict

_ns = Quantity.not_searched


class ConsequenceAnnotation(Strict):
    consequence: str | None = None
    severity: str | None = None
    coding_status: str | None = None
    exon: str | None = None
    intron: str | None = None
    splice_region: bool | None = None
    canonical_transcript_effect: str | None = None
    position_in_transcript: Quantity = Field(default_factory=_ns)
    nmd_predicted: Quantity = Field(default_factory=_ns)


class PopulationFrequency(Strict):
    """Spec 9.2. Not collapsed prematurely -- max-population AF and homozygote
    counts carry information that global AF hides."""

    global_af: Quantity = Field(default_factory=_ns)
    max_population_af: Quantity = Field(default_factory=_ns)
    homozygote_count: Quantity = Field(default_factory=_ns)
    hemizygote_count: Quantity = Field(default_factory=_ns)
    ancestry_breakdown: dict[str, float] = Field(default_factory=dict)
    gnomad_version: str | None = None


class ClinVarAnnotation(Strict):
    """Spec 9.3. Deliberately NOT a single clinvar_count.

    Multiple weak submissions are not automatically stronger than an expert-panel
    review, and condition-specific classifications are what decide whether this is
    evidence for the QUERY disease or for a different one.
    """

    variation_id: str | None = None
    vcv_id: str | None = None
    rcv_ids: list[str] = Field(default_factory=list)
    germline_classification: str | None = None
    condition_specific_classifications: dict[str, str] = Field(default_factory=dict)
    review_status: str | None = None
    review_stars: Quantity = Field(default_factory=_ns)
    submission_count: int | None = None
    independent_submitters: int | None = None
    disease_matched_submitters: int | None = None
    conflict_status: str | None = None
    last_evaluated: str | None = None
    citations: list[str] = Field(default_factory=list)
    clinvar_release: str | None = None


class GwasAnnotation(Strict):
    associated_trait: str | None = None
    lead_variant: str | None = None
    p_value: Quantity = Field(default_factory=_ns)
    log10_p: Quantity = Field(default_factory=_ns)
    beta: Quantity = Field(default_factory=_ns)
    odds_ratio: Quantity = Field(default_factory=_ns)
    credible_set_member: bool | None = None
    posterior_inclusion_probability: Quantity = Field(default_factory=_ns)
    locus_id: str | None = None
    mapped_gene: str | None = None
    gene_mapping_method: str | None = None


class AlphaGenomeAnnotation(Strict):
    """Spec 9.5.1. Local lookup table, not an API.

    AVI is NOT a pure AlphaGenome output. It is a supervised 18-feature model whose
    inputs are: aggregated AlphaGenome ISM scores across modalities, AlphaMissense,
    three VEP protein-termination annotations, two conservation scores (PhastCons
    470-way, Zoonomia Cactus 241-way), and two indel indicators. It is trained to
    separate variants with filtering allele frequency above vs below 0.1%.

    Two consequences the scorer has to respect:

    1. Using `avi_score` in the same branch as `alphamissense_score` double-counts
       AlphaMissense. Coding branches use `avi_regulatory_shap` instead -- the
       SHAP-isolated AlphaGenome-modality contribution, i.e. the increment AVI
       carries that AlphaMissense does not.

    2. THE BULK DOWNLOAD IS SNV-ONLY. The paper scores >100M observed indels, but
       those are not in the downloadable tables -- only in the Atlas portal/API.
       The distributed assets are all named SNV:

         AVI SNV scores (raw + PHRED), tabix          ~88.5 GB
         AVI SNV feature importance (SHAP), tabix    ~283.9 GB
         AlphaGenome SNV merged splicing, tabix       ~20.6 GB

       The Ensembl VEP AVI plugin makes this explicit in code: it returns empty
       unless start == end and both ref and alt match /^[ACGT]$/.

       For a prime-editing panel this is not an edge case. Every non-SNV allele
       misses the table, so `avi_variant_in_table` is False and the term is
       NOT_APPLICABLE -- excluded with weights renormalized, never a low score.
    """

    avi_score: Quantity = Field(
        default_factory=_ns,
        description="AVI_PHRED, 0 to ~80. PHRED 20 = top 1% most deleterious "
        "substitutions genome-wide. Column 6 of the tabix TSV.",
    )
    avi_raw: Quantity = Field(
        default_factory=_ns,
        description="AVI_RAW, approx -1.3 to +4.0. Column 5 of the tabix TSV. The "
        "scale the SHAP contributions sum to.",
    )
    avi_table_release: str | None = Field(
        default=None, description="Required in the run manifest when avi_score is present"
    )
    avi_regulatory_shap: Quantity = Field(
        default_factory=_ns,
        description="Sum of the AlphaGenome-modality SHAP contributions (splicing, "
        "ATAC, DNase, ChIP-TF, ChIP-histone, CAGE, PRO-cap, RNA-seq, polyadenylation, "
        "contact maps), EXCLUDING the AlphaMissense and protein-termination "
        "contributions. Signed, on the raw pre-PHRED AVI scale.",
    )
    avi_coding_shap: Quantity = Field(
        default_factory=_ns,
        description="AlphaMissense + protein-termination SHAP contributions. Stored "
        "for audit and for checking redundancy against alphamissense_score; NOT used "
        "in any score, because that value is already carried by ProteinEffectAnnotation.",
    )
    avi_variant_in_table: bool | None = Field(
        default=None,
        description="Whether this variant was found in the table at all. False for an "
        "unobserved indel -- a coverage gap, not a low impact prediction.",
    )
    rna_score: Quantity = Field(default_factory=_ns)
    atac_score: Quantity = Field(default_factory=_ns)
    dnase_score: Quantity = Field(default_factory=_ns)
    tf_binding_score: Quantity = Field(default_factory=_ns)
    histone_score: Quantity = Field(default_factory=_ns)
    splice_score: Quantity = Field(default_factory=_ns)
    contact_score: Quantity = Field(default_factory=_ns)
    strongest_modality: str | None = None
    strongest_context: str | None = None
    modality_relevance: Quantity = Field(
        default_factory=_ns,
        description="How well strongest_context matches the screen cell type",
    )
    affected_gene: str | None = None
    coverage_note: str | None = Field(
        default=None,
        description="e.g. 'non-SNV: absent from the SNV-only bulk table' -- an "
        "uncovered variant is NOT_APPLICABLE, never a low score",
    )


class ProteinEffectAnnotation(Strict):
    """Spec 9.5.2. `label_leakage` is not optional.

    Benchmarking a ClinVar-supervised predictor against a ClinVar-derived answer
    key is circular. Metrics are reported per leakage class, never pooled.
    """

    alphamissense_score: Quantity = Field(default_factory=_ns)
    alphamissense_class: Quantity = Field(default_factory=_ns)
    esm1b_llr: Quantity = Field(default_factory=_ns)
    cadd_phred: Quantity = Field(default_factory=_ns)
    revel_score: Quantity = Field(default_factory=_ns)
    mpc_score: Quantity = Field(default_factory=_ns)
    regional_missense_constraint: Quantity = Field(default_factory=_ns)
    label_leakage: dict[str, LabelLeakage] = Field(
        default_factory=lambda: {
            # Not trained on ClinVar labels, but benchmarked on DMS -- only weakly
            # independent of a DMS answer key.
            "alphamissense_score": LabelLeakage.DMS,
            # Unsupervised: independent of both ClinVar and DMS keys.
            "esm1b_llr": LabelLeakage.NONE,
            # Supervised on curated pathogenic sets: a ClinVar-keyed evaluation of
            # these is circular.
            "cadd_phred": LabelLeakage.CLINVAR,
            "revel_score": LabelLeakage.CLINVAR,
            "mpc_score": LabelLeakage.CLINVAR,
            # AVI: trained on a population-frequency proxy label (FAF 0.1%), so
            # relatively fair against a ClinVar key. But it takes AlphaMissense as
            # an input feature, so the DMS caveat propagates, and it is NOT
            # independent of the rarity term in D, which also uses allele frequency.
            "alphagenome_avi": LabelLeakage.ALLELE_FREQUENCY,
            "avi_regulatory_shap": LabelLeakage.ALLELE_FREQUENCY,
        }
    )


class SpliceAnnotation(Strict):
    spliceai_ds_ag: Quantity = Field(default_factory=_ns)
    spliceai_ds_al: Quantity = Field(default_factory=_ns)
    spliceai_ds_dg: Quantity = Field(default_factory=_ns)
    spliceai_ds_dl: Quantity = Field(default_factory=_ns)
    spliceai_max_delta: Quantity = Field(default_factory=_ns)
    pangolin_score: Quantity = Field(default_factory=_ns)
    distance_to_nearest_junction: Quantity = Field(default_factory=_ns)
    predicted_event: str | None = None


class ProteinContextAnnotation(Strict):
    uniprot_id: str | None = None
    residue: str | None = None
    domain: str | None = None
    domain_start: int | None = None
    domain_end: int | None = None
    motif: str | None = None
    active_site: bool | None = None
    binding_site: bool | None = None
    disordered_region: bool | None = None
    domain_context: Quantity = Field(default_factory=_ns)
    conservation_phylop: Quantity = Field(default_factory=_ns)
    alphafold_plddt: Quantity = Field(default_factory=_ns)
    structural_context: str | None = None


class MaveAnnotation(Strict):
    """Spec 9.6b. A measurement, not a prediction -- routes into F, not M."""

    urn: str | None = None
    assay_type: str | None = None
    model_system: str | None = None
    cell_type: str | None = None
    score: Quantity = Field(default_factory=_ns)
    score_normalized: Quantity = Field(default_factory=_ns)
    functional_class: Quantity = Field(default_factory=_ns)
    n_replicates: int | None = None
    score_set_version: str | None = None
    retrieved_at: datetime | None = None


class EditabilityAnnotation(Strict):
    """Spec 17.1. The design tool must be named: the E anchor maps THAT tool's
    output distribution, so swapping tools requires re-anchoring."""

    design_tool: str | None = None
    design_tool_version: str | None = None

    pam_available: bool | None = None
    pegrna_count: Quantity = Field(default_factory=_ns)
    predicted_efficiency: Quantity = Field(default_factory=_ns)
    bystander_risk: Quantity = Field(default_factory=_ns)
    off_target_risk: Quantity = Field(default_factory=_ns)
    indel_length: Quantity = Field(default_factory=_ns)
    mutation_type: str | None = None

    host_line: str | None = None
    host_line_genotype_verified: bool = False
    host_line_ref_allele_confirmed: bool | None = None
    host_line_conflict: str | None = Field(
        default=None, description="null | heterozygous | homozygous_alt | cnv"
    )

    variant_status: VariantStatus = VariantStatus.EDITABILITY_UNKNOWN

    def resolve_status(self) -> VariantStatus:
        """Spec 17.1.2. Host-line conflict is a hard filter, not a penalty.

        A biologically excellent variant at a position where the host line is not
        reference cannot be built, and nothing in a purely biological ranking would
        surface that. Discovering it after a panel is ordered is expensive.
        """
        if self.host_line_conflict:
            return VariantStatus.NOT_EDITABLE_IN_HOST_LINE
        if self.pam_available is False:
            return VariantStatus.NOT_EDITABLE_NO_PAM
        if not self.host_line_genotype_verified:
            return VariantStatus.EDITABILITY_UNKNOWN
        if self.pegrna_count.is_present and (self.pegrna_count.value or 0) <= 0:
            return VariantStatus.NOT_EDITABLE_NO_PAM
        return VariantStatus.EDITABLE


class GeneConstraintAnnotation(Strict):
    gnomad_loeuf: Quantity = Field(default_factory=_ns)
    gnomad_pli: Quantity = Field(default_factory=_ns)
    gnomad_mis_z: Quantity = Field(default_factory=_ns)
    gnomad_syn_z: Quantity = Field(default_factory=_ns)
    clingen_haploinsufficiency_score: Quantity = Field(default_factory=_ns)
    clingen_triplosensitivity_score: Quantity = Field(default_factory=_ns)
    gnomad_version: str | None = None


class VariantAnnotations(Strict):
    """Everything deterministic known about one variant."""

    vrs_id: str
    consequence: ConsequenceAnnotation = Field(default_factory=ConsequenceAnnotation)
    frequency: PopulationFrequency = Field(default_factory=PopulationFrequency)
    clinvar: ClinVarAnnotation = Field(default_factory=ClinVarAnnotation)
    gwas: GwasAnnotation = Field(default_factory=GwasAnnotation)
    alphagenome: AlphaGenomeAnnotation = Field(default_factory=AlphaGenomeAnnotation)
    protein_effect: ProteinEffectAnnotation = Field(default_factory=ProteinEffectAnnotation)
    splice: SpliceAnnotation = Field(default_factory=SpliceAnnotation)
    protein_context: ProteinContextAnnotation = Field(default_factory=ProteinContextAnnotation)
    mave: MaveAnnotation = Field(default_factory=MaveAnnotation)
    editability: EditabilityAnnotation = Field(default_factory=EditabilityAnnotation)
    gene_constraint: GeneConstraintAnnotation = Field(default_factory=GeneConstraintAnnotation)

    annotation_version: str | None = None
    snapshots: list[SourceSnapshot] = Field(default_factory=list)
    retrieved_at: datetime | None = None

    def raw_input(self, name: str) -> Quantity:
        """Flat lookup by the raw-input name used in anchors.yaml and the M branches.

        One dispatch table so that a config naming a raw input nothing produces is
        an immediate, locatable failure rather than a silently missing term.
        """
        table: dict[str, Quantity] = {
            # mechanism, coding
            "alphamissense_score": self.protein_effect.alphamissense_score,
            "alphamissense_class": self.protein_effect.alphamissense_class,
            "esm1b_llr": self.protein_effect.esm1b_llr,
            "cadd_phred": self.protein_effect.cadd_phred,
            "revel_score": self.protein_effect.revel_score,
            "mpc_score": self.protein_effect.mpc_score,
            "regional_missense_constraint": self.protein_effect.regional_missense_constraint,
            # mechanism, splice
            "spliceai_max_delta": self.splice.spliceai_max_delta,
            "pangolin_score": self.splice.pangolin_score,
            "distance_to_nearest_junction": self.splice.distance_to_nearest_junction,
            # mechanism, regulatory
            "alphagenome_avi": self.alphagenome.avi_score,
            "avi_raw": self.alphagenome.avi_raw,
            "avi_regulatory_shap": self.alphagenome.avi_regulatory_shap,
            "alphagenome_modality_relevance": self.alphagenome.modality_relevance,
            # mechanism, ptv / position
            "position_in_transcript": self.consequence.position_in_transcript,
            "nmd_predicted": self.consequence.nmd_predicted,
            # protein / conservation
            "protein_domain_context": self.protein_context.domain_context,
            "conservation_phylop": self.protein_context.conservation_phylop,
            # disease
            "gwas_pip": self.gwas.posterior_inclusion_probability,
            "gwas_log10_p": self.gwas.log10_p,
            "clinvar_review_stars": self.clinvar.review_stars,
            "gnomad_af": self.frequency.global_af,
            "gnomad_homozygote_count": self.frequency.homozygote_count,
            # gene-level constraint
            "gnomad_loeuf": self.gene_constraint.gnomad_loeuf,
            "gnomad_pli": self.gene_constraint.gnomad_pli,
            "clingen_haploinsufficiency_score": (
                self.gene_constraint.clingen_haploinsufficiency_score
            ),
            # functional measurement
            "mave_functional_class": self.mave.functional_class,
            "mave_score_normalized": self.mave.score_normalized,
            # editability
            "pegrna_predicted_efficiency": self.editability.predicted_efficiency,
            "pegrna_design_count": self.editability.pegrna_count,
            "bystander_risk": self.editability.bystander_risk,
            "off_target_risk": self.editability.off_target_risk,
            "indel_length": self.editability.indel_length,
        }
        if name not in table:
            raise KeyError(
                f"raw input {name!r} is referenced by config but not produced by any "
                f"annotation adapter. Either add it to VariantAnnotations.raw_input or "
                f"remove it from the config. Available: {sorted(table)}"
            )
        return table[name]

    def blocking_flags(self) -> list[str]:
        """Missingness states that must keep a variant out of a final panel."""
        flags: list[str] = []
        for name in (
            "alphamissense_score",
            "spliceai_max_delta",
            "alphagenome_avi",
            "gnomad_af",
            "pegrna_predicted_efficiency",
        ):
            try:
                q = self.raw_input(name)
            except KeyError:
                continue
            if q.state is Missingness.FAILED_RETRIEVAL:
                flags.append(f"retrieval_failure:{name}")
        return flags
