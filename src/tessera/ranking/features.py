"""Per-allele feature derivation for nomination_v0.1.

Pure per allele: nothing here reads other candidates, so adding an allele to a
run cannot change another allele's features. Only accepted links are read
(INV-09); only exact-allele links can characterize an allele or set O (INV-04);
the gene mechanism reaches K only through the concordance table (INV-17).
"""

from __future__ import annotations

from dataclasses import dataclass

from tessera.ranking.policy import NominationPolicy
from tessera.schemas.annotation import GeneMechanismAssessment, VariantAnnotation
from tessera.schemas.evidence import Claim, EvidenceKind, EvidenceLink, OutcomeType
from tessera.schemas.features import (
    Channel,
    Characterization,
    Concordance,
    FeatureRow,
    ObservationTier,
)
from tessera.schemas.identity import CanonicalVariant
from tessera.schemas.observed import Observed


@dataclass(frozen=True)
class AcceptedEvidence:
    link: EvidenceLink
    claim: Claim


def mechanism_concordance(
    policy: NominationPolicy,
    mechanism: GeneMechanismAssessment,
    stratum: str,
    annotation: VariantAnnotation | None,
) -> tuple[Concordance, str]:
    mech = mechanism.mechanism
    rule = policy.subclass[mech].get(stratum)
    if rule is None:
        c = policy.stratum_level[mech][stratum]
        return c, f"{mech.value}/{stratum}: class-level only ({c.value})"
    feature = getattr(annotation, rule.feature, None) if annotation else None
    if feature is None or not feature.is_present:
        return rule.missing, f"{mech.value}/{stratum}: {rule.feature} not available"
    key = str(feature.value).lower()
    c = rule.values.get(key, rule.missing)
    return c, f"{mech.value}/{stratum}: {rule.feature}={key}"


def derive_features(
    variant: CanonicalVariant,
    roles: set[str],
    annotation: VariantAnnotation | None,
    accepted: list[AcceptedEvidence],
    mechanism: GeneMechanismAssessment,
    observation_search_status: str,
    policy: NominationPolicy,
    evidence_digest: str,
) -> FeatureRow:
    flags: list[str] = []
    mine = [e for e in accepted if e.link.variant_id == variant.variant_id]
    exact = [e for e in mine if e.link.identity_match.is_exact]
    contributing: dict[str, list[str]] = {}

    # --- stratum
    effect = annotation.selected_effect() if annotation else None
    if annotation is None:
        flags.append("no_annotation_record")
    elif effect is None:
        flags.append("no_single_selected_transcript")
    terms = list(effect.consequences) if effect else []
    stratum = policy.stratum_for(terms)

    # --- characterization (DEC-02)
    tested = [
        e
        for e in exact
        if e.claim.evidence_kind is EvidenceKind.FUNCTIONAL_EXPERIMENT
        and e.claim.functional_assay is not None
        and e.claim.functional_assay.perturbation_type.installs_allele
    ]
    if any(e.claim.functional_assay.outcome_type is OutcomeType.ALTERED for e in tested):
        characterization = Characterization.CHARACTERIZED_ALTERED
    elif tested:
        characterization = Characterization.CHARACTERIZED_OTHER
    else:
        characterization = Characterization.UNCHARACTERIZED
    if tested:
        contributing["characterization"] = sorted(e.link.link_id for e in tested)
    if any(
        e.claim.functional_assay and not e.claim.functional_assay.perturbation_type.installs_allele
        for e in exact
        if e.claim.evidence_kind is EvidenceKind.FUNCTIONAL_EXPERIMENT
    ):
        flags.append("carrier_or_non_installation_experiment_only")

    # --- O (DEC-07)
    o = ObservationTier.NONE
    o_links: list[str] = []
    for e in exact:
        if e.claim.evidence_kind not in (EvidenceKind.PATIENT_OBSERVATION, EvidenceKind.ASSOCIATION):
            continue
        assoc = e.claim.association
        if assoc is None or not assoc.cases_with_allele.is_present or assoc.cases_with_allele.value < 1:
            continue
        if e.link.disease_relationship.is_query_disease:
            tier = (
                ObservationTier.CASE_CONTROLLED
                if assoc.has_control_denominator
                else ObservationTier.CASE_UNCONTROLLED
            )
        elif e.link.disease_relationship.value == "related":
            tier = ObservationTier.RELATED_CONDITION
        else:
            continue
        o_links.append(e.link.link_id)
        if tier.order > o.order:
            o = tier
    if o_links:
        contributing["observation"] = sorted(o_links)

    # Non-exact accepted links stay visible as context and do nothing else.
    context = [e.link.link_id for e in mine if not e.link.identity_match.is_exact]
    if context:
        contributing["context_only"] = sorted(context)

    # --- K
    k, k_basis = mechanism_concordance(policy, mechanism, stratum, annotation)

    # --- P
    predictor = policy.predictors.get(stratum)
    p: Observed[float] | None = None
    if predictor:
        p = (annotation.predictions.get(predictor) if annotation else None) or Observed[float].not_searched()
        if not p.is_present:
            flags.append(f"p_{p.status.value}")

    # --- channel and roles
    out_roles = set(roles)
    if out_roles & policy.control_roles:
        channel = Channel.CONTROL
    elif characterization is not Characterization.UNCHARACTERIZED and policy.characterized_policy == "route_to_control":
        channel = Channel.CONTROL
    else:
        channel = Channel.DISCOVERY
    if characterization is Characterization.CHARACTERIZED_ALTERED:
        out_roles.add("positive_control_candidate")
    elif characterization is Characterization.CHARACTERIZED_OTHER:
        out_roles.add("characterized_other")
    if not out_roles:
        out_roles.add("discovery")

    return FeatureRow(
        variant_id=variant.variant_id,
        label=variant.label,
        policy_id=policy.policy_id,
        policy_digest=policy.digest,
        policy_status=policy.status,
        evidence_digest=evidence_digest,
        primary_channel=channel,
        roles=sorted(out_roles),
        stratum=stratum,
        consequence_terms=terms,
        characterization=characterization,
        k=k,
        k_basis=k_basis,
        o=o,
        o_search_status=observation_search_status,
        p_predictor=predictor,
        p=p,
        contributing_link_ids=contributing,
        flags=sorted(flags),
    )
