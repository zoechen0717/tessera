"""screening_policy_v0.1: priority is a host rule over model answers (no network)."""

from tessera.screening import PaperScreen, priority, rule_fields


def _doc(**kw):
    d = {"pub_types": ["research-article"], "is_open_access": True, "in_pmc": True, "has_suppl": True,
         "abstract": "x", "litvar_variants": []}
    d.update(kw)
    d.update(rule_fields(d))
    return d


def _s(**kw):
    base = dict(pmid="1", about_human_gene="yes", variant_reporting="none", variant_class="rare_coding_or_splice",
                allele_functional_experiment="no",
                systems=["human_cells"], case_or_cohort_study="no", is_review="no", evidence_quote="q", rationale="r")
    base.update(kw)
    return PaperScreen(**base)


def test_allele_functional_experiment_is_p1():
    assert priority(_doc(), _s(allele_functional_experiment="yes"))[0] == "P1"


def test_unnamed_variants_in_cohort_is_p1_because_supplements_may_list_them():
    assert priority(_doc(), _s(variant_reporting="variants_reported_not_named", case_or_cohort_study="yes"))[0] == "P1"


def test_litvar_named_variant_overrides_llm_none():
    assert priority(_doc(litvar_variants=["p.R913C"]), _s())[0] == "P2"
    assert priority(_doc(litvar_variants=["rs4889603"]), _s())[0] == "P3"  # rsIDs alone are not named alleles


def test_common_snp_and_somatic_papers_are_not_prioritised():
    assert priority(_doc(), _s(variant_class="common_snp_association", case_or_cohort_study="yes",
                               variant_reporting="specific_alleles_named"))[0] == "P4"
    assert priority(_doc(), _s(variant_class="somatic_cancer", variant_reporting="specific_alleles_named"))[0] == "P4"
    assert priority(_doc(), _s(variant_class="engineered_mutation", allele_functional_experiment="yes"))[0] == "P3"


def test_litvar_only_link_needs_full_text_confirmation():
    body_only = dict(litvar_variants=["p.Y1052F"], body_alleles=[], body_gene_mentions=1)
    assert priority(_doc(**body_only), None)[0] == "P4"                       # OA, scan found nothing
    assert priority(_doc(is_open_access=False, **body_only), None)[0] == "P3"  # cannot confirm
    assert priority(_doc(body_alleles=["c.4582-2delAG"]), None)[0] == "P2"


def test_reviews_and_off_target_papers_are_not_deep_read():
    assert priority(_doc(pub_types=["Review"]), _s(allele_functional_experiment="yes"))[0] == "P4"
    assert priority(_doc(), _s(about_human_gene="no", allele_functional_experiment="yes"))[0] == "P5"


def test_access_rules():
    assert _doc(is_open_access=False)["access_status"] == "abstract_only"
    assert _doc(is_open_access=False)["supplement_status"] == "supplement_unavailable"
    assert _doc(abstract="", is_open_access=False)["access_status"] == "metadata_only"
