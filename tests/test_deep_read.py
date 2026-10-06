"""Passage selection for deep reading (no network, no LLM)."""

import re

from tessera.deep_read import article_units, select_passages

XML = """<article><body>
<sec><title>Results</title>
<p>We studied <italic>SETD1A</italic> in neurons carrying c.4596_4597insG (p.Leu1533fs).</p>
<p>Unrelated paragraph about culture conditions.</p>
<p>The c.4596_4597insG line showed reduced synaptic density.</p>
<table-wrap><table><tr><td>KMT2D</td><td>c.100A&gt;G</td></tr><tr><td>SETD 1A</td><td>c.4582-2delAG</td></tr></table></table-wrap>
<fig><caption><p>Western blot of SETD1A protein.</p></caption></fig>
</sec></body></article>"""

ALIAS = re.compile(r"\bSETD\s?1A\b", re.I)


def test_units_split_paragraphs_rows_captions():
    us = article_units(XML)
    assert "We studied SETD1A in neurons carrying c.4596_4597insG (p.Leu1533fs)." in us
    assert any(u.startswith("SETD 1A | c.4582-2delAG") for u in us)


def test_selection_keeps_gene_units_and_cooccurring_variant_units_only():
    us = article_units(XML)
    picked = [us[i] for i in select_passages(us, ALIAS, max_units=10, max_chars=10_000)]
    assert any("reduced synaptic density" in u for u in picked)       # same allele as a gene unit
    assert not any("culture conditions" in u for u in picked)
    assert not any(u.startswith("KMT2D") for u in picked)              # another gene's table row


def test_selection_respects_caps():
    us = article_units(XML)
    assert len(select_passages(us, ALIAS, max_units=1, max_chars=10_000)) == 1
