"""Resolved disease. Synonym expansion is bounded, never free-form."""

from __future__ import annotations

from datetime import datetime

from pydantic import Field

from tessera.schemas.common import Provenance, Strict


class ResolvedDisease(Strict):
    """Spec 7.2.

    `approved_synonyms` is the complete set an agent may search with. Anything
    outside it is an unconstrained disease expansion, which silently widens the
    query and produces evidence for a different disease than the one asked about.
    `related_not_identical` exists so that such concepts can be recorded without
    being treated as the query disease.
    """

    input: str
    canonical_name: str
    mondo_id: str | None = None
    efo_id: str | None = None
    mesh_terms: list[str] = Field(default_factory=list)
    umls_cuis: list[str] = Field(default_factory=list)
    approved_synonyms: list[str] = Field(default_factory=list)
    abbreviations: list[str] = Field(default_factory=list)
    historical_synonyms: list[str] = Field(default_factory=list)
    related_not_identical: list[str] = Field(default_factory=list)
    resolved_at: datetime | None = None
    provenance: list[Provenance] = Field(default_factory=list)

    def search_terms(self) -> list[str]:
        out = [self.canonical_name, *self.approved_synonyms, *self.abbreviations]
        seen, uniq = set(), []
        for term in out:
            key = term.strip().lower()
            if key and key not in seen:
                seen.add(key)
                uniq.append(term)
        return uniq
