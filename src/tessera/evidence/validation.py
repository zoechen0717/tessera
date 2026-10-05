"""Host-side mechanical validation and effective-state resolution.

Mechanical checks (locator resolves in the pinned snapshot, identity link is
justified, scope permits the link) run on every link regardless of what a model
or importer claimed. They cannot prove semantic entailment (R-11); they can
only reject. A human rejection overrides automatic acceptance (ARCH §8).
"""

from __future__ import annotations

import json

from tessera.schemas.base import ValidationState, canonical_json, sha256_hex
from tessera.schemas.evidence import (
    CheckResult,
    Claim,
    ClaimScope,
    EvidenceKind,
    EvidenceLink,
    ValidationDecision,
)
from tessera.schemas.identity import IdentityDecision, IdentityStatus
from tessera.schemas.source import JsonFieldLocator, SourceSnapshot, TextSpanLocator

POLICY_VERSION = "host-mechanical-1"


def resolve_pointer(doc, pointer: str):
    """RFC 6901."""
    if pointer == "":
        return doc
    if not pointer.startswith("/"):
        raise KeyError(pointer)
    cur = doc
    for raw in pointer[1:].split("/"):
        tok = raw.replace("~1", "/").replace("~0", "~")
        cur = cur[int(tok)] if isinstance(cur, list) else cur[tok]
    return cur


def check_locator(loc, snapshots: dict[str, SourceSnapshot], contents: dict[str, bytes]) -> str | None:
    """Return None if the locator resolves exactly, else a reason."""
    snap = snapshots.get(loc.snapshot_id)
    raw = contents.get(loc.snapshot_id)
    if snap is None or raw is None:
        return f"snapshot {loc.snapshot_id} missing"
    if sha256_hex(raw) != snap.raw_content_digest:
        return f"snapshot {loc.snapshot_id} content altered"
    if isinstance(loc, JsonFieldLocator):
        if snap.media_type != "application/json":
            return "json_field locator on non-JSON snapshot"
        try:
            value = resolve_pointer(json.loads(raw), loc.pointer)
        except (KeyError, IndexError, ValueError):
            return f"pointer {loc.pointer} does not resolve"
        if sha256_hex(canonical_json(value)) != loc.expected_value_digest:
            return f"value at {loc.pointer} does not match expected digest"
        return None
    if isinstance(loc, TextSpanLocator):
        text = raw.decode("utf-8")
        if sha256_hex(text) != loc.text_digest:
            return "text digest mismatch"
        if loc.end > len(text) or text[loc.start : loc.end] != loc.excerpt:
            return f"excerpt not found at offsets {loc.start}:{loc.end}"
        return None
    return "unknown locator kind"


def mechanical_checks(
    link: EvidenceLink,
    claim: Claim,
    identity_decisions: dict[str, IdentityDecision],
    snapshots: dict[str, SourceSnapshot],
    contents: dict[str, bytes],
) -> tuple[dict[str, CheckResult], list[str]]:
    checks: dict[str, CheckResult] = {}
    reasons: list[str] = []

    loc_failures = [r for r in (check_locator(l, snapshots, contents) for l in claim.locators) if r]
    checks["citation_location"] = CheckResult.FAIL if loc_failures else CheckResult.PASS
    reasons += loc_failures

    if link.identity_match.is_exact:
        problems: list[str] = []
        # INV-04: only variant-scoped claims can be exact-allele evidence.
        if claim.scope is not ClaimScope.VARIANT:
            problems.append(f"exact link from a {claim.scope.value}-scoped claim")
        # INV-02: the identity decision must resolve to this very allele.
        decs = [identity_decisions.get(d) for d in link.identity_decision_ids]
        if not decs or any(
            d is None or d.status is not IdentityStatus.RESOLVED or d.variant_ids != [link.variant_id]
            for d in decs
        ):
            problems.append("no resolved identity decision for this allele")
        # INV-02: the source must report this allele, and for an experiment the
        # tested allele -- not merely mention it elsewhere in the paper.
        if link.mention_id not in claim.source_reported_mention_ids:
            problems.append("allele not among the claim's source-reported variants")
        if claim.evidence_kind is EvidenceKind.FUNCTIONAL_EXPERIMENT:
            assay = claim.functional_assay
            if assay is None or link.mention_id not in assay.tested_allele_mention_ids:
                problems.append("allele is not the tested allele of this experiment")
        checks["identity"] = CheckResult.FAIL if problems else CheckResult.PASS
        reasons += problems
    else:
        checks["identity"] = CheckResult.NOT_APPLICABLE

    return checks, reasons


def host_decision(
    link: EvidenceLink,
    claim: Claim,
    identity_decisions: dict[str, IdentityDecision],
    snapshots: dict[str, SourceSnapshot],
    contents: dict[str, bytes],
    sequence: int,
    out_of_scope: bool = False,
) -> ValidationDecision:
    if out_of_scope:
        # The allele resolved but lies outside the declared consequence scope:
        # the claim is kept for audit, the link cannot contribute.
        checks, reasons = {"scope": CheckResult.FAIL}, ["allele outside declared consequence scope"]
    else:
        checks, reasons = mechanical_checks(link, claim, identity_decisions, snapshots, contents)
    failed = any(v is CheckResult.FAIL for v in checks.values())
    return ValidationDecision(
        decision_id=f"vd_host_{link.link_id}_{claim.revision}",
        link_id=link.link_id,
        claim_id=claim.claim_id,
        claim_revision=claim.revision,
        policy_version=POLICY_VERSION,
        checks=checks,
        reasons=reasons,
        reviewer_type="host",
        reviewer_id="host",
        effective_state=ValidationState.REJECTED if failed else ValidationState.MECHANICALLY_VALID,
        sequence=sequence,
    )


def effective_state(decisions: list[ValidationDecision]) -> ValidationState:
    """Resolve the state of one link revision from its append-only decisions.

    - Any host mechanical rejection is final for that revision.
    - A human rejection overrides any acceptance.
    - Otherwise the latest human decision, else the latest automatic one.
    - Mechanical validity alone is never acceptance.
    """
    if not decisions:
        return ValidationState.DRAFT
    ordered = sorted(decisions, key=lambda d: d.sequence)
    if any(d.reviewer_type == "host" and d.effective_state is ValidationState.REJECTED for d in ordered):
        return ValidationState.REJECTED
    human = [d for d in ordered if d.reviewer_type == "human"]
    if any(d.human_review == "rejected" or d.effective_state is ValidationState.REJECTED for d in human):
        return ValidationState.REJECTED
    if human:
        return human[-1].effective_state
    non_host = [d for d in ordered if d.reviewer_type != "host"]
    if not non_host:
        return ValidationState.MECHANICALLY_VALID
    state = non_host[-1].effective_state
    if state is ValidationState.ACCEPTED_HUMAN:
        # Only a human reviewer can produce human acceptance.
        return ValidationState.QUARANTINED
    return state
