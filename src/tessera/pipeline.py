"""M0 offline vertical slice: freeze → derive → export, and replay.

    freeze(input_dir)  : mentions → identity → registry; import annotations,
                         snapshots, claims, links, decisions; host mechanical
                         validation; write frozen objects + evidence digest.
    derive(bundle)     : frozen objects → features → stratified ranks → exports.

`run` is freeze followed by derive on the same bundle. `replay` runs derive on
an existing bundle into a new directory and compares derived digests. Replay
imports nothing that can reach a network or a model.
"""

from __future__ import annotations

import csv
import json
import shutil
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

import yaml

import tessera
from tessera import bundle as B
from tessera.evidence.validation import effective_state, host_decision
from tessera.identity.normalize import Capabilities
from tessera.identity.reference import ReferenceBundle
from tessera.identity.registry import build_registry
from tessera.ranking.features import AcceptedEvidence, derive_features
from tessera.ranking.policy import NominationPolicy
from tessera.ranking.rank import rank
from tessera.reporting import render_csv, render_report
from tessera.schemas.annotation import GeneEntity, GeneMechanismAssessment, VariantAnnotation
from tessera.schemas.base import SCHEMA_VERSION, sha256_hex
from tessera.schemas.evidence import Claim, EvidenceLink, ValidationDecision
from tessera.schemas.features import FeatureRow
from tessera.schemas.identity import CanonicalVariant, Exclusion, VariantMention
from tessera.schemas.source import SourceSnapshot


def _code_digest() -> str:
    root = Path(tessera.__file__).parent
    files = sorted(p for p in root.rglob("*.py"))
    return sha256_hex("".join(f"{p.relative_to(root)}\0{sha256_hex(p.read_bytes())}\n" for p in files))


def _read_variants_tsv(path: Path, route: str) -> list[VariantMention]:
    out = []
    with path.open() as fh:
        for row in csv.DictReader((l for l in fh if not l.startswith("##")), delimiter="\t"):
            def opt(k):
                v = (row.get(k) or "").strip()
                return v or None
            pos = opt("pos")
            out.append(
                VariantMention(
                    mention_id=row["mention_id"],
                    raw_text=row["raw_text"],
                    discovery_route=opt("route") or route,
                    reported_assembly=opt("assembly"),
                    reported_chrom=opt("chrom"),
                    reported_pos=int(pos) if pos else None,
                    reported_ref=opt("ref"),
                    reported_alt=opt("alt"),
                    roles=[r for r in (opt("roles") or "").split(",") if r],
                    resolution_method=opt("resolution_method"),
                    resolution_status=opt("resolution_status"),
                    resolution_note=opt("resolution_note"),
                    source_refs=[x for x in (opt("source_refs") or "").split(";") if x],
                )
            )
    return out


def _vcf_variant_id(reg_bundle: ReferenceBundle, rec: dict) -> str:
    from tessera.identity.normalize import normalize

    res = normalize(reg_bundle, rec["chrom"], rec["pos"], rec["ref"], rec["alt"])
    if res.allele is None:
        raise ValueError(f"annotation for unnormalizable allele {rec}: {res.errors}")
    return res.allele.variant_id


# --------------------------------------------------------------------- freeze


def freeze(input_path: Path, out: Path) -> dict:
    input_path = input_path.resolve()
    base = input_path.parent
    cfg = yaml.safe_load(input_path.read_text())
    out.mkdir(parents=True, exist_ok=True)

    ref = ReferenceBundle.load(base / cfg["reference"])
    gene = GeneEntity.model_validate(cfg["gene"])
    policy_path = (base / cfg["policy"]).resolve()
    policy = NominationPolicy.load(policy_path)

    # identity
    mentions = _read_variants_tsv(base / cfg["variants"], route="manual")
    reg = build_registry(mentions, ref, gene.gene_id, Capabilities())

    # annotations, keyed by normalized allele
    annotations: list[VariantAnnotation] = []
    for line in (base / cfg["annotations"]).read_text().splitlines():
        if not line.strip():
            continue
        rec = json.loads(line)
        vid = _vcf_variant_id(ref, rec.pop("allele"))
        if vid in reg.variants:
            annotations.append(VariantAnnotation.model_validate({**rec, "variant_id": vid}))

    # consequence scope: resolved alleles outside it leave the registry but stay
    # visible as exclusions (R-03, SPEC §8.1)
    allowed = set((cfg.get("scope") or {}).get("allowed_consequences") or [])
    if allowed:
        by_vid = {a.variant_id: a for a in annotations}
        for vid in sorted(reg.variants):
            if vid not in by_vid:  # annotation failed: keep, flagged downstream
                continue
            effect = by_vid[vid].selected_effect()
            terms = list(effect.consequences) if effect else []
            if allowed & set(terms):
                continue
            for mid in [m for m, v in reg.mention_to_variant.items() if v == vid]:
                reg.exclusions.append(
                    Exclusion(
                        mention_id=mid,
                        decision_id=reg.mention_to_decision[mid],
                        status="out_of_scope",
                        reason=("consequence " + ",".join(terms) if terms else "no consequence on the selected transcript")
                        + " outside declared scope",
                        variant_id=vid,
                    )
                )
                del reg.mention_to_variant[mid]
            del reg.variants[vid]
            reg.roles.pop(vid, None)
        annotations = [a for a in annotations if a.variant_id in reg.variants]

    # snapshots: copy raw content, verify content addressing
    snap_dir = out / "snapshots"
    snap_dir.mkdir(exist_ok=True)
    snapshots: list[SourceSnapshot] = []
    contents: dict[str, bytes] = {}
    for line in (base / cfg["snapshots"] / "snapshots.jsonl").read_text().splitlines():
        if not line.strip():
            continue
        meta = json.loads(line)
        raw = (base / cfg["snapshots"] / meta.pop("file")).read_bytes()
        snap = SourceSnapshot.model_validate(meta)
        if sha256_hex(raw) != snap.raw_content_digest:
            raise ValueError(f"snapshot {snap.source_record_id}: content does not match its digest")
        contents[snap.snapshot_id] = raw
        (snap_dir / snap.snapshot_id).write_bytes(raw)
        snapshots.append(snap)
    snap_by_id = {s.snapshot_id: s for s in snapshots}
    if any(s.is_synthetic != ref.is_synthetic for s in snapshots):
        raise ValueError("synthetic and real sources cannot be mixed in one run (INV-12)")

    claims = [Claim.model_validate_json(l) for l in (base / cfg["claims"]).read_text().splitlines() if l.strip()]
    claim_by_key = {c.key: c for c in claims}

    # links: resolve the source mention through the identity registry
    links: list[EvidenceLink] = []
    for l in (base / cfg["links"]).read_text().splitlines():
        if not l.strip():
            continue
        rec = json.loads(l)
        mid = rec.get("mention_id")
        rec["variant_id"] = reg.mention_to_variant.get(mid)
        rec["identity_decision_ids"] = [reg.mention_to_decision[mid]] if mid in reg.mention_to_decision else []
        links.append(EvidenceLink.model_validate(rec))

    # validation: imported decisions first, then host mechanical decision per link
    decisions = [
        ValidationDecision.model_validate_json(l)
        for l in (base / cfg["decisions"]).read_text().splitlines()
        if l.strip()
    ]
    id_decisions = {d.decision_id: d for d in reg.decisions}
    out_of_scope = {e.mention_id for e in reg.exclusions if e.status == "out_of_scope"}
    seq = max((d.sequence for d in decisions), default=0)
    for link in links:
        claim = claim_by_key[(link.claim_id, link.claim_revision)]
        seq += 1
        decisions.append(host_decision(link, claim, id_decisions, snap_by_id, contents, seq,
                                       out_of_scope=link.mention_id in out_of_scope))

    mechanism = GeneMechanismAssessment.model_validate(json.loads((base / cfg["mechanism"]).read_text()))
    if mechanism.gene_id != gene.gene_id:
        raise ValueError("mechanism assessment is for a different gene")

    # write frozen objects
    shutil.copy(input_path, out / "input.yaml")
    pol_dir = out / "policy"
    pol_dir.mkdir(exist_ok=True)
    B.write_json(pol_dir / "policy.json", policy.raw)
    B.write_json(
        out / "resolved_scope.json",
        {
            "schema_version": SCHEMA_VERSION,
            "gene": gene.model_dump(mode="json"),
            "disease": cfg["disease"],
            "screen_context": cfg.get("screen_context", {}),
            "reference_bundle_id": ref.bundle_id,
            "assembly": ref.assembly,
            "fasta_sha256": ref.fasta_sha256,
            "is_synthetic": ref.is_synthetic,
            "policy_id": policy.policy_id,
            "policy_digest": policy.digest,
            "requested_panel_size": cfg.get("requested_panel_size"),
        },
    )
    B.write_jsonl(out / "candidate_mentions.jsonl", reg.mentions)
    B.write_jsonl(out / "identity_decisions.jsonl", reg.decisions)
    B.write_jsonl(out / "exclusions.jsonl", reg.exclusions)
    B.write_jsonl(out / "candidate_variants.jsonl", sorted(reg.variants.values(), key=lambda v: v.variant_id))
    B.write_jsonl(
        out / "memberships.jsonl",
        [{"variant_id": v, "roles": sorted(r)} for v, r in sorted(reg.roles.items())],
    )
    B.write_jsonl(out / "annotations.jsonl", sorted(annotations, key=lambda a: a.variant_id))
    B.write_jsonl(out / "source_snapshots.jsonl", sorted(snapshots, key=lambda s: s.snapshot_id))
    B.write_jsonl(out / "claims.jsonl", sorted(claims, key=lambda c: c.key))
    B.write_jsonl(out / "evidence_links.jsonl", sorted(links, key=lambda l: l.link_id))
    B.write_jsonl(out / "validation.jsonl", sorted(decisions, key=lambda d: d.sequence))
    B.write_json(out / "gene_mechanism.json", mechanism)
    B.write_json(out / "observation_search.json", cfg["observation_search"])
    (out / "datasets.jsonl").write_text("")  # GEO is M1b; empty by scope, recorded in manifest

    digests = B.digest_files(out, B.FROZEN_FILES)
    frozen = {
        "evidence_digest": B.combined_digest(digests),
        "frozen_files": digests,
        "policy_digest": policy.digest,
        "is_synthetic": ref.is_synthetic,
    }
    B.write_json(out / "frozen.json", frozen)
    return frozen


# --------------------------------------------------------------------- derive


@dataclass
class Derived:
    features: list[FeatureRow]
    content_digest: str
    artifact_digests: dict[str, str]


def derive(bundle_dir: Path, out: Path, *, mode: str) -> Derived:
    out.mkdir(parents=True, exist_ok=True)
    frozen = B.read_json(bundle_dir / "frozen.json")
    current = B.digest_files(bundle_dir, B.FROZEN_FILES)
    if B.combined_digest(current) != frozen["evidence_digest"]:
        raise ValueError("frozen bundle was modified after freezing (INV-13)")

    pol = B.read_json(bundle_dir / "policy" / "policy.json")
    policy = NominationPolicy.from_dicts(pol["policy"], pol["table"])
    if policy.digest != frozen["policy_digest"]:
        raise ValueError("policy in bundle does not match frozen policy digest")

    variants = B.read_jsonl(bundle_dir / "candidate_variants.jsonl", CanonicalVariant)
    roles = {m["variant_id"]: set(m["roles"]) for m in map(json.loads, (bundle_dir / "memberships.jsonl").read_text().splitlines())}
    annotations = {a.variant_id: a for a in B.read_jsonl(bundle_dir / "annotations.jsonl", VariantAnnotation)}
    claims = {c.key: c for c in B.read_jsonl(bundle_dir / "claims.jsonl", Claim)}
    links = B.read_jsonl(bundle_dir / "evidence_links.jsonl", EvidenceLink)
    decisions = B.read_jsonl(bundle_dir / "validation.jsonl", ValidationDecision)
    mechanism = GeneMechanismAssessment.model_validate(B.read_json(bundle_dir / "gene_mechanism.json"))
    obs_search = B.read_json(bundle_dir / "observation_search.json")

    by_link: dict[tuple[str, int], list[ValidationDecision]] = defaultdict(list)
    for d in decisions:
        by_link[(d.link_id, d.claim_revision)].append(d)
    states = {l.link_id: effective_state(by_link[(l.link_id, l.claim_revision)]) for l in links}
    accepted = [
        AcceptedEvidence(link=l, claim=claims[(l.claim_id, l.claim_revision)])
        for l in links
        if l.variant_id is not None and states[l.link_id].is_accepted
    ]

    features = [
        derive_features(
            v,
            roles.get(v.variant_id, set()),
            annotations.get(v.variant_id),
            accepted,
            mechanism,
            obs_search["status"],
            policy,
            frozen["evidence_digest"],
        )
        for v in sorted(variants, key=lambda v: v.variant_id)
    ]
    coords = {v.variant_id: v.vcf_key for v in variants}
    ranked, strata = rank(features, coords, policy, mechanism)

    # provenance per allele: discovery routes, source references, accepted literature evidence
    decided = {d["mention_id"]: d for d in map(json.loads, (bundle_dir / "identity_decisions.jsonl").read_text().splitlines()) if d}
    prov: dict[str, dict] = defaultdict(lambda: {"routes": set(), "refs": set(), "lit": set()})
    for m in map(json.loads, (bundle_dir / "candidate_mentions.jsonl").read_text().splitlines()):
        d = decided.get(m["mention_id"])
        if d and d["status"] == "resolved":
            p = prov[d["variant_ids"][0]]
            p["routes"].add(m["discovery_route"])
            p["refs"].update(m.get("source_refs") or [])
    for e in accepted:
        c = e.claim
        if c.publication_group_id and c.publication_group_id.startswith("PMID:"):
            a = c.functional_assay
            detail = f"{a.perturbation_type.value}/{a.outcome_type.value}" if a else c.evidence_kind.value
            prov[e.link.variant_id]["lit"].add(f"{c.publication_group_id} {c.evidence_kind.value}:{detail}")
    provenance = {v: {"routes": sorted(p["routes"]), "refs": sorted(p["refs"]), "lit": sorted(p["lit"])}
                  for v, p in prov.items()}
    for vid, a in annotations.items():
        eff = a.selected_effect()
        if eff:
            provenance.setdefault(vid, {"routes": [], "refs": [], "lit": []})
            provenance[vid]["hgvs_c"] = (eff.hgvs_c or "").split(":")[-1]
            provenance[vid]["hgvs_p"] = (eff.hgvs_p or "").split(":")[-1]

    B.write_jsonl(out / "features.jsonl", features)
    feat_by_id = {f.variant_id: f for f in features}
    (out / "ranked_variants.csv").write_text(render_csv(ranked, feat_by_id, provenance), encoding="utf-8")
    B.write_json(out / "strata_report.json", [s.model_dump(mode="json") for s in strata])
    scope = B.read_json(bundle_dir / "resolved_scope.json")
    source_failures = (yaml.safe_load((bundle_dir / "input.yaml").read_text()).get("provenance") or {}).get(
        "source_failures") or []
    exclusions = [json.loads(l) for l in (bundle_dir / "exclusions.jsonl").read_text().splitlines() if l]
    mentions = [json.loads(l) for l in (bundle_dir / "candidate_mentions.jsonl").read_text().splitlines() if l]
    routes = dict(sorted(Counter(m["discovery_route"] for m in mentions).items()))
    scope_links = sorted({d.link_id for d in decisions if d.checks.get("scope") == "fail"})
    report = render_report(
        scope=scope,
        policy=policy,
        mechanism=mechanism,
        features=features,
        ranked=ranked,
        strata=strata,
        exclusions=exclusions,
        n_mentions=len(mentions),
        routes=routes,
        n_scope_links=len(scope_links),
        link_states={k: v.value for k, v in sorted(states.items())},
        obs_search=obs_search,
        source_failures=source_failures,
        provenance=provenance,
    )
    (out / "report.md").write_text(report, encoding="utf-8")

    artifact_digests = B.digest_files(out, B.DERIVED_FILES)
    content = B.combined_digest(artifact_digests)
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "mode": mode,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "bundle": str(bundle_dir.resolve()),
        "code_digest": _code_digest(),
        "policy_id": policy.policy_id,
        "policy_status": policy.status,
        "policy_digest": policy.digest,
        "evidence_digest": frozen["evidence_digest"],
        "is_synthetic": frozen["is_synthetic"],
        "derived_artifacts": artifact_digests,
        "content_digest": content,
        "network_calls": 0,
        "llm_calls": 0,
        "not_in_scope": {"datasets": "GEO research is M1b", "editing": "not_assessed", "panel": "panel step not implemented"},
        "source_failures": source_failures,
        "run_status": "partial" if source_failures else "complete",
    }
    B.write_json(out / "manifest.json", manifest)
    return Derived(features=features, content_digest=content, artifact_digests=artifact_digests)


def run(input_path: Path, out: Path) -> Derived:
    freeze(input_path, out)
    return derive(out, out, mode="run")


def replay(bundle_dir: Path, out: Path) -> tuple[Derived, bool]:
    """Re-derive from frozen objects only and compare with the original run."""
    original = B.read_json(bundle_dir / "manifest.json")
    result = derive(bundle_dir, out, mode="replay")
    return result, result.content_digest == original["content_digest"]
