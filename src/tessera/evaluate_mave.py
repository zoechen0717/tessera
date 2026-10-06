"""Evaluate a tessera ranking against MAVE readouts (EVAL §8.1a; pre-registered per gene).

MAVE scores are fetched from MaveDB and used ONLY here, never as ranking inputs.
Orders are compared on the same allele set: tessera (tier rank), AlphaMissense alone,
domain membership alone (K), a ClinVar lookup baseline, and random (expected values).
Ties are handled by expectation (average rank / expected hits within a tie group).
"""

from __future__ import annotations

import csv
import io
import json
import random
import re
import urllib.parse
from dataclasses import dataclass
from pathlib import Path

from tessera.sources.http import Fetcher

MAVEDB = "https://api.mavedb.org/api/v1"
_SUB = re.compile(r"^p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$")


def mave_scores(f: Fetcher, urn: str) -> dict[str, float]:
    """{'p.Ala121Gly': score} for single amino-acid substitutions (mean if repeated)."""
    r = f.get("mavedb", f"scores/{urn}", "live", f"{MAVEDB}/score-sets/{urllib.parse.quote(urn)}/scores",
              media="text/plain", ttl_days=float("inf"))
    acc: dict[str, list[float]] = {}
    for row in csv.DictReader(io.StringIO(r.text())):
        hp, sc = (row.get("hgvs_pro") or "").strip(), row.get("score")
        if not _SUB.match(hp) or sc in (None, "", "NA"):
            continue
        try:
            acc.setdefault(hp, []).append(float(sc))
        except ValueError:
            continue
    return {k: sum(v) / len(v) for k, v in acc.items()}


def avg_ranks(keys: list[tuple]) -> list[float]:
    """Average rank (1 = best) for 'higher key is better', ties share their mean rank."""
    order = sorted(range(len(keys)), key=lambda i: keys[i], reverse=True)
    ranks = [0.0] * len(keys)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and keys[order[j + 1]] == keys[order[i]]:
            j += 1
        for t in range(i, j + 1):
            ranks[order[t]] = (i + j) / 2 + 1
        i = j + 1
    return ranks


def spearman(x: list[float], y: list[float]) -> float:
    rx, ry = avg_ranks([(-v,) for v in x]), avg_ranks([(-v,) for v in y])
    n = len(x)
    mx, my = sum(rx) / n, sum(ry) / n
    cov = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    vx = sum((a - mx) ** 2 for a in rx) ** 0.5
    vy = sum((b - my) ** 2 for b in ry) ** 0.5
    return cov / (vx * vy) if vx and vy else float("nan")


def auroc(keys: list[tuple], labels: list[int]) -> float:
    pos = [k for k, l in zip(keys, labels) if l == 1]
    neg = [k for k, l in zip(keys, labels) if l == 0]
    if not pos or not neg:
        return float("nan")
    s = 0.0
    for p in pos:
        for q in neg:
            s += 1.0 if p > q else 0.5 if p == q else 0.0
    return s / (len(pos) * len(neg))


def expected_hits_at_k(keys: list[tuple], labels: list[int], k: int) -> float:
    groups: dict[tuple, list[int]] = {}
    for key, l in zip(keys, labels):
        groups.setdefault(key, []).append(l)
    taken, hits = 0, 0.0
    for key in sorted(groups, reverse=True):
        g = groups[key]
        if taken + len(g) <= k:
            hits += sum(g)
            taken += len(g)
        else:
            hits += (k - taken) * sum(g) / len(g)
            break
    return hits


@dataclass
class Row:
    hgvs_p: str
    tier: int
    k_order: int
    am: float | None
    clinvar: str


_K = {"concordant": 3, "partial": 2, "not_assessable": 1, "discordant": 0}


def clinvar_key(c: str) -> int:
    c = c.lower()
    if "pathogenic" in c and "benign" not in c and "conflicting" not in c:
        return 2
    if "uncertain" in c or "conflicting" in c:
        return 1
    if "benign" in c:
        return 0
    return 1


def load_rows(run_dir: Path) -> list[Row]:
    ann = {a["variant_id"]: a for a in map(json.loads, (run_dir / "annotations.jsonl").read_text().splitlines())}
    rows = []
    for r in csv.DictReader(open(run_dir / "ranked_variants.csv")):
        if r["primary_channel"] != "discovery" or r["stratum"] != "missense" or not r["hgvs_p"]:
            continue
        a = ann.get(r["variant_id"], {})
        cv = (a.get("clinvar_classification") or {})
        rows.append(Row(r["hgvs_p"], int(r["tier_rank"]), _K.get(r["K"], 1),
                        float(r["P"]) if re.match(r"^-?[\d.]+(e-?\d+)?$", r["P"] or "") else None,
                        cv.get("value") or ""))
    return rows


def evaluate(run_dir: Path, urns: dict[str, dict], cache_dir: Path, out: Path, seed: int = 42,
             n_boot: int = 2000) -> dict:
    f = Fetcher(out / "snapshots", cache_dir=cache_dir)
    rows = load_rows(run_dir)
    orders = {
        "tessera": lambda r: (-r.tier,),
        "alphamissense_only": lambda r: (r.am if r.am is not None else -1.0,),
        "domain_only": lambda r: (r.k_order,),
        "clinvar_lookup": lambda r: (clinvar_key(r.clinvar),),
    }
    results = {"run": str(run_dir), "n_ranked_missense": len(rows), "score_sets": {}}
    rng = random.Random(seed)
    for urn, spec in urns.items():
        scores = mave_scores(f, urn)
        ev = [(r, scores[r.hgvs_p]) for r in rows if r.hgvs_p in scores]
        res = {"n_scored": len(ev), "n_mave_substitutions": len(scores), "role": spec["role"], "orders": {}}
        if len(ev) < 5:
            results["score_sets"][urn] = res
            continue
        labels = None
        if "abnormal_max" in spec:
            lab = [1 if s <= spec["abnormal_max"] else 0 if s > spec["normal_min"] else None for _, s in ev]
            keep = [i for i, l in enumerate(lab) if l is not None]
            labels = (keep, [lab[i] for i in keep])
            res["n_abnormal"] = sum(labels[1])
            res["n_normal"] = len(labels[1]) - res["n_abnormal"]
        for name, fn in orders.items():
            keys = [fn(r) for r, _ in ev]
            # higher key = nominated earlier; low MAVE score = loss of function → expect negative
            # correlation between key and score; report as positive "agreement".
            o = {"spearman_agreement": -spearman([k[0] for k in keys], [s for _, s in ev]),
                 "distinct_keys": len(set(keys))}
            if labels:
                kk = [keys[i] for i in labels[0]]
                ll = labels[1]
                o["auroc"] = auroc(kk, ll)
                for k in (10, 25, 50):
                    if k <= len(ll):
                        h = expected_hits_at_k(kk, ll, k)
                        o[f"precision@{k}"] = h / k
                        o[f"recall@{k}"] = h / max(1, sum(ll))
                boots = []
                for _ in range(n_boot):
                    idx = [rng.randrange(len(ll)) for _ in ll]
                    boots.append(auroc([kk[i] for i in idx], [ll[i] for i in idx]))
                boots = sorted(b for b in boots if b == b)
                o["auroc_ci95"] = [boots[int(0.025 * len(boots))], boots[int(0.975 * len(boots)) - 1]]
            res["orders"][name] = o
        if labels:
            res["random_expected"] = {"auroc": 0.5, "precision": res["n_abnormal"] / len(labels[1])}
        results["score_sets"][urn] = res
    f.write_index()
    out.mkdir(parents=True, exist_ok=True)
    (out / "results.json").write_text(json.dumps(results, indent=1))
    return results
