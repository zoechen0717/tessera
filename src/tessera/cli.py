"""tessera CLI: run | replay (M0, offline)."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from tessera.pipeline import replay, run


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="tessera")
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run", help="freeze inputs and rank (offline)")
    r.add_argument("--input", required=True, type=Path)
    r.add_argument("--output", required=True, type=Path)
    fe = sub.add_parser("fetch", help="M1a: fetch structured sources and assemble an input dir (network)")
    fe.add_argument("--request", required=True, type=Path)
    fe.add_argument("--output", required=True, type=Path)
    sc = sub.add_parser("screen", help="M1b: build the literature pool and triage it (network + LLM)")
    sc.add_argument("--request", required=True, type=Path)
    sc.add_argument("--output", required=True, type=Path)
    sc.add_argument("--llm", choices=["live", "replay", "off"], default="live",
                    help="off: rules + full-text scan only (recorded LLM outputs are still reused)")
    sc.add_argument("--refresh-pool", action="store_true",
                    help="re-run the literature searches; only papers not screened before cost LLM calls")
    dr = sub.add_parser("deepread", help="M1b layer 3: supplements scan + passage-based claim extraction")
    dr.add_argument("--request", required=True, type=Path)
    dr.add_argument("--literature", required=True, type=Path, help="directory written by `screen`")
    dr.add_argument("--llm", choices=["live", "replay", "off"], default="live")
    ev = sub.add_parser("evaluate-mave", help="score a ranked run against MAVE readouts (pre-registered spec)")
    ev.add_argument("--spec", required=True, type=Path)
    ev.add_argument("--output", required=True, type=Path)
    p = sub.add_parser("replay", help="re-derive ranks from a frozen bundle; no network, no model")
    p.add_argument("--bundle", required=True, type=Path)
    p.add_argument("--output", required=True, type=Path)
    args = ap.parse_args(argv)

    if args.cmd == "fetch":
        from tessera.m1a import fetch

        inp = fetch(args.request, args.output)
        print(f"input assembled: {inp}")
        return 0
    if args.cmd == "screen":
        from tessera.screening import screen

        print(f"documents table: {screen(args.request, args.output, llm=args.llm, refresh_pool=args.refresh_pool)}")
        return 0
    if args.cmd == "deepread":
        from tessera.deep_read import deep_read

        deep_read(args.request, args.literature, llm=args.llm)
        return 0
    if args.cmd == "evaluate-mave":
        import json as _json
        import yaml as _yaml

        from tessera.evaluate_mave import evaluate

        spec = _yaml.safe_load(args.spec.read_text())
        root = Path(__file__).resolve().parents[2]
        res = evaluate((args.spec.parent / spec["run"]).resolve(), spec["score_sets"], root / "cache" / "http", args.output)
        print(_json.dumps(res, indent=1))
        return 0
    if args.cmd == "run":
        res = run(args.input, args.output)
        print(f"run complete: {len(res.features)} alleles; content digest {res.content_digest[:16]}")
        print(f"report: {args.output / 'report.md'}")
        return 0
    res, same = replay(args.bundle, args.output)
    print(f"replay content digest {res.content_digest[:16]} — {'IDENTICAL' if same else 'DIFFERS'}")
    return 0 if same else 1


if __name__ == "__main__":
    sys.exit(main())
