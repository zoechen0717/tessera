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
    p = sub.add_parser("replay", help="re-derive ranks from a frozen bundle; no network, no model")
    p.add_argument("--bundle", required=True, type=Path)
    p.add_argument("--output", required=True, type=Path)
    args = ap.parse_args(argv)

    if args.cmd == "fetch":
        from tessera.m1a import fetch

        inp = fetch(args.request, args.output)
        print(f"input assembled: {inp}")
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
