"""Build the synthetic M0 fixture, run, and replay. Offline; no credentials.

    ./.venv/bin/python examples/synthetic_demo.py [output_dir]
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from tessera.pipeline import replay, run  # noqa: E402
from tessera.synthetic import build  # noqa: E402

out = Path(sys.argv[1] if len(sys.argv) > 1 else "runs/synthetic_demo").resolve()
inp = build(out / "input")
res = run(inp, out / "run")
_, same = replay(out / "run", out / "replay")
print(f"{len(res.features)} alleles ranked; replay {'identical' if same else 'DIFFERS'}")
print(f"report: {out / 'run' / 'report.md'}")
