#!/usr/bin/env python3
"""Run the declarative column-collapse aspect sweep and verify its result plot.

The physical cases are all positive experimental checks.  A non-zero exit from
one is preserved as a non-zero sweep result; the driver still runs every case
and regenerates the figure so a failure is quantitative rather than hidden.
Only config_negctl is deliberately inverted: it must be rejected by the same
external band.
"""

from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
CASES = [
    "config_sweep_a0p5.toml",
    "config_sweep_a1.toml",
    "config_a.toml",
    "config_b.toml",
    "config_c.toml",
    "config_sweep_a3.toml",
    "config_sweep_a6.toml",
    "config_negctl.toml",
]


def run_case(binary: Path, config: str) -> bool:
    print(f"=== {config} ===", flush=True)
    completed = subprocess.run([str(binary), str(HERE / config)], cwd=ROOT)
    # The executable owns the control inversion: the wrong-physics config exits
    # zero only when the external band rejects it.  The driver must therefore
    # require zero from every listed case, rather than inverting it a second time.
    passed = completed.returncode == 0
    if not passed:
        print(f"FAIL: {config} exit={completed.returncode}; expected 0")
    return passed


def verify_plot_summary() -> bool:
    summary = HERE / "plots" / "column_collapse_results.csv"
    with summary.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != 6:
        raise RuntimeError(f"expected six plotted aspect/control results, found {len(rows)}")
    if not (HERE / "plots" / "column_collapse_reference_bands.png").is_file():
        raise RuntimeError("aspect-sweep reference-band figure was not generated")
    failed = [row["case"] for row in rows if row["run_result"] != "PASS"]
    if failed:
        print("FAIL: external-reference misses in " + ", ".join(failed))
        return False
    return True


def main() -> int:
    try:
        subprocess.run(
            ["cargo", "build", "--release", "--example", "column_collapse"],
            cwd=ROOT,
            check=True,
        )
        binary = ROOT / "target" / "release" / "examples" / "column_collapse"
        case_ok = [run_case(binary, config) for config in CASES]
        subprocess.run([sys.executable, str(HERE / "plot_results.py")], cwd=ROOT, check=True)
        plot_ok = verify_plot_summary()
    except (OSError, RuntimeError, subprocess.CalledProcessError) as error:
        print(f"CHECKS FAILED: {error}")
        return 1
    if all(case_ok) and plot_ok:
        print("ALL CHECKS PASSED: column-collapse aspect sweep and negative control")
        return 0
    print("CHECKS FAILED: positive cases outside their external-reference bands")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
