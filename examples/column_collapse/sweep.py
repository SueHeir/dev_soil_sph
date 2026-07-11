#!/usr/bin/env python3
"""Run the declarative column-collapse aspect sweep and verify its result plot.

Each configuration owns its geometry, material, integration window, and expected
reference verdict.  The Rust example enforces the Lube/Lajeunesse/LSP bands;
this driver only orchestrates that checked-in case matrix and regenerates the
reviewer-facing measured-versus-reference figure.
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


def run_case(config: str) -> None:
    print(f"=== {config} ===", flush=True)
    command = ["cargo", "run", "--release", "--example", "column_collapse", "--", str(HERE / config)]
    completed = subprocess.run(command, cwd=ROOT)
    if completed.returncode:
        raise RuntimeError(f"{config} exited {completed.returncode}")


def verify_plot_summary() -> None:
    summary = HERE / "plots" / "column_collapse_results.csv"
    with summary.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != 6:
        raise RuntimeError(f"expected six plotted aspect/control results, found {len(rows)}")
    failed = [row["case"] for row in rows if row["run_result"] != "PASS"]
    if failed:
        raise RuntimeError(f"plot summary reports failed cases: {', '.join(failed)}")
    if not (HERE / "plots" / "column_collapse_reference_bands.png").is_file():
        raise RuntimeError("aspect-sweep reference-band figure was not generated")


def main() -> int:
    try:
        for config in CASES:
            run_case(config)
        subprocess.run([sys.executable, str(HERE / "plot_results.py")], cwd=ROOT, check=True)
        verify_plot_summary()
    except (OSError, RuntimeError, subprocess.CalledProcessError) as error:
        print(f"CHECKS FAILED: {error}")
        return 1
    print("ALL CHECKS PASSED: column-collapse aspect sweep and negative control")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
