#!/usr/bin/env python3
"""Plot column-collapse measured results against cited reference bands.

The input profiles are emitted by `column_collapse/main.rs` at the end of each
run. This script recomputes the same normalized run-out and deposit-height
metrics from those profiles and writes the reviewer-facing CSV + PNG.
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parent
PLOTS = ROOT / "plots"
R0 = 0.025
H_TOE = 0.005
CASES = [
    ("a", "res_a", "accept", "5.0 mm"),
    ("b", "res_b", "accept", "3.3 mm"),
    ("c", "res_c", "accept", "2.5 mm"),
    ("negctl", "res_negctl", "reject", "5.0 mm"),
]
RUNOUT_BAND = (2.40, 3.60)
HEIGHT_BAND = (0.80, 1.70)
LSP_RUNOUT = 4.40


def metrics(profile_path: Path) -> tuple[float, float]:
    rows = []
    with profile_path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            rows.append((float(row["x"]), float(row["h"])))
    toe = max((abs(x) for x, h in rows if h >= H_TOE), default=0.0)
    h_max = max((h for _, h in rows), default=0.0)
    return (toe - R0) / R0, h_max / R0


def verdict(value: float, band: tuple[float, float]) -> bool:
    lo, hi = band
    return lo <= value <= hi


def main() -> None:
    PLOTS.mkdir(exist_ok=True)
    rows = []
    for case, result_dir, expectation, spacing in CASES:
        runout, height = metrics(ROOT / result_dir / "profile.csv")
        in_band = verdict(runout, RUNOUT_BAND) and verdict(height, HEIGHT_BAND)
        passed = in_band if expectation == "accept" else not in_band
        rows.append(
            {
                "case": case,
                "spacing": spacing,
                "expectation": expectation,
                "runout_n": f"{runout:.2f}",
                "height_n": f"{height:.2f}",
                "band_verdict": "accept" if in_band else "reject",
                "run_result": "PASS" if passed else "FAIL",
            }
        )

    csv_path = PLOTS / "column_collapse_results.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    labels = [row["case"] for row in rows]
    x = range(len(rows))
    colors = ["#2b8cbe" if row["expectation"] == "accept" else "#c75146" for row in rows]

    fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.8), constrained_layout=True)
    for ax, key, band, ylabel in (
        (axes[0], "runout_n", RUNOUT_BAND, "normalized run-out (Linf-L0)/L0"),
        (axes[1], "height_n", HEIGHT_BAND, "normalized deposit height Hinf/L0"),
    ):
        lo, hi = band
        ax.axhspan(lo, hi, color="#d7ead2", alpha=0.9, label="pass band")
        ax.axhline(lo, color="#367c39", linewidth=1.2)
        ax.axhline(hi, color="#367c39", linewidth=1.2)
        ax.scatter(list(x), [float(row[key]) for row in rows], s=78, c=colors, zorder=3)
        ax.set_xticks(list(x), labels)
        ax.set_ylabel(ylabel)
        ax.grid(axis="y", color="#cccccc", linewidth=0.6, alpha=0.7)
        ax.set_axisbelow(True)
        ax.margins(x=0.12)

    axes[0].axhline(
        LSP_RUNOUT,
        color="#5c5c5c",
        linewidth=1.0,
        linestyle="--",
        label="LSP continuum",
    )
    axes[0].annotate("LSP 2011 continuum 4.40", (2.05, LSP_RUNOUT + 0.06), fontsize=8)
    axes[0].set_ylim(-0.35, 4.75)
    axes[1].set_ylim(0.55, 2.08)

    for ax in axes:
        ax.legend(loc="upper left", frameon=False, fontsize=8)

    fig.suptitle("Column collapse: measured output vs Lube/Lajeunesse/LSP bands", fontsize=11)
    fig.savefig(PLOTS / "column_collapse_reference_bands.png", dpi=180)


if __name__ == "__main__":
    main()
