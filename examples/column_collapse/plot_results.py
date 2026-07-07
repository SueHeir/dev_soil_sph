#!/usr/bin/env python3
"""Plot column-collapse aspect-sweep results against cited reference bands.

Profiles are emitted by `column_collapse/main.rs`. This script recomputes the
same normalized run-out and deposit-height metrics, writes the reviewer-facing
CSV, and draws the pass/fail bands from Lube/Lajeunesse/LSP Eqs. 3.1-3.2.
"""

from __future__ import annotations

import csv
import tomllib
from pathlib import Path

import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parent
PLOTS = ROOT / "plots"
H_TOE = 0.005
CASES = [
    ("a=0.5", "config_sweep_a0p5.toml", "res_sweep_a0p5", "outside_reference"),
    ("a=1", "config_sweep_a1.toml", "res_sweep_a1", "outside_reference"),
    ("a=2", "config_a.toml", "res_a", "accept"),
    ("a=3", "config_sweep_a3.toml", "res_sweep_a3", "accept"),
    ("a=6", "config_sweep_a6.toml", "res_sweep_a6", "accept"),
    ("a=2 neg", "config_negctl.toml", "res_negctl", "reject"),
]


def column_geometry(config_path: Path) -> tuple[float, float]:
    with config_path.open("rb") as handle:
        data = tomllib.load(handle)
    for insert in data["sph"]["insert"]:
        if insert.get("frozen", False):
            continue
        xmin, _, zmin = insert["region_min"]
        xmax, _, zmax = insert["region_max"]
        return max(abs(xmin), abs(xmax)), zmax - zmin
    raise RuntimeError(f"no fluid insert in {config_path}")


def runout_band(a: float) -> tuple[float, float]:
    if abs(a - 2.0) < 1.0e-9:
        return 2.40, 3.60
    if a < 2.0:
        return 1.2 * a, 2.2 * a
    return 1.9 * a ** (2.0 / 3.0), max(2.3 * a ** (2.0 / 3.0), 2.2 * a)


def height_band(a: float) -> tuple[float, float]:
    if abs(a - 2.0) < 1.0e-9:
        return 0.80, 1.70
    if a < 2.0:
        return 0.75 * a, 1.25 * a
    return 0.65 * a**0.35, 1.30 * a**0.40


def lsp_continuum_runout(a: float) -> float:
    if a < 7.0:
        return 2.2 * a
    return 3.9 * a**0.7


def metrics(profile_path: Path, r0: float) -> tuple[float, float]:
    rows = []
    with profile_path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            rows.append((float(row["x"]), float(row["h"])))
    toe = max((abs(x) for x, h in rows if h >= H_TOE), default=0.0)
    h_max = max((h for _, h in rows), default=0.0)
    return (toe - r0) / r0, h_max / r0


def in_band(value: float, band: tuple[float, float]) -> bool:
    lo, hi = band
    return lo <= value <= hi


def main() -> None:
    PLOTS.mkdir(exist_ok=True)
    rows = []
    for label, config_name, result_dir, expectation in CASES:
        r0, h0 = column_geometry(ROOT / config_name)
        aspect = h0 / r0
        runout, height = metrics(ROOT / result_dir / "profile.csv", r0)
        run_band = runout_band(aspect)
        h_band = height_band(aspect)
        band_accept = in_band(runout, run_band) and in_band(height, h_band)
        passed = band_accept if expectation == "accept" else not band_accept
        rows.append(
            {
                "case": label,
                "config": config_name,
                "aspect": f"{aspect:.2f}",
                "expectation": expectation,
                "runout_n": f"{runout:.2f}",
                "runout_lo": f"{run_band[0]:.2f}",
                "runout_hi": f"{run_band[1]:.2f}",
                "height_n": f"{height:.2f}",
                "height_lo": f"{h_band[0]:.2f}",
                "height_hi": f"{h_band[1]:.2f}",
                "band_verdict": "accept" if band_accept else "reject",
                "run_result": "PASS" if passed else "FAIL",
            }
        )

    csv_path = PLOTS / "column_collapse_results.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    aspects = [float(row["aspect"]) for row in rows]
    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.0), constrained_layout=True)
    for ax, key, lo_key, hi_key, ylabel in (
        (axes[0], "runout_n", "runout_lo", "runout_hi", "normalized run-out (Linf-L0)/L0"),
        (axes[1], "height_n", "height_lo", "height_hi", "normalized height Hinf/L0"),
    ):
        band_lo = [float(row[lo_key]) for row in rows]
        band_hi = [float(row[hi_key]) for row in rows]
        band_a = [float(row["aspect"]) for row in rows]
        ax.fill_between(band_a, band_lo, band_hi, color="#d7ead2", alpha=0.9, label="pass band")
        ax.plot(band_a, band_lo, color="#367c39", linewidth=1.2)
        ax.plot(band_a, band_hi, color="#367c39", linewidth=1.2)
        colors = ["#2b8cbe" if row["expectation"] == "accept" else "#c75146" for row in rows]
        ax.scatter(aspects, [float(row[key]) for row in rows], s=78, c=colors, zorder=3)
        for row in rows:
            ax.annotate(row["case"], (float(row["aspect"]), float(row[key])), xytext=(4, 4),
                        textcoords="offset points", fontsize=8)
        ax.set_xlabel("initial aspect ratio a = H0/L0")
        ax.set_ylabel(ylabel)
        ax.grid(color="#cccccc", linewidth=0.6, alpha=0.7)
        ax.set_axisbelow(True)

    a_curve = [0.5, 1.0, 2.0, 3.0, 6.0]
    axes[0].plot(
        a_curve,
        [lsp_continuum_runout(a) for a in a_curve],
        color="#5c5c5c",
        linewidth=1.0,
        linestyle="--",
        label="LSP continuum",
    )
    axes[0].set_ylim(0.0, max(14.0, axes[0].get_ylim()[1]))
    axes[1].set_ylim(0.0, max(2.2, axes[1].get_ylim()[1]))
    for ax in axes:
        ax.legend(loc="upper left", frameon=False, fontsize=8)

    fig.suptitle("Column-collapse aspect sweep vs Lube/Lajeunesse/LSP Eqs. 3.1-3.2", fontsize=11)
    fig.savefig(PLOTS / "column_collapse_reference_bands.png", dpi=180)


if __name__ == "__main__":
    main()
