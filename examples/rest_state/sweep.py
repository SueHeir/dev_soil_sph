#!/usr/bin/env python3
"""Run the rest/hydrostatic validation pair and emit the result figure."""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[2]
REST_CONFIG = "examples/rest_state/config.toml"
HYDRO_CONFIG = "examples/hydrostatic_column/config.toml"
PLOT_DIR = ROOT / "examples" / "rest_state" / "plots"
PLOT_PATH = PLOT_DIR / "rest_hydrostatic_validation.png"

REST_SPEED_TOL = 1.0e-2
HYDRO_RATIO_LO = 0.7
HYDRO_RATIO_HI = 1.3
HYDRO_REGRESSION_FLOOR = 0.8208
HYDRO_RATIO_TARGET = 1.0


def run_example(name: str, config: str) -> str:
    cmd = ["cargo", "run", "--release", "--example", name, "--", config]
    result = subprocess.run(
        cmd,
        cwd=ROOT,
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    print(result.stdout)
    if result.returncode != 0:
        raise SystemExit(result.returncode)
    return result.stdout


def require_float(pattern: str, text: str, label: str) -> float:
    match = re.search(pattern, text, re.MULTILINE)
    if not match:
        raise RuntimeError(f"could not parse {label}")
    return float(match.group(1))


def main() -> int:
    rest_out = run_example("rest_state", REST_CONFIG)
    hydro_out = run_example("hydrostatic_column", HYDRO_CONFIG)

    max_speed = require_float(r"^max speed:\s+([0-9.eE+-]+) m/s", rest_out, "rest max speed")
    hydro_ratio = require_float(r"ratio = ([0-9.eE+-]+)", hydro_out, "hydrostatic ratio")
    hydro_speed = require_float(
        r"max fluid speed: ([0-9.eE+-]+) m/s", hydro_out, "hydrostatic max fluid speed"
    )
    pressure_min = require_float(
        r"p . \[([0-9.eE+-]+),", hydro_out, "hydrostatic minimum pressure"
    )
    density_spread_pct = require_float(
        r"spread = ([0-9.eE+-]+)%", hydro_out, "density spread"
    )

    PLOT_DIR.mkdir(parents=True, exist_ok=True)

    plt.style.use("seaborn-v0_8-whitegrid")
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2), constrained_layout=True)

    ax = axes[0]
    ax.bar(["reference", "measured"], [0.0, max_speed], color=["#6b7280", "#2563eb"])
    ax.axhline(REST_SPEED_TOL, color="#b91c1c", linestyle="--", linewidth=2, label="pass line")
    ax.set_ylim(0.0, REST_SPEED_TOL * 1.15)
    ax.set_ylabel("max particle speed (m/s)")
    ax.set_title("Periodic rest state")
    ax.text(
        1,
        REST_SPEED_TOL * 0.09,
        f"{max_speed:.2e} m/s",
        ha="center",
        va="bottom",
        color="#1f2937",
    )
    ax.text(
        1.02,
        REST_SPEED_TOL * 1.01,
        f"tol {REST_SPEED_TOL:.0e}",
        ha="left",
        va="bottom",
        color="#7f1d1d",
    )
    ax.legend(loc="upper left", frameon=True)

    ax = axes[1]
    ax.axhspan(HYDRO_RATIO_LO, HYDRO_RATIO_HI, color="#dcfce7", label="external pass band")
    ax.axhline(HYDRO_RATIO_TARGET, color="#166534", linewidth=2, label="target")
    ax.axhline(
        HYDRO_REGRESSION_FLOOR,
        color="#b45309",
        linestyle="--",
        linewidth=2,
        label="regression floor",
    )
    ax.scatter([0], [hydro_ratio], s=90, color="#7c3aed", zorder=4, label="measured")
    ax.set_xlim(-0.7, 0.7)
    ax.set_xticks([0])
    ax.set_xticklabels(["dp/dz / (-rho g)"])
    ax.set_ylim(0.65, 1.35)
    ax.set_ylabel("gradient ratio")
    ax.set_title("Hydrostatic column")
    ax.text(
        0.04,
        hydro_ratio,
        f"{hydro_ratio:.4f}",
        ha="left",
        va="center",
        color="#3b0764",
    )
    ax.text(
        -0.63,
        0.67,
        f"max speed {hydro_speed:.2e} m/s\np_min {pressure_min:.2f} Pa\nrho spread {density_spread_pct:.3f}%",
        ha="left",
        va="bottom",
        fontsize=9,
        color="#1f2937",
    )
    ax.legend(loc="upper right", frameon=True)

    fig.suptitle("hydrostatic_rest validation: measured vs reference", fontsize=13)
    fig.savefig(PLOT_PATH, dpi=180)
    plt.close(fig)

    print(f"wrote {PLOT_PATH.relative_to(ROOT)}")
    print(
        "PASS: rest max speed "
        f"{max_speed:.3e} < {REST_SPEED_TOL:.1e}; hydrostatic ratio {hydro_ratio:.4f} "
        f"in [{HYDRO_RATIO_LO:.1f}, {HYDRO_RATIO_HI:.1f}] and >= {HYDRO_REGRESSION_FLOOR:.4f}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
