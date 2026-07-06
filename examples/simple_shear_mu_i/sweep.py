#!/usr/bin/env python3
"""Run the simple-shear mu(I) gate and plot the measured-vs-reference result."""

from __future__ import annotations

import math
import re
import subprocess
import sys
from pathlib import Path

import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[2]
EXAMPLE_DIR = ROOT / "examples" / "simple_shear_mu_i"
PLOT_DIR = EXAMPLE_DIR / "plots"
PLOT_PATH = PLOT_DIR / "mu_i_recovery.png"

POINT_TOL = 5.0e-3
PARAM_TOL_FRAC = 0.02
RMS_TOL = 2.0e-3


def jop(i_value: float, mu_s: float, mu_2: float, i0: float) -> float:
    return mu_s + (mu_2 - mu_s) * i_value / (i_value + i0)


def parse_output(output: str) -> dict[str, object]:
    point_re = re.compile(
        r"^\s*([0-9.]+e[+-][0-9]+)\s+"
        r"([0-9.]+)\s+"
        r"([0-9.]+)\s+"
        r"([0-9.]+)\s+"
        r"([0-9.]+e[+-][0-9]+)\s*$",
        re.IGNORECASE,
    )
    fit_re = re.compile(
        r"fitted:\s+.*?=([0-9.]+)\s+.*?=([0-9.]+)\s+.*?=([0-9.]+)"
        r"\s+\(RMS residual ([0-9.]+e[+-][0-9]+)\)",
        re.IGNORECASE,
    )
    target_re = re.compile(
        r"target:\s+.*?=([0-9.]+)\s+.*?=([0-9.]+)\s+.*?=([0-9.]+)"
    )
    collapse_re = re.compile(r"I-collapse check.*\|Delta\|=([0-9.]+e[+-][0-9]+)")

    points: list[tuple[float, float, float, float]] = []
    fitted = None
    target = None
    rms = None
    collapse = None
    passed = "ALL CHECKS PASSED" in output

    for line in output.splitlines():
        if match := point_re.match(line):
            i_value = float(match.group(1))
            mu_measured = float(match.group(3))
            mu_target = float(match.group(4))
            residual = mu_measured - mu_target
            points.append((i_value, mu_measured, mu_target, residual))
        elif match := fit_re.search(line):
            fitted = tuple(float(match.group(k)) for k in range(1, 4))
            rms = float(match.group(4))
        elif match := target_re.search(line):
            target = tuple(float(match.group(k)) for k in range(1, 4))
        elif match := collapse_re.search(line.replace("Δ", "Delta")):
            collapse = float(match.group(1))

    if not points or fitted is None or target is None or rms is None or collapse is None:
        raise ValueError("could not parse simple_shear_mu_i output")

    return {
        "points": points,
        "fitted": fitted,
        "target": target,
        "rms": rms,
        "collapse": collapse,
        "passed": passed,
    }


def plot(parsed: dict[str, object]) -> None:
    points = parsed["points"]
    fitted = parsed["fitted"]
    target = parsed["target"]
    rms = parsed["rms"]
    collapse = parsed["collapse"]
    passed = parsed["passed"]
    assert isinstance(points, list)
    assert isinstance(fitted, tuple)
    assert isinstance(target, tuple)
    assert isinstance(rms, float)
    assert isinstance(collapse, float)
    assert isinstance(passed, bool)

    i_values = [p[0] for p in points]
    mu_measured = [p[1] for p in points]
    residuals = [p[3] for p in points]
    xs = [
        math.exp(math.log(min(i_values)) + k * (math.log(max(i_values)) - math.log(min(i_values))) / 300)
        for k in range(301)
    ]
    target_curve = [jop(x, *target) for x in xs]
    fitted_curve = [jop(x, *fitted) for x in xs]

    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "figure.dpi": 130,
        }
    )
    fig = plt.figure(figsize=(9.5, 7.0), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=[1.45, 1.0])
    ax_curve = fig.add_subplot(gs[0, :])
    ax_resid = fig.add_subplot(gs[1, 0])
    ax_params = fig.add_subplot(gs[1, 1])

    ax_curve.fill_between(
        xs,
        [y - POINT_TOL for y in target_curve],
        [y + POINT_TOL for y in target_curve],
        color="#d9ead3",
        alpha=0.9,
        label="point pass band: target +/- 5e-3",
    )
    ax_curve.plot(xs, target_curve, color="#2f5d50", linewidth=2.2, label="JFP/GDR MiDi target")
    ax_curve.plot(xs, fitted_curve, color="#a6422b", linewidth=2.0, linestyle="--", label="fit to measured")
    ax_curve.scatter(i_values, mu_measured, color="#1d4e89", edgecolor="white", linewidth=0.7, zorder=4, label="measured")
    ax_curve.set_xscale("log")
    ax_curve.set_xlabel("inertial number I")
    ax_curve.set_ylabel("stress ratio mu")
    ax_curve.set_title("simple_shear_mu_i: return-map recovery of the Jop mu(I) law")
    status = "PASS" if passed else "FAIL"
    ax_curve.text(
        0.02,
        0.95,
        (
            f"{status}: fitted mu_s={fitted[0]:.5f}, mu_2={fitted[1]:.5f}, I0={fitted[2]:.5f}\n"
            f"target=({target[0]:.3f}, {target[1]:.3f}, {target[2]:.3f}); "
            f"I-collapse {collapse:.2e} <= 5e-3"
        ),
        transform=ax_curve.transAxes,
        ha="left",
        va="top",
        fontsize=9,
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.82, "pad": 3},
    )
    ax_curve.legend(loc="lower right", frameon=False)
    ax_curve.grid(True, which="both", color="#e5e5e5", linewidth=0.7)

    ax_resid.axhspan(-POINT_TOL, POINT_TOL, color="#d9ead3", alpha=0.9, label="+/- 5e-3 point gate")
    ax_resid.axhline(0.0, color="#666666", linewidth=1.0)
    ax_resid.plot(i_values, residuals, color="#1d4e89", marker="o", linewidth=1.4)
    ax_resid.set_xscale("log")
    ax_resid.set_xlabel("inertial number I")
    ax_resid.set_ylabel("measured - target")
    ax_resid.set_title(f"residuals: RMS {rms:.2e} <= {RMS_TOL:.1e}")
    ax_resid.legend(loc="lower right", frameon=False)
    ax_resid.grid(True, which="both", color="#e5e5e5", linewidth=0.7)

    labels = ["mu_s", "mu_2", "I0"]
    pct_errors = [(fitted[k] - target[k]) / target[k] * 100.0 for k in range(3)]
    colors = ["#1d4e89" if abs(e) <= PARAM_TOL_FRAC * 100.0 else "#a6422b" for e in pct_errors]
    ax_params.axhspan(-PARAM_TOL_FRAC * 100.0, PARAM_TOL_FRAC * 100.0, color="#d9ead3", alpha=0.9)
    ax_params.axhline(0.0, color="#666666", linewidth=1.0)
    ax_params.bar(labels, pct_errors, color=colors, width=0.58)
    ax_params.set_ylabel("fit error vs target [%]")
    ax_params.set_title("fitted constants within +/- 2% gate")
    ax_params.grid(True, axis="y", color="#e5e5e5", linewidth=0.7)
    for idx, error in enumerate(pct_errors):
        vertical = "bottom" if error >= 0 else "top"
        offset = 0.05 if error >= 0 else -0.05
        ax_params.text(idx, error + offset, f"{error:+.3f}%", ha="center", va=vertical, fontsize=9)

    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(PLOT_PATH, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    result = subprocess.run(
        ["cargo", "run", "--release", "--example", "simple_shear_mu_i"],
        cwd=ROOT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    sys.stdout.write(result.stdout)
    if result.returncode != 0:
        return result.returncode
    parsed = parse_output(result.stdout)
    plot(parsed)
    print(f"\nwrote {PLOT_PATH.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
