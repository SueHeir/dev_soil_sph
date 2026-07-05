#!/usr/bin/env python3
"""Run hydrostatic_column and plot the measured validation gates."""

from __future__ import annotations

import argparse
import re
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[2]
EXAMPLE = ROOT / "examples" / "hydrostatic_column"
CONFIG = EXAMPLE / "config.toml"
PLOTS = EXAMPLE / "plots"
FIGURE = PLOTS / "hydrostatic_column_validation.svg"

RHO_REF = 1500.0
G = 9.81
GRADIENT_LO = 0.7
GRADIENT_HI = 1.3


def run_example() -> str:
    cmd = [
        "cargo",
        "run",
        "--release",
        "--example",
        "hydrostatic_column",
        "--",
        str(CONFIG.relative_to(ROOT)),
    ]
    proc = subprocess.run(
        cmd,
        cwd=ROOT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    if proc.returncode != 0:
        print(proc.stdout)
        raise SystemExit(proc.returncode)
    return proc.stdout


def parse_output(output: str) -> dict[str, object]:
    bins: list[tuple[float, int, float, float, float]] = []
    for line in output.splitlines():
        match = re.match(
            r"^\s*([0-9.]+)\s+(\d+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s*$",
            line,
        )
        if match:
            bins.append(
                (
                    float(match.group(1)),
                    int(match.group(2)),
                    float(match.group(3)),
                    float(match.group(4)),
                    float(match.group(5)),
                )
            )

    slope_match = re.search(
        r"dp/dz =\s*([-+0-9.]+) Pa/m\s+expected .*? =\s*([-+0-9.]+) Pa/m\s+ratio =\s*([-+0-9.]+)",
        output,
    )
    floor_match = re.search(r"ratio >=\s*([-+0-9.]+)", output)
    if floor_match is None:
        floor_match = re.search(r"ratio .?=\s*([-+0-9.]+)", output)
    tension_match = re.search(
        r"p\s+[^[]+\[([-+0-9.]+),\s*([-+0-9.]+)\] Pa\s+\(want p_min [^0-9-]*([-+0-9.]+)\)",
        output,
    )
    rho_match = re.search(r"spread =\s*([-+0-9.]+)%\s+\(want <\s*([-+0-9.]+)%\)", output)
    speed_match = re.search(r"settled: true \(max speed\s+([-+0-9.eE]+)\)", output)
    pass_match = "PASS: hydrostatic gradient" in output

    if not (bins and slope_match and floor_match and tension_match and rho_match and speed_match and pass_match):
        raise SystemExit("could not parse hydrostatic_column PASS output")

    return {
        "bins": bins,
        "slope": float(slope_match.group(1)),
        "expected": float(slope_match.group(2)),
        "ratio": float(slope_match.group(3)),
        "floor": float(floor_match.group(1)),
        "p_min": float(tension_match.group(1)),
        "p_max": float(tension_match.group(2)),
        "p_limit": float(tension_match.group(3)),
        "rho_spread_pct": float(rho_match.group(1)),
        "rho_limit_pct": float(rho_match.group(2)),
        "max_speed": float(speed_match.group(1)),
    }


def save_plot(data: dict[str, object], path: Path) -> None:
    bins = sorted(data["bins"])  # type: ignore[arg-type]
    z = [row[0] for row in bins]
    p = [row[2] for row in bins]
    p_hydro = [row[3] for row in bins]
    rho = [row[4] for row in bins]

    z_surface = sum(zi + ph / (RHO_REF * G) for zi, ph in zip(z, p_hydro)) / len(z)
    p_lo = [GRADIENT_LO * RHO_REF * G * (z_surface - zi) for zi in z]
    p_hi = [GRADIENT_HI * RHO_REF * G * (z_surface - zi) for zi in z]
    p_floor = [float(data["floor"]) * RHO_REF * G * (z_surface - zi) for zi in z]

    z_bar = sum(z) / len(z)
    p_bar = sum(p) / len(p)
    p_fit = [p_bar + float(data["slope"]) * (zi - z_bar) for zi in z]

    fig, (ax_p, ax_gate) = plt.subplots(
        1,
        2,
        figsize=(10.0, 4.8),
        gridspec_kw={"width_ratios": [1.55, 1.0]},
        constrained_layout=True,
    )

    ax_p.fill_betweenx(z, p_lo, p_hi, color="#d7e8f7", label="0.7-1.3 x -rho g")
    ax_p.plot(p_hydro, z, color="#1b6f5f", linewidth=2.0, label="-rho g reference")
    ax_p.plot(p_floor, z, color="#bb6b00", linestyle="--", linewidth=1.8, label="regression floor")
    ax_p.plot(p_fit, z, color="#4b2e83", linewidth=2.0, label="fitted dp/dz")
    ax_p.scatter(p, z, color="#202124", s=32, zorder=4, label="measured slab mean")
    ax_p.set_xlabel("Pressure (Pa)")
    ax_p.set_ylabel("Height z (m)")
    ax_p.set_title(
        f"dp/dz ratio {float(data['ratio']):.4f}; PASS >= {float(data['floor']):.4f}",
        fontsize=10,
    )
    ax_p.grid(True, alpha=0.25)
    ax_p.legend(loc="lower right", fontsize=8)

    labels = ["p_min (Pa)", "rho spread (%)"]
    values = [float(data["p_min"]), float(data["rho_spread_pct"])]
    limits = [float(data["p_limit"]), float(data["rho_limit_pct"])]
    colors = ["#2f7d32", "#2f7d32"]
    y = range(len(labels))
    ax_gate.barh(y, values, color=colors, height=0.45, label="measured")
    ax_gate.axvline(limits[0], color="#9f2f27", linestyle="--", linewidth=1.8)
    ax_gate.plot([limits[1]], [1], marker="|", color="#9f2f27", markersize=22, markeredgewidth=2)
    ax_gate.text(limits[0], -0.32, f"p_min limit {limits[0]:.2f}", color="#9f2f27", fontsize=8)
    ax_gate.text(limits[1], 1.22, f"density limit {limits[1]:.1f}%", color="#9f2f27", fontsize=8, ha="center")
    ax_gate.set_yticks(list(y), labels)
    ax_gate.set_xlabel("Gate value")
    ax_gate.set_title(
        f"tensile/density gates PASS; max speed {float(data['max_speed']):.2e} m/s",
        fontsize=10,
    )
    ax_gate.grid(True, axis="x", alpha=0.25)

    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, metadata={"Date": None})
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=FIGURE)
    args = parser.parse_args()

    output = run_example()
    data = parse_output(output)
    save_plot(data, args.output)
    print(
        "PASS hydrostatic_column: "
        f"ratio={float(data['ratio']):.4f}, "
        f"p_min={float(data['p_min']):.2f} Pa, "
        f"rho_spread={float(data['rho_spread_pct']):.3f}%, "
        f"figure={args.output}"
    )


if __name__ == "__main__":
    main()
