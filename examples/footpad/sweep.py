#!/usr/bin/env python3
"""Footpad bearing/sinkage validation against an independent Bekker/DEM oracle.

The Rust example writes force vs sinkage. This driver fits only the seated
loading branch, subtracting the force and sinkage at the seating depth, then
checks the excess pressure follows the Bekker/Wong pressure-sinkage form
validated in DIRT's independent DEM plate-sinkage benchmark.
"""

import csv
import math
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HERE = ROOT / "examples" / "footpad"
DATA = HERE / "data"
PLOTS = HERE / "plots"
POSITIVE_CFG = HERE / "config.toml"
NEGATIVE_CFG = HERE / "config_negctl.toml"
POSITIVE_CSV = DATA / "positive" / "sinkage.csv"
NEGATIVE_CSV = DATA / "zero_g" / "sinkage.csv"
REFERENCE_CSV = DATA / "dirt_bekker_reference.csv"

# Geometry in config.toml: x footprint [-0.02, 0.02], y period [0, 0.015].
PLATE_AREA = 0.04 * 0.015

# Ignore first 2 mm of plate travel: this SPH plate is built from three frozen
# layers, so early samples include contact seating and force-offset pickup rather
# than a Bekker loading branch.
SEATING_DEPTH = 0.002
MAX_FIT_DEPTH = 0.0085
R2_MIN = 0.90
RISE_MIN_PA = 500.0

# DIRT's DEM plate-sinkage benchmark validates the same Bekker form with
# granular-soil exponent band 0.4 <= n <= 1.6. For this SPH gate we also require
# the fitted exponent to land within a modest absolute margin around the
# DIRT-DEM cases in data/dirt_bekker_reference.csv.
BEKKER_N_MIN = 0.4
BEKKER_N_MAX = 1.6
DIRT_N_MARGIN = 0.25


def run_example(config: Path, log: Path) -> int:
    cmd = ["cargo", "run", "--release", "--example", "footpad", "--", str(config)]
    with log.open("w") as f:
        proc = subprocess.run(cmd, cwd=ROOT, stdout=f, stderr=subprocess.STDOUT)
    return proc.returncode


def load_curve(path: Path):
    rows = []
    with path.open() as f:
        for r in csv.DictReader(f):
            z = float(r["sinkage"])
            force = float(r["force_z"])
            if math.isfinite(z) and math.isfinite(force):
                rows.append((z, force / PLATE_AREA))
    return rows


def reference_band():
    exponents = []
    with REFERENCE_CSV.open() as f:
        for r in csv.DictReader(f):
            if r["monotone"] == "True" and r["p_increases"] == "True":
                exponents.append(float(r["n"]))
    lo = max(BEKKER_N_MIN, min(exponents) - DIRT_N_MARGIN)
    hi = min(BEKKER_N_MAX, max(exponents) + DIRT_N_MARGIN)
    return lo, hi, min(exponents), max(exponents)


def fit_seated_bekker(rows):
    window = [(z, p) for z, p in rows if SEATING_DEPTH <= z <= MAX_FIT_DEPTH]
    if len(window) < 8:
        return {"ok": False, "reason": "insufficient seated samples"}

    z0, p0 = window[0]
    excess = [(z - z0, p - p0) for z, p in window if z > z0 and p > p0]
    if len(excess) < 8:
        return {"ok": False, "reason": "insufficient rising excess-pressure samples"}

    xs = [math.log(z) for z, _ in excess]
    ys = [math.log(p) for _, p in excess]
    mx = sum(xs) / len(xs)
    my = sum(ys) / len(ys)
    denom = sum((x - mx) ** 2 for x in xs)
    n = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / denom
    log_a = my - n * mx
    ss = sum((y - my) ** 2 for y in ys)
    res = sum((y - (log_a + n * x)) ** 2 for x, y in zip(xs, ys))
    r2 = 1.0 - res / ss if ss > 0.0 else 0.0
    rise = excess[-1][1]
    binned = []
    nbin = 10
    for i in range(nbin):
        chunk = excess[i * len(excess) // nbin : (i + 1) * len(excess) // nbin]
        if chunk:
            binned.append(sum(p for _, p in chunk) / len(chunk))
    monotone = all(binned[i + 1] >= binned[i] * 0.98 for i in range(len(binned) - 1))
    lo, hi, dirt_lo, dirt_hi = reference_band()
    checks = {
        "n_in_reference_band": lo <= n <= hi,
        "r2": r2 >= R2_MIN,
        "pressure_rises": rise >= RISE_MIN_PA,
        "monotone": monotone,
    }
    return {
        "ok": all(checks.values()),
        "reason": "",
        "n": n,
        "A": math.exp(log_a),
        "r2": r2,
        "rise": rise,
        "z0": z0,
        "p0": p0,
        "count": len(excess),
        "checks": checks,
        "band_lo": lo,
        "band_hi": hi,
        "dirt_lo": dirt_lo,
        "dirt_hi": dirt_hi,
        "excess": excess,
        "binned": binned,
    }


def write_summary(pos, neg, pos_rc, neg_rc):
    DATA.mkdir(parents=True, exist_ok=True)
    with (DATA / "validation_summary.csv").open("w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["case", "expect", "exit_code", "n", "r2", "rise_pa", "accepted"])
        for name, expect, rc, result in [
            ("positive", "accept", pos_rc, pos),
            ("zero_g", "reject", neg_rc, neg),
        ]:
            w.writerow([
                name,
                expect,
                rc,
                f"{result.get('n', float('nan')):.8g}" if "n" in result else "",
                f"{result.get('r2', float('nan')):.8g}" if "r2" in result else "",
                f"{result.get('rise', float('nan')):.8g}" if "rise" in result else "",
                result.get("ok", False),
            ])


def cleanup_transient_outputs():
    for path in [DATA / "positive.log", DATA / "zero_g.log"]:
        path.unlink(missing_ok=True)
    for path in [DATA / "positive" / "dump", DATA / "zero_g" / "dump"]:
        if path.exists():
            shutil.rmtree(path)


def plot(pos, neg):
    PLOTS.mkdir(parents=True, exist_ok=True)
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7.0, 5.0))
    for result, label, color, marker in [
        (pos, "SPH footpad seated branch", "#1b7f5a", "o"),
        (neg, "zero-g control (must reject)", "#a33c2f", "x"),
    ]:
        if "excess" not in result:
            continue
        z = [dz * 1000.0 for dz, _ in result["excess"]]
        p = [p for _, p in result["excess"]]
        ax.loglog(z, p, marker, ms=4, linestyle="none", color=color, label=label)
        if result.get("ok"):
            zz = [min(z) * (max(z) / min(z)) ** (i / 80) for i in range(81)]
            pp = [result["A"] * (mm / 1000.0) ** result["n"] for mm in zz]
            ax.loglog(zz, pp, "-", color=color, lw=1.4, label=f"fit n={result['n']:.2f}")

    lo, hi, _, _ = reference_band()
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    ax.text(
        0.02,
        0.98,
        f"Bekker/DIRT DEM exponent gate: {lo:.2f} <= n <= {hi:.2f}\n"
        f"positive: n={pos.get('n', float('nan')):.2f}, R2={pos.get('r2', 0):.3f} -> "
        f"{'PASS' if pos.get('ok') else 'FAIL'}\n"
        f"zero-g control -> {'REJECTED' if not neg.get('ok') else 'ACCEPTED'}",
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=9,
        bbox={"facecolor": "white", "edgecolor": "0.7", "alpha": 0.9},
    )
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)
    ax.set_xlabel("seated sinkage increment dz (mm)")
    ax.set_ylabel("excess plate pressure dp (Pa)")
    ax.set_title("Footpad pressure-sinkage vs Bekker/DIRT DEM reference gate")
    ax.grid(True, which="both", ls=":", alpha=0.4)
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(PLOTS / "footpad_bekker_validation.png")
    plt.close(fig)


def main():
    DATA.mkdir(parents=True, exist_ok=True)
    pos_rc = run_example(POSITIVE_CFG, DATA / "positive.log")
    neg_rc = run_example(NEGATIVE_CFG, DATA / "zero_g.log")

    pos = fit_seated_bekker(load_curve(POSITIVE_CSV))
    neg = fit_seated_bekker(load_curve(NEGATIVE_CSV)) if NEGATIVE_CSV.exists() else {"ok": False}
    write_summary(pos, neg, pos_rc, neg_rc)
    plot(pos, neg)
    cleanup_transient_outputs()

    print("=== footpad bearing/sinkage validation ===")
    print(
        f"reference: Bekker/Wong p=(kc/b+kphi) z^n; DIRT DEM n range "
        f"{pos.get('dirt_lo', float('nan')):.3f}-{pos.get('dirt_hi', float('nan')):.3f}, "
        f"gate {pos.get('band_lo', float('nan')):.3f}-{pos.get('band_hi', float('nan')):.3f}"
    )
    if "n" in pos:
        print(
            f"positive: n={pos['n']:.3f}, R2={pos['r2']:.3f}, "
            f"rise={pos['rise']:.1f} Pa, samples={pos['count']}"
        )
    else:
        print(f"positive: {pos.get('reason', 'no fit')}")
    print(f"positive mechanics exit code: {pos_rc}")
    print(f"zero-g control mechanics exit code: {neg_rc} (non-zero is expected)")
    if "n" in neg:
        print(f"zero-g control oracle fit: n={neg['n']:.3f}, R2={neg['r2']:.3f}, accepted={neg['ok']}")
    else:
        print(f"zero-g control oracle fit: {neg.get('reason', 'rejected before fit')}")

    ok = (pos_rc == 0) and pos.get("ok", False) and not neg.get("ok", False)
    if ok:
        print("PASS: footpad force-sinkage curve matches the independent Bekker/DIRT DEM gate; zero-g control rejected")
        return 0
    print("FAIL: footpad bearing/sinkage validation failed")
    return 1


if __name__ == "__main__":
    sys.exit(main())
