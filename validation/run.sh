#!/usr/bin/env bash
# dev_soil_sph validation set — one-command, skeptic-facing runner.
#
# Builds and runs ONLY the validation examples (rest_state, hydrostatic_column,
# column_collapse aspect sweep + a/b/c resolution check, and the μ(I) return-map recovery), each of which asserts
# a NUMERIC pass/fail against an external reference (Bui 2008 tensile stability;
# Lube 2005 / Lagrée-Staron-Popinet 2011 run-out & deposit scalings; Jop-Forterre-
# Pouliquen 2006 / GDR MiDi 2004 μ(I) constants). Any FAIL exits non-zero.
#
# The demoted demos (haff_cooling, shear_heating, conduction_test, footpad,
# defluidization) are intentionally NOT run here — see validation/manifest.toml.
#
# Usage:  source ~/projects/.build-env && validation/run.sh
set -uo pipefail
cd "$(dirname "$0")/.."   # repo root

pass=0; fail=0
run() { # <name> <config> <extra-run-args...>
  local name="$1" cfg="$2"; shift 2
  echo "── $name  ($cfg) ───────────────────────────────────────────"
  local log; log="$(mktemp)"
  # Run ONCE; the example itself exits non-zero on FAIL (its own numeric checks).
  if cargo run --release --example "$name" -- "$cfg" "$@" >"$log" 2>/dev/null; then
    grep -E "^(PASS|FAIL)" "$log" || echo "  (ran, no PASS line)"
    pass=$((pass+1))
  else
    grep -E "^(PASS|FAIL)" "$log" || echo "  (no result)"
    fail=$((fail+1)); echo "  -> FAILED: $name ($cfg)"
  fi
  rm -f "$log"
}

run_noarg() { # <name>   self-contained example (target constants baked in, no config file)
  local name="$1"
  echo "── $name  (self-contained) ─────────────────────────────────"
  local log; log="$(mktemp)"
  # Exit code is the verdict (0 = PASS); the example prints its own checks.
  if cargo run --release --example "$name" >"$log" 2>/dev/null; then
    grep -E "^(PASS|FAIL|ALL CHECKS|checks:)" "$log" || echo "  (ran, exit 0 = PASS)"
    pass=$((pass+1))
  else
    grep -E "^(PASS|FAIL|CHECKS FAILED)" "$log" || echo "  (no result)"
    fail=$((fail+1)); echo "  -> FAILED: $name"
  fi
  rm -f "$log"
}

echo "=== dev_sph validation set ==="
run rest_state         examples/rest_state/config.toml
run hydrostatic_column examples/hydrostatic_column/config.toml
run column_collapse    examples/column_collapse/config_sweep_a0p5.toml
run column_collapse    examples/column_collapse/config_sweep_a1.toml
run column_collapse    examples/column_collapse/config_a.toml
run column_collapse    examples/column_collapse/config_b.toml
run column_collapse    examples/column_collapse/config_c.toml
run column_collapse    examples/column_collapse/config_sweep_a3.toml
run column_collapse    examples/column_collapse/config_sweep_a6.toml
# Falsifiability control: an over-frictional material MUST be rejected by the same
# reference band that a/b/c sit inside. The config declares [validation] expect =
# "reject"; the example inverts its verdict, so a green here proves the band can
# fail. (A green would be impossible if the gate were vacuous — see README.)
run column_collapse    examples/column_collapse/config_negctl.toml
# μ(I) return-map recovery: shear a pressurized sample across a range of inertial
# numbers I and fit (μ_s, μ_2, I_0); the example asserts the fit reproduces the
# Jop-Forterre-Pouliquen 2006 glass-bead constants within tolerance (exit != 0 on
# fail). Self-contained, so no config file.
run_noarg simple_shear_mu_i

echo "=========================================="
echo "validation set: $pass passed, $fail failed"
[ "$fail" -eq 0 ] || exit 1
echo "ALL VALIDATIONS PASSED"
