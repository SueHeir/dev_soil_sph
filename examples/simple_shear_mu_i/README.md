# `simple_shear_mu_i` — μ(I) flow-law recovery (return-map element test)

An **isolation** validation of the Dunatunga–Kamrin / μ(I) Drucker–Prager stress
update ([`sph_constitutive::update_stress`]). No SPH neighbors, no `App`, no
substrate — just the constitutive update driven by a homogeneous simple-shear
velocity gradient. It answers one question: *does the return map, in a
controlled homogeneous flow, reproduce the μ(I) inertial rheology it is built
to encode?* This is the standalone gate that must hold before the same update is
trusted inside the full SPH solver.

## Method

For a log-spaced sweep of inertial numbers `I ∈ [2×10⁻³, 0.8]`:

1. Fix the density (`ρ = 1600 kg/m³ > ρ_c`) → a well-defined confining pressure
   `p`. Choose the shear rate `γ̇ = I·√p/(d√ρ_s)` so the steady inertial number
   equals the target `I`.
2. Drive `update_stress` with `L_xy = γ̇` to steady state (fixed total shear
   strain; convergence checked — a 10× finer step leaves μ unchanged to 6 s.f.).
3. Read the steady stress ratio `μ = τ̄/p`.

Then fit the three-parameter Jop form
`μ(I) = μ_s + (μ_2−μ_s)·I/(I+I_0)` to the recovered `(I, μ)` points by
Levenberg–Marquardt (initial guess taken from the data, not the target) and
compare the fitted constants to the material's target values.

## Reference

- P. Jop, Y. Forterre & O. Pouliquen, *A constitutive law for dense granular
  flows*, **Nature 441**, 727–730 (2006). Glass-bead constants
  μ_s = 0.38 (≈ tan 20.9°), μ_2 = 0.64 (≈ tan 32.6°), I_0 = 0.28
  (`MaterialParams::glass_beads_v0`).
- GDR MiDi, *On dense granular flows*, **Eur. Phys. J. E 14**, 341–365 (2004) —
  the μ(I) inertial-number rheology.

The target constants are fixed by the citation; the test does not adjust them.

## Pass criteria

- Fitted `(μ_s, μ_2, I_0)` each within **2%** of the JFP-2006 constants
  (observed error ≈ 0.1%, a ≳20× margin).
- Every swept point within `5×10⁻³` in stress ratio of the exact μ(I); RMS
  residual of the fit `≤ 2×10⁻³`.
- **I-collapse:** the same `I` reached from a different `(γ̇, p)` (via a 4× higher
  pressure) yields the same μ within `5×10⁻³`. This one bound is sized to a named
  numerical artifact — the return map's weak dependence on the dimensionless
  elastic stiffness `G/p`, an O(few×10⁻³) effect — not to the answer.

## Run

```sh
cargo run --release --example simple_shear_mu_i
```

Exit 0 = PASS, nonzero = FAIL. Latest recovery: fitted
`μ_s = 0.38001, μ_2 = 0.63985, I_0 = 0.28042` (RMS 5.2×10⁻⁶).
