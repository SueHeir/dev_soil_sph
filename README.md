# dev_soil_sph

> **Development-stage research tier.** dev_soil_sph is a proof-of-concept
> granular-SPH physics tier built on the `grass` scheduler and `soil` particle
> substrate. Treat it as research code: read the docs, run the validation gates,
> and do not assume production readiness just because an example runs.

dev_soil_sph implements a meshfree granular-continuum model for dry glass-bead-like
materials: SPH density/gradient operators, a density-based tension-free
separation criterion, and a Dunatunga-Kamrin-style elasto-viscoplastic stress
update with a Jop/GDR MiDi μ(I) Drucker-Prager flow law. It is a sibling physics
tier to DIRT and dev_soil_peri, not a DEM solver.

## From an SPH App to a coupled participant

SOIL defines and stores the particle data model; SPH plugins in this repo own
the continuum operators and constitutive state. GRASS holds those resources in
the App and schedules the systems, but does not define particles or SPH fields.
A cross-substrate coupling therefore lives outside both solvers. The runnable
consumer is
[`dev_couple_sph_cfd`](https://github.com/SueHeir/dev_couple_sph_cfd): its
[`seam.rs`](https://github.com/SueHeir/dev_couple_sph_cfd/blob/main/crates/sph_cfd/src/seam.rs)
reads and writes solver resources from the parent between child ticks. It does
not place cross-App access inside the SPH child scheduler.

Start here:

| Document | Purpose |
|---|---|
| [`DISCLAIMER.md`](DISCLAIMER.md) | Provenance and dev-tier caution. |
| [`docs/sph-primer.md`](docs/sph-primer.md) | SPH method primer for this problem. |
| [`docs/literature-review.md`](docs/literature-review.md) | Source papers and model-family choices. |
| [`docs/physics-design.md`](docs/physics-design.md) | Governing equations, constitutive law, and parameters. |
| [`docs/architecture.md`](docs/architecture.md) | How dev_soil_sph rides the GRASS -> SOIL stack. |
| [`docs/dem-campaign.md`](docs/dem-campaign.md) | DEM calibration/validation campaign spec. |
| [`docs/dem-campaign-dirt.md`](docs/dem-campaign-dirt.md) | Mapping the campaign onto DIRT/LAMMPS capabilities. |
| [`docs/dem-lebc-kt-spec.md`](docs/dem-lebc-kt-spec.md) | LEBC rheometer and kinetic-theory validation brief. |
| [`validation/README.md`](validation/README.md) | Skeptic-facing validation status and caveats. |

## Validation Status

The current validated set is deliberately small and numeric. Each gate exits
non-zero when it leaves the documented pass band.

| Gate | Example(s) | Independent reference | Current status |
|---|---|---|---|
| μ(I) return-map recovery | [`examples/simple_shear_mu_i`](examples/simple_shear_mu_i/README.md) | Jop, Forterre & Pouliquen 2006; GDR MiDi 2004 | PASS: fitted `mu_s = 0.38001`, `mu_2 = 0.63985`, `I_0 = 0.28042` against glass-bead constants. |
| Hydrostatic/rest stability | [`examples/rest_state`](examples/rest_state), [`examples/hydrostatic_column`](examples/hydrostatic_column) | Bui, Fukagawa, Sako & Ohno 2008 tensile-instability criterion | PASS: rest state preserved; hydrostatic gradient ratio `0.821`; compressive pressure and no clumping. |
| Column collapse aspect sweep | [`examples/column_collapse`](examples/column_collapse) | Lube/Lajeunesse 2005 via Lagree, Staron & Popinet 2011 | **FAIL:** `a = 0.5, 1, 6` miss the cited experimental run-out bands; `a = 2, 3` pass and the wrong-physics negative control is rejected. The generated graph keeps all misses visible. |
| Footpad bearing/sinkage | [`examples/footpad`](examples/footpad/README.md) | Bekker/Wong pressure-sinkage; independent DIRT DEM plate-sinkage benchmark | PASS: seated branch `n = 0.957`, `R^2 = 0.960`; zero-gravity control is rejected. |

Run the full checked-in validation harness:

```bash
source ~/projects/.build-env
validation/run.sh
```

The μ(I) isolation gate can also be run directly:

```bash
cargo run --release --example simple_shear_mu_i
```

## Demoted Demos

These examples remain useful smoke tests or coupling seeds, but they are **not**
part of the validation claim because they do not recover an independent
experimental or cross-code reference:

| Demo | Role |
|---|---|
| [`examples/haff_cooling`](examples/haff_cooling) | Kinetic-theory granular-temperature decay smoke check. |
| [`examples/shear_heating`](examples/shear_heating) | Shear-production / Bagnold steady-state operator showcase. |
| [`examples/conduction_test`](examples/conduction_test) | SPH Laplacian conduction operator check. |
| [`examples/defluidization`](examples/defluidization) | KT-to-contact stress hand-off transient; SPH-CFD seed. |

The distinction matters: demos can prove that machinery executes, while
validations prove that dev_soil_sph reproduces a published physical target within a
declared numerical band.
