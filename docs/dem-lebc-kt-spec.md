# LEBC Rheometer + Kinetic-Theory Validation — implementation brief

**For the `dirt`/`soil`-side work** (a separate effort). This is the Tier-1 rheometer of `docs/dem-campaign.md`: homogeneous simple shear of glass-bead spheres under Lees–Edwards boundary conditions (LEBC), producing the μ(I) and Φ(I) closure the dev_soil_sph SPH solver consumes — and validated, in the regime where it is exact, against granular kinetic theory.

**Division of labor.** This brief defines the DEM/substrate side (build LEBC + measure stress/temperature + KT cross-check + fit μ(I),Φ(I)). The dev_soil_sph/SPH side consumes the resulting calibration file (`docs/dem-campaign.md` §6) and is unaffected by how it's produced. Keep all changes in `soil`/`dirt`; dev_soil_sph does not depend on them at build time.

---

## Part A — the LEBC simple-shear primitive (the missing piece)

`soil_deform` today only stretches the diagonal box lengths (engineering strain rate on x/y/z); it has **no xy-shear / Lees–Edwards** primitive. That is the gap. Two paths:

- **Quick path (recommended for first data):** run Tier-1 in **LAMMPS** via Liz's existing `../lammps_shear_cell` (`fix deform xy erate … remap v` already does true LEBC, stress already wired to `stress_tensors.dat`). Two edits: turn on interparticle friction (`xmu` → μ_p = 0.5; currently 0) and set restitution e ≈ 0.7 with a softened E; and either add σ_yy pressure control or keep the fixed-volume φ-sweep as the fallback.
- **Native path (the build):** add a Lees–Edwards homogeneous-shear capability to `soil`/`dirt` so the whole calibration→SPH pipeline lives in one stack.

### A1. Lees–Edwards homogeneous shear (native build)
Standard LEBC for shear rate γ̇ with flow = x, gradient = y, vorticity = z:
- Box periodic in all three directions; the y-periodic images are offset in x by `Δx(t) = γ̇ · L_y · t (mod L_x)`.
- A particle crossing the **upper** y-boundary has its position wrapped by `−Δx` in x and its velocity by `−γ̇ L_y` in x (and vice-versa at the lower boundary).
- **Comm/ghost implication (the real work):** ghost atoms built across the y-boundary must carry the same `±Δx` position shift and `±γ̇ L_y` velocity offset. This touches `soil_core`'s `borders`/`forward_comm` (the periodic-offset machinery already exists for ordinary PBC — `CommTopology.periodic_swap` — so this is an extension of that offset to a time-dependent, shear-coupled x-offset on the y-swaps, not a new subsystem). Honor the no-physics-in-soil rule: LEBC is a boundary/domain concern, so it belongs in `soil` (e.g. `soil_deform` or a new `soil` module), not in a physics tier.
- Initialize with a linear velocity profile `v_x = γ̇ y` so the system starts near steady.

### A2. Stress control (constant pressure) vs fixed volume
The campaign wants **pressure-controlled** shear (fixed σ_yy = P, Φ an output). Options:
- **Fixed-volume φ-sweep (do this first — simplest):** hold the box fixed, prescribe Φ, shear, measure the resulting P and τ. Sweep Φ to trace μ(I), Φ(I). This is what `lammps_shear_cell` already does and is fully sufficient to fit the closure.
- **Pressure control (add later):** a barostat on `L_y` that rescales the gradient direction to drive measured σ_yy → P (Berendsen-style relaxation toward the target). Gives Φ(P) directly and matches the spec, but is not required for a valid μ(I)/Φ(I) fit.

---

## Part B — measurement (record more than τ and P)

Per run, after reaching steady state, time-average:

1. **Full stress tensor σ_ij** (not just τ=σ_xy and P). DIRT already accumulates the complete Love–Weber tensor — `VirialStress` in `soil_core/src/virial.rs`, summed (normal + tangential) by the Hertz–Mindlin contact loop. It has **no output path**; add a ~20-line recorder (read `Res<VirialStress>`, divide by box volume, write CSV — same shape as the plate-sinkage recorder). Report `p = ⅓ tr σ`, `τ = √(½ σ′:σ′)`, and the **normal-stress differences** `N₁ = σ_xx − σ_yy`, `N₂ = σ_yy − σ_zz`.
2. **Granular temperature** `T = ⅓⟨|v_i − v̄(y_i)|²⟩` — **subtract the mean shear profile** `v̄(y) = γ̇ y` before computing fluctuations, or T is swamped by the mean flow. This is the key new measurement that unlocks the KT comparison.
3. **Solid fraction Φ** (and ρ = Φ ρ_s).
4. The control inputs γ̇ and (P or Φ).

Form the dimensionless groups: `I = γ̇ d / √(P/ρ_s)`, `μ_eff = τ/P`, and `Φ(I)`.

---

## Part C — kinetic-theory validation (frictionless sub-sweep)

KT is exact for **smooth (frictionless) inelastic spheres** in the dilute-to-moderate-Φ (collisional) regime, and breaks near jamming. Use it to **validate the stress-measurement pipeline** before trusting it in the dense regime where there is no analytic check.

**Protocol:** run a dedicated **frictionless sub-sweep** (μ_p = 0, fixed e, e.g. 0.7 and 0.9) across several Φ and γ̇, and compare the measured stress tensor to the KT prediction below. Expect quantitative agreement at moderate Φ (say Φ ≲ 0.55); growing deviation as Φ → Φ_c ≈ 0.64 (enduring contacts) marks the dense regime where DEM is irreplaceable.

**KT closure (3D smooth inelastic spheres).** With `g₀(Φ) = (2 − Φ)/(2(1 − Φ)³)` (Carnahan–Starling):

- **Pressure:** `p = ρ_s Φ T [ 1 + 2(1 + e) Φ g₀ ]`  (kinetic `ρΦT` + collisional term).
- **Shear stress:** `σ_xy = η γ̇`, with `η = η*(Φ, e) · ρ_s d √T` (the dimensionless `η*` from your chosen closure — Lun et al. 1984 or Garzó–Dufty 1999; both standard).
- **Dissipation rate** (collisional cooling): `Γ = γ*(Φ, e) · ρ_s Φ² g₀ (1 − e²) T^{3/2} / d`.
- **Steady-shear energy balance** (production = dissipation): `σ_xy γ̇ = Γ` ⇒ closes `T`, giving `T ∝ (γ̇ d)² · F(Φ, e)` (**Bagnold scaling**, σ ∝ ρ_s d² γ̇²). Substituting back makes `μ = σ_xy/p` and `Φ` **functions of I alone** — i.e. KT *derives* a μ(I)/Φ(I) law in the collisional regime.

Use the exact `η*`, `γ*` from one canonical reference and cite it (Lun, Savage, Jeffrey & Chepurniy 1984; Garzó & Dufty 1999; Brilliantov & Pöschel 2004). The validation passes if measured `p`, `σ_xy`, `T` match KT within a few % in the moderate-Φ collisional band.

**Why this matters downstream:** (i) it validates the virial-stress recorder you just built; (ii) KT is the natural closure for the **collisional branch (I ≳ 0.3)** that μ(I) cannot represent — the fast-impact/ejecta regime flagged in `docs/literature-review.md` — and a future dev_soil_sph constitutive may add it; (iii) KT predicts the normal-stress differences `N₁, N₂` that the dev_soil_sph μ(I) Drucker–Prager omits, telling us quantitatively whether that omission matters at our operating point.

---

## Part D — production sweep → the calibration the SPH solver consumes

With the pipeline validated:
- Run the **frictional** production sweep (μ_p = 0.5, e ≈ 0.7) across `I ∈ [1e-4, 0.5]` (≥5 points/decade) per `docs/dem-campaign.md` §3.2.
- Fit `μ(I) = μ_s + (μ_2 − μ_s)/(I_0/I + 1)` and `Φ(I)`; extract `ρ_c = Φ_max ρ_s` and an effective `K` (separate isotropic-compression run).
- Emit the calibration file (`docs/dem-campaign.md` §6 YAML). That file is the entire interface to dev_soil_sph — `MaterialParams` in `mud_constitutive` is populated from it (replacing the `glass_beads_v0()` literature anchors).

---

## Acceptance criteria
1. LEBC produces a homogeneous, steady linear shear profile (no shear banding, plateau in T and stress).
2. Frictionless sub-sweep matches KT (`p`, `σ_xy`, `T`) within a few % at moderate Φ.
3. Frictional sweep yields a clean μ(I), Φ(I) collapse; μ_s, μ_2, I_0 in plausible glass-bead ranges (≈0.38 / 0.64 / 0.28 as anchors, possibly lower for spheres).
4. Calibration YAML written and round-trips into `mud_constitutive::MaterialParams`.

## References
- GDR MiDi 2004; da Cruz et al. 2005 — the I-scaling from DEM simple shear (`docs/literature-review.md`).
- Lun, Savage, Jeffrey & Chepurniy 1984; Garzó & Dufty 1999; Brilliantov & Pöschel 2004 — KT transport coefficients.
- Jenkins 2007 / Berzi–Jenkins — extended (dense) kinetic theory near jamming.
- Dunatunga & Kamrin 2015 (our SPH core) — uses Jenkins–Berzi KT to justify the stress-free separation.

*Drafted 2026-06-16. Companions: `docs/dem-campaign.md` (spec), `docs/dem-campaign-dirt.md` (DIRT/LAMMPS mapping), `docs/physics-design.md` (what dev_soil_sph does with the fit).*
