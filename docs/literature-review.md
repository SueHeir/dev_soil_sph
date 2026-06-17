# Literature Review — DEM-Informed SPH for Lunar Regolith (Lander-Leg Touchdown)

**Goal of the codebase:** a continuum SPH model of dry lunar regolith for lander-leg touchdown, in which **DEM micro-scale simulations are the ground truth that calibrates (or derives) the continuum constitutive law** the SPH solver integrates. The hard physics is the **quasi-static-solid ↔ dense-fluid transition**, and the target application — a footpad penetrating, decelerating, and arresting in a granular bed — *is* a granular-intrusion problem.

This review covers the 10 papers in `./papers`. Each entry is condensed from a full read; equation numbers refer to the original papers. The cross-cutting synthesis (§III) is the part that drives architecture decisions; the per-paper entries (§II) are reference.

> Companion: `docs/sph-primer.md` (the SPH method itself). This file is about *what constitutive physics to put inside it and how to calibrate it from DEM*.

---

## I. The big picture in one diagram

```
            MICRO (grain scale)                 MACRO (continuum)              SOLVER
        ┌────────────────────────┐       ┌───────────────────────────┐   ┌──────────────┐
  DEM   │ plane-shear / triaxial │  -->  │  μ(I), Φ(I), μ_s, μ_2,     │-->│  SPH momentum│
 (ours) │ controlled-P or -V     │  fit  │  I_0, ρ_c, K, cohesion c  │   │  + return map│
        └────────────────────────┘       └───────────────────────────┘   └──────────────┘
              da Cruz '05                  Jop '06 / GDR MiDi '04            this codebase
              (origin of μ(I))             Dunatunga-Kamrin '15
                                                                          ▲ validate against ▲
                                           Lagrée-Staron-Popinet '11  ──  DEM column collapse
                                           Szewc '16 (SPH vs DEM)
                                                                          ▲ cross-check forces ▲
                                           Agarwal '21/'23 DRFT / 3D-RFT  — footpad oracle
```

**Three constitutive families appear in these papers** — this is the central design fork:

| Family | Representative papers | Core idea | Solid↔fluid transition handled by |
|---|---|---|---|
| **A. μ(I) effective-viscosity (fluid-first)** | Minatti–Paris '15, Szewc '16, Salehizadeh–Shafiei '19, Lagrée–Staron–Popinet '11 | Treat dense granular as a non-Newtonian fluid, viscosity η = μ(I)·p/\|γ̇\| | A **viscosity cap** (Papanastasiou / η_max) — quasi-static = very-high-viscosity creep |
| **B. Elasto-viscoplastic + density separation (Kamrin)** | Dunatunga–Kamrin '15, Agarwal '21, Agarwal '23 | Hypoelastic solid below yield; μ(I) Drucker–Prager flow above; **stress-free below ρ_c** | A genuine **elastic solid branch** + density-triggered tension-free separation |
| **C. Elastic-plastic Drucker–Prager (soil-mechanics)** | Bui '08 | Classic elasto-perfectly-plastic soil; pressure from the constitutive model, not an EOS | Yield-surface return mapping; flow = plastic flow on the DP cone |

Families **B** and **C** share a return-mapping skeleton and a real solid branch; **A** is the simplest to implement but has no true solid (only capped creep). The foundations (GDR MiDi '04, Jop '06) supply the μ(I) law that **A** and **B** both consume.

---

## II. Per-paper entries

### Foundations — where μ(I) comes from (and its DEM provenance)

#### 1. GDR MiDi (2004), *Eur. Phys. J. E* 14:341 — "On dense granular flows" · DOI 10.1140/epje/i2003-10153-0
The meta-analysis that established the **inertial number** `I = γ̇ d / √(P/ρ)` as the single governing parameter of dense flow, by collating six geometries (plane shear, annular, chute, incline, heap, drum). Result: `μ_eff = τ/P = μ(I)` and volume fraction `Φ(I)` are functions of I alone; `Φ_max ≈ 0.85 → ~0.80` by I≈0.1.
- **DEM role — central.** Plane shear (the cleanest "rheometer") is available *only* in DEM, not experiment. Reported DEM (MD and contact-dynamics) setups use fixed-pressure or fixed-volume control, periodic BCs, rigid-grain limit, restitution `e` and contact friction `μ_p` swept. **Key finding: the dense rheology is independent of `e` and `μ_p`** (for μ_p ≳ 0.1) and of material — only the dense→collisional transition I depends on `e`. This is what makes DEM calibration robust: a few contact parameters can be set loosely.
- **For us:** the I-number + dual scaling μ(I), Φ(I) is the baseline; the plane-shear DEM protocol (fixed-P, periodic, rigid grains) is our calibration rig. Caveats: cohesionless, local rheology only (fails for surface flows / near thresholds → nonlocal needed), mostly 2D.

#### 2. Jop, Forterre & Pouliquen (2006), *Nature* 441:727 — "A constitutive law for dense granular flows" · DOI 10.1038/nature04801
The 3D tensorial μ(I) closure: `σ = −p I + η γ̇`, `η = μ(I) p / |γ̇|`, with
`μ(I) = μ_s + (μ_2 − μ_s)/(I_0/I + 1)`. Built-in yield: η → ∞ as γ̇ → 0 gives a Drucker–Prager threshold `|τ| > μ_s p`. **Parameters (glass beads): μ_s ≈ 0.38, μ_2 ≈ 0.64, I_0 = 0.279.**
- **DEM role:** none run here (calibrated on inclined-plane *experiments*), but the τ∝P, μ(I) form descends from da Cruz et al. 2005 DEM.
- **For us:** this *is* the target continuum law for families A/B. Regolith enters only through (μ_s, μ_2, I_0, d, ρ_s). **Limitations to fix:** incompressible (no dilatancy Φ(I)), cohesionless, steady-flow only, viscosity singularity needs regularization.

### The solid↔fluid transition models (the core fork)

#### 3. Dunatunga & Kamrin (2015), *J. Fluid Mech.* 779:483 — "Continuum modelling ... through their many phases" · DOI 10.1017/jfm.2015.383
**The canonical solid↔fluid↔separated model (Family B).** One framework spanning elastic solid, μ(I) viscoplastic flow, and a **disconnected, tension-free gas state**. State variable = density.
- **Separation / EOS (Eq. 2.7):** `p = 0 if ρ<ρ_c`, else `p = (K_c/ρ)(ρ−ρ_c)`.
- **Three-phase master relation (Eqs. 2.14–2.15):** `σ = 0 if ρ<ρ_c`; else hypoelastic `σ̌ = ℂ:(D − D̂ᵖ)` with codirectional isochoric plastic flow and μ(I) Drucker–Prager yield.
- **Implicit return map (Eqs. 3.9–3.16):** density via exponential volume map; elastic trial stress; if `p_tr ≤ 0` → σ=0; else solve a **quadratic for the equivalent shear stress** with a cancellation-safe root `τ̄ = 2H/(B + √(B²−4H))`, then rescale the deviator. **Fully point-local — no mesh/neighbor dependence.**
- **Parameters:** E=1 GPa, ν=0.3, μ_s=0.3819, μ_2=0.6435, ρ_s=2450, ρ_c=1500 (Φ_c≈60%).
- **DEM role:** μ(I) and the codirectionality assumption are DEM-justified (Silbert, da Cruz, Koval); column-collapse runout matches Staron–Hinch 2005 **2D DEM**.
- **For us (high relevance):** the granular-slug-hits-wall case (loose body densifies, transmits stress, arrests into a load-bearing heap) is kinematically a touchdown. **The constitutive engine ports MPM→SPH essentially unchanged** — it consumes (σⁿ, D, W, ρ, Δt) and returns σⁿ⁺¹; SPH supplies ρ and ∇v natively. The only SPH-specific work is the disconnection/free-surface numerics (tensile instability) that MPM sidesteps via its grid. **Limitations:** no cohesion, local (no length scale → no nonlocal creep), isochoric/critical-state (no dilatancy transients), 2D demos.

#### 4. Bui, Fukagawa, Sako & Ohno (2008), *IJNAMG* 32:1537 — SPH elastic-plastic geomaterial · DOI 10.1002/nag.688
**The canonical SPH-soil backbone (Family C).** First Drucker–Prager elastic-plastic soil in SPH.
- **Two architectural choices we want:** (i) **pressure from the constitutive model**, `p = −(σˣˣ+σʸʸ+σᶻᶻ)/3`, *not* an EOS — the natural slot for a DEM-calibrated law; (ii) **artificial-stress** fix for SPH tensile instability.
- **Yield (Eq. 16):** `√J₂ + α_φ I₁ − k_c = 0` with DP↔Mohr–Coulomb match `α_φ = tanφ/√(9+12tan²φ)`, `k_c = 3c/√(9+12tan²φ)`. Non-associated flow with ψ=0 (plastically incompressible) for the granular runs. Jaumann objective rate; return mapping by tension-cracking + deviator stress-scaling.
- **Artificial stress (Eqs. 63–69):** repulsive factor `f_ij = (W_ij/W(Δd,h))ⁿ`, **n≈2.55, ε≈0.5**; activates only when particles clump in tension. Cubic spline, h=1.2Δd, leapfrog, C_cour=0.2, sound speed c≈600 m/s.
- **DEM role:** none (calibrated on an aluminum-bar shear-box experiment). **This is exactly the DEM→continuum gap our project fills** — we replace the shear box with virtual DEM shear/triaxial tests.
- **For us:** borrow `p = −I₁/3`, the DP return-map skeleton, and the artificial-stress fix (or tension-cracking alone for cohesionless regolith). **Limitations:** elastic-perfectly-plastic only (no hardening/critical-state/rate dependence), small-deformation theory, 2D, quasi-static loading (no impact dynamics).

### μ(I)-in-SPH implementations & DEM validation (Family A)

#### 5. Minatti & Paris (2015), *Appl. Math. Modelling* 39:363 — SPH dense free-surface granular · DOI 10.1016/j.apm.2014.05.034
The primary **μ(I)-as-SPH-effective-viscosity** template. WCSPH; `η = μ(i) p / γ̇` (Eq. 11) with yield stress `τ_y = μ_s p`.
- **The viscosity cap we need (Eqs. 46–55):** a yield-stress-dependent **Papanastasiou** regularization, `η̂ = η_0(τ_y,γ̇) + τ_y(1−e^{−m_p γ̇})/γ̇`, with `m_p` adapted to the local yield stress and a cap `η̂_max ≥ 10× residual viscosity`. This is what keeps the quasi-static limit from collapsing the timestep.
- **Numerics:** Wendland C4 kernel, **conservative (uncorrected) gradients** (deliberately no MLS, to preserve momentum), Cleary harmonic-mean viscous operator, Morris ghost-particle walls, symplectic Verlet, dual sound/viscous CFL (Eq. 39, CFL=0.5). Optional artificial grain-compressibility EOS (cubic-per-particle) for pressure indeterminacy at rest.
- **DEM role:** none (μ(I) from da Cruz/GDR MiDi DEM historically; validated on experiments).
- **For us:** the cleanest effective-viscosity + cap recipe. **Limitations:** valid only dense (fails i>0.3, the collisional/impact regime); μ(I) has **no genuine solid branch** — at rest it's a high-viscosity fluid, not an elastoplastic solid.

#### 6. Szewc (2016/17), *Granular Matter* 19:3 — SPH granular column collapse · DOI 10.1007/s10035-016-0684-3
WCSPH column collapse in **2D and 3D**, but with a **Mohr–Coulomb/Bingham** effective viscosity (`τ_y = c + p tanφ`, μ = μ∞ + τ_y/γ̇, capped at μ_solid) — *not* a true μ(I). Wendland C2 kernel (kills pairing), Tait EOS (γ=7, c≥10 v_max), three CFL constraints.
- **DEM role — direct benchmark.** Compares SPH against **Utili et al. DEM** on deposit profiles, **energy partition vs. time**, and dissipated-energy-vs-aspect-ratio. Verdict: SPH matches DEM/experiment to first order but **drifts at high aspect ratio**, and the author attributes all error to the **rheology parameters** (μ∞, μ_solid are non-physical knobs) — i.e. the exact gap DEM calibration closes. Candid note: for *dry* granular, DEM is more accurate/efficient; SPH earns its keep in multiphase/continuum-coupled problems.
- **For us:** a fully-specified **column-collapse validation harness** (aspect-ratio sweep, 5 validation axes, runout scaling laws) and a clean **SPH-vs-DEM comparison methodology**. Lesson: **double precision is mandatory** (viscosity spans >5 orders). WCSPH free-surface pressure noise needs filtering.

#### 7. Salehizadeh & Shafiei (2019), *Granular Matter* 21:32 — μ(I) SPH column collapses · DOI 10.1007/s10035-019-0886-6
The **newest, most complete μ(I)-in-SPH** column-collapse/dam-break implementation. `η = μ(I)p/|γ̇|` (Eq. 38) with an explicit **regularization (Eq. 39): α_r=0.01, α_s=1e-6** (exponential Papanastasiou cap + singularity removal), plus a **velocity–pressure-coupling pressure stabilizer** (∇·v ↔ ∇²p) to kill WCSPH free-surface oscillations. Schwaiger corrected Laplacian + normalized kernel gradient for free-surface accuracy. EOS `K = 50 ρ_0 g H_max`, CFL=0.2.
- **Parameters (sand column):** μ_s=tan(30.5°), μ_2=tan(51.3°), I_0=2.65.
- **DEM role:** none (parameters from terrestrial experiments).
- **For us:** the most directly clonable Family-A skeleton, with the best-documented regularization. Same limitations as Minatti: local μ(I), no genuine solid, 2D, terrestrial g, no cohesion.

#### 8. Lagrée, Staron & Popinet (2011), *J. Fluid Mech.* 686:378 — μ(I) continuum vs DEM · DOI 10.1017/jfm.2011.335
**The reference DEM-vs-continuum validation** (continuum Navier–Stokes in Gerris, *not* SPH, but the closure transfers). `η = max(μ(I)p/(√2 D₂), 0)` with a high-viscosity **cap η_M ≈ 250 ρ√(gH³)** (arrest ≈ slow creep; results insensitive to η_M down to 1×).
- **DEM role — the whole point.** Ground truth = 2D **contact-dynamics** DEM (rigid disks, contact friction μ=0.5, restitution e=0.5), aspect ratios a=0.5–67.9. μ(I) reproduces outer shape, **inner deformation, and the static triangular core** to within **~1–2 grain diameters**, and recovers runout/height scalings. **Failure mode (critical):** μ(I) **systematically under-predicts runout in the deceleration/arrest phase (~10%, worse at high a)** because the thin, dilute, energetic front reaches **I = O(1)** — outside μ(I)'s validated regime; a kinetic-theory front would be needed.
- **For us:** defines the calibration/validation loop exactly — DEM ground truth, **error-map parameter optimization** (sweep μ_s, Δμ, I_0; minimize shape-averaged height error; pick one optimum across aspect ratios). Default fit: μ_s=0.32, Δμ=0.28, I_0=0.4. And a standing warning: expect the same under-spread at a footpad's dilute leading edge.

### Intrusion / landing-force models (the application)

#### 9. Agarwal, Karsai, Goldman & Kamrin (2021), *Sci. Adv.* 7:eabe0631 — Dynamic granular intrusion · DOI 10.1126/sciadv.abe0631
**The payoff paper for touchdown.** A **rate-independent** Drucker–Prager + tension-free separation continuum (Eqs. 1–2), with **macro-inertia (ρv̇) retained**, reproduces *dynamic* intrusion across speeds. Surprising result: **all rate dependence comes from macro-inertia** — no μ(I), no dynamic-friction drop needed.
- **Force decomposition (one static + two dynamic) — DRFT (Eq. 4):** `t = α(β,γ) H(−z̃)|z̃| − n λρ v_n²`, where
  - **static** `α(β,γ)|z|` = classic depth-linear RFT (the bearing term `K|z|`),
  - **dynamic inertial** `−n λρ v_n²` (λ ≈ O(1), ≈1–1.1 for plates) = velocity-squared momentum reaction,
  - **dynamic structural** `z̃ = z + δh`, `δh = r(rω²/g)` = free-surface remolding (matters for wheels/legs, **negligible for a symmetric footpad**).
- **For a vertical footpad:** `F(z, v_n) ≈ K|z̃| + λ ρ A v_n²`. The **ρv_n² term is what produces rapid-penetration → strong deceleration → arrest** (force ∝ v² brakes the leg; as v_n decays, force relaxes to the static K|z| term). **Ring-down is NOT modeled** — separation gives unloading, but the elastic↔plastic coupling that produces damped oscillation is a gap we must add.
- **DEM role:** none (MPM continuum vs. poppy-seed experiments + PIV).
- **For us:** the **DEM-calibration target** — design a footpad-intrusion DEM campaign to separately extract (i) static K at quasi-static rate, (ii) inertial λ from the v² scaling, (iii) structural δh. λ≈O(1) is a strong sanity check.

#### 10. Agarwal, Goldman & Kamrin (2023), *PNAS* 120:e2214017120 — 3D granular intrusion · DOI 10.1073/pnas.2214017120
Derives a **3D Resistive Force Theory** from the same Eq. 2 constitutive law (separation + DP yield + codirectional flow) via localization + dimensional analysis + isotropic representation theory. `F = ∫_surf α(n̂,v̂,ĝ)|z| ds`, with α = ρ_c g f̃(μ_int)·α^gen and α^gen = f_1 n̂ + f_2 v̂ + f_3 ĝ (three cubic-polynomial scalar functions). Fit from **~3000 full-field runs**.
- **Method note (important):** the full-field solver is **MPM, not SPH** — the "SPH" subject tag is topical. But the constitutive law is solver-agnostic and is what an SPH solver would integrate.
- **DEM role:** independent validation (drilling tests, ~6×10⁵–2.1×10⁶ grains); μ_int=0.21 extracted from **simple-shear DEM** — a direct example of the DEM→parameter step.
- **Validity bound (quote for touchdown):** quasi-static RFT requires `I ≲ 0.01` and **macro-inertial `I_mac = v/√(P/ρ) ≲ 0.15`**. Touchdown speeds may exceed this — exactly where the '21 dynamic ρv² term is required.
- **For us:** 3D-RFT is a **cheap analytic oracle** to cross-check SPH footpad forces (net force + moment + surface traction) after a tiny recalibration (the f_i and ξ_n = ρ_c g f̃). g enters explicitly → easy to rescale to lunar gravity.

---

## III. Synthesis — what this means for the codebase

### III.1 What the literature agrees you NEED to predict footpad touchdown force

1. **A pressure-proportional frictional shear strength** (Drucker–Prager / μ-P). Universal across all 10 papers. Sets the static bearing/penetration resistance.
2. **A tension-free / separation mechanism.** Mandatory in the Kamrin line (Eq. 1/2.7) and physically essential behind a fast footpad (ejecta, wake, unloading). Family A approximates it via "stress vanishes as material dilates"; Family B/C make it explicit.
3. **Macro-inertia in the momentum balance (ρv̇).** Per Agarwal '21, *this single term is the source of all rate dependence* — the `λ ρ A v_n²` contribution that brakes the leg. SPH retains it natively.
4. **A static depth term `K|z|`** (the RFT α-integral) — sufficient alone only in the quasi-static settling limit.

**μ(I) itself is sufficient-but-not-necessary for the force magnitude:** at intrusion pressures `I ≪ 0.1`, so a *constant*-μ Drucker–Prager already captures footpad force (Agarwal '21/'23). μ(I) becomes necessary only if the footpad reaches the dense-inertial regime (high-speed impact, large I) — then it sets velocity-strengthening of *shear* resistance, distinct from the bulk ρv² inertia.

### III.2 The constitutive fork — pick the core, then bolt on

| | **Family A: μ(I) viscosity** | **Family B: Kamrin elasto-VP + separation** | **Family C: Bui DP elastic-plastic** |
|---|---|---|---|
| Solid branch | ✗ (capped creep only) | ✓ real hypoelastic solid | ✓ real elastic solid |
| Separation/tension-free | weak (cap) | ✓ explicit (ρ_c) | via tension-cracking |
| Implementation cost | lowest (effective viscosity) | medium (return map) | medium (return map) |
| Free-surface/wake fidelity | poor at dilute front | best (separation built in) | needs artificial stress |
| Ports to SPH | trivially (viscous term) | "essentially unchanged" (point-local) | yes (canonical SPH-soil) |
| Data model | one σ via η | one σ + density scalar ρ | one σ, p from −I₁/3 |

**Recommendation: build the core on Family B (Dunatunga–Kamrin).** Reasons specific to *this* problem:
- It is the only family with a **genuine solid branch AND explicit tension-free separation** — exactly the two things footpad bearing + wake/ejecta demand.
- Its constitutive engine is **point-local and proven to port MPM→SPH unchanged** (consumes σ, D, W, ρ, Δt → returns σ).
- It already *uses* the μ(I) law, so a DEM-calibrated μ(I) drops straight in, and you can degrade to constant-μ for the quasi-static regime per Agarwal.
- Borrow Bui's **`p = −I₁/3` pressure-from-constitutive** convention and **artificial-stress / tension-cracking** for the SPH-specific tensile-instability work that DK's MPM grid sidesteps.
- Borrow the **Papanastasiou/η_max cap** (Minatti, Salehizadeh, Lagrée) as a *fallback* regularizer if the elastic branch chatters near ρ_c in WCSPH.

This is consistent with the earlier deep-read fork (Kamrin density-switch vs Wang two-stress): the 10-paper corpus leans the same way for touchdown, because the dominant dynamic effect is the inertial ρv² term and tension-free separation is non-negotiable.

### III.3 The DEM → SPH calibration workflow (concrete)

From GDR MiDi, Lagrée, and Agarwal, the pipeline is:

1. **Element-test DEM (the rheometer).** Run homogeneous **plane-shear** (and triaxial) DEM on your regolith grains, fixed-pressure control, periodic BCs, rigid-grain limit. Vary γ̇ and P.
2. **Form dimensionless groups & fit.** `I = γ̇ d/√(P/ρ_s)`, `μ_eff = τ/P`, `Φ = ρ/ρ_s`. Fit **μ(I) = μ_s + (μ_2−μ_s)/(I_0/I+1)** and **Φ(I)**; extract `ρ_c` (density at p→0) and bulk modulus `K`. Add a **cohesion intercept c** for regolith (a Drucker–Prager k_c term, absent from all 10 papers).
3. **Validation-loop DEM (the proving ground).** Run **granular column collapse** in DEM and in SPH; tune via the **Lagrée error-map** (minimize shape-averaged height error across aspect ratios). This is your first acceptance gate (Szewc's 5 axes: deposit shape, runout scaling, energy partition vs time, failure angle, pressure field).
4. **Footpad DEM campaign (the application calibration).** Per Agarwal '21, design DEM footpad-intrusion runs to *separately* extract the static `K`, inertial `λ` (from v² scaling), and (if relevant) structural `δh`. λ≈O(1) sanity-checks the continuum forces.
5. **3D-RFT oracle.** Fit Agarwal '23 3D-RFT (f_1,f_2,f_3, ξ_n) from the same data as a near-real-time cross-check on SPH footpad force/moment.

### III.4 Numerics consensus across the corpus

- **Kernel:** Wendland (C2/C4) — suppresses pairing/tensile instability (Minatti, Szewc). Cubic spline works with artificial stress (Bui).
- **Pressure:** WCSPH + Tait/linear EOS is the field default (Minatti, Szewc, Salehizadeh); `c ≥ 10 v_max`, γ=7. **OR** `p = −I₁/3` from the constitutive model (Bui) — preferred for Family B/C.
- **Low-γ̇ regularization:** mandatory. Papanastasiou cap (Minatti α_r=0.01/α_s=1e-6 in Salehizadeh) or η_max≈250ρ√(gH³) (Lagrée), or — for Family B — the elastic branch makes it moot.
- **Tensile instability:** artificial stress (Bui, n≈2.55, ε≈0.5) or tension-cracking (cohesionless) or particle shifting / stress regularization.
- **Free-surface pressure noise:** velocity–pressure coupling stabilizer (Salehizadeh) or filtering (Szewc); a known WCSPH headache — consider δ-SPH.
- **Time integration:** explicit leapfrog / symplectic Verlet; CFL ≈ 0.2–0.5; dual sound-speed + viscous constraint. **Double precision required.**
- **Boundaries:** dummy/ghost-particle walls (Morris/Adami); the **footpad is a moving rigid frictional contact** and deserves first-class treatment.

### III.5 Gaps the corpus does NOT cover (your contribution space)

These appear in *none* (or almost none) of the 10 papers and are exactly where lunar-regolith touchdown departs from the literature:

1. **Cohesion** — every paper is cohesionless. Regolith has van der Waals / vacuum cohesion → add a Drucker–Prager cohesion intercept `k_c` (Bui has the slot) calibrated from DEM with cohesive contacts.
2. **Low gravity** — all terrestrial g. g enters explicitly in I_mac, δh=r(rω²/g), and ξ_n=ρ_c g f̃ — rescalable but must be done deliberately.
3. **3D / axisymmetric** — nearly all demos are 2D plane-strain. Touchdown is 3D.
4. **Dynamic ring-down / rebound / repeated contact** — Agarwal models monotonic intrusion and steady locomotion only; the damped oscillatory leg response (your "landing stability") is unmodeled. This is the elastic-leg ↔ plastic-bed coupling you must add on top of the bed model.
5. **The dilute-front under-spread** (Lagrée) and the **i>0.3 collisional regime** (Minatti) — high-speed impact and ejecta exit μ(I)'s validated range; may need a kinetic-theory branch or accept the known ~10% bias.
6. **Nonlocal / quasi-static creep** — local μ(I) lacks a length scale; genuine slow quasi-static loading may need Kamrin–Koval/Henann–Kamrin nonlocal terms.

---

## IV. Consolidated parameter table (μ(I), as reported)

| Source | Material | μ_s | μ_2 | I_0 | ρ_s (kg/m³) | Notes |
|---|---|---|---|---|---|---|
| Jop '06 | glass beads | 0.382 (20.9°) | 0.643 (32.8°) | 0.279 | — | inclined-plane calibration |
| Dunatunga–Kamrin '15 | glass beads | 0.3819 | 0.6435 | (in ξ) | 2450 | ρ_c=1500, E=1GPa, ν=0.3 |
| Minatti–Paris '15 | sand | tan(repose) | μ_2/μ_s≈1.5 | 0.279 / 1.09 / 2.65 | 2600 | I_0 grows for finer sand |
| Salehizadeh '19 (beads) | glass beads | 0.382 | 0.643 | 0.279 | 2500 | φ=0.6 |
| Salehizadeh '19 (sand) | sand | tan(30.5°) | tan(51.3°) | 2.65 | 2600 | column collapse |
| Lagrée '11 (best fit) | DEM disks | 0.32 | 0.60 | 0.40 | — | Δμ=0.28; from DEM error-map |
| Szewc '16 | sand (M-C) | tan(φ), φ=30° | — | — (Bingham) | 2600 | not μ(I); c=0 |
| Bui '08 | Al bars / soil | φ=19.8°, c=0 | — | — (DP) | 1850–2650 | not μ(I); K≈0.7 MPa |
| Agarwal '23 (DEM) | poppy/grains | μ_int=0.21 | — | — (rate-indep) | 2470 | from simple-shear DEM; ρ_c=1310 |

*Regolith values to be obtained from our own DEM — these are anchors/sanity ranges, not regolith parameters.*

---

## V. Citation index

| # | Short | DOI | Family / role |
|---|---|---|---|
| 1 | GDR MiDi 2004 | 10.1140/epje/i2003-10153-0 | Foundation (I-number, DEM rheometry) |
| 2 | Jop, Forterre, Pouliquen 2006 | 10.1038/nature04801 | Foundation (μ(I) closure) |
| 3 | Dunatunga & Kamrin 2015 | 10.1017/jfm.2015.383 | **B — recommended core** |
| 4 | Bui et al. 2008 | 10.1002/nag.688 | C — SPH-soil backbone, artificial stress |
| 5 | Minatti & Paris 2015 | 10.1016/j.apm.2014.05.034 | A — μ(I)→SPH viscosity + Papanastasiou cap |
| 6 | Szewc 2016/17 | 10.1007/s10035-016-0684-3 | A — SPH-vs-DEM validation harness |
| 7 | Salehizadeh & Shafiei 2019 | 10.1007/s10035-019-0886-6 | A — newest μ(I)-SPH + regularization |
| 8 | Lagrée, Staron, Popinet 2011 | 10.1017/jfm.2011.335 | DEM-vs-continuum validation + error-map |
| 9 | Agarwal et al. 2021 | 10.1126/sciadv.abe0631 | Intrusion force decomposition (DRFT) |
| 10 | Agarwal et al. 2023 | 10.1073/pnas.2214017120 | 3D-RFT oracle |

*Generated 2026-06-16 from full reads of all 10 PDFs in `./papers`.*
