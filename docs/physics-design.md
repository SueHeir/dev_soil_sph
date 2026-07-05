# Solver Physics Design — DEM-Informed Granular SPH (v0, glass beads)

**Purpose.** Specify the *physics* the SPH solver integrates: governing equations, the constitutive model, the stress-update algorithm, and the regularizations — at enough precision to implement and unit-test. Software architecture (data structures, neighbor search, parallelism) is a separate doc; this one defines what each particle must compute.

**Scope of v0.** Dry, cohesionless **glass beads**, terrestrial gravity, the constitutive **Family B** core (Dunatunga–Kamrin elasto-viscoplastic + density-based tension-free separation) cast for SPH, validated on granular column collapse. Deliberately deferred (see §10): cohesion, lunar gravity, full 3D production, leg ring-down, the dilute/collisional front.

**Provenance.** Every choice below traces to `docs/literature-review.md`. Key sources: Dunatunga–Kamrin 2015 (core), Jop 2006 / GDR MiDi 2004 (μ(I)), Bui 2008 (SPH-soil stress treatment, tensile instability), Minatti–Paris 2015 / Salehizadeh 2019 / Lagrée 2011 (regularization), Agarwal 2021 (intrusion force, validity bound).

---

## 1. Modeling decisions (the commitments)

| Decision | Choice (v0) | Why |
|---|---|---|
| Constitutive family | **B — elasto-viscoplastic + density separation** | Only family with a *real solid branch* AND *explicit tension-free separation*; engine is point-local, ports MPM→SPH unchanged. |
| Pressure | **Weakly-compressible, p from a granular EOS p(ρ)** | The EOS *is* the separation criterion (p=0 below ρ_c). Decouples volumetric (density) from deviatoric (shear) response — clean for SPH. |
| Deviatoric stress | **Hypoelastic rate + μ(I) Drucker–Prager return map** | Gives the quasi-static elastic solid AND the dense flow in one local update. |
| Shear rheology | **μ(I)** (Jop form), reducible to constant-μ | μ(I) needed only at high I; intrusion runs at I≪0.1 (Agarwal). |
| Kernel | **Wendland C2** | Suppresses pairing/tensile instability without artificial stress. |
| Tensile-instability backup | **Bui artificial stress / tension-cracking** | Cohesionless beads: tension-cracking usually enough; keep artificial stress available. |
| Time integration | **Explicit, symplectic (kick–drift–kick)** | Standard, momentum-friendly, GPU-friendly. |
| Precision | **float64 mandatory** | Viscosity/stiffness span >5 orders (Szewc). |

---

## 2. Governing equations (Lagrangian continuum)

Per material particle, with material derivative $D/Dt$ (advection is free in SPH):

**Mass (continuity form):**
$$\frac{D\rho}{Dt} = -\rho\,\nabla\!\cdot\!\mathbf{v}$$

**Momentum:**
$$\frac{D\mathbf{v}}{Dt} = \frac{1}{\rho}\nabla\!\cdot\!\boldsymbol{\sigma} + \mathbf{g}, \qquad \boldsymbol{\sigma} = -p\,\mathbf{I} + \mathbf{s}$$

where $\mathbf{s}$ is the deviatoric stress. **Macro-inertia (the $\rho\,D\mathbf v/Dt$ term) is retained** — per Agarwal 2021 it is the sole source of rate dependence in intrusion (the $\lambda\rho v_n^2$ braking force emerges from it automatically; we do not add it by hand).

Sign convention: **compression positive** for $p$ (i.e. $p=-\tfrac13\,\mathrm{tr}\,\boldsymbol\sigma$); tension is $p<0$ and is not sustained (separation).

Kinematics from the SPH velocity gradient $\mathbf{L}=\nabla\mathbf{v}$:
$$\mathbf{D} = \tfrac12(\mathbf{L}+\mathbf{L}^\top)\ \text{(strain rate)}, \qquad \mathbf{W} = \tfrac12(\mathbf{L}-\mathbf{L}^\top)\ \text{(spin)}, \qquad \mathbf{D}' = \mathbf{D} - \tfrac13(\mathrm{tr}\,\mathbf{D})\mathbf{I}$$

---

## 3. Constitutive model

The state carried per particle is $(\rho,\ \mathbf{s})$ — density and deviatoric stress. Pressure is derived from $\rho$; full stress is reassembled each step.

### 3.1 Three phases (the master switch)

$$\boldsymbol{\sigma} =
\begin{cases}
\mathbf{0} & \rho < \rho_c \quad\text{(disconnected / separated — stress-free)}\\[4pt]
-p(\rho)\,\mathbf{I} + \mathbf{s} & \rho \ge \rho_c \quad\text{(dense: elastic solid OR μ(I) flow)}
\end{cases}$$

Below the critical (close-packed) density $\rho_c$ the grains lose contact and carry **no stress** — the tension-free separation that gives ejecta, the wake behind a footpad, and clean unloading. This is one `if`.

### 3.2 Pressure — granular EOS (and the separation criterion)

Weakly-compressible, tied to density (Dunatunga–Kamrin Eq. 2.7 form):
$$p(\rho) = \begin{cases} 0 & \rho < \rho_c \\[3pt] K\left(\dfrac{\rho}{\rho_c}-1\right) & \rho \ge \rho_c \end{cases}$$

$K$ = effective bulk modulus, chosen for weak compressibility: sound speed $c_s=\sqrt{K/\rho_c}\ge 10\,v_\text{max}$ (Mach ≲ 0.1), so $K \gtrsim 100\,\rho_c v_\text{max}^2$. **Do not use the true glass modulus** — it would force a crippling timestep; use the smallest $K$ that keeps density fluctuations ≲1% (rigid-grain-limit insensitivity, GDR MiDi). The $\rho<\rho_c\Rightarrow p=0$ branch is the separation trigger and must be evaluated *before* the deviatoric update.

### 3.3 Deviatoric stress — hypoelastic predictor + μ(I) return map

**Inertial number** (rate sensitivity; $d$ = grain diameter, $\rho_s$ = solid grain density):
$$I = \dot\gamma^{p}\,\frac{d\sqrt{\rho_s}}{\sqrt{p}}$$

**Friction law (Jop μ(I)):**
$$\mu(I) = \mu_s + \frac{\mu_2-\mu_s}{I_0/I + 1}, \qquad \tau_\text{yield} = \mu(I)\,p$$

with equivalent shear stress $\bar\tau = \sqrt{\tfrac12\,\mathbf{s}\!:\!\mathbf{s}}$. No plastic flow while $\bar\tau \le \mu_s p$ (the elastic, quasi-static solid). The yield is Drucker–Prager (pressure-dependent).

**Return-mapping update** (per particle, per step — Dunatunga–Kamrin §3.2, the cancellation-safe form). Let $G$ = elastic shear modulus, $\Delta t$ = step:

1. **Density & pressure first.** Update $\rho$ (continuity), get $p$ from §3.2. If $\rho<\rho_c$ or $p\le 0$: set $\mathbf{s}=\mathbf 0$, $\boldsymbol\sigma=\mathbf 0$, **return** (disconnected).
2. **Jaumann elastic trial deviator** (objective rate keeps it frame-indifferent):
$$\mathbf{s}^{\text{tr}} = \mathbf{s}^{n} + \Delta t\left(2G\,\mathbf{D}' + \mathbf{s}^{n}\mathbf{W} - \mathbf{W}\mathbf{s}^{n}\right), \qquad \bar\tau^{\text{tr}}=\sqrt{\tfrac12\mathbf{s}^{\text{tr}}\!:\!\mathbf{s}^{\text{tr}}}$$
3. **Yield check.** If $\bar\tau^{\text{tr}} \le \mu_s\,p$: **elastic**, $\mathbf{s}^{n+1}=\mathbf{s}^{\text{tr}}$.
4. **Else plastic.** With $\xi = I_0/(d\sqrt{\rho_s})$, define
$$S_0=\mu_s p,\quad S_2=\mu_2 p,\quad \alpha=\xi G\,\Delta t\sqrt{p},\quad B=S_2+\bar\tau^{\text{tr}}+\alpha,\quad H=S_2\,\bar\tau^{\text{tr}}+S_0\,\alpha$$
$$\boxed{\ \bar\tau^{n+1} = \frac{2H}{B+\sqrt{B^2-4H}}\ }\qquad \dot\gamma^{p}=\frac{\bar\tau^{\text{tr}}-\bar\tau^{n+1}}{G\,\Delta t}$$
   then **radially rescale** the deviator: $\mathbf{s}^{n+1} = \dfrac{\bar\tau^{n+1}}{\bar\tau^{\text{tr}}}\,\mathbf{s}^{\text{tr}}$.
5. **Reassemble** $\boldsymbol\sigma = -p\,\mathbf{I} + \mathbf{s}^{n+1}$.

The boxed root is the numerically stable form (avoids catastrophic cancellation near $\mu\to\mu_2$). This update is **fully local** — no neighbor coupling — so it parallelizes trivially and is unit-testable in isolation (see §8).

### 3.4 Low-I regularization (fallback)

With a genuine elastic branch (§3.3), the $\eta\to\infty$ singularity of pure-viscosity μ(I) models **does not arise** — below yield the material is an elastic solid, not an infinite-viscosity fluid. If, in WCSPH, the deviator chatters near $\rho_c$ or at vanishing $\bar\tau$, apply the **Papanastasiou cap** as a fallback (Minatti Eq. 48 / Salehizadeh Eq. 39, $\alpha_r=0.01,\ \alpha_s=10^{-6}$), or cap an effective viscosity at $\eta_M\approx 250\,\rho\sqrt{gH^3}$ (Lagrée). Treat this as a numerical safety net, not the primary mechanism.

---

## 4. Stress into the SPH momentum equation

Symmetric, momentum-conserving form (full stress tensor):
$$\frac{D\mathbf{v}_i}{Dt} = \sum_j m_j\left(\frac{\boldsymbol\sigma_i}{\rho_i^2}+\frac{\boldsymbol\sigma_j}{\rho_j^2}\right)\!\cdot\nabla_i W_{ij} + \mathbf{g} \;+\; (\text{artificial viscosity } \Pi_{ij}) \;+\; (\text{artificial stress } R_{ij})$$

- **Artificial viscosity** $\Pi_{ij}$ (Monaghan): small coefficient for numerical stability/shock damping; physical dissipation is mostly carried by friction. Start $\alpha_\Pi\!\approx\!0.1$.
- **Artificial stress** $R_{ij}$ (Bui Eqs. 63–69): activate only if pairing appears. Repulsive factor $f_{ij}=(W_{ij}/W(\Delta d,h))^n$, $n\approx2.55$, $\varepsilon\approx0.5$. For cohesionless beads, **tension-cracking** (force $p\to0$ where the trial pressure goes tensile) is usually sufficient and is already implied by §3.1.

---

## 5. Boundaries and the footpad

- **Walls / container:** dummy/ghost particles (Morris/Adami), mirrored pressure, no-slip or free-slip as needed.
- **Footpad (the load path):** a **moving rigid frictional contact**, not a default no-slip wall — the touchdown force is read off here. v0: rigid intruder represented by boundary particles with prescribed kinematics; Coulomb surface friction $\mu_\text{surf}\le\mu_\text{int}$. Record net force/moment on the footpad each step (this is the deliverable observable).
- **Validity flag:** compute the macro-inertial number $I_\text{mac}=v/\sqrt{p/\rho}$ at the footpad. Agarwal: quasi-static RFT holds for $I_\text{mac}\lesssim0.15$; above it the $\rho v_n^2$ inertial term dominates (the solver captures it natively, but flag the regime in output).

---

## 6. Time integration & stability

- **Scheme:** kick–drift–kick (leapfrog) / symplectic Verlet.
- **Timestep** = min of: CFL/sound $\Delta t\le 0.25\,h/c_s$; body-force $\Delta t\le 0.25\sqrt{h/|\mathbf a|}$; viscous (if cap active) $\Delta t\le 0.125\,\rho h^2/\eta_\text{max}$.
- **CFL coefficient** 0.2–0.5 (Bui 0.2; Minatti 0.5).
- Free-surface WCSPH pressure noise: velocity–pressure-coupling stabilizer (Salehizadeh) or δ-SPH diffusion; clip tensile $p$ at the surface.

---

## 7. Glass-bead parameter set (v0 starting values)

These are **literature anchors** to bring the solver up; v1 replaces them with our DEM-fitted values (see `docs/dem-campaign.md`).

| Symbol | Meaning | v0 value | Source |
|---|---|---|---|
| $\rho_s$ | solid grain density | 2500 kg/m³ | glass |
| $d$ | grain diameter | 0.5 mm | typical |
| $\Phi_c$ | critical packing fraction | 0.60 | Jop/DK |
| $\rho_c$ | critical density $=\Phi_c\rho_s$ | 1500 kg/m³ | — |
| $\mu_s$ | static friction coeff | 0.38 (≈20.9°) | Jop |
| $\mu_2$ | limiting friction coeff | 0.64 (≈32.8°) | Jop |
| $I_0$ | inertial-number scale | 0.28 | Jop |
| $\nu$ | Poisson ratio | 0.3 | — |
| $K$ | effective bulk modulus | set via $c_s\ge10v_\text{max}$ | WCSPH |
| $G$ | shear modulus $=\tfrac{3(1-2\nu)}{2(1+\nu)}K$ | derived | — |
| $g$ | gravity | 9.81 m/s² (v0) | terrestrial first |

---

## 8. Validation plan (acceptance gates)

1. **Unit test — single-particle stress update.** Drive the §3.3 return map with a prescribed $\mathbf{L}(t)$ (shear, then compression, then extension) and compare against a 4th-order Runge–Kutta integration of the same ODE (Dunatunga–Kamrin App. B). Must show 1st-order convergence in $\Delta t$. *This is the first thing to build and the cheapest to verify.*
2. **Bagnold profile** — steady flow down an incline reproduces the analytic $u(z)\propto[1-(1-z/H)^{3/2}]$; no flow below $\arctan\mu_s$, no steady state above $\arctan\mu_2$.
3. **Granular column collapse** — the main gate, against DEM + experiment (Szewc's 5 axes): deposit shape, runout scaling $(r_\infty-r_0)/r_0$ vs aspect ratio $a$, energy partition vs time, failure-plane angle, pressure field. Calibrate via the **Lagrée error-map** (minimize shape-averaged height error across $a$). Expect ~1–2 grain-diameter agreement in the bulk and a **known ~10% runout under-spread at the dilute front** ($I=O(1)$) — accept it in v0.
4. **Footpad intrusion** — force–depth–velocity vs DEM; cross-check against 3D-RFT (Agarwal 2023) once fit.

---

## 9. Data-model implication (one line for the architecture doc)

Each particle stores $(\mathbf{x},\mathbf{v},\rho,\mathbf{s})$ + mass. Pressure and full stress are derived. The constitutive update (§3.3) is a pure function `(s_n, D, W, ρ, Δt, params) → (s_{n+1}, σ)` behind an interface — so swapping Family B for a Family-A μ(I)-viscosity or adding cohesion later touches *only that function*, not the solver core.

---

## 10. Deferred (v1+), with the hook for each

| Deferred | Hook already in place |
|---|---|
| **Cohesion** (regolith) | Add intercept $k_c$ to the yield: $\bar\tau\le\mu_s p + k_c$; calibrate from cohesive-contact DEM. |
| **Lunar gravity** | $\mathbf g$ is a parameter; note $g$ enters $I_\text{mac}$ and any structural $\delta h$. |
| **Full 3D production** | Tensors are written 3D; v0 demos may be 2D plane-strain. |
| **Leg ring-down / rebound** | Couple footpad boundary to an external elastic-leg model; the bed already unloads via separation (§3.1). |
| **Dilute/collisional front ($I>0.3$)** | Add a kinetic-theory branch above an $I$ threshold; until then, document the front bias. |
| **Φ(I) dilatancy transients** | v0 is critical-state/isochoric; add a $\Phi(I)$ evolution law when DEM provides it. |

---

## 11. Planned v1 — two-branch stress + granular temperature (de-fluidization)

**Why.** The landing-critical physics is *memory*: regolith arrives plume-fluidized
(hot, dilated, high granular temperature `T`, low `Φ`) and **consolidates** —
bearing capacity and sinkage are time-dependent as `T → 0` and the contact stress
ramps up. The v0 μ(I)/Drucker–Prager model is memoryless and cannot produce this
transient. Granular temperature `T` is the state variable that remembers prior
agitation. (Project direction, 2026-06-16; KT closure already in
`docs/dem-lebc-kt-spec.md`.)

### 11.1 Two-branch stress

$$\boldsymbol{\sigma} = \boldsymbol{\sigma}_{\mathrm{KT}}(T,\Phi) + \boldsymbol{\sigma}_{\mathrm{contact}}(\Phi)$$

- **Collisional (kinetic-theory) branch** — agitated, T-dependent:
  $$p_{\mathrm{KT}} = \rho_s\,\Phi\,T\big[\,1 + 2(1+e)\,\Phi\,g_0(\Phi)\,\big], \qquad
  \boldsymbol{\tau}_{\mathrm{KT}} = \eta_{\mathrm{KT}}(\Phi,T)\,\dot{\boldsymbol{\gamma}}$$
  with `g_0` the Carnahan–Starling pair correlation. Known from KT (Lun/Garzó–Dufty),
  validated on the LEBC rig.
- **Enduring-contact branch** — the calibrated residual `σ_contact = p_DEM − p_KT`
  vs `Φ`, switching on above a crossover `Φ_c`. This is the v0 contact stress (μ(I)
  Drucker–Prager) recast as a function of `Φ`. **Needs the vf-grid DEM campaign.**

`Φ = ρ/ρ_s` is already available from the carried (continuity) density — no new
volumetric field.

### 11.2 Granular-temperature balance

$$\frac{DT}{Dt} = \underbrace{\mathcal{P}}_{\text{production}} - \underbrace{\Gamma}_{\text{dissipation}} + \underbrace{\nabla\!\cdot(\kappa\,\nabla T)}_{\text{conduction}}$$

- **Production** `P` — collisional stress working against shear, local
  (`∝ σ_KT : D`); uses the velocity gradient already computed in pass 1.
- **Dissipation** `Γ` — inelastic cooling `∝ γ*(Φ,e)\,ρ_s Φ² g_0 (1-e²) T^{3/2}/d`,
  local. **Known/validated** (LEBC `T*` + `bench_*_haff_cooling`).
- **Conduction** `∇·(κ∇T)` — an SPH Laplacian (neighbor gather, Brookshaw/Cleary
  form with harmonic-mean `κ`). The coefficient `κ(Φ,e)` is the **one missing DEM
  measurement** — needs the inhomogeneous (boundary-heated) DIRT rig.

### 11.3 dev_soil_sph code impact (modular)

- **New `MudAtom` column** `temperature` (`#[forward]` — neighbors read `T` for the
  conduction Laplacian and for their `σ_KT`), plus a `#[zero]` accumulator for the
  conduction gather. `T` is persistent state, integrated by §11.2.
- **One new neighbor pass** (conduction Laplacian) + **one per-particle T
  integration** (production − dissipation + conduction), slotted into `PreForce`
  between the density pass and the constitutive update. `T` joins the forwarded
  columns at the mid-step halo.
- **`mud_constitutive::update_stress` becomes two-branch** — `σ_KT(T,Φ) +
  σ_contact(Φ)` — behind the same pure-function interface (the architecture's payoff:
  this swap touches only the constitutive crate + adds the T systems).

### 11.4 Sequencing & first dev_soil_sph milestone

- **Can build now** (theory-based, no new DEM): `σ_KT`, production, dissipation. The
  **homogeneous Haff-cooling test** (uniform `T_0`, no shear/gravity → `T(t)` decays
  as Haff's law) needs *neither* `κ` *nor* `σ_contact` → the natural first dev_soil_sph-side
  `T` milestone, analogous to rest_state/hydrostatic.
- **Waits on DEM:** `σ_contact(Φ)` (vf campaign, running) and `κ(Φ,e)` (inhomogeneous
  rig, to be built). These swap in like `glass_beads_v0` did for μ(I).
- **Headline validation:** a DEM de-fluidization transient (hot+dilated bed cools and
  consolidates) reproduced by dev_soil_sph — `T`-decay, `Φ`-rise, contact-stress hand-off.

Active-plume / two-phase gas coupling is a later, larger step; dry de-fluidization
*aftermath* is the tractable, landing-relevant first target.

---

*Drafted 2026-06-16. §11 added 2026-06-16. Companion: `docs/literature-review.md`, `docs/sph-primer.md`, `docs/dem-campaign.md`, `docs/dem-lebc-kt-spec.md`.*
