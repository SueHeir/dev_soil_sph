# DEM Campaign Spec — Glass-Bead Spheres (v0, to calibrate & validate the SPH model)

**Purpose.** Generate the micro-scale ground truth that (1) **calibrates** the continuum constitutive law the SPH solver integrates, (2) **validates** the SPH solver, and (3) **anchors** the footpad force law. First pass uses **monodisperse-ish glass-bead spheres** — deliberately the simplest, best-characterized material, and the one the μ(I) literature (Jop, GDR MiDi, Lagrée) is calibrated on, so every result is cross-checkable.

**Why spheres / glass beads first.** Spheres roll, so DEM spheres give a *lower* macroscopic friction $\mu_s$ than angular regolith — that is a feature here, not a bug: it keeps us inside the validated μ(I) regime and lets us confirm the whole DEM→SPH pipeline before adding the complications (angularity, cohesion, polydispersity, low-g) that real regolith needs. See `docs/physics-design.md` §10 for the deferred list.

**Output the solver needs (close the loop):** $\mu_s,\ \mu_2,\ I_0$ (the μ(I) fit), $\Phi(I)$, $\rho_c$, effective bulk modulus $K$ — plus column-collapse and footpad reference data. These map one-to-one onto `physics-design.md` §7.

---

## 1. The three tiers

| Tier | Purpose | DEM configuration | Feeds |
|---|---|---|---|
| **1. Element tests** | Derive the constitutive law | Homogeneous **simple shear**, controlled pressure | $\mu(I),\ \Phi(I),\ \rho_c,\ K$ |
| **2. Column collapse** | Validate the SPH solver | Free-surface collapse, aspect-ratio sweep | runout/shape/energy gates |
| **3. Footpad intrusion** | Anchor the application | Rigid intruder, vertical penetration | static $K_\text{bear}$, inertial $\lambda$ |

Tier 1 is the rheometer; Tier 2 is the proving ground; Tier 3 is the actual touchdown observable. Build them in that order.

---

## 2. DEM model specification (shared by all tiers)

### 2.1 Particles
- **Shape:** spheres.
- **Polydispersity:** uniform in radius, **±10–20%** about $d=0.5$ mm (i.e. $d\in[0.45,0.55]$ mm). *Mild polydispersity is mandatory* — perfectly monodisperse spheres crystallize and give spurious rheology.
- **Solid density:** $\rho_s = 2500\ \text{kg/m}^3$.

### 2.2 Contact model
Hertz–Mindlin (nonlinear elastic) with tangential history and Coulomb cap — the standard for spheres. (Linear spring-dashpot is an acceptable cheaper alternative; results in the dense regime are insensitive per GDR MiDi.)

| Parameter | Symbol | Value (v0) | Notes |
|---|---|---|---|
| Grain Young's modulus | $E$ | **$5\times10^{7}$–$10^{8}$ Pa** (reduced) | NOT real glass (70 GPa) — reduced for timestep. Verify rigid-grain-limit insensitivity (§5). |
| Poisson ratio | $\nu$ | 0.25 | glass-like |
| Restitution | $e$ | 0.7 | dense rheology ~insensitive (GDR MiDi); affects only dense→collisional I |
| **Interparticle friction** | $\mu_p$ | **0.5** | THIS sets macroscopic $\mu_s$ — the one contact param that matters most |
| Rolling friction | $\mu_r$ | 0.0 (v0) | spheres; add small $\mu_r$ later to mimic angularity |
| Cohesion | — | **none (v0)** | glass beads dry; cohesion is a v1 regolith addition |

### 2.3 Numerics
- **Timestep:** $\Delta t \le 0.15\,t_R$, Rayleigh time $t_R = \pi R\sqrt{\rho_s/G_\text{grain}}/(0.1631\nu+0.8766)$. Re-derive when $E$ changes.
- **Gravity:** $g=9.81\ \text{m/s}^2$ (terrestrial, v0) — keeps results comparable to the glass-bead literature; swap to $1.62\ \text{m/s}^2$ only in the lunar phase.
- **Suggested codes:** LIGGGHTS / LAMMPS-granular, MercuryDPM, Yade, or Project Chrono. Any with Hertz–Mindlin + Lees–Edwards (Tier 1) + arbitrary rigid intruders (Tier 3). *(Pick one; LIGGGHTS and MercuryDPM are the common granular choices — flag your preference and I'll tailor the input-deck details.)*

---

## 3. Tier 1 — Element tests (the rheometer)

**Goal:** measure $\mu_\text{eff}(I)$ and $\Phi(I)$, extract $\mu_s,\mu_2,I_0,\rho_c,K$.

### 3.1 Configuration — homogeneous simple shear, controlled pressure
The canonical da Cruz / GDR MiDi rheometer:
- **Box:** fully periodic in flow ($x$) and vorticity ($z$); shear imposed via **Lees–Edwards** boundaries (or moving rough walls in $y$).
- **Control mode:** **fixed normal stress (pressure-controlled)** — the top boundary moves to hold $\sigma_{yy}=P$ constant, so the sample dilates/compacts freely and $\Phi$ is an *output*. (This is preferred over fixed-volume because it directly gives $\Phi(I)$ and $\rho_c$.)
- **Size:** $\ge 20\,d$ in each direction; target $\sim 10^4$ particles (2D-thin) to $\sim 10^5$ (fully 3D). Confirm size-independence by doubling.
- **Measurement:** after steady state, time-average the stress tensor (Love–Weber) → $\tau=\sigma_{xy}$, $P=\tfrac12(\sigma_{xx}+\sigma_{yy})$ (or $\sigma_{yy}$), and $\Phi$.

### 3.2 The sweep (build the μ(I) curve)
Each run is one $(P,\dot\gamma)$ point giving one $I$:
$$I = \dot\gamma\,d\big/\sqrt{P/\rho_s}, \qquad \mu_\text{eff}=\tau/P, \qquad \Phi=\rho/\rho_s$$

Target **$I$ from $\sim10^{-4}$ (quasi-static) to $\sim0.5$ (dense-inertial)** — at least **5 points per decade**, ~20–25 runs. Reach low $I$ by raising $P$ and/or lowering $\dot\gamma$. Suggested grid:

| $P$ (Pa) | $\dot\gamma$ (s⁻¹) sweep | $I$ range covered |
|---|---|---|
| 1000 | 1 → 300 | ~$3\times10^{-4}$ → ~0.1 |
| 100 | 1 → 300 | ~$10^{-3}$ → ~0.3 |
| 10 | 3 → 300 | ~$10^{-2}$ → ~0.5 |

*(Exact $\dot\gamma$ values set after a pilot run fixes the $I(P,\dot\gamma)$ mapping for the chosen $d,\rho_s$.)*

### 3.3 Fits / deliverables
- **μ(I):** fit $\mu(I)=\mu_s+(\mu_2-\mu_s)/(I_0/I+1)$ → report $\mu_s,\mu_2,I_0$ (+CIs). Sanity anchor: glass beads ≈ (0.38, 0.64, 0.28); spheres may read a bit lower.
- **Φ(I):** fit $\Phi(I)=\Phi_\text{max}-(\Phi_\text{max}-\Phi_\text{min})\,I$ (or a saturating form). $\Phi_\text{max}=\Phi(I\!\to\!0)$.
- **$\rho_c$:** $\rho_c=\Phi_c\rho_s$ with $\Phi_c=\Phi_\text{max}$ (the $P\to0$ / loosest sustaining packing).
- **$K$:** separate **isotropic compression test** (no shear, ramp $P$, measure $\Delta\Phi$) → effective bulk modulus for the WCSPH EOS. (Only an *effective* $K$ is needed; the WCSPH $K$ will be chosen for Mach number, not grain stiffness — this DEM $K$ is a sanity bound.)

### 3.4 Robustness checks (cheap, high-value)
- **Insensitivity sweep:** repeat 2–3 $I$ points at $e\in\{0.5,0.9\}$ and $\mu_p\in\{0.3,0.7\}$. Expect $\mu(I),\Phi(I)$ ~unchanged except the dense→collisional knee (confirms GDR MiDi → loose tolerance on $e$, tight on $\mu_p$).
- **Stiffness check:** one $I$ point at $10\times E$ → confirms rigid-grain limit (rheology independent of grain stiffness).

### 3.5 Kinetic-theory validation (frictionless sub-sweep)

For **smooth (frictionless) inelastic spheres** the LEBC stress tensor is analytic from granular kinetic theory (collisional regime), so it validates the DEM stress-measurement pipeline before we trust it in the dense regime where there is no analytic check. Additions to Tier 1:
- **Record granular temperature** $T = \tfrac13\langle |\mathbf v_i - \bar{\mathbf v}(y_i)|^2\rangle$ (subtract the mean shear profile $\bar v_x = \dot\gamma y$) and the **full stress tensor** (including normal-stress differences $N_1, N_2$), not just $\tau$ and $P$.
- **Frictionless sub-sweep** ($\mu_p = 0$, $e \in \{0.7, 0.9\}$) compared to the KT closure $p = \rho_s\Phi T[1 + 2(1+e)\Phi g_0]$, $\sigma_{xy} = \eta\dot\gamma$, with the steady-shear energy balance $\sigma_{xy}\dot\gamma = \Gamma$ fixing $T$ (Bagnold scaling). Agreement at moderate $\Phi$; deviation as $\Phi\to\Phi_c$ marks the dense regime.
- KT also anchors the **collisional branch** ($I\gtrsim 0.3$) that $\mu(I)$ cannot represent (the fast-impact/ejecta regime, `docs/literature-review.md`) and predicts the $N_1,N_2$ that the MUD $\mu(I)$ model omits.

Full implementation brief (LEBC primitive, measurement, KT formulas): **`docs/dem-lebc-kt-spec.md`**.

---

## 4. Tier 2 — Column collapse (SPH validation)

**Goal:** the reference dataset the SPH solver must reproduce (`physics-design.md` §8 gate 3).

- **Setup:** rectangular column of beads, width $r_0$, height $H_0$, released onto a rough rigid floor (glued-grain roughness). Side wall removed instantaneously (or gate lifted ≥ collapse speed).
- **Aspect-ratio sweep:** $a=H_0/r_0 \in \{0.5,\ 1,\ 2,\ 3,\ 6\}$ — spans both runout regimes (transition near $a\approx2$–3). Add $a=0.5,6$ to bracket.
- **Size:** $\gtrsim 10^4$–$10^5$ grains so the deposit is many grains thick.
- **Measure (the 5 SPH gates):**
  1. deposit profile $h(x)$ vs time and final;
  2. **runout scaling** $(r_\infty-r_0)/r_0$ vs $a$ (expect $\sim1.2\,a$ then $\sim1.6\,a^{1/2}$, Lube/Staron);
  3. **energy partition** (kinetic/potential/dissipated) vs dimensionless time (Szewc/Utili template);
  4. failure-plane inclination;
  5. internal velocity field (for the static-core / Bagnold-front structure).
- **Use:** feed (1)+(2) into the **Lagrée error-map** to lock the SPH μ(I) parameters (minimize shape-averaged height error across all $a$). The DEM here is ground truth; the SPH must land within ~1–2 grain diameters in the bulk.

---

## 5. Tier 3 — Footpad intrusion (application anchor)

**Goal:** the touchdown force law $F(z,v_n)\approx K_\text{bear}|z| + \lambda\rho A v_n^2$ (Agarwal 2021), and its two coefficients.

- **Intruder:** rigid, vertically driven. v0 shapes: **flat circular plate** and **sphere** (the two cleanest; a footpad is between them). Diameter $\gg d$ (≥ 20–40 grains) so it's a continuum-scale intruder.
- **Bed:** deep relative to intruder ($\ge 10\times$ intruder size; below that, lift saturates — Agarwal depth bound), prepared at a known $\Phi$ (loose and dense states).

### 5.1 Two coefficient-extraction campaigns
1. **Quasi-static (static term $K_\text{bear}$):** drive at low speed (target $I_\text{mac}=v/\sqrt{P/\rho}\lesssim 0.05$). Measure force vs depth $F(z)$; the slope gives $K_\text{bear}$ (the $\int\alpha H(-z)|z|$ RFT term).
2. **Velocity sweep (inertial term $\lambda$):** fixed shallow depth band, sweep impact speed $v$ across $I_\text{mac}\sim 0.05\to 0.5+$. Fit the *excess* force above quasi-static to $\lambda\rho A v^2$ → extract $\lambda$. **Expect $\lambda\approx O(1)$** (≈1–1.1 for plates) — a strong sanity check on the whole chain.

### 5.2 Deliverables
- $K_\text{bear}(\Phi)$, $\lambda$, and raw $F(z,v)$ tables.
- The **two-stage signature** (rapid penetration → oscillatory attenuation): record full force–time and penetration–time so the SPH (and later the coupled leg model) can be checked against the *transient*, not just steady values.
- These also seed the Agarwal-2023 **3D-RFT oracle** fit ($\xi_n=\rho_c g\tilde f(\mu_\text{int})$, $f_1,f_2,f_3$).

---

## 6. Data products (what to hand the SPH side)

A single calibration file (JSON/YAML) the solver reads at startup:
```yaml
material: glass_beads_v0
grain: {d_mean: 0.5e-3, d_spread: 0.2, rho_s: 2500}
mu_I:  {mu_s: ___, mu_2: ___, I_0: ___}      # Tier 1
phi_I: {phi_max: ___, phi_min: ___}          # Tier 1
rho_c: ___                                    # Tier 1  (= phi_max * rho_s)
K_eff: ___                                    # Tier 1 compression (sanity bound)
contact_used: {E: ___, nu: 0.25, e: 0.7, mu_p: 0.5}
footpad: {K_bear: ___, lambda: ___}           # Tier 3
provenance: {code: ___, n_particles: ___, date: ___}
```
Plus reference datasets: Tier 2 column-collapse profiles/scalings, Tier 3 $F(z,v)$ tables — versioned alongside.

---

## 7. Execution order & rough cost

1. **Pilot** (1 shear run + 1 small collapse) — fix $I(P,\dot\gamma)$ mapping, timestep, steady-state time. *Low cost.*
2. **Tier 1** (~25 shear runs + 1 compression + robustness checks) — the critical path; gives the constitutive law. *Cheap each; embarrassingly parallel.*
3. **Tier 2** (5 collapses) — moderate ($10^4$–$10^5$ grains).
4. **Tier 3** (2 shapes × {quasi-static + ~6 speeds} × 2 packings ≈ 28 runs) — most expensive (large beds), but only after Tiers 1–2 validate.

---

## 8. Open choices to confirm

- **DEM code** (LIGGGHTS vs MercuryDPM vs Yade vs Chrono) — sets input-deck specifics. *Recommend LIGGGHTS or MercuryDPM.*
- **2D-thin vs full 3D** for Tier 1 — full 3D is more faithful; 2D-thin is far cheaper and matches much of the literature. *Recommend a thin 3D slab (periodic in vorticity) as the compromise.*
- **Grain size $d$** — 0.5 mm assumed; confirm against the specific glass-bead lot you want to mirror.

*Drafted 2026-06-16. Companion: `docs/physics-design.md` (consumes these outputs), `docs/literature-review.md`.*
