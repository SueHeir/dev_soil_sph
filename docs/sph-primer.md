# SPH Primer — for the dry-granular landing-stability codebase

A working mental model of Smoothed Particle Hydrodynamics, pitched for someone fluent in
continuum mechanics and constitutive modeling but new to the *particle-discretization*
machinery. Focus: the SPH-specific moves, where the method bites, and where granular
material makes it worse.

---

## 1. The one idea everything hangs on

SPH represents a continuous field and its spatial derivatives using a scattered cloud of
points — **no mesh, no connectivity**. Write any field $A(\mathbf{x})$ as a convolution
against a smoothing kernel $W$, then approximate that integral as a sum over neighbor
particles:

$$
A(\mathbf{x}) \approx \int A(\mathbf{x}')\, W(\mathbf{x}-\mathbf{x}', h)\, d\mathbf{x}'
\approx \sum_j \frac{m_j}{\rho_j}\, A_j\, W(\mathbf{x}-\mathbf{x}_j, h)
$$

That second step **is** the method: the integral $d\mathbf{x}'$ becomes $\sum_j (m_j/\rho_j)$,
i.e. each particle carries a volume

$$
V_j = \frac{m_j}{\rho_j}.
$$

Particles are quadrature points that move with the material (Lagrangian). $h$ is the
**smoothing length** — the kernel's support-radius scale. As $h \to 0$ and neighbor count
$\to \infty$, $W \to \delta$ and the approximation becomes exact.

---

## 2. Why this fits *your* problem

- **Lagrangian + meshfree** → arbitrarily large deformation, free surfaces, and material
  separation with no remeshing. Exactly why granular collapse, cratering, and a footpad
  opening a wake behind it are natural in SPH and painful in FEM.
  (MPM — Kamrin's tool — is the mesh-based cousin: particles carry state, but momentum is
  solved on a background grid.)
- **Advection is free.** Particles *are* the material, so there is no convective term
  $\mathbf{v}\cdot\nabla(\cdot)$ to discretize. The material derivative

  $$
  \frac{D A}{Dt} = \frac{\partial A}{\partial t} + \mathbf{v}\cdot\nabla A
  $$

  reduces to just the rate of change following a particle, $\dot{A}_i$. A real
  simplification for a velocity-strengthening granular rheology.

---

## 3. The derivative trick — the actual "magic"

You never differentiate your field data. You differentiate the **kernel** — a smooth
analytic function you chose. Integration by parts moves the derivative off $A$ and onto $W$:

$$
\nabla A(\mathbf{x}_i) \approx \sum_j \frac{m_j}{\rho_j}\, A_j\, \nabla_i W_{ij},
\qquad W_{ij} \equiv W(\mathbf{x}_i - \mathbf{x}_j, h).
$$

A gradient is a weighted sum of neighbor values times the known kernel gradient.

**But this naïve form is bad** — it does not even return zero for the gradient of a constant,
because the particle distribution is never perfectly regular. In practice always use a
**symmetrized** form. The two standard ones:

$$
\nabla A_i \approx \rho_i \sum_j m_j
\left(\frac{A_i}{\rho_i^2} + \frac{A_j}{\rho_j^2}\right)\nabla_i W_{ij}
\qquad\text{(symmetric)}
$$

$$
\nabla A_i \approx \frac{1}{\rho_i}\sum_j m_j\,(A_j - A_i)\,\nabla_i W_{ij}
\qquad\text{(difference)}
$$

This choice is not cosmetic:

- The **symmetric** form gives pairwise antisymmetric forces, $\mathbf{f}_{ij} = -\mathbf{f}_{ji}$,
  which makes discrete linear/angular momentum *exactly* conserved. Use it for momentum.
- The **difference** form vanishes for constant $A$ (better consistency) but is not
  conservative. Use it for things like velocity divergence in the continuity equation.

---

## 4. The kernel $W$

Requirements: normalizes to 1, compact support (radius $\kappa h$, typically $2h$ or $3h$),
positive, symmetric, $\to \delta$ as $h\to 0$, and smooth enough that $\nabla W$
(and, for some viscosity forms, $\nabla^2 W$) is well-behaved.

| Kernel | Support | Notes |
|---|---|---|
| **Cubic spline** | $2h$ | Classic, cheap. Its $\nabla W$ has an inflection that triggers the **pairing/tensile instability** (particles clump in pairs). |
| **Wendland C2** | $2h$ | Modern default. Non-negative Fourier transform → suppresses pairing; tolerates more neighbors. **Start here.** |
| **Quintic spline** | $3h$ | Smoother, more accurate, more expensive. |

For granular/solid SPH the kernel's *stability* properties matter more than raw accuracy,
because you run near the failure modes in §8.

---

## 5. Discretizing the governing equations

### Density — a real design choice

**(a) Summation density**

$$
\rho_i = \sum_j m_j\, W_{ij}
$$

Exactly conserves mass, but **underestimates density at free surfaces** (kernel deficiency —
missing neighbors on the open side). Fatal for a granular free surface.

**(b) Continuity / density-rate** — integrates $\dot\rho = -\rho\,\nabla\cdot\mathbf{v}$:

$$
\frac{d\rho_i}{dt} = \sum_j m_j\,(\mathbf{v}_i - \mathbf{v}_j)\cdot\nabla_i W_{ij}
$$

No surface deficit. **The standard choice for free-surface granular flow — you will almost
certainly use this one.**

### Momentum — where your constitutive model plugs in

With the full Cauchy stress $\boldsymbol{\sigma}$:

$$
\frac{d\mathbf{v}_i}{dt} = \sum_j m_j
\left(\frac{\boldsymbol{\sigma}_i}{\rho_i^2} + \frac{\boldsymbol{\sigma}_j}{\rho_j^2}\right)
\cdot \nabla_i W_{ij} \;+\; \mathbf{g}
$$

For a fluid, $\boldsymbol{\sigma} = -p\,\mathbf{I} + \boldsymbol{\tau}$. For granular,
$\boldsymbol{\sigma}$ is whatever your µ(I) / Drucker–Prager / hypoplastic model returns.

> **The SPH momentum equation does not care where the stress comes from.** That is exactly
> why the Wang-vs-Kamrin choice is a *constitutive* decision layered cleanly on top of this
> same momentum sum — and why it should sit behind a clean interface in the code (§11).

---

## 6. Pressure and the weakly-compressible trick (WCSPH)

Real grains are nearly incompressible, but true incompressibility (ISPH) means a global
pressure-Poisson solve every step. WCSPH avoids this with an **artificial equation of state** —
treat the material as slightly compressible and get pressure algebraically from density
(Tait form, $\gamma \approx 7$):

$$
p_i = \frac{\rho_0\, c_s^2}{\gamma}
\left[\left(\frac{\rho_i}{\rho_0}\right)^{\gamma} - 1\right]
$$

Pick a **numerical** sound speed $c_s$ large enough that density fluctuations stay under
~1%. Rule of thumb:

$$
c_s \gtrsim 10\, v_{\max}.
$$

The price is the **CFL timestep**, set by $c_s$:

$$
\Delta t \lesssim 0.25\,\frac{h}{c_s}.
$$

A stiffer EOS → smaller steps. For you, $v_{\max}$ includes the footpad impact velocity, so
impact speed directly sets your timestep budget. Most dry-granular SPH papers (Bui, Minatti,
Szewc) are WCSPH — **start there over ISPH.**

---

## 7. Neighbor search and smoothing length

Each step, find every particle's neighbors within $\kappa h$. Use a **cell-linked list /
hash grid** (bin particles into cells of size $\kappa h$, check only adjacent cells) → $O(N)$.
This is usually the performance-critical kernel of the whole code.

- Typical neighbor count: ~30–50 in 2D, ~50–80 in 3D. Too few → noise/instability; too many
  → cost and over-smoothing.
- Smoothing length can be constant or adaptive, $h \propto (m/\rho)^{1/d}$. Variable $h$
  helps under the strong compaction of an impact but complicates the symmetry of
  $\nabla W_{ij}$ (symmetrize via $\bar h = \tfrac12(h_i+h_j)$ or average the two kernels).

---

## 8. Where it bites — the pathologies you *will* hit

This is what separates "ran the demo" from "trust the footpad force."

- **Tensile instability / particle pairing.** When the stress state is tensile (or
  $\nabla W$ enters its wrong-sign region), particles attract and clump into unphysical
  clusters/voids. Granular materials sit near tensile states at the free surface and behind
  an intruder, so this is the dominant headache.
  *Fixes:* Wendland kernels; **artificial stress** (Monaghan–Gray — what Bui 2008 uses);
  **stress-point / particle-shifting** schemes (the 2021 *Computers & Geotechnics* and
  Wang 2023 papers).

- **Zeroth/first-order inconsistency.** Disordered particles mean
  $\sum_j V_j \nabla W_{ij} \ne 0$, so even a constant field shows spurious gradients.
  *Cures:* kernel-gradient correction (a per-particle renormalization matrix), or MLS/CSPH
  corrections. Improves accuracy but can hurt conservation — a tradeoff.

- **Stress noise / "checkerboarding."** High-frequency oscillation in the stress field from
  SPH's collocated nature; granular constitutive models amplify it. Hence **stress
  regularization** (MLS filtering, particle shifting) is load-bearing infrastructure, not
  polish.

- **Artificial viscosity.** A pairwise dissipative term $\Pi_{ij}$ added to the momentum
  equation to damp post-shock oscillation and stabilize. Tune its coefficients: too much
  smears penetration depth, too little lets noise blow up. For granular you often let
  physical frictional dissipation do this job, but usually still need a little for stability.

---

## 9. Boundary conditions — the genuinely hard part

No mesh boundary, so walls and the rigid footpad are awkward. Common approaches:

- **Dummy / ghost particles** — fill the wall with frozen SPH particles so kernels stay
  complete near boundaries; assign mirrored pressure/velocity. Most common.
- **Boundary force particles** (Monaghan repulsive forces) — simpler, less accurate.
- **Mirror / semi-analytic** boundaries.

For you, the **footpad is a moving rigid boundary with frictional contact.** That contact
law (and the Wang 2023 "generalized frictional boundary" paper) is where the touchdown force
is actually computed — it deserves first-class attention, not a default no-slip wall.

---

## 10. Time integration

Explicit, low-storage schemes: **kick–drift–kick / leapfrog**, or predictor–corrector
(Verlet). Symplectic-ish, cheap, good energy behavior. The timestep is the minimum of three
constraints:

$$
\Delta t \le \min\!\left(
\underbrace{0.25\,\frac{h}{c_s}}_{\text{CFL / sound}},\;\;
\underbrace{0.25\,\sqrt{\frac{h}{|\mathbf{a}|}}}_{\text{body force}},\;\;
\underbrace{0.125\,\frac{h^2}{\nu}}_{\text{viscous}}
\right)
$$

Everything explicit → trivially parallel, GPU-friendly. That is why these codes scale.

---

## 11. The architecture takeaway

The layering this implies for the codebase:

1. **Neighbor search** (hash grid) — performance core, physics-independent.
2. **Kernel module** ($W$, $\nabla W$, corrections) — swappable.
3. **Field operators** (SPH gradient/divergence forms) — the discrete calculus.
4. **Constitutive model** (returns $\boldsymbol{\sigma}$ from deformation/history) — *this*
   is where Kamrin-vs-Wang lives; it plugs into one line of the momentum sum.
5. **Boundary / contact** (footpad) — where the answer comes from.
6. **Time integrator** — orchestrates the above.

Keep (4) behind a clean interface and you can prototype both the **density-switched
single-stress** model (Dunatunga–Kamrin) and the **two-stress additive** model (Wang/Wu)
against the same column-collapse validation without touching the solver core.

---

## Symbol reference

| Symbol | Meaning |
|---|---|
| $W_{ij}$, $\nabla_i W_{ij}$ | kernel and its gradient between particles $i,j$ |
| $h$ | smoothing length; $\kappa h$ = support radius |
| $m_j,\ \rho_j,\ V_j=m_j/\rho_j$ | particle mass, density, volume |
| $\boldsymbol{\sigma},\ p,\ \boldsymbol{\tau}$ | Cauchy stress, pressure, deviatoric stress |
| $c_s$ | numerical sound speed (WCSPH) |
| $\rho_0,\ \gamma$ | reference density and Tait exponent |
| $\Pi_{ij}$ | artificial viscosity term |
| $\mathbf{g},\ \mathbf{a}$ | body force, acceleration |
