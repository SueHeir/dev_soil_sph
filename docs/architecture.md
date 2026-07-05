# SPH Solver Architecture — riding the GRASS → SOIL stack

How the SPH solver from `docs/physics-design.md` maps onto Liz's `grass`/`soil`/`dirt` Rust stack. This doc is the bridge between the *physics* (what each particle computes) and the *code* (how it's organized, scheduled, and communicated). It is written from a full read of `grass`, `soil`, `dirt`, and `dev_soil_peri`; file references are concrete so they can be opened directly.

> **Headline:** the bottom half of the SPH primer's architecture (neighbor search, domain decomposition, halo exchange, migration, integration, IO, MPI) **already exists in SOIL**. Our SPH tier is a thin physics layer — a *sibling of DIRT and dev_soil_peri* — that adds kernels, density, and the µ(I) constitutive update. The single real divergence from DEM is that **SPH needs two neighbor passes per step with a halo exchange between them** (§4).

---

## 1. Where SPH sits

```
GRASS   framework: App / Plugin / Scheduler / IO / MPI / coupling      (no particles)
  └─ SOIL   substrate: Atom, AtomData, domain decomp, halo comm, neighbor lists   (no physics)
       ├─ DIRT            DEM physics  (pairwise contact)        ← template for structure
       ├─ dev_soil_peri   peridynamics (bond forces, fracture)   ← template for "new method on SOIL"
       └─ SPH             ← us: continuum granular (kernel + µ(I))
```

`dev_soil_peri` is **peridynamics, not fluids** — but it is the cleanest precedent for *adding a fresh non-DEM particle method to SOIL*, so we mirror its crate shape. We are a sibling tier, not built on DIRT or dev_soil_peri.

---

## 2. Crate layout (mirror `dirt` / `dev_soil_peri`)

| Our crate | Role | Mirrors |
|---|---|---|
| `sph_atom` | `SphAtom` per-particle data column (ρ, deviatoric stress, h, accumulators) + `SphMaterialTable` + particle insertion | `dirt_atom` / `peri_atom` |
| `sph_kernel` | Wendland C2 kernel `W`, `∇W`, support radius helpers (pure functions) | (new; trivial, no DEM analog) |
| `sph_constitutive` | the Dunatunga–Kamrin stress update (`update_stress`), `MaterialParams`, EOS, μ(I) — **pure, substrate-free** so gate #1 is unit-testable in isolation | (new; the §3.3 math) |
| `sph_physics` | the per-step SOIL systems: density/velocity-gradient pass, calls `sph_constitutive::update_stress` per particle, momentum-force pass | `dirt_granular` / `peri_bond` |
| `sph_core` | umbrella: re-export, `SphDefaultPlugins` PluginGroup, `prelude` | `dirt_core` / `peri_core` |

**Reused as-is** (method-agnostic, no fork needed): `grass_app`, `grass_scheduler`, `grass_io`, `soil_core`, `soil_verlet`, `soil_print`, `soil_derive`, `dirt_fixes` (gravity / add-force / viscous damping), and `dirt_wall` (walls + the footpad rig — §6). `CorePlugins` (from `dirt_core`, DEM-free) gives us App/IO/domain/neighbor/run/print.

Dependency shape (per crate `Cargo.toml`, copying `dirt`):
```toml
soil_core    = { git = "...SueHeir/soil", branch = "main", default-features = false }
soil_verlet  = { git = "...SueHeir/soil", branch = "main" }
grass_app    = { git = "...SueHeir/grass", branch = "main" }
grass_scheduler = { git = "...SueHeir/grass", branch = "main" }
sph_atom     = { path = "../sph_atom" }
```

> **Naming:** the tier is presented as **dev_soil_sph** (the `dev_` prefix marks a non-DEM method not personally validated by the author). The crates keep the friendly `sph_*` prefix (`sph_core` / `sph_atom` / `sph_kernel` / `sph_physics`) and the `SphAtom` data type — the same convention as `dev_field_efvm`, whose crates stay `cfd_*`.

---

## 3. The `SphAtom` data column

Per-particle state, declared as a SOIL `AtomData` column. The attribute on each field is a *communication contract*, not decoration:
- `#[forward]` — owner→ghost each step, **overwrite** on unpack. For values a neighbor must *read* (h, and — see §4 — ρ, stress).
- `#[reverse] #[zero]` — ghost→owner, **summed** (`+=`); zeroed each step. For values *accumulated* over neighbors (only needed if we use a half neighbor list — see §4).
- *(no attr)* — migrates with the atom, never ghosted. Persistent per-particle state.

```rust
use soil_derive::AtomData;

#[derive(AtomData)]
pub struct SphAtom {
    /// Smoothing length. Neighbors need it to evaluate the kernel → replicated to ghosts.
    #[forward] pub h: Vec<f64>,

    /// Density ρ — PERSISTENT state (continuity-integrated), AND read by neighbors in the
    /// momentum pass → forward to ghosts. (Mid-step refresh: see §4.)
    #[forward] pub density: Vec<f64>,

    /// Pressure p(ρ). Read by neighbors in the momentum pass → forward to ghosts.
    #[forward] pub pressure: Vec<f64>,

    /// Deviatoric stress (symmetric 3x3 → 6 unique). PERSISTENT state (hypoelastic, evolved
    /// by the return map) AND read by neighbors → forward. Store as [s_xx,s_yy,s_zz,s_xy,s_xz,s_yz].
    #[forward] pub dev_stress: Vec<[f64; 6]>,

    /// Velocity gradient L accumulated in pass 1 (9 comps). Reset each step.
    #[zero] pub velgrad: Vec<[f64; 9]>,

    /// dρ/dt from the continuity sum in pass 1. Reset each step.
    #[zero] pub drho_dt: Vec<f64>,

    /// Rest mass — neighbors need it for kernel sums (ρ_i = Σ m_j W, volume m_j/ρ_j), and
    /// the base `Atom.mass` is NOT forward-comm'd (FORWARD_PACK_SIZE=6 = pos+vel only), so
    /// mass must be a #[forward] SphAtom column to be visible on ghosts.
    #[forward] pub particle_mass: Vec<f64>,
}
```
Register once in the plugin's `build()`:
```rust
register_atom_data!(app, SphAtom::new());   // requires AtomPlugin already added
```
The substrate then carries these columns through every migration, halo exchange, spatial-sort permutation, and restart **with zero comm code from us**.

`density` and `dev_stress` are *both* persistent state *and* `#[forward]` — they are integrated on owners and must be visible to ghost-neighbors. The `#[zero]` accumulators (`velgrad`, `drho_dt`) hold this step's neighbor sums. (Whether they also need `#[reverse]` depends on the list choice in §4.)

---

## 4. The per-step pipeline — where SPH diverges from DEM

DEM is **one** neighbor pass (pairwise contact → force). SPH is **two** passes — density/velocity-gradient, then momentum — with a **per-particle constitutive update and a halo exchange of stress in between**. This is the crux of the whole design.

### 4.1 Neighbor-list choice: gather (full list) vs scatter (half list)

| | Half list (`newton=true`) — DEM default | **Full list (`newton=false`) — recommended for SPH v0** |
|---|---|---|
| Each pair | once; write both `i` and ghost `j`, then reverse-comm | twice; each owner `i` gathers all neighbors, writes only `i` |
| Reverse comm | needed (auto at `PostForce`) | **none** |
| Timing risk | reverse-comm only fires once/step (at `PostForce`) — awkward for a *pre*-force density sum | each owner's sums are self-contained at any phase |
| Cost | ~½ the pair evaluations | ~2× pair evaluations |

**Recommendation: full list (`[neighbor] newton = false`) for v0.** SPH's two passes each need a *complete* per-owner sum *before* the next phase; a half list's reverse comm only folds ghost contributions at `PostForce`, which doesn't align with a density sum that must finish at `PreForce`. The gather formulation makes each pass self-contained (owner reads ghosts, writes only itself), eliminating reverse-comm timing entirely. Revisit the half-list optimization later if pair cost dominates.

With a full list the `#[reverse]` attributes on `velgrad`/`drho_dt` are unnecessary — `#[zero]` alone suffices.

### 4.2 The pipeline (mapped to `ParticleSimScheduleSet` phases)

SOIL's per-step phase order: `Setup → PreInitialIntegration → InitialIntegration → PostInitialIntegration → PreExchange → Exchange → PreNeighbor → Neighbor → PreForce → Force → PostForce → PreFinalIntegration → FinalIntegration → PostFinalIntegration`.

```
InitialIntegration   [soil_verlet]   half-kick v, drift x                     (positions move)
PreNeighbor          [soil]          borders / forward-comm ghosts (pos,vel,h,ρ,p,s from prev step)
Neighbor             [soil]          (re)build full neighbor list
PreForce  ┌ pass 1   [sph_physics]   GATHER: each owner i sums over neighbors:
          │                            velgrad L_i  +=  (m_j/ρ_j)(v_j−v_i) ⊗ ∇W_ij
          │                            drho_dt_i    +=  ρ_i (v_i−v_j)·∇W_ij      (continuity)
          ├ integrate ρ [sph_physics] ρ_i  +=  drho_dt_i · dt                    (continuity state)
          ├ constitutive[sph_physics] p_i = EOS(ρ_i)  [separation: 0 if ρ<ρ_c]   ← physics-design §3.2
          │                            D,W from L_i; hypoelastic trial; µ(I) return map → dev_stress_i   ← §3.3
          └ HALO EXCHANGE ★ [sph]      forward-comm density, pressure, dev_stress  owner→ghost
Force     ┌ pass 2   [sph_physics]   GATHER: each owner i sums neighbor stress divergence:
          │                            force_i += Σ_j m_j(σ_i/ρ_i² + σ_j/ρ_j²)·∇W_ij   ← physics-design §4
          └ body force [dirt_fixes]   + gravity ; + artificial viscosity / stress if enabled
FinalIntegration     [soil_verlet]   final half-kick v
(PreExchange/print)  [soil_print]    thermo + dumps
```

★ **The one genuinely SPH-specific piece of plumbing.** The stress computed at `PreForce` on owners must reach ghost-neighbors before the `Force` pass reads it. SOIL's automatic forward-comm of `#[forward]` columns runs at `PreNeighbor` — i.e. it carries *previous-step* stress, not the stress we just computed. So we need an **explicit intra-step forward-comm of `density`/`pressure`/`dev_stress` between `PreForce` and `Force`**.

**RESOLVED (verified against `soil_core/src/comm.rs`, 2026-06-16).** SOIL already exposes exactly the hook we need: `forward_comm_borders` is a **`pub fn` system** (`comm.rs:553`) that forward-comms ghost pos/vel **plus every registry `#[forward]` column** (`registry.pack_forward_all`, `comm.rs:392`). The `CommunicationPlugin` registers it at `PreNeighbor` gated on `CommState::CommunicateOnly`, but nothing stops a tier from **registering a second, ungated instance at `PreForce`**, after the constitutive system:

```rust
use soil_core::comm::forward_comm_borders;
app.add_update_system(
    forward_comm_borders.after("sph_const"),   // push freshly-computed ρ, p, σ to ghosts
    ParticleSimScheduleSet::PreForce,
);
```

It reuses the swap topology recorded at the last rebuild (`CommTopology.swap_data`), so it is correct on both FullRebuild and CommunicateOnly steps. **No substrate change is needed.** The only cost is mild redundancy — it re-sends pos/vel/h and *all* `#[forward]` columns, not just stress (can't scope to a subset) — negligible for v0. Everything else is a direct copy of the DIRT idiom.

### 4.3 The gather force loop (idiom, copied from `dirt_granular`)

```rust
fn sph_momentum(atoms: ResMut<Atom>, neighbor: Res<Neighbor>, registry: Res<AtomDataRegistry>) {
    let sph = registry.expect::<SphAtom>("SphAtom");
    let mut atoms = atoms;
    let nlocal = atoms.nlocal as usize;
    for (i, j) in neighbor.pairs(nlocal) {        // full list: every neighbor of i, ghost or local
        let gradw = grad_w(&atoms, &sph, i, j);   // ∇W_ij from sph_kernel
        // σ/ρ² for both endpoints; pressure & dev_stress on ghost j are valid (forward-comm'd ★)
        let term = stress_term(&sph, i, j);        // m_j (σ_i/ρ_i² + σ_j/ρ_j²)
        for d in 0..3 { atoms.force[i][d] += dot(term, gradw, d); }
        // full list → write ONLY i (owner gathers; no ghost write, no reverse comm)
    }
}
```
Register the systems into phases, ordering within `PreForce` by label:
```rust
impl Plugin for SphPhysicsPlugin {
    fn build(&self, app: &mut App) {
        app.add_update_system(sph_density_velgrad.label("sph_density"), ParticleSimScheduleSet::PreForce)
           .add_update_system(sph_constitutive.label("sph_const").after("sph_density"), ParticleSimScheduleSet::PreForce)
           .add_update_system(sph_forward_stress.after("sph_const"), ParticleSimScheduleSet::PreForce) // ★
           .add_update_system(sph_momentum.label("sph_force"), ParticleSimScheduleSet::Force);
    }
}
```

---

## 5. The constitutive update lives here

`physics-design.md` §3.3 (the cancellation-safe Dunatunga–Kamrin return map) is the body of `sph_constitutive` — a **per-particle, no-neighbor** system at `PreForce`:
- in: `velgrad` L_i (from pass 1), `density` ρ_i (just integrated), `dev_stress` s_i^n (persistent state), `dt`, material params (a `Res<SphMaterialTable>`).
- out: `pressure` p_i, updated `dev_stress` s_i^{n+1}; assemble σ_i = −p_i I + s_i for pass 2.
- It is the pure function `(s_n, D, W, ρ, Δt, params) → (s_{n+1}, p)` behind a clean boundary — swapping in cohesion or a different rheology touches only this system. Its **unit test** (vs RK4, no SPH) is the first acceptance gate and needs none of the substrate.

Density is carried as integrated state (continuity, for the free surface) — note this means `sph_density_velgrad` produces `drho_dt` and a small step integrates ρ; `soil_verlet` integrates only x and v, so **ρ integration is our system**, not the substrate's. (Alternative: summation density `ρ_i = Σ m_j W_ij` avoids ρ-state but underestimates at free surfaces — deferred.)

---

## 6. Boundaries and the footpad — largely already built

`dirt_wall` gives plane/cylinder/sphere/region walls with a `force_accumulator` and prescribed `WallMotion`. **The footpad touchdown rig already exists as a DEM example: `bench_plate_sinkage`** — a downward-facing plane wall with a finite footprint (`bound_x_low/high`…), driven at `ConstantVelocity`, with the net vertical load read straight off `WallPlane::force_accumulator`:

```toml
[[wall]]
name = "footpad"
normal_z = -1.0 ; point_z = 0.17
velocity = [0.0, 0.0, -0.08]      # prescribed push-in (or use Servo for force-control)
bound_x_low = -0.02 ; bound_x_high = 0.02
```
```rust
let plate = walls.planes.iter().find(|w| w.name == "footpad").unwrap();
let load  = plate.force_accumulator.abs();   // N on the footpad, this step
let depth = z_contact - plate.point_z;       // point_z advances by v·dt in wall_move
```
- **Prescribed-motion footpad** (impose v(t)): `WallMotion::ConstantVelocity` — verbatim from `bench_plate_sinkage`.
- **Force-controlled / floating footpad** (mass m, settles under load + gravity): `WallMotion::Servo` is the template P-controller; or add a `Dynamic { mass }` motion variant that integrates `v += (force_accumulator + applied)/m · dt`. The accumulator plumbing (zeroed `PreForce`, summed `Force`, consumed `PreInitialIntegration`) is already in place. **This is the hook for the leg ring-down deferred in `physics-design.md` §10.**

**What we replace:** `dirt_wall`'s Hertz contact *force kernel* is DEM-specific. For SPH we swap it for a boundary scheme (dummy/ghost boundary particles à la Bui/Morris, or a penalty/repulsive wall) — but **reuse the `Walls` resource, geometry types, active-flags, motion, and `force_accumulator` wholesale**. Runtime control like `Walls::deactivate_by_name("gate")` gives us the column-collapse gate release for free.

---

## 7. Assembling a run — config + `main()` + examples

Thin driver; geometry/materials/stages live in `config.toml` (the DIRT convention):
```rust
use sph_core::prelude::*;
fn main() {
    let mut app = App::new();
    app.add_plugins(CorePlugins)          // GRASS+SOIL: app, IO, domain, neighbor, run, print
       .add_plugins(SphDefaultPlugins)    // our umbrella: SphAtom + insert + verlet + ρ-integrate + physics
       .add_plugins(GravityPlugin)        // dirt_fixes
       .add_plugins(WallPlugin);          // dirt_wall (footpad, container, gate)
    app.start();                          // builds [domain], inserts [[particles.insert]], runs [[run]] stages
    dump_results(&app);
}
```
`SphDefaultPlugins` is a `PluginGroup` bundling, in phase order: `SphAtomPlugin` (column + materials) → `SphInsertPlugin` → `VelocityVerletPlugin::new()` → `SphDensityIntegratePlugin` → `SphPhysicsPlugin` (the §4 systems).

**Validation gate #3 — column collapse** is a near-clone of `bench_column_collapse`: settle stage, then a `collapse` stage whose `on_enter` opens the gate wall via `Walls::deactivate_by_name`; the multi-stage state machine uses `#[derive(StageEnum)]` + `StatesPlugin` + `StageAdvancePlugin`, gated by `in_stage("collapse")`. Post-run, read final `(pos, ρ, σ)` off the `Atom` + `AtomDataRegistry` resources and write the deposit CSV for the Lagrée error-map fit.

Config sections we consume: `[comm]`, `[domain]`, `[neighbor]` (set `bin_size ≥ κh`, `newton = false`), `[gravity]`, `[sph]` (model + `[[sph.materials]]` with µ_s, µ_2, I_0, ρ_c, K from the DEM fit), `[[particles.insert]]`, `[[wall]]`, `[output]`, `[vtp]`, `[[run]]` stages.

Dump columns via `DumpRegistry::register_scalar("density", …)` / `register_vector("velocity", …)` in `build()`.

---

## 8. Diagnostics

`dirt_measure_plane` is the template for profiles: a config-driven array of plane definitions, per-tag signed-distance tracking, windowed averaging into `Thermo`. For **stress/velocity profiles** (Bagnold check, column-collapse internal fields) we follow the same plugin shape but bin particles into slabs and accumulate per-bin momentum/stress instead of counting crossings. Simple bench-style recorders (stream `pos`/`vel`/`σ`/`force_accumulator` to CSV from `PostFinalIntegration`) cover most validation needs.

---

## 9. Open decisions / risks (resolve before/at first code)

1. ~~**Mid-step halo exchange of stress (★, §4.2)**~~ — **RESOLVED.** Register a second ungated `soil_core::comm::forward_comm_borders` at `PreForce` after the constitutive system; it forward-comms all `#[forward]` columns (incl. ρ, p, σ) with no substrate change. The two-pass design is de-risked.
2. **Full vs half neighbor list** — recommending full (`newton=false`) gather for v0 simplicity/correctness; half-list is a later perf optimization.
3. **Continuity vs summation density** — recommending continuity (free surface), which means a tier-side ρ-integration system (soil_verlet won't do it).
4. **Boundary scheme** — dummy boundary particles vs penalty wall for the SPH–wall/footpad contact; affects column-collapse basal friction (DIRT flags runout sensitivity to exactly this).
5. ~~**Crate naming**~~ — DECIDED: tier presented as **dev_soil_sph**; crates keep the `sph_*` prefix.

---

## 10. Suggested build order (first commits)

1. **`sph_kernel`** — Wendland C2 `W`, `∇W` + unit tests (partition of unity, gradient sum). Pure, no substrate. *Trivial, high-confidence start.*
2. **Constitutive unit test** — implement `physics-design.md` §3.3 as a standalone function + the RK4 comparison test (gate #1). No substrate needed.
3. **`sph_atom`** — `SphAtom` column + `register_atom_data!` + materials + insertion (mirror `dirt_atom`).
4. **`sph_physics` pass 1 + constitutive + pass 2** + the ★ halo exchange; wire `SphDefaultPlugins`.
5. **Hydrostatic-column test** — a settled bed under gravity holds still with correct pressure (shakes out EOS, BCs, the halo exchange).
6. **Column-collapse example** — gate #3, against the DEM data from `docs/dem-campaign.md`.

*Drafted 2026-06-16 from full reads of grass/soil/dirt/dev_soil_peri. Companions: `docs/physics-design.md`, `docs/dem-campaign.md`, `docs/sph-primer.md`, `docs/literature-review.md`.*
