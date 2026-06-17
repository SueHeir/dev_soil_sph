# Executing the DEM Calibration Campaign on DIRT — Bridging Doc

**Purpose.** `docs/dem-campaign.md` is the *spec* (3 tiers, glass beads, target outputs). This doc maps that spec onto code that exists today — the DIRT DEM crates (`/Users/suehr/Documents/GitHub/dirt/`), the SOIL base (`/Users/suehr/Documents/GitHub/soil/`), and Liz's existing LAMMPS shear cell (`/Users/suehr/Documents/GitHub/lammps_shear_cell/`). It says, per tier, *which existing example/capability to run, what config to write, and what to change*. It also identifies the one real gap: a pressure-controlled homogeneous-shear rheometer for Tier 1.

All file paths are absolute. Read alongside `docs/dem-campaign.md`.

---

## 0. How a DIRT run is configured and launched (ground rules)

A DIRT example is a thin Rust binary (`main.rs`) that assembles plugins and a recorder, plus a `config.toml` that holds all the physics. Run with:

```bash
cd /Users/suehr/Documents/GitHub/dirt
cargo run --release --example <name> --no-default-features -- examples/<name>/config.toml
```

(confirmed in `/Users/suehr/Documents/GitHub/dirt/README.md:54` and every example's header comment.)

The TOML sections used across the granular examples:

| Section | Role | Source crate |
|---|---|---|
| `[comm]` | MPI process grid (`processors_x/y/z`) | grass/soil core |
| `[domain]` | box bounds `x_low..z_high`, `boundary_x/y/z = "fixed"|"periodic"` | soil_core |
| `[neighbor]` | `skin_fraction`, `bin_size`, `every` | soil_verlet |
| `[gravity]` | `gx/gy/gz` | soil_core |
| `[dem]` | `contact_model = "hertz"` (or `"hooke"`) | dirt_granular |
| `[[dem.materials]]` | per-material contact params (see §below) | `dirt_atom/src/lib.rs:77` (`MaterialConfig`) |
| `[[particles.insert]]` | fill a region with spheres | `dirt_atom/src/insert.rs:73` (`InsertConfig`) |
| `[[wall]]` | plane/cylinder/sphere/region walls + motion/servo | `dirt_wall/src/lib.rs` |
| `[contact_analysis]` | coordination / rattlers / fabric tensor | dirt_contact_analysis |
| `[[measure_plane]]` | particle/mass flux across a plane | dirt_measure_plane |
| `[deform]` | continuous box deformation (x/y/z axis lengths only) | `soil_deform/src/lib.rs` |
| `[output]`, `[vtp]`, `[[run]]` | output dir, VTP interval, staged run blocks | soil core |

### Material (contact) keys — `[[dem.materials]]`

From `dirt_atom/src/lib.rs:77` (`MaterialConfig`) and confirmed in every example config. **The serde key is `youngs_mod`, not `youngs_modulus`** (the dirt_granular doc-comment uses `youngs_modulus`, but the live struct and all real configs use `youngs_mod`):

```toml
[dem]
contact_model = "hertz"        # Hertz-Mindlin nonlinear elastic (default); "hooke" = linear spring

[[dem.materials]]
name           = "glass"
youngs_mod     = 7.0e7         # Pa  (E)  — softened, see §timestep
poisson_ratio  = 0.25          # ν
restitution    = 0.7           # e
friction       = 0.5           # μ_p  (Coulomb sliding cap on the Mindlin tangential spring)
rolling_friction = 0.0         # μ_r  (default 0)
```

The contact model is **Hertz normal + Mindlin incremental tangential spring with Coulomb cap `μ|F_n|`** (`dirt_granular/src/lib.rs:10-15`) — exactly the Hertz–Mindlin the spec asks for. Cross-material params mix by geometric mean (`friction_ij`, etc.).

### Particle insertion + polydispersity — `[[particles.insert]]`

From `dirt_atom/src/insert.rs:73` and `dirt_atom/src/radius.rs`. Random rejection-sampled insertion into a region; **polydispersity is a first-class feature** via a uniform radius distribution:

```toml
[[particles.insert]]
material = "glass"
count    = 10000
radius   = { distribution = "uniform", min = 0.000225, max = 0.000275 }  # d = 0.45–0.55 mm
density  = 2500.0
region   = { type = "block", min = [..], max = [..] }
# seed = 12345   # deterministic, reproducible across ranks (insert.rs:127)
```

This directly satisfies the spec's "uniform in radius, ±10–20 % about d = 0.5 mm" (§2.1). `radius` also accepts a plain number for monodisperse, or a `lognormal` distribution.

### Walls + the servo (pressure control) — `[[wall]]`

From `dirt_wall/src/lib.rs`. Wall types: `plane`, `cylinder`, `sphere`, `region`. **Plane walls** support four motion modes (`dirt_wall/src/lib.rs:17-26`): static, constant `velocity`, `oscillate`, and **`servo`** — a proportional controller that moves the wall to hold a target contact force:

```toml
[[wall]]
type = "plane"
point_z = 0.1
normal_z = -1.0
material = "glass"
servo = { target_force = 100.0, max_velocity = 0.1, gain = 0.001 }   # ServoDef, lib.rs:157
```

Each step: `error = target_force − measured_force`, `v = clamp(gain·error, ±max_velocity)` along the normal (`lib.rs:144-164, 810`). The measured force is `WallPlane::force_accumulator` (`lib.rs:339`), the summed scalar normal contact force on that wall this step. **This is the pressure-control mechanism** — set `target_force = P · A_wall`.

**Wall friction:** `dirt_wall` *does* apply Mindlin tangential (sliding) friction on plane walls via per-contact `tangential_springs` (`dirt_wall/src/lib.rs:474-478`), using the material `friction` (μ) through `friction_ij`. The angle-of-repose example relies on exactly this for basal arrest (`bench_angle_of_repose/main.rs:14-18`, `config.toml:9-12`). (Note: `bench_column_collapse/config.toml` carries an older "no particle–wall sliding friction" comment — that is stale relative to the angle-of-repose example, which demonstrates working wall friction. Verify on first run via the repose angle.)

### Reading forces out of a run

The plate-sinkage recorder is the template for any in-run measurement: it pulls `Res<Walls>`, finds the moving plane, and reads `plate.force_accumulator.abs()` each step into a CSV (`bench_plate_sinkage/main.rs:53-88`). The same `force_accumulator` exists on every wall, so a Tier-1 cell can read normal force off the top/bottom walls with a ~20-line recorder modeled on `record_sinkage`.

**Bulk stress is also available** (key for Tier 1): `VirialStress` (`/Users/suehr/Documents/GitHub/soil/crates/soil_core/src/virial.rs`) is the full Love–Weber tensor, auto-registered by the Hertz–Mindlin contact plugin and accumulated (normal + tangential) every step. It is an in-memory `Res<VirialStress>` with **no built-in output** — a recorder must read it and divide by box volume (σ_ij = −VirialStress_ij/V). See §4.

---

## 1. Tier → capability map

| Tier | Spec config | DIRT capability today | Status |
|---|---|---|---|
| **1. Element tests** (pressure-controlled homogeneous simple shear → μ(I), Φ(I), ρ_c, K) | Lees–Edwards or moving rough walls, fixed σ_yy | **No direct example.** `soil_deform` does box stretch only (no shear/tilt). Servo wall gives pressure control. Moving-wall shear is buildable; Lees–Edwards is not. | **GAP — see §2** |
| **2. Column collapse** (aspect sweep → runout/shape/energy) | Free-surface collapse on rough floor, gate removal | **`bench_column_collapse`** — runs as-is; quasi-2D slab, gate wall removed at stage 2, dumps deposit CSV. | **Ready** (§5) |
| **3. Footpad intrusion** (rigid intruder, F(z,v) → K_bear, λ) | Rigid plate/sphere driven into a deep bed, read reaction force | **`bench_plate_sinkage`** — flat plate (clipped plane wall) driven down, reads `force_accumulator` vs sinkage to CSV. Sphere intruder via `type="sphere"` wall. | **Ready** (§5) |

Supporting examples that confirm the building blocks: `bench_angle_of_repose` (Mindlin wall friction, cylinder-lift insertion), `hopper` (basic granular: insert + multiple plane walls + staged run), `bench_sliding_friction`/`bench_hertz_rebound`/`bench_oblique_impact` (contact-law validation — use these to confirm e, μ_p calibration before the campaign).

---

## 2. Tier 1 rheometer — the gap and the recommendation

### Can DIRT do pressure-controlled homogeneous simple shear *today*?

**Not the canonical Lees–Edwards version, no.** The shear-driving primitive is missing:

- `soil_deform` (`/Users/suehr/Documents/GitHub/soil/crates/soil_deform/src/lib.rs`, 925 lines) is DIRT's analog of LAMMPS `fix deform`, but it only deforms the three **diagonal** box lengths (`x`, `y`, `z`) with styles `erate` / `vel` / `final`. There is **no triclinic tilt / xy-shear style and no Lees–Edwards remapping** (grep for `shear|triclinic|tilt|xy|lees` returns nothing). So you cannot impose `fix deform xy erate` — the exact primitive Liz's LAMMPS cell uses (`lammps_shear_cell/shear_cell.lam`: `fix shear all deform 1 xy erate 100.0 remap v`).
- DIRT *does* have the other half — pressure control — via the servo wall (§0).

**Two ways to close it, in increasing effort:**

**(A) Moving-rough-wall shear cell in DIRT (no code change).** Confine the sample between two rough plane walls in `y`; drive the top wall in `+x` at constant `velocity = [V, 0, 0]` (γ̇ = V/L_y), make both walls "rough" by giving them high `friction` (or, more robustly, gluing a layer of frozen grains as the wall — `dirt_fixes` `[[freeze]]` immobilizes a grain group). Periodic in `x` and `z`. For pressure control, replace the top wall with a `servo` wall **and** drive it in x — but the servo only actuates *along the wall normal* (`y`), so you need a y-servo wall that is *also* given a constant x-velocity. The WallDef allows `velocity` and `servo` to coexist as fields; whether the integrator superposes the servo-normal motion with a constant tangential velocity needs a one-line check/confirmation in `wall_move` (`dirt_wall/src/lib.rs:766-832`). If it does, this is a **zero-code-change** Couette rheometer. If not, it's a few-line addition to `wall_move`.
  - *Caveat:* moving-wall shear is not perfectly homogeneous (shear bands / wall layering near the plates), which is exactly why the literature prefers Lees–Edwards. Mitigate with a tall cell (≥ 30 d in y) and measure stress in the central bulk only.

**(B) Add Lees–Edwards / xy-deform to `soil_deform` (small, well-scoped code addition).** Add a 4th style that updates the box `xy` tilt factor and applies the Lees–Edwards image shift in the neighbor wrap, mirroring `fix deform xy erate ... remap v`. This is the "right" homogeneous-shear primitive and would make DIRT self-sufficient for Tier 1. It is a contained change in one crate (the remap and box-tilt plumbing already exist for the diagonal styles). Pressure control then comes from a y-servo wall *or* from a stress-controlled box (harder). **Estimated scope: moderate — one new `DeformStyle` variant + tilt bookkeeping + image-shift in the periodic wrap.**

### Recommendation: **run Tier 1 in LAMMPS, not DIRT (for v0).**

Liz already has a working, validated `fix deform xy erate` shear cell with the exact glass-bead material and the stress-tensor post-processing in place (`/Users/suehr/Documents/GitHub/lammps_shear_cell/`). Reasons:

1. **The primitive exists there and not in DIRT.** LAMMPS does true Lees–Edwards homogeneous shear with `temp/deform` streaming-velocity removal — the GDR MiDi / da Cruz standard the spec (§2.3, §3.1) explicitly calls for.
2. **The stress measurement is already wired** (`compute pressure`, collisional + kinetic split, `fix ave/time` → `stress_tensors.dat`). DIRT has *no* virial/stress-tensor output today (see §4) — only fabric/coordination and wall forces.
3. **It's the critical path and embarrassingly parallel** (spec §7); using the ready tool gets the constitutive law fastest.

**Two changes are needed to Liz's existing deck to make it the spec's rheometer** — it is currently *frictionless and fixed-volume*:

1. **Turn on tangential friction.** Current deck: `tangential mindlin NULL 0.0 0.0` (kt=0, xmu=0) → μ_p = 0. The spec needs **μ_p = 0.5**, the one param that sets μ_s. Change to a Hertz–Mindlin tangential with `xmu = 0.5` and a physical `kt` (LAMMPS `tangential mindlin` derives kt from G*, ν; set the friction coefficient to 0.5). Also set restitution to **0.7** (deck currently 0.9/0.95) and E to the softened **5e7–1e8 Pa** (deck currently 8.7e9 — real glass; fine but stiffer ⇒ smaller dt).
2. **Switch volume-controlled → pressure-controlled.** The deck runs fixed-φ (`generate_runs.py` sweeps φ at fixed box). The spec wants **fixed σ_yy = P** so Φ is an *output* (§3.1). In LAMMPS: keep `fix deform xy erate` for shear, but make the y-dimension stress-controlled with `fix deform y ... 1 controller` is not available; instead use **`fix nph`-style barostat on y only** (`fix press/berendsen` or `fix npt/sphere`-equivalent is not granular-clean) — the standard granular route is a **servo wall in y**: replace the periodic-y box face with two `fix wall/gran` walls and drive the top one with `fix move`/`fix smd` (constant-force) or a `fix indent`/`fix aveforce` PID to hold σ_yy. This mirrors DIRT option (A). *This is the one real piece of LAMMPS plumbing to add.* If pressure control proves fiddly, fall back to **fixed-volume runs across a φ-grid** (which the deck already does) and convert φ→P post-hoc via the measured σ_yy(φ); you recover μ(I) and Φ(I) either way, you just don't get ρ_c as cleanly.

Use the LAMMPS deck for Tier 1; **build option (B) in DIRT later** if/when you want the whole pipeline in one code (e.g., for the lunar/regolith phase). Tiers 2 and 3 run natively in DIRT (§5) — so the campaign is LAMMPS-for-Tier-1, DIRT-for-Tiers-2-and-3.

---

## 3. Concrete Tier-1 input

### 3a. LAMMPS deck (recommended path) — pressure-controlled glass-bead shear

Sketch building on `lammps_shear_cell/generate_runs.py`. One run = one (P, γ̇) point. Spec params: d ∈ [0.45,0.55] mm, ρ_s = 2500, E ≈ 7e7, ν = 0.25, e = 0.7, μ_p = 0.5.

```lammps
# === Tier-1 rheometer: pressure-controlled simple shear, glass beads ===
units           si
atom_style      sphere
boundary        p f p          # periodic x (flow), z (vorticity); fixed y (walls for pressure)
comm_modify     vel yes

variable        d   equal 0.0005          # mean diameter 0.5 mm
variable        P   equal 1000.0          # target normal stress sigma_yy (Pa)  -> sweep {1000,100,10}
variable        gdot equal 30.0           # shear rate (1/s)                    -> sweep per table
variable        Lx  equal 30*${d}
variable        Lz  equal 30*${d}
variable        Ly  equal 40*${d}         # tall in shear-gradient dir to fit a bulk core

region          box block 0 ${Lx} 0 ${Ly} 0 ${Lz}
create_box      2 box                      # type 1 = grains, type 2 = wall grains
# ... create_atoms 1 random N seed box, with N from target initial phi ~0.55 ...
set             type 1 diameter ${d} density 2500.0
# mild polydispersity: assign per-atom diameters uniform in [0.45d_lo, 0.55d_hi]

# Hertz-Mindlin WITH friction (the key change vs Liz's frictionless deck):
#   hertz/material  E  COR  nu  ;  tangential mindlin NULL kt_ratio xmu(=mu_p)
pair_style      granular
pair_coeff      * * hertz/material 7.0e7 0.7 0.25 &
                    tangential mindlin NULL 1.0 0.5 &
                    damping coeff_restitution

fix             integrate all nve/sphere

# --- equilibrate random overlaps (Liz's pattern) ---
fix             damp all viscous 0.1
timestep        1e-7
run             100000
unfix           damp
velocity        all set 0.0 0.0 0.0

# --- pressure control on y: rough servo walls hold sigma_yy = P ---
# Two layers of frozen grains as rough walls; top wall driven to hold force = P * Lx * Lz.
# (Use fix wall/gran + a constant-force / PID move on the top wall, OR
#  a frozen-grain wall group moved by fix aveforce to target P*A.)

# --- shear: drive top wall block in +x at V = gdot * Ly  (Couette) ---
#     OR, if staying fully periodic, fix deform xy erate ${gdot} remap v (Lees-Edwards, fixed volume).

# --- stress (Liz's stress block, verbatim) ---
compute         Tdef   all temp/deform
compute         p_coll all pressure NULL pair
compute         p_total all pressure Tdef
# kin = total - coll for each component (see shear_cell.lam)
fix             stress_out all ave/time 100 10 1000 c_p_coll[*] ... file stress_tensors.dat
timestep        1e-7
run             10000000        # to steady state (strain >> 1)
```

The **(P, γ̇) sweep** (spec §3.2 table) targeting I ∈ [1e-4, 0.5]:

| P (Pa) | γ̇ (s⁻¹) | I = γ̇·d·√(ρ_s/P) |
|---|---|---|
| 1000 | 1 → 300 | ~3e-4 → ~0.1 |
| 100 | 1 → 300 | ~1e-3 → ~0.3 |
| 10 | 3 → 300 | ~1e-2 → ~0.5 |

≥5 points/decade, ~20–25 runs, generated by extending `generate_runs.py` to loop over (P, γ̇) instead of φ. Pilot one (P,γ̇) first to fix the I(P,γ̇) map and steady-state time. Add the **isotropic compression test** (no shear, ramp P, measure ΔΦ → K) as 1 extra deck, and the robustness checks (e ∈ {0.5,0.9}, μ_p ∈ {0.3,0.7}, 10×E) as ~7 extra runs.

### 3b. DIRT Couette alternative (if you want it native — option (A))

```toml
[domain]
x_low = 0; x_high = 0.015     # 30 d, periodic (flow)
y_low = 0; y_high = 0.020     # 40 d, walls
z_low = 0; z_high = 0.015     # 30 d, periodic (vorticity)
boundary_x = "periodic"
boundary_y = "fixed"
boundary_z = "periodic"

[gravity]                      # OFF for the rheometer (homogeneous shear, no body force)
gx = 0.0; gy = 0.0; gz = 0.0

[dem]
contact_model = "hertz"
[[dem.materials]]
name = "glass"
youngs_mod = 7.0e7
poisson_ratio = 0.25
restitution = 0.7
friction = 0.5
rolling_friction = 0.0

[[particles.insert]]
material = "glass"
count = 12000
radius = { distribution = "uniform", min = 0.000225, max = 0.000275 }
density = 2500.0
region = { type = "block", min = [0.0, 0.002, 0.0], max = [0.015, 0.018, 0.015] }

# bottom rough wall (static), top rough wall: servo-normal (hold P) + tangential drive (shear)
[[wall]]
type = "plane"; point_y = 0.002; normal_y = 1.0; material = "glass"; name = "bottom"
[[wall]]
type = "plane"; point_y = 0.018; normal_y = -1.0; material = "glass"; name = "top"
velocity = [0.30, 0.0, 0.0]                                   # V = gdot * Ly ; gdot = 0.30/0.016 ≈ 19 s^-1
servo = { target_force = 0.225, max_velocity = 0.05, gain = 1e-6 }  # P*A = 1000 * (0.015*0.015)

[contact_analysis]              # single table (not array); diagnostics only — stress comes from VirialStress (§4)
coordination = true
fabric_tensor = true

[[run]]
name = "shear"; steps = 2000000; thermo = 20000; dt = 2e-7
```

Then read τ and P off `top.force_accumulator` (tangential vs normal components) with a recorder modeled on `bench_plate_sinkage/main.rs:record_sinkage`. **Prerequisite check:** confirm `wall_move` superposes the servo (normal) velocity with the constant tangential `velocity` — see §2(A). Rough walls: high `friction` on the wall material, or a frozen grain layer via `[[freeze]]` (`dirt_fixes`).

---

## 4. Measurement — extracting τ, P, φ, γ̇, then I, μ_eff, Φ

**γ̇** is imposed (γ̇ = V/L_y for Couette, or the `erate` for Lees–Edwards) — known exactly.

**LAMMPS (Tier-1 recommended path):** stress is **already computed** in Liz's deck. `compute pressure` gives the full 6-component virial stress, split into **collisional** (`p_coll`, pair virial = Love–Weber) and **kinetic** (`p_total − p_coll`, fluctuating-velocity contribution via `temp/deform`). Time-averaged to `stress_tensors.dat` by `fix ave/time 100 10 1000` (cols: `coll_xx..coll_yz kin_xx..kin_yz`). From these:
- **τ = σ_xy = coll_xy + kin_xy**
- **P = σ_yy** (pressure-controlled) **or** P = ½(σ_xx+σ_yy) [or ⅓ trace] = ½(coll+kin diagonal)
- **φ** = N·(π/6)·⟨d³⟩ / V_box (from box volume; under pressure control read the live L_y)
- **I = γ̇·d·√(ρ_s/P)**, **μ_eff = τ/P**, **Φ = φ**

This is the closed loop, no new code. (Note: stress signs — LAMMPS pressure is the *negative* of the stress tensor; take magnitudes consistently.)

**DIRT (if running 3b):** DIRT **does compute a bulk virial (Love–Weber) stress tensor** — the physics is present; only the *output path* is missing.
- `VirialStress` lives in `/Users/suehr/Documents/GitHub/soil/crates/soil_core/src/virial.rs` (symmetric 3×3: `xx,yy,zz,xy,xz,yz`). It is auto-registered by `HertzMindlinContactPlugin` (`dirt_granular/src/contact.rs:66`), zeroed each step at `PreForce`, and accumulated by the contact loop as the **full pairwise virial** — normal *and* tangential contributions, with the `newton=false` half-count correction (`contact.rs:339-343, 559-567`, `add_pair(dx,dy,dz, fx,fy,fz)` with `d = pos[j]−pos[i]`). Documented convention: `P = NkT/V − trace/(3V)`. So **σ_ij = −VirialStress_ij / V_box** (plus a kinetic term if you want it).
- **The gap is narrow:** nothing in DIRT *reads `VirialStress` back out* — there is no thermo key or CSV exporter, only the in-memory `Res<VirialStress>`. Closing it is a **~20-line recorder system** modeled on `bench_plate_sinkage/main.rs:record_sinkage`: pull `Res<VirialStress>`, divide by domain volume, stream `σ_xy, σ_yy, ...` to CSV each thermo step. No core/physics change — purely an output add in the Tier-1 `main.rs`. This makes **DIRT-native Tier-1 stress feasible** if you take the option-(B)/Couette route.
- `dirt_contact_analysis` (`[contact_analysis]`) gives **coordination number** (`coord_avg/max/min`), **rattler fraction**, and the **fabric tensor** `F_ij = (1/N_c)Σ n_i n_j` (`fabric_xx..fabric_yz`, `contacts` to thermo), plus an optional per-contact CSV (`i_tag, j_tag, overlap, cx,cy,cz, nx,ny,nz`). Fabric is *directional structure*, not stress, and the per-contact CSV lacks contact force — so use `VirialStress`, not this, for stress. (Coordination/fabric are still useful diagnostics: jamming, anisotropy.)
- Alternatively, **wall-based** stress: read `force_accumulator` off the shear walls → P = F_n/A. The accumulator is currently the scalar *normal* force only (`lib.rs:339`), so wall-based τ would need the tangential wall force exposed — the `VirialStress` route avoids that.
- `dirt_measure_plane` only counts crossings/mass flux — not stress; not useful here.
- **Conclusion:** stress can come from DIRT (`VirialStress` + a small recorder) *or* from LAMMPS (already wired). For v0 the LAMMPS path is still recommended because the *shear primitive* (Lees–Edwards) is the harder missing piece, not the stress — but DIRT-native stress is no longer a blocker.

**Fits (both paths, spec §3.3):** μ(I) = μ_s + (μ_2−μ_s)/(I_0/I + 1); Φ(I) = Φ_max − (Φ_max−Φ_min)I; ρ_c = Φ_max·ρ_s; K from the separate compression test (ΔP/(ΔΦ/Φ)).

---

## 5. Tier 2 & Tier 3 — directly runnable in DIRT

### Tier 2 — Column collapse (`bench_column_collapse`)

`/Users/suehr/Documents/GitHub/dirt/examples/bench_column_collapse/{main.rs,config.toml}` runs the spec setup as-is: quasi-2D slab, loose insert → settle (stage 1) → gate wall removed (stage 2) → spread; dumps the final deposit to `data/column_collapse_results.csv`. Run:

```bash
cargo run --release --example bench_column_collapse --no-default-features -- examples/bench_column_collapse/config.toml
```

**Changes to map onto glass beads + the aspect sweep:**
1. **Grain size:** config uses `radius = 0.0015` (d = 3 mm). For the spec's d = 0.5 mm with mild polydispersity, set `radius = { distribution = "uniform", min = 0.000225, max = 0.000275 }`. This raises particle count sharply (∝ d⁻³ at fixed volume) — keep the column small or accept ~10⁴–10⁵ grains.
2. **Material:** set `restitution = 0.7`, `friction = 0.5`, `youngs_mod = 7.0e7`, `poisson_ratio = 0.25` (config currently 0.5 / 0.5 / 7e7 / 0.25 — only e and the radius distribution differ).
3. **Aspect sweep a = H₀/r₀ ∈ {0.5, 1, 2, 3, 6}** (config is the a = 2 case): vary the insert `region` height and `count`, and the gate position `point_x = r₀`. Lengthen `domain.x_high` for the long-runout (a = 6) cases. A `sweep.py` looping these (like the other benches' sweep scripts) is the clean way.
4. **Basal friction:** confirm the floor wall arrests the deposit (it should — Mindlin wall friction works, per angle-of-repose). Watch for the stale "no wall friction" comment in the config; validate the deposit comes to rest and measure the repose-consistent runout.

Measure deposit profile h(x), runout (r∞−r₀)/r₀ vs a, energy partition, failure plane, internal velocity from the CSV + VTP dumps (spec §4 five gates).

### Tier 3 — Footpad intrusion (`bench_plate_sinkage`)

`/Users/suehr/Documents/GitHub/dirt/examples/bench_plate_sinkage/{main.rs,config.toml}` is the flat-plate intrusion: a deep bed settles under enhanced gravity, a downward-facing plane wall clipped to a finite footprint (`bound_x_low/high`) drives down at constant velocity, and the recorder streams (sinkage z, reaction force F = `plate.force_accumulator`) to `data/plate_sinkage_results.csv`. Run:

```bash
cargo run --release --example bench_plate_sinkage --no-default-features -- examples/bench_plate_sinkage/config.toml
```

**Changes to map onto the spec (§5):**
1. **Grain size / material:** config uses d ≈ 2.5 mm, `friction = 0.5`, `restitution = 0.3`, `youngs_mod = 5.0e6`. Set d = 0.5 mm polydisperse, `restitution = 0.7`, and `youngs_mod ≈ 7e7` (stiffer than 5e6 → smaller dt; verify rigid-limit insensitivity). Deep bed (≥10× intruder size).
2. **Two campaigns:** (a) *quasi-static K_bear* — slow plate `velocity` (low I_mac), F(z) slope; the existing config is already in this regime (enhanced g, slow descent). (b) *velocity sweep λ* — sweep the plate `velocity` z-component across I_mac ~0.05→0.5, fit excess force to λρAv². The config's single-velocity plate is one point; `sweep.py` should loop the plate velocity (and `bound_x_*` footprint width b, already swept).
3. **Sphere intruder** (spec's second shape): replace the plate plane wall with `type = "sphere"` (`center`, `radius`, `velocity`) — supported by `dirt_wall`. Read its `force_accumulator` the same way.
4. **Two packing states** (loose/dense): prepare via the settle stage (loose insert) vs a pre-compaction stage (servo top wall, or higher restitution-damping settle).

Force is read natively — `force_accumulator` is already the mechanism (`main.rs:53-88`), no code change for the plate; sphere intruder needs the recorder to look at `walls.spheres` instead of `walls.planes`.

---

## 6. Recommended execution order, compute notes, gaps

**Order (mirrors spec §7):**
1. **Contact-law sanity** (cheap, do first): run `bench_hertz_rebound`, `bench_sliding_friction`, `bench_oblique_impact` with the glass-bead material to confirm e = 0.7 and μ_p = 0.5 reproduce target rebound/friction — this validates the contact params before any campaign run.
2. **Pilot:** one LAMMPS shear run (fix the I(P,γ̇) map, steady-state strain, dt) + one small DIRT collapse. Low cost.
3. **Tier 1 in LAMMPS** (~25 shear + 1 compression + ~7 robustness): the critical path → μ_s, μ_2, I_0, Φ(I), ρ_c, K. Embarrassingly parallel; extend `generate_runs.py`.
4. **Tier 2 in DIRT** (5 collapses): moderate (10⁴–10⁵ grains).
5. **Tier 3 in DIRT** (~28 runs): most expensive (large beds); only after 1–2 validate.

**Compute / numerics:**
- **Timestep:** spec §2.3 — Δt ≤ 0.15 t_R, t_R = πR√(ρ_s/G)/(0.1631ν+0.8766). For E = 7e7, ν = 0.25, d = 0.5 mm, R = 0.25 mm: G ≈ E/(2(1+ν)) ≈ 2.8e7, t_R ≈ π·2.5e-4·√(2500/2.8e7)/0.9 ≈ 2.6e-6 s → **Δt ≲ 4e-7 s**. The DIRT example configs use dt = 4e-6 for d = 1.5–3 mm; scale down ~10× for the smaller, stiffer glass beads. Re-derive whenever E changes (the stiffness robustness run at 10×E needs ~√10 smaller dt).
- **Softened E (5e7–1e8 Pa)** is mandatory for tractable dt and is exactly what the examples already do (`youngs_mod` 5e6–7e7). The 10×E run confirms rigid-grain insensitivity.
- **Gravity off** for the Tier-1 rheometer (homogeneous shear); g = 9.81 for Tiers 2–3 (swap to 1.62 only in the lunar phase).
- Particle counts: spec wants ~10⁴ (2D-thin) to 10⁵ (3D) for Tier 1; the DIRT Couette config above uses 12 k. At d = 0.5 mm the collapse/intrusion beds balloon — keep cells minimal and rely on periodicity.

**Gaps / risks (ranked):**
1. **No Lees–Edwards / xy-shear in `soil_deform`** — the core reason Tier 1 goes to LAMMPS. Fix = add a `DeformStyle` tilt variant (option B, §2). Medium effort.
2. **DIRT computes virial stress but never exports it** — `VirialStress` (`soil_core/src/virial.rs`) is the full Love–Weber tensor, auto-accumulated by the Hertz–Mindlin contact loop, but no thermo/CSV consumer exists. Closing it is a ~20-line recorder (read `Res<VirialStress>`/V_box → CSV), *not* a physics change. So DIRT-native Tier-1 stress is feasible; only the output plumbing is missing. (`dirt_contact_analysis` gives fabric/coordination, not stress.) LAMMPS already exports stress.
3. **Servo + tangential-velocity superposition on one wall** is unverified — needed for the DIRT-native Couette (§2A, §3b). One-line check in `dirt_wall::wall_move`.
4. **Liz's LAMMPS deck is frictionless + fixed-volume** — must add μ_p = 0.5 (Hertz–Mindlin tangential) and pressure control (y servo wall) to be the spec rheometer (§2). The fixed-volume φ-sweep is a working fallback that still yields μ(I), Φ(I).
5. **Stale "no wall friction" comment** in `bench_column_collapse/config.toml` contradicts the working Mindlin wall friction in `bench_angle_of_repose` — validate basal arrest empirically on the first collapse run.

**Bottom line:** Tiers 2 and 3 are *runnable today in DIRT* with config-only edits (grain size, e, polydispersity, sweeps). Tier 1 — the pressure-controlled homogeneous-shear rheometer — has no DIRT example and a missing shear primitive, so run it in **LAMMPS** using Liz's existing shear cell, after adding interparticle friction (μ_p = 0.5) and pressure control. Build the DIRT Lees–Edwards + stress-output path later if a single-code pipeline is wanted for the regolith/lunar phase.

*Drafted 2026-06-16. Companion to `docs/dem-campaign.md` (the spec) and `docs/physics-design.md` (the consumer).*
