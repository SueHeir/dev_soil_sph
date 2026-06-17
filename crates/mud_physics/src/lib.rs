//! # mud_physics — the per-step SPH pipeline
//!
//! The MUD physics tier's per-timestep work, riding the SOIL substrate. Unlike
//! DEM's single pairwise-contact pass, SPH needs **two neighbor passes** with a
//! per-particle constitutive update and a halo exchange between them
//! (`docs/architecture.md` §4):
//!
//! ```text
//! PreForce: mud_density_velgrad   (gather: dρ/dt and L = ∇v from neighbors)
//!         → mud_integrate_density (ρ += dρ/dt · dt)
//!         → mud_constitutive_update (p = EOS(ρ); μ(I) return map → dev_stress)
//!         → forward_comm_borders   ★ push fresh ρ, p, σ to ghosts
//! Force  : mud_momentum           (gather: dv/dt = Σ mⱼ(σᵢ/ρᵢ² + σⱼ/ρⱼ²)·∇W)
//! ```
//!
//! ## Gather formulation
//! Uses the **full neighbor list** (`[neighbor] newton = false`): each owner `i`
//! sums over all its neighbors `j` (local or ghost) and writes only `i`, so no
//! reverse communication is needed and each pass is self-contained. The systems
//! assert `!neighbor.newton`.
//!
//! ## v0 simplifications
//! - 3D Wendland C2 kernel, uniform smoothing length per particle.
//! - Continuity density (ρ carried as integrated state).
//! - Body force (gravity) is added separately (e.g. `dirt_fixes::GravityPlugin`).

// Index-based loops are intentional for the small fixed-size tensor algebra.
#![allow(clippy::needless_range_loop)]

use grass_app::prelude::*;
use grass_scheduler::prelude::*;

use serde::Deserialize;
use soil_core::{
    forward_comm_borders, Atom, AtomDataRegistry, Config, Neighbor, ParticleSimScheduleSet,
    RunConfig, ScheduleSetupSet,
};

use mud_atom::{MudAtom, MudAtomPlugin, MudMaterialTable};
use mud_constitutive::{kt_cooling_rate, kt_production_rate, two_branch_stress};
use mud_kernel::Kernel;

/// Spatial dimension of the kernel (v0: 3D).
const KERNEL: Kernel = Kernel::Dim3;

// ── Small tensor helpers ─────────────────────────────────────────────────────

/// Reconstruct the full Cauchy stress `σ = −p I + s` (as a `Sym3` `[xx,yy,zz,xy,xz,yz]`)
/// for particle `k` from its stored pressure and deviatoric stress.
#[inline]
fn sigma(sph: &MudAtom, k: usize) -> [f64; 6] {
    let p = sph.pressure[k];
    let s = sph.dev_stress[k];
    [s[0] - p, s[1] - p, s[2] - p, s[3], s[4], s[5]]
}

/// Symmetric-3×3 times vector: `(a · g)` with `a = [xx,yy,zz,xy,xz,yz]`.
#[inline]
fn sym3_mat_vec(a: &[f64; 6], g: &[f64; 3]) -> [f64; 3] {
    [
        a[0] * g[0] + a[3] * g[1] + a[4] * g[2],
        a[3] * g[0] + a[1] * g[1] + a[5] * g[2],
        a[4] * g[0] + a[5] * g[1] + a[2] * g[2],
    ]
}

// ── Pass 1: density rate + velocity gradient (PreForce) ──────────────────────

/// Gather pass: for each owner `i`, accumulate the continuity density rate
/// `dρ/dt_i = Σⱼ mⱼ (vᵢ − vⱼ)·∇Wᵢⱼ` and the velocity gradient
/// `Lᵢ = Σⱼ (mⱼ/ρⱼ)(vⱼ − vᵢ) ⊗ ∇Wᵢⱼ`. Both accumulators are zeroed by SOIL each step.
pub fn mud_density_velgrad(atoms: Res<Atom>, neighbor: Res<Neighbor>, registry: Res<AtomDataRegistry>) {
    assert!(
        !neighbor.newton,
        "MUD uses a gather formulation; set [neighbor] newton = false"
    );
    let mut sph = registry.expect_mut::<MudAtom>("mud_density_velgrad");
    let nlocal = atoms.nlocal as usize;

    for (i, j) in neighbor.pairs(nlocal) {
        let dx = [
            atoms.pos[i][0] - atoms.pos[j][0],
            atoms.pos[i][1] - atoms.pos[j][1],
            atoms.pos[i][2] - atoms.pos[j][2],
        ];
        let gradw = KERNEL.grad_w(dx, sph.h[i]);
        if gradw == [0.0, 0.0, 0.0] {
            continue;
        }
        let vol_j = sph.particle_mass[j] / sph.density[j];
        let mj = sph.particle_mass[j];
        // v_j − v_i
        let dvji = [
            atoms.vel[j][0] - atoms.vel[i][0],
            atoms.vel[j][1] - atoms.vel[i][1],
            atoms.vel[j][2] - atoms.vel[i][2],
        ];
        // L_i += vol_j (v_j − v_i) ⊗ ∇W   (row-major a*3+b)
        for a in 0..3 {
            for b in 0..3 {
                sph.velgrad[i][a * 3 + b] += vol_j * dvji[a] * gradw[b];
            }
        }
        // dρ/dt_i += m_j (v_i − v_j)·∇W = −m_j (v_j − v_i)·∇W
        let dot = dvji[0] * gradw[0] + dvji[1] * gradw[1] + dvji[2] * gradw[2];
        sph.drho_dt[i] -= mj * dot;
    }
}

// ── Density integration (PreForce, after pass 1) ─────────────────────────────

/// Integrate the continuity density: `ρ_i += dρ/dt_i · dt`.
pub fn mud_integrate_density(atoms: Res<Atom>, registry: Res<AtomDataRegistry>) {
    let mut sph = registry.expect_mut::<MudAtom>("mud_integrate_density");
    let dt = atoms.dt;
    for i in 0..atoms.nlocal as usize {
        sph.density[i] += sph.drho_dt[i] * dt;
    }
}

// ── Granular-temperature update (PreForce, after density integration) ────────

/// Evolve the granular temperature `T` (`physics-design.md` §11.2). v0 scaffolding:
/// the **homogeneous** balance `dT/dt = −Γ` (inelastic collisional cooling →
/// Haff's law). Production (`σ_KT : D`) and conduction (`∇·κ∇T`) are later
/// increments — flagged below. Per-particle; clamps `T ≥ 0`.
pub fn mud_temperature_update(
    atoms: Res<Atom>,
    registry: Res<AtomDataRegistry>,
    table: Res<MudMaterialTable>,
) {
    let mut sph = registry.expect_mut::<MudAtom>("mud_temperature_update");
    let dt = atoms.dt;
    for i in 0..atoms.nlocal as usize {
        let mat = &table.params[atoms.atom_type[i] as usize];
        let rho = sph.density[i];
        let t = sph.temperature[i];
        let l = sph.velgrad[i];
        // dT/dt = production (collisional shear heating) − dissipation.
        // TODO(v1): + conduction ∇·(κ∇T) (needs κ from the inhomogeneous DEM rig).
        let dt_dt = kt_production_rate(rho, t, &l, mat) + kt_cooling_rate(rho, t, mat);
        sph.temperature[i] = (t + dt_dt * dt).max(0.0);
    }
}

// ── Constitutive update (PreForce, after integration) ────────────────────────

/// Per-particle stress update: `p = EOS(ρ)` and the μ(I) return map evolving the
/// deviatoric stress (delegates to `mud_constitutive::update_stress`).
pub fn mud_constitutive_update(
    atoms: Res<Atom>,
    registry: Res<AtomDataRegistry>,
    table: Res<MudMaterialTable>,
) {
    let mut sph = registry.expect_mut::<MudAtom>("mud_constitutive_update");
    let dt = atoms.dt;
    for i in 0..atoms.nlocal as usize {
        let mat = &table.params[atoms.atom_type[i] as usize];
        let s_n = sph.dev_stress_elastic[i];
        let l = sph.velgrad[i];
        let rho = sph.density[i];
        let t = sph.temperature[i];
        let out = two_branch_stress(&s_n, &l, rho, t, dt, mat);
        sph.pressure[i] = out.pressure; // total: p_contact + p_KT
        sph.dev_stress[i] = out.dev_total; // total: s_contact + 2η_KT D' (for force)
        sph.dev_stress_elastic[i] = out.dev_elastic; // persistent s_contact (next s_n)
    }
}

// ── Pass 2: momentum (Force) ─────────────────────────────────────────────────

/// Monaghan artificial-viscosity coefficients (linear, quadratic). Provides the
/// dissipation that stabilizes SPH (a perfect lattice is otherwise an unstable
/// equilibrium). `docs/physics-design.md` §4. v0: hard-coded; configurable later.
const AV_ALPHA: f64 = 1.0;
const AV_BETA: f64 = 2.0;

/// Gather pass: stress-divergence force plus Monaghan artificial viscosity.
/// For each owner `i`,
/// `force_i += Σⱼ mᵢ mⱼ (σᵢ/ρᵢ² + σⱼ/ρⱼ² − Πᵢⱼ I)·∇Wᵢⱼ`. Gravity is added separately.
pub fn mud_momentum(
    mut atoms: ResMut<Atom>,
    neighbor: Res<Neighbor>,
    registry: Res<AtomDataRegistry>,
    table: Res<MudMaterialTable>,
) {
    let sph = registry.expect::<MudAtom>("mud_momentum");
    let nlocal = atoms.nlocal as usize;

    for (i, j) in neighbor.pairs(nlocal) {
        let dx = [
            atoms.pos[i][0] - atoms.pos[j][0],
            atoms.pos[i][1] - atoms.pos[j][1],
            atoms.pos[i][2] - atoms.pos[j][2],
        ];
        let h = sph.h[i];
        let gradw = KERNEL.grad_w(dx, h);
        if gradw == [0.0, 0.0, 0.0] {
            continue;
        }
        let mi = sph.particle_mass[i];
        let mj = sph.particle_mass[j];
        let rho_i = sph.density[i];
        let rho_j = sph.density[j];
        let rhoi2 = rho_i * rho_i;
        let rhoj2 = rho_j * rho_j;
        let sig_i = sigma(&sph, i);
        let sig_j = sigma(&sph, j);
        let mut term: [f64; 6] = std::array::from_fn(|k| sig_i[k] / rhoi2 + sig_j[k] / rhoj2);

        // Monaghan artificial viscosity Π_ij: active only for approaching pairs
        // (v_ij·r_ij < 0). Adds an isotropic dissipative pressure (−Π onto the
        // stress-tensor diagonal, since σ = −p I + s here uses compression-positive p).
        let vij = [
            atoms.vel[i][0] - atoms.vel[j][0],
            atoms.vel[i][1] - atoms.vel[j][1],
            atoms.vel[i][2] - atoms.vel[j][2],
        ];
        let vdotr = vij[0] * dx[0] + vij[1] * dx[1] + vij[2] * dx[2];
        if vdotr < 0.0 {
            let r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
            let mu = h * vdotr / (r2 + 0.01 * h * h);
            // Use the owner's material sound speed. `i` is always local; base
            // `atom_type` is NOT forward-communicated, so `atom_type[j]` is invalid
            // for ghost `j`. Single-material v0 → c_i == c_j anyway.
            let c_bar = table.params[atoms.atom_type[i] as usize].sound_speed();
            let rho_bar = 0.5 * (rho_i + rho_j);
            let pi_ij = (-AV_ALPHA * c_bar * mu + AV_BETA * mu * mu) / rho_bar;
            // a_i += Σ_j m_j (−Π_ij I)·∇W  ⇒  subtract Π_ij from the diagonal.
            term[0] -= pi_ij;
            term[1] -= pi_ij;
            term[2] -= pi_ij;
        }

        let tg = sym3_mat_vec(&term, &gradw);
        for d in 0..3 {
            atoms.force[i][d] += mi * mj * tg[d];
        }
    }
}

// ── Plugin ───────────────────────────────────────────────────────────────────

// ── Gravity (body force) ─────────────────────────────────────────────────────

#[derive(Deserialize, Clone, Default)]
#[serde(deny_unknown_fields)]
struct GravityConfig {
    #[serde(default)]
    gx: f64,
    #[serde(default)]
    gy: f64,
    #[serde(default)]
    gz: f64,
}

/// The gravity acceleration vector resource.
pub struct MudGravity {
    pub g: [f64; 3],
}

/// Adds the gravity body force `m·g` to each particle, at the Force phase.
pub fn mud_gravity(
    mut atoms: ResMut<Atom>,
    gravity: Res<MudGravity>,
    registry: Res<AtomDataRegistry>,
) {
    let sph = registry.expect::<MudAtom>("mud_gravity");
    let g = gravity.g;
    for i in 0..atoms.nlocal as usize {
        let m = sph.particle_mass[i];
        atoms.force[i][0] += m * g[0];
        atoms.force[i][1] += m * g[1];
        atoms.force[i][2] += m * g[2];
    }
}

/// Registers the gravity body force from `[gravity]` (gx, gy, gz).
pub struct MudGravityPlugin;

impl Plugin for MudGravityPlugin {
    fn default_config(&self) -> Option<&str> {
        Some("[gravity]\ngx = 0.0\ngy = 0.0\ngz = -9.81")
    }

    fn build(&self, app: &mut App) {
        let cfg = Config::load::<GravityConfig>(app, "gravity");
        app.add_resource(MudGravity {
            g: [cfg.gx, cfg.gy, cfg.gz],
        });
        app.add_update_system(mud_gravity.label("mud_gravity"), ParticleSimScheduleSet::Force);
    }
}

// ── Boundary freeze ──────────────────────────────────────────────────────────

/// Freeze boundary particles: zero their force and velocity each step so they
/// stay fixed in place, while still participating in the SPH sums (they develop
/// pressure via continuity/EOS and support the fluid). Runs at PostForce, after
/// all force contributions (SPH stress + gravity). No-op when there are no
/// boundary particles.
pub fn mud_freeze_boundary(mut atoms: ResMut<Atom>, registry: Res<AtomDataRegistry>) {
    let sph = registry.expect::<MudAtom>("mud_freeze_boundary");
    for i in 0..atoms.nlocal as usize {
        if sph.is_boundary[i] > 0.5 {
            atoms.force[i] = [0.0; 3];
            atoms.vel[i] = [0.0; 3];
        }
    }
}

/// Copy the `[[run]] dt` into `Atom::dt` (the Verlet integrator reads `Atom::dt`,
/// which otherwise defaults to 1.0 — a catastrophic CFL violation). Mirrors POND's
/// `set_timestep` / DIRT's `calculate_delta_time`.
fn mud_set_timestep(mut atoms: ResMut<Atom>, run_config: Res<RunConfig>) {
    let dt = run_config.current_stage(0).dt;
    if dt > 0.0 {
        atoms.dt = dt;
    }
}

/// Registers the MUD per-step SPH systems in their schedule phases.
pub struct MudPhysicsPlugin;

impl Plugin for MudPhysicsPlugin {
    fn dependencies(&self) -> Vec<std::any::TypeId> {
        grass_app::type_ids![MudAtomPlugin]
    }

    fn provides(&self) -> Vec<&str> {
        vec!["mud_forces"]
    }

    fn requires(&self) -> Vec<&str> {
        vec!["mud_particles", "neighbor_list"]
    }

    fn build(&self, app: &mut App) {
        app.add_setup_system(mud_set_timestep, ScheduleSetupSet::PostSetup);
        app.add_update_system(
            mud_density_velgrad.label("mud_density"),
            ParticleSimScheduleSet::PreForce,
        )
        .add_update_system(
            mud_integrate_density.label("mud_integ_rho").after("mud_density"),
            ParticleSimScheduleSet::PreForce,
        )
        .add_update_system(
            mud_temperature_update.label("mud_temp").after("mud_integ_rho"),
            ParticleSimScheduleSet::PreForce,
        )
        .add_update_system(
            mud_constitutive_update.label("mud_const").after("mud_temp"),
            ParticleSimScheduleSet::PreForce,
        )
        // ★ mid-step halo: push freshly-computed ρ, p, σ to ghosts before the
        // momentum pass reads them (soil_core::forward_comm_borders, re-registered
        // here in addition to its PreNeighbor registration).
        .add_update_system(
            forward_comm_borders.label("mud_halo").after("mud_const"),
            ParticleSimScheduleSet::PreForce,
        )
        .add_update_system(
            mud_momentum.label("mud_force"),
            ParticleSimScheduleSet::Force,
        )
        .add_update_system(
            mud_freeze_boundary.label("mud_freeze"),
            ParticleSimScheduleSet::PostForce,
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sym3_mat_vec_matches_dense() {
        // a = [[1,4,5],[4,2,6],[5,6,3]], g = [1,1,1] → row sums [10,12,14]
        let a = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let g = [1.0, 1.0, 1.0];
        assert_eq!(sym3_mat_vec(&a, &g), [10.0, 12.0, 14.0]);
    }

    #[test]
    fn sigma_reconstructs_minus_p_i_plus_s() {
        let mut sph = MudAtom::new();
        sph.pressure.push(10.0);
        sph.dev_stress.push([1.0, 2.0, 3.0, 0.5, 0.6, 0.7]);
        // σ_xx = s_xx − p = 1 − 10 = −9; off-diagonals unchanged
        assert_eq!(sigma(&sph, 0), [-9.0, -8.0, -7.0, 0.5, 0.6, 0.7]);
    }
}
