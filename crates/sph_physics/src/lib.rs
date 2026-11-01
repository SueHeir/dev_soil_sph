//! # sph_physics — the per-step SPH pipeline
//!
//! The dev_soil_sph physics tier's per-timestep work, riding the SOIL substrate. Unlike
//! DEM's single pairwise-contact pass, SPH needs **two neighbor passes** with a
//! per-particle constitutive update and a halo exchange between them
//! (`docs/architecture.md` §4):
//!
//! ```text
//! PreForce: sph_density_velgrad   (gather: dρ/dt and L = ∇v from neighbors)
//!         → sph_integrate_density (ρ += dρ/dt · dt)
//!         → sph_constitutive_update (p = EOS(ρ); μ(I) return map → dev_stress)
//!         → forward_comm_borders   ★ push fresh ρ, p, σ to ghosts
//! Force  : sph_momentum           (gather: dv/dt = Σ mⱼ(σᵢ/ρᵢ² + σⱼ/ρⱼ²)·∇W)
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

use soil_core::{forward_comm_borders, Atom, AtomDataRegistry, Neighbor, ParticleSimScheduleSet};

use sph_atom::{SphAtom, SphAtomPlugin, SphMaterialTable};
use sph_constitutive::{
    kt_conductivity, kt_cooling_rate, kt_production_rate, two_branch_stress,
};
use sph_kernel::Kernel;

/// Spatial dimension of the kernel (v0: 3D).
const KERNEL: Kernel = Kernel::Dim3;

// ── Small tensor helpers ─────────────────────────────────────────────────────

/// Reconstruct the full Cauchy stress `σ = −p I + s` (as a `Sym3` `[xx,yy,zz,xy,xz,yz]`)
/// for particle `k` from its stored pressure and deviatoric stress.
#[inline]
fn sigma(sph: &SphAtom, k: usize) -> [f64; 6] {
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
pub fn sph_density_velgrad(atoms: Res<Atom>, neighbor: Res<Neighbor>, registry: Res<AtomDataRegistry>) {
    assert!(
        !neighbor.newton,
        "dev_soil_sph uses a gather formulation; set [neighbor] newton = false"
    );
    let mut sph = registry.expect_mut::<SphAtom>("sph_density_velgrad");
    let nlocal = atoms.nlocal() as usize;

    // Row-wise gather: everything that depends only on `i` is read once, and the
    // accumulators are written once per particle rather than once per pair.
    for i in 0..nlocal {
        let xi = atoms.pos[i];
        let vi = atoms.vel[i];
        let hi = sph.h[i];
        let mut l = sph.velgrad[i];
        let mut drho = sph.drho_dt[i];
        for &j in neighbor.neighbors_of(i) {
            let j = j as usize;
            let dx = [xi[0] - atoms.pos[j][0], xi[1] - atoms.pos[j][1], xi[2] - atoms.pos[j][2]];
            let gradw = KERNEL.grad_w(dx, hi);
            if gradw == [0.0, 0.0, 0.0] {
                continue;
            }
            let mj = sph.particle_mass[j];
            let vol_j = mj / sph.density[j];
            // v_j − v_i
            let vj = atoms.vel[j];
            let dvji = [vj[0] - vi[0], vj[1] - vi[1], vj[2] - vi[2]];
            // L_i += vol_j (v_j − v_i) ⊗ ∇W   (row-major a*3+b)
            for a in 0..3 {
                for b in 0..3 {
                    l[a * 3 + b] += vol_j * dvji[a] * gradw[b];
                }
            }
            // dρ/dt_i += m_j (v_i − v_j)·∇W = −m_j (v_j − v_i)·∇W
            let dot = dvji[0] * gradw[0] + dvji[1] * gradw[1] + dvji[2] * gradw[2];
            drho -= mj * dot;
        }
        sph.velgrad[i] = l;
        sph.drho_dt[i] = drho;
    }
}

// ── Density integration (PreForce, after pass 1) ─────────────────────────────

/// Integrate the continuity density: `ρ_i += dρ/dt_i · dt`.
pub fn sph_integrate_density(atoms: Res<Atom>, registry: Res<AtomDataRegistry>) {
    let mut sph = registry.expect_mut::<SphAtom>("sph_integrate_density");
    let dt = atoms.dt;
    for i in 0..atoms.nlocal() as usize {
        sph.density[i] += sph.drho_dt[i] * dt;
    }
}

// ── Granular-temperature conduction (PreForce gather, before the T update) ───

/// SPH-Laplacian gather for the conduction term `∇·(κ∇T)` (Cleary/Brookshaw
/// form), accumulated into the `lap_t` column. For each owner `i`:
/// `∇·(κ∇T)_i = Σⱼ (mⱼ/ρⱼ)(κᵢ+κⱼ)(Tᵢ−Tⱼ)(rᵢⱼ·∇Wᵢⱼ)/(‖rᵢⱼ‖²+εh²)`.
/// Reduces to `κ∇²T` for uniform κ; the sign diffuses heat hot→cold. κ uses the
/// owner's material (base `atom_type` isn't forwarded — single-material v0).
pub fn sph_conduction(atoms: Res<Atom>, neighbor: Res<Neighbor>, registry: Res<AtomDataRegistry>, table: Res<SphMaterialTable>) {
    let mut sph = registry.expect_mut::<SphAtom>("sph_conduction");
    let nlocal = atoms.nlocal() as usize;
    for (i, j) in neighbor.pairs(nlocal) {
        let ti = sph.temperature[i];
        let tj = sph.temperature[j];
        // No flux if both cold (κ = 0 there) — cheap skip.
        if ti <= 0.0 && tj <= 0.0 {
            continue;
        }
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
        let mat = &table.params[atoms.atom_type[i] as usize];
        let ki = kt_conductivity(sph.density[i], ti, mat);
        let kj = kt_conductivity(sph.density[j], tj, mat);
        let r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
        let rdotgw = dx[0] * gradw[0] + dx[1] * gradw[1] + dx[2] * gradw[2];
        let vol_j = sph.particle_mass[j] / sph.density[j];
        sph.lap_t[i] += vol_j * (ki + kj) * (ti - tj) * rdotgw / (r2 + 0.01 * h * h);
    }
}

// ── Granular-temperature update (PreForce, after density integration) ────────

/// Evolve the granular temperature `T` (`physics-design.md` §11.2). v0 scaffolding:
/// the **homogeneous** balance `dT/dt = −Γ` (inelastic collisional cooling →
/// Haff's law). Production (`σ_KT : D`) and conduction (`∇·κ∇T`) are later
/// increments — flagged below. Per-particle; clamps `T ≥ 0`.
pub fn sph_temperature_update(
    atoms: Res<Atom>,
    registry: Res<AtomDataRegistry>,
    table: Res<SphMaterialTable>,
) {
    let mut sph = registry.expect_mut::<SphAtom>("sph_temperature_update");
    let dt = atoms.dt;
    for i in 0..atoms.nlocal() as usize {
        let mat = &table.params[atoms.atom_type[i] as usize];
        let rho = sph.density[i];
        let t = sph.temperature[i];
        let l = sph.velgrad[i];
        // dT/dt = production (shear heating) − dissipation + (2/3ρ)∇·(κ∇T).
        let conduction = if rho > 0.0 { (2.0 / (3.0 * rho)) * sph.lap_t[i] } else { 0.0 };
        let dt_dt = kt_production_rate(rho, t, &l, mat) + kt_cooling_rate(rho, t, mat) + conduction;
        sph.temperature[i] = (t + dt_dt * dt).max(0.0);
    }
}

// ── Constitutive update (PreForce, after integration) ────────────────────────

/// Per-particle stress update: `p = EOS(ρ)` and the μ(I) return map evolving the
/// deviatoric stress (delegates to `sph_constitutive::update_stress`).
pub fn sph_constitutive_update(
    atoms: Res<Atom>,
    registry: Res<AtomDataRegistry>,
    table: Res<SphMaterialTable>,
) {
    let mut sph = registry.expect_mut::<SphAtom>("sph_constitutive_update");
    let dt = atoms.dt;
    for i in 0..atoms.nlocal() as usize {
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
pub fn sph_momentum(
    mut atoms: ResMut<Atom>,
    neighbor: Res<Neighbor>,
    registry: Res<AtomDataRegistry>,
    table: Res<SphMaterialTable>,
) {
    let sph = registry.expect::<SphAtom>("sph_momentum");
    let nlocal = atoms.nlocal() as usize;

    // Row-wise gather: σᵢ/ρᵢ², the sound speed and the other `i`-only terms are
    // hoisted out of the neighbor loop; the force is written once per particle.
    for i in 0..nlocal {
        let xi = atoms.pos[i];
        let vi = atoms.vel[i];
        let h = sph.h[i];
        let mi = sph.particle_mass[i];
        let rho_i = sph.density[i];
        let sig_i = sigma(&sph, i);
        let rhoi2 = rho_i * rho_i;
        let s_i: [f64; 6] = std::array::from_fn(|k| sig_i[k] / rhoi2);
        // Owner's material sound speed. `i` is always local; base `atom_type` is
        // NOT forward-communicated, so `atom_type[j]` is invalid for ghost `j`.
        // Single-material v0 → c_i == c_j anyway.
        let c_bar = table.params[atoms.atom_type[i] as usize].sound_speed();
        let mut f = [0.0f64; 3];

        for &j in neighbor.neighbors_of(i) {
            let j = j as usize;
            let dx = [xi[0] - atoms.pos[j][0], xi[1] - atoms.pos[j][1], xi[2] - atoms.pos[j][2]];
            let gradw = KERNEL.grad_w(dx, h);
            if gradw == [0.0, 0.0, 0.0] {
                continue;
            }
            let mj = sph.particle_mass[j];
            let rho_j = sph.density[j];
            let rhoj2 = rho_j * rho_j;
            let sig_j = sigma(&sph, j);
            let mut term: [f64; 6] = std::array::from_fn(|k| s_i[k] + sig_j[k] / rhoj2);

            // Monaghan artificial viscosity Π_ij: active only for approaching pairs
            // (v_ij·r_ij < 0). Adds an isotropic dissipative pressure (−Π onto the
            // stress-tensor diagonal, since σ = −p I + s here uses compression-positive p).
            let vj = atoms.vel[j];
            let vij = [vi[0] - vj[0], vi[1] - vj[1], vi[2] - vj[2]];
            let vdotr = vij[0] * dx[0] + vij[1] * dx[1] + vij[2] * dx[2];
            if vdotr < 0.0 {
                let r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
                let mu = h * vdotr / (r2 + 0.01 * h * h);
                let rho_bar = 0.5 * (rho_i + rho_j);
                let pi_ij = (-AV_ALPHA * c_bar * mu + AV_BETA * mu * mu) / rho_bar;
                // a_i += Σ_j m_j (−Π_ij I)·∇W  ⇒  subtract Π_ij from the diagonal.
                term[0] -= pi_ij;
                term[1] -= pi_ij;
                term[2] -= pi_ij;
            }

            let tg = sym3_mat_vec(&term, &gradw);
            for d in 0..3 {
                f[d] += mi * mj * tg[d];
            }
        }
        for d in 0..3 {
            atoms.force[i][d] += f[d];
        }
    }
}

// ── Plugin ───────────────────────────────────────────────────────────────────

// ── Gravity (body force) ─────────────────────────────────────────────────────

// Gravity is SOIL's shared `soil_fixes::GravityPlugin` (`[gravity]` gx/gy/gz),
// the same plugin DIRT uses. Rigid boundary particles receive `m g` like every
// other particle, but `sph_freeze_boundary` zeroes their force, and it removes
// the driven set's own weight from the measured plate reaction.
pub use soil_fixes::{GravityConfig, GravityPlugin};

// ── Boundary freeze ──────────────────────────────────────────────────────────

/// Net reaction of the bed on the **driven** rigid boundary (footpad), plus its
/// mean height — updated each step by [`sph_freeze_boundary`]. `force` is the
/// SPH bed force (gravity excluded); `z` is the footpad's mean height (sinkage =
/// initial − current).
#[derive(Default, Clone, Copy)]
pub struct SphPlateForce {
    pub force: [f64; 3],
    pub z: f64,
    pub n: usize,
}

/// Freeze boundary particles: pin each to its prescribed velocity (`boundary_vel`
/// — `[0,0,0]` static floor/wall, or a driven footpad velocity) and zero its
/// force, so it moves rigidly while still participating in the SPH sums. For the
/// **driven** set, first accumulate the net bed reaction (the force the bed exerts
/// before it is zeroed) into [`SphPlateForce`]. Runs at PostForce, after all force
/// contributions. No-op when there are no boundary particles.
pub fn sph_freeze_boundary(
    mut atoms: ResMut<Atom>,
    registry: Res<AtomDataRegistry>,
    mut plate: ResMut<SphPlateForce>,
    gravity: Option<Res<GravityConfig>>,
) {
    let sph = registry.expect::<SphAtom>("sph_freeze_boundary");
    let g = gravity.map_or([0.0; 3], |g| [g.gx, g.gy, g.gz]);
    let mut fsum = [0.0f64; 3];
    let mut zsum = 0.0f64;
    let mut npl = 0usize;
    for i in 0..atoms.nlocal() as usize {
        if sph.is_boundary[i] <= 0.5 {
            continue;
        }
        let bv = sph.boundary_vel[i];
        let driven = bv[0] != 0.0 || bv[1] != 0.0 || bv[2] != 0.0;
        if driven {
            // Reaction = bed force on the footpad this step (before zeroing),
            // excluding the footpad's own weight added by the gravity plugin.
            let m = sph.particle_mass[i];
            for d in 0..3 {
                fsum[d] += atoms.force[i][d] - m * g[d];
            }
            zsum += atoms.pos[i][2];
            npl += 1;
        }
        atoms.force[i] = [0.0; 3];
        atoms.vel[i] = bv; // pin to prescribed velocity (static or driven)
    }
    if npl > 0 {
        *plate = SphPlateForce { force: fsum, z: zsum / npl as f64, n: npl };
    }
}

/// Registers the dev_soil_sph per-step SPH systems in their schedule phases.
pub struct SphPhysicsPlugin;

impl Plugin for SphPhysicsPlugin {
    fn dependencies(&self) -> Vec<std::any::TypeId> {
        grass_app::type_ids![SphAtomPlugin]
    }

    fn provides(&self) -> Vec<&str> {
        vec!["sph_forces"]
    }

    fn requires(&self) -> Vec<&str> {
        vec!["sph_particles", "neighbor_list"]
    }

    fn build(&self, app: &mut App) {
        app.add_resource(SphPlateForce::default());
        app.add_update_system(
            sph_density_velgrad.label("sph_density"),
            ParticleSimScheduleSet::PreForce,
        )
        .add_update_system(
            sph_integrate_density.label("sph_integ_rho").after("sph_density"),
            ParticleSimScheduleSet::PreForce,
        )
        .add_update_system(
            sph_conduction.label("sph_conduction").after("sph_integ_rho"),
            ParticleSimScheduleSet::PreForce,
        )
        .add_update_system(
            sph_temperature_update.label("sph_temp").after("sph_conduction"),
            ParticleSimScheduleSet::PreForce,
        )
        .add_update_system(
            sph_constitutive_update.label("sph_const").after("sph_temp"),
            ParticleSimScheduleSet::PreForce,
        )
        // ★ mid-step halo: push freshly-computed ρ, p, σ to ghosts before the
        // momentum pass reads them (soil_core::forward_comm_borders, re-registered
        // here in addition to its PreNeighbor registration).
        .add_update_system(
            forward_comm_borders.label("sph_halo").after("sph_const"),
            ParticleSimScheduleSet::PreForce,
        )
        .add_update_system(
            sph_momentum.label("sph_force"),
            ParticleSimScheduleSet::Force,
        )
        .add_update_system(
            sph_freeze_boundary.label("sph_freeze"),
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
        let mut sph = SphAtom::new();
        sph.pressure.push(10.0);
        sph.dev_stress.push([1.0, 2.0, 3.0, 0.5, 0.6, 0.7]);
        // σ_xx = s_xx − p = 1 − 10 = −9; off-diagonals unchanged
        assert_eq!(sigma(&sph, 0), [-9.0, -8.0, -7.0, 0.5, 0.6, 0.7]);
    }
}
