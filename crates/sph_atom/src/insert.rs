//! Particle insertion for dev_soil_sph — regular-lattice block fill.
//!
//! SPH wants an *ordered* initial packing (not the random placement DEM uses), so
//! each `[[sph.insert]]` block fills its region on a regular cubic lattice of the
//! given spacing. Per-particle mass is `ρ₀ · Δ³` (the lattice cell volume) and the
//! neighbor `cutoff_radius` is set to the kernel support `2h`.

use grass_app::prelude::*;
use grass_scheduler::prelude::*;

use soil_core::{Atom, AtomDataRegistry, Domain, ParticleStore, ScheduleSetupSet};

use sph_constitutive::pressure;

use crate::{SphAtom, SphConfig, SphMaterialTable};

/// Inserts dev_soil_sph particles at setup time from `[[sph.insert]]` blocks.
pub struct SphAtomInsertPlugin;

impl Plugin for SphAtomInsertPlugin {
    fn build(&self, app: &mut App) {
        // Run after domain decomposition so each rank's subdomain bounds exist.
        app.add_setup_system(
            sph_insert_atoms.after("domain_read_input"),
            ScheduleSetupSet::Setup,
        );
    }
}

/// Lattice coordinates along one axis: cell-centered points spaced by `spacing`
/// within `[min, max]`. A too-thin axis (extent < spacing) gets a single centered
/// layer, so quasi-2D slabs work.
fn lattice_axis(min: f64, max: f64, spacing: f64) -> Vec<f64> {
    let n = ((max - min) / spacing).floor() as i64;
    if n < 1 {
        vec![0.5 * (min + max)]
    } else {
        (0..n).map(|k| min + (k as f64 + 0.5) * spacing).collect()
    }
}

/// Half-open subdomain ownership `[low, high)` per axis (matches `exchange()`).
fn owns_position(domain: &Domain, pos: &[f64; 3]) -> bool {
    (0..3).all(|d| pos[d] >= domain.sub_domain_low[d] && pos[d] < domain.sub_domain_high[d])
}

/// Append one SPH particle to the core `Atom` and every registered extension
/// through SOIL's `ParticleStore` (keeping every column row-aligned), then
/// fill in the core and `SphAtom` rows.
#[allow(clippy::too_many_arguments)]
fn insert_particle(
    atom: &mut Atom,
    registry: &AtomDataRegistry,
    pos: [f64; 3],
    vel: [f64; 3],
    mass: f64,
    density: f64,
    pressure0: f64,
    temperature0: f64,
    h: f64,
    is_boundary: f64,
    boundary_vel: [f64; 3],
    mat_idx: u32,
    tag: u32,
) {
    let global_natoms = atom
        .natoms()
        .checked_add(1)
        .expect("global particle count overflow during SPH insertion");
    ParticleStore::new(atom, registry)
        .push_default_local(global_natoms)
        .expect("registered SPH rows must accept transactional insertion");
    let i = atom.len() - 1;
    atom.tag[i] = tag;
    atom.origin_index[i] = 0;
    // Neighbor pair cutoff is (r_i + r_j)·skin; with uniform h this gives the
    // kernel support 2h when cutoff_radius = h.
    atom.cutoff_radius[i] = h;
    atom.image[i] = [0, 0, 0];
    atom.pos[i] = pos;
    atom.vel[i] = vel;
    atom.force[i] = [0.0; 3];
    atom.mass[i] = mass;
    atom.inv_mass[i] = 1.0 / mass;
    atom.atom_type[i] = mat_idx;

    let mut sph = registry.expect_mut::<SphAtom>("sph_insert_atoms");
    sph.h[i] = h;
    sph.density[i] = density;
    sph.pressure[i] = pressure0;
    sph.dev_stress[i] = [0.0; 6];
    sph.dev_stress_elastic[i] = [0.0; 6];
    sph.temperature[i] = temperature0;
    sph.velgrad[i] = [0.0; 9];
    sph.drho_dt[i] = 0.0;
    sph.lap_t[i] = 0.0;
    sph.particle_mass[i] = mass;
    sph.is_boundary[i] = is_boundary;
    sph.boundary_vel[i] = boundary_vel;
}

/// Setup system: fill each `[[sph.insert]]` block with a lattice of particles.
pub fn sph_insert_atoms(
    domain: Res<Domain>,
    mut atom: ResMut<Atom>,
    registry: Res<AtomDataRegistry>,
    table: Res<SphMaterialTable>,
    cfg: Res<SphConfig>,
) {
    let inserts = match &cfg.insert {
        Some(v) if !v.is_empty() => v,
        _ => return,
    };
    let mut tag = atom.natoms() as u32;

    for ins in inserts {
        let mat_idx = table.index_of(&ins.material).unwrap_or_else(|| {
            eprintln!(
                "ERROR: unknown material '{}' in [[sph.insert]]. Available: {:?}",
                ins.material, table.names
            );
            std::process::exit(1);
        });
        let params = &table.params[mat_idx];
        let rho0 = ins.rest_density.unwrap_or(params.rho_c);
        let p0 = pressure(rho0, params);
        let h = ins.h_factor.unwrap_or(1.3) * ins.spacing;
        let mass = rho0 * ins.spacing.powi(3);
        let vel = ins.velocity.unwrap_or([0.0, 0.0, 0.0]);
        let t0 = ins.initial_temperature.unwrap_or(0.0);
        let is_boundary = if ins.frozen.unwrap_or(false) {
            1.0
        } else {
            0.0
        };
        let bvel = ins.boundary_velocity.unwrap_or([0.0, 0.0, 0.0]);

        let xs = lattice_axis(ins.region_min[0], ins.region_max[0], ins.spacing);
        let ys = lattice_axis(ins.region_min[1], ins.region_max[1], ins.spacing);
        let zs = lattice_axis(ins.region_min[2], ins.region_max[2], ins.spacing);

        for &x in &xs {
            for &y in &ys {
                for &z in &zs {
                    let pos = [x, y, z];
                    if owns_position(&domain, &pos) {
                        tag += 1;
                        insert_particle(
                            &mut atom,
                            &registry,
                            pos,
                            vel,
                            mass,
                            rho0,
                            p0,
                            t0,
                            h,
                            is_boundary,
                            bvel,
                            mat_idx as u32,
                            tag,
                        );
                    }
                }
            }
        }
    }

    println!("dev_soil_sph: inserted {} local particles", atom.nlocal());
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lattice_axis_cell_centered() {
        let xs = lattice_axis(0.0, 1.0, 0.25);
        assert_eq!(xs.len(), 4);
        assert!((xs[0] - 0.125).abs() < 1e-12);
        assert!((xs[3] - 0.875).abs() < 1e-12);
    }

    #[test]
    fn lattice_axis_thin_slab_single_layer() {
        let ys = lattice_axis(0.0, 0.1, 0.25); // extent < spacing
        assert_eq!(ys.len(), 1);
        assert!((ys[0] - 0.05).abs() < 1e-12);
    }
}
