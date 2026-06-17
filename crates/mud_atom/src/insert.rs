//! Particle insertion for MUD — regular-lattice block fill.
//!
//! SPH wants an *ordered* initial packing (not the random placement DEM uses), so
//! each `[[mud.insert]]` block fills its region on a regular cubic lattice of the
//! given spacing. Per-particle mass is `ρ₀ · Δ³` (the lattice cell volume) and the
//! neighbor `cutoff_radius` is set to the kernel support `2h`.

use grass_app::prelude::*;
use grass_scheduler::prelude::*;

use soil_core::{Atom, AtomDataRegistry, Domain, ScheduleSetupSet};

use mud_constitutive::pressure;

use crate::{MudAtom, MudConfig, MudMaterialTable};

/// Inserts MUD particles at setup time from `[[mud.insert]]` blocks.
pub struct MudAtomInsertPlugin;

impl Plugin for MudAtomInsertPlugin {
    fn build(&self, app: &mut App) {
        // Run after domain decomposition so each rank's subdomain bounds exist.
        app.add_setup_system(
            mud_insert_atoms.after("domain_read_input"),
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

/// Append one SPH particle to the shared `Atom` arrays and the `MudAtom` column,
/// keeping every column length-synchronized.
#[allow(clippy::too_many_arguments)]
fn insert_particle(
    atom: &mut Atom,
    sph: &mut MudAtom,
    pos: [f64; 3],
    vel: [f64; 3],
    mass: f64,
    density: f64,
    pressure0: f64,
    h: f64,
    is_boundary: f64,
    mat_idx: u32,
    tag: u32,
) {
    atom.natoms += 1;
    atom.nlocal += 1;
    atom.tag.push(tag);
    atom.origin_index.push(0);
    // Neighbor pair cutoff is (r_i + r_j)·skin; with uniform h this gives the
    // kernel support 2h when cutoff_radius = h.
    atom.cutoff_radius.push(h);
    atom.image.push([0, 0, 0]);
    atom.is_ghost.push(false);
    atom.pos.push(pos);
    atom.vel.push(vel);
    atom.force.push([0.0; 3]);
    atom.mass.push(mass);
    atom.inv_mass.push(1.0 / mass);
    atom.atom_type.push(mat_idx);

    sph.h.push(h);
    sph.density.push(density);
    sph.pressure.push(pressure0);
    sph.dev_stress.push([0.0; 6]);
    sph.velgrad.push([0.0; 9]);
    sph.drho_dt.push(0.0);
    sph.particle_mass.push(mass);
    sph.is_boundary.push(is_boundary);
}

/// Setup system: fill each `[[mud.insert]]` block with a lattice of particles.
pub fn mud_insert_atoms(
    domain: Res<Domain>,
    mut atom: ResMut<Atom>,
    registry: Res<AtomDataRegistry>,
    table: Res<MudMaterialTable>,
    cfg: Res<MudConfig>,
) {
    let inserts = match &cfg.insert {
        Some(v) if !v.is_empty() => v,
        _ => return,
    };
    let mut sph = registry.expect_mut::<MudAtom>("mud_insert_atoms");
    let mut tag = atom.natoms as u32;

    for ins in inserts {
        let mat_idx = table.index_of(&ins.material).unwrap_or_else(|| {
            eprintln!(
                "ERROR: unknown material '{}' in [[mud.insert]]. Available: {:?}",
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
        let is_boundary = if ins.frozen.unwrap_or(false) { 1.0 } else { 0.0 };

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
                            &mut atom, &mut sph, pos, vel, mass, rho0, p0, h, is_boundary,
                            mat_idx as u32, tag,
                        );
                    }
                }
            }
        }
    }

    println!("MUD: inserted {} local particles", atom.nlocal);
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
