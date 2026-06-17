//! # mud_atom — per-particle SPH data and granular material table
//!
//! The MUD physics tier's per-particle extension on the SOIL substrate, mirroring
//! the role `dirt_atom` plays for DEM. Provides:
//!
//! - [`MudAtom`] — the per-particle SPH column (smoothing length, density,
//!   pressure, deviatoric stress, and the per-step accumulators), registered via
//!   SOIL's `AtomData` derive. The `#[forward]`/`#[zero]` attributes encode the
//!   communication contract (see `docs/architecture.md` §3).
//! - [`MudMaterialTable`] — named granular materials, each a
//!   [`mud_constitutive::MaterialParams`], loaded from `[[mud.materials]]`.
//! - [`MudAtomPlugin`] — registers the column + material table.
//!
//! Insertion (`[[particles.insert]]`) is handled by SOIL/`dirt_atom`-style code
//! and added in a later increment.
//!
//! ## TOML
//! ```toml
//! [[mud.materials]]
//! name = "glass"
//! mu_s = 0.38
//! mu_2 = 0.64
//! i0 = 0.28
//! rho_s = 2500.0
//! rho_c = 1500.0
//! bulk_modulus = 3.75e6
//! poisson = 0.3
//! d = 0.5e-3
//! ```

use grass_app::prelude::*;
use grass_scheduler::prelude::*;
use serde::Deserialize;
use soil_derive::AtomData;

use soil_core::{register_atom_data, Atom, AtomData, AtomPlugin, Config, ScheduleSetupSet};

use mud_constitutive::MaterialParams;

pub mod insert;
pub use insert::*;

// ── Per-particle SPH column ──────────────────────────────────────────────────

/// Per-particle SPH extension data (`docs/architecture.md` §3).
///
/// Communication contract (the field attributes are load-bearing, not cosmetic):
/// - `#[forward]` — replicated owner→ghost each step; values a neighbor must read
///   (smoothing length, density, pressure, deviatoric stress, mass).
/// - `#[zero]` — reset each step; the per-step neighbor-sum accumulators.
/// - *(no attribute)* — migrates with the atom but never ghosts.
#[derive(AtomData)]
pub struct MudAtom {
    /// Smoothing length `h` (m). Neighbors need it to evaluate the kernel.
    #[forward]
    pub h: Vec<f64>,
    /// Density `ρ` (kg/m³) — persistent state (continuity-integrated), and read by
    /// neighbors in the momentum pass.
    #[forward]
    pub density: Vec<f64>,
    /// Pressure `p` (Pa) from the EOS. Read by neighbors in the momentum pass.
    #[forward]
    pub pressure: Vec<f64>,
    /// Deviatoric stress `s`, symmetric 3×3 as `[xx, yy, zz, xy, xz, yz]` (Pa).
    /// Persistent state (hypoelastic, evolved by the return map) and read by
    /// neighbors.
    #[forward]
    pub dev_stress: Vec<[f64; 6]>,
    /// Granular temperature `T = ⅓⟨δv²⟩` (m²/s²) — the collisional-branch state
    /// variable (`physics-design.md` §11). Persistent; read by neighbors (for
    /// their `σ_KT` and the future conduction Laplacian). Defaults to 0 → the
    /// collisional branch is off and behaviour reduces to the v0 contact model.
    #[forward]
    pub temperature: Vec<f64>,
    /// Velocity gradient `L = ∇v` (row-major `[f64; 9]`), accumulated in the
    /// density/velocity-gradient pass; reset each step.
    #[zero]
    pub velgrad: Vec<[f64; 9]>,
    /// `dρ/dt` from the continuity sum; reset each step.
    #[zero]
    pub drho_dt: Vec<f64>,
    /// Rest mass (kg). Neighbors need it for kernel sums; base `Atom.mass` is not
    /// forward-communicated, so it is carried here as a `#[forward]` column.
    #[forward]
    pub particle_mass: Vec<f64>,
    /// Boundary flag: 1.0 = frozen boundary particle (fixed in place, but still
    /// participates in the SPH sums to support the fluid), 0.0 = free. Migrates
    /// with the atom; not forwarded (only the freeze system reads it, on owners).
    pub is_boundary: Vec<f64>,
}

impl Default for MudAtom {
    fn default() -> Self {
        Self::new()
    }
}

impl MudAtom {
    /// An empty column with no particles; data is appended at insertion.
    pub fn new() -> Self {
        MudAtom {
            h: Vec::new(),
            density: Vec::new(),
            pressure: Vec::new(),
            dev_stress: Vec::new(),
            temperature: Vec::new(),
            velgrad: Vec::new(),
            drho_dt: Vec::new(),
            particle_mass: Vec::new(),
            is_boundary: Vec::new(),
        }
    }
}

// ── Material table ───────────────────────────────────────────────────────────

/// A single granular material from `[[mud.materials]]`.
#[derive(Deserialize, Clone)]
#[serde(deny_unknown_fields)]
pub struct MudMaterialConfig {
    /// Material name, referenced by particle insert blocks.
    pub name: String,
    /// Static friction coefficient μ_s.
    pub mu_s: f64,
    /// Limiting friction coefficient μ_2.
    pub mu_2: f64,
    /// Inertial-number scale I_0.
    pub i0: f64,
    /// Solid-grain density ρ_s (kg/m³).
    pub rho_s: f64,
    /// Critical (close-packed) density ρ_c (kg/m³).
    pub rho_c: f64,
    /// Effective bulk modulus K (Pa) — weakly compressible.
    pub bulk_modulus: f64,
    /// Poisson ratio ν (used to derive the shear modulus G from K).
    pub poisson: f64,
    /// Grain diameter d (m).
    pub d: f64,
    /// Coefficient of restitution e (kinetic-theory branch); default 0.7.
    #[serde(default = "default_restitution")]
    pub restitution: f64,
}

fn default_restitution() -> f64 {
    0.7
}

impl MudMaterialConfig {
    /// Convert to the constitutive [`MaterialParams`], deriving `G` from `K, ν`.
    pub fn to_params(&self) -> MaterialParams {
        MaterialParams {
            mu_s: self.mu_s,
            mu_2: self.mu_2,
            i0: self.i0,
            rho_s: self.rho_s,
            rho_c: self.rho_c,
            k_bulk: self.bulk_modulus,
            g_shear: MaterialParams::shear_from_bulk_poisson(self.bulk_modulus, self.poisson),
            d: self.d,
            restitution: self.restitution,
        }
    }
}

/// A single lattice-fill insertion block from `[[mud.insert]]`.
#[derive(Deserialize, Clone)]
#[serde(deny_unknown_fields)]
pub struct MudInsertConfig {
    /// Material name (must match a `[[mud.materials]]` name).
    pub material: String,
    /// Lower corner of the fill block (m).
    pub region_min: [f64; 3],
    /// Upper corner of the fill block (m).
    pub region_max: [f64; 3],
    /// Lattice spacing Δ (m); particle volume is Δ³.
    pub spacing: f64,
    /// Initial velocity (m/s); defaults to zero.
    #[serde(default)]
    pub velocity: Option<[f64; 3]>,
    /// Initial (rest) density (kg/m³); defaults to the material's ρ_c.
    #[serde(default)]
    pub rest_density: Option<f64>,
    /// Smoothing length as a multiple of spacing, `h = h_factor · Δ`; default 1.3.
    #[serde(default)]
    pub h_factor: Option<f64>,
    /// If true, these are frozen boundary particles (fixed in place, support the
    /// fluid). Default false.
    #[serde(default)]
    pub frozen: Option<bool>,
    /// Initial granular temperature T (m²/s²); default 0 (collisional branch off).
    #[serde(default)]
    pub initial_temperature: Option<f64>,
}

/// The `[mud]` config section.
#[derive(Deserialize, Clone, Default)]
#[serde(deny_unknown_fields)]
pub struct MudConfig {
    /// Granular material definitions.
    #[serde(default)]
    pub materials: Option<Vec<MudMaterialConfig>>,
    /// Particle insertion blocks.
    #[serde(default)]
    pub insert: Option<Vec<MudInsertConfig>>,
}

/// Named granular materials, each resolved to constitutive [`MaterialParams`].
///
/// Type index `i` (the SOIL `Atom::atom_type`) maps to `params[i]` / `names[i]`.
#[derive(Default)]
pub struct MudMaterialTable {
    /// Material names, in type-index order.
    pub names: Vec<String>,
    /// Constitutive parameters, in type-index order.
    pub params: Vec<MaterialParams>,
}

impl MudMaterialTable {
    /// An empty table.
    pub fn new() -> Self {
        Self::default()
    }

    /// Append a material; returns its type index.
    pub fn add(&mut self, name: &str, params: MaterialParams) -> usize {
        let idx = self.names.len();
        self.names.push(name.to_string());
        self.params.push(params);
        idx
    }

    /// Look up a material's type index by name.
    pub fn index_of(&self, name: &str) -> Option<usize> {
        self.names.iter().position(|n| n == name)
    }
}

// ── Plugin ───────────────────────────────────────────────────────────────────

/// Registers the [`MudAtom`] column and the [`MudMaterialTable`] from
/// `[[mud.materials]]`.
pub struct MudAtomPlugin;

impl Plugin for MudAtomPlugin {
    fn provides(&self) -> Vec<&str> {
        vec!["mud_particles"]
    }

    fn default_config(&self) -> Option<&str> {
        Some(
            r#"# Granular material definitions for MUD (SPH) particles
[[mud.materials]]
name = "glass"
mu_s = 0.38
mu_2 = 0.64
i0 = 0.28
rho_s = 2500.0
rho_c = 1500.0
bulk_modulus = 3.75e6
poisson = 0.3
d = 0.5e-3"#,
        )
    }

    fn build(&self, app: &mut App) {
        app.add_plugins(AtomPlugin);

        register_atom_data!(app, MudAtom::new());

        let cfg = Config::load::<MudConfig>(app, "mud");

        let mut table = MudMaterialTable::new();
        if let Some(ref materials) = cfg.materials {
            for mat in materials {
                table.add(&mat.name, mat.to_params());
            }
        }
        app.add_resource(table);

        app.add_setup_system(set_mud_ntypes, ScheduleSetupSet::Setup);
    }
}

/// Setup system: set `Atom::ntypes` from the number of registered materials.
fn set_mud_ntypes(mut atoms: ResMut<Atom>, table: Res<MudMaterialTable>) {
    if !table.names.is_empty() {
        atoms.ntypes = table.names.len();
    }
}

// ── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empty_column_is_empty() {
        let a = MudAtom::new();
        assert_eq!(a.h.len(), 0);
        assert_eq!(a.dev_stress.len(), 0);
    }

    #[test]
    fn material_config_to_params_derives_shear() {
        let cfg = MudMaterialConfig {
            name: "glass".into(),
            mu_s: 0.38,
            mu_2: 0.64,
            i0: 0.28,
            rho_s: 2500.0,
            rho_c: 1500.0,
            bulk_modulus: 3.75e6,
            poisson: 0.3,
            d: 0.5e-3,
            restitution: 0.7,
        };
        let p = cfg.to_params();
        assert_eq!(p.mu_s, 0.38);
        assert_eq!(p.k_bulk, 3.75e6);
        // G = 3(1−2ν)/(2(1+ν)) K = (1.2/2.6)·K
        let expect_g = (1.2 / 2.6) * 3.75e6;
        assert!((p.g_shear - expect_g).abs() < 1.0);
    }

    #[test]
    fn material_table_index_lookup() {
        let mut t = MudMaterialTable::new();
        let i = t.add("glass", MaterialParams::glass_beads_v0());
        assert_eq!(i, 0);
        assert_eq!(t.index_of("glass"), Some(0));
        assert_eq!(t.index_of("steel"), None);
    }
}
