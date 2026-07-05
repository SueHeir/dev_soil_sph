//! # sph_core — dev_soil_sph umbrella crate
//!
//! The batteries-included entry point for the dev_soil_sph SPH tier, mirroring
//! `dirt_core`. Re-exports the sub-crates and provides two plugin groups plus a
//! prelude so a simulation is a few lines:
//!
//! ```rust,ignore
//! use sph_core::prelude::*;
//!
//! fn main() {
//!     let mut app = App::new();
//!     app.add_plugins(CorePlugins)        // GRASS+SOIL: input, comm, domain, neighbor, run, print
//!        .add_plugins(SphDefaultPlugins); // SPH: atom column, insertion, Verlet, the two-pass physics
//!     app.start();
//! }
//! ```
//!
//! `CorePlugins` is a verbatim replica of the generic GRASS+SOIL infrastructure
//! group (it contains no physics), so dev_soil_sph does not depend on DIRT.

// ── Sub-crate re-exports ─────────────────────────────────────────────────────
pub use sph_atom;
pub use sph_constitutive;
pub use sph_kernel;
pub use sph_physics;
pub use soil_core;
pub use soil_print;
pub use soil_verlet;

use grass_app::prelude::*;

/// Core simulation infrastructure (input, communication, domain decomposition,
/// neighbor lists, groups, run loop, output). Contains no physics.
///
/// Velocity Verlet integration is **not** included here — `SphDefaultPlugins`
/// adds it.
pub struct CorePlugins;

impl PluginGroup for CorePlugins {
    fn build(self) -> PluginGroupBuilder {
        PluginGroupBuilder::start::<Self>()
            .add(|app: &mut App| {
                app.set_warning_fn(soil_core::verlet_schedule_warnings);
            })
            .add(soil_core::InputPlugin)
            .add(soil_core::CommunicationPlugin)
            .add(soil_core::DomainPlugin)
            .add(soil_core::NeighborPlugin)
            .add(soil_core::GroupPlugin)
            .add(soil_core::RunPlugin)
            .add(soil_print::PrintPlugin)
    }
}

/// The default dev_soil_sph SPH physics group, in registration order:
/// - [`SphAtomPlugin`](sph_atom::SphAtomPlugin) — the `SphAtom` column + material table
/// - [`SphAtomInsertPlugin`](sph_atom::SphAtomInsertPlugin) — lattice particle insertion
/// - [`VelocityVerletPlugin`](soil_verlet::VelocityVerletPlugin) — translational integration
/// - [`SphPhysicsPlugin`](sph_physics::SphPhysicsPlugin) — the two-pass SPH pipeline
pub struct SphDefaultPlugins;

impl PluginGroup for SphDefaultPlugins {
    fn build(self) -> PluginGroupBuilder {
        PluginGroupBuilder::start::<Self>()
            .add(sph_atom::SphAtomPlugin)
            .add(sph_atom::SphAtomInsertPlugin)
            .add(soil_verlet::VelocityVerletPlugin::new())
            .add(sph_physics::SphPhysicsPlugin)
    }
}

/// The dev_soil_sph prelude — `use sph_core::prelude::*;` for everything a simulation needs.
pub mod prelude {
    pub use crate::{CorePlugins, SphDefaultPlugins};

    // dev_soil_sph types
    pub use sph_atom::{
        SphAtom, SphAtomInsertPlugin, SphAtomPlugin, SphConfig, SphInsertConfig, SphMaterialConfig,
        SphMaterialTable,
    };
    pub use sph_constitutive::{
        kt_conductivity, kt_cooling_rate, kt_pressure, kt_production_rate, kt_shear_viscosity,
        pair_correlation, pressure, two_branch_stress, update_stress, MaterialParams, StressOut,
        TwoBranchStress,
    };
    pub use sph_kernel::Kernel;
    pub use sph_physics::{SphGravity, SphGravityPlugin, SphPhysicsPlugin, SphPlateForce};

    // Derive macros (multi-stage runs, etc.)
    pub use grass_derive::{ScheduleSet, StageEnum};
    pub use soil_derive::AtomData;

    // Core framework / substrate (glob)
    pub use grass_app::prelude::*;
    pub use grass_scheduler::prelude::*;
    pub use soil_core::ParticleSimScheduleSet;
    pub use soil_core::*;
    pub use soil_print::*;
    pub use soil_verlet::*;
}
