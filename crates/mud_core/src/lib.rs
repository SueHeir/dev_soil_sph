//! # mud_core — MUD umbrella crate
//!
//! The batteries-included entry point for the MUD SPH tier, mirroring
//! `dirt_core`. Re-exports the sub-crates and provides two plugin groups plus a
//! prelude so a simulation is a few lines:
//!
//! ```rust,ignore
//! use mud_core::prelude::*;
//!
//! fn main() {
//!     let mut app = App::new();
//!     app.add_plugins(CorePlugins)        // GRASS+SOIL: input, comm, domain, neighbor, run, print
//!        .add_plugins(MudDefaultPlugins); // SPH: atom column, insertion, Verlet, the two-pass physics
//!     app.start();
//! }
//! ```
//!
//! `CorePlugins` is a verbatim replica of the generic GRASS+SOIL infrastructure
//! group (it contains no physics), so MUD does not depend on DIRT.

// ── Sub-crate re-exports ─────────────────────────────────────────────────────
pub use mud_atom;
pub use mud_constitutive;
pub use mud_kernel;
pub use mud_physics;
pub use soil_core;
pub use soil_print;
pub use soil_verlet;

use grass_app::prelude::*;

/// Core simulation infrastructure (input, communication, domain decomposition,
/// neighbor lists, groups, run loop, output). Contains no physics.
///
/// Velocity Verlet integration is **not** included here — `MudDefaultPlugins`
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

/// The default MUD SPH physics group, in registration order:
/// - [`MudAtomPlugin`](mud_atom::MudAtomPlugin) — the `MudAtom` column + material table
/// - [`MudAtomInsertPlugin`](mud_atom::MudAtomInsertPlugin) — lattice particle insertion
/// - [`VelocityVerletPlugin`](soil_verlet::VelocityVerletPlugin) — translational integration
/// - [`MudPhysicsPlugin`](mud_physics::MudPhysicsPlugin) — the two-pass SPH pipeline
pub struct MudDefaultPlugins;

impl PluginGroup for MudDefaultPlugins {
    fn build(self) -> PluginGroupBuilder {
        PluginGroupBuilder::start::<Self>()
            .add(mud_atom::MudAtomPlugin)
            .add(mud_atom::MudAtomInsertPlugin)
            .add(soil_verlet::VelocityVerletPlugin::new())
            .add(mud_physics::MudPhysicsPlugin)
    }
}

/// The MUD prelude — `use mud_core::prelude::*;` for everything a simulation needs.
pub mod prelude {
    pub use crate::{CorePlugins, MudDefaultPlugins};

    // MUD types
    pub use mud_atom::{
        MudAtom, MudAtomInsertPlugin, MudAtomPlugin, MudConfig, MudInsertConfig, MudMaterialConfig,
        MudMaterialTable,
    };
    pub use mud_constitutive::{
        kt_conductivity, kt_cooling_rate, kt_pressure, kt_production_rate, kt_shear_viscosity,
        pair_correlation, pressure, two_branch_stress, update_stress, MaterialParams, StressOut,
        TwoBranchStress,
    };
    pub use mud_kernel::Kernel;
    pub use mud_physics::{MudGravity, MudGravityPlugin, MudPhysicsPlugin};

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
