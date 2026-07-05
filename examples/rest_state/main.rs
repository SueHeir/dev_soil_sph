//! Rest-state smoke test — the minimal end-to-end dev_soil_sph run.
//!
//! A uniform lattice block in a fully periodic box with no gravity should stay at
//! rest: the density is uniform, so the pressure is uniform, so the symmetric SPH
//! stress-gradient force cancels on the lattice (Σⱼ ∇Wᵢⱼ = 0). If the particles
//! stay put, the whole pipeline is self-consistent — insertion, neighbor list
//! (with PBC ghosts), the density/velocity-gradient pass, the constitutive
//! update, the mid-step halo exchange, the momentum pass, and integration.
//!
//! Run:
//!   cargo run --release --example rest_state -- examples/rest_state/config.toml

use sph_core::prelude::*;

fn main() {
    let mut app = App::new();
    app.add_plugins(CorePlugins).add_plugins(SphDefaultPlugins);
    app.start();

    // ── Post-run check ───────────────────────────────────────────────────────
    let atoms = app.get_resource_ref::<Atom>().expect("Atom resource");
    let registry = app
        .get_resource_ref::<AtomDataRegistry>()
        .expect("AtomDataRegistry resource");
    let sph = registry.expect::<SphAtom>("rest_state post-check");
    let n = atoms.nlocal as usize;

    let mut max_speed = 0.0f64;
    for i in 0..n {
        let v = atoms.vel[i];
        max_speed = max_speed.max((v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt());
    }
    let mut rho_min = f64::INFINITY;
    let mut rho_max = f64::NEG_INFINITY;
    for i in 0..n {
        rho_min = rho_min.min(sph.density[i]);
        rho_max = rho_max.max(sph.density[i]);
    }

    println!("\n=== rest_state result ===");
    println!("particles:     {n}");
    println!("max speed:     {max_speed:.3e} m/s");
    println!("density range: [{rho_min:.4}, {rho_max:.4}] kg/m^3");

    // Sound speed ~50 m/s; a stable rest state keeps speeds many orders below it.
    if max_speed < 1.0e-2 {
        println!("PASS: rest state preserved");
    } else {
        eprintln!("FAIL: particles drifted (max speed {max_speed:.3e} m/s)");
        std::process::exit(1);
    }
}
