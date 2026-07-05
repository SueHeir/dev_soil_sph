//! DEMO (not a validation) — self-consistent operator unit check, excluded from
//! the dev_sph validation set (see validation/README.md). It checks the SPH
//! Laplacian against its own analytic form, not an independent reference.
//!
//! Conduction operator test — checks the SPH Laplacian for `∇·(κ∇T)`.
//!
//! Impose `T(x) = T0(1 + A sin(2πx/Lx))`, run one step, and check the gathered
//! conduction term `lap_t` matches the analytical `∇·(κ∇T) = κ ∇²T = −κ k² (T−T0)`
//! (for small A so κ ≈ const). Least-squares slope of `lap_t` vs `(T−T0)` should
//! equal `−κ(T0) k²`. Granular conduction is collisional (short range), so this
//! isolates the operator rather than a hot/cold relaxation (where dissipation
//! dominates at SPH scales).
//!
//! Run:
//!   cargo run --release --example conduction_test -- examples/conduction_test/config.toml

use sph_core::prelude::*;
use std::f64::consts::PI;

const L_X: f64 = 0.05;
const T0: f64 = 0.01;
const AMP: f64 = 0.1;

/// Setup system: overwrite T with the sinusoidal field after insertion.
fn init_sinusoid(atoms: Res<Atom>, registry: Res<AtomDataRegistry>) {
    let mut sph = registry.expect_mut::<SphAtom>("init_sinusoid");
    let k = 2.0 * PI / L_X;
    for i in 0..atoms.nlocal as usize {
        sph.temperature[i] = T0 * (1.0 + AMP * (k * atoms.pos[i][0]).sin());
    }
}

fn main() {
    let mut app = App::new();
    app.add_plugins(CorePlugins).add_plugins(SphDefaultPlugins);
    app.add_setup_system(init_sinusoid, ScheduleSetupSet::PostSetup);
    app.start();

    let atoms = app.get_resource_ref::<Atom>().expect("Atom");
    let registry = app.get_resource_ref::<AtomDataRegistry>().expect("registry");
    let sph = registry.expect::<SphAtom>("conduction post-check");
    let table = app.get_resource_ref::<SphMaterialTable>().expect("materials");
    let n = atoms.nlocal as usize;
    let mat = &table.params[0];

    // Least-squares slope of lap_t vs (T − T0), through the origin.
    let mut sxy = 0.0;
    let mut sxx = 0.0;
    let mut sum_rho = 0.0;
    for i in 0..n {
        let dt = sph.temperature[i] - T0;
        sxy += sph.lap_t[i] * dt;
        sxx += dt * dt;
        sum_rho += sph.density[i];
    }
    let slope = sxy / sxx; // measured ∇·(κ∇T)/(T−T0) ≈ −κ k²
    let mean_rho = sum_rho / n as f64;

    let k = 2.0 * PI / L_X;
    let kappa = kt_conductivity(mean_rho, T0, mat);
    let predicted = -kappa * k * k;
    let rel_err = (slope - predicted).abs() / predicted.abs();

    println!("\n=== conduction_test result ===");
    println!("particles: {n}   κ(T0) = {kappa:.4e} Pa·s   k = {k:.1} 1/m");
    println!("measured ∇·(κ∇T)/(T−T0) = {slope:.4e}");
    println!("analytical −κ k²        = {predicted:.4e}   (rel err {rel_err:.3})");

    if rel_err < 0.20 {
        println!("PASS: SPH conduction operator matches ∇·(κ∇T) = −κ k² (T−T0)");
    } else {
        eprintln!("FAIL: conduction operator rel err {rel_err:.3} (want <0.20)");
        std::process::exit(1);
    }
}
