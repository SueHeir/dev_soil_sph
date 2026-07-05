//! DEMO (not a validation) — self-consistent showcase, excluded from the dev_sph
//! validation set (see validation/README.md). Its PASS/FAIL is a smoke check
//! against dev_soil_sph itself, not an independent reference.
//!
//! Haff cooling — the first granular-temperature milestone (`physics-design.md` §11).
//!
//! A homogeneous dilute granular gas (Φ = 0.4, below ρ_c so only the collisional
//! KT branch is active) with no gravity and no shear. The granular temperature
//! `T` should decay by inelastic dissipation following Haff's law
//! `T(t) = T0 / (1 + t/τ)²`, while the bed stays at rest (uniform `p_KT` → no
//! pressure gradient → no force). Validates the `T` field, the dissipation term,
//! the integration, and that the collisional pressure does not spuriously move
//! the bed.
//!
//! Run:
//!   cargo run --release --example haff_cooling -- examples/haff_cooling/config.toml

use sph_core::prelude::*;

fn main() {
    let mut app = App::new();
    app.add_plugins(CorePlugins).add_plugins(SphDefaultPlugins);
    app.start();

    // ── Analyze ──────────────────────────────────────────────────────────────
    let atoms = app.get_resource_ref::<Atom>().expect("Atom");
    let registry = app
        .get_resource_ref::<AtomDataRegistry>()
        .expect("registry");
    let sph = registry.expect::<SphAtom>("haff post-check");
    let table = app
        .get_resource_ref::<SphMaterialTable>()
        .expect("material table");
    let run = app.get_resource_ref::<RunConfig>().expect("run config");
    let n = atoms.nlocal as usize;

    let mat = &table.params[0];
    let stage = run.current_stage(0);
    let t_end = stage.steps as f64 * stage.dt;

    let mut sum_rho = 0.0;
    let mut sum_t = 0.0;
    let mut max_speed = 0.0f64;
    for i in 0..n {
        sum_rho += sph.density[i];
        sum_t += sph.temperature[i];
        let v = atoms.vel[i];
        max_speed = max_speed.max((v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt());
    }
    let mean_rho = sum_rho / n as f64;
    let mean_t = sum_t / n as f64;

    // Analytical Haff law from the run's actual material parameters.
    const T0: f64 = 0.01; // must match config initial_temperature
    let phi = mean_rho / mat.rho_s;
    let g0 = pair_correlation(phi);
    let zeta = 12.0 / std::f64::consts::PI.sqrt();
    let a = 2.0 * zeta * phi * g0 * (1.0 - mat.restitution * mat.restitution) / (3.0 * mat.d);
    let tau = 2.0 / (a * T0.sqrt());
    let analytic = T0 / (1.0 + t_end / tau).powi(2);
    let rel_err = (mean_t - analytic).abs() / analytic;

    println!("\n=== haff_cooling result ===");
    println!("particles: {n}   Φ = {phi:.3}   τ = {tau:.3e} s   t_end = {t_end:.3e} s ({:.2} τ)", t_end / tau);
    println!("mean T:      {mean_t:.4e} m²/s²");
    println!("Haff T(t):   {analytic:.4e} m²/s²   (rel err {rel_err:.3})");
    println!("max speed:   {max_speed:.3e} m/s  (bed should stay at rest)");

    if rel_err < 0.05 && max_speed < 1.0e-3 {
        println!("PASS: granular temperature decays by Haff's law, bed at rest");
    } else {
        eprintln!("FAIL: rel_err={rel_err:.3} (want <0.05), max_speed={max_speed:.3e} (want <1e-3)");
        std::process::exit(1);
    }
}
