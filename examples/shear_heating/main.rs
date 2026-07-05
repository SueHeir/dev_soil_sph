//! DEMO (not a validation) — self-consistent showcase, excluded from the dev_sph
//! validation set (see validation/README.md). Its PASS/FAIL is a smoke check
//! against dev_soil_sph itself, not an independent reference.
//!
//! Shear heating — end-to-end demo of the KT shear-production term.
//!
//! A homogeneous dilute granular gas (pure KT branch) is sheared by Lees-Edwards
//! simple shear (`soil_deform`'s xy deform — the same rig as the DEM LEBC). The
//! imposed shear `γ̇` produces granular temperature via `(4η_KT/3ρ)(D':D')`, which
//! heats `T` from cold until it balances dissipation at the Bagnold steady value
//! `T ∝ γ̇²`. Validates production in a real SPH shear flow.
//!
//! Run:
//!   cargo run --release --example shear_heating -- examples/shear_heating/config.toml

use sph_core::prelude::*;
use soil_deform::DeformPlugin;

const GAMMA_DOT: f64 = 50.0; // must match [deform] xy rate

fn main() {
    let mut app = App::new();
    app.add_plugins(CorePlugins)
        .add_plugins(SphDefaultPlugins)
        .add_plugins(DeformPlugin);
    app.start();

    let atoms = app.get_resource_ref::<Atom>().expect("Atom");
    let registry = app.get_resource_ref::<AtomDataRegistry>().expect("registry");
    let sph = registry.expect::<SphAtom>("shear post-check");
    let table = app.get_resource_ref::<SphMaterialTable>().expect("materials");
    let n = atoms.nlocal as usize;
    let mat = &table.params[0];

    // Measured: mean density, mean T, mean velocity gradient L_xy (the imposed shear).
    let mut sum_rho = 0.0;
    let mut sum_t = 0.0;
    let mut sum_lxy = 0.0;
    for i in 0..n {
        sum_rho += sph.density[i];
        sum_t += sph.temperature[i];
        sum_lxy += sph.velgrad[i][1]; // L_xy (row-major index 0*3+1)
    }
    let mean_rho = sum_rho / n as f64;
    let mean_t = sum_t / n as f64;
    let mean_lxy = sum_lxy / n as f64;

    // Analytical Bagnold steady T: solve production(T) + cooling(T) = 0 by bisection,
    // using the *imposed* shear (L_xy = γ̇) and the measured mean density.
    let l = [0.0, GAMMA_DOT, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    let balance = |t: f64| kt_production_rate(mean_rho, t, &l, mat) + kt_cooling_rate(mean_rho, t, mat);
    let (mut lo, mut hi) = (1.0e-12, 1.0);
    for _ in 0..100 {
        let mid = 0.5 * (lo + hi);
        if balance(mid) > 0.0 {
            lo = mid; // production still wins → hotter steady state
        } else {
            hi = mid;
        }
    }
    let t_steady = 0.5 * (lo + hi);
    let rel_err = (mean_t - t_steady).abs() / t_steady;

    println!("\n=== shear_heating result ===");
    println!("particles: {n}   Φ = {:.3}", mean_rho / mat.rho_s);
    println!("imposed γ̇ = {GAMMA_DOT:.1} 1/s   measured mean L_xy = {mean_lxy:.2} 1/s");
    println!("mean T:        {mean_t:.4e} m²/s²");
    println!("Bagnold T_ss:  {t_steady:.4e} m²/s²   (rel err {rel_err:.3})");

    // The shear must actually be imposed (L_xy ≈ γ̇), and T must reach the Bagnold
    // steady value. Tolerances are generous (v0: SPH gradient noise + the dilute
    // free-streaming regime).
    let shear_ok = (mean_lxy / GAMMA_DOT - 1.0).abs() < 0.25;
    let temp_ok = rel_err < 0.20;
    if shear_ok && temp_ok {
        println!("PASS: shear heats T to the Bagnold steady state");
    } else {
        eprintln!("FAIL: shear_ok={shear_ok} (L_xy ratio {:.3}), temp_ok={temp_ok} (rel err {rel_err:.3})", mean_lxy / GAMMA_DOT);
        std::process::exit(1);
    }
}
