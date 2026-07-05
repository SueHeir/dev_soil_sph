//! Hydrostatic column — the first non-trivial physics test.
//!
//! A fluid column settles under gravity onto frozen boundary layers (periodic in
//! x,y). At equilibrium the pressure should follow the hydrostatic law
//! `p(z) ≈ ρ g (z_top − z)`, rising linearly with depth, and the column should be
//! at rest. This exercises gravity, the boundary-freeze, and the SPH pressure
//! response together.
//!
//! Run:
//!   cargo run --release --example hydrostatic_column -- examples/hydrostatic_column/config.toml

use sph_core::prelude::*;

fn main() {
    let mut app = App::new();
    app.add_plugins(CorePlugins)
        .add_plugins(SphDefaultPlugins)
        .add_plugins(SphGravityPlugin);
    app.start();

    // ── Analyze the settled fluid column ─────────────────────────────────────
    let atoms = app.get_resource_ref::<Atom>().expect("Atom");
    let registry = app
        .get_resource_ref::<AtomDataRegistry>()
        .expect("registry");
    let sph = registry.expect::<SphAtom>("hydrostatic post-check");
    let n = atoms.nlocal as usize;

    const G: f64 = 9.81;
    const RHO_REF: f64 = 1500.0;
    const GRADIENT_RATIO_REGRESSION_FLOOR: f64 = 0.8208;
    const GRADIENT_RATIO_TARGET: f64 = 1.0;

    // Collect fluid particles (non-boundary): (z, pressure, density, speed).
    let mut fluid: Vec<(f64, f64, f64, f64)> = Vec::new();
    let mut max_speed_fluid = 0.0f64;
    for i in 0..n {
        if sph.is_boundary[i] > 0.5 {
            continue;
        }
        let z = atoms.pos[i][2];
        let v = atoms.vel[i];
        let speed = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
        max_speed_fluid = max_speed_fluid.max(speed);
        fluid.push((z, sph.pressure[i], sph.density[i], speed));
    }
    assert!(!fluid.is_empty(), "no fluid particles");

    let z_top = fluid.iter().fold(f64::NEG_INFINITY, |a, &(z, ..)| a.max(z));
    let z_bot = fluid.iter().fold(f64::INFINITY, |a, &(z, ..)| a.min(z));

    // Bin by z into slabs; mean pressure/density per slab.
    let nbins = 10usize;
    let dz = (z_top - z_bot) / nbins as f64;
    let mut sum_p = vec![0.0f64; nbins];
    let mut sum_rho = vec![0.0f64; nbins];
    let mut cnt = vec![0usize; nbins];
    for &(z, p, rho, _) in &fluid {
        let mut b = ((z - z_bot) / dz) as usize;
        if b >= nbins {
            b = nbins - 1;
        }
        sum_p[b] += p;
        sum_rho[b] += rho;
        cnt[b] += 1;
    }

    println!("\n=== hydrostatic_column result ===");
    println!("fluid particles: {}", fluid.len());
    println!("z range: [{z_bot:.4}, {z_top:.4}]  max fluid speed: {max_speed_fluid:.3e} m/s");
    println!(
        "{:>8} {:>8} {:>12} {:>12} {:>10}",
        "z", "n", "p_mean", "p_hydro", "rho_mean"
    );
    let z_surface = z_top + 0.5 * dz; // free surface ≈ half a slab above the top centroid
    for b in (0..nbins).rev() {
        if cnt[b] == 0 {
            continue;
        }
        let zc = z_bot + (b as f64 + 0.5) * dz;
        let pm = sum_p[b] / cnt[b] as f64;
        let rm = sum_rho[b] / cnt[b] as f64;
        let p_hydro = RHO_REF * G * (z_surface - zc);
        println!("{zc:8.4} {:8} {pm:12.2} {p_hydro:12.2} {rm:10.3}", cnt[b]);
    }

    // Physics check: the pressure gradient. Least-squares slope of p vs z over all
    // fluid particles should match dp/dz = −ρg (free of binning/staircase artifacts).
    let m = fluid.len() as f64;
    let z_bar = fluid.iter().map(|&(z, ..)| z).sum::<f64>() / m;
    let p_bar = fluid.iter().map(|&(_, p, ..)| p).sum::<f64>() / m;
    let mut sxy = 0.0;
    let mut sxx = 0.0;
    for &(z, p, ..) in &fluid {
        sxy += (z - z_bar) * (p - p_bar);
        sxx += (z - z_bar) * (z - z_bar);
    }
    let slope = sxy / sxx; // dp/dz [Pa/m]
    let expected = -RHO_REF * G; // −ρg
    let slope_ratio = slope / expected;

    let settled = max_speed_fluid < 0.5; // << sound speed (50 m/s)

    // ── Skeptic check #2: NO SPH tensile instability (Bui et al. 2008) ────────
    // Bui, Fukagawa, Sako & Ohno (2008), IJNAMG 32:1537 (DOI 10.1002/nag.688):
    // "the so-called SPH tensile instability ... may result in unrealistic
    // fracture and particles clustering ... when the material is in tension."
    // A physical hydrostatic column is compressive everywhere (p ≥ 0); a run that
    // develops spurious *tension* (p < 0) is exhibiting the tensile instability
    // Bui's artificial-stress fix exists to suppress. We assert the column stays
    // compressive to a small tolerance (−0.5% of the base hydrostatic pressure).
    let (mut p_min, mut p_max) = (f64::INFINITY, f64::NEG_INFINITY);
    for &(_, p, ..) in &fluid {
        p_min = p_min.min(p);
        p_max = p_max.max(p);
    }
    let p_base_hydro = RHO_REF * G * (z_top - z_bot); // ρg·H over the fluid span
    let p_tol = 0.005 * p_base_hydro; // allow only a tiny negative excursion
    let no_tension = p_min >= -p_tol;

    // ── Skeptic check #3: NO particle clumping / voids (Bui et al. 2008) ──────
    // The other signature of the tensile instability is density clustering
    // ("particles forming clumps"). A settled hydrostatic column is only weakly
    // compressible (K = 3.75 MPa, ρg·H ≈ 100 Pa ⇒ Δρ/ρ ~ 1e-4), so the density
    // field must stay tight. A clumped/voided field would show a wide spread.
    let (mut rho_min, mut rho_max) = (f64::INFINITY, f64::NEG_INFINITY);
    for &(_, _, rho, _) in &fluid {
        rho_min = rho_min.min(rho);
        rho_max = rho_max.max(rho);
    }
    let rho_spread = (rho_max - rho_min) / RHO_REF;
    let no_clumping = rho_spread < 0.02; // < 2 % — no clustering/voids

    println!(
        "\nsettled: {settled} (max speed {max_speed_fluid:.3e})\n\
         dp/dz = {slope:.1} Pa/m   expected −ρg = {expected:.1} Pa/m   ratio = {slope_ratio:.4}\n\
         regression: dp/dz ratio >= {GRADIENT_RATIO_REGRESSION_FLOOR:.4} (target {GRADIENT_RATIO_TARGET:.1})\n\
         tensile stability (Bui 2008): p ∈ [{p_min:.2}, {p_max:.2}] Pa  (want p_min ≥ {:.2})\n\
         clumping (Bui 2008):          ρ spread = {:.3}%  (want < 2%)",
        -p_tol,
        rho_spread * 100.0
    );

    let grad_ok = (0.7..=1.3).contains(&slope_ratio);
    let grad_regression_ok = slope_ratio >= GRADIENT_RATIO_REGRESSION_FLOOR;
    if settled && grad_ok && grad_regression_ok && no_tension && no_clumping {
        println!(
            "PASS: hydrostatic gradient dp/dz ≈ −ρg, compressive throughout (no tensile\n\
             instability), density field tight (no clumping), gradient regression floor held\n\
             — Bui 2008 skeptic checks green"
        );
    } else {
        eprintln!(
            "FAIL: settled={settled}, grad_ok={grad_ok} (ratio {slope_ratio:.4}, want 0.7–1.3), \
             grad_regression_ok={grad_regression_ok} (want >= {GRADIENT_RATIO_REGRESSION_FLOOR:.4}), \
             no_tension={no_tension} (p_min {p_min:.2}), no_clumping={no_clumping} (spread {:.3}%)",
            rho_spread * 100.0
        );
        std::process::exit(1);
    }
}
