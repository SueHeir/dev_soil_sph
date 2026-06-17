//! Granular column collapse — the classic μ(I) validation, with an OVITO dump.
//!
//! A dense column slumps under gravity onto a frozen floor and spreads into a
//! deposit. The contact (μ(I)) branch drives the dense flow; a small seeded
//! granular temperature lets the collisional branch light up where it shears.
//! Writes OVITO-native `.lammpstrj` frames (load the directory in OVITO and color
//! by `temperature`, `speed`, `pressure`, or `is_boundary`).
//!
//! Run:
//!   cargo run --release --example column_collapse -- examples/column_collapse/config.toml
//! Then open examples/column_collapse/dump/ in OVITO.

use mud_core::prelude::*;

const R0: f64 = 0.025; // initial column half-width

fn main() {
    let mut app = App::new();
    app.add_plugins(CorePlugins)
        .add_plugins(MudDefaultPlugins)
        .add_plugins(MudGravityPlugin);

    // Register MUD fields for the dump (OVITO can color by any of these).
    {
        let dump = app
            .get_resource_ref::<DumpRegistry>()
            .expect("DumpRegistry (PrintPlugin)");
        dump.register_scalar("density", |a, r| {
            r.expect::<MudAtom>("dump").density[..a.nlocal as usize].to_vec()
        });
        dump.register_scalar("pressure", |a, r| {
            r.expect::<MudAtom>("dump").pressure[..a.nlocal as usize].to_vec()
        });
        dump.register_scalar("temperature", |a, r| {
            r.expect::<MudAtom>("dump").temperature[..a.nlocal as usize].to_vec()
        });
        dump.register_scalar("is_boundary", |a, r| {
            r.expect::<MudAtom>("dump").is_boundary[..a.nlocal as usize].to_vec()
        });
        dump.register_scalar("speed", |a, _r| {
            (0..a.nlocal as usize)
                .map(|i| {
                    let v = a.vel[i];
                    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
                })
                .collect()
        });
        dump.register_vector("velocity", |a, _r| a.vel[..a.nlocal as usize].to_vec());
    }

    app.start();

    // ── Analyze the deposit ──────────────────────────────────────────────────
    let atoms = app.get_resource_ref::<Atom>().expect("Atom");
    let registry = app.get_resource_ref::<AtomDataRegistry>().expect("registry");
    let sph = registry.expect::<MudAtom>("collapse post-check");
    let n = atoms.nlocal as usize;

    let mut runout = 0.0f64; // max |x| of fluid
    let mut max_speed = 0.0f64;
    let mut max_t = 0.0f64;
    let mut n_fluid = 0;
    for i in 0..n {
        if sph.is_boundary[i] > 0.5 {
            continue;
        }
        n_fluid += 1;
        runout = runout.max(atoms.pos[i][0].abs());
        let v = atoms.vel[i];
        max_speed = max_speed.max((v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt());
        max_t = max_t.max(sph.temperature[i]);
    }

    // Lube planar scaling (rough): (r∞ − r0)/r0 ≈ 1.6 √a for a ≳ 1.7 (a = 2 here).
    let spread = (runout - R0) / R0;
    let lube = 1.6 * 2.0f64.sqrt();

    println!("\n=== column_collapse result ===");
    println!("fluid particles: {n_fluid}");
    println!("runout (max |x|): {runout:.4} m   (initial r0 = {R0})");
    println!("spread (r∞−r0)/r0: {spread:.2}   (Lube ~{lube:.2} for a=2)");
    println!("max fluid speed:  {max_speed:.3e} m/s");
    println!("max granular T:   {max_t:.3e} m²/s²  (seed was 1e-5 → grew where it sheared)");
    println!("OVITO frames:     examples/column_collapse/dump/*.lammpstrj");

    // Acceptance (v0, lenient): it collapsed/spread, did not blow up, and the
    // seeded temperature stayed bounded (didn't run away).
    let collapsed = spread > 0.5;
    let bounded = max_speed < 5.0 && max_t < 1.0;
    if collapsed && bounded {
        println!("PASS: column collapsed and spread into a deposit");
    } else {
        eprintln!("FAIL: collapsed={collapsed} (spread {spread:.2}), bounded={bounded} (vmax {max_speed:.2e}, Tmax {max_t:.2e})");
        std::process::exit(1);
    }
}
