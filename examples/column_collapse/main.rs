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

    // Deposit profile on a RESOLUTION-INDEPENDENT grid: fixed bin width over a
    // fixed span, surface height (max z) per bin. This is the shape we compare.
    const DX_BIN: f64 = 0.005; // fixed physical bin width
    const SPAN: f64 = 0.14; // half-width of the binned region (covers the floor)
    const H_TOE: f64 = 0.005; // deposit-edge height threshold (robust runout)
    let nbins = (2.0 * SPAN / DX_BIN) as usize;

    let mut profile = vec![0.0f64; nbins];
    let mut front_reach = 0.0f64; // max |x| of any fluid particle (incl. stragglers)
    let mut max_speed = 0.0f64;
    let mut max_t = 0.0f64;
    let mut n_fluid = 0;
    for i in 0..n {
        if sph.is_boundary[i] > 0.5 {
            continue;
        }
        n_fluid += 1;
        let x = atoms.pos[i][0];
        front_reach = front_reach.max(x.abs());
        let b = (((x + SPAN) / DX_BIN) as usize).min(nbins - 1);
        profile[b] = profile[b].max(atoms.pos[i][2]);
        let v = atoms.vel[i];
        max_speed = max_speed.max((v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt());
        max_t = max_t.max(sph.temperature[i]);
    }
    let bin_x = |b: usize| -> f64 { -SPAN + (b as f64 + 0.5) * DX_BIN };
    let h_max = profile.iter().cloned().fold(0.0, f64::max);
    // Robust runout = furthest |x| where the deposit is still ≥ H_TOE thick.
    let toe = (0..nbins)
        .filter(|&b| profile[b] >= H_TOE)
        .map(|b| bin_x(b).abs())
        .fold(0.0, f64::max);
    let h_at = |frac: f64| -> f64 {
        let b = (((frac * toe + SPAN) / DX_BIN) as usize).min(nbins - 1);
        profile[b]
    };

    // Write the profile CSV (x, h) for overlay vs DEM / between resolutions.
    if let Some(dir) = app.get_resource_ref::<Input>().and_then(|i| i.output_dir.clone()) {
        let mut s = String::from("x,h\n");
        for b in 0..nbins {
            s.push_str(&format!("{:.5},{:.5}\n", bin_x(b), profile[b]));
        }
        let _ = std::fs::write(format!("{dir}/profile.csv"), s);
        println!("(deposit profile written to {dir}/profile.csv)");
    }

    // Optional DEM/reference overlay: set MUD_DEM_PROFILE=path.csv (x,h rows) to
    // get a one-number normalized L2 shape error vs that profile (e.g. the 100k-DEM
    // deposit). Binned onto the same grid as the SPH profile.
    if let Ok(path) = std::env::var("MUD_DEM_PROFILE") {
        if let Ok(txt) = std::fs::read_to_string(&path) {
            let mut refp = vec![0.0f64; nbins];
            for line in txt.lines().skip(1) {
                let mut it = line.split(',');
                if let (Some(xs), Some(hs)) = (it.next(), it.next()) {
                    if let (Ok(x), Ok(h)) = (xs.trim().parse::<f64>(), hs.trim().parse::<f64>()) {
                        let b = (((x + SPAN) / DX_BIN) as usize).min(nbins - 1);
                        refp[b] = refp[b].max(h);
                    }
                }
            }
            let (mut num, mut den) = (0.0, 0.0);
            for b in 0..nbins {
                num += (profile[b] - refp[b]).powi(2);
                den += refp[b].powi(2);
            }
            let l2 = (num / den.max(1e-12)).sqrt();
            println!("DEM overlay ({path}): normalized L2 shape error = {l2:.3}");
        }
    }

    // Lube planar scaling (rough): (r∞ − r0)/r0 ≈ 1.6 √a for a ≳ 1.7 (a = 2 here).
    let spread = (toe - R0) / R0;
    let lube = 1.6 * 2.0f64.sqrt();

    println!("\n=== column_collapse result ===");
    println!("fluid particles: {n_fluid}");
    println!("runout toe (h≥{H_TOE}): {toe:.4} m   front reach (max|x|): {front_reach:.4} m");
    println!("spread (toe−r0)/r0: {spread:.2}   (Lube ~{lube:.2} for a=2)");
    println!("deposit height:   h(0)={:.4}  h(.5r∞)={:.4}  h(.8r∞)={:.4}  h_max={h_max:.4}",
        h_at(0.0), h_at(0.5), h_at(0.8));
    println!("max fluid speed:  {max_speed:.3e} m/s");
    println!("max granular T:   {max_t:.3e} m²/s²  (seed was 1e-5 → grew where it sheared)");
    println!("OVITO frames:     <output>/dump/*.lammpstrj");

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
