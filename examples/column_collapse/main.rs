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

use sph_core::prelude::*;

const DEFAULT_R0: f64 = 0.025; // initial column half-width L0 used by checked-in configs

fn value_f64(value: &toml::Value) -> Option<f64> {
    value
        .as_float()
        .or_else(|| value.as_integer().map(|i| i as f64))
}

fn table_f64(table: &toml::Table, key: &str) -> Option<f64> {
    table.get(key).and_then(value_f64)
}

fn array3(table: &toml::Table, key: &str) -> Option<[f64; 3]> {
    let arr = table.get(key)?.as_array()?;
    if arr.len() != 3 {
        return None;
    }
    Some([
        value_f64(&arr[0])?,
        value_f64(&arr[1])?,
        value_f64(&arr[2])?,
    ])
}

fn column_geometry(config: &Config) -> (f64, f64) {
    let Some(inserts) = config
        .table
        .get("sph")
        .and_then(|v| v.as_table())
        .and_then(|sph| sph.get("insert"))
        .and_then(|v| v.as_array())
    else {
        return (DEFAULT_R0, 2.0 * DEFAULT_R0);
    };

    for insert in inserts {
        let Some(t) = insert.as_table() else { continue };
        if t.get("frozen").and_then(|v| v.as_bool()) == Some(true) {
            continue;
        }
        let Some(min) = array3(t, "region_min") else {
            continue;
        };
        let Some(max) = array3(t, "region_max") else {
            continue;
        };
        let r0 = min[0].abs().max(max[0].abs());
        let h0 = max[2] - min[2];
        if r0 > 0.0 && h0 > 0.0 {
            return (r0, h0);
        }
    }
    (DEFAULT_R0, 2.0 * DEFAULT_R0)
}

fn domain_span(config: &Config) -> f64 {
    config
        .table
        .get("domain")
        .and_then(|v| v.as_table())
        .map(|domain| {
            let x_low = table_f64(domain, "x_low").unwrap_or(-0.15);
            let x_high = table_f64(domain, "x_high").unwrap_or(0.15);
            x_low.abs().max(x_high.abs())
        })
        .unwrap_or(0.15)
}

fn expect_reject(config: &Config) -> bool {
    config
        .table
        .get("validation")
        .and_then(|v| v.as_table())
        .and_then(|t| t.get("expect"))
        .and_then(|e| e.as_str())
        == Some("reject")
}

fn runout_band(a: f64) -> (f64, f64) {
    if (a - 2.0).abs() < 1.0e-9 {
        return (2.40, 3.60);
    }
    if a < 2.0 {
        (1.2 * a, 2.2 * a)
    } else {
        // The experimental high-a envelope is the Lube/Lajeunesse range.
        // LSP's 2.2*a continuum curve is reported separately below; it is not
        // an experimental tolerance that can turn an experimental miss green.
        (1.9 * a.powf(2.0 / 3.0), 2.3 * a.powf(2.0 / 3.0))
    }
}

fn height_band(a: f64) -> (f64, f64) {
    if (a - 2.0).abs() < 1.0e-9 {
        return (0.80, 1.70);
    }
    if a < 2.0 {
        (0.75 * a, 1.25 * a)
    } else {
        (0.65 * a.powf(0.35), 1.30 * a.powf(0.40))
    }
}

fn main() {
    let mut app = App::new();
    app.add_plugins(CorePlugins)
        .add_plugins(SphDefaultPlugins)
        .add_plugins(SphGravityPlugin);

    // Register dev_soil_sph fields for the dump (OVITO can color by any of these).
    {
        let dump = app
            .get_resource_ref::<DumpRegistry>()
            .expect("DumpRegistry (PrintPlugin)");
        dump.register_scalar("density", |a, r| {
            r.expect::<SphAtom>("dump").density[..a.nlocal as usize].to_vec()
        });
        dump.register_scalar("pressure", |a, r| {
            r.expect::<SphAtom>("dump").pressure[..a.nlocal as usize].to_vec()
        });
        dump.register_scalar("temperature", |a, r| {
            r.expect::<SphAtom>("dump").temperature[..a.nlocal as usize].to_vec()
        });
        dump.register_scalar("is_boundary", |a, r| {
            r.expect::<SphAtom>("dump").is_boundary[..a.nlocal as usize].to_vec()
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
    let registry = app
        .get_resource_ref::<AtomDataRegistry>()
        .expect("registry");
    let sph = registry.expect::<SphAtom>("collapse post-check");
    let n = atoms.nlocal as usize;

    // Deposit profile on a RESOLUTION-INDEPENDENT grid: fixed bin width over the
    // declared domain span, surface height (max z) per bin.
    const DX_BIN: f64 = 0.005; // fixed physical bin width
    const H_TOE: f64 = 0.005; // deposit-edge height threshold (robust runout)
    let config = app.get_resource_ref::<Config>().expect("Config");
    let (r0, h0) = column_geometry(&config);
    let a = h0 / r0;
    let span = domain_span(&config) - DX_BIN; // keep profile bins away from fixed x walls
    let nbins = (2.0 * span / DX_BIN) as usize;

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
        let b = (((x + span) / DX_BIN) as usize).min(nbins - 1);
        profile[b] = profile[b].max(atoms.pos[i][2]);
        let v = atoms.vel[i];
        max_speed = max_speed.max((v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt());
        max_t = max_t.max(sph.temperature[i]);
    }
    let bin_x = |b: usize| -> f64 { -span + (b as f64 + 0.5) * DX_BIN };
    let h_max = profile.iter().cloned().fold(0.0, f64::max);
    // Robust runout = furthest |x| where the deposit is still ≥ H_TOE thick.
    let toe = (0..nbins)
        .filter(|&b| profile[b] >= H_TOE)
        .map(|b| bin_x(b).abs())
        .fold(0.0, f64::max);
    let h_at = |frac: f64| -> f64 {
        let b = (((frac * toe + span) / DX_BIN) as usize).min(nbins - 1);
        profile[b]
    };

    // Write the profile CSV (x, h) for overlay vs DEM / between resolutions.
    if let Some(dir) = app
        .get_resource_ref::<Input>()
        .and_then(|i| i.output_dir.clone())
    {
        let mut s = String::from("x,h\n");
        for b in 0..nbins {
            s.push_str(&format!("{:.5},{:.5}\n", bin_x(b), profile[b]));
        }
        let _ = std::fs::write(format!("{dir}/profile.csv"), s);
        println!("(deposit profile written to {dir}/profile.csv)");
    }

    // Optional DEM/reference overlay: set SPH_DEM_PROFILE=path.csv (x,h rows) to
    // get a one-number normalized L2 shape error vs that profile (e.g. the 100k-DEM
    // deposit). Binned onto the same grid as the SPH profile.
    if let Ok(path) = std::env::var("SPH_DEM_PROFILE") {
        if let Ok(txt) = std::fs::read_to_string(&path) {
            let mut refp = vec![0.0f64; nbins];
            for line in txt.lines().skip(1) {
                let mut it = line.split(',');
                if let (Some(xs), Some(hs)) = (it.next(), it.next()) {
                    if let (Ok(x), Ok(h)) = (xs.trim().parse::<f64>(), hs.trim().parse::<f64>()) {
                        let b = (((x + span) / DX_BIN) as usize).min(nbins - 1);
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

    // ── Skeptic reference: Lube 2005 experiment / Lagrée-Staron-Popinet 2011 ──
    // Aspect ratio a = H0/L0. Lagrée, Staron & Popinet (2011), JFM 686:378
    // (DOI 10.1017/jfm.2011.335), Eqs. (3.1)-(3.2), consolidate the Lube and
    // Lajeunesse experimental scalings plus the LSP continuum fit. Low a uses
    // λ1*a; high a uses λ2*a^(2/3), while the LSP continuum gives 2.2*a up to
    // a≈7. Height uses λ4*a^α. The a=2 gate keeps its original experimental
    // envelope exactly: run-out [2.40,3.60], height [0.80,1.70].
    let runout_n = (toe - r0) / r0; // (L∞−L0)/L0
    let (runout_lo, runout_hi) = runout_band(a);
    let hinf_n = h_max / r0; // H∞/L0
    let (hinf_lo, hinf_hi) = height_band(a);

    println!("\n=== column_collapse result ===");
    println!("fluid particles: {n_fluid}   aspect ratio a = H0/L0 = {a:.2}");
    println!("runout toe (h≥{H_TOE}): {toe:.4} m   front reach (max|x|): {front_reach:.4} m");
    println!(
        "normalized runout (L∞−L0)/L0: {runout_n:.2}   \
         (Lube/Lajeunesse band at this a: [{runout_lo:.2}, {runout_hi:.2}], LSP 2011 Eq.3.1)"
    );
    println!(
        "normalized height  H∞/L0:     {hinf_n:.2}   (LSP 2011 Eq.3.2 band: [{hinf_lo:.2}, {hinf_hi:.2}])"
    );
    println!(
        "deposit height:   h(0)={:.4}  h(.5r∞)={:.4}  h(.8r∞)={:.4}  h_max={h_max:.4}",
        h_at(0.0),
        h_at(0.5),
        h_at(0.8)
    );
    println!("max fluid speed:  {max_speed:.3e} m/s");
    println!("max granular T:   {max_t:.3e} m²/s²  (seed was 1e-5 → grew where it sheared)");
    println!("OVITO frames:     <output>/dump/*.lammpstrj");

    // Numeric acceptance vs the cited references: run-out AND deposit height fall
    // in the experimental 2-D envelope, the flow arrested (did not blow up), and
    // the seeded granular temperature stayed bounded (no runaway).
    let runout_ok = (runout_lo..=runout_hi).contains(&runout_n);
    let height_ok = (hinf_lo..=hinf_hi).contains(&hinf_n);
    let arrested = max_speed < 1.0; // deposit at rest (<< sound speed ~50 m/s)
    let bounded = max_speed < 5.0 && max_t < 1.0;
    // The reference-band verdict: does this deposit reproduce the Lube/Lajeunesse
    // experimental run-out AND deposit-height scaling at this aspect ratio?
    let matches_scaling = runout_ok && height_ok;

    // ── Negative control (declarative) ───────────────────────────────────────
    // A validation is only trustworthy if it is *capable of failing*. The config
    // may declare `[validation] expect = "reject"` — a deliberately wrong
    // material (e.g. an over-frictional / cohesive column, μ ≫ real granular)
    // that should NOT reproduce the experimental scaling. In that mode we INVERT
    // the verdict: this run PASSES iff the reference band correctly REJECTS it
    // (run-out or height leaves the cited envelope), and FAILS iff the wrong
    // physics slipped through the band (which would prove the gate is vacuous).
    // Anything other than "reject" (incl. absent) is a normal positive check.
    let expect_reject = expect_reject(&config);
    if expect_reject {
        println!("\n=== NEGATIVE CONTROL (config declares [validation] expect = \"reject\") ===");
        if !matches_scaling {
            let why = if !runout_ok {
                "run-out"
            } else {
                "deposit height"
            };
            println!(
                "PASS: reference band correctly REJECTED the deliberately-wrong material \
                 ({why} out of envelope: runout {runout_n:.2} want [{runout_lo:.2},{runout_hi:.2}], \
                 height {hinf_n:.2} want [{hinf_lo:.2},{hinf_hi:.2}]).\n\
                 The Lube/LSP gate is falsifiable — it fails on non-granular physics."
            );
        } else {
            eprintln!(
                "FAIL: negative control was NOT rejected — wrong physics landed INSIDE the \
                 cited band (runout {runout_n:.2} in [{runout_lo:.2},{runout_hi:.2}], \
                 height {hinf_n:.2} in [{hinf_lo:.2},{hinf_hi:.2}]). The gate is vacuous."
            );
            std::process::exit(1);
        }
    } else if matches_scaling && arrested && bounded {
        println!(
            "PASS: run-out {runout_n:.2} and height {hinf_n:.2} match the Lube 2005 / LSP 2011\n\
             2-D column-collapse scalings; deposit arrested and bounded"
        );
    } else {
        eprintln!(
            "FAIL: runout_ok={runout_ok} ({runout_n:.2} want [{runout_lo:.2},{runout_hi:.2}]), \
             height_ok={height_ok} ({hinf_n:.2} want [{hinf_lo:.2},{hinf_hi:.2}]), \
             arrested={arrested} (vmax {max_speed:.2e}), bounded={bounded} (Tmax {max_t:.2e})"
        );
        std::process::exit(1);
    }
}
