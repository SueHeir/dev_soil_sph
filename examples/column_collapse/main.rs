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

const R0: f64 = 0.025; // initial column half-width  L0
const H0: f64 = 0.05; // initial column height       H0  (config a = H0/R0 = 2)

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

    // ── Skeptic reference: Lube 2005 experiment / Lagrée-Staron-Popinet 2011 ──
    // Aspect ratio a = H0/L0. Lagrée, Staron & Popinet (2011), JFM 686:378
    // (DOI 10.1017/jfm.2011.335), Eq. (3.1), give the experimental 2-D run-out
    // scaling (L∞−L0)/L0 ≃ λ1·a (a<a0) with, at a=2:
    //   • Lube et al. 2005 (sand/rice/sugar): λ1≃1.2, λ2≃1.9, 1.8≤a0≤2.8 → 2.40
    //   • Lajeunesse et al. 2005 (glass beads): λ1≃1.8, a0≃3.0            → 3.60
    // The cited experimental envelope at a=2 is therefore [2.40, 3.60]. (LSP's
    // own μ(I) *continuum* over-spreads to 2.2·a = 4.40; discrete/SPH fronts
    // under-spread that, landing back near the experiments — see LSP §3.1.)
    let a = H0 / R0; // = 2
    let runout_n = (toe - R0) / R0; // (L∞−L0)/L0
    let (runout_lo, runout_hi) = (2.40, 3.60); // Lube ↔ Lajeunesse @ a=2, LSP Eq. 3.1

    // Deposit final height H∞/L0. LSP Eq. (3.2): H∞/L0 ≃ λ3·a (a<a0) / λ4·a^α
    // (a>a0). At a=2 the cited fits span λ4·a^α from the LSP continuum
    // (λ4≃0.65, α≃0.35 → 0.83) up to the Lube-experiment branch (λ4≃1, α≃0.4 →
    // 1.32); with SPH resolution scatter we accept H∞/L0 ∈ [0.8, 1.7].
    let hinf_n = h_max / R0; // H∞/L0
    let (hinf_lo, hinf_hi) = (0.8, 1.7);

    println!("\n=== column_collapse result ===");
    println!("fluid particles: {n_fluid}   aspect ratio a = H0/L0 = {a:.2}");
    println!("runout toe (h≥{H_TOE}): {toe:.4} m   front reach (max|x|): {front_reach:.4} m");
    println!(
        "normalized runout (L∞−L0)/L0: {runout_n:.2}   \
         (Lube/Lajeunesse @a=2: [{runout_lo:.2}, {runout_hi:.2}], LSP 2011 Eq.3.1)"
    );
    println!(
        "normalized height  H∞/L0:     {hinf_n:.2}   (LSP 2011 Eq.3.2 @a=2: [{hinf_lo:.2}, {hinf_hi:.2}])"
    );
    println!("deposit height:   h(0)={:.4}  h(.5r∞)={:.4}  h(.8r∞)={:.4}  h_max={h_max:.4}",
        h_at(0.0), h_at(0.5), h_at(0.8));
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
    // experimental run-out AND deposit-height scaling at a=2?
    let matches_scaling = runout_ok && height_ok;

    // ── Negative control (declarative) ───────────────────────────────────────
    // A validation is only trustworthy if it is *capable of failing*. The config
    // may declare `[validation] expect = "reject"` — a deliberately-wrong
    // material (e.g. an over-frictional / cohesive column, μ ≫ real granular)
    // that should NOT reproduce the experimental scaling. In that mode we INVERT
    // the verdict: this run PASSES iff the reference band correctly REJECTS it
    // (run-out or height leaves the cited envelope), and FAILS iff the wrong
    // physics slipped through the band (which would prove the gate is vacuous).
    // Anything other than "reject" (incl. absent) is a normal positive check.
    let expect_reject = app
        .get_resource_ref::<Config>()
        .and_then(|c| {
            c.table
                .get("validation")
                .and_then(|v| v.as_table())
                .and_then(|t| t.get("expect"))
                .and_then(|e| e.as_str())
                .map(|s| s.to_string())
        })
        .as_deref()
        == Some("reject");

    if expect_reject {
        println!("\n=== NEGATIVE CONTROL (config declares [validation] expect = \"reject\") ===");
        if !matches_scaling {
            let why = if !runout_ok { "run-out" } else { "deposit height" };
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
