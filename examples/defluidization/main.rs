//! De-fluidization — the landing-relevant transient (`physics-design.md` §11).
//!
//! A hot, dilated bed (Φ < Φ_c, so the contact branch is off) is held up purely by
//! the collisional pressure `p_KT(T)`. As `T` decays by dissipation, `p_KT` fades,
//! gravity consolidates the bed (ρ rises past ρ_c), and the **contact branch takes
//! over** — the stress hands off from collisional to enduring-contact, exactly the
//! plume-fluidized → load-bearing transition under a lander. Uses only the
//! two-branch stress + granular temperature + continuity already in the model.
//!
//! Run:
//!   cargo run --release --example defluidization -- examples/defluidization/config.toml
//! Then open examples/defluidization/dump/ in OVITO (color by temperature / density).

use mud_core::prelude::*;

const T0: f64 = 0.05; // initial granular temperature
const RHO0: f64 = 1375.0; // initial bed density (Φ = 0.55)
const RHO_C: f64 = 1500.0; // contact-branch onset

fn main() {
    let mut app = App::new();
    app.add_plugins(CorePlugins)
        .add_plugins(MudDefaultPlugins)
        .add_plugins(MudGravityPlugin);

    {
        let dump = app.get_resource_ref::<DumpRegistry>().expect("DumpRegistry");
        dump.register_scalar("density", |a, r| {
            r.expect::<MudAtom>("d").density[..a.nlocal as usize].to_vec()
        });
        dump.register_scalar("temperature", |a, r| {
            r.expect::<MudAtom>("d").temperature[..a.nlocal as usize].to_vec()
        });
        dump.register_scalar("pressure", |a, r| {
            r.expect::<MudAtom>("d").pressure[..a.nlocal as usize].to_vec()
        });
        dump.register_scalar("is_boundary", |a, r| {
            r.expect::<MudAtom>("d").is_boundary[..a.nlocal as usize].to_vec()
        });
    }

    app.start();

    let atoms = app.get_resource_ref::<Atom>().expect("Atom");
    let registry = app.get_resource_ref::<AtomDataRegistry>().expect("registry");
    let sph = registry.expect::<MudAtom>("defluidization");
    let n = atoms.nlocal as usize;

    let mut sum_t = 0.0;
    let mut sum_rho = 0.0;
    let mut max_rho = 0.0f64;
    let mut bed_height = 0.0f64;
    let mut max_speed = 0.0f64;
    let mut n_connected = 0; // ρ ≥ ρ_c → contact branch active
    let mut n_fluid = 0;
    for i in 0..n {
        if sph.is_boundary[i] > 0.5 {
            continue;
        }
        n_fluid += 1;
        sum_t += sph.temperature[i];
        sum_rho += sph.density[i];
        max_rho = max_rho.max(sph.density[i]);
        bed_height = bed_height.max(atoms.pos[i][2]);
        if sph.density[i] >= RHO_C {
            n_connected += 1;
        }
        let v = atoms.vel[i];
        max_speed = max_speed.max((v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt());
    }
    let mean_t = sum_t / n_fluid as f64;
    let mean_rho = sum_rho / n_fluid as f64;
    let connected = n_connected as f64 / n_fluid as f64;

    println!("\n=== defluidization result ===");
    println!("fluid particles: {n_fluid}");
    println!("granular T:   {T0:.4} → {mean_t:.4e} m²/s²   ({:.1}% of initial)", 100.0 * mean_t / T0);
    println!("mean density: {RHO0:.0} → {mean_rho:.1} kg/m³   (max ρ = {max_rho:.1})");
    println!("contact-supported (ρ≥ρ_c): {:.0}% of bed   (was 0% — KT-supported)", 100.0 * connected);
    println!("bed height:   0.050 → {bed_height:.4} m   max speed {max_speed:.3e} m/s");
    println!("OVITO:        examples/defluidization/dump/*.lammpstrj (color by temperature)");

    // The de-fluidization signature: cooled, consolidated past ρ_c (stress
    // handed off to the contact branch), and at rest.
    let cooled = mean_t < 0.3 * T0;
    let consolidated = mean_rho > RHO0 && max_rho >= RHO_C;
    let bounded = max_speed < 2.0;
    if cooled && consolidated && bounded {
        println!("PASS: bed cooled and consolidated — stress handed off KT → contact");
    } else {
        eprintln!("FAIL: cooled={cooled} (T {mean_t:.2e}), consolidated={consolidated} (ρ̄ {mean_rho:.0}, ρmax {max_rho:.0}), bounded={bounded} ({max_speed:.2e})");
        std::process::exit(1);
    }
}
