//! DEMO (not a validation) — machinery showcase, excluded from the dev_sph
//! validation set (see validation/README.md). Load-bearing seed for the future
//! SPH-CFD plume-surface coupling; no independent oracle yet.
//!
//! Footpad (foundation) — a driven rigid plate penetrating a granular bed.
//!
//! Validates the moving-rigid-boundary + reaction-force machinery: the plate is
//! pinned to a prescribed downward velocity (a footpad penetrating), and we read
//! the net bed reaction on it via [`MudPlateForce`]. Logs the force–sinkage curve.
//! Cold bed for now; the de-fluidization-coupled time-dependent bearing (the real
//! landing observable) builds on this.
//!
//! Run:
//!   cargo run --release --example footpad -- examples/footpad/config.toml
//! Then open examples/footpad/dump/ in OVITO (color by speed / pressure).

use mud_core::prelude::*;

const PLATE_Z0: f64 = 0.0575; // initial mean height of the 3-layer plate

/// Force–sinkage log: (sinkage [m], vertical bed reaction [N]).
#[derive(Default)]
struct SinkageLog(Vec<(f64, f64)>);

/// Sample the plate reaction every 100 steps (runs after the freeze updates it).
fn record(plate: Res<MudPlateForce>, mut log: ResMut<SinkageLog>, mut step: Local<u64>) {
    *step += 1;
    if *step % 100 == 0 && plate.n > 0 {
        log.0.push((PLATE_Z0 - plate.z, plate.force[2]));
    }
}

fn main() {
    let mut app = App::new();
    app.add_plugins(CorePlugins)
        .add_plugins(MudDefaultPlugins)
        .add_plugins(MudGravityPlugin);
    app.add_resource(SinkageLog::default());
    // After the freeze (PostForce) has updated MudPlateForce.
    app.add_update_system(record.after("mud_freeze"), ParticleSimScheduleSet::PostForce);

    {
        let dump = app.get_resource_ref::<DumpRegistry>().expect("DumpRegistry");
        dump.register_scalar("pressure", |a, r| {
            r.expect::<MudAtom>("d").pressure[..a.nlocal as usize].to_vec()
        });
        dump.register_scalar("is_boundary", |a, r| {
            r.expect::<MudAtom>("d").is_boundary[..a.nlocal as usize].to_vec()
        });
        dump.register_scalar("speed", |a, _r| {
            (0..a.nlocal as usize)
                .map(|i| {
                    let v = a.vel[i];
                    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
                })
                .collect()
        });
    }

    app.start();

    let plate = app.get_resource_ref::<MudPlateForce>().expect("MudPlateForce");
    let log = app.get_resource_ref::<SinkageLog>().expect("SinkageLog");
    let sinkage = PLATE_Z0 - plate.z;
    let f_z = plate.force[2];

    // Write the force–sinkage curve.
    if let Some(dir) = app.get_resource_ref::<Input>().and_then(|i| i.output_dir.clone()) {
        let mut s = String::from("sinkage,force_z\n");
        for &(d, f) in &log.0 {
            s.push_str(&format!("{d:.5},{f:.4}\n"));
        }
        let _ = std::fs::write(format!("{dir}/sinkage.csv"), s);
        println!("(force–sinkage curve written to {dir}/sinkage.csv)");
    }

    // Did the reaction grow with depth? Compare early vs late samples.
    let early = log.0.iter().take(5).map(|&(_, f)| f).sum::<f64>() / 5.0_f64.max(1.0);
    let late = log.0.iter().rev().take(5).map(|&(_, f)| f).sum::<f64>() / 5.0_f64.max(1.0);

    println!("\n=== footpad result ===");
    println!("plate particles: {}", plate.n);
    println!("final sinkage:   {sinkage:.4} m");
    println!("bed reaction Fz: {f_z:.3} N  (early ~{early:.3}, late ~{late:.3})");
    println!("OVITO:           examples/footpad/dump/*.lammpstrj");

    // Mechanism check: the footpad penetrated, the bed pushed back upward (Fz > 0),
    // the reaction grew as it sank, and nothing blew up.
    let penetrated = sinkage > 0.004;
    let resisting = f_z > 0.0;
    let grew = late > early;
    let bounded = f_z.is_finite() && f_z < 1.0e6;
    if penetrated && resisting && grew && bounded {
        println!("PASS: footpad penetrated; bed reaction grows with sinkage");
    } else {
        eprintln!("FAIL: penetrated={penetrated} (sink {sinkage:.4}), resisting={resisting}, grew={grew} ({early:.3}→{late:.3}), bounded={bounded}");
        std::process::exit(1);
    }
}
