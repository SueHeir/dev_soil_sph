//! Footpad (foundation) — a driven rigid plate penetrating a granular bed.
//!
//! This binary remains the direct mechanics run: it drives the frozen plate,
//! records the bed reaction through [`SphPlateForce`], writes `sinkage.csv`, and
//! exits non-zero if the moving-boundary machinery does not produce a bounded,
//! growing upward reaction. The external bearing/sinkage validation lives in
//! `sweep.py`, which compares the emitted force-sinkage curve against an
//! independent Bekker/DIRT-DEM reference band and runs a zero-gravity control
//! that must be rejected.
//!
//! Cold bed for now; the de-fluidization-coupled time-dependent bearing (the real
//! landing observable) builds on this.
//!
//! Run:
//!   cargo run --release --example footpad -- examples/footpad/config.toml
//!   $BENCH_PYTHON examples/footpad/sweep.py
//! Then open examples/footpad/dump/ in OVITO (color by speed / pressure).

use sph_core::prelude::*;

const PLATE_Z0: f64 = 0.0575; // initial mean height of the 3-layer plate

/// Force–sinkage log: (sinkage [m], vertical bed reaction [N]).
#[derive(Default)]
struct SinkageLog(Vec<(f64, f64)>);

/// Sample the plate reaction every 100 steps (runs after the freeze updates it).
fn record(plate: Res<SphPlateForce>, mut log: ResMut<SinkageLog>, mut step: Local<u64>) {
    *step += 1;
    if *step % 100 == 0 && plate.n > 0 {
        log.0.push((PLATE_Z0 - plate.z, plate.force[2]));
    }
}

fn main() {
    let mut app = App::new();
    app.add_plugins(CorePlugins)
        .add_plugins(SphDefaultPlugins)
        .add_plugins(SphGravityPlugin);
    app.add_resource(SinkageLog::default());
    // After the freeze (PostForce) has updated SphPlateForce.
    app.add_update_system(
        record.after("sph_freeze"),
        ParticleSimScheduleSet::PostForce,
    );

    {
        let dump = app
            .get_resource_ref::<DumpRegistry>()
            .expect("DumpRegistry");
        dump.register_scalar("pressure", |a, r| {
            r.expect::<SphAtom>("d").pressure[..a.nlocal as usize].to_vec()
        });
        dump.register_scalar("is_boundary", |a, r| {
            r.expect::<SphAtom>("d").is_boundary[..a.nlocal as usize].to_vec()
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

    let plate = app
        .get_resource_ref::<SphPlateForce>()
        .expect("SphPlateForce");
    let log = app.get_resource_ref::<SinkageLog>().expect("SinkageLog");
    let sinkage = PLATE_Z0 - plate.z;
    let f_z = plate.force[2];

    // Write the force–sinkage curve.
    if let Some(dir) = app
        .get_resource_ref::<Input>()
        .and_then(|i| i.output_dir.clone())
    {
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
