//! Unit tests for the granular constitutive update.
//!
//! Gate #1 (`docs/physics-design.md` §8.1, the Dunatunga–Kamrin App. B
//! verification): drive a single material point with a prescribed velocity-
//! gradient history `L(t)` (simple shear → isotropic compression → extension),
//! integrate the *same* constitutive ODE with an independent 4th-order
//! Runge–Kutta reference at very small `dt`, and assert the forward
//! [`update_stress`] scheme converges to the RK4 reference at first order in
//! `dt` (error halves when `dt` halves).
//!
//! Plus focused tests: separation (`ρ < ρ_c → σ = 0`), pure elastic step below
//! yield (stress equals the trial), and pure shear steady state
//! (`τ̄/p → μ(I)`).

use super::*;

// ---------------------------------------------------------------------------
// Continuous constitutive ODE — the reference law the forward scheme discretizes
// ---------------------------------------------------------------------------

/// The continuous rate `ṡ` of the deviatoric stress at fixed density (the dense
/// branch, `ρ ≥ ρ_c`, `p > 0`). This is the ODE that [`update_stress`] is a
/// first-order (return-map / backward-projected forward-Euler) discretization
/// of, and that the RK4 reference below integrates.
///
/// Rate form (§3.3, continuous limit):
/// - Elastic trial rate `ṡ_tr = 2G D' + sW − Ws`.
/// - Below yield (`τ̄ ≤ μ_s p`): `ṡ = ṡ_tr`.
/// - At/above yield: viscoplastic flow removes the excess shear stress along the
///   current deviator direction at rate `G γ̇ᵖ`, where `γ̇ᵖ` is fixed by the
///   inertial consistency condition `τ̄ = μ(I) p`, `I = γ̇ᵖ d √ρ_s / √p`.
///   Inverting μ(I) for `I` and writing `γ̇ᵖ = I √p /(d √ρ_s)` gives a closed
///   form (see below). The deviatoric flow direction is `s/τ̄` and
///   `ṡ_plastic = − G γ̇ᵖ s / τ̄`.
fn dev_stress_rate(s: &Sym3, l: &[f64; 9], rho: f64, params: &MaterialParams) -> Sym3 {
    let p = pressure(rho, params);
    let g = params.g_shear;
    let (d, w) = decompose_velocity_gradient(l);
    let d_dev = deviator(&d);
    let jt = jaumann_term(s, &w);

    // Elastic trial rate.
    let mut rate: Sym3 = [
        2.0 * g * d_dev[XX] + jt[XX],
        2.0 * g * d_dev[YY] + jt[YY],
        2.0 * g * d_dev[ZZ] + jt[ZZ],
        2.0 * g * d_dev[XY] + jt[XY],
        2.0 * g * d_dev[XZ] + jt[XZ],
        2.0 * g * d_dev[YZ] + jt[YZ],
    ];

    let tau = equiv_shear_stress(s);
    let s0 = params.mu_s * p;
    if tau <= s0 || tau == 0.0 {
        return rate; // elastic
    }

    // Viscoplastic: relax τ̄ toward the rate-dependent yield μ(I) p.
    // The current ratio μ_now = τ̄/p (≥ μ_s). Invert μ(I) = μ_s + (μ_2−μ_s)/(I0/I+1):
    //   for μ in (μ_s, μ_2),  I = I0 (μ−μ_s)/(μ_2−μ). For μ ≥ μ_2, I → ∞ (cap).
    let mu_now = (tau / p).min(params.mu_2 - 1e-12);
    let denom = params.mu_2 - mu_now;
    let i = if denom <= 0.0 {
        f64::INFINITY
    } else {
        params.i0 * (mu_now - params.mu_s) / denom
    };
    let gamma_dot_p = i * p.sqrt() / (params.d * params.rho_s.sqrt());

    // Plastic correction along the deviator direction s/τ̄, magnitude G γ̇ᵖ.
    let coef = g * gamma_dot_p / tau;
    rate[XX] -= coef * s[XX];
    rate[YY] -= coef * s[YY];
    rate[ZZ] -= coef * s[ZZ];
    rate[XY] -= coef * s[XY];
    rate[XZ] -= coef * s[XZ];
    rate[YZ] -= coef * s[YZ];
    rate
}

/// One RK4 step of the continuous ODE `ds/dt = dev_stress_rate(s, L, ρ)` at
/// fixed `L`, `ρ`.
fn rk4_step(s: &Sym3, l: &[f64; 9], rho: f64, dt: f64, params: &MaterialParams) -> Sym3 {
    let add = |a: &Sym3, b: &Sym3, c: f64| -> Sym3 {
        [
            a[0] + c * b[0],
            a[1] + c * b[1],
            a[2] + c * b[2],
            a[3] + c * b[3],
            a[4] + c * b[4],
            a[5] + c * b[5],
        ]
    };
    let k1 = dev_stress_rate(s, l, rho, params);
    let k2 = dev_stress_rate(&add(s, &k1, 0.5 * dt), l, rho, params);
    let k3 = dev_stress_rate(&add(s, &k2, 0.5 * dt), l, rho, params);
    let k4 = dev_stress_rate(&add(s, &k3, dt), l, rho, params);
    [
        s[0] + dt / 6.0 * (k1[0] + 2.0 * k2[0] + 2.0 * k3[0] + k4[0]),
        s[1] + dt / 6.0 * (k1[1] + 2.0 * k2[1] + 2.0 * k3[1] + k4[1]),
        s[2] + dt / 6.0 * (k1[2] + 2.0 * k2[2] + 2.0 * k3[2] + k4[2]),
        s[3] + dt / 6.0 * (k1[3] + 2.0 * k2[3] + 2.0 * k3[3] + k4[3]),
        s[4] + dt / 6.0 * (k1[4] + 2.0 * k2[4] + 2.0 * k3[4] + k4[4]),
        s[5] + dt / 6.0 * (k1[5] + 2.0 * k2[5] + 2.0 * k3[5] + k4[5]),
    ]
}

// ---------------------------------------------------------------------------
// Prescribed L(t) history: shear → compression → extension
// ---------------------------------------------------------------------------

/// Velocity gradient at time `t` and the density rate `ρ̇ = −ρ tr(D)`-consistent
/// evolution handled by the driver. Three phases, each `T_PHASE` long.
///
/// Phase 1 (shear): pure simple shear `L_xy = γ̇`, exercises the plastic branch.
/// Phase 2 (compression): isotropic compression `L_xx=L_yy=L_zz = −ε̇` raises ρ
/// well above ρ_c (more pressure, more yield headroom).
/// Phase 3 (extension): isotropic extension `+ε̇` lowers ρ back toward/below ρ_c,
/// crossing the separation threshold and exercising the elastic + disconnect.
const T_PHASE: f64 = 1.0e-3; // s
const SHEAR_RATE: f64 = 200.0; // 1/s  (strong → plastic flow)
const VOL_RATE: f64 = 80.0; // 1/s

fn l_history(t: f64) -> [f64; 9] {
    let mut l = [0.0; 9];
    if t < T_PHASE {
        // simple shear: v_x = γ̇ y → L_xy = γ̇
        l[1] = SHEAR_RATE;
    } else if t < 2.0 * T_PHASE {
        // isotropic compression: diag = −ε̇
        l[0] = -VOL_RATE;
        l[4] = -VOL_RATE;
        l[8] = -VOL_RATE;
    } else {
        // isotropic extension: diag = +ε̇
        l[0] = VOL_RATE;
        l[4] = VOL_RATE;
        l[8] = VOL_RATE;
    }
    l
}

/// Advance density by the continuity equation `ρ̇ = −ρ tr(D)` over one step
/// (forward Euler — shared by both integrators so it is not the variable under
/// test; the deviatoric stress is what we compare).
fn advance_density(rho: f64, l: &[f64; 9], dt: f64) -> f64 {
    let (d, _) = decompose_velocity_gradient(l);
    let tr_d = trace(&d);
    rho * (1.0 - dt * tr_d)
}

// ---------------------------------------------------------------------------
// Gate #1: convergence of the forward scheme to the RK4 reference
// ---------------------------------------------------------------------------

/// Run the forward `update_stress` scheme over the full L(t) history at step
/// `dt`, returning the final deviatoric stress.
fn run_forward(dt: f64, params: &MaterialParams) -> Sym3 {
    let t_end = 3.0 * T_PHASE;
    let n = (t_end / dt).round() as usize;
    let mut s = [0.0; 6];
    let mut rho = 1.55 * 1000.0; // 1550 kg/m³, just above ρ_c=1500
    let mut t = 0.0;
    for _ in 0..n {
        let l = l_history(t + 0.5 * dt); // midpoint sample of the history
        rho = advance_density(rho, &l, dt);
        let out = update_stress(&s, &l, rho, dt, params);
        s = out.dev_stress;
        t += dt;
    }
    s
}

/// Run the RK4 reference (very small internal substep) over the same history,
/// returning the final deviatoric stress. The density is advanced with the same
/// continuity update at the *fine* substep so it is effectively exact.
fn run_reference(dt: f64, params: &MaterialParams) -> Sym3 {
    let t_end = 3.0 * T_PHASE;
    let n = (t_end / dt).round() as usize;
    let mut s = [0.0; 6];
    let mut rho = 1.55 * 1000.0;
    let mut t = 0.0;
    for _ in 0..n {
        let l = l_history(t + 0.5 * dt);
        rho = advance_density(rho, &l, dt);
        // Reference must also honor separation (the ODE only holds for ρ≥ρ_c).
        let p = pressure(rho, params);
        if rho < params.rho_c || p <= 0.0 {
            s = [0.0; 6];
        } else {
            // If we just reconnected with τ̄ above yield, the continuous law
            // would instantly relax it; integrate with RK4 either way.
            s = rk4_step(&s, &l, rho, dt, params);
        }
        t += dt;
    }
    s
}

fn dist(a: &Sym3, b: &Sym3) -> f64 {
    let mut acc = 0.0;
    for i in 0..6 {
        acc += (a[i] - b[i]).powi(2);
    }
    acc.sqrt()
}

#[test]
fn gate1_first_order_convergence_to_rk4() {
    let params = MaterialParams::glass_beads_v0();

    // A very fine RK4 solution is the "truth".
    let dt_truth = 1.0e-9;
    let truth = run_reference(dt_truth, &params);

    // Forward scheme at a sequence of halving timesteps; measure the error to
    // the truth and verify it ~halves each time (first order).
    let dts = [4.0e-7, 2.0e-7, 1.0e-7, 0.5e-7];
    let mut errs = Vec::new();
    for &dt in &dts {
        let s = run_forward(dt, &params);
        errs.push(dist(&s, &truth));
    }

    eprintln!("gate1 forward-vs-RK4 errors:");
    for (dt, e) in dts.iter().zip(&errs) {
        eprintln!("  dt = {:.2e}  err = {:.6e}", dt, e);
    }

    // Convergence rate r = log2(err[k]/err[k+1]) should be ≈ 1 (first order).
    let mut rates = Vec::new();
    for k in 0..errs.len() - 1 {
        let r = (errs[k] / errs[k + 1]).log2();
        rates.push(r);
        eprintln!("  rate (dt {:.1e}→{:.1e}) = {:.3}", dts[k], dts[k + 1], r);
    }

    // Each refinement at least roughly halves the error...
    for k in 0..errs.len() - 1 {
        assert!(
            errs[k + 1] < errs[k],
            "error did not decrease on refinement: {:.3e} -> {:.3e}",
            errs[k],
            errs[k + 1]
        );
    }
    // ...and the observed order is first-order (allow a band around 1.0 since the
    // history mixes elastic/plastic/separation phases).
    let avg_rate: f64 = rates.iter().sum::<f64>() / rates.len() as f64;
    assert!(
        avg_rate > 0.7 && avg_rate < 1.4,
        "average convergence rate {:.3} not first-order (expected ≈1)",
        avg_rate
    );
    // Errors must actually be small at the finest dt (sanity on magnitude).
    let rel = errs[errs.len() - 1] / equiv_shear_stress(&truth).max(1.0);
    assert!(rel < 1e-2, "final relative error too large: {:.3e}", rel);
}

// ---------------------------------------------------------------------------
// Focused tests
// ---------------------------------------------------------------------------

#[test]
fn separation_below_rho_c_is_stress_free() {
    let params = MaterialParams::glass_beads_v0();
    // Some nonzero prior deviator and an arbitrary L: must be wiped to zero.
    let s_n: Sym3 = [1.0e4, -0.5e4, -0.5e4, 2.0e3, 0.0, 0.0];
    let l = [0.0, 50.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    let out = update_stress(&s_n, &l, params.rho_c - 1.0, 1e-6, &params);
    assert!(out.disconnected);
    assert_eq!(out.pressure, 0.0);
    assert_eq!(out.dev_stress, [0.0; 6]);
    assert_eq!(out.sigma, [0.0; 6]);
}

#[test]
fn pressure_eos_and_separation_threshold() {
    let params = MaterialParams::glass_beads_v0();
    assert_eq!(pressure(params.rho_c - 10.0, &params), 0.0);
    assert_eq!(pressure(params.rho_c, &params), 0.0); // exactly at ρ_c → p=0
    let p = pressure(1.01 * params.rho_c, &params);
    assert!((p - params.k_bulk * 0.01).abs() < 1e-6);
    // Sound speed of the chosen K must satisfy the documented 50 m/s.
    assert!((params.sound_speed() - 50.0).abs() < 1e-9);
}

#[test]
fn elastic_step_below_yield_equals_trial() {
    let params = MaterialParams::glass_beads_v0();
    // Pressure high (well-packed), shear small and short → stays below μ_s p.
    let rho = 1.6 * 1000.0;
    let p = pressure(rho, &params);
    let s_n = [0.0; 6];
    // Tiny shear rate and tiny dt so τ̄_tr ≪ μ_s p.
    let gdot = 1.0e-3;
    let dt = 1.0e-6;
    let l = [0.0, gdot, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    let out = update_stress(&s_n, &l, rho, dt, &params);
    assert!(!out.disconnected);

    // Trial deviator computed independently.
    let (d, w) = decompose_velocity_gradient(&l);
    let d_dev = deviator(&d);
    let jt = jaumann_term(&s_n, &w);
    let s_tr: Sym3 = [
        s_n[0] + dt * (2.0 * params.g_shear * d_dev[0] + jt[0]),
        s_n[1] + dt * (2.0 * params.g_shear * d_dev[1] + jt[1]),
        s_n[2] + dt * (2.0 * params.g_shear * d_dev[2] + jt[2]),
        s_n[3] + dt * (2.0 * params.g_shear * d_dev[3] + jt[3]),
        s_n[4] + dt * (2.0 * params.g_shear * d_dev[4] + jt[4]),
        s_n[5] + dt * (2.0 * params.g_shear * d_dev[5] + jt[5]),
    ];
    let tau_tr = equiv_shear_stress(&s_tr);
    assert!(tau_tr <= params.mu_s * p, "test setup must stay below yield");
    assert_eq!(out.gamma_dot_p, 0.0);
    assert!(dist(&out.dev_stress, &s_tr) < 1e-20);
    // σ = −p I + s.
    assert!((out.sigma[0] - (-p + s_tr[0])).abs() < 1e-9);
}

#[test]
fn pure_shear_steady_state_gives_mu_of_i() {
    // Drive simple shear at fixed ρ until the stress reaches steady state; then
    // τ̄/p must equal μ(I) for the imposed inertial number I = γ̇ d√ρ_s/√p.
    let params = MaterialParams::glass_beads_v0();
    let rho = 1.6 * 1000.0;
    let p = pressure(rho, &params);
    let gdot = 500.0; // strong shear → clearly plastic, finite I
    let dt = 2.0e-8;
    let l = [0.0, gdot, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];

    let mut s = [0.0; 6];
    // Integrate long enough (relative to G/τ time scale) to reach steady state.
    for _ in 0..2_000_000 {
        let out = update_stress(&s, &l, rho, dt, &params);
        s = out.dev_stress;
    }
    let tau = equiv_shear_stress(&s);
    let mu_measured = tau / p;

    // The plastic shear rate at steady state equals the imposed γ̇ (all imposed
    // strain rate is plastic when the elastic stress is no longer changing).
    let i_imposed = gdot * params.d * params.rho_s.sqrt() / p.sqrt();
    let mu_target = mu_of_i(i_imposed, &params);

    eprintln!(
        "steady μ = {:.5}, μ(I={:.4e}) = {:.5}",
        mu_measured, i_imposed, mu_target
    );
    assert!(
        (mu_measured - mu_target).abs() < 1e-3,
        "steady-state μ {:.5} ≠ μ(I) {:.5}",
        mu_measured,
        mu_target
    );
    assert!(mu_measured > params.mu_s && mu_measured < params.mu_2);
}

#[test]
fn plastic_caps_stress_at_yield_surface() {
    // Strong, sustained shear must never let τ̄ exceed μ_2 p (the friction cap).
    let params = MaterialParams::glass_beads_v0();
    let rho = 1.7 * 1000.0;
    let p = pressure(rho, &params);
    let l = [0.0, 1000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    let dt = 1.0e-8;
    let mut s = [0.0; 6];
    for _ in 0..500_000 {
        let out = update_stress(&s, &l, rho, dt, &params);
        s = out.dev_stress;
        let tau = equiv_shear_stress(&s);
        // Always within (slightly above due to high-I) the cap μ_2 p.
        assert!(tau <= params.mu_2 * p + 1e-6 * p);
    }
}

#[test]
fn jaumann_term_is_traceless_and_symmetric_consistent() {
    // s·W − W·s must be deviatoric-preserving: its trace is zero.
    let s: Sym3 = [3.0, -1.0, -2.0, 0.7, 0.4, -0.9];
    let w = [0.3, -0.2, 0.5];
    let r = jaumann_term(&s, &w);
    assert!(trace(&r).abs() < 1e-12);
    // Cross-check against a brute-force 3×3 matrix multiply.
    let sm = [
        [s[XX], s[XY], s[XZ]],
        [s[XY], s[YY], s[YZ]],
        [s[XZ], s[YZ], s[ZZ]],
    ];
    let wm = [
        [0.0, w[0], w[1]],
        [-w[0], 0.0, w[2]],
        [-w[1], -w[2], 0.0],
    ];
    // R = s·w − w·s
    let mut rm = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            let mut acc = 0.0;
            for k in 0..3 {
                acc += sm[i][k] * wm[k][j] - wm[i][k] * sm[k][j];
            }
            rm[i][j] = acc;
        }
    }
    let r_ref: Sym3 = [rm[0][0], rm[1][1], rm[2][2], rm[0][1], rm[0][2], rm[1][2]];
    assert!(dist(&r, &r_ref) < 1e-12);
}

// ── Kinetic-theory branch + granular temperature (§11) ───────────────────────

#[test]
fn kt_pressure_zero_at_zero_temperature() {
    let p = MaterialParams::glass_beads_v0();
    assert_eq!(kt_pressure(1200.0, 0.0, &p), 0.0);
    // positive and ∝ T for T > 0
    let p1 = kt_pressure(1200.0, 0.01, &p);
    let p2 = kt_pressure(1200.0, 0.02, &p);
    assert!(p1 > 0.0);
    assert!((p2 / p1 - 2.0).abs() < 1e-12, "p_KT linear in T");
}

#[test]
fn two_branch_reduces_to_contact_at_zero_temperature() {
    let p = MaterialParams::glass_beads_v0();
    let l = [0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    let s_n = [0.0; 6];
    let rho = 1.05 * p.rho_c;
    let a = update_stress(&s_n, &l, rho, 1e-5, &p);
    let b = two_branch_stress(&s_n, &l, rho, 0.0, 1e-5, &p);
    assert_eq!(a.pressure, b.pressure);
    assert_eq!(a.dev_stress, b.dev_elastic);
    assert_eq!(a.dev_stress, b.dev_total); // no KT viscous stress at T = 0
}

#[test]
fn two_branch_adds_kt_pressure_and_viscous_stress() {
    let p = MaterialParams::glass_beads_v0();
    let rho = 1.05 * p.rho_c;
    let t = 0.02;
    // No shear: pressure gains p_KT, deviatoric unchanged.
    let contact = update_stress(&[0.0; 6], &[0.0; 9], rho, 1e-5, &p);
    let two = two_branch_stress(&[0.0; 6], &[0.0; 9], rho, t, 1e-5, &p);
    let pkt = kt_pressure(rho, t, &p);
    assert!(pkt > 0.0);
    assert!((two.pressure - (contact.pressure + pkt)).abs() < 1e-6);
    assert_eq!(two.dev_total, contact.dev_stress); // no shear → no τ_KT

    // With shear: τ_KT = 2 η_KT D' is added to the total deviatoric only.
    let l = [0.0, 2.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0]; // D_xy = 2
    let two_sh = two_branch_stress(&[0.0; 6], &l, rho, t, 1e-5, &p);
    let eta = kt_shear_viscosity(rho, t, &p);
    assert!(eta > 0.0);
    let added = two_sh.dev_total[XY] - two_sh.dev_elastic[XY];
    assert!((added - 2.0 * eta * 2.0).abs() < 1e-6 * (eta + 1.0), "τ_KT,xy = 2 η_KT D'_xy");
}

#[test]
fn kt_shear_viscosity_scales_with_sqrt_t() {
    let p = MaterialParams::glass_beads_v0();
    let rho = 0.5 * p.rho_s;
    assert_eq!(kt_shear_viscosity(rho, 0.0, &p), 0.0);
    let e1 = kt_shear_viscosity(rho, 0.01, &p);
    let e4 = kt_shear_viscosity(rho, 0.04, &p);
    assert!(e1 > 0.0);
    assert!((e4 / e1 - 2.0).abs() < 1e-9, "η_KT ∝ √T");
}

#[test]
fn steady_shear_gives_bagnold_temperature() {
    // Under constant shear, dT/dt = production + cooling drives T to a steady
    // value where they balance; that T ∝ γ̇² (Bagnold). Integrate two shear rates.
    let p = MaterialParams::glass_beads_v0();
    let rho = 0.5 * p.rho_s;
    let to_steady = |gdot: f64| -> f64 {
        let l = [0.0, gdot, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]; // simple shear
        let mut t = 1.0e-4_f64;
        let dt = 1.0e-6;
        for _ in 0..200_000 {
            let rate = kt_production_rate(rho, t, &l, &p) + kt_cooling_rate(rho, t, &p);
            t = (t + rate * dt).max(1.0e-12);
        }
        t
    };
    let t1 = to_steady(50.0);
    let t2 = to_steady(100.0); // 2× shear rate
    assert!(((t2 / t1) - 4.0).abs() < 0.2, "Bagnold T∝γ̇²: ratio {:.3}", t2 / t1);
    // balance holds at steady
    let l1 = [0.0, 50.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    let prod = kt_production_rate(rho, t1, &l1, &p);
    let bal = prod + kt_cooling_rate(rho, t1, &p);
    assert!(bal.abs() < 1.0e-2 * prod.abs(), "production ≈ dissipation at steady");
}

#[test]
fn homogeneous_cooling_follows_haff_law() {
    // Integrate dT/dt = kt_cooling_rate (= −A T^{3/2}) and compare to the
    // analytical Haff solution T(t) = T0 / (1 + t/τ)², τ = 2/(A√T0).
    let p = MaterialParams::glass_beads_v0();
    let rho = 0.4 * p.rho_s; // dilute granular gas (Φ = 0.4)
    let t0 = 0.01_f64;

    // A from the same formula the code uses (so τ is exact for this material).
    let phi = rho / p.rho_s;
    let g0 = pair_correlation(phi);
    let zeta = 12.0 / std::f64::consts::PI.sqrt();
    let a = 2.0 * zeta * phi * g0 * (1.0 - p.restitution * p.restitution) / (3.0 * p.d);
    let tau = 2.0 / (a * t0.sqrt());

    let dt = tau / 2000.0; // well-resolved
    let mut t = t0;
    let mut time = 0.0;
    for _ in 0..6000 {
        // ~3τ
        t = (t + kt_cooling_rate(rho, t, &p) * dt).max(0.0);
        time += dt;
    }
    let analytic = t0 / (1.0 + time / tau).powi(2);
    let rel_err = (t - analytic).abs() / analytic;
    assert!(
        rel_err < 1e-2,
        "Haff cooling: numeric {t:.6e} vs analytic {analytic:.6e} (rel err {rel_err:.2e})"
    );
}

#[test]
fn kt_conductivity_scales_with_sqrt_t_and_is_positive() {
    let p = MaterialParams::glass_beads_v0();
    let rho = 0.5 * p.rho_s;
    assert_eq!(kt_conductivity(rho, 0.0, &p), 0.0);
    let k1 = kt_conductivity(rho, 0.01, &p);
    let k4 = kt_conductivity(rho, 0.04, &p);
    assert!(k1 > 0.0);
    assert!((k4 / k1 - 2.0).abs() < 1e-9, "κ ∝ √T"); // 4× T → 2× κ
}
