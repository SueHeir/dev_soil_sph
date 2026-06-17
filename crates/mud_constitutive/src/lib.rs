//! MUD — granular constitutive core.
//!
//! Pure, substrate-free. Implements the Dunatunga–Kamrin elasto-viscoplastic
//! stress update with a μ(I) Drucker–Prager yield and density-based tension-free
//! separation, as specified in `docs/physics-design.md` §3. The public entry
//! point is the per-particle update [`update_stress`]:
//! `(s_n, L, ρ, dt, params) → StressOut { s_{n+1}, p, σ }`.
//!
//! The update is a **pure function** — no I/O, no globals, no neighbor coupling
//! — so `mud_physics` can call it per particle inside a SOIL system and so it is
//! unit-testable in isolation (the App. B verification, `docs/physics-design.md`
//! §8.1, lives in `tests` below).
//!
//! ## Conventions
//! - **float64 everywhere** (mandated, §1: viscosity/stiffness span >5 orders).
//! - **Compression-positive pressure**: `p = -⅓ tr σ`; tension is `p < 0` and is
//!   not sustained (separation, §3.1).
//! - Tensors are full 3D. Symmetric tensors (deviatoric stress, strain rate) are
//!   stored as `[f64; 6]` in the order `[xx, yy, zz, xy, xz, yz]`. The velocity
//!   gradient `L = ∇v` is a general `[f64; 9]` in row-major order
//!   `[xx, xy, xz, yx, yy, yz, zx, zy, zz]`.
//!
//! All symmetric-tensor algebra is hand-rolled (std-only, no external crates),
//! matching the sibling stack's `[f64; 3]` array-math style.

#![forbid(unsafe_code)]

/// Symmetric 3×3 tensor stored as the six independent components
/// `[xx, yy, zz, xy, xz, yz]`.
///
/// Used for the deviatoric stress `s` and the strain rate `D`. Off-diagonal
/// entries appear once: e.g. `xy` *is* both `σ_xy` and `σ_yx`.
pub type Sym3 = [f64; 6];

/// Index helpers into a [`Sym3`] for readability.
const XX: usize = 0;
const YY: usize = 1;
const ZZ: usize = 2;
const XY: usize = 3;
const XZ: usize = 4;
const YZ: usize = 5;

// ---------------------------------------------------------------------------
// Material parameters
// ---------------------------------------------------------------------------

/// Material parameters for the granular constitutive update (§3, §7).
///
/// All in SI. `mu_s`, `mu_2`, `i0` are the μ(I) friction law; `rho_s` is the
/// solid-grain density and `rho_c` the critical (close-packed) density; `k_bulk`
/// (`K`) and `g_shear` (`G`) are the weakly-compressible elastic moduli; `d` is
/// the grain diameter.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialParams {
    /// Static friction coefficient μ_s (yield onset).
    pub mu_s: f64,
    /// Limiting friction coefficient μ_2 (high-I plateau).
    pub mu_2: f64,
    /// Inertial-number scale I_0.
    pub i0: f64,
    /// Solid-grain density ρ_s [kg/m³].
    pub rho_s: f64,
    /// Critical (close-packed) density ρ_c [kg/m³]; below it, stress-free.
    pub rho_c: f64,
    /// Effective bulk modulus K [Pa] (weakly compressible, §3.2).
    pub k_bulk: f64,
    /// Elastic shear modulus G [Pa].
    pub g_shear: f64,
    /// Grain diameter d [m].
    pub d: f64,
}

impl MaterialParams {
    /// Derive the shear modulus from bulk modulus `K` and Poisson ratio `ν`:
    /// `G = 3(1−2ν)/(2(1+ν)) · K` (§7).
    pub fn shear_from_bulk_poisson(k_bulk: f64, nu: f64) -> f64 {
        3.0 * (1.0 - 2.0 * nu) / (2.0 * (1.0 + nu)) * k_bulk
    }

    /// Glass-bead v0 parameter set (`docs/physics-design.md` §7).
    ///
    /// μ_s = 0.38, μ_2 = 0.64, I_0 = 0.28, ρ_s = 2500, ρ_c = 1500, ν = 0.3,
    /// d = 0.5 mm.
    ///
    /// **Choice of K.** The true glass bulk modulus (~40 GPa) would cripple the
    /// explicit timestep, so we use the *smallest* K consistent with weak
    /// compressibility (§3.2): pick a target sound speed and set
    /// `K = ρ_c · c_s²`. For a v0 column-collapse regime with `v_max ≈ 5 m/s`
    /// we want `c_s ≥ 10 v_max = 50 m/s` (Mach ≲ 0.1, density fluctuations
    /// ≲ 1%). Taking `c_s = 50 m/s` gives
    /// `K = 1500 · 50² = 3.75 × 10⁶ Pa`, and
    /// `G = 3(1−2ν)/(2(1+ν)) K = (1.2/2.6) K ≈ 1.7308 × 10⁶ Pa`.
    pub fn glass_beads_v0() -> Self {
        let nu = 0.3;
        let rho_c = 1500.0;
        let c_s = 50.0; // target sound speed [m/s] = 10 · v_max (v_max ≈ 5 m/s)
        let k_bulk = rho_c * c_s * c_s; // = 3.75e6 Pa
        MaterialParams {
            mu_s: 0.38,
            mu_2: 0.64,
            i0: 0.28,
            rho_s: 2500.0,
            rho_c,
            k_bulk,
            g_shear: Self::shear_from_bulk_poisson(k_bulk, nu),
            d: 0.5e-3,
        }
    }

    /// Sound speed `c_s = sqrt(K/ρ_c)` implied by the EOS (§3.2).
    pub fn sound_speed(&self) -> f64 {
        (self.k_bulk / self.rho_c).sqrt()
    }

    /// The plastic-flow prefactor `ξ = I_0 / (d √ρ_s)` used in the return map
    /// (§3.3 step 4).
    #[inline]
    pub fn xi(&self) -> f64 {
        self.i0 / (self.d * self.rho_s.sqrt())
    }
}

// ---------------------------------------------------------------------------
// Symmetric-tensor algebra (hand-rolled, std-only)
// ---------------------------------------------------------------------------

/// Trace of a [`Sym3`].
#[inline]
pub fn trace(a: &Sym3) -> f64 {
    a[XX] + a[YY] + a[ZZ]
}

/// Deviatoric part `a' = a − ⅓(tr a) I`.
#[inline]
pub fn deviator(a: &Sym3) -> Sym3 {
    let m = trace(a) / 3.0;
    [a[XX] - m, a[YY] - m, a[ZZ] - m, a[XY], a[XZ], a[YZ]]
}

/// Full tensor contraction `a : b = Σ_ij a_ij b_ij` for symmetric tensors.
/// Off-diagonal components count twice (they represent two entries each).
#[inline]
pub fn double_dot(a: &Sym3, b: &Sym3) -> f64 {
    a[XX] * b[XX]
        + a[YY] * b[YY]
        + a[ZZ] * b[ZZ]
        + 2.0 * (a[XY] * b[XY] + a[XZ] * b[XZ] + a[YZ] * b[YZ])
}

/// Equivalent (von-Mises-like) shear stress `τ̄ = sqrt(½ s : s)` (§3.3).
#[inline]
pub fn equiv_shear_stress(s: &Sym3) -> f64 {
    (0.5 * double_dot(s, s)).sqrt()
}

/// The Jaumann co-rotational term `s·W − W·s` for a symmetric `s` and an
/// antisymmetric spin `W`.
///
/// `w` packs the three independent spin components `[w_xy, w_xz, w_yz]`, i.e.
/// `W = [[0, w_xy, w_xz], [−w_xy, 0, w_yz], [−w_xz, −w_yz, 0]]`. The result
/// `s·W − W·s` is symmetric, returned as a [`Sym3`].
#[inline]
pub fn jaumann_term(s: &Sym3, w: &[f64; 3]) -> Sym3 {
    // Spin matrix W (antisymmetric):
    //   W = [[ 0,    wxy,  wxz],
    //        [-wxy,  0,    wyz],
    //        [-wxz, -wyz,  0  ]]
    let (wxy, wxz, wyz) = (w[0], w[1], w[2]);
    let (sxx, syy, szz, sxy, sxz, syz) = (s[XX], s[YY], s[ZZ], s[XY], s[XZ], s[YZ]);

    // R = s·W − W·s, with
    //   W = [[ 0,    wxy,  wxz],
    //        [-wxy,  0,    wyz],
    //        [-wxz, -wyz,  0  ]].
    // (s·W)_ik = Σ_j s_ij W_jk ; (W·s)_ik = Σ_j W_ij s_jk. R is symmetric.
    //
    // R_xx = -2(sxy·wxy + sxz·wxz)
    // R_yy =  2(sxy·wxy - syz·wyz)
    // R_zz =  2(sxz·wxz + syz·wyz)
    // R_xy = (sxx - syy)·wxy - sxz·wyz - syz·wxz
    // R_xz = (sxx - szz)·wxz + sxy·wyz - syz·wxy
    // R_yz = (syy - szz)·wyz + sxy·wxz + sxz·wxy
    let r_xx = -2.0 * (sxy * wxy + sxz * wxz);
    let r_yy = 2.0 * (sxy * wxy - syz * wyz);
    let r_zz = 2.0 * (sxz * wxz + syz * wyz);
    let r_xy = (sxx - syy) * wxy - sxz * wyz - syz * wxz;
    let r_xz = (sxx - szz) * wxz + sxy * wyz - syz * wxy;
    let r_yz = (syy - szz) * wyz + sxy * wxz + sxz * wxy;

    [r_xx, r_yy, r_zz, r_xy, r_xz, r_yz]
}

/// Decompose the velocity gradient `L = ∇v` (row-major `[f64; 9]`:
/// `[L_xx, L_xy, L_xz, L_yx, L_yy, L_yz, L_zx, L_zy, L_zz]`) into the symmetric
/// strain rate `D = ½(L + Lᵀ)` (as [`Sym3`]) and the spin `W = ½(L − Lᵀ)`
/// (packed `[w_xy, w_xz, w_yz]`). See §2.
#[inline]
pub fn decompose_velocity_gradient(l: &[f64; 9]) -> (Sym3, [f64; 3]) {
    let (lxx, lxy, lxz) = (l[0], l[1], l[2]);
    let (lyx, lyy, lyz) = (l[3], l[4], l[5]);
    let (lzx, lzy, lzz) = (l[6], l[7], l[8]);

    let d: Sym3 = [
        lxx,
        lyy,
        lzz,
        0.5 * (lxy + lyx),
        0.5 * (lxz + lzx),
        0.5 * (lyz + lzy),
    ];
    // W_ij = ½(L_ij − L_ji): w_xy = ½(lxy − lyx), etc.
    let w = [
        0.5 * (lxy - lyx),
        0.5 * (lxz - lzx),
        0.5 * (lyz - lzy),
    ];
    (d, w)
}

// ---------------------------------------------------------------------------
// Constitutive scalar laws
// ---------------------------------------------------------------------------

/// Granular EOS with tension-free separation (§3.2):
/// `p(ρ) = 0` if `ρ < ρ_c`, else `K (ρ/ρ_c − 1)`.
///
/// Compression-positive. Note this can still return `p ≤ 0` exactly at
/// `ρ = ρ_c`; the caller treats `p ≤ 0` as disconnected (§3.3 step 1).
#[inline]
pub fn pressure(rho: f64, params: &MaterialParams) -> f64 {
    if rho < params.rho_c {
        0.0
    } else {
        params.k_bulk * (rho / params.rho_c - 1.0)
    }
}

/// μ(I) friction law (Jop form, §3.3):
/// `μ(I) = μ_s + (μ_2 − μ_s) / (I_0/I + 1)`.
///
/// At `I = 0` returns `μ_s`; as `I → ∞` it approaches `μ_2`.
#[inline]
pub fn mu_of_i(i: f64, params: &MaterialParams) -> f64 {
    if i <= 0.0 {
        params.mu_s
    } else {
        params.mu_s + (params.mu_2 - params.mu_s) / (params.i0 / i + 1.0)
    }
}

// ---------------------------------------------------------------------------
// The per-particle stress update (§3.3)
// ---------------------------------------------------------------------------

/// Result of [`update_stress`]: the updated deviatoric stress, the pressure, the
/// full Cauchy stress, and the accumulated plastic shear rate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StressOut {
    /// Updated deviatoric stress `s_{n+1}` ([`Sym3`]).
    pub dev_stress: Sym3,
    /// Pressure `p` (compression-positive). Zero when disconnected.
    pub pressure: f64,
    /// Full Cauchy stress `σ = −p I + s` ([`Sym3`]).
    pub sigma: Sym3,
    /// Plastic shear strain-rate `γ̇ᵖ` realized this step (0 if elastic).
    pub gamma_dot_p: f64,
    /// True if the point was disconnected this step (`ρ < ρ_c` or `p ≤ 0`):
    /// stress-free.
    pub disconnected: bool,
}

impl StressOut {
    /// The all-zero, disconnected state (§3.1 first branch).
    #[inline]
    fn disconnected() -> Self {
        StressOut {
            dev_stress: [0.0; 6],
            pressure: 0.0,
            sigma: [0.0; 6],
            gamma_dot_p: 0.0,
            disconnected: true,
        }
    }
}

/// Per-particle granular constitutive update — Dunatunga–Kamrin elasto-
/// viscoplastic return map cast for weakly-compressible SPH (§3.3).
///
/// Pure function: given the old deviatoric stress `s_n`, the velocity gradient
/// `l = ∇v` (row-major `[f64; 9]`), the *already-updated* density `rho`, the
/// step `dt`, and the material `params`, returns the new state ([`StressOut`]).
///
/// # Algorithm (exactly §3.3)
/// 1. **Density & pressure first.** `p = pressure(rho)`. If `rho < ρ_c` or
///    `p ≤ 0`, return stress-free (disconnected).
/// 2. **Jaumann elastic trial.** With `D, W` from `L`,
///    `s_tr = s_n + dt (2G D' + s_n·W − W·s_n)`; `τ̄_tr = sqrt(½ s_tr:s_tr)`.
/// 3. **Yield check.** If `τ̄_tr ≤ μ_s p`: elastic, `s_{n+1} = s_tr`.
/// 4. **Else plastic.** With `ξ = I_0/(d√ρ_s)`,
///    `S0 = μ_s p`, `S2 = μ_2 p`, `α = ξ G dt √p`,
///    `B = S2 + τ̄_tr + α`, `H = S2 τ̄_tr + S0 α`,
///    the cancellation-safe root `τ̄_{n+1} = 2H/(B + √(B² − 4H))`,
///    `γ̇ᵖ = (τ̄_tr − τ̄_{n+1})/(G dt)`,
///    radial return `s_{n+1} = (τ̄_{n+1}/τ̄_tr) s_tr`.
/// 5. **Reassemble** `σ = −p I + s_{n+1}`.
pub fn update_stress(
    s_n: &Sym3,
    l: &[f64; 9],
    rho: f64,
    dt: f64,
    params: &MaterialParams,
) -> StressOut {
    // --- Step 1: density & pressure first; separation check. ---
    let p = pressure(rho, params);
    if rho < params.rho_c || p <= 0.0 {
        return StressOut::disconnected();
    }

    let g = params.g_shear;
    let (d, w) = decompose_velocity_gradient(l);
    let d_dev = deviator(&d);

    // --- Step 2: Jaumann elastic trial deviator. ---
    let jt = jaumann_term(s_n, &w);
    let s_tr: Sym3 = [
        s_n[XX] + dt * (2.0 * g * d_dev[XX] + jt[XX]),
        s_n[YY] + dt * (2.0 * g * d_dev[YY] + jt[YY]),
        s_n[ZZ] + dt * (2.0 * g * d_dev[ZZ] + jt[ZZ]),
        s_n[XY] + dt * (2.0 * g * d_dev[XY] + jt[XY]),
        s_n[XZ] + dt * (2.0 * g * d_dev[XZ] + jt[XZ]),
        s_n[YZ] + dt * (2.0 * g * d_dev[YZ] + jt[YZ]),
    ];
    let tau_tr = equiv_shear_stress(&s_tr);

    let s0 = params.mu_s * p; // yield threshold μ_s p

    // --- Step 3: yield check. ---
    let (s_next, gamma_dot_p) = if tau_tr <= s0 || tau_tr == 0.0 {
        // Elastic (also the trivial τ̄_tr = 0 case avoids a 0/0 radial return).
        (s_tr, 0.0)
    } else {
        // --- Step 4: plastic return map (cancellation-safe root). ---
        let s2 = params.mu_2 * p;
        let alpha = params.xi() * g * dt * p.sqrt();
        let b = s2 + tau_tr + alpha;
        let h = s2 * tau_tr + s0 * alpha;
        // τ̄_{n+1} = 2H / (B + √(B² − 4H)). B > 0 and H > 0 here, and the
        // discriminant is non-negative for physical inputs; clamp to be safe.
        let disc = (b * b - 4.0 * h).max(0.0);
        let tau_next = 2.0 * h / (b + disc.sqrt());
        let gamma_dot_p = (tau_tr - tau_next) / (g * dt);
        let scale = tau_next / tau_tr; // radial return
        let s_next: Sym3 = [
            scale * s_tr[XX],
            scale * s_tr[YY],
            scale * s_tr[ZZ],
            scale * s_tr[XY],
            scale * s_tr[XZ],
            scale * s_tr[YZ],
        ];
        (s_next, gamma_dot_p)
    };

    // --- Step 5: reassemble σ = −p I + s. ---
    let sigma: Sym3 = [
        -p + s_next[XX],
        -p + s_next[YY],
        -p + s_next[ZZ],
        s_next[XY],
        s_next[XZ],
        s_next[YZ],
    ];

    StressOut {
        dev_stress: s_next,
        pressure: p,
        sigma,
        gamma_dot_p,
        disconnected: false,
    }
}

#[cfg(test)]
mod tests;
