//! MUD — Wendland smoothing kernels for SPH.
//!
//! Pure, substrate-free. Provides the Wendland C2 smoothing kernel `W(r, h)` and
//! its gradient `∇W` in 2D and 3D. See `docs/sph-primer.md` §4 and
//! `docs/architecture.md` (`mud_kernel`).
//!
//! # The kernel
//!
//! The Wendland C2 kernel is the modern SPH default: its Fourier transform is
//! non-negative, which suppresses the pairing/tensile instability that plagues
//! the cubic spline. Support radius is `2h`; the kernel vanishes for `r ≥ 2h`.
//!
//! With `q = r/h` and `0 ≤ q < 2`:
//!
//! - **3D**: `W(q) = (21 / (16 π h³)) · (1 − q/2)⁴ · (2q + 1)`
//! - **2D**: `W(q) = (7  / (4  π h²)) · (1 − q/2)⁴ · (2q + 1)`
//!
//! The gradient (with respect to the position of particle `i`, displacement
//! `dx = x_i − x_j`) is derived analytically below and is well-behaved as
//! `r → 0`, where it tends to the zero vector.
//!
//! # Usage in a pair loop
//!
//! ```
//! use mud_kernel::Kernel;
//!
//! let k = Kernel::Dim3;
//! let h = 0.1;
//! let dx = [0.05, 0.0, 0.0];          // x_i − x_j
//! let _w = k.w(0.05, h);              // scalar kernel value
//! let _grad = k.grad_w(dx, h);        // ∇_i W_ij
//! assert!(k.support_radius(h) == 2.0 * h);
//! ```

use std::f64::consts::PI;

/// A Wendland C2 SPH smoothing kernel, parameterized by spatial dimension.
///
/// All methods are allocation-free pure functions. Construct one of the two
/// variants and call [`w`](Kernel::w), [`grad_w`](Kernel::grad_w),
/// [`support_radius`](Kernel::support_radius), or [`w_zero`](Kernel::w_zero).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Kernel {
    /// Two-dimensional Wendland C2 kernel (normalization `7 / (4 π h²)`).
    Dim2,
    /// Three-dimensional Wendland C2 kernel (normalization `21 / (16 π h³)`).
    Dim3,
}

impl Kernel {
    /// The compact-support radius of the kernel: `2h`.
    ///
    /// `W(r, h) = 0` and `grad_w(dx, h) = 0` for `r ≥ 2h`.
    #[inline]
    pub fn support_radius(self, h: f64) -> f64 {
        2.0 * h
    }

    /// Dimension-dependent normalization constant `α_d` such that
    /// `W(q) = (α_d / h^d) · (1 − q/2)⁴ (2q + 1)`.
    ///
    /// `α_2 = 7 / (4π)`, `α_3 = 21 / (16π)`.
    #[inline]
    fn alpha(self) -> f64 {
        match self {
            Kernel::Dim2 => 7.0 / (4.0 * PI),
            Kernel::Dim3 => 21.0 / (16.0 * PI),
        }
    }

    /// The full per-`h` normalization `α_d / h^d`.
    #[inline]
    fn norm(self, h: f64) -> f64 {
        match self {
            Kernel::Dim2 => self.alpha() / (h * h),
            Kernel::Dim3 => self.alpha() / (h * h * h),
        }
    }

    /// The smoothing kernel value `W(r, h)` for a non-negative scalar distance `r`.
    ///
    /// Returns `0.0` for `r ≥ 2h` (compact support). `r` is expected to be
    /// non-negative; negative inputs are treated as their magnitude.
    #[inline]
    pub fn w(self, r: f64, h: f64) -> f64 {
        let q = r.abs() / h;
        if q >= 2.0 {
            return 0.0;
        }
        // (1 − q/2)⁴ (2q + 1)
        let t = 1.0 - 0.5 * q;
        let t4 = t * t * t * t;
        self.norm(h) * t4 * (2.0 * q + 1.0)
    }

    /// The kernel value at zero separation, `W(0, h)` — useful for self-contribution
    /// terms in SPH sums. Equal to `α_d / h^d` since `(1)⁴ (1) = 1` at `q = 0`.
    #[inline]
    pub fn w_zero(self, h: f64) -> f64 {
        self.norm(h)
    }

    /// The radial derivative `dW/dr` at scalar distance `r`.
    ///
    /// Derivation: with `q = r/h` and `f(q) = (1 − q/2)⁴ (2q + 1)`,
    ///
    /// ```text
    /// f'(q) = −2(1 − q/2)³ (2q + 1) + 2(1 − q/2)⁴
    ///       = 2(1 − q/2)³ [ −(2q + 1) + (1 − q/2) ]
    ///       = −5 q (1 − q/2)³.
    /// ```
    ///
    /// Hence `dW/dr = (1/h) · (α_d / h^d) · f'(q) = −5 q (1 − q/2)³ · α_d / h^{d+1}`,
    /// which is `≤ 0` on `(0, 2h)` (the kernel decreases monotonically) and is
    /// `0` at `r = 0` and `r ≥ 2h`.
    #[inline]
    pub fn dw_dr(self, r: f64, h: f64) -> f64 {
        let q = r.abs() / h;
        if q >= 2.0 || q == 0.0 {
            return 0.0;
        }
        let t = 1.0 - 0.5 * q;
        let t3 = t * t * t;
        // (α_d / h^d) · (1/h) · (−5 q t³)
        self.norm(h) * (-5.0 * q * t3) / h
    }

    /// The kernel gradient `∇_i W_ij`, with respect to the position of particle `i`.
    ///
    /// `dx = x_i − x_j` is the displacement vector (always pass all three
    /// components; for 2D leave the third as `0.0`). The result is
    ///
    /// ```text
    /// ∇_i W = (dW/dr) · (dx / r) = [ −5 (1 − q/2)³ · α_d / h^{d+2} ] · dx,
    /// ```
    ///
    /// where the bracketed scalar is finite at `r = 0`, so the gradient tends
    /// smoothly to the zero vector there (it is `dx` times a finite scalar, and
    /// `dx = 0`). Returns the zero vector for `r ≥ 2h`.
    ///
    /// Because `dW/dr ≤ 0`, the gradient points along `−r̂` (i.e. toward the
    /// kernel center, opposite to `dx`), as expected for a monotonically
    /// decreasing kernel.
    #[inline]
    pub fn grad_w(self, dx: [f64; 3], h: f64) -> [f64; 3] {
        let r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
        let r = r2.sqrt();
        let q = r / h;
        if q >= 2.0 || r == 0.0 {
            return [0.0, 0.0, 0.0];
        }
        // factor = (dW/dr) / r = −5 (1 − q/2)³ · α_d / h^{d+2}
        let t = 1.0 - 0.5 * q;
        let t3 = t * t * t;
        let factor = self.norm(h) * (-5.0 * t3) / (h * h);
        [factor * dx[0], factor * dx[1], factor * dx[2]]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const TWO_PI: f64 = 2.0 * PI;

    // ── Normalization: ∫ W dV ≈ 1 (numerical quadrature) ─────────────────────

    #[test]
    fn normalization_3d() {
        let k = Kernel::Dim3;
        let h = 1.0;
        let rmax = 2.0 * h;
        // Integrate radially: ∫₀^{2h} W(r) · 4π r² dr (spherical shells).
        let n = 200_000;
        let dr = rmax / n as f64;
        let mut integral = 0.0;
        for i in 0..n {
            let r = (i as f64 + 0.5) * dr; // midpoint rule
            integral += k.w(r, h) * 4.0 * PI * r * r * dr;
        }
        assert!(
            (integral - 1.0).abs() < 1e-3,
            "3D normalization integral = {integral}, expected ≈ 1"
        );
    }

    #[test]
    fn normalization_2d() {
        let k = Kernel::Dim2;
        let h = 1.0;
        let rmax = 2.0 * h;
        // Integrate radially: ∫₀^{2h} W(r) · 2π r dr (annular rings).
        let n = 200_000;
        let dr = rmax / n as f64;
        let mut integral = 0.0;
        for i in 0..n {
            let r = (i as f64 + 0.5) * dr;
            integral += k.w(r, h) * TWO_PI * r * dr;
        }
        assert!(
            (integral - 1.0).abs() < 1e-3,
            "2D normalization integral = {integral}, expected ≈ 1"
        );
    }

    // ── Compact support: W and grad_w vanish for r ≥ 2h ──────────────────────

    #[test]
    fn compact_support() {
        for k in [Kernel::Dim2, Kernel::Dim3] {
            let h = 0.37;
            let support = k.support_radius(h);
            assert_eq!(support, 2.0 * h);
            // exactly at the support radius and beyond → zero
            for &r in &[support, support + 1e-9, support * 1.5, 10.0 * h] {
                assert_eq!(k.w(r, h), 0.0, "W should be 0 at r = {r} (≥ 2h)");
                let grad = k.grad_w([r, 0.0, 0.0], h);
                assert_eq!(grad, [0.0, 0.0, 0.0], "grad_w should be 0 at r = {r}");
            }
            // just inside support → strictly positive
            assert!(k.w(support - 1e-6, h) > 0.0);
        }
    }

    // ── Positivity and monotonic decrease on (0, 2h) ─────────────────────────

    #[test]
    fn positive_and_monotonic() {
        for k in [Kernel::Dim2, Kernel::Dim3] {
            let h = 0.5;
            let n = 1000;
            let mut prev = f64::INFINITY;
            for i in 0..n {
                let r = (i as f64 / n as f64) * 2.0 * h; // [0, 2h)
                let w = k.w(r, h);
                assert!(w > 0.0, "W must be positive on [0, 2h): W({r}) = {w}");
                assert!(
                    w <= prev + 1e-15,
                    "W must be non-increasing: W({r}) = {w} > previous {prev}"
                );
                prev = w;
            }
            // strict: the value at 0 exceeds the value near the edge
            assert!(k.w(0.0, h) > k.w(2.0 * h - 1e-3, h));
        }
    }

    // ── Gradient consistency: sign and finite-difference magnitude ───────────

    #[test]
    fn gradient_sign_points_inward() {
        // grad_w(dx) should point along −dx (toward center) because dW/dr < 0.
        for k in [Kernel::Dim2, Kernel::Dim3] {
            let h = 0.8;
            let dx = [0.6, 0.0, 0.0]; // r = 0.6 < 2h = 1.6
            let grad = k.grad_w(dx, h);
            assert!(grad[0] < 0.0, "grad_w should oppose +x displacement");
            assert_eq!(grad[1], 0.0);
            assert_eq!(grad[2], 0.0);
            // and dW/dr must be negative there
            assert!(k.dw_dr(0.6, h) < 0.0);
        }
    }

    #[test]
    fn gradient_magnitude_matches_finite_difference() {
        for k in [Kernel::Dim2, Kernel::Dim3] {
            let h = 0.6;
            let eps = 1e-6;
            // sample a set of distances strictly inside the support
            for i in 1..30 {
                let r = i as f64 / 30.0 * (2.0 * h) * 0.98 + 1e-3;
                // central finite difference of W(r) w.r.t. r
                let fd = (k.w(r + eps, h) - k.w(r - eps, h)) / (2.0 * eps);
                // analytic |grad_w| along a radial displacement equals |dW/dr|
                let grad = k.grad_w([r, 0.0, 0.0], h);
                let gmag = (grad[0] * grad[0] + grad[1] * grad[1] + grad[2] * grad[2]).sqrt();
                // |grad_w| should equal |dW/dr| = |fd|
                assert!(
                    (gmag - fd.abs()).abs() < 1e-4,
                    "at r={r}: |grad_w|={gmag} vs |dW/dr|≈{}",
                    fd.abs()
                );
                // and the analytic dw_dr should match the finite difference too
                assert!(
                    (k.dw_dr(r, h) - fd).abs() < 1e-4,
                    "at r={r}: dw_dr={} vs fd={fd}",
                    k.dw_dr(r, h)
                );
            }
        }
    }

    // ── Symmetry: W is radial; grad_w is odd in dx ───────────────────────────

    #[test]
    fn symmetry() {
        let k = Kernel::Dim3;
        let h = 0.45;
        // W depends only on |r|: equal-magnitude displacements give equal W.
        let r = 0.5;
        assert_eq!(k.w(r, h), k.w(-r, h));

        // grad_w(dx) = −grad_w(−dx)
        for dx in [
            [0.3, -0.2, 0.1],
            [0.7, 0.0, 0.0],
            [0.0, 0.0, 0.55],
            [-0.4, 0.1, -0.3],
        ] {
            let g_pos = k.grad_w(dx, h);
            let g_neg = k.grad_w([-dx[0], -dx[1], -dx[2]], h);
            for d in 0..3 {
                assert!(
                    (g_pos[d] + g_neg[d]).abs() < 1e-15,
                    "grad_w should be odd: {g_pos:?} vs {g_neg:?}"
                );
            }
        }

        // also: W from grad over equal-|dx| in different directions is consistent
        let a = k.grad_w([0.5, 0.0, 0.0], h);
        let b = k.grad_w([0.0, 0.5, 0.0], h);
        let amag = (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt();
        let bmag = (b[0] * b[0] + b[1] * b[1] + b[2] * b[2]).sqrt();
        assert!((amag - bmag).abs() < 1e-15);
    }

    // ── Partition of unity on a regular lattice ──────────────────────────────
    //
    // For a uniform lattice with spacing Δ, each particle carries volume
    // V = Δ^d (= m/ρ). The SPH interpolation of the constant field 1 is
    // Σ_j V_j W(r_ij) ≈ ∫ W dV = 1 — exact in the continuum limit, with a
    // lattice-quadrature residual that shrinks as Δ/h → 0.
    //
    // Achieved accuracy (documented), summing a particle's full 2h support on a
    // cubic lattice:
    //   spacing = h     (coarse, ~33 neighbors 3D):  err ≈ 3.4% (3D), 3.8% (2D)
    //   spacing = h/2   (realistic SPH resolution):  err ≈ 0.06% (3D), 0.3% (2D)
    //   spacing = h/4   (fine):                      err < 5e-4 in both
    // We assert at the realistic h/2 resolution (the regime SPH actually runs at,
    // ~250 neighbors in 3D within 2h) and require sub-percent accuracy.

    #[test]
    fn partition_of_unity_lattice_3d() {
        let k = Kernel::Dim3;
        let h: f64 = 1.0;
        let spacing = 0.5 * h; // realistic SPH resolution (Δ/h = 0.5)
        let vol = spacing * spacing * spacing; // m/ρ per particle
        let reach = (2.0 * h / spacing).ceil() as i64 + 1;

        let mut sum = 0.0;
        for ix in -reach..=reach {
            for iy in -reach..=reach {
                for iz in -reach..=reach {
                    let dx = ix as f64 * spacing;
                    let dy = iy as f64 * spacing;
                    let dz = iz as f64 * spacing;
                    let r = (dx * dx + dy * dy + dz * dz).sqrt();
                    sum += vol * k.w(r, h);
                }
            }
        }
        let err = (sum - 1.0).abs();
        assert!(
            err < 5e-3,
            "3D lattice partition-of-unity Σ V W = {sum} (err {err})"
        );
        eprintln!("3D lattice partition-of-unity (Δ/h=0.5): Σ V W = {sum}, err = {err:e}");
    }

    #[test]
    fn partition_of_unity_lattice_2d() {
        let k = Kernel::Dim2;
        let h: f64 = 1.0;
        let spacing = 0.5 * h; // realistic SPH resolution (Δ/h = 0.5)
        let vol = spacing * spacing; // m/ρ per particle (area)
        let reach = (2.0 * h / spacing).ceil() as i64 + 1;

        let mut sum = 0.0;
        for ix in -reach..=reach {
            for iy in -reach..=reach {
                let dx = ix as f64 * spacing;
                let dy = iy as f64 * spacing;
                let r = (dx * dx + dy * dy).sqrt();
                sum += vol * k.w(r, h);
            }
        }
        let err = (sum - 1.0).abs();
        assert!(
            err < 5e-3,
            "2D lattice partition-of-unity Σ V W = {sum} (err {err})"
        );
        eprintln!("2D lattice partition-of-unity (Δ/h=0.5): Σ V W = {sum}, err = {err:e}");
    }

    // ── w_zero and support_radius helpers ────────────────────────────────────

    #[test]
    fn w_zero_matches_w_at_origin() {
        for k in [Kernel::Dim2, Kernel::Dim3] {
            let h = 0.9;
            assert!((k.w_zero(h) - k.w(0.0, h)).abs() < 1e-15);
            assert!(k.w_zero(h) > 0.0);
        }
    }
}
