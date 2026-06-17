//! MUD — granular constitutive core.
//!
//! Pure, substrate-free. Implements the Dunatunga–Kamrin elasto-viscoplastic
//! stress update with a μ(I) Drucker–Prager yield and density-based tension-free
//! separation, as specified in `docs/physics-design.md` §3. The public entry
//! point is the per-particle update `(s_n, L, ρ, dt, params) → (s_{n+1}, p, σ)`.
//!
//! TODO(stub): implemented by the constitutive coding pass.
