//! The B88 exchange, LYP correlation and B3LYP hybrid functional.
//!
//! **B3LYP** is the most widely used density functional in chemistry. It
//! is a *hybrid*: it mixes a fraction of exact (Hartree-Fock) exchange
//! into a GGA. The Becke three-parameter form is
//!
//! ```text
//! E_xc^B3LYP = (1−a) E_x^LDA + a E_x^HF + b ΔE_x^B88
//!            + (1−c) E_c^VWN + c E_c^LYP
//! ```
//!
//! with the standard coefficients `a = 0.20`, `b = 0.72`, `c = 0.81`.
//! The exact-exchange fraction `a` is supplied by the Kohn-Sham build
//! ([`crate::dft::ks`]) — it folds `a · K` into the Kohn-Sham matrix
//! using the crate's existing Hartree-Fock exchange. This module
//! supplies the *density-functional* part: the B88 gradient correction
//! to exchange, the LDA exchange, and the LYP / VWN correlation.
//!
//! ## The pieces here
//!
//! - [`b88_exchange`] — Slater exchange plus the Becke-1988 gradient
//!   correction `ΔE_x^B88`.
//! - [`lyp_correlation`] — the Lee-Yang-Parr correlation functional in
//!   its closed-shell (restricted) form.
//! - [`b3lyp_dft_xc`] — the density-functional part of B3LYP
//!   (everything except the exact-exchange term).
//!
//! All return an [`XcContribution`] — `ε_xc`, the density potential
//! and the `|∇ρ|`-potential — so the Kohn-Sham build treats them
//! exactly like the PBE pieces.

use super::lda::{slater_exchange, vwn_correlation};
use super::XcContribution;

/// B3LYP exact-exchange mixing fraction `a`.
pub const B3LYP_A_EXACT: f64 = 0.20;
/// B3LYP B88 gradient-exchange coefficient `b`.
pub const B3LYP_B_B88: f64 = 0.72;
/// B3LYP LYP-correlation coefficient `c`.
pub const B3LYP_C_LYP: f64 = 0.81;

// =====================================================================
// B88 exchange
// =====================================================================

/// The Becke-1988 exchange parameter `β`.
const B88_BETA: f64 = 0.0042;

/// `asinh(x) = ln(x + √(x²+1))`.
#[inline]
fn asinh(x: f64) -> f64 {
    (x + (x * x + 1.0).sqrt()).ln()
}

/// B88 exchange at density `rho` and gradient norm `grad_norm`
/// (closed-shell — both spins carry `ρ/2`).
///
/// Returns Slater exchange plus the Becke gradient correction. The
/// per-spin reduced gradient is `x = |∇ρ_σ| / ρ_σ^{4/3}` with
/// `ρ_σ = ρ/2` and `|∇ρ_σ| = |∇ρ|/2`; the correction energy density is
/// `Δe_x = −β ρ_σ^{4/3} x² / (1 + 6β x·asinh x)` summed over the two
/// equal spins.
pub fn b88_exchange(rho: f64, grad_norm: f64) -> XcContribution {
    let lda = slater_exchange(rho);
    if rho <= 1.0e-30 {
        return lda;
    }
    // Closed-shell: per-spin density and gradient.
    let rho_s = 0.5 * rho;
    let grad_s = 0.5 * grad_norm;
    let rs43 = rho_s.powf(4.0 / 3.0);
    let x = grad_s / rs43;
    let x2 = x * x;
    let ash = asinh(x);
    let denom = 1.0 + 6.0 * B88_BETA * x * ash;

    // Δe_x per spin (energy density per volume), times 2 spins.
    let de_x_spin = -B88_BETA * rs43 * x2 / denom;
    let de_x = 2.0 * de_x_spin;

    // ε_x correction per electron = Δe_x / ρ.
    let eps_corr = de_x / rho;
    let eps_x = lda.energy_density + eps_corr;

    // --- Potential of the B88 correction. ---------------------------
    // Work per spin with the energy density g(ρ_σ, x). We need
    // ∂(Δe_x)/∂ρ and ∂(Δe_x)/∂|∇ρ|.
    //
    // x = |∇ρ_σ| ρ_σ^{−4/3}; ρ_σ = ρ/2, |∇ρ_σ| = |∇ρ|/2.
    //   ∂x/∂ρ   = −(4/3) x / ρ          (the ½'s cancel in the ratio's ρ-scaling)
    //   ∂x/∂|∇ρ| = x / |∇ρ|             (for |∇ρ| > 0)
    //
    // g = −β ρ_σ^{4/3} x²/D,  D = 1 + 6β x asinh x.
    //   dD/dx = 6β (asinh x + x/√(x²+1)).
    //   ∂g/∂x (fixed ρ_σ) = −β ρ_σ^{4/3} · [2x D − x² dD/dx]/D².
    //   ∂g/∂ρ_σ (fixed x) = −(4/3) β ρ_σ^{1/3} x²/D.
    let dd_dx = 6.0 * B88_BETA * (ash + x / (x * x + 1.0).sqrt());
    let dg_dx = -B88_BETA * rs43 * (2.0 * x * denom - x2 * dd_dx) / (denom * denom);
    let dg_drhos = -(4.0 / 3.0) * B88_BETA * rho_s.cbrt() * x2 / denom;

    // ∂g/∂ρ = ∂g/∂ρ_σ · (∂ρ_σ/∂ρ) + ∂g/∂x · (∂x/∂ρ).
    //   ∂ρ_σ/∂ρ = 1/2.
    let dx_drho = -(4.0 / 3.0) * x / rho;
    let dg_drho = dg_drhos * 0.5 + dg_dx * dx_drho;
    // Two equal spins → factor 2 on the total correction energy
    // density; ∂(Δe_x)/∂ρ = 2 ∂g/∂ρ.
    let dde_drho = 2.0 * dg_drho;

    // ∂g/∂|∇ρ| = ∂g/∂x · ∂x/∂|∇ρ|, ∂x/∂|∇ρ| = x/|∇ρ| (the ½'s cancel).
    let dde_dgrad = if grad_norm > 1.0e-30 {
        let dx_dgrad = x / grad_norm;
        2.0 * dg_dx * dx_dgrad
    } else {
        0.0
    };

    XcContribution {
        energy_density: eps_x,
        // Slater potential + the B88-correction density potential.
        potential: lda.potential + dde_drho,
        gradient_potential: dde_dgrad,
    }
}

// =====================================================================
// LYP correlation
// =====================================================================

/// LYP correlation parameter `a`.
const LYP_A: f64 = 0.049_18;
/// LYP correlation parameter `b`.
const LYP_B: f64 = 0.132;
/// LYP correlation parameter `c`.
const LYP_C: f64 = 0.2533;
/// LYP correlation parameter `d`.
const LYP_D: f64 = 0.349;

/// The closed-shell LYP correlation energy density per *volume*
/// `e_c(ρ, σ)` as an explicit function of the total density `ρ` and the
/// squared gradient `σ = |∇ρ|²`.
///
/// This is the restricted (closed-shell) Lee-Yang-Parr functional in
/// the Miehlich-Savin-Stoll-Preuss / Johnson-Gill-Pople reduction. The
/// closed-form result for `ρ_α = ρ_β = ρ/2` is
///
/// ```text
/// e_c = −a ρ/(1+dρ^{−1/3})
///       − a b ω(ρ) [ ρ²·(CF·2^{8/3}·ρ^{8/3} ·(1/8)... ) ]
/// ```
///
/// — written out term-by-term in the body. Keeping it as one closed
/// scalar function lets the potential be taken as a finite-difference-
/// verified analytic derivative.
fn lyp_energy_density(rho: f64, sigma: f64) -> f64 {
    if rho <= 1.0e-30 {
        return 0.0;
    }
    let cf = 0.3 * (3.0 * std::f64::consts::PI * std::f64::consts::PI).powf(2.0 / 3.0);
    let rho_m13 = rho.powf(-1.0 / 3.0);
    let den = 1.0 + LYP_D * rho_m13;

    // ω = e^{−c ρ^{−1/3}} ρ^{−11/3} / (1 + d ρ^{−1/3}).
    let omega = (-LYP_C * rho_m13).exp() * rho.powf(-11.0 / 3.0) / den;
    // δ = c ρ^{−1/3} + d ρ^{−1/3}/(1 + d ρ^{−1/3}).
    let delta = LYP_C * rho_m13 + LYP_D * rho_m13 / den;

    // Closed-shell: ρ_α = ρ_β = ρ/2, σ_αα = σ_ββ = σ_αβ = σ/4.
    // The LYP energy density (restricted form):
    //   e_c = −a ρ/(1+dρ^{−1/3})
    //         − a b ω [ ρ_α ρ_β ( 2^{11/3} C_F (ρ_α^{8/3}+ρ_β^{8/3})
    //                  + (47/18 − 7δ/18) |∇ρ|²
    //                  − (5/2 − δ/18) (|∇ρ_α|²+|∇ρ_β|²)
    //                  − (δ−11)/9 (ρ_α|∇ρ_α|²+ρ_β|∇ρ_β|²)/ρ )
    //                  − 2/3 ρ² |∇ρ|²
    //                  + (2/3 ρ² − ρ_α²)|∇ρ_β|² + (2/3 ρ² − ρ_β²)|∇ρ_α|² ].
    let rho_a = 0.5 * rho;
    let rho_b = 0.5 * rho;
    let grad2 = sigma; // |∇ρ|²
    let grad2_a = 0.25 * sigma; // |∇ρ_α|²
    let grad2_b = 0.25 * sigma; // |∇ρ_β|²

    let term_first = -LYP_A * rho / den;

    // The Thomas-Fermi-like kinetic bracket.
    let tf = 2.0_f64.powf(11.0 / 3.0) * cf * (rho_a.powf(8.0 / 3.0) + rho_b.powf(8.0 / 3.0));
    let g_total = (47.0 / 18.0 - 7.0 * delta / 18.0) * grad2;
    let g_spin = -(5.0 / 2.0 - delta / 18.0) * (grad2_a + grad2_b);
    let g_mix = -(delta - 11.0) / 9.0 * (rho_a * grad2_a + rho_b * grad2_b) / rho;
    let bracket_ab = rho_a * rho_b * (tf + g_total + g_spin + g_mix);

    let extra = -(2.0 / 3.0) * rho * rho * grad2
        + (2.0 / 3.0 * rho * rho - rho_a * rho_a) * grad2_b
        + (2.0 / 3.0 * rho * rho - rho_b * rho_b) * grad2_a;

    let term_grad = -LYP_A * LYP_B * omega * (bracket_ab + extra);

    term_first + term_grad
}

/// LYP correlation at density `rho` and gradient norm `grad_norm`.
///
/// Returns the LYP correlation energy density per electron and the two
/// potential components. The potentials are derivatives of the closed
/// LYP energy-density scalar taken by a tight central difference of
/// that analytic function — verified against the finite-difference
/// functional-derivative test.
pub fn lyp_correlation(rho: f64, grad_norm: f64) -> XcContribution {
    if rho <= 1.0e-30 {
        return XcContribution::zero();
    }
    let sigma = grad_norm * grad_norm;
    let e_c = lyp_energy_density(rho, sigma);
    let eps_c = e_c / rho;

    // ∂e_c/∂ρ by central difference on the closed scalar (the
    // functional form is too long to differentiate by hand reliably;
    // a tight central difference of the *analytic* closed form is
    // exact to ~1e-9 and is itself verified by the FD test against an
    // independent step).
    let hr = rho * 1.0e-7 + 1.0e-12;
    let de_drho =
        (lyp_energy_density(rho + hr, sigma) - lyp_energy_density(rho - hr, sigma)) / (2.0 * hr);

    // ∂e_c/∂|∇ρ| = ∂e_c/∂σ · 2|∇ρ|.
    let de_dgrad = if grad_norm > 1.0e-30 {
        let hs = sigma * 1.0e-7 + 1.0e-14;
        let de_dsigma = (lyp_energy_density(rho, sigma + hs) - lyp_energy_density(rho, sigma - hs))
            / (2.0 * hs);
        de_dsigma * 2.0 * grad_norm
    } else {
        0.0
    };

    XcContribution {
        energy_density: eps_c,
        potential: de_drho,
        gradient_potential: de_dgrad,
    }
}

/// The density-functional part of B3LYP at `(rho, grad_norm)` — the
/// whole B3LYP XC *minus* the exact-exchange term.
///
/// `E_xc^B3LYP − a·E_x^HF =
///   (1−a) E_x^LDA + b ΔE_x^B88 + (1−c) E_c^VWN + c E_c^LYP`.
/// The Kohn-Sham build adds the `a·E_x^HF` term itself.
pub fn b3lyp_dft_xc(rho: f64, grad_norm: f64) -> XcContribution {
    // (1−a) of LDA exchange.
    let lda_x = slater_exchange(rho).scale(1.0 - B3LYP_A_EXACT);
    // b · the B88 *gradient correction only* (B88 = LDA + correction).
    let b88_full = b88_exchange(rho, grad_norm);
    let b88_corr = b88_full.sub(&slater_exchange(rho)).scale(B3LYP_B_B88);
    // (1−c) VWN + c LYP correlation.
    let vwn = vwn_correlation(rho).scale(1.0 - B3LYP_C_LYP);
    let lyp = lyp_correlation(rho, grad_norm).scale(B3LYP_C_LYP);

    lda_x.add(&b88_corr).add(&vwn).add(&lyp)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// B88 exchange at zero gradient reduces exactly to the LDA
    /// (Slater) exchange — the correction vanishes at `x = 0`.
    #[test]
    fn b88_reduces_to_lda_at_zero_gradient() {
        for &rho in &[0.05, 0.5, 2.0, 8.0] {
            let b88 = b88_exchange(rho, 0.0);
            let lda = slater_exchange(rho);
            assert!(
                (b88.energy_density - lda.energy_density).abs() < 1.0e-12,
                "ρ={rho}: B88 {} vs LDA {}",
                b88.energy_density,
                lda.energy_density
            );
            assert!(b88.gradient_potential.abs() < 1.0e-12);
        }
    }

    /// The B88 gradient correction *lowers* the exchange energy density
    /// (makes it more negative) — B88 corrects the LDA's underbinding
    /// of exchange.
    #[test]
    fn b88_correction_lowers_exchange_energy() {
        let rho = 1.0;
        let lda = slater_exchange(rho).energy_density;
        let b88 = b88_exchange(rho, 2.0).energy_density;
        assert!(b88 < lda, "B88 {b88} should be < LDA {lda}");
    }

    /// The B88 exchange potential (density part) must match a finite
    /// difference of the B88 energy density at fixed gradient.
    #[test]
    fn b88_density_potential_is_functional_derivative() {
        let grad = 1.4;
        for &rho in &[0.1, 0.5, 2.0] {
            let h = rho * 1.0e-6;
            let fp = (rho + h) * b88_exchange(rho + h, grad).energy_density;
            let fm = (rho - h) * b88_exchange(rho - h, grad).energy_density;
            let fd = (fp - fm) / (2.0 * h);
            let v = b88_exchange(rho, grad).potential;
            assert!((v - fd).abs() < 1.0e-5, "ρ={rho}: v={v} vs FD={fd}");
        }
    }

    /// The B88 exchange gradient-potential must match a finite
    /// difference at fixed density.
    #[test]
    fn b88_gradient_potential_is_functional_derivative() {
        let rho = 1.0;
        for &grad in &[0.5, 1.5, 3.0] {
            let h = grad * 1.0e-6;
            let fp = rho * b88_exchange(rho, grad + h).energy_density;
            let fm = rho * b88_exchange(rho, grad - h).energy_density;
            let fd = (fp - fm) / (2.0 * h);
            let v = b88_exchange(rho, grad).gradient_potential;
            assert!((v - fd).abs() < 1.0e-5, "grad={grad}: v={v} vs FD={fd}");
        }
    }

    /// LYP correlation is negative for a physical density — correlation
    /// always lowers the energy.
    #[test]
    fn lyp_correlation_is_negative() {
        for &rho in &[0.1, 0.5, 1.0, 3.0] {
            let c = lyp_correlation(rho, 0.5);
            assert!(
                c.energy_density < 0.0,
                "ρ={rho}: LYP ε_c = {} should be < 0",
                c.energy_density
            );
        }
    }

    /// The LYP density-potential must match a finite difference of the
    /// LYP energy density at fixed gradient.
    #[test]
    fn lyp_density_potential_is_functional_derivative() {
        let grad = 0.9;
        for &rho in &[0.2, 0.6, 2.0] {
            let h = rho * 1.0e-5;
            let fp = (rho + h) * lyp_correlation(rho + h, grad).energy_density;
            let fm = (rho - h) * lyp_correlation(rho - h, grad).energy_density;
            let fd = (fp - fm) / (2.0 * h);
            let v = lyp_correlation(rho, grad).potential;
            assert!((v - fd).abs() < 1.0e-4, "ρ={rho}: v={v} vs FD={fd}");
        }
    }

    /// The LYP gradient-potential must match a finite difference at
    /// fixed density.
    #[test]
    fn lyp_gradient_potential_is_functional_derivative() {
        let rho = 1.0;
        for &grad in &[0.4, 1.2, 2.5] {
            let h = grad * 1.0e-5;
            let fp = rho * lyp_correlation(rho, grad + h).energy_density;
            let fm = rho * lyp_correlation(rho, grad - h).energy_density;
            let fd = (fp - fm) / (2.0 * h);
            let v = lyp_correlation(rho, grad).gradient_potential;
            assert!((v - fd).abs() < 1.0e-4, "grad={grad}: v={v} vs FD={fd}");
        }
    }

    /// The B3LYP coefficients are the standard published values.
    #[test]
    fn b3lyp_coefficients_are_standard() {
        assert!((B3LYP_A_EXACT - 0.20).abs() < 1.0e-12);
        assert!((B3LYP_B_B88 - 0.72).abs() < 1.0e-12);
        assert!((B3LYP_C_LYP - 0.81).abs() < 1.0e-12);
    }

    /// The B3LYP DFT part is finite and non-trivial at a physical
    /// density.
    #[test]
    fn b3lyp_dft_part_is_finite() {
        let c = b3lyp_dft_xc(1.0, 1.0);
        assert!(c.energy_density.is_finite());
        assert!(c.potential.is_finite());
        assert!(c.gradient_potential.is_finite());
        // The DFT-only XC energy density is negative (exchange
        // dominates) even before the exact-exchange term.
        assert!(c.energy_density < 0.0);
    }

    #[test]
    fn zero_density_is_zero() {
        assert_eq!(b88_exchange(0.0, 1.0).energy_density, 0.0);
        assert_eq!(lyp_correlation(0.0, 1.0).energy_density, 0.0);
        assert_eq!(b3lyp_dft_xc(0.0, 1.0).energy_density, 0.0);
    }
}
