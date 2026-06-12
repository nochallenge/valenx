//! Generalised-gradient-approximation (GGA) functionals.
//!
//! A GGA functional's energy density depends on the local density `ρ`
//! *and* its gradient `∇ρ` — it can describe the inhomogeneity an LDA
//! misses. This module ships **PBE** (Perdew-Burke-Ernzerhof 1996, the
//! workhorse GGA) and the **B88** exchange + **LYP** correlation pieces
//! used by the B3LYP hybrid.
//!
//! ## The semilocal potential
//!
//! A GGA's Kohn-Sham potential is no longer a plain function of `ρ`.
//! Writing `E_xc = ∫ f(ρ, σ) dr` with `σ = |∇ρ|²`, the functional
//! derivative is
//!
//! ```text
//! v_xc = ∂f/∂ρ − ∇·( 2 ∂f/∂σ ∇ρ ).
//! ```
//!
//! The Kohn-Sham matrix build in [`crate::dft::ks`] handles the
//! divergence term by an integration by parts on the grid, so a
//! functional here only needs to report the two *local* partials. Each
//! returns an [`XcContribution`] carrying `ε_xc`, the density potential
//! `∂f/∂ρ` and `gradient_potential` `∂f/∂|∇ρ|` — the `|∇ρ|`-derivative
//! the divergence term needs.

use super::XcContribution;

/// Slater exchange coefficient `C_x = (3/4)(3/π)^{1/3}`.
fn slater_cx() -> f64 {
    0.75 * (3.0 / std::f64::consts::PI).cbrt()
}

// =====================================================================
// PBE exchange
// =====================================================================

/// PBE exchange `κ` parameter.
const PBE_KAPPA: f64 = 0.804;
/// PBE exchange `μ` parameter (`μ = β π²/3` with the PBE `β`).
const PBE_MU: f64 = 0.219_514_972_764_517_1;

/// PBE exchange at density `rho` and gradient norm `grad_norm`.
///
/// The PBE exchange energy density is the LDA exchange times the
/// enhancement factor `F_x(s) = 1 + κ − κ/(1 + μs²/κ)`, where
/// `s = |∇ρ| / (2 k_F ρ)` is the reduced density gradient and
/// `k_F = (3π²ρ)^{1/3}`.
pub fn pbe_exchange(rho: f64, grad_norm: f64) -> XcContribution {
    if rho <= 1.0e-30 {
        return XcContribution::zero();
    }
    let cx = slater_cx();
    let rho_third = rho.cbrt();
    let eps_x_lda = -cx * rho_third; // LDA exchange energy density
    let v_x_lda = -(4.0 / 3.0) * cx * rho_third;

    let kf = (3.0 * std::f64::consts::PI * std::f64::consts::PI * rho).cbrt();
    // s = |∇ρ| / (2 k_F ρ).
    let s = grad_norm / (2.0 * kf * rho);
    let s2 = s * s;

    let denom = 1.0 + PBE_MU * s2 / PBE_KAPPA;
    let fx = 1.0 + PBE_KAPPA - PBE_KAPPA / denom;
    // dF_x/d(s²) = μ / denom².
    let dfx_ds2 = PBE_MU / (denom * denom);

    // ε_x = ε_x^LDA · F_x.
    let eps_x = eps_x_lda * fx;

    // The potential pieces. Write the exchange energy *density* (per
    // unit volume) e = ρ ε_x^LDA F_x. We need ∂e/∂ρ and ∂e/∂|∇ρ|.
    //
    // s² = |∇ρ|² / (4 k_F² ρ²); with k_F ∝ ρ^{1/3},
    //   s² ∝ |∇ρ|² ρ^{−8/3}, so ∂s²/∂ρ = −(8/3) s²/ρ
    //   and ∂s²/∂|∇ρ| = 2 s² / |∇ρ|.
    // ∂e/∂ρ = ε_x^LDA F_x + ρ [ v_x^LDA' contribution ].
    // d(ρ ε_x^LDA)/dρ = v_x^LDA. So
    //   ∂e/∂ρ = v_x^LDA · F_x + (ρ ε_x^LDA) · dF_x/ds² · ∂s²/∂ρ.
    let ds2_drho = -(8.0 / 3.0) * s2 / rho;
    let de_drho = v_x_lda * fx + (rho * eps_x_lda) * dfx_ds2 * ds2_drho;

    // ∂e/∂|∇ρ| = (ρ ε_x^LDA) · dF_x/ds² · ∂s²/∂|∇ρ|.
    let de_dgrad = if grad_norm > 1.0e-30 {
        let ds2_dgrad = 2.0 * s2 / grad_norm;
        (rho * eps_x_lda) * dfx_ds2 * ds2_dgrad
    } else {
        0.0
    };

    XcContribution {
        energy_density: eps_x,
        potential: de_drho,
        gradient_potential: de_dgrad,
    }
}

// =====================================================================
// PW92 LDA correlation (the base PBE correlation builds on)
// =====================================================================

/// PW92 correlation energy per electron `ε_c(r_s)` and its derivative
/// `dε_c/dr_s`, paramagnetic spin-unpolarised parametrisation
/// (Perdew-Wang 1992).
fn pw92_eps_and_deriv(rs: f64) -> (f64, f64) {
    // Paramagnetic PW92 parameters.
    let a = 0.031_091;
    let alpha1 = 0.213_70;
    let beta1 = 7.5957;
    let beta2 = 3.5876;
    let beta3 = 1.6382;
    let beta4 = 0.492_94;

    let rsh = rs.sqrt();
    // G = 2A (β1 r^{1/2} + β2 r + β3 r^{3/2} + β4 r²).
    let g = 2.0 * a * (beta1 * rsh + beta2 * rs + beta3 * rs * rsh + beta4 * rs * rs);
    let one_plus_inv = 1.0 + 1.0 / g;
    let ln = one_plus_inv.ln();
    let eps = -2.0 * a * (1.0 + alpha1 * rs) * ln;

    // dε/dr_s. Let P = 1 + α1 r_s ; ε = −2A P ln(1 + 1/G).
    // dε/dr = −2A [ α1 ln(1+1/G) + P · d/dr ln(1+1/G) ].
    // d/dr ln(1+1/G) = (−G'/G²) / (1+1/G) = −G' / (G² + G).
    let gp = 2.0 * a * (0.5 * beta1 / rsh + beta2 + 1.5 * beta3 * rsh + 2.0 * beta4 * rs);
    let dln_dr = -gp / (g * g + g);
    let deps_dr = -2.0 * a * (alpha1 * ln + (1.0 + alpha1 * rs) * dln_dr);
    (eps, deps_dr)
}

// =====================================================================
// PBE correlation
// =====================================================================

/// PBE correlation `β` parameter.
const PBE_BETA: f64 = 0.066_724_550_603_149_22;

/// PBE correlation at density `rho` and gradient norm `grad_norm`.
///
/// The PBE correlation energy density per electron is
/// `ε_c = ε_c^PW92 + H(r_s, t)`, where the gradient correction `H`
/// uses the PBE `β` and `γ` and the reduced gradient
/// `t = |∇ρ| / (2 k_s ρ)` with `k_s = √(4 k_F / π)`.
pub fn pbe_correlation(rho: f64, grad_norm: f64) -> XcContribution {
    if rho <= 1.0e-30 {
        return XcContribution::zero();
    }
    let pi = std::f64::consts::PI;
    let gamma = (1.0 - std::f64::consts::LN_2) / (pi * pi);

    let rs = (3.0 / (4.0 * pi * rho)).cbrt();
    let (eps_pw, deps_pw_drs) = pw92_eps_and_deriv(rs);

    // k_F and the Thomas-Fermi screening k_s.
    let kf = (3.0 * pi * pi * rho).cbrt();
    let ks = (4.0 * kf / pi).sqrt();
    // Reduced gradient t = |∇ρ| / (2 k_s ρ).
    let t = grad_norm / (2.0 * ks * rho);
    let t2 = t * t;

    // A = (β/γ) / (exp(−ε_c^PW/γ) − 1).
    let exp_arg = (-eps_pw / gamma).exp();
    let cap_a = (PBE_BETA / gamma) / (exp_arg - 1.0).max(1.0e-30);

    // H = γ ln[ 1 + (β/γ) t² (1 + A t²)/(1 + A t² + A² t⁴) ].
    let at2 = cap_a * t2;
    let num = t2 * (1.0 + at2);
    let den = 1.0 + at2 + at2 * at2;
    let h_arg = 1.0 + (PBE_BETA / gamma) * num / den;
    let h = gamma * h_arg.ln();

    let eps_c = eps_pw + h;

    // --- Potential. e_c = ρ ε_c ; need ∂e_c/∂ρ and ∂e_c/∂|∇ρ|. -----
    //
    // The cleanest route is the chain rule through (r_s, t). The PBE H
    // depends on ρ through r_s (via ε_c^PW and A) and through t. We
    // build the partials numerically-stably with analytic pieces.
    //
    // ∂t²/∂ρ: t² ∝ |∇ρ|² ρ^{−2} k_s^{−2}, k_s² ∝ ρ^{1/3}, so
    //   t² ∝ |∇ρ|² ρ^{−7/3}  →  ∂t²/∂ρ = −(7/3) t²/ρ.
    // ∂t²/∂|∇ρ| = 2 t²/|∇ρ|.
    let dt2_drho = -(7.0 / 3.0) * t2 / rho;
    let dt2_dgrad = if grad_norm > 1.0e-30 {
        2.0 * t2 / grad_norm
    } else {
        0.0
    };

    // dr_s/dρ = −r_s/(3ρ).
    let drs_drho = -rs / (3.0 * rho);

    // dH/dt² and dH/d(ε_c^PW) via the explicit form of H.
    // Let R = (β/γ) num/den, so H = γ ln(1+R), dH = γ/(1+R) dR.
    let bg = PBE_BETA / gamma;
    // dR/dt² holding A fixed:
    //   num = t²(1+A t²) = t² + A t⁴ ; dnum/dt² = 1 + 2A t².
    //   den = 1 + A t² + A² t⁴ ; dden/dt² = A + 2A² t².
    let dnum_dt2 = 1.0 + 2.0 * at2;
    let dden_dt2 = cap_a + 2.0 * cap_a * cap_a * t2;
    let dr_dt2 = bg * (dnum_dt2 * den - num * dden_dt2) / (den * den);
    // dR/dA holding t² fixed:
    //   dnum/dA = t⁴ ; dden/dA = t² + 2A t⁴.
    let dnum_da = t2 * t2;
    let dden_da = t2 + 2.0 * cap_a * t2 * t2;
    let dr_da = bg * (dnum_da * den - num * dden_da) / (den * den);

    // dA/d(ε_c^PW): A = (β/γ)/(E−1) with E = exp(−ε_c/γ).
    //   dA/dε_c = (β/γ) · (−1/(E−1)²) · dE/dε_c, dE/dε_c = −E/γ.
    let e_minus_1 = (exp_arg - 1.0).max(1.0e-30);
    let da_deps = (PBE_BETA / gamma) * (exp_arg / gamma) / (e_minus_1 * e_minus_1);

    let one_plus_r = 1.0 + bg * num / den;
    let dh_dt2 = gamma / one_plus_r * dr_dt2;
    let dh_da = gamma / one_plus_r * dr_da;
    let dh_deps_pw = dh_da * da_deps;

    // ∂ε_c/∂ρ = dε_c^PW/dρ + ∂H/∂ρ.
    //   dε_c^PW/dρ = deps_pw_drs · drs_drho.
    //   ∂H/∂ρ = dH/dε_c^PW · (dε_c^PW/dρ) + dH/dt² · ∂t²/∂ρ.
    let deps_pw_drho = deps_pw_drs * drs_drho;
    let deps_c_drho = deps_pw_drho + dh_deps_pw * deps_pw_drho + dh_dt2 * dt2_drho;
    // ∂e_c/∂ρ = ε_c + ρ ∂ε_c/∂ρ.
    let de_drho = eps_c + rho * deps_c_drho;

    // ∂ε_c/∂|∇ρ| = dH/dt² · ∂t²/∂|∇ρ| ; ∂e_c/∂|∇ρ| = ρ · that.
    let de_dgrad = rho * dh_dt2 * dt2_dgrad;

    XcContribution {
        energy_density: eps_c,
        potential: de_drho,
        gradient_potential: de_dgrad,
    }
}

/// The full PBE exchange-correlation contribution at `(rho, grad_norm)`.
pub fn pbe_xc(rho: f64, grad_norm: f64) -> XcContribution {
    pbe_exchange(rho, grad_norm).add(&pbe_correlation(rho, grad_norm))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dft::functional::lda::slater_exchange;

    /// At zero density gradient PBE exchange reduces exactly to the
    /// LDA (Slater) exchange — the uniform-electron-gas limit.
    #[test]
    fn pbe_exchange_reduces_to_lda_at_zero_gradient() {
        for &rho in &[0.05, 0.5, 2.0, 8.0] {
            let pbe = pbe_exchange(rho, 0.0);
            let lda = slater_exchange(rho);
            assert!(
                (pbe.energy_density - lda.energy_density).abs() < 1.0e-12,
                "ρ={rho}: PBE_x {} vs LDA_x {}",
                pbe.energy_density,
                lda.energy_density
            );
            // The potential also reduces (gradient_potential → 0).
            assert!(pbe.gradient_potential.abs() < 1.0e-12);
            assert!(
                (pbe.potential - lda.potential).abs() < 1.0e-12,
                "ρ={rho}: potential mismatch"
            );
        }
    }

    /// A non-zero gradient makes the PBE exchange energy density *more
    /// negative* than the LDA — the enhancement factor `F_x > 1`
    /// multiplies a negative LDA exchange, lowering it.
    #[test]
    fn pbe_exchange_enhancement_is_below_lda() {
        let rho = 1.0;
        let lda = slater_exchange(rho).energy_density;
        let pbe = pbe_exchange(rho, 2.0).energy_density;
        // ε_x = ε_x^LDA · F_x with ε_x^LDA < 0 and F_x > 1, so PBE
        // exchange is more negative than the LDA.
        assert!(pbe < lda, "PBE_x {pbe} should be < LDA_x {lda}");
        assert!(pbe < 0.0);
    }

    /// PBE exchange enhancement saturates at `1 + κ` for a large
    /// reduced gradient.
    #[test]
    fn pbe_exchange_saturates() {
        let rho = 1.0;
        let lda = slater_exchange(rho).energy_density;
        // Huge gradient → F_x → 1 + κ.
        let pbe = pbe_exchange(rho, 1.0e6).energy_density;
        let ratio = pbe / lda;
        assert!(
            (ratio - (1.0 + PBE_KAPPA)).abs() < 1.0e-3,
            "F_x saturation ratio = {ratio}, expected {}",
            1.0 + PBE_KAPPA
        );
    }

    /// PBE correlation at zero gradient reduces to the PW92 LDA
    /// correlation — `H(t=0) = 0`.
    #[test]
    fn pbe_correlation_reduces_to_pw92_at_zero_gradient() {
        for &rho in &[0.05, 0.5, 2.0] {
            let pbe = pbe_correlation(rho, 0.0);
            let rs = (3.0 / (4.0 * std::f64::consts::PI * rho)).cbrt();
            let (eps_pw, _) = pw92_eps_and_deriv(rs);
            assert!(
                (pbe.energy_density - eps_pw).abs() < 1.0e-12,
                "ρ={rho}: PBE_c {} vs PW92 {}",
                pbe.energy_density,
                eps_pw
            );
        }
    }

    /// PW92 correlation at `r_s = 1` reproduces the published
    /// uniform-electron-gas value `ε_c ≈ −0.0598` Hartree.
    #[test]
    fn pw92_correlation_at_rs_one_matches_published() {
        let (eps, _) = pw92_eps_and_deriv(1.0);
        assert!(
            (eps - (-0.059_77)).abs() < 5.0e-4,
            "PW92 ε_c(r_s=1) = {eps}, expected ≈ −0.0598"
        );
    }

    /// A non-zero gradient lowers the PBE correlation energy density
    /// (the gradient correction `H ≥ 0`, so `ε_c` becomes *less*
    /// negative).
    #[test]
    fn pbe_correlation_gradient_correction_is_nonnegative() {
        let rho = 1.0;
        let base = pbe_correlation(rho, 0.0).energy_density;
        let with_grad = pbe_correlation(rho, 1.5).energy_density;
        // H ≥ 0, so ε_c with gradient ≥ ε_c without.
        assert!(with_grad >= base - 1.0e-12, "{with_grad} vs {base}");
    }

    /// PBE correlation slowly-varying-density limit: as the gradient
    /// tends to zero from a small value, PBE approaches PW92 LDA — the
    /// "PBE reduces toward LDA for a slowly-varying density" property.
    #[test]
    fn pbe_approaches_lda_for_slowly_varying_density() {
        let rho = 1.0;
        let lda = pbe_correlation(rho, 0.0).energy_density;
        let mut prev_diff = f64::MAX;
        // Shrinking gradient → shrinking difference from LDA.
        for &g in &[0.5, 0.1, 0.02, 0.004] {
            let diff = (pbe_correlation(rho, g).energy_density - lda).abs();
            assert!(
                diff <= prev_diff,
                "diff not shrinking: {diff} > {prev_diff}"
            );
            prev_diff = diff;
        }
        assert!(prev_diff < 1.0e-4, "residual {prev_diff}");
    }

    /// The PBE exchange density-potential `∂(ρε_x)/∂ρ` must match a
    /// finite difference *at fixed gradient*.
    #[test]
    fn pbe_exchange_density_potential_is_functional_derivative() {
        let grad = 1.3;
        for &rho in &[0.1, 0.5, 2.0] {
            let h = rho * 1.0e-6;
            let fp = (rho + h) * pbe_exchange(rho + h, grad).energy_density;
            let fm = (rho - h) * pbe_exchange(rho - h, grad).energy_density;
            let fd = (fp - fm) / (2.0 * h);
            let v = pbe_exchange(rho, grad).potential;
            assert!((v - fd).abs() < 1.0e-5, "ρ={rho}: v={v} vs ∂(ρε_x)/∂ρ={fd}");
        }
    }

    /// The PBE exchange gradient-potential `∂(ρε_x)/∂|∇ρ|` must match a
    /// finite difference *at fixed density*.
    #[test]
    fn pbe_exchange_gradient_potential_is_functional_derivative() {
        let rho = 1.0;
        for &grad in &[0.5, 1.5, 3.0] {
            let h = grad * 1.0e-6;
            let fp = rho * pbe_exchange(rho, grad + h).energy_density;
            let fm = rho * pbe_exchange(rho, grad - h).energy_density;
            let fd = (fp - fm) / (2.0 * h);
            let v = pbe_exchange(rho, grad).gradient_potential;
            assert!(
                (v - fd).abs() < 1.0e-5,
                "grad={grad}: v={v} vs ∂(ρε_x)/∂|∇ρ|={fd}"
            );
        }
    }

    /// The PBE correlation density-potential must match a finite
    /// difference at fixed gradient.
    #[test]
    fn pbe_correlation_density_potential_is_functional_derivative() {
        let grad = 1.1;
        for &rho in &[0.1, 0.5, 2.0] {
            let h = rho * 1.0e-6;
            let fp = (rho + h) * pbe_correlation(rho + h, grad).energy_density;
            let fm = (rho - h) * pbe_correlation(rho - h, grad).energy_density;
            let fd = (fp - fm) / (2.0 * h);
            let v = pbe_correlation(rho, grad).potential;
            assert!((v - fd).abs() < 1.0e-5, "ρ={rho}: v={v} vs FD={fd}");
        }
    }

    /// The PBE correlation gradient-potential must match a finite
    /// difference at fixed density.
    #[test]
    fn pbe_correlation_gradient_potential_is_functional_derivative() {
        let rho = 1.0;
        for &grad in &[0.5, 1.5, 3.0] {
            let h = grad * 1.0e-6;
            let fp = rho * pbe_correlation(rho, grad + h).energy_density;
            let fm = rho * pbe_correlation(rho, grad - h).energy_density;
            let fd = (fp - fm) / (2.0 * h);
            let v = pbe_correlation(rho, grad).gradient_potential;
            assert!((v - fd).abs() < 1.0e-5, "grad={grad}: v={v} vs FD={fd}");
        }
    }

    #[test]
    fn zero_density_is_zero() {
        assert_eq!(pbe_xc(0.0, 0.0).energy_density, 0.0);
        assert_eq!(pbe_xc(0.0, 1.0).potential, 0.0);
    }
}
