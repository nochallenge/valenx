//! The local-density-approximation (LDA) exchange-correlation
//! functional.
//!
//! In the LDA the exchange-correlation energy density at a point
//! depends only on the local density `ρ(r)` — the value the uniform
//! electron gas would have at that density. This module ships the
//! standard spin-unpolarised LDA used for closed-shell molecules:
//!
//! - **Slater (Dirac) exchange** — the exact exchange energy of the
//!   uniform electron gas,
//!   `ε_x(ρ) = −C_x ρ^{1/3}`, `C_x = (3/4)(3/π)^{1/3}`.
//! - **VWN correlation** — the Vosko-Wilk-Nusair (1980) parametrisation
//!   of the correlation energy of the uniform electron gas, fitted to
//!   the Ceperley-Alder quantum Monte Carlo data. This module uses the
//!   **VWN5** paramagnetic functional (their formula 4.4 with the
//!   functional-V parameters), the form `"SVWN"` / `"LDA"` denotes in
//!   every quantum-chemistry program.
//!
//! ## What a functional returns
//!
//! Every functional here returns an [`XcContribution`] at a point: the
//! energy density per electron `ε_xc`, so the XC energy is
//! `E_xc = ∫ ρ ε_xc dr`, and the **potential** `v_xc = ∂(ρ ε_xc)/∂ρ`,
//! which enters the Kohn-Sham matrix. Splitting it this way lets the
//! Kohn-Sham build add `v_xc` to the Fock matrix directly.

use super::XcContribution;

/// The Slater exchange coefficient `C_x = (3/4)(3/π)^{1/3}`.
fn slater_cx() -> f64 {
    0.75 * (3.0 / std::f64::consts::PI).cbrt()
}

/// Slater (Dirac) exchange at density `rho`.
///
/// Returns the exchange energy density per electron `ε_x = −C_x ρ^{1/3}`
/// and the exchange potential `v_x = −(4/3) C_x ρ^{1/3}` (the
/// functional derivative `∂(ρ ε_x)/∂ρ`).
pub fn slater_exchange(rho: f64) -> XcContribution {
    if rho <= 0.0 {
        return XcContribution::zero();
    }
    let cx = slater_cx();
    let rho_third = rho.cbrt();
    let eps_x = -cx * rho_third;
    let v_x = -(4.0 / 3.0) * cx * rho_third;
    XcContribution {
        energy_density: eps_x,
        potential: v_x,
        // LDA: no density-gradient dependence.
        gradient_potential: 0.0,
    }
}

// --- VWN5 correlation -------------------------------------------------

/// VWN5 fit parameters for the paramagnetic (spin-unpolarised)
/// correlation energy.
struct VwnParams {
    a: f64,
    x0: f64,
    b: f64,
    c: f64,
}

/// The VWN5 paramagnetic parameter set (Vosko-Wilk-Nusair 1980,
/// their functional V). `A = 0.0310907` is the paramagnetic value used
/// by every quantum-chemistry program (Gaussian, NWChem, libxc) — half
/// the `0.0621814` quoted in some forms of the original paper.
const VWN5_PARA: VwnParams = VwnParams {
    a: 0.031_090_7,
    x0: -0.104_98,
    b: 3.727_44,
    c: 12.935_2,
};

/// The VWN correlation energy per electron `ε_c(x)` and its derivative
/// with respect to `x = √r_s`, for one parameter set.
///
/// The VWN functional form is
///
/// ```text
/// ε_c = A { ln(x²/X) + 2b/Q · atan(Q/(2x+b))
///           − b x0/X(x0) [ ln((x−x0)²/X) + 2(b+2x0)/Q · atan(Q/(2x+b)) ] }
/// ```
///
/// with `X(x) = x² + b x + c` and `Q = √(4c − b²)`.
fn vwn_eps_and_deriv(x: f64, p: &VwnParams) -> (f64, f64) {
    let q = (4.0 * p.c - p.b * p.b).sqrt();
    let xx = |y: f64| y * y + p.b * y + p.c;
    let big_x = xx(x);
    let big_x0 = xx(p.x0);
    let atan_arg = q / (2.0 * x + p.b);
    let atan = atan_arg.atan();

    let term1 = (x * x / big_x).ln();
    let term2 = 2.0 * p.b / q * atan;
    let term3 = -p.b * p.x0 / big_x0 * ((x - p.x0).powi(2) / big_x).ln();
    let term4 = -p.b * p.x0 / big_x0 * (2.0 * (p.b + 2.0 * p.x0) / q * atan);
    let eps = p.a * (term1 + term2 + term3 + term4);

    // dε_c/dx — analytic derivative of each term.
    // d/dx ln(x²/X) = 2/x − (2x+b)/X.
    let d1 = 2.0 / x - (2.0 * x + p.b) / big_x;
    // d/dx [2b/Q atan(Q/(2x+b))] = 2b/Q · (−2Q/((2x+b)²+Q²)).
    let denom = (2.0 * x + p.b).powi(2) + q * q;
    let d2 = 2.0 * p.b / q * (-2.0 * q / denom);
    // d/dx ln((x−x0)²/X) = 2/(x−x0) − (2x+b)/X.
    let d3_inner = 2.0 / (x - p.x0) - (2.0 * x + p.b) / big_x;
    let d3 = -p.b * p.x0 / big_x0 * d3_inner;
    let d4 = -p.b * p.x0 / big_x0 * (2.0 * (p.b + 2.0 * p.x0) / q * (-2.0 * q / denom));
    let deps_dx = p.a * (d1 + d2 + d3 + d4);
    (eps, deps_dx)
}

/// VWN5 correlation at density `rho` (spin-unpolarised).
///
/// Returns the correlation energy density per electron `ε_c` and the
/// correlation potential `v_c = ε_c − (r_s/3) dε_c/dr_s`, the standard
/// LDA functional derivative.
pub fn vwn_correlation(rho: f64) -> XcContribution {
    if rho <= 1.0e-30 {
        return XcContribution::zero();
    }
    // Wigner-Seitz radius r_s from ρ = 3/(4π r_s³).
    let rs = (3.0 / (4.0 * std::f64::consts::PI * rho)).cbrt();
    let x = rs.sqrt();
    let (eps_c, deps_dx) = vwn_eps_and_deriv(x, &VWN5_PARA);

    // v_c = ε_c − (r_s/3) dε_c/dr_s, and dε_c/dr_s = dε_c/dx · 1/(2x).
    let deps_drs = deps_dx / (2.0 * x);
    let v_c = eps_c - (rs / 3.0) * deps_drs;
    XcContribution {
        energy_density: eps_c,
        potential: v_c,
        gradient_potential: 0.0,
    }
}

/// The full LDA (SVWN) exchange-correlation contribution at `rho` —
/// Slater exchange plus VWN5 correlation.
pub fn lda_xc(rho: f64) -> XcContribution {
    slater_exchange(rho).add(&vwn_correlation(rho))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Slater exchange potential is exactly `4/3` of the exchange
    /// energy density — a defining LDA identity.
    #[test]
    fn slater_potential_is_four_thirds_energy() {
        for &rho in &[0.01, 0.1, 1.0, 5.0] {
            let c = slater_exchange(rho);
            let ratio = c.potential / c.energy_density;
            assert!(
                (ratio - 4.0 / 3.0).abs() < 1.0e-12,
                "ρ={rho}: v_x/ε_x = {ratio}"
            );
        }
    }

    /// The Slater exchange energy density at unit density is exactly
    /// `−C_x` with `C_x = (3/4)(3/π)^{1/3} ≈ 0.7385588`.
    #[test]
    fn slater_exchange_unit_density_reference() {
        let c = slater_exchange(1.0);
        let cx = 0.75 * (3.0 / std::f64::consts::PI).cbrt();
        assert!((c.energy_density - (-cx)).abs() < 1.0e-12);
        // Numerically C_x ≈ 0.738558766.
        assert!((cx - 0.738_558_766_4).abs() < 1.0e-8, "C_x = {cx}");
    }

    /// Zero (and negative) density returns a zero contribution — no
    /// `NaN` from `ρ^{1/3}` or `ln`.
    #[test]
    fn zero_density_is_zero_and_finite() {
        let x = lda_xc(0.0);
        assert_eq!(x.energy_density, 0.0);
        assert_eq!(x.potential, 0.0);
        let x = lda_xc(-1.0);
        assert!(x.energy_density.is_finite());
    }

    /// VWN5 correlation at the well-known density `r_s = 1` (i.e.
    /// `ρ = 3/4π`). The VWN5 paramagnetic correlation energy per
    /// electron at `r_s = 1` is the published uniform-electron-gas
    /// value `ε_c ≈ −0.0600` Hartree.
    #[test]
    fn vwn_correlation_at_rs_one_matches_published() {
        let rho = 3.0 / (4.0 * std::f64::consts::PI); // r_s = 1
        let c = vwn_correlation(rho);
        // VWN5 paramagnetic ε_c(r_s=1) ≈ −0.06002 Ha.
        assert!(
            (c.energy_density - (-0.060_02)).abs() < 5.0e-4,
            "ε_c(r_s=1) = {} (expected ≈ −0.0600)",
            c.energy_density
        );
        // Correlation energy is always negative.
        assert!(c.energy_density < 0.0);
    }

    /// VWN5 correlation at `r_s = 2` against the published value
    /// `ε_c ≈ −0.0448` Hartree.
    #[test]
    fn vwn_correlation_at_rs_two_matches_published() {
        let rho = 3.0 / (4.0 * std::f64::consts::PI * 8.0); // r_s = 2
        let c = vwn_correlation(rho);
        assert!(
            (c.energy_density - (-0.044_78)).abs() < 5.0e-4,
            "ε_c(r_s=2) = {} (expected ≈ −0.0448)",
            c.energy_density
        );
    }

    /// The VWN correlation potential `v_c` must be the functional
    /// derivative `∂(ρ ε_c)/∂ρ`. Check it against a finite difference
    /// of `ρ ε_c`.
    #[test]
    fn vwn_potential_is_functional_derivative() {
        for &rho in &[0.05, 0.2, 1.0, 3.0] {
            let h = rho * 1.0e-6;
            let fp = (rho + h) * vwn_correlation(rho + h).energy_density;
            let fm = (rho - h) * vwn_correlation(rho - h).energy_density;
            let fd = (fp - fm) / (2.0 * h);
            let v = vwn_correlation(rho).potential;
            assert!(
                (v - fd).abs() < 1.0e-6,
                "ρ={rho}: v_c={v} vs d(ρε_c)/dρ={fd}"
            );
        }
    }

    /// The Slater exchange potential must likewise be `∂(ρ ε_x)/∂ρ`.
    #[test]
    fn slater_potential_is_functional_derivative() {
        for &rho in &[0.05, 0.2, 1.0, 3.0] {
            let h = rho * 1.0e-6;
            let fp = (rho + h) * slater_exchange(rho + h).energy_density;
            let fm = (rho - h) * slater_exchange(rho - h).energy_density;
            let fd = (fp - fm) / (2.0 * h);
            let v = slater_exchange(rho).potential;
            assert!((v - fd).abs() < 1.0e-7, "ρ={rho}: v_x={v} vs FD={fd}");
        }
    }
}
