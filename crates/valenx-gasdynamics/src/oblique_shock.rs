//! Oblique-shock relations via the θ–β–M equation.
//!
//! An attached oblique shock turns supersonic flow through a deflection angle
//! `θ` at a shock angle `β` (measured from the upstream flow). The three are
//! tied by the **θ–β–M relation**:
//!
//! ```text
//! tan θ = 2 cot β · (M₁² sin²β − 1) / (M₁² (γ + cos 2β) + 2)
//! ```
//!
//! - [`deflection_angle`] is the forward direction: `(M₁, β) → θ`.
//! - For a given `(M₁, θ)` there are **two** shock angles — a *weak* (small `β`,
//!   usually leaving the flow supersonic) and a *strong* (large `β`) solution —
//!   returned by [`shock_angle`]. They coincide at the maximum deflection
//!   [`max_deflection_angle`]; beyond it no attached shock exists and the shock
//!   detaches.
//! - [`oblique_shock`] bundles the downstream state by applying the normal-shock
//!   jump to the shock-normal Mach component `M₁ sin β`.
//!
//! The shock-normal component must be supersonic (`M₁ sin β ≥ 1`) for a shock to
//! exist. Reference: Anderson, *Modern Compressible Flow*; NACA Report 1135.
//! Same perfect-gas scope and caveats as the rest of the crate.

use crate::error::{check_gamma, GasError, Result};
use crate::normal_shock::normal_shock;

/// Raw θ–β–M relation (radians), without input validation.
fn theta_of_beta(m1: f64, beta: f64, gamma: f64) -> f64 {
    let (sin_b, cos_b) = beta.sin_cos();
    let num = m1 * m1 * sin_b * sin_b - 1.0;
    let den = m1 * m1 * (gamma + (2.0 * beta).cos()) + 2.0;
    (2.0 * (cos_b / sin_b) * num / den).atan()
}

fn check_supersonic(m1: f64) -> Result<()> {
    if !m1.is_finite() || m1 <= 1.0 {
        return Err(GasError::bad_mach(
            m1,
            "oblique_shock",
            "must be finite and > 1 (supersonic)",
        ));
    }
    Ok(())
}

/// Maximum deflection and the shock angle at which it occurs, by a coarse scan
/// over `β ∈ (μ, π/2)` followed by a local refinement (`μ` is the Mach angle).
fn max_deflection_internal(m1: f64, gamma: f64) -> (f64, f64) {
    let mu = (1.0 / m1).asin();
    let hi = std::f64::consts::FRAC_PI_2;
    let scan = |lo: f64, hi: f64, n: usize| {
        let mut best_theta = f64::NEG_INFINITY;
        let mut best_beta = lo;
        for i in 0..=n {
            let beta = lo + (hi - lo) * (i as f64 / n as f64);
            let theta = theta_of_beta(m1, beta, gamma);
            if theta > best_theta {
                best_theta = theta;
                best_beta = beta;
            }
        }
        (best_theta, best_beta)
    };
    let (_, b0) = scan(mu, hi, 2000);
    let span = (hi - mu) / 2000.0;
    scan((b0 - span).max(mu), (b0 + span).min(hi), 2000)
}

/// The flow deflection angle `θ` (radians) produced by an oblique shock at shock
/// angle `beta` (radians) in a Mach `m1` flow.
///
/// # Errors
///
/// [`GasError::BadGamma`] for invalid `gamma`; [`GasError::BadMach`] if `m1` is
/// not supersonic, if `beta` is outside `(0, π/2]`, or if the shock-normal Mach
/// component `m1 sin β` is below 1 (no shock).
pub fn deflection_angle(m1: f64, beta: f64, gamma: f64) -> Result<f64> {
    let gamma = check_gamma(gamma)?;
    check_supersonic(m1)?;
    if !beta.is_finite() || beta <= 0.0 || beta > std::f64::consts::FRAC_PI_2 {
        return Err(GasError::bad_mach(
            beta,
            "shock_angle",
            "beta must be finite and in (0, pi/2]",
        ));
    }
    let mn1 = m1 * beta.sin();
    if mn1 < 1.0 {
        return Err(GasError::bad_mach(
            mn1,
            "oblique_shock",
            "normal Mach component m1*sin(beta) must be >= 1",
        ));
    }
    Ok(theta_of_beta(m1, beta, gamma))
}

/// The maximum flow deflection (radians) an attached oblique shock can produce
/// at Mach `m1`. Beyond this the shock detaches.
///
/// # Errors
///
/// [`GasError::BadGamma`] / [`GasError::BadMach`] for invalid inputs.
pub fn max_deflection_angle(m1: f64, gamma: f64) -> Result<f64> {
    let gamma = check_gamma(gamma)?;
    check_supersonic(m1)?;
    Ok(max_deflection_internal(m1, gamma).0)
}

/// The shock angle `β` (radians) for an oblique shock that turns Mach `m1` flow
/// through deflection `theta`. `strong = false` returns the weak solution (the
/// physically common one), `true` the strong solution.
///
/// # Errors
///
/// [`GasError::BadGamma`] / [`GasError::BadMach`] for invalid inputs, or if
/// `theta` is not in `(0, θ_max)` — at or beyond the maximum deflection the
/// shock detaches and no attached solution exists.
pub fn shock_angle(m1: f64, theta: f64, gamma: f64, strong: bool) -> Result<f64> {
    let gamma = check_gamma(gamma)?;
    check_supersonic(m1)?;
    let (theta_max, beta_max) = max_deflection_internal(m1, gamma);
    if !theta.is_finite() || theta <= 0.0 || theta >= theta_max {
        return Err(GasError::bad_mach(
            theta,
            "oblique_shock_deflection",
            "theta must be in (0, theta_max); beyond it the shock detaches",
        ));
    }
    let mu = (1.0 / m1).asin();
    // θ(β) rises 0→θ_max on (μ, β_max) [weak] and falls θ_max→0 on (β_max, π/2)
    // [strong]; bisect the appropriate monotone branch.
    let (mut lo, mut hi) = if strong {
        (beta_max, std::f64::consts::FRAC_PI_2)
    } else {
        (mu, beta_max)
    };
    for _ in 0..80 {
        let mid = 0.5 * (lo + hi);
        let t = theta_of_beta(m1, mid, gamma);
        // On the weak branch θ increases with β; on the strong branch it
        // decreases. `below` tracks whether we are short of the target.
        let below = if strong { t > theta } else { t < theta };
        if below {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    Ok(0.5 * (lo + hi))
}

/// The complete oblique-shock state for a `(m1, theta)` turn.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ObliqueShock {
    /// Shock angle `β` (radians).
    pub shock_angle: f64,
    /// Flow deflection `θ` (radians).
    pub deflection: f64,
    /// Upstream shock-normal Mach component `M₁ sin β` (`≥ 1`).
    pub upstream_normal_mach: f64,
    /// Downstream Mach number `M₂`.
    pub downstream_mach: f64,
    /// Static-pressure ratio `p₂/p₁` (`≥ 1`).
    pub pressure_ratio: f64,
    /// Static-density ratio `ρ₂/ρ₁` (`≥ 1`).
    pub density_ratio: f64,
    /// Static-temperature ratio `T₂/T₁` (`≥ 1`).
    pub temperature_ratio: f64,
}

/// Solve the oblique shock that turns Mach `m1` flow through `theta` (radians),
/// taking the weak (`strong = false`) or strong solution, and apply the
/// normal-shock jump to the shock-normal component.
///
/// # Errors
///
/// As [`shock_angle`]; also propagates the [`normal_shock`] error for the
/// shock-normal component.
pub fn oblique_shock(m1: f64, theta: f64, gamma: f64, strong: bool) -> Result<ObliqueShock> {
    let beta = shock_angle(m1, theta, gamma, strong)?;
    let mn1 = m1 * beta.sin();
    let ns = normal_shock(mn1, gamma)?;
    let m2 = ns.downstream_mach / (beta - theta).sin();
    Ok(ObliqueShock {
        shock_angle: beta,
        deflection: theta,
        upstream_normal_mach: mn1,
        downstream_mach: m2,
        pressure_ratio: ns.pressure_ratio,
        density_ratio: ns.density_ratio,
        temperature_ratio: ns.temperature_ratio,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    const G: f64 = 1.4;

    #[test]
    fn deflection_matches_naca_1135() {
        // NACA 1135 θ–β–M chart: M=2, β=40° → θ ≈ 10.6°.
        let theta = deflection_angle(2.0, 40.0_f64.to_radians(), G)
            .unwrap()
            .to_degrees();
        assert!((theta - 10.6).abs() < 0.1, "θ = {theta:.3}° vs ~10.6°");
    }

    #[test]
    fn limits_give_zero_deflection() {
        // A normal shock (β = 90°) and a Mach wave (β = μ = asin(1/M)) both turn
        // the flow through zero angle.
        assert!(
            deflection_angle(2.0, std::f64::consts::FRAC_PI_2, G)
                .unwrap()
                .abs()
                < 1e-9
        );
        let mu = (1.0_f64 / 2.0).asin();
        assert!(deflection_angle(2.0, mu, G).unwrap().abs() < 1e-9);
    }

    #[test]
    fn max_deflection_for_mach_2() {
        // NACA 1135: at M=2 the maximum deflection is ≈ 22.97°.
        let theta_max = max_deflection_angle(2.0, G).unwrap().to_degrees();
        assert!(
            (theta_max - 22.97).abs() < 0.1,
            "θ_max = {theta_max:.3}° vs ~22.97°"
        );
    }

    #[test]
    fn weak_and_strong_roots_round_trip() {
        let theta = 15.0_f64.to_radians();
        let weak = shock_angle(2.0, theta, G, false).unwrap();
        let strong = shock_angle(2.0, theta, G, true).unwrap();
        assert!(weak < strong, "weak {weak} should be < strong {strong}");
        // Both reproduce the requested deflection.
        assert!((deflection_angle(2.0, weak, G).unwrap() - theta).abs() < 1e-6);
        assert!((deflection_angle(2.0, strong, G).unwrap() - theta).abs() < 1e-6);
        // The weak root for ~15° at M=2 is ≈ 45.3° (NACA chart).
        assert!(
            (weak.to_degrees() - 45.3).abs() < 0.5,
            "weak β = {:.2}°",
            weak.to_degrees()
        );
    }

    #[test]
    fn detached_deflection_is_an_error() {
        // 25° > θ_max(22.97°) at M=2 → no attached shock.
        assert!(shock_angle(2.0, 25.0_f64.to_radians(), G, false).is_err());
        assert!(deflection_angle(0.8, 40.0_f64.to_radians(), G).is_err()); // subsonic
    }

    #[test]
    fn weak_oblique_shock_keeps_flow_supersonic() {
        // A weak oblique shock at M=2, 10° deflection: downstream still
        // supersonic, with a pressure rise.
        let s = oblique_shock(2.0, 10.0_f64.to_radians(), G, false).unwrap();
        assert!(s.upstream_normal_mach >= 1.0);
        assert!(
            s.downstream_mach < 2.0 && s.downstream_mach > 1.0,
            "M2 = {}",
            s.downstream_mach
        );
        assert!(s.pressure_ratio > 1.0 && s.density_ratio > 1.0 && s.temperature_ratio > 1.0);
        // The strong solution drops the flow subsonic.
        let strong = oblique_shock(2.0, 10.0_f64.to_radians(), G, true).unwrap();
        assert!(
            strong.downstream_mach < 1.0,
            "strong M2 = {}",
            strong.downstream_mach
        );
    }
}
