//! Stationary normal-shock jump relations for a calorically-perfect
//! ideal gas.
//!
//! A normal shock is a thin, stationary discontinuity standing
//! perpendicular to a supersonic stream. Across it the flow jumps from a
//! supersonic upstream state `1` to a subsonic downstream state `2`. For
//! a perfect gas the jumps are exact algebraic functions of the upstream
//! Mach number `M1` and `gamma` (Anderson, *Modern Compressible Flow*,
//! ch. 3; NACA Report 1135):
//!
//! ```text
//! M2^2     = ( 1 + (gamma-1)/2 * M1^2 ) / ( gamma * M1^2 - (gamma-1)/2 )
//! p2 / p1  = 1 + 2 gamma / (gamma + 1) * (M1^2 - 1)
//! rho2/rho1= (gamma + 1) M1^2 / ( (gamma - 1) M1^2 + 2 )
//! T2 / T1  = (p2/p1) * (rho1/rho2)                     [perfect-gas EOS]
//! p02/p01  = [ (gamma+1)M1^2 / ((gamma-1)M1^2 + 2) ] ^ ( gamma/(gamma-1) )
//!            * [ (gamma+1) / (2 gamma M1^2 - (gamma-1)) ] ^ ( 1/(gamma-1) )
//! ```
//!
//! The second law forces `M1 >= 1`: a stationary shock can only sit in a
//! supersonic flow, and the entropy rise makes the stagnation-pressure
//! ratio `p02/p01 <= 1`. At exactly `M1 = 1` every ratio is unity (the
//! vanishingly weak shock).
//!
//! ## Honest scope
//!
//! These relations are the textbook *normal*-shock jumps for a single
//! perfect gas with constant specific heats. They model neither oblique
//! shocks, shock thickness / structure, real-gas chemistry, nor wall
//! effects. They are for study and first-order estimates, not certified
//! design.

use serde::{Deserialize, Serialize};

use crate::error::{check_gamma, GasError, Result};

/// The complete set of normal-shock property ratios across a stationary
/// shock, as a function of the upstream Mach number and `gamma`.
///
/// Build it with [`normal_shock`]. Every ratio is downstream `2` over
/// upstream `1`, except [`downstream_mach`](Self::downstream_mach) which
/// is the absolute Mach number behind the shock.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct NormalShock {
    /// Upstream (supersonic) Mach number `M1` that produced this shock.
    pub upstream_mach: f64,
    /// Downstream (subsonic) Mach number `M2`.
    pub downstream_mach: f64,
    /// Static-pressure ratio `p2 / p1` (`>= 1`).
    pub pressure_ratio: f64,
    /// Static-temperature ratio `T2 / T1` (`>= 1`).
    pub temperature_ratio: f64,
    /// Static-density ratio `rho2 / rho1` (`>= 1`).
    pub density_ratio: f64,
    /// Stagnation-pressure ratio `p02 / p01` (`<= 1`; the entropy-rise
    /// signature of the shock).
    pub stagnation_pressure_ratio: f64,
}

/// Validate the upstream Mach number for a normal-shock evaluation: it
/// must be finite and `>= 1`. Returns it unchanged on success.
fn check_upstream_mach(m1: f64) -> Result<f64> {
    if !m1.is_finite() || m1 < 1.0 {
        return Err(GasError::subsonic_shock(m1));
    }
    Ok(m1)
}

/// Downstream Mach number `M2` behind a stationary normal shock.
///
/// ```text
/// M2 = sqrt( ( 1 + (gamma-1)/2 * M1^2 ) / ( gamma * M1^2 - (gamma-1)/2 ) )
/// ```
///
/// Subsonic for every `M1 > 1`; exactly one at `M1 = 1`.
///
/// # Errors
///
/// Returns [`GasError::SubsonicShock`] if `M1` is non-finite or `< 1`,
/// and [`GasError::BadGamma`] if `gamma` is non-finite or `<= 1`.
pub fn downstream_mach(m1: f64, gamma: f64) -> Result<f64> {
    let g = check_gamma(gamma)?;
    let m1 = check_upstream_mach(m1)?;
    let num = 1.0 + (g - 1.0) / 2.0 * m1 * m1;
    let den = g * m1 * m1 - (g - 1.0) / 2.0;
    Ok((num / den).sqrt())
}

/// Static-pressure ratio `p2 / p1 = 1 + 2 gamma/(gamma+1) * (M1^2 - 1)`.
///
/// `>= 1`, increasing without bound as `M1` grows; equal to one at
/// `M1 = 1`.
///
/// # Errors
///
/// As [`downstream_mach`].
pub fn pressure_ratio(m1: f64, gamma: f64) -> Result<f64> {
    let g = check_gamma(gamma)?;
    let m1 = check_upstream_mach(m1)?;
    Ok(1.0 + 2.0 * g / (g + 1.0) * (m1 * m1 - 1.0))
}

/// Static-density ratio
/// `rho2 / rho1 = (gamma+1) M1^2 / ((gamma-1) M1^2 + 2)`.
///
/// `>= 1`, approaching the finite limit `(gamma+1)/(gamma-1)` as
/// `M1 -> infinity`; equal to one at `M1 = 1`.
///
/// # Errors
///
/// As [`downstream_mach`].
pub fn density_ratio(m1: f64, gamma: f64) -> Result<f64> {
    let g = check_gamma(gamma)?;
    let m1 = check_upstream_mach(m1)?;
    let m2 = m1 * m1;
    Ok((g + 1.0) * m2 / ((g - 1.0) * m2 + 2.0))
}

/// Static-temperature ratio `T2 / T1`.
///
/// From the perfect-gas equation of state the temperature ratio is the
/// pressure ratio divided by the density ratio,
/// `T2/T1 = (p2/p1) / (rho2/rho1)`. `>= 1`, increasing without bound with
/// `M1`; equal to one at `M1 = 1`.
///
/// # Errors
///
/// As [`downstream_mach`].
pub fn temperature_ratio(m1: f64, gamma: f64) -> Result<f64> {
    let g = check_gamma(gamma)?;
    let p = pressure_ratio(m1, g)?;
    let rho = density_ratio(m1, g)?;
    Ok(p / rho)
}

/// Stagnation-pressure ratio `p02 / p01` across the shock.
///
/// ```text
/// p02/p01 = [ (gamma+1)M1^2 / ((gamma-1)M1^2 + 2) ] ^ ( gamma/(gamma-1) )
///           * [ (gamma+1) / (2 gamma M1^2 - (gamma-1)) ] ^ ( 1/(gamma-1) )
/// ```
///
/// Always `<= 1` — the irreversible entropy rise destroys stagnation
/// pressure — and equal to one at `M1 = 1`.
///
/// # Errors
///
/// As [`downstream_mach`].
pub fn stagnation_pressure_ratio(m1: f64, gamma: f64) -> Result<f64> {
    let g = check_gamma(gamma)?;
    let m1 = check_upstream_mach(m1)?;
    let m2 = m1 * m1;
    let term1 = ((g + 1.0) * m2 / ((g - 1.0) * m2 + 2.0)).powf(g / (g - 1.0));
    let term2 = ((g + 1.0) / (2.0 * g * m2 - (g - 1.0))).powf(1.0 / (g - 1.0));
    Ok(term1 * term2)
}

/// All normal-shock property ratios at once.
///
/// Convenience wrapper that returns the full [`NormalShock`] bundle. The
/// density ratio is computed once and reused for the temperature ratio.
///
/// # Errors
///
/// As [`downstream_mach`].
pub fn normal_shock(m1: f64, gamma: f64) -> Result<NormalShock> {
    let g = check_gamma(gamma)?;
    let m1 = check_upstream_mach(m1)?;
    let p = pressure_ratio(m1, g)?;
    let rho = density_ratio(m1, g)?;
    Ok(NormalShock {
        upstream_mach: m1,
        downstream_mach: downstream_mach(m1, g)?,
        pressure_ratio: p,
        temperature_ratio: p / rho,
        density_ratio: rho,
        stagnation_pressure_ratio: stagnation_pressure_ratio(m1, g)?,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Absolute tolerance for analytic comparisons (closed-form, so the
    /// only error is f64 round-off).
    const EPS: f64 = 1e-9;
    /// Looser bound for the four-significant-figure NACA-1135 table
    /// entries (`M2`, `p02/p01`).
    const TABLE_EPS: f64 = 5e-5;

    fn close(a: f64, b: f64, eps: f64) -> bool {
        (a - b).abs() < eps
    }

    #[test]
    fn all_ratios_are_unity_at_mach_one() {
        let s = normal_shock(1.0, 1.4).unwrap();
        assert!(
            close(s.downstream_mach, 1.0, EPS),
            "M2 = {}",
            s.downstream_mach
        );
        assert!(
            close(s.pressure_ratio, 1.0, EPS),
            "p2/p1 = {}",
            s.pressure_ratio
        );
        assert!(
            close(s.temperature_ratio, 1.0, EPS),
            "T2/T1 = {}",
            s.temperature_ratio
        );
        assert!(
            close(s.density_ratio, 1.0, EPS),
            "rho2/rho1 = {}",
            s.density_ratio
        );
        assert!(
            close(s.stagnation_pressure_ratio, 1.0, EPS),
            "p02/p01 = {}",
            s.stagnation_pressure_ratio
        );
    }

    #[test]
    fn qualitative_shock_signs_for_supersonic_upstream() {
        // For every M1 > 1: M2 < 1, p2 > p1, T2 > T1, rho2 > rho1, and
        // the stagnation pressure is destroyed (p02/p01 < 1).
        for &m1 in &[1.5_f64, 2.0, 3.0, 5.0] {
            let s = normal_shock(m1, 1.4).unwrap();
            assert!(
                s.downstream_mach < 1.0,
                "M1={m1}: M2 = {}",
                s.downstream_mach
            );
            assert!(
                s.pressure_ratio > 1.0,
                "M1={m1}: p2/p1 = {}",
                s.pressure_ratio
            );
            assert!(
                s.temperature_ratio > 1.0,
                "M1={m1}: T2/T1 = {}",
                s.temperature_ratio
            );
            assert!(
                s.density_ratio > 1.0,
                "M1={m1}: rho2/rho1 = {}",
                s.density_ratio
            );
            assert!(
                s.stagnation_pressure_ratio < 1.0,
                "M1={m1}: p02/p01 = {}",
                s.stagnation_pressure_ratio
            );
        }
    }

    #[test]
    fn known_values_m1_2_gamma14() {
        // Ground truth (NACA Report 1135, gamma = 1.4, M1 = 2.0):
        //   M2       = 0.57735   (= 1/sqrt(3), exact)
        //   p2/p1    = 4.5       (exact)
        //   rho2/rho1= 8/3 = 2.66667 (exact)
        //   T2/T1    = 1.6875    (exact)
        //   p02/p01  = 0.72087   (table, 5 sig figs)
        let s = normal_shock(2.0, 1.4).unwrap();
        assert!(
            close(s.downstream_mach, 1.0 / 3.0_f64.sqrt(), EPS),
            "M2 = {}",
            s.downstream_mach
        );
        assert!(
            close(s.pressure_ratio, 4.5, EPS),
            "p2/p1 = {}",
            s.pressure_ratio
        );
        assert!(
            close(s.density_ratio, 8.0 / 3.0, EPS),
            "rho2/rho1 = {}",
            s.density_ratio
        );
        assert!(
            close(s.temperature_ratio, 1.6875, EPS),
            "T2/T1 = {}",
            s.temperature_ratio
        );
        assert!(
            close(s.stagnation_pressure_ratio, 0.720_874, TABLE_EPS),
            "p02/p01 = {}",
            s.stagnation_pressure_ratio
        );
    }

    #[test]
    fn known_values_m1_3_gamma14() {
        // Second independent table point (NACA 1135, gamma = 1.4,
        // M1 = 3.0): M2 = 0.475191, p2/p1 = 10.33333, T2/T1 = 2.679012,
        // rho2/rho1 = 3.857143, p02/p01 = 0.328344.
        let s = normal_shock(3.0, 1.4).unwrap();
        assert!(
            close(s.pressure_ratio, 31.0 / 3.0, EPS),
            "p2/p1 = {}",
            s.pressure_ratio
        );
        assert!(
            close(s.density_ratio, 27.0 / 7.0, EPS),
            "rho2/rho1 = {}",
            s.density_ratio
        );
        assert!(
            close(s.downstream_mach, 0.475_191_792_359_017, 1e-6),
            "M2 = {}",
            s.downstream_mach
        );
        assert!(
            close(s.temperature_ratio, 2.679_012_345_679_012, 1e-9),
            "T2/T1 = {}",
            s.temperature_ratio
        );
        assert!(
            close(s.stagnation_pressure_ratio, 0.328_344, TABLE_EPS),
            "p02/p01 = {}",
            s.stagnation_pressure_ratio
        );
    }

    #[test]
    fn temperature_ratio_equals_pressure_over_density() {
        // The EOS identity must hold at an arbitrary off-table state.
        let m1 = 2.73;
        let g = 1.3;
        let t = temperature_ratio(m1, g).unwrap();
        let p = pressure_ratio(m1, g).unwrap();
        let rho = density_ratio(m1, g).unwrap();
        assert!(close(t, p / rho, EPS), "T={t}, p/rho={}", p / rho);
    }

    #[test]
    fn density_ratio_saturates_at_strong_shock_limit() {
        // As M1 -> infinity, rho2/rho1 -> (gamma+1)/(gamma-1) = 6 for
        // gamma = 1.4. A very strong shock should sit just under it.
        let limit = (1.4 + 1.0) / (1.4 - 1.0); // = 6.0
        let s = density_ratio(1.0e4, 1.4).unwrap();
        assert!(s < limit, "rho2/rho1 = {s} should stay below {limit}");
        assert!(
            close(s, limit, 1e-3),
            "rho2/rho1 = {s} should approach {limit}"
        );
    }

    #[test]
    fn rankine_hugoniot_entropy_is_nonnegative() {
        // The non-dimensional entropy change across the shock,
        // (s2 - s1)/R = (gamma/(gamma-1)) ln(T2/T1) - ln(p2/p1), must be
        // >= 0 for every admissible shock (second law). It is exactly
        // zero at M1 = 1 and positive beyond.
        let g = 1.4;
        for &m1 in &[1.0_f64, 1.2, 2.0, 4.0] {
            let t = temperature_ratio(m1, g).unwrap();
            let p = pressure_ratio(m1, g).unwrap();
            let ds = g / (g - 1.0) * t.ln() - p.ln();
            assert!(ds >= -EPS, "M1={m1}: (s2-s1)/R = {ds} must be >= 0");
            if m1 > 1.0 {
                assert!(ds > 0.0, "M1={m1}: entropy must strictly rise");
            }
        }
    }

    #[test]
    fn domain_errors() {
        // Subsonic upstream Mach is rejected by every relation.
        assert!(normal_shock(0.8, 1.4).is_err());
        assert!(downstream_mach(0.99, 1.4).is_err());
        assert!(pressure_ratio(0.5, 1.4).is_err());
        assert!(stagnation_pressure_ratio(f64::NAN, 1.4).is_err());
        // gamma must be > 1.
        assert!(normal_shock(2.0, 1.0).is_err());
        assert!(pressure_ratio(2.0, 0.95).is_err());
    }
}
