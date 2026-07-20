//! Fin parameter and characteristic length.
//!
//! The 1-D constant-area fin equation is
//!
//! ```text
//! d^2 theta / dx^2 - m^2 * theta = 0,   theta(x) = T(x) - T_inf,
//! ```
//!
//! whose single grouping of the physical inputs is the **fin parameter**
//!
//! ```text
//! m = sqrt( h * P / (k * A) )      [units: 1 / length]
//! ```
//!
//! with convection coefficient `h` (W/m^2 K), wetted perimeter `P` (m),
//! thermal conductivity `k` (W/m K), and cross-sectional area `A` (m^2).
//! The dimensionless product `m * L` (the fin length scaled by `1/m`) is the
//! sole argument that governs efficiency and effectiveness.

use crate::error::{FinError, Result};

/// Fin parameter `m = sqrt(h P / (k A))`, in inverse metres.
///
/// All four inputs must be finite and strictly positive — a fin with zero
/// area, conductivity, perimeter, or convection coefficient is unphysical.
///
/// # Errors
///
/// Returns [`FinError::NonPositive`] (or [`FinError::NotFinite`]) if any of
/// `h`, `perimeter`, `k`, or `area` is non-positive or non-finite.
///
/// # Example
///
/// ```
/// use valenx_fin::fin_parameter;
/// // h=100, P=0.1, k=200, A=1e-4  ->  hP/kA = 10/0.02 = 500, m = sqrt(500).
/// let m = fin_parameter(100.0, 0.1, 200.0, 1.0e-4).unwrap();
/// assert!((m - 500.0_f64.sqrt()).abs() < 1e-9);
/// ```
pub fn fin_parameter(h: f64, perimeter: f64, k: f64, area: f64) -> Result<f64> {
    let h = FinError::positive("h", h)?;
    let perimeter = FinError::positive("perimeter", perimeter)?;
    let k = FinError::positive("k", k)?;
    let area = FinError::positive("area", area)?;
    Ok((h * perimeter / (k * area)).sqrt())
}

/// Dimensionless fin length `mL` from physical inputs.
///
/// Convenience wrapper that computes [`fin_parameter`] and multiplies by the
/// fin length `length` (m). This `mL` is exactly the argument consumed by
/// [`crate::efficiency()`] and friends.
///
/// # Errors
///
/// Propagates the validation errors of [`fin_parameter`], and additionally
/// rejects a non-positive or non-finite `length`.
///
/// # Example
///
/// ```
/// use valenx_fin::dimensionless_length;
/// let mlength = dimensionless_length(100.0, 0.1, 200.0, 1.0e-4, 0.05).unwrap();
/// // m = sqrt(500) ~= 22.3607, L = 0.05  ->  mL ~= 1.118034.
/// assert!((mlength - 0.05 * 500.0_f64.sqrt()).abs() < 1e-9);
/// ```
pub fn dimensionless_length(h: f64, perimeter: f64, k: f64, area: f64, length: f64) -> Result<f64> {
    let m = fin_parameter(h, perimeter, k, area)?;
    let length = FinError::positive("length", length)?;
    Ok(m * length)
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    #[test]
    fn fin_parameter_matches_hand_worked_value() {
        // hP/kA = (100*0.1)/(200*1e-4) = 10/0.02 = 500.
        let m = fin_parameter(100.0, 0.1, 200.0, 1.0e-4).unwrap();
        assert!((m - 500.0_f64.sqrt()).abs() < EPS);
    }

    #[test]
    fn fin_parameter_unit_inputs_give_one() {
        // h=P=k=A=1  ->  m = sqrt(1) = 1.
        let m = fin_parameter(1.0, 1.0, 1.0, 1.0).unwrap();
        assert!((m - 1.0).abs() < EPS);
    }

    #[test]
    fn fin_parameter_scales_as_sqrt_of_h() {
        // Quadrupling h must double m (sqrt dependence).
        let m1 = fin_parameter(50.0, 0.2, 150.0, 5.0e-4).unwrap();
        let m4 = fin_parameter(200.0, 0.2, 150.0, 5.0e-4).unwrap();
        assert!((m4 - 2.0 * m1).abs() < 1e-7);
    }

    #[test]
    fn fin_parameter_rejects_nonpositive_inputs() {
        assert_eq!(
            fin_parameter(0.0, 0.1, 200.0, 1e-4).unwrap_err().code(),
            "fin.non-positive"
        );
        assert!(fin_parameter(100.0, -0.1, 200.0, 1e-4).is_err());
        assert!(fin_parameter(100.0, 0.1, 0.0, 1e-4).is_err());
        assert!(fin_parameter(100.0, 0.1, 200.0, 0.0).is_err());
    }

    #[test]
    fn fin_parameter_rejects_non_finite() {
        assert_eq!(
            fin_parameter(f64::NAN, 0.1, 200.0, 1e-4)
                .unwrap_err()
                .code(),
            "fin.not-finite"
        );
    }

    #[test]
    fn dimensionless_length_matches_m_times_l() {
        let m = fin_parameter(100.0, 0.1, 200.0, 1.0e-4).unwrap();
        let mlength = dimensionless_length(100.0, 0.1, 200.0, 1.0e-4, 0.05).unwrap();
        assert!((mlength - m * 0.05).abs() < EPS);
    }

    #[test]
    fn dimensionless_length_rejects_bad_length() {
        assert!(dimensionless_length(100.0, 0.1, 200.0, 1e-4, 0.0).is_err());
        assert!(dimensionless_length(100.0, 0.1, 200.0, 1e-4, -0.05).is_err());
    }
}
