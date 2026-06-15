//! Transmitted hydraulic power.
//!
//! The power carried by a hydraulic stream is the product of pressure
//! and volumetric flow:
//!
//! ```text
//! Power = p * Q
//! ```
//!
//! In SI, pascals times cubic-metres-per-second gives watts. This is
//! the ideal fluid power crossing a section; it does not account for
//! pump or motor efficiency, so the mechanical shaft power needed to
//! supply it is larger in practice.

use crate::error::HydraulicsError;

/// Hydraulic power transmitted at pressure `pressure_pa` carrying
/// volumetric flow `flow_m3_s`, watts, via `Power = p Q`.
///
/// With pascals and cubic metres per second the result is watts. A
/// zero flow yields zero power.
///
/// # Errors
///
/// - [`HydraulicsError::NonPositive`] if `pressure_pa` is not strictly
///   positive (or not finite).
/// - [`HydraulicsError::Negative`] if `flow_m3_s` is negative (or not
///   finite). Zero is allowed and gives zero power.
pub fn hydraulic_power(pressure_pa: f64, flow_m3_s: f64) -> Result<f64, HydraulicsError> {
    let p = HydraulicsError::require_positive("pressure_pa", pressure_pa)?;
    let q = HydraulicsError::require_non_negative("flow_m3_s", flow_m3_s)?;
    Ok(p * q)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Tolerance for floating-point comparisons.
    const EPS: f64 = 1e-6;

    #[test]
    fn power_equals_pressure_times_flow() {
        // 200 bar = 20e6 Pa, 60 L/min = 1e-3 m^3/s  ->  20000 W.
        let p = 20.0e6;
        let q = 1.0e-3;
        let power = hydraulic_power(p, q).unwrap();
        assert!((power - p * q).abs() < EPS);
        assert!((power - 20_000.0).abs() < 1e-3);
    }

    #[test]
    fn power_is_linear_in_each_argument() {
        let base = hydraulic_power(10.0e6, 2.0e-3).unwrap();
        let double_p = hydraulic_power(20.0e6, 2.0e-3).unwrap();
        let double_q = hydraulic_power(10.0e6, 4.0e-3).unwrap();
        assert!((double_p - 2.0 * base).abs() < EPS);
        assert!((double_q - 2.0 * base).abs() < EPS);
    }

    #[test]
    fn zero_flow_gives_zero_power() {
        let power = hydraulic_power(35.0e6, 0.0).unwrap();
        assert!(power.abs() < EPS);
    }

    #[test]
    fn rejects_bad_inputs() {
        assert_eq!(
            hydraulic_power(0.0, 1.0e-3).unwrap_err().code(),
            "hydraulics.non_positive"
        );
        assert_eq!(
            hydraulic_power(1.0e6, -1.0e-3).unwrap_err().code(),
            "hydraulics.negative"
        );
    }
}
