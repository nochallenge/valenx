//! Wheatstone-bridge balance condition.
//!
//! ## Model
//!
//! A Wheatstone bridge has four arms `R1, R2, R3, R4`. Label the two
//! series legs so that `R1` is above `R2` in the first leg and `R3`
//! is above `R4` in the second, both legs across the same excitation.
//! The bridge is *balanced* (zero galvanometer / detector voltage)
//! when the two voltage-divider ratios match:
//!
//! `R1 / R2 = R3 / R4`,
//!
//! equivalently the cross-products are equal, `R1 * R4 = R2 * R3`.
//! The cross-product form is the one used here because it has no
//! division and so stays exact for the comparison even when the
//! ratios themselves are awkward fractions.
//!
//! For an out-of-balance bridge the detector voltage (relative to
//! the excitation `Vex`) follows from the two dividers:
//!
//! `Vdet = Vex * (R2 / (R1 + R2) - R4 / (R3 + R4))`.
//!
//! At balance this difference is exactly zero, which the tests check
//! both directly and via the cross-product predicate.

use crate::error::{check_finite, check_positive, ResistorError};

/// Default relative tolerance for [`is_balanced`].
///
/// Two bridge ratios are treated as equal when the normalised
/// cross-product mismatch is within this fraction. It is exposed so
/// callers can compare against the same default they would get from
/// the no-tolerance helper.
pub const DEFAULT_BALANCE_TOL: f64 = 1e-9;

/// Detector (bridge-output) voltage of a Wheatstone bridge.
///
/// Returns `Vdet = Vex * (R2/(R1+R2) - R4/(R3+R4))`, the voltage
/// across the detector arm relative to the excitation `vex`. It is
/// exactly zero at balance.
///
/// All four resistances must be finite and strictly positive; `vex`
/// must be finite.
///
/// # Errors
///
/// - [`ResistorError::NonPositive`] / [`ResistorError::NonFinite`]
///   for an out-of-domain arm (named `"r1"`..`"r4"`).
/// - [`ResistorError::NonFinite`] if `vex` is `NaN`/infinite.
///
/// # Examples
///
/// ```
/// use valenx_resistor_network::bridge::detector_voltage;
/// // A balanced bridge has zero output.
/// let v = detector_voltage(5.0, 100.0, 100.0, 100.0, 100.0).unwrap();
/// assert!(v.abs() < 1e-12);
/// ```
pub fn detector_voltage(
    vex: f64,
    r1: f64,
    r2: f64,
    r3: f64,
    r4: f64,
) -> Result<f64, ResistorError> {
    let vex = check_finite("vex", vex)?;
    let r1 = check_positive("r1", r1)?;
    let r2 = check_positive("r2", r2)?;
    let r3 = check_positive("r3", r3)?;
    let r4 = check_positive("r4", r4)?;
    let left = r2 / (r1 + r2);
    let right = r4 / (r3 + r4);
    Ok(vex * (left - right))
}

/// Test whether a Wheatstone bridge is balanced within `tol`.
///
/// Balanced means `R1 / R2 = R3 / R4`. The comparison is done on the
/// equivalent cross-product form `R1 * R4 == R2 * R3`, normalised by
/// the product magnitude so `tol` is a relative tolerance.
///
/// All four resistances must be finite and strictly positive.
///
/// # Errors
///
/// [`ResistorError::NonPositive`] / [`ResistorError::NonFinite`] for
/// an out-of-domain arm, or [`ResistorError::NonFinite`] if `tol`
/// itself is `NaN`/infinite.
///
/// # Examples
///
/// ```
/// use valenx_resistor_network::bridge::{is_balanced, DEFAULT_BALANCE_TOL};
/// // 100/200 == 150/300, so this bridge is balanced.
/// assert!(is_balanced(100.0, 200.0, 150.0, 300.0, DEFAULT_BALANCE_TOL).unwrap());
/// // 100/200 != 100/300, so this one is not.
/// assert!(!is_balanced(100.0, 200.0, 100.0, 300.0, DEFAULT_BALANCE_TOL).unwrap());
/// ```
pub fn is_balanced(r1: f64, r2: f64, r3: f64, r4: f64, tol: f64) -> Result<bool, ResistorError> {
    let r1 = check_positive("r1", r1)?;
    let r2 = check_positive("r2", r2)?;
    let r3 = check_positive("r3", r3)?;
    let r4 = check_positive("r4", r4)?;
    let tol = check_finite("tol", tol)?;
    let cross_a = r1 * r4;
    let cross_b = r2 * r3;
    let scale = cross_a.abs().max(cross_b.abs()).max(1.0);
    Ok((cross_a - cross_b).abs() <= tol * scale)
}

/// Solve for the unknown arm `R4` that balances the bridge.
///
/// From `R1 / R2 = R3 / R4` the balancing value is
/// `R4 = R2 * R3 / R1`. Useful for the classic "measure an unknown
/// resistor by nulling the bridge" arrangement where `R4` is the
/// unknown and `R1, R2, R3` are known.
///
/// All three inputs must be finite and strictly positive.
///
/// # Errors
///
/// [`ResistorError::NonPositive`] / [`ResistorError::NonFinite`] for
/// an out-of-domain input.
///
/// # Examples
///
/// ```
/// use valenx_resistor_network::bridge::{balancing_r4, is_balanced, DEFAULT_BALANCE_TOL};
/// let r4 = balancing_r4(100.0, 200.0, 150.0).unwrap();
/// assert!((r4 - 300.0).abs() < 1e-9);
/// // The solved value does in fact balance the bridge.
/// assert!(is_balanced(100.0, 200.0, 150.0, r4, DEFAULT_BALANCE_TOL).unwrap());
/// ```
pub fn balancing_r4(r1: f64, r2: f64, r3: f64) -> Result<f64, ResistorError> {
    let r1 = check_positive("r1", r1)?;
    let r2 = check_positive("r2", r2)?;
    let r3 = check_positive("r3", r3)?;
    Ok(r2 * r3 / r1)
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    #[test]
    fn equal_arms_are_balanced() {
        assert!(is_balanced(100.0, 100.0, 100.0, 100.0, DEFAULT_BALANCE_TOL).expect("valid"));
    }

    #[test]
    fn matching_ratios_are_balanced() {
        // 100/200 == 150/300 == 1/2.
        assert!(is_balanced(100.0, 200.0, 150.0, 300.0, DEFAULT_BALANCE_TOL).expect("valid"));
    }

    #[test]
    fn mismatched_ratios_are_unbalanced() {
        // 100/200 = 1/2 but 100/300 = 1/3.
        assert!(!is_balanced(100.0, 200.0, 100.0, 300.0, DEFAULT_BALANCE_TOL).expect("valid"));
    }

    #[test]
    fn balanced_bridge_has_zero_detector_voltage() {
        let v = detector_voltage(5.0, 100.0, 200.0, 150.0, 300.0).expect("valid");
        assert!(
            v.abs() < EPS,
            "balanced bridge output should be ~0, got {v}"
        );
    }

    #[test]
    fn balanced_predicate_agrees_with_zero_output() {
        // For a range of balanced arms, the predicate is true exactly
        // when the detector voltage is (numerically) zero.
        let cases = [
            (10.0, 10.0, 10.0, 10.0),
            (100.0, 200.0, 150.0, 300.0),
            (47.0, 94.0, 33.0, 66.0),
            (1000.0, 1.0, 5000.0, 5.0),
        ];
        for (r1, r2, r3, r4) in cases {
            let balanced = is_balanced(r1, r2, r3, r4, DEFAULT_BALANCE_TOL).expect("valid");
            let v = detector_voltage(12.0, r1, r2, r3, r4).expect("valid");
            assert!(balanced, "expected balanced for {r1},{r2},{r3},{r4}");
            assert!(v.abs() < 1e-9, "output should be ~0, got {v}");
        }
    }

    #[test]
    fn unbalanced_bridge_has_nonzero_output_of_expected_sign() {
        // R2/(R1+R2) > R4/(R3+R4) -> positive output for positive Vex.
        // Left ratio 200/300 ~ 0.667, right ratio 100/400 = 0.25.
        let v = detector_voltage(10.0, 100.0, 200.0, 300.0, 100.0).expect("valid");
        assert!(v > 0.0, "expected positive output, got {v}");
        // Magnitude: 10 * (2/3 - 1/4) = 10 * 5/12 = 4.16666...
        assert!((v - 10.0 * (2.0 / 3.0 - 0.25)).abs() < EPS, "got {v}");
    }

    #[test]
    fn balancing_r4_known_value() {
        // R4 = R2 * R3 / R1 = 200 * 150 / 100 = 300.
        let r4 = balancing_r4(100.0, 200.0, 150.0).expect("valid");
        assert!((r4 - 300.0).abs() < EPS, "got {r4}");
    }

    #[test]
    fn balancing_r4_actually_balances() {
        let (r1, r2, r3) = (330.0, 470.0, 220.0);
        let r4 = balancing_r4(r1, r2, r3).expect("valid");
        assert!(is_balanced(r1, r2, r3, r4, DEFAULT_BALANCE_TOL).expect("valid"));
        let v = detector_voltage(5.0, r1, r2, r3, r4).expect("valid");
        assert!(v.abs() < 1e-9, "solved bridge output should be ~0, got {v}");
    }

    #[test]
    fn bridge_rejects_bad_arm() {
        assert_eq!(
            detector_voltage(5.0, 0.0, 100.0, 100.0, 100.0),
            Err(ResistorError::non_positive("r1", 0.0))
        );
        assert_eq!(
            is_balanced(100.0, 100.0, -1.0, 100.0, DEFAULT_BALANCE_TOL),
            Err(ResistorError::non_positive("r3", -1.0))
        );
        assert_eq!(
            balancing_r4(100.0, 100.0, f64::NEG_INFINITY),
            Err(ResistorError::non_finite("r3", f64::NEG_INFINITY))
        );
    }

    #[test]
    fn bridge_rejects_non_finite_excitation() {
        match detector_voltage(f64::NAN, 100.0, 100.0, 100.0, 100.0) {
            Err(ResistorError::NonFinite { name, .. }) => assert_eq!(name, "vex"),
            other => panic!("expected NonFinite, got {other:?}"),
        }
    }
}
