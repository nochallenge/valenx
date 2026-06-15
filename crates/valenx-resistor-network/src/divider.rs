//! Voltage- and current-divider closed forms.
//!
//! ## Model
//!
//! A voltage divider is two resistors `R1` (top, source side) and
//! `R2` (bottom, output side) in series across a source `Vin`; the
//! voltage across `R2` is
//!
//! `Vout = Vin * R2 / (R1 + R2)`.
//!
//! A current divider splits a total current `I_in` entering two
//! parallel branches `R1` and `R2`; the current through `R1` is the
//! *opposite* resistor over the sum,
//!
//! `I1 = I_in * R2 / (R1 + R2)`,
//!
//! and symmetrically `I2 = I_in * R1 / (R1 + R2)`. The two branch
//! currents always sum back to `I_in` (Kirchhoff's current law).
//!
//! The resistances must be finite and strictly positive; the
//! source `Vin` / `I_in` may be any finite real (including zero or
//! negative, for a reversed source).

use crate::error::{check_finite, check_positive, ResistorError};

/// Output voltage of a two-resistor voltage divider.
///
/// Returns `Vout = Vin * R2 / (R1 + R2)`, the voltage measured
/// across `r2` (the bottom / output-side resistor).
///
/// `r1` and `r2` must be finite and strictly positive; `vin` must be
/// finite (it may be zero or negative).
///
/// # Errors
///
/// - [`ResistorError::NonPositive`] / [`ResistorError::NonFinite`]
///   if a resistance is out of domain.
/// - [`ResistorError::NonFinite`] if `vin` is `NaN`/infinite.
///
/// # Examples
///
/// ```
/// use valenx_resistor_network::divider::voltage_divider;
/// // Equal legs halve the source.
/// let v = voltage_divider(12.0, 1000.0, 1000.0).unwrap();
/// assert!((v - 6.0).abs() < 1e-9);
/// ```
pub fn voltage_divider(vin: f64, r1: f64, r2: f64) -> Result<f64, ResistorError> {
    let vin = check_finite("vin", vin)?;
    let r1 = check_positive("r1", r1)?;
    let r2 = check_positive("r2", r2)?;
    Ok(vin * r2 / (r1 + r2))
}

/// Current through branch `r1` of a two-branch current divider.
///
/// Returns `I1 = i_in * r2 / (r1 + r2)` — the *opposite* resistor
/// over the sum, since the larger branch resistance carries less of
/// the shared current.
///
/// `r1` and `r2` must be finite and strictly positive; `i_in` must
/// be finite.
///
/// # Errors
///
/// As [`voltage_divider`], with the resistance names `"r1"` / `"r2"`
/// and the source named `"i_in"`.
///
/// # Examples
///
/// ```
/// use valenx_resistor_network::divider::current_divider_i1;
/// // 1 A into 100 ohm || 300 ohm: the 100-ohm branch takes 0.75 A.
/// let i1 = current_divider_i1(1.0, 100.0, 300.0).unwrap();
/// assert!((i1 - 0.75).abs() < 1e-9);
/// ```
pub fn current_divider_i1(i_in: f64, r1: f64, r2: f64) -> Result<f64, ResistorError> {
    let i_in = check_finite("i_in", i_in)?;
    let r1 = check_positive("r1", r1)?;
    let r2 = check_positive("r2", r2)?;
    Ok(i_in * r2 / (r1 + r2))
}

/// Current through branch `r2` of a two-branch current divider.
///
/// Returns `I2 = i_in * r1 / (r1 + r2)`. Together with
/// [`current_divider_i1`] the two branch currents sum to `i_in`.
///
/// # Errors
///
/// As [`current_divider_i1`].
///
/// # Examples
///
/// ```
/// use valenx_resistor_network::divider::{current_divider_i1, current_divider_i2};
/// let i1 = current_divider_i1(1.0, 100.0, 300.0).unwrap();
/// let i2 = current_divider_i2(1.0, 100.0, 300.0).unwrap();
/// assert!((i1 + i2 - 1.0).abs() < 1e-9);
/// ```
pub fn current_divider_i2(i_in: f64, r1: f64, r2: f64) -> Result<f64, ResistorError> {
    let i_in = check_finite("i_in", i_in)?;
    let r1 = check_positive("r1", r1)?;
    let r2 = check_positive("r2", r2)?;
    Ok(i_in * r1 / (r1 + r2))
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    #[test]
    fn voltage_divider_equal_legs_halves() {
        let v = voltage_divider(12.0, 1000.0, 1000.0).expect("valid");
        assert!((v - 6.0).abs() < EPS, "got {v}");
    }

    #[test]
    fn voltage_divider_known_value() {
        // Ground truth: 9 V * 2k / (1k + 2k) = 6 V.
        let v = voltage_divider(9.0, 1000.0, 2000.0).expect("valid");
        assert!((v - 6.0).abs() < EPS, "got {v}");
    }

    #[test]
    fn voltage_divider_full_at_open_top() {
        // As R1 -> 0 the output approaches the full source. Use a
        // tiny top resistor as a sanity check on the limit.
        let v = voltage_divider(5.0, 1e-6, 1000.0).expect("valid");
        assert!((v - 5.0).abs() < 1e-3, "got {v}");
    }

    #[test]
    fn voltage_divider_zero_source_is_zero() {
        let v = voltage_divider(0.0, 100.0, 200.0).expect("valid");
        assert!(v.abs() < EPS, "got {v}");
    }

    #[test]
    fn voltage_divider_negative_source() {
        // A reversed source flips the output sign but keeps the ratio.
        let v = voltage_divider(-10.0, 1000.0, 1000.0).expect("valid");
        assert!((v - (-5.0)).abs() < EPS, "got {v}");
    }

    #[test]
    fn current_divider_known_split() {
        // Ground truth: 1 A into 100 || 300 -> 0.75 A and 0.25 A.
        let i1 = current_divider_i1(1.0, 100.0, 300.0).expect("valid");
        let i2 = current_divider_i2(1.0, 100.0, 300.0).expect("valid");
        assert!((i1 - 0.75).abs() < EPS, "i1 = {i1}");
        assert!((i2 - 0.25).abs() < EPS, "i2 = {i2}");
    }

    #[test]
    fn current_divider_branches_sum_to_input() {
        let i_in = 2.5;
        let i1 = current_divider_i1(i_in, 47.0, 220.0).expect("valid");
        let i2 = current_divider_i2(i_in, 47.0, 220.0).expect("valid");
        assert!((i1 + i2 - i_in).abs() < EPS, "sum = {}", i1 + i2);
    }

    #[test]
    fn current_divider_equal_branches_split_evenly() {
        let i1 = current_divider_i1(1.0, 500.0, 500.0).expect("valid");
        let i2 = current_divider_i2(1.0, 500.0, 500.0).expect("valid");
        assert!((i1 - 0.5).abs() < EPS, "i1 = {i1}");
        assert!((i2 - 0.5).abs() < EPS, "i2 = {i2}");
    }

    #[test]
    fn current_divider_consistent_with_ohms_law() {
        // The branch currents must equal V_parallel / R_branch, where
        // V_parallel = I_in * (R1 || R2).
        let i_in = 1.3;
        let (r1, r2) = (120.0, 270.0);
        let r_par = r1 * r2 / (r1 + r2);
        let v_par = i_in * r_par;
        let i1 = current_divider_i1(i_in, r1, r2).expect("valid");
        let i2 = current_divider_i2(i_in, r1, r2).expect("valid");
        assert!((i1 - v_par / r1).abs() < EPS, "i1 = {i1}");
        assert!((i2 - v_par / r2).abs() < EPS, "i2 = {i2}");
    }

    #[test]
    fn divider_rejects_bad_resistance() {
        assert_eq!(
            voltage_divider(5.0, 0.0, 100.0),
            Err(ResistorError::non_positive("r1", 0.0))
        );
        assert_eq!(
            current_divider_i1(1.0, 100.0, -5.0),
            Err(ResistorError::non_positive("r2", -5.0))
        );
    }

    #[test]
    fn divider_rejects_non_finite_source() {
        match voltage_divider(f64::NAN, 100.0, 100.0) {
            Err(ResistorError::NonFinite { name, .. }) => assert_eq!(name, "vin"),
            other => panic!("expected NonFinite, got {other:?}"),
        }
        match current_divider_i1(f64::INFINITY, 100.0, 100.0) {
            Err(ResistorError::NonFinite { name, .. }) => assert_eq!(name, "i_in"),
            other => panic!("expected NonFinite, got {other:?}"),
        }
    }
}
