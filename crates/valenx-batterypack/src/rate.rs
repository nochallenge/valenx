//! C-rate ↔ current conversions.
//!
//! The **C-rate** expresses charge / discharge current as a multiple of
//! a capacity. For a capacity `Q` (in ampere-hours) the textbook
//! relation is:
//!
//! - Current from C-rate: `I = C · Q` (amperes). A `1C` rate on a
//!   3.0 Ah cell is 3.0 A and empties it in one hour; `2C` is 6.0 A in
//!   half an hour; `0.5C` is 1.5 A over two hours.
//! - C-rate from current: `C = I / Q`.
//! - Nominal run-time at a constant C-rate: `t = 1 / C` hours,
//!   independent of `Q` (the capacity cancels), so a `2C` draw lasts
//!   0.5 h and a `0.5C` draw lasts 2 h.
//!
//! These are *nominal* figures: real cells deliver less usable capacity
//! at higher C-rates (the Peukert effect) and sag under load. This
//! module is the idealised textbook relation only.

use crate::error::{require_non_negative, require_positive, BatteryPackError};

/// Constant-current draw for a given C-rate against a capacity.
///
/// Computes `I = c_rate · capacity_ah`, returning amperes. Works for a
/// single cell's capacity or a whole pack's capacity — the relation is
/// the same; pass whichever `Ah` figure the C-rate is referenced to.
///
/// # Errors
///
/// Returns [`BatteryPackError::BadParameter`] if `c_rate` is negative
/// or non-finite, or if `capacity_ah` is not strictly positive.
pub fn current_from_c_rate(c_rate: f64, capacity_ah: f64) -> Result<f64, BatteryPackError> {
    let c_rate = require_non_negative("c_rate", c_rate)?;
    let capacity_ah = require_positive("capacity_ah", capacity_ah)?;
    Ok(c_rate * capacity_ah)
}

/// C-rate that a given constant current represents against a capacity.
///
/// Computes `C = current_a / capacity_ah`, the inverse of
/// [`current_from_c_rate`].
///
/// # Errors
///
/// Returns [`BatteryPackError::BadParameter`] if `current_a` is negative
/// or non-finite, or if `capacity_ah` is not strictly positive.
pub fn c_rate_from_current(current_a: f64, capacity_ah: f64) -> Result<f64, BatteryPackError> {
    let current_a = require_non_negative("current_a", current_a)?;
    let capacity_ah = require_positive("capacity_ah", capacity_ah)?;
    Ok(current_a / capacity_ah)
}

/// Nominal constant-current run-time at a given C-rate, in hours.
///
/// Computes `t = 1 / c_rate`. The capacity cancels out of `Q / (C·Q)`,
/// so the run-time depends only on the C-rate: `1C` lasts 1 h, `2C`
/// lasts 0.5 h, `0.5C` lasts 2 h.
///
/// # Errors
///
/// Returns [`BatteryPackError::BadParameter`] if `c_rate` is not
/// strictly positive (a zero or negative rate has no finite run-time).
pub fn runtime_hours_at_c_rate(c_rate: f64) -> Result<f64, BatteryPackError> {
    let c_rate = require_positive("c_rate", c_rate)?;
    Ok(1.0 / c_rate)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn one_c_equals_capacity_in_amps() {
        // 1C on 3.0 Ah = 3.0 A.
        let i = current_from_c_rate(1.0, 3.0).unwrap();
        assert!((i - 3.0).abs() < 1e-12);
    }

    #[test]
    fn two_c_doubles_current() {
        // 2C on 3.0 Ah = 6.0 A.
        let i = current_from_c_rate(2.0, 3.0).unwrap();
        assert!((i - 6.0).abs() < 1e-12);
    }

    #[test]
    fn half_c_halves_current() {
        // 0.5C on 3.0 Ah = 1.5 A.
        let i = current_from_c_rate(0.5, 3.0).unwrap();
        assert!((i - 1.5).abs() < 1e-12);
    }

    #[test]
    fn zero_c_is_zero_current() {
        let i = current_from_c_rate(0.0, 3.0).unwrap();
        assert!(i.abs() < 1e-12);
    }

    #[test]
    fn c_rate_and_current_are_inverses() {
        // Round-trip: current -> C-rate -> current.
        let cap = 12.0;
        let c = c_rate_from_current(6.0, cap).unwrap();
        assert!((c - 0.5).abs() < 1e-12);
        let i = current_from_c_rate(c, cap).unwrap();
        assert!((i - 6.0).abs() < 1e-12);
    }

    #[test]
    fn runtime_is_reciprocal_of_c_rate() {
        assert!((runtime_hours_at_c_rate(1.0).unwrap() - 1.0).abs() < 1e-12);
        assert!((runtime_hours_at_c_rate(2.0).unwrap() - 0.5).abs() < 1e-12);
        assert!((runtime_hours_at_c_rate(0.5).unwrap() - 2.0).abs() < 1e-12);
    }

    #[test]
    fn rejects_bad_inputs() {
        assert!(current_from_c_rate(-1.0, 3.0).is_err());
        assert!(current_from_c_rate(1.0, 0.0).is_err());
        assert!(c_rate_from_current(-1.0, 3.0).is_err());
        assert!(c_rate_from_current(1.0, 0.0).is_err());
        assert!(runtime_hours_at_c_rate(0.0).is_err());
        assert!(current_from_c_rate(f64::NAN, 3.0).is_err());
    }
}
