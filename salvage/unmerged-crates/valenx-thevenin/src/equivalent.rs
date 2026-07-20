//! Thevenin and Norton equivalent two-terminal sources and the source
//! transformation that converts between them.
//!
//! A linear two-terminal DC network is fully described, as seen from its
//! port, by either:
//!
//! Thevenin form: an ideal voltage source `v_th` in series with a
//! resistance `r_th`. The open-circuit terminal voltage is `v_th`.
//!
//! Norton form: an ideal current source `i_n` in parallel with the same
//! resistance `r_n = r_th`. The short-circuit terminal current is `i_n`.
//!
//! ## Source transformation
//!
//! The two forms describe the *same* external behaviour and are related by
//! Ohm's law applied to the internal resistance:
//!
//! ```text
//! i_n  = v_th / r_th        r_n  = r_th
//! v_th = i_n  * r_n         r_th = r_n
//! ```
//!
//! Converting Thevenin to Norton and back is the identity (a round-trip),
//! which is one of the analytic ground truths exercised by the tests.

use crate::error::{finite, positive_resistance, TheveninError};

/// A Thevenin equivalent: ideal voltage source `v_th` in series with
/// internal resistance `r_th` (> 0).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Thevenin {
    /// Open-circuit voltage `V_th` (volts). May be any finite value
    /// (sign denotes polarity).
    pub v_th: f64,
    /// Internal / equivalent resistance `R_th` (ohms), strictly positive.
    pub r_th: f64,
}

/// A Norton equivalent: ideal current source `i_n` in parallel with
/// internal resistance `r_n` (> 0).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Norton {
    /// Short-circuit current `I_n` (amperes). May be any finite value
    /// (sign denotes direction).
    pub i_n: f64,
    /// Internal / equivalent resistance `R_n` (ohms), strictly positive.
    pub r_n: f64,
}

impl Thevenin {
    /// Construct a validated Thevenin source.
    ///
    /// `v_th` must be finite; `r_th` must be finite and strictly positive.
    ///
    /// # Errors
    ///
    /// Returns [`TheveninError::NonFinite`] for a non-finite `v_th`, and
    /// [`TheveninError::NonFinite`] or
    /// [`TheveninError::NonPositiveResistance`] for a bad `r_th`.
    pub fn new(v_th: f64, r_th: f64) -> Result<Self, TheveninError> {
        let v_th = finite("v_th", v_th)?;
        let r_th = positive_resistance("r_th", r_th)?;
        Ok(Self { v_th, r_th })
    }

    /// Convert to the equivalent Norton source via the source
    /// transformation `i_n = v_th / r_th`, `r_n = r_th`.
    pub fn to_norton(self) -> Norton {
        Norton {
            i_n: self.v_th / self.r_th,
            r_n: self.r_th,
        }
    }

    /// Terminal (load) current delivered into a resistive load `r_load`
    /// (>= 0 ohms), from the voltage divider `i = v_th / (r_th + r_load)`.
    ///
    /// # Errors
    ///
    /// Returns an error if `r_load` is non-finite or negative.
    pub fn load_current(self, r_load: f64) -> Result<f64, TheveninError> {
        let r_load = crate::error::non_negative("r_load", r_load)?;
        Ok(self.v_th / (self.r_th + r_load))
    }

    /// Terminal voltage across a resistive load `r_load` (>= 0 ohms),
    /// from the voltage divider `v = v_th * r_load / (r_th + r_load)`.
    ///
    /// # Errors
    ///
    /// Returns an error if `r_load` is non-finite or negative.
    pub fn load_voltage(self, r_load: f64) -> Result<f64, TheveninError> {
        let r_load = crate::error::non_negative("r_load", r_load)?;
        Ok(self.v_th * r_load / (self.r_th + r_load))
    }
}

impl Norton {
    /// Construct a validated Norton source.
    ///
    /// `i_n` must be finite; `r_n` must be finite and strictly positive.
    ///
    /// # Errors
    ///
    /// Returns [`TheveninError::NonFinite`] for a non-finite `i_n`, and
    /// [`TheveninError::NonFinite`] or
    /// [`TheveninError::NonPositiveResistance`] for a bad `r_n`.
    pub fn new(i_n: f64, r_n: f64) -> Result<Self, TheveninError> {
        let i_n = finite("i_n", i_n)?;
        let r_n = positive_resistance("r_n", r_n)?;
        Ok(Self { i_n, r_n })
    }

    /// Convert to the equivalent Thevenin source via the source
    /// transformation `v_th = i_n * r_n`, `r_th = r_n`.
    pub fn to_thevenin(self) -> Thevenin {
        Thevenin {
            v_th: self.i_n * self.r_n,
            r_th: self.r_n,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() < EPS
    }

    #[test]
    fn thevenin_constructor_validates() {
        let t = Thevenin::new(12.0, 4.0).unwrap();
        assert!(close(t.v_th, 12.0));
        assert!(close(t.r_th, 4.0));

        // Non-finite voltage rejected.
        assert!(matches!(
            Thevenin::new(f64::NAN, 4.0),
            Err(TheveninError::NonFinite { name: "v_th", .. })
        ));
        // Zero resistance rejected.
        assert!(matches!(
            Thevenin::new(12.0, 0.0),
            Err(TheveninError::NonPositiveResistance { .. })
        ));
        // Negative resistance rejected.
        assert!(matches!(
            Thevenin::new(12.0, -3.0),
            Err(TheveninError::NonPositiveResistance { .. })
        ));
    }

    #[test]
    fn norton_constructor_validates() {
        let n = Norton::new(3.0, 4.0).unwrap();
        assert!(close(n.i_n, 3.0));
        assert!(close(n.r_n, 4.0));

        assert!(matches!(
            Norton::new(f64::INFINITY, 4.0),
            Err(TheveninError::NonFinite { name: "i_n", .. })
        ));
        assert!(matches!(
            Norton::new(3.0, 0.0),
            Err(TheveninError::NonPositiveResistance { .. })
        ));
    }

    #[test]
    fn thevenin_to_norton_textbook_value() {
        // Textbook: V_th = 12 V, R_th = 4 ohm  =>  I_n = 12/4 = 3 A.
        let t = Thevenin::new(12.0, 4.0).unwrap();
        let n = t.to_norton();
        assert!(close(n.i_n, 3.0));
        assert!(close(n.r_n, 4.0));
    }

    #[test]
    fn norton_to_thevenin_textbook_value() {
        // I_n = 2 A, R_n = 5 ohm  =>  V_th = 2 * 5 = 10 V.
        let n = Norton::new(2.0, 5.0).unwrap();
        let t = n.to_thevenin();
        assert!(close(t.v_th, 10.0));
        assert!(close(t.r_th, 5.0));
    }

    #[test]
    fn source_transform_round_trip_is_identity() {
        // Ground truth: Thevenin -> Norton -> Thevenin recovers the original.
        let t0 = Thevenin::new(9.0, 3.0).unwrap();
        let t1 = t0.to_norton().to_thevenin();
        assert!(close(t1.v_th, t0.v_th));
        assert!(close(t1.r_th, t0.r_th));

        // And the other direction: Norton -> Thevenin -> Norton.
        let n0 = Norton::new(1.5, 8.0).unwrap();
        let n1 = n0.to_thevenin().to_norton();
        assert!(close(n1.i_n, n0.i_n));
        assert!(close(n1.r_n, n0.r_n));
    }

    #[test]
    fn round_trip_with_negative_polarity() {
        // Negative source values are legal (polarity / direction).
        let t0 = Thevenin::new(-6.0, 2.0).unwrap();
        let t1 = t0.to_norton().to_thevenin();
        assert!(close(t1.v_th, -6.0));
        assert!(close(t1.r_th, 2.0));
        assert!(close(t0.to_norton().i_n, -3.0));
    }

    #[test]
    fn load_current_voltage_divider() {
        // V_th = 10, R_th = 5. Into R_load = 5: i = 10/(5+5) = 1 A.
        let t = Thevenin::new(10.0, 5.0).unwrap();
        assert!(close(t.load_current(5.0).unwrap(), 1.0));
        // Open circuit (very large load) -> current -> 0; here exact limit
        // at infinite load is not allowed, so check a large finite load.
        assert!(t.load_current(1.0e12).unwrap() < 1.0e-9);
        // Short circuit (R_load = 0) -> i = V_th / R_th = 2 A.
        assert!(close(t.load_current(0.0).unwrap(), 2.0));
        // Negative load rejected.
        assert!(matches!(
            t.load_current(-1.0),
            Err(TheveninError::Negative { .. })
        ));
    }

    #[test]
    fn load_voltage_divider_and_open_circuit_limit() {
        // V_th = 10, R_th = 5.
        let t = Thevenin::new(10.0, 5.0).unwrap();
        // R_load = R_th = 5: v = 10 * 5/(5+5) = 5 V (half of V_th).
        assert!(close(t.load_voltage(5.0).unwrap(), 5.0));
        // R_load = 0 (short): terminal voltage collapses to 0.
        assert!(close(t.load_voltage(0.0).unwrap(), 0.0));
        // Very large load: terminal voltage approaches V_th (open circuit).
        assert!((t.load_voltage(1.0e12).unwrap() - 10.0).abs() < 1.0e-6);
        assert!(matches!(
            t.load_voltage(-2.0),
            Err(TheveninError::Negative { .. })
        ));
    }
}
