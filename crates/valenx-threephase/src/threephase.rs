//! Balanced three-phase line/phase relations and real-power formulas.
//!
//! All quantities are RMS magnitudes of a *balanced* three-phase
//! system, where the three phases are equal in magnitude and spaced
//! 120 degrees apart. Under that assumption the line-to-phase relations
//! reduce to the standard `sqrt(3)` factors and the total real power is
//! `P = sqrt(3) * V_line * I_line * cos(phi)`.
//!
//! # Conventions
//!
//! `V_phase` / `I_phase` are the voltage across, and current through, a
//! single load element (one arm of the wye or one side of the delta).
//! `V_line` / `I_line` are the line-to-line voltage and the current in
//! a supply conductor. The angle `phi` is the phase angle between phase
//! voltage and phase current, so `cos(phi)` is the power factor.

use serde::{Deserialize, Serialize};

use crate::error::{require_positive, require_power_factor, ThreePhaseError};

/// `sqrt(3)`, the conversion factor between line and phase quantities in
/// a balanced three-phase system.
///
/// Exposed as a named constant so callers can reuse the exact value the
/// crate validates against rather than recomputing it.
pub const SQRT_3: f64 = 1.732_050_807_568_877_2;

/// How the three identical single-phase load elements are wired.
///
/// In a [`Connection::Wye`] (star, "Y") connection the three elements
/// share a common neutral point; in a [`Connection::Delta`] (mesh)
/// connection they are joined end-to-end in a closed loop. The wiring
/// determines which line/phase relation gets the `sqrt(3)` factor.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub enum Connection {
    /// Star / "Y" connection. Line voltage leads: `V_line = sqrt(3) *
    /// V_phase`, while `I_line = I_phase`.
    Wye,
    /// Mesh / triangle connection. Line current leads: `I_line =
    /// sqrt(3) * I_phase`, while `V_line = V_phase`.
    Delta,
}

impl Connection {
    /// Line-to-line voltage for this connection given the per-element
    /// phase voltage.
    ///
    /// Wye: `V_line = sqrt(3) * V_phase`. Delta: `V_line = V_phase`.
    #[must_use]
    pub fn line_voltage(self, v_phase: f64) -> f64 {
        match self {
            Connection::Wye => SQRT_3 * v_phase,
            Connection::Delta => v_phase,
        }
    }

    /// Per-element phase voltage for this connection given the
    /// line-to-line voltage. Inverse of [`Connection::line_voltage`].
    ///
    /// Wye: `V_phase = V_line / sqrt(3)`. Delta: `V_phase = V_line`.
    #[must_use]
    pub fn phase_voltage(self, v_line: f64) -> f64 {
        match self {
            Connection::Wye => v_line / SQRT_3,
            Connection::Delta => v_line,
        }
    }

    /// Line (conductor) current for this connection given the
    /// per-element phase current.
    ///
    /// Wye: `I_line = I_phase`. Delta: `I_line = sqrt(3) * I_phase`.
    #[must_use]
    pub fn line_current(self, i_phase: f64) -> f64 {
        match self {
            Connection::Wye => i_phase,
            Connection::Delta => SQRT_3 * i_phase,
        }
    }

    /// Per-element phase current for this connection given the line
    /// current. Inverse of [`Connection::line_current`].
    ///
    /// Wye: `I_phase = I_line`. Delta: `I_phase = I_line / sqrt(3)`.
    #[must_use]
    pub fn phase_current(self, i_line: f64) -> f64 {
        match self {
            Connection::Wye => i_line,
            Connection::Delta => i_line / SQRT_3,
        }
    }
}

/// A balanced three-phase load specified by its wiring and the RMS
/// magnitudes seen by a single load element.
///
/// Construct via [`BalancedLoad::new`], which validates that both
/// magnitudes are finite and strictly positive. The accessor methods
/// then derive line quantities, per-phase and total power.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct BalancedLoad {
    connection: Connection,
    v_phase: f64,
    i_phase: f64,
    power_factor: f64,
}

impl BalancedLoad {
    /// Build a balanced load from its connection, per-element phase
    /// voltage and current (RMS magnitudes), and power factor
    /// `cos(phi)`.
    ///
    /// # Errors
    ///
    /// Returns [`ThreePhaseError::NonPositive`] or
    /// [`ThreePhaseError::NotFinite`] if either magnitude is not a
    /// finite strictly-positive number, and
    /// [`ThreePhaseError::PowerFactorOutOfRange`] if `power_factor` is
    /// outside `[-1, 1]`.
    pub fn new(
        connection: Connection,
        v_phase: f64,
        i_phase: f64,
        power_factor: f64,
    ) -> Result<Self, ThreePhaseError> {
        let v_phase = require_positive("v_phase", v_phase)?;
        let i_phase = require_positive("i_phase", i_phase)?;
        let power_factor = require_power_factor(power_factor)?;
        Ok(Self {
            connection,
            v_phase,
            i_phase,
            power_factor,
        })
    }

    /// The wiring of this load.
    #[must_use]
    pub fn connection(&self) -> Connection {
        self.connection
    }

    /// Per-element phase voltage (RMS).
    #[must_use]
    pub fn phase_voltage(&self) -> f64 {
        self.v_phase
    }

    /// Per-element phase current (RMS).
    #[must_use]
    pub fn phase_current(&self) -> f64 {
        self.i_phase
    }

    /// Power factor `cos(phi)`.
    #[must_use]
    pub fn power_factor(&self) -> f64 {
        self.power_factor
    }

    /// Line-to-line voltage `V_line` for this load.
    #[must_use]
    pub fn line_voltage(&self) -> f64 {
        self.connection.line_voltage(self.v_phase)
    }

    /// Line (conductor) current `I_line` for this load.
    #[must_use]
    pub fn line_current(&self) -> f64 {
        self.connection.line_current(self.i_phase)
    }

    /// Real power dissipated in a single load element,
    /// `P_phase = V_phase * I_phase * cos(phi)`.
    #[must_use]
    pub fn power_per_phase(&self) -> f64 {
        self.v_phase * self.i_phase * self.power_factor
    }

    /// Total real power of the balanced three-phase load, computed from
    /// the line quantities as `P = sqrt(3) * V_line * I_line *
    /// cos(phi)`.
    ///
    /// For a balanced system this is identically three times
    /// [`BalancedLoad::power_per_phase`]; see the module tests for the
    /// algebraic cross-check.
    #[must_use]
    pub fn total_power(&self) -> f64 {
        SQRT_3 * self.line_voltage() * self.line_current() * self.power_factor
    }
}

/// Total real power of a balanced three-phase load directly from line
/// quantities: `P = sqrt(3) * V_line * I_line * cos(phi)`.
///
/// This is the form used when only the measured line-to-line voltage,
/// line current, and power factor are known.
///
/// # Errors
///
/// Returns [`ThreePhaseError::NonPositive`] / [`ThreePhaseError::NotFinite`]
/// if `v_line` or `i_line` is not finite-positive, and
/// [`ThreePhaseError::PowerFactorOutOfRange`] if `power_factor` is
/// outside `[-1, 1]`.
pub fn power_from_line(
    v_line: f64,
    i_line: f64,
    power_factor: f64,
) -> Result<f64, ThreePhaseError> {
    let v_line = require_positive("v_line", v_line)?;
    let i_line = require_positive("i_line", i_line)?;
    let power_factor = require_power_factor(power_factor)?;
    Ok(SQRT_3 * v_line * i_line * power_factor)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::error::ThreePhaseError;

    /// Absolute tolerance for floating-point ground-truth comparisons.
    const EPS: f64 = 1e-9;

    #[test]
    fn sqrt_3_constant_matches_library() {
        let diff = (SQRT_3 - 3.0_f64.sqrt()).abs();
        assert!(diff < EPS, "SQRT_3 deviates from 3f64.sqrt(): {diff}");
    }

    #[test]
    fn wye_line_voltage_is_sqrt3_times_phase() {
        // Canonical 230 V phase / 400 V line low-voltage example.
        let v_phase = 230.0;
        let v_line = Connection::Wye.line_voltage(v_phase);
        let expected = SQRT_3 * v_phase;
        assert!((v_line - expected).abs() < EPS, "got {v_line}");
        // Ground-truth numeric value: 230 * sqrt(3) ~= 398.371686...
        assert!((v_line - 398.371_685_7).abs() < 1e-6, "got {v_line}");
    }

    #[test]
    fn wye_line_current_equals_phase_current() {
        let i_phase = 12.5;
        let i_line = Connection::Wye.line_current(i_phase);
        assert!((i_line - i_phase).abs() < EPS, "got {i_line}");
    }

    #[test]
    fn delta_line_current_is_sqrt3_times_phase() {
        let i_phase = 10.0;
        let i_line = Connection::Delta.line_current(i_phase);
        let expected = SQRT_3 * i_phase;
        assert!((i_line - expected).abs() < EPS, "got {i_line}");
        // Ground-truth numeric value: 10 * sqrt(3) ~= 17.320508...
        assert!((i_line - 17.320_508_07).abs() < 1e-6, "got {i_line}");
    }

    #[test]
    fn delta_line_voltage_equals_phase_voltage() {
        let v_phase = 400.0;
        let v_line = Connection::Delta.line_voltage(v_phase);
        assert!((v_line - v_phase).abs() < EPS, "got {v_line}");
    }

    #[test]
    fn voltage_round_trips_through_phase_and_line() {
        for conn in [Connection::Wye, Connection::Delta] {
            let v_phase = 277.0;
            let back = conn.phase_voltage(conn.line_voltage(v_phase));
            assert!((back - v_phase).abs() < EPS, "{conn:?}: got {back}");
        }
    }

    #[test]
    fn current_round_trips_through_phase_and_line() {
        for conn in [Connection::Wye, Connection::Delta] {
            let i_phase = 8.4;
            let back = conn.phase_current(conn.line_current(i_phase));
            assert!((back - i_phase).abs() < EPS, "{conn:?}: got {back}");
        }
    }

    #[test]
    fn power_formula_matches_known_value() {
        // 400 V line, 10 A line, unity PF.
        // P = sqrt(3) * 400 * 10 * 1 = 6928.203... W.
        let p = power_from_line(400.0, 10.0, 1.0).unwrap();
        assert!((p - SQRT_3 * 400.0 * 10.0).abs() < EPS, "got {p}");
        assert!((p - 6_928.203_230_3).abs() < 1e-6, "got {p}");
    }

    #[test]
    fn power_formula_applies_power_factor() {
        // PF 0.8 lagging: P = sqrt(3) * 400 * 10 * 0.8 = 5542.562... W.
        let p = power_from_line(400.0, 10.0, 0.8).unwrap();
        let expected = SQRT_3 * 400.0 * 10.0 * 0.8;
        assert!((p - expected).abs() < EPS, "got {p}");
        assert!((p - 5_542.562_584_2).abs() < 1e-6, "got {p}");
    }

    #[test]
    fn total_power_is_three_times_per_phase_wye() {
        // Wye: per-phase P = V_phase * I_phase * PF.
        let load = BalancedLoad::new(Connection::Wye, 230.0, 12.0, 0.9).unwrap();
        let per_phase = load.power_per_phase();
        let total = load.total_power();
        assert!(
            (total - 3.0 * per_phase).abs() < 1e-6,
            "total {total} != 3 * {per_phase}"
        );
        // Explicit ground-truth: 3 * 230 * 12 * 0.9 = 7452 W.
        assert!((total - 7452.0).abs() < 1e-6, "got {total}");
    }

    #[test]
    fn total_power_is_three_times_per_phase_delta() {
        // Delta: V_line == V_phase, I_line == sqrt(3) * I_phase, so
        // sqrt(3) * V_line * I_line == 3 * V_phase * I_phase.
        let load = BalancedLoad::new(Connection::Delta, 400.0, 5.0, 0.85).unwrap();
        let per_phase = load.power_per_phase();
        let total = load.total_power();
        assert!(
            (total - 3.0 * per_phase).abs() < 1e-6,
            "total {total} != 3 * {per_phase}"
        );
        // Explicit ground-truth: 3 * 400 * 5 * 0.85 = 5100 W.
        assert!((total - 5100.0).abs() < 1e-6, "got {total}");
    }

    #[test]
    fn balanced_load_derives_line_quantities_wye() {
        let load = BalancedLoad::new(Connection::Wye, 230.0, 12.5, 1.0).unwrap();
        assert!((load.line_voltage() - SQRT_3 * 230.0).abs() < EPS);
        assert!((load.line_current() - 12.5).abs() < EPS);
    }

    #[test]
    fn balanced_load_derives_line_quantities_delta() {
        let load = BalancedLoad::new(Connection::Delta, 400.0, 10.0, 1.0).unwrap();
        assert!((load.line_voltage() - 400.0).abs() < EPS);
        assert!((load.line_current() - SQRT_3 * 10.0).abs() < EPS);
    }

    #[test]
    fn line_based_power_matches_balanced_load_total() {
        // The two independent power paths must agree for any wiring.
        for conn in [Connection::Wye, Connection::Delta] {
            let load = BalancedLoad::new(conn, 230.0, 9.0, 0.95).unwrap();
            let via_line = power_from_line(
                load.line_voltage(),
                load.line_current(),
                load.power_factor(),
            )
            .unwrap();
            assert!(
                (via_line - load.total_power()).abs() < 1e-6,
                "{conn:?}: {via_line} != {}",
                load.total_power()
            );
        }
    }

    #[test]
    fn rejects_non_positive_voltage() {
        let err = BalancedLoad::new(Connection::Wye, 0.0, 10.0, 1.0).unwrap_err();
        assert_eq!(
            err,
            ThreePhaseError::NonPositive {
                name: "v_phase",
                value: 0.0
            }
        );
        assert_eq!(err.code(), "threephase.non-positive");
    }

    #[test]
    fn rejects_non_positive_current() {
        let err = BalancedLoad::new(Connection::Delta, 400.0, -3.0, 1.0).unwrap_err();
        assert_eq!(
            err,
            ThreePhaseError::NonPositive {
                name: "i_phase",
                value: -3.0
            }
        );
    }

    #[test]
    fn rejects_non_finite_voltage() {
        let err = BalancedLoad::new(Connection::Wye, f64::NAN, 10.0, 1.0).unwrap_err();
        assert!(matches!(
            err,
            ThreePhaseError::NotFinite {
                name: "v_phase",
                ..
            }
        ));
        assert_eq!(err.code(), "threephase.not-finite");
    }

    #[test]
    fn rejects_power_factor_out_of_range() {
        let err = power_from_line(400.0, 10.0, 1.5).unwrap_err();
        assert_eq!(err, ThreePhaseError::PowerFactorOutOfRange { value: 1.5 });
        assert_eq!(err.code(), "threephase.power-factor-out-of-range");
    }

    #[test]
    fn power_factor_endpoints_are_accepted() {
        assert!(power_from_line(400.0, 10.0, 1.0).is_ok());
        assert!(power_from_line(400.0, 10.0, -1.0).is_ok());
        assert!(BalancedLoad::new(Connection::Wye, 230.0, 10.0, 0.0).is_ok());
    }
}
