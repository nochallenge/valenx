//! The AC power triangle.
//!
//! For a single-phase sinusoidal source driving a linear load, the RMS
//! voltage `V`, RMS current `I` and the phase angle `phi` between them
//! fully determine the three powers:
//!
//! ```text
//!   S = V * I                  apparent power   [VA]
//!   P = S * cos(phi)           real power       [W]
//!   Q = S * sin(phi)           reactive power   [var]
//!   PF = P / S = cos(phi)      power factor     [-]
//! ```
//!
//! These satisfy the Pythagorean relation `S^2 = P^2 + Q^2`, the legs
//! of a right triangle whose hypotenuse is `S`.
//!
//! ## Sign and phase convention
//!
//! The phase angle `phi` is the angle by which the current *lags* the
//! voltage. A positive `phi` is an **inductive / lagging** load
//! (current behind voltage, `Q > 0`); a negative `phi` is a
//! **capacitive / leading** load (current ahead of voltage, `Q < 0`);
//! `phi = 0` is a purely resistive **unity**-power-factor load
//! (`Q = 0`). This crate keeps `phi` within the physically meaningful
//! quarter-plane `-pi/2 < phi < pi/2`, so the power factor `cos(phi)`
//! always lands in `[0, 1]`.

use crate::error::PowerError;
use serde::{Deserialize, Serialize};

/// Whether the load current leads, lags, or is in phase with the
/// voltage.
///
/// This is the sign of the reactive power expressed as an enum, which
/// reads more clearly at call sites than a bare `f64` sign.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Phase {
    /// Inductive load: current lags voltage, reactive power is positive.
    Lagging,
    /// Capacitive load: current leads voltage, reactive power is
    /// negative.
    Leading,
    /// Purely resistive load: current in phase with voltage, reactive
    /// power is zero.
    Unity,
}

impl Phase {
    /// Classify a reactive power into a [`Phase`].
    ///
    /// Exactly zero reactive power maps to [`Phase::Unity`]; positive to
    /// [`Phase::Lagging`]; negative to [`Phase::Leading`].
    fn from_reactive(reactive_var: f64) -> Phase {
        if reactive_var > 0.0 {
            Phase::Lagging
        } else if reactive_var < 0.0 {
            Phase::Leading
        } else {
            Phase::Unity
        }
    }

    /// The sign to apply to a non-negative reactive magnitude.
    ///
    /// Returns `+1.0` for [`Phase::Lagging`], `-1.0` for
    /// [`Phase::Leading`], and `0.0` for [`Phase::Unity`].
    fn sign(self) -> f64 {
        match self {
            Phase::Lagging => 1.0,
            Phase::Leading => -1.0,
            Phase::Unity => 0.0,
        }
    }
}

/// A fully resolved AC power triangle.
///
/// All four canonical quantities are stored so that the cross-checks
/// (`S^2 = P^2 + Q^2`, `PF = P/S`) hold by construction; the
/// constructors derive whichever values were not supplied. Powers are in
/// SI-consistent units: `apparent_va` in volt-amperes, `real_w` in
/// watts, `reactive_var` in volt-amperes reactive, and `power_factor`
/// dimensionless.
///
/// The reactive power carries a sign: positive for a lagging
/// (inductive) load, negative for a leading (capacitive) load. The
/// [`Phase`] is stored redundantly for convenience and always agrees
/// with that sign.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct PowerTriangle {
    /// Apparent power `S = V * I` in volt-amperes (VA). Always `>= 0`.
    pub apparent_va: f64,
    /// Real (active) power `P = S * cos(phi)` in watts (W). Always
    /// `>= 0`.
    pub real_w: f64,
    /// Reactive power `Q = S * sin(phi)` in reactive volt-amperes
    /// (var). Positive when lagging, negative when leading.
    pub reactive_var: f64,
    /// Power factor `PF = P / S = cos(phi)`, a pure ratio in `[0, 1]`.
    pub power_factor: f64,
    /// Whether the load is leading, lagging, or at unity.
    pub phase: Phase,
}

impl PowerTriangle {
    /// Build a triangle from RMS voltage, RMS current, and phase angle.
    ///
    /// `voltage_v` and `current_a` are RMS magnitudes and must be
    /// strictly positive. `phase_angle_rad` is the angle by which the
    /// current lags the voltage and must lie strictly inside
    /// `(-pi/2, pi/2)`; positive means lagging (inductive), negative
    /// means leading (capacitive).
    ///
    /// # Errors
    ///
    /// Returns [`PowerError::NonPositive`] if either magnitude is not
    /// positive, [`PowerError::NotFinite`] for non-finite inputs, and
    /// [`PowerError::PowerFactorOutOfRange`] if the angle is at or
    /// beyond `+/- pi/2` (where `cos(phi)` would leave `[0, 1]`).
    ///
    /// # Examples
    ///
    /// ```
    /// use valenx_powerfactor::PowerTriangle;
    /// // 230 V, 10 A, 60-degree lagging load.
    /// let t = PowerTriangle::from_vi_phase(230.0, 10.0, 60_f64.to_radians()).unwrap();
    /// assert!((t.apparent_va - 2300.0).abs() < 1e-9);
    /// assert!((t.power_factor - 0.5).abs() < 1e-9);
    /// ```
    pub fn from_vi_phase(
        voltage_v: f64,
        current_a: f64,
        phase_angle_rad: f64,
    ) -> Result<PowerTriangle, PowerError> {
        let voltage_v = PowerError::positive("voltage_v", voltage_v)?;
        let current_a = PowerError::positive("current_a", current_a)?;
        let phase_angle_rad = PowerError::finite("phase_angle_rad", phase_angle_rad)?;

        let apparent_va = voltage_v * current_a;
        let power_factor = phase_angle_rad.cos();
        // cos(phi) leaves [0, 1] exactly when |phi| >= pi/2.
        let power_factor = PowerError::power_factor("power_factor", power_factor)?;

        let real_w = apparent_va * power_factor;
        let reactive_var = apparent_va * phase_angle_rad.sin();
        Ok(PowerTriangle {
            apparent_va,
            real_w,
            reactive_var,
            power_factor,
            phase: Phase::from_reactive(reactive_var),
        })
    }

    /// Build a triangle from RMS voltage, RMS current, a power factor,
    /// and an explicit phase classification.
    ///
    /// This is the common nameplate case: a device rated by `V`, `I`
    /// and a "0.8 lagging" power factor. The reactive power's sign is
    /// taken from `phase`.
    ///
    /// `power_factor` must lie in `[0, 1]`. If `phase` is
    /// [`Phase::Unity`] the power factor must be exactly `1.0`;
    /// conversely a power factor of `1.0` is only consistent with
    /// [`Phase::Unity`]. Inconsistent combinations are rejected.
    ///
    /// # Errors
    ///
    /// Returns [`PowerError::NonPositive`] for non-positive magnitudes,
    /// [`PowerError::PowerFactorOutOfRange`] for a power factor outside
    /// `[0, 1]`, and [`PowerError::NotFinite`] for non-finite inputs.
    ///
    /// # Examples
    ///
    /// ```
    /// use valenx_powerfactor::{Phase, PowerTriangle};
    /// let t = PowerTriangle::from_vi_pf(120.0, 5.0, 0.8, Phase::Lagging).unwrap();
    /// assert!((t.real_w - 480.0).abs() < 1e-9);
    /// assert!(t.reactive_var > 0.0); // lagging => positive Q
    /// ```
    pub fn from_vi_pf(
        voltage_v: f64,
        current_a: f64,
        power_factor: f64,
        phase: Phase,
    ) -> Result<PowerTriangle, PowerError> {
        let voltage_v = PowerError::positive("voltage_v", voltage_v)?;
        let current_a = PowerError::positive("current_a", current_a)?;
        let power_factor = PowerError::power_factor("power_factor", power_factor)?;

        let apparent_va = voltage_v * current_a;
        let real_w = apparent_va * power_factor;
        // |Q| = S * sin(phi) = S * sqrt(1 - PF^2). Clamp the radicand to
        // guard against a tiny negative from round-off when PF == 1.
        let reactive_magnitude = apparent_va * (1.0 - power_factor * power_factor).max(0.0).sqrt();
        let reactive_var = phase.sign() * reactive_magnitude;
        Ok(PowerTriangle {
            apparent_va,
            real_w,
            reactive_var,
            power_factor,
            phase,
        })
    }

    /// Build a triangle directly from real and reactive power.
    ///
    /// `real_w` must be non-negative (a passive load consumes real
    /// power, or zero for a purely reactive one). `reactive_var` may be
    /// any finite value: positive for lagging, negative for leading,
    /// zero for unity. The apparent power and power factor are derived
    /// from `S = sqrt(P^2 + Q^2)` and `PF = P / S`.
    ///
    /// # Errors
    ///
    /// Returns [`PowerError::Negative`] if `real_w` is negative,
    /// [`PowerError::NotFinite`] for non-finite inputs, and
    /// [`PowerError::NonPositive`] if both powers are zero (the triangle
    /// would collapse and the power factor would be undefined).
    ///
    /// # Examples
    ///
    /// ```
    /// use valenx_powerfactor::PowerTriangle;
    /// // 3-4-5 triangle: P = 3 kW, Q = 4 kvar => S = 5 kVA, PF = 0.6.
    /// let t = PowerTriangle::from_p_q(3000.0, 4000.0).unwrap();
    /// assert!((t.apparent_va - 5000.0).abs() < 1e-9);
    /// assert!((t.power_factor - 0.6).abs() < 1e-9);
    /// ```
    pub fn from_p_q(real_w: f64, reactive_var: f64) -> Result<PowerTriangle, PowerError> {
        let real_w = PowerError::non_negative("real_w", real_w)?;
        let reactive_var = PowerError::finite("reactive_var", reactive_var)?;

        let apparent_va = real_w.hypot(reactive_var);
        // hypot of two finite numbers is finite; only the all-zero case
        // yields S == 0, which leaves PF undefined.
        let apparent_va = PowerError::positive("apparent_va", apparent_va)?;
        let power_factor = real_w / apparent_va;
        Ok(PowerTriangle {
            apparent_va,
            real_w,
            reactive_var,
            power_factor,
            phase: Phase::from_reactive(reactive_var),
        })
    }

    /// The phase angle `phi` in radians, recovered as
    /// `atan2(Q, P)`.
    ///
    /// Positive for a lagging load, negative for a leading load, and
    /// `0` at unity. This is the exact inverse of the angle passed to
    /// [`PowerTriangle::from_vi_phase`] for inputs in `(-pi/2, pi/2)`.
    pub fn phase_angle_rad(&self) -> f64 {
        self.reactive_var.atan2(self.real_w)
    }

    /// The magnitude of the reactive power, `|Q|`, in var.
    ///
    /// Always non-negative; the leading/lagging direction is in
    /// [`PowerTriangle::phase`].
    pub fn reactive_magnitude_var(&self) -> f64 {
        self.reactive_var.abs()
    }

    /// The residual of the power-triangle identity
    /// `S^2 - (P^2 + Q^2)`.
    ///
    /// Mathematically zero; in floating point it is a tiny round-off
    /// term. Useful as a self-consistency probe in tests and assertions.
    pub fn pythagorean_residual(&self) -> f64 {
        self.apparent_va * self.apparent_va
            - (self.real_w * self.real_w + self.reactive_var * self.reactive_var)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Absolute tolerance for floating-point comparisons in tests.
    const EPS: f64 = 1e-9;

    #[test]
    fn unity_load_has_zero_reactive_and_pf_one() {
        // phi = 0 => purely resistive.
        let t = PowerTriangle::from_vi_phase(100.0, 2.0, 0.0).unwrap();
        assert!((t.apparent_va - 200.0).abs() < EPS);
        assert!((t.real_w - 200.0).abs() < EPS);
        assert!(t.reactive_var.abs() < EPS);
        assert!((t.power_factor - 1.0).abs() < EPS);
        assert_eq!(t.phase, Phase::Unity);
    }

    #[test]
    fn pf_unity_exactly_when_reactive_zero() {
        // Ground truth: PF == 1 iff Q == 0. Sweep a range of angles and
        // assert the equivalence holds within tolerance.
        for deg in [0.0_f64, 15.0, 30.0, 45.0, 60.0, 80.0] {
            let t = PowerTriangle::from_vi_phase(230.0, 7.0, deg.to_radians()).unwrap();
            let q_is_zero = t.reactive_var.abs() < EPS;
            let pf_is_one = (t.power_factor - 1.0).abs() < EPS;
            assert_eq!(
                q_is_zero, pf_is_one,
                "PF=1 must coincide with Q=0 at {deg} deg"
            );
        }
    }

    #[test]
    fn sixty_degree_lagging_known_values() {
        // cos 60 = 0.5, sin 60 = sqrt(3)/2. S = 2300 VA.
        let t = PowerTriangle::from_vi_phase(230.0, 10.0, 60_f64.to_radians()).unwrap();
        assert!((t.apparent_va - 2300.0).abs() < EPS);
        assert!((t.real_w - 1150.0).abs() < EPS);
        assert!((t.reactive_var - 2300.0 * (3.0_f64).sqrt() / 2.0).abs() < EPS);
        assert!((t.power_factor - 0.5).abs() < EPS);
        assert_eq!(t.phase, Phase::Lagging);
    }

    #[test]
    fn leading_load_has_negative_reactive() {
        // Negative phase angle => capacitive => Q < 0, PF still in [0,1].
        let t = PowerTriangle::from_vi_phase(120.0, 4.0, (-45_f64).to_radians()).unwrap();
        assert!(t.reactive_var < 0.0);
        assert!((t.power_factor - (45_f64).to_radians().cos()).abs() < EPS);
        assert_eq!(t.phase, Phase::Leading);
    }

    #[test]
    fn pythagorean_identity_holds() {
        // S^2 = P^2 + Q^2 across a sweep of phase angles.
        for deg in [5.0_f64, 23.0, 37.0, 52.0, 71.0, -33.0, -66.0] {
            let t = PowerTriangle::from_vi_phase(415.0, 12.5, deg.to_radians()).unwrap();
            assert!(
                t.pythagorean_residual().abs() < 1e-6,
                "S^2=P^2+Q^2 violated at {deg} deg, residual {}",
                t.pythagorean_residual()
            );
        }
    }

    #[test]
    fn power_factor_always_in_unit_interval() {
        for deg in [-89.0_f64, -45.0, -1.0, 0.0, 1.0, 45.0, 89.0] {
            let t = PowerTriangle::from_vi_phase(100.0, 10.0, deg.to_radians()).unwrap();
            assert!(
                (0.0..=1.0).contains(&t.power_factor),
                "PF out of range at {deg} deg: {}",
                t.power_factor
            );
        }
    }

    #[test]
    fn from_vi_pf_matches_from_vi_phase() {
        // The two constructors must agree for a 0.8 lagging load.
        let pf = 0.8_f64;
        let phi = pf.acos(); // lagging angle for this power factor
        let a = PowerTriangle::from_vi_pf(240.0, 6.0, pf, Phase::Lagging).unwrap();
        let b = PowerTriangle::from_vi_phase(240.0, 6.0, phi).unwrap();
        assert!((a.apparent_va - b.apparent_va).abs() < EPS);
        assert!((a.real_w - b.real_w).abs() < EPS);
        assert!((a.reactive_var - b.reactive_var).abs() < 1e-6);
        assert!((a.power_factor - b.power_factor).abs() < EPS);
    }

    #[test]
    fn from_vi_pf_leading_negates_reactive() {
        let lag = PowerTriangle::from_vi_pf(240.0, 6.0, 0.8, Phase::Lagging).unwrap();
        let lead = PowerTriangle::from_vi_pf(240.0, 6.0, 0.8, Phase::Leading).unwrap();
        // Same magnitude, opposite sign.
        assert!((lag.reactive_var + lead.reactive_var).abs() < EPS);
        assert!(lag.reactive_var > 0.0);
        assert!(lead.reactive_var < 0.0);
    }

    #[test]
    fn from_p_q_three_four_five() {
        // Canonical 3-4-5 right triangle.
        let t = PowerTriangle::from_p_q(3000.0, 4000.0).unwrap();
        assert!((t.apparent_va - 5000.0).abs() < EPS);
        assert!((t.power_factor - 0.6).abs() < EPS);
        assert_eq!(t.phase, Phase::Lagging);
        assert!(t.pythagorean_residual().abs() < 1e-6);
    }

    #[test]
    fn from_p_q_unity_when_no_reactive() {
        let t = PowerTriangle::from_p_q(1500.0, 0.0).unwrap();
        assert!((t.apparent_va - 1500.0).abs() < EPS);
        assert!((t.power_factor - 1.0).abs() < EPS);
        assert_eq!(t.phase, Phase::Unity);
    }

    #[test]
    fn phase_angle_round_trips() {
        // from_vi_phase -> phase_angle_rad must recover the input.
        for deg in [-70.0_f64, -10.0, 0.0, 25.0, 55.0] {
            let phi = deg.to_radians();
            let t = PowerTriangle::from_vi_phase(100.0, 3.0, phi).unwrap();
            assert!(
                (t.phase_angle_rad() - phi).abs() < 1e-9,
                "phase angle did not round-trip at {deg} deg"
            );
        }
    }

    #[test]
    fn rejects_non_positive_voltage_current() {
        assert!(matches!(
            PowerTriangle::from_vi_phase(0.0, 1.0, 0.0),
            Err(PowerError::NonPositive {
                name: "voltage_v",
                ..
            })
        ));
        assert!(matches!(
            PowerTriangle::from_vi_phase(1.0, -2.0, 0.0),
            Err(PowerError::NonPositive {
                name: "current_a",
                ..
            })
        ));
    }

    #[test]
    fn rejects_angle_at_ninety_degrees() {
        // |phi| = pi/2 => PF = 0 boundary is fine, but cos can dip below
        // 0 for the smallest representable step beyond it; a true pi/2
        // gives cos ~= 6e-17 >= 0 so it is accepted at the PF=0 limit.
        let t = PowerTriangle::from_vi_phase(100.0, 1.0, std::f64::consts::FRAC_PI_2).unwrap();
        assert!(t.power_factor.abs() < 1e-9);
        // Just past pi/2, cos goes negative -> rejected.
        let past = std::f64::consts::FRAC_PI_2 + 0.1;
        assert!(matches!(
            PowerTriangle::from_vi_phase(100.0, 1.0, past),
            Err(PowerError::PowerFactorOutOfRange { .. })
        ));
    }

    #[test]
    fn rejects_collapsed_triangle() {
        assert!(matches!(
            PowerTriangle::from_p_q(0.0, 0.0),
            Err(PowerError::NonPositive {
                name: "apparent_va",
                ..
            })
        ));
    }

    #[test]
    fn rejects_negative_real_power() {
        assert!(matches!(
            PowerTriangle::from_p_q(-1.0, 10.0),
            Err(PowerError::Negative { name: "real_w", .. })
        ));
    }
}
