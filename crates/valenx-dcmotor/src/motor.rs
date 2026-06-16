//! Steady-state brushed permanent-magnet DC-motor model.
//!
//! All relations are the textbook armature equations evaluated in
//! steady state (the armature-inductance term `L di/dt` is zero, so it
//! drops out). SI units throughout unless a method name says otherwise:
//! speed `omega` in rad/s, current `i` in A, torque `T` in N*m, voltage
//! `V` in V, and power in W.
//!
//! In a coherent SI motor the back-EMF constant `Ke` (V*s/rad) and the
//! torque constant `Kt` (N*m/A) are numerically equal; the model keeps
//! them as separate fields so a caller may supply catalogue values in
//! mixed unit systems, but the analytic identities below assume the
//! ideal lossless coupling `T*omega = E*i` only when `Kt == Ke`.

use serde::{Deserialize, Serialize};

use crate::error::DcMotorError;

/// Parameters of an ideal steady-state brushed PM DC motor.
///
/// Construct with [`DcMotor::new`], which validates every field, rather
/// than building the struct literally — the analytic methods assume the
/// invariants the constructor enforces (`r > 0`, `ke > 0`, `kt > 0`).
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct DcMotor {
    /// Armature (terminal-to-terminal) resistance `R`, ohms. `R > 0`.
    pub resistance_ohm: f64,
    /// Back-EMF constant `Ke`, V*s/rad (volts per rad/s). `Ke > 0`.
    pub ke_v_s_per_rad: f64,
    /// Torque constant `Kt`, N*m/A. `Kt > 0`.
    pub kt_nm_per_a: f64,
}

impl DcMotor {
    /// Build a validated motor.
    ///
    /// # Errors
    ///
    /// Returns [`DcMotorError::NotFinite`] if any argument is `NaN` or
    /// infinite, or [`DcMotorError::OutOfRange`] if any argument is not
    /// strictly positive.
    pub fn new(
        resistance_ohm: f64,
        ke_v_s_per_rad: f64,
        kt_nm_per_a: f64,
    ) -> Result<Self, DcMotorError> {
        DcMotorError::require_positive("resistance_ohm", resistance_ohm)?;
        DcMotorError::require_positive("ke_v_s_per_rad", ke_v_s_per_rad)?;
        DcMotorError::require_positive("kt_nm_per_a", kt_nm_per_a)?;
        Ok(Self {
            resistance_ohm,
            ke_v_s_per_rad,
            kt_nm_per_a,
        })
    }

    /// Build the common coherent-SI motor where `Kt == Ke == k`.
    ///
    /// # Errors
    ///
    /// As [`DcMotor::new`], with `k` validated for both constants.
    pub fn from_constant(resistance_ohm: f64, k: f64) -> Result<Self, DcMotorError> {
        DcMotor::new(resistance_ohm, k, k)
    }

    /// Back-EMF `E = Ke * omega` (volts) at shaft speed `omega` (rad/s).
    ///
    /// `omega` may be negative (reverse rotation); the EMF follows sign.
    pub fn back_emf(&self, omega_rad_s: f64) -> f64 {
        self.ke_v_s_per_rad * omega_rad_s
    }

    /// Electromagnetic torque `T = Kt * I` (N*m) at armature current
    /// `i` (A). Sign follows the current.
    pub fn torque(&self, current_a: f64) -> f64 {
        self.kt_nm_per_a * current_a
    }

    /// Terminal voltage required to hold `(i, omega)` in steady state:
    /// `V = I*R + Ke*omega`.
    pub fn terminal_voltage(&self, current_a: f64, omega_rad_s: f64) -> f64 {
        current_a * self.resistance_ohm + self.back_emf(omega_rad_s)
    }

    /// Armature current drawn at terminal voltage `v` and speed
    /// `omega`, from `I = (V - Ke*omega) / R`.
    pub fn current(&self, voltage_v: f64, omega_rad_s: f64) -> f64 {
        (voltage_v - self.back_emf(omega_rad_s)) / self.resistance_ohm
    }

    /// Stall current `I = V / R` (rotor held, `omega = 0`), A.
    ///
    /// # Errors
    ///
    /// [`DcMotorError::NotFinite`] if `voltage_v` is not finite.
    pub fn stall_current(&self, voltage_v: f64) -> Result<f64, DcMotorError> {
        DcMotorError::require_finite("voltage_v", voltage_v)?;
        Ok(voltage_v / self.resistance_ohm)
    }

    /// Stall torque `T = Kt * V / R` (rotor held), N*m.
    ///
    /// # Errors
    ///
    /// [`DcMotorError::NotFinite`] if `voltage_v` is not finite.
    pub fn stall_torque(&self, voltage_v: f64) -> Result<f64, DcMotorError> {
        Ok(self.kt_nm_per_a * self.stall_current(voltage_v)?)
    }

    /// Ideal no-load speed `omega = V / Ke`, rad/s.
    ///
    /// This is the speed at which the torque-producing current is zero
    /// (`T = 0`), found by setting `I = 0` in the voltage balance so
    /// `V = Ke*omega`. A real motor settles slightly below this because
    /// friction draws a small no-load current; the ideal model ignores
    /// that (see crate-level honest-scope note).
    ///
    /// # Errors
    ///
    /// [`DcMotorError::NotFinite`] if `voltage_v` is not finite.
    pub fn no_load_speed(&self, voltage_v: f64) -> Result<f64, DcMotorError> {
        DcMotorError::require_finite("voltage_v", voltage_v)?;
        Ok(voltage_v / self.ke_v_s_per_rad)
    }

    /// Speed on the linear torque-speed line for a demanded shaft
    /// torque `t` at terminal voltage `v`.
    ///
    /// Eliminating `I` between `T = Kt*I` and `V = I*R + Ke*omega`
    /// gives the straight line
    ///
    /// `omega = V/Ke - (R / (Kt*Ke)) * T`
    ///
    /// which runs from the no-load point `(T = 0, omega = V/Ke)` down to
    /// the stall point `(T = T_stall, omega = 0)`.
    ///
    /// # Errors
    ///
    /// [`DcMotorError::NotFinite`] if `voltage_v` or `torque_nm` is not
    /// finite.
    pub fn speed_at_torque(&self, voltage_v: f64, torque_nm: f64) -> Result<f64, DcMotorError> {
        DcMotorError::require_finite("voltage_v", voltage_v)?;
        DcMotorError::require_finite("torque_nm", torque_nm)?;
        let slope = self.resistance_ohm / (self.kt_nm_per_a * self.ke_v_s_per_rad);
        Ok(self.no_load_speed(voltage_v)? - slope * torque_nm)
    }

    /// Inverse of [`DcMotor::speed_at_torque`]: the shaft torque produced
    /// at terminal voltage `v` and speed `omega`, equal to
    /// `Kt * (V - Ke*omega) / R`.
    ///
    /// # Errors
    ///
    /// [`DcMotorError::NotFinite`] if `voltage_v` or `omega_rad_s` is not
    /// finite.
    pub fn torque_at_speed(&self, voltage_v: f64, omega_rad_s: f64) -> Result<f64, DcMotorError> {
        DcMotorError::require_finite("voltage_v", voltage_v)?;
        DcMotorError::require_finite("omega_rad_s", omega_rad_s)?;
        Ok(self.torque(self.current(voltage_v, omega_rad_s)))
    }

    /// A full steady-state operating point at terminal voltage `v` and
    /// speed `omega`, packaging current, torque, the power split and
    /// efficiency.
    ///
    /// # Errors
    ///
    /// [`DcMotorError::NotFinite`] if `voltage_v` or `omega_rad_s` is not
    /// finite.
    pub fn operating_point(
        &self,
        voltage_v: f64,
        omega_rad_s: f64,
    ) -> Result<OperatingPoint, DcMotorError> {
        DcMotorError::require_finite("voltage_v", voltage_v)?;
        DcMotorError::require_finite("omega_rad_s", omega_rad_s)?;
        let current_a = self.current(voltage_v, omega_rad_s);
        let torque_nm = self.torque(current_a);
        let electrical_power_w = voltage_v * current_a;
        let mechanical_power_w = torque_nm * omega_rad_s;
        let copper_loss_w = current_a * current_a * self.resistance_ohm;
        let efficiency = if electrical_power_w.abs() > 0.0 {
            mechanical_power_w / electrical_power_w
        } else {
            0.0
        };
        Ok(OperatingPoint {
            voltage_v,
            omega_rad_s,
            current_a,
            torque_nm,
            electrical_power_w,
            mechanical_power_w,
            copper_loss_w,
            efficiency,
        })
    }

    /// Shaft speed at which the motor delivers its maximum mechanical
    /// output power: exactly half the no-load speed, `omega = V/(2*Ke)`,
    /// rad/s.
    ///
    /// On the straight torque-speed line the output power `P = T*omega`
    /// is a downward parabola, so it peaks at the midpoint of the line —
    /// half the [`no_load_speed`](DcMotor::no_load_speed) and half the
    /// stall torque.
    ///
    /// # Errors
    ///
    /// [`DcMotorError::NotFinite`] if `voltage_v` is not finite.
    pub fn max_power_speed(&self, voltage_v: f64) -> Result<f64, DcMotorError> {
        DcMotorError::require_finite("voltage_v", voltage_v)?;
        Ok(voltage_v / (2.0 * self.ke_v_s_per_rad))
    }

    /// Shaft torque at the maximum-power operating point: half the stall
    /// torque, `T = Kt*V/(2*R)`, N*m.
    ///
    /// # Errors
    ///
    /// [`DcMotorError::NotFinite`] if `voltage_v` is not finite.
    pub fn max_power_torque(&self, voltage_v: f64) -> Result<f64, DcMotorError> {
        DcMotorError::require_finite("voltage_v", voltage_v)?;
        Ok(self.kt_nm_per_a * voltage_v / (2.0 * self.resistance_ohm))
    }

    /// Maximum mechanical output power at terminal voltage `v`, watts.
    ///
    /// It occurs at the midpoint of the torque-speed line, where
    /// `T = T_stall/2` and `omega = omega_nl/2`, giving
    ///
    /// `P_max = T_stall * omega_nl / 4 = Kt * V^2 / (4 * R * Ke)`.
    ///
    /// For a coherent motor (`Kt == Ke`) this reduces to the
    /// maximum-power-transfer value `V^2 / (4 R)`, drawn at half the
    /// stall current with the back-EMF equal to half the terminal
    /// voltage (so the ideal motor is exactly 50% efficient there).
    ///
    /// # Errors
    ///
    /// [`DcMotorError::NotFinite`] if `voltage_v` is not finite.
    pub fn max_output_power(&self, voltage_v: f64) -> Result<f64, DcMotorError> {
        DcMotorError::require_finite("voltage_v", voltage_v)?;
        Ok(self.kt_nm_per_a * voltage_v * voltage_v
            / (4.0 * self.resistance_ohm * self.ke_v_s_per_rad))
    }

    /// The full steady-state [`OperatingPoint`] at the maximum-power
    /// speed.
    ///
    /// Convenience wrapper that evaluates
    /// [`operating_point`](DcMotor::operating_point) at
    /// [`max_power_speed`](DcMotor::max_power_speed); its
    /// `mechanical_power_w` equals
    /// [`max_output_power`](DcMotor::max_output_power), and for a
    /// coherent motor its efficiency is exactly `0.5`.
    ///
    /// # Errors
    ///
    /// [`DcMotorError::NotFinite`] if `voltage_v` is not finite.
    pub fn max_power_point(&self, voltage_v: f64) -> Result<OperatingPoint, DcMotorError> {
        let omega = self.max_power_speed(voltage_v)?;
        self.operating_point(voltage_v, omega)
    }
}

/// A solved steady-state operating point.
///
/// Produced by [`DcMotor::operating_point`]. All powers are in watts.
/// The model's only loss is armature copper loss `I^2 R`; with `Kt == Ke`
/// the identity `P_electrical = P_mechanical + I^2 R` holds exactly
/// (see [`OperatingPoint::power_residual`]).
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct OperatingPoint {
    /// Terminal voltage `V` applied, volts.
    pub voltage_v: f64,
    /// Shaft speed `omega`, rad/s.
    pub omega_rad_s: f64,
    /// Armature current `I`, amperes.
    pub current_a: f64,
    /// Shaft torque `T = Kt*I`, N*m.
    pub torque_nm: f64,
    /// Electrical input power `V*I`, watts.
    pub electrical_power_w: f64,
    /// Mechanical output power `T*omega`, watts.
    pub mechanical_power_w: f64,
    /// Armature copper loss `I^2 * R`, watts.
    pub copper_loss_w: f64,
    /// Efficiency `T*omega / (V*I)` (dimensionless). Zero when the
    /// electrical power is exactly zero.
    pub efficiency: f64,
}

impl OperatingPoint {
    /// Power-balance residual `P_electrical - P_mechanical - I^2 R`.
    ///
    /// For an ideal motor with `Kt == Ke` this is zero to floating-point
    /// round-off, because the back-EMF power `E*I = Ke*omega*I` equals
    /// the mechanical power `Kt*I*omega`. A non-zero residual flags that
    /// the supplied `Kt` and `Ke` were not equal.
    pub fn power_residual(&self) -> f64 {
        self.electrical_power_w - self.mechanical_power_w - self.copper_loss_w
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::error::ErrorCategory;

    const EPS: f64 = 1e-12;

    /// A simple coherent-SI motor: R = 2 ohm, Ke = Kt = 0.05.
    fn motor() -> DcMotor {
        DcMotor::from_constant(2.0, 0.05).unwrap()
    }

    #[test]
    fn constructor_rejects_non_positive_and_non_finite() {
        assert_eq!(
            DcMotor::new(0.0, 0.05, 0.05).unwrap_err().code(),
            "dcmotor.out_of_range"
        );
        assert_eq!(
            DcMotor::new(2.0, -0.05, 0.05).unwrap_err().code(),
            "dcmotor.out_of_range"
        );
        assert_eq!(
            DcMotor::new(2.0, 0.05, f64::NAN).unwrap_err().code(),
            "dcmotor.not_finite"
        );
        assert_eq!(
            DcMotor::new(f64::INFINITY, 0.05, 0.05)
                .unwrap_err()
                .category(),
            ErrorCategory::Domain
        );
        assert_eq!(
            DcMotor::new(-1.0, 0.05, 0.05).unwrap_err().category(),
            ErrorCategory::Range
        );
    }

    #[test]
    fn from_constant_sets_both_constants_equal() {
        let m = DcMotor::from_constant(2.0, 0.05).unwrap();
        assert!((m.ke_v_s_per_rad - m.kt_nm_per_a).abs() < EPS);
        assert!((m.ke_v_s_per_rad - 0.05).abs() < EPS);
    }

    #[test]
    fn back_emf_is_ke_times_omega() {
        let m = motor();
        // Ke = 0.05, omega = 100 -> E = 5.0 V.
        assert!((m.back_emf(100.0) - 5.0).abs() < EPS);
        // Linear and odd in omega.
        assert!((m.back_emf(0.0) - 0.0).abs() < EPS);
        assert!((m.back_emf(-40.0) - (-2.0)).abs() < EPS);
        // Doubling omega doubles E.
        assert!((m.back_emf(200.0) - 2.0 * m.back_emf(100.0)).abs() < EPS);
    }

    #[test]
    fn torque_is_kt_times_current() {
        let m = motor();
        // Kt = 0.05, I = 3 -> T = 0.15 N*m.
        assert!((m.torque(3.0) - 0.15).abs() < EPS);
        assert!((m.torque(0.0) - 0.0).abs() < EPS);
        assert!((m.torque(-2.0) - (-0.1)).abs() < EPS);
    }

    #[test]
    fn voltage_balance_v_equals_ir_plus_emf() {
        let m = motor();
        // I = 3 A, omega = 100 rad/s, R = 2:
        // V = 3*2 + 0.05*100 = 6 + 5 = 11 V.
        let v = m.terminal_voltage(3.0, 100.0);
        assert!((v - 11.0).abs() < EPS);
        // Check the decomposition exactly equals I*R + E.
        let expected = 3.0 * m.resistance_ohm + m.back_emf(100.0);
        assert!((v - expected).abs() < EPS);
    }

    #[test]
    fn current_is_inverse_of_voltage_balance() {
        let m = motor();
        // Round-trip: pick (I, omega), compute V, recover I.
        let i = 4.5;
        let omega = 75.0;
        let v = m.terminal_voltage(i, omega);
        let i_back = m.current(v, omega);
        assert!((i_back - i).abs() < 1e-10);
    }

    #[test]
    fn stall_current_is_v_over_r() {
        let m = motor();
        // V = 12, R = 2 -> 6 A.
        assert!((m.stall_current(12.0).unwrap() - 6.0).abs() < EPS);
        // At stall the back-EMF is zero, so current(V, 0) must agree.
        assert!((m.current(12.0, 0.0) - 6.0).abs() < EPS);
    }

    #[test]
    fn stall_torque_is_kt_v_over_r() {
        let m = motor();
        // Kt * V / R = 0.05 * 12 / 2 = 0.3 N*m.
        assert!((m.stall_torque(12.0).unwrap() - 0.3).abs() < EPS);
        // Consistent with torque(stall_current).
        let t = m.torque(m.stall_current(12.0).unwrap());
        assert!((m.stall_torque(12.0).unwrap() - t).abs() < EPS);
    }

    #[test]
    fn no_load_speed_is_v_over_ke_and_gives_zero_torque() {
        let m = motor();
        // V = 12, Ke = 0.05 -> 240 rad/s.
        let omega_nl = m.no_load_speed(12.0).unwrap();
        assert!((omega_nl - 240.0).abs() < EPS);
        // At no-load speed the current (hence torque) is zero.
        assert!(m.current(12.0, omega_nl).abs() < 1e-10);
        assert!(m.torque_at_speed(12.0, omega_nl).unwrap().abs() < 1e-10);
    }

    #[test]
    fn torque_speed_line_passes_through_both_endpoints() {
        let m = motor();
        let v = 12.0;
        let stall_t = m.stall_torque(v).unwrap();
        let omega_nl = m.no_load_speed(v).unwrap();
        // At T = 0 the line is at the no-load speed.
        assert!((m.speed_at_torque(v, 0.0).unwrap() - omega_nl).abs() < 1e-10);
        // At T = T_stall the line is at zero speed.
        assert!(m.speed_at_torque(v, stall_t).unwrap().abs() < 1e-10);
    }

    #[test]
    fn torque_speed_line_is_straight() {
        let m = motor();
        let v = 12.0;
        // Equal torque steps must give equal, opposite speed steps.
        let w0 = m.speed_at_torque(v, 0.0).unwrap();
        let w1 = m.speed_at_torque(v, 0.1).unwrap();
        let w2 = m.speed_at_torque(v, 0.2).unwrap();
        let d01 = w0 - w1;
        let d12 = w1 - w2;
        assert!((d01 - d12).abs() < 1e-10);
        // Slope equals R / (Kt*Ke): d(omega)/dT.
        let slope = m.resistance_ohm / (m.kt_nm_per_a * m.ke_v_s_per_rad);
        assert!((d01 / 0.1 - slope).abs() < 1e-9);
    }

    #[test]
    fn speed_and_torque_maps_are_mutual_inverses() {
        let m = motor();
        let v = 12.0;
        // omega -> T -> omega.
        let omega = 150.0;
        let t = m.torque_at_speed(v, omega).unwrap();
        let omega_back = m.speed_at_torque(v, t).unwrap();
        assert!((omega_back - omega).abs() < 1e-9);
    }

    #[test]
    fn operating_point_obeys_power_balance_for_ideal_motor() {
        let m = motor();
        let op = m.operating_point(12.0, 150.0).unwrap();
        // I = (12 - 0.05*150) / 2 = (12 - 7.5)/2 = 2.25 A.
        assert!((op.current_a - 2.25).abs() < 1e-12);
        // T = 0.05 * 2.25 = 0.1125 N*m.
        assert!((op.torque_nm - 0.1125).abs() < 1e-12);
        // P_in = 12 * 2.25 = 27 W.
        assert!((op.electrical_power_w - 27.0).abs() < 1e-12);
        // P_mech = 0.1125 * 150 = 16.875 W.
        assert!((op.mechanical_power_w - 16.875).abs() < 1e-12);
        // I^2 R = 2.25^2 * 2 = 10.125 W.
        assert!((op.copper_loss_w - 10.125).abs() < 1e-12);
        // Balance residual is zero for Kt == Ke.
        assert!(op.power_residual().abs() < 1e-12);
        // Efficiency = P_mech / P_in = 16.875 / 27 = 0.625.
        assert!((op.efficiency - 0.625).abs() < 1e-12);
    }

    #[test]
    fn efficiency_is_zero_at_stall_and_at_no_load() {
        let m = motor();
        let v = 12.0;
        // Stall: omega = 0 -> P_mech = 0 -> eta = 0.
        let stall = m.operating_point(v, 0.0).unwrap();
        assert!(stall.efficiency.abs() < EPS);
        assert!(stall.mechanical_power_w.abs() < EPS);
        // No-load: I = 0 -> P_in = 0 -> eta defined as 0.
        let omega_nl = m.no_load_speed(v).unwrap();
        let nl = m.operating_point(v, omega_nl).unwrap();
        assert!(nl.current_a.abs() < 1e-10);
        assert!(nl.efficiency.abs() < EPS);
    }

    #[test]
    fn efficiency_stays_in_unit_interval_across_the_line() {
        let m = motor();
        let v = 12.0;
        let omega_nl = m.no_load_speed(v).unwrap();
        // Sample the operating line between stall and no-load.
        for k in 0..=20 {
            let omega = omega_nl * (k as f64 / 20.0);
            let op = m.operating_point(v, omega).unwrap();
            assert!(
                op.efficiency >= -EPS && op.efficiency <= 1.0 + EPS,
                "efficiency {eff} out of [0,1] at omega {omega}",
                eff = op.efficiency
            );
        }
    }

    #[test]
    fn unequal_kt_ke_breaks_power_balance() {
        // A non-coherent motor: Kt != Ke. The back-EMF power no longer
        // equals the mechanical power, so the residual is non-zero.
        let m = DcMotor::new(2.0, 0.05, 0.06).unwrap();
        let op = m.operating_point(12.0, 150.0).unwrap();
        assert!(op.power_residual().abs() > 1e-6);
    }

    #[test]
    fn finite_guards_reject_nan_inputs() {
        let m = motor();
        assert_eq!(
            m.stall_current(f64::NAN).unwrap_err().code(),
            "dcmotor.not_finite"
        );
        assert!(m.stall_torque(f64::INFINITY).is_err());
        assert!(m.no_load_speed(f64::NAN).is_err());
        assert!(m.speed_at_torque(12.0, f64::NAN).is_err());
        assert!(m.torque_at_speed(f64::NAN, 10.0).is_err());
        assert!(m.operating_point(12.0, f64::NAN).is_err());
    }

    #[test]
    fn serde_round_trips_via_json() {
        let m = motor();
        let json = serde_json::to_string(&m).unwrap();
        let back: DcMotor = serde_json::from_str(&json).unwrap();
        assert_eq!(m, back);
    }

    #[test]
    fn max_power_is_quarter_stall_torque_times_no_load_speed() {
        let m = motor();
        let v = 12.0;
        let p_max = m.max_output_power(v).unwrap();
        let stall_t = m.stall_torque(v).unwrap();
        let omega_nl = m.no_load_speed(v).unwrap();
        // P_max = T_stall * omega_nl / 4.
        assert!(
            (p_max - stall_t * omega_nl / 4.0).abs() < EPS,
            "got {p_max}"
        );
        // For this motor: 0.3 * 240 / 4 = 18 W.
        assert!((p_max - 18.0).abs() < 1e-9, "got {p_max}");
    }

    #[test]
    fn max_power_equals_v_squared_over_4r_for_coherent_motor() {
        // Maximum-power-transfer theorem: with Kt == Ke, P_max = V^2/(4R).
        let m = motor(); // Kt == Ke
        let v = 12.0;
        let p_max = m.max_output_power(v).unwrap();
        assert!(
            (p_max - v * v / (4.0 * m.resistance_ohm)).abs() < EPS,
            "got {p_max}"
        );
    }

    #[test]
    fn max_power_occurs_at_line_midpoint() {
        let m = motor();
        let v = 12.0;
        // Speed = half no-load, torque = half stall.
        let w_mp = m.max_power_speed(v).unwrap();
        let t_mp = m.max_power_torque(v).unwrap();
        assert!((w_mp - m.no_load_speed(v).unwrap() / 2.0).abs() < EPS);
        assert!((t_mp - m.stall_torque(v).unwrap() / 2.0).abs() < EPS);
        // And the (T, omega) pair lies on the torque-speed line.
        assert!((m.speed_at_torque(v, t_mp).unwrap() - w_mp).abs() < 1e-9);
    }

    #[test]
    fn max_power_point_matches_and_is_half_efficient() {
        let m = motor(); // coherent
        let v = 12.0;
        let op = m.max_power_point(v).unwrap();
        // Its mechanical power equals max_output_power.
        assert!((op.mechanical_power_w - m.max_output_power(v).unwrap()).abs() < 1e-9);
        // Speed and torque are the midpoint values.
        assert!((op.omega_rad_s - m.max_power_speed(v).unwrap()).abs() < EPS);
        assert!((op.torque_nm - m.max_power_torque(v).unwrap()).abs() < EPS);
        // At maximum power transfer the ideal motor is exactly 50% efficient.
        assert!((op.efficiency - 0.5).abs() < 1e-12, "eff {}", op.efficiency);
    }

    #[test]
    fn max_power_is_the_true_maximum_along_the_line() {
        let m = motor();
        let v = 12.0;
        let p_max = m.max_output_power(v).unwrap();
        let omega_nl = m.no_load_speed(v).unwrap();
        // No sampled point on the line may exceed the claimed peak power.
        for k in 0..=100 {
            let omega = omega_nl * (k as f64 / 100.0);
            let p = m.operating_point(v, omega).unwrap().mechanical_power_w;
            assert!(
                p <= p_max + 1e-9,
                "p {p} exceeds p_max {p_max} at omega {omega}"
            );
        }
    }

    #[test]
    fn max_power_methods_reject_non_finite_voltage() {
        let m = motor();
        assert!(m.max_power_speed(f64::NAN).is_err());
        assert!(m.max_power_torque(f64::INFINITY).is_err());
        assert!(m.max_output_power(f64::NAN).is_err());
        assert!(m.max_power_point(f64::NAN).is_err());
    }
}
