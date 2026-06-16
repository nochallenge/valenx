//! Steady-state rotating-field kinematics of a poly-phase induction motor.
//!
//! All quantities are the standard textbook closed-form relations for
//! the rotating magnetic field of a three-phase asynchronous machine.
//! Speeds are in revolutions per minute (rev/min), frequencies in
//! hertz (Hz), and slip is a dimensionless fraction.

use serde::{Deserialize, Serialize};

use crate::error::InductionMotorError;

/// Numerator constant in the synchronous-speed relation
/// `Ns = 120 f / poles`.
///
/// It is `60 s/min` (Hz to rev/min) times `2` (two poles per pole
/// pair), i.e. `Ns = (60 * f) / pole_pairs = 120 f / poles`.
pub const SYNC_SPEED_CONSTANT: f64 = 120.0;

/// Synchronous speed of the rotating stator field, in rev/min.
///
/// ```text
/// Ns = 120 * f / poles
/// ```
///
/// `f` is the line (supply) frequency in hertz and `poles` is the
/// number of magnetic poles (an even integer).
///
/// Returns [`InductionMotorError::InvalidFrequency`] if `freq_hz` is
/// not finite and `> 0`, or [`InductionMotorError::InvalidPoles`] if
/// `poles` is zero or odd.
///
/// # Examples
///
/// ```
/// use valenx_inductionmotor::sync_speed_rpm;
/// // 4-pole machine on a 60 Hz line spins its field at 1800 rev/min.
/// let ns = sync_speed_rpm(60.0, 4).unwrap();
/// assert!((ns - 1800.0).abs() < 1e-9);
/// ```
pub fn sync_speed_rpm(freq_hz: f64, poles: u32) -> Result<f64, InductionMotorError> {
    validate_frequency(freq_hz)?;
    validate_poles(poles)?;
    Ok(SYNC_SPEED_CONSTANT * freq_hz / poles as f64)
}

/// Fractional slip of a rotor turning at `rotor_rpm` against a field
/// turning at `sync_rpm`.
///
/// ```text
/// s = (Ns - Nr) / Ns
/// ```
///
/// The result is dimensionless: `0` at synchronous speed and `1` at
/// standstill (`Nr = 0`). This is the *raw* arithmetic slip and may
/// fall outside `[0, 1]` for hyper-synchronous (generating) or
/// reverse (plugging) operation; use [`InductionMotor::slip`] for the
/// range-checked motoring value.
///
/// Returns [`InductionMotorError::InvalidFrequency`] if `sync_rpm` is
/// not finite and `> 0` (a non-positive synchronous speed has no
/// physical meaning and would divide by zero), or
/// [`InductionMotorError::InvalidRotorSpeed`] if `rotor_rpm` is not
/// finite.
///
/// # Examples
///
/// ```
/// use valenx_inductionmotor::slip;
/// // Field at 1800 rev/min, rotor at 1746 rev/min -> 3 % slip.
/// let s = slip(1800.0, 1746.0).unwrap();
/// assert!((s - 0.03).abs() < 1e-12);
/// ```
pub fn slip(sync_rpm: f64, rotor_rpm: f64) -> Result<f64, InductionMotorError> {
    if !sync_rpm.is_finite() || sync_rpm <= 0.0 {
        // Reuse the frequency variant: a synchronous speed is just a
        // scaled frequency, and both must be finite and strictly
        // positive for the slip ratio to be defined.
        return Err(InductionMotorError::InvalidFrequency { hz: sync_rpm });
    }
    if !rotor_rpm.is_finite() {
        return Err(InductionMotorError::InvalidRotorSpeed { rpm: rotor_rpm });
    }
    Ok((sync_rpm - rotor_rpm) / sync_rpm)
}

/// Rotor electrical frequency (slip frequency), in hertz.
///
/// ```text
/// f_r = s * f
/// ```
///
/// The frequency of the currents induced in the rotor bars equals the
/// slip times the supply frequency: `f` at standstill (`s = 1`) and
/// `0` at synchronous speed (`s = 0`).
///
/// `slip` must be finite; `freq_hz` must be finite and `> 0`.
///
/// # Examples
///
/// ```
/// use valenx_inductionmotor::rotor_frequency_hz;
/// // 3 % slip on a 60 Hz supply -> 1.8 Hz rotor frequency.
/// let fr = rotor_frequency_hz(0.03, 60.0).unwrap();
/// assert!((fr - 1.8).abs() < 1e-12);
/// ```
pub fn rotor_frequency_hz(slip: f64, freq_hz: f64) -> Result<f64, InductionMotorError> {
    if !slip.is_finite() {
        return Err(InductionMotorError::InvalidSlip { slip });
    }
    validate_frequency(freq_hz)?;
    Ok(slip * freq_hz)
}

/// Rotor mechanical speed implied by a synchronous speed and slip, in
/// rev/min.
///
/// ```text
/// Nr = Ns * (1 - s)
/// ```
///
/// This is the algebraic inverse of [`slip`]: at `s = 0` the rotor
/// runs at the field speed `Ns`, and at `s = 1` it is at rest.
///
/// `slip` must be finite; `sync_rpm` must be finite and `> 0`.
///
/// # Examples
///
/// ```
/// use valenx_inductionmotor::rotor_speed_rpm;
/// // 1800 rev/min field at 3 % slip -> 1746 rev/min rotor.
/// let nr = rotor_speed_rpm(1800.0, 0.03).unwrap();
/// assert!((nr - 1746.0).abs() < 1e-9);
/// ```
pub fn rotor_speed_rpm(sync_rpm: f64, slip: f64) -> Result<f64, InductionMotorError> {
    if !sync_rpm.is_finite() || sync_rpm <= 0.0 {
        return Err(InductionMotorError::InvalidFrequency { hz: sync_rpm });
    }
    if !slip.is_finite() {
        return Err(InductionMotorError::InvalidSlip { slip });
    }
    Ok(sync_rpm * (1.0 - slip))
}

/// Rotor copper (`I^2 R`) loss as a fraction of the air-gap power:
/// `P_rotor_cu = s * P_airgap`.
///
/// The air-gap power `P_airgap` is the power transferred across the air gap
/// into the rotor; energy conservation splits it into the rotor copper loss
/// `s * P_airgap` and the developed mechanical power `(1 - s) * P_airgap`
/// ([`developed_mechanical_power_w`]). This split is a pure consequence of
/// the slip — it needs no equivalent-circuit parameters — and is the
/// relation that gives slip its energetic meaning: at standstill (`s = 1`)
/// all air-gap power is lost in the rotor; at synchronous speed (`s = 0`)
/// none is.
///
/// # Errors
///
/// Returns [`InductionMotorError::InvalidSlip`] if `slip` is not finite or
/// outside the motoring range `[0, 1]`, or
/// [`InductionMotorError::InvalidPower`] if `air_gap_power_w` is not finite
/// or negative.
///
/// # Examples
///
/// ```
/// use valenx_inductionmotor::rotor_copper_loss_w;
/// // 10 kW across the gap at 3 % slip wastes 300 W in the rotor.
/// let p_cu = rotor_copper_loss_w(0.03, 10_000.0).unwrap();
/// assert!((p_cu - 300.0).abs() < 1e-9);
/// ```
pub fn rotor_copper_loss_w(slip: f64, air_gap_power_w: f64) -> Result<f64, InductionMotorError> {
    let s = validate_motoring_slip(slip)?;
    let pag = validate_power(air_gap_power_w)?;
    Ok(s * pag)
}

/// Developed (gross) mechanical power as a fraction of the air-gap power:
/// `P_mech = (1 - s) * P_airgap`.
///
/// The complement of [`rotor_copper_loss_w`] in the air-gap power split
/// `P_airgap = P_rotor_cu + P_mech`. This is the *developed* (gross)
/// mechanical power — before friction and windage — and is not the shaft
/// output; it does not account for stator or iron losses.
///
/// # Errors
///
/// Returns [`InductionMotorError::InvalidSlip`] if `slip` is not finite or
/// outside `[0, 1]`, or [`InductionMotorError::InvalidPower`] if
/// `air_gap_power_w` is not finite or negative.
///
/// # Examples
///
/// ```
/// use valenx_inductionmotor::developed_mechanical_power_w;
/// // 10 kW across the gap at 3 % slip develops 9.7 kW of mechanical power.
/// let p_mech = developed_mechanical_power_w(0.03, 10_000.0).unwrap();
/// assert!((p_mech - 9_700.0).abs() < 1e-9);
/// ```
pub fn developed_mechanical_power_w(
    slip: f64,
    air_gap_power_w: f64,
) -> Result<f64, InductionMotorError> {
    let s = validate_motoring_slip(slip)?;
    let pag = validate_power(air_gap_power_w)?;
    Ok((1.0 - s) * pag)
}

fn validate_frequency(freq_hz: f64) -> Result<(), InductionMotorError> {
    if !freq_hz.is_finite() || freq_hz <= 0.0 {
        return Err(InductionMotorError::InvalidFrequency { hz: freq_hz });
    }
    Ok(())
}

fn validate_poles(poles: u32) -> Result<(), InductionMotorError> {
    if poles == 0 || poles % 2 != 0 {
        return Err(InductionMotorError::InvalidPoles { poles });
    }
    Ok(())
}

/// Validate a motoring slip: finite and within `[0, 1]`.
fn validate_motoring_slip(slip: f64) -> Result<f64, InductionMotorError> {
    if !slip.is_finite() || !(0.0..=1.0).contains(&slip) {
        return Err(InductionMotorError::InvalidSlip { slip });
    }
    Ok(slip)
}

/// Validate an air-gap power: finite and non-negative.
fn validate_power(air_gap_power_w: f64) -> Result<f64, InductionMotorError> {
    if !air_gap_power_w.is_finite() || air_gap_power_w < 0.0 {
        return Err(InductionMotorError::InvalidPower {
            watts: air_gap_power_w,
        });
    }
    Ok(air_gap_power_w)
}

/// A fully specified three-phase induction-motor operating point.
///
/// Construct one with [`InductionMotor::new`] (validated against the
/// motoring slip range) or [`InductionMotor::from_slip`]. Once built,
/// every derived kinematic quantity is available infallibly through the
/// accessor methods, since all inputs were validated at construction.
///
/// Speeds are in rev/min, the supply frequency in Hz, and slip is a
/// dimensionless fraction in `[0, 1]`.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct InductionMotor {
    /// Number of magnetic poles (a positive even integer).
    poles: u32,
    /// Line / supply frequency in hertz (finite, `> 0`).
    freq_hz: f64,
    /// Synchronous speed in rev/min (`120 f / poles`).
    sync_rpm: f64,
    /// Rotor mechanical speed in rev/min (`Ns (1 - s)`).
    rotor_rpm: f64,
    /// Fractional slip in `[0, 1]` (`(Ns - Nr) / Ns`).
    slip: f64,
}

impl InductionMotor {
    /// Build a model from a supply frequency, pole count, and a
    /// measured rotor speed.
    ///
    /// The synchronous speed is computed as `120 f / poles`, the slip
    /// as `(Ns - Nr) / Ns`, and the result is **rejected** if that slip
    /// falls outside the motoring range `[0, 1]`.
    ///
    /// # Errors
    ///
    /// - [`InductionMotorError::InvalidFrequency`] if `freq_hz` is not
    ///   finite and `> 0`.
    /// - [`InductionMotorError::InvalidPoles`] if `poles` is zero or odd.
    /// - [`InductionMotorError::InvalidRotorSpeed`] if `rotor_rpm` is
    ///   not finite.
    /// - [`InductionMotorError::SlipOutOfRange`] if the implied slip is
    ///   `< 0` (rotor faster than the field) or `> 1` (rotor driven
    ///   backwards).
    ///
    /// # Examples
    ///
    /// ```
    /// use valenx_inductionmotor::InductionMotor;
    /// let m = InductionMotor::new(60.0, 4, 1746.0).unwrap();
    /// assert!((m.sync_speed_rpm() - 1800.0).abs() < 1e-9);
    /// assert!((m.slip() - 0.03).abs() < 1e-12);
    /// ```
    pub fn new(freq_hz: f64, poles: u32, rotor_rpm: f64) -> Result<Self, InductionMotorError> {
        let sync_rpm = sync_speed_rpm(freq_hz, poles)?;
        if !rotor_rpm.is_finite() {
            return Err(InductionMotorError::InvalidRotorSpeed { rpm: rotor_rpm });
        }
        let s = (sync_rpm - rotor_rpm) / sync_rpm;
        if !(0.0..=1.0).contains(&s) {
            return Err(InductionMotorError::SlipOutOfRange {
                slip: s,
                sync_rpm,
                rotor_rpm,
            });
        }
        Ok(Self {
            poles,
            freq_hz,
            sync_rpm,
            rotor_rpm,
            slip: s,
        })
    }

    /// Build a model from a supply frequency, pole count, and a slip
    /// fraction, deriving the rotor speed as `Nr = Ns (1 - s)`.
    ///
    /// # Errors
    ///
    /// - [`InductionMotorError::InvalidFrequency`] if `freq_hz` is not
    ///   finite and `> 0`.
    /// - [`InductionMotorError::InvalidPoles`] if `poles` is zero or odd.
    /// - [`InductionMotorError::InvalidSlip`] if `slip` is not finite or
    ///   lies outside `[0, 1]`.
    ///
    /// # Examples
    ///
    /// ```
    /// use valenx_inductionmotor::InductionMotor;
    /// // Standstill: slip = 1 -> rotor at rest.
    /// let locked = InductionMotor::from_slip(50.0, 2, 1.0).unwrap();
    /// assert!(locked.rotor_speed_rpm().abs() < 1e-9);
    /// ```
    pub fn from_slip(freq_hz: f64, poles: u32, slip: f64) -> Result<Self, InductionMotorError> {
        let sync_rpm = sync_speed_rpm(freq_hz, poles)?;
        if !slip.is_finite() || !(0.0..=1.0).contains(&slip) {
            return Err(InductionMotorError::InvalidSlip { slip });
        }
        let rotor_rpm = sync_rpm * (1.0 - slip);
        Ok(Self {
            poles,
            freq_hz,
            sync_rpm,
            rotor_rpm,
            slip,
        })
    }

    /// Number of magnetic poles (a positive even integer).
    pub fn poles(&self) -> u32 {
        self.poles
    }

    /// Line / supply frequency in hertz.
    pub fn supply_frequency_hz(&self) -> f64 {
        self.freq_hz
    }

    /// Synchronous speed of the rotating field, in rev/min
    /// (`Ns = 120 f / poles`).
    pub fn sync_speed_rpm(&self) -> f64 {
        self.sync_rpm
    }

    /// Rotor mechanical speed, in rev/min (`Nr = Ns (1 - s)`).
    pub fn rotor_speed_rpm(&self) -> f64 {
        self.rotor_rpm
    }

    /// Fractional slip, dimensionless and within `[0, 1]`
    /// (`s = (Ns - Nr) / Ns`).
    pub fn slip(&self) -> f64 {
        self.slip
    }

    /// Slip expressed as a percentage (`100 s`).
    pub fn slip_percent(&self) -> f64 {
        self.slip * 100.0
    }

    /// Rotor electrical frequency (slip frequency), in hertz
    /// (`f_r = s f`).
    pub fn rotor_frequency_hz(&self) -> f64 {
        self.slip * self.freq_hz
    }

    /// Slip speed: the difference between synchronous and rotor speed,
    /// in rev/min (`Ns - Nr = s Ns`).
    pub fn slip_speed_rpm(&self) -> f64 {
        self.sync_rpm - self.rotor_rpm
    }

    /// Rotor copper (`I^2 R`) loss for a given air-gap power,
    /// `P_rotor_cu = s * P_airgap` (see [`rotor_copper_loss_w`]).
    ///
    /// # Errors
    ///
    /// Returns [`InductionMotorError::InvalidPower`] if `air_gap_power_w`
    /// is not finite or negative. (The slip was validated at construction.)
    pub fn rotor_copper_loss_w(&self, air_gap_power_w: f64) -> Result<f64, InductionMotorError> {
        let pag = validate_power(air_gap_power_w)?;
        Ok(self.slip * pag)
    }

    /// Developed (gross) mechanical power for a given air-gap power,
    /// `P_mech = (1 - s) * P_airgap` (see [`developed_mechanical_power_w`]).
    ///
    /// # Errors
    ///
    /// Returns [`InductionMotorError::InvalidPower`] if `air_gap_power_w`
    /// is not finite or negative.
    pub fn developed_mechanical_power_w(
        &self,
        air_gap_power_w: f64,
    ) -> Result<f64, InductionMotorError> {
        let pag = validate_power(air_gap_power_w)?;
        Ok((1.0 - self.slip) * pag)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::error::ErrorCategory;

    const EPS: f64 = 1e-9;

    // -- sync_speed_rpm: Ns = 120 f / poles --------------------------

    #[test]
    fn sync_speed_60hz_4pole_is_1800() {
        let ns = sync_speed_rpm(60.0, 4).unwrap();
        assert!((ns - 1800.0).abs() < EPS);
    }

    #[test]
    fn sync_speed_50hz_4pole_is_1500() {
        let ns = sync_speed_rpm(50.0, 4).unwrap();
        assert!((ns - 1500.0).abs() < EPS);
    }

    #[test]
    fn sync_speed_60hz_2pole_is_3600() {
        let ns = sync_speed_rpm(60.0, 2).unwrap();
        assert!((ns - 3600.0).abs() < EPS);
    }

    #[test]
    fn sync_speed_50hz_6pole_is_1000() {
        let ns = sync_speed_rpm(50.0, 6).unwrap();
        assert!((ns - 1000.0).abs() < EPS);
    }

    #[test]
    fn sync_speed_matches_closed_form_over_grid() {
        for &poles in &[2u32, 4, 6, 8, 12] {
            for &f in &[50.0_f64, 60.0, 25.0, 400.0] {
                let ns = sync_speed_rpm(f, poles).unwrap();
                let expected = 120.0 * f / poles as f64;
                assert!((ns - expected).abs() < EPS, "poles={poles} f={f}");
            }
        }
    }

    #[test]
    fn sync_speed_rejects_zero_poles() {
        let err = sync_speed_rpm(60.0, 0).unwrap_err();
        assert_eq!(err, InductionMotorError::InvalidPoles { poles: 0 });
        assert_eq!(err.category(), ErrorCategory::Input);
    }

    #[test]
    fn sync_speed_rejects_odd_poles() {
        assert!(matches!(
            sync_speed_rpm(60.0, 3),
            Err(InductionMotorError::InvalidPoles { poles: 3 })
        ));
    }

    #[test]
    fn sync_speed_rejects_nonpositive_and_nonfinite_frequency() {
        assert!(matches!(
            sync_speed_rpm(0.0, 4),
            Err(InductionMotorError::InvalidFrequency { .. })
        ));
        assert!(matches!(
            sync_speed_rpm(-60.0, 4),
            Err(InductionMotorError::InvalidFrequency { .. })
        ));
        assert!(matches!(
            sync_speed_rpm(f64::NAN, 4),
            Err(InductionMotorError::InvalidFrequency { .. })
        ));
        assert!(matches!(
            sync_speed_rpm(f64::INFINITY, 4),
            Err(InductionMotorError::InvalidFrequency { .. })
        ));
    }

    // -- slip: s = (Ns - Nr) / Ns ------------------------------------

    #[test]
    fn slip_is_zero_at_synchronous_speed() {
        let s = slip(1800.0, 1800.0).unwrap();
        assert!(s.abs() < EPS);
    }

    #[test]
    fn slip_is_one_at_standstill() {
        let s = slip(1800.0, 0.0).unwrap();
        assert!((s - 1.0).abs() < EPS);
    }

    #[test]
    fn slip_three_percent_operating_point() {
        let s = slip(1800.0, 1746.0).unwrap();
        assert!((s - 0.03).abs() < 1e-12);
    }

    #[test]
    fn slip_is_half_at_half_synchronous_speed() {
        let s = slip(1500.0, 750.0).unwrap();
        assert!((s - 0.5).abs() < EPS);
    }

    #[test]
    fn slip_rejects_nonpositive_sync_speed() {
        assert!(matches!(
            slip(0.0, 0.0),
            Err(InductionMotorError::InvalidFrequency { .. })
        ));
    }

    #[test]
    fn slip_rejects_nonfinite_rotor_speed() {
        assert!(matches!(
            slip(1800.0, f64::NAN),
            Err(InductionMotorError::InvalidRotorSpeed { .. })
        ));
    }

    // -- rotor_frequency_hz: f_r = s f -------------------------------

    #[test]
    fn rotor_frequency_zero_at_sync() {
        let fr = rotor_frequency_hz(0.0, 60.0).unwrap();
        assert!(fr.abs() < EPS);
    }

    #[test]
    fn rotor_frequency_equals_supply_at_standstill() {
        let fr = rotor_frequency_hz(1.0, 60.0).unwrap();
        assert!((fr - 60.0).abs() < EPS);
    }

    #[test]
    fn rotor_frequency_three_percent() {
        let fr = rotor_frequency_hz(0.03, 60.0).unwrap();
        assert!((fr - 1.8).abs() < 1e-12);
    }

    #[test]
    fn rotor_frequency_is_slip_times_supply_over_grid() {
        for &s in &[0.0_f64, 0.01, 0.05, 0.2, 0.5, 1.0] {
            for &f in &[50.0_f64, 60.0, 400.0] {
                let fr = rotor_frequency_hz(s, f).unwrap();
                assert!((fr - s * f).abs() < EPS, "s={s} f={f}");
            }
        }
    }

    #[test]
    fn rotor_frequency_rejects_nonfinite_slip() {
        assert!(matches!(
            rotor_frequency_hz(f64::NAN, 60.0),
            Err(InductionMotorError::InvalidSlip { .. })
        ));
    }

    // -- rotor_speed_rpm: Nr = Ns (1 - s) ----------------------------

    #[test]
    fn rotor_speed_equals_sync_at_zero_slip() {
        let nr = rotor_speed_rpm(1800.0, 0.0).unwrap();
        assert!((nr - 1800.0).abs() < EPS);
    }

    #[test]
    fn rotor_speed_is_zero_at_unit_slip() {
        let nr = rotor_speed_rpm(1800.0, 1.0).unwrap();
        assert!(nr.abs() < EPS);
    }

    #[test]
    fn rotor_speed_three_percent_slip() {
        let nr = rotor_speed_rpm(1800.0, 0.03).unwrap();
        assert!((nr - 1746.0).abs() < EPS);
    }

    // -- round-trip identities ---------------------------------------

    #[test]
    fn slip_and_rotor_speed_are_inverse() {
        // Ns (1 - (Ns - Nr)/Ns) == Nr for arbitrary valid inputs.
        let sync_rpm = 1500.0;
        for &rotor in &[0.0_f64, 375.0, 750.0, 1425.0, 1500.0] {
            let s = slip(sync_rpm, rotor).unwrap();
            let back = rotor_speed_rpm(sync_rpm, s).unwrap();
            assert!((back - rotor).abs() < EPS, "rotor={rotor}");
        }
    }

    // -- InductionMotor aggregate ------------------------------------

    #[test]
    fn motor_new_matches_textbook_quantities() {
        // 4-pole, 60 Hz, 1746 rev/min: Ns=1800, s=0.03, fr=1.8 Hz.
        let m = InductionMotor::new(60.0, 4, 1746.0).unwrap();
        assert_eq!(m.poles(), 4);
        assert!((m.supply_frequency_hz() - 60.0).abs() < EPS);
        assert!((m.sync_speed_rpm() - 1800.0).abs() < EPS);
        assert!((m.rotor_speed_rpm() - 1746.0).abs() < EPS);
        assert!((m.slip() - 0.03).abs() < 1e-12);
        assert!((m.slip_percent() - 3.0).abs() < 1e-9);
        assert!((m.rotor_frequency_hz() - 1.8).abs() < 1e-12);
        assert!((m.slip_speed_rpm() - 54.0).abs() < EPS);
    }

    #[test]
    fn motor_from_slip_derives_rotor_speed() {
        let m = InductionMotor::from_slip(60.0, 4, 0.03).unwrap();
        assert!((m.rotor_speed_rpm() - 1746.0).abs() < EPS);
        assert!((m.rotor_frequency_hz() - 1.8).abs() < 1e-12);
    }

    #[test]
    fn motor_at_synchronous_speed_has_zero_slip() {
        let m = InductionMotor::new(50.0, 4, 1500.0).unwrap();
        assert!(m.slip().abs() < EPS);
        assert!(m.rotor_frequency_hz().abs() < EPS);
        assert!(m.slip_speed_rpm().abs() < EPS);
    }

    #[test]
    fn motor_at_standstill_has_unit_slip() {
        let m = InductionMotor::new(50.0, 2, 0.0).unwrap();
        assert!((m.slip() - 1.0).abs() < EPS);
        assert!((m.rotor_frequency_hz() - 50.0).abs() < EPS);
    }

    #[test]
    fn motor_new_rejects_hypersynchronous_negative_slip() {
        // Rotor faster than the 1800 rev/min field -> slip < 0.
        let err = InductionMotor::new(60.0, 4, 1900.0).unwrap_err();
        assert!(matches!(err, InductionMotorError::SlipOutOfRange { .. }));
        assert_eq!(err.category(), ErrorCategory::Domain);
    }

    #[test]
    fn motor_new_rejects_reverse_rotation_slip_above_one() {
        // Rotor driven backwards -> slip > 1.
        assert!(matches!(
            InductionMotor::new(60.0, 4, -100.0),
            Err(InductionMotorError::SlipOutOfRange { .. })
        ));
    }

    #[test]
    fn motor_from_slip_rejects_out_of_range_slip() {
        assert!(matches!(
            InductionMotor::from_slip(60.0, 4, 1.5),
            Err(InductionMotorError::InvalidSlip { .. })
        ));
        assert!(matches!(
            InductionMotor::from_slip(60.0, 4, -0.1),
            Err(InductionMotorError::InvalidSlip { .. })
        ));
    }

    #[test]
    fn motor_constructors_agree() {
        // new(.., Nr) and from_slip(.., s) must yield the same point.
        let a = InductionMotor::new(60.0, 6, 1140.0).unwrap();
        let b = InductionMotor::from_slip(60.0, 6, a.slip()).unwrap();
        assert!((a.sync_speed_rpm() - b.sync_speed_rpm()).abs() < EPS);
        assert!((a.rotor_speed_rpm() - b.rotor_speed_rpm()).abs() < EPS);
        assert!((a.slip() - b.slip()).abs() < EPS);
    }

    #[test]
    fn motor_serde_roundtrips() {
        let m = InductionMotor::new(60.0, 4, 1746.0).unwrap();
        let json = serde_json::to_string(&m).unwrap();
        let back: InductionMotor = serde_json::from_str(&json).unwrap();
        assert_eq!(m, back);
    }

    // -- air-gap power split: Pag = P_cu + P_mech, P_cu = s Pag -------

    #[test]
    fn power_split_hand_values() {
        // 10 kW air-gap power at 3 % slip: 300 W rotor loss, 9.7 kW mech.
        let p_cu = rotor_copper_loss_w(0.03, 10_000.0).unwrap();
        let p_mech = developed_mechanical_power_w(0.03, 10_000.0).unwrap();
        assert!((p_cu - 300.0).abs() < EPS, "p_cu = {p_cu}");
        assert!((p_mech - 9_700.0).abs() < EPS, "p_mech = {p_mech}");
    }

    #[test]
    fn power_split_sums_to_air_gap_power() {
        // Energy conservation: P_cu + P_mech == Pag for any s in [0, 1].
        let pag = 7_500.0;
        for &s in &[0.0_f64, 0.02, 0.1, 0.5, 1.0] {
            let p_cu = rotor_copper_loss_w(s, pag).unwrap();
            let p_mech = developed_mechanical_power_w(s, pag).unwrap();
            assert!((p_cu + p_mech - pag).abs() < EPS, "sum at s = {s}");
        }
    }

    #[test]
    fn power_split_limits() {
        let pag = 5_000.0;
        // Synchronous speed (s = 0): no rotor loss, all power mechanical.
        assert!(rotor_copper_loss_w(0.0, pag).unwrap().abs() < EPS);
        assert!((developed_mechanical_power_w(0.0, pag).unwrap() - pag).abs() < EPS);
        // Standstill (s = 1): all power lost in the rotor, none mechanical.
        assert!((rotor_copper_loss_w(1.0, pag).unwrap() - pag).abs() < EPS);
        assert!(developed_mechanical_power_w(1.0, pag).unwrap().abs() < EPS);
    }

    #[test]
    fn power_split_fractions_track_slip() {
        // P_cu/Pag = s and P_mech/Pag = 1 - s; more slip -> more rotor loss.
        let pag = 4_000.0;
        let lo = rotor_copper_loss_w(0.02, pag).unwrap();
        let hi = rotor_copper_loss_w(0.08, pag).unwrap();
        assert!(hi > lo, "more slip should waste more: {hi} vs {lo}");
        for &s in &[0.01_f64, 0.05, 0.25] {
            assert!((rotor_copper_loss_w(s, pag).unwrap() / pag - s).abs() < EPS);
            assert!((developed_mechanical_power_w(s, pag).unwrap() / pag - (1.0 - s)).abs() < EPS);
        }
    }

    #[test]
    fn power_split_zero_air_gap_power_is_zero() {
        assert!(rotor_copper_loss_w(0.3, 0.0).unwrap().abs() < EPS);
        assert!(developed_mechanical_power_w(0.3, 0.0).unwrap().abs() < EPS);
    }

    #[test]
    fn motor_power_methods_match_free_functions() {
        let m = InductionMotor::new(60.0, 4, 1746.0).unwrap(); // s = 0.03
        let pag = 12_000.0;
        assert!(
            (m.rotor_copper_loss_w(pag).unwrap() - rotor_copper_loss_w(m.slip(), pag).unwrap())
                .abs()
                < EPS
        );
        assert!(
            (m.developed_mechanical_power_w(pag).unwrap()
                - developed_mechanical_power_w(m.slip(), pag).unwrap())
            .abs()
                < EPS
        );
        // On the motor, the two pieces still sum to the air-gap power.
        assert!(
            (m.rotor_copper_loss_w(pag).unwrap() + m.developed_mechanical_power_w(pag).unwrap()
                - pag)
                .abs()
                < EPS
        );
    }

    #[test]
    fn power_split_rejects_bad_inputs() {
        // Bad slip.
        assert!(matches!(
            rotor_copper_loss_w(-0.1, 1000.0),
            Err(InductionMotorError::InvalidSlip { .. })
        ));
        assert!(matches!(
            developed_mechanical_power_w(1.5, 1000.0),
            Err(InductionMotorError::InvalidSlip { .. })
        ));
        assert!(matches!(
            rotor_copper_loss_w(f64::NAN, 1000.0),
            Err(InductionMotorError::InvalidSlip { .. })
        ));
        // Bad power.
        assert!(matches!(
            rotor_copper_loss_w(0.03, -1.0),
            Err(InductionMotorError::InvalidPower { .. })
        ));
        assert!(matches!(
            developed_mechanical_power_w(0.03, f64::INFINITY),
            Err(InductionMotorError::InvalidPower { .. })
        ));
        let m = InductionMotor::new(60.0, 4, 1746.0).unwrap();
        assert!(matches!(
            m.rotor_copper_loss_w(f64::NAN),
            Err(InductionMotorError::InvalidPower { .. })
        ));
    }
}
