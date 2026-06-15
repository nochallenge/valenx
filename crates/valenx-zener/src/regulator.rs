//! Ideal shunt-regulator design and operating-point analysis.
//!
//! All quantities are SI base units: volts (V), amperes (A), ohms (Ω),
//! watts (W). The diode is modelled as an ideal constant-voltage shunt
//! once it is in reverse breakdown: terminal voltage is exactly the
//! nominal zener voltage `Vz` and is independent of the current through
//! it (no dynamic impedance `r_z`). See the crate-level honest-scope
//! note for the consequences of that idealisation.

use serde::{Deserialize, Serialize};

use crate::error::{require_non_negative, require_positive, ZenerError};

/// A zener-diode shunt regulator with a chosen series resistor.
///
/// Topology: the unregulated supply `Vin` drives a series resistor
/// `Rs`; the zener diode and the load hang in parallel from the
/// resistor's output node to ground. The diode pins that node at its
/// reverse-breakdown voltage `Vz`, so the load sees a regulated `Vz`
/// while the resistor absorbs the difference `Vin - Vz`.
///
/// ```text
///   Vin ──[ Rs ]──┬────────── Vout = Vz
///                 │        │
///               [ Dz ]   [ load ]
///                 │        │
///                gnd      gnd
/// ```
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct ZenerRegulator {
    /// Nominal zener (reverse-breakdown) voltage `Vz`, in volts.
    pub zener_voltage_v: f64,
    /// Series resistor `Rs` between the supply and the output node,
    /// in ohms.
    pub series_resistor_ohm: f64,
}

impl ZenerRegulator {
    /// Build a regulator from a known zener voltage and series resistor.
    ///
    /// # Errors
    /// Returns [`ZenerError::BadParameter`] if either argument is not a
    /// finite, strictly-positive number.
    pub fn new(zener_voltage_v: f64, series_resistor_ohm: f64) -> Result<Self, ZenerError> {
        let zener_voltage_v = require_positive("zener_voltage_v", zener_voltage_v)?;
        let series_resistor_ohm = require_positive("series_resistor_ohm", series_resistor_ohm)?;
        Ok(Self {
            zener_voltage_v,
            series_resistor_ohm,
        })
    }

    /// Total current drawn from the supply through the series resistor
    /// when the output node sits at `Vz`:
    ///
    /// `I_Rs = (Vin - Vz) / Rs`.
    ///
    /// This is the current the resistor delivers; it splits between the
    /// diode and the load. Returns `0.0` (clamped, never negative) when
    /// `Vin <= Vz`, i.e. when there is no headroom and the diode is out
    /// of regulation.
    ///
    /// # Errors
    /// Returns [`ZenerError::BadParameter`] if `vin_v` is not finite.
    pub fn resistor_current_a(&self, vin_v: f64) -> Result<f64, ZenerError> {
        if !vin_v.is_finite() {
            return Err(ZenerError::BadParameter {
                name: "vin_v",
                reason: format!("must be finite, got {vin_v}"),
            });
        }
        let headroom = vin_v - self.zener_voltage_v;
        if headroom <= 0.0 {
            return Ok(0.0);
        }
        Ok(headroom / self.series_resistor_ohm)
    }

    /// Regulated output voltage seen by the load.
    ///
    /// While the diode is in breakdown the output is clamped to `Vz`,
    /// independent of `Vin` and the load. If the supply cannot push the
    /// node up to `Vz` (`Vin <= Vz`, so `I_Rs == 0`), the diode is dark
    /// and the output simply follows the supply through the resistor;
    /// with the load drawing `i_load_a` that unregulated level is
    /// `Vin - i_load_a * Rs` (floored at zero).
    ///
    /// # Errors
    /// Returns [`ZenerError::BadParameter`] if `vin_v` is not finite or
    /// `i_load_a` is negative / non-finite.
    pub fn output_voltage_v(&self, vin_v: f64, i_load_a: f64) -> Result<f64, ZenerError> {
        let i_load_a = require_non_negative("i_load_a", i_load_a)?;
        let i_rs = self.resistor_current_a(vin_v)?;
        if i_rs > 0.0 {
            // Diode is in breakdown: hard clamp at Vz.
            Ok(self.zener_voltage_v)
        } else {
            // No headroom: unregulated, output droops below the supply
            // by the resistor drop from the load current.
            Ok((vin_v - i_load_a * self.series_resistor_ohm).max(0.0))
        }
    }

    /// Current actually flowing through the zener diode at this
    /// operating point:
    ///
    /// `Iz = I_Rs - I_load`.
    ///
    /// The resistor current splits between diode and load by KCL; what
    /// the load does not take, the diode shunts to ground. A result of
    /// `0.0` means the load has consumed the entire resistor current and
    /// the diode is on the edge of dropping out of regulation; a value
    /// at or below zero is reported as the [`OperatingPoint::in_regulation`]
    /// flag being `false`.
    ///
    /// # Errors
    /// Returns [`ZenerError::BadParameter`] if `vin_v` is not finite or
    /// `i_load_a` is negative / non-finite.
    pub fn zener_current_a(&self, vin_v: f64, i_load_a: f64) -> Result<f64, ZenerError> {
        let i_load_a = require_non_negative("i_load_a", i_load_a)?;
        let i_rs = self.resistor_current_a(vin_v)?;
        Ok(i_rs - i_load_a)
    }

    /// Power dissipated in the zener diode:
    ///
    /// `Pz = Vz * Iz`.
    ///
    /// Uses the clamped (non-negative) diode current, so this never
    /// reports a spurious negative dissipation when the load has
    /// stolen all of the resistor current.
    ///
    /// # Errors
    /// Returns [`ZenerError::BadParameter`] if `vin_v` is not finite or
    /// `i_load_a` is negative / non-finite.
    pub fn zener_power_w(&self, vin_v: f64, i_load_a: f64) -> Result<f64, ZenerError> {
        let iz = self.zener_current_a(vin_v, i_load_a)?.max(0.0);
        Ok(self.zener_voltage_v * iz)
    }

    /// Power dissipated in the series resistor: `P_Rs = (Vin - Vz) * I_Rs`
    /// while in regulation (equivalently `I_Rs² * Rs`).
    ///
    /// # Errors
    /// Returns [`ZenerError::BadParameter`] if `vin_v` is not finite.
    pub fn resistor_power_w(&self, vin_v: f64) -> Result<f64, ZenerError> {
        let i_rs = self.resistor_current_a(vin_v)?;
        Ok(i_rs * i_rs * self.series_resistor_ohm)
    }

    /// Minimum supply voltage that still keeps the diode in regulation
    /// at the given worst-case operating point.
    ///
    /// To stay in breakdown the diode must retain at least its minimum
    /// holding current `iz_min_a` while the load simultaneously draws
    /// its maximum `i_load_max_a`. The resistor must therefore carry
    /// `iz_min_a + i_load_max_a`, which across `Rs` demands a supply of
    ///
    /// `Vin_min = Vz + Rs * (Iz_min + I_load_max)`.
    ///
    /// # Errors
    /// Returns [`ZenerError::BadParameter`] if either current argument
    /// is negative / non-finite.
    pub fn min_input_voltage_v(&self, iz_min_a: f64, i_load_max_a: f64) -> Result<f64, ZenerError> {
        let iz_min_a = require_non_negative("iz_min_a", iz_min_a)?;
        let i_load_max_a = require_non_negative("i_load_max_a", i_load_max_a)?;
        Ok(self.zener_voltage_v + self.series_resistor_ohm * (iz_min_a + i_load_max_a))
    }

    /// Full operating-point summary at one `(Vin, I_load)` working point.
    ///
    /// Bundles the regulated output, the resistor and diode currents,
    /// and the diode and resistor dissipations, plus an `in_regulation`
    /// flag (`true` iff the diode still carries strictly-positive
    /// current). Convenient for tables and sweeps.
    ///
    /// # Errors
    /// Returns [`ZenerError::BadParameter`] if `vin_v` is not finite or
    /// `i_load_a` is negative / non-finite.
    pub fn operating_point(&self, vin_v: f64, i_load_a: f64) -> Result<OperatingPoint, ZenerError> {
        let i_load_a = require_non_negative("i_load_a", i_load_a)?;
        let i_rs = self.resistor_current_a(vin_v)?;
        let iz = i_rs - i_load_a;
        let in_regulation = i_rs > 0.0 && iz > 0.0;
        let v_out = if i_rs > 0.0 {
            self.zener_voltage_v
        } else {
            (vin_v - i_load_a * self.series_resistor_ohm).max(0.0)
        };
        let iz_clamped = iz.max(0.0);
        Ok(OperatingPoint {
            output_voltage_v: v_out,
            resistor_current_a: i_rs,
            zener_current_a: iz,
            zener_power_w: self.zener_voltage_v * iz_clamped,
            resistor_power_w: i_rs * i_rs * self.series_resistor_ohm,
            in_regulation,
        })
    }
}

/// Resolved electrical state of a regulator at a single operating point.
///
/// Returned by [`ZenerRegulator::operating_point`]. All fields are SI
/// units (V, A, W). `zener_current_a` is reported signed so callers can
/// detect drop-out (a non-positive value), while the power figures use
/// the clamped current and are therefore never negative.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct OperatingPoint {
    /// Output voltage delivered to the load, in volts. Equals `Vz`
    /// while in regulation.
    pub output_voltage_v: f64,
    /// Current through the series resistor `I_Rs`, in amperes.
    pub resistor_current_a: f64,
    /// Signed current through the zener diode `Iz = I_Rs - I_load`, in
    /// amperes. A non-positive value means the diode has dropped out of
    /// regulation.
    pub zener_current_a: f64,
    /// Power dissipated in the diode `Pz = Vz * Iz` (clamped current),
    /// in watts.
    pub zener_power_w: f64,
    /// Power dissipated in the series resistor `I_Rs² * Rs`, in watts.
    pub resistor_power_w: f64,
    /// `true` iff the diode is carrying strictly-positive current and
    /// the supply has headroom, i.e. the output is genuinely clamped.
    pub in_regulation: bool,
}

/// Size the series resistor `Rs` for a target operating point.
///
/// Given the unregulated supply `vin_v`, the zener voltage `zener_voltage_v`,
/// the desired (nominal) zener current `iz_nominal_a`, and the load
/// current `i_load_a`, the resistor must carry the sum `Iz + Il` across
/// the headroom `Vin - Vz`:
///
/// `Rs = (Vin - Vz) / (Iz + Il)`.
///
/// # Errors
/// Returns [`ZenerError::BadParameter`] for non-finite / non-positive
/// `vin_v` or `zener_voltage_v`, for negative / non-finite currents, or
/// when both currents are zero (the divisor would vanish). Returns
/// [`ZenerError::CannotRegulate`] when `vin_v <= zener_voltage_v`, since
/// no positive resistor can drop a non-positive headroom.
pub fn size_series_resistor(
    vin_v: f64,
    zener_voltage_v: f64,
    iz_nominal_a: f64,
    i_load_a: f64,
) -> Result<f64, ZenerError> {
    let vin_v = require_positive("vin_v", vin_v)?;
    let zener_voltage_v = require_positive("zener_voltage_v", zener_voltage_v)?;
    let iz_nominal_a = require_non_negative("iz_nominal_a", iz_nominal_a)?;
    let i_load_a = require_non_negative("i_load_a", i_load_a)?;

    let total_current = iz_nominal_a + i_load_a;
    if total_current <= 0.0 {
        return Err(ZenerError::BadParameter {
            name: "iz_nominal_a + i_load_a",
            reason: "total branch current must be > 0".to_string(),
        });
    }
    if vin_v <= zener_voltage_v {
        return Err(ZenerError::CannotRegulate(format!(
            "Vin ({vin_v} V) must exceed Vz ({zener_voltage_v} V) for headroom across Rs"
        )));
    }
    Ok((vin_v - zener_voltage_v) / total_current)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::error::ErrorCategory;

    const EPS: f64 = 1e-12;

    // --- size_series_resistor: Rs = (Vin - Vz) / (Iz + Il) ----------

    #[test]
    fn sizes_resistor_from_canonical_formula() {
        // Vin=12, Vz=5.1, Iz=5 mA, Il=20 mA -> Rs = 6.9 / 0.025 = 276 Ω
        let rs = size_series_resistor(12.0, 5.1, 0.005, 0.020).unwrap();
        assert!((rs - 276.0).abs() < 1e-9, "got {rs}");
    }

    #[test]
    fn sizes_resistor_no_load() {
        // Pure reference: Il=0, Iz=10 mA. Rs = (15-5)/0.010 = 1000 Ω
        let rs = size_series_resistor(15.0, 5.0, 0.010, 0.0).unwrap();
        assert!((rs - 1000.0).abs() < 1e-9, "got {rs}");
    }

    #[test]
    fn sizing_rejects_zero_total_current() {
        let err = size_series_resistor(12.0, 5.0, 0.0, 0.0).unwrap_err();
        assert_eq!(err.code(), "zener.bad_parameter");
        assert_eq!(err.category(), ErrorCategory::Input);
    }

    #[test]
    fn sizing_rejects_vin_not_above_vz() {
        let err = size_series_resistor(5.0, 5.0, 0.005, 0.010).unwrap_err();
        assert_eq!(err.code(), "zener.cannot_regulate");
        assert_eq!(err.category(), ErrorCategory::Operating);
    }

    #[test]
    fn sizing_rejects_nonfinite() {
        assert!(size_series_resistor(f64::NAN, 5.0, 0.005, 0.010).is_err());
        assert!(size_series_resistor(12.0, f64::INFINITY, 0.005, 0.010).is_err());
        assert!(size_series_resistor(12.0, 5.0, -0.001, 0.010).is_err());
    }

    // --- round trip: size Rs, then recompute the operating point ----

    #[test]
    fn sized_resistor_reproduces_target_currents() {
        // Design point: Vin=12, Vz=5.1, Iz=5 mA, Il=20 mA.
        let vin = 12.0;
        let vz = 5.1;
        let iz_target = 0.005;
        let il = 0.020;
        let rs = size_series_resistor(vin, vz, iz_target, il).unwrap();
        let reg = ZenerRegulator::new(vz, rs).unwrap();

        // I_Rs should equal Iz + Il = 25 mA.
        let i_rs = reg.resistor_current_a(vin).unwrap();
        assert!((i_rs - 0.025).abs() < EPS, "I_Rs = {i_rs}");

        // Iz = I_Rs - Il = 5 mA (the design current).
        let iz = reg.zener_current_a(vin, il).unwrap();
        assert!((iz - iz_target).abs() < EPS, "Iz = {iz}");
    }

    // --- Vout clamped at Vz -----------------------------------------

    #[test]
    fn output_is_clamped_at_vz_in_regulation() {
        let reg = ZenerRegulator::new(5.1, 276.0).unwrap();
        // Across a wide supply range the output stays pinned at Vz.
        for &vin in &[8.0, 10.0, 12.0, 15.0, 20.0] {
            let v_out = reg.output_voltage_v(vin, 0.020).unwrap();
            assert!((v_out - 5.1).abs() < EPS, "Vin={vin} -> Vout={v_out}");
        }
    }

    #[test]
    fn output_droops_unregulated_when_no_headroom() {
        let reg = ZenerRegulator::new(5.1, 276.0).unwrap();
        // Vin below Vz: diode dark, output follows supply minus IR drop.
        // Vin=4, Il=10 mA -> 4 - 0.010*276 = 1.24 V.
        let v_out = reg.output_voltage_v(4.0, 0.010).unwrap();
        assert!((v_out - 1.24).abs() < 1e-9, "Vout = {v_out}");
        // And it never goes negative even for a heavy load.
        let v_out2 = reg.output_voltage_v(4.0, 1.0).unwrap();
        assert!(v_out2 >= 0.0, "Vout = {v_out2}");
    }

    // --- Iz = I_Rs - I_load -----------------------------------------

    #[test]
    fn zener_current_is_resistor_minus_load() {
        let reg = ZenerRegulator::new(5.0, 500.0).unwrap();
        // Vin=15 -> I_Rs = 10/500 = 20 mA.
        let i_rs = reg.resistor_current_a(15.0).unwrap();
        assert!((i_rs - 0.020).abs() < EPS, "I_Rs = {i_rs}");
        // Il = 8 mA -> Iz = 12 mA.
        let iz = reg.zener_current_a(15.0, 0.008).unwrap();
        assert!((iz - 0.012).abs() < EPS, "Iz = {iz}");
        // Identity Iz == I_Rs - Il holds for an arbitrary load.
        let il = 0.013;
        let iz2 = reg.zener_current_a(15.0, il).unwrap();
        assert!((iz2 - (i_rs - il)).abs() < EPS, "Iz = {iz2}");
    }

    #[test]
    fn zener_current_goes_negative_on_overload() {
        let reg = ZenerRegulator::new(5.0, 500.0).unwrap();
        // I_Rs = 20 mA but load wants 30 mA -> signed Iz = -10 mA.
        let iz = reg.zener_current_a(15.0, 0.030).unwrap();
        assert!((iz - (-0.010)).abs() < EPS, "Iz = {iz}");
        // And the operating point flags drop-out.
        let op = reg.operating_point(15.0, 0.030).unwrap();
        assert!(!op.in_regulation);
    }

    // --- min input voltage for regulation ---------------------------

    #[test]
    fn min_input_voltage_matches_formula() {
        let reg = ZenerRegulator::new(5.1, 276.0).unwrap();
        // Vin_min = Vz + Rs*(Iz_min + Il_max)
        //         = 5.1 + 276*(0.001 + 0.020) = 5.1 + 5.796 = 10.896 V
        let vin_min = reg.min_input_voltage_v(0.001, 0.020).unwrap();
        assert!((vin_min - 10.896).abs() < 1e-9, "Vin_min = {vin_min}");
    }

    #[test]
    fn at_min_input_voltage_diode_holds_iz_min() {
        // Self-consistency: drive the regulator at exactly Vin_min with
        // the worst-case load and the diode current should land on the
        // requested minimum holding current.
        let reg = ZenerRegulator::new(5.1, 276.0).unwrap();
        let iz_min = 0.001;
        let il_max = 0.020;
        let vin_min = reg.min_input_voltage_v(iz_min, il_max).unwrap();
        let iz = reg.zener_current_a(vin_min, il_max).unwrap();
        assert!((iz - iz_min).abs() < 1e-9, "Iz = {iz}");
    }

    #[test]
    fn min_input_voltage_rejects_negative_currents() {
        let reg = ZenerRegulator::new(5.0, 500.0).unwrap();
        assert!(reg.min_input_voltage_v(-0.001, 0.020).is_err());
        assert!(reg.min_input_voltage_v(0.001, -0.020).is_err());
    }

    // --- power Pz = Vz * Iz -----------------------------------------

    #[test]
    fn zener_power_matches_formula() {
        let reg = ZenerRegulator::new(5.0, 500.0).unwrap();
        // Iz = 12 mA, Pz = 5.0 * 0.012 = 0.060 W.
        let pz = reg.zener_power_w(15.0, 0.008).unwrap();
        assert!((pz - 0.060).abs() < EPS, "Pz = {pz}");
    }

    #[test]
    fn zener_power_worst_case_is_no_load() {
        // The diode dissipates most when the load draws nothing, so all
        // of I_Rs flows through it: Pz = Vz * I_Rs.
        let reg = ZenerRegulator::new(5.0, 500.0).unwrap();
        let i_rs = reg.resistor_current_a(15.0).unwrap();
        let pz_noload = reg.zener_power_w(15.0, 0.0).unwrap();
        assert!((pz_noload - 5.0 * i_rs).abs() < EPS, "Pz = {pz_noload}");
        // Adding load can only reduce diode dissipation.
        let pz_loaded = reg.zener_power_w(15.0, 0.010).unwrap();
        assert!(pz_loaded < pz_noload);
    }

    #[test]
    fn zener_power_never_negative_on_overload() {
        let reg = ZenerRegulator::new(5.0, 500.0).unwrap();
        // Overloaded: signed Iz < 0 but reported power clamps to 0.
        let pz = reg.zener_power_w(15.0, 0.030).unwrap();
        assert!(pz.abs() < EPS, "Pz = {pz}");
    }

    // --- resistor power + energy balance ----------------------------

    #[test]
    fn resistor_power_matches_i_squared_r() {
        let reg = ZenerRegulator::new(5.0, 500.0).unwrap();
        // I_Rs = 20 mA -> P_Rs = 0.020^2 * 500 = 0.2 W.
        let p_rs = reg.resistor_power_w(15.0).unwrap();
        assert!((p_rs - 0.2).abs() < EPS, "P_Rs = {p_rs}");
    }

    #[test]
    fn energy_balance_source_equals_resistor_plus_zener_plus_load() {
        // Conservation: Vin * I_Rs = P_Rs + Pz + P_load (in regulation).
        let reg = ZenerRegulator::new(5.0, 500.0).unwrap();
        let vin = 15.0;
        let il = 0.008;
        let i_rs = reg.resistor_current_a(vin).unwrap();
        let p_source = vin * i_rs;
        let p_rs = reg.resistor_power_w(vin).unwrap();
        let pz = reg.zener_power_w(vin, il).unwrap();
        let p_load = reg.output_voltage_v(vin, il).unwrap() * il;
        let sum = p_rs + pz + p_load;
        assert!((p_source - sum).abs() < 1e-9, "source={p_source} sum={sum}");
    }

    // --- operating_point bundle -------------------------------------

    #[test]
    fn operating_point_bundles_consistent_values() {
        let reg = ZenerRegulator::new(5.1, 276.0).unwrap();
        let op = reg.operating_point(12.0, 0.020).unwrap();
        assert!(op.in_regulation);
        assert!((op.output_voltage_v - 5.1).abs() < EPS);
        // Cross-check each field against the standalone methods.
        let i_rs = reg.resistor_current_a(12.0).unwrap();
        assert!((op.resistor_current_a - i_rs).abs() < EPS);
        assert!((op.zener_current_a - reg.zener_current_a(12.0, 0.020).unwrap()).abs() < EPS);
        assert!((op.zener_power_w - reg.zener_power_w(12.0, 0.020).unwrap()).abs() < EPS);
        assert!((op.resistor_power_w - reg.resistor_power_w(12.0).unwrap()).abs() < EPS);
    }

    // --- constructor validation -------------------------------------

    #[test]
    fn constructor_rejects_bad_values() {
        assert!(ZenerRegulator::new(0.0, 500.0).is_err());
        assert!(ZenerRegulator::new(5.0, -1.0).is_err());
        assert!(ZenerRegulator::new(f64::NAN, 500.0).is_err());
        let err = ZenerRegulator::new(-5.0, 500.0).unwrap_err();
        assert_eq!(err.code(), "zener.bad_parameter");
    }

    #[test]
    fn methods_reject_nonfinite_inputs() {
        let reg = ZenerRegulator::new(5.0, 500.0).unwrap();
        assert!(reg.resistor_current_a(f64::NAN).is_err());
        assert!(reg.output_voltage_v(12.0, -0.001).is_err());
        assert!(reg.zener_current_a(12.0, f64::INFINITY).is_err());
    }
}
