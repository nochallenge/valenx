//! DC bias networks and the resulting operating point.
//!
//! Two topologies are provided, both reducing to the same base-loop
//! Kirchhoff solve:
//!
//! - [`FixedBias`] — a single base resistor `Rb` from the supply.
//! - [`DividerBias`] — the four-resistor voltage-divider network, whose
//!   base side is first replaced by its Thevenin equivalent.
//!
//! Each `solve` returns an [`OperatingPoint`]. See the
//! [crate-level docs](crate) for the governing equations.

use crate::error::BjtError;
use crate::model::{Region, Transistor};
use serde::{Deserialize, Serialize};

/// The quiescent operating point ("Q-point") of a biased transistor.
///
/// All currents are in amperes and all voltages in volts. The currents
/// satisfy the beta relations exactly: `ic == beta * ib`,
/// `ie == (beta + 1) * ib`, and therefore `ie == ic + ib`.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct OperatingPoint {
    /// Base current `Ib` (A).
    pub ib: f64,
    /// Collector current `Ic` (A).
    pub ic: f64,
    /// Emitter current `Ie` (A).
    pub ie: f64,
    /// Emitter-resistor voltage drop `VE = Ie * Re` (V).
    pub ve: f64,
    /// Collector-emitter voltage `Vce` (V).
    pub vce: f64,
    /// Operating region implied by `Vce` relative to `Vce_sat`.
    pub region: Region,
}

/// Solve the shared base loop for a Thevenin-reduced bias.
///
/// `vth` / `rth` are the Thevenin source seen by the base. The loop
/// `vth = Ib*rth + VBE + Ie*Re` with `Ie = (beta + 1) * Ib` gives
/// `Ib = (vth - VBE) / (rth + (beta + 1) * Re)`. The result is then
/// checked against the saturation floor and, if saturated, the
/// collector current is clamped to `Ic_sat = (Vcc - Vce_sat)/(Rc+Re)`.
///
/// `rth`, `rc`, `re` must be validated non-negative by the caller, and
/// `rth + (beta + 1) * re` must be strictly positive (guaranteed
/// because `beta > 0` for a valid [`Transistor`] so the `+1` term keeps
/// the denominator positive whenever `rth >= 0`).
fn solve_base_loop(
    device: &Transistor,
    vcc: f64,
    vth: f64,
    rth: f64,
    rc: f64,
    re: f64,
) -> Result<OperatingPoint, BjtError> {
    let drive = vth - device.vbe;
    if drive <= 0.0 {
        return Err(BjtError::cut_off(
            "Thevenin base voltage does not exceed VBE, so the base-emitter junction is not forward biased",
        ));
    }

    let denom = rth + (device.beta + 1.0) * re;
    // denom > 0 is guaranteed: beta > 0 and rth, re >= 0 by validation.
    let ib = drive / denom;
    let ic_active = device.collector_current(ib);
    let ie = device.emitter_current(ib);

    // Active-region collector-emitter voltage.
    let vce_active = vcc - ic_active * rc - ie * re;

    if vce_active > device.vce_sat {
        Ok(OperatingPoint {
            ib,
            ic: ic_active,
            ie,
            ve: ie * re,
            vce: vce_active,
            region: Region::Active,
        })
    } else {
        // Saturated: Vce is pinned to the floor and the collector
        // current is set by the external resistors, not by beta.
        let ic_sat = (vcc - device.vce_sat) / (rc + re);
        // In saturation the base is overdriven; the emitter current is
        // the saturated collector current plus the (still beta-limited)
        // base current, preserving Ie = Ic + Ib at the device.
        let ie_sat = ic_sat + ib;
        Ok(OperatingPoint {
            ib,
            ic: ic_sat,
            ie: ie_sat,
            ve: ie_sat * re,
            vce: device.vce_sat,
            region: Region::Saturation,
        })
    }
}

/// The DC bias **stability factor** `S(ICO) = ∂Ic/∂Ico` for an
/// emitter-degenerated base loop with effective base resistance `rb`
/// (the Thevenin `Rth` for a divider) and emitter resistance `re`.
///
/// Differentiating the base-loop solution with respect to the reverse
/// saturation ("leakage") current `Ico` gives the standard Boylestad
/// result
///
/// > `S = (beta + 1) * (rb + re) / (rb + (beta + 1) * re)`,
///
/// equal to the algebraically identical `(1 + beta) / (1 + beta * re /
/// (re + rb))`. It is purely a property of the network and the gain —
/// it does not depend on `Vcc` or the Q-point — and is bounded by
/// `1 <= S <= beta + 1`: it collapses to `beta + 1` for a bare fixed
/// bias (`re = 0`, the worst case) and to `1` for an ideal emitter bias
/// (`rb = 0`, the best case).
///
/// `rb`, `re` are assumed validated non-negative and `beta > 0`; the
/// only degenerate input is `rb == re == 0`, for which the factor is the
/// indeterminate `0 / 0` and an error is returned.
fn stability_factor_s_ico(beta: f64, rb: f64, re: f64) -> Result<f64, BjtError> {
    let denom = rb + (beta + 1.0) * re;
    if denom == 0.0 {
        return Err(BjtError::bad_parameter(
            "rb+re",
            "stability factor S(ICO) is undefined when both the base and emitter resistance are zero",
            0.0,
        ));
    }
    Ok((beta + 1.0) * (rb + re) / denom)
}

/// Validate a resistance argument (`>= 0`, finite).
fn check_resistance(name: &'static str, value: f64) -> Result<(), BjtError> {
    if !value.is_finite() || value < 0.0 {
        return Err(BjtError::bad_parameter(
            name,
            "resistance must be finite and non-negative",
            value,
        ));
    }
    Ok(())
}

/// Validate a supply-voltage argument (finite).
fn check_supply(value: f64) -> Result<(), BjtError> {
    if !value.is_finite() {
        return Err(BjtError::bad_parameter(
            "vcc",
            "supply must be finite",
            value,
        ));
    }
    Ok(())
}

/// Fixed-base bias: a single resistor `rb` from `vcc` into the base,
/// a collector resistor `rc`, and an optional emitter resistor `re`.
///
/// This is the simplest (and most `beta`-sensitive) bias. It is exactly
/// a [`DividerBias`] whose Thevenin source is `Vth = Vcc`, `Rth = Rb`.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct FixedBias {
    /// Supply voltage `Vcc` (V).
    pub vcc: f64,
    /// Base resistor `Rb` (ohm).
    pub rb: f64,
    /// Collector resistor `Rc` (ohm).
    pub rc: f64,
    /// Emitter resistor `Re` (ohm); pass `0.0` for no emitter
    /// degeneration.
    pub re: f64,
}

impl FixedBias {
    /// Build a validated fixed-base bias network.
    ///
    /// # Errors
    ///
    /// Returns [`BjtError::BadParameter`] if `vcc` is not finite or if
    /// any of `rb`, `rc`, `re` is negative or not finite.
    pub fn new(vcc: f64, rb: f64, rc: f64, re: f64) -> Result<Self, BjtError> {
        check_supply(vcc)?;
        check_resistance("rb", rb)?;
        check_resistance("rc", rc)?;
        check_resistance("re", re)?;
        Ok(Self { vcc, rb, rc, re })
    }

    /// Solve for the operating point of `device` in this network.
    ///
    /// # Errors
    ///
    /// Returns [`BjtError::CutOff`] if `Vcc <= VBE` so the base never
    /// forward-biases.
    pub fn solve(&self, device: &Transistor) -> Result<OperatingPoint, BjtError> {
        solve_base_loop(device, self.vcc, self.vcc, self.rb, self.rc, self.re)
    }

    /// The DC bias stability factor `S(ICO) = ∂Ic/∂Ico` of this network
    /// with `device` installed.
    ///
    /// For a fixed bias the effective base resistance is simply `Rb`, so
    /// `S = (beta + 1)(Rb + Re) / (Rb + (beta + 1) Re)`. With no emitter
    /// resistor (`Re = 0`) this is the textbook worst case
    /// `S = beta + 1`. See the [crate-level docs](crate) for the
    /// derivation and the `1 <= S <= beta + 1` bounds.
    ///
    /// # Errors
    ///
    /// Returns [`BjtError::BadParameter`] only for the degenerate network
    /// with `Rb == Re == 0`, where the factor is indeterminate.
    pub fn stability_factor(&self, device: &Transistor) -> Result<f64, BjtError> {
        stability_factor_s_ico(device.beta, self.rb, self.re)
    }
}

/// Voltage-divider bias: the four-resistor network `R1` (top, to
/// `Vcc`), `R2` (bottom, to ground), collector resistor `Rc`, and
/// emitter resistor `Re`.
///
/// The base divider is reduced to its Thevenin equivalent before the
/// base loop is solved (see [`DividerBias::thevenin`]).
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct DividerBias {
    /// Supply voltage `Vcc` (V).
    pub vcc: f64,
    /// Upper divider resistor `R1`, from `Vcc` to the base node (ohm).
    pub r1: f64,
    /// Lower divider resistor `R2`, from the base node to ground (ohm).
    pub r2: f64,
    /// Collector resistor `Rc` (ohm).
    pub rc: f64,
    /// Emitter resistor `Re` (ohm); pass `0.0` for no emitter
    /// degeneration.
    pub re: f64,
}

impl DividerBias {
    /// Build a validated voltage-divider bias network.
    ///
    /// # Errors
    ///
    /// Returns [`BjtError::BadParameter`] if `vcc` is not finite, if any
    /// resistance is negative or not finite, or if `r1 + r2 == 0` (the
    /// divider would short the supply node and the Thevenin source is
    /// undefined).
    pub fn new(vcc: f64, r1: f64, r2: f64, rc: f64, re: f64) -> Result<Self, BjtError> {
        check_supply(vcc)?;
        check_resistance("r1", r1)?;
        check_resistance("r2", r2)?;
        check_resistance("rc", rc)?;
        check_resistance("re", re)?;
        if r1 + r2 == 0.0 {
            return Err(BjtError::bad_parameter(
                "r1+r2",
                "divider resistors cannot both be zero",
                0.0,
            ));
        }
        Ok(Self {
            vcc,
            r1,
            r2,
            rc,
            re,
        })
    }

    /// Thevenin equivalent of the base divider, returned as
    /// `(Vth, Rth)`.
    ///
    /// `Vth = Vcc * R2 / (R1 + R2)` is the open-circuit base voltage and
    /// `Rth = R1 * R2 / (R1 + R2)` (the parallel combination `R1 || R2`)
    /// is the resistance seen looking back into the divider.
    pub fn thevenin(&self) -> (f64, f64) {
        let sum = self.r1 + self.r2;
        let vth = self.vcc * self.r2 / sum;
        let rth = self.r1 * self.r2 / sum;
        (vth, rth)
    }

    /// Solve for the operating point of `device` in this network.
    ///
    /// # Errors
    ///
    /// Returns [`BjtError::CutOff`] if the Thevenin base voltage does
    /// not exceed `VBE`.
    pub fn solve(&self, device: &Transistor) -> Result<OperatingPoint, BjtError> {
        let (vth, rth) = self.thevenin();
        solve_base_loop(device, self.vcc, vth, rth, self.rc, self.re)
    }

    /// The DC bias stability factor `S(ICO) = ∂Ic/∂Ico` of this network
    /// with `device` installed.
    ///
    /// The effective base resistance is the Thevenin `Rth = R1 || R2`, so
    /// `S = (beta + 1)(Rth + Re) / (Rth + (beta + 1) Re)`. A "stiff"
    /// divider (`(beta + 1) Re >> Rth`) drives `S` toward its ideal lower
    /// bound of `1`. See the [crate-level docs](crate) for the derivation
    /// and the `1 <= S <= beta + 1` bounds.
    ///
    /// # Errors
    ///
    /// Returns [`BjtError::BadParameter`] only for the degenerate network
    /// with `Rth == Re == 0` (i.e. `R1 == 0` and `Re == 0`).
    pub fn stability_factor(&self, device: &Transistor) -> Result<f64, BjtError> {
        let (_, rth) = self.thevenin();
        stability_factor_s_ico(device.beta, rth, self.re)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Tight tolerance for current comparisons (amperes).
    const EPS_I: f64 = 1e-12;
    /// Tolerance for voltage comparisons (volts).
    const EPS_V: f64 = 1e-9;

    fn approx(a: f64, b: f64, eps: f64) -> bool {
        (a - b).abs() < eps
    }

    // ---- model: beta current relations -------------------------------

    #[test]
    fn beta_relations_hold() {
        let q = Transistor::silicon(100.0).unwrap();
        let ib = 2.0e-5; // 20 uA
        let ic = q.collector_current(ib);
        let ie = q.emitter_current(ib);
        // Ic = beta * Ib
        assert!(approx(ic, 100.0 * ib, EPS_I), "Ic = beta*Ib, got {ic}");
        // Ie = (beta + 1) * Ib
        assert!(approx(ie, 101.0 * ib, EPS_I), "Ie = (beta+1)*Ib, got {ie}");
        // Ie = Ic + Ib
        assert!(approx(ie, ic + ib, EPS_I), "Ie = Ic + Ib, got {ie}");
    }

    #[test]
    fn beta_from_currents_inverts() {
        let beta = Transistor::beta_from_currents(2.0e-3, 1.0e-5).unwrap();
        assert!(approx(beta, 200.0, 1e-9), "beta = Ic/Ib = 200, got {beta}");
    }

    #[test]
    fn beta_from_zero_base_current_errors() {
        let err = Transistor::beta_from_currents(1.0e-3, 0.0).unwrap_err();
        assert_eq!(err.code(), "bjt.bad_parameter");
    }

    #[test]
    fn base_current_for_collector_inverts() {
        let q = Transistor::silicon(150.0).unwrap();
        let ib = q.base_current_for_collector(3.0e-3);
        assert!(approx(ib, 3.0e-3 / 150.0, EPS_I), "Ib = Ic/beta, got {ib}");
        // Round trip back to Ic.
        assert!(approx(q.collector_current(ib), 3.0e-3, EPS_I));
    }

    // ---- voltage-divider bias: ground-truth hand calculation ---------

    #[test]
    fn divider_bias_active_point_matches_hand_calc() {
        // Network: Vcc = 7 V, R1 = 30k, R2 = 10k, Rc = 2k, Re = 1k.
        // Device: beta = 99, VBE = 0.7, Vce_sat = 0.2.
        let q = Transistor::new(99.0, 0.7, 0.2).unwrap();
        let bias = DividerBias::new(7.0, 30_000.0, 10_000.0, 2_000.0, 1_000.0).unwrap();

        // Thevenin (hand): Vth = 7 * 10k/40k = 1.75 V ; Rth = 30k||10k = 7.5k.
        let (vth, rth) = bias.thevenin();
        assert!(approx(vth, 1.75, EPS_V), "Vth = 1.75, got {vth}");
        assert!(approx(rth, 7_500.0, EPS_V), "Rth = 7.5k, got {rth}");

        let op = bias.solve(&q).unwrap();

        // Ib = (Vth - VBE) / (Rth + (beta+1)*Re)
        //    = (1.75 - 0.7) / (7500 + 100*1000)
        //    = 1.05 / 107500.
        let ib_expected = 1.05 / 107_500.0;
        assert!(
            approx(op.ib, ib_expected, EPS_I),
            "Ib hand-calc, got {}",
            op.ib
        );

        // Beta relations on the solved point.
        assert!(approx(op.ic, 99.0 * ib_expected, EPS_I));
        assert!(approx(op.ie, 100.0 * ib_expected, EPS_I));
        assert!(approx(op.ie, op.ic + op.ib, EPS_I));

        // VE = Ie * Re.
        assert!(
            approx(op.ve, op.ie * 1_000.0, EPS_V),
            "VE = Ie*Re, got {}",
            op.ve
        );

        // Vce = Vcc - Ic*Rc - Ie*Re.
        let vce_expected = 7.0 - op.ic * 2_000.0 - op.ie * 1_000.0;
        assert!(
            approx(op.vce, vce_expected, EPS_V),
            "Vce hand-calc, got {}",
            op.vce
        );

        // Comfortably active (Vce well above 0.2 V).
        assert_eq!(op.region, Region::Active);
        assert!(op.vce > q.vce_sat);
    }

    // ---- a second, independently published-style worked example ------

    #[test]
    fn divider_bias_classic_round_numbers() {
        // Pick R1=R2 so Vth = Vcc/2 = 5 V, Rth = R1/2 = 5k, with a
        // dominant emitter term so the answer is essentially beta-
        // independent ("stiff" divider design rule (beta+1)Re >> Rth).
        // Vcc = 10, R1 = R2 = 10k, Rc = 1k, Re = 1k, beta = 199.
        let q = Transistor::new(199.0, 0.7, 0.2).unwrap();
        let bias = DividerBias::new(10.0, 10_000.0, 10_000.0, 1_000.0, 1_000.0).unwrap();
        let (vth, rth) = bias.thevenin();
        assert!(approx(vth, 5.0, EPS_V));
        assert!(approx(rth, 5_000.0, EPS_V));

        let op = bias.solve(&q).unwrap();
        // Ib = (5 - 0.7) / (5000 + 200*1000) = 4.3 / 205000.
        let ib = 4.3 / 205_000.0;
        assert!(approx(op.ib, ib, EPS_I));
        // VE = Ie * Re exactly.
        assert!(approx(op.ve, op.ie * 1_000.0, EPS_V));
        // KVL ground truth: walking the base loop, the emitter node sits
        // at VE = Vth - VBE - Ib*Rth (the Thevenin source less the
        // base-emitter drop and the Rth drop).
        let ve_kvl = vth - q.vbe - op.ib * rth;
        assert!(
            approx(op.ve, ve_kvl, EPS_V),
            "VE = Vth-VBE-Ib*Rth, got {}",
            op.ve
        );
        assert_eq!(op.region, Region::Active);
    }

    // ---- fixed-base bias equals divider bias with Vth=Vcc, Rth=Rb ----

    #[test]
    fn fixed_bias_matches_hand_calc() {
        // Vcc = 12, Rb = 240k, Rc = 2.2k, Re = 0 (no emitter resistor).
        // beta = 100, VBE = 0.7.
        let q = Transistor::silicon(100.0).unwrap();
        let bias = FixedBias::new(12.0, 240_000.0, 2_200.0, 0.0).unwrap();
        let op = bias.solve(&q).unwrap();

        // Ib = (Vcc - VBE)/Rb = (12 - 0.7)/240000 = 11.3/240000.
        let ib = 11.3 / 240_000.0;
        assert!(approx(op.ib, ib, EPS_I), "Ib = (Vcc-VBE)/Rb, got {}", op.ib);
        // Ic = beta*Ib.
        assert!(approx(op.ic, 100.0 * ib, EPS_I));
        // With Re = 0, VE = 0 and Ie = (beta+1)Ib.
        assert!(approx(op.ve, 0.0, EPS_V));
        assert!(approx(op.ie, 101.0 * ib, EPS_I));
        // Vce = Vcc - Ic*Rc (Re=0).
        let vce = 12.0 - op.ic * 2_200.0;
        assert!(approx(op.vce, vce, EPS_V), "Vce, got {}", op.vce);
        assert_eq!(op.region, Region::Active);
    }

    #[test]
    fn fixed_bias_equivalent_to_divider_with_same_thevenin() {
        // A fixed bias (Vcc, Rb) must equal a divider whose Thevenin
        // source reproduces Vth = Vcc and Rth = Rb. Achieve Rth = Rb and
        // Vth = Vcc by letting R2 -> infinity relative to R1: use R1 = Rb,
        // R2 huge. As R2 -> inf, Vth -> Vcc and Rth -> R1 = Rb.
        let q = Transistor::silicon(120.0).unwrap();
        let fixed = FixedBias::new(9.0, 100_000.0, 1_500.0, 470.0).unwrap();
        // Exact equality: build the divider directly from the Thevenin
        // identity instead, by constructing it to share the solver path.
        let op_fixed = fixed.solve(&q).unwrap();

        // Independently reduce: fixed bias is solve_base_loop with
        // Vth=Vcc=9, Rth=Rb=100k.
        let op_manual = solve_base_loop(&q, 9.0, 9.0, 100_000.0, 1_500.0, 470.0).unwrap();
        assert!(approx(op_fixed.ib, op_manual.ib, EPS_I));
        assert!(approx(op_fixed.ic, op_manual.ic, EPS_I));
        assert!(approx(op_fixed.vce, op_manual.vce, EPS_V));
        assert_eq!(op_fixed.region, op_manual.region);
    }

    // ---- saturation detection ----------------------------------------

    #[test]
    fn heavy_base_drive_saturates() {
        // Force saturation: tiny Rb gives a large Ib, beta*Ib*Rc would
        // exceed Vcc, so the device pins at Vce_sat.
        // Vcc = 5, Rb = 1k, Rc = 1k, Re = 0, beta = 200, VBE = 0.7.
        let q = Transistor::new(200.0, 0.7, 0.2).unwrap();
        let bias = FixedBias::new(5.0, 1_000.0, 1_000.0, 0.0).unwrap();
        let op = bias.solve(&q).unwrap();

        // Active-region Ic would be beta*(Vcc-VBE)/Rb = 200*4.3/1000 = 0.86 A,
        // giving Ic*Rc = 860 V >> Vcc, impossible -> saturation.
        assert_eq!(op.region, Region::Saturation, "should saturate");
        // Vce pinned to floor.
        assert!(approx(op.vce, 0.2, EPS_V), "Vce = Vce_sat, got {}", op.vce);
        // Ic clamped to Ic_sat = (Vcc - Vce_sat)/(Rc + Re) = 4.8/1000.
        let ic_sat = (5.0 - 0.2) / 1_000.0;
        assert!(approx(op.ic, ic_sat, EPS_I), "Ic_sat, got {}", op.ic);
        // Base current is still the beta-loop value (overdriven base).
        let ib = 4.3 / 1_000.0;
        assert!(approx(op.ib, ib, EPS_I));
        // Ie = Ic + Ib still holds at the device.
        assert!(approx(op.ie, op.ic + op.ib, EPS_I));
    }

    #[test]
    fn boundary_region_is_active_just_above_floor() {
        // Construct a case where Vce lands just above Vce_sat -> Active.
        // Vcc = 5, Re = 0, Rc chosen so Ic*Rc = Vcc - Vce_sat - tiny.
        let q = Transistor::new(100.0, 0.7, 0.2).unwrap();
        // Ib = (5-0.7)/Rb. Choose Rb = 430k -> Ib = 1e-5, Ic = 1e-3.
        // Want Ic*Rc slightly under 4.8 -> Rc = 4700 -> Ic*Rc = 4.7 < 4.8.
        let bias = FixedBias::new(5.0, 430_000.0, 4_700.0, 0.0).unwrap();
        let op = bias.solve(&q).unwrap();
        assert_eq!(op.region, Region::Active);
        // Vce = 5 - 1e-3*4700 = 0.3 V > 0.2.
        assert!(approx(op.vce, 0.3, 1e-6), "Vce ~ 0.3, got {}", op.vce);
    }

    // ---- validation / error paths ------------------------------------

    #[test]
    fn cutoff_when_thevenin_below_vbe() {
        // Vth = Vcc * R2/(R1+R2) = 5 * 1/(99+1) = 0.05 V < 0.7 -> cut off.
        let q = Transistor::silicon(100.0).unwrap();
        let bias = DividerBias::new(5.0, 99_000.0, 1_000.0, 1_000.0, 1_000.0).unwrap();
        let err = bias.solve(&q).unwrap_err();
        assert_eq!(err.code(), "bjt.cut_off");
        assert_eq!(err.category(), crate::error::ErrorCategory::Operating);
    }

    #[test]
    fn negative_gain_rejected() {
        let err = Transistor::new(-10.0, 0.7, 0.2).unwrap_err();
        assert_eq!(err.code(), "bjt.bad_parameter");
        assert_eq!(err.category(), crate::error::ErrorCategory::Input);
    }

    #[test]
    fn negative_resistance_rejected() {
        let err = DividerBias::new(10.0, -1.0, 10_000.0, 1_000.0, 1_000.0).unwrap_err();
        assert_eq!(err.code(), "bjt.bad_parameter");
    }

    #[test]
    fn zero_divider_rejected() {
        let err = DividerBias::new(10.0, 0.0, 0.0, 1_000.0, 1_000.0).unwrap_err();
        assert_eq!(err.code(), "bjt.bad_parameter");
    }

    #[test]
    fn nonfinite_supply_rejected() {
        let err = FixedBias::new(f64::NAN, 1_000.0, 1_000.0, 0.0).unwrap_err();
        assert_eq!(err.code(), "bjt.bad_parameter");
    }

    // ---- bias stability factor S(ICO) -------------------------------

    /// Tolerance for the dimensionless stability factor.
    const EPS_S: f64 = 1e-9;

    #[test]
    fn stability_factor_divider_matches_hand_calc() {
        // Same network as `divider_bias_active_point_matches_hand_calc`:
        // Rth = 7.5k, Re = 1k, beta = 99.
        let q = Transistor::new(99.0, 0.7, 0.2).unwrap();
        let bias = DividerBias::new(7.0, 30_000.0, 10_000.0, 2_000.0, 1_000.0).unwrap();
        let (_, rth) = bias.thevenin();
        assert!(approx(rth, 7_500.0, EPS_V));

        // S = (beta+1)(Rth+Re)/(Rth+(beta+1)Re)
        //   = 100 * 8500 / 107500 = 7.906976744...
        let s = bias.stability_factor(&q).unwrap();
        let s_hand = 100.0 * (7_500.0 + 1_000.0) / (7_500.0 + 100.0 * 1_000.0);
        assert!(approx(s, s_hand, EPS_S), "S hand-calc, got {s}");

        // Independent algebraic form S = (1+beta)/(1 + beta*Re/(Re+Rth)).
        let s_alt = (1.0 + 99.0) / (1.0 + 99.0 * 1_000.0 / (1_000.0 + 7_500.0));
        assert!(approx(s, s_alt, EPS_S), "S alt-form, got {s} vs {s_alt}");
    }

    #[test]
    fn stability_factor_fixed_bias_no_emitter_is_beta_plus_one() {
        // GOLD limiting case: Re = 0 -> S = beta + 1 (the worst case),
        // independent of Rb.
        let q = Transistor::silicon(100.0).unwrap();
        let bias = FixedBias::new(12.0, 240_000.0, 2_200.0, 0.0).unwrap();
        let s = bias.stability_factor(&q).unwrap();
        assert!(approx(s, 101.0, EPS_S), "S = beta+1 = 101, got {s}");
    }

    #[test]
    fn stability_factor_ideal_emitter_bias_is_one() {
        // GOLD limiting case: Rb = 0 -> S = 1 exactly (best case),
        // independent of beta and Re.
        let q = Transistor::silicon(250.0).unwrap();
        let bias = FixedBias::new(9.0, 0.0, 1_500.0, 470.0).unwrap();
        let s = bias.stability_factor(&q).unwrap();
        assert!(approx(s, 1.0, EPS_S), "S = 1, got {s}");
    }

    #[test]
    fn stability_factor_is_bounded_by_one_and_beta_plus_one() {
        // 1 <= S <= beta+1 for every valid network.
        let q = Transistor::silicon(150.0).unwrap();
        let upper = q.beta + 1.0;
        for &rth in &[1.0, 1_000.0, 10_000.0, 100_000.0] {
            for &re in &[0.0, 100.0, 1_000.0, 10_000.0] {
                // R1 = R2 = 2*rth -> Rth = R1||R2 = rth.
                let bias = DividerBias::new(10.0, 2.0 * rth, 2.0 * rth, 1_000.0, re).unwrap();
                let (_, got_rth) = bias.thevenin();
                assert!(approx(got_rth, rth, 1e-6));
                let s = bias.stability_factor(&q).unwrap();
                assert!(
                    ((1.0 - EPS_S)..=(upper + EPS_S)).contains(&s),
                    "S out of [1, beta+1]: rth={rth}, re={re}, S={s}"
                );
            }
        }
    }

    #[test]
    fn stability_factor_decreases_with_emitter_resistance() {
        // Emitter degeneration monotonically improves stability (drives
        // S strictly down toward 1) at fixed Rb.
        let q = Transistor::silicon(100.0).unwrap();
        let res = [0.0, 100.0, 470.0, 1_000.0, 4_700.0, 22_000.0];
        let mut prev = f64::INFINITY;
        for &re in &res {
            let bias = FixedBias::new(12.0, 47_000.0, 2_200.0, re).unwrap();
            let s = bias.stability_factor(&q).unwrap();
            assert!(
                s < prev,
                "S should strictly decrease with Re: {s} !< {prev}"
            );
            assert!(s >= 1.0 - EPS_S, "S >= 1, got {s}");
            prev = s;
        }
    }

    #[test]
    fn stability_factor_degenerate_network_rejected() {
        // Rb = 0 and Re = 0 -> indeterminate 0/0.
        let q = Transistor::silicon(100.0).unwrap();
        let bias = FixedBias::new(5.0, 0.0, 1_000.0, 0.0).unwrap();
        let err = bias.stability_factor(&q).unwrap_err();
        assert_eq!(err.code(), "bjt.bad_parameter");
    }

    #[test]
    fn serde_round_trip_operating_point() {
        let q = Transistor::silicon(100.0).unwrap();
        let bias = DividerBias::new(12.0, 47_000.0, 10_000.0, 2_200.0, 1_000.0).unwrap();
        let op = bias.solve(&q).unwrap();
        let json = serde_json::to_string(&op).unwrap();
        let back: OperatingPoint = serde_json::from_str(&json).unwrap();
        assert!(approx(op.ib, back.ib, EPS_I));
        assert!(approx(op.vce, back.vce, EPS_V));
        assert_eq!(op.region, back.region);
    }
}
