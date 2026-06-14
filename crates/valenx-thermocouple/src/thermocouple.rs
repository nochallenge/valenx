//! Linear-Seebeck thermocouple model: EMF from a junction pair, cold-junction
//! compensation, and the inverse temperature-from-voltage map.

use serde::{Deserialize, Serialize};

use crate::error::ThermocoupleError;

/// A named standard thermocouple type.
///
/// Each variant selects a representative near-room-temperature Seebeck
/// sensitivity via [`TcType::sensitivity_v_per_c`]. These are textbook
/// nominal values, not the full NIST ITS-90 reference functions (see the
/// crate-level *Honest scope* note).
#[derive(Copy, Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub enum TcType {
    /// Type K (chromel/alumel). The workhorse general-purpose
    /// thermocouple, ~41 uV/C near room temperature.
    K,
    /// Type J (iron/constantan), ~52 uV/C near room temperature.
    J,
    /// Type T (copper/constantan), ~41 uV/C near room temperature.
    T,
    /// Type E (chromel/constantan), the highest common output at
    /// ~61 uV/C near room temperature.
    E,
}

impl TcType {
    /// Representative Seebeck sensitivity `S` for this type, in volts per
    /// degree-Celsius (equivalently volts per kelvin).
    ///
    /// Values are nominal near-room-temperature figures: type K and T
    /// ~41 uV/C, type J ~52 uV/C, type E ~61 uV/C.
    pub fn sensitivity_v_per_c(self) -> f64 {
        match self {
            TcType::K => 41.0e-6,
            TcType::J => 52.0e-6,
            TcType::T => 41.0e-6,
            TcType::E => 61.0e-6,
        }
    }
}

/// A thermocouple characterised by a single (constant) Seebeck
/// sensitivity.
///
/// Construct one with [`Thermocouple::new`] from an explicit
/// sensitivity, or with [`Thermocouple::of_type`] from a named
/// [`TcType`]. All temperatures are in degrees Celsius and all voltages
/// in volts.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Thermocouple {
    /// Seebeck sensitivity `S` in volts per degree-Celsius. Strictly
    /// positive (enforced by the constructors).
    sensitivity_v_per_c: f64,
}

impl Thermocouple {
    /// Build a thermocouple from an explicit Seebeck sensitivity `S`
    /// (volts per degree-Celsius).
    ///
    /// # Errors
    ///
    /// Returns [`ThermocoupleError::NonFinite`] if `sensitivity_v_per_c`
    /// is `NaN` or infinite, and
    /// [`ThermocoupleError::NonPositiveSensitivity`] if it is zero or
    /// negative.
    pub fn new(sensitivity_v_per_c: f64) -> Result<Self, ThermocoupleError> {
        if !sensitivity_v_per_c.is_finite() {
            return Err(ThermocoupleError::NonFinite {
                name: "sensitivity_v_per_c",
                value: sensitivity_v_per_c,
            });
        }
        if sensitivity_v_per_c <= 0.0 {
            return Err(ThermocoupleError::NonPositiveSensitivity(
                sensitivity_v_per_c,
            ));
        }
        Ok(Self {
            sensitivity_v_per_c,
        })
    }

    /// Build a thermocouple from a named [`TcType`], using that type's
    /// representative sensitivity.
    ///
    /// This never fails because every preset sensitivity is a positive
    /// finite constant.
    pub fn of_type(kind: TcType) -> Self {
        // Safe to unwrap: all preset sensitivities are positive finite
        // constants, so `new` cannot reject them.
        Self::new(kind.sensitivity_v_per_c())
            .expect("preset sensitivities are always positive and finite")
    }

    /// The Seebeck sensitivity `S` in volts per degree-Celsius.
    pub fn sensitivity_v_per_c(&self) -> f64 {
        self.sensitivity_v_per_c
    }

    /// Open-circuit thermoelectric EMF of the junction pair,
    /// `EMF = S * (t_hot_c - t_cold_c)`, in volts.
    ///
    /// A hotter measurement junction than reference junction yields a
    /// positive EMF; equal junction temperatures yield exactly `0`.
    ///
    /// # Errors
    ///
    /// Returns [`ThermocoupleError::NonFinite`] if either temperature is
    /// `NaN` or infinite.
    pub fn emf(&self, t_hot_c: f64, t_cold_c: f64) -> Result<f64, ThermocoupleError> {
        Self::check_finite("t_hot_c", t_hot_c)?;
        Self::check_finite("t_cold_c", t_cold_c)?;
        Ok(self.sensitivity_v_per_c * (t_hot_c - t_cold_c))
    }

    /// Cold-junction-compensated EMF, i.e. the voltage an instrument
    /// reports relative to a `0` reference junction.
    ///
    /// This is the raw EMF across the leads plus the EMF the cold
    /// junction itself produces against `0`, which collapses to
    /// `S * t_hot_c`:
    ///
    /// ```text
    /// V_comp = S * (t_hot_c - t_cold_c) + S * (t_cold_c - 0) = S * t_hot_c
    /// ```
    ///
    /// # Errors
    ///
    /// Returns [`ThermocoupleError::NonFinite`] if either temperature is
    /// `NaN` or infinite.
    pub fn emf_compensated(&self, t_hot_c: f64, t_cold_c: f64) -> Result<f64, ThermocoupleError> {
        let measured = self.emf(t_hot_c, t_cold_c)?;
        // Add back the reference-junction contribution against 0.
        let reference = self.emf(t_cold_c, 0.0)?;
        Ok(measured + reference)
    }

    /// Invert a measured EMF back to the hot-junction temperature,
    /// `t_hot_c = t_cold_c + emf_volts / S`, in degrees Celsius.
    ///
    /// This is the exact inverse of [`Thermocouple::emf`] for the same
    /// `t_cold_c`.
    ///
    /// # Errors
    ///
    /// Returns [`ThermocoupleError::NonFinite`] if `emf_volts` or
    /// `t_cold_c` is `NaN` or infinite.
    pub fn temperature_from_emf(
        &self,
        emf_volts: f64,
        t_cold_c: f64,
    ) -> Result<f64, ThermocoupleError> {
        Self::check_finite("emf_volts", emf_volts)?;
        Self::check_finite("t_cold_c", t_cold_c)?;
        Ok(t_cold_c + emf_volts / self.sensitivity_v_per_c)
    }

    /// Internal guard rejecting non-finite measurement inputs.
    fn check_finite(name: &'static str, value: f64) -> Result<(), ThermocoupleError> {
        if value.is_finite() {
            Ok(())
        } else {
            Err(ThermocoupleError::NonFinite { name, value })
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::error::ErrorCategory;

    /// Tight tolerance for exact closed-form comparisons.
    const EPS: f64 = 1e-12;

    fn type_k() -> Thermocouple {
        Thermocouple::of_type(TcType::K)
    }

    #[test]
    fn type_k_sensitivity_is_about_41_microvolts() {
        // Ground truth: nominal type-K sensitivity ~41 uV/C.
        let s = TcType::K.sensitivity_v_per_c();
        assert!((s - 41.0e-6).abs() < EPS, "type-K S = {s}");
        assert!((type_k().sensitivity_v_per_c() - 41.0e-6).abs() < EPS);
    }

    #[test]
    fn emf_equals_sensitivity_times_delta_t() {
        // VALIDATE: EMF = S * dT against a hand-computed value.
        let tc = type_k();
        let dt = 100.0;
        let emf = tc.emf(125.0, 25.0).expect("valid junctions");
        let expected = 41.0e-6 * dt; // 4.1 mV
        assert!(
            (emf - expected).abs() < EPS,
            "emf = {emf}, expected {expected}"
        );
        assert!((emf - 0.0041).abs() < EPS, "emf = {emf}");
    }

    #[test]
    fn zero_delta_t_gives_zero_emf() {
        // VALIDATE: zero dT -> 0 EMF, at several reference temperatures.
        let tc = type_k();
        for &t in &[-40.0, 0.0, 25.0, 300.0, 1000.0] {
            let emf = tc.emf(t, t).expect("valid junctions");
            assert!(emf.abs() < EPS, "emf at equal junctions ({t}) = {emf}");
        }
    }

    #[test]
    fn cold_junction_compensation_adds_reference() {
        // VALIDATE: cold-junction compensation adds the reference EMF, and
        // the compensated reading collapses to S * t_hot.
        let tc = type_k();
        let t_hot = 200.0;
        let t_cold = 25.0;
        let measured = tc.emf(t_hot, t_cold).expect("valid");
        let comp = tc.emf_compensated(t_hot, t_cold).expect("valid");

        // The compensated value must exceed the raw value by exactly the
        // reference-junction EMF against 0 (here positive, t_cold > 0).
        let reference = tc.emf(t_cold, 0.0).expect("valid");
        assert!(reference > 0.0, "reference EMF = {reference}");
        assert!(
            ((comp - measured) - reference).abs() < EPS,
            "comp - measured = {}, reference = {reference}",
            comp - measured
        );

        // And the closed form V_comp = S * t_hot.
        let expected = tc.sensitivity_v_per_c() * t_hot;
        assert!(
            (comp - expected).abs() < EPS,
            "comp = {comp}, expected {expected}"
        );
    }

    #[test]
    fn compensation_is_identity_when_reference_is_zero() {
        // With a 0 C reference there is nothing to add back.
        let tc = type_k();
        let measured = tc.emf(150.0, 0.0).expect("valid");
        let comp = tc.emf_compensated(150.0, 0.0).expect("valid");
        assert!(
            (comp - measured).abs() < EPS,
            "comp = {comp}, measured = {measured}"
        );
    }

    #[test]
    fn temperature_from_emf_inverts_emf() {
        // VALIDATE: T from V inverts EMF over a sweep of hot/cold pairs.
        let tc = type_k();
        for &(t_hot, t_cold) in &[
            (125.0, 25.0),
            (0.0, 0.0),
            (-50.0, 10.0),
            (980.0, 23.0),
            (37.0, 37.0),
        ] {
            let emf = tc.emf(t_hot, t_cold).expect("valid");
            let recovered = tc.temperature_from_emf(emf, t_cold).expect("valid emf");
            assert!(
                (recovered - t_hot).abs() < 1e-9,
                "recovered {recovered} from ({t_hot}, {t_cold})"
            );
        }
    }

    #[test]
    fn temperature_from_compensated_emf_recovers_hot() {
        // The compensated EMF inverts against a 0 reference back to t_hot.
        let tc = type_k();
        let t_hot = 314.0;
        let comp = tc.emf_compensated(t_hot, 25.0).expect("valid");
        let recovered = tc.temperature_from_emf(comp, 0.0).expect("valid");
        assert!((recovered - t_hot).abs() < 1e-9, "recovered {recovered}");
    }

    #[test]
    fn higher_delta_t_gives_higher_emf() {
        // VALIDATE: higher dT -> higher EMF (strict monotonicity).
        let tc = type_k();
        let mut prev = f64::NEG_INFINITY;
        for &dt in &[0.0, 10.0, 50.0, 100.0, 500.0, 1000.0] {
            let emf = tc.emf(dt, 0.0).expect("valid");
            assert!(emf > prev, "emf {emf} at dt {dt} not above prev {prev}");
            prev = emf;
        }
    }

    #[test]
    fn negative_delta_t_gives_negative_emf() {
        // Sign convention: colder hot junction than reference -> negative.
        let tc = type_k();
        let emf = tc.emf(0.0, 100.0).expect("valid");
        assert!(emf < 0.0, "emf = {emf}");
        assert!((emf + 0.0041).abs() < EPS, "emf = {emf}");
    }

    #[test]
    fn emf_scales_linearly_with_sensitivity() {
        // Double the sensitivity -> double the EMF for the same dT.
        let lo = Thermocouple::new(20.0e-6).expect("valid S");
        let hi = Thermocouple::new(40.0e-6).expect("valid S");
        let a = lo.emf(100.0, 0.0).expect("valid");
        let b = hi.emf(100.0, 0.0).expect("valid");
        assert!((b - 2.0 * a).abs() < EPS, "a = {a}, b = {b}");
    }

    #[test]
    fn distinct_types_have_expected_sensitivities() {
        assert!((TcType::J.sensitivity_v_per_c() - 52.0e-6).abs() < EPS);
        assert!((TcType::T.sensitivity_v_per_c() - 41.0e-6).abs() < EPS);
        assert!((TcType::E.sensitivity_v_per_c() - 61.0e-6).abs() < EPS);
        // Type E has the largest nominal output, type K/T the smallest.
        assert!(TcType::E.sensitivity_v_per_c() > TcType::J.sensitivity_v_per_c());
        assert!(TcType::J.sensitivity_v_per_c() > TcType::K.sensitivity_v_per_c());
    }

    #[test]
    fn new_rejects_non_positive_sensitivity() {
        for &s in &[0.0, -1.0e-6, -1.0] {
            let err = Thermocouple::new(s).expect_err("should reject");
            assert_eq!(err.code(), "thermocouple.non_positive_sensitivity");
            assert_eq!(err.category(), ErrorCategory::Config);
        }
    }

    #[test]
    fn new_rejects_non_finite_sensitivity() {
        for s in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let err = Thermocouple::new(s).expect_err("should reject");
            assert_eq!(err.code(), "thermocouple.non_finite");
            assert_eq!(err.category(), ErrorCategory::Input);
        }
    }

    #[test]
    fn emf_rejects_non_finite_temperatures() {
        let tc = type_k();
        assert!(tc.emf(f64::NAN, 25.0).is_err());
        assert!(tc.emf(25.0, f64::INFINITY).is_err());
        let err = tc.emf(f64::NAN, 0.0).expect_err("should reject");
        assert_eq!(err.code(), "thermocouple.non_finite");
    }

    #[test]
    fn temperature_from_emf_rejects_non_finite_inputs() {
        let tc = type_k();
        assert!(tc.temperature_from_emf(f64::NAN, 25.0).is_err());
        assert!(tc.temperature_from_emf(0.001, f64::NEG_INFINITY).is_err());
    }

    #[test]
    fn serde_round_trips_thermocouple() {
        let tc = Thermocouple::of_type(TcType::E);
        let json = serde_json::to_string(&tc).expect("serialize");
        let back: Thermocouple = serde_json::from_str(&json).expect("deserialize");
        assert_eq!(tc, back);
        assert!((back.sensitivity_v_per_c() - 61.0e-6).abs() < EPS);
    }

    #[test]
    fn kelvin_difference_matches_celsius_difference() {
        // The model is linear in temperature difference, so shifting both
        // junctions by the 273.15 C->K offset leaves the EMF unchanged.
        let tc = type_k();
        let celsius = tc.emf(125.0, 25.0).expect("valid");
        let kelvin = tc.emf(125.0 + 273.15, 25.0 + 273.15).expect("valid");
        assert!(
            (celsius - kelvin).abs() < EPS,
            "C = {celsius}, K = {kelvin}"
        );
    }
}
