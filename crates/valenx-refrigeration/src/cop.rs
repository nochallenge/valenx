//! Coefficient-of-performance (COP) definitions.
//!
//! The COP of a vapor-compression machine is the useful energy effect
//! divided by the work that drives the cycle. Which effect is "useful"
//! depends on the application:
//!
//! When the machine is run as a **refrigerator / air conditioner** the
//! useful effect is the heat absorbed in the evaporator, the
//! refrigerating effect `Q_evap`, so the cooling coefficient of
//! performance is
//!
//! ```text
//! COP_cool = Q_evap / W_comp
//! ```
//!
//! When the same machine is run as a **heat pump** the useful effect is
//! the heat rejected at the condenser `Q_cond`, so the heating
//! coefficient of performance is
//!
//! ```text
//! COP_heat = Q_cond / W_comp
//! ```
//!
//! For an idealised closed cycle exchanging heat with only two
//! reservoirs the steady-flow energy balance is
//! `Q_cond = Q_evap + W_comp`, from which the two coefficients differ by
//! exactly one:
//!
//! ```text
//! COP_heat = (Q_evap + W_comp) / W_comp = COP_cool + 1
//! ```
//!
//! The same numbers can be formed directly from the four cycle-corner
//! specific enthalpies: with state 1 the evaporator outlet (compressor
//! inlet), state 2 the compressor outlet, state 3 the condenser outlet
//! and state 4 the throttle outlet, an isenthalpic expansion gives
//! `h4 = h3`, so
//!
//! ```text
//! COP_cool = (h1 - h4) / (h2 - h1)
//! ```
//!
//! All enthalpies and works here are specific quantities (per unit mass,
//! conventionally kJ/kg), and all heats are non-negative magnitudes.

use crate::error::{RefrigError, Result};

/// Cooling coefficient of performance from the refrigerating effect and
/// the compressor work, `COP_cool = Q_evap / W_comp`.
///
/// `q_evap` is the heat absorbed in the evaporator (non-negative) and
/// `w_comp` is the compressor work input (strictly positive). Both are
/// specific quantities in the same energy-per-mass unit.
///
/// # Errors
///
/// Returns [`RefrigError::BadParameter`] if `q_evap` is negative or
/// non-finite, and [`RefrigError::NonPositiveWork`] if `w_comp` is not
/// strictly positive.
pub fn cop_cool(q_evap: f64, w_comp: f64) -> Result<f64> {
    let q = RefrigError::require_non_negative("q_evap", q_evap)?;
    let w = require_work(w_comp)?;
    Ok(q / w)
}

/// Heating coefficient of performance from the condenser heat rejection
/// and the compressor work, `COP_heat = Q_cond / W_comp`.
///
/// `q_cond` is the heat rejected at the condenser (non-negative) and
/// `w_comp` is the compressor work input (strictly positive).
///
/// # Errors
///
/// Returns [`RefrigError::BadParameter`] if `q_cond` is negative or
/// non-finite, and [`RefrigError::NonPositiveWork`] if `w_comp` is not
/// strictly positive.
pub fn cop_heat(q_cond: f64, w_comp: f64) -> Result<f64> {
    let q = RefrigError::require_non_negative("q_cond", q_cond)?;
    let w = require_work(w_comp)?;
    Ok(q / w)
}

/// Heating COP derived from the cooling COP via the energy-balance
/// identity `COP_heat = COP_cool + 1`.
///
/// This holds for any single-stage cycle exchanging heat with two
/// reservoirs because the condenser rejects exactly the refrigerating
/// effect plus the compressor work. It is a thin, total helper used by
/// callers that already hold a cooling COP.
///
/// `cop_cool` must be finite and non-negative (a physical cooling COP is
/// never negative).
///
/// # Errors
///
/// Returns [`RefrigError::BadParameter`] if `cop_cool` is negative or
/// non-finite.
pub fn cop_heat_from_cool(cop_cool: f64) -> Result<f64> {
    let c = RefrigError::require_non_negative("cop_cool", cop_cool)?;
    Ok(c + 1.0)
}

/// Cooling COP from the cycle-corner specific enthalpies,
/// `COP_cool = (h1 - h4) / (h2 - h1)`.
///
/// The arguments are the specific enthalpies at the four cycle corners:
///
/// - `h1`: evaporator outlet / compressor inlet (saturated or superheated vapor),
/// - `h2`: compressor outlet (high-pressure superheated vapor),
/// - `h4`: throttle (expansion-valve) outlet entering the evaporator.
///
/// The refrigerating effect is `h1 - h4` and the compression work is
/// `h2 - h1`. For a physically valid cooling cycle both differences must
/// be strictly positive.
///
/// # Errors
///
/// Returns [`RefrigError::BadParameter`] if any enthalpy is non-finite,
/// and [`RefrigError::InconsistentCycle`] if the refrigerating effect or
/// the compression work is not strictly positive.
pub fn cop_cool_from_enthalpies(h1: f64, h2: f64, h4: f64) -> Result<f64> {
    let h1 = RefrigError::require_finite("h1", h1)?;
    let h2 = RefrigError::require_finite("h2", h2)?;
    let h4 = RefrigError::require_finite("h4", h4)?;

    let refrigerating_effect = h1 - h4;
    let compression_work = h2 - h1;

    if compression_work <= 0.0 {
        return Err(RefrigError::InconsistentCycle(
            "compression work h2 - h1 must be strictly positive (h2 > h1)",
        ));
    }
    if refrigerating_effect <= 0.0 {
        return Err(RefrigError::InconsistentCycle(
            "refrigerating effect h1 - h4 must be strictly positive (h1 > h4)",
        ));
    }
    Ok(refrigerating_effect / compression_work)
}

/// Heating COP from the cycle-corner specific enthalpies,
/// `COP_heat = (h2 - h3) / (h2 - h1)`.
///
/// The arguments are the specific enthalpies at the relevant corners:
///
/// - `h1`: compressor inlet,
/// - `h2`: compressor outlet (condenser inlet),
/// - `h3`: condenser outlet (saturated or sub-cooled liquid).
///
/// The condenser heat rejection is `h2 - h3` and the compression work is
/// `h2 - h1`. Because the throttle is isenthalpic (`h4 = h3`), this is
/// numerically equal to [`cop_cool_from_enthalpies`] plus one.
///
/// # Errors
///
/// Returns [`RefrigError::BadParameter`] if any enthalpy is non-finite,
/// and [`RefrigError::InconsistentCycle`] if the heat rejection or the
/// compression work is not strictly positive.
pub fn cop_heat_from_enthalpies(h1: f64, h2: f64, h3: f64) -> Result<f64> {
    let h1 = RefrigError::require_finite("h1", h1)?;
    let h2 = RefrigError::require_finite("h2", h2)?;
    let h3 = RefrigError::require_finite("h3", h3)?;

    let heat_rejected = h2 - h3;
    let compression_work = h2 - h1;

    if compression_work <= 0.0 {
        return Err(RefrigError::InconsistentCycle(
            "compression work h2 - h1 must be strictly positive (h2 > h1)",
        ));
    }
    if heat_rejected <= 0.0 {
        return Err(RefrigError::InconsistentCycle(
            "condenser heat rejection h2 - h3 must be strictly positive (h2 > h3)",
        ));
    }
    Ok(heat_rejected / compression_work)
}

/// Shared guard turning a non-positive compressor work into the
/// dedicated [`RefrigError::NonPositiveWork`] variant.
fn require_work(w_comp: f64) -> Result<f64> {
    if w_comp.is_finite() && w_comp > 0.0 {
        Ok(w_comp)
    } else {
        Err(RefrigError::NonPositiveWork(w_comp))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Absolute tolerance for floating-point comparisons in this module.
    const EPS: f64 = 1e-9;

    #[test]
    fn cop_cool_basic_ratio() {
        // 12 kW refrigerating effect on 4 kW of work → COP 3.
        let cop = cop_cool(12.0, 4.0).unwrap();
        assert!((cop - 3.0).abs() < EPS, "got {cop}");
    }

    #[test]
    fn cop_heat_basic_ratio() {
        // Condenser rejects 16 kW on 4 kW of work → COP 4.
        let cop = cop_heat(16.0, 4.0).unwrap();
        assert!((cop - 4.0).abs() < EPS, "got {cop}");
    }

    #[test]
    fn heating_cop_is_cooling_cop_plus_one() {
        // Energy balance: with Q_evap = 12, W = 4, the condenser must
        // reject Q_cond = 16, so COP_heat = COP_cool + 1 exactly.
        let q_evap = 12.0;
        let w = 4.0;
        let q_cond = q_evap + w;

        let cc = cop_cool(q_evap, w).unwrap();
        let ch = cop_heat(q_cond, w).unwrap();
        assert!((ch - (cc + 1.0)).abs() < EPS, "cc={cc} ch={ch}");

        // And the dedicated helper agrees.
        let ch2 = cop_heat_from_cool(cc).unwrap();
        assert!((ch2 - ch).abs() < EPS, "ch={ch} ch2={ch2}");
    }

    #[test]
    fn cop_from_enthalpies_matches_textbook() {
        // A classic R-134a worked example (Cengel & Boles, "Thermodynamics",
        // ideal vapor-compression cycle, -20 C / 0.8 MPa):
        //   h1 = 239.16, h2 = 275.39, h3 = h4 = 95.47 kJ/kg.
        // Refrigerating effect 143.69, work 36.23 → COP_cool ≈ 3.966.
        let h1 = 239.16;
        let h2 = 275.39;
        let h3 = 95.47;
        let h4 = h3; // isenthalpic throttle

        let cc = cop_cool_from_enthalpies(h1, h2, h4).unwrap();
        let expected = (h1 - h4) / (h2 - h1);
        assert!((cc - expected).abs() < EPS, "cc={cc}");
        // Sanity against the hand-computed textbook figure.
        assert!((cc - 3.9660_f64).abs() < 1e-3, "cc={cc}");
    }

    #[test]
    fn enthalpy_heating_equals_cooling_plus_one() {
        // With an isenthalpic throttle (h4 = h3), the enthalpy form must
        // reproduce the +1 identity to machine precision.
        let h1 = 239.16;
        let h2 = 275.39;
        let h3 = 95.47;
        let h4 = h3;

        let cc = cop_cool_from_enthalpies(h1, h2, h4).unwrap();
        let ch = cop_heat_from_enthalpies(h1, h2, h3).unwrap();
        assert!((ch - (cc + 1.0)).abs() < EPS, "cc={cc} ch={ch}");
    }

    #[test]
    fn rejects_non_positive_work() {
        assert!(matches!(
            cop_cool(10.0, 0.0),
            Err(RefrigError::NonPositiveWork(_))
        ));
        assert!(matches!(
            cop_heat(10.0, -2.0),
            Err(RefrigError::NonPositiveWork(_))
        ));
    }

    #[test]
    fn rejects_negative_heat() {
        assert!(matches!(
            cop_cool(-1.0, 4.0),
            Err(RefrigError::BadParameter { .. })
        ));
    }

    #[test]
    fn enthalpy_rejects_unordered_states() {
        // h2 <= h1 → no compression work.
        assert!(matches!(
            cop_cool_from_enthalpies(250.0, 250.0, 95.0),
            Err(RefrigError::InconsistentCycle(_))
        ));
        // h1 <= h4 → no refrigerating effect.
        assert!(matches!(
            cop_cool_from_enthalpies(95.0, 275.0, 95.0),
            Err(RefrigError::InconsistentCycle(_))
        ));
        // NaN enthalpy.
        assert!(matches!(
            cop_cool_from_enthalpies(f64::NAN, 275.0, 95.0),
            Err(RefrigError::BadParameter { .. })
        ));
    }
}
