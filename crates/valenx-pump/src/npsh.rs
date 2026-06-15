//! **Net positive suction head (NPSH).**
//!
//! Cavitation is avoided when the head available at the pump suction
//! exceeds the head the pump requires. The *available* NPSH is
//!
//! ```text
//!   NPSHa = (P_atm - P_vap) / (rho * g) + H_s - H_loss      [m]
//! ```
//!
//! where `P_atm` is the absolute pressure on the source liquid surface,
//! `P_vap` the fluid's vapour pressure, `H_s` the static suction head
//! (positive for a flooded suction, **negative** for a suction lift), and
//! `H_loss ≥ 0` the friction loss in the suction line. Raising the
//! suction lift (making `H_s` more negative) lowers `NPSHa`, pushing the
//! pump towards cavitation.
//!
//! The pump is safe when the **margin** `NPSHa − NPSHr` is positive,
//! where `NPSHr` is the pump's required NPSH datum.

use serde::{Deserialize, Serialize};

use crate::error::{require_finite, require_non_negative, require_positive, PumpError};
use crate::G;

/// The suction-side conditions feeding the NPSH calculation.
///
/// Pressures are absolute, in pascals. Heads are in metres of the pumped
/// fluid. `static_suction_head_m` carries the sign convention: positive
/// when the source surface is **above** the pump (flooded), negative for
/// a lift.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct SuctionConditions {
    /// Absolute pressure on the source liquid surface `P_atm`, in pascals.
    pub atmospheric_pa: f64,
    /// Fluid vapour pressure `P_vap` at the operating temperature, in
    /// pascals.
    pub vapor_pressure_pa: f64,
    /// Fluid density `ρ`, in kg/m³.
    pub density_kg_m3: f64,
    /// Static suction head `H_s`, in metres. Positive = flooded suction
    /// (surface above the pump); negative = suction lift.
    pub static_suction_head_m: f64,
    /// Suction-line friction loss `H_loss ≥ 0`, in metres.
    pub suction_loss_m: f64,
}

impl SuctionConditions {
    /// Build a validated set of suction conditions.
    ///
    /// `atmospheric_pa`, `vapor_pressure_pa` and `suction_loss_m` must be
    /// finite and non-negative; `density_kg_m3` must be finite and
    /// strictly positive (it appears in a denominator);
    /// `static_suction_head_m` must be finite (any sign). The vapour
    /// pressure must not exceed the surface pressure, otherwise the
    /// liquid is already boiling at the source.
    ///
    /// # Errors
    ///
    /// Returns [`PumpError::BadParameter`] for an out-of-domain scalar,
    /// or [`PumpError::Inconsistent`] if `vapor_pressure_pa >
    /// atmospheric_pa`.
    pub fn new(
        atmospheric_pa: f64,
        vapor_pressure_pa: f64,
        density_kg_m3: f64,
        static_suction_head_m: f64,
        suction_loss_m: f64,
    ) -> Result<Self, PumpError> {
        let atmospheric_pa = require_non_negative("atmospheric_pa", atmospheric_pa)?;
        let vapor_pressure_pa = require_non_negative("vapor_pressure_pa", vapor_pressure_pa)?;
        let density_kg_m3 = require_positive("density_kg_m3", density_kg_m3)?;
        let static_suction_head_m = require_finite("static_suction_head_m", static_suction_head_m)?;
        let suction_loss_m = require_non_negative("suction_loss_m", suction_loss_m)?;
        if vapor_pressure_pa > atmospheric_pa {
            return Err(PumpError::Inconsistent(format!(
                "vapour pressure {vapor_pressure_pa} Pa exceeds surface pressure {atmospheric_pa} Pa (liquid is boiling)"
            )));
        }
        Ok(Self {
            atmospheric_pa,
            vapor_pressure_pa,
            density_kg_m3,
            static_suction_head_m,
            suction_loss_m,
        })
    }

    /// The pressure margin above vapour expressed as a head, in metres:
    /// `(P_atm − P_vap) / (ρ·g)`. This is the barometric term of NPSHa,
    /// before the static head and line loss are applied.
    pub fn pressure_head_m(&self) -> f64 {
        (self.atmospheric_pa - self.vapor_pressure_pa) / (self.density_kg_m3 * G)
    }
}

/// Available net positive suction head, in metres of fluid:
/// `(P_atm − P_vap)/(ρ·g) + H_s − H_loss`.
///
/// # Examples
///
/// ```
/// use valenx_pump::npsh::{available_npsh_m, SuctionConditions};
///
/// // Cold water (vp ≈ 2.34 kPa) at sea level, flooded 2 m, 0.5 m loss.
/// let c = SuctionConditions::new(101_325.0, 2_340.0, 1000.0, 2.0, 0.5).unwrap();
/// let npsha = available_npsh_m(&c);
/// let barometric = (101_325.0 - 2_340.0) / (1000.0 * 9.806_65);
/// assert!((npsha - (barometric + 2.0 - 0.5)).abs() < 1e-9);
/// ```
pub fn available_npsh_m(conditions: &SuctionConditions) -> f64 {
    conditions.pressure_head_m() + conditions.static_suction_head_m - conditions.suction_loss_m
}

/// The cavitation margin `NPSHa − NPSHr`, in metres.
///
/// A positive margin means the available suction head exceeds what the
/// pump requires, so it will not cavitate at this duty.
///
/// # Errors
///
/// Returns [`PumpError::BadParameter`] if `required_npsh_m` is not finite
/// and non-negative.
pub fn npsh_margin_m(
    conditions: &SuctionConditions,
    required_npsh_m: f64,
) -> Result<f64, PumpError> {
    let required_npsh_m = require_non_negative("required_npsh_m", required_npsh_m)?;
    Ok(available_npsh_m(conditions) - required_npsh_m)
}

/// Whether the pump is free of cavitation at this duty, i.e. the
/// [`npsh_margin_m`] is greater than or equal to zero.
///
/// # Errors
///
/// Returns [`PumpError::BadParameter`] if `required_npsh_m` is not finite
/// and non-negative.
pub fn is_cavitation_free(
    conditions: &SuctionConditions,
    required_npsh_m: f64,
) -> Result<bool, PumpError> {
    Ok(npsh_margin_m(conditions, required_npsh_m)? >= 0.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    /// Sea-level barometric head over cold-water vapour pressure, the
    /// recurring constant in these tests.
    fn barometric() -> f64 {
        (101_325.0 - 2_340.0) / (1000.0 * G)
    }

    #[test]
    fn npsha_matches_closed_form() {
        let c = SuctionConditions::new(101_325.0, 2_340.0, 1000.0, 2.0, 0.5).unwrap();
        let expected = barometric() + 2.0 - 0.5;
        assert!((available_npsh_m(&c) - expected).abs() < EPS);
        // Sanity: barometric head of water at 1 atm is ~10.07 m.
        assert!((barometric() - 10.094).abs() < 1e-2);
    }

    #[test]
    fn suction_lift_lowers_npsha() {
        // Flooded +2 m vs lift -3 m, all else equal.
        let flooded = SuctionConditions::new(101_325.0, 2_340.0, 1000.0, 2.0, 0.5).unwrap();
        let lift = SuctionConditions::new(101_325.0, 2_340.0, 1000.0, -3.0, 0.5).unwrap();
        let drop = available_npsh_m(&flooded) - available_npsh_m(&lift);
        // The only difference is the static head term: 2 - (-3) = 5 m.
        assert!((drop - 5.0).abs() < EPS);
        assert!(available_npsh_m(&lift) < available_npsh_m(&flooded));
    }

    #[test]
    fn deeper_lift_falls_monotonically() {
        let mut prev = f64::INFINITY;
        for lift in [-1.0, -2.0, -3.0, -4.0, -5.0] {
            let c = SuctionConditions::new(101_325.0, 2_340.0, 1000.0, lift, 0.3).unwrap();
            let v = available_npsh_m(&c);
            assert!(v < prev, "NPSHa must fall as the lift deepens");
            prev = v;
        }
    }

    #[test]
    fn line_loss_lowers_npsha_one_for_one() {
        let low = SuctionConditions::new(101_325.0, 2_340.0, 1000.0, 0.0, 0.2).unwrap();
        let high = SuctionConditions::new(101_325.0, 2_340.0, 1000.0, 0.0, 1.7).unwrap();
        // 1.5 m more loss -> 1.5 m less NPSHa.
        assert!((available_npsh_m(&low) - available_npsh_m(&high) - 1.5).abs() < EPS);
    }

    #[test]
    fn margin_and_cavitation_flag_agree() {
        let c = SuctionConditions::new(101_325.0, 2_340.0, 1000.0, 1.0, 0.5).unwrap();
        let npsha = available_npsh_m(&c);
        // Require comfortably less than available -> positive margin, safe.
        let safe_req = npsha - 3.0;
        assert!(
            (npsh_margin_m(&c, safe_req.max(0.0)).unwrap() - (npsha - safe_req.max(0.0))).abs()
                < EPS
        );
        assert!(is_cavitation_free(&c, 2.0).unwrap());
        // Require more than available -> cavitates.
        assert!(!is_cavitation_free(&c, npsha + 1.0).unwrap());
    }

    #[test]
    fn margin_exactly_zero_is_cavitation_free() {
        let c = SuctionConditions::new(101_325.0, 2_340.0, 1000.0, 0.0, 0.0).unwrap();
        let npsha = available_npsh_m(&c);
        assert!((npsh_margin_m(&c, npsha).unwrap()).abs() < EPS);
        assert!(is_cavitation_free(&c, npsha).unwrap());
    }

    #[test]
    fn boiling_source_is_rejected() {
        let err = SuctionConditions::new(2_000.0, 3_000.0, 1000.0, 0.0, 0.0).unwrap_err();
        assert_eq!(err.code(), "pump.inconsistent");
    }

    #[test]
    fn rejects_non_positive_density() {
        assert!(SuctionConditions::new(101_325.0, 2_340.0, 0.0, 0.0, 0.0).is_err());
    }
}
