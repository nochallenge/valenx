//! Reynolds number and flow-regime classification for internal
//! (pipe / duct) flow.
//!
//! The Reynolds number is the dimensionless ratio of inertial to
//! viscous forces,
//!
//! ```text
//!   Re = rho * v * D / mu
//! ```
//!
//! where `rho` is the fluid density (kg/m^3), `v` the bulk (mean)
//! velocity (m/s), `D` the internal pipe diameter (m), and `mu` the
//! dynamic viscosity (Pa.s). For a circular pipe the characteristic
//! length is the diameter; for non-circular ducts substitute the
//! hydraulic diameter `D_h = 4 A / P`.

use crate::error::{require_positive, PipeFlowError};
use serde::{Deserialize, Serialize};

/// The classical pipe-flow lower critical Reynolds number. Below this,
/// pipe flow is laminar in practice (textbook value, e.g. White,
/// *Fluid Mechanics*).
pub const RE_LAMINAR_UPPER: f64 = 2300.0;

/// The Reynolds number above which fully turbulent pipe flow is the
/// usual engineering assumption. Between [`RE_LAMINAR_UPPER`] and this
/// value the flow is transitional and intermittent.
pub const RE_TURBULENT_LOWER: f64 = 4000.0;

/// Qualitative pipe-flow regime, classified by Reynolds number.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum FlowRegime {
    /// `Re < 2300`: smooth, ordered, viscous-dominated flow. The
    /// friction factor follows the exact analytic law `f = 64/Re`.
    Laminar,
    /// `2300 <= Re < 4000`: intermittent, hard-to-predict transitional
    /// flow. No single friction correlation is reliable here.
    Transitional,
    /// `Re >= 4000`: fully turbulent flow, well described by the
    /// Colebrook-White / Haaland correlations.
    Turbulent,
}

impl FlowRegime {
    /// Classify a Reynolds number into a [`FlowRegime`] using the
    /// standard `2300` / `4000` thresholds.
    ///
    /// The thresholds are treated as: laminar for `Re < 2300`,
    /// transitional for `2300 <= Re < 4000`, turbulent for `Re >= 4000`.
    /// A non-positive or non-finite `re` is clamped to laminar, since
    /// the only physically meaningful Reynolds numbers are positive and
    /// the smallest of them are laminar.
    pub fn classify(re: f64) -> FlowRegime {
        // A non-finite (NaN/Inf) or non-positive Reynolds number is not a
        // physical flow; treat it as the (viscous, smallest-Re) laminar
        // case rather than panicking. Screen those first so the regime
        // boundaries below only ever see well-ordered positive values.
        if !re.is_finite() || re < RE_LAMINAR_UPPER {
            FlowRegime::Laminar
        } else if re < RE_TURBULENT_LOWER {
            FlowRegime::Transitional
        } else {
            FlowRegime::Turbulent
        }
    }

    /// `true` if the regime is [`FlowRegime::Laminar`].
    pub fn is_laminar(self) -> bool {
        matches!(self, FlowRegime::Laminar)
    }

    /// `true` if the regime is [`FlowRegime::Turbulent`].
    pub fn is_turbulent(self) -> bool {
        matches!(self, FlowRegime::Turbulent)
    }
}

/// Compute the Reynolds number `Re = rho * v * D / mu`.
///
/// All inputs are SI: `rho` in kg/m^3, `velocity` in m/s, `diameter` in
/// m, `viscosity` (dynamic) in Pa.s. Each must be finite and strictly
/// positive, otherwise a [`PipeFlowError`] is returned.
///
/// # Examples
///
/// ```
/// use valenx_pipeflow::reynolds::reynolds_number;
///
/// // Water at 20 C (rho = 998, mu = 1.002e-3) at 2 m/s in a 50 mm pipe.
/// let re = reynolds_number(998.0, 2.0, 0.05, 1.002e-3).unwrap();
/// assert!((re - 99_600.8).abs() < 1.0);
/// ```
pub fn reynolds_number(
    rho: f64,
    velocity: f64,
    diameter: f64,
    viscosity: f64,
) -> Result<f64, PipeFlowError> {
    let rho = require_positive("rho", rho)?;
    let velocity = require_positive("velocity", velocity)?;
    let diameter = require_positive("diameter", diameter)?;
    let viscosity = require_positive("viscosity", viscosity)?;
    Ok(rho * velocity * diameter / viscosity)
}

/// Compute the Reynolds number from the *kinematic* viscosity
/// `nu = mu / rho` (m^2/s): `Re = v * D / nu`.
///
/// Convenient when the fluid property to hand is the kinematic
/// viscosity (as is common for air and water tables).
pub fn reynolds_number_kinematic(
    velocity: f64,
    diameter: f64,
    kinematic_viscosity: f64,
) -> Result<f64, PipeFlowError> {
    let velocity = require_positive("velocity", velocity)?;
    let diameter = require_positive("diameter", diameter)?;
    let nu = require_positive("kinematic_viscosity", kinematic_viscosity)?;
    Ok(velocity * diameter / nu)
}

/// Classify the flow regime directly from the flow conditions, computing
/// the Reynolds number internally. Equivalent to
/// `FlowRegime::classify(reynolds_number(..)?)`.
pub fn classify_flow(
    rho: f64,
    velocity: f64,
    diameter: f64,
    viscosity: f64,
) -> Result<FlowRegime, PipeFlowError> {
    let re = reynolds_number(rho, velocity, diameter, viscosity)?;
    Ok(FlowRegime::classify(re))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `Re = rho v D / mu` against a hand-computed value.
    #[test]
    fn reynolds_matches_hand_calculation() {
        // rho=1000, v=1, D=0.1, mu=1e-3  =>  Re = 1000*1*0.1/1e-3 = 1e5.
        let re = reynolds_number(1000.0, 1.0, 0.1, 1.0e-3).unwrap();
        assert!((re - 1.0e5).abs() < 1e-6);
    }

    /// Kinematic and dynamic forms agree when `nu = mu/rho`.
    #[test]
    fn kinematic_and_dynamic_agree() {
        let rho = 998.0;
        let mu = 1.002e-3;
        let nu = mu / rho;
        let re_dyn = reynolds_number(rho, 2.0, 0.05, mu).unwrap();
        let re_kin = reynolds_number_kinematic(2.0, 0.05, nu).unwrap();
        assert!((re_dyn - re_kin).abs() < 1e-6);
    }

    /// Re scales linearly with velocity and with diameter, inversely
    /// with viscosity.
    #[test]
    fn reynolds_scaling_is_linear() {
        let base = reynolds_number(1000.0, 1.0, 0.1, 1.0e-3).unwrap();
        let double_v = reynolds_number(1000.0, 2.0, 0.1, 1.0e-3).unwrap();
        let double_d = reynolds_number(1000.0, 1.0, 0.2, 1.0e-3).unwrap();
        let half_mu = reynolds_number(1000.0, 1.0, 0.1, 0.5e-3).unwrap();
        assert!((double_v - 2.0 * base).abs() < 1e-6);
        assert!((double_d - 2.0 * base).abs() < 1e-6);
        assert!((half_mu - 2.0 * base).abs() < 1e-6);
    }

    /// Non-positive and non-finite inputs are rejected.
    #[test]
    fn reynolds_rejects_bad_inputs() {
        assert!(reynolds_number(0.0, 1.0, 0.1, 1e-3).is_err());
        assert!(reynolds_number(1000.0, -1.0, 0.1, 1e-3).is_err());
        assert!(reynolds_number(1000.0, 1.0, 0.0, 1e-3).is_err());
        assert!(reynolds_number(1000.0, 1.0, 0.1, f64::NAN).is_err());
    }

    /// Regime thresholds at exactly 2300 and 4000.
    #[test]
    fn regime_thresholds_are_correct() {
        assert_eq!(FlowRegime::classify(1000.0), FlowRegime::Laminar);
        assert_eq!(FlowRegime::classify(2299.0), FlowRegime::Laminar);
        // 2300 is the first transitional value.
        assert_eq!(FlowRegime::classify(2300.0), FlowRegime::Transitional);
        assert_eq!(FlowRegime::classify(3000.0), FlowRegime::Transitional);
        assert_eq!(FlowRegime::classify(3999.0), FlowRegime::Transitional);
        // 4000 is the first turbulent value.
        assert_eq!(FlowRegime::classify(4000.0), FlowRegime::Turbulent);
        assert_eq!(FlowRegime::classify(1.0e6), FlowRegime::Turbulent);
    }

    /// Just below / just above each boundary, to lock the open/closed
    /// interval convention.
    #[test]
    fn regime_boundaries_are_half_open() {
        // [0, 2300) laminar, [2300, 4000) transitional, [4000, inf) turbulent.
        assert!(FlowRegime::classify(2299.999).is_laminar());
        assert!(!FlowRegime::classify(2300.001).is_laminar());
        assert!(!FlowRegime::classify(3999.999).is_turbulent());
        assert!(FlowRegime::classify(4000.001).is_turbulent());
    }

    /// Degenerate Reynolds numbers fall back to laminar.
    #[test]
    fn regime_handles_degenerate_reynolds() {
        assert_eq!(FlowRegime::classify(0.0), FlowRegime::Laminar);
        assert_eq!(FlowRegime::classify(-5.0), FlowRegime::Laminar);
        assert_eq!(FlowRegime::classify(f64::NAN), FlowRegime::Laminar);
    }

    /// `classify_flow` agrees with classifying the explicit Reynolds.
    #[test]
    fn classify_flow_matches_two_step() {
        let re = reynolds_number(1000.0, 5.0, 0.1, 1e-3).unwrap();
        let direct = classify_flow(1000.0, 5.0, 0.1, 1e-3).unwrap();
        assert_eq!(direct, FlowRegime::classify(re));
        assert_eq!(direct, FlowRegime::Turbulent);
    }
}
