//! The Terzaghi bearing-capacity equation and the allowable pressure.
//!
//! This module assembles the soil properties, footing geometry, and
//! bearing-capacity factors into the three-term Terzaghi ultimate
//! bearing capacity and divides by a factor of safety to obtain the
//! allowable bearing pressure.
//!
//! # Equation (general shear, strip footing)
//!
//! `qult = c * Nc + q * Nq + 0.5 * gamma * B * Ngamma`
//!
//! where `q = gamma * Df` is the surcharge at the founding plane and
//! `Nc`, `Nq`, `Ngamma` come from [`crate::factors`].
//!
//! `qall = qult / FS`

use serde::{Deserialize, Serialize};

use crate::error::SoilBearingError;
use crate::factors::BearingFactors;
use crate::footing::Footing;
use crate::soil::SoilProperties;

/// Fully decomposed result of a bearing-capacity calculation.
///
/// Carries the individual term contributions alongside the totals so a
/// caller can inspect which mechanism (cohesion, surcharge, or
/// self-weight) dominates.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct BearingResult {
    /// The bearing-capacity factors used (`Nc`, `Nq`, `Ngamma`).
    pub factors: BearingFactors,
    /// Cohesion contribution `c * Nc` (stress units).
    pub cohesion_term: f64,
    /// Surcharge contribution `q * Nq` (stress units).
    pub surcharge_term: f64,
    /// Self-weight contribution `0.5 * gamma * B * Ngamma` (stress units).
    pub self_weight_term: f64,
    /// Ultimate bearing capacity `qult`, the sum of the three terms.
    pub q_ultimate: f64,
    /// Allowable bearing pressure `qall = qult / FS`.
    pub q_allowable: f64,
    /// The factor of safety applied to obtain `q_allowable`.
    pub factor_of_safety: f64,
}

/// Compute the ultimate Terzaghi bearing capacity `qult` (no factor of
/// safety applied).
///
/// `qult = c * Nc + q * Nq + 0.5 * gamma * B * Ngamma`.
///
/// This is infallible: both `soil` and `footing` are already validated
/// by their constructors, and the factors are finite for any
/// admissible friction angle.
///
/// # Examples
///
/// ```
/// use valenx_soilbearing::{ultimate_bearing_capacity, Footing, SoilProperties};
///
/// // Undrained clay, phi = 0: qult = c*Nc + q (Nq = 1, Ngamma = 0).
/// let clay = SoilProperties::new(0.0, 50.0, 18.0).unwrap();
/// let surface = Footing::new(2.0, 0.0).unwrap();
/// let qult = ultimate_bearing_capacity(&clay, &surface);
/// // c * (pi + 2) with q = 0 and the self-weight term = 0.
/// assert!((qult - 50.0 * (std::f64::consts::PI + 2.0)).abs() < 1e-9);
/// ```
pub fn ultimate_bearing_capacity(soil: &SoilProperties, footing: &Footing) -> f64 {
    ultimate_terms(soil, footing).3
}

/// Internal: return `(cohesion_term, surcharge_term, self_weight_term, qult)`.
fn ultimate_terms(soil: &SoilProperties, footing: &Footing) -> (f64, f64, f64, f64) {
    let factors = BearingFactors::from_friction_angle(soil);
    let q = footing.surcharge(soil.unit_weight());

    let cohesion_term = soil.cohesion() * factors.nc;
    let surcharge_term = q * factors.nq;
    let self_weight_term = 0.5 * soil.unit_weight() * footing.width() * factors.ngamma;
    let qult = cohesion_term + surcharge_term + self_weight_term;

    (cohesion_term, surcharge_term, self_weight_term, qult)
}

/// Compute the allowable bearing pressure `qall = qult / FS`, returning
/// the full term-by-term [`BearingResult`].
///
/// # Errors
///
/// Returns [`SoilBearingError::NotFinite`] if `factor_of_safety` is not
/// finite, or [`SoilBearingError::InvalidParameter`] if it is not
/// strictly greater than `1` (a factor of safety at or below unity
/// provides no margin and is rejected).
///
/// # Examples
///
/// ```
/// use valenx_soilbearing::{bearing_capacity, Footing, SoilProperties};
///
/// let clay = SoilProperties::new(0.0, 50.0, 18.0).unwrap();
/// let surface = Footing::new(2.0, 0.0).unwrap();
/// let r = bearing_capacity(&clay, &surface, 3.0).unwrap();
/// // qall is exactly qult / FS.
/// assert!((r.q_allowable - r.q_ultimate / 3.0).abs() < 1e-9);
/// assert!(bearing_capacity(&clay, &surface, 1.0).is_err()); // FS must exceed 1
/// ```
pub fn bearing_capacity(
    soil: &SoilProperties,
    footing: &Footing,
    factor_of_safety: f64,
) -> Result<BearingResult, SoilBearingError> {
    let factor_of_safety = SoilBearingError::require_finite("factor_of_safety", factor_of_safety)?;
    if factor_of_safety <= 1.0 {
        return Err(SoilBearingError::invalid(
            "factor_of_safety",
            factor_of_safety,
            "factor of safety must be strictly greater than 1",
        ));
    }

    let factors = BearingFactors::from_friction_angle(soil);
    let (cohesion_term, surcharge_term, self_weight_term, q_ultimate) =
        ultimate_terms(soil, footing);
    let q_allowable = q_ultimate / factor_of_safety;

    Ok(BearingResult {
        factors,
        cohesion_term,
        surcharge_term,
        self_weight_term,
        q_ultimate,
        q_allowable,
        factor_of_safety,
    })
}

/// Apply a factor of safety to a known ultimate capacity:
/// `qall = qult / FS`.
///
/// A small convenience for callers that already hold a `qult` value and
/// only need the allowable pressure.
///
/// # Errors
///
/// Same validation as [`bearing_capacity`]: `factor_of_safety` must be
/// finite and strictly greater than `1`. `q_ultimate` must also be
/// finite.
///
/// # Examples
///
/// ```
/// use valenx_soilbearing::allowable_from_ultimate;
///
/// let qall = allowable_from_ultimate(300.0, 3.0).unwrap();
/// assert!((qall - 100.0).abs() < 1e-12);
/// ```
pub fn allowable_from_ultimate(
    q_ultimate: f64,
    factor_of_safety: f64,
) -> Result<f64, SoilBearingError> {
    let q_ultimate = SoilBearingError::require_finite("q_ultimate", q_ultimate)?;
    let factor_of_safety = SoilBearingError::require_finite("factor_of_safety", factor_of_safety)?;
    if factor_of_safety <= 1.0 {
        return Err(SoilBearingError::invalid(
            "factor_of_safety",
            factor_of_safety,
            "factor of safety must be strictly greater than 1",
        ));
    }
    Ok(q_ultimate / factor_of_safety)
}
