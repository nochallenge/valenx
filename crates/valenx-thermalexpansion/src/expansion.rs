//! Linear, areal and volumetric thermal expansion.
//!
//! For a small temperature change `dT` an isotropic solid with a constant
//! linear coefficient of thermal expansion `alpha` changes its dimensions
//! to first order as
//!
//! ```text
//! dL = alpha * L0 * dT                  (length)
//! dA = (2 * alpha) * A0 * dT            (area)
//! dV = (3 * alpha) * V0 * dT            (volume)
//! ```
//!
//! The areal coefficient is therefore `~ 2 * alpha` and the volumetric
//! coefficient `~ 3 * alpha`. These factors are the leading terms of
//! `(1 + alpha dT)^2` and `(1 + alpha dT)^3`; the cross terms are
//! `O((alpha dT)^2)` and are dropped, which is the standard textbook
//! linearisation valid for the small `alpha dT` of ordinary materials.
//!
//! `alpha` carries units of inverse kelvin (1/K), `dT` is a temperature
//! change in kelvin (equivalently degrees Celsius — only the *difference*
//! matters), and lengths / areas / volumes are returned in whatever unit
//! the reference dimension was given in.

use crate::error::{require_finite, require_positive, ThermalError};
use serde::{Deserialize, Serialize};

/// The areal expansion coefficient as a multiple of the linear coefficient
/// `alpha`. To first order the areal coefficient is `AREA_FACTOR * alpha`.
pub const AREA_FACTOR: f64 = 2.0;

/// The volumetric expansion coefficient as a multiple of the linear
/// coefficient `alpha`. To first order the volumetric coefficient is
/// `VOLUME_FACTOR * alpha`.
pub const VOLUME_FACTOR: f64 = 3.0;

/// Linear coefficient of thermal expansion, in inverse kelvin (1/K).
///
/// A newtype so a coefficient cannot be silently confused with a length or
/// a temperature. Construct one with [`LinearCoefficient::new`], which
/// rejects non-positive and non-finite values.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct LinearCoefficient(f64);

impl LinearCoefficient {
    /// Build a validated linear coefficient of thermal expansion.
    ///
    /// # Errors
    ///
    /// Returns [`ThermalError::NonPositive`] if `alpha <= 0` and
    /// [`ThermalError::NonFinite`] if `alpha` is `NaN` or infinite. A
    /// physically meaningful coefficient for a solid that expands on
    /// heating is strictly positive.
    pub fn new(alpha: f64) -> Result<Self, ThermalError> {
        Ok(Self(require_positive("alpha", alpha)?))
    }

    /// The linear coefficient value in 1/K.
    pub fn per_kelvin(self) -> f64 {
        self.0
    }

    /// The corresponding areal coefficient `~ 2 * alpha`, in 1/K.
    pub fn areal(self) -> f64 {
        AREA_FACTOR * self.0
    }

    /// The corresponding volumetric coefficient `~ 3 * alpha`, in 1/K.
    pub fn volumetric(self) -> f64 {
        VOLUME_FACTOR * self.0
    }
}

/// The change in length `dL = alpha * L0 * dT`.
///
/// The result has the same units as `length`. A positive `delta_t`
/// (heating) with a positive `alpha` gives a positive `dL` (expansion); a
/// negative `delta_t` (cooling) gives a contraction; a zero `delta_t`
/// gives exactly `0.0`.
///
/// # Errors
///
/// Returns [`ThermalError::NonPositive`] if `length` is not strictly
/// positive, and [`ThermalError::NonFinite`] if `length` or `delta_t` is
/// not finite. `delta_t` itself may be any finite sign.
pub fn linear_expansion(
    alpha: LinearCoefficient,
    length: f64,
    delta_t: f64,
) -> Result<f64, ThermalError> {
    let length = require_positive("length", length)?;
    let delta_t = require_finite("delta_t", delta_t)?;
    Ok(alpha.per_kelvin() * length * delta_t)
}

/// The new length after heating / cooling, `L = L0 * (1 + alpha * dT)`.
///
/// This is simply `length + linear_expansion(..)`; it is provided so
/// callers do not have to remember to add the delta back to the original
/// dimension.
///
/// # Errors
///
/// Same conditions as [`linear_expansion`].
pub fn linear_final_length(
    alpha: LinearCoefficient,
    length: f64,
    delta_t: f64,
) -> Result<f64, ThermalError> {
    Ok(length + linear_expansion(alpha, length, delta_t)?)
}

/// The change in area `dA = (2 * alpha) * A0 * dT`, to first order.
///
/// Uses the areal coefficient `2 * alpha`. The result has the same units
/// as `area`.
///
/// # Errors
///
/// Returns [`ThermalError::NonPositive`] if `area` is not strictly
/// positive, and [`ThermalError::NonFinite`] if `area` or `delta_t` is not
/// finite.
pub fn area_expansion(
    alpha: LinearCoefficient,
    area: f64,
    delta_t: f64,
) -> Result<f64, ThermalError> {
    let area = require_positive("area", area)?;
    let delta_t = require_finite("delta_t", delta_t)?;
    Ok(alpha.areal() * area * delta_t)
}

/// The new area after heating / cooling, `A = A0 * (1 + 2 * alpha * dT)`,
/// to first order.
///
/// # Errors
///
/// Same conditions as [`area_expansion`].
pub fn area_final(alpha: LinearCoefficient, area: f64, delta_t: f64) -> Result<f64, ThermalError> {
    Ok(area + area_expansion(alpha, area, delta_t)?)
}

/// The change in volume `dV = (3 * alpha) * V0 * dT`, to first order.
///
/// Uses the volumetric coefficient `3 * alpha`. The result has the same
/// units as `volume`.
///
/// # Errors
///
/// Returns [`ThermalError::NonPositive`] if `volume` is not strictly
/// positive, and [`ThermalError::NonFinite`] if `volume` or `delta_t` is
/// not finite.
pub fn volume_expansion(
    alpha: LinearCoefficient,
    volume: f64,
    delta_t: f64,
) -> Result<f64, ThermalError> {
    let volume = require_positive("volume", volume)?;
    let delta_t = require_finite("delta_t", delta_t)?;
    Ok(alpha.volumetric() * volume * delta_t)
}

/// The new volume after heating / cooling,
/// `V = V0 * (1 + 3 * alpha * dT)`, to first order.
///
/// # Errors
///
/// Same conditions as [`volume_expansion`].
pub fn volume_final(
    alpha: LinearCoefficient,
    volume: f64,
    delta_t: f64,
) -> Result<f64, ThermalError> {
    Ok(volume + volume_expansion(alpha, volume, delta_t)?)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Absolute tolerance for floating-point comparisons. The numbers here
    /// are O(1)-to-O(1e3); 1e-9 is far tighter than the model error yet
    /// loose enough to absorb IEEE-754 rounding.
    const EPS: f64 = 1e-9;

    fn alpha_aluminium() -> LinearCoefficient {
        // Representative aluminium CTE, 23.1e-6 / K.
        LinearCoefficient::new(23.1e-6).unwrap()
    }

    #[test]
    fn linear_expansion_matches_alpha_l_dt() {
        let alpha = alpha_aluminium();
        let l0 = 2.0; // metres
        let dt = 100.0; // kelvin
        let got = linear_expansion(alpha, l0, dt).unwrap();
        let expected = 23.1e-6 * l0 * dt;
        assert!(
            (got - expected).abs() < EPS,
            "dL: got {got}, expected {expected}"
        );
        // Hand-checked magnitude: 23.1e-6 * 2 * 100 = 4.62e-3 m.
        assert!((got - 4.62e-3).abs() < EPS, "dL magnitude wrong: {got}");
    }

    #[test]
    fn final_length_is_l0_plus_delta() {
        let alpha = alpha_aluminium();
        let l0 = 2.0;
        let dt = 100.0;
        let lf = linear_final_length(alpha, l0, dt).unwrap();
        let expected = l0 * (1.0 + 23.1e-6 * dt);
        assert!(
            (lf - expected).abs() < EPS,
            "Lf: got {lf}, expected {expected}"
        );
    }

    #[test]
    fn area_coefficient_is_two_alpha() {
        let alpha = alpha_aluminium();
        // The areal coefficient must be exactly 2 * alpha to first order.
        assert!(
            (alpha.areal() - 2.0 * alpha.per_kelvin()).abs() < EPS,
            "areal coeff != 2 alpha"
        );
        // And dA = 2 alpha A dT must use it.
        let a0 = 5.0;
        let dt = 50.0;
        let got = area_expansion(alpha, a0, dt).unwrap();
        let expected = 2.0 * 23.1e-6 * a0 * dt;
        assert!(
            (got - expected).abs() < EPS,
            "dA: got {got}, expected {expected}"
        );
    }

    #[test]
    fn volume_coefficient_is_three_alpha() {
        let alpha = alpha_aluminium();
        // The volumetric coefficient must be exactly 3 * alpha.
        assert!(
            (alpha.volumetric() - 3.0 * alpha.per_kelvin()).abs() < EPS,
            "volumetric coeff != 3 alpha"
        );
        let v0 = 1.5;
        let dt = 80.0;
        let got = volume_expansion(alpha, v0, dt).unwrap();
        let expected = 3.0 * 23.1e-6 * v0 * dt;
        assert!(
            (got - expected).abs() < EPS,
            "dV: got {got}, expected {expected}"
        );
    }

    #[test]
    fn coefficient_ratios_are_one_two_three() {
        let alpha = alpha_aluminium();
        // areal / linear == 2, volumetric / linear == 3, exactly.
        assert!((alpha.areal() / alpha.per_kelvin() - 2.0).abs() < EPS);
        assert!((alpha.volumetric() / alpha.per_kelvin() - 3.0).abs() < EPS);
    }

    #[test]
    fn positive_alpha_expands_on_heating() {
        let alpha = alpha_aluminium();
        let dl = linear_expansion(alpha, 1.0, 50.0).unwrap();
        let da = area_expansion(alpha, 1.0, 50.0).unwrap();
        let dv = volume_expansion(alpha, 1.0, 50.0).unwrap();
        assert!(dl > 0.0, "length should grow: {dl}");
        assert!(da > 0.0, "area should grow: {da}");
        assert!(dv > 0.0, "volume should grow: {dv}");
    }

    #[test]
    fn cooling_contracts() {
        let alpha = alpha_aluminium();
        let dl = linear_expansion(alpha, 1.0, -50.0).unwrap();
        assert!(dl < 0.0, "cooling should contract: {dl}");
        // And it must be the exact negative of the heating case.
        let dl_heat = linear_expansion(alpha, 1.0, 50.0).unwrap();
        assert!((dl + dl_heat).abs() < EPS, "asymmetric heating/cooling");
    }

    #[test]
    fn zero_delta_t_means_no_change() {
        let alpha = alpha_aluminium();
        assert!(linear_expansion(alpha, 3.3, 0.0).unwrap().abs() < EPS);
        assert!(area_expansion(alpha, 3.3, 0.0).unwrap().abs() < EPS);
        assert!(volume_expansion(alpha, 3.3, 0.0).unwrap().abs() < EPS);
        // Final dimensions equal the originals when dT = 0.
        assert!((linear_final_length(alpha, 3.3, 0.0).unwrap() - 3.3).abs() < EPS);
        assert!((area_final(alpha, 3.3, 0.0).unwrap() - 3.3).abs() < EPS);
        assert!((volume_final(alpha, 3.3, 0.0).unwrap() - 3.3).abs() < EPS);
    }

    #[test]
    fn rejects_bad_inputs() {
        let alpha = alpha_aluminium();
        assert!(matches!(
            linear_expansion(alpha, 0.0, 10.0),
            Err(ThermalError::NonPositive { name: "length", .. })
        ));
        assert!(matches!(
            area_expansion(alpha, -1.0, 10.0),
            Err(ThermalError::NonPositive { name: "area", .. })
        ));
        assert!(matches!(
            volume_expansion(alpha, 1.0, f64::NAN),
            Err(ThermalError::NonFinite {
                name: "delta_t",
                ..
            })
        ));
        assert!(matches!(
            LinearCoefficient::new(-1.0),
            Err(ThermalError::NonPositive { name: "alpha", .. })
        ));
    }
}
