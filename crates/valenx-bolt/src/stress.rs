//! Bolt cross-section and stress relations.
//!
//! A threaded bolt does not carry tension on its full nominal-diameter
//! circle; it carries it on the smaller **tensile-stress area** `A_t`,
//! an effective area at the mean of the pitch and minor (root) diameters.
//! For ISO metric threads the standard fit is
//!
//! ```text
//! A_t = (pi / 4) * (d - 0.938_194 * P)^2
//! ```
//!
//! where `d` is the nominal diameter and `P` the thread pitch (both in
//! the same length unit). Axial stress is then simply force over `A_t`,
//! and the proof / tensile *loads* are the material strengths times
//! `A_t`.
//!
//! All functions are unit-agnostic as long as the inputs are
//! consistent: pass `d` and `P` in metres and a strength in pascals to
//! get an area in m² and a load in newtons.

use crate::error::BoltError;
use crate::material::BoltMaterial;

/// The ISO-metric coefficient in `A_t = (pi/4)(d - C*P)^2`.
///
/// `C = (3/2) * sqrt(3) / something` historically; the standard
/// tabulated constant is `0.9382` (the mean of the basic pitch and minor
/// diameter offsets), here carried to six places.
pub const ISO_TENSILE_AREA_PITCH_COEFF: f64 = 0.938_194;

/// Tensile-stress area `A_t` of an ISO metric thread.
///
/// `nominal_diameter` and `pitch` must share a length unit; the result
/// is in that unit squared.
///
/// # Errors
///
/// Returns [`BoltError`] if either input is not strictly positive and
/// finite, or if the pitch is so coarse relative to the diameter that
/// the effective diameter `d - C*P` is non-positive (a degenerate
/// thread).
pub fn tensile_stress_area(nominal_diameter: f64, pitch: f64) -> Result<f64, BoltError> {
    let d = BoltError::require_positive("nominal_diameter", nominal_diameter)?;
    let p = BoltError::require_positive("pitch", pitch)?;
    let effective = d - ISO_TENSILE_AREA_PITCH_COEFF * p;
    // Guard the degenerate "pitch larger than diameter" case.
    let effective = BoltError::require_positive("effective_diameter", effective)?;
    Ok(std::f64::consts::FRAC_PI_4 * effective * effective)
}

/// Axial tensile stress in a bolt carrying force `force` over a
/// tensile-stress area `area`.
///
/// `stress = force / area`. Units follow the inputs (N over m² gives
/// Pa).
///
/// # Errors
///
/// Returns [`BoltError`] if `area` is not strictly positive and finite,
/// or if `force` is not finite.
pub fn axial_stress(force: f64, area: f64) -> Result<f64, BoltError> {
    let a = BoltError::require_positive("area", area)?;
    if !force.is_finite() {
        return Err(BoltError::NotFinite {
            name: "force",
            value: force,
        });
    }
    Ok(force / a)
}

/// The **proof load** of a bolt — the largest tensile force it can carry
/// with no permanent set: `F_p = S_p * A_t`.
///
/// # Errors
///
/// Returns [`BoltError`] if `area` is not strictly positive and finite.
pub fn proof_load(material: &BoltMaterial, area: f64) -> Result<f64, BoltError> {
    let a = BoltError::require_positive("area", area)?;
    Ok(material.proof_strength_pa * a)
}

/// The **tensile (proof-to-failure) load** of a bolt:
/// `F_u = S_u * A_t`.
///
/// # Errors
///
/// Returns [`BoltError`] if `area` is not strictly positive and finite.
pub fn tensile_load(material: &BoltMaterial, area: f64) -> Result<f64, BoltError> {
    let a = BoltError::require_positive("area", area)?;
    Ok(material.tensile_strength_pa * a)
}

/// The recommended preload force for a reused joint, `0.75 S_p A_t`.
///
/// # Errors
///
/// Returns [`BoltError`] if `area` is not strictly positive and finite.
pub fn recommended_preload(material: &BoltMaterial, area: f64) -> Result<f64, BoltError> {
    let a = BoltError::require_positive("area", area)?;
    Ok(material.recommended_preload_stress_pa() * a)
}
