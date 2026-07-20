//! Lame stress-field evaluators and the thin-wall comparison.
//!
//! Given a [`Cylinder`] and the solved [`LameConstants`], the radial and
//! hoop (circumferential) stresses through the wall are
//!
//! ```text
//! sigma_r(r)     = A - B / r^2
//! sigma_theta(r) = A + B / r^2
//! ```
//!
//! Because `B >= 0` whenever the net pressure `p_i - p_o >= 0`, the hoop
//! stress is largest at the bore (`r = a`) and the radial stress equals
//! `-p_i` there, recovering the bore boundary condition exactly.

use crate::cylinder::{Cylinder, LameConstants};
use crate::error::ThickCylinderError;

/// Radial stress `sigma_r(r) = A - B/r^2` at radius `r`.
///
/// At `r = a` this returns `-p_i` and at `r = b` it returns `-p_o`,
/// reproducing the pressure boundary conditions (compression is negative
/// under the tension-positive convention).
///
/// # Errors
///
/// [`ThickCylinderError::OutsideWall`] if `r` is not in `[a, b]`.
pub fn radial_stress(cyl: &Cylinder, k: &LameConstants, r: f64) -> Result<f64, ThickCylinderError> {
    let r = cyl.check_radius(r)?;
    Ok(k.a - k.b / (r * r))
}

/// Hoop (circumferential) stress `sigma_theta(r) = A + B/r^2` at radius `r`.
///
/// # Errors
///
/// [`ThickCylinderError::OutsideWall`] if `r` is not in `[a, b]`.
pub fn hoop_stress(cyl: &Cylinder, k: &LameConstants, r: f64) -> Result<f64, ThickCylinderError> {
    let r = cyl.check_radius(r)?;
    Ok(k.a + k.b / (r * r))
}

/// Maximum hoop stress, which for a pressurised cylinder occurs at the
/// bore `r = a`.
///
/// For internal pressure only (`p_o = 0`) this equals the classic
/// `sigma_theta,max = p_i (b^2 + a^2) / (b^2 - a^2)`.
///
/// # Errors
///
/// Propagates evaluation errors (none are expected for `r = a`, which is
/// always inside the wall).
pub fn max_hoop_stress(cyl: &Cylinder, k: &LameConstants) -> Result<f64, ThickCylinderError> {
    hoop_stress(cyl, k, cyl.inner())
}

/// Thin-wall ("membrane") hoop-stress estimate `sigma = p r_m / t`.
///
/// This is the elementary pressure-vessel formula obtained by ignoring the
/// radial stress gradient. It uses the mean radius `r_m = (a + b)/2` and
/// the wall thickness `t = b - a`. The function exists so the exact Lame
/// result can be validated against the limiting case: as `t/r -> 0` the
/// bore hoop stress converges to this value.
///
/// Only the *net* pressure `p_i - p_o` drives membrane hoop stress, matching
/// the thin-shell derivation.
///
/// # Errors
///
/// [`ThickCylinderError::NonPositive`] if a pressure is negative and
/// [`ThickCylinderError::NotFinite`] for NaN/Inf.
pub fn thin_wall_hoop_stress(
    cyl: &Cylinder,
    p_internal: f64,
    p_external: f64,
) -> Result<f64, ThickCylinderError> {
    let p_i = crate::error::non_negative("p_internal", p_internal)?;
    let p_o = crate::error::non_negative("p_external", p_external)?;
    let net = p_i - p_o;
    Ok(net * cyl.mean_radius() / cyl.thickness())
}
