//! The Lewis bending-stress equation and its supporting kinematics.
//!
//! ## Model
//!
//! Wilfred Lewis (1892) modelled a gear tooth as a cantilever beam of
//! uniform-strength (parabolic) outline loaded by the transmitted
//! tangential force `Wt`. The result, in SI metric-module form, is the
//! root bending stress
//!
//! ```text
//! sigma = Wt / (F * m * Y)
//! ```
//!
//! where
//!
//! - `Wt` is the tangential (transmitted) load at the pitch line, in N,
//! - `F` is the face width, in mm,
//! - `m` is the module, in mm,
//! - `Y` is the dimensionless [Lewis form factor](crate::lewis_factor).
//!
//! With `Wt` in newtons and `F`, `m` in millimetres the stress comes
//! out in newtons per square millimetre, i.e. **megapascals (MPa)** —
//! the convention used throughout this crate.
//!
//! ## Honest scope
//!
//! The bare Lewis equation ignores stress concentration at the fillet,
//! dynamic (impact) effects, and load sharing between teeth. It is a
//! first-pass sizing estimate, not a rating. The AGMA refinement in
//! [`crate::agma`] layers the standard correction factors on top.

use crate::error::{require_positive, GearToothError};
use crate::lewis_factor::lewis_form_factor;
use crate::spec::ToothLoad;
use serde::{Deserialize, Serialize};

/// Result of a Lewis bending-stress evaluation.
///
/// All stresses are in megapascals (MPa); the form factor is
/// dimensionless.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct LewisResult {
    /// Lewis form factor `Y` used (dimensionless).
    pub form_factor_y: f64,
    /// Root bending stress `sigma`, in MPa.
    pub bending_stress_mpa: f64,
}

/// Compute the Lewis root bending stress `sigma = Wt / (F m Y)`.
///
/// `tangential_load_n` is `Wt` in newtons, `face_width_mm` is `F`,
/// `module_mm` is `m`, and `form_factor_y` is the dimensionless Lewis
/// factor `Y` (look it up with
/// [`crate::lewis_factor::lewis_form_factor`]). The returned stress is
/// in MPa.
///
/// This is the low-level form that takes `Y` directly; see
/// [`lewis_bending_stress_for_teeth`] for the variant that derives `Y`
/// from the tooth count.
///
/// # Errors
///
/// Returns [`GearToothError::BadParameter`] if any of `Wt`, `F`, `m`,
/// or `Y` is non-finite or non-positive.
pub fn lewis_bending_stress(
    tangential_load_n: f64,
    face_width_mm: f64,
    module_mm: f64,
    form_factor_y: f64,
) -> Result<f64, GearToothError> {
    let wt = require_positive("tangential_load_n", tangential_load_n)?;
    let f = require_positive("face_width_mm", face_width_mm)?;
    let m = require_positive("module_mm", module_mm)?;
    let y = require_positive("form_factor_y", form_factor_y)?;
    Ok(wt / (f * m * y))
}

/// Compute the Lewis root bending stress for a tooth count, looking up
/// the Lewis form factor `Y` from the embedded 20-degree full-depth
/// table.
///
/// Convenience wrapper around [`lewis_bending_stress`] that first
/// resolves `Y(N)` via [`crate::lewis_factor::lewis_form_factor`], then
/// applies the Lewis equation. The returned [`LewisResult`] reports both
/// the `Y` that was used and the resulting stress so callers can show
/// their work.
///
/// # Errors
///
/// Returns [`GearToothError::OutOfDomain`] if `teeth` is below the
/// minimum tabulated count, or [`GearToothError::BadParameter`] if any
/// load / geometry input is non-finite or non-positive.
pub fn lewis_bending_stress_for_teeth(
    tangential_load_n: f64,
    face_width_mm: f64,
    module_mm: f64,
    teeth: u32,
) -> Result<LewisResult, GearToothError> {
    let y = lewis_form_factor(teeth)?;
    let sigma = lewis_bending_stress(tangential_load_n, face_width_mm, module_mm, y)?;
    Ok(LewisResult {
        form_factor_y: y,
        bending_stress_mpa: sigma,
    })
}

/// Compute the Lewis root bending stress from a [`ToothLoad`] bundle.
///
/// Equivalent to [`lewis_bending_stress_for_teeth`] but reads its
/// geometry and load from a validated [`ToothLoad`] value, which is the
/// ergonomic entry point for callers that already hold one.
///
/// # Errors
///
/// Propagates the same errors as [`lewis_bending_stress_for_teeth`].
pub fn lewis_bending_stress_of(load: &ToothLoad) -> Result<LewisResult, GearToothError> {
    lewis_bending_stress_for_teeth(
        load.tangential_load_n,
        load.face_width_mm,
        load.module_mm,
        load.teeth,
    )
}

/// Pitch-line (pitch-circle tangential) velocity, in metres per second.
///
/// Given the pitch diameter `d` (mm) and the rotational speed `n`
/// (rev/min), the pitch-line velocity is
///
/// ```text
/// V = pi * d * n / 60000
/// ```
///
/// The `60000` folds together the per-minute-to-per-second factor (60)
/// and the millimetre-to-metre factor (1000).
///
/// # Errors
///
/// Returns [`GearToothError::BadParameter`] if either input is
/// non-finite or non-positive.
pub fn pitch_line_velocity_m_per_s(
    pitch_diameter_mm: f64,
    speed_rpm: f64,
) -> Result<f64, GearToothError> {
    let d = require_positive("pitch_diameter_mm", pitch_diameter_mm)?;
    let n = require_positive("speed_rpm", speed_rpm)?;
    Ok(std::f64::consts::PI * d * n / 60_000.0)
}

/// Transmitted tangential load `Wt`, in newtons, from transmitted power
/// and pitch-line velocity.
///
/// Mechanical power equals force times velocity, so
///
/// ```text
/// Wt = P / V
/// ```
///
/// with `P` in watts and `V` in metres per second, giving `Wt` in
/// newtons.
///
/// # Errors
///
/// Returns [`GearToothError::BadParameter`] if either input is
/// non-finite or non-positive.
pub fn tangential_load_from_power_n(
    power_w: f64,
    pitch_line_velocity_m_per_s: f64,
) -> Result<f64, GearToothError> {
    let p = require_positive("power_w", power_w)?;
    let v = require_positive("pitch_line_velocity_m_per_s", pitch_line_velocity_m_per_s)?;
    Ok(p / v)
}

/// Transmitted tangential load `Wt`, in newtons, from transmitted
/// torque and pitch radius.
///
/// Torque is force times radius, so `Wt = T / r`. Here `T` is in
/// newton-metres and the pitch *diameter* is supplied in millimetres,
/// so the radius in metres is `d / 2000`:
///
/// ```text
/// Wt = T / (d / 2000) = 2000 * T / d
/// ```
///
/// # Errors
///
/// Returns [`GearToothError::BadParameter`] if either input is
/// non-finite or non-positive.
pub fn tangential_load_from_torque_n(
    torque_nm: f64,
    pitch_diameter_mm: f64,
) -> Result<f64, GearToothError> {
    let t = require_positive("torque_nm", torque_nm)?;
    let d = require_positive("pitch_diameter_mm", pitch_diameter_mm)?;
    Ok(2000.0 * t / d)
}
