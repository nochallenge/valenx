//! AGMA bending-stress basics (AGMA 2001-D04 / Shigley, SI form).
//!
//! ## Model
//!
//! The American Gear Manufacturers Association refines the bare Lewis
//! estimate by multiplying in a family of correction factors. In SI
//! (metric-module) form the AGMA bending stress is
//!
//! ```text
//! sigma = Wt * Ko * Kv * Ks * (1 / (b * m_t)) * (Kh * Kb / Y_j)
//! ```
//!
//! where (Shigley notation, equation 14-15):
//!
//! - `Wt` transmitted tangential load, N,
//! - `Ko` overload factor (driven/driving machine roughness),
//! - `Kv` dynamic (velocity) factor (>= 1; impact from tooth errors),
//! - `Ks` size factor (>= 1),
//! - `b`  net face width, mm,
//! - `m_t` transverse module, mm,
//! - `Kh` load-distribution factor (face misalignment),
//! - `Kb` rim-thickness factor,
//! - `Y_j` (or `J`) the AGMA bending-strength geometry factor.
//!
//! With load in newtons and lengths in millimetres the stress is in
//! MPa, matching [`crate::lewis`].
//!
//! The geometry factor `J = Y / Kf` folds the Lewis form factor `Y`
//! together with a fillet stress-concentration factor `Kf` and a
//! load-sharing ratio; this module exposes the simple definition
//! `J = Y / Kf` plus the standard `Kv` curve-fit by transmission
//! quality number `Qv`.
//!
//! ## Honest scope
//!
//! These are the textbook correction-factor *forms*. The actual factor
//! values for a real gearset come from AGMA charts, the application,
//! and measured mounting tolerances. The numbers here are illustrative
//! defaults, not a substitute for a full AGMA 2001 rating.

use crate::error::{require_in_range, require_positive, GearToothError};
use serde::{Deserialize, Serialize};

/// AGMA transmission accuracy-level (quality) number `Qv`.
///
/// Higher numbers mean more precise teeth and therefore a dynamic
/// factor closer to unity. The standard curve-fit covers `Qv` in
/// `6..=11`; commercial gears sit around `6..=8`, precision gears
/// `10..=12`.
pub const QV_MIN: f64 = 6.0;

/// Upper bound of the `Qv` curve-fit domain.
pub const QV_MAX: f64 = 11.0;

/// AGMA dynamic (velocity) factor `Kv` for the given transmission
/// quality number `Qv` and pitch-line velocity `V` (m/s).
///
/// Shigley equations 14-27..14-29 give the curve-fit
///
/// ```text
/// B  = 0.25 * (12 - Qv)^(2/3)
/// A  = 50 + 56 * (1 - B)
/// Kv = ((A + sqrt(200 V)) / A)^B
/// ```
///
/// with `V` in metres per second. `Kv` is always at least 1; it grows
/// with speed and shrinks toward 1 as quality improves.
///
/// # Errors
///
/// Returns [`GearToothError::BadParameter`] if `Qv` is outside
/// `[QV_MIN, QV_MAX]` or if `V` is non-finite or non-positive.
pub fn dynamic_factor_kv(qv: f64, velocity_m_per_s: f64) -> Result<f64, GearToothError> {
    let qv = require_in_range("qv", qv, QV_MIN, QV_MAX)?;
    let v = require_positive("velocity_m_per_s", velocity_m_per_s)?;
    let b = 0.25 * (12.0 - qv).powf(2.0 / 3.0);
    let a = 50.0 + 56.0 * (1.0 - b);
    Ok(((a + (200.0 * v).sqrt()) / a).powf(b))
}

/// AGMA bending geometry factor `J` from the Lewis form factor `Y` and
/// a fillet stress-concentration factor `Kf`.
///
/// Definition: `J = Y / Kf`. The form factor `Y` is the dimensionless
/// Lewis profile factor (see [`crate::lewis_factor`]); `Kf >= 1`
/// accounts for the stress riser at the root fillet, so `J <= Y`.
///
/// # Errors
///
/// Returns [`GearToothError::BadParameter`] if `Y` is non-positive or
/// if `Kf` is below 1 (a concentration factor cannot relieve stress).
pub fn geometry_factor_j(form_factor_y: f64, fillet_kf: f64) -> Result<f64, GearToothError> {
    let y = require_positive("form_factor_y", form_factor_y)?;
    let kf = require_in_range("fillet_kf", fillet_kf, 1.0, f64::INFINITY)?;
    Ok(y / kf)
}

/// The full set of AGMA bending correction factors.
///
/// Each is dimensionless. Construct with [`AgmaFactors::new`], which
/// enforces the physical lower bounds (`Kv`, `Ks`, `Kh`, `Kb`, `Ko`
/// are all at least the values a benign gearset would see).
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct AgmaFactors {
    /// Overload factor `Ko` (>= 1).
    pub overload_ko: f64,
    /// Dynamic / velocity factor `Kv` (>= 1).
    pub dynamic_kv: f64,
    /// Size factor `Ks` (>= 1).
    pub size_ks: f64,
    /// Load-distribution factor `Kh` (>= 1).
    pub load_dist_kh: f64,
    /// Rim-thickness factor `Kb` (>= 1; 1 for a solid-blank gear).
    pub rim_kb: f64,
}

impl AgmaFactors {
    /// Build a validated [`AgmaFactors`].
    ///
    /// # Errors
    ///
    /// Returns [`GearToothError::BadParameter`] if any factor is
    /// non-finite or less than 1.
    pub fn new(
        overload_ko: f64,
        dynamic_kv: f64,
        size_ks: f64,
        load_dist_kh: f64,
        rim_kb: f64,
    ) -> Result<Self, GearToothError> {
        let overload_ko = require_in_range("overload_ko", overload_ko, 1.0, f64::INFINITY)?;
        let dynamic_kv = require_in_range("dynamic_kv", dynamic_kv, 1.0, f64::INFINITY)?;
        let size_ks = require_in_range("size_ks", size_ks, 1.0, f64::INFINITY)?;
        let load_dist_kh = require_in_range("load_dist_kh", load_dist_kh, 1.0, f64::INFINITY)?;
        let rim_kb = require_in_range("rim_kb", rim_kb, 1.0, f64::INFINITY)?;
        Ok(Self {
            overload_ko,
            dynamic_kv,
            size_ks,
            load_dist_kh,
            rim_kb,
        })
    }

    /// A benign baseline: every correction factor equal to 1, i.e. the
    /// AGMA equation collapses to `sigma = Wt / (b m J)` — the Lewis
    /// equation with `Y` replaced by the geometry factor `J`.
    pub fn unity() -> Self {
        Self {
            overload_ko: 1.0,
            dynamic_kv: 1.0,
            size_ks: 1.0,
            load_dist_kh: 1.0,
            rim_kb: 1.0,
        }
    }
}

/// AGMA root bending stress `sigma`, in MPa.
///
/// Evaluates
///
/// ```text
/// sigma = Wt * Ko * Kv * Ks * (1 / (b m)) * (Kh * Kb / J)
/// ```
///
/// with `Wt` in newtons, face width `b` and module `m` in millimetres,
/// the geometry factor `J` dimensionless, and the correction
/// [`AgmaFactors`] dimensionless.
///
/// When every correction factor is 1 (see [`AgmaFactors::unity`]) this
/// reduces to `Wt / (b m J)`, the Lewis equation in which the form
/// factor `Y` has been replaced by the AGMA geometry factor `J`.
///
/// # Errors
///
/// Returns [`GearToothError::BadParameter`] if `Wt`, `b`, `m`, or `J`
/// is non-finite or non-positive. The factors are validated at
/// [`AgmaFactors`] construction time.
pub fn agma_bending_stress(
    tangential_load_n: f64,
    face_width_mm: f64,
    module_mm: f64,
    geometry_factor_j: f64,
    factors: &AgmaFactors,
) -> Result<f64, GearToothError> {
    let wt = require_positive("tangential_load_n", tangential_load_n)?;
    let b = require_positive("face_width_mm", face_width_mm)?;
    let m = require_positive("module_mm", module_mm)?;
    let j = require_positive("geometry_factor_j", geometry_factor_j)?;
    let sigma = wt * factors.overload_ko * factors.dynamic_kv * factors.size_ks / (b * m)
        * (factors.load_dist_kh * factors.rim_kb / j);
    Ok(sigma)
}
