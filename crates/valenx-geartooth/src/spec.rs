//! Validated input bundles for gear-tooth bending calculations.
//!
//! [`ToothLoad`] groups the four quantities the Lewis equation needs —
//! tangential load, face width, module, and tooth count — behind a
//! checked constructor so downstream calculators never see a
//! non-positive or non-finite value.

use crate::error::{require_positive, GearToothError};
use serde::{Deserialize, Serialize};

/// A validated bundle of the inputs to the Lewis bending equation.
///
/// Construct one with [`ToothLoad::new`], which rejects non-finite or
/// non-positive values. The fields are public for read access and
/// serialization, but the constructor is the only checked way to build
/// a value with confidence in its invariants.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ToothLoad {
    /// Tangential (transmitted) load at the pitch line `Wt`, in newtons.
    pub tangential_load_n: f64,
    /// Face width `F`, in millimetres.
    pub face_width_mm: f64,
    /// Module `m`, in millimetres. Pitch diameter = `module * teeth`.
    pub module_mm: f64,
    /// Number of teeth `N` on the gear (drives the Lewis form factor).
    pub teeth: u32,
}

impl ToothLoad {
    /// Build a validated [`ToothLoad`].
    ///
    /// # Errors
    ///
    /// Returns [`GearToothError::BadParameter`] if any of the floating
    /// inputs is non-finite or non-positive, or if `teeth` is zero.
    pub fn new(
        tangential_load_n: f64,
        face_width_mm: f64,
        module_mm: f64,
        teeth: u32,
    ) -> Result<Self, GearToothError> {
        let tangential_load_n = require_positive("tangential_load_n", tangential_load_n)?;
        let face_width_mm = require_positive("face_width_mm", face_width_mm)?;
        let module_mm = require_positive("module_mm", module_mm)?;
        if teeth == 0 {
            return Err(GearToothError::bad_parameter("teeth", "must be at least 1"));
        }
        Ok(Self {
            tangential_load_n,
            face_width_mm,
            module_mm,
            teeth,
        })
    }

    /// Pitch (reference) circle diameter, in millimetres: `module * teeth`.
    pub fn pitch_diameter_mm(&self) -> f64 {
        self.module_mm * self.teeth as f64
    }
}
