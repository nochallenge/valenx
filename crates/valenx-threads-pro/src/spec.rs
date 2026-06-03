//! [`ThreadSpecPro`] — a single thread row in any of the extended
//! tables.

use serde::{Deserialize, Serialize};

use crate::standard::{ProfileShape, ThreadStandardPro};

/// A single extended thread spec.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ThreadSpecPro {
    /// Family.
    pub standard: ThreadStandardPro,
    /// Display designation (`"M8"`, `"1/4-20 UNC"`, `"G 1/2"`,
    /// `"Tr 30x6"`, `"3/8-12 Acme"`).
    pub designation: String,
    /// Nominal diameter in millimetres.
    pub nominal_diameter: f64,
    /// Thread pitch in millimetres.
    pub pitch: f64,
    /// Thread profile shape — V, Acme, Trapezoidal, or Buttress.
    pub profile: ProfileShape,
}

impl ThreadSpecPro {
    /// Construct directly.
    pub fn new(
        standard: ThreadStandardPro,
        designation: impl Into<String>,
        nominal_diameter: f64,
        pitch: f64,
    ) -> Self {
        let profile = standard.profile_shape();
        Self {
            standard,
            designation: designation.into(),
            nominal_diameter,
            pitch,
            profile,
        }
    }

    /// Outer (major) diameter — equal to nominal for V profiles.
    pub fn major_diameter(&self) -> f64 {
        self.nominal_diameter
    }

    /// Minor diameter — `nominal - 2 * H`, where `H` is the canonical
    /// profile height factor times pitch (0.61343 for V, 0.5 for
    /// Acme / Trapezoidal, ~0.75 for Buttress' tall side).
    pub fn minor_diameter(&self) -> f64 {
        let h = match self.profile {
            ProfileShape::V => 0.61343,
            ProfileShape::Acme => 0.5,
            ProfileShape::Trapezoidal => 0.5,
            ProfileShape::Buttress => 0.75,
        };
        (self.nominal_diameter - 2.0 * h * self.pitch).max(0.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn major_equals_nominal() {
        let s = ThreadSpecPro::new(ThreadStandardPro::IsoMetric, "M8", 8.0, 1.25);
        assert!((s.major_diameter() - 8.0).abs() < 1e-9);
    }

    #[test]
    fn minor_smaller_than_major_for_v_thread() {
        let s = ThreadSpecPro::new(ThreadStandardPro::IsoMetric, "M8", 8.0, 1.25);
        assert!(s.minor_diameter() < s.major_diameter());
        assert!(s.minor_diameter() > 0.0);
    }

    #[test]
    fn acme_minor_uses_half_pitch_height() {
        let s = ThreadSpecPro::new(ThreadStandardPro::Acme, "5/8-8 Acme", 15.875, 3.175);
        // h = 0.5 * 3.175 → minor = 15.875 - 3.175 = 12.7
        assert!((s.minor_diameter() - 12.7).abs() < 1e-9);
    }
}
