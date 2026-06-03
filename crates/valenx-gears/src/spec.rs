//! [`GearSpec`] + [`GearKind`].

use serde::{Deserialize, Serialize};

/// Type of gear to generate.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum GearKind {
    /// Parallel-axis straight-tooth (spur).
    Spur,
    /// Parallel-axis with helix angle.
    Helical,
    /// Intersecting-axis truncated-cone bevel.
    Bevel,
    /// Crossed-axis worm (helical thread + worm gear pinion).
    Worm,
}

impl GearKind {
    /// Short UI label.
    pub fn label(self) -> &'static str {
        match self {
            Self::Spur => "Spur",
            Self::Helical => "Helical",
            Self::Bevel => "Bevel",
            Self::Worm => "Worm",
        }
    }
}

/// Parametric description of a single gear.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct GearSpec {
    /// Which family.
    pub kind: GearKind,
    /// Module (mm). Pitch diameter = module × teeth.
    pub module_mm: f64,
    /// Tooth count.
    pub teeth: u32,
    /// Pressure angle, degrees. Standard = 20°.
    pub pressure_angle_deg: f64,
    /// Helix angle, degrees. 0 for spur. ~20-30 for helical.
    pub helix_angle_deg: f64,
    /// Face width, mm.
    pub face_width_mm: f64,
}

impl GearSpec {
    /// Convenience: a standard 1-module, 20° spur gear.
    pub fn standard_spur(teeth: u32) -> Self {
        Self {
            kind: GearKind::Spur,
            module_mm: 1.0,
            teeth,
            pressure_angle_deg: 20.0,
            helix_angle_deg: 0.0,
            face_width_mm: 10.0,
        }
    }

    /// Pitch (reference) circle diameter — `module × teeth`.
    pub fn pitch_diameter_mm(&self) -> f64 {
        self.module_mm * self.teeth as f64
    }

    /// Base circle diameter — `pitch × cos(pressure_angle)`. The
    /// involute curve is generated from this circle.
    pub fn base_diameter_mm(&self) -> f64 {
        self.pitch_diameter_mm() * self.pressure_angle_deg.to_radians().cos()
    }

    /// Addendum diameter — pitch + 2 × module (standard).
    pub fn addendum_diameter_mm(&self) -> f64 {
        self.pitch_diameter_mm() + 2.0 * self.module_mm
    }

    /// Dedendum diameter — pitch − 2.5 × module (clearance 0.25 m).
    pub fn dedendum_diameter_mm(&self) -> f64 {
        (self.pitch_diameter_mm() - 2.5 * self.module_mm).max(0.0)
    }
}

impl Default for GearSpec {
    fn default() -> Self {
        Self::standard_spur(20)
    }
}
