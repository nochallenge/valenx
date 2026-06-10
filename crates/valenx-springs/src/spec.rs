//! [`SpringSpec`] + [`SpringKind`] + end-treatment enum.

use serde::{Deserialize, Serialize};

/// Type of helical spring.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum SpringKind {
    /// Compression spring (loaded along the helix axis).
    Compression,
    /// Extension spring (with hook ends).
    Extension,
    /// Torsion spring (with radial leg ends).
    Torsion,
}

impl SpringKind {
    /// UI label.
    pub fn label(self) -> &'static str {
        match self {
            Self::Compression => "Compression",
            Self::Extension => "Extension",
            Self::Torsion => "Torsion",
        }
    }
}

/// End-coil treatment for compression springs.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum EndTreatment {
    /// Closed ends — last coil pressed flat.
    Closed,
    /// Closed and ground — last coil pressed flat and ground.
    ClosedGround,
    /// Open ends — no end treatment.
    Open,
    /// Open and ground.
    OpenGround,
}

impl EndTreatment {
    /// Inactive coils added by the treatment (per pair).
    pub fn inactive_coils(self) -> f64 {
        match self {
            Self::Closed | Self::ClosedGround => 2.0,
            Self::Open => 0.0,
            Self::OpenGround => 1.0,
        }
    }
}

/// Parametric helical spring description.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct SpringSpec {
    /// Spring kind.
    pub kind: SpringKind,
    /// Wire diameter (mm).
    pub wire_diameter_mm: f64,
    /// Mean coil diameter (mm) — measured to the centreline of the
    /// wire.
    pub mean_coil_diameter_mm: f64,
    /// Free length (mm) — unloaded overall length.
    pub free_length_mm: f64,
    /// Number of active coils.
    pub n_active_coils: f64,
    /// End treatment.
    pub end_treatment: EndTreatment,
    /// Shear modulus G (MPa). Steel default = 79_300 MPa.
    pub shear_modulus_mpa: f64,
}

impl SpringSpec {
    /// Convenience: a standard compression spring (steel, 1 mm wire,
    /// 10 mm coil, 30 mm free length, 8 active coils).
    pub fn default_compression() -> Self {
        Self {
            kind: SpringKind::Compression,
            wire_diameter_mm: 1.0,
            mean_coil_diameter_mm: 10.0,
            free_length_mm: 30.0,
            n_active_coils: 8.0,
            end_treatment: EndTreatment::Closed,
            shear_modulus_mpa: 79_300.0,
        }
    }

    /// Outer diameter — mean + wire.
    pub fn outer_diameter_mm(&self) -> f64 {
        self.mean_coil_diameter_mm + self.wire_diameter_mm
    }

    /// Inner diameter — mean − wire.
    pub fn inner_diameter_mm(&self) -> f64 {
        self.mean_coil_diameter_mm - self.wire_diameter_mm
    }

    /// Pitch (axial advance per coil) for the active section.
    pub fn pitch_mm(&self) -> f64 {
        if self.n_active_coils <= 0.0 {
            return self.free_length_mm;
        }
        self.free_length_mm / self.n_active_coils
    }

    /// Developed wire length of the active coils (mm) — the arc length of the
    /// helical centerline, `ℓ = n·√((π·D)² + p²)`, where `n` is the active-coil
    /// count, `D` the mean coil diameter, and `p` the per-coil [`pitch`](Self::pitch_mm)
    /// (one turn unrolls to a right triangle with base `π·D` and height `p`). A
    /// geometric length, distinct from the diameters and the mechanical scalars;
    /// the primary input to wire mass / cost. Returns `0.0` for a non-positive or
    /// non-finite active-coil count.
    pub fn helix_length_mm(&self) -> f64 {
        if !self.n_active_coils.is_finite() || self.n_active_coils <= 0.0 {
            return 0.0;
        }
        let circumference = std::f64::consts::PI * self.mean_coil_diameter_mm;
        let p = self.pitch_mm();
        self.n_active_coils * (circumference * circumference + p * p).sqrt()
    }
}

impl Default for SpringSpec {
    fn default() -> Self {
        Self::default_compression()
    }
}
