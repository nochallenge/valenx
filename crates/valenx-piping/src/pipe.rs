//! Pipe section + CAD solid emit.

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

use valenx_cad::primitives::cylinder;
use valenx_cad::Solid;

use crate::dims::{nominal_to_od_in, Schedule};
use crate::error::PipingError;

/// Material the pipe is fabricated from. Affects roughness +
/// downstream BOM but not v1 geometry.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Material {
    /// ASTM A53 / A106 carbon steel.
    CarbonSteel,
    /// 304 / 316 stainless steel.
    StainlessSteel,
    /// Copper Type L / M.
    Copper,
    /// PVC (Schedule 40 / 80).
    Pvc,
    /// Cross-linked polyethylene.
    Pex,
}

/// A straight pipe section between two world-space points.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct PipeSection {
    /// Start point in world coordinates (mm).
    pub start: Vector3<f64>,
    /// End point in world coordinates (mm).
    pub end: Vector3<f64>,
    /// Nominal pipe size designation (e.g. `"2"` for NPS 2).
    pub nominal_size: String,
    /// Wall schedule.
    pub schedule: Schedule,
    /// Material.
    pub material: Material,
}

impl PipeSection {
    /// Construct a new section.
    pub fn new(
        start: Vector3<f64>,
        end: Vector3<f64>,
        nominal_size: impl Into<String>,
        schedule: Schedule,
        material: Material,
    ) -> Self {
        Self {
            start,
            end,
            nominal_size: nominal_size.into(),
            schedule,
            material,
        }
    }

    /// Outside diameter in millimetres.
    pub fn outer_diameter_mm(&self) -> Result<f64, PipingError> {
        nominal_to_od_in(&self.nominal_size)
            .map(|inches| inches * 25.4)
            .ok_or_else(|| PipingError::UnknownNps(self.nominal_size.clone()))
    }

    /// Section length in millimetres.
    pub fn length_mm(&self) -> f64 {
        (self.end - self.start).norm()
    }
}

/// Convert a [`PipeSection`] to a [`Solid`].
///
/// v1: emits a cylinder of length = section.length() and radius =
/// OD/2 oriented along +Z (callers translate / rotate as needed).
/// A future v2 will sweep the cross-section along the section axis
/// directly via `valenx_cad::sweep`.
pub fn to_solid(section: &PipeSection) -> Result<Solid, PipingError> {
    let od = section.outer_diameter_mm()?;
    let len = section.length_mm();
    if len <= 0.0 {
        return Err(PipingError::BadParameter {
            name: "length",
            reason: "section start == end".into(),
        });
    }
    cylinder(od / 2.0, len).map_err(|e| PipingError::Cad(e.to_string()))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn outer_diameter_matches_nps_table() {
        let s = PipeSection::new(
            Vector3::zeros(),
            Vector3::new(0.0, 0.0, 100.0),
            "2",
            Schedule::Sch40,
            Material::CarbonSteel,
        );
        let od = s.outer_diameter_mm().unwrap();
        // 2.375 in × 25.4 = 60.325 mm
        assert!((od - 60.325).abs() < 1e-6);
    }

    #[test]
    fn unknown_nps_surfaces_error() {
        let s = PipeSection::new(
            Vector3::zeros(),
            Vector3::new(1.0, 0.0, 0.0),
            "99",
            Schedule::Sch40,
            Material::Pvc,
        );
        let err = s.outer_diameter_mm().unwrap_err();
        assert!(matches!(err, PipingError::UnknownNps(_)));
    }

    #[test]
    fn to_solid_rejects_zero_length() {
        let s = PipeSection::new(
            Vector3::zeros(),
            Vector3::zeros(),
            "1",
            Schedule::Sch40,
            Material::Copper,
        );
        let err = to_solid(&s).unwrap_err();
        assert!(matches!(err, PipingError::BadParameter { .. }));
    }
}
