//! Tool library: the cutter that physically removes stock.
//!
//! v1 ships six [`ToolKind`] flavors. A [`Tool`] bundles an `id`
//! (the integer the postprocessor emits in `T{id} M6`), a display
//! `name`, the kind, diameter, length, flute count, and a free-form
//! `material` string.
//!
//! Validation happens in [`Tool::new`] so any tool that makes it into
//! an operation already has a positive diameter / non-zero flute count.

use serde::{Deserialize, Serialize};

use crate::error::CamError;

/// Catalog of tool geometries the v1 CAM workbench knows how to plan
/// for.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, Serialize, Deserialize)]
pub enum ToolKind {
    /// Flat-end mill — the default tool for Profile, Pocket, Face.
    #[default]
    EndMill,
    /// Ball-end mill — used for 3D contouring; v1 plans 2.5D paths
    /// only but the kind is captured for postprocessor labelling.
    BallMill,
    /// Twist drill — the only kind valid for `Drill` operations.
    Drill,
    /// Face mill — large-diameter cutter for facing operations.
    FaceMill,
    /// Tap — threads pre-drilled holes; v1 captures the kind but
    /// does not yet plan tapping cycles.
    Tap,
    /// Reamer — finishes pre-drilled holes to tight tolerance.
    Reamer,
}

impl ToolKind {
    /// Short human-readable label for the panel + audit log.
    pub fn label(self) -> &'static str {
        match self {
            ToolKind::EndMill => "EndMill",
            ToolKind::BallMill => "BallMill",
            ToolKind::Drill => "Drill",
            ToolKind::FaceMill => "FaceMill",
            ToolKind::Tap => "Tap",
            ToolKind::Reamer => "Reamer",
        }
    }
}

impl std::fmt::Display for ToolKind {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.label())
    }
}

/// A physical cutter the post will reference by `id`.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Tool {
    /// Postprocessor tool number — emitted as `T{id}` in tool change
    /// commands. Must be unique within a [`crate::persist::CamFile`]
    /// (callers enforce uniqueness; the struct itself does not).
    pub id: u32,
    /// Human-readable name (panel label).
    pub name: String,
    /// Tool geometry / category.
    pub kind: ToolKind,
    /// Cutting diameter in millimetres.
    pub diameter_mm: f64,
    /// Overall flute / cutting length in millimetres.
    pub length_mm: f64,
    /// Flute count (1, 2, 3, 4, …).
    pub flutes: u32,
    /// Free-form material descriptor (e.g. `"HSS"`, `"carbide"`).
    pub material: String,
}

impl Tool {
    /// Construct a validated tool. Returns
    /// [`CamError::BadTool`] if `diameter_mm`/`length_mm` is not
    /// strictly positive or `flutes == 0`.
    pub fn new(
        id: u32,
        name: impl Into<String>,
        kind: ToolKind,
        diameter_mm: f64,
        length_mm: f64,
        flutes: u32,
        material: impl Into<String>,
    ) -> Result<Self, CamError> {
        if !(diameter_mm > 0.0) {
            return Err(CamError::BadTool {
                reason: format!("diameter must be > 0 (got {diameter_mm})"),
            });
        }
        if !(length_mm > 0.0) {
            return Err(CamError::BadTool {
                reason: format!("length must be > 0 (got {length_mm})"),
            });
        }
        if flutes == 0 {
            return Err(CamError::BadTool {
                reason: "flutes must be >= 1".into(),
            });
        }
        Ok(Self {
            id,
            name: name.into(),
            kind,
            diameter_mm,
            length_mm,
            flutes,
            material: material.into(),
        })
    }

    /// Cutting radius (half the diameter).
    pub fn radius_mm(&self) -> f64 {
        self.diameter_mm * 0.5
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn new_validates() {
        assert!(Tool::new(1, "EM6", ToolKind::EndMill, 6.0, 25.0, 2, "carbide").is_ok());
        let bad_dia = Tool::new(2, "x", ToolKind::EndMill, 0.0, 1.0, 2, "");
        assert!(matches!(bad_dia, Err(CamError::BadTool { .. })));
        let bad_len = Tool::new(3, "x", ToolKind::EndMill, 6.0, 0.0, 2, "");
        assert!(matches!(bad_len, Err(CamError::BadTool { .. })));
        let bad_flutes = Tool::new(4, "x", ToolKind::EndMill, 6.0, 25.0, 0, "");
        assert!(matches!(bad_flutes, Err(CamError::BadTool { .. })));
    }

    #[test]
    fn radius_is_half_diameter() {
        let t = Tool::new(1, "EM10", ToolKind::EndMill, 10.0, 30.0, 4, "carbide").unwrap();
        assert!((t.radius_mm() - 5.0).abs() < 1e-12);
    }

    #[test]
    fn kind_display() {
        assert_eq!(format!("{}", ToolKind::EndMill), "EndMill");
        assert_eq!(ToolKind::Drill.label(), "Drill");
    }
}
