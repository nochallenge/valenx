//! Printer + Part data model.

use nalgebra::UnitQuaternion;
use serde::{Deserialize, Serialize};

use valenx_mesh::Mesh;

/// Bed heating capability.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub enum BedType {
    /// Heated.
    Heated,
    /// Unheated (PLA-only printers).
    Unheated,
}

/// Bed surface material.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub enum BedMaterial {
    /// Tempered glass.
    Glass,
    /// PEI sheet.
    Pei,
    /// BuildTak surface.
    BuildTak,
}

/// A 3D printer.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Printer {
    /// Bed bounding-volume size in mm (x_width, y_depth, z_height).
    pub bed_size: (f64, f64, f64),
    /// Bed type.
    pub bed_type: BedType,
    /// Bed surface material.
    pub bed_material: BedMaterial,
}

impl Printer {
    /// New Printer.
    pub fn new(bed_size: (f64, f64, f64), bed_type: BedType, bed_material: BedMaterial) -> Self {
        Self {
            bed_size,
            bed_type,
            bed_material,
        }
    }
}

/// A part ready to print: its mesh, an orientation quaternion, and a
/// 2D bed position.
#[derive(Clone, Debug)]
pub struct Part {
    /// Short identifier for UI / error display.
    pub name: String,
    /// Triangle mesh in the part's local frame.
    pub mesh: Mesh,
    /// Orientation that will be applied before laying on the bed.
    pub orientation: UnitQuaternion<f64>,
    /// Center position on the bed `(x, y)` in mm.
    pub bed_position: [f64; 2],
}

impl Part {
    /// New Part with identity orientation at the bed origin.
    pub fn new(name: impl Into<String>, mesh: Mesh) -> Self {
        Self {
            name: name.into(),
            mesh,
            orientation: UnitQuaternion::identity(),
            bed_position: [0.0, 0.0],
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn printer_holds_bed_size() {
        let p = Printer::new((220.0, 220.0, 250.0), BedType::Heated, BedMaterial::Pei);
        assert_eq!(p.bed_size.0, 220.0);
    }

    #[test]
    fn default_part_has_identity_orientation() {
        let m = Mesh::new("empty");
        let p = Part::new("part1", m);
        assert!(p.orientation == UnitQuaternion::identity());
        assert_eq!(p.bed_position, [0.0, 0.0]);
    }
}
