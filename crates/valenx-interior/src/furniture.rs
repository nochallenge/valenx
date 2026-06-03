//! Furniture catalog — 12 item kinds with parametric block geometry.
//!
//! Each kind has a default bounding-box size in metres (length,
//! width, height). [`to_solid`] returns a single axis-aligned box
//! representation suitable for the v1 viewport (a placeholder
//! suitable until the BRep import pipeline handles SH3D's `.obj`
//! catalog).

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

/// One catalog entry kind.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash, Serialize, Deserialize)]
pub enum Furniture {
    /// Chair.
    Chair,
    /// Table (dining or coffee).
    Table,
    /// Bed.
    Bed,
    /// Sofa.
    Sofa,
    /// Desk.
    Desk,
    /// Wardrobe.
    Wardrobe,
    /// Bookcase.
    Bookcase,
    /// Stove (kitchen).
    Stove,
    /// Refrigerator.
    Fridge,
    /// Toilet (bathroom).
    Toilet,
    /// Sink (bathroom or kitchen).
    Sink,
    /// Bathtub.
    Bathtub,
}

impl Furniture {
    /// Default size in metres `(length_x, width_y, height_z)`.
    pub fn default_size(self) -> Vector3<f64> {
        match self {
            Self::Chair => Vector3::new(0.50, 0.50, 0.90),
            Self::Table => Vector3::new(1.60, 0.90, 0.75),
            Self::Bed => Vector3::new(2.00, 1.60, 0.55),
            Self::Sofa => Vector3::new(2.10, 0.90, 0.85),
            Self::Desk => Vector3::new(1.40, 0.70, 0.75),
            Self::Wardrobe => Vector3::new(1.20, 0.60, 2.10),
            Self::Bookcase => Vector3::new(0.90, 0.30, 2.10),
            Self::Stove => Vector3::new(0.60, 0.60, 0.90),
            Self::Fridge => Vector3::new(0.65, 0.65, 1.80),
            Self::Toilet => Vector3::new(0.40, 0.65, 0.80),
            Self::Sink => Vector3::new(0.60, 0.45, 0.20),
            Self::Bathtub => Vector3::new(1.70, 0.75, 0.55),
        }
    }

    /// Catalog of all 12 kinds (handy for the palette picker).
    pub fn all() -> &'static [Furniture] {
        &[
            Self::Chair,
            Self::Table,
            Self::Bed,
            Self::Sofa,
            Self::Desk,
            Self::Wardrobe,
            Self::Bookcase,
            Self::Stove,
            Self::Fridge,
            Self::Toilet,
            Self::Sink,
            Self::Bathtub,
        ]
    }

    /// Stable kebab-style name (matches Sweet Home 3D's lowercase
    /// catalog references).
    pub fn name(self) -> &'static str {
        match self {
            Self::Chair => "chair",
            Self::Table => "table",
            Self::Bed => "bed",
            Self::Sofa => "sofa",
            Self::Desk => "desk",
            Self::Wardrobe => "wardrobe",
            Self::Bookcase => "bookcase",
            Self::Stove => "stove",
            Self::Fridge => "fridge",
            Self::Toilet => "toilet",
            Self::Sink => "sink",
            Self::Bathtub => "bathtub",
        }
    }
}

/// Lightweight axis-aligned box used as a v1 stand-in for the BRep
/// `valenx_cad::Solid`. The caller can pipe this into the
/// `truck-modeling` box primitive when they need an actual Solid.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Solid {
    /// Origin of the box (lowest corner).
    pub origin: Vector3<f64>,
    /// Box size `(length_x, width_y, height_z)`.
    pub size: Vector3<f64>,
    /// Tagged catalog kind for round-trip identification.
    pub kind: Furniture,
}

/// Build a v1 axis-aligned box [`Solid`] for `kind` at the origin.
pub fn to_solid(kind: Furniture, size: Vector3<f64>) -> Solid {
    Solid {
        origin: Vector3::zeros(),
        size,
        kind,
    }
}

/// Placement of a furniture item in world coordinates.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Placement {
    /// Furniture kind.
    pub kind: Furniture,
    /// World position of the box's `origin` corner.
    pub position: Vector3<f64>,
    /// Yaw rotation about +Z axis in radians.
    pub rotation_rad: f64,
    /// Room id this placement belongs to.
    pub room_id: String,
}

impl Furniture {
    /// Construct a [`Placement`] of `self` at the given position +
    /// rotation in the named room.
    pub fn place(
        self,
        position: Vector3<f64>,
        rotation_rad: f64,
        room_id: impl Into<String>,
    ) -> Placement {
        Placement {
            kind: self,
            position,
            rotation_rad,
            room_id: room_id.into(),
        }
    }
}
