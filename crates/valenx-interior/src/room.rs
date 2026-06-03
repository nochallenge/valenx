//! Room + wall reference + floor polygon helpers.
//!
//! The interior workflow keeps walls as opaque `WallRef`s — the
//! actual wall geometry lives in `valenx-arch`. This module avoids
//! a hard dependency on the arch crate so the interior workbench
//! can compile standalone in tests.

use nalgebra::Vector2;
use serde::{Deserialize, Serialize};

/// Opaque reference to a `valenx-arch` wall. Just a string id — the
/// arch crate is the source of truth for wall geometry.
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct WallRef {
    /// String id mapping to a `valenx-arch::Wall`.
    pub id: String,
}

impl WallRef {
    /// Construct from an id-like string.
    pub fn new(id: impl Into<String>) -> Self {
        Self { id: id.into() }
    }
}

/// A single room — a closed floor polygon, ceiling height, label,
/// and the walls bounding it.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct Room {
    /// Room id (unique within an [`crate::panel::InteriorPanelState`]).
    pub id: String,
    /// Bounding walls (zero-or-more).
    pub walls: Vec<WallRef>,
    /// Floor polygon (CCW for outward normal facing +Z).
    pub floor_polygon: Vec<Vector2<f64>>,
    /// Ceiling height in metres.
    pub ceiling_height: f64,
    /// Display label (e.g. `"Kitchen"`).
    pub label: String,
    /// Cached area in m^2 (recomputed via [`compute_area`]).
    pub area_m2: f64,
}

impl Room {
    /// New empty room with `id` + `label`.
    pub fn new(id: impl Into<String>, label: impl Into<String>, ceiling_height: f64) -> Self {
        Self {
            id: id.into(),
            walls: Vec::new(),
            floor_polygon: Vec::new(),
            ceiling_height,
            label: label.into(),
            area_m2: 0.0,
        }
    }
}

/// Compute the signed area of a CCW polygon (m^2 if input is in metres).
pub fn compute_area(polygon: &[Vector2<f64>]) -> f64 {
    let n = polygon.len();
    if n < 3 {
        return 0.0;
    }
    let mut acc = 0.0;
    for i in 0..n {
        let p = polygon[i];
        let q = polygon[(i + 1) % n];
        acc += p.x * q.y - q.x * p.y;
    }
    (acc * 0.5).abs()
}
