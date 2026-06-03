//! Fixture / vise / clamp collision-check (Phase 17F).
//!
//! v1 represents a fixture as a list of axis-aligned bounding boxes
//! (AABBs) in world space. The tool body is approximated as a
//! cylinder of `tool.diameter * tool.length`. The collision check
//! samples the toolpath at every move and asks "does the tool
//! cylinder intersect any fixture AABB?".
//!
//! v1 limitations:
//!
//! - **AABB-only fixtures** — no mesh-level collision; complex vise
//!   geometry (T-slot bars, swivel jaws) is out of scope.
//! - **Cylindrical tool approximation** — shank vs flute distinction
//!   is ignored (the tool is one homogeneous cylinder).

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

use crate::tool::Tool;
use crate::toolpath::Toolpath;

/// A single fixture AABB.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct FixtureAabb {
    /// Min XYZ corner.
    pub min: Vector3<f64>,
    /// Max XYZ corner.
    pub max: Vector3<f64>,
    /// Free-form fixture-component label (`"vise jaw"`, `"clamp"`).
    pub label: String,
}

/// Fixture description for collision checks.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct Fixture {
    /// AABB list. v1 unions all AABBs as the keep-out region.
    pub aabbs: Vec<FixtureAabb>,
}

impl Fixture {
    /// Empty fixture (no keep-out region).
    pub fn new() -> Self {
        Self::default()
    }

    /// Append an AABB.
    pub fn push_aabb(&mut self, min: Vector3<f64>, max: Vector3<f64>, label: impl Into<String>) {
        self.aabbs.push(FixtureAabb {
            min,
            max,
            label: label.into(),
        });
    }
}

/// Collision report for one detected violation.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Collision {
    /// Index of the offending [`crate::toolpath::Move`] in the toolpath.
    pub move_index: usize,
    /// Label of the fixture AABB the tool hit.
    pub fixture_label: String,
    /// Tool centre at the moment of collision.
    pub tool_position: Vector3<f64>,
}

/// Test the tool cylinder (centred at `centre`, radius `r`, axial
/// length `len` extending downward in -Z) for overlap with `aabb`.
fn cylinder_aabb_overlap(centre: Vector3<f64>, r: f64, len: f64, aabb: &FixtureAabb) -> bool {
    // Tool extends from `centre.z` down to `centre.z - len`.
    let tool_z_max = centre.z;
    let tool_z_min = centre.z - len;
    if tool_z_max < aabb.min.z || tool_z_min > aabb.max.z {
        return false;
    }
    // XY: closest-point AABB-to-circle.
    let cx = centre.x.clamp(aabb.min.x, aabb.max.x);
    let cy = centre.y.clamp(aabb.min.y, aabb.max.y);
    let dx = cx - centre.x;
    let dy = cy - centre.y;
    (dx * dx + dy * dy) <= r * r
}

/// Sample the toolpath and report every move whose tool cylinder
/// overlaps a fixture AABB.
pub fn collision_check(toolpath: &Toolpath, tool: &Tool, fixture: &Fixture) -> Vec<Collision> {
    let mut out = Vec::new();
    let r = tool.radius_mm();
    let len = tool.length_mm;
    for (i, m) in toolpath.moves.iter().enumerate() {
        for aabb in &fixture.aabbs {
            if cylinder_aabb_overlap(m.position, r, len, aabb) {
                out.push(Collision {
                    move_index: i,
                    fixture_label: aabb.label.clone(),
                    tool_position: m.position,
                });
            }
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tool::ToolKind;
    use crate::toolpath::{Move, MoveKind};

    #[test]
    fn empty_fixture_yields_no_collisions() {
        let tool = Tool::new(1, "EM6", ToolKind::EndMill, 6.0, 25.0, 2, "").unwrap();
        let mut tp = Toolpath::new();
        tp.push(Move::new(MoveKind::Rapid, Vector3::new(0.0, 0.0, 5.0), 0.0));
        let fix = Fixture::new();
        assert!(collision_check(&tp, &tool, &fix).is_empty());
    }

    #[test]
    fn collision_when_tool_passes_through_aabb() {
        let tool = Tool::new(1, "EM6", ToolKind::EndMill, 6.0, 25.0, 2, "").unwrap();
        let mut tp = Toolpath::new();
        tp.push(Move::new(
            MoveKind::Rapid,
            Vector3::new(50.0, 50.0, 10.0),
            0.0,
        ));
        let mut fix = Fixture::new();
        // Tool centre at z=10 extends from 10 down to -15. AABB sits
        // in that span.
        fix.push_aabb(
            Vector3::new(45.0, 45.0, 0.0),
            Vector3::new(55.0, 55.0, 5.0),
            "vise jaw",
        );
        let hits = collision_check(&tp, &tool, &fix);
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].fixture_label, "vise jaw");
    }

    #[test]
    fn no_collision_when_tool_above_fixture() {
        let tool = Tool::new(1, "EM6", ToolKind::EndMill, 6.0, 25.0, 2, "").unwrap();
        let mut tp = Toolpath::new();
        // Tool centre at z=200, extends from 175 down — well above.
        tp.push(Move::new(
            MoveKind::Rapid,
            Vector3::new(50.0, 50.0, 200.0),
            0.0,
        ));
        let mut fix = Fixture::new();
        fix.push_aabb(
            Vector3::new(45.0, 45.0, 0.0),
            Vector3::new(55.0, 55.0, 5.0),
            "vise jaw",
        );
        assert!(collision_check(&tp, &tool, &fix).is_empty());
    }
}
