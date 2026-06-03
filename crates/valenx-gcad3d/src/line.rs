//! Line3d primitive — infinite line through a point with a direction.

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

use crate::error::Gcad3dError;

/// Infinite line.
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
pub struct Line3d {
    /// Origin point.
    pub origin: Vector3<f64>,
    /// Unit direction.
    pub direction: Vector3<f64>,
}

impl Line3d {
    /// Line between two points.
    pub fn between(p1: Vector3<f64>, p2: Vector3<f64>) -> Result<Self, Gcad3dError> {
        let d = p2 - p1;
        let l = d.norm();
        if l < 1e-12 {
            return Err(Gcad3dError::Degenerate(
                "between: p1 and p2 coincide".into(),
            ));
        }
        Ok(Self {
            origin: p1,
            direction: d / l,
        })
    }

    /// Line parallel to `line` offset by `offset` along an arbitrary
    /// perpendicular (the lexicographically-first standard basis
    /// vector that's not colinear with the line direction).
    pub fn parallel_to(line: &Line3d, offset: f64) -> Result<Self, Gcad3dError> {
        if !offset.is_finite() {
            return Err(Gcad3dError::BadParameter {
                name: "offset",
                reason: "must be finite".into(),
            });
        }
        // Pick a basis vector not parallel to direction.
        let candidates = [Vector3::x(), Vector3::y(), Vector3::z()];
        let perp = candidates
            .iter()
            .map(|c| line.direction.cross(c))
            .max_by(|a, b| a.norm().partial_cmp(&b.norm()).unwrap())
            .unwrap();
        let perp_n = perp.normalize();
        Ok(Self {
            origin: line.origin + perp_n * offset,
            direction: line.direction,
        })
    }

    /// Line perpendicular to `line` passing through `point`. The
    /// perpendicular direction is the rejection of `point - origin`
    /// from `line.direction`.
    pub fn perpendicular_at(line: &Line3d, point: Vector3<f64>) -> Result<Self, Gcad3dError> {
        let v = point - line.origin;
        let proj = line.direction * line.direction.dot(&v);
        let perp = v - proj;
        let l = perp.norm();
        if l < 1e-12 {
            return Err(Gcad3dError::Degenerate(
                "perpendicular_at: point lies on the line".into(),
            ));
        }
        Ok(Self {
            origin: point,
            direction: perp / l,
        })
    }
}
