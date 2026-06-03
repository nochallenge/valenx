//! Arc3d primitive — 3-point and tangent-to-line constructors.

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

use crate::error::Gcad3dError;
use crate::line::Line3d;

/// 3D arc — centre + normal + two axes + radius + angle sweep.
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
pub struct Arc3d {
    /// Centre of the arc's circle.
    pub centre: Vector3<f64>,
    /// Plane normal (right-hand sweep along this).
    pub normal: Vector3<f64>,
    /// First in-plane axis (unit) — corresponds to angle 0.
    pub u: Vector3<f64>,
    /// Second in-plane axis (unit) — corresponds to angle pi/2.
    pub v: Vector3<f64>,
    /// Radius.
    pub radius: f64,
    /// Start angle in radians.
    pub start_angle: f64,
    /// Sweep angle in radians (positive = CCW about normal).
    pub sweep: f64,
}

impl Arc3d {
    /// Arc through three points. Fails if colinear.
    pub fn three_point(
        p1: Vector3<f64>,
        p2: Vector3<f64>,
        p3: Vector3<f64>,
    ) -> Result<Self, Gcad3dError> {
        let a = p2 - p1;
        let b = p3 - p1;
        let n = a.cross(&b);
        let n_len = n.norm();
        if n_len < 1e-12 {
            return Err(Gcad3dError::Degenerate(
                "three_point: points are colinear".into(),
            ));
        }
        let normal = n / n_len;

        // Circumcentre of the triangle, written relative to p1 with
        // a = p2−p1, b = p3−p1, n = a×b:
        //   centre − p1 = ( |a|²·(b×n) + |b|²·(n×a) ) / (2·|n|²)
        // (the standard closed form — see e.g. the circumscribed-circle
        // identity in Schneider & Eberly, *Geometric Tools*).
        let a2 = a.dot(&a);
        let b2 = b.dot(&b);
        let denom = 2.0 * n.dot(&n);
        if denom < 1e-18 {
            return Err(Gcad3dError::Degenerate(
                "three_point: degenerate triangle".into(),
            ));
        }
        let centre =
            p1 + (b.cross(&n) * a2 + n.cross(&a) * b2) / denom;
        let radius = (p1 - centre).norm();
        let u_axis = (p1 - centre) / radius;
        let v_axis = normal.cross(&u_axis);
        let p2_u = (p2 - centre).dot(&u_axis);
        let p2_v = (p2 - centre).dot(&v_axis);
        let p3_u = (p3 - centre).dot(&u_axis);
        let p3_v = (p3 - centre).dot(&v_axis);
        let theta2 = p2_v.atan2(p2_u);
        let theta3 = p3_v.atan2(p3_u);
        // Sweep from 0 (p1) to theta3, going through theta2 — pick
        // sign by the sign of theta2.
        let sweep = if theta2 >= 0.0 { theta3.rem_euclid(2.0 * std::f64::consts::PI) } else { -((-theta3).rem_euclid(2.0 * std::f64::consts::PI)) };
        Ok(Self {
            centre,
            normal,
            u: u_axis,
            v: v_axis,
            radius,
            start_angle: 0.0,
            sweep,
        })
    }

    /// Arc tangent to `line` with given `radius`, centred at the
    /// foot of the perpendicular from `line.origin + line.direction * 0`
    /// (i.e. at `line.origin`), in the plane spanned by `line.direction`
    /// and a stable perpendicular.
    pub fn tangent_to(line: &Line3d, radius: f64) -> Result<Self, Gcad3dError> {
        if !radius.is_finite() || radius <= 0.0 {
            return Err(Gcad3dError::BadParameter {
                name: "radius",
                reason: format!("must be > 0 (got {radius})"),
            });
        }
        // Stable perpendicular — same trick as Line3d::parallel_to.
        let candidates = [Vector3::x(), Vector3::y(), Vector3::z()];
        let perp = candidates
            .iter()
            .map(|c| line.direction.cross(c))
            .max_by(|a, b| a.norm().partial_cmp(&b.norm()).unwrap())
            .unwrap();
        let perp_n = perp.normalize();
        let normal = line.direction.cross(&perp_n).normalize();
        let centre = line.origin + perp_n * radius;
        Ok(Self {
            centre,
            normal,
            u: -perp_n, // start at the tangent point
            v: normal.cross(&-perp_n),
            radius,
            start_angle: 0.0,
            sweep: std::f64::consts::PI,
        })
    }
}
