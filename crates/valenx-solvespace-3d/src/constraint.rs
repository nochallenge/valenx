//! Constraint variants and their residual + Jacobian contributions.
//!
//! Each constraint contributes one or more scalar residual rows to the
//! global `r` vector and a corresponding rectangular Jacobian block.
//! All residuals are written so that `r == 0` ⇔ constraint satisfied.

use serde::{Deserialize, Serialize};

use crate::entity::{Entity3D, EntityId};
use crate::sketch::Sketch3D;

/// Triplet emitted into the global Jacobian matrix: `(row_offset_into_constraint, col, value)`.
pub type Triplet = (usize, usize, f64);

/// 3D constraint zoo.
#[derive(Copy, Clone, Debug, Serialize, Deserialize)]
pub enum Constraint3D {
    /// Two 3D points coincide. 3 residuals (Δx, Δy, Δz).
    Coincident3 {
        /// First point.
        a: EntityId,
        /// Second point.
        b: EntityId,
    },
    /// `target` distance between two points. 1 residual.
    PointDistance3 {
        /// First point.
        a: EntityId,
        /// Second point.
        b: EntityId,
        /// Target distance.
        target: f64,
    },
    /// A line's length equals `target`. 1 residual.
    LineLength3 {
        /// Line entity.
        line: EntityId,
        /// Target length.
        target: f64,
    },
    /// Angle (radians) between two lines equals `target`. 1 residual.
    /// Computed via the angle of the cross / dot of their direction
    /// vectors — see [`Constraint3D::residuals`] for the exact form.
    LineAngle3 {
        /// First line entity.
        a: EntityId,
        /// Second line entity.
        b: EntityId,
        /// Target angle (radians).
        target: f64,
    },
    /// Two lines parallel — direction cross product = 0. 3 residuals.
    LineParallel3 {
        /// First line entity.
        a: EntityId,
        /// Second line entity.
        b: EntityId,
    },
    /// Two lines perpendicular — direction dot product = 0. 1 residual.
    LinePerpendicular3 {
        /// First line entity.
        a: EntityId,
        /// Second line entity.
        b: EntityId,
    },
    /// Point lies on a plane: `n · (p - origin) = 0`. 1 residual.
    PointInPlane {
        /// Point id.
        point: EntityId,
        /// Plane id.
        plane: EntityId,
    },
    /// Point lies on a workplane — alias of `PointInPlane` against a
    /// `Workplane`. 1 residual.
    OnPlane {
        /// Point id.
        point: EntityId,
        /// Workplane id.
        workplane: EntityId,
    },
    /// Point lies on a line: `(p - a) × dir = 0` where `dir = b - a`.
    /// 3 residuals (one per component of the cross product).
    OnLine {
        /// Point id.
        point: EntityId,
        /// Line id.
        line: EntityId,
    },
    /// Coincident3 alias kept for symmetry — same scalar count as
    /// `Coincident3` but emphasises "this point is the same as that
    /// point even though they may have been added at different times".
    /// 3 residuals.
    SamePoint {
        /// First point.
        a: EntityId,
        /// Second point.
        b: EntityId,
    },
    /// Pin a plane as a fixed datum — holds its origin point and its
    /// normal direction at the captured values. 6 residuals (Δox, Δoy,
    /// Δoz, Δnx, Δny, Δnz). Without this a plane added with a free
    /// normal is itself a free body, so a single `PointInPlane` against
    /// it is under-determined (the plane can tilt to meet the point).
    /// Construct via [`Sketch3D::lock_plane`], which captures the
    /// current origin + normal.
    PlaneFixed {
        /// Plane (or workplane) id.
        plane: EntityId,
        /// Captured origin (ox, oy, oz).
        origin: [f64; 3],
        /// Captured normal (nx, ny, nz).
        normal: [f64; 3],
    },
}

impl Constraint3D {
    /// Number of scalar residual equations contributed.
    pub fn n_residuals(&self) -> usize {
        match self {
            Constraint3D::Coincident3 { .. } | Constraint3D::SamePoint { .. } => 3,
            Constraint3D::PointDistance3 { .. }
            | Constraint3D::LineLength3 { .. }
            | Constraint3D::LineAngle3 { .. }
            | Constraint3D::LinePerpendicular3 { .. }
            | Constraint3D::PointInPlane { .. }
            | Constraint3D::OnPlane { .. } => 1,
            Constraint3D::LineParallel3 { .. } | Constraint3D::OnLine { .. } => 3,
            Constraint3D::PlaneFixed { .. } => 6,
        }
    }

    /// Write residual values into `out` (length = `self.n_residuals()`).
    pub fn residuals(&self, s: &Sketch3D, out: &mut [f64]) {
        match self {
            Constraint3D::Coincident3 { a, b } | Constraint3D::SamePoint { a, b } => {
                let (ax, ay, az) = s.point_xyz(*a);
                let (bx, by, bz) = s.point_xyz(*b);
                out[0] = bx - ax;
                out[1] = by - ay;
                out[2] = bz - az;
            }
            Constraint3D::PointDistance3 { a, b, target } => {
                let (ax, ay, az) = s.point_xyz(*a);
                let (bx, by, bz) = s.point_xyz(*b);
                let dx = bx - ax;
                let dy = by - ay;
                let dz = bz - az;
                let len = (dx * dx + dy * dy + dz * dz).sqrt();
                out[0] = len - target;
            }
            Constraint3D::LineLength3 { line, target } => {
                let (a, b) = s.line_endpoints(*line);
                let (ax, ay, az) = s.point_xyz(a);
                let (bx, by, bz) = s.point_xyz(b);
                let dx = bx - ax;
                let dy = by - ay;
                let dz = bz - az;
                let len = (dx * dx + dy * dy + dz * dz).sqrt();
                out[0] = len - target;
            }
            Constraint3D::LineAngle3 { a, b, target } => {
                let (ua, ub) = (s.line_endpoints(*a), s.line_endpoints(*b));
                let (a1x, a1y, a1z) = s.point_xyz(ua.0);
                let (a2x, a2y, a2z) = s.point_xyz(ua.1);
                let (b1x, b1y, b1z) = s.point_xyz(ub.0);
                let (b2x, b2y, b2z) = s.point_xyz(ub.1);
                let (ax, ay, az) = (a2x - a1x, a2y - a1y, a2z - a1z);
                let (bx, by, bz) = (b2x - b1x, b2y - b1y, b2z - b1z);
                let la = (ax * ax + ay * ay + az * az).sqrt().max(1e-18);
                let lb = (bx * bx + by * by + bz * bz).sqrt().max(1e-18);
                let dot = ax * bx + ay * by + az * bz;
                let cos = (dot / (la * lb)).clamp(-1.0, 1.0);
                let ang = cos.acos();
                out[0] = ang - target;
            }
            Constraint3D::LineParallel3 { a, b } => {
                let (ua, ub) = (s.line_endpoints(*a), s.line_endpoints(*b));
                let (a1x, a1y, a1z) = s.point_xyz(ua.0);
                let (a2x, a2y, a2z) = s.point_xyz(ua.1);
                let (b1x, b1y, b1z) = s.point_xyz(ub.0);
                let (b2x, b2y, b2z) = s.point_xyz(ub.1);
                let (ax, ay, az) = (a2x - a1x, a2y - a1y, a2z - a1z);
                let (bx, by, bz) = (b2x - b1x, b2y - b1y, b2z - b1z);
                // Cross product components.
                out[0] = ay * bz - az * by;
                out[1] = az * bx - ax * bz;
                out[2] = ax * by - ay * bx;
            }
            Constraint3D::LinePerpendicular3 { a, b } => {
                let (ua, ub) = (s.line_endpoints(*a), s.line_endpoints(*b));
                let (a1x, a1y, a1z) = s.point_xyz(ua.0);
                let (a2x, a2y, a2z) = s.point_xyz(ua.1);
                let (b1x, b1y, b1z) = s.point_xyz(ub.0);
                let (b2x, b2y, b2z) = s.point_xyz(ub.1);
                let (ax, ay, az) = (a2x - a1x, a2y - a1y, a2z - a1z);
                let (bx, by, bz) = (b2x - b1x, b2y - b1y, b2z - b1z);
                out[0] = ax * bx + ay * by + az * bz;
            }
            Constraint3D::PointInPlane { point, plane } => {
                let (px, py, pz) = s.point_xyz(*point);
                let (ox, oy, oz, nx, ny, nz) = s.plane_data(*plane);
                // True geometric point-to-plane distance: n·(p-o)/‖n‖.
                // Dividing by ‖n‖ makes the residual scale-invariant in
                // the normal, so the solver cannot satisfy the
                // constraint by collapsing the (free) normal toward
                // zero — it must move the point.
                let nlen = (nx * nx + ny * ny + nz * nz).sqrt().max(1e-18);
                out[0] = (nx * (px - ox) + ny * (py - oy) + nz * (pz - oz)) / nlen;
            }
            Constraint3D::OnPlane { point, workplane } => {
                let (px, py, pz) = s.point_xyz(*point);
                let (ox, oy, oz, nx, ny, nz) = s.workplane_data(*workplane);
                let nlen = (nx * nx + ny * ny + nz * nz).sqrt().max(1e-18);
                out[0] = (nx * (px - ox) + ny * (py - oy) + nz * (pz - oz)) / nlen;
            }
            Constraint3D::OnLine { point, line } => {
                let (px, py, pz) = s.point_xyz(*point);
                let (a, b) = s.line_endpoints(*line);
                let (ax, ay, az) = s.point_xyz(a);
                let (bx, by, bz) = s.point_xyz(b);
                let (dx, dy, dz) = (bx - ax, by - ay, bz - az);
                let (rx, ry, rz) = (px - ax, py - ay, pz - az);
                out[0] = dy * rz - dz * ry;
                out[1] = dz * rx - dx * rz;
                out[2] = dx * ry - dy * rx;
            }
            Constraint3D::PlaneFixed {
                plane,
                origin,
                normal,
            } => {
                let (ox, oy, oz, nx, ny, nz) = s.plane_or_workplane_data(*plane);
                out[0] = ox - origin[0];
                out[1] = oy - origin[1];
                out[2] = oz - origin[2];
                out[3] = nx - normal[0];
                out[4] = ny - normal[1];
                out[5] = nz - normal[2];
            }
        }
    }

    /// Append Jacobian triplets `(row_within_constraint, column, value)`
    /// to `tri`. Uses analytic derivatives for the common shapes and
    /// numerical differentiation (central differences) for the harder
    /// trigonometric ones — accurate enough for Newton damping.
    pub fn jacobian_triplets(&self, s: &Sketch3D, tri: &mut Vec<Triplet>) {
        // Numerical Jacobian: cheap and robust for all variants. Each
        // constraint's residual count is small (≤3) so the cost is
        // negligible compared to the linear solve.
        let n = self.n_residuals();
        let cols = self.referenced_vars(s);
        let eps = 1e-7;
        for &c in &cols {
            let mut s_plus = s.clone();
            let mut s_minus = s.clone();
            s_plus.vars[c] += eps;
            s_minus.vars[c] -= eps;
            let mut r_plus = vec![0.0; n];
            let mut r_minus = vec![0.0; n];
            self.residuals(&s_plus, &mut r_plus);
            self.residuals(&s_minus, &mut r_minus);
            for row in 0..n {
                let d = (r_plus[row] - r_minus[row]) / (2.0 * eps);
                if d.abs() > 1e-18 {
                    tri.push((row, c, d));
                }
            }
        }
    }

    /// Set of variable indices this constraint reads. Used both by the
    /// numerical Jacobian (which only perturbs touched variables) and
    /// by downstream consumers that want to render dependency arrows.
    pub fn referenced_vars(&self, s: &Sketch3D) -> Vec<usize> {
        let mut out = Vec::new();
        let push_point = |id: EntityId, out: &mut Vec<usize>| {
            if let Some(Entity3D::Point3(p)) = s.entities.get(id.0) {
                out.push(p.x_var);
                out.push(p.y_var);
                out.push(p.z_var);
            }
        };
        let push_line = |id: EntityId, out: &mut Vec<usize>| {
            if let Some(Entity3D::Line3(l)) = s.entities.get(id.0) {
                push_point(l.a, out);
                push_point(l.b, out);
            }
        };
        match self {
            Constraint3D::Coincident3 { a, b }
            | Constraint3D::SamePoint { a, b }
            | Constraint3D::PointDistance3 { a, b, .. } => {
                push_point(*a, &mut out);
                push_point(*b, &mut out);
            }
            Constraint3D::LineLength3 { line, .. } => push_line(*line, &mut out),
            Constraint3D::LineAngle3 { a, b, .. }
            | Constraint3D::LineParallel3 { a, b }
            | Constraint3D::LinePerpendicular3 { a, b } => {
                push_line(*a, &mut out);
                push_line(*b, &mut out);
            }
            Constraint3D::PointInPlane { point, plane } => {
                push_point(*point, &mut out);
                if let Some(Entity3D::Plane3(p)) = s.entities.get(plane.0) {
                    push_point(p.origin, &mut out);
                    out.push(p.nx_var);
                    out.push(p.ny_var);
                    out.push(p.nz_var);
                }
            }
            Constraint3D::OnPlane { point, workplane } => {
                push_point(*point, &mut out);
                if let Some(Entity3D::Workplane(w)) = s.entities.get(workplane.0) {
                    push_point(w.origin, &mut out);
                    out.push(w.nx_var);
                    out.push(w.ny_var);
                    out.push(w.nz_var);
                }
            }
            Constraint3D::OnLine { point, line } => {
                push_point(*point, &mut out);
                push_line(*line, &mut out);
            }
            Constraint3D::PlaneFixed { plane, .. } => match s.entities.get(plane.0) {
                Some(Entity3D::Plane3(p)) => {
                    push_point(p.origin, &mut out);
                    out.push(p.nx_var);
                    out.push(p.ny_var);
                    out.push(p.nz_var);
                }
                Some(Entity3D::Workplane(w)) => {
                    push_point(w.origin, &mut out);
                    out.push(w.nx_var);
                    out.push(w.ny_var);
                    out.push(w.nz_var);
                }
                _ => {}
            },
        }
        out.sort_unstable();
        out.dedup();
        out
    }
}
