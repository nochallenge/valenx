//! 2-D geometric primitives for the sketcher.
//!
//! Each primitive holds *variable indices*, not coordinates directly.
//! The actual variable values live in [`crate::sketch::Sketch::vars`],
//! a flat `Vec<f64>` that the solver mutates. This layout lets the
//! solver pack the parameter vector tightly for Newton-Raphson.

use serde::{Deserialize, Serialize};

/// Stable handle to an entity inside a sketch. 1-based for human
/// readability (entity #0 is reserved for "none").
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash, Serialize, Deserialize)]
pub struct EntityId(pub usize);

/// 2-D point — holds the variable indices of its (x, y) coordinates
/// in the parent sketch's variable vector.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash, Serialize, Deserialize)]
pub struct Point2 {
    /// Index into [`crate::sketch::Sketch::vars`] for the x coordinate.
    pub x_var: usize,
    /// Index into [`crate::sketch::Sketch::vars`] for the y coordinate.
    pub y_var: usize,
}

impl Point2 {
    /// Read the current (x, y) from a variable vector.
    ///
    /// Out-of-range handles read as `0.0` rather than panicking. The
    /// primary defence against a corrupt handle is
    /// [`crate::sketch::Sketch::validate`] (reject-early at load); this
    /// `.get()` is belt-and-suspenders so any un-validated path
    /// (e.g. a future caller) can never panic with "index out of
    /// bounds". R33 H1.
    pub fn read(&self, vars: &[f64]) -> (f64, f64) {
        (
            vars.get(self.x_var).copied().unwrap_or(0.0),
            vars.get(self.y_var).copied().unwrap_or(0.0),
        )
    }
}

/// Straight line segment between two [`Point2`]s.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash, Serialize, Deserialize)]
pub struct Line2 {
    /// Start endpoint.
    pub start: Point2,
    /// End endpoint.
    pub end: Point2,
}

impl Line2 {
    /// Both endpoint coordinates as `((sx, sy), (ex, ey))`.
    pub fn endpoints(&self, vars: &[f64]) -> ((f64, f64), (f64, f64)) {
        (self.start.read(vars), self.end.read(vars))
    }

    /// Euclidean length of the segment.
    pub fn length(&self, vars: &[f64]) -> f64 {
        let ((sx, sy), (ex, ey)) = self.endpoints(vars);
        let dx = ex - sx;
        let dy = ey - sy;
        (dx * dx + dy * dy).sqrt()
    }

    /// Direction vector (unnormalized) from start to end.
    pub fn direction(&self, vars: &[f64]) -> (f64, f64) {
        let ((sx, sy), (ex, ey)) = self.endpoints(vars);
        (ex - sx, ey - sy)
    }
}

/// Circle defined by its centre point and a radius variable.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash, Serialize, Deserialize)]
pub struct Circle2 {
    /// Centre point.
    pub center: Point2,
    /// Index into [`crate::sketch::Sketch::vars`] for the radius.
    pub radius_var: usize,
}

impl Circle2 {
    /// Current radius from a variable vector. Out-of-range handle reads
    /// as `0.0` (R33 H1 defense-in-depth — see [`Point2::read`]).
    pub fn radius(&self, vars: &[f64]) -> f64 {
        vars.get(self.radius_var).copied().unwrap_or(0.0)
    }
}

/// Circular arc — same as [`Circle2`] but with start and end angles
/// (in radians, measured CCW from +x). Sweep is always end - start.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash, Serialize, Deserialize)]
pub struct Arc2 {
    /// Centre point.
    pub center: Point2,
    /// Radius variable index.
    pub radius_var: usize,
    /// Start angle variable index (radians).
    pub start_angle_var: usize,
    /// End angle variable index (radians).
    pub end_angle_var: usize,
}

impl Arc2 {
    /// Current radius from a variable vector. Out-of-range handle reads
    /// as `0.0` (R33 H1 defense-in-depth — see [`Point2::read`]).
    pub fn radius(&self, vars: &[f64]) -> f64 {
        vars.get(self.radius_var).copied().unwrap_or(0.0)
    }

    /// (start_angle, end_angle) in radians. Out-of-range handles read
    /// as `0.0` (R33 H1 defense-in-depth — see [`Point2::read`]).
    pub fn angles(&self, vars: &[f64]) -> (f64, f64) {
        (
            vars.get(self.start_angle_var).copied().unwrap_or(0.0),
            vars.get(self.end_angle_var).copied().unwrap_or(0.0),
        )
    }

    /// Sweep = end - start (radians; may be negative for CW arcs).
    pub fn sweep(&self, vars: &[f64]) -> f64 {
        let (s, e) = self.angles(vars);
        e - s
    }
}

/// Sum type over all primitive kinds — the entity table stores these.
///
/// `PartialEq` is derived structurally — two entities compare equal
/// when their variant kind matches and every nested field
/// (variable index, embedded knot vector, etc.) compares equal.
/// The float fields (`BSpline2::knots`, `BSpline2::weights`)
/// use IEEE 754 semantics — NaNs never compare equal, which means an
/// `Entity::BSpline` whose params contain NaNs will fail to dedupe
/// in `History` snapshots. In practice the sketcher never lets a NaN
/// reach a knot / weight, so this is a documented degradation
/// rather than a correctness issue.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum Entity {
    /// 2-D point.
    Point(Point2),
    /// Line segment.
    Line(Line2),
    /// Full circle.
    Circle(Circle2),
    /// Circular arc.
    Arc(Arc2),
    /// B-spline curve (Phase 12A).
    BSpline(crate::geom_bspline::BSpline2),
    /// Full ellipse (Phase 12A).
    Ellipse(crate::geom_ellipse::Ellipse2),
    /// Elliptical arc (Phase 12A).
    EllipticalArc(crate::geom_ellipse::EllipticalArc2),
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn entity_id_is_copy_and_eq() {
        let a = EntityId(7);
        let b = a;
        assert_eq!(a, b);
        assert_eq!(a.0, 7);
    }

    #[test]
    fn point_reads_from_var_vector() {
        let p = Point2 { x_var: 2, y_var: 3 };
        let vars = vec![0.0, 0.0, 1.5, -2.0];
        assert_eq!(p.read(&vars), (1.5, -2.0));
    }

    #[test]
    fn line_returns_endpoint_coords() {
        let line = Line2 {
            start: Point2 { x_var: 0, y_var: 1 },
            end: Point2 { x_var: 2, y_var: 3 },
        };
        let vars = vec![0.0, 0.0, 3.0, 4.0];
        let ((sx, sy), (ex, ey)) = line.endpoints(&vars);
        assert_eq!((sx, sy), (0.0, 0.0));
        assert_eq!((ex, ey), (3.0, 4.0));
        // Length = sqrt(9 + 16) = 5.
        assert!((line.length(&vars) - 5.0).abs() < 1e-12);
    }

    #[test]
    fn circle_returns_center_and_radius() {
        let c = Circle2 {
            center: Point2 { x_var: 0, y_var: 1 },
            radius_var: 2,
        };
        let vars = vec![10.0, 20.0, 5.0];
        let (cx, cy) = c.center.read(&vars);
        assert_eq!((cx, cy), (10.0, 20.0));
        assert_eq!(c.radius(&vars), 5.0);
    }

    // R33 H1 defense-in-depth: an out-of-range handle on the hot read
    // path must degrade to 0.0 rather than panic with "index out of
    // bounds". `validate()` is the primary reject-early gate; this keeps
    // any un-validated path (e.g. a future caller) panic-free too.
    #[test]
    fn point_read_with_out_of_range_handle_does_not_panic() {
        let p = Point2 {
            x_var: 999,
            y_var: 0,
        };
        let vars = vec![1.5];
        // Pre-fix this indexed vars[999] and panicked.
        let (x, y) = p.read(&vars);
        assert_eq!(x, 0.0, "out-of-range x_var reads as 0.0");
        assert_eq!(y, 1.5);
    }

    #[test]
    fn circle_radius_with_out_of_range_handle_does_not_panic() {
        let c = Circle2 {
            center: Point2 { x_var: 0, y_var: 0 },
            radius_var: 999,
        };
        let vars = vec![0.0];
        assert_eq!(c.radius(&vars), 0.0);
    }

    #[test]
    fn arc_angles_with_out_of_range_handle_does_not_panic() {
        let a = Arc2 {
            center: Point2 { x_var: 0, y_var: 0 },
            radius_var: 0,
            start_angle_var: 999,
            end_angle_var: 999,
        };
        let vars = vec![0.0];
        let (s, e) = a.angles(&vars);
        assert_eq!((s, e), (0.0, 0.0));
        assert_eq!(a.sweep(&vars), 0.0);
    }

    #[test]
    fn arc_returns_center_radius_and_sweep() {
        let a = Arc2 {
            center: Point2 { x_var: 0, y_var: 1 },
            radius_var: 2,
            start_angle_var: 3,
            end_angle_var: 4,
        };
        let vars = vec![0.0, 0.0, 1.0, 0.0, std::f64::consts::PI];
        assert_eq!(a.radius(&vars), 1.0);
        let (start, end) = a.angles(&vars);
        assert!((start - 0.0).abs() < 1e-12);
        assert!((end - std::f64::consts::PI).abs() < 1e-12);
        let sweep = a.sweep(&vars);
        assert!((sweep - std::f64::consts::PI).abs() < 1e-12);
    }
}
