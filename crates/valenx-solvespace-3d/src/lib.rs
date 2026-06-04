//! # valenx-solvespace-3d
//!
//! SolveSpace-style 3D mechanical constraint solver — Phase 53 of the
//! FreeCAD-parity roadmap. This crate is the **3D
//! analogue** of `valenx_sketch`: same structure-of-arrays variable
//! layout, same Newton-Raphson + Levenberg-Marquardt solver, but with
//! 3D entities ([`entity::Point3`], [`entity::Line3`], [`entity::Plane3`],
//! [`entity::Workplane`]) and a 3D constraint zoo
//! ([`constraint::Constraint3D`]).
//!
//! # Why a separate crate
//!
//! Phase 1's `valenx_sketch` is hard-coded to a 2D sketch plane: every
//! point owns exactly two variables (x, y), every constraint emits 1-2
//! scalar residuals, and the Jacobian assembly walks a 2D-shaped variable
//! vector. The 3D analogue needs three variables per point, plane
//! entities with their own (free) normal vector, and constraints whose
//! residual count varies (1 for distance, 3 for parallel, etc.). Rather
//! than retrofit a generic `N`-dim sketcher we keep the 2D crate clean
//! and reuse only the **algorithmic** pattern here.
//!
//! # Quick example
//!
//! ```
//! use valenx_solvespace_3d::{Constraint3D, Sketch3D};
//!
//! let mut s = Sketch3D::new();
//! let a = s.add_point(0.0, 0.0, 0.0);
//! let b = s.add_point(3.0, 4.0, 12.0);
//! s.add_constraint(Constraint3D::Coincident3 { a, b });
//! let report = s.solve().unwrap();
//! assert!(matches!(
//!     report.status,
//!     valenx_solvespace_3d::SolverStatus::Converged
//! ));
//! ```
//!
//! # Module map
//!
//! - [`entity`]   — `Point3` / `Line3` / `Plane3` / `Workplane`.
//! - [`constraint`] — `Constraint3D` enum + residual / Jacobian impls.
//! - [`sketch`]   — `Sketch3D` owning vars + entities + constraints.
//! - [`solver`]   — Newton-LM driver (configurable `SolverConfig`).
//! - [`panel`]    — UI panel state envelope.
//! - [`persist`]  — RON round-trip.
//! - [`error`]    — typed `Solve3DError`.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod constraint;
pub mod entity;
pub mod error;
pub mod panel;
pub mod parameters;
pub mod persist;
pub mod sketch;
pub mod solver;
pub mod timeline;

pub use constraint::Constraint3D;
pub use entity::{Arc3, Circle3, Entity3D, EntityId, Line3, Plane3, Point3, Spline3, Workplane};
pub use error::{ErrorCategory, Solve3DError};
pub use panel::SolveSpace3DPanelState;
pub use parameters::{ParamError, ParameterTable};
pub use persist::{
    from_ron_str, timeline_from_ron, timeline_to_ron, to_ron_string, PanelFile, TimelineFile,
    TIMELINE_VERSION, VERSION,
};
pub use sketch::Sketch3D;
pub use solver::{SolverConfig, SolverDiagnostics, SolverReport, SolverStatus};
pub use timeline::{Feature, FeatureTimeline, Op, RebuiltModel, Step, TimelineError};

#[cfg(test)]
mod tests_integration {
    use super::*;

    /// Two free 3D points + a Coincident3 — should land them on top of
    /// each other.
    #[test]
    fn two_coincident_points_converge() {
        let mut s = Sketch3D::new();
        let a = s.add_point(0.0, 0.0, 0.0);
        let b = s.add_point(5.0, 6.0, 7.0);
        s.add_constraint(Constraint3D::Coincident3 { a, b });
        let rep = s.solve().expect("ok");
        assert_eq!(rep.status, SolverStatus::Converged);
        let (ax, ay, az) = s.point_xyz(a);
        let (bx, by, bz) = s.point_xyz(b);
        assert!((ax - bx).abs() < 1e-6);
        assert!((ay - by).abs() < 1e-6);
        assert!((az - bz).abs() < 1e-6);
    }

    /// Point + plane + PointInPlane — the point should land on the
    /// plane. The Z = 0 plane is pinned as a datum first (a plane with
    /// a free normal is itself a free body, so an unpinned plane would
    /// tilt to meet the point instead).
    #[test]
    fn point_lands_on_plane() {
        let mut s = Sketch3D::new();
        let origin = s.add_point(0.0, 0.0, 0.0);
        // Z = 0 plane (normal = +Z), pinned as a fixed datum.
        let plane = s.add_plane(origin, 0.0, 0.0, 1.0).unwrap();
        s.lock_plane(plane).unwrap();
        let p = s.add_point(1.0, 2.0, 5.0);
        s.add_constraint(Constraint3D::PointInPlane { point: p, plane });
        let rep = s.solve().expect("ok");
        assert_eq!(rep.status, SolverStatus::Converged);
        let (px, py, pz) = s.point_xyz(p);
        // The point drops straight onto z = 0; x and y are unconstrained
        // and stay essentially put (a sub-1e-4 drift is numerical-
        // Jacobian noise — ∂r/∂px is analytically 0 for a +Z normal).
        assert!(pz.abs() < 1e-6, "z = {pz}");
        assert!((px - 1.0).abs() < 1e-4, "x moved: {px}");
        assert!((py - 2.0).abs() < 1e-4, "y moved: {py}");
    }

    /// `lock_plane` holds a plane's normal direction — an `OnPlane` of
    /// a point against a tilted, *unpinned* workplane is satisfiable by
    /// rotating the plane, but once locked the plane stays put.
    #[test]
    fn lock_plane_holds_the_normal() {
        let mut s = Sketch3D::new();
        let o = s.add_point(0.0, 0.0, 0.0);
        let plane = s.add_plane(o, 0.0, 0.0, 1.0).unwrap();
        s.lock_plane(plane).unwrap();
        let p = s.add_point(3.0, -1.0, 7.0);
        s.add_constraint(Constraint3D::PointInPlane { point: p, plane });
        let rep = s.solve().expect("ok");
        assert_eq!(rep.status, SolverStatus::Converged);
        // The plane normal is unchanged by the solve.
        let (_, _, _, nx, ny, nz) = s.plane_data(plane);
        assert!((nx).abs() < 1e-6 && (ny).abs() < 1e-6 && (nz - 1.0).abs() < 1e-6);
        // The origin point is unchanged too.
        let (ox, oy, oz) = s.point_xyz(o);
        assert!(ox.abs() < 1e-6 && oy.abs() < 1e-6 && oz.abs() < 1e-6);
    }

    /// PointDistance3 — pin a point at the origin (no constraint, just
    /// near-zero start) and pull a second point to distance 10.
    #[test]
    fn point_distance_converges() {
        let mut s = Sketch3D::new();
        let a = s.add_point(0.0, 0.0, 0.0);
        let b = s.add_point(1.0, 1.0, 1.0);
        s.add_constraint(Constraint3D::PointDistance3 {
            a,
            b,
            target: 10.0,
        });
        let rep = s.solve().expect("ok");
        assert_eq!(rep.status, SolverStatus::Converged);
        let (ax, ay, az) = s.point_xyz(a);
        let (bx, by, bz) = s.point_xyz(b);
        let d = ((bx - ax).powi(2) + (by - ay).powi(2) + (bz - az).powi(2)).sqrt();
        assert!((d - 10.0).abs() < 1e-5, "distance = {d}");
    }

    /// LineAngle3 — start two lines at ~30° apart and pull them to 90°.
    #[test]
    fn line_angle_converges() {
        let mut s = Sketch3D::new();
        let o = s.add_point(0.0, 0.0, 0.0);
        let a = s.add_point(1.0, 0.0, 0.0);
        let b = s.add_point(0.7, 0.3, 0.0);
        let la = s.add_line(o, a).unwrap();
        let lb = s.add_line(o, b).unwrap();
        s.add_constraint(Constraint3D::LineAngle3 {
            a: la,
            b: lb,
            target: std::f64::consts::FRAC_PI_2,
        });
        let _ = s.solve().expect("ok");
        let (ox, oy, oz) = s.point_xyz(o);
        let (ax2, ay2, az2) = s.point_xyz(a);
        let (bx2, by2, bz2) = s.point_xyz(b);
        let ua = (ax2 - ox, ay2 - oy, az2 - oz);
        let ub = (bx2 - ox, by2 - oy, bz2 - oz);
        let dot = ua.0 * ub.0 + ua.1 * ub.1 + ua.2 * ub.2;
        let na = (ua.0.powi(2) + ua.1.powi(2) + ua.2.powi(2)).sqrt();
        let nb = (ub.0.powi(2) + ub.1.powi(2) + ub.2.powi(2)).sqrt();
        let cos = dot / (na * nb);
        assert!(cos.abs() < 1e-4, "lines not perpendicular: cos = {cos}");
    }

    /// CircleRadius — a free-radius circle pulled to a target radius.
    #[test]
    fn circle_radius_converges() {
        let mut s = Sketch3D::new();
        let c = s.add_point(0.0, 0.0, 0.0);
        let circle = s.add_circle(c, 1.0, 0.0, 0.0, 1.0).unwrap();
        s.add_constraint(Constraint3D::CircleRadius { circle, target: 5.0 });
        let rep = s.solve().expect("ok");
        assert_eq!(rep.status, SolverStatus::Converged);
        assert!(
            (s.circle_radius(circle) - 5.0).abs() < 1e-6,
            "radius = {}",
            s.circle_radius(circle)
        );
    }

    /// PointOnCircle — a free point dropped onto a radius-5 circle ends up
    /// in the circle's plane and at the radius distance from its centre.
    #[test]
    fn point_lands_on_circle() {
        let mut s = Sketch3D::new();
        let center = s.add_point(0.0, 0.0, 0.0);
        let circle = s.add_circle(center, 5.0, 0.0, 0.0, 1.0).unwrap();
        s.add_constraint(Constraint3D::CircleRadius { circle, target: 5.0 });
        let p = s.add_point(2.0, 1.0, 3.0);
        s.add_constraint(Constraint3D::PointOnCircle { point: p, circle });
        let rep = s.solve().expect("ok");
        assert_eq!(rep.status, SolverStatus::Converged);
        let (cx, cy, cz, r, nx, ny, nz) = s.circle_data(circle);
        let (px, py, pz) = s.point_xyz(p);
        let (rx, ry, rz) = (px - cx, py - cy, pz - cz);
        let dist = (rx * rx + ry * ry + rz * rz).sqrt();
        let nlen = (nx * nx + ny * ny + nz * nz).sqrt();
        let plane_resid = (nx * rx + ny * ry + nz * rz) / nlen;
        assert!((dist - r).abs() < 1e-5, "point not at radius: dist={dist} r={r}");
        assert!(plane_resid.abs() < 1e-5, "point not in circle plane: {plane_resid}");
        assert!((r - 5.0).abs() < 1e-5, "radius drifted: {r}");
    }

    /// EqualRadius — two circles forced to share a radius.
    #[test]
    fn equal_radius_matches_two_circles() {
        let mut s = Sketch3D::new();
        let ca = s.add_point(0.0, 0.0, 0.0);
        let cb = s.add_point(10.0, 0.0, 0.0);
        let a = s.add_circle(ca, 3.0, 0.0, 0.0, 1.0).unwrap();
        let b = s.add_circle(cb, 1.0, 0.0, 0.0, 1.0).unwrap();
        s.add_constraint(Constraint3D::CircleRadius { circle: a, target: 3.0 });
        s.add_constraint(Constraint3D::EqualRadius { a, b });
        let rep = s.solve().expect("ok");
        assert_eq!(rep.status, SolverStatus::Converged);
        assert!(
            (s.circle_radius(b) - 3.0).abs() < 1e-6,
            "b radius = {}",
            s.circle_radius(b)
        );
    }

    /// The parametric-modelling loop: a named parameter (defined by an
    /// expression referencing another) drives a constraint target, which
    /// drives the geometry. Editing `base` would re-drive the radius.
    #[test]
    fn parameter_drives_a_constraint_target() {
        let mut params = ParameterTable::new();
        params.set("base", "4");
        params.set("radius", "base + 1"); // = 5
        let r = params.value("radius").unwrap();
        let mut s = Sketch3D::new();
        let c = s.add_point(0.0, 0.0, 0.0);
        let circle = s.add_circle(c, 1.0, 0.0, 0.0, 1.0).unwrap();
        s.add_constraint(Constraint3D::CircleRadius { circle, target: r });
        let rep = s.solve().expect("ok");
        assert_eq!(rep.status, SolverStatus::Converged);
        assert!(
            (s.circle_radius(circle) - 5.0).abs() < 1e-6,
            "radius = {}",
            s.circle_radius(circle)
        );
    }

    /// ArcRadius — a free-radius arc pulled to a target radius.
    #[test]
    fn arc_radius_converges() {
        let mut s = Sketch3D::new();
        let c = s.add_point(0.0, 0.0, 0.0);
        let start = s.add_point(5.0, 0.0, 0.0);
        let end = s.add_point(0.0, 5.0, 0.0);
        let arc = s.add_arc(c, 1.0, 0.0, 0.0, 1.0, start, end).unwrap();
        s.add_constraint(Constraint3D::ArcRadius { arc, target: 5.0 });
        let rep = s.solve().expect("ok");
        assert_eq!(rep.status, SolverStatus::Converged);
        assert!((s.arc_radius(arc) - 5.0).abs() < 1e-6, "radius = {}", s.arc_radius(arc));
    }

    /// ArcEndpointsOnArc — the start + end points are pulled onto the arc's
    /// circle: in its plane and at the radius distance from the centre.
    #[test]
    fn arc_endpoints_land_on_circle() {
        let mut s = Sketch3D::new();
        let c = s.add_point(0.0, 0.0, 0.0);
        let start = s.add_point(4.0, 1.0, 2.0);
        let end = s.add_point(1.0, 4.0, -1.0);
        let arc = s.add_arc(c, 5.0, 0.0, 0.0, 1.0, start, end).unwrap();
        s.add_constraint(Constraint3D::ArcRadius { arc, target: 5.0 });
        s.add_constraint(Constraint3D::ArcEndpointsOnArc { arc });
        let rep = s.solve().expect("ok");
        assert_eq!(rep.status, SolverStatus::Converged);
        let (cx, cy, cz, r, nx, ny, nz) = s.arc_circle(arc);
        let nlen = (nx * nx + ny * ny + nz * nz).sqrt();
        let (es, ee) = s.arc_endpoints(arc);
        for pt in [es, ee] {
            let (px, py, pz) = s.point_xyz(pt);
            let (rx, ry, rz) = (px - cx, py - cy, pz - cz);
            let dist = (rx * rx + ry * ry + rz * rz).sqrt();
            let plane = (nx * rx + ny * ry + nz * rz) / nlen;
            assert!((dist - r).abs() < 1e-5, "endpoint off radius: {dist} vs {r}");
            assert!(plane.abs() < 1e-5, "endpoint off plane: {plane}");
        }
    }

    /// Spline — a cubic Bézier evaluates to its endpoints at t=0/1, and its
    /// control points are solvable like any other points.
    #[test]
    fn spline_evaluates_and_tracks_its_control_points() {
        let mut s = Sketch3D::new();
        let p0 = s.add_point(0.0, 0.0, 0.0);
        let p1 = s.add_point(1.0, 2.0, 0.0);
        let p2 = s.add_point(2.0, 2.0, 0.0);
        let p3 = s.add_point(3.0, 0.0, 0.0);
        let spline = s.add_spline(p0, p1, p2, p3).unwrap();
        let start = s.spline_point_at(spline, 0.0);
        let end = s.spline_point_at(spline, 1.0);
        assert!(start.0.abs() < 1e-9 && start.1.abs() < 1e-9, "B(0) = p0");
        assert!((end.0 - 3.0).abs() < 1e-9 && end.1.abs() < 1e-9, "B(1) = p3");
        // The control points solve under ordinary constraints.
        s.add_constraint(Constraint3D::PointDistance3 { a: p0, b: p3, target: 10.0 });
        let rep = s.solve().expect("ok");
        assert_eq!(rep.status, SolverStatus::Converged);
    }
}
