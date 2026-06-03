//! Geometric constraints applied to sketch entities.
//!
//! Each [`Constraint`] variant defines one or more *residual functions*
//! `r(vars) -> f64` that the solver drives to zero. The variant also
//! defines a *Jacobian row* `dr/d(vars)` so the solver can step.

use serde::{Deserialize, Serialize};

use crate::geom::EntityId;

pub mod angle;
pub mod arc_angle;
pub mod arc_radius;
pub mod block;
pub mod coincident;
pub mod distance;
pub mod distance_x;
pub mod distance_y;
pub mod ellipse_radius_a;
pub mod ellipse_radius_b;
pub mod equal_length;
pub mod equal_length_scaled;
pub mod horizontal;
pub mod internal_alignment;
pub mod line_length;
pub mod parallel;
pub mod perpendicular;
pub mod point_on_arc;
pub mod point_on_bspline;
pub mod point_on_circle;
pub mod point_on_ellipse;
pub mod point_on_line;
pub mod radius;
pub mod snells_law;
pub mod symmetric;
pub mod tangent;
pub mod vertical;

pub use internal_alignment::AlignmentKind;

/// One geometric constraint between sketch entities.
///
/// `PartialEq` is derived structurally — variant kind plus every
/// nested field (including float targets / pitches / refractive
/// indices). Float fields use IEEE 754 semantics so a NaN target
/// never compares equal to itself; the practical effect is that an
/// undo snapshot containing NaN-valued constraints fails to dedupe.
/// The sketcher prevents NaN entry through its drag-value widgets.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum Constraint {
    /// Two points are coincident (occupy the same xy).
    Coincident {
        /// First point.
        a: EntityId,
        /// Second point — driven to share `a`'s location.
        b: EntityId,
    },
    /// A line is horizontal (slope = 0).
    Horizontal(EntityId),
    /// A line is vertical.
    Vertical(EntityId),
    /// Two lines are parallel.
    Parallel {
        /// First line.
        a: EntityId,
        /// Second line.
        b: EntityId,
    },
    /// Two lines are perpendicular.
    Perpendicular {
        /// First line.
        a: EntityId,
        /// Second line.
        b: EntityId,
    },
    /// A line is tangent to a circle (or two circles tangent).
    Tangent {
        /// Either a line or a circle.
        line_or_circle_a: EntityId,
        /// A circle that the first entity must be tangent to.
        circle_b: EntityId,
    },
    /// Two lines have equal length.
    EqualLength {
        /// First line.
        a: EntityId,
        /// Second line — driven to match `a`'s length.
        b: EntityId,
    },
    /// Two points are at a fixed Euclidean distance.
    Distance {
        /// First point.
        a: EntityId,
        /// Second point.
        b: EntityId,
        /// Target Euclidean distance between `a` and `b`.
        target: f64,
    },
    /// Two lines meet at a fixed angle (radians).
    Angle {
        /// First line.
        a: EntityId,
        /// Second line.
        b: EntityId,
        /// Target angle in radians.
        target: f64,
    },
    /// A circle / arc has a fixed radius.
    Radius {
        /// The circle or arc whose radius is constrained.
        circle_or_arc: EntityId,
        /// Target radius (sketch units).
        target: f64,
    },
    // ===== Phase 12B — additional constraints =====
    /// Two points symmetric about a third (the midpoint).
    Symmetric {
        /// First point.
        a: EntityId,
        /// Second point.
        b: EntityId,
        /// Midpoint that `a` and `b` are mirrored across.
        midpoint: EntityId,
    },
    /// Lock every coordinate variable of an entity to its current value.
    Block {
        /// `(var_idx, target_value)` pairs the solver must hold fixed.
        frozen: Vec<(usize, f64)>,
    },
    /// A point lies on a line's infinite extension.
    PointOnLine {
        /// Point that should lie on the line.
        point: EntityId,
        /// Line whose extension the point lies on.
        line: EntityId,
    },
    /// A point lies on a circle's perimeter.
    PointOnCircle {
        /// Point that lies on the circle.
        point: EntityId,
        /// Circle whose perimeter the point lies on.
        circle: EntityId,
    },
    /// A point lies on an arc's circular path.
    PointOnArc {
        /// Point on the arc.
        point: EntityId,
        /// Arc the point lies on.
        arc: EntityId,
    },
    /// A point lies on a B-spline curve.
    PointOnBSpline {
        /// Point projected onto the curve.
        point: EntityId,
        /// Target B-spline.
        bspline: EntityId,
    },
    /// A point lies on an ellipse perimeter.
    PointOnEllipse {
        /// Point on the ellipse.
        point: EntityId,
        /// Target ellipse.
        ellipse: EntityId,
    },
    /// A point is aligned with an ellipse's major/minor-axis endpoint.
    InternalAlignment {
        /// Point bound to the axis endpoint.
        point: EntityId,
        /// Parent ellipse.
        ellipse: EntityId,
        /// Which axis endpoint the point tracks.
        kind: AlignmentKind,
    },
    /// Horizontal signed distance between two points.
    DistanceX {
        /// First point.
        a: EntityId,
        /// Second point.
        b: EntityId,
        /// Target b.x - a.x.
        target: f64,
    },
    /// Vertical signed distance between two points.
    DistanceY {
        /// First point.
        a: EntityId,
        /// Second point.
        b: EntityId,
        /// Target b.y - a.y.
        target: f64,
    },
    /// A line has a fixed length (UI convenience for Distance on endpoints).
    LineLength {
        /// Line whose length is constrained.
        line: EntityId,
        /// Target length.
        target: f64,
    },
    /// An arc has a fixed radius.
    ArcRadius {
        /// Arc whose radius is constrained.
        arc: EntityId,
        /// Target radius.
        target: f64,
    },
    /// An arc sweeps a fixed angle.
    ArcAngle {
        /// Arc whose sweep is constrained.
        arc: EntityId,
        /// Target sweep (radians).
        target: f64,
    },
    /// Ellipse semi-major length = target.
    EllipseRadiusA {
        /// Ellipse whose semi-major is constrained.
        ellipse: EntityId,
        /// Target semi-major length.
        target: f64,
    },
    /// Ellipse semi-minor length = target.
    EllipseRadiusB {
        /// Ellipse whose semi-minor is constrained.
        ellipse: EntityId,
        /// Target semi-minor length.
        target: f64,
    },
    /// |a| = |b| * factor.
    EqualLengthScaled {
        /// First line.
        a: EntityId,
        /// Second line.
        b: EntityId,
        /// Scale factor: target |a| = factor * |b|.
        factor: f64,
    },
    /// Snell's law refraction at an interface line.
    SnellsLaw {
        /// Incoming ray (line).
        ray_in: EntityId,
        /// Outgoing (refracted) ray (line).
        ray_out: EntityId,
        /// Interface line (surface).
        interface_line: EntityId,
        /// Refractive index of the incoming medium.
        n1: f64,
        /// Refractive index of the outgoing medium.
        n2: f64,
    },
}

impl Constraint {
    /// Number of scalar residual equations this constraint contributes.
    /// Coincident = 2 (Δx + Δy), most others = 1.
    pub fn n_residuals(&self) -> usize {
        match self {
            Constraint::Coincident { .. } => 2,
            Constraint::Symmetric { .. } => 2,
            Constraint::InternalAlignment { .. } => 2,
            Constraint::Block { frozen } => frozen.len(),
            _ => 1,
        }
    }
}

use crate::sketch::Sketch;

impl Constraint {
    /// Compute residuals for this constraint into the provided slice.
    /// The slice length must equal [`Self::n_residuals`].
    pub fn residuals(&self, sketch: &Sketch, out: &mut [f64]) {
        assert_eq!(
            out.len(),
            self.n_residuals(),
            "residual slice length mismatch"
        );
        match self {
            Constraint::Coincident { a, b } => coincident::residuals(sketch, *a, *b, out),
            Constraint::Horizontal(line) => horizontal::residuals(sketch, *line, out),
            Constraint::Vertical(line) => vertical::residuals(sketch, *line, out),
            Constraint::Parallel { a, b } => parallel::residuals(sketch, *a, *b, out),
            Constraint::Perpendicular { a, b } => perpendicular::residuals(sketch, *a, *b, out),
            Constraint::Tangent {
                line_or_circle_a,
                circle_b,
            } => tangent::residuals(sketch, *line_or_circle_a, *circle_b, out),
            Constraint::EqualLength { a, b } => equal_length::residuals(sketch, *a, *b, out),
            Constraint::Distance { a, b, target } => {
                distance::residuals(sketch, *a, *b, *target, out)
            }
            Constraint::Angle { a, b, target } => angle::residuals(sketch, *a, *b, *target, out),
            Constraint::Radius {
                circle_or_arc,
                target,
            } => radius::residuals(sketch, *circle_or_arc, *target, out),
            Constraint::Symmetric { a, b, midpoint } => {
                symmetric::residuals(sketch, *a, *b, *midpoint, out)
            }
            Constraint::Block { frozen } => block::residuals(sketch, frozen, out),
            Constraint::PointOnLine { point, line } => {
                point_on_line::residuals(sketch, *point, *line, out)
            }
            Constraint::PointOnCircle { point, circle } => {
                point_on_circle::residuals(sketch, *point, *circle, out)
            }
            Constraint::PointOnArc { point, arc } => {
                point_on_arc::residuals(sketch, *point, *arc, out)
            }
            Constraint::PointOnBSpline { point, bspline } => {
                point_on_bspline::residuals(sketch, *point, *bspline, out)
            }
            Constraint::PointOnEllipse { point, ellipse } => {
                point_on_ellipse::residuals(sketch, *point, *ellipse, out)
            }
            Constraint::InternalAlignment {
                point,
                ellipse,
                kind,
            } => internal_alignment::residuals(sketch, *point, *ellipse, *kind, out),
            Constraint::DistanceX { a, b, target } => {
                distance_x::residuals(sketch, *a, *b, *target, out)
            }
            Constraint::DistanceY { a, b, target } => {
                distance_y::residuals(sketch, *a, *b, *target, out)
            }
            Constraint::LineLength { line, target } => {
                line_length::residuals(sketch, *line, *target, out)
            }
            Constraint::ArcRadius { arc, target } => {
                arc_radius::residuals(sketch, *arc, *target, out)
            }
            Constraint::ArcAngle { arc, target } => {
                arc_angle::residuals(sketch, *arc, *target, out)
            }
            Constraint::EllipseRadiusA { ellipse, target } => {
                ellipse_radius_a::residuals(sketch, *ellipse, *target, out)
            }
            Constraint::EllipseRadiusB { ellipse, target } => {
                ellipse_radius_b::residuals(sketch, *ellipse, *target, out)
            }
            Constraint::EqualLengthScaled { a, b, factor } => {
                equal_length_scaled::residuals(sketch, *a, *b, *factor, out)
            }
            Constraint::SnellsLaw {
                ray_in,
                ray_out,
                interface_line,
                n1,
                n2,
            } => snells_law::residuals(sketch, *ray_in, *ray_out, *interface_line, *n1, *n2, out),
        }
    }

    /// Write the non-zero Jacobian entries into `triplets` as
    /// `(row_offset, var_index, derivative)`. `row_offset` is the
    /// per-constraint row index in `0..n_residuals`.
    pub fn jacobian_triplets(&self, sketch: &Sketch, triplets: &mut Vec<(usize, usize, f64)>) {
        match self {
            Constraint::Coincident { a, b } => coincident::jacobian(sketch, *a, *b, triplets),
            Constraint::Horizontal(line) => horizontal::jacobian(sketch, *line, triplets),
            Constraint::Vertical(line) => vertical::jacobian(sketch, *line, triplets),
            Constraint::Parallel { a, b } => parallel::jacobian(sketch, *a, *b, triplets),
            Constraint::Perpendicular { a, b } => perpendicular::jacobian(sketch, *a, *b, triplets),
            Constraint::Tangent {
                line_or_circle_a,
                circle_b,
            } => tangent::jacobian(sketch, *line_or_circle_a, *circle_b, triplets),
            Constraint::EqualLength { a, b } => equal_length::jacobian(sketch, *a, *b, triplets),
            Constraint::Distance { a, b, target } => {
                distance::jacobian(sketch, *a, *b, *target, triplets)
            }
            Constraint::Angle { a, b, target } => {
                angle::jacobian(sketch, *a, *b, *target, triplets)
            }
            Constraint::Radius {
                circle_or_arc,
                target,
            } => radius::jacobian(sketch, *circle_or_arc, *target, triplets),
            Constraint::Symmetric { a, b, midpoint } => {
                symmetric::jacobian(sketch, *a, *b, *midpoint, triplets)
            }
            Constraint::Block { frozen } => block::jacobian(sketch, frozen, triplets),
            Constraint::PointOnLine { point, line } => {
                point_on_line::jacobian(sketch, *point, *line, triplets)
            }
            Constraint::PointOnCircle { point, circle } => {
                point_on_circle::jacobian(sketch, *point, *circle, triplets)
            }
            Constraint::PointOnArc { point, arc } => {
                point_on_arc::jacobian(sketch, *point, *arc, triplets)
            }
            Constraint::PointOnBSpline { point, bspline } => {
                point_on_bspline::jacobian(sketch, *point, *bspline, triplets)
            }
            Constraint::PointOnEllipse { point, ellipse } => {
                point_on_ellipse::jacobian(sketch, *point, *ellipse, triplets)
            }
            Constraint::InternalAlignment {
                point,
                ellipse,
                kind,
            } => internal_alignment::jacobian(sketch, *point, *ellipse, *kind, triplets),
            Constraint::DistanceX { a, b, target } => {
                distance_x::jacobian(sketch, *a, *b, *target, triplets)
            }
            Constraint::DistanceY { a, b, target } => {
                distance_y::jacobian(sketch, *a, *b, *target, triplets)
            }
            Constraint::LineLength { line, target } => {
                line_length::jacobian(sketch, *line, *target, triplets)
            }
            Constraint::ArcRadius { arc, target } => {
                arc_radius::jacobian(sketch, *arc, *target, triplets)
            }
            Constraint::ArcAngle { arc, target } => {
                arc_angle::jacobian(sketch, *arc, *target, triplets)
            }
            Constraint::EllipseRadiusA { ellipse, target } => {
                ellipse_radius_a::jacobian(sketch, *ellipse, *target, triplets)
            }
            Constraint::EllipseRadiusB { ellipse, target } => {
                ellipse_radius_b::jacobian(sketch, *ellipse, *target, triplets)
            }
            Constraint::EqualLengthScaled { a, b, factor } => {
                equal_length_scaled::jacobian(sketch, *a, *b, *factor, triplets)
            }
            Constraint::SnellsLaw {
                ray_in,
                ray_out,
                interface_line,
                n1,
                n2,
            } => snells_law::jacobian(
                sketch,
                *ray_in,
                *ray_out,
                *interface_line,
                *n1,
                *n2,
                triplets,
            ),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn coincident_contributes_two_residuals() {
        let c = Constraint::Coincident {
            a: EntityId(1),
            b: EntityId(2),
        };
        assert_eq!(c.n_residuals(), 2);
    }

    #[test]
    fn other_constraints_contribute_one_residual() {
        let constraints = vec![
            Constraint::Horizontal(EntityId(1)),
            Constraint::Vertical(EntityId(2)),
            Constraint::Parallel {
                a: EntityId(1),
                b: EntityId(2),
            },
            Constraint::Perpendicular {
                a: EntityId(1),
                b: EntityId(2),
            },
            Constraint::Tangent {
                line_or_circle_a: EntityId(1),
                circle_b: EntityId(2),
            },
            Constraint::EqualLength {
                a: EntityId(1),
                b: EntityId(2),
            },
            Constraint::Distance {
                a: EntityId(1),
                b: EntityId(2),
                target: 5.0,
            },
            Constraint::Angle {
                a: EntityId(1),
                b: EntityId(2),
                target: std::f64::consts::PI / 2.0,
            },
            Constraint::Radius {
                circle_or_arc: EntityId(1),
                target: 3.0,
            },
        ];
        for c in constraints {
            assert_eq!(c.n_residuals(), 1, "wrong count for {c:?}");
        }
    }

    /// Phase 12B Task 28: DOF balance regression.
    /// Sketch with one Symmetric constraint should report
    /// `n_residuals == 2`.
    #[test]
    fn dof_balance_with_symmetric() {
        let mut s = crate::sketch::Sketch::new();
        let a = s.add_point(0.0, 0.0);
        let b = s.add_point(2.0, 0.0);
        let c = s.add_point(1.0, 0.0);
        s.add_constraint(Constraint::Symmetric { a, b, midpoint: c });
        assert_eq!(s.total_residuals(), 2);
    }

    /// Phase 12B Task 28: Block contributes n_residuals = frozen.len().
    #[test]
    fn dof_balance_with_block() {
        let mut s = crate::sketch::Sketch::new();
        let _p = s.add_point(3.0, 4.0);
        s.add_constraint(Constraint::Block {
            frozen: vec![(0, 3.0), (1, 4.0)],
        });
        assert_eq!(s.total_residuals(), 2);
    }

    /// Phase 12B Task 28: each new single-residual constraint adds 1.
    #[test]
    fn dof_balance_with_distance_x_y() {
        let mut s = crate::sketch::Sketch::new();
        let a = s.add_point(0.0, 0.0);
        let b = s.add_point(5.0, 3.0);
        s.add_constraint(Constraint::DistanceX { a, b, target: 5.0 });
        s.add_constraint(Constraint::DistanceY { a, b, target: 3.0 });
        assert_eq!(s.total_residuals(), 2);
    }

    use crate::sketch::Sketch;

    #[test]
    fn residual_dispatch_returns_correct_count() {
        let mut s = Sketch::new();
        let a = s.add_point(0.0, 0.0);
        let b = s.add_point(1.0, 1.0);
        let c = Constraint::Coincident { a, b };
        let mut out = vec![0.0; 2];
        c.residuals(&s, &mut out);
        assert_eq!(out.len(), 2);
    }
}
