//! [`Sketch`] holds the parameter vector + the entity table + the
//! constraint set. The solver mutates `vars` in place; entities and
//! constraints index into `vars` by [`crate::geom::EntityId`].

use std::collections::HashSet;

use serde::{Deserialize, Serialize};

use crate::constraint::Constraint;
use crate::external_geom::{ExternalGeomRef, FeatureTreeLookup};
use crate::geom::{Arc2, Circle2, Entity, EntityId, Line2, Point2};
use crate::geom_bspline::BSpline2;
use crate::geom_ellipse::{Ellipse2, EllipticalArc2};

/// A 2-D parametric sketch on a single plane.
///
/// `PartialEq` is derived structurally so the host app can wrap a
/// `Sketch` in an `undo::History<Sketch>` for the Sketcher panel's
/// `↶ ↷` controls. The float-valued `vars` use IEEE 754 semantics:
/// a NaN coordinate never compares equal to itself, so a snapshot
/// containing NaN coordinates fails to dedupe in the undo stack.
/// The sketcher's `add_point` / drag-value widgets never let a NaN
/// reach `vars` in practice — this is a documented degradation
/// rather than a correctness issue.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct Sketch {
    /// Flat parameter vector that all entities index into.
    pub vars: Vec<f64>,
    /// Entity table — addressed by 1-based [`EntityId`].
    pub entities: Vec<Entity>,
    /// Constraint list applied to entities.
    pub constraints: Vec<Constraint>,
    /// Construction-geometry flag per entity (Phase 12C). Parallel to
    /// `entities`; default `false` (regular geometry).
    #[serde(default)]
    pub construction: Vec<bool>,
    /// External-geometry source ref per entity (Phase 12D). `Some(_)`
    /// means the entity was projected from another feature and the
    /// solver should freeze its variables.
    #[serde(default)]
    pub external: Vec<Option<ExternalGeomRef>>,
    /// Set of variable indices the solver must keep frozen. Populated
    /// by [`Sketch::add_external_edge`] (and friends) — these are the
    /// x/y indices of external-geometry points.
    #[serde(default)]
    pub frozen_vars: HashSet<usize>,
}

impl Sketch {
    /// Empty sketch.
    pub fn new() -> Self {
        Self::default()
    }

    /// Allocate a new variable slot initialised to `value`. Returns
    /// the variable's index.
    pub fn alloc_var(&mut self, value: f64) -> usize {
        self.vars.push(value);
        self.vars.len() - 1
    }

    /// Push an entity into the table along with the parallel
    /// construction/external metadata.
    fn push_entity(&mut self, e: Entity) -> EntityId {
        self.entities.push(e);
        self.construction.push(false);
        self.external.push(None);
        EntityId(self.entities.len())
    }

    /// Add a point at (x, y). Returns its [`EntityId`].
    pub fn add_point(&mut self, x: f64, y: f64) -> EntityId {
        let x_var = self.alloc_var(x);
        let y_var = self.alloc_var(y);
        self.push_entity(Entity::Point(Point2 { x_var, y_var }))
    }

    /// Add a line between the given existing points.
    pub fn add_line(
        &mut self,
        start: EntityId,
        end: EntityId,
    ) -> Result<EntityId, crate::SketchError> {
        let s = self.point_at(start)?;
        let e = self.point_at(end)?;
        Ok(self.push_entity(Entity::Line(Line2 { start: s, end: e })))
    }

    /// Add a circle at the given centre with the given radius.
    pub fn add_circle(
        &mut self,
        center: EntityId,
        radius: f64,
    ) -> Result<EntityId, crate::SketchError> {
        let c = self.point_at(center)?;
        let r_var = self.alloc_var(radius);
        Ok(self.push_entity(Entity::Circle(Circle2 {
            center: c,
            radius_var: r_var,
        })))
    }

    /// Add an arc.
    pub fn add_arc(
        &mut self,
        center: EntityId,
        radius: f64,
        start_angle: f64,
        end_angle: f64,
    ) -> Result<EntityId, crate::SketchError> {
        let c = self.point_at(center)?;
        let r_var = self.alloc_var(radius);
        let sa = self.alloc_var(start_angle);
        let ea = self.alloc_var(end_angle);
        Ok(self.push_entity(Entity::Arc(Arc2 {
            center: c,
            radius_var: r_var,
            start_angle_var: sa,
            end_angle_var: ea,
        })))
    }

    /// Add a B-spline curve (Phase 12A). `control_points` are existing
    /// point entity ids — variables are reused.
    pub fn add_bspline(
        &mut self,
        degree: usize,
        knots: Vec<f64>,
        control_points: &[EntityId],
        weights: Vec<f64>,
    ) -> Result<EntityId, crate::SketchError> {
        let cps: Vec<Point2> = control_points
            .iter()
            .map(|id| self.point_at(*id))
            .collect::<Result<_, _>>()?;
        Ok(self.push_entity(Entity::BSpline(BSpline2 {
            degree,
            knots,
            control_points: cps,
            weights,
        })))
    }

    /// Add a full ellipse (Phase 12A). `center` is an existing point
    /// id; `major_axis` and `minor_radius` allocate new variables.
    pub fn add_ellipse(
        &mut self,
        center: EntityId,
        major_axis: (f64, f64),
        minor_radius: f64,
    ) -> Result<EntityId, crate::SketchError> {
        let c = self.point_at(center)?;
        let mx = self.alloc_var(major_axis.0);
        let my = self.alloc_var(major_axis.1);
        let mr = self.alloc_var(minor_radius);
        Ok(self.push_entity(Entity::Ellipse(Ellipse2 {
            center: c,
            major_x_var: mx,
            major_y_var: my,
            minor_radius_var: mr,
        })))
    }

    /// Add an elliptical arc (Phase 12A). Adds the underlying ellipse
    /// inline (centre is a point, major/minor/angles get fresh vars).
    pub fn add_elliptical_arc(
        &mut self,
        center: EntityId,
        major_axis: (f64, f64),
        minor_radius: f64,
        start_angle: f64,
        end_angle: f64,
    ) -> Result<EntityId, crate::SketchError> {
        let c = self.point_at(center)?;
        let mx = self.alloc_var(major_axis.0);
        let my = self.alloc_var(major_axis.1);
        let mr = self.alloc_var(minor_radius);
        let sa = self.alloc_var(start_angle);
        let ea = self.alloc_var(end_angle);
        let ellipse = Ellipse2 {
            center: c,
            major_x_var: mx,
            major_y_var: my,
            minor_radius_var: mr,
        };
        Ok(self.push_entity(Entity::EllipticalArc(EllipticalArc2 {
            ellipse,
            start_angle_var: sa,
            end_angle_var: ea,
        })))
    }

    /// Append a constraint. Does NOT validate that referenced entities
    /// exist or are of the right kind — validation happens during
    /// [`crate::solver::solve`].
    pub fn add_constraint(&mut self, c: Constraint) {
        self.constraints.push(c);
    }

    /// Total scalar residual equations across all constraints.
    pub fn total_residuals(&self) -> usize {
        self.constraints.iter().map(|c| c.n_residuals()).sum()
    }

    /// Look up an entity by id. Returns the underlying [`Point2`] if it
    /// is a point, else `ConstraintTypeMismatch`.
    pub fn point_at(&self, id: EntityId) -> Result<Point2, crate::SketchError> {
        let entity = self
            .entities
            .get(id.0.wrapping_sub(1))
            .ok_or(crate::SketchError::UnknownEntity(id.0))?;
        match entity {
            Entity::Point(p) => Ok(*p),
            other => Err(crate::SketchError::ConstraintTypeMismatch(format!(
                "expected Point, got {other:?}"
            ))),
        }
    }

    /// Look up a line by id.
    pub fn line_at(&self, id: EntityId) -> Result<Line2, crate::SketchError> {
        match self.entities.get(id.0.wrapping_sub(1)) {
            Some(Entity::Line(l)) => Ok(*l),
            Some(other) => Err(crate::SketchError::ConstraintTypeMismatch(format!(
                "expected Line, got {other:?}"
            ))),
            None => Err(crate::SketchError::UnknownEntity(id.0)),
        }
    }

    /// Look up a circle by id.
    pub fn circle_at(&self, id: EntityId) -> Result<Circle2, crate::SketchError> {
        match self.entities.get(id.0.wrapping_sub(1)) {
            Some(Entity::Circle(c)) => Ok(*c),
            Some(other) => Err(crate::SketchError::ConstraintTypeMismatch(format!(
                "expected Circle, got {other:?}"
            ))),
            None => Err(crate::SketchError::UnknownEntity(id.0)),
        }
    }

    /// Look up an arc by id.
    pub fn arc_at(&self, id: EntityId) -> Result<Arc2, crate::SketchError> {
        match self.entities.get(id.0.wrapping_sub(1)) {
            Some(Entity::Arc(a)) => Ok(*a),
            Some(other) => Err(crate::SketchError::ConstraintTypeMismatch(format!(
                "expected Arc, got {other:?}"
            ))),
            None => Err(crate::SketchError::UnknownEntity(id.0)),
        }
    }

    /// Look up a B-spline by id.
    pub fn bspline_at(&self, id: EntityId) -> Result<&BSpline2, crate::SketchError> {
        match self.entities.get(id.0.wrapping_sub(1)) {
            Some(Entity::BSpline(b)) => Ok(b),
            Some(other) => Err(crate::SketchError::ConstraintTypeMismatch(format!(
                "expected BSpline, got {other:?}"
            ))),
            None => Err(crate::SketchError::UnknownEntity(id.0)),
        }
    }

    /// Look up an ellipse by id.
    pub fn ellipse_at(&self, id: EntityId) -> Result<Ellipse2, crate::SketchError> {
        match self.entities.get(id.0.wrapping_sub(1)) {
            Some(Entity::Ellipse(e)) => Ok(*e),
            Some(other) => Err(crate::SketchError::ConstraintTypeMismatch(format!(
                "expected Ellipse, got {other:?}"
            ))),
            None => Err(crate::SketchError::UnknownEntity(id.0)),
        }
    }

    /// Look up an elliptical arc by id.
    pub fn elliptical_arc_at(&self, id: EntityId) -> Result<EllipticalArc2, crate::SketchError> {
        match self.entities.get(id.0.wrapping_sub(1)) {
            Some(Entity::EllipticalArc(a)) => Ok(*a),
            Some(other) => Err(crate::SketchError::ConstraintTypeMismatch(format!(
                "expected EllipticalArc, got {other:?}"
            ))),
            None => Err(crate::SketchError::UnknownEntity(id.0)),
        }
    }

    /// Phase 12C: toggle the construction flag for an entity. Out-of-
    /// range ids are silently ignored.
    pub fn toggle_construction(&mut self, id: EntityId) {
        if let Some(slot) = self.construction.get_mut(id.0.wrapping_sub(1)) {
            *slot = !*slot;
        }
    }

    /// Phase 12C: mark `id` as construction (set flag to true).
    pub fn mark_construction(&mut self, id: EntityId, on: bool) {
        if let Some(slot) = self.construction.get_mut(id.0.wrapping_sub(1)) {
            *slot = on;
        }
    }

    /// Phase 12C: returns true if `id` is flagged construction.
    pub fn is_construction(&self, id: EntityId) -> bool {
        self.construction
            .get(id.0.wrapping_sub(1))
            .copied()
            .unwrap_or(false)
    }

    /// Phase 12D: add an external edge reference. Project to the
    /// sketch plane via `lookup`, create a fixed Line with two new
    /// frozen points at the projected endpoints, and record the source
    /// ref on the new line entity.
    pub fn add_external_edge(
        &mut self,
        ext: ExternalGeomRef,
        lookup: &dyn FeatureTreeLookup,
    ) -> Result<EntityId, crate::SketchError> {
        let ((sx, sy), (ex, ey)) =
            lookup
                .edge_endpoints(&ext)
                .ok_or(crate::SketchError::ConstraintTypeMismatch(format!(
                    "external geom {ext:?} returned no endpoints"
                )))?;
        let a = self.add_point(sx, sy);
        let b = self.add_point(ex, ey);
        // Freeze the new points' variables.
        let pa = self.point_at(a).unwrap();
        let pb = self.point_at(b).unwrap();
        self.frozen_vars.insert(pa.x_var);
        self.frozen_vars.insert(pa.y_var);
        self.frozen_vars.insert(pb.x_var);
        self.frozen_vars.insert(pb.y_var);
        let line = self.add_line(a, b)?;
        // Mark the line itself as external (the endpoints are
        // implicitly external by being frozen).
        let idx = line.0 - 1;
        self.external[idx] = Some(ext);
        Ok(line)
    }

    /// Phase 12D: re-fetch external sources from `lookup` and update
    /// the frozen variable values to match. Call after the source
    /// features have been edited.
    pub fn resolve_externals(&mut self, lookup: &dyn FeatureTreeLookup) {
        for idx in 0..self.entities.len() {
            let Some(ext) = self.external[idx] else {
                continue;
            };
            if let Entity::Line(l) = self.entities[idx] {
                if let Some(((sx, sy), (ex, ey))) = lookup.edge_endpoints(&ext) {
                    self.vars[l.start.x_var] = sx;
                    self.vars[l.start.y_var] = sy;
                    self.vars[l.end.x_var] = ex;
                    self.vars[l.end.y_var] = ey;
                }
            }
        }
    }

    /// Returns true if variable `var` is frozen by an external-geom
    /// reference and must be skipped by the solver.
    pub fn is_var_frozen(&self, var: usize) -> bool {
        self.frozen_vars.contains(&var)
    }

    /// Convenience wrapper: extract the closed profile and extrude
    /// it along +Z by `depth`. Defaults to a 1e-6 profile tolerance.
    pub fn extrude(&self, depth: f64) -> Result<valenx_cad::Solid, crate::SketchError> {
        crate::extrude::extrude(self, depth, 1e-6)
    }

    /// Validate that every entity's variable handles index inside
    /// [`Self::vars`]. R33 H1.
    ///
    /// `vars` and `entities` are independent public fields, so a
    /// hand-edited, corrupt, or version-skewed `.valenx` (RON) can
    /// deserialize into a `Sketch` whose entity carries a var-handle
    /// past the end of `vars`. Consuming that sketch — feature-tree
    /// replay (`hole`, `pad`, `revolve`, …), the solver's Jacobian
    /// assembly, profile extraction — would otherwise panic with
    /// "index out of bounds". This walks every entity variant and every
    /// handle field it carries, returning
    /// [`crate::SketchError::CorruptHandle`] on the first violation so a
    /// bad document is rejected at load instead of crashing the app.
    ///
    /// Called by [`crate::persist::SketchFile::from_ron`] (standalone
    /// sketch load) and
    /// [`crate::persist`]'s project envelope counterpart in
    /// `valenx-feature-tree` (project load validates every sketch).
    pub fn validate(&self) -> Result<(), crate::SketchError> {
        let len = self.vars.len();
        let check = |entity: usize,
                     field: &'static str,
                     var: usize|
         -> Result<(), crate::SketchError> {
            if var >= len {
                Err(crate::SketchError::CorruptHandle {
                    entity,
                    field,
                    var,
                    len,
                })
            } else {
                Ok(())
            }
        };
        // Helper: every Point2 carries an (x_var, y_var) pair. The
        // `kind` prefix disambiguates nested points (line start/end,
        // ellipse centre, etc.) in the error message.
        let check_point = |entity: usize,
                           x_field: &'static str,
                           y_field: &'static str,
                           p: &Point2|
         -> Result<(), crate::SketchError> {
            check(entity, x_field, p.x_var)?;
            check(entity, y_field, p.y_var)?;
            Ok(())
        };

        for (idx, entity) in self.entities.iter().enumerate() {
            match entity {
                Entity::Point(p) => {
                    check_point(idx, "x_var", "y_var", p)?;
                }
                Entity::Line(l) => {
                    check_point(idx, "start.x_var", "start.y_var", &l.start)?;
                    check_point(idx, "end.x_var", "end.y_var", &l.end)?;
                }
                Entity::Circle(c) => {
                    check_point(idx, "center.x_var", "center.y_var", &c.center)?;
                    check(idx, "radius_var", c.radius_var)?;
                }
                Entity::Arc(a) => {
                    check_point(idx, "center.x_var", "center.y_var", &a.center)?;
                    check(idx, "radius_var", a.radius_var)?;
                    check(idx, "start_angle_var", a.start_angle_var)?;
                    check(idx, "end_angle_var", a.end_angle_var)?;
                }
                Entity::BSpline(b) => {
                    // Structural invariants. A corrupt/hand-edited document can
                    // violate these, and curve evaluation indexes `knots[n_cp]`
                    // / `weights[cp_idx]` and computes `span - degree` without
                    // re-checking (see geom_bspline), so a malformed curve must
                    // be rejected here rather than panic during replay.
                    let n_cp = b.control_points.len();
                    let bad =
                        |reason: String| crate::SketchError::CorruptBSpline { entity: idx, reason };
                    if b.degree < 1 {
                        return Err(bad(format!("degree {} must be >= 1", b.degree)));
                    }
                    if n_cp < b.degree + 1 {
                        return Err(bad(format!(
                            "needs at least degree+1 = {} control points, has {n_cp}",
                            b.degree + 1
                        )));
                    }
                    if b.weights.len() != n_cp {
                        return Err(bad(format!(
                            "weights length {} != control-point count {n_cp}",
                            b.weights.len()
                        )));
                    }
                    let expected_knots = n_cp + b.degree + 1;
                    if b.knots.len() != expected_knots {
                        return Err(bad(format!(
                            "knot length {} != control_points + degree + 1 = {expected_knots}",
                            b.knots.len()
                        )));
                    }
                    if b.knots.windows(2).any(|w| w[1] < w[0]) {
                        return Err(bad("knot vector must be non-decreasing".to_string()));
                    }
                    for cp in &b.control_points {
                        check_point(idx, "control_point.x_var", "control_point.y_var", cp)?;
                    }
                }
                Entity::Ellipse(e) => {
                    check_point(idx, "center.x_var", "center.y_var", &e.center)?;
                    check(idx, "major_x_var", e.major_x_var)?;
                    check(idx, "major_y_var", e.major_y_var)?;
                    check(idx, "minor_radius_var", e.minor_radius_var)?;
                }
                Entity::EllipticalArc(a) => {
                    check_point(idx, "ellipse.center.x_var", "ellipse.center.y_var", &a.ellipse.center)?;
                    check(idx, "ellipse.major_x_var", a.ellipse.major_x_var)?;
                    check(idx, "ellipse.major_y_var", a.ellipse.major_y_var)?;
                    check(idx, "ellipse.minor_radius_var", a.ellipse.minor_radius_var)?;
                    check(idx, "start_angle_var", a.start_angle_var)?;
                    check(idx, "end_angle_var", a.end_angle_var)?;
                }
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empty_sketch_has_no_vars_no_entities() {
        let s = Sketch::new();
        assert!(s.vars.is_empty());
        assert!(s.entities.is_empty());
    }

    #[test]
    fn add_point_allocates_two_vars_returns_id_1() {
        let mut s = Sketch::new();
        let id = s.add_point(1.0, 2.0);
        assert_eq!(id, EntityId(1));
        assert_eq!(s.vars, vec![1.0, 2.0]);
        assert_eq!(s.entities.len(), 1);
    }

    #[test]
    fn add_line_between_two_points_succeeds() {
        let mut s = Sketch::new();
        let a = s.add_point(0.0, 0.0);
        let b = s.add_point(3.0, 4.0);
        let line = s.add_line(a, b).unwrap();
        assert_eq!(line, EntityId(3));
    }

    #[test]
    fn add_line_with_circle_endpoint_errors() {
        let mut s = Sketch::new();
        let a = s.add_point(0.0, 0.0);
        let c_center = s.add_point(1.0, 0.0);
        let circle = s.add_circle(c_center, 0.5).unwrap();
        let err = s.add_line(a, circle).unwrap_err();
        assert_eq!(err.code(), "sketch.constraint_type_mismatch");
    }

    #[test]
    fn add_circle_allocates_one_extra_var() {
        let mut s = Sketch::new();
        let c = s.add_point(0.0, 0.0);
        let circle = s.add_circle(c, 5.0).unwrap();
        assert_eq!(circle, EntityId(2));
        assert_eq!(s.vars, vec![0.0, 0.0, 5.0]);
    }

    #[test]
    fn line_at_returns_line() {
        let mut s = Sketch::new();
        let a = s.add_point(0.0, 0.0);
        let b = s.add_point(1.0, 0.0);
        let l = s.add_line(a, b).unwrap();
        let line = s.line_at(l).unwrap();
        assert_eq!(line.length(&s.vars), 1.0);
    }

    #[test]
    fn circle_at_returns_circle() {
        let mut s = Sketch::new();
        let c = s.add_point(0.0, 0.0);
        let cid = s.add_circle(c, 2.5).unwrap();
        assert_eq!(s.circle_at(cid).unwrap().radius(&s.vars), 2.5);
    }

    #[test]
    fn wrong_kind_lookup_errors() {
        let mut s = Sketch::new();
        let p = s.add_point(0.0, 0.0);
        let err = s.line_at(p).unwrap_err();
        assert_eq!(err.code(), "sketch.constraint_type_mismatch");
    }

    use crate::constraint::Constraint;

    #[test]
    fn add_constraint_appends_to_list() {
        let mut s = Sketch::new();
        let a = s.add_point(0.0, 0.0);
        let b = s.add_point(1.0, 1.0);
        s.add_constraint(Constraint::Coincident { a, b });
        assert_eq!(s.constraints.len(), 1);
    }

    #[test]
    fn total_residuals_sums_per_constraint() {
        let mut s = Sketch::new();
        let a = s.add_point(0.0, 0.0);
        let b = s.add_point(1.0, 0.0);
        let line = s.add_line(a, b).unwrap();
        s.add_constraint(Constraint::Coincident { a, b }); // 2 residuals
        s.add_constraint(Constraint::Horizontal(line)); // 1 residual
        assert_eq!(s.total_residuals(), 3);
    }

    #[test]
    fn sketch_extrude_method_exists() {
        let mut s = Sketch::new();
        let a = s.add_point(0.0, 0.0);
        let b = s.add_point(1.0, 0.0);
        let c = s.add_point(0.5, 1.0);
        s.add_line(a, b).unwrap();
        s.add_line(b, c).unwrap();
        s.add_line(c, a).unwrap();
        // Either succeeds (truck API works) or fails with a documented error
        // — both are acceptable, the method exists.
        let _ = s.extrude(1.0);
    }

    // Phase 12A — new primitives

    #[test]
    fn add_bspline_returns_id_and_lookup_works() {
        let mut s = Sketch::new();
        let p0 = s.add_point(0.0, 0.0);
        let p1 = s.add_point(1.0, 1.0);
        let p2 = s.add_point(2.0, 1.0);
        let p3 = s.add_point(3.0, 0.0);
        let id = s
            .add_bspline(
                3,
                vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
                &[p0, p1, p2, p3],
                vec![1.0; 4],
            )
            .unwrap();
        let b = s.bspline_at(id).unwrap();
        assert_eq!(b.degree, 3);
        assert_eq!(b.n_control_points(), 4);
    }

    #[test]
    fn add_ellipse_returns_id_and_lookup_works() {
        let mut s = Sketch::new();
        let c = s.add_point(0.0, 0.0);
        let id = s.add_ellipse(c, (2.0, 0.0), 1.0).unwrap();
        let e = s.ellipse_at(id).unwrap();
        assert!((e.major_radius(&s.vars) - 2.0).abs() < 1e-12);
        assert!((e.minor_radius(&s.vars) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn add_elliptical_arc_returns_id_and_lookup_works() {
        let mut s = Sketch::new();
        let c = s.add_point(0.0, 0.0);
        let id = s
            .add_elliptical_arc(c, (2.0, 0.0), 1.0, 0.0, std::f64::consts::PI)
            .unwrap();
        let arc = s.elliptical_arc_at(id).unwrap();
        assert!((arc.sweep(&s.vars) - std::f64::consts::PI).abs() < 1e-12);
    }

    // R33 H1 — validated-deserialize: a hand-edited / version-skewed
    // `.valenx` (RON) whose entity carries a var-handle past the end of
    // `vars` must be rejected at load via `validate()`, not panic with
    // "index out of bounds" during feature-tree replay.

    #[test]
    fn validate_rejects_point_with_out_of_range_x_var() {
        // One var slot, but a Point referencing index 999.
        let s = Sketch {
            vars: vec![0.0],
            entities: vec![Entity::Point(Point2 {
                x_var: 999,
                y_var: 0,
            })],
            ..Sketch::default()
        };
        let err = s.validate().expect_err("out-of-range x_var must be rejected");
        assert_eq!(err.code(), "sketch.corrupt_handle");
    }

    #[test]
    fn validate_rejects_point_with_out_of_range_y_var() {
        let s = Sketch {
            vars: vec![0.0],
            entities: vec![Entity::Point(Point2 {
                x_var: 0,
                y_var: 999,
            })],
            ..Sketch::default()
        };
        assert!(s.validate().is_err());
    }

    #[test]
    fn validate_rejects_line_endpoint_handle_out_of_range() {
        let s = Sketch {
            vars: vec![0.0, 0.0],
            entities: vec![Entity::Line(Line2 {
                start: Point2 { x_var: 0, y_var: 1 },
                end: Point2 {
                    x_var: 2, // == len → out of range
                    y_var: 1,
                },
            })],
            ..Sketch::default()
        };
        assert!(s.validate().is_err());
    }

    #[test]
    fn validate_rejects_circle_radius_handle_out_of_range() {
        let s = Sketch {
            vars: vec![0.0, 0.0],
            entities: vec![Entity::Circle(Circle2 {
                center: Point2 { x_var: 0, y_var: 1 },
                radius_var: 7,
            })],
            ..Sketch::default()
        };
        assert!(s.validate().is_err());
    }

    #[test]
    fn validate_rejects_arc_angle_handle_out_of_range() {
        let s = Sketch {
            vars: vec![0.0, 0.0, 1.0, 0.0],
            entities: vec![Entity::Arc(Arc2 {
                center: Point2 { x_var: 0, y_var: 1 },
                radius_var: 2,
                start_angle_var: 3,
                end_angle_var: 42, // out of range
            })],
            ..Sketch::default()
        };
        assert!(s.validate().is_err());
    }

    #[test]
    fn validate_rejects_ellipse_major_handle_out_of_range() {
        use crate::geom_ellipse::Ellipse2;
        let s = Sketch {
            vars: vec![0.0, 0.0, 1.0, 0.0],
            entities: vec![Entity::Ellipse(Ellipse2 {
                center: Point2 { x_var: 0, y_var: 1 },
                major_x_var: 2,
                major_y_var: 99, // out of range
                minor_radius_var: 3,
            })],
            ..Sketch::default()
        };
        assert!(s.validate().is_err());
    }

    #[test]
    fn validate_rejects_elliptical_arc_angle_handle_out_of_range() {
        use crate::geom_ellipse::{Ellipse2, EllipticalArc2};
        let s = Sketch {
            vars: vec![0.0, 0.0, 1.0, 1.0, 0.0],
            entities: vec![Entity::EllipticalArc(EllipticalArc2 {
                ellipse: Ellipse2 {
                    center: Point2 { x_var: 0, y_var: 1 },
                    major_x_var: 2,
                    major_y_var: 3,
                    minor_radius_var: 4,
                },
                start_angle_var: 0,
                end_angle_var: 500, // out of range
            })],
            ..Sketch::default()
        };
        assert!(s.validate().is_err());
    }

    #[test]
    fn validate_rejects_bspline_control_point_handle_out_of_range() {
        use crate::geom_bspline::BSpline2;
        let s = Sketch {
            vars: vec![0.0, 0.0],
            entities: vec![Entity::BSpline(BSpline2 {
                degree: 1,
                knots: vec![0.0, 0.0, 1.0, 1.0],
                control_points: vec![
                    Point2 { x_var: 0, y_var: 1 },
                    Point2 {
                        x_var: 2, // out of range
                        y_var: 1,
                    },
                ],
                weights: vec![1.0, 1.0],
            })],
            ..Sketch::default()
        };
        assert!(s.validate().is_err());
    }

    #[test]
    fn validate_rejects_bspline_with_mismatched_knot_vector() {
        use crate::geom_bspline::BSpline2;
        // Valid control-point handles, but a degree-3 / 4-CP curve needs 8
        // knots and here has none — the kind of corruption a hand-edited or
        // version-skewed document carries. Pre-fix, evaluation indexed
        // `knots[4]` on an empty vec and panicked during replay.
        let s = Sketch {
            vars: vec![0.0, 0.0, 1.0, 1.0, 2.0, 1.0, 3.0, 0.0],
            entities: vec![Entity::BSpline(BSpline2 {
                degree: 3,
                knots: vec![],
                control_points: vec![
                    Point2 { x_var: 0, y_var: 1 },
                    Point2 { x_var: 2, y_var: 3 },
                    Point2 { x_var: 4, y_var: 5 },
                    Point2 { x_var: 6, y_var: 7 },
                ],
                weights: vec![1.0; 4],
            })],
            ..Sketch::default()
        };
        match s.validate() {
            Err(crate::SketchError::CorruptBSpline { reason, .. }) => {
                assert!(reason.contains("knot"), "reason names the knot mismatch: {reason}");
            }
            other => panic!("expected CorruptBSpline, got {other:?}"),
        }
    }

    #[test]
    fn validate_rejects_bspline_with_zero_degree() {
        use crate::geom_bspline::BSpline2;
        let s = Sketch {
            vars: vec![0.0, 0.0],
            entities: vec![Entity::BSpline(BSpline2 {
                degree: 0,
                knots: vec![0.0, 1.0],
                control_points: vec![Point2 { x_var: 0, y_var: 1 }],
                weights: vec![1.0],
            })],
            ..Sketch::default()
        };
        assert!(matches!(
            s.validate(),
            Err(crate::SketchError::CorruptBSpline { .. })
        ));
    }

    #[test]
    fn validate_accepts_a_well_formed_sketch_with_every_primitive() {
        let mut s = Sketch::new();
        let a = s.add_point(0.0, 0.0);
        let b = s.add_point(1.0, 0.0);
        let c = s.add_point(0.5, 1.0);
        s.add_line(a, b).unwrap();
        s.add_line(b, c).unwrap();
        s.add_circle(a, 0.5).unwrap();
        s.add_arc(a, 0.5, 0.0, std::f64::consts::PI).unwrap();
        s.add_ellipse(a, (2.0, 0.0), 1.0).unwrap();
        s.add_elliptical_arc(a, (2.0, 0.0), 1.0, 0.0, std::f64::consts::PI)
            .unwrap();
        s.add_bspline(
            1,
            vec![0.0, 0.0, 1.0, 1.0],
            &[a, b],
            vec![1.0, 1.0],
        )
        .unwrap();
        assert!(
            s.validate().is_ok(),
            "a sketch built through the normal API must always validate"
        );
    }
}
