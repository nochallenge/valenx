//! [`Sketch3D`] — the 3D analogue of `valenx_sketch::Sketch`.
//!
//! Owns the global `vars: Vec<f64>` variable vector plus the heterogeneous
//! `entities: Vec<Entity3D>` table. Constraints reference entities by
//! [`EntityId`], and entities reference variables by index — exactly the
//! same structure-of-arrays pattern Phase 1 uses.

use serde::{Deserialize, Serialize};

use crate::constraint::Constraint3D;
use crate::entity::{Arc3, Circle3, Entity3D, EntityId, Line3, Plane3, Point3, Spline3, Workplane};
use crate::error::Solve3DError;

/// Variable-backed 3D sketch.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct Sketch3D {
    /// Global variable vector. Entities hold indices into it.
    pub vars: Vec<f64>,
    /// Entity table.
    pub entities: Vec<Entity3D>,
    /// Constraint list. Order matters: residuals are emitted in this
    /// order so callers can correlate solver-report rows back to source
    /// constraints.
    pub constraints: Vec<Constraint3D>,
}

impl Sketch3D {
    /// Empty sketch.
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a free 3D point at `(x, y, z)`. Returns its handle.
    pub fn add_point(&mut self, x: f64, y: f64, z: f64) -> EntityId {
        let xv = self.vars.len();
        self.vars.push(x);
        let yv = self.vars.len();
        self.vars.push(y);
        let zv = self.vars.len();
        self.vars.push(z);
        let id = EntityId(self.entities.len());
        self.entities.push(Entity3D::Point3(Point3 {
            x_var: xv,
            y_var: yv,
            z_var: zv,
        }));
        id
    }

    /// Add a line between two existing points.
    pub fn add_line(&mut self, a: EntityId, b: EntityId) -> Result<EntityId, Solve3DError> {
        self.expect_point(a)?;
        self.expect_point(b)?;
        let id = EntityId(self.entities.len());
        self.entities.push(Entity3D::Line3(Line3 { a, b }));
        Ok(id)
    }

    /// Add a plane through `origin` with normal `(nx, ny, nz)`. The
    /// normal is stored as three free variables; callers who want a
    /// fixed normal should add an `equality` constraint themselves.
    pub fn add_plane(
        &mut self,
        origin: EntityId,
        nx: f64,
        ny: f64,
        nz: f64,
    ) -> Result<EntityId, Solve3DError> {
        self.expect_point(origin)?;
        let nxv = self.vars.len();
        self.vars.push(nx);
        let nyv = self.vars.len();
        self.vars.push(ny);
        let nzv = self.vars.len();
        self.vars.push(nz);
        let id = EntityId(self.entities.len());
        self.entities.push(Entity3D::Plane3(Plane3 {
            origin,
            nx_var: nxv,
            ny_var: nyv,
            nz_var: nzv,
        }));
        Ok(id)
    }

    /// Add a workplane (same shape as a plane).
    pub fn add_workplane(
        &mut self,
        origin: EntityId,
        nx: f64,
        ny: f64,
        nz: f64,
    ) -> Result<EntityId, Solve3DError> {
        self.expect_point(origin)?;
        let nxv = self.vars.len();
        self.vars.push(nx);
        let nyv = self.vars.len();
        self.vars.push(ny);
        let nzv = self.vars.len();
        self.vars.push(nz);
        let id = EntityId(self.entities.len());
        self.entities.push(Entity3D::Workplane(Workplane {
            origin,
            nx_var: nxv,
            ny_var: nyv,
            nz_var: nzv,
        }));
        Ok(id)
    }

    /// Add a circle centred on `center` with `radius`, lying in the plane
    /// with normal `(nx, ny, nz)`. Radius and normal are free variables.
    pub fn add_circle(
        &mut self,
        center: EntityId,
        radius: f64,
        nx: f64,
        ny: f64,
        nz: f64,
    ) -> Result<EntityId, Solve3DError> {
        self.expect_point(center)?;
        let rv = self.vars.len();
        self.vars.push(radius);
        let nxv = self.vars.len();
        self.vars.push(nx);
        let nyv = self.vars.len();
        self.vars.push(ny);
        let nzv = self.vars.len();
        self.vars.push(nz);
        let id = EntityId(self.entities.len());
        self.entities.push(Entity3D::Circle3(Circle3 {
            center,
            radius_var: rv,
            nx_var: nxv,
            ny_var: nyv,
            nz_var: nzv,
        }));
        Ok(id)
    }

    /// (cx, cy, cz, radius, nx, ny, nz) for a circle entity.
    pub fn circle_data(&self, id: EntityId) -> (f64, f64, f64, f64, f64, f64, f64) {
        match &self.entities[id.0] {
            Entity3D::Circle3(c) => {
                let (cx, cy, cz) = self.point_xyz(c.center);
                (
                    cx,
                    cy,
                    cz,
                    self.vars[c.radius_var],
                    self.vars[c.nx_var],
                    self.vars[c.ny_var],
                    self.vars[c.nz_var],
                )
            }
            other => panic!("expected Circle3 at {id:?}, got {}", other.kind()),
        }
    }

    /// A circle's current radius.
    pub fn circle_radius(&self, id: EntityId) -> f64 {
        match &self.entities[id.0] {
            Entity3D::Circle3(c) => self.vars[c.radius_var],
            other => panic!("expected Circle3 at {id:?}, got {}", other.kind()),
        }
    }

    /// Add a circular arc centred on `center` with `radius`, in the plane
    /// with normal `(nx, ny, nz)`, bounded by the `start` and `end` points.
    #[allow(clippy::too_many_arguments)]
    pub fn add_arc(
        &mut self,
        center: EntityId,
        radius: f64,
        nx: f64,
        ny: f64,
        nz: f64,
        start: EntityId,
        end: EntityId,
    ) -> Result<EntityId, Solve3DError> {
        self.expect_point(center)?;
        self.expect_point(start)?;
        self.expect_point(end)?;
        let rv = self.vars.len();
        self.vars.push(radius);
        let nxv = self.vars.len();
        self.vars.push(nx);
        let nyv = self.vars.len();
        self.vars.push(ny);
        let nzv = self.vars.len();
        self.vars.push(nz);
        let id = EntityId(self.entities.len());
        self.entities.push(Entity3D::Arc3(Arc3 {
            center,
            radius_var: rv,
            nx_var: nxv,
            ny_var: nyv,
            nz_var: nzv,
            start,
            end,
        }));
        Ok(id)
    }

    /// (cx, cy, cz, radius, nx, ny, nz) for an arc's underlying circle.
    pub fn arc_circle(&self, id: EntityId) -> (f64, f64, f64, f64, f64, f64, f64) {
        match &self.entities[id.0] {
            Entity3D::Arc3(a) => {
                let (cx, cy, cz) = self.point_xyz(a.center);
                (
                    cx,
                    cy,
                    cz,
                    self.vars[a.radius_var],
                    self.vars[a.nx_var],
                    self.vars[a.ny_var],
                    self.vars[a.nz_var],
                )
            }
            other => panic!("expected Arc3 at {id:?}, got {}", other.kind()),
        }
    }

    /// The (start, end) endpoint ids of an arc.
    pub fn arc_endpoints(&self, id: EntityId) -> (EntityId, EntityId) {
        match &self.entities[id.0] {
            Entity3D::Arc3(a) => (a.start, a.end),
            other => panic!("expected Arc3 at {id:?}, got {}", other.kind()),
        }
    }

    /// An arc's current radius.
    pub fn arc_radius(&self, id: EntityId) -> f64 {
        match &self.entities[id.0] {
            Entity3D::Arc3(a) => self.vars[a.radius_var],
            other => panic!("expected Arc3 at {id:?}, got {}", other.kind()),
        }
    }

    /// Add a cubic Bézier spline through the four existing control points
    /// `p0` (start), `p1`, `p2`, `p3` (end).
    pub fn add_spline(
        &mut self,
        p0: EntityId,
        p1: EntityId,
        p2: EntityId,
        p3: EntityId,
    ) -> Result<EntityId, Solve3DError> {
        self.expect_point(p0)?;
        self.expect_point(p1)?;
        self.expect_point(p2)?;
        self.expect_point(p3)?;
        let id = EntityId(self.entities.len());
        self.entities
            .push(Entity3D::Spline3(Spline3 { p0, p1, p2, p3 }));
        Ok(id)
    }

    /// The four control-point ids of a spline `(p0, p1, p2, p3)`.
    pub fn spline_control_points(&self, id: EntityId) -> (EntityId, EntityId, EntityId, EntityId) {
        match &self.entities[id.0] {
            Entity3D::Spline3(sp) => (sp.p0, sp.p1, sp.p2, sp.p3),
            other => panic!("expected Spline3 at {id:?}, got {}", other.kind()),
        }
    }

    /// Evaluate a spline's cubic Bézier at parameter `t` ∈ [0, 1].
    pub fn spline_point_at(&self, id: EntityId, t: f64) -> (f64, f64, f64) {
        let (p0, p1, p2, p3) = self.spline_control_points(id);
        let (x0, y0, z0) = self.point_xyz(p0);
        let (x1, y1, z1) = self.point_xyz(p1);
        let (x2, y2, z2) = self.point_xyz(p2);
        let (x3, y3, z3) = self.point_xyz(p3);
        let u = 1.0 - t;
        let (b0, b1, b2, b3) = (u * u * u, 3.0 * u * u * t, 3.0 * u * t * t, t * t * t);
        (
            b0 * x0 + b1 * x1 + b2 * x2 + b3 * x3,
            b0 * y0 + b1 * y1 + b2 * y2 + b3 * y3,
            b0 * z0 + b1 * z1 + b2 * z2 + b3 * z3,
        )
    }

    /// Append a constraint.
    pub fn add_constraint(&mut self, c: Constraint3D) {
        self.constraints.push(c);
    }

    /// Sum of `n_residuals()` over all constraints.
    pub fn total_residuals(&self) -> usize {
        self.constraints.iter().map(|c| c.n_residuals()).sum()
    }

    /// (x, y, z) for a point entity. Panics if the id is wrong kind —
    /// callers should pre-validate with [`Self::expect_point`].
    pub fn point_xyz(&self, id: EntityId) -> (f64, f64, f64) {
        match &self.entities[id.0] {
            Entity3D::Point3(p) => p.read(&self.vars),
            other => panic!("expected Point3 at {id:?}, got {}", other.kind()),
        }
    }

    /// Endpoints of a line entity.
    pub fn line_endpoints(&self, id: EntityId) -> (EntityId, EntityId) {
        match &self.entities[id.0] {
            Entity3D::Line3(l) => (l.a, l.b),
            other => panic!("expected Line3 at {id:?}, got {}", other.kind()),
        }
    }

    /// (ox, oy, oz, nx, ny, nz) for a plane entity.
    pub fn plane_data(&self, id: EntityId) -> (f64, f64, f64, f64, f64, f64) {
        match &self.entities[id.0] {
            Entity3D::Plane3(p) => {
                let (ox, oy, oz) = self.point_xyz(p.origin);
                (
                    ox,
                    oy,
                    oz,
                    self.vars[p.nx_var],
                    self.vars[p.ny_var],
                    self.vars[p.nz_var],
                )
            }
            other => panic!("expected Plane3 at {id:?}, got {}", other.kind()),
        }
    }

    /// (ox, oy, oz, nx, ny, nz) for a workplane.
    pub fn workplane_data(&self, id: EntityId) -> (f64, f64, f64, f64, f64, f64) {
        match &self.entities[id.0] {
            Entity3D::Workplane(w) => {
                let (ox, oy, oz) = self.point_xyz(w.origin);
                (
                    ox,
                    oy,
                    oz,
                    self.vars[w.nx_var],
                    self.vars[w.ny_var],
                    self.vars[w.nz_var],
                )
            }
            other => panic!("expected Workplane at {id:?}, got {}", other.kind()),
        }
    }

    /// (ox, oy, oz, nx, ny, nz) for a plane *or* workplane — the two
    /// share an identical origin + normal layout. Panics on any other
    /// entity kind.
    pub fn plane_or_workplane_data(&self, id: EntityId) -> (f64, f64, f64, f64, f64, f64) {
        match &self.entities[id.0] {
            Entity3D::Plane3(_) => self.plane_data(id),
            Entity3D::Workplane(_) => self.workplane_data(id),
            other => panic!(
                "expected Plane3 or Workplane at {id:?}, got {}",
                other.kind()
            ),
        }
    }

    /// Pin a plane (or workplane) as a fixed datum: captures its
    /// current origin + normal and adds a [`Constraint3D::PlaneFixed`]
    /// holding them. A plane added with a free normal is otherwise a
    /// free body — pin it before constraining points against it as a
    /// reference plane.
    pub fn lock_plane(&mut self, plane: EntityId) -> Result<(), Solve3DError> {
        let (ox, oy, oz, nx, ny, nz) = match self.entities.get(plane.0) {
            Some(Entity3D::Plane3(_)) | Some(Entity3D::Workplane(_)) => {
                self.plane_or_workplane_data(plane)
            }
            Some(_) => {
                return Err(Solve3DError::KindMismatch {
                    id: plane,
                    expected: "Plane3 or Workplane",
                })
            }
            None => return Err(Solve3DError::UnknownEntity(plane)),
        };
        self.constraints.push(Constraint3D::PlaneFixed {
            plane,
            origin: [ox, oy, oz],
            normal: [nx, ny, nz],
        });
        Ok(())
    }

    /// Solve via the [`crate::solver`].
    pub fn solve(&mut self) -> Result<crate::solver::SolverReport, Solve3DError> {
        crate::solver::solve(self, crate::solver::SolverConfig::default())
    }

    /// Quick validity check: the id must resolve to a `Point3`.
    pub fn expect_point(&self, id: EntityId) -> Result<(), Solve3DError> {
        match self.entities.get(id.0) {
            Some(Entity3D::Point3(_)) => Ok(()),
            Some(_) => Err(Solve3DError::KindMismatch {
                id,
                expected: "Point3",
            }),
            None => Err(Solve3DError::UnknownEntity(id)),
        }
    }

    /// Like [`Self::expect_point`] but for an arbitrary entity kind:
    /// `id` must resolve in range and satisfy `ok`, else `UnknownEntity`
    /// / `KindMismatch`.
    fn expect_kind(
        &self,
        id: EntityId,
        expected: &'static str,
        ok: impl Fn(&Entity3D) -> bool,
    ) -> Result<(), Solve3DError> {
        match self.entities.get(id.0) {
            Some(e) if ok(e) => Ok(()),
            Some(_) => Err(Solve3DError::KindMismatch { id, expected }),
            None => Err(Solve3DError::UnknownEntity(id)),
        }
    }

    /// Check that every entity- and variable-index reference resolves in
    /// range and to the kind the solver's accessors will demand.
    ///
    /// The `add_*` builders enforce this at construction (via
    /// [`Self::expect_point`]), but serde deserialisation builds a
    /// `Sketch3D` field-by-field, bypassing them. A sketch loaded from a
    /// crafted/corrupt `.ron` could otherwise drive the residual
    /// accessors (`point_xyz`/`line_endpoints`/`circle_data`/…), which
    /// raw-index `entities[id.0]` / `vars[var]` and `panic!` on a wrong
    /// kind, into a crash on the first re-solve.
    /// [`crate::persist::from_ron_str`] calls this on load.
    pub fn validate(&self) -> Result<(), Solve3DError> {
        let n_vars = self.vars.len();
        let var = |v: usize| -> Result<(), Solve3DError> {
            if v < n_vars {
                Ok(())
            } else {
                Err(Solve3DError::BadParameter {
                    name: "var",
                    reason: format!("variable index {v} out of range (have {n_vars})"),
                })
            }
        };
        // Entities: every variable index in range; every entity->entity
        // reference resolves to a `Point3` (all of them are point refs).
        // Validating the table first lets the constraint pass below trust
        // that a referenced `Line3`/`Circle3`/… has sound sub-references.
        for e in &self.entities {
            match e {
                Entity3D::Point3(p) => {
                    var(p.x_var)?;
                    var(p.y_var)?;
                    var(p.z_var)?;
                }
                Entity3D::Line3(l) => {
                    self.expect_point(l.a)?;
                    self.expect_point(l.b)?;
                }
                Entity3D::Plane3(p) => {
                    self.expect_point(p.origin)?;
                    var(p.nx_var)?;
                    var(p.ny_var)?;
                    var(p.nz_var)?;
                }
                Entity3D::Workplane(w) => {
                    self.expect_point(w.origin)?;
                    var(w.nx_var)?;
                    var(w.ny_var)?;
                    var(w.nz_var)?;
                }
                Entity3D::Circle3(c) => {
                    self.expect_point(c.center)?;
                    var(c.radius_var)?;
                    var(c.nx_var)?;
                    var(c.ny_var)?;
                    var(c.nz_var)?;
                }
                Entity3D::Arc3(a) => {
                    self.expect_point(a.center)?;
                    self.expect_point(a.start)?;
                    self.expect_point(a.end)?;
                    var(a.radius_var)?;
                    var(a.nx_var)?;
                    var(a.ny_var)?;
                    var(a.nz_var)?;
                }
                Entity3D::Spline3(sp) => {
                    self.expect_point(sp.p0)?;
                    self.expect_point(sp.p1)?;
                    self.expect_point(sp.p2)?;
                    self.expect_point(sp.p3)?;
                }
            }
        }
        // Constraints: every directly-referenced entity resolves in range
        // and to the kind its residual accessor demands (a referenced
        // `Line3`/`Circle3`/`Arc3`'s own point sub-refs were checked above).
        for c in &self.constraints {
            match c {
                Constraint3D::Coincident3 { a, b }
                | Constraint3D::SamePoint { a, b }
                | Constraint3D::PointDistance3 { a, b, .. } => {
                    self.expect_point(*a)?;
                    self.expect_point(*b)?;
                }
                Constraint3D::LineLength3 { line, .. } => {
                    self.expect_kind(*line, "Line3", |e| matches!(e, Entity3D::Line3(_)))?;
                }
                Constraint3D::LineAngle3 { a, b, .. }
                | Constraint3D::LineParallel3 { a, b }
                | Constraint3D::LinePerpendicular3 { a, b } => {
                    self.expect_kind(*a, "Line3", |e| matches!(e, Entity3D::Line3(_)))?;
                    self.expect_kind(*b, "Line3", |e| matches!(e, Entity3D::Line3(_)))?;
                }
                Constraint3D::PointInPlane { point, plane } => {
                    self.expect_point(*point)?;
                    self.expect_kind(*plane, "Plane3", |e| matches!(e, Entity3D::Plane3(_)))?;
                }
                Constraint3D::OnPlane { point, workplane } => {
                    self.expect_point(*point)?;
                    self.expect_kind(*workplane, "Workplane", |e| {
                        matches!(e, Entity3D::Workplane(_))
                    })?;
                }
                Constraint3D::OnLine { point, line } => {
                    self.expect_point(*point)?;
                    self.expect_kind(*line, "Line3", |e| matches!(e, Entity3D::Line3(_)))?;
                }
                Constraint3D::PlaneFixed { plane, .. } => {
                    self.expect_kind(*plane, "Plane3 or Workplane", |e| {
                        matches!(e, Entity3D::Plane3(_) | Entity3D::Workplane(_))
                    })?;
                }
                Constraint3D::CircleRadius { circle, .. } => {
                    self.expect_kind(*circle, "Circle3", |e| matches!(e, Entity3D::Circle3(_)))?;
                }
                Constraint3D::PointOnCircle { point, circle } => {
                    self.expect_point(*point)?;
                    self.expect_kind(*circle, "Circle3", |e| matches!(e, Entity3D::Circle3(_)))?;
                }
                Constraint3D::EqualRadius { a, b } => {
                    self.expect_kind(*a, "Circle3", |e| matches!(e, Entity3D::Circle3(_)))?;
                    self.expect_kind(*b, "Circle3", |e| matches!(e, Entity3D::Circle3(_)))?;
                }
                Constraint3D::ArcRadius { arc, .. } | Constraint3D::ArcEndpointsOnArc { arc } => {
                    self.expect_kind(*arc, "Arc3", |e| matches!(e, Entity3D::Arc3(_)))?;
                }
            }
        }
        Ok(())
    }
}
