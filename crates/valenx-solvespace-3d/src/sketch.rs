//! [`Sketch3D`] — the 3D analogue of `valenx_sketch::Sketch`.
//!
//! Owns the global `vars: Vec<f64>` variable vector plus the heterogeneous
//! `entities: Vec<Entity3D>` table. Constraints reference entities by
//! [`EntityId`], and entities reference variables by index — exactly the
//! same structure-of-arrays pattern Phase 1 uses.

use serde::{Deserialize, Serialize};

use crate::constraint::Constraint3D;
use crate::entity::{Circle3, Entity3D, EntityId, Line3, Plane3, Point3, Workplane};
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
            other => panic!("expected Plane3 or Workplane at {id:?}, got {}", other.kind()),
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
}
