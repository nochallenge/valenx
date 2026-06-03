//! Geometric entities living in a [`crate::sketch::Sketch3D`].
//!
//! Like `valenx_sketch`'s 2D entities, every primitive only stores
//! variable indices — never raw coordinates — so the solver can mutate
//! the global variable vector in place.

use serde::{Deserialize, Serialize};

/// Stable handle into [`crate::sketch::Sketch3D::entities`].
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash, Serialize, Deserialize)]
pub struct EntityId(pub usize);

/// 3D point — three variable indices (x, y, z).
#[derive(Copy, Clone, Debug, Serialize, Deserialize)]
pub struct Point3 {
    /// Variable index for X.
    pub x_var: usize,
    /// Variable index for Y.
    pub y_var: usize,
    /// Variable index for Z.
    pub z_var: usize,
}

impl Point3 {
    /// Read the current (x, y, z) from the variable vector.
    pub fn read(&self, vars: &[f64]) -> (f64, f64, f64) {
        (vars[self.x_var], vars[self.y_var], vars[self.z_var])
    }
}

/// Line segment between two [`Point3`]s.
#[derive(Copy, Clone, Debug, Serialize, Deserialize)]
pub struct Line3 {
    /// First endpoint (id into entity table — must resolve to `Point3`).
    pub a: EntityId,
    /// Second endpoint.
    pub b: EntityId,
}

/// Implicit plane — origin point + normal direction. The normal is
/// stored as three free variables so constraints can rotate it; an
/// auxiliary unit-length constraint keeps it normalised.
#[derive(Copy, Clone, Debug, Serialize, Deserialize)]
pub struct Plane3 {
    /// Origin point id (must resolve to `Point3`).
    pub origin: EntityId,
    /// Normal X variable.
    pub nx_var: usize,
    /// Normal Y variable.
    pub ny_var: usize,
    /// Normal Z variable.
    pub nz_var: usize,
}

/// Workplane = origin + normal (no in-plane axes; we treat it like a
/// plane for the v1 constraint set).
#[derive(Copy, Clone, Debug, Serialize, Deserialize)]
pub struct Workplane {
    /// Origin point id.
    pub origin: EntityId,
    /// Normal X variable.
    pub nx_var: usize,
    /// Normal Y variable.
    pub ny_var: usize,
    /// Normal Z variable.
    pub nz_var: usize,
}

/// Circle — a centre point, a free radius variable, and a free normal
/// vector defining the circle's plane (three free variables, like a plane
/// normal).
#[derive(Copy, Clone, Debug, Serialize, Deserialize)]
pub struct Circle3 {
    /// Centre point id (must resolve to `Point3`).
    pub center: EntityId,
    /// Variable index for the radius.
    pub radius_var: usize,
    /// Normal X variable.
    pub nx_var: usize,
    /// Normal Y variable.
    pub ny_var: usize,
    /// Normal Z variable.
    pub nz_var: usize,
}

/// Sum-type wrapper so a sketch can hold heterogeneous entities in a
/// single `Vec`.
#[derive(Copy, Clone, Debug, Serialize, Deserialize)]
pub enum Entity3D {
    /// 3D point.
    Point3(Point3),
    /// 3D line segment.
    Line3(Line3),
    /// 3D plane.
    Plane3(Plane3),
    /// 3D workplane.
    Workplane(Workplane),
    /// 3D circle.
    Circle3(Circle3),
}

impl Entity3D {
    /// Discriminator label, used by [`crate::error::Solve3DError::KindMismatch`].
    pub fn kind(&self) -> &'static str {
        match self {
            Entity3D::Point3(_) => "Point3",
            Entity3D::Line3(_) => "Line3",
            Entity3D::Plane3(_) => "Plane3",
            Entity3D::Workplane(_) => "Workplane",
            Entity3D::Circle3(_) => "Circle3",
        }
    }
}
