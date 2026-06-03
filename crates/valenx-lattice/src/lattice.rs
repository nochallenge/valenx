//! Lattice recipes.
//!
//! See the crate docs for the full enumeration of variants.

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

/// Which features of a mesh become placement positions.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub enum MeshSamplingMode {
    /// Place one instance at each mesh node.
    Vertices,
    /// Place one instance at the centroid of each triangle.
    FaceCentroids,
}

/// One lattice recipe. Each variant carries its own parameters; the
/// dispatch lives in [`crate::generate::generate`].
#[derive(Clone, Debug, Serialize, Deserialize)]
pub enum Lattice {
    /// Box grid: `rows` × `cols` × `levels` with per-axis spacing.
    Grid {
        /// Number of rows (x).
        rows: usize,
        /// Number of columns (y).
        cols: usize,
        /// Number of levels (z).
        levels: usize,
        /// Spacing in each axis.
        spacing: Vector3<f64>,
    },
    /// Polar (circular) array around an axis.
    Polar {
        /// Centre of rotation.
        center: Vector3<f64>,
        /// Axis to rotate around.
        axis: Vector3<f64>,
        /// Number of placements.
        count: usize,
        /// Total sweep angle (radians).
        total_angle: f64,
    },
    /// `n_samples` placements along a 3D Bezier curve.
    Bezier {
        /// Control points (degree = control_points.len() - 1).
        control_points: Vec<Vector3<f64>>,
        /// Number of placements.
        n_samples: usize,
    },
    /// Placements along a NURBS curve from `valenx-surface`.
    OnCurve {
        /// The source curve (serialised inline).
        curve: valenx_surface::NurbsCurve,
        /// Sample count.
        n_samples: usize,
    },
    /// Placements at iso-parameter crosses on a NURBS surface.
    OnSurface {
        /// The source surface (serialised inline).
        surface: valenx_surface::NurbsSurface,
        /// Number of u samples.
        n_u: usize,
        /// Number of v samples.
        n_v: usize,
    },
    /// Placements at mesh features (vertices or face centroids).
    OnMesh {
        /// The mesh (serialised inline).
        mesh: valenx_mesh::Mesh,
        /// Whether to sample vertices or face centroids.
        mode: MeshSamplingMode,
    },
}

impl Lattice {
    /// Short label used in error messages and the UI dropdown.
    pub fn label(&self) -> &'static str {
        match self {
            Self::Grid { .. } => "Grid",
            Self::Polar { .. } => "Polar",
            Self::Bezier { .. } => "Bezier",
            Self::OnCurve { .. } => "OnCurve",
            Self::OnSurface { .. } => "OnSurface",
            Self::OnMesh { .. } => "OnMesh",
        }
    }
}
