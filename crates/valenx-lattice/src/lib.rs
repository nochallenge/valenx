//! # valenx-lattice
//!
//! Lattice2 workbench — advanced parametric arrays beyond Part
//! Design's basic Linear / Circular patterns. The FreeCAD Lattice2
//! community workbench equivalent.
//!
//! Phase 28 of the FreeCAD-parity roadmap.
//!
//! # Concept
//!
//! A [`Placement`] = `(position, orientation)`. A
//! [`Lattice`] is a recipe for generating a list of placements
//! (`generate(&Lattice) -> Vec<Placement>`). The desktop shell
//! instances a source `valenx_cad::Solid` at each placement to
//! produce a `Vec<Solid>` — that step is left to the caller because
//! solid-cloning lives downstream.
//!
//! # Generators
//!
//! - [`Lattice::Grid`] — n_x × n_y × n_z box grid with per-axis spacing.
//! - [`Lattice::Polar`] — `count` placements around `axis` over
//!   `total_angle` radians.
//! - [`Lattice::Bezier`] — `n_samples` placements along a 3D Bezier
//!   curve defined by `control_points`.
//! - [`Lattice::OnCurve`] — `n_samples` placements along a
//!   [`valenx_surface::NurbsCurve`].
//! - [`Lattice::OnSurface`] — `n_u × n_v` placements at iso-parameter
//!   crosses on a [`valenx_surface::NurbsSurface`].
//! - [`Lattice::OnMesh`] — placements at each mesh vertex (mode
//!   `Vertices`) or triangle centroid (mode `FaceCentroids`).

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod error;
pub mod generate;
pub mod lattice;
pub mod persist;
pub mod placement;

pub use error::{ErrorCategory, LatticeError};
pub use generate::generate;
pub use lattice::{Lattice, MeshSamplingMode};
pub use persist::LatticeFile;
pub use placement::Placement;
