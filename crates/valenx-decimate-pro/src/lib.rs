//! # valenx-decimate-pro
//!
//! Curvature-aware, UV-preserving, and feature-aware mesh decimation
//! layered on top of valenx-mesh's Phase 7 QEM (Garland & Heckbert).
//!
//! Phase 47 of the FreeCAD-parity roadmap.  FreeCAD `MeshDecimation
//! Pro` community workbench equivalent.
//!
//! # Surface
//!
//! - [`curvature::per_vertex`] — discrete mean curvature via the
//!   cotangent Laplacian.
//! - [`appearance::QuadricMatrix`] + [`appearance::uv_aware_quadric`] —
//!   per-vertex 2x2 UV quadric for UV-preserving cost.
//! - [`decimate_pro::weighted_qem`] — curvature-aware QEM.
//! - [`decimate_pro::uv_preserving`] — UV-preserving QEM (returns the
//!   re-mapped UV array alongside the mesh).
//! - [`decimate_pro::feature_aware`] — protects every endpoint of a
//!   feature-edge list from collapse.
//! - [`panel::DecimateProPanelState`] — UI state envelope.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod appearance;
pub mod curvature;
pub mod decimate_pro;
pub mod error;
pub mod panel;

pub use appearance::{uv_aware_quadric, uv_stretch_weight, QuadricMatrix};
pub use curvature::per_vertex;
pub use decimate_pro::{feature_aware, uv_preserving, weighted_qem};
pub use error::{DecimateProError, ErrorCategory};
pub use panel::{DecimateProPanelState, Mode};
