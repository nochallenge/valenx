//! # valenx-curves
//!
//! Curves workbench — higher-level operations on the
//! [`valenx_surface::NurbsCurve`] type that Phase 9 / Phase 19
//! introduced. This is the FreeCAD `Curves` community workbench
//! equivalent.
//!
//! Phase 27 of the FreeCAD-parity roadmap.
//!
//! # Operations
//!
//! - **Offset** — parallel-offset a planar curve by `d` mm.
//! - **BlendCorner** — fillet two curves with a circular-arc blend.
//! - **Approximate** — fit a NURBS curve through a polyline cloud.
//! - **Project** — project a curve onto a NURBS surface.
//! - **Discretize** — sample a curve into an even / chordal /
//!   curvature-adapted polyline.
//! - **Reverse** — flip parameter direction.
//! - **Trim** — restrict a curve to `[t_start, t_end]`.
//! - **Extend** — extend a curve past its end by a tangent
//!   continuation.
//! - **IsoCurve** — extract an iso-u or iso-v curve from a surface.
//!
//! All operations are pure functions returning either a new
//! [`valenx_surface::NurbsCurve`] (or polyline `Vec<Vector3<f64>>`)
//! plus a [`CurvesError`] on bad inputs.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod error;
pub mod iso;
pub mod ops;
pub mod persist;

pub use error::{CurvesError, ErrorCategory};
pub use iso::{extract_iso, IsoKind};
pub use ops::{
    approximate, blend_corner, discretize, extend, offset_planar, project_curve, reverse, trim,
    DiscretizeMode, ExtendEnd,
};
pub use persist::CurvesFile;
