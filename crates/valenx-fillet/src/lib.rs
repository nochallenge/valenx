//! # valenx-fillet
//!
//! Mesh-domain fillet and chamfer for Valenx.
//!
//! `valenx-fillet` takes a triangle mesh, classifies its edges, and
//! replaces every sharp convex edge with either:
//!
//! - a smooth cylindrical strip (`apply_fillet`), or
//! - a flat bevel strip (`apply_chamfer`).
//!
//! The result is a new [`valenx_mesh::Mesh`] with the corners rounded
//! or beveled.
//!
//! # v1 honest scope
//!
//! v1 is a **mesh-domain** operation, not a true BRep fillet:
//!
//! - The output is a triangle mesh; topology of the original BRep is
//!   discarded.
//! - Strips are added on top of the original mesh — the original
//!   sharp triangles are **not** clipped to make room for the strip
//!   (a known overlap of size `radius`).
//! - Corners where 3+ filleted edges meet are not blended with a
//!   spherical patch; strips simply overlap each other.
//!
//! For visualization and 3D-printing, the overlap is invisible from
//! the outside and the print is correct. For STEP/IGES export with
//! parametric fillet preservation, a true BRep fillet
//! (Phase 3.5) is required.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod chamfer;
pub mod cyl_strip;
pub mod edge_graph;
pub mod error;
pub mod fillet;

pub use chamfer::apply_chamfer;
pub use error::{ErrorCategory, FilletError};
pub use fillet::apply_fillet;
