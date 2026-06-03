//! # valenx-manipulator
//!
//! Free-form direct push/pull editing of BRep solids.
//!
//! Phase 41 of the FreeCAD-parity roadmap. FreeCAD `Manipulator` /
//! `Direct Modelling` community workbench analogue.
//!
//! v1 pipeline: tessellate the input BRep (default chord error
//! 0.5 mm — caller's choice if it already has a finer-resolution
//! mesh), mutate the resulting triangle mesh in place, and return
//! a [`valenx_cad::Solid::Mesh`].
//!
//! # Surface
//!
//! - [`ManipulateOp`] — six push/pull primitives: MoveFace,
//!   RotateFace, MoveEdge, MoveVertex, ExtrudeFace, OffsetFace.
//! - [`apply::apply`] — single op.
//! - [`apply::apply_sequence`] — chain of ops.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod apply;
pub mod error;
pub mod op;

pub use apply::{apply, apply_sequence, DEFAULT_TOLERANCE_MM};
pub use error::{ErrorCategory, ManipulatorError};
pub use op::ManipulateOp;
