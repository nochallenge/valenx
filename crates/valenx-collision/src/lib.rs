//! # valenx-collision
//!
//! Advanced assembly collision detection — AABB prune + triangle-
//! triangle Möller separating-axis test. Pairwise check across
//! every distinct part in a [`valenx_assembly::Assembly`].
//!
//! Phase 44 of the FreeCAD-parity roadmap. Sister to Phase 6
//! Assembly mate solver.
//!
//! # Surface
//!
//! - [`Aabb`] / [`aabb::intersect`] / [`aabb::distance`] — fast
//!   bounding-box overlap test.
//! - [`mesh_pair::collide`] — AABB-pruned per-triangle SAT test;
//!   returns the first [`mesh_pair::CollisionInfo`] hit.
//! - [`assembly::check_collisions`] — pairwise check across an
//!   assembly with single tessellation per part.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod aabb;
pub mod assembly;
pub mod error;
pub mod mesh_pair;

pub use aabb::{distance, distance_squared, intersect, Aabb};
pub use assembly::check_collisions;
pub use error::{CollisionError, ErrorCategory};
pub use mesh_pair::{collide, CollisionInfo};
