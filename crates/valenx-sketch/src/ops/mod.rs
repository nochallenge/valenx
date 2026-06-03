//! Sketch-level operations — mirror, copy, move, rotate, arrays.
//!
//! Phase 12E.
//!
//! These operations *create new entities* (mirror/copy/array) or
//! *mutate existing variable values* (move/rotate). All ops operate
//! on a list of [`crate::geom::EntityId`]s and return the ids of any
//! newly created entities.

pub mod copy;
pub mod linear_array;
pub mod mirror;
#[allow(clippy::module_inception)]
pub mod r#move;
pub mod polar_array;
#[cfg(test)]
mod regression;
pub mod rotate;
