//! Per-operation evaluators called from [`crate::replay::replay`].
//!
//! Each submodule exposes a single `evaluate(...)` function that takes
//! the relevant [`crate::tree::FeatureTree`] context plus its operation
//! parameters and produces a [`valenx_cad::Solid`]. Tasks 9 / 14 / 19 /
//! 24 / 28 / 32 fill these in one by one; until then they return
//! [`crate::FeatureError::BadParameter`] with a `"not implemented
//! (stub)"` reason so the replay pipeline can be exercised end-to-end.

pub mod boolean_history;
pub mod chamfer;
pub mod draft_angle;
pub mod fillet;
pub mod helix;
pub mod hole;
pub mod imported_advanced;
pub mod imported_solid;
pub mod loft;
pub mod mirror;
pub mod multi_transform;
pub mod pad;
pub mod pattern_circular;
pub mod pattern_linear;
pub mod pipe;
pub mod pocket;
pub mod revolve;
pub mod shell;
pub mod sweep;
pub mod thickness;
