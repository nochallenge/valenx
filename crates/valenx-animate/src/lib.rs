//! # valenx-animate
//!
//! Animation workbench — keyframe animation of assembly joint
//! parameters (revolute angle / prismatic distance / etc.). Builds on
//! the Phase 6 kinematics pipeline.
//!
//! Phase 29 of the FreeCAD-parity roadmap.
//!
//! # Data model
//!
//! - [`Keyframe`] — `time: f64` (seconds) + a sparse
//!   `joint_parameters: HashMap<JointId, f64>` mapping joint id →
//!   target parameter.
//! - [`Animation`] — `Vec<Keyframe>` sorted by `time`. Sampling at any
//!   `t` interpolates the parameter for each joint between its
//!   bracketing keyframes using the keyframe's [`TweenMode`].
//! - [`Animation::frames`] — sample the animation across
//!   `[start_t, end_t]` at `fps` and return the per-frame
//!   joint-parameter snapshots. [`Animation::apply_frame`] then
//!   applies one snapshot to a clone of the assembly via
//!   `valenx_assembly::kinematics::apply_all_joints`.
//!
//! # Tween modes
//!
//! Linear / EaseIn / EaseOut / EaseInOut / Cubic. Selected on the
//! *outgoing* keyframe (i.e. when interpolating from keyframe `i` to
//! keyframe `i+1`, the [`TweenMode`] from `i+1.tween` is used).

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod animation;
pub mod error;
pub mod persist;
pub mod tween;

pub use animation::{Animation, FrameSnapshot, Keyframe};
pub use error::{AnimateError, ErrorCategory};
pub use persist::AnimateFile;
pub use tween::TweenMode;
