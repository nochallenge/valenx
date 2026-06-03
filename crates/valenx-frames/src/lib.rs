//! # valenx-frames
//!
//! Structural-frame workbench — sweep a [`Profile`] cross-section
//! along a polyline path to build mesh-backed members; collect
//! members into a [`Frame`] with auto-detected joints.
//!
//! Phase 38 of the FreeCAD-parity roadmap. FreeCAD `Frames` /
//! `Reinforcement Members` community workbench analogue.
//!
//! # Surface
//!
//! - [`Profile`] — 6 parametric section catalogues (I-beam,
//!   C-channel, L-angle, rect HSS, round CHS, T-beam).
//! - [`profile::cross_section_polygon`] — closed CCW outline in the
//!   local (u, v) frame.
//! - [`Member`] + [`member::to_solid`] — linear-segment sweep that
//!   returns a [`valenx_cad::Solid::Mesh`].
//! - [`Frame`] + [`Frame::auto_joints`] — collection of members
//!   with coincident-endpoint joint detection.
//! - [`frame::to_ron_string`] / [`frame::from_ron_str`] — round-trip
//!   through a versioned [`frame::FrameFile`] envelope.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod error;
pub mod frame;
pub mod member;
pub mod profile;

pub use error::{ErrorCategory, FramesError};
pub use frame::{Frame, FrameFile, Joint, VERSION};
pub use member::Member;
pub use profile::{cross_section_polygon, Profile};
