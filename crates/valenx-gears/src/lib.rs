//! # valenx-gears
//!
//! Involute spur / helical / bevel / worm gear generator.
//!
//! Phase 39 of the FreeCAD-parity roadmap. FreeCAD `Gears` /
//! `GearWB` community workbench analogue.
//!
//! # Surface
//!
//! - [`GearKind`] / [`GearSpec`] — gear family + parametric inputs
//!   (module, teeth, pressure angle, helix angle, face width).
//! - [`profile::involute_point`] — canonical involute curve point.
//! - [`profile::tooth_profile`] / [`profile::full_profile`] — 2D
//!   tooth + complete outline.
//! - [`solid::to_solid_spur`] / [`solid::to_solid_helical`] /
//!   [`solid::to_solid_bevel`] / [`solid::to_solid_worm`] —
//!   tessellated 3D output.
//! - [`solid::to_solid`] — dispatcher.
//! - [`spec::base_pitch_mm`](GearSpec::base_pitch_mm) /
//!   [`spec::contact_ratio`] — meshing geometry: the base pitch and the
//!   transverse contact ratio `mₚ` of a gear pair (`> 1` for continuous
//!   transmission).
//! - [`spec::lewis_bending_stress_mpa`] — Lewis tooth-root bending
//!   stress `σ = W_t / (F·m·Y)` (first-order strength screen).

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod error;
pub mod profile;
pub mod solid;
pub mod spec;

pub use error::{ErrorCategory, GearsError};
pub use profile::{full_profile, involute_point, tooth_profile};
pub use solid::{to_solid, to_solid_bevel, to_solid_helical, to_solid_spur, to_solid_worm};
pub use spec::{
    circular_pitch_mm, contact_ratio, gear_ratio, lewis_bending_stress_mpa, GearKind, GearSpec,
};
