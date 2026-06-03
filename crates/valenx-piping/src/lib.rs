//! # valenx-piping
//!
//! P&ID-aware piping system — pipe sections, fittings, valves, and an
//! NPS-to-OD lookup. Each [`PipeSection`] knows its NPS designation,
//! schedule, and material, and can emit a CAD [`valenx_cad::Solid`]
//! for visualisation.
//!
//! Phase 49 of the FreeCAD-parity roadmap. FreeCAD `Pipes & Tubing`
//! community workbench equivalent.
//!
//! # Surface
//!
//! - [`pipe::PipeSection`] + [`pipe::Material`] + [`pipe::to_solid`].
//! - [`fitting::PipeFitting`] (Elbow90/45, Tee, Reducer, Cap,
//!   Coupling, Union) + [`fitting::Valve`] (Gate, Globe, Ball,
//!   Check, Butterfly, Needle, Diaphragm).
//! - [`dims::nominal_to_od_in`] / [`dims::nominal_to_od_mm`] +
//!   [`dims::wall_thickness_in`] + [`dims::Schedule`].
//! - [`system::Piping`] / [`system::Junction`] /
//!   [`system::JunctionKind`].
//! - [`panel::PipingPanelState`] — UI state envelope.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod dims;
pub mod error;
pub mod fitting;
pub mod panel;
pub mod pipe;
pub mod system;

pub use dims::{nominal_to_od_in, nominal_to_od_mm, wall_thickness_in, Schedule};
pub use error::{ErrorCategory, PipingError};
pub use fitting::{PipeFitting, Valve};
pub use panel::PipingPanelState;
pub use pipe::{to_solid, Material, PipeSection};
pub use system::{Junction, JunctionKind, Piping};
