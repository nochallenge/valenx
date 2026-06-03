//! # valenx-hvac
//!
//! HVAC ducting + equipment placement.  Defines rectangular and round
//! ducts, their per-section CAD solids, an equipment kind enumeration
//! with parametric solids, plus sizing and pressure-drop helpers.
//!
//! Phase 50 of the FreeCAD-parity roadmap. FreeCAD `HVAC` community
//! workbench equivalent.
//!
//! # Surface
//!
//! - [`duct::Duct`] + [`duct::CrossSection`] (Rect / Round) +
//!   [`duct::to_solid`].
//! - [`equipment::Equipment`] (AHU / VAV / Diffuser / Grille / Fan /
//!   Heater / Chiller / Damper) + [`equipment::to_solid`].
//! - [`flow::cfm_to_duct_size`] — CFM + max-velocity → duct in/in.
//! - [`pressure_drop::darcy_weisbach`] — ΔP across a duct length.
//! - [`panel::HvacPanelState`] — UI state envelope.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod duct;
pub mod equipment;
pub mod error;
pub mod flow;
pub mod panel;
pub mod pressure_drop;

pub use duct::{to_solid as duct_to_solid, CrossSection, Duct};
pub use equipment::{to_solid as equipment_to_solid, Equipment};
pub use error::{ErrorCategory, HvacError};
pub use flow::cfm_to_duct_size;
pub use panel::HvacPanelState;
pub use pressure_drop::{darcy_weisbach, darcy_weisbach_rho};
