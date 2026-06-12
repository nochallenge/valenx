//! Model analysis — features 18 through 21.
//!
//! Higher-level tasks built on top of the simulation layer: they ask
//! questions *about* a model rather than just running it.
//!
//! - [`param`] — [`ParamTarget`], the shared "address one knob of a
//!   model" abstraction every analysis depends on.
//! - [`scan`] — 1-D and 2-D parameter scans (feature 18).
//! - [`sensitivity`] — local finite-difference sensitivity and the
//!   Morris global elementary-effects screen (feature 19).
//! - [`conservation`] — conserved-moiety detection from the
//!   stoichiometry-matrix null space (feature 20).
//! - [`bifurcation`] — a steady-state continuation bifurcation scan
//!   with stability classification (feature 21).

pub mod bifurcation;
pub mod conservation;
pub mod estimation;
pub mod param;
pub mod scan;
pub mod sensitivity;

pub use bifurcation::{bifurcation_scan, BifurcationDiagram, BifurcationPoint, Stability};
pub use conservation::{conservation_laws, ConservationLaw};
pub use estimation::{
    estimate_parameters, EstimationOptions, EstimationReport, EstimationTarget, ObservedPoint,
};
pub use param::ParamTarget;
pub use scan::{linspace, scan_1d, scan_2d, Scan1d, Scan2d};
pub use sensitivity::{
    global_sensitivity, local_sensitivity, relative_ranges, GlobalSensitivity, LocalSensitivity,
};
