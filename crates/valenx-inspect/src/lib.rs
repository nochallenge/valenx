//! # valenx-inspect
//!
//! Inspection workbench — geometric measurement and GD&T tolerance
//! verification. This is the FreeCAD `Inspection` community workbench
//! re-imagined as a pure-data Rust crate.
//!
//! Phase 25 of the FreeCAD-parity roadmap.
//!
//! # Two layers
//!
//! - [`measurement`] — compute scalar geometric quantities from
//!   primitives: distance between two points, angle at a vertex, radius
//!   of a circle, polyline length, polygon area, mesh volume, mesh
//!   bounding-box diagonal.
//! - [`tolerance`] / [`gdt`] — verify a measured value against a
//!   nominal-with-deviations [`Tolerance`] band, or against an ASME
//!   Y14.5 [`valenx_techdraw::gdt::GdtSymbol`] feature-control frame.
//!
//! Aggregate via [`InspectReport`] (a list of
//! (measurement, tolerance, result) rows) which serializes to RON and
//! exports to CSV for hand-off to a metrology lab.
//!
//! # Example
//!
//! ```
//! use valenx_inspect::{Measurement, Tolerance, CheckResult, InspectReport};
//! use nalgebra::Vector3;
//!
//! let m = Measurement::Distance {
//!     from: Vector3::new(0.0, 0.0, 0.0),
//!     to:   Vector3::new(10.0, 0.0, 0.0),
//! };
//! let actual = valenx_inspect::measurement::compute(&m).unwrap();
//! let tol = Tolerance::symmetric(10.0, 0.05);
//! assert_eq!(actual, 10.0);
//! assert_eq!(tol.evaluate(actual), CheckResult::Pass);
//!
//! let mut report = InspectReport::new();
//! report.add_row(m, tol, CheckResult::Pass);
//! assert!(report.to_csv().contains("Pass"));
//! ```

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod error;
pub mod gdt;
pub mod measurement;
pub mod persist;
pub mod report;
pub mod tolerance;

pub use error::{ErrorCategory, InspectError};
pub use gdt::{GdtCheck, GdtRule};
pub use measurement::Measurement;
pub use persist::InspectFile;
pub use report::{CheckResult, InspectReport, ReportRow};
pub use tolerance::Tolerance;
