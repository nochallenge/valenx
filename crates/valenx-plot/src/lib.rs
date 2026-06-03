//! # valenx-plot
//!
//! Plot workbench — 2D plotting + charting. The FreeCAD `Plot`
//! community workbench equivalent.
//!
//! Phase 31 of the FreeCAD-parity roadmap.
//!
//! # Pipeline
//!
//! 1. Build one or more [`Series`] — `name + Vec<(f64, f64)>` data +
//!    a [`SeriesStyle`] (Line / Scatter / Bar / Area).
//! 2. Wrap them in a [`Plot`] — title + x_label + y_label + optional
//!    explicit axis ranges (defaults to auto-fit).
//! 3. Call [`to_svg`] to render — returns a `String` containing a
//!    standalone `<svg>` document with axes / ticks / labels / series.
//! 4. [`to_png`] returns [`PlotError::PngDeferred`] in v1 — would need
//!    the `image` crate; Phase 31.5 follow-up.
//!
//! Persistence via [`PlotFile`] RON envelope.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod error;
pub mod persist;
pub mod plot;
pub mod series;
pub mod svg;

pub use error::{ErrorCategory, PlotError};
pub use persist::PlotFile;
pub use plot::{AxisRange, Plot};
pub use series::{Series, SeriesStyle};
pub use svg::{to_png, to_svg};
