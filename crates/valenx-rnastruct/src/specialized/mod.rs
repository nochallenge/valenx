//! Specialized tools and structure-analysis utilities.
//!
//! The remaining single-purpose tools:
//!
//! - [`trna`] — tRNA cloverleaf structure scan / detection.
//! - [`stats`] — structure statistics: hairpin / bulge / internal /
//!   multiloop / stem counts and loop sizes.
//! - [`mountain`] — mountain-plot data (single structure and
//!   ensemble-averaged).
//! - [`dotplot`] — energy / probability dot plots (mfold-class).
//! - [`report`] — batch folding utilities and the [`report::FoldingReport`]
//!   bundle (MFE, ensemble free energy, centroid, MEA, MFE frequency).
//!
//! The 2-D drawing layout lives one level up in [`mod@crate::layout`].

pub mod dotplot;
pub mod mountain;
pub mod report;
pub mod stats;
pub mod trna;
