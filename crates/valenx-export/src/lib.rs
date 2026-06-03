//! # valenx-export
//!
//! Export Valenx run results to formats other tools consume. First
//! concrete chunk of [RFC 0012](../../../rfcs/0012-ml-training-data-export.md).
//!
//! Today: **CSV writer** for any `Results.scalars` catalog (one row
//! per ScalarRecord, columns: name / value / units / time-key /
//! source). Useful as both a debug tool and a quick handoff to
//! pandas / Excel.
//!
//! Also today: **npy writer** ([`write_npy_f64`]) — emits the
//! standard NumPy `.npy` format that pandas / NumPy / PyTorch /
//! TensorFlow can load directly with no Valenx-side dependency.
//! `.npy` per-array beats `.npz` (zip-of-arrays) for the tight-
//! scope first cut: no zip dep, callers can group files however
//! they want, numpy reads each one with a one-liner.
//!
//! Sketched but not implemented (per RFC 0012):
//!
//! - `.npz` zip wrapping multiple `.npy` arrays into one file.
//! - JSON-manifest schema with provenance + train/val/test split
//!   (needs the sweep harness to land first).
//! - Cross-sample mesh interpolation for geometry sweeps.
//! - Per-framework converters (`.pt`, `.tfrecord`).

#![forbid(unsafe_code)]
#![allow(missing_docs)]

use thiserror::Error;

pub mod csv;
pub mod html;
pub mod manifest;
pub mod npy;
pub mod remap;
pub mod sweep;

/// Canonical I/O error type for the export crate. Lives in `lib.rs`
/// because every submodule's writer returns it.
#[derive(Debug, Error)]
pub enum ExportError {
    #[error("export to {path}: {source}")]
    Io {
        path: std::path::PathBuf,
        #[source]
        source: std::io::Error,
    },
}

// Re-exports so the existing `valenx_export::<name>` public surface
// keeps working unchanged after the split.

pub use crate::csv::{format_csv_row, write_scalars_csv};
pub use crate::html::{
    render_html_report, render_markdown_report, write_html_report, write_markdown_report,
};
pub use crate::manifest::{write_manifest, ExportManifest, ExportSchemaSection};
pub use crate::npy::{
    build_npy_bytes_f64, build_npy_bytes_f64_nd, write_npy_f64, write_npy_f64_nd, write_npz_f64,
    NpzEntry,
};
pub use crate::remap::{remap_sample_fields, FieldSample};
pub use crate::sweep::{export_sweep_dataset, DatasetExportConfig, DatasetExportError, Sample};
