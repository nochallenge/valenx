//! # valenx-partlib
//!
//! Offline parts library — installs CAD files (STEP / IGES / STL /
//! NURBS) into a local root directory, computes a SHA-256 checksum at
//! install time, and persists a [`PartLibrary`] index file
//! (`index.ron`).
//!
//! Phase 46 of the FreeCAD-parity roadmap.  FreeCAD `PartLibrary`
//! community workbench equivalent.
//!
//! # Surface
//!
//! - [`PartEntry`] / [`PartKind`] — single library item.
//! - [`PartLibrary`] — name → entry HashMap rooted at a directory.
//! - [`load_index`] — read `<root>/index.ron` (empty on missing).
//! - [`save_index`] — write `<root>/index.ron`.
//! - [`install_local`] — copy a file into the library + register.
//! - [`fetch_remote`] — v1 stub returning
//!   [`PartLibError::FetchRequiresNetwork`]; v2 will route through
//!   the Phase 22 Add-on Manager pipeline.
//! - [`PartLibPanelState`] — UI state for the Valenx panel.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod entry;
pub mod error;
pub mod library;
pub mod panel;

pub use entry::{PartEntry, PartKind};
pub use error::{ErrorCategory, PartLibError};
pub use library::{
    fetch_remote, install_local, load_index, save_index, LibraryFile, PartLibrary,
    INDEX_FILENAME, VERSION,
};
pub use panel::PartLibPanelState;
