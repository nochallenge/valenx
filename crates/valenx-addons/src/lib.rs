//! # valenx-addons
//!
//! In-app add-on manager for Valenx.
//!
//! Lets users install community-authored workbenches at runtime without
//! recompiling Valenx. Each add-on is a directory under
//! `~/.valenx/addons/{name}/` with a manifest (`valenx-addon.toml`) and
//! either a Python entry-point (shipped with v1) or a WASM module
//! (deferred). The Python entry-point gets imported into the Phase 11
//! `valenx-py` interpreter on next session start.
//!
//! ## v1 deliberate limitations
//!
//! - **No GitHub-search wiring.** The plan called for a `reqwest`-
//!   backed search of `topic:valenx-addon` repos. Pulling `reqwest`
//!   in across every workspace target (Windows MSVC, macOS, Linux,
//!   the in-tree MSRV pin) needed multi-platform verification that
//!   wasn't on the Phase 22 budget. v1 ships **manual install by
//!   directory** (`addons::install_from_dir`) + a documented surface
//!   the GitHub-search path can plug into in a follow-up phase.
//! - **WASM entry-points deferred.** The manifest already has a
//!   `entry_point` enum that holds a `Wasm` variant for future use;
//!   v1 only dispatches `Python` entry-points to the embedded
//!   interpreter.
//! - **No version pinning.** v1 lets the user upgrade or downgrade
//!   manually by re-installing; semver-aware update is a v1.5
//!   follow-up.
//!
//! ## Quick start
//!
//! ```no_run
//! use std::path::PathBuf;
//! use valenx_addons::AddonRegistry;
//!
//! let mut registry = AddonRegistry::new(PathBuf::from("/tmp/addons"));
//! let _ = registry.refresh();
//! for a in registry.installed() {
//!     println!("{} v{}", a.manifest.name, a.manifest.version);
//! }
//! ```

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod entrypoint;
pub mod error;
pub mod install;
pub mod manifest;
pub mod registry;

pub use entrypoint::EntryPoint;
pub use error::AddonError;
pub use install::{install_from_dir, uninstall, update_from_dir};
pub use manifest::{parse_manifest, AddonManifest};
pub use registry::{AddonRegistry, LocalAddon, RemoteAddon};
