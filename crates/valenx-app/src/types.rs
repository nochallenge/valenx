//! Plain-data types used across the app shell.
//!
//! These structs and enums were extracted from `lib.rs` to keep the
//! public-facing surface area small and to reduce the size of the
//! root module. They're re-exported from `lib.rs` so the public API
//! is unchanged: `use valenx_app::LoadedStl;` still resolves.
//!
//! Nothing here owns behaviour — the operating logic lives in the
//! `impl ValenxApp { ... }` blocks back in `lib.rs`.

use std::path::PathBuf;

use valenx_mesh::{Mesh, QualityReport};
use valenx_viz::TriangleMesh;

/// Which tab is visible in the bottom dock panel.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum BottomTab {
    #[default]
    Residuals,
    Log,
}

/// An STL file the user dropped into the viewport — source path plus
/// the parsed triangle mesh.
pub struct LoadedStl {
    /// Disk path the STL was loaded from.
    pub path: PathBuf,
    /// Parsed triangle mesh.
    pub mesh: TriangleMesh,
}

/// One entry in `ValenxApp::run_history` — the outcome of the
/// most recent run for a given case. Lightweight on purpose; the
/// full RunReport / Results live in `last_run_report` /
/// `last_run_results`, which only carry the LAST run's data. The
/// history map is the "I ran this case ten minutes ago and it
/// converged" memory the case browser needs to show a tick mark.
///
/// Serializable so the map can be persisted to
/// `<state_dir>/run-history.json` and survive app restarts.
#[derive(Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct RunHistoryEntry {
    /// Whether the run finished cleanly (exit code 0).
    pub succeeded: bool,
    /// Wall time the run took.
    pub wall_time: std::time::Duration,
    /// `Some(true)` if residuals dropped below the target,
    /// `Some(false)` if not, `None` for transient runs (no notion
    /// of convergence). Mirrors `RunReport.converged`.
    pub converged: Option<bool>,
}

/// One entry in `ValenxApp::sweep_history` — the outcome of the
/// most recent sweep for a given case. Recorded by both the sync
/// and async sweep runners so the case browser can show "swept
/// 32 cases (24 succeeded) 5 minutes ago" without needing to
/// keep the full per-derived-case state in memory.
///
/// Serializable so the map persists to
/// `<state_dir>/sweep-history.json` across app restarts.
#[derive(Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct SweepHistoryEntry {
    /// Total derived cases the sweep planned.
    pub planned: usize,
    /// Number that finished cleanly (exit 0 + Results::collect
    /// succeeded).
    pub succeeded: usize,
    /// Number that failed any pipeline stage.
    pub failed: usize,
    /// Parent sweep workdir — useful for the "open in file
    /// browser" affordance after restart.
    pub workdir: PathBuf,
    /// ISO 8601 UTC timestamp of when the sweep finished.
    pub completed_at: String,
}

/// A canonical `valenx_mesh::Mesh` loaded from disk, plus a
/// pre-computed quality report so the browser pane can render
/// stats without recomputing every frame.
pub struct LoadedMesh {
    pub path: PathBuf,
    pub mesh: Mesh,
    pub quality: QualityReport,
    /// Aspect-ratio histogram on the default buckets — computed
    /// once at load time, walked by the Quality panel for the
    /// "distribution" section.
    pub aspect_hist: valenx_mesh::AspectRatioHistogram,
    /// Skewness histogram on the default quality bands.
    pub skew_hist: valenx_mesh::SkewnessHistogram,
}
