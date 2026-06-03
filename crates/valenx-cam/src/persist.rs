//! RON envelope for round-tripping CAM workbench files.
//!
//! Same pattern as `valenx-surface::persist::SurfaceFile`,
//! `valenx-sketch::persist::SketchFile`, etc.: a thin envelope
//! wrapping the live in-memory state with a `version` field so we
//! can evolve the schema without breaking older files.
//!
//! The generated [`crate::Toolpath`] is **not** persisted — it's
//! regenerated on demand from the tools + stock + operations.

use std::path::Path;

use serde::{Deserialize, Serialize};

use crate::error::CamError;
use crate::fixture::Fixture;
use crate::operation::Operation;
use crate::setup::SetupSet;
use crate::stock::Stock;
use crate::tool::Tool;

/// On-disk envelope wrapping a CAM-workbench session.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct CamFile {
    /// Format version — bumped when the on-disk schema changes.
    pub version: u32,
    /// Stock block in this session.
    pub stock: Stock,
    /// All tools in the session's tool table.
    pub tools: Vec<Tool>,
    /// Ordered list of operations the user has set up.
    pub operations: Vec<Operation>,
    /// Phase 17F — optional fixture geometry. Empty by default for
    /// backwards-compat with Phase 10 files.
    #[serde(default)]
    pub fixture: Fixture,
    /// Phase 17F — optional multi-setup list. Empty by default.
    #[serde(default)]
    pub setups: SetupSet,
}

impl CamFile {
    /// Current on-disk format version. Bump and add a migration
    /// step in `from_ron` when the schema changes.
    pub const VERSION: u32 = 2;

    /// Construct an empty CAM file at the current `VERSION`.
    pub fn new() -> Self {
        Self {
            version: Self::VERSION,
            stock: Stock::default(),
            tools: Vec::new(),
            operations: Vec::new(),
            fixture: Fixture::new(),
            setups: SetupSet::new(),
        }
    }

    /// Serialize to a pretty-printed RON string.
    pub fn to_ron(&self) -> Result<String, CamError> {
        ron::ser::to_string_pretty(self, ron::ser::PrettyConfig::new())
            .map_err(|e| CamError::Ron(e.to_string()))
    }

    /// Write to a file. Overwrites if the file exists. Round-28 H2:
    /// goes through the canonical
    /// `valenx_core::io_caps::atomic_write_str` for sidecar-based
    /// atomic publication (O_NOFOLLOW, fsync-before-rename, parent
    /// dir fsync on Unix). Pre-fix this was a bare `std::fs::write`
    /// which silently followed leaf symlinks and was non-atomic.
    pub fn write_to(&self, path: &Path) -> Result<(), CamError> {
        valenx_core::io_caps::atomic_write_str(path, &self.to_ron()?)?;
        Ok(())
    }

    /// Parse from a RON string.
    pub fn from_ron(s: &str) -> Result<Self, CamError> {
        ron::from_str(s).map_err(|e| CamError::Ron(e.to_string()))
    }

    /// Read from a file.
    ///
    /// R29 D: bounded at [`valenx_core::io_caps::MAX_DOC_FILE_BYTES`]
    /// (16 MiB) via the canonical helper, replacing the per-crate
    /// private `read_capped_to_string` duplicate.
    pub fn read_from(path: &Path) -> Result<Self, CamError> {
        Self::from_ron(&valenx_core::io_caps::read_capped_to_string(
            path,
            valenx_core::io_caps::MAX_DOC_FILE_BYTES,
        )?)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::op::adaptive_clearing::AdaptiveParams;
    use crate::operation::Operation;

    #[test]
    fn round_trip_v2_includes_phase17_ops() {
        let mut f = CamFile::new();
        assert_eq!(f.version, CamFile::VERSION);
        f.operations
            .push(Operation::AdaptiveClearing(AdaptiveParams::default()));
        let ron = f.to_ron().unwrap();
        let g = CamFile::from_ron(&ron).unwrap();
        assert_eq!(g.operations.len(), 1);
        assert_eq!(g.operations[0].label(), "Adaptive Clearing");
        assert_eq!(g.version, CamFile::VERSION);
    }

    #[test]
    fn round_trip_fixture_and_setups_default_when_absent() {
        // A Phase 10 file that lacks fixture / setups deserialises
        // cleanly thanks to `#[serde(default)]`.
        let v1_ron = "(version: 1, stock: (origin: (0.0, 0.0, 0.0), size: (10.0, 10.0, 10.0), material: \"al\"), tools: [], operations: [])";
        let f = CamFile::from_ron(v1_ron).unwrap();
        assert!(f.fixture.aabbs.is_empty());
        assert!(f.setups.is_empty());
    }
}
