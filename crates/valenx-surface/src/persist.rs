//! RON envelope for round-tripping surface workbench files.
//!
//! Same pattern as `valenx-sketch::persist::SketchFile`,
//! `valenx-assembly::persist::AssemblyFile`, etc.: a thin envelope
//! wrapping the live in-memory state with a `version` field so we
//! can evolve the schema without breaking older files.

use std::path::Path;

use serde::{Deserialize, Serialize};

use crate::error::SurfaceError;
use crate::nurbs_curve::NurbsCurve;
use crate::nurbs_surface::NurbsSurface;

/// On-disk envelope wrapping a surface-workbench session.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct SurfaceFile {
    /// Format version — bumped when the on-disk schema changes.
    pub version: u32,
    /// All NURBS curves in the session.
    pub curves: Vec<NurbsCurve>,
    /// All NURBS surfaces in the session.
    pub surfaces: Vec<NurbsSurface>,
}

impl SurfaceFile {
    /// Current on-disk format version. Bump and add a migration
    /// step in `from_ron` when the schema changes.
    ///
    /// - **v1** (Phase 9): curves + surfaces with rational CP grids.
    /// - **v2** (Phase 19): same wire format — bumped to flag files
    ///   that *may* contain surfaces produced by knot insertion /
    ///   removal / degree elevation / G2 sew / SSI curve fits /
    ///   ruled surfaces / point-cloud fits. v1 files still load
    ///   transparently (we don't gate any field on the version);
    ///   v2 just lets newer readers display a "produced with Phase
    ///   19 features" hint when relevant.
    pub const VERSION: u32 = 2;

    /// Construct an empty surface file at the current `VERSION`.
    pub fn new() -> Self {
        Self {
            version: Self::VERSION,
            curves: Vec::new(),
            surfaces: Vec::new(),
        }
    }

    /// Serialize to a pretty-printed RON string.
    pub fn to_ron(&self) -> Result<String, SurfaceError> {
        ron::ser::to_string_pretty(self, ron::ser::PrettyConfig::new())
            .map_err(|e| SurfaceError::Ron(e.to_string()))
    }

    /// Write to a file. Overwrites if the file exists. Round-28 H2:
    /// routed through the canonical
    /// `valenx_core::io_caps::atomic_write_str` (sidecar O_NOFOLLOW,
    /// fsync-before-rename, parent-dir fsync on Unix). Pre-fix this
    /// was a bare `std::fs::write` which silently followed leaf
    /// symlinks and was non-atomic.
    pub fn write_to(&self, path: &Path) -> Result<(), SurfaceError> {
        valenx_core::io_caps::atomic_write_str(path, &self.to_ron()?)?;
        Ok(())
    }

    /// Parse from a RON string.
    pub fn from_ron(s: &str) -> Result<Self, SurfaceError> {
        ron::from_str(s).map_err(|e| SurfaceError::Ron(e.to_string()))
    }

    /// Read from a file.
    pub fn read_from(path: &Path) -> Result<Self, SurfaceError> {
        // R29 D: canonical valenx_core::io_caps::read_capped_to_string at
        // MAX_DOC_FILE_BYTES (16 MiB), replacing the private dupe.
        Self::from_ron(&valenx_core::io_caps::read_capped_to_string(
            path,
            valenx_core::io_caps::MAX_DOC_FILE_BYTES,
        )?)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Vector3;

    #[test]
    fn round_trips_empty() {
        let f = SurfaceFile::new();
        let ron = f.to_ron().unwrap();
        let parsed = SurfaceFile::from_ron(&ron).unwrap();
        assert_eq!(parsed.version, SurfaceFile::VERSION);
        assert!(parsed.curves.is_empty());
        assert!(parsed.surfaces.is_empty());
    }

    #[test]
    fn v1_files_load_transparently() {
        // A v1-style file (just version: 1) should still round-trip:
        // the schema for curves + surfaces hasn't changed in Phase 19,
        // so we don't even need a migration.
        let ron = r#"(version: 1, curves: [], surfaces: [])"#;
        let parsed = SurfaceFile::from_ron(ron).unwrap();
        assert_eq!(parsed.version, 1);
        assert!(parsed.curves.is_empty());
        assert!(parsed.surfaces.is_empty());
    }

    #[test]
    fn round_trips_phase19_artefacts() {
        // Build a curve via knot insertion, a ruled surface, and a
        // G2-stitched surface; persist them and read back.
        let knots = vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0];
        let cps = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 2.0, 0.0),
            Vector3::new(2.0, 2.0, 0.0),
            Vector3::new(3.0, 0.0, 0.0),
        ];
        let curve = NurbsCurve::new(3, knots, cps, vec![1.0; 4]).unwrap();
        let inserted = curve.insert_knot(0.5).unwrap();
        let ruled =
            crate::ruled::extrude_along_vector(&curve, Vector3::new(0.0, 0.0, 1.0)).unwrap();
        let mut f = SurfaceFile::new();
        f.curves.push(curve);
        f.curves.push(inserted);
        f.surfaces.push(ruled);
        let ron = f.to_ron().unwrap();
        let parsed = SurfaceFile::from_ron(&ron).unwrap();
        assert_eq!(parsed.curves.len(), 2);
        assert_eq!(parsed.surfaces.len(), 1);
        // Inserted curve has one more CP than the original.
        assert_eq!(
            parsed.curves[1].n_control_points(),
            parsed.curves[0].n_control_points() + 1
        );
    }

    #[test]
    fn round_trips_one_curve() {
        let mut f = SurfaceFile::new();
        f.curves.push(
            NurbsCurve::new(
                3,
                vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
                vec![
                    Vector3::new(0.0, 0.0, 0.0),
                    Vector3::new(1.0, 1.0, 0.0),
                    Vector3::new(2.0, 1.0, 0.0),
                    Vector3::new(3.0, 0.0, 0.0),
                ],
                vec![1.0; 4],
            )
            .unwrap(),
        );
        let ron = f.to_ron().unwrap();
        let parsed = SurfaceFile::from_ron(&ron).unwrap();
        assert_eq!(parsed.curves.len(), 1);
        assert_eq!(parsed.curves[0].degree, 3);
        assert_eq!(parsed.curves[0].n_control_points(), 4);
    }

    #[test]
    fn writes_to_file_round_trips() {
        let mut f = SurfaceFile::new();
        f.curves.push(
            NurbsCurve::new(
                3,
                vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
                vec![Vector3::zeros(); 4],
                vec![1.0; 4],
            )
            .unwrap(),
        );
        let tmp = std::env::temp_dir().join("valenx_surface_persist_test.ron");
        f.write_to(&tmp).unwrap();
        let parsed = SurfaceFile::read_from(&tmp).unwrap();
        assert_eq!(parsed.curves.len(), 1);
        let _ = std::fs::remove_file(&tmp);
    }
}
