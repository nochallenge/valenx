//! RON-based persistence for sketches.

use std::path::Path;

use serde::{Deserialize, Serialize};

use crate::sketch::Sketch;
use crate::SketchError;

/// On-disk envelope wrapping a sketch with format-version metadata.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SketchFile {
    /// Format version — bumped when on-disk schema changes.
    pub version: u32,
    /// The sketch payload.
    pub sketch: Sketch,
}

impl SketchFile {
    /// Current on-disk format version.
    pub const VERSION: u32 = 1;

    /// Wrap a sketch.
    pub fn from_sketch(sketch: &Sketch) -> Self {
        Self {
            version: Self::VERSION,
            sketch: sketch.clone(),
        }
    }

    /// Serialize to a RON string.
    pub fn to_ron(&self) -> Result<String, SketchError> {
        ron::ser::to_string_pretty(self, ron::ser::PrettyConfig::new())
            .map_err(|e| SketchError::Ron(e.to_string()))
    }

    /// Write to a file. Round-28 H2: routed through the canonical
    /// `valenx_core::io_caps::atomic_write_str` (sidecar O_NOFOLLOW,
    /// fsync-before-rename, parent-dir fsync on Unix). Pre-fix this
    /// was a bare `std::fs::write` which silently followed leaf
    /// symlinks and was non-atomic.
    pub fn write_to(&self, path: &Path) -> Result<(), SketchError> {
        let ron = self.to_ron()?;
        valenx_core::io_caps::atomic_write_str(path, &ron)?;
        Ok(())
    }

    /// Parse from a RON string.
    ///
    /// R33 H1: after structural deserialization the embedded sketch is
    /// run through [`Sketch::validate`] so a hand-edited / corrupt /
    /// version-skewed document carrying an out-of-range variable handle
    /// is rejected here, rather than panicking with "index out of
    /// bounds" the first time the sketch is consumed (solver, extrude,
    /// feature-tree replay).
    pub fn from_ron(s: &str) -> Result<Self, SketchError> {
        let file: Self = ron::from_str(s).map_err(|e| SketchError::Ron(e.to_string()))?;
        file.sketch.validate()?;
        Ok(file)
    }

    /// Read from a file.
    ///
    /// R29 D: bounded at [`valenx_core::io_caps::MAX_DOC_FILE_BYTES`]
    /// (16 MiB) via the canonical helper, replacing the per-crate
    /// private `read_capped_to_string` duplicate.
    pub fn read_from(path: &Path) -> Result<Self, SketchError> {
        let s = valenx_core::io_caps::read_capped_to_string(
            path,
            valenx_core::io_caps::MAX_DOC_FILE_BYTES,
        )?;
        Self::from_ron(&s)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn round_trips_empty_sketch() {
        let s = Sketch::new();
        let f = SketchFile::from_sketch(&s);
        let ron = f.to_ron().unwrap();
        assert!(ron.contains("version: 1"));
    }

    #[test]
    fn writes_to_file() {
        let s = Sketch::new();
        let f = SketchFile::from_sketch(&s);
        let tmp = std::env::temp_dir().join("valenx_sketch_test.ron");
        f.write_to(&tmp).unwrap();
        assert!(tmp.exists());
        let _ = std::fs::remove_file(&tmp);
    }

    #[test]
    fn round_trip_one_point() {
        let mut s = Sketch::new();
        s.add_point(3.0, 4.0);
        let f = SketchFile::from_sketch(&s);
        let ron = f.to_ron().unwrap();
        let parsed = SketchFile::from_ron(&ron).unwrap();
        assert_eq!(parsed.version, 1);
        assert_eq!(parsed.sketch.vars, vec![3.0, 4.0]);
    }

    #[test]
    fn round_trips_with_constraint() {
        use crate::constraint::Constraint;
        use crate::geom::EntityId;
        let mut s = Sketch::new();
        let a = s.add_point(0.0, 0.0);
        let b = s.add_point(3.0, 4.0);
        s.add_constraint(Constraint::Distance { a, b, target: 5.0 });
        let ron = SketchFile::from_sketch(&s).to_ron().unwrap();
        let parsed = SketchFile::from_ron(&ron).unwrap();
        assert_eq!(parsed.sketch.constraints.len(), 1);
        match &parsed.sketch.constraints[0] {
            Constraint::Distance {
                a: pa,
                b: pb,
                target,
            } => {
                assert_eq!(*pa, EntityId(1));
                assert_eq!(*pb, EntityId(2));
                assert_eq!(*target, 5.0);
            }
            other => panic!("wrong constraint variant: {other:?}"),
        }
    }

    /// Round-12 M1: a file larger than MAX_DOC_FILE_BYTES must be
    /// rejected at the read-cap layer rather than slurped into
    /// memory and then handed to the RON parser. This is the
    /// representative test for the 11-crate persist.rs sweep —
    /// each crate has its own copy of the cap helper, so a single
    /// regression test here pins the contract for the family.
    #[test]
    fn read_from_rejects_oversize_file() {
        let tmp = std::env::temp_dir().join("valenx_sketch_oversize_test.ron");
        // 17 MiB of zeros — over the 16 MiB cap, well past anything
        // a realistic sketch could ever produce.
        let oversize = valenx_core::io_caps::MAX_DOC_FILE_BYTES + 1024;
        std::fs::write(&tmp, vec![0u8; oversize]).unwrap();
        let err = SketchFile::read_from(&tmp).expect_err("must reject oversize file");
        match err {
            SketchError::Io(io_err) => {
                assert_eq!(io_err.kind(), std::io::ErrorKind::InvalidData);
            }
            other => panic!("expected Io(InvalidData), got: {other:?}"),
        }
        let _ = std::fs::remove_file(&tmp);
    }

    /// R33 H1: a standalone sketch RON whose entity carries a variable
    /// handle past the end of `vars` must be rejected by `from_ron`
    /// (sketch.corrupt_handle), not parsed and later panicked on.
    #[test]
    fn from_ron_rejects_out_of_range_var_handle() {
        let ron = r#"(
    version: 1,
    sketch: (
        vars: [0.0],
        entities: [
            Point((x_var: 999, y_var: 0)),
        ],
        constraints: [],
    ),
)"#;
        let err =
            SketchFile::from_ron(ron).expect_err("corrupt handle must be rejected at sketch load");
        assert_eq!(err.code(), "sketch.corrupt_handle");
    }

    /// Phase 12A Task 7: persistence regression covering the new
    /// BSpline / Ellipse / EllipticalArc primitives.
    #[test]
    fn round_trip_with_bspline_ellipse_and_arc() {
        let mut s = Sketch::new();
        let c = s.add_point(0.0, 0.0);
        let _e = s.add_ellipse(c, (2.0, 0.0), 1.0).unwrap();
        let _ea = s
            .add_elliptical_arc(c, (2.0, 0.0), 1.0, 0.0, std::f64::consts::PI)
            .unwrap();
        let p0 = s.add_point(0.0, 0.0);
        let p1 = s.add_point(1.0, 1.0);
        let p2 = s.add_point(2.0, 1.0);
        let p3 = s.add_point(3.0, 0.0);
        let _b = s
            .add_bspline(
                3,
                vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
                &[p0, p1, p2, p3],
                vec![1.0; 4],
            )
            .unwrap();
        let ron = SketchFile::from_sketch(&s).to_ron().unwrap();
        let parsed = SketchFile::from_ron(&ron).unwrap();
        assert_eq!(parsed.sketch.entities.len(), s.entities.len());
    }
}
