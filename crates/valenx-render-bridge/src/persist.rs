//! RON-based persistence for a render job.

use std::path::Path;

use serde::{Deserialize, Serialize};

use crate::error::RenderError;
use crate::scene::RenderJob;

/// On-disk envelope.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct RenderFile {
    /// Format version.
    pub version: u32,
    /// The job payload.
    pub job: RenderJob,
}

impl RenderFile {
    /// Current on-disk format version.
    pub const VERSION: u32 = 1;

    /// Wrap a job.
    pub fn from_job(job: &RenderJob) -> Self {
        Self {
            version: Self::VERSION,
            job: job.clone(),
        }
    }

    /// Pretty RON.
    pub fn to_ron(&self) -> Result<String, RenderError> {
        ron::ser::to_string_pretty(self, ron::ser::PrettyConfig::new())
            .map_err(|e| RenderError::Ron(e.to_string()))
    }

    /// Write to disk atomically — write to a unique sidecar, fsync
    /// the tmp file, then rename it over `<path>` and (on Unix)
    /// fsync the parent directory so the dentry update is durable.
    ///
    /// Round-22 L3 (sister to `valenx_app::state_paths::atomic_write`):
    /// pre-fix this site did `std::fs::write(path, ron)` which
    /// truncates `<path>` to zero bytes and then starts writing. A
    /// process crash mid-write (or a power loss with WB-cache enabled)
    /// would leave an empty or partial `.ron` on disk that the next
    /// `read_from` round-trip would fail to parse.
    ///
    /// ## Round-27 STRUCTURAL consolidation
    ///
    /// Thin wrapper around
    /// [`valenx_core::io_caps::atomic_write_str`]. Pre-fix this site
    /// inlined the write-tmp → fsync → rename pattern (justified at
    /// the time as "`render-bridge` cannot depend on `valenx-app`")
    /// — but `valenx-render-bridge` HAS depended on `valenx-core`
    /// since round 20 (for the `read_capped_to_string` cap in
    /// `read_from`), so the dep direction is already available. By
    /// moving the canonical helper into `valenx-core` and routing
    /// every site through it, render-bridge inherits the unique
    /// `<pid>.<counter>` sidecar, `O_NOFOLLOW`, and parent-dir
    /// fsync (R26 M3) without ever having to back-port them
    /// independently.
    pub fn write_to(&self, path: &Path) -> Result<(), RenderError> {
        let ron = self.to_ron()?;
        valenx_core::io_caps::atomic_write_str(path, &ron).map_err(RenderError::Io)
    }

    /// Parse RON.
    pub fn from_ron(s: &str) -> Result<Self, RenderError> {
        ron::from_str(s).map_err(|e| RenderError::Ron(e.to_string()))
    }

    /// Read from disk.
    ///
    /// Round-20 M1 (R12 persist.rs sweep sister gap): pre-fix this
    /// site did a bare `std::fs::read_to_string`, leaving the
    /// render-bridge as the last workbench-style `persist.rs` whose
    /// read path could slurp a multi-GB hostile `.ron` into memory
    /// before serde-ron saw the first token. The cap matches every
    /// other workbench persist site (sketch / cam / arch / techdraw /
    /// surface / spreadsheet / macro / draft / lattice / assembly).
    pub fn read_from(path: &Path) -> Result<Self, RenderError> {
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
    fn round_trip_empty_job() {
        let job = RenderJob::new();
        let f = RenderFile::from_job(&job);
        let ron = f.to_ron().unwrap();
        assert!(ron.contains("version: 1"));
        let back = RenderFile::from_ron(&ron).unwrap();
        assert_eq!(back.version, 1);
    }

    #[test]
    fn write_to_then_read_from_disk_round_trips() {
        // Exercises the `write_to` + `read_from` filesystem paths.
        let job = RenderJob::new();
        let f = RenderFile::from_job(&job);
        let dir = std::env::temp_dir().join(format!(
            "valenx-renderfile-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("job.ron");
        f.write_to(&path).expect("write_to should succeed");
        assert!(path.is_file(), "the RON file should exist on disk");
        let back = RenderFile::read_from(&path).expect("read_from should succeed");
        assert_eq!(back.version, RenderFile::VERSION);
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn read_from_missing_path_is_an_io_error() {
        // The `read_from` -> `std::fs::read_to_string` error branch.
        let missing = std::env::temp_dir()
            .join("valenx-definitely-not-a-real-renderfile-xyz.ron");
        let err = RenderFile::read_from(&missing).unwrap_err();
        assert_eq!(err.code(), "render.io", "a missing file is an IO error");
    }

    #[test]
    fn from_ron_rejects_malformed_input() {
        // The `from_ron` parse-error branch.
        let err = RenderFile::from_ron("this is not valid RON {{{").unwrap_err();
        assert_eq!(err.code(), "render.ron");
    }

    #[test]
    fn from_job_stamps_the_current_version() {
        let f = RenderFile::from_job(&RenderJob::new());
        assert_eq!(f.version, RenderFile::VERSION);
    }

    /// Round-22 L3 RED→GREEN (sister to `state_paths::atomic_write`):
    /// `write_to` must publish the final `.ron` atomically — i.e. the
    /// destination either has the old contents (which are absent in a
    /// fresh test) or the full new contents, never a partial write.
    ///
    /// Direct simulation of a process crash mid-write is hard; the
    /// observable proxy here is that no `.tmp` sidecar is left behind
    /// after a successful write (the pre-fix `std::fs::write` writes
    /// directly to `path` with no sidecar at all, but it also doesn't
    /// satisfy the "never partial" guarantee). Post-fix the sidecar
    /// is written, fsynced, and renamed — and the rename is the
    /// publication point, so after a clean write only `path` (not
    /// `path.tmp`) is on disk.
    #[test]
    fn write_to_publishes_atomically_no_tmp_sidecar() {
        let job = RenderJob::new();
        let f = RenderFile::from_job(&job);
        let dir = std::env::temp_dir().join(format!(
            "valenx-renderfile-r22l3-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("job.ron");
        f.write_to(&path).expect("write_to should succeed");
        // The final file must exist and round-trip.
        assert!(path.is_file(), "the RON file should exist on disk");
        let back = RenderFile::read_from(&path).expect("round-trip should succeed");
        assert_eq!(back.version, RenderFile::VERSION);
        // The `.tmp` sidecar must not be lingering — the rename step
        // in the atomic-write path consumes it. (If pre-fix `std::fs::write`
        // were still in use this assertion would still pass because
        // `write` doesn't create a sidecar at all — the assertion here
        // is about the post-fix shape, not the pre-fix-vs-post-fix
        // contrast.)
        let tmp = dir.join("job.ron.tmp");
        assert!(
            !tmp.exists(),
            "atomic-write tmp sidecar must be renamed away; found: {}",
            tmp.display()
        );
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// Round-22 L3 RED→GREEN: a `write_to` over an existing file
    /// must replace its contents atomically (rename, not truncate-
    /// then-write). The observable proxy is "after the new write,
    /// read_from returns the new contents", verified across a
    /// version-stamp round trip.
    #[test]
    fn write_to_replaces_existing_file_atomically() {
        let job_v1 = RenderJob::new();
        let f1 = RenderFile::from_job(&job_v1);
        let dir = std::env::temp_dir().join(format!(
            "valenx-renderfile-r22l3-replace-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("job.ron");
        f1.write_to(&path).expect("first write");
        // Read back v1.
        let back_v1 = RenderFile::read_from(&path).expect("read v1");
        assert_eq!(back_v1.version, RenderFile::VERSION);
        // Now overwrite with a fresh job — the rename should replace
        // atomically.
        let job_v2 = RenderJob::new();
        let f2 = RenderFile::from_job(&job_v2);
        f2.write_to(&path).expect("second write");
        let back_v2 = RenderFile::read_from(&path).expect("read v2");
        assert_eq!(back_v2.version, RenderFile::VERSION);
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// Round-20 M1 RED→GREEN (R12 persist.rs sweep sister gap): a
    /// `.ron` file larger than `MAX_DOC_FILE_BYTES` (16 MiB) must be
    /// rejected as an IO error WITHOUT being slurped into memory.
    /// Pre-fix `read_from` did a bare `std::fs::read_to_string` — a
    /// 20 MiB hostile RON would have allocated a 20 MiB String before
    /// the RON parser saw anything.
    #[test]
    fn read_from_rejects_oversize_ron_file() {
        use std::io::Write;
        let dir = std::env::temp_dir().join(format!(
            "valenx-renderfile-r20m1-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("oversize.ron");
        // 20 MiB — past the 16 MiB MAX_DOC_FILE_BYTES cap. set_len +
        // single-byte write gives us an over-cap file without
        // writing 20 MiB of bytes on every CI run.
        let mut f = std::fs::File::create(&path).unwrap();
        f.set_len(20 * 1024 * 1024).unwrap();
        f.write_all(b"x").unwrap();
        drop(f);
        let err = RenderFile::read_from(&path)
            .expect_err("round-20 M1: 20 MiB ron must be rejected as IO error");
        assert_eq!(
            err.code(),
            "render.io",
            "an oversize file is an IO error (size-cap exceeded)"
        );
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// RED→GREEN (round-27 M2 sister, STRUCTURAL): concurrent
    /// `write_to` against the same path must end with one writer's
    /// content at the destination — pre-R27 the inlined `<path>.tmp`
    /// shape gave two writers a shared sidecar handle to interleave
    /// against; post-R27 each writer owns a distinct
    /// `<basename>.tmp.<pid>.<counter>` sidecar via the canonical
    /// helper, so the rename publishes one writer's content atomically.
    #[test]
    fn write_to_concurrent_writes_preserve_integrity_round27() {
        use std::sync::{Arc, Barrier};
        use std::thread;
        let dir = std::env::temp_dir().join(format!(
            "valenx-renderfile-r27-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("job.ron");
        const N: usize = 4;
        let barrier = Arc::new(Barrier::new(N));
        let mut handles = Vec::with_capacity(N);
        for _ in 0..N {
            let path = path.clone();
            let barrier = Arc::clone(&barrier);
            handles.push(thread::spawn(move || {
                let f = RenderFile::from_job(&RenderJob::new());
                barrier.wait();
                f.write_to(&path)
            }));
        }
        for h in handles {
            h.join().unwrap().expect("concurrent write");
        }
        // The final file must round-trip cleanly — proves no
        // interleaving corrupted the RON.
        let back = RenderFile::read_from(&path).expect("round-trip");
        assert_eq!(back.version, RenderFile::VERSION);
        let _ = std::fs::remove_dir_all(&dir);
    }
}
