//! Zero-copy memory-mapped reads for very large files.
//!
//! The `io_caps` module enforces a *stat-and-bounded-take* pattern that
//! rejects multi-GB payloads before they are slurped into RAM — the
//! right default for untrusted inputs that must be parsed in full. This
//! module covers the complementary case: a file that is *legitimately*
//! multi-GB (a CAD tessellation dump, a FASTQ/BAM sequencing run, a
//! telemetry capture) and that we want to read **without** allocating a
//! heap buffer the size of the file. Memory-mapping lets the OS page the
//! bytes in lazily and back the [`&[u8]`] view directly with the file's
//! pages, so a random-access reader touches only the pages it needs.
//!
//! Only a read-only, private mapping is exposed; nothing here writes
//! through the map. The single audited `unsafe` site is
//! [`MappedFile::open`], where `memmap2::Mmap::map` is called.
//!
//! # Safety contract
//!
//! Memory-mapping is `unsafe` in `memmap2` for one reason: the borrow
//! checker cannot police the file's *backing storage*, which lives
//! outside the process. If another process (or another handle in this
//! process) **truncates or shrinks** the file while it is mapped, then
//! touching a page past the new end of file faults with `SIGBUS`
//! (Unix) / an access violation (Windows) — undefined behaviour from
//! Rust's point of view, and not something a `Result` can catch.
//!
//! By opening the file read-only and keeping the [`File`] handle alive
//! inside [`MappedFile`] for the lifetime of the mapping, this type
//! removes the *append/grow* and *use-after-free-of-fd* hazards. The
//! remaining obligation is on the **caller / operator**:
//!
//! > While a [`MappedFile`] is alive, the underlying file MUST NOT be
//! > truncated, shrunk, or otherwise have bytes within `[0, len())`
//! > made unreadable by any other process or handle.
//!
//! Appending to the file is harmless (the mapping simply does not see
//! the new tail); only making already-mapped bytes disappear is unsound.
//! For files under our own control (we created them, or they are
//! immutable inputs in a read-only case directory) this holds trivially.

use std::fs::File;
use std::path::Path;

use crate::error::AdapterError;

/// A read-only memory-mapped view of a file on disk.
///
/// Construct one with [`MappedFile::open`]; read its contents as a byte
/// slice with [`as_slice`](MappedFile::as_slice), a bounds-checked
/// sub-range with [`chunk`](MappedFile::chunk), or the mapped length
/// with [`len`](MappedFile::len). The mapping (and the file handle that
/// backs it) is released when the value is dropped.
///
/// See the [module docs](self) for the safety contract: the file must
/// not be truncated while this value is alive.
#[derive(Debug)]
pub struct MappedFile {
    // Field order matters for Drop: `mmap` is dropped before `_file`,
    // so the mapping is torn down while its backing fd is still open.
    mmap: memmap2::Mmap,
    // Held purely to keep the descriptor alive for the lifetime of the
    // mapping; never read after construction.
    _file: File,
}

impl MappedFile {
    /// Open `path` read-only and memory-map its entire contents.
    ///
    /// An empty file maps to a zero-length view (no bytes, no error):
    /// [`is_empty`](MappedFile::is_empty) returns `true` and
    /// [`as_slice`](MappedFile::as_slice) returns `&[]`.
    ///
    /// # Errors
    ///
    /// Returns [`AdapterError::Io`] if the path cannot be opened (missing,
    /// permission denied, is a directory, …) or if the OS refuses the
    /// mapping.
    ///
    /// # Safety
    ///
    /// This is a safe function, but it relies on the module-level safety
    /// contract: the caller must ensure the file is **not truncated or
    /// shrunk** by any other process or handle while the returned
    /// [`MappedFile`] is alive. Violating that can fault on page access
    /// (`SIGBUS` / access violation), which is undefined behaviour.
    pub fn open(path: impl AsRef<Path>) -> Result<Self, AdapterError> {
        let path = path.as_ref();
        let file = File::open(path).map_err(|e| {
            AdapterError::Io(std::io::Error::new(
                e.kind(),
                format!("mmap: cannot open {}: {e}", path.display()),
            ))
        })?;
        // SAFETY: `file` is opened read-only just above and is moved into
        // the returned struct (kept alive for the whole lifetime of the
        // mapping). The only remaining requirement memmap2 places on the
        // caller — that the file is not shrunk/truncated underneath the
        // mapping — is documented as this type's safety contract (see the
        // module docs and `open`'s `# Safety` section) and is the
        // operator's responsibility. We never write through the mapping.
        #[allow(unsafe_code)]
        let mmap = unsafe { memmap2::Mmap::map(&file) }.map_err(|e| {
            AdapterError::Io(std::io::Error::new(
                e.kind(),
                format!("mmap: cannot map {}: {e}", path.display()),
            ))
        })?;
        Ok(Self { mmap, _file: file })
    }

    /// The number of mapped bytes (the file's length at map time).
    #[inline]
    pub fn len(&self) -> usize {
        self.mmap.len()
    }

    /// Whether the mapping is empty (the file was zero-length).
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.mmap.is_empty()
    }

    /// The full mapped contents as a byte slice.
    ///
    /// The slice borrows the mapping, so it cannot outlive `self`. For a
    /// zero-length file this is `&[]`.
    #[inline]
    pub fn as_slice(&self) -> &[u8] {
        &self.mmap
    }

    /// A bounds-checked sub-range `[offset, offset + len)` of the mapping.
    ///
    /// Returns `None` (never panics) if the range is out of bounds —
    /// i.e. if `offset > self.len()` or `offset + len > self.len()`, or
    /// if `offset + len` overflows `usize`. A zero-length request at any
    /// `offset <= self.len()` yields `Some(&[])`.
    #[inline]
    pub fn chunk(&self, offset: usize, len: usize) -> Option<&[u8]> {
        let end = offset.checked_add(len)?;
        self.mmap.get(offset..end)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    /// Write `bytes` to a uniquely-named temp file and return its path.
    /// We use the process id plus a monotonic counter so parallel test
    /// threads never collide, and we do not rely on any external
    /// tempfile crate (keeping `valenx-core`'s dep set lean).
    fn temp_file_with(bytes: &[u8], tag: &str) -> std::path::PathBuf {
        use std::sync::atomic::{AtomicU64, Ordering};
        static COUNTER: AtomicU64 = AtomicU64::new(0);
        let n = COUNTER.fetch_add(1, Ordering::Relaxed);
        let mut path = std::env::temp_dir();
        path.push(format!("valenx-mmap-{tag}-{}-{n}.bin", std::process::id()));
        let mut f = File::create(&path).expect("create temp file");
        f.write_all(bytes).expect("write temp file");
        f.flush().expect("flush temp file");
        path
    }

    #[test]
    fn round_trips_bytes() {
        let data: Vec<u8> = (0..=255u8).cycle().take(4096).collect();
        let path = temp_file_with(&data, "roundtrip");
        let m = MappedFile::open(&path).expect("open mmap");
        assert_eq!(m.len(), data.len());
        assert!(!m.is_empty());
        assert_eq!(m.as_slice(), &data[..]);
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn sub_range_matches() {
        let data: Vec<u8> = (0..100u8).collect();
        let path = temp_file_with(&data, "subrange");
        let m = MappedFile::open(&path).expect("open mmap");
        // A middle slice.
        assert_eq!(m.chunk(10, 20), Some(&data[10..30]));
        // A slice that ends exactly at EOF.
        assert_eq!(m.chunk(90, 10), Some(&data[90..100]));
        // A zero-length slice at EOF is valid and empty.
        let empty: &[u8] = &[];
        assert_eq!(m.chunk(100, 0), Some(empty));
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn out_of_bounds_chunk_is_none_not_panic() {
        let data: Vec<u8> = (0..10u8).collect();
        let path = temp_file_with(&data, "oob");
        let m = MappedFile::open(&path).expect("open mmap");
        // Range extends past EOF.
        assert_eq!(m.chunk(5, 10), None);
        // Offset itself past EOF.
        assert_eq!(m.chunk(11, 0), None);
        // Overflowing offset + len must not panic; returns None.
        assert_eq!(m.chunk(usize::MAX, 1), None);
        assert_eq!(m.chunk(1, usize::MAX), None);
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn empty_file_is_handled() {
        let path = temp_file_with(&[], "empty");
        let m = MappedFile::open(&path).expect("open mmap of empty file");
        assert_eq!(m.len(), 0);
        assert!(m.is_empty());
        let empty: &[u8] = &[];
        assert_eq!(m.as_slice(), empty);
        // Zero-length chunk at offset 0 is the only in-bounds request.
        assert_eq!(m.chunk(0, 0), Some(empty));
        assert_eq!(m.chunk(0, 1), None);
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn missing_path_errors() {
        let mut path = std::env::temp_dir();
        path.push(format!(
            "valenx-mmap-does-not-exist-{}.bin",
            std::process::id()
        ));
        let _ = std::fs::remove_file(&path);
        assert!(MappedFile::open(&path).is_err());
    }
}
