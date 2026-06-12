//! [`PartLibrary`] — index file + filesystem ops.

use std::collections::HashMap;
use std::fs;
use std::io::Read;
use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};

use crate::entry::{PartEntry, PartKind};
use crate::error::PartLibError;

/// Filename of the on-disk index inside the library root.
pub const INDEX_FILENAME: &str = "index.ron";

/// Versioned envelope for the index file.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct LibraryFile {
    /// Format version.
    pub version: u32,
    /// Wrapped library.
    pub library: PartLibrary,
}

/// Current version.
pub const VERSION: u32 = 1;

/// In-memory parts library — root path + the name → entry map.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct PartLibrary {
    /// Root directory the index file + parts live under.
    pub root: PathBuf,
    /// name → PartEntry.
    pub entries: HashMap<String, PartEntry>,
}

impl PartLibrary {
    /// Empty library rooted at `root`.
    pub fn new(root: impl Into<PathBuf>) -> Self {
        Self {
            root: root.into(),
            entries: HashMap::new(),
        }
    }

    /// Count of installed parts.
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// True when the library has no parts.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }
}

/// Read the on-disk index file at `root/index.ron`. Returns an empty
/// library (rooted at `root`) when the file doesn't exist — that's
/// the bootstrap case for a fresh library.
pub fn load_index(root: impl AsRef<Path>) -> Result<PartLibrary, PartLibError> {
    let root = root.as_ref();
    let path = root.join(INDEX_FILENAME);
    if !path.exists() {
        return Ok(PartLibrary::new(root.to_path_buf()));
    }
    // Round-23 sweep: bound the index read at MAX_PARTLIB_INDEX_BYTES
    // (16 MiB) — sister to MAX_DOC_FILE_BYTES. Even a large org's
    // parts library is in the low MiB; 16 MiB refuses pathological
    // inputs.
    let text = valenx_core::io_caps::read_capped_to_string(
        &path,
        valenx_core::io_caps::MAX_PARTLIB_INDEX_BYTES as usize,
    )
    .map_err(|e| PartLibError::Io {
        path: path.display().to_string(),
        reason: e.to_string(),
    })?;
    let file: LibraryFile =
        ron::de::from_str(&text).map_err(|e| PartLibError::Ron(e.to_string()))?;
    if file.version != VERSION {
        return Err(PartLibError::Ron(format!(
            "index version mismatch: file = {}, expected = {}",
            file.version, VERSION
        )));
    }
    let mut library = file.library;
    library.root = root.to_path_buf();
    Ok(library)
}

/// Save the index file at `root/index.ron`.
pub fn save_index(library: &PartLibrary) -> Result<(), PartLibError> {
    let path = library.root.join(INDEX_FILENAME);
    fs::create_dir_all(&library.root).map_err(|e| PartLibError::Io {
        path: library.root.display().to_string(),
        reason: e.to_string(),
    })?;
    let file = LibraryFile {
        version: VERSION,
        library: library.clone(),
    };
    let text = ron::ser::to_string_pretty(&file, ron::ser::PrettyConfig::default())
        .map_err(|e| PartLibError::Ron(e.to_string()))?;
    valenx_core::io_caps::atomic_write_str(&path, &text).map_err(|e| PartLibError::Io {
        path: path.display().to_string(),
        reason: e.to_string(),
    })?;
    Ok(())
}

/// Copy `file_path` into `library.root` and register the resulting
/// [`PartEntry`] under `name`. Computes a SHA-256 of the file at
/// install time so callers can detect on-disk tampering later.
///
/// Returns [`PartLibError::DuplicatePart`] when `name` is already
/// present — overwriting is intentional opt-in via `remove` + re-add.
pub fn install_local(
    library: &mut PartLibrary,
    name: &str,
    file_path: impl AsRef<Path>,
) -> Result<PartEntry, PartLibError> {
    if name.is_empty() {
        return Err(PartLibError::BadParameter {
            name: "name",
            reason: "must not be empty".into(),
        });
    }
    if library.entries.contains_key(name) {
        return Err(PartLibError::DuplicatePart {
            name: name.to_string(),
        });
    }
    let src = file_path.as_ref();
    if !src.exists() {
        return Err(PartLibError::Io {
            path: src.display().to_string(),
            reason: "source file does not exist".into(),
        });
    }

    fs::create_dir_all(&library.root).map_err(|e| PartLibError::Io {
        path: library.root.display().to_string(),
        reason: e.to_string(),
    })?;

    let ext = src
        .extension()
        .and_then(|s| s.to_str())
        .unwrap_or("")
        .to_string();
    let kind = PartKind::from_extension(&ext).ok_or_else(|| PartLibError::BadParameter {
        name: "file_path",
        reason: format!("unknown extension `{ext}`"),
    })?;

    let dst_name = format!("{name}.{ext}");
    let dst = library.root.join(&dst_name);
    fs::copy(src, &dst).map_err(|e| PartLibError::Io {
        path: dst.display().to_string(),
        reason: e.to_string(),
    })?;

    let checksum = sha256_of_file(&dst)?;
    let entry = PartEntry {
        name: name.to_string(),
        kind,
        source_url: None,
        local_path: dst,
        checksum,
    };
    library.entries.insert(name.to_string(), entry.clone());
    Ok(entry)
}

/// Stub for remote fetch — returns
/// [`PartLibError::FetchRequiresNetwork`] in v1. v2 will route through
/// the Phase 22 Add-on Manager network pipeline (gh-release).
pub fn fetch_remote(
    _library: &mut PartLibrary,
    url: &str,
    _name: &str,
) -> Result<PartEntry, PartLibError> {
    Err(PartLibError::FetchRequiresNetwork {
        url: url.to_string(),
    })
}

/// Compute the SHA-256 hex digest of a file.
fn sha256_of_file(path: &Path) -> Result<String, PartLibError> {
    let mut f = fs::File::open(path).map_err(|e| PartLibError::Io {
        path: path.display().to_string(),
        reason: e.to_string(),
    })?;
    let mut hasher = Sha256::new();
    let mut buf = [0u8; 4096];
    loop {
        let n = f.read(&mut buf).map_err(|e| PartLibError::Io {
            path: path.display().to_string(),
            reason: e.to_string(),
        })?;
        if n == 0 {
            break;
        }
        hasher.update(&buf[..n]);
    }
    let out = hasher.finalize();
    Ok(out.iter().map(|b| format!("{b:02x}")).collect())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::env;

    fn tmpdir(tag: &str) -> PathBuf {
        let dir = env::temp_dir().join(format!(
            "valenx_partlib_{}_{}_{}",
            tag,
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .map(|d| d.as_nanos())
                .unwrap_or(0)
        ));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir).unwrap();
        dir
    }

    #[test]
    fn load_missing_index_returns_empty() {
        let dir = tmpdir("load_missing");
        let lib = load_index(&dir).expect("ok");
        assert!(lib.is_empty());
        assert_eq!(lib.root, dir);
    }

    #[test]
    fn install_local_copies_and_indexes() {
        let dir = tmpdir("install");
        let src = dir.join("cube.stl");
        fs::write(&src, b"solid cube\nendsolid").unwrap();
        let mut lib = PartLibrary::new(dir.join("library"));
        let entry = install_local(&mut lib, "cube", &src).expect("ok");
        assert_eq!(entry.name, "cube");
        assert_eq!(entry.kind, PartKind::StlMesh);
        assert!(entry.local_path.exists());
        assert_eq!(entry.checksum.len(), 64);
        assert_eq!(lib.len(), 1);
    }

    #[test]
    fn install_local_rejects_duplicate() {
        let dir = tmpdir("dup");
        let src = dir.join("a.stl");
        fs::write(&src, b"x").unwrap();
        let mut lib = PartLibrary::new(dir.join("library"));
        install_local(&mut lib, "a", &src).unwrap();
        let err = install_local(&mut lib, "a", &src).unwrap_err();
        assert!(matches!(err, PartLibError::DuplicatePart { .. }));
    }

    #[test]
    fn fetch_remote_returns_capability_error() {
        let mut lib = PartLibrary::new("/tmp/no");
        let err = fetch_remote(&mut lib, "https://x/y.step", "y").unwrap_err();
        assert!(matches!(err, PartLibError::FetchRequiresNetwork { .. }));
    }

    #[test]
    fn save_and_reload_round_trips() {
        let dir = tmpdir("rt");
        let src = dir.join("plate.step");
        fs::write(&src, b"ISO-10303-21;").unwrap();
        let mut lib = PartLibrary::new(dir.join("lib"));
        install_local(&mut lib, "plate", &src).unwrap();
        save_index(&lib).unwrap();
        let reloaded = load_index(&lib.root).unwrap();
        assert_eq!(reloaded.len(), 1);
        assert!(reloaded.entries.contains_key("plate"));
    }

    /// Round-23 RED→GREEN: `load_index` rejects an over-cap
    /// `index.ron` at the read-cap layer (MAX_PARTLIB_INDEX_BYTES,
    /// 16 MiB) rather than slurping unbounded.
    #[test]
    fn load_index_rejects_oversize() {
        let dir = tmpdir("oversize");
        let index = dir.join(INDEX_FILENAME);
        let oversize = (valenx_core::io_caps::MAX_PARTLIB_INDEX_BYTES as usize) + 1024;
        fs::write(&index, vec![b'x'; oversize]).unwrap();
        let err = load_index(&dir).expect_err("must reject oversize");
        match err {
            PartLibError::Io { reason, .. } => {
                assert!(
                    reason.contains("exceeds") || reason.contains("cap"),
                    "expected cap-exceeded reason, got: {reason}"
                );
            }
            other => panic!("expected Io, got: {other:?}"),
        }
    }
}
