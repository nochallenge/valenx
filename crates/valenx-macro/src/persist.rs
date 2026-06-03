//! Filesystem helpers for the macro library.
//!
//! Macros live in `~/.valenx/macros/{name}.ron`. This module exposes
//! load / save / list helpers that the Macro Library UI calls.

use std::fs;
use std::io;
use std::path::{Path, PathBuf};

use thiserror::Error;

use crate::Macro;

/// Errors writing / loading macros.
#[derive(Debug, Error)]
pub enum MacroPersistError {
    /// I/O failure (file system, etc.).
    #[error("io: {0}")]
    Io(#[from] io::Error),
    /// RON encoding failure.
    #[error("ron encode: {0}")]
    Encode(#[from] ron::Error),
    /// RON parse failure.
    #[error("ron parse: {0}")]
    Decode(#[from] ron::error::SpannedError),
}

/// Default macro directory under `~/.valenx/macros/`. Returns `None`
/// if the home directory can't be resolved (rare; CI environments).
pub fn default_macro_dir() -> Option<PathBuf> {
    let home = std::env::var_os("USERPROFILE").or_else(|| std::env::var_os("HOME"))?;
    Some(PathBuf::from(home).join(".valenx").join("macros"))
}

/// Save a macro to `dir/{macro.name}.ron`. Creates `dir` if missing.
///
/// # Errors
///
/// - [`MacroPersistError::Io`] on filesystem failure (permission
///   denied, no space, etc.).
/// - [`MacroPersistError::Encode`] on serialization failure.
pub fn save(m: &Macro, dir: &Path) -> Result<PathBuf, MacroPersistError> {
    fs::create_dir_all(dir)?;
    let path = dir.join(format!("{}.ron", sanitise_filename(&m.name)));
    let text = m.to_ron()?;
    // Round-28 H2: routed through the canonical
    // `valenx_core::io_caps::atomic_write_str` (sidecar O_NOFOLLOW,
    // fsync-before-rename, parent-dir fsync on Unix). Pre-fix this
    // was a bare `fs::write` which silently followed leaf symlinks
    // and was non-atomic.
    valenx_core::io_caps::atomic_write_str(&path, &text)?;
    Ok(path)
}

/// Load a macro from `path`.
///
/// # Errors
///
/// - [`MacroPersistError::Io`] on read failure.
/// - [`MacroPersistError::Decode`] on RON parse failure.
pub fn load(path: &Path) -> Result<Macro, MacroPersistError> {
    // R29 D: bounded at valenx_core::io_caps::MAX_DOC_FILE_BYTES (16 MiB)
    // via the canonical helper, replacing the per-crate private dupe.
    let text = valenx_core::io_caps::read_capped_to_string(
        path,
        valenx_core::io_caps::MAX_DOC_FILE_BYTES,
    )?;
    let m = Macro::from_ron(&text)?;
    Ok(m)
}

/// List `.ron` files in `dir`, returning their full paths sorted
/// lexicographically. Returns an empty vec if `dir` doesn't exist.
pub fn list_files(dir: &Path) -> Result<Vec<PathBuf>, io::Error> {
    if !dir.exists() {
        return Ok(Vec::new());
    }
    let mut paths: Vec<PathBuf> = fs::read_dir(dir)?
        .filter_map(|e| e.ok())
        .map(|e| e.path())
        .filter(|p| p.extension().and_then(|s| s.to_str()) == Some("ron"))
        .collect();
    paths.sort();
    Ok(paths)
}

/// Delete a macro file. Returns `Ok(false)` if the file didn't exist.
///
/// # Errors
///
/// - [`MacroPersistError::Io`] for any deletion error other than
///   `NotFound`.
pub fn delete(path: &Path) -> Result<bool, MacroPersistError> {
    match fs::remove_file(path) {
        Ok(()) => Ok(true),
        Err(e) if e.kind() == io::ErrorKind::NotFound => Ok(false),
        Err(e) => Err(MacroPersistError::Io(e)),
    }
}

/// Replace any character that's bad in a filename with `_`.
fn sanitise_filename(name: &str) -> String {
    name.chars()
        .map(|c| {
            if c.is_alphanumeric() || c == '-' || c == '_' || c == '.' {
                c
            } else {
                '_'
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::action::{MacroAction, PanelId};

    fn tmpdir(name: &str) -> PathBuf {
        let p = std::env::temp_dir().join(format!("valenx_macro_test_{name}"));
        let _ = fs::remove_dir_all(&p);
        p
    }

    #[test]
    fn sanitise_replaces_bad_chars() {
        assert_eq!(sanitise_filename("ok-name"), "ok-name");
        assert_eq!(sanitise_filename("bad/name"), "bad_name");
        assert_eq!(sanitise_filename("space here"), "space_here");
    }

    #[test]
    fn save_then_load_round_trip() {
        let dir = tmpdir("save_load");
        let mut m = Macro::new("Test");
        m.push(MacroAction::SwitchPanel {
            panel_id: PanelId::Sketcher,
        });
        let path = save(&m, &dir).unwrap();
        let m2 = load(&path).unwrap();
        let _ = fs::remove_dir_all(&dir);
        assert_eq!(m.name, m2.name);
        assert_eq!(m.actions.len(), m2.actions.len());
    }

    #[test]
    fn list_files_returns_sorted_ron_files() {
        let dir = tmpdir("list_files");
        fs::create_dir_all(&dir).unwrap();
        fs::write(dir.join("b.ron"), "()").unwrap();
        fs::write(dir.join("a.ron"), "()").unwrap();
        fs::write(dir.join("ignored.txt"), "x").unwrap();
        let files = list_files(&dir).unwrap();
        let _ = fs::remove_dir_all(&dir);
        assert_eq!(files.len(), 2);
        assert!(files[0].ends_with("a.ron"));
        assert!(files[1].ends_with("b.ron"));
    }

    #[test]
    fn delete_missing_file_returns_false() {
        let path = std::env::temp_dir().join("valenx_macro_doesnt_exist.ron");
        let _ = fs::remove_file(&path);
        assert!(!delete(&path).unwrap());
    }

    #[test]
    fn list_files_on_missing_dir_is_empty() {
        let dir = std::env::temp_dir().join("valenx_macro_no_such_dir");
        let _ = fs::remove_dir_all(&dir);
        let files = list_files(&dir).unwrap();
        assert!(files.is_empty());
    }
}
