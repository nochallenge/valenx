//! Loader for `.valenx` projects on disk.
//!
//! Per RFC 0001 a project is a **directory** named `foo.valenx/`
//! containing `project.toml`, optional `tools.lock`, a `cases/`
//! subdirectory, and assets. The loader:
//!
//! 1. Checks `project.toml` exists and parses as TOML
//! 2. Validates format-version compatibility (SemVer)
//! 3. Loads `tools.lock` if present (optional per RFC 0001)
//! 4. Loads every case referenced in `project.toml [cases].order`
//! 5. Rejects absolute paths and paths escaping the project root
//! 6. Returns a `LoadedProject` with the root pinned

use std::collections::BTreeMap;
use std::io::Read;
use std::path::{Component, Path, PathBuf};

use thiserror::Error;

/// Per-file size cap for `project.toml`, `tools.lock`, and each
/// `case.toml`. Round-11 hardening: pre-fix every site went through a
/// bare `std::fs::read_to_string` which would slurp an arbitrarily
/// large file before the TOML parser ever saw it. 1 MiB is well above
/// any realistic project file (the schema is a few dozen fields per
/// case at most); the cap kills a hostile or accidental multi-GB
/// payload before allocation. Mirrors the cap in
/// `valenx-app::rbac_io::rbac_override_from_project_toml`.
pub const MAX_PROJECT_FILE_BYTES: u64 = 1024 * 1024;

use super::case_def::CaseDef;
use super::manifest::Project;
use super::tools_lock::ToolsLock;

/// Supported project format major version. RFC 0001 promises
/// forward-compat within the same major; newer majors require
/// migration (not yet implemented).
pub const SUPPORTED_MAJOR: u64 = 1;

/// Errors raised by [`LoadedProject::load`] and helpers.
#[derive(Debug, Error)]
pub enum ProjectLoadError {
    /// The supplied path is not a directory at all.
    #[error("project root {path} does not exist or is not a directory")]
    NotADirectory {
        /// The offending path.
        path: PathBuf,
    },

    /// The directory exists but contains no `project.toml`.
    #[error("missing project.toml at {path}")]
    MissingManifest {
        /// Path to the expected manifest.
        path: PathBuf,
    },

    /// Underlying filesystem IO error while reading a project file.
    #[error("failed to read {path}: {source}")]
    Io {
        /// Path that failed to read.
        path: PathBuf,
        /// The underlying [`std::io::Error`].
        #[source]
        source: std::io::Error,
    },

    /// A `.toml` file is syntactically invalid.
    #[error("{path} is not valid TOML: {source}")]
    Parse {
        /// Path to the offending file.
        path: PathBuf,
        // Boxed because `toml::de::Error` is ~280 bytes, which makes
        // every `Result<_, ProjectLoadError>` expensive to pass
        // around. See clippy's `result_large_err`.
        /// The underlying [`toml::de::Error`] (boxed for size).
        #[source]
        source: Box<toml::de::Error>,
    },

    /// The project format major version is newer than what this build
    /// of Valenx supports.
    #[error("project format {found} is newer than supported major v{supported}; upgrade Valenx")]
    UnsupportedFormatMajor {
        /// The `format` field as seen in the manifest.
        found: String,
        /// The highest major this binary understands.
        supported: u64,
    },

    /// The `format` field couldn't be parsed as a SemVer string.
    #[error("project.toml format field {0:?} is not a valid SemVer")]
    InvalidFormatVersion(String),

    /// `[cases].order` referenced a case directory that does not
    /// exist on disk.
    #[error("case {name} listed in [cases].order but {path} is missing")]
    MissingCase {
        /// Case name from `[cases].order`.
        name: String,
        /// Path the loader expected to find it at.
        path: PathBuf,
    },

    /// A relative path inside the project escapes the root (uses
    /// `..` traversal or an absolute prefix).
    #[error("path {0:?} is absolute or escapes the project root")]
    UnsafePath(PathBuf),

    /// A case name in `[cases].order` is not a single path component —
    /// it contains `..`, a path separator (`/` or `\`), or an absolute
    /// prefix. Such a name would let `cases/<name>` escape the
    /// `cases/` subdir (e.g. `abc/../def` or `../escape`), and at
    /// `save()` time the write could land outside the project root if
    /// a traversed component is a symlink.
    #[error("case name {0:?} is not a single path component (contains .., separators, or is absolute)")]
    UnsafeCaseName(String),

    /// A project file (`project.toml`, `tools.lock`, or a per-case
    /// `case.toml`) exceeds [`MAX_PROJECT_FILE_BYTES`]. Round-11
    /// hardening: pre-fix a multi-GB file would force the loader to
    /// allocate the whole payload before the TOML parser rejected it.
    #[error("project file {path} is {size} bytes (cap {cap}); refusing to load")]
    FileTooLarge {
        /// Path to the oversized file.
        path: PathBuf,
        /// Observed file size in bytes.
        size: u64,
        /// The hard cap in bytes.
        cap: u64,
    },
}

/// A fully-loaded project, pinned to its on-disk location.
#[derive(Debug)]
pub struct LoadedProject {
    /// Absolute path to the directory containing `project.toml`.
    pub root: PathBuf,
    pub project: Project,
    pub tools_lock: Option<ToolsLock>,
    pub cases: BTreeMap<String, CaseDef>,
}

impl LoadedProject {
    /// Open a project at `root` (the `.valenx` directory).
    pub fn load(root: impl Into<PathBuf>) -> Result<Self, ProjectLoadError> {
        let root = root.into();
        let root = root.canonicalize().map_err(|e| ProjectLoadError::Io {
            path: root.clone(),
            source: e,
        })?;

        if !root.is_dir() {
            return Err(ProjectLoadError::NotADirectory { path: root });
        }

        // 1. project.toml
        let manifest_path = root.join("project.toml");
        if !manifest_path.is_file() {
            return Err(ProjectLoadError::MissingManifest {
                path: manifest_path,
            });
        }
        let raw = read_capped(&manifest_path)?;
        let project: Project = toml::from_str(&raw).map_err(|e| ProjectLoadError::Parse {
            path: manifest_path.clone(),
            source: Box::new(e),
        })?;

        // 2. Format-version gate
        let format_major = parse_format_major(&project.project.format)?;
        if format_major > SUPPORTED_MAJOR {
            return Err(ProjectLoadError::UnsupportedFormatMajor {
                found: project.project.format.clone(),
                supported: SUPPORTED_MAJOR,
            });
        }

        // 3. Path safety — every geometry source must stay inside the root.
        // Round-24 L4: lexical check first (cheap rejection of `..` and
        // absolute paths), then canonical check (refuses symlinks
        // pointing outside the project root). The canonical check is
        // best-effort: if the source doesn't exist yet (geometry stub),
        // we skip it; the lexical check still protected against the
        // path traversal case.
        for entry in &project.geometry.entries {
            ensure_within(&entry.source)?;
            ensure_canonical_within(&root, &entry.source)?;
        }
        for mesh in project.mesh.values() {
            ensure_within(&mesh.source)?;
            ensure_canonical_within(&root, &mesh.source)?;
        }

        // 4. tools.lock (optional)
        let lock_path = root.join("tools.lock");
        let tools_lock = if lock_path.is_file() {
            let raw = read_capped(&lock_path)?;
            let lock: ToolsLock = toml::from_str(&raw).map_err(|e| ProjectLoadError::Parse {
                path: lock_path.clone(),
                source: Box::new(e),
            })?;
            Some(lock)
        } else {
            None
        };

        // 5. Cases — every name in `order` must have a readable case.toml.
        let mut cases: BTreeMap<String, CaseDef> = BTreeMap::new();
        for name in &project.cases.order {
            // M1: a case name must be a single path component. A
            // hostile name like `../escape` or `abc/../def` would let
            // `cases/<name>` resolve outside the `cases/` subdir; the
            // old `is_file()` check below would happily follow it, and
            // `save()` would then write through it. Reject before the
            // join so neither read nor (later) write can escape.
            if !is_safe_case_name(name) {
                return Err(ProjectLoadError::UnsafeCaseName(name.clone()));
            }
            let case_path = root.join("cases").join(name).join("case.toml");
            if !case_path.is_file() {
                return Err(ProjectLoadError::MissingCase {
                    name: name.clone(),
                    path: case_path,
                });
            }
            let raw = read_capped(&case_path)?;
            let case: CaseDef = toml::from_str(&raw).map_err(|e| ProjectLoadError::Parse {
                path: case_path.clone(),
                source: Box::new(e),
            })?;
            cases.insert(name.clone(), case);
        }

        Ok(LoadedProject {
            root,
            project,
            tools_lock,
            cases,
        })
    }

    /// Resolve a project-relative path into an absolute disk path.
    /// Every caller that touches the filesystem should go through
    /// this — `ensure_within` has already been enforced on load.
    pub fn absolute_path(&self, relative: &Path) -> PathBuf {
        self.root.join(relative)
    }

    /// Ordered list of case names per `project.toml`.
    pub fn case_names(&self) -> &[String] {
        &self.project.cases.order
    }
}

/// Read a project file (TOML) with the [`MAX_PROJECT_FILE_BYTES`]
/// cap enforced both by stat AND by a bounded `take(...)` on the
/// reader (TOCTOU defence-in-depth — the stat might lie if the file
/// grew between the metadata call and the open).
fn read_capped(path: &Path) -> Result<String, ProjectLoadError> {
    let md = std::fs::metadata(path).map_err(|e| ProjectLoadError::Io {
        path: path.to_path_buf(),
        source: e,
    })?;
    if md.len() > MAX_PROJECT_FILE_BYTES {
        return Err(ProjectLoadError::FileTooLarge {
            path: path.to_path_buf(),
            size: md.len(),
            cap: MAX_PROJECT_FILE_BYTES,
        });
    }
    let mut buf = Vec::new();
    std::fs::File::open(path)
        .map_err(|e| ProjectLoadError::Io {
            path: path.to_path_buf(),
            source: e,
        })?
        .take(MAX_PROJECT_FILE_BYTES)
        .read_to_end(&mut buf)
        .map_err(|e| ProjectLoadError::Io {
            path: path.to_path_buf(),
            source: e,
        })?;
    String::from_utf8(buf).map_err(|e| ProjectLoadError::Io {
        path: path.to_path_buf(),
        source: std::io::Error::new(std::io::ErrorKind::InvalidData, e),
    })
}

fn parse_format_major(s: &str) -> Result<u64, ProjectLoadError> {
    // Accept "1", "1.0", "1.2", "1.2.3" — first component is major.
    let first = s.split('.').next().unwrap_or("");
    first
        .parse::<u64>()
        .map_err(|_| ProjectLoadError::InvalidFormatVersion(s.to_string()))
}

/// M1: a case name must resolve to exactly one `Component::Normal` —
/// no `..`, no `.`, no `/` or `\` separators, no absolute / drive
/// prefix. This is stricter than [`ensure_within`] (which permits
/// multi-segment relative paths like `mesh/a/b`) because a case name
/// is joined as a single directory under `cases/`; anything with a
/// separator or traversal could escape that subdir. Shared with
/// `writer.rs::save` as a belt-and-suspenders defensive check.
pub(crate) fn is_safe_case_name(name: &str) -> bool {
    // Fast lexical rejections first: empty, separators, traversal.
    if name.is_empty()
        || name.contains('/')
        || name.contains('\\')
        || name.contains("..")
    {
        return false;
    }
    // R30 L (cosmetic, Windows-only): Windows silently strips trailing
    // dots and spaces from filenames, so `foo.` / `foo ` alias `foo` on
    // disk and two declared case names would collide into one directory.
    // Reject the trailing-dot/space forms on Windows so the on-disk name
    // always matches the declared name. (This is de-aliasing, not a
    // traversal guard — `..` and separators are already handled above.)
    // On Unix a trailing dot/space is a legitimate, distinct filename, so
    // it is left untouched.
    #[cfg(windows)]
    if name.ends_with('.') || name.ends_with(' ') {
        return false;
    }
    let path = Path::new(name);
    if path.is_absolute() {
        return false;
    }
    // Exactly one Normal component, nothing else (no CurDir, ParentDir,
    // RootDir, or Prefix).
    let mut comps = path.components();
    matches!(
        (comps.next(), comps.next()),
        (Some(Component::Normal(_)), None)
    )
}

/// Reject absolute paths and any path whose resolution leaves the
/// project root via `..`.
fn ensure_within(rel: &Path) -> Result<(), ProjectLoadError> {
    if rel.is_absolute() {
        return Err(ProjectLoadError::UnsafePath(rel.to_path_buf()));
    }
    let mut depth: i32 = 0;
    for c in rel.components() {
        match c {
            Component::ParentDir => depth -= 1,
            Component::Normal(_) => depth += 1,
            Component::CurDir => {}
            Component::Prefix(_) | Component::RootDir => {
                return Err(ProjectLoadError::UnsafePath(rel.to_path_buf()));
            }
        }
        if depth < 0 {
            return Err(ProjectLoadError::UnsafePath(rel.to_path_buf()));
        }
    }
    Ok(())
}

/// Round-24 L4: defence-in-depth check on top of [`ensure_within`].
/// Resolves the relative source path against the project root and
/// canonicalises both sides; refuses any source that — after
/// symlink resolution — lives outside the canonical project root.
/// The lexical `ensure_within` blocks `../../etc/passwd`, but a
/// `geometry/case.step` symlinked to `/etc/passwd` would pass the
/// lexical check while still escaping. This check catches that.
///
/// Best-effort: a source path that doesn't exist on disk (geometry
/// stub, mesh placeholder) skips the canonical check — the lexical
/// `ensure_within` already protected against the traversal case.
fn ensure_canonical_within(root: &Path, rel: &Path) -> Result<(), ProjectLoadError> {
    let joined = root.join(rel);
    // If the source doesn't exist (yet), the lexical check was the
    // only guard available. canonicalize would error with
    // NotFound — skip to preserve the load-stub behaviour.
    let source_canon = match joined.canonicalize() {
        Ok(c) => c,
        Err(_) => return Ok(()),
    };
    let root_canon = match root.canonicalize() {
        Ok(c) => c,
        // If the root itself can't be canonicalised the whole load
        // is suspect anyway; fall through to the lexical safety net.
        Err(_) => return Ok(()),
    };
    if !source_canon.starts_with(&root_canon) {
        return Err(ProjectLoadError::UnsafePath(rel.to_path_buf()));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn format_major_parsing() {
        assert_eq!(parse_format_major("1").unwrap(), 1);
        assert_eq!(parse_format_major("1.0").unwrap(), 1);
        assert_eq!(parse_format_major("2.3.1").unwrap(), 2);
        assert!(parse_format_major("abc").is_err());
    }

    #[test]
    fn within_root_allows_normal_paths() {
        assert!(ensure_within(Path::new("geometry/a.step")).is_ok());
        assert!(ensure_within(Path::new("mesh/abc/mesh.msh")).is_ok());
        assert!(ensure_within(Path::new("a/b/./c")).is_ok());
    }

    #[test]
    fn within_root_rejects_absolute() {
        #[cfg(unix)]
        {
            assert!(ensure_within(Path::new("/etc/passwd")).is_err());
        }
        #[cfg(windows)]
        {
            assert!(ensure_within(Path::new("C:\\Windows\\system.ini")).is_err());
        }
    }

    #[test]
    fn within_root_rejects_escape() {
        assert!(ensure_within(Path::new("../secret.txt")).is_err());
        assert!(ensure_within(Path::new("a/../../etc")).is_err());
    }

    #[test]
    fn safe_case_name_accepts_single_component() {
        assert!(is_safe_case_name("smoke"));
        assert!(is_safe_case_name("case-01"));
        assert!(is_safe_case_name("Case_With.Dots"));
    }

    #[test]
    fn safe_case_name_rejects_traversal_and_separators() {
        assert!(!is_safe_case_name(""));
        assert!(!is_safe_case_name(".."));
        assert!(!is_safe_case_name("../escape"));
        assert!(!is_safe_case_name("abc/../def"));
        assert!(!is_safe_case_name("a/b"));
        assert!(!is_safe_case_name("a\\b"));
        assert!(!is_safe_case_name("."));
        #[cfg(unix)]
        assert!(!is_safe_case_name("/etc"));
        #[cfg(windows)]
        {
            assert!(!is_safe_case_name("C:\\Windows"));
            assert!(!is_safe_case_name("\\\\server\\share"));
        }
    }

    /// R30 L (cosmetic): Windows silently strips trailing dots/spaces
    /// from filenames, so `foo.` and `foo ` alias `foo` on disk — two
    /// distinct case names collide to one directory. Reject them on
    /// Windows so the on-disk name always matches the declared name.
    /// (Not a traversal escape — purely de-aliasing.) On Unix a trailing
    /// dot/space is a legitimate distinct filename, so it stays allowed.
    #[test]
    #[cfg(windows)]
    fn safe_case_name_rejects_windows_trailing_dot_or_space() {
        assert!(!is_safe_case_name("foo."));
        assert!(!is_safe_case_name("foo "));
        assert!(!is_safe_case_name("foo.."));
        assert!(!is_safe_case_name("foo  "));
        assert!(!is_safe_case_name("foo. "));
        // Internal dots/spaces remain fine (the existing accept test
        // covers `Case_With.Dots`); only the trailing position aliases.
        assert!(is_safe_case_name("a.b"));
        assert!(is_safe_case_name("a b"));
    }

    /// The trailing-dot/space rejection is Windows-only: on Unix these
    /// are valid, distinct filenames and must stay accepted.
    #[test]
    #[cfg(unix)]
    fn safe_case_name_allows_unix_trailing_dot_or_space() {
        assert!(is_safe_case_name("foo."));
        assert!(is_safe_case_name("foo "));
    }

    /// R29 M1 RED→GREEN — a `project.toml` whose `[cases].order` names
    /// a case that traverses out of `cases/` (`../escape`) must be
    /// rejected at load with `UnsafeCaseName`. Pre-fix the loader only
    /// checked `case_path.is_file()`, so a name like `../escape`
    /// survived ingestion and `save()` could later write through it.
    #[test]
    fn load_rejects_traversing_case_name() {
        let mut tmp = std::env::temp_dir();
        let nanos = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map(|d| d.as_nanos())
            .unwrap_or(0);
        tmp.push(format!("valenx-r29-m1-loader-{nanos}"));
        std::fs::create_dir_all(&tmp).expect("create tmp dir");

        // Minimal manifest with a traversing case name in order.
        let manifest = "\
[project]
name = \"m1\"
format = \"1.0\"

[cases]
order = [\"../escape\"]
";
        std::fs::write(tmp.join("project.toml"), manifest).expect("write manifest");

        // Plant a real case.toml at the escaped location so that the
        // OLD `is_file()` check would have passed (proving the lexical
        // guard, not a missing file, is what rejects it).
        let escaped_dir = tmp.join("escape");
        std::fs::create_dir_all(&escaped_dir).expect("create escaped dir");
        std::fs::write(
            escaped_dir.join("case.toml"),
            "[case]\nformat = \"1.0\"\nname = \"escape\"\nphysics = \"cfd\"\nsolver = \"x.y\"\nmesh = \"(none)\"\n",
        )
        .expect("write escaped case.toml");

        let err = LoadedProject::load(&tmp)
            .expect_err("traversing case name must be rejected at load");
        match err {
            ProjectLoadError::UnsafeCaseName(n) => {
                assert_eq!(n, "../escape", "error must name the offending case");
            }
            other => panic!("expected UnsafeCaseName, got {other:?}"),
        }

        let _ = std::fs::remove_dir_all(&tmp);
    }

    /// Round-11 RED→GREEN — R11-2. Pre-fix `LoadedProject::load`
    /// went through unbounded `std::fs::read_to_string` at three
    /// sites (`project.toml`, `tools.lock`, each `case.toml`). A
    /// hostile or malformed multi-MiB `project.toml` would force
    /// the loader to allocate the whole payload before the TOML
    /// parser rejected it. The cap (`MAX_PROJECT_FILE_BYTES = 1
    /// MiB`) hard-rejects with `FileTooLarge` before allocation.
    #[test]
    fn load_rejects_project_toml_larger_than_cap() {
        // Build a unique temp directory by hand — no `tempfile`
        // dev-dep in valenx-core.
        let mut tmp = std::env::temp_dir();
        let nanos = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map(|d| d.as_nanos())
            .unwrap_or(0);
        tmp.push(format!("valenx-r11-loader-{nanos}"));
        std::fs::create_dir_all(&tmp).expect("create tmp dir");

        let manifest_path = tmp.join("project.toml");
        // Pad past the 1 MiB cap with a syntactically-valid TOML
        // comment so the size guard fires first (the cap must
        // reject before parsing).
        let mut padded = String::with_capacity((MAX_PROJECT_FILE_BYTES as usize) + 4096);
        padded.push_str(
            "[project]\nname = \"oversized\"\nformat = \"1.0\"\n# pad below\n# ",
        );
        // Roughly 5 MiB of '#' comment bytes so we comfortably exceed
        // the 1 MiB cap regardless of OS line-ending rewrites.
        padded.extend(std::iter::repeat_n('A', 5 * 1024 * 1024));
        padded.push('\n');
        std::fs::write(&manifest_path, &padded).expect("write oversized manifest");

        let err = LoadedProject::load(&tmp)
            .expect_err("oversized project.toml must be rejected before alloc");
        match err {
            ProjectLoadError::FileTooLarge { size, cap, .. } => {
                assert!(
                    size > cap,
                    "FileTooLarge size ({size}) must exceed cap ({cap})"
                );
                assert_eq!(cap, MAX_PROJECT_FILE_BYTES);
            }
            other => panic!("expected ProjectLoadError::FileTooLarge, got {other:?}"),
        }

        // Best-effort cleanup; the OS will sweep eventually.
        let _ = std::fs::remove_dir_all(&tmp);
    }
}
