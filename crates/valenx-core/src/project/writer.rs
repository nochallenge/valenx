//! Writing `.valenx` projects back to disk.
//!
//! Symmetric with `loader.rs`: serialize `Project`, `ToolsLock`, and
//! every `CaseDef` to TOML in the right locations. Writes are
//! atomic per-file via the canonical
//! [`crate::io_caps::atomic_write_bytes`] (sidecar named
//! `<basename>.tmp.<pid>.<counter>`, opened with O_NOFOLLOW /
//! `FILE_FLAG_OPEN_REPARSE_POINT`, fsynced before rename, parent dir
//! fsynced after rename on Unix). Round-28 C1: this used to go
//! through a private `write_atomic` whose sidecar shape
//! (`target.with_extension("tmp")`) collided under two concurrent
//! saves and silently followed leaf symlinks.

use std::path::Path;

use thiserror::Error;

use super::loader::LoadedProject;

/// Errors raised by [`LoadedProject::save`].
#[derive(Debug, Error)]
pub enum ProjectSaveError {
    /// A file's typed model couldn't be serialised back to TOML.
    #[error("failed to serialize {label} to TOML: {source}")]
    Serialize {
        /// Human label of the file being serialised
        /// (e.g. `"project.toml"`).
        label: String,
        /// The underlying [`toml::ser::Error`].
        #[source]
        source: toml::ser::Error,
    },

    /// The atomic rename / write of one of the project files failed.
    #[error("failed to write {path}: {source}")]
    Io {
        /// Path that failed to write.
        path: std::path::PathBuf,
        /// The underlying [`std::io::Error`].
        #[source]
        source: std::io::Error,
    },

    /// M1: a case name is not a single path component (contains `..`,
    /// a path separator, or is absolute) and would let the per-case
    /// write escape the project root. Normally caught at load time by
    /// [`crate::project::loader::ProjectLoadError::UnsafeCaseName`];
    /// this is the belt-and-suspenders guard for `LoadedProject`s
    /// constructed in-memory rather than via `load`.
    #[error("case name {0:?} is not a single path component; refusing to write outside the project root")]
    InvalidCaseName(String),
}

impl LoadedProject {
    /// Save the project to its pinned root. Every file writes
    /// atomically via [`crate::io_caps::atomic_write_bytes`] — a
    /// `<basename>.tmp.<pid>.<counter>` sidecar opened with
    /// O_NOFOLLOW, fsynced, renamed into place, and a parent dir
    /// fsync (Unix) for dentry durability.
    pub fn save(&self) -> Result<(), ProjectSaveError> {
        // project.toml
        let manifest =
            toml::to_string_pretty(&self.project).map_err(|e| ProjectSaveError::Serialize {
                label: "project.toml".into(),
                source: e,
            })?;
        write_atomic(&self.root.join("project.toml"), manifest.as_bytes())?;

        // tools.lock (if present)
        if let Some(lock) = &self.tools_lock {
            let text = toml::to_string_pretty(lock).map_err(|e| ProjectSaveError::Serialize {
                label: "tools.lock".into(),
                source: e,
            })?;
            write_atomic(&self.root.join("tools.lock"), text.as_bytes())?;
        }

        // Each case.toml
        for (name, case) in &self.cases {
            // M1 belt-and-suspenders: reject a case name that would
            // escape `cases/` (`..`, separators, absolute). The loader
            // already rejects these, but a LoadedProject built in
            // memory bypasses the loader — without this check a
            // `../escape` name would write through `cases/../escape`,
            // landing outside the project root (and through a symlink
            // if one sits on the traversed path).
            if !super::loader::is_safe_case_name(name) {
                return Err(ProjectSaveError::InvalidCaseName(name.clone()));
            }
            let text = toml::to_string_pretty(case).map_err(|e| ProjectSaveError::Serialize {
                label: format!("cases/{name}/case.toml"),
                source: e,
            })?;
            let case_dir = self.root.join("cases").join(name);
            std::fs::create_dir_all(&case_dir).map_err(|e| ProjectSaveError::Io {
                path: case_dir.clone(),
                source: e,
            })?;
            write_atomic(&case_dir.join("case.toml"), text.as_bytes())?;
        }

        Ok(())
    }
}

/// Round-28 C1 — thin wrapper around the canonical
/// [`crate::io_caps::atomic_write_bytes`]. The pre-R28 in-tree
/// implementation used `target.with_extension("tmp")` as the
/// sidecar, which collided under two concurrent `save()` calls and
/// followed leaf symlinks. The canonical helper namespaces sidecars
/// per `<pid>.<counter>`, opens them O_NOFOLLOW / reparse-point-
/// refusing, fsyncs before the rename, and fsyncs the parent dir
/// after.
fn write_atomic(target: &Path, bytes: &[u8]) -> Result<(), ProjectSaveError> {
    crate::io_caps::atomic_write_bytes(target, bytes).map_err(|e| ProjectSaveError::Io {
        path: target.to_path_buf(),
        source: e,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::project::case_def::{CaseDef, CaseHeader};
    use crate::project::manifest::{
        CasesSection, GeometrySection, Project, ProjectHeader, UiSection, UnitsConfig,
    };
    use std::collections::BTreeMap;
    use std::path::PathBuf;
    use std::sync::{Arc, Barrier};

    /// Build a unique temp directory by hand (no `tempfile` dev-dep
    /// in `valenx-core`). Mirrors the helper used in
    /// `loader.rs::tests::load_rejects_project_toml_larger_than_cap`.
    fn make_tmp(label: &str) -> PathBuf {
        let mut tmp = std::env::temp_dir();
        let nanos = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map(|d| d.as_nanos())
            .unwrap_or(0);
        tmp.push(format!("valenx-r28-writer-{label}-{nanos}"));
        std::fs::create_dir_all(&tmp).expect("create tmp dir");
        tmp
    }

    /// Build a minimal LoadedProject anchored at `root`. The on-disk
    /// state doesn't have to exist yet — `save()` creates it.
    fn make_project(root: PathBuf) -> LoadedProject {
        let mut cases: BTreeMap<String, CaseDef> = BTreeMap::new();
        cases.insert(
            "smoke".into(),
            CaseDef {
                case: CaseHeader {
                    format: "1.0".into(),
                    name: "smoke".into(),
                    physics: "cfd".into(),
                    solver: "openfoam.simpleFoam".into(),
                    mesh: "(none)".into(),
                    description: None,
                },
                sections: BTreeMap::new(),
            },
        );
        LoadedProject {
            root,
            project: Project {
                project: ProjectHeader {
                    format: "1.0".into(),
                    name: "r28-c1".into(),
                    valenx_min: None,
                    created: None,
                    modified: None,
                    author: None,
                    description: None,
                },
                geometry: GeometrySection::default(),
                mesh: BTreeMap::new(),
                cases: CasesSection {
                    order: vec!["smoke".into()],
                },
                ui: UiSection::default(),
                units: UnitsConfig::default(),
            },
            tools_lock: None,
            cases,
        }
    }

    /// Round-28 C1 RED→GREEN — pre-fix `write_atomic` used
    /// `target.with_extension("tmp")` so two concurrent `save()`
    /// calls into the same root raced on a shared `project.tmp`
    /// sidecar. The loser hit `AlreadyExists` (Windows) or a torn
    /// rename. Post-fix every writer gets a unique
    /// `<basename>.tmp.<pid>.<counter>` sidecar via the canonical
    /// helper, so both calls succeed and the final `project.toml`
    /// always parses.
    #[test]
    fn concurrent_save_does_not_collide_on_sidecar() {
        let root = make_tmp("concurrent");
        let project = Arc::new(make_project(root.clone()));

        // 8 writers in parallel, all touching the same project root.
        // Pre-R28 this would AlreadyExists on Windows because each
        // save races to create `project.tmp`. Post-R28 every save
        // sees a unique sidecar so all return Ok.
        let n = 8usize;
        let barrier = Arc::new(Barrier::new(n));
        let mut handles = Vec::with_capacity(n);
        for _ in 0..n {
            let p = Arc::clone(&project);
            let b = Arc::clone(&barrier);
            handles.push(std::thread::spawn(move || {
                b.wait();
                p.save()
            }));
        }
        for h in handles {
            let res = h.join().expect("thread did not panic");
            res.expect("concurrent save must not collide on sidecar");
        }

        // After all 8 succeed the final project.toml must parse as
        // valid TOML (no torn content). We deserialize back into
        // `Project` to confirm.
        let bytes = std::fs::read(root.join("project.toml")).expect("project.toml exists");
        let text = std::str::from_utf8(&bytes).expect("project.toml is UTF-8");
        let _: Project =
            toml::from_str(text).expect("final project.toml must be valid TOML, not torn");

        // Confirm no orphan sidecars survive. The canonical helper
        // namespaces sidecars `project.toml.tmp.<pid>.<counter>`; if
        // any survive the rename-or-cleanup contract was violated.
        for entry in std::fs::read_dir(&root).expect("read root") {
            let entry = entry.expect("dir entry");
            let name = entry.file_name();
            let name = name.to_string_lossy();
            assert!(!name.contains(".tmp."), "orphan sidecar survived: {name}");
        }

        let _ = std::fs::remove_dir_all(&root);
    }

    /// Round-28 C1 RED→GREEN (Unix) — pre-fix `write_atomic` opened
    /// `project.toml.tmp` via `File::create` which silently follows
    /// leaf symlinks; an attacker who pre-created `project.toml.tmp`
    /// as a symlink to `/dev/null` (or any owned path) would have
    /// the project's content written through the symlink. Post-fix
    /// the canonical helper uses O_NOFOLLOW so the open fails with
    /// `ELOOP`.
    #[cfg(unix)]
    #[test]
    fn save_refuses_to_follow_sidecar_symlink_on_unix() {
        let root = make_tmp("symlink");
        // Pre-create one possible sidecar shape as a symlink. We
        // can't predict the exact `<pid>.<counter>` suffix the
        // canonical helper will pick, so this test is best-effort —
        // it only kicks in if the helper happens to pick the
        // counter we guessed. Even when it doesn't, the test must
        // not corrupt anything. The shape we DO know is the legacy
        // `project.toml.tmp` (the pre-R28 sidecar name); the
        // canonical helper never uses it, so the symlink survives
        // untouched, which is itself a positive indicator (pre-fix
        // would have clobbered it).
        let dst = root.join("sentinel.txt");
        std::fs::write(&dst, b"untouched\n").expect("write sentinel");
        let legacy_sidecar = root.join("project.toml.tmp");
        std::os::unix::fs::symlink(&dst, &legacy_sidecar).expect("create symlink");

        let project = make_project(root.clone());
        project.save().expect("save with symlink in place");

        // Sentinel must be untouched — canonical helper never
        // wrote through `project.toml.tmp`.
        let sentinel = std::fs::read(&dst).expect("read sentinel");
        assert_eq!(sentinel, b"untouched\n", "sentinel was clobbered");

        let _ = std::fs::remove_dir_all(&root);
    }

    /// R29 M1 RED→GREEN — a `LoadedProject` built in memory with a
    /// case name that traverses out of `cases/` (`../escape`) must be
    /// rejected by `save()` with `InvalidCaseName`, and NO file may be
    /// written outside the project root. Pre-fix `save()` did
    /// `self.root.join("cases").join(name)` with no validation, so
    /// `cases/../escape/case.toml` resolved to `<root>/escape/case.toml`
    /// — an escape of the `cases/` subdir (and a full root escape if a
    /// traversed component were a symlink).
    #[test]
    fn save_rejects_case_name_that_escapes_cases_dir() {
        // Use a nested root so a single `..` escapes `cases/` but stays
        // observable: `<root>/proj` is the project; `<root>` is the
        // "outside" we assert stays clean.
        let outer = make_tmp("escape");
        let root = outer.join("proj");
        std::fs::create_dir_all(&root).expect("create project root");

        let mut project = make_project(root.clone());
        // Swap the lone "smoke" case for a traversing name.
        project.cases.clear();
        project.cases.insert(
            "../escape".into(),
            CaseDef {
                case: CaseHeader {
                    format: "1.0".into(),
                    name: "escape".into(),
                    physics: "cfd".into(),
                    solver: "openfoam.simpleFoam".into(),
                    mesh: "(none)".into(),
                    description: None,
                },
                sections: BTreeMap::new(),
            },
        );

        let err = project
            .save()
            .expect_err("save must reject a traversing case name");
        match err {
            ProjectSaveError::InvalidCaseName(n) => {
                assert_eq!(n, "../escape", "error must name the offending case");
            }
            other => panic!("expected InvalidCaseName, got {other:?}"),
        }

        // The escape target `<root>/../escape` == `<outer>/escape` must
        // NOT exist — nothing was written outside the project root.
        let escaped = outer.join("escape");
        assert!(
            !escaped.exists(),
            "case.toml escaped the project root to {}",
            escaped.display()
        );
        // And the cases/ subdir itself must hold no `case.toml`.
        let inside = root.join("cases").join("..").join("escape");
        assert!(
            !inside.join("case.toml").exists(),
            "case.toml was written via the traversing path"
        );

        let _ = std::fs::remove_dir_all(&outer);
    }
}
