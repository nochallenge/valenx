//! # valenx-adapter-biopython
//!
//! Adapter for [Biopython](https://biopython.org/) — the Python
//! library for bioinformatics: sequence I/O, BLAST parsing, motif
//! analysis, structural biology helpers.
//!
//! **Phase 17 — subprocess wrapper for user-provided scripts.**
//! The adapter doesn't generate Python; the user supplies an
//! `analyse.py` (or whatever filename) referenced from
//! `[bio.biopython].script` in `case.toml`. `prepare()` stages the
//! script into the workdir and `run()` invokes
//! `python <script>` via the shared subprocess runner.
//!
//! Future work: bundled-recipe variants (read FASTA -> ORF
//! prediction, multiple-sequence-alignment scoring, etc.) that
//! generate the Python the way `valenx-adapter-cantera` does. The
//! plumbing here is the same shape the recipes will graduate into.

#![forbid(unsafe_code)]
#![allow(missing_docs)]

pub mod case_input;

use std::ffi::OsString;
use std::fs;
use std::path::Path;
use std::time::Duration;

use semver::Version;

use valenx_core::{
    adapter_helpers::{confined_join, find_on_path, live_provenance},
    error::RunPhase,
    subprocess, Adapter, AdapterError, AdapterInfo, Capabilities, Case, LicenseMode, Physics,
    PreparedJob, ProbeReport, RunContext, RunReport, VersionRange,
};
use valenx_fields::{
    artifact::{Artifact, ArtifactKind},
    Results,
};

use crate::case_input::BiopythonInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(BiopythonAdapter::new())
}

pub struct BiopythonAdapter;

impl BiopythonAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for BiopythonAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "biopython";
/// Python interpreter candidates. `python3` first because on Linux
/// `python` may still be Python 2 on legacy distros; on Windows
/// `python` typically resolves to the Windows Store / 3.x install.
const PYTHON_BINARIES: &[&str] = &["python3", "python"];

impl Adapter for BiopythonAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "Biopython",
            version_range: VersionRange {
                // Biopython 1.80 is the first release that fully
                // supports Python 3.9+ and ships the modern
                // `Bio.PDB` / `Bio.SeqIO` API surface we lean on.
                min_inclusive: Version::new(1, 80, 0),
                max_exclusive: Version::new(2, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            tool_license: "BSD-3-Clause",
            docs_url: "https://biopython.org/docs/",
            homepage_url: "https://biopython.org/",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(PYTHON_BINARIES) {
            Some(binary_path) => {
                // Try `import Bio; print(Bio.__version__)` first —
                // that's the only string that confirms Biopython
                // itself is installed (vs. just having Python on
                // PATH). Fall back to `python --version` if the
                // import fails so we can still surface a useful
                // probe state.
                let found_version = detect_biopython_version(&binary_path);
                let mut warnings = Vec::new();
                if found_version.is_none() {
                    warnings.push(
                        "probe found `python` on PATH but could not import \
                         `Bio` — install Biopython with `pip install \
                         biopython` for runs to succeed"
                            .into(),
                    );
                }
                Ok(ProbeReport {
                    ok: true,
                    found_version,
                    binary_path: Some(binary_path),
                    warnings,
                    required_env: Vec::new(),
                })
            }
            None => Err(AdapterError::ToolNotInstalled {
                name: INFO_ID,
                hint: "Python 3.9+ with Biopython installed; `pip install \
                       biopython` after ensuring python3 is on PATH"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = BiopythonInput::from_case_dir(&case.path)?;

        fs::create_dir_all(workdir)?;

        // Resolve the script path against the case directory.
        // `confined_join` rejects absolute paths and `..` traversal so
        // the staged copy stays confined to the case directory.
        let source_script = confined_join(&case.path, &input.script)?;
        if !source_script.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.biopython].script `{}` not found (resolved {})",
                    input.script.display(),
                    source_script.display()
                ),
            });
        }
        let script_filename =
            input
                .script
                .file_name()
                .ok_or_else(|| AdapterError::InvalidCase {
                    case_path: case.path.join("case.toml"),
                    reason: format!(
                        "[bio.biopython].script path `{}` has no filename",
                        input.script.display()
                    ),
                })?;
        let dest_script = workdir.join(script_filename);
        if source_script != dest_script {
            fs::copy(&source_script, &dest_script)?;
        }

        // Resolve the Python binary. If the user pinned a specific
        // interpreter via `python = "..."`, honor it; otherwise
        // walk PATH.
        // Round-4 security: validate python interpreter spec
        // against the allow-list AND resolve to a real binary
        // in one step. Closes the arbitrary-binary-exec class
        // that round-3 only patched in 8 of the 48 affected
        // adapters.
        let binary_path = valenx_core::adapter_helpers::resolve_python_binary(
            &input.python,
            PYTHON_BINARIES,
        )
        // Round-5: do NOT rewrap as ToolNotInstalled — the resolver
        // returns InvalidCase for allow-list rejections (which a hint
        // string would have hidden) and a clear Other for PATH lookup
        // failures. Pass the error through unchanged.
        ?;

        let native_command: Vec<OsString> = vec![
            binary_path.into_os_string(),
            OsString::from(script_filename),
        ];

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            // Bio scripts can run for a while — sequence alignment,
            // structure prediction, BLAST searches all routinely
            // hit minutes-to-tens-of-minutes. Pick 30 minutes as a
            // reasonable default; long-running cases override via
            // their own progress reporting.
            estimated_runtime: Some(Duration::from_secs(30 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting Biopython", |line| {
            let mut hint = subprocess::Hint::default();
            // Convention: the user-supplied script can emit a
            // sentinel line `[valenx] biopython done` to signal
            // completion before exit; we lift that to a 95% progress
            // tick so the UI doesn't sit at "indeterminate" until
            // the process actually exits.
            if line.contains("[valenx] biopython done") {
                hint.progress = Some((95.0, line.to_string()));
            } else if line.contains("Traceback") || line.contains("Error") {
                hint.warning = Some(line.trim().to_string());
            }
            hint
        })?;
        Ok(RunReport {
            exit_code: report.exit_code,
            wall_time: report.wall_time,
            converged: Some(true),
            residual_history: Vec::new(),
            warnings: report.warnings,
            final_phase: Some(RunPhase::Shutdown),
        })
    }

    fn collect(&self, job: &PreparedJob) -> Result<Results, AdapterError> {
        // Provenance: hash the staged script (the canonical
        // "this case is configured this way" input). We don't
        // know the user's mesh / lock files so leave those empty.
        let script_path = first_script_in_workdir(&job.workdir);
        let case_hash_input = script_path
            .clone()
            .unwrap_or_else(|| job.workdir.join("case.toml"));
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "Biopython",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);

        // Walk the workdir top level and classify any output files
        // the user's script left behind. We deliberately don't
        // recurse — Biopython scripts that need nested output
        // directories will surface their key artefacts at the top
        // level by convention.
        let entries = match fs::read_dir(&job.workdir) {
            Ok(e) => e,
            Err(e) => {
                tracing::warn!(target: "valenx-biopython", ?e, "workdir read failed");
                return Ok(results);
            }
        };
        let mut artefacts: Vec<Artifact> = Vec::new();
        for entry in entries.flatten() {
            let path = entry.path();
            if !path.is_file() {
                continue;
            }
            let ext = path
                .extension()
                .and_then(|s| s.to_str())
                .map(|s| s.to_ascii_lowercase());
            let (kind, label) = match ext.as_deref() {
                Some("fasta") | Some("fa") | Some("faa") | Some("fna") => {
                    // Soft-validate FASTA but never fail — a script
                    // that emits a partial fasta during a long run
                    // shouldn't make `collect()` blow up.
                    //
                    // Round-22 M2: cap the read at MAX_PDB_FILE_BYTES
                    // (256 MiB) so a poisoned workdir with a multi-GB
                    // `.fasta` can't OOM `collect()` before the parser
                    // runs.
                    let label = match valenx_core::io_caps::read_capped_to_string(
                        &path,
                        valenx_core::io_caps::MAX_PDB_FILE_BYTES as usize,
                    ) {
                        Ok(text) => match valenx_bio::format::fasta::read(
                            &text,
                            valenx_bio::Alphabet::Protein,
                        ) {
                            Ok(seqs) => format!("FASTA ({} sequences)", seqs.len()),
                            Err(e) => format!(
                                "FASTA (parse warning: {})",
                                e.to_string().lines().next().unwrap_or("invalid")
                            ),
                        },
                        Err(_) => "FASTA".to_string(),
                    };
                    (ArtifactKind::Other, label)
                }
                Some("json") => (ArtifactKind::Other, "JSON output".to_string()),
                Some("csv") | Some("tsv") => (ArtifactKind::Tabular, "Tabular output".to_string()),
                Some("py") => (ArtifactKind::Other, "Biopython script".to_string()),
                _ => continue,
            };
            artefacts.push(Artifact {
                path,
                kind,
                checksum: None,
                label,
            });
        }
        artefacts.sort_by(|a, b| a.path.cmp(&b.path));
        results.artifacts = artefacts;
        Ok(results)
    }

    fn capabilities(&self) -> Capabilities {
        // The bio-specific Capability variants land in a follow-up
        // task; for now we publish an empty capability vector and
        // a single ribbon contribution so the registry can wire the
        // adapter in without crashing the UI's capability-index
        // builder.
        Capabilities {
            capabilities: Vec::new(),
            ribbon_contributions: vec!["bio.biopython.script"],
        }
    }
}

/// Lift the staged Python script out of a workdir for provenance
/// hashing. Returns the lexicographically-first `.py` file at the
/// top level, or `None` if none exists yet.
fn first_script_in_workdir(workdir: &Path) -> Option<std::path::PathBuf> {
    let entries = fs::read_dir(workdir).ok()?;
    let mut hits: Vec<std::path::PathBuf> = entries
        .flatten()
        .map(|e| e.path())
        .filter(|p| {
            p.extension()
                .and_then(|s| s.to_str())
                .map(|s| s.eq_ignore_ascii_case("py"))
                .unwrap_or(false)
        })
        .collect();
    hits.sort();
    hits.into_iter().next()
}

/// Run `python -c "import Bio; print(Bio.__version__)"` and parse a
/// `semver::Version` out of stdout. Returns `None` on any failure
/// (interpreter unusable, Biopython not importable, version string
/// malformed) — `probe()` falls back to a "Biopython not importable"
/// warning in that case.
fn detect_biopython_version(python_binary: &Path) -> Option<Version> {
    let output = std::process::Command::new(python_binary)
        .arg("-c")
        .arg("import Bio; print(Bio.__version__)")
        .stdin(std::process::Stdio::null())
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::piped())
        .output()
        .ok()?;
    if !output.status.success() {
        return None;
    }
    let stdout = String::from_utf8(output.stdout).ok()?;
    let raw = valenx_core::adapter_helpers::extract_semver(&stdout)?;
    let dots = raw.chars().filter(|c| *c == '.').count();
    let normalised: String = match dots {
        0 => format!("{raw}.0.0"),
        1 => format!("{raw}.0"),
        _ => raw,
    };
    Version::parse(&normalised).ok()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn info_is_bio_domain() {
        let info = BiopythonAdapter::new().info();
        assert_eq!(info.id, "biopython");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "BSD-3-Clause");
        assert_eq!(info.display_name, "Biopython");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = BiopythonAdapter::new().info();
        // We support Biopython >= 1.80 (modern Bio.PDB API) and
        // expect to revisit the upper bound when a major bump lands.
        assert_eq!(info.version_range.min_inclusive, Version::new(1, 80, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(2, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = BiopythonAdapter::new().capabilities();
        // Capability variants land in a future task; ribbon
        // contributions are already enough for the registry to
        // surface the adapter.
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.biopython.script"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = BiopythonAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }

    #[test]
    fn first_script_in_workdir_picks_lexicographic_first() {
        let tmp = std::env::temp_dir().join(format!(
            "valenx-biopython-script-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(&tmp).unwrap();
        fs::write(tmp.join("z_late.py"), b"# placeholder").unwrap();
        fs::write(tmp.join("a_first.py"), b"# placeholder").unwrap();
        fs::write(tmp.join("not_python.txt"), b"placeholder").unwrap();
        let f = first_script_in_workdir(&tmp).expect("found");
        assert!(f.ends_with("a_first.py"));
        let _ = fs::remove_dir_all(&tmp);
    }

    #[test]
    fn first_script_in_workdir_returns_none_when_empty() {
        let tmp = std::env::temp_dir().join(format!(
            "valenx-biopython-empty-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(&tmp).unwrap();
        fs::write(tmp.join("notes.md"), b"placeholder").unwrap();
        assert!(first_script_in_workdir(&tmp).is_none());
        let _ = fs::remove_dir_all(&tmp);
    }
}
