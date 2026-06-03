//! # valenx-adapter-pymol
//!
//! Adapter for [PyMOL](https://pymol.org/) — Schrödinger's open-source
//! molecular visualisation system, the go-to companion to ChimeraX
//! for cartoon / surface / ray-traced rendering of proteins,
//! ligands, and small molecules. PyMOL reads PDB / CIF / many other
//! structural formats and executes `.pml` command scripts that pose
//! the camera, set styles, and produce publication-quality output.
//!
//! **Phase 23 — subprocess wrapper for user-provided PyMOL scripts.**
//! Sister adapter to ChimeraX: same shape, same script-driven flow.
//! The user supplies `render.pml` (or whatever filename) referenced
//! from `[bio.pymol].script` in `case.toml`. `prepare()` stages the
//! script into the workdir and `run()` invokes
//! `pymol -c -q <script> ...extras` via the shared subprocess
//! runner. `-c` keeps the run headless (no GUI window); `-q`
//! suppresses the startup banner. Both flags default to true so
//! headless CI takes the happy path; either flips to false on the
//! rare interactive / banner-on workstation run.
//!
//! On `collect()` we walk the workdir for PyMOL's customary outputs:
//! `.png` rendered images, `.pse` session files, and any `.pdb` /
//! `.cif` structures the script wrote out via `save`. Images surface
//! as `Image`; sessions and structure files surface as `Native` so
//! the user can re-open them in PyMOL (or hand them off to ChimeraX
//! / VMD).
//!
//! The license name we publish for `tool_license` is
//! `BSD-3-Clause-Open-Source` — the open-source build at
//! `github.com/schrodinger/pymol-open-source` ships under BSD-3
//! whereas Schrödinger's commercial PyMOL distribution carries
//! proprietary terms. We surface the open-source identifier
//! verbatim so license-aware tooling can distinguish the two
//! distributions.

#![forbid(unsafe_code)]
#![allow(missing_docs)]

pub mod case_input;

use std::ffi::OsString;
use std::fs;
use std::path::{Path, PathBuf};
use std::time::Duration;

use semver::Version;

use valenx_core::{
    adapter_helpers::{confined_join, detect_tool_version_semver, find_on_path, live_provenance},
    error::RunPhase,
    subprocess, Adapter, AdapterError, AdapterInfo, Capabilities, Case, LicenseMode, Physics,
    PreparedJob, ProbeReport, RunContext, RunReport, VersionRange,
};
use valenx_fields::{
    artifact::{Artifact, ArtifactKind},
    Results,
};

use crate::case_input::PymolInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(PymolAdapter::new())
}

pub struct PymolAdapter;

impl PymolAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for PymolAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "pymol";
/// PyMOL binary candidates. Bioconda, the open-source GitHub build,
/// and the Schrödinger commercial installer all ship the launcher
/// under the canonical lowercase name.
const BINARIES: &[&str] = &["pymol"];

impl Adapter for PymolAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "PyMOL",
            // PyMOL 2.5 (2021) is the floor we test against — it
            // carries the modernised CLI flags we lean on (`-c` / `-q`
            // semantics, the open-source build's stable plugin
            // surface). The upper bound 3.0 reserves room for the
            // long-rumoured 3.x line; bump when it lands.
            version_range: VersionRange {
                min_inclusive: Version::new(2, 5, 0),
                max_exclusive: Version::new(3, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            // The open-source build at github.com/schrodinger/pymol-open-source
            // ships under BSD-3-Clause. We surface a custom identifier
            // (BSD-3-Clause-Open-Source) so license-aware tooling can
            // distinguish this from Schrödinger's proprietary
            // commercial PyMOL distribution.
            tool_license: "BSD-3-Clause-Open-Source",
            docs_url: "https://pymol.org/dokuwiki/",
            homepage_url: "https://github.com/schrodinger/pymol-open-source",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(BINARIES) {
            Some(binary_path) => {
                // PyMOL prints its version on `--version` (stdout)
                // and on `-c -V` (rare). The combined scanner picks
                // the SemVer up cleanly from `--version`.
                let found_version = detect_tool_version_semver(&binary_path, &["--version"]);
                Ok(ProbeReport {
                    ok: true,
                    found_version,
                    binary_path: Some(binary_path),
                    warnings: Vec::new(),
                    required_env: Vec::new(),
                })
            }
            None => Err(AdapterError::ToolNotInstalled {
                name: INFO_ID,
                hint: "PyMOL 2.5+ required; install the open-source build via \
                       `conda install -c conda-forge pymol-open-source` or download \
                       from https://github.com/schrodinger/pymol-open-source"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = PymolInput::from_case_dir(&case.path)?;

        fs::create_dir_all(workdir)?;

        // Resolve the script path against the case directory. Mirrors
        // ChimeraX's `script = "render.cxc"` next to `case.toml`
        // convention. `confined_join` rejects absolute paths and `..`
        // traversal so a malicious case bundle can't smuggle arbitrary
        // host files into the workdir.
        let source_script = confined_join(&case.path, &input.script)?;
        if !source_script.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.pymol].script `{}` not found (resolved {})",
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
                        "[bio.pymol].script path `{}` has no filename",
                        input.script.display()
                    ),
                })?;
        let dest_script = workdir.join(script_filename);
        if source_script != dest_script {
            fs::copy(&source_script, &dest_script)?;
        }

        let binary_path = find_on_path(BINARIES).ok_or_else(|| AdapterError::ToolNotInstalled {
            name: INFO_ID,
            hint: "PyMOL 2.5+ required; install the open-source build via \
                       `conda install -c conda-forge pymol-open-source` or download \
                       from https://github.com/schrodinger/pymol-open-source"
                .into(),
        })?;

        // Build the command. `-c` runs PyMOL headlessly; `-q`
        // suppresses the startup banner. The script filename is
        // passed as a positional argument; PyMOL detects the `.pml`
        // extension and executes it as a command file.
        let mut native_command: Vec<OsString> = vec![binary_path.into_os_string()];
        if input.nogui {
            native_command.push(OsString::from("-c"));
        }
        if input.quiet {
            native_command.push(OsString::from("-q"));
        }
        native_command.push(OsString::from(script_filename));
        for arg in &input.extra_args {
            native_command.push(OsString::from(arg));
        }

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            // PyMOL rendering jobs typically finish in seconds to a
            // couple of minutes; ray-traced animations can run for
            // 30+ minutes. Same generous 30-minute default as
            // ChimeraX.
            estimated_runtime: Some(Duration::from_secs(30 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting PyMOL", |line| {
            let mut hint = subprocess::Hint::default();
            // PyMOL banners are loose — the script can emit arbitrary
            // text via `print` / `cmd.feedback`. We pick three weak
            // signals as best-effort progress hints:
            //   * "Ray:" / "ray:" — the ray-tracer started a render
            //   * "Saving" / "PNG" — an output was written
            //   * "PyMOL>" / "exit" — the session is shutting down
            // These are heuristics; mismatches just leave the spinner
            // alone.
            if line.contains("exit") || line.contains("PyMOL: normal program termination") {
                hint.progress = Some((95.0, line.to_string()));
            } else if line.contains("Saving") || line.contains("PNG") {
                hint.progress = Some((75.0, line.to_string()));
            } else if line.contains("Ray:") || line.contains("ray:") {
                hint.progress = Some((50.0, line.to_string()));
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
        // Provenance: hash the staged .pml script (the canonical
        // "this case is configured this way" input).
        let script_path = first_script_in_workdir(&job.workdir);
        let case_hash_input = script_path
            .clone()
            .unwrap_or_else(|| job.workdir.join("case.toml"));
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "PyMOL",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);

        // Walk the workdir top level. PyMOL scripts conventionally
        // write outputs to the working directory; deeply nested
        // outputs surface only via the script's own explicit `cd` /
        // `save` paths.
        let entries = match fs::read_dir(&job.workdir) {
            Ok(e) => e,
            Err(e) => {
                tracing::warn!(target: "valenx-pymol", ?e, "workdir read failed");
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
                // Rendered images — `png` / `ray` outputs from PyMOL.
                Some("png") => (ArtifactKind::Image, "PyMOL render".to_string()),
                // .pse — PyMOL session file. Re-openable via
                // `load session.pse` in another PyMOL run.
                Some("pse") => (ArtifactKind::Native, "PyMOL session".to_string()),
                // Structure outputs the script wrote (typically via
                // `save out.pdb` or `save out.cif`). We surface as
                // Native so other viewers (ChimeraX, VMD) can pick
                // them up.
                Some("pdb") | Some("cif") => {
                    (ArtifactKind::Native, "PyMOL exported structure".to_string())
                }
                Some("pml") => (ArtifactKind::Other, "PyMOL command script".to_string()),
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
        // Bio-specific Capability variants land in a follow-up task;
        // ribbon contributions are already enough for the registry to
        // surface the adapter without crashing the UI's
        // capability-index builder.
        Capabilities {
            capabilities: Vec::new(),
            ribbon_contributions: vec!["bio.pymol.render"],
        }
    }
}

/// Lift the staged PyMOL script out of a workdir for provenance
/// hashing. Returns the lexicographically-first `.pml` file at the
/// top level, or `None` if none exists yet.
fn first_script_in_workdir(workdir: &Path) -> Option<PathBuf> {
    let entries = fs::read_dir(workdir).ok()?;
    let mut hits: Vec<PathBuf> = entries
        .flatten()
        .map(|e| e.path())
        .filter(|p| {
            p.extension()
                .and_then(|s| s.to_str())
                .map(|s| s.eq_ignore_ascii_case("pml"))
                .unwrap_or(false)
        })
        .collect();
    hits.sort();
    hits.into_iter().next()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn info_is_bio_domain() {
        let info = PymolAdapter::new().info();
        assert_eq!(info.id, "pymol");
        assert_eq!(info.physics, &[Physics::Bio]);
        // The license identifier must surface the open-source
        // distinguisher so license-aware tooling can tell this apart
        // from Schrödinger's commercial PyMOL.
        assert_eq!(info.tool_license, "BSD-3-Clause-Open-Source");
        assert_eq!(info.display_name, "PyMOL");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = PymolAdapter::new().info();
        // PyMOL >= 2.5 (stable open-source CLI surface); upper bound
        // 3.0 reserves room for the next major.
        assert_eq!(info.version_range.min_inclusive, Version::new(2, 5, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(3, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = PymolAdapter::new().capabilities();
        // Capability variants land in a future task; ribbon
        // contributions are already enough for the registry to
        // surface the adapter.
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.pymol.render"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = PymolAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }
}
