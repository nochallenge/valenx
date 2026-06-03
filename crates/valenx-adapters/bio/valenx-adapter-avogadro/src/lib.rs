//! # valenx-adapter-avogadro
//!
//! Adapter for [Avogadro 2](https://two.avogadro.cc/) — the
//! cross-platform, Python-scriptable molecular editor and renderer
//! that succeeded the original Avogadro. Renders chemistry structures
//! (CML / MOL / XYZ / PDB) to publication-quality images, edits
//! molecules from scripts, and runs computational chemistry workflows
//! against MMFF94 / UFF / GFN-FF backends.
//!
//! **Phase 24 — script-driven headless renderer.** The user supplies a
//! `render.py` (or whatever filename) referenced from
//! `[bio.avogadro].script` in `case.toml`. `prepare()` stages the
//! script (and optional structure file) into the workdir; `run()`
//! invokes `avogadro2 --script <script>` with `--no-gui` when running
//! headlessly via the shared subprocess runner.
//!
//! On `collect()` we walk the workdir for `.png` (rendered images) and
//! `.cml` / `.mol` / `.xyz` (exported structures the script wrote
//! out).

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

use crate::case_input::AvogadroInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(AvogadroAdapter::new())
}

pub struct AvogadroAdapter;

impl AvogadroAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for AvogadroAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "avogadro";
/// Avogadro 2's CLI binary. `avogadro2` is the canonical install name
/// (the `2` distinguishes it from the legacy Avogadro 1 line, which
/// some distros still ship as `avogadro`).
const BINARIES: &[&str] = &["avogadro2"];

impl Adapter for AvogadroAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "Avogadro 2",
            // Avogadro 2 ships on a 1.97.x calendar-ish track; the
            // 1.97 line is the first with the modern `--script` /
            // `--no-gui` headless surface this adapter targets.
            // Upper bound 2.0 reserves room for the eventual major
            // bump.
            version_range: VersionRange {
                min_inclusive: Version::new(1, 97, 0),
                max_exclusive: Version::new(2, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            tool_license: "GPL-2.0-or-later",
            docs_url: "https://two.avogadro.cc/docs/",
            homepage_url: "https://two.avogadro.cc/",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(BINARIES) {
            Some(binary_path) => {
                let found_version = detect_tool_version_semver(&binary_path, &["--version", "-v"]);
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
                hint: "Avogadro 2 (1.97+) required; download from \
                       https://two.avogadro.cc/install.html"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = AvogadroInput::from_case_dir(&case.path)?;

        fs::create_dir_all(workdir)?;

        // Stage the user-supplied Python script. Resolved against the
        // case directory when relative — same convention as every
        // other Phase 17/24 bio adapter. `confined_join` rejects
        // absolute paths and `..` traversal that would otherwise let a
        // malicious case bundle stage arbitrary files into the workdir.
        let source_script = confined_join(&case.path, &input.script)?;
        if !source_script.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.avogadro].script `{}` not found (resolved {})",
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
                        "[bio.avogadro].script path `{}` has no filename",
                        input.script.display()
                    ),
                })?;
        let dest_script = workdir.join(script_filename);
        if source_script != dest_script {
            fs::copy(&source_script, &dest_script)?;
        }

        // Stage the optional structure file. Avogadro can load CML /
        // MOL / XYZ / PDB and many more.
        let structure_filename: Option<OsString> = if let Some(ref structure) = input.structure {
            let source_structure = confined_join(&case.path, structure)?;
            if !source_structure.is_file() {
                return Err(AdapterError::InvalidCase {
                    case_path: case.path.join("case.toml"),
                    reason: format!(
                        "[bio.avogadro].structure `{}` not found (resolved {})",
                        structure.display(),
                        source_structure.display()
                    ),
                });
            }
            let fname = structure
                .file_name()
                .ok_or_else(|| AdapterError::InvalidCase {
                    case_path: case.path.join("case.toml"),
                    reason: format!(
                        "[bio.avogadro].structure path `{}` has no filename",
                        structure.display()
                    ),
                })?
                .to_os_string();
            let dest = workdir.join(&fname);
            if source_structure != dest {
                fs::copy(&source_structure, &dest)?;
            }
            Some(fname)
        } else {
            None
        };

        let binary_path = find_on_path(BINARIES).ok_or_else(|| AdapterError::ToolNotInstalled {
            name: INFO_ID,
            hint: "Avogadro 2 (1.97+) required; download from \
                       https://two.avogadro.cc/install.html"
                .into(),
        })?;

        // Compose the invocation:
        //   avogadro2 --script <script> [<structure>] [--no-gui] [extras...]
        let mut native_command: Vec<OsString> = vec![
            binary_path.into_os_string(),
            OsString::from("--script"),
            OsString::from(script_filename),
        ];
        if let Some(structure) = structure_filename {
            native_command.push(structure);
        }
        if input.headless {
            native_command.push(OsString::from("--no-gui"));
        }
        for arg in &input.extra_args {
            native_command.push(OsString::from(arg));
        }

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            // Rendering jobs typically finish in seconds to a couple of
            // minutes; complex animations stretch to tens of minutes.
            // 30 minutes covers the long tail without being absurd.
            estimated_runtime: Some(Duration::from_secs(30 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting Avogadro 2", |line| {
            let mut hint = subprocess::Hint::default();
            // Avogadro Python scripts can emit arbitrary text; pick a
            // few weak signals as best-effort progress hints. A
            // mismatch just leaves the spinner alone.
            if line.contains("Saving") || line.contains("Saved") {
                hint.progress = Some((75.0, line.to_string()));
            } else if line.contains("Rendering") {
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
        // Provenance: hash the staged script (the canonical
        // "this case is configured this way" input).
        let script_path = first_script_in_workdir(&job.workdir);
        let case_hash_input = script_path
            .clone()
            .unwrap_or_else(|| job.workdir.join("case.toml"));
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "Avogadro 2",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);

        let entries = match fs::read_dir(&job.workdir) {
            Ok(e) => e,
            Err(e) => {
                tracing::warn!(target: "valenx-avogadro", ?e, "workdir read failed");
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
                // Rendered images.
                Some("png") => (ArtifactKind::Image, "Avogadro 2 render".to_string()),
                // Exported structure files. CML is Avogadro's native
                // format; MOL / XYZ are common alternatives.
                Some("cml") | Some("mol") | Some("xyz") => (
                    ArtifactKind::Native,
                    "Avogadro 2 exported structure".to_string(),
                ),
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
        Capabilities {
            capabilities: Vec::new(),
            ribbon_contributions: vec!["bio.avogadro.render"],
        }
    }
}

/// Lift the staged Python script out of the workdir for provenance
/// hashing. Returns the lexicographically-first `.py` file at the top
/// level, or `None` if none exists yet.
fn first_script_in_workdir(workdir: &Path) -> Option<PathBuf> {
    let entries = fs::read_dir(workdir).ok()?;
    let mut hits: Vec<PathBuf> = entries
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn info_is_bio_domain() {
        let info = AvogadroAdapter::new().info();
        assert_eq!(info.id, "avogadro");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "GPL-2.0-or-later");
        assert_eq!(info.display_name, "Avogadro 2");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = AvogadroAdapter::new().info();
        // Avogadro 2 1.97+ for the modern --script / --no-gui surface;
        // upper bound 2.0 reserves room for the major bump.
        assert_eq!(info.version_range.min_inclusive, Version::new(1, 97, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(2, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = AvogadroAdapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.avogadro.render"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = AvogadroAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }
}
