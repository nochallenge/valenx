//! # valenx-adapter-vmd
//!
//! Adapter for [VMD](https://www.ks.uiuc.edu/Research/vmd/) — UIUC's
//! Visual Molecular Dynamics renderer, the canonical viewer for
//! large MD trajectories from NAMD / GROMACS / LAMMPS / AMBER.
//! VMD's strength over PyMOL / ChimeraX is trajectory handling: it
//! streams multi-million-atom DCDs frame-by-frame and exposes a
//! Tcl scripting surface that drives loading, representation
//! styling, frame iteration, and image / data export.
//!
//! **Phase 23 — subprocess wrapper for user-provided VMD scripts.**
//! Sister adapter to ChimeraX / PyMOL, with one structural addition:
//! VMD setups commonly need a topology file (`.psf`, `.gro`,
//! `.parm7`) loaded before the script runs so the script can
//! `mol addfile` the trajectory. The case-input shape carries an
//! optional `structure` path that the adapter passes as a positional
//! argument before `-e <script>`.
//!
//! `prepare()` stages the script (and optional structure file) into
//! the workdir and `run()` invokes
//! `vmd -dispdev text -e <script>` via the shared subprocess
//! runner. `-dispdev text` selects the headless (no-OpenGL)
//! renderer; defaults to true so CI takes the happy path. Flip
//! `headless = false` to drive the OpenGL renderer on a workstation.
//!
//! On `collect()` we walk the workdir for VMD's customary outputs:
//! `.png` / `.tga` / `.bmp` rendered frames, `.pdb` / `.gro`
//! exported structures, and `.dat` / `.csv` analysis tables (RMSD,
//! RDF, etc.) the script wrote out.
//!
//! ## License flag
//!
//! VMD ships under a custom non-OSS license that restricts use to
//! non-commercial / academic contexts. We surface this accurately
//! via a `tool_license` value of `VMD-License` and emit a probe
//! warning when the binary is found so downstream tooling and
//! end-users get a clear "check your license" signal before
//! redistributing renders or derived data. The probe-warning text
//! contains the literal string `"academic"` as a stable anchor for
//! tests and downstream filters.

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

use crate::case_input::VmdInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(VmdAdapter::new())
}

pub struct VmdAdapter;

impl VmdAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for VmdAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "vmd";
/// VMD binary candidates. The launcher script installs as the
/// canonical lowercase `vmd` on every supported platform.
const BINARIES: &[&str] = &["vmd"];

/// The probe-warning surfaced whenever VMD is detected. Anchors a
/// stable "academic / non-commercial only" reminder for downstream
/// tooling and tests; the literal string `"academic"` is part of
/// the asserted contract.
const LICENSE_WARNING: &str = "VMD is licensed for non-commercial / academic use only. \
     Confirm your use case complies with the upstream license \
     before redistributing renders or derived data.";

impl Adapter for VmdAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "VMD",
            // VMD 1.9.x has been the stable line for over a decade;
            // 1.9 (2014) is the floor we test against and covers
            // every recent release through 1.9.4. The upper bound
            // 2.0 reserves room for the long-rumoured 2.x line.
            version_range: VersionRange {
                min_inclusive: Version::new(1, 9, 0),
                max_exclusive: Version::new(2, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            // VMD's terms aren't a recognised SPDX identifier; the
            // closest accurate label is the project's own
            // "VMD-License" name. Surfacing it here (instead of
            // mislabeling as MIT / BSD) keeps license-aware tooling
            // honest.
            tool_license: "VMD-License",
            docs_url: "https://www.ks.uiuc.edu/Research/vmd/current/ug/",
            homepage_url: "https://www.ks.uiuc.edu/Research/vmd/",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(BINARIES) {
            Some(binary_path) => {
                // VMD prints its version on `--version` (some
                // distributions) and via the banner on a quick
                // headless start; the combined scanner picks it up
                // from `--version` when available.
                let found_version = detect_tool_version_semver(&binary_path, &["--version"]);
                Ok(ProbeReport {
                    ok: true,
                    found_version,
                    binary_path: Some(binary_path),
                    // Always surface the license reminder when VMD
                    // is detected — it's a custom non-OSS license
                    // and we'd rather over-warn than have a user
                    // ship commercial output without checking.
                    warnings: vec![LICENSE_WARNING.to_string()],
                    required_env: Vec::new(),
                })
            }
            None => Err(AdapterError::ToolNotInstalled {
                name: INFO_ID,
                hint: "VMD 1.9+ required; download from \
                       https://www.ks.uiuc.edu/Research/vmd/ \
                       (registration required, academic-use license)"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = VmdInput::from_case_dir(&case.path)?;

        fs::create_dir_all(workdir)?;

        // Resolve the script path against the case directory.
        // `confined_join` rejects absolute paths and `..` traversal so
        // a malicious case bundle can't smuggle arbitrary host files
        // into the workdir.
        let source_script = confined_join(&case.path, &input.script)?;
        if !source_script.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.vmd].script `{}` not found (resolved {})",
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
                        "[bio.vmd].script path `{}` has no filename",
                        input.script.display()
                    ),
                })?;
        let dest_script = workdir.join(script_filename);
        if source_script != dest_script {
            fs::copy(&source_script, &dest_script)?;
        }

        // Optional topology file: same staging contract as the
        // script. We only stage it if the user supplied one; bare
        // trajectory-only setups leave `structure = None` and skip
        // this entire branch.
        let staged_structure = if let Some(struct_path) = &input.structure {
            let source_struct = confined_join(&case.path, struct_path)?;
            if !source_struct.is_file() {
                return Err(AdapterError::InvalidCase {
                    case_path: case.path.join("case.toml"),
                    reason: format!(
                        "[bio.vmd].structure `{}` not found (resolved {})",
                        struct_path.display(),
                        source_struct.display()
                    ),
                });
            }
            let struct_filename =
                struct_path
                    .file_name()
                    .ok_or_else(|| AdapterError::InvalidCase {
                        case_path: case.path.join("case.toml"),
                        reason: format!(
                            "[bio.vmd].structure path `{}` has no filename",
                            struct_path.display()
                        ),
                    })?;
            let dest_struct = workdir.join(struct_filename);
            if source_struct != dest_struct {
                fs::copy(&source_struct, &dest_struct)?;
            }
            Some(OsString::from(struct_filename))
        } else {
            None
        };

        let binary_path = find_on_path(BINARIES).ok_or_else(|| AdapterError::ToolNotInstalled {
            name: INFO_ID,
            hint: "VMD 1.9+ required; download from \
                       https://www.ks.uiuc.edu/Research/vmd/ \
                       (registration required, academic-use license)"
                .into(),
        })?;

        // Build the command. `-dispdev text` selects the headless
        // renderer (no GUI / OpenGL window). `-e <script>` runs the
        // staged Tcl driver. Optional structure path lands as a
        // positional argument before `-e` so VMD loads it as the
        // initial molecule before the script runs.
        let mut native_command: Vec<OsString> = vec![binary_path.into_os_string()];
        if input.headless {
            native_command.push(OsString::from("-dispdev"));
            native_command.push(OsString::from("text"));
        }
        if let Some(struct_filename) = staged_structure {
            native_command.push(struct_filename);
        }
        native_command.push(OsString::from("-e"));
        native_command.push(OsString::from(script_filename));
        for arg in &input.extra_args {
            native_command.push(OsString::from(arg));
        }

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            // VMD trajectory rendering can run from a few seconds
            // (single-frame snapshot) to multiple hours (full
            // trajectory animation with Tachyon ray-tracing). 1
            // hour is a generous default; the user can extend via
            // executor-side limits for the long tail.
            estimated_runtime: Some(Duration::from_secs(60 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting VMD", |line| {
            let mut hint = subprocess::Hint::default();
            // VMD's banner / progress output is loose. Best-effort
            // hints from common Tcl-script idioms:
            //   * "Loading" / "Reading" — input load progress
            //   * "frame" — trajectory frame iteration
            //   * "rendering" / "Tachyon" — render step started
            //   * "Info)" with "exit" — shutting down
            // Heuristics; mismatches just leave the spinner alone.
            if line.contains("normal exit") || line.contains("VMD exit") {
                hint.progress = Some((95.0, line.to_string()));
            } else if line.contains("rendering") || line.contains("Tachyon") {
                hint.progress = Some((75.0, line.to_string()));
            } else if line.contains("Loading") || line.contains("Reading") {
                hint.progress = Some((25.0, line.to_string()));
            } else if line.contains("ERROR") || line.contains("Tcl Error") {
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
        // Provenance: hash the staged .tcl driver script.
        let script_path = first_script_in_workdir(&job.workdir);
        let case_hash_input = script_path
            .clone()
            .unwrap_or_else(|| job.workdir.join("case.toml"));
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "VMD",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);

        // Walk the workdir top level. VMD scripts conventionally
        // write outputs to the working directory; deeply nested
        // outputs surface only via the script's own explicit `cd` /
        // `render` paths.
        let entries = match fs::read_dir(&job.workdir) {
            Ok(e) => e,
            Err(e) => {
                tracing::warn!(target: "valenx-vmd", ?e, "workdir read failed");
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
                // Rendered images — VMD's `render snapshot ...` /
                // `render Tachyon ...` writes to PNG (via Tachyon),
                // TGA (the native renderer's default), or BMP.
                Some("png") | Some("tga") | Some("bmp") => {
                    (ArtifactKind::Image, "VMD render".to_string())
                }
                // Structure outputs the script wrote (typically via
                // `mol writepdb out.pdb` or `mol writegro out.gro`).
                Some("pdb") | Some("gro") => {
                    (ArtifactKind::Native, "VMD exported structure".to_string())
                }
                // Analysis outputs — RMSD / RDF / contact-map
                // tables that user scripts commonly write via
                // Tcl `puts`.
                Some("dat") | Some("csv") => {
                    (ArtifactKind::Tabular, "VMD analysis data".to_string())
                }
                Some("tcl") => (ArtifactKind::Other, "VMD command script".to_string()),
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
            ribbon_contributions: vec!["bio.vmd.render"],
        }
    }
}

/// Lift the staged VMD driver script out of a workdir for provenance
/// hashing. Returns the lexicographically-first `.tcl` file at the
/// top level, or `None` if none exists yet.
fn first_script_in_workdir(workdir: &Path) -> Option<PathBuf> {
    let entries = fs::read_dir(workdir).ok()?;
    let mut hits: Vec<PathBuf> = entries
        .flatten()
        .map(|e| e.path())
        .filter(|p| {
            p.extension()
                .and_then(|s| s.to_str())
                .map(|s| s.eq_ignore_ascii_case("tcl"))
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
        let info = VmdAdapter::new().info();
        assert_eq!(info.id, "vmd");
        assert_eq!(info.physics, &[Physics::Bio]);
        // The license identifier must surface VMD's custom non-OSS
        // license rather than mislabel as MIT / BSD.
        assert_eq!(info.tool_license, "VMD-License");
        assert_eq!(info.display_name, "VMD");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = VmdAdapter::new().info();
        // VMD >= 1.9 (the stable line for over a decade); upper
        // bound 2.0 reserves room for the next major.
        assert_eq!(info.version_range.min_inclusive, Version::new(1, 9, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(2, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = VmdAdapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.vmd.render"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = VmdAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }

    #[test]
    fn probe_warning_mentions_academic() {
        // The license-flag warning is mandatory: VMD is non-OSS
        // academic-use, and we surface that on every successful
        // probe. The literal "academic" anchor is what downstream
        // tooling and license-aware filters key off — pin it.
        assert!(
            LICENSE_WARNING.contains("academic"),
            "probe warning must contain `academic` anchor; got: {LICENSE_WARNING}"
        );
    }
}
