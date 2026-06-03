//! # valenx-adapter-ilastik
//!
//! Adapter for [Ilastik](https://www.ilastik.org/) — the interactive
//! ML-based pixel / object classification suite (Hamprecht et al,
//! GPL-3.0). Ilastik is the de-facto reference tool for "draw a few
//! brushstrokes, train a Random Forest, segment an entire image
//! stack" workflows in bioimage analysis: users author a trained
//! classifier in the GUI (saved as a `.ilp` project) and then
//! re-run that project headlessly across new image batches.
//!
//! **Phase 40 — sister adapter to Fiji and CellProfiler for the
//! microscopy / bioimage-analysis surface.** Same app-launcher
//! pattern as Fiji: the user supplies the absolute path to the
//! platform-appropriate launcher binary
//! (`run_ilastik.sh` on Linux / macOS, `ilastik.exe` on Windows)
//! via `[bio.ilastik].ilastik_app`, plus the trained `.ilp` project
//! and one or more input images.
//!
//! `prepare()` composes
//! `<ilastik_app> --headless --project=<project>
//!  --output_filename_format=<output_basename>_{nickname}.h5
//!  <input_images...> [extras...]`.
//! The literal `{nickname}` token is Ilastik's per-input substitution
//! placeholder — it must reach Ilastik unmodified so each input
//! image gets a unique disambiguated output filename. `run()` streams
//! the run via the shared subprocess runner. `collect()` walks the
//! workdir for the canonical `<output_basename>*.h5` probability
//! maps, `<output_basename>*.tif` segmentations, plus any `*.log`
//! files Ilastik emits.

#![forbid(unsafe_code)]
#![allow(missing_docs)]

pub mod case_input;

use std::ffi::OsString;
use std::fs;
use std::path::Path;
use std::time::Duration;

use semver::Version;

use valenx_core::{
    adapter_helpers::{find_on_path, live_provenance},
    error::RunPhase,
    subprocess, Adapter, AdapterError, AdapterInfo, Capabilities, Case, LicenseMode, Physics,
    PreparedJob, ProbeReport, RunContext, RunReport, VersionRange,
};
use valenx_fields::{
    artifact::{Artifact, ArtifactKind},
    Results,
};

use crate::case_input::IlastikInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(IlastikAdapter::new())
}

pub struct IlastikAdapter;

impl IlastikAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for IlastikAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "ilastik";
/// Probe binary candidates. `ilastik` is the bare symlink some
/// distros (Linux package managers, conda) create; `run_ilastik.sh`
/// is the canonical launcher in the portable Linux / macOS bundle;
/// `ilastik.exe` is the Windows launcher.
const BINARIES: &[&str] = &["ilastik", "run_ilastik.sh", "ilastik.exe"];

impl Adapter for IlastikAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "Ilastik",
            // Ilastik 1.4 (2022) is the modern stable line —
            // first to ship the unified Workflow API + the
            // tiktorch / neural-network workflow integration.
            // 2.0 reserves room for the next major bump.
            version_range: VersionRange {
                min_inclusive: Version::new(1, 4, 0),
                max_exclusive: Version::new(2, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            tool_license: "GPL-3.0",
            docs_url: "https://www.ilastik.org/documentation/",
            homepage_url: "https://www.ilastik.org/",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(BINARIES) {
            Some(binary_path) => Ok(ProbeReport {
                ok: true,
                // Ilastik's `--version` output format has shifted
                // across releases; we surface no version here —
                // the user pins the Ilastik release implicitly by
                // the launcher path they point at via
                // `[bio.ilastik].ilastik_app`.
                found_version: None,
                binary_path: Some(binary_path),
                warnings: Vec::new(),
                required_env: Vec::new(),
            }),
            None => {
                // Ilastik isn't on PATH, but the user can still
                // supply the launcher path via
                // `[bio.ilastik].ilastik_app` and the run will
                // succeed. Report `ok: true` with an actionable
                // warning rather than `ToolNotInstalled` so case
                // execution isn't blocked when the user has
                // opted into the explicit-path workflow.
                Ok(ProbeReport {
                    ok: true,
                    found_version: None,
                    binary_path: None,
                    warnings: vec!["Ilastik not found on PATH; download from \
                         https://www.ilastik.org/download.html and \
                         add bin to PATH"
                        .into()],
                    required_env: Vec::new(),
                })
            }
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = IlastikInput::from_case_dir(&case.path)?;

        // Round-4 security: reject `output_basename = "../etc/passwd"`
        // and friends before the value flows into any path join.
        // Same pattern as the round-3 fix in bionetgen/iqtree/art/fasttree.
        valenx_core::adapter_helpers::validate_output_basename(
            &input.output_basename,
            "[bio.ilastik].output_basename",
        )
        .map_err(|e| AdapterError::InvalidCase {
            case_path: case.path.join("case.toml"),
            reason: format!("{e}"),
        })?;

        fs::create_dir_all(workdir)?;

        // Resolve the launcher path against the case directory if
        // relative. Almost always absolute (Ilastik installs land
        // under /opt or C:\Program Files), but support the relative
        // form too.
        // Round-9 classification: KEEP `case.path.join` — the
        // Ilastik launcher is an admin-managed system binary. The
        // project file and input images below go through `confined_join`.
        let resolved_app = if input.ilastik_app.is_absolute() {
            input.ilastik_app.clone()
        } else {
            case.path.join(&input.ilastik_app)
        };
        if !resolved_app.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.ilastik].ilastik_app `{}` not found (resolved {})",
                    input.ilastik_app.display(),
                    resolved_app.display()
                ),
            });
        }

        // Resolve the trained project file against the case dir.
        // Round-9 hardening: `project` is user-supplied data and
        // flows into `--project=<path>`; wrap with `confined_join`.
        let resolved_project = if input.project.is_absolute() {
            input.project.clone()
        } else {
            valenx_core::adapter_helpers::confined_join(&case.path, &input.project)?
        };
        if !resolved_project.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.ilastik].project `{}` not found (resolved {})",
                    input.project.display(),
                    resolved_project.display()
                ),
            });
        }

        // Resolve every input image against the case dir and
        // validate it exists. Ilastik will fail fast on its own,
        // but surfacing the missing-file error at validation time
        // is much friendlier than a Python traceback at run time.
        // Round-9 hardening: per-image relative paths flow into the
        // Ilastik command line; wrap with `confined_join`.
        let mut resolved_inputs: Vec<std::path::PathBuf> =
            Vec::with_capacity(input.input_images.len());
        for img in &input.input_images {
            let resolved = if img.is_absolute() {
                img.clone()
            } else {
                valenx_core::adapter_helpers::confined_join(&case.path, img)?
            };
            if !resolved.is_file() {
                return Err(AdapterError::InvalidCase {
                    case_path: case.path.join("case.toml"),
                    reason: format!(
                        "[bio.ilastik].input_images entry `{}` not found \
                         (resolved {})",
                        img.display(),
                        resolved.display()
                    ),
                });
            }
            resolved_inputs.push(resolved);
        }

        // Compose
        // `<ilastik_app> --headless --project=<project>
        //  --output_filename_format=<basename>_{nickname}.h5
        //  <input_images...> [extras...]`.
        // The `{nickname}` literal is Ilastik's per-input
        // substitution placeholder — Ilastik replaces it at run
        // time with each input file's stem so each input gets a
        // unique disambiguated output filename. The literal must
        // reach Ilastik unmodified.
        let output_filename_format = format!(
            "--output_filename_format={}_{{nickname}}.h5",
            input.output_basename
        );
        let mut native_command: Vec<OsString> = vec![
            resolved_app.into_os_string(),
            OsString::from("--headless"),
            OsString::from(format!("--project={}", resolved_project.display())),
            OsString::from(output_filename_format),
        ];
        for img in resolved_inputs {
            native_command.push(img.into_os_string());
        }
        for arg in &input.extra_args {
            native_command.push(OsString::from(arg));
        }

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            // Headless Ilastik runtime scales with image-set count,
            // image dimensions (3D / 4D stacks add up fast), and
            // classifier complexity. Single-plane batches finish in
            // minutes; volumetric high-content screens run for
            // hours. 4 hours is a generous default that covers most
            // workflows without artificially capping batch jobs.
            estimated_runtime: Some(Duration::from_secs(4 * 60 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting Ilastik", |line| {
            let mut hint = subprocess::Hint::default();
            // Ilastik prints per-input progress markers ("Saving
            // results", "Project file saved") and a closing
            // "Done" / "Successfully exported" message. We pattern
            // match conservatively — log formatting has shifted
            // across the 1.3.x → 1.4.x line.
            if line.contains("Successfully exported")
                || line.contains("Done.")
                || line.contains("Project file saved")
            {
                hint.progress = Some((95.0, line.to_string()));
            } else if line.contains("Saving results")
                || line.contains("Loading")
                || line.contains("Predicting")
            {
                hint.progress = Some((50.0, line.to_string()));
            } else if line.contains("Traceback") || line.contains("Error") || line.contains("ERROR")
            {
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
        // Provenance: hash the case.toml as the canonical input
        // descriptor. Ilastik output filenames are
        // `{output_basename}_{nickname}.h5` — driven by the
        // case.toml's `output_basename` field — so the case.toml
        // is the right anchor for run identity.
        let case_hash_input = job.workdir.join("case.toml");
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "ilastik",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);
        let mut artefacts: Vec<Artifact> = Vec::new();

        // Restrict format-specific outputs to those whose stem
        // starts with the configured `output_basename` so any
        // case-staged input image that happens to share an
        // extension doesn't pollute the artefact list. Logs
        // (`*.log`) are accepted regardless of stem since
        // Ilastik / the underlying Python logger may name them
        // on their own.
        let basename = read_output_basename(&job.workdir);

        let entries = match fs::read_dir(&job.workdir) {
            Ok(e) => e,
            Err(e) => {
                tracing::warn!(target: "valenx-ilastik", ?e, "workdir read failed");
                return Ok(results);
            }
        };
        for entry in entries.flatten() {
            let path = entry.path();
            if !path.is_file() {
                continue;
            }
            let ext = path
                .extension()
                .and_then(|s| s.to_str())
                .map(|s| s.to_ascii_lowercase());
            let stem = path
                .file_stem()
                .and_then(|s| s.to_str())
                .unwrap_or("")
                .to_string();
            let stem_matches_basename = match basename.as_deref() {
                Some(b) => stem.starts_with(b),
                None => true,
            };
            let (kind, label) = match ext.as_deref() {
                Some("h5") if stem_matches_basename => (
                    ArtifactKind::Native,
                    "Ilastik probability map (HDF5)".to_string(),
                ),
                Some("tif") if stem_matches_basename => {
                    (ArtifactKind::Native, "Ilastik segmentation".to_string())
                }
                Some("log") => (ArtifactKind::Log, "Ilastik log".to_string()),
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
            ribbon_contributions: vec!["bio.ilastik.classify"],
        }
    }
}

/// Re-read the `[bio.ilastik].output_basename` from a staged
/// `case.toml` for collect()-time output filtering. Returns
/// `None` when the case.toml is missing or unparseable —
/// collect() then accepts every recognised file in the workdir
/// (best-effort).
fn read_output_basename(workdir: &Path) -> Option<String> {
    // Round-23 sweep: bound staged case.toml at MAX_PROJECT_FILE_BYTES.
    let text = valenx_core::io_caps::read_capped_to_string(
        &workdir.join("case.toml"),
        valenx_core::project::loader::MAX_PROJECT_FILE_BYTES as usize,
    )
    .ok()?;
    let parsed: toml::Value = toml::from_str(&text).ok()?;
    parsed
        .get("bio")?
        .get("ilastik")?
        .get("output_basename")?
        .as_str()
        .map(str::to_string)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn info_is_bio_domain() {
        let info = IlastikAdapter::new().info();
        assert_eq!(info.id, "ilastik");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "GPL-3.0");
        assert_eq!(info.display_name, "Ilastik");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = IlastikAdapter::new().info();
        // Ilastik 1.4 (2022) is the modern stable line; 2.0
        // reserves room for the next major bump.
        assert_eq!(info.version_range.min_inclusive, Version::new(1, 4, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(2, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = IlastikAdapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.ilastik.classify"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = IlastikAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }

    /// Round-9 RED→GREEN: `[bio.ilastik].project` used to be joined
    /// with bare `case.path.join`. Wrap with `confined_join`.
    #[test]
    fn prepare_rejects_project_traversing_outside_case_dir() {
        use valenx_test_utils::tempdir;
        let d = tempdir("ilastik-project-trav");
        std::fs::write(d.join("ilastik-fake"), b"").unwrap();
        std::fs::write(d.join("img.tif"), b"x").unwrap();
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "ilastik.classify"

[bio.ilastik]
ilastik_app     = "ilastik-fake"
project         = "../../etc/passwd"
input_images    = ["img.tif"]
output_basename = "out"
"#,
        )
        .unwrap();
        let case = Case {
            id: "ilastik-project-trav".into(),
            path: d.clone(),
        };
        let workdir = d.join("workdir");
        let err = IlastikAdapter::new().prepare(&case, &workdir).unwrap_err();
        let msg = format!("{err}");
        assert!(
            msg.contains("..") || msg.contains("stay within") || msg.contains("escape"),
            "expected confined_join rejection, got: {msg}"
        );
        let _ = std::fs::remove_dir_all(&d);
    }
}
