//! # valenx-adapter-freecad
//!
//! Adapter that drives FreeCAD's headless `FreeCADCmd` as a
//! subprocess to import STEP / IGES / STL / BREP / FCStd documents,
//! bake shapes, emit a feature-tree summary, and export into the
//! canonical geometry formats downstream adapters consume.
//!
//! **Phase 2 — live for the import path.** `prepare()` generates a
//! deterministic Python script + command line; `run()` spawns
//! FreeCADCmd via the shared [`valenx_core::subprocess`] runner;
//! `collect()` parses the `summary.json` the script wrote back and
//! attaches every produced artifact to `Results`.
//!
//! Parametric rebuilds (edit a pattern count, re-execute the
//! feature tree in under 500 ms) land once the adapter has a real
//! parameter-override surface.

#![forbid(unsafe_code)]
#![allow(missing_docs)]

pub mod case_input;
pub mod python_script;
pub mod summary_parser;

use std::ffi::OsString;
use std::fs;
use std::path::Path;
use std::time::Duration;

use semver::Version;

use valenx_core::{
    adapter_helpers::{find_on_path, first_workdir_match},
    subprocess, Adapter, AdapterError, AdapterInfo, Capabilities, Capability, Case, LicenseMode,
    Physics, PreparedJob, ProbeReport, RunContext, RunReport, VersionRange,
};
use valenx_fields::{
    artifact::{Artifact, ArtifactKind},
    Results,
};

use crate::case_input::GeometryImportInput;
use crate::python_script::{SCRIPT_FILENAME, SUMMARY_FILENAME};

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(FreeCadAdapter::new())
}

pub struct FreeCadAdapter;

impl FreeCadAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for FreeCadAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "freecad";
const BINARIES: &[&str] = &["FreeCADCmd", "freecadcmd", "FreeCAD"];

impl Adapter for FreeCadAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "FreeCAD",
            version_range: VersionRange {
                min_inclusive: Version::new(0, 21, 0),
                max_exclusive: Version::new(2, 0, 0),
            },
            physics: &[Physics::Geometry],
            license_mode: LicenseMode::Subprocess,
            tool_license: "LGPL-2.1-or-later",
            docs_url: "https://wiki.freecad.org/",
            homepage_url: "https://www.freecad.org/",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(BINARIES) {
            Some(binary_path) => {
                // FreeCAD's CLI is `FreeCADCmd --version`. The
                // helper tries each candidate flag and returns the
                // first parseable semver — None if FreeCAD is a
                // wrapper script that ate the version output.
                let found_version = valenx_core::adapter_helpers::detect_tool_version_semver(
                    &binary_path,
                    &["--version", "-v"],
                );
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
                hint: "FreeCAD 0.21+ required; install from freecad.org or your \
                       distribution (FreeCADCmd / freecadcmd on PATH)"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let (_header, input) = GeometryImportInput::from_case_dir(&case.path)?;

        fs::create_dir_all(workdir)?;

        // Source file resolution — we accept an absolute path, or a
        // path relative to the case dir.
        let source_abs = if input.source.is_absolute() {
            input.source.clone()
        } else {
            valenx_core::adapter_helpers::confined_join(
            &case.path,
            &input.source,
        )?
        };
        if !source_abs.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[geometry] source {} not found (resolved {})",
                    input.source.display(),
                    source_abs.display()
                ),
            });
        }

        // Stage the source alongside the script — FreeCADCmd's cwd is
        // the workdir, and relative paths are fragile across versions.
        let source_in_wd = workdir.join(
            source_abs
                .file_name()
                .expect("source must have a filename")
                .to_string_lossy()
                .to_string(),
        );
        if source_abs != source_in_wd {
            fs::copy(&source_abs, &source_in_wd)?;
        }

        // Generate + write the Python script.
        let script_path = workdir.join(SCRIPT_FILENAME);
        python_script::write_to_file(&input, &script_path)?;

        let binary_path = find_on_path(BINARIES).ok_or_else(|| AdapterError::ToolNotInstalled {
            name: INFO_ID,
            hint: "FreeCADCmd / freecadcmd not found on PATH; install \
                       FreeCAD 0.21+"
                .into(),
        })?;

        // Command: FreeCADCmd <script> -- <source>
        // FreeCADCmd passes everything after `--` through argv, which
        // is where the script reads `sys.argv[-1]`.
        let native_command: Vec<OsString> = vec![
            binary_path.into_os_string(),
            OsString::from(SCRIPT_FILENAME),
            OsString::from("--"),
            OsString::from(
                source_in_wd
                    .file_name()
                    .expect("staged file must have a name"),
            ),
        ];

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            // Imports are usually fast — under 30 s even for
            // moderately complex STEP files. Upper bound the estimate
            // generously; the UI's progress bar prefers over-
            // estimating to under-estimating.
            estimated_runtime: Some(Duration::from_secs(60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting FreeCAD", |line| {
            let mut hint = subprocess::Hint::default();
            if let Some(pct) = freecad_progress_hint(line) {
                hint.progress = Some((pct, line.to_string()));
            }
            if line.contains("<Error>") || line.contains("Traceback") {
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
            final_phase: Some(valenx_core::error::RunPhase::Shutdown),
        })
    }

    fn collect(&self, job: &PreparedJob) -> Result<Results, AdapterError> {
        // Real provenance: hash whichever source CAD format was
        // imported (.step / .iges / .brep / .fcstd / .stl) plus
        // the generated FreeCAD Python script. The first match
        // becomes the canonical case_hash; mesh_path stays None
        // since FreeCAD is upstream of any meshing.
        let case_path = first_workdir_match(
            &job.workdir,
            &["step", "stp", "iges", "igs", "brep", "fcstd", "stl"],
        )
        .or_else(|| {
            let py = job.workdir.join(SCRIPT_FILENAME);
            py.exists().then_some(py)
        })
        .unwrap_or_else(|| job.workdir.join("(no-input-found)"));
        let prov = valenx_core::adapter_helpers::live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "FreeCAD",
            "unknown",
            &case_path,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);

        // Parse the summary the script wrote, if present, and lift
        // its headline figures into a `description` on the Results
        // meta. The raw summary.json also goes on the artifact list.
        let summary_path = job.workdir.join(SUMMARY_FILENAME);
        if summary_path.is_file() {
            match summary_parser::parse_file(&summary_path) {
                Ok(summary) => {
                    let desc = format!(
                        "FreeCAD import · {} part(s) · vol {:.3e} · area {:.3e}",
                        summary.parts.len(),
                        summary.volume.unwrap_or(0.0),
                        summary.area.unwrap_or(0.0),
                    );
                    results.meta.description = Some(desc);
                    results.artifacts.push(Artifact {
                        path: summary_path.clone(),
                        kind: ArtifactKind::Other,
                        checksum: None,
                        label: format!("FreeCAD summary · {} parts", summary.parts.len()),
                    });
                }
                Err(e) => {
                    tracing::warn!(target: "valenx-freecad", ?e, "summary parse failed");
                    results.artifacts.push(Artifact {
                        path: summary_path,
                        kind: ArtifactKind::Other,
                        checksum: None,
                        label: format!("FreeCAD summary (parse error: {e})"),
                    });
                }
            }
        }

        // Exported geometry files.
        for (name, kind, label) in [
            ("output.stl", ArtifactKind::VizData, "FreeCAD STL export"),
            ("output.brep", ArtifactKind::Native, "FreeCAD BREP export"),
            ("output.step", ArtifactKind::Native, "FreeCAD STEP export"),
            ("output.iges", ArtifactKind::Native, "FreeCAD IGES export"),
        ] {
            let path = job.workdir.join(name);
            if path.is_file() {
                results.artifacts.push(Artifact {
                    path,
                    kind,
                    checksum: None,
                    label: label.to_string(),
                });
            }
        }

        // The generated script itself — keep it on the artifact list
        // for transparency.
        let script_path = job.workdir.join(SCRIPT_FILENAME);
        if script_path.is_file() {
            results.artifacts.push(Artifact {
                path: script_path,
                kind: ArtifactKind::Other,
                checksum: None,
                label: "FreeCAD script (generated)".into(),
            });
        }

        results.artifacts.sort_by(|a, b| a.path.cmp(&b.path));
        Ok(results)
    }

    fn capabilities(&self) -> Capabilities {
        Capabilities {
            capabilities: vec![
                Capability::GeoStep,
                Capability::GeoIges,
                Capability::GeoBRep,
                Capability::GeoStl,
                Capability::GeoSketch,
            ],
            ribbon_contributions: vec![
                "cad.freecad.import",
                "cad.freecad.export",
                "cad.freecad.parameters",
            ],
        }
    }
}

/// Coarse progress hints derived from the messages our generated
/// script prints + FreeCAD's own startup banners.
fn freecad_progress_hint(line: &str) -> Option<f32> {
    if line.contains("[valenx] opening") {
        Some(15.0)
    } else if line.contains("Open document") || line.contains("Import:") {
        Some(35.0)
    } else if line.contains("Mesh.export") || line.contains("exportBrep") {
        Some(65.0)
    } else if line.contains(&format!("wrote {SUMMARY_FILENAME}")) {
        Some(95.0)
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn info_declares_subprocess_mode() {
        let info = FreeCadAdapter::new().info();
        assert_eq!(info.id, "freecad");
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
        assert_eq!(info.physics, &[Physics::Geometry]);
    }

    #[test]
    fn progress_hints_are_monotonic() {
        let pts = [
            freecad_progress_hint("[valenx] opening bracket.step"),
            freecad_progress_hint("Open document bracket"),
            freecad_progress_hint("Mesh.export finished"),
            freecad_progress_hint(&format!("[valenx] wrote {SUMMARY_FILENAME}")),
        ];
        let mut last = 0.0f32;
        for (i, p) in pts.iter().enumerate() {
            let v = p.expect("known token");
            assert!(v >= last, "step {i}: {last} -> {v}");
            last = v;
        }
    }
}
