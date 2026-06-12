//! # valenx-adapter-cellprofiler
//!
//! Adapter for [CellProfiler](https://cellprofiler.org/) — the Broad
//! Institute's pipeline-driven cell-segmentation + measurement suite
//! (Carpenter / Lamprecht / Sabatini et al). CellProfiler pipelines
//! are authored in the desktop GUI and saved as `.cppipe` (text) or
//! `.cpproj` (binary project) files; in production they're typically
//! re-run headlessly across new image batches via the `cellprofiler`
//! CLI's `-c` (no-GUI) + `-r` (run-immediately) flags.
//!
//! **Phase 40 — sister adapter to Fiji and Ilastik for the
//! microscopy / bioimage-analysis surface.** The adapter composes
//!
//! ```text
//! cellprofiler -c -r -p <pipeline> -i <input_dir> -o <output_basename> [extras...]
//! ```
//!
//! and falls back to `<python> -m cellprofiler ...` when the
//! standalone `cellprofiler` shim isn't on PATH but the package was
//! installed into a Python environment (typical conda layout).
//!
//! On `collect()` we walk **one level deep** into the
//! `<output_basename>/` subdir for the canonical CellProfiler
//! outputs (`.csv` measurement tables, `.tif` / `.tiff` segmented
//! images, `.png` plots) and surface any `*.log` files at the top
//! level of the workdir as logs.

#![forbid(unsafe_code)]
#![allow(missing_docs)]

pub mod case_input;

use std::ffi::OsString;
use std::fs;
use std::path::{Path, PathBuf};
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

use crate::case_input::CellProfilerInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(CellProfilerAdapter::new())
}

pub struct CellProfilerAdapter;

impl CellProfilerAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for CellProfilerAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "cellprofiler";
/// Probe binary candidates. `cellprofiler` first because that's the
/// canonical install-name shim; `python3` / `python` cover the
/// fallback path (`python -m cellprofiler ...`) when only the package
/// is reachable.
const PROBE_BINARIES: &[&str] = &["cellprofiler", "python3", "python"];
/// Python interpreter candidates for the `python -m cellprofiler`
/// fallback.
const PYTHON_BINARIES: &[&str] = &["python3", "python"];

impl Adapter for CellProfilerAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "CellProfiler",
            // CellProfiler 4.0 (2020) is the modern stable line —
            // first to ship the unified DeepProfiler + skimage-based
            // image-processing modules; 5.0 reserves the next major.
            version_range: VersionRange {
                min_inclusive: Version::new(4, 0, 0),
                max_exclusive: Version::new(5, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            tool_license: "BSD-3-Clause",
            docs_url: "https://cellprofiler.org/manuals",
            homepage_url: "https://cellprofiler.org/",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        // Prefer the standalone `cellprofiler` shim; fall back to
        // Python so the `python -m cellprofiler` path still works.
        // The combined PATH walk matches the spec convention; the
        // per-binary lookup below isolates the "Python is here but
        // CellProfiler isn't" warning condition.
        match find_on_path(PROBE_BINARIES) {
            Some(binary_path) => {
                let cellprofiler_present = find_on_path(&["cellprofiler"]).is_some();
                let mut warnings = Vec::new();
                if !cellprofiler_present {
                    // Python is on PATH (some hit) but the dedicated
                    // CLI shim isn't — the `python -m cellprofiler`
                    // fallback will kick in at run time. Surface a
                    // warning so the user knows why their PATH lookup
                    // didn't find the binary.
                    warnings.push(
                        "CellProfiler not found on PATH; install via \
                         `pip install cellprofiler` or download from \
                         https://cellprofiler.org/releases"
                            .into(),
                    );
                }
                Ok(ProbeReport {
                    ok: true,
                    found_version: None,
                    binary_path: Some(binary_path),
                    warnings,
                    required_env: Vec::new(),
                })
            }
            None => Err(AdapterError::ToolNotInstalled {
                name: INFO_ID,
                hint: "CellProfiler 4.0+ required; install via \
                       `pip install cellprofiler` or download from \
                       https://cellprofiler.org/releases and ensure \
                       `cellprofiler` (or `python3`) is on PATH"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = CellProfilerInput::from_case_dir(&case.path)?;

        // Round-4 security: reject `output_basename = "../etc/passwd"`
        // and friends before the value flows into any path join.
        // Same pattern as the round-3 fix in bionetgen/iqtree/art/fasttree.
        valenx_core::adapter_helpers::validate_output_basename(
            &input.output_basename,
            "[bio.cellprofiler].output_basename",
        )
        .map_err(|e| AdapterError::InvalidCase {
            case_path: case.path.join("case.toml"),
            reason: format!("{e}"),
        })?;

        fs::create_dir_all(workdir)?;

        // Resolve the pipeline path against the case directory if
        // relative. CellProfiler reads the pipeline by absolute path
        // — we don't stage it (some `.cpproj` projects are large
        // binaries and re-running with the original path keeps the
        // pipeline-author's relative-path references intact).
        let resolved_pipeline = if input.pipeline.is_absolute() {
            input.pipeline.clone()
        } else {
            valenx_core::adapter_helpers::confined_join(&case.path, &input.pipeline)?
        };
        if !resolved_pipeline.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.cellprofiler].pipeline `{}` not found (resolved {})",
                    input.pipeline.display(),
                    resolved_pipeline.display()
                ),
            });
        }

        // Resolve the input image directory against the case dir.
        // CellProfiler's `Images` module walks this directory at run
        // time, so we just need it to exist as a directory (the
        // pipeline itself decides which files to consume).
        // Round-9 hardening: `input_dir` is user-supplied data and
        // flows into CellProfiler's `-i` flag where it walks the dir;
        // wrap the relative branch with `confined_join` so a hostile
        // case can't aim it at `..`-traversal targets.
        let resolved_input_dir = if input.input_dir.is_absolute() {
            input.input_dir.clone()
        } else {
            valenx_core::adapter_helpers::confined_join(&case.path, &input.input_dir)?
        };
        if !resolved_input_dir.is_dir() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.cellprofiler].input_dir `{}` not found or not a directory \
                     (resolved {})",
                    input.input_dir.display(),
                    resolved_input_dir.display()
                ),
            });
        }

        // Output directory is workdir-relative — CellProfiler
        // creates `<workdir>/<output_basename>` and writes there.
        // Pass the bare basename; the subprocess runner's cwd is
        // the workdir, so CellProfiler resolves it correctly.
        let output_dir_arg = OsString::from(&input.output_basename);

        // Prefer the standalone `cellprofiler` binary; fall back to
        // the `python -m cellprofiler` invocation so the adapter
        // still works when only the package is installed (typical
        // conda layout). Same fallback shape as OmegaFold.
        let mut native_command: Vec<OsString> = if let Some(bin) = find_on_path(&["cellprofiler"]) {
            vec![
                bin.into_os_string(),
                OsString::from("-c"),
                OsString::from("-r"),
                OsString::from("-p"),
                resolved_pipeline.into_os_string(),
                OsString::from("-i"),
                resolved_input_dir.into_os_string(),
                OsString::from("-o"),
                output_dir_arg,
            ]
        } else {
            // Resolve the Python binary. Same logic as every other
            // Phase 17 / 5.7 Python-script adapter: bare `python` /
            // `python3` walks PATH; absolute paths or pinned
            // interpreters are honored verbatim if present, with a
            // final PATH fallback.
            // Round-4 security: validate python interpreter spec
            // against the allow-list AND resolve to a real binary
            // in one step. Closes the arbitrary-binary-exec class
            // that round-3 only patched in 8 of the 48 affected
            // adapters.
            let binary_path = valenx_core::adapter_helpers::resolve_python_binary(
                &input.python,
                PYTHON_BINARIES,
            )
            .map_err(|e| AdapterError::ToolNotInstalled {
                name: INFO_ID,
                hint: format!("neither `cellprofiler` nor Python interpreter the configured value reachable on PATH: {e}"),
            })?;
            vec![
                binary_path.into_os_string(),
                OsString::from("-m"),
                OsString::from("cellprofiler"),
                OsString::from("-c"),
                OsString::from("-r"),
                OsString::from("-p"),
                resolved_pipeline.into_os_string(),
                OsString::from("-i"),
                resolved_input_dir.into_os_string(),
                OsString::from("-o"),
                output_dir_arg,
            ]
        };

        for arg in &input.extra_args {
            native_command.push(OsString::from(arg));
        }

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            // CellProfiler runtime scales with image-set count and
            // pipeline complexity — single-plate batches finish in
            // minutes, multi-plate high-content screens run for
            // hours. 4 hours is a generous default that covers most
            // workflows without artificially capping batch jobs.
            estimated_runtime: Some(Duration::from_secs(4 * 60 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting CellProfiler", |line| {
            let mut hint = subprocess::Hint::default();
            // CellProfiler prints per-image-set progress markers
            // (`Image # of M:` style) and a closing "Complete" /
            // "All cycles processed" message. We pattern-match
            // conservatively — log formatting has shifted across
            // CellProfiler 4.x point releases.
            if line.contains("All cycles processed")
                || line.contains("Pipeline completed")
                || line.contains("Complete")
            {
                hint.progress = Some((95.0, line.to_string()));
            } else if line.contains("Image #") || line.contains("Processing") {
                hint.progress = Some((50.0, line.to_string()));
            } else if line.contains("Traceback")
                || line.contains("Error")
                || line.contains("ModuleError")
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
        // descriptor. CellProfiler output filenames are pipeline-
        // defined and a single run produces dozens of per-object
        // CSVs + per-image masks; the case.toml captures the
        // pipeline path + basename that drove the run.
        let case_hash_input = job.workdir.join("case.toml");
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "CellProfiler",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);
        let mut artefacts: Vec<Artifact> = Vec::new();

        // Top-level: Python's `logging` output drops `*.log` files
        // into the workdir cwd. We surface those regardless of name
        // (CellProfiler doesn't enforce a fixed log basename).
        let entries = match fs::read_dir(&job.workdir) {
            Ok(e) => e,
            Err(e) => {
                tracing::warn!(target: "valenx-cellprofiler", ?e, "workdir read failed");
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
            if let Some("log") = ext.as_deref() {
                artefacts.push(Artifact {
                    path,
                    kind: ArtifactKind::Log,
                    checksum: None,
                    label: "CellProfiler log".to_string(),
                });
            }
        }

        // One level deep: walk into the configured `output_basename`
        // subdir for the canonical CellProfiler outputs. We don't
        // know the user's basename without reparsing case.toml, so
        // we read it back and walk that single subdir (rather than
        // every subdir in the workdir, which would scoop up any
        // pipeline-staged scratch directories).
        let basename = read_output_basename(&job.workdir);
        if let Some(b) = basename {
            let output_dir = job.workdir.join(&b);
            if output_dir.is_dir() {
                let inner = match fs::read_dir(&output_dir) {
                    Ok(e) => e,
                    Err(e) => {
                        tracing::warn!(
                            target: "valenx-cellprofiler",
                            ?e,
                            output_dir = %output_dir.display(),
                            "output dir read failed"
                        );
                        artefacts.sort_by(|a, b| a.path.cmp(&b.path));
                        results.artifacts = artefacts;
                        return Ok(results);
                    }
                };
                let mut csv_paths: Vec<PathBuf> = Vec::new();
                let mut tif_paths: Vec<PathBuf> = Vec::new();
                let mut png_paths: Vec<PathBuf> = Vec::new();
                for entry in inner.flatten() {
                    let path = entry.path();
                    if !path.is_file() {
                        continue;
                    }
                    let ext = path
                        .extension()
                        .and_then(|s| s.to_str())
                        .map(|s| s.to_ascii_lowercase());
                    match ext.as_deref() {
                        Some("csv") => csv_paths.push(path),
                        Some("tif") | Some("tiff") => tif_paths.push(path),
                        Some("png") => png_paths.push(path),
                        _ => continue,
                    }
                }
                csv_paths.sort();
                tif_paths.sort();
                png_paths.sort();
                for path in csv_paths {
                    artefacts.push(Artifact {
                        path,
                        kind: ArtifactKind::Tabular,
                        checksum: None,
                        label: "CellProfiler measurements".to_string(),
                    });
                }
                for path in tif_paths {
                    artefacts.push(Artifact {
                        path,
                        kind: ArtifactKind::Native,
                        checksum: None,
                        label: "CellProfiler segmented image".to_string(),
                    });
                }
                for path in png_paths {
                    artefacts.push(Artifact {
                        path,
                        kind: ArtifactKind::Native,
                        checksum: None,
                        label: "CellProfiler plot".to_string(),
                    });
                }
            }
        }

        artefacts.sort_by(|a, b| a.path.cmp(&b.path));
        results.artifacts = artefacts;
        Ok(results)
    }

    fn capabilities(&self) -> Capabilities {
        Capabilities {
            capabilities: Vec::new(),
            ribbon_contributions: vec!["bio.cellprofiler.segment"],
        }
    }
}

/// Re-read the `[bio.cellprofiler].output_basename` from a staged
/// `case.toml` for collect()-time output filtering. Returns `None`
/// when the case.toml is missing or unparseable — collect() then
/// returns just the top-level log artefacts.
fn read_output_basename(workdir: &Path) -> Option<String> {
    // Round-23 sweep: bound the staged case.toml read at
    // MAX_PROJECT_FILE_BYTES (1 MiB) — sister to the round-21 L3
    // adapter-params sweep. A poisoned workdir with a multi-GB
    // case.toml would slurp into memory before toml::from_str saw
    // it. Returns None on cap-failure so collect() degrades gracefully.
    let text = valenx_core::io_caps::read_capped_to_string(
        &workdir.join("case.toml"),
        valenx_core::project::loader::MAX_PROJECT_FILE_BYTES as usize,
    )
    .ok()?;
    let parsed: toml::Value = toml::from_str(&text).ok()?;
    parsed
        .get("bio")?
        .get("cellprofiler")?
        .get("output_basename")?
        .as_str()
        .map(str::to_string)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn info_is_bio_domain() {
        let info = CellProfilerAdapter::new().info();
        assert_eq!(info.id, "cellprofiler");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "BSD-3-Clause");
        assert_eq!(info.display_name, "CellProfiler");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = CellProfilerAdapter::new().info();
        // CellProfiler 4.0 (2020) shipped the modern unified
        // skimage-based image-processing modules; 5.0 reserves
        // room for the next major bump.
        assert_eq!(info.version_range.min_inclusive, Version::new(4, 0, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(5, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = CellProfilerAdapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.cellprofiler.segment"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = CellProfilerAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }

    /// Round-9 RED→GREEN: `[bio.cellprofiler].input_dir` used to be
    /// joined with bare `case.path.join`, letting a hostile case
    /// supply `input_dir = "../../etc"` and have CellProfiler walk
    /// arbitrary directories. The fix wraps the relative branch with
    /// `confined_join`.
    #[test]
    fn prepare_rejects_input_dir_traversing_outside_case_dir() {
        use valenx_test_utils::tempdir;
        let d = tempdir("cellprofiler-input-trav");
        std::fs::write(d.join("pipe.cppipe"), b"# placeholder\n").unwrap();
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "cellprofiler.segment"

[bio.cellprofiler]
pipeline        = "pipe.cppipe"
input_dir       = "../../etc"
output_basename = "out"
"#,
        )
        .unwrap();
        let case = Case {
            id: "cellprofiler-input-trav".into(),
            path: d.clone(),
        };
        let workdir = d.join("workdir");
        let err = CellProfilerAdapter::new()
            .prepare(&case, &workdir)
            .unwrap_err();
        let msg = format!("{err}");
        assert!(
            msg.contains("..") || msg.contains("stay within") || msg.contains("escape"),
            "expected confined_join rejection, got: {msg}"
        );
        let _ = std::fs::remove_dir_all(&d);
    }
}
