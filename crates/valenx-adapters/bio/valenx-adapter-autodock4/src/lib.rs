//! # valenx-adapter-autodock4
//!
//! Adapter for [AutoDock 4](https://autodock.scripps.edu/) — Scripps'
//! original genetic-algorithm small-molecule docking engine. Predates
//! Vina by a decade and is still the reference implementation in the
//! published literature: explicit force field, fully scriptable
//! docking parameter file, the canonical pose-clustering + binding-
//! energy report. Two binaries cooperate:
//!
//! 1. `autogrid4` — pre-computes the per-atom-type affinity grids
//!    around the receptor (one `.map` per ligand atom type, plus a
//!    `.fld` describing the grid geometry). Slow, but deterministic;
//!    if you dock many ligands against the same pocket the grids are
//!    re-used.
//! 2. `autodock4` — runs the Lamarckian genetic-algorithm search
//!    against those grids and writes a `.dlg` log with embedded
//!    poses and the clustering summary.
//!
//! **Phase 34 — two-stage subprocess wrapper.** Mirrors the BWA
//! adapter's `index → mem` pattern:
//!
//! * `prepare()` runs `autogrid4` synchronously (skippable via
//!   `skip_grid = true` if maps already exist), then composes the
//!   `autodock4 -p <dpf> -l <dock_log>` invocation as the
//!   `PreparedJob`. Doing the grid step inline keeps the subsequent
//!   `run()` call as a single subprocess so the shared runner can
//!   stream `autodock4`'s stderr line-by-line.
//! * `run()` executes only the `autodock4` step.
//! * `collect()` walks the workdir for the `.dlg` log + any `.pdbqt`
//!   poses, surfacing them as `Log` and `Native` artifacts.
//!
//! License: AutoDock 4 ships under GPL-2.0-or-later; we run it as a
//! subprocess (`LicenseMode::Subprocess`) so no source is linked into
//! the Valenx binary.

#![forbid(unsafe_code)]
#![allow(missing_docs)]

pub mod case_input;

use std::ffi::OsString;
use std::path::Path;
use std::time::Duration;

use semver::Version;

use valenx_core::{
    adapter_helpers::{detect_tool_version_semver, find_on_path, live_provenance},
    error::RunPhase,
    subprocess, Adapter, AdapterError, AdapterInfo, Capabilities, Case, LicenseMode, Physics,
    PreparedJob, ProbeReport, RunContext, RunReport, VersionRange,
};
use valenx_fields::{
    artifact::{Artifact, ArtifactKind},
    Results,
};

use crate::case_input::AutoDock4Input;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(AutoDock4Adapter::new())
}

pub struct AutoDock4Adapter;

impl AutoDock4Adapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for AutoDock4Adapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "autodock4";
const DOCK_BINARIES: &[&str] = &["autodock4"];
const GRID_BINARIES: &[&str] = &["autogrid4"];

/// Probe-time warning emitted when `autodock4` is on PATH but the
/// companion `autogrid4` binary is not. We do *not* fail probe — the
/// user might legitimately rely on pre-computed grids and
/// `skip_grid = true`. Pinned as a constant so a unit test can assert
/// on the literal substring without re-stating it.
pub const AUTOGRID_MISSING_WARNING: &str =
    "autogrid4 not found on PATH; AutoDock 4 docking requires both \
     autogrid4 + autodock4. Skip the grid step via skip_grid=true \
     only if grid maps already exist.";

impl Adapter for AutoDock4Adapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "AutoDock 4",
            // 4.2.6 (2014) is the current frozen version; 4.2.x is the
            // long-stable line that every distro and the Scripps mirror
            // ship. 5.0 reserves room for an eventual major bump
            // (the AutoDock-GPU effort lives at a different binary).
            version_range: VersionRange {
                min_inclusive: Version::new(4, 2, 0),
                max_exclusive: Version::new(5, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            tool_license: "GPL-2.0-or-later",
            docs_url:
                "https://autodock.scripps.edu/wp-content/uploads/sites/56/2021/10/AutoDock4.2.6_UserGuide.pdf",
            homepage_url: "https://autodock.scripps.edu/",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(DOCK_BINARIES) {
            Some(binary_path) => {
                // `autodock4 --version` (or no args) prints a banner
                // including the version string to stdout / stderr.
                let found_version = detect_tool_version_semver(&binary_path, &["--version", ""]);
                let mut warnings: Vec<String> = Vec::new();
                // The grid step is mandatory unless the user explicitly
                // re-uses pre-computed maps. If autogrid4 is missing
                // we still report ok=true (the user might be running a
                // skip_grid case) but surface a clear warning.
                if find_on_path(GRID_BINARIES).is_none() {
                    warnings.push(AUTOGRID_MISSING_WARNING.to_string());
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
                hint: "AutoDock 4.2+ required; install via \
                       `apt install autodock`, \
                       `conda install -c bioconda autodock`, or build from \
                       https://autodock.scripps.edu/"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = AutoDock4Input::from_case_dir(&case.path)?;

        std::fs::create_dir_all(workdir)?;

        // Stage `case.toml` into the workdir so collect() can recover
        // the configured `output_basename` for prefix-filtering output
        // artifacts. Without this stage, the basename filter silently
        // degrades to "match everything".
        let staged_case_toml = workdir.join("case.toml");
        let source_case_toml = case.path.join("case.toml");
        if source_case_toml.is_file() {
            std::fs::copy(&source_case_toml, &staged_case_toml)
                .map_err(|e| AdapterError::Other(anyhow::anyhow!("stage case.toml: {e}")))?;
        }

        // Resolve the GPF + DPF parameter files relative to the case
        // directory (the convention every other bio adapter follows
        // for case-relative paths).
        let resolved_gpf = if input.gpf.is_absolute() {
            input.gpf.clone()
        } else {
            valenx_core::adapter_helpers::confined_join(
            &case.path,
            &input.gpf,
        )?
        };
        let resolved_dpf = if input.dpf.is_absolute() {
            input.dpf.clone()
        } else {
            valenx_core::adapter_helpers::confined_join(
            &case.path,
            &input.dpf,
        )?
        };
        if !resolved_gpf.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.autodock4].gpf `{}` not found (resolved {})",
                    input.gpf.display(),
                    resolved_gpf.display()
                ),
            });
        }
        if !resolved_dpf.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.autodock4].dpf `{}` not found (resolved {})",
                    input.dpf.display(),
                    resolved_dpf.display()
                ),
            });
        }

        let dock_binary =
            find_on_path(DOCK_BINARIES).ok_or_else(|| AdapterError::ToolNotInstalled {
                name: INFO_ID,
                hint: "AutoDock 4.2+ required; install via \
                       `apt install autodock`, \
                       `conda install -c bioconda autodock`, or build from \
                       https://autodock.scripps.edu/"
                    .into(),
            })?;

        // Stage 1: autogrid4 -p <gpf> -l <grid_log> [extras...].
        // Mirrors the BWA `bwa index <reference>` pre-stage at
        // valenx-adapter-bwa/src/lib.rs:178-201 — synchronous run in
        // prepare() so the subsequent run() call has only one
        // subprocess to stream from.
        if !input.skip_grid {
            let grid_binary =
                find_on_path(GRID_BINARIES).ok_or_else(|| AdapterError::ToolNotInstalled {
                    name: "autogrid4",
                    hint: "autogrid4 (sister tool to autodock4) is required \
                           unless skip_grid=true and the grid maps already exist."
                        .into(),
                })?;
            let mut cmd = std::process::Command::new(&grid_binary);
            cmd.arg("-p")
                .arg(&resolved_gpf)
                .arg("-l")
                .arg(&input.grid_log);
            for extra in &input.extra_grid_args {
                cmd.arg(extra);
            }
            cmd.current_dir(workdir)
                .stdin(std::process::Stdio::null())
                .stdout(std::process::Stdio::piped())
                .stderr(std::process::Stdio::piped());
            match cmd.output() {
                Ok(out) if out.status.success() => {}
                Ok(out) => {
                    let stderr = String::from_utf8_lossy(&out.stderr);
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "`autogrid4 -p {}` failed (exit {}): {}",
                        resolved_gpf.display(),
                        out.status.code().unwrap_or(-1),
                        stderr.lines().next().unwrap_or("(no stderr)")
                    )));
                }
                Err(e) => {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "spawning `autogrid4 -p {}` failed: {e}",
                        resolved_gpf.display()
                    )));
                }
            }
        }

        // Stage 2: build the autodock4 command. The DPF references its
        // grid maps by relative path, so the working directory matters;
        // the runner sets `current_dir` to `workdir` for us.
        let mut native_command: Vec<OsString> = vec![
            dock_binary.into_os_string(),
            OsString::from("-p"),
            resolved_dpf.into_os_string(),
            OsString::from("-l"),
            OsString::from(&input.dock_log),
        ];
        for arg in &input.extra_dock_args {
            native_command.push(OsString::from(arg));
        }

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            // GA docking with default params (50 runs, ~25k evaluations
            // each) takes 10–60 minutes on a single CPU; 4 hours is
            // generous headroom for high-precision parameter sets.
            estimated_runtime: Some(Duration::from_secs(4 * 60 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting AutoDock 4", |line| {
            let mut hint = subprocess::Hint::default();
            // autodock4 prints "Run:   N / M" lines roughly per GA run
            // (default M=50). Lift those to mid-progress markers and
            // pin "Successful Completion" / "Real time:" at 95%.
            if line.contains("Successful Completion") || line.contains("Real time:") {
                hint.progress = Some((95.0, line.to_string()));
            } else if line.contains("Run:") || line.contains("DOCKED:") {
                hint.progress = Some((50.0, line.to_string()));
            } else if line.contains("ERROR") || line.contains("FATAL") {
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
        // Re-read the case so we know what dock_log filename to look
        // for. Falling back to the canonical default when re-parse
        // fails keeps collect() resilient against post-run edits.
        let dock_log_name = AutoDock4Input::from_case_dir(&job.workdir)
            .map(|i| i.dock_log)
            .unwrap_or_else(|_| "autodock4.dlg".to_string());

        // Provenance: hash the dock log when present, or fall back to
        // case.toml so the provenance block is well-formed even on
        // partial / failed runs.
        let case_hash_input = {
            let dlg = job.workdir.join(&dock_log_name);
            if dlg.is_file() {
                dlg
            } else {
                job.workdir.join("case.toml")
            }
        };
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "AutoDock 4",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);
        let mut artefacts: Vec<Artifact> = Vec::new();

        // Walk the workdir top-level. We surface:
        //   - the configured dock_log (Log)
        //   - any `.dlg` (Native; the docked-poses + clustering log)
        //   - any `.pdbqt` (Native; some workflows emit ranked poses)
        let entries = match std::fs::read_dir(&job.workdir) {
            Ok(e) => e,
            Err(e) => {
                tracing::warn!(target: "valenx-autodock4", ?e, "workdir read failed");
                return Ok(results);
            }
        };
        for entry in entries.flatten() {
            let path = entry.path();
            if !path.is_file() {
                continue;
            }
            let file_name = path
                .file_name()
                .and_then(|s| s.to_str())
                .unwrap_or_default()
                .to_string();
            let ext = path
                .extension()
                .and_then(|s| s.to_str())
                .map(|s| s.to_ascii_lowercase());
            let (kind, label) = if file_name == dock_log_name {
                (ArtifactKind::Log, "AutoDock 4 docking log".to_string())
            } else {
                match ext.as_deref() {
                    // `.dlg` outputs that aren't the configured dock_log
                    // name still count — some users rename them.
                    Some("dlg") => (ArtifactKind::Native, "AutoDock 4 docked poses".to_string()),
                    Some("pdbqt") => (ArtifactKind::Native, "AutoDock 4 docked poses".to_string()),
                    _ => continue,
                }
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
            ribbon_contributions: vec!["bio.autodock4.dock"],
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn info_is_bio_domain() {
        let info = AutoDock4Adapter::new().info();
        assert_eq!(info.id, "autodock4");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "GPL-2.0-or-later");
        assert_eq!(info.display_name, "AutoDock 4");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = AutoDock4Adapter::new().info();
        // AutoDock 4.2.x is the long-stable line; 5.0 reserves room
        // for an eventual major bump.
        assert_eq!(info.version_range.min_inclusive, Version::new(4, 2, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(5, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = AutoDock4Adapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.autodock4.dock"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = AutoDock4Adapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }

    #[test]
    fn probe_warning_mentions_autogrid4_missing() {
        // Lock the warning text so a future refactor that drops the
        // `autogrid4 not found` substring trips a test failure — the
        // probe ribbon UI greps for this exact phrase to render a
        // missing-tool hint.
        assert!(
            AUTOGRID_MISSING_WARNING.contains("autogrid4 not found"),
            "warning constant must mention `autogrid4 not found`: {AUTOGRID_MISSING_WARNING}"
        );
    }
}
