//! # valenx-adapter-planemo
//!
//! Adapter for [Planemo](https://planemo.readthedocs.io/) — the Galaxy
//! project's official command-line companion for tool development and
//! workflow execution outside a full Galaxy server. The same binary
//! lints tool wrappers, runs Galaxy workflow tests, and executes
//! `.ga` / `.gxwf.yml` workflows; the `action` knob picks which
//! subcommand to invoke.
//!
//! **Phase 22.5 — workflow-CLI wrapper.** Sister adapter to
//! Snakemake's: composes `planemo <action> <workflow> [inputs]
//! [extras...]` from `[bio.planemo]` in `case.toml`, with optional
//! inputs JSON and pass-through extras. `action` is constrained to
//! `run` / `test` / `lint` at parse time so the adapter doesn't
//! forward unsupported subcommands.
//!
//! `collect()` walks the workdir and surfaces `<output_basename>*.html`
//! reports as `Native`, `*.json` files as `Tabular`, and `*.log` files
//! as `Log` — the standard set of artefacts Planemo drops next to the
//! workflow during a run / test invocation.

#![forbid(unsafe_code)]
#![allow(missing_docs)]

pub mod case_input;

use std::ffi::OsString;
use std::fs;
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

use crate::case_input::PlanemoInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(PlanemoAdapter::new())
}

pub struct PlanemoAdapter;

impl PlanemoAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for PlanemoAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "planemo";
/// Planemo's binary candidates. Bioconda + pip both install the
/// canonical lowercase entry-point.
const BINARIES: &[&str] = &["planemo"];

impl Adapter for PlanemoAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "Planemo (Galaxy)",
            // 0.75 (early 2023) is the floor where the modern
            // `workflow_run` semantics and Galaxy 23.0+ compatibility
            // stabilised; 1.0 reserves room for the next major bump.
            version_range: VersionRange {
                min_inclusive: Version::new(0, 75, 0),
                max_exclusive: Version::new(1, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            tool_license: "AFL-3.0",
            docs_url: "https://planemo.readthedocs.io/",
            homepage_url: "https://github.com/galaxyproject/planemo",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(BINARIES) {
            Some(binary_path) => {
                // `planemo --version` prints the bare semver on
                // stdout — same convention as every other Python
                // entry-point.
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
                hint: "Planemo 0.75+ required; install via \
                       `pip install planemo`, `conda install -c bioconda planemo`, \
                       or see https://planemo.readthedocs.io/en/latest/installation.html"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = PlanemoInput::from_case_dir(&case.path)?;

        // Round-4 security: reject `output_basename = "../etc/passwd"`
        // and friends before the value flows into any path join.
        // Same pattern as the round-3 fix in bionetgen/iqtree/art/fasttree.
        valenx_core::adapter_helpers::validate_output_basename(
            &input.output_basename,
            "[bio.planemo].output_basename",
        )
        .map_err(|e| AdapterError::InvalidCase {
            case_path: case.path.join("case.toml"),
            reason: format!("{e}"),
        })?;

        fs::create_dir_all(workdir)?;

        // Stage `case.toml` into the workdir so collect() can recover
        // the configured `output_basename` for prefix-filtering output
        // artifacts. Without this stage, the basename filter silently
        // degrades to "match everything".
        let staged_case_toml = workdir.join("case.toml");
        let source_case_toml = case.path.join("case.toml");
        if source_case_toml.is_file() {
            fs::copy(&source_case_toml, &staged_case_toml)
                .map_err(|e| AdapterError::Other(anyhow::anyhow!("stage case.toml: {e}")))?;
        }

        // Resolve workflow against case dir if relative. Same
        // convention as every other Phase 17/18/22 bio adapter.
        let source_workflow = if input.workflow.is_absolute() {
            input.workflow.clone()
        } else {
            valenx_core::adapter_helpers::confined_join(
            &case.path,
            &input.workflow,
        )?
        };
        if !source_workflow.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.planemo].workflow `{}` not found (resolved {})",
                    input.workflow.display(),
                    source_workflow.display()
                ),
            });
        }

        // Resolve optional inputs path against case dir via
        // `confined_join` — sibling-field sweep round-8: `workflow` got
        // it in round-6, `inputs` needs the same sandboxing.
        let source_inputs = match &input.inputs {
            Some(p) => {
                let resolved = valenx_core::adapter_helpers::confined_join(&case.path, p)?;
                if !resolved.is_file() {
                    return Err(AdapterError::InvalidCase {
                        case_path: case.path.join("case.toml"),
                        reason: format!(
                            "[bio.planemo].inputs `{}` not found (resolved {})",
                            p.display(),
                            resolved.display()
                        ),
                    });
                }
                Some(resolved)
            }
            None => None,
        };

        let binary_path = find_on_path(BINARIES).ok_or_else(|| AdapterError::ToolNotInstalled {
            name: INFO_ID,
            hint: "Planemo 0.75+ required; install via \
                       `pip install planemo`, `conda install -c bioconda planemo`, \
                       or see https://planemo.readthedocs.io/en/latest/installation.html"
                .into(),
        })?;

        // Compose `planemo <action> <workflow> [inputs] [extras...]`.
        let mut native_command: Vec<OsString> = vec![
            binary_path.into_os_string(),
            OsString::from(&input.action),
            source_workflow.into_os_string(),
        ];
        if let Some(p) = source_inputs {
            native_command.push(p.into_os_string());
        }
        for arg in &input.extra_args {
            native_command.push(OsString::from(arg));
        }

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            // Planemo runs vary widely — a `lint` is sub-second, a
            // full `test` of a multi-tool workflow can take an hour
            // plus while Galaxy spins up. 4 hours is a generous default.
            estimated_runtime: Some(Duration::from_secs(4 * 60 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting Planemo", |line| {
            let mut hint = subprocess::Hint::default();
            // Planemo's stderr / stdout banners include "Galaxy" /
            // "Tool tests" / per-step status. Surface the
            // canonical "All tests passed" / "FAIL" markers as
            // progress / warning hints.
            if line.contains("All tests passed") || line.contains("Workflow finished") {
                hint.progress = Some((95.0, line.to_string()));
            } else if line.contains("Running tool")
                || line.contains("Executing job")
                || line.contains("Step ")
            {
                hint.progress = Some((50.0, line.to_string()));
            } else if line.contains("Traceback") || line.contains("FAIL") || line.contains("Error:")
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
        // Provenance: hash the case.toml — Planemo workdirs are a
        // mix of tool outputs, Galaxy state, and logs rather than a
        // single canonical file.
        let case_hash_input = job.workdir.join("case.toml");
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "Planemo",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);
        let mut artefacts: Vec<Artifact> = Vec::new();

        // Lift the configured output_basename out of case.toml so
        // collect() filters HTML reports by prefix; falls back to the
        // empty prefix when the file isn't present yet.
        let output_basename = PlanemoInput::from_case_dir(&job.workdir)
            .map(|i| i.output_basename)
            .unwrap_or_default();

        // Walk the workdir and classify outputs by extension. HTML
        // reports use the `<output_basename>*` filter; JSON and log
        // files are always collected regardless of name.
        let entries = match fs::read_dir(&job.workdir) {
            Ok(e) => e,
            Err(e) => {
                tracing::warn!(target: "valenx-planemo", ?e, "workdir read failed");
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
                .file_name()
                .and_then(|s| s.to_str())
                .unwrap_or_default();
            match ext.as_deref() {
                Some("html") if stem.starts_with(&output_basename) => artefacts.push(Artifact {
                    path,
                    kind: ArtifactKind::Native,
                    checksum: None,
                    label: "Planemo report".to_string(),
                }),
                Some("json") => artefacts.push(Artifact {
                    path,
                    kind: ArtifactKind::Tabular,
                    checksum: None,
                    label: "Planemo run JSON".to_string(),
                }),
                Some("log") => artefacts.push(Artifact {
                    path,
                    kind: ArtifactKind::Log,
                    checksum: None,
                    label: "Planemo log".to_string(),
                }),
                _ => continue,
            }
        }

        artefacts.sort_by(|a, b| a.path.cmp(&b.path));
        results.artifacts = artefacts;
        Ok(results)
    }

    fn capabilities(&self) -> Capabilities {
        Capabilities {
            capabilities: Vec::new(),
            ribbon_contributions: vec!["bio.planemo.run"],
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn info_is_bio_domain() {
        let info = PlanemoAdapter::new().info();
        assert_eq!(info.id, "planemo");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "AFL-3.0");
        assert_eq!(info.display_name, "Planemo (Galaxy)");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = PlanemoAdapter::new().info();
        // 0.75 is the modern-flag floor; 1.0 reserves the next major
        // bump.
        assert_eq!(info.version_range.min_inclusive, Version::new(0, 75, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(1, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = PlanemoAdapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.planemo.run"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = PlanemoAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }

    #[test]
    fn prepare_rejects_traversal_inputs_path() {
        // Round-8 RED→GREEN: `[bio.planemo].inputs` (test data /
        // inputs file) now routes through `confined_join` — sibling
        // to `workflow`, which got it in round-6.
        //
        // We use `../etc/passwd` (relative traversal) for cross-platform
        // portability — see the bcftools test for the rationale.
        use valenx_test_utils::tempdir;
        let d = tempdir("planemo-traversal");
        std::fs::write(d.join("workflow.ga"), b"{}").unwrap();
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "planemo.run"

[bio.planemo]
workflow        = "workflow.ga"
inputs          = "../etc/passwd"
output_basename = "report"
"#,
        )
        .unwrap();
        let case = Case {
            id: "planemo-traversal".into(),
            path: d.clone(),
        };
        let workdir = d.join("workdir");
        let err = PlanemoAdapter::new().prepare(&case, &workdir).unwrap_err();
        let msg = format!("{err}");
        assert!(
            msg.contains("absolute") || msg.contains("escape") || msg.contains("`..`") || msg.contains("traversal"),
            "expected confined_join rejection on inputs, got: {msg}"
        );
        let _ = std::fs::remove_dir_all(&d);
    }
}
