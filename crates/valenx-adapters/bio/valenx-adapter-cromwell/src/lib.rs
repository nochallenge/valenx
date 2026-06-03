//! # valenx-adapter-cromwell
//!
//! Adapter for [Cromwell](https://github.com/broadinstitute/cromwell) —
//! the Broad Institute's WDL workflow execution engine. Cromwell ships
//! as a single Java JAR (`cromwell-<version>.jar`); the user supplies
//! the path to that JAR via `[bio.cromwell].jar` in `case.toml` and
//! the adapter composes a `java -jar <jar> <action> <workflow>
//! [-i <inputs>] [extras...]` invocation.
//!
//! **Phase 22.5 — JAR-distributed workflow-CLI wrapper.** Sister to
//! the Jalview / j5 / Cello adapters (probe `java`, jar via case
//! input). `action` is constrained to `run` / `submit` / `validate`
//! at parse time so the adapter never forwards unsupported
//! subcommands. The `-i <inputs>` flag pair is only emitted when
//! `inputs` is `Some(_)`.
//!
//! `collect()` walks the **top level** of the workdir for
//! `<output_basename>*.json` (Tabular, "Cromwell metadata") and
//! `*.log` (Log, "Cromwell log") — Cromwell drops both alongside
//! the workflow run by default.

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

use crate::case_input::CromwellInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(CromwellAdapter::new())
}

pub struct CromwellAdapter;

impl CromwellAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for CromwellAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "cromwell";
/// Cromwell is JAR-distributed — we probe `java` itself, not a
/// `cromwell` launcher. The user supplies the jar path via case input.
const BINARIES: &[&str] = &["java"];

impl Adapter for CromwellAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "Cromwell",
            // Cromwell's modern release line started at 80 (2023);
            // 100 reserves room for the next decade of point releases
            // without pinning a too-narrow band.
            version_range: VersionRange {
                min_inclusive: Version::new(80, 0, 0),
                max_exclusive: Version::new(100, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            tool_license: "BSD-3-Clause",
            docs_url: "https://cromwell.readthedocs.io/",
            homepage_url: "https://github.com/broadinstitute/cromwell",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(BINARIES) {
            Some(binary_path) => Ok(ProbeReport {
                ok: true,
                // Cromwell's version comes from the jar filename / manifest,
                // not from `java`; we surface no version here. The user
                // pins the Cromwell release implicitly by the jar they
                // point at.
                found_version: None,
                binary_path: Some(binary_path),
                warnings: vec![
                    "probe found `java` on PATH but cannot verify the cromwell-*.jar \
                     release without invoking it; ensure `[bio.cromwell].jar` points \
                     at a valid Cromwell distribution"
                        .into(),
                ],
                required_env: Vec::new(),
            }),
            None => Err(AdapterError::ToolNotInstalled {
                name: INFO_ID,
                hint: "Java 11+ JRE required to run Cromwell; install via your \
                       package manager (`apt install default-jre`, \
                       `brew install openjdk`, etc.) and ensure `java` is \
                       on PATH"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = CromwellInput::from_case_dir(&case.path)?;

        // Round-4 security: reject `output_basename = "../etc/passwd"`
        // and friends before the value flows into any path join.
        // Same pattern as the round-3 fix in bionetgen/iqtree/art/fasttree.
        valenx_core::adapter_helpers::validate_output_basename(
            &input.output_basename,
            "[bio.cromwell].output_basename",
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

        // Resolve the jar path against the case directory if relative.
        // Almost always absolute (jars live under /opt or similar),
        // but support the relative form too for self-contained cases.
        // Round-9 hardening: the jar value flows straight into
        // `java -jar <jar>`, so a hostile case.toml could otherwise
        // point a *relative* form at `../../some-arbitrary.jar` and
        // turn "Run case" into arbitrary jar exec. Wrap the relative
        // branch with `confined_join`; absolute paths remain explicit
        // because the standard install layout puts Cromwell under
        // `/opt/cromwell/cromwell.jar` (admin-managed, not user data).
        let source_jar = if input.jar.is_absolute() {
            input.jar.clone()
        } else {
            valenx_core::adapter_helpers::confined_join(&case.path, &input.jar)?
        };
        if !source_jar.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.cromwell].jar `{}` not found (resolved {})",
                    input.jar.display(),
                    source_jar.display()
                ),
            });
        }

        // Resolve workflow against case dir if relative.
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
                    "[bio.cromwell].workflow `{}` not found (resolved {})",
                    input.workflow.display(),
                    source_workflow.display()
                ),
            });
        }

        // Resolve optional inputs JSON against case dir via
        // `confined_join` — sibling-field sweep round-8: `workflow` got
        // it in round-6, `inputs` needs the same sandboxing.
        let source_inputs = match &input.inputs {
            Some(p) => {
                let resolved = valenx_core::adapter_helpers::confined_join(&case.path, p)?;
                if !resolved.is_file() {
                    return Err(AdapterError::InvalidCase {
                        case_path: case.path.join("case.toml"),
                        reason: format!(
                            "[bio.cromwell].inputs `{}` not found (resolved {})",
                            p.display(),
                            resolved.display()
                        ),
                    });
                }
                Some(resolved)
            }
            None => None,
        };

        let java_path = find_on_path(BINARIES).ok_or_else(|| AdapterError::ToolNotInstalled {
            name: INFO_ID,
            hint: "Java 11+ JRE required to run Cromwell; install via your \
                   package manager and ensure `java` is on PATH"
                .into(),
        })?;

        // Compose `java -jar <jar> <action> <workflow> [-i <inputs>]
        // [extras...]`. The `-i <inputs>` flag is two separate args
        // and only emitted when an inputs JSON is configured.
        let mut native_command: Vec<OsString> = vec![
            java_path.into_os_string(),
            OsString::from("-jar"),
            source_jar.into_os_string(),
            OsString::from(&input.action),
            source_workflow.into_os_string(),
        ];
        if let Some(p) = source_inputs {
            native_command.push(OsString::from("-i"));
            native_command.push(p.into_os_string());
        }
        for arg in &input.extra_args {
            native_command.push(OsString::from(arg));
        }

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            // Cromwell runs vary widely — `validate` is sub-second,
            // a `run` of a multi-step WDL pipeline can take many
            // hours while individual tools execute. 4 hours is a
            // generous default (matches Planemo / Snakemake).
            estimated_runtime: Some(Duration::from_secs(4 * 60 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting Cromwell", |line| {
            let mut hint = subprocess::Hint::default();
            // Cromwell's stdout banners follow a predictable shape:
            // "Workflow ... succeeded" / "Workflow ... failed" at end
            // of run, "Starting calls" / "job id" mid-run, and the
            // usual JVM Exception / Error markers when something
            // blows up.
            if line.contains("succeeded") || line.contains("Workflow finished") {
                hint.progress = Some((95.0, line.to_string()));
            } else if line.contains("Starting calls")
                || line.contains("job id")
                || line.contains("WorkflowExecutionActor")
            {
                hint.progress = Some((50.0, line.to_string()));
            } else if line.contains("Exception")
                || line.contains("ERROR")
                || line.contains("failed")
                || line.contains("java.lang.")
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
        // descriptor. Cromwell workdirs are a mix of metadata JSON,
        // logs, and task outputs rather than a single canonical file.
        let case_hash_input = job.workdir.join("case.toml");
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "Cromwell",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);
        let mut artefacts: Vec<Artifact> = Vec::new();

        // Lift `output_basename` out of the staged case.toml so
        // collect() can filter metadata JSON files by stem prefix —
        // falls back to the empty prefix when the file isn't there yet
        // (best-effort) and accepts every recognised file.
        let output_basename = CromwellInput::from_case_dir(&job.workdir)
            .map(|i| i.output_basename)
            .unwrap_or_default();

        // Walk only the top level of the workdir — Cromwell drops the
        // metadata JSON and engine log alongside the workflow root.
        // Per-task subdirectories under `cromwell-executions/` are
        // out of scope for this adapter (too deep to enumerate
        // safely without a recursive walker).
        let entries = match fs::read_dir(&job.workdir) {
            Ok(e) => e,
            Err(e) => {
                tracing::warn!(target: "valenx-cromwell", ?e, "workdir read failed");
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
                Some("json") if stem.starts_with(&output_basename) => artefacts.push(Artifact {
                    path,
                    kind: ArtifactKind::Tabular,
                    checksum: None,
                    label: "Cromwell metadata".to_string(),
                }),
                Some("log") => artefacts.push(Artifact {
                    path,
                    kind: ArtifactKind::Log,
                    checksum: None,
                    label: "Cromwell log".to_string(),
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
            ribbon_contributions: vec!["bio.cromwell.run"],
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn info_is_bio_domain() {
        let info = CromwellAdapter::new().info();
        assert_eq!(info.id, "cromwell");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "BSD-3-Clause");
        assert_eq!(info.display_name, "Cromwell");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = CromwellAdapter::new().info();
        // Modern Cromwell line started at 80 (2023); 100 reserves
        // room for the next decade of point releases.
        assert_eq!(info.version_range.min_inclusive, Version::new(80, 0, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(100, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = CromwellAdapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.cromwell.run"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = CromwellAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }

    #[test]
    fn prepare_rejects_traversal_inputs_path() {
        // Round-8 RED→GREEN: `[bio.cromwell].inputs` (the workflow
        // inputs JSON) now routes through `confined_join` — sibling to
        // `workflow`, which got it in round-6.
        //
        // Round-9 follow-up: `jar` ALSO routes through `confined_join`
        // for its relative form now — see the sibling test
        // `prepare_rejects_jar_traversing_outside_case_dir`. Absolute
        // jar paths still flow through unchanged because the typical
        // install lives under `/opt/cromwell/cromwell-86.jar`.
        //
        // We use `../etc/passwd` (relative traversal) for cross-platform
        // portability — see the bcftools test for the rationale.
        use valenx_test_utils::tempdir;
        let d = tempdir("cromwell-traversal");
        std::fs::write(d.join("workflow.wdl"), b"workflow {}").unwrap();
        // Plausible jar path the case_input parser will accept (the
        // file-exists check fires only inside prepare()), but
        // confined_join on `inputs` should fail first.
        std::fs::write(d.join("dummy.jar"), b"PK\x03\x04").unwrap();
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "cromwell.run"

[bio.cromwell]
jar             = "dummy.jar"
workflow        = "workflow.wdl"
inputs          = "../etc/passwd"
output_basename = "metadata"
"#,
        )
        .unwrap();
        let case = Case {
            id: "cromwell-traversal".into(),
            path: d.clone(),
        };
        let workdir = d.join("workdir");
        let err = CromwellAdapter::new()
            .prepare(&case, &workdir)
            .unwrap_err();
        let msg = format!("{err}");
        assert!(
            msg.contains("absolute") || msg.contains("escape") || msg.contains("`..`") || msg.contains("traversal"),
            "expected confined_join rejection on inputs, got: {msg}"
        );
        let _ = std::fs::remove_dir_all(&d);
    }

    /// Round-9 RED→GREEN: `[bio.cromwell].jar`, in its relative form,
    /// used to be joined with bare `case.path.join` — letting a hostile
    /// case supply `jar = "../../my-evil.jar"` and have `java -jar`
    /// execute whatever was placed there. Wrap with `confined_join`;
    /// absolute paths (the documented `/opt/cromwell/cromwell-86.jar`
    /// install layout) are left to the existing branch.
    #[test]
    fn prepare_rejects_jar_traversing_outside_case_dir() {
        use valenx_test_utils::tempdir;
        let d = tempdir("cromwell-jar-trav");
        std::fs::write(d.join("workflow.wdl"), b"workflow {}").unwrap();
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "cromwell.run"

[bio.cromwell]
jar             = "../../my-evil.jar"
workflow        = "workflow.wdl"
output_basename = "metadata"
"#,
        )
        .unwrap();
        let case = Case {
            id: "cromwell-jar-trav".into(),
            path: d.clone(),
        };
        let workdir = d.join("workdir");
        let err = CromwellAdapter::new()
            .prepare(&case, &workdir)
            .unwrap_err();
        let msg = format!("{err}");
        assert!(
            msg.contains("absolute") || msg.contains("escape") || msg.contains("`..`") || msg.contains("traversal") || msg.contains("stay within"),
            "expected confined_join rejection on jar, got: {msg}"
        );
        let _ = std::fs::remove_dir_all(&d);
    }
}
