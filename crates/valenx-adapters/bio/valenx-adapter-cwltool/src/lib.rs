//! # valenx-adapter-cwltool
//!
//! Adapter for [cwltool](https://github.com/common-workflow-language/cwltool) —
//! the reference implementation of the
//! [Common Workflow Language](https://www.commonwl.org/) (CWL).
//! cwltool runs CWL `.cwl` tool / workflow documents against an
//! optional input-object document, staging tools either in-process or
//! through a container runtime (Docker / Singularity / podman) per the
//! workflow's `DockerRequirement` hint.
//!
//! **Phase 22.5 — workflow-managers expansion.** Sister adapter to
//! Snakemake / Nextflow / Planemo / Cromwell on the workflow-runner
//! surface. The adapter composes
//!
//! ```text
//! cwltool --outdir <output_dir> [extras...] <workflow> [inputs]
//! ```
//!
//! and resolves both the workflow document and the optional input
//! object against `case.path`.
//!
//! `collect()` walks **one level deep** into the configured
//! `<output_dir>/` subdir for any file (cwltool's outputs are
//! workflow-defined and span every conceivable extension — BAM, VCF,
//! tabular CSV, plots, structured JSON, etc.) and surfaces top-level
//! `*.log` files as `Log` artefacts.

#![forbid(unsafe_code)]
#![allow(missing_docs)]

pub mod case_input;

use std::ffi::OsString;
use std::fs;
use std::path::Path;
use std::time::Duration;

use semver::Version;

use valenx_core::{
    adapter_helpers::{
        detect_tool_version_semver, find_on_path, live_provenance, validate_output_dir,
    },
    error::RunPhase,
    subprocess, Adapter, AdapterError, AdapterInfo, Capabilities, Case, LicenseMode, Physics,
    PreparedJob, ProbeReport, RunContext, RunReport, VersionRange,
};
use valenx_fields::{
    artifact::{Artifact, ArtifactKind},
    Results,
};

use crate::case_input::CwltoolInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(CwltoolAdapter::new())
}

pub struct CwltoolAdapter;

impl CwltoolAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for CwltoolAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "cwltool";
/// cwltool ships a single console-script entry-point (`cwltool`); the
/// canonical lowercase form is what `pip install cwltool` and the
/// bioconda recipe both expose.
const BINARIES: &[&str] = &["cwltool"];
/// Python interpreter probe — used to surface a more helpful warning
/// when Python is reachable but the cwltool entry-point isn't (the
/// "you forgot to `pip install cwltool`" case).
const PYTHON_BINARIES: &[&str] = &["python3", "python"];
/// Probe-warning surfaced when Python is on PATH but the `cwltool`
/// console-script entry-point isn't. The probe ribbon UI greps for
/// the `cwltool not found` substring to render a missing-tool hint;
/// keep that phrase verbatim.
const PYTHON_ONLY_WARNING: &str = "cwltool not found on PATH; install via `pip install cwltool`";

impl Adapter for CwltoolAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "cwltool",
            // cwltool 3.1 (2020) is the modern release line tracking
            // CWL v1.2; 4.0 reserves room for the next major bump.
            version_range: VersionRange {
                min_inclusive: Version::new(3, 1, 0),
                max_exclusive: Version::new(4, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            tool_license: "Apache-2.0",
            docs_url: "https://cwltool.readthedocs.io/",
            homepage_url: "https://github.com/common-workflow-language/cwltool",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        // Prefer the `cwltool` console-script entry-point. If only
        // Python is reachable we surface a warning (rather than
        // erroring) so users with a Python environment ready but no
        // `cwltool` package see the targeted install hint without
        // failing the probe — same shape as the CellProfiler adapter.
        if let Some(binary_path) = find_on_path(BINARIES) {
            // `cwltool --version` prints `<path> <version>` (e.g.
            // `/usr/local/bin/cwltool 3.1.20240508115724`); the
            // semver helper extracts the version segment.
            let found_version = detect_tool_version_semver(&binary_path, &["--version", "-v"]);
            return Ok(ProbeReport {
                ok: true,
                found_version,
                binary_path: Some(binary_path),
                warnings: Vec::new(),
                required_env: Vec::new(),
            });
        }
        if let Some(python_path) = find_on_path(PYTHON_BINARIES) {
            // Python is on PATH but the dedicated `cwltool` console
            // script isn't — surface the targeted install hint as a
            // warning. The probe still reports `ok` so downstream
            // tooling can act on it.
            return Ok(ProbeReport {
                ok: true,
                found_version: None,
                binary_path: Some(python_path),
                warnings: vec![PYTHON_ONLY_WARNING.into()],
                required_env: Vec::new(),
            });
        }
        Err(AdapterError::ToolNotInstalled {
            name: INFO_ID,
            hint: "cwltool 3.1+ required; install via \
                   `pip install cwltool`, \
                   `conda install -c bioconda cwltool`, \
                   or see https://cwltool.readthedocs.io/en/latest/#install"
                .into(),
        })
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = CwltoolInput::from_case_dir(&case.path)?;

        fs::create_dir_all(workdir)?;

        // Round-5: validate the `output_dir` field BEFORE we resolve
        // any paths. A hostile case.toml could otherwise set
        // `output_dir = "../../etc/cron.d"` and cwltool would write
        // workflow outputs anywhere it has permission for.
        validate_output_dir(
            std::path::Path::new(&input.output_dir),
            "[bio.cwltool].output_dir",
        )?;

        // Resolve workflow against case dir if relative. Same
        // convention as every other Phase 17 / 18 / 22 / 22.5 bio
        // adapter.
        let source_workflow = if input.workflow.is_absolute() {
            input.workflow.clone()
        } else {
            valenx_core::adapter_helpers::confined_join(&case.path, &input.workflow)?
        };
        if !source_workflow.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.cwltool].workflow `{}` not found (resolved {})",
                    input.workflow.display(),
                    source_workflow.display()
                ),
            });
        }

        // Resolve optional inputs path against case dir. cwltool
        // accepts either YAML or JSON for the input object; we just
        // pass the path through. Round-5: wrap with confined_join so
        // a hostile `inputs = "../../etc/passwd.yml"` is rejected
        // exactly like the sister `workflow` field.
        let source_inputs = match &input.inputs {
            Some(p) => {
                let resolved = if p.is_absolute() {
                    p.clone()
                } else {
                    valenx_core::adapter_helpers::confined_join(&case.path, p)?
                };
                if !resolved.is_file() {
                    return Err(AdapterError::InvalidCase {
                        case_path: case.path.join("case.toml"),
                        reason: format!(
                            "[bio.cwltool].inputs `{}` not found (resolved {})",
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
            hint: "cwltool 3.1+ required; install via \
                       `pip install cwltool`, \
                       `conda install -c bioconda cwltool`, \
                       or see https://cwltool.readthedocs.io/en/latest/#install"
                .into(),
        })?;

        // Compose
        // `cwltool --outdir <output_dir> [extras...] <workflow> [inputs]`.
        // `--outdir` is workdir-relative — the subprocess runner's cwd
        // is the workdir, so cwltool resolves the basename correctly
        // and writes into `<workdir>/<output_dir>/`.
        let mut native_command: Vec<OsString> = vec![
            binary_path.into_os_string(),
            OsString::from("--outdir"),
            OsString::from(&input.output_dir),
        ];

        // Round-4 fix: extra_args after positionals — see
        // security/code-review.md. cwltool's CLI grammar requires the
        // workflow first and the input object second; pushing extras
        // before them would let a hostile case.toml supply something
        // like `extra_args = ["--print-deps"]` that swallows the
        // workflow path positional.
        native_command.push(source_workflow.into_os_string());
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
            // CWL workflows vary as widely as Nextflow / Snakemake
            // pipelines; cap the estimate at 24 hours for the long
            // tail of multi-tool genomics workflows.
            estimated_runtime: Some(Duration::from_secs(24 * 60 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting cwltool", |line| {
            let mut hint = subprocess::Hint::default();
            // cwltool's standard run banner sequence:
            //   "[job <name>] /tmp/...$ <cmd>"           (per-step start)
            //   "[step <name>] start"                    (workflow step)
            //   "[workflow <id>] completed success"      (success tail)
            //   "Final process status is success"        (run tail)
            // Errors surface as "WorkflowException" / "permanentFailure".
            if line.contains("Final process status is success")
                || line.contains("completed success")
            {
                hint.progress = Some((95.0, line.to_string()));
            } else if line.contains("[step ")
                || line.contains("[job ")
                || line.contains("[workflow ")
            {
                hint.progress = Some((50.0, line.to_string()));
            } else if line.contains("WorkflowException")
                || line.contains("permanentFailure")
                || line.contains("Final process status is permanentFail")
                || line.contains("ERROR")
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
        // Provenance: hash the case.toml — like Snakemake / Nextflow,
        // a CWL run produces a tree of workflow-defined outputs rather
        // than a single canonical file.
        let case_hash_input = job.workdir.join("case.toml");
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "cwltool",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);
        let mut artefacts: Vec<Artifact> = Vec::new();

        // Top-level: cwltool's logging defaults to stderr but users
        // routinely tee into a `*.log` file via `--log-stderr` or
        // shell redirection. Surface any top-level `*.log` regardless
        // of basename — same convention as the CellProfiler adapter.
        let entries = match fs::read_dir(&job.workdir) {
            Ok(e) => e,
            Err(e) => {
                tracing::warn!(target: "valenx-cwltool", ?e, "workdir read failed");
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
                    label: "cwltool log".to_string(),
                });
            }
        }

        // One level deep: walk into the configured `output_dir`
        // subdir for cwltool's workflow outputs. The output-document
        // shapes are workflow-author-defined (BAM / VCF / CSV / JSON /
        // PNG / arbitrary tool stdout), so we surface every regular
        // file rather than filtering by extension. We re-read the
        // case.toml to recover the configured `output_dir` without
        // needing to thread the input struct into collect().
        let output_dir_name = read_output_dir(&job.workdir);
        if let Some(name) = output_dir_name {
            let output_dir = job.workdir.join(&name);
            if output_dir.is_dir() {
                let inner = match fs::read_dir(&output_dir) {
                    Ok(e) => e,
                    Err(e) => {
                        tracing::warn!(
                            target: "valenx-cwltool",
                            ?e,
                            output_dir = %output_dir.display(),
                            "output dir read failed"
                        );
                        artefacts.sort_by(|a, b| a.path.cmp(&b.path));
                        results.artifacts = artefacts;
                        return Ok(results);
                    }
                };
                for entry in inner.flatten() {
                    let path = entry.path();
                    if !path.is_file() {
                        continue;
                    }
                    artefacts.push(Artifact {
                        path,
                        kind: ArtifactKind::Native,
                        checksum: None,
                        label: "cwltool output".to_string(),
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
            ribbon_contributions: vec!["bio.cwltool.run"],
        }
    }
}

/// Re-read the `[bio.cwltool].output_dir` from a staged `case.toml`
/// for collect()-time output-tree filtering. Returns `None` when the
/// case.toml is missing or unparseable — collect() then returns just
/// the top-level log artefacts (a sensible degraded mode).
fn read_output_dir(workdir: &Path) -> Option<String> {
    // Round-23 sweep: bound staged case.toml at MAX_PROJECT_FILE_BYTES.
    let text = valenx_core::io_caps::read_capped_to_string(
        &workdir.join("case.toml"),
        valenx_core::project::loader::MAX_PROJECT_FILE_BYTES as usize,
    )
    .ok()?;
    let parsed: toml::Value = toml::from_str(&text).ok()?;
    parsed
        .get("bio")?
        .get("cwltool")?
        .get("output_dir")?
        .as_str()
        .map(str::to_string)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn info_is_bio_domain() {
        let info = CwltoolAdapter::new().info();
        assert_eq!(info.id, "cwltool");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "Apache-2.0");
        assert_eq!(info.display_name, "cwltool");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = CwltoolAdapter::new().info();
        // cwltool 3.1 (2020) tracks CWL v1.2 — the first stable
        // release line that's persisted across the bioconda /
        // common-workflow-language ecosystem; 4.0 reserves the next
        // major bump.
        assert_eq!(info.version_range.min_inclusive, Version::new(3, 1, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(4, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = CwltoolAdapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.cwltool.run"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = CwltoolAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }

    /// Round-5 RED→GREEN: the `inputs` field of [bio.cwltool] used to
    /// be joined with bare `case.path.join(p)`, letting a hostile case
    /// supply `inputs = "../../etc/passwd.yml"` and have cwltool
    /// happily read whatever file the user has access to. The fix
    /// wraps the join with `confined_join` (same as `workflow`).
    #[test]
    fn prepare_rejects_inputs_traversing_outside_case_dir() {
        let case_dir = std::env::temp_dir().join(format!(
            "valenx-cwltool-trav-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&case_dir).unwrap();
        // Write a case.toml with a `..`-traversing inputs path. The
        // workflow path is a plain relative one so the workflow check
        // doesn't fire first.
        std::fs::write(case_dir.join("workflow.cwl"), b"# placeholder").unwrap();
        std::fs::write(
            case_dir.join("case.toml"),
            br#"[case]
physics = "bio"
solver  = "cwltool.run"

[bio.cwltool]
workflow   = "workflow.cwl"
inputs     = "../../etc/passwd.yml"
output_dir = "results"
"#,
        )
        .unwrap();
        let workdir = case_dir.join("workdir");
        let case = Case {
            id: "test".to_string(),
            path: case_dir.clone(),
        };
        let adapter = CwltoolAdapter::new();
        let err = adapter
            .prepare(&case, &workdir)
            .expect_err("must reject ../../etc/passwd.yml inputs");
        let msg = format!("{err}");
        assert!(
            msg.contains("..") || msg.contains("stay within"),
            "msg: {msg}"
        );
        let _ = std::fs::remove_dir_all(&case_dir);
    }

    /// Round-5 RED→GREEN: same idea for `output_dir`, which used to
    /// flow straight into cwltool's `--outdir` without validation.
    /// The fix calls `validate_output_dir` BEFORE resolving paths so
    /// the rejection lands as `Other` with a clear `..` message.
    #[test]
    fn prepare_rejects_output_dir_traversing_outside_workdir() {
        let case_dir = std::env::temp_dir().join(format!(
            "valenx-cwltool-out-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&case_dir).unwrap();
        std::fs::write(case_dir.join("workflow.cwl"), b"# placeholder").unwrap();
        std::fs::write(
            case_dir.join("case.toml"),
            br#"[case]
physics = "bio"
solver  = "cwltool.run"

[bio.cwltool]
workflow   = "workflow.cwl"
output_dir = "../escape"
"#,
        )
        .unwrap();
        let workdir = case_dir.join("workdir");
        let case = Case {
            id: "test".to_string(),
            path: case_dir.clone(),
        };
        let adapter = CwltoolAdapter::new();
        let err = adapter
            .prepare(&case, &workdir)
            .expect_err("must reject ../escape output_dir");
        let msg = format!("{err}");
        assert!(msg.contains(".."), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&case_dir);
    }

    #[test]
    fn python_only_warning_anchors_spec_phrase() {
        // The Phase 22.5 spec requires probe() to surface the
        // verbatim phrase "cwltool not found" (and a `pip install
        // cwltool` install hint) as a *warning* — not an error —
        // when Python is on PATH but the cwltool console-script
        // entry-point isn't. The probe ribbon UI greps for that
        // substring to render the missing-tool hint, so a future
        // refactor that drops the phrase or downgrades the
        // construction to an error path must trip a test failure.
        assert!(
            PYTHON_ONLY_WARNING.contains("cwltool not found"),
            "warning must contain `cwltool not found` anchor; got: {PYTHON_ONLY_WARNING}"
        );
        assert!(
            PYTHON_ONLY_WARNING.contains("pip install cwltool"),
            "warning must contain `pip install cwltool` install hint; \
             got: {PYTHON_ONLY_WARNING}"
        );
    }
}
