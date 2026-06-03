//! # valenx-adapter-nextflow
//!
//! Adapter for [Nextflow](https://www.nextflow.io/) — Paolo Di
//! Tommaso's reactive workflow language and runtime, the de-facto
//! orchestrator for nf-core's curated bioinformatics pipelines (rnaseq,
//! sarek, ampliseq, ...). Nextflow processes form a DAG of tasks
//! shipped through cloud / cluster / containerised executors; the
//! Valenx adapter just runs `nextflow run <pipeline>` and reports the
//! workdir as the produced artifact bundle.
//!
//! **Phase 22 — pipeline-orchestrator wrapper.** Unlike single-tool
//! adapters (BWA, samtools), the wrapped CLI runs an *internal*
//! pipeline that itself dispatches other bio tools. The adapter's job
//! is composition: stitch a `nextflow run` invocation from
//! `[bio.nextflow]` in `case.toml` (pipeline ref, profile, config,
//! `--<key> <value>` params, `-resume`) and surface Nextflow's
//! standard observability outputs — `report.html`, `timeline.html`,
//! `dag.svg` — alongside the workdir.
//!
//! `params` lands in deterministic alphabetical order via a
//! `BTreeMap`, so two runs of the same case produce byte-identical
//! command lines (matters for reproducibility and for diff stability).

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

use crate::case_input::NextflowInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(NextflowAdapter::new())
}

pub struct NextflowAdapter;

impl NextflowAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for NextflowAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "nextflow";
/// Nextflow's binary candidates. The launcher is the canonical
/// install name from Bioconda, the standalone bundle, and the
/// `nextflow self-update` distribution; `nf` is a community alias
/// that some installers ship as a convenience symlink.
const BINARIES: &[&str] = &["nextflow"];

impl Adapter for NextflowAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "Nextflow",
            // 23.10.0 is the long-term support release (Oct 2023);
            // newer 24.x's are stable — we cap below 25.0 to reserve
            // room for the next major bump.
            version_range: VersionRange {
                min_inclusive: Version::new(23, 10, 0),
                max_exclusive: Version::new(25, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            tool_license: "Apache-2.0",
            docs_url: "https://www.nextflow.io/docs/latest/",
            homepage_url: "https://www.nextflow.io/",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(BINARIES) {
            Some(binary_path) => {
                // `nextflow -version` is the canonical incantation
                // (single-dash); newer builds also accept `--version`.
                let found_version =
                    detect_tool_version_semver(&binary_path, &["-version", "--version"]);
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
                hint: "Nextflow 23.10+ required; install via \
                       `conda install -c bioconda nextflow`, \
                       `curl -s https://get.nextflow.io | bash`, or \
                       see https://www.nextflow.io/docs/latest/install.html"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = NextflowInput::from_case_dir(&case.path)?;

        fs::create_dir_all(workdir)?;

        // Pipeline path resolution: a value containing `/`, `\`, or
        // ending in `.nf` is treated as a path and resolved against
        // the case dir if relative; everything else is passed through
        // verbatim so registry slugs like `nf-core/rnaseq` reach the
        // CLI unchanged.
        //
        // Round-8 sibling-field sweep: when the path-shaped branch
        // fires, relative paths go through `confined_join` to refuse
        // `..` traversal escapes; absolute paths are forwarded
        // verbatim (Nextflow accepts absolute `.nf` paths for system
        // installs). Path validation runs *before* the binary probe
        // so a hostile case bundle gets a sandbox-rejection error
        // even on hosts where Nextflow isn't installed.
        //
        // Note we deliberately don't error on a non-existent `.nf`
        // path here — Nextflow itself produces a clearer error than
        // we could ("could not find pipeline ..."). The job-prep
        // layer just hands the resolved string to the CLI.
        let pipeline_arg = if looks_like_path(&input.pipeline) {
            let p = PathBuf::from(&input.pipeline);
            let resolved = if p.is_absolute() {
                p
            } else {
                confined_join(&case.path, &p)?
            };
            resolved.into_os_string()
        } else {
            OsString::from(&input.pipeline)
        };

        // Resolve config-path early too so traversal-escape rejection
        // doesn't depend on the binary being installed.
        let resolved_config = match &input.config {
            Some(config) => {
                let path = if config.is_absolute() {
                    config.clone()
                } else {
                    confined_join(&case.path, config)?
                };
                Some(path)
            }
            None => None,
        };

        let binary_path = find_on_path(BINARIES).ok_or_else(|| AdapterError::ToolNotInstalled {
            name: INFO_ID,
            hint: "Nextflow 23.10+ required; install via \
                       `conda install -c bioconda nextflow`, \
                       `curl -s https://get.nextflow.io | bash`, or \
                       see https://www.nextflow.io/docs/latest/install.html"
                .into(),
        })?;

        // Compose `nextflow run <pipeline> [-c <config>]
        //                       [-profile <profile>] [-resume]
        //                       [--<key> <value> for each param]
        //                       [extras...]`.
        //
        // BTreeMap iteration ensures `--<key> <value>` pairs land in
        // deterministic alphabetical order regardless of TOML source
        // order — important for reproducibility + diff stability.
        let mut native_command: Vec<OsString> = vec![
            binary_path.into_os_string(),
            OsString::from("run"),
            pipeline_arg,
        ];

        if let Some(config_path) = resolved_config {
            native_command.push(OsString::from("-c"));
            native_command.push(config_path.into_os_string());
        }

        if let Some(profile) = &input.profile {
            native_command.push(OsString::from("-profile"));
            native_command.push(OsString::from(profile));
        }

        if input.resume {
            native_command.push(OsString::from("-resume"));
        }

        // Sorted-by-key thanks to BTreeMap; `--<key>` flags use the
        // double-dash form (Nextflow's pipeline-param convention),
        // distinct from the single-dash CLI control flags above.
        for (key, value) in &input.params {
            native_command.push(OsString::from(format!("--{key}")));
            native_command.push(OsString::from(value));
        }

        for arg in &input.extra_args {
            native_command.push(OsString::from(arg));
        }

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            // Nextflow pipelines vary wildly: a smoke-test
            // `hello-world` finishes in seconds, while a whole-genome
            // sarek run can take days. Cap the estimate at 24 hours
            // — the supervisor will simply not pre-empt within this
            // budget; long-runners pay only the cancellation latency.
            estimated_runtime: Some(Duration::from_secs(24 * 60 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting Nextflow", |line| {
            let mut hint = subprocess::Hint::default();
            // Nextflow's progress output is process-oriented; the
            // most reliable markers are the executor banner at start
            // (`executor >  local`), per-process status lines
            // (`[xx/yyyyyy] process > ...`), and the closing
            // `Completed at:` summary. Lift the obvious ones to
            // progress hints so the UI doesn't sit at 0%.
            if line.contains("Completed at:") || line.contains("Workflow finished") {
                hint.progress = Some((95.0, line.to_string()));
            } else if line.contains("] process >") || line.contains("executor >") {
                hint.progress = Some((50.0, line.to_string()));
            } else if line.contains("ERROR") || line.contains("Error executing process") {
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
        // Provenance: hash the case.toml — a Nextflow workdir is a
        // tree of per-process subdirectories rather than a single
        // canonical output, so we don't pick one as the "primary".
        let case_hash_input = job.workdir.join("case.toml");
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "Nextflow",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);
        let mut artefacts: Vec<Artifact> = Vec::new();

        // Surface the workdir itself as the primary `Native` artifact
        // — Nextflow's per-process output tree lives under
        // `<workdir>/work/`, alongside the user's pipeline outputs in
        // whatever `--outdir` they supplied. The viewer / power user
        // navigates from this anchor.
        if job.workdir.is_dir() {
            artefacts.push(Artifact {
                path: job.workdir.clone(),
                kind: ArtifactKind::Native,
                checksum: None,
                label: "Nextflow run workdir".to_string(),
            });
        }

        // Walk the workdir top-level for Nextflow's standard
        // observability outputs:
        //   - report.html   — execution report (process stats, CPU, RAM)
        //   - timeline.html — per-process Gantt chart
        //   - dag.svg       — workflow DAG render
        // These are off by default but commonly enabled via
        // `-with-report` / `-with-timeline` / `-with-dag` extras, so
        // surfacing them when present makes the post-run UX feel
        // first-class.
        let entries = match fs::read_dir(&job.workdir) {
            Ok(e) => e,
            Err(e) => {
                tracing::warn!(target: "valenx-nextflow", ?e, "workdir read failed");
                results.artifacts = artefacts;
                return Ok(results);
            }
        };
        for entry in entries.flatten() {
            let path = entry.path();
            if !path.is_file() {
                continue;
            }
            let name = match path.file_name().and_then(|s| s.to_str()) {
                Some(n) => n,
                None => continue,
            };
            let (kind, label) = match name {
                "report.html" => (ArtifactKind::Log, "Nextflow execution report".to_string()),
                "timeline.html" => (ArtifactKind::Log, "Nextflow timeline report".to_string()),
                "dag.svg" => (ArtifactKind::Native, "Nextflow workflow DAG".to_string()),
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
            ribbon_contributions: vec!["bio.nextflow.run"],
        }
    }
}

/// Returns true when the pipeline string looks like a filesystem path
/// rather than a registry slug. We treat anything containing a
/// path separator or ending in `.nf` as a path; everything else
/// (e.g. `nf-core/rnaseq`) goes through verbatim.
///
/// Note that `nf-core/rnaseq` *does* contain a forward slash, but
/// that's the registry's `<org>/<pipeline>` convention rather than a
/// filesystem path. To disambiguate we look at the leading character
/// — a path that starts with `/`, `\`, `.`, `~`, or a drive letter
/// (Windows `C:`) is a path; an `<org>/<pipeline>` slug starts with
/// an alphanumeric and contains no further structure suggesting
/// filesystem semantics. Using the `.nf` suffix as the second
/// indicator keeps relative-path local pipelines (e.g.
/// `subdir/main.nf`) on the path side of the fence.
fn looks_like_path(s: &str) -> bool {
    if s.starts_with('/') || s.starts_with('\\') || s.starts_with('.') || s.starts_with('~') {
        return true;
    }
    // Windows drive letter prefix like `C:\foo` or `C:/foo`.
    if s.len() >= 3
        && s.as_bytes()[0].is_ascii_alphabetic()
        && s.as_bytes()[1] == b':'
        && (s.as_bytes()[2] == b'\\' || s.as_bytes()[2] == b'/')
    {
        return true;
    }
    s.ends_with(".nf")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn info_is_bio_domain() {
        let info = NextflowAdapter::new().info();
        assert_eq!(info.id, "nextflow");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "Apache-2.0");
        assert_eq!(info.display_name, "Nextflow");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = NextflowAdapter::new().info();
        // Nextflow 23.10 is the LTS floor; 25.0 reserves the next
        // major bump.
        assert_eq!(info.version_range.min_inclusive, Version::new(23, 10, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(25, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = NextflowAdapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.nextflow.run"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = NextflowAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }

    #[test]
    fn prepare_rejects_relative_traversal_config() {
        // Round-8 RED→GREEN: relative `config` entries now route
        // through `confined_join` (absolute config paths still pass
        // through verbatim for the shared-config use-case). A
        // `../etc/passwd` traversal escape is rejected.
        use valenx_test_utils::tempdir;
        let d = tempdir("nextflow-config-traversal");
        // Use a registry-slug pipeline so the pipeline branch is not
        // exercised; the test targets `config` specifically.
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "nextflow.run"

[bio.nextflow]
pipeline = "nf-core/rnaseq"
config   = "../etc/passwd"
"#,
        )
        .unwrap();
        let case = Case {
            id: "nextflow-config-traversal".into(),
            path: d.clone(),
        };
        let workdir = d.join("workdir");
        let err = NextflowAdapter::new()
            .prepare(&case, &workdir)
            .unwrap_err();
        let msg = format!("{err}");
        assert!(
            msg.contains("absolute") || msg.contains("escape") || msg.contains("traversal") || msg.contains(".."),
            "expected confined_join rejection on relative-traversal config, got: {msg}"
        );
        let _ = std::fs::remove_dir_all(&d);
    }
}
