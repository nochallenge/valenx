//! # valenx-adapter-snakemake
//!
//! Adapter for [Snakemake](https://snakemake.github.io/) — Johannes
//! Köster's rule-based / data-flow workflow engine, the long-standing
//! peer to Nextflow in the bioinformatics-pipeline space. Snakemake
//! pipelines describe a DAG of rules over file targets — Snakemake
//! resolves which rules need rerunning by comparing input / output
//! mtimes. The Valenx adapter just runs `snakemake -s <Snakefile>`
//! and reports the workdir.
//!
//! **Phase 22 — pipeline-orchestrator wrapper.** Sister adapter to
//! Nextflow's: composes `snakemake -s <snakefile> --cores N` from
//! `[bio.snakemake]` in `case.toml`, with optional `--use-conda`
//! environment isolation, `-n` dry-run mode, an outboard
//! `--configfile`, and arbitrary positional rule / file targets.
//!
//! `collect()` surfaces the workdir as the primary `Native`
//! artifact; if `.snakemake/log/` exists (Snakemake's standard log
//! location), the most-recent log file is surfaced as a `Log`
//! artifact too.

#![forbid(unsafe_code)]
#![allow(missing_docs)]

pub mod case_input;

use std::ffi::OsString;
use std::fs;
use std::path::{Path, PathBuf};
use std::time::{Duration, SystemTime};

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

use crate::case_input::SnakemakeInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(SnakemakeAdapter::new())
}

pub struct SnakemakeAdapter;

impl SnakemakeAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for SnakemakeAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "snakemake";
/// Snakemake's binary candidates. Bioconda + pip both install the
/// canonical lowercase entry-point.
const BINARIES: &[&str] = &["snakemake"];

impl Adapter for SnakemakeAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "Snakemake",
            // 7.0 (Feb 2022) is where the modern `--cores` flag and
            // `--use-conda` semantics stabilised; 9.0 reserves room
            // for the next major bump.
            version_range: VersionRange {
                min_inclusive: Version::new(7, 0, 0),
                max_exclusive: Version::new(9, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            tool_license: "MIT",
            docs_url: "https://snakemake.readthedocs.io/",
            homepage_url: "https://snakemake.github.io/",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(BINARIES) {
            Some(binary_path) => {
                // `snakemake --version` prints the bare semver on
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
                hint:
                    "Snakemake 7.0+ required; install via \
                       `conda install -c bioconda snakemake`, \
                       `pip install snakemake`, or see \
                       https://snakemake.readthedocs.io/en/stable/getting_started/installation.html"
                        .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = SnakemakeInput::from_case_dir(&case.path)?;

        fs::create_dir_all(workdir)?;

        // Resolve Snakefile against case dir if relative. Same
        // convention as every other Phase 17/18 bio adapter.
        let source_snakefile = if input.snakefile.is_absolute() {
            input.snakefile.clone()
        } else {
            valenx_core::adapter_helpers::confined_join(
            &case.path,
            &input.snakefile,
        )?
        };
        if !source_snakefile.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.snakemake].snakefile `{}` not found (resolved {})",
                    input.snakefile.display(),
                    source_snakefile.display()
                ),
            });
        }

        // Resolve config-file path against case dir via
        // `confined_join` — sibling-field sweep round-8: `snakefile`
        // got it in round-6, `config_file` needs the same sandboxing.
        // Path validation runs *before* the binary probe so a hostile
        // case bundle gets a sandbox-rejection error even on hosts
        // where Snakemake isn't installed.
        let resolved_config_file = match &input.config_file {
            Some(config) => {
                Some(valenx_core::adapter_helpers::confined_join(&case.path, config)?)
            }
            None => None,
        };

        let binary_path = find_on_path(BINARIES).ok_or_else(|| AdapterError::ToolNotInstalled {
            name: INFO_ID,
            hint:
                "Snakemake 7.0+ required; install via \
                       `conda install -c bioconda snakemake`, \
                       `pip install snakemake`, or see \
                       https://snakemake.readthedocs.io/en/stable/getting_started/installation.html"
                    .into(),
        })?;

        // Compose `snakemake -s <snakefile> --cores N
        //                    [--use-conda] [-n]
        //                    [--configfile <path>]
        //                    [<targets>...] [extras...]`.
        let mut native_command: Vec<OsString> = vec![
            binary_path.into_os_string(),
            OsString::from("-s"),
            source_snakefile.into_os_string(),
            OsString::from("--cores"),
            OsString::from(input.cores.to_string()),
        ];

        if input.use_conda {
            native_command.push(OsString::from("--use-conda"));
        }

        if input.dry_run {
            // Snakemake accepts `-n` (short) or `--dry-run` (long);
            // pick the short form to mirror its docs / examples.
            native_command.push(OsString::from("-n"));
        }

        if let Some(config_path) = resolved_config_file {
            native_command.push(OsString::from("--configfile"));
            native_command.push(config_path.into_os_string());
        }

        // Targets are positional after the flags. Snakemake accepts
        // both rule names (`align`) and file targets
        // (`results/all.bam`); we just pass them through.
        for target in &input.targets {
            native_command.push(OsString::from(target));
        }

        for arg in &input.extra_args {
            native_command.push(OsString::from(arg));
        }

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            // Snakemake pipelines vary as widely as Nextflow's; cap
            // the estimate at 24 hours for the long tail.
            estimated_runtime: Some(Duration::from_secs(24 * 60 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting Snakemake", |line| {
            let mut hint = subprocess::Hint::default();
            // Snakemake's progress markers on stderr include
            // per-rule banners ("rule align:") and a closing
            // "Finished job N." / "(100%) done" status line. The
            // rule banner is a useful 50% tick; "(100%)" / "Complete
            // log:" pin the tail at 95%.
            if line.contains("Complete log:") || line.contains("(100%) done") {
                hint.progress = Some((95.0, line.to_string()));
            } else if line.contains("rule ") || line.contains("Job counts:") {
                hint.progress = Some((50.0, line.to_string()));
            } else if line.contains("Error in rule")
                || line.contains("MissingInputException")
                || line.contains("WorkflowError")
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
        // Provenance: hash the case.toml — like Nextflow, a
        // Snakemake workdir is a tree of per-rule outputs rather
        // than a single canonical file.
        let case_hash_input = job.workdir.join("case.toml");
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "Snakemake",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);
        let mut artefacts: Vec<Artifact> = Vec::new();

        // Surface the workdir itself as the primary `Native`
        // artifact — Snakemake's per-rule outputs land throughout
        // the user's directory tree (plus the `.snakemake/`
        // bookkeeping subdir). The viewer / power user navigates
        // from this anchor.
        if job.workdir.is_dir() {
            artefacts.push(Artifact {
                path: job.workdir.clone(),
                kind: ArtifactKind::Native,
                checksum: None,
                label: "Snakemake run workdir".to_string(),
            });
        }

        // Surface the most-recent log under `.snakemake/log/`. That
        // directory accumulates one log file per `snakemake` run
        // (timestamped), so picking the newest by mtime gets us
        // *this* run's log post-mortem. Skipping silently when the
        // directory doesn't exist keeps `collect()` robust against
        // partial / pre-init runs.
        let log_dir = job.workdir.join(".snakemake").join("log");
        if let Some(latest) = newest_file_in(&log_dir) {
            artefacts.push(Artifact {
                path: latest,
                kind: ArtifactKind::Log,
                checksum: None,
                label: "Snakemake run log".to_string(),
            });
        }

        artefacts.sort_by(|a, b| a.path.cmp(&b.path));
        results.artifacts = artefacts;
        Ok(results)
    }

    fn capabilities(&self) -> Capabilities {
        Capabilities {
            capabilities: Vec::new(),
            ribbon_contributions: vec!["bio.snakemake.run"],
        }
    }
}

/// Return the most-recently-modified regular file in `dir`, or
/// `None` if the directory doesn't exist / is empty / can't be
/// read. Used to pick the latest Snakemake log out of
/// `.snakemake/log/`.
fn newest_file_in(dir: &Path) -> Option<PathBuf> {
    let entries = fs::read_dir(dir).ok()?;
    let mut newest: Option<(SystemTime, PathBuf)> = None;
    for entry in entries.flatten() {
        let path = entry.path();
        if !path.is_file() {
            continue;
        }
        let mtime = match entry.metadata().and_then(|m| m.modified()) {
            Ok(t) => t,
            Err(_) => continue,
        };
        match &newest {
            Some((cur, _)) if *cur >= mtime => {}
            _ => newest = Some((mtime, path)),
        }
    }
    newest.map(|(_, p)| p)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn info_is_bio_domain() {
        let info = SnakemakeAdapter::new().info();
        assert_eq!(info.id, "snakemake");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "MIT");
        assert_eq!(info.display_name, "Snakemake");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = SnakemakeAdapter::new().info();
        // Snakemake 7.0 is the modern-flag floor; 9.0 reserves the
        // next major bump.
        assert_eq!(info.version_range.min_inclusive, Version::new(7, 0, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(9, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = SnakemakeAdapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.snakemake.run"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = SnakemakeAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }

    #[test]
    fn prepare_rejects_traversal_config_file_path() {
        // Round-8 RED→GREEN: `[bio.snakemake].config_file` now routes
        // through `confined_join` — sibling to `snakefile`, which got
        // it in round-6.
        //
        // We use `../etc/passwd` (relative traversal) for cross-platform
        // portability — see the bcftools test for the rationale.
        use valenx_test_utils::tempdir;
        let d = tempdir("snakemake-traversal");
        std::fs::write(d.join("Snakefile"), b"rule all:\n    input: []\n").unwrap();
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "snakemake.run"

[bio.snakemake]
snakefile   = "Snakefile"
config_file = "../etc/passwd"
cores       = 1
"#,
        )
        .unwrap();
        let case = Case {
            id: "snakemake-traversal".into(),
            path: d.clone(),
        };
        let workdir = d.join("workdir");
        let err = SnakemakeAdapter::new()
            .prepare(&case, &workdir)
            .unwrap_err();
        let msg = format!("{err}");
        assert!(
            msg.contains("absolute") || msg.contains("escape") || msg.contains("`..`") || msg.contains("traversal"),
            "expected confined_join rejection on config_file, got: {msg}"
        );
        let _ = std::fs::remove_dir_all(&d);
    }
}
