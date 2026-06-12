//! # valenx-adapter-mrbayes
//!
//! Adapter for [MrBayes](http://nbisweden.github.io/MrBayes/) — the
//! long-standing Bayesian MCMC phylogenetic inference engine. MrBayes
//! is the historic workhorse for Bayesian phylogenetics: alongside
//! BEAST 2 it remains the de-facto choice for posterior tree sampling
//! across nucleotide / amino-acid / morphological datasets, with its
//! own NEXUS-embedded model-and-mcmc command language and built-in
//! Metropolis-coupled MCMC ("MC^3") swapping.
//!
//! **Phase 30.5 — subprocess wrapper around the `mb` binary.** The
//! user supplies a NEXUS file (DATA block + MrBayes block embedding
//! the model and `mcmc` command) via `[bio.mrbayes].nexus` in
//! `case.toml`, optionally a batch-mode flag. `prepare()` composes a
//! `mb [-i] <nexus> [extras...]` invocation; `run()` streams the run
//! via the shared subprocess runner.
//!
//! On `collect()` we walk the workdir for the canonical MrBayes
//! output families: `*.t` (sampled trees), `*.p` (parameter
//! samples), and `*.con.tre` (consensus tree).

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

use crate::case_input::MrBayesInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(MrBayesAdapter::new())
}

pub struct MrBayesAdapter;

impl MrBayesAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for MrBayesAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "mrbayes";
/// MrBayes' binary candidate. Both source and Bioconda installs use
/// the canonical short `mb` name (the project's own convention; the
/// long form `mrbayes` is not the canonical entry point).
const BINARIES: &[&str] = &["mb"];

impl Adapter for MrBayesAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "MrBayes",
            // MrBayes 3.2.x is the long-running stable line that
            // every distro ships; 3.2 (2012) is the floor we test
            // against and covers every release through 3.2.7. Upper
            // bound 4.0 reserves room for an eventual major bump.
            version_range: VersionRange {
                min_inclusive: Version::new(3, 2, 0),
                max_exclusive: Version::new(4, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            tool_license: "GPL-3.0",
            docs_url: "http://nbisweden.github.io/MrBayes/",
            homepage_url: "http://nbisweden.github.io/MrBayes/",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(BINARIES) {
            Some(binary_path) => {
                // `mb` (no args) prints the version banner on stdout
                // before dropping to its prompt; the generic detector
                // tries `--version` first (some packaged builds learn
                // it for free) and falls back to a bare-name scan.
                let found_version = detect_tool_version_semver(&binary_path, &["--version", ""]);
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
                hint: "MrBayes 3.2+ required; install via \
                       `conda install -c bioconda mrbayes`, \
                       `apt install mrbayes`, or build from source"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = MrBayesInput::from_case_dir(&case.path)?;

        fs::create_dir_all(workdir)?;

        // Resolve the NEXUS path against the case directory if
        // relative. Same convention as every other Phase 17/18 bio
        // adapter — `nexus = "data.nex"` next to `case.toml`.
        let source_nexus = if input.nexus.is_absolute() {
            input.nexus.clone()
        } else {
            valenx_core::adapter_helpers::confined_join(&case.path, &input.nexus)?
        };
        if !source_nexus.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.mrbayes].nexus `{}` not found (resolved {})",
                    input.nexus.display(),
                    source_nexus.display()
                ),
            });
        }

        let binary_path = find_on_path(BINARIES).ok_or_else(|| AdapterError::ToolNotInstalled {
            name: INFO_ID,
            hint: "MrBayes 3.2+ required; install via \
                       `conda install -c bioconda mrbayes`, \
                       `apt install mrbayes`, or build from source"
                .into(),
        })?;

        // Compose `mb [-i] <nexus> [extras...]`.
        // The NEXUS is the positional input — must come last so MrBayes
        // treats it as the model file rather than as the value of an
        // earlier flag.
        //
        // Round-3 fix: extras MUST come after the positional `<nexus>`.
        // Pre-fix they were appended between `-i` and the NEXUS, which
        // let a hostile case.toml slip an extra positional in via
        // `extra_args = ["phantom"]` and shift `<nexus>` onto a
        // different argument slot.
        let mut native_command: Vec<OsString> = vec![binary_path.into_os_string()];
        if input.batch {
            native_command.push(OsString::from("-i"));
        }
        native_command.push(source_nexus.into_os_string());
        for arg in &input.extra_args {
            native_command.push(OsString::from(arg));
        }

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            // Toy MCMC chains finish in seconds; production
            // posterior sampling on a serious dataset (multi-locus
            // partition, MC^3 across 4 chains) regularly runs for
            // days. 24 hours is a generous default that covers the
            // typical long tail without being absurd.
            estimated_runtime: Some(Duration::from_secs(24 * 60 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting MrBayes", |line| {
            let mut hint = subprocess::Hint::default();
            // MrBayes' progress chatter on stdout: a "MrBayes vX.Y.Z"
            // banner at startup, periodic "Generation NNNN" / "Avg
            // standard deviation of split frequencies" lines as the
            // MCMC progresses, and a "Analysis completed" sentinel
            // at end-of-run.
            if line.contains("Analysis completed") || line.contains("Continue with analysis") {
                hint.progress = Some((95.0, line.to_string()));
            } else if line.contains("Generation ") || line.contains("Avg standard deviation") {
                hint.progress = Some((50.0, line.to_string()));
            } else if line.contains("MrBayes v") || line.contains("Initializing") {
                hint.progress = Some((5.0, line.to_string()));
            }
            if line.contains("Error in command")
                || line.contains("ERROR")
                || line.contains("MrBayes terminated abnormally")
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
        // descriptor. MrBayes writes outputs alongside the input
        // NEXUS with `.run<N>.t` / `.run<N>.p` / `.con.tre`
        // filenames derived from the NEXUS basename, so we can't
        // assume a single fixed-name artifact for the prov hash.
        let case_hash_input = job.workdir.join("case.toml");
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "MrBayes",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);
        let mut artefacts: Vec<Artifact> = Vec::new();

        // Walk the workdir top-level. MrBayes writes outputs
        // alongside the input NEXUS:
        //   * `*.con.tre` — consensus tree (TreeAnnotator-style summary)
        //   * `*.t`        — sampled tree posterior (per chain)
        //   * `*.p`        — parameter samples (per chain)
        // The `.con.tre` check has to come before the bare `.tre` /
        // generic suffix walk because MrBayes nests two extensions.
        let entries = match fs::read_dir(&job.workdir) {
            Ok(e) => e,
            Err(e) => {
                tracing::warn!(target: "valenx-mrbayes", ?e, "workdir read failed");
                return Ok(results);
            }
        };
        for entry in entries.flatten() {
            let path = entry.path();
            if !path.is_file() {
                continue;
            }
            let name = path
                .file_name()
                .and_then(|s| s.to_str())
                .map(|s| s.to_ascii_lowercase());
            // Consensus tree — pick up first by exact suffix match.
            if let Some(n) = &name {
                if n.ends_with(".con.tre") {
                    artefacts.push(Artifact {
                        path,
                        kind: ArtifactKind::Native,
                        checksum: None,
                        label: "MrBayes consensus tree".to_string(),
                    });
                    continue;
                }
            }
            let ext = path
                .extension()
                .and_then(|s| s.to_str())
                .map(|s| s.to_ascii_lowercase());
            let (kind, label) = match ext.as_deref() {
                Some("t") => (ArtifactKind::Native, "MrBayes tree samples".to_string()),
                Some("p") => (
                    ArtifactKind::Tabular,
                    "MrBayes parameter samples".to_string(),
                ),
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
            ribbon_contributions: vec!["bio.mrbayes.mcmc"],
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn info_is_bio_domain() {
        let info = MrBayesAdapter::new().info();
        assert_eq!(info.id, "mrbayes");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "GPL-3.0");
        assert_eq!(info.display_name, "MrBayes");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = MrBayesAdapter::new().info();
        // MrBayes 3.2 (2012) is the floor; upper bound 4.0 reserves
        // room for the next major.
        assert_eq!(info.version_range.min_inclusive, Version::new(3, 2, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(4, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = MrBayesAdapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.mrbayes.mcmc"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = MrBayesAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }
}
