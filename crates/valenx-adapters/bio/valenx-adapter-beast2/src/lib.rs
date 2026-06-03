//! # valenx-adapter-beast2
//!
//! Adapter for [BEAST 2](https://beast2.org/) — the cross-platform
//! Bayesian Evolutionary Analysis by Sampling Trees engine. BEAST 2 is
//! the canonical Bayesian MCMC framework for time-calibrated
//! phylogenetics: tip-dated trees, relaxed molecular clocks, coalescent
//! demographic models, birth-death speciation models, and the
//! ever-growing universe of BEAST 2 packages (BDSKY, MASCOT, BEASTling,
//! StarBEAST3, ...). It complements the maximum-likelihood Phase 30
//! tools (IQ-TREE, RAxML-NG, FastTree) with a full posterior over
//! tree topologies and parameters.
//!
//! **Phase 30.5 — subprocess wrapper around the `beast` binary.** The
//! user supplies a BEAUti-generated `.xml` model via
//! `[bio.beast2].xml` in `case.toml`, optionally a seed, a thread
//! count, and an overwrite flag. `prepare()` composes a
//! `beast [-seed N] -threads N [-overwrite] <xml> [extras...]`
//! invocation; `run()` streams the run via the shared subprocess
//! runner.
//!
//! On `collect()` we walk the workdir for `*.log` (the trace log
//! Tracer reads) and `*.trees` (the sampled tree posterior
//! TreeAnnotator / DensiTree consumes).

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

use crate::case_input::Beast2Input;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(Beast2Adapter::new())
}

pub struct Beast2Adapter;

impl Beast2Adapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for Beast2Adapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "beast2";
/// BEAST 2's binary candidate. Both source and Bioconda installs
/// expose the canonical lowercase `beast` launcher script.
const BINARIES: &[&str] = &["beast"];

impl Adapter for Beast2Adapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "BEAST 2",
            // BEAST 2's modern stable line is the 2.7.x series (2022+);
            // 2.7 introduced the modern threading + package manager.
            // Upper bound 3.0 reserves room for an eventual major bump.
            version_range: VersionRange {
                min_inclusive: Version::new(2, 7, 0),
                max_exclusive: Version::new(3, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            tool_license: "LGPL-2.1",
            docs_url: "https://beast2.org/",
            homepage_url: "https://beast2.org/",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(BINARIES) {
            Some(binary_path) => {
                // `beast -version` (one dash) prints a banner with
                // the BEAST 2 release on stdout. The generic detector
                // tries both the conventional `--version` and BEAST's
                // own `-version` form.
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
                hint: "BEAST 2.7+ required; install via \
                       `conda install -c bioconda beast2`, the standalone \
                       installer from https://beast2.org/, or build from \
                       source"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = Beast2Input::from_case_dir(&case.path)?;

        fs::create_dir_all(workdir)?;

        // Resolve the .xml model against the case directory if
        // relative. Same convention as every other Phase 17/18 bio
        // adapter — `xml = "model.xml"` next to `case.toml`.
        let source_xml = if input.xml.is_absolute() {
            input.xml.clone()
        } else {
            valenx_core::adapter_helpers::confined_join(
            &case.path,
            &input.xml,
        )?
        };
        if !source_xml.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.beast2].xml `{}` not found (resolved {})",
                    input.xml.display(),
                    source_xml.display()
                ),
            });
        }

        let binary_path = find_on_path(BINARIES).ok_or_else(|| AdapterError::ToolNotInstalled {
            name: INFO_ID,
            hint: "BEAST 2.7+ required; install via \
                       `conda install -c bioconda beast2`, the standalone \
                       installer from https://beast2.org/, or build from \
                       source"
                .into(),
        })?;

        // Compose `beast [-seed N] -threads N [-overwrite] <xml> [extras...]`.
        // The XML is the positional model file — must come last so
        // BEAST treats it as the model rather than as the value of an
        // earlier flag.
        let mut native_command: Vec<OsString> = vec![binary_path.into_os_string()];
        if let Some(seed) = input.seed {
            native_command.push(OsString::from("-seed"));
            native_command.push(OsString::from(seed.to_string()));
        }
        native_command.push(OsString::from("-threads"));
        native_command.push(OsString::from(input.threads.to_string()));
        if input.overwrite {
            native_command.push(OsString::from("-overwrite"));
        }
        // Round-4 fix: extra_args after positionals — see
        // security/code-review.md. Otherwise a hostile case.toml could
        // inject options that shadow the positional source_xml path.
        native_command.push(source_xml.into_os_string());
        for arg in &input.extra_args {
            native_command.push(OsString::from(arg));
        }

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            // Toy MCMC chains finish in seconds; production
            // posterior sampling on a serious dataset (multi-locus,
            // tip-dated, relaxed clock) regularly runs for days.
            // 24 hours is a generous default that covers the typical
            // long tail without being absurd.
            estimated_runtime: Some(Duration::from_secs(24 * 60 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting BEAST 2", |line| {
            let mut hint = subprocess::Hint::default();
            // BEAST 2's progress chatter on stdout: a "Random number
            // seed" banner at startup, periodic "Sample" / "Chain"
            // status lines as the MCMC progresses, and a "End
            // likelihood" / "Total calculation time" sentinel at
            // end-of-run.
            if line.contains("End likelihood") || line.contains("Total calculation time") {
                hint.progress = Some((95.0, line.to_string()));
            } else if line.contains("Sample") || line.contains("posterior") {
                hint.progress = Some((50.0, line.to_string()));
            } else if line.contains("Random number seed") || line.contains("BEAST v2") {
                hint.progress = Some((5.0, line.to_string()));
            }
            if line.contains("Error:") || line.contains("ERROR") || line.contains("java.lang.") {
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
        // descriptor. BEAST's output filenames are baked into the
        // XML at <log fileName="..."> sites, so we can't assume a
        // single fixed-name artifact for the prov hash.
        let case_hash_input = job.workdir.join("case.toml");
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "BEAST 2",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);
        let mut artefacts: Vec<Artifact> = Vec::new();

        // Walk the workdir top-level. BEAST's two canonical output
        // families are `.log` (parameter trace, consumed by Tracer)
        // and `.trees` (sampled tree posterior, consumed by
        // TreeAnnotator / DensiTree). Both filenames are baked into
        // the XML model so we accept any file with those extensions.
        let entries = match fs::read_dir(&job.workdir) {
            Ok(e) => e,
            Err(e) => {
                tracing::warn!(target: "valenx-beast2", ?e, "workdir read failed");
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
            let (kind, label) = match ext.as_deref() {
                Some("log") => (ArtifactKind::Log, "BEAST 2 trace log".to_string()),
                Some("trees") => (ArtifactKind::Native, "BEAST 2 sampled trees".to_string()),
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
            ribbon_contributions: vec!["bio.beast2.mcmc"],
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn info_is_bio_domain() {
        let info = Beast2Adapter::new().info();
        assert_eq!(info.id, "beast2");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "LGPL-2.1");
        assert_eq!(info.display_name, "BEAST 2");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = Beast2Adapter::new().info();
        // BEAST 2.7+ is the modern stable line; upper bound 3.0
        // reserves room for the next major.
        assert_eq!(info.version_range.min_inclusive, Version::new(2, 7, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(3, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = Beast2Adapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.beast2.mcmc"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = Beast2Adapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }
}
