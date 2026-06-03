//! # valenx-adapter-slim
//!
//! Adapter for [SLiM](https://messerlab.org/slim/) — Philipp Messer's
//! forward-time population-genetics simulator. SLiM evolves a
//! finite-population model generation by generation under a
//! user-defined Eidos script: mutation rates, selection
//! coefficients, recombination maps, demographic events,
//! migrations, mating systems. The state is sampled at any
//! generation the script asks for, and tree-sequence recording
//! (the `treeSeqOutput()` family) feeds straight into
//! tskit / msprime downstream.
//!
//! **Phase 29 — subprocess wrapper around the `slim` binary.** The
//! user supplies a `.slim` Eidos script via `[bio.slim].script` in
//! `case.toml`, optionally a seed and Eidos-constant overrides.
//! `prepare()` composes a `slim [-s <seed>] <script> [extras...]`
//! invocation; `run()` streams the run via the shared subprocess
//! runner.
//!
//! Output paths are determined by the script (typically via
//! `treeSeqOutput("sim.trees")` or `writeFile(...)` calls), so the
//! adapter doesn't try to predict them — `collect()` walks the
//! workdir for any `.trees` (tree-sequence) and `.log` (run log)
//! the script left behind.

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

use crate::case_input::SlimInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(SlimAdapter::new())
}

pub struct SlimAdapter;

impl SlimAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for SlimAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "slim";
/// SLiM's binary candidate. Conda-forge / source / Homebrew all
/// install under the canonical lowercase `slim` name.
const BINARIES: &[&str] = &["slim"];

impl Adapter for SlimAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "SLiM",
            // SLiM's modern release line is the 4.x series (2022+);
            // the 4.0 release introduced the streamlined Eidos
            // type system and the `treeSeqOutput()` helpers we
            // rely on for tskit interop. Upper bound 5.0 reserves
            // room for an eventual major bump.
            version_range: VersionRange {
                min_inclusive: Version::new(4, 0, 0),
                max_exclusive: Version::new(5, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            tool_license: "GPL-3.0",
            docs_url: "https://messerlab.org/slim/",
            homepage_url: "https://messerlab.org/slim/",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(BINARIES) {
            Some(binary_path) => {
                // `slim -version` (one dash) prints a banner with
                // "SLiM version X.Y.Z" on stdout. The generic
                // detector tries both the conventional `--version`
                // and the SLiM-native `-version` form.
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
                hint: "SLiM 4.0+ required; install via \
                       `conda install -c conda-forge slim`, `brew install \
                       messerlab/slim/slim`, or build from \
                       https://messerlab.org/slim/"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = SlimInput::from_case_dir(&case.path)?;

        // Round-4 security: reject `output_basename = "../etc/passwd"`
        // and friends before the value flows into any path join.
        // Same pattern as the round-3 fix in bionetgen/iqtree/art/fasttree.
        valenx_core::adapter_helpers::validate_output_basename(
            &input.output_basename,
            "[bio.slim].output_basename",
        )
        .map_err(|e| AdapterError::InvalidCase {
            case_path: case.path.join("case.toml"),
            reason: format!("{e}"),
        })?;

        fs::create_dir_all(workdir)?;

        // Resolve the .slim script against the case directory if
        // relative. Same convention as every other Phase 17/18 bio
        // adapter — `script = "model.slim"` next to `case.toml`.
        let source_script = if input.script.is_absolute() {
            input.script.clone()
        } else {
            valenx_core::adapter_helpers::confined_join(
            &case.path,
            &input.script,
        )?
        };
        if !source_script.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.slim].script `{}` not found (resolved {})",
                    input.script.display(),
                    source_script.display()
                ),
            });
        }

        let binary_path = find_on_path(BINARIES).ok_or_else(|| AdapterError::ToolNotInstalled {
            name: INFO_ID,
            hint: "SLiM 4.0+ required; install via \
                       `conda install -c conda-forge slim`, `brew install \
                       messerlab/slim/slim`, or build from \
                       https://messerlab.org/slim/"
                .into(),
        })?;

        // Compose `slim [-s <seed>] <script> [extras...]`.
        // SLiM's `-s` flag pins the PRNG seed; without it SLiM
        // picks its own seed and prints it on the run banner.
        let mut native_command: Vec<OsString> = vec![binary_path.into_os_string()];
        if let Some(seed) = input.seed {
            native_command.push(OsString::from("-s"));
            native_command.push(OsString::from(seed.to_string()));
        }
        // The script is positional — must come last so SLiM treats
        // it as the model file rather than as the value of an
        // earlier flag. Round-4 fix: extra_args go AFTER positionals
        // (was before) — see security/code-review.md.
        native_command.push(source_script.into_os_string());
        for arg in &input.extra_args {
            native_command.push(OsString::from(arg));
        }

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            // Tiny test models finish in seconds; whole-genome
            // forward simulations with selection sweeps run for
            // hours. 8 hours is a generous default that covers
            // the typical long tail.
            estimated_runtime: Some(Duration::from_secs(8 * 60 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting SLiM", |line| {
            let mut hint = subprocess::Hint::default();
            // SLiM's progress chatter on stdout: a "// Initial
            // random seed" banner at startup, periodic
            // "// generation N" lines as the simulation
            // progresses, and a "// Run finished" sentinel at
            // end-of-run.
            if line.contains("Run finished") || line.contains("// Run finished") {
                hint.progress = Some((95.0, line.to_string()));
            } else if line.contains("// generation") {
                hint.progress = Some((50.0, line.to_string()));
            } else if line.contains("Initial random seed") || line.contains("// Initial") {
                hint.progress = Some((5.0, line.to_string()));
            }
            if line.contains("ERROR") || line.contains("Eidos error") {
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
        // descriptor. SLiM scripts choose their own output paths
        // so we can't assume a fixed-name artifact for the prov
        // hash.
        let case_hash_input = job.workdir.join("case.toml");
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "SLiM",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);
        let mut artefacts: Vec<Artifact> = Vec::new();

        // Walk the workdir top-level. SLiM scripts conventionally
        // write outputs to the working directory; tree sequences
        // (`.trees`) feed tskit / msprime, and ad-hoc per-run logs
        // (`.log`) capture stochastic-sample summaries.
        let entries = match fs::read_dir(&job.workdir) {
            Ok(e) => e,
            Err(e) => {
                tracing::warn!(target: "valenx-slim", ?e, "workdir read failed");
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
                Some("trees") => (ArtifactKind::Native, "SLiM tree sequence".to_string()),
                Some("log") => (ArtifactKind::Log, "SLiM run log".to_string()),
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
            ribbon_contributions: vec!["bio.slim.simulate"],
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn info_is_bio_domain() {
        let info = SlimAdapter::new().info();
        assert_eq!(info.id, "slim");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "GPL-3.0");
        assert_eq!(info.display_name, "SLiM");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = SlimAdapter::new().info();
        // SLiM 4.x is the modern stable line (2022+); upper bound
        // 5.0 reserves room for the next major.
        assert_eq!(info.version_range.min_inclusive, Version::new(4, 0, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(5, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = SlimAdapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.slim.simulate"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = SlimAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }
}
