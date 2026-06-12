//! # valenx-adapter-simrna
//!
//! Adapter for [SimRNA](https://genesilico.pl/SimRNAweb/) — the
//! Bujnicki-lab Monte Carlo engine for predicting three-dimensional
//! RNA tertiary structure from primary sequence. SimRNA samples a
//! coarse-grained RNA model under a knowledge-based statistical
//! potential, optionally using replica-exchange Monte Carlo to span
//! a temperature ladder; the engine writes candidate PDB models,
//! replica trajectories (`*.trafl`), and an energy log into the
//! working directory.
//!
//! **Phase 45 — subprocess wrapper around the SimRNA CLI.** The user
//! supplies the config and sequence files via `[bio.simrna]` in
//! `case.toml`; `prepare()` composes
//! `SimRNA -c <config> -s <sequence> -o <basename> -R <n_replicas> [extras...]`,
//! `run()` streams progress via the shared subprocess runner, and
//! `collect()` walks the workdir for `<basename>*.pdb`,
//! `<basename>*.trafl`, `<basename>*.txt`, and `*.log` outputs.

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

use crate::case_input::SimRnaInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(SimRnaAdapter::new())
}

pub struct SimRnaAdapter;

impl SimRnaAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for SimRnaAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "simrna";
/// SimRNA's binary candidates. Upstream tarballs ship the binary
/// capitalised (`SimRNA`); some package managers normalise to the
/// lowercase `simrna`. Probing both covers either install style.
const BINARIES: &[&str] = &["SimRNA", "simrna"];

impl Adapter for SimRnaAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "SimRNA",
            // SimRNA 3.20 is the modern stable line shipped from
            // genesilico.pl; upper bound 4.0 reserves room for the
            // next major bump.
            version_range: VersionRange {
                min_inclusive: Version::new(3, 20, 0),
                max_exclusive: Version::new(4, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            tool_license: "GPL-3.0",
            docs_url: "https://genesilico.pl/software/stand-alone/simrna",
            homepage_url: "https://genesilico.pl/SimRNAweb/",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(BINARIES) {
            Some(binary_path) => {
                // SimRNA prints its version on `--version` (and the
                // shorter `-v` on older builds). Try both.
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
                hint: "SimRNA 3.20+ required; download from \
                       https://genesilico.pl/software/stand-alone/simrna \
                       and place the binary on PATH"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = SimRnaInput::from_case_dir(&case.path)?;

        // Round-4 security: reject `output_basename = "../etc/passwd"`
        // and friends before the value flows into any path join.
        // Same pattern as the round-3 fix in bionetgen/iqtree/art/fasttree.
        valenx_core::adapter_helpers::validate_output_basename(
            &input.output_basename,
            "[bio.simrna].output_basename",
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

        // Resolve config + sequence against the case directory if
        // relative. Both are read in place by SimRNA — no staging,
        // no `confined_join`. Plain `case.path.join` is correct.
        let source_config = if input.config.is_absolute() {
            input.config.clone()
        } else {
            valenx_core::adapter_helpers::confined_join(&case.path, &input.config)?
        };
        if !source_config.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.simrna].config `{}` not found (resolved {})",
                    input.config.display(),
                    source_config.display()
                ),
            });
        }

        let source_sequence = if input.sequence.is_absolute() {
            input.sequence.clone()
        } else {
            valenx_core::adapter_helpers::confined_join(&case.path, &input.sequence)?
        };
        if !source_sequence.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.simrna].sequence `{}` not found (resolved {})",
                    input.sequence.display(),
                    source_sequence.display()
                ),
            });
        }

        let binary_path = find_on_path(BINARIES).ok_or_else(|| AdapterError::ToolNotInstalled {
            name: INFO_ID,
            hint: "SimRNA 3.20+ required; download from \
                       https://genesilico.pl/software/stand-alone/simrna \
                       and place the binary on PATH"
                .into(),
        })?;

        // Compose `SimRNA -c <config> -s <sequence> -o <basename>
        //          -R <n_replicas> [extras...]`. SimRNA writes its
        // outputs into the cwd (workdir) under the supplied basename.
        let mut native_command: Vec<OsString> = vec![
            binary_path.into_os_string(),
            OsString::from("-c"),
            source_config.into_os_string(),
            OsString::from("-s"),
            source_sequence.into_os_string(),
            OsString::from("-o"),
            OsString::from(&input.output_basename),
            OsString::from("-R"),
            OsString::from(input.n_replicas.to_string()),
        ];
        for arg in &input.extra_args {
            native_command.push(OsString::from(arg));
        }

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            // Short RNAs (~30 nt) finish in minutes; long-loop /
            // riboswitch sampling with replica exchange routinely
            // runs for many hours. 4 hours is the canonical adapter
            // ceiling.
            estimated_runtime: Some(Duration::from_secs(4 * 60 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting SimRNA", |line| {
            let mut hint = subprocess::Hint::default();
            // SimRNA emits a startup banner ("SimRNA"), per-step MC
            // energy lines, and a "finished" sentinel near the end.
            if line.contains("Simulation finished") || line.contains("simulation finished") {
                hint.progress = Some((95.0, line.to_string()));
            } else if line.starts_with("SimRNA") {
                hint.progress = Some((5.0, line.to_string()));
            } else if line.contains("step ") || line.contains("energy") {
                hint.progress = Some((50.0, line.to_string()));
            } else if line.contains("ERROR") || line.contains("Error:") {
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
        // descriptor. SimRNA reads config + sequence in place, so
        // case.toml is the closest stable input fingerprint.
        let case_hash_input = job.workdir.join("case.toml");
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "SimRNA",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);
        let mut artefacts: Vec<Artifact> = Vec::new();

        // Read the output basename back from case.toml so we can
        // prefix-filter the canonical output families.
        let basename = case_input::SimRnaInput::from_case_dir(&job.workdir)
            .ok()
            .map(|i| i.output_basename);

        // Walk the workdir top-level. SimRNA writes
        // `<basename>*.pdb` (tertiary structures),
        // `<basename>*.trafl` (replica trajectories),
        // `<basename>*.txt` (energy log) plus optional `*.log` files.
        let entries = match fs::read_dir(&job.workdir) {
            Ok(e) => e,
            Err(e) => {
                tracing::warn!(target: "valenx-simrna", ?e, "workdir read failed");
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
                .file_stem()
                .and_then(|s| s.to_str())
                .unwrap_or("")
                .to_string();
            let stem_matches_basename = match basename.as_deref() {
                Some(b) => stem.starts_with(b),
                None => true,
            };
            match ext.as_deref() {
                Some("pdb") => {
                    if !stem_matches_basename {
                        continue;
                    }
                    artefacts.push(Artifact {
                        path,
                        kind: ArtifactKind::Native,
                        checksum: None,
                        label: "SimRNA tertiary structure".to_string(),
                    });
                }
                Some("trafl") => {
                    if !stem_matches_basename {
                        continue;
                    }
                    artefacts.push(Artifact {
                        path,
                        kind: ArtifactKind::Native,
                        checksum: None,
                        label: "SimRNA trajectory".to_string(),
                    });
                }
                Some("txt") => {
                    if !stem_matches_basename {
                        continue;
                    }
                    artefacts.push(Artifact {
                        path,
                        kind: ArtifactKind::Tabular,
                        checksum: None,
                        label: "SimRNA energy log".to_string(),
                    });
                }
                Some("log") => {
                    artefacts.push(Artifact {
                        path,
                        kind: ArtifactKind::Log,
                        checksum: None,
                        label: "SimRNA log".to_string(),
                    });
                }
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
            ribbon_contributions: vec!["bio.simrna.fold"],
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn info_is_bio_domain() {
        let info = SimRnaAdapter::new().info();
        assert_eq!(info.id, "simrna");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "GPL-3.0");
        assert_eq!(info.display_name, "SimRNA");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = SimRnaAdapter::new().info();
        // SimRNA 3.20 is the modern stable line; 4.0 reserves room
        // for the next major bump.
        assert_eq!(info.version_range.min_inclusive, Version::new(3, 20, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(4, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = SimRnaAdapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.simrna.fold"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = SimRnaAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }
}
