//! # valenx-adapter-iqtree
//!
//! Adapter for [IQ-TREE](http://www.iqtree.org/) — the de-facto
//! ML-phylogenetics workhorse. IQ-TREE searches for the maximum-
//! likelihood phylogeny under a chosen substitution model, with
//! ultrafast bootstrap (UFBoot) for branch support and ModelFinder
//! Plus (`MFP`) for automatic model selection from a catalogue of
//! 88 nucleotide and 22 amino-acid models. Single-binary CLI;
//! every Bioconda / Homebrew install ships either `iqtree` (1.x
//! legacy line) or `iqtree2` (2.x current line).
//!
//! **Phase 30 — subprocess wrapper around `iqtree2`.** The user
//! supplies a multi-FASTA / PHYLIP / NEXUS alignment via
//! `[bio.iqtree]` in `case.toml`. `prepare()` resolves the input
//! against the case directory and composes the IQ-TREE invocation:
//! `-s <alignment> -m <model> [-B <bootstrap>] -T <threads>
//! --prefix <prefix> [extras...]`. `run()` streams via the shared
//! subprocess runner — IQ-TREE writes a chatty progress log on
//! stdout that we lift to UI ticks. `collect()` walks for the
//! canonical `<prefix>.treefile` (Newick ML tree),
//! `<prefix>.iqtree` (model summary), and `<prefix>.log`.

#![forbid(unsafe_code)]
#![allow(missing_docs)]

pub mod case_input;

use std::ffi::OsString;
use std::fs;
use std::path::{Path, PathBuf};
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

use crate::case_input::IqTreeInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(IqTreeAdapter::new())
}

pub struct IqTreeAdapter;

impl IqTreeAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for IqTreeAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "iqtree";
/// IQ-TREE binary candidates. Newer 2.x ships as `iqtree2` on every
/// distro that supports parallel installs of the 1.x line; some
/// minimal images install only `iqtree`. Probe for both, prefer
/// `iqtree2` (the current generation).
const BINARIES: &[&str] = &["iqtree2", "iqtree"];

impl Adapter for IqTreeAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "IQ-TREE",
            // IQ-TREE 2.2.x is the current stable line (2.2.0 landed
            // in 2022). Floor at 2.2 covers UFBoot2 behaviour the
            // adapter relies on; upper bound 3.0 reserves room for
            // an eventual major bump.
            version_range: VersionRange {
                min_inclusive: Version::new(2, 2, 0),
                max_exclusive: Version::new(3, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            tool_license: "GPL-2.0",
            docs_url: "http://www.iqtree.org/doc/",
            homepage_url: "http://www.iqtree.org/",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(BINARIES) {
            Some(binary_path) => {
                // `iqtree2 --version` prints "IQ-TREE multicore version
                // 2.2.0 ..." on stdout; the generic detector picks up
                // the leading SemVer.
                let found_version = detect_tool_version_semver(&binary_path, &["--version"]);
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
                hint: "IQ-TREE 2.2+ required; install via `apt install iqtree`, \
                       `brew install iqtree`, or `conda install -c bioconda iqtree`"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = IqTreeInput::from_case_dir(&case.path)?;

        // Round-3 security fix: `prefix` is fed to `--prefix` and
        // becomes the basename of every output file. A hostile value
        // like `"../../etc/cron.d/x"` would otherwise let IQ-TREE
        // write outside the workdir.
        valenx_core::adapter_helpers::validate_output_basename(
            &input.prefix,
            "[bio.iqtree].prefix",
        )
        .map_err(|e| AdapterError::InvalidCase {
            case_path: case.path.join("case.toml"),
            reason: format!("{e}"),
        })?;

        fs::create_dir_all(workdir)?;

        // Resolve the alignment path against the case directory if
        // relative. Same convention as every other Phase 17/18 bio
        // adapter — `alignment = "aln.fa"` next to `case.toml`.
        let source_alignment = if input.alignment.is_absolute() {
            input.alignment.clone()
        } else {
            valenx_core::adapter_helpers::confined_join(&case.path, &input.alignment)?
        };
        if !source_alignment.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.iqtree].alignment `{}` not found (resolved {})",
                    input.alignment.display(),
                    source_alignment.display()
                ),
            });
        }

        let binary_path = find_on_path(BINARIES).ok_or_else(|| AdapterError::ToolNotInstalled {
            name: INFO_ID,
            hint: "IQ-TREE 2.2+ required; install via `apt install iqtree`, \
                       `brew install iqtree`, or `conda install -c bioconda iqtree`"
                .into(),
        })?;

        // Compose the IQ-TREE invocation. `-s` for the alignment,
        // `-m` for the substitution model, `-B` for UFBoot replicates
        // (only when > 0 — passing 0 errors out), `-T` for threads,
        // `--prefix` to namespace every output file, then user extras.
        // IQ-TREE writes outputs into the working directory with the
        // chosen prefix as the basename.
        let mut native_command: Vec<OsString> = vec![
            binary_path.into_os_string(),
            OsString::from("-s"),
            source_alignment.into_os_string(),
            OsString::from("-m"),
            OsString::from(&input.model),
        ];
        if input.bootstrap > 0 {
            native_command.push(OsString::from("-B"));
            native_command.push(OsString::from(input.bootstrap.to_string()));
        }
        native_command.push(OsString::from("-T"));
        native_command.push(OsString::from(&input.threads));
        native_command.push(OsString::from("--prefix"));
        native_command.push(OsString::from(&input.prefix));
        for arg in &input.extra_args {
            native_command.push(OsString::from(arg));
        }

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            // Small alignments (~50 sequences, ~1 kbp) finish in
            // minutes; large ones (~thousands of sequences with
            // model-finder + UFBoot) can run for many hours. 6 hours
            // covers the long tail without being absurd.
            estimated_runtime: Some(Duration::from_secs(6 * 60 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting IQ-TREE", |line| {
            let mut hint = subprocess::Hint::default();
            // IQ-TREE writes a chatty stdout log: model-test progress,
            // tree-search iterations, bootstrap replicates. Lift the
            // major sentinels to UI ticks so progress is visible.
            if line.contains("Performing ModelFinder") || line.contains("ModelFinder will test") {
                hint.progress = Some((20.0, line.to_string()));
            } else if line.contains("Tree search") || line.contains("BEST SCORE FOUND") {
                hint.progress = Some((60.0, line.to_string()));
            } else if line.contains("UFBoot trees printed") || line.contains("CONSENSUS TREE") {
                hint.progress = Some((85.0, line.to_string()));
            } else if line.contains("Total CPU time used") || line.contains("Date and Time:") {
                hint.progress = Some((95.0, line.to_string()));
            } else if line.contains("ERROR") || line.contains("WARNING:") {
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
        // Provenance: hash the staged treefile if present (the
        // canonical run output). Falls back to case.toml when the
        // run hasn't produced a tree yet — keeps the provenance
        // block well-formed for partial / failed runs.
        let case_hash_input = {
            // Try to find any *.treefile in the workdir.
            let mut tree_path: Option<PathBuf> = None;
            if let Ok(entries) = fs::read_dir(&job.workdir) {
                for entry in entries.flatten() {
                    let p = entry.path();
                    if p.extension().and_then(|s| s.to_str()) == Some("treefile") {
                        tree_path = Some(p);
                        break;
                    }
                }
            }
            tree_path.unwrap_or_else(|| job.workdir.join("case.toml"))
        };
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "IQ-TREE",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);
        let mut artefacts: Vec<Artifact> = Vec::new();

        // Walk the workdir top-level. IQ-TREE writes many files keyed
        // by prefix: `.treefile` is the ML tree (Newick), `.iqtree`
        // is the human-readable summary, `.log` is the run log,
        // `.contree` is the consensus tree (UFBoot), `.mldist` is the
        // ML distance matrix, etc. We surface the three the user is
        // most likely to need.
        let entries = match fs::read_dir(&job.workdir) {
            Ok(e) => e,
            Err(e) => {
                tracing::warn!(target: "valenx-iqtree", ?e, "workdir read failed");
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
                // `<prefix>.treefile` — the ML tree in Newick.
                Some("treefile") => (ArtifactKind::Native, "IQ-TREE ML tree".to_string()),
                // `<prefix>.iqtree` — the human-readable run summary
                // with model parameters, likelihood, AIC/BIC.
                Some("iqtree") => (ArtifactKind::Log, "IQ-TREE summary".to_string()),
                // `<prefix>.log` — the verbose run log.
                Some("log") => (ArtifactKind::Log, "IQ-TREE log".to_string()),
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
        // Bio-specific Capability variants land in a follow-up task;
        // ribbon contributions are already enough for the registry
        // to surface the adapter.
        Capabilities {
            capabilities: Vec::new(),
            ribbon_contributions: vec!["bio.iqtree.tree"],
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn info_is_bio_domain() {
        let info = IqTreeAdapter::new().info();
        assert_eq!(info.id, "iqtree");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "GPL-2.0");
        assert_eq!(info.display_name, "IQ-TREE");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = IqTreeAdapter::new().info();
        // 2.2.x is the modern stable line; 3.0 reserves room for an
        // eventual major bump.
        assert_eq!(info.version_range.min_inclusive, Version::new(2, 2, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(3, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = IqTreeAdapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.iqtree.tree"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = IqTreeAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }
}
