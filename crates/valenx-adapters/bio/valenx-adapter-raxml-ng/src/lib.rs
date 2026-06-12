//! # valenx-adapter-raxml-ng
//!
//! Adapter for [RAxML-NG](https://github.com/amkozlov/raxml-ng) —
//! the next-generation rewrite of the classic RAxML / ExaML
//! phylogenetics codebase. Uses the same likelihood machinery
//! (PThreads + SSE / AVX) but with a streamlined CLI: every output
//! is namespaced with a user-chosen prefix (`<prefix>.raxml.bestTree`,
//! `<prefix>.raxml.support`, `<prefix>.raxml.log`), and the modes
//! map cleanly onto run intent (`--search`, `--bootstrap`, `--all`).
//!
//! **Phase 30 — subprocess wrapper around `raxml-ng`.** The user
//! supplies an alignment plus a substitution model via
//! `[bio.raxml-ng]` in `case.toml`. `prepare()` resolves the input
//! against the case directory and composes the invocation:
//! `raxml-ng --<mode> --msa <alignment> --model <model>
//! --threads N --prefix <prefix> [--bs-trees <bootstrap>] [extras...]`.
//! `run()` streams via the shared subprocess runner. `collect()`
//! walks for the canonical `<prefix>.raxml.bestTree` (Newick ML
//! tree), `<prefix>.raxml.support` (bootstrap support tree), and
//! `<prefix>.raxml.log`.

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

use crate::case_input::RaxmlNgInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(RaxmlNgAdapter::new())
}

pub struct RaxmlNgAdapter;

impl RaxmlNgAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for RaxmlNgAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "raxml-ng";
/// RAxML-NG ships under exactly one binary name on every distro
/// (`raxml-ng`); the legacy hyphen-rich variants (`raxmlHPC`,
/// `raxmlHPC-PTHREADS-SSE3`, etc.) belong to the original RAxML
/// codebase, which is a separate adapter.
const BINARIES: &[&str] = &["raxml-ng"];

impl Adapter for RaxmlNgAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "RAxML-NG",
            // RAxML-NG 1.2.x is the current stable line. Floor at
            // 1.2 covers the `--bs-trees` / `--prefix` CLI surface
            // the adapter relies on; upper bound 2.0 reserves room
            // for an eventual major bump.
            version_range: VersionRange {
                min_inclusive: Version::new(1, 2, 0),
                max_exclusive: Version::new(2, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            tool_license: "AGPL-3.0",
            docs_url: "https://github.com/amkozlov/raxml-ng/wiki",
            homepage_url: "https://github.com/amkozlov/raxml-ng",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(BINARIES) {
            Some(binary_path) => {
                // `raxml-ng -v` and `raxml-ng --version` both print a
                // banner with "RAxML-NG v. 1.2.0 ..." on stdout.
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
                hint: "RAxML-NG 1.2+ required; install via `conda install -c bioconda raxml-ng` \
                       or build from https://github.com/amkozlov/raxml-ng"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = RaxmlNgInput::from_case_dir(&case.path)?;

        // Round-4 security: reject `prefix = "../etc/passwd"`
        // and friends before the value flows into any path join.
        // Same pattern as the round-3 fix in bionetgen/iqtree/art/fasttree.
        valenx_core::adapter_helpers::validate_output_basename(
            &input.prefix,
            "[bio.raxml-ng].prefix",
        )
        .map_err(|e| AdapterError::InvalidCase {
            case_path: case.path.join("case.toml"),
            reason: format!("{e}"),
        })?;

        fs::create_dir_all(workdir)?;

        // Resolve alignment against the case directory if relative.
        let source_alignment = if input.alignment.is_absolute() {
            input.alignment.clone()
        } else {
            valenx_core::adapter_helpers::confined_join(&case.path, &input.alignment)?
        };
        if !source_alignment.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.raxml-ng].alignment `{}` not found (resolved {})",
                    input.alignment.display(),
                    source_alignment.display()
                ),
            });
        }

        let binary_path = find_on_path(BINARIES).ok_or_else(|| AdapterError::ToolNotInstalled {
            name: INFO_ID,
            hint: "RAxML-NG 1.2+ required; install via `conda install -c bioconda raxml-ng` \
                       or build from https://github.com/amkozlov/raxml-ng"
                .into(),
        })?;

        // Compose `raxml-ng --<mode> --msa <alignment> --model <model>
        // --threads N --prefix <prefix> [--bs-trees N] [extras...]`.
        // Mode comes first so subsequent flags are unambiguous; the
        // mode-flag form (`--search`, `--all`, `--bootstrap`) is what
        // RAxML-NG documents as canonical.
        let mut native_command: Vec<OsString> = vec![
            binary_path.into_os_string(),
            OsString::from(format!("--{}", input.mode)),
            OsString::from("--msa"),
            source_alignment.into_os_string(),
            OsString::from("--model"),
            OsString::from(&input.model),
            OsString::from("--threads"),
            OsString::from(input.threads.to_string()),
            OsString::from("--prefix"),
            OsString::from(&input.prefix),
        ];
        // Bootstrap replicates only matter for `all` and `bootstrap`
        // modes; case-input validation already enforces the count is
        // > 0 in those modes, so this is a straight pass-through.
        if input.mode == "all" || input.mode == "bootstrap" {
            native_command.push(OsString::from("--bs-trees"));
            native_command.push(OsString::from(input.bootstrap.to_string()));
        }
        for arg in &input.extra_args {
            native_command.push(OsString::from(arg));
        }

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            // Search-only on a small alignment finishes in minutes;
            // all-mode on a thousand-taxon alignment with hundreds
            // of bootstrap replicates can run for a day. 8 hours is
            // a generous default that still fails fast on stuck runs.
            estimated_runtime: Some(Duration::from_secs(8 * 60 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting RAxML-NG", |line| {
            let mut hint = subprocess::Hint::default();
            // RAxML-NG writes a chatty stdout log: an "Analysis
            // started" banner, "Initial LogLikelihood" markers per
            // tree-search start, "Bootstrap tree" lines per replicate,
            // and a "Final LogLikelihood" when search converges.
            if line.contains("Analysis started") || line.contains("Loading binary MSA") {
                hint.progress = Some((10.0, line.to_string()));
            } else if line.contains("Initial LogLikelihood") || line.contains("ML tree search") {
                hint.progress = Some((40.0, line.to_string()));
            } else if line.contains("Bootstrap tree #") {
                hint.progress = Some((75.0, line.to_string()));
            } else if line.contains("Final LogLikelihood") || line.contains("Analysis completed") {
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
        // Provenance: hash the bestTree if present (the canonical
        // run output). Falls back to case.toml otherwise.
        let case_hash_input = {
            let mut tree_path: Option<PathBuf> = None;
            if let Ok(entries) = fs::read_dir(&job.workdir) {
                for entry in entries.flatten() {
                    let p = entry.path();
                    if let Some(name) = p.file_name().and_then(|s| s.to_str()) {
                        if name.ends_with(".raxml.bestTree") {
                            tree_path = Some(p);
                            break;
                        }
                    }
                }
            }
            tree_path.unwrap_or_else(|| job.workdir.join("case.toml"))
        };
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "RAxML-NG",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);
        let mut artefacts: Vec<Artifact> = Vec::new();

        // Walk the workdir top-level. RAxML-NG writes outputs with
        // the user-chosen prefix and a `.raxml.<kind>` suffix, so
        // matching by extension alone won't cut it; instead, match by
        // the trailing `.raxml.<kind>` substring of the file name.
        let entries = match fs::read_dir(&job.workdir) {
            Ok(e) => e,
            Err(e) => {
                tracing::warn!(target: "valenx-raxml-ng", ?e, "workdir read failed");
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
            // RAxML-NG writes many files: bestTree, support, log,
            // bestModel, mlTrees, bootstraps, rba (binary MSA), etc.
            // Surface the three the user is most likely to need.
            let (kind, label) = if name.ends_with(".raxml.bestTree") {
                (ArtifactKind::Native, "RAxML-NG ML tree".to_string())
            } else if name.ends_with(".raxml.support") {
                (
                    ArtifactKind::Native,
                    "RAxML-NG bootstrap support".to_string(),
                )
            } else if name.ends_with(".raxml.log") {
                (ArtifactKind::Log, "RAxML-NG log".to_string())
            } else {
                continue;
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
            ribbon_contributions: vec!["bio.raxml-ng.tree"],
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn info_is_bio_domain() {
        let info = RaxmlNgAdapter::new().info();
        assert_eq!(info.id, "raxml-ng");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "AGPL-3.0");
        assert_eq!(info.display_name, "RAxML-NG");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = RaxmlNgAdapter::new().info();
        // 1.2.x is the modern stable line; 2.0 reserves room for an
        // eventual major bump.
        assert_eq!(info.version_range.min_inclusive, Version::new(1, 2, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(2, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = RaxmlNgAdapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.raxml-ng.tree"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = RaxmlNgAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }
}
