//! # valenx-adapter-minimap2
//!
//! Adapter for [minimap2](https://github.com/lh3/minimap2) — Heng
//! Li's "fast aligner for noisy long reads and assemblies." minimap2
//! is the workhorse mapper for Oxford Nanopore + PacBio reads, for
//! genome-vs-genome alignment, and for spliced RNA-seq mapping. It
//! also competes with BWA-MEM on short reads when invoked with the
//! `sr` preset.
//!
//! **Phase 18 — subprocess wrapper around `minimap2`.** Unlike BWA,
//! minimap2 builds its index on the fly from the reference, so there
//! is no separate index step. `prepare()` resolves the reference +
//! reads against the case directory, picks the preset (`map-ont` by
//! default), and composes a single `minimap2 -a -x <preset> -t N -o
//! out.sam <reference> <reads...>` invocation. `run()` streams that
//! through the shared subprocess runner; minimap2 prints chatty
//! progress to stderr ("[M::mm_idx_gen::...]") so the line handler
//! can lift those into UI progress hints.
//!
//! On `collect()` we surface the canonical `out.sam` plus any `.log`
//! files the user redirected stderr to, the same shape every other
//! Phase 18 aligner adapter publishes.

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

use crate::case_input::Minimap2Input;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(Minimap2Adapter::new())
}

pub struct Minimap2Adapter;

impl Minimap2Adapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for Minimap2Adapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "minimap2";
/// minimap2's binary candidates. The canonical name is `minimap2`
/// regardless of install path — Bioconda, Homebrew, and source
/// builds all produce the same name.
const BINARIES: &[&str] = &["minimap2"];

/// The aligned-reads filename written by `minimap2 -o`. Pinned so
/// `prepare()`, `collect()`, and the artifact label all agree.
const OUT_SAM: &str = "out.sam";

impl Adapter for Minimap2Adapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "minimap2",
            // 2.24 is the long-running stable line that every distro
            // ships and is the floor we test against; the upper bound
            // 3.0 reserves room for an eventual major bump that would
            // likely change CLI semantics enough to need a re-test.
            version_range: VersionRange {
                min_inclusive: Version::new(2, 24, 0),
                max_exclusive: Version::new(3, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            tool_license: "MIT",
            docs_url: "https://lh3.github.io/minimap2/minimap2.html",
            homepage_url: "https://github.com/lh3/minimap2",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(BINARIES) {
            Some(binary_path) => {
                // `minimap2 --version` prints just `2.x.y` on stdout
                // — exactly the shape `extract_semver` expects.
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
                hint: "minimap2 2.24+ required; install via `apt install minimap2`, \
                       `brew install minimap2`, or `conda install -c bioconda minimap2`"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = Minimap2Input::from_case_dir(&case.path)?;

        fs::create_dir_all(workdir)?;

        // Resolve the reference path against the case directory if
        // relative. Same convention as every other Phase 17/18 bio
        // adapter — `reference = "ref.fa"` next to `case.toml`.
        let source_reference = if input.reference.is_absolute() {
            input.reference.clone()
        } else {
            valenx_core::adapter_helpers::confined_join(
            &case.path,
            &input.reference,
        )?
        };
        if !source_reference.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.minimap2].reference `{}` not found (resolved {})",
                    input.reference.display(),
                    source_reference.display()
                ),
            });
        }

        // Resolve each read file against the case directory via
        // `confined_join` — rejects absolute paths and `..` traversal
        // out of the case sandbox (round-6 hardening). We *don't*
        // copy them into the workdir — long-read FASTQs run to many
        // GB and minimap2 reads them by path, so duplicating them is
        // pure waste.
        let mut resolved_reads: Vec<PathBuf> = Vec::with_capacity(input.reads.len());
        for read in &input.reads {
            let resolved = confined_join(&case.path, read)?;
            if !resolved.is_file() {
                return Err(AdapterError::InvalidCase {
                    case_path: case.path.join("case.toml"),
                    reason: format!(
                        "[bio.minimap2].reads entry `{}` not found (resolved {})",
                        read.display(),
                        resolved.display()
                    ),
                });
            }
            resolved_reads.push(resolved);
        }

        let binary_path = find_on_path(BINARIES).ok_or_else(|| AdapterError::ToolNotInstalled {
            name: INFO_ID,
            hint: "minimap2 2.24+ required; install via `apt install minimap2`, \
                       `brew install minimap2`, or `conda install -c bioconda minimap2`"
                .into(),
        })?;

        // Compose `minimap2 -a -x <preset> -t <N> -o out.sam
        //                  <reference> <reads...>`.
        // `-a` requests SAM output (the long-read default is PAF);
        // `-x <preset>` picks the scoring profile; `-o <file>` keeps
        // stdout clean so the subprocess runner can stream it
        // line-by-line if the user later swaps presets.
        let mut native_command: Vec<OsString> = vec![
            binary_path.into_os_string(),
            OsString::from("-a"),
            OsString::from("-x"),
            OsString::from(&input.preset),
            OsString::from("-t"),
            OsString::from(input.threads.to_string()),
            OsString::from("-o"),
            OsString::from(OUT_SAM),
        ];
        // Round-4 fix: extra_args after positionals — see
        // security/code-review.md.
        native_command.push(source_reference.into_os_string());
        for read in resolved_reads {
            native_command.push(read.into_os_string());
        }
        for arg in &input.extra_args {
            native_command.push(OsString::from(arg));
        }

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            // Whole-genome long-read alignment runs in tens of
            // minutes on a workstation; large all-vs-all overlaps
            // run for hours. 4 hours is a generous default that
            // covers the long tail without being absurd.
            estimated_runtime: Some(Duration::from_secs(4 * 60 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting minimap2", |line| {
            let mut hint = subprocess::Hint::default();
            // minimap2 emits "[M::mm_idx_gen::...]" lines while
            // building the on-the-fly index, then "[M::worker_pipeline::...]"
            // lines once mapping starts. Lift them to coarse progress
            // ticks so the UI shows forward motion.
            if line.contains("[M::main]") || line.contains("Real time:") {
                hint.progress = Some((95.0, line.to_string()));
            } else if line.contains("[M::worker_pipeline") || line.contains("[M::mm_idx_gen") {
                hint.progress = Some((50.0, line.to_string()));
            } else if line.contains("[E::") || line.contains("ERROR") {
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
        // Provenance: hash the staged out.sam if present; fall back
        // to case.toml when the alignment hasn't produced a SAM yet.
        let case_hash_input = {
            let sam = job.workdir.join(OUT_SAM);
            if sam.is_file() {
                sam
            } else {
                job.workdir.join("case.toml")
            }
        };
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "minimap2",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);
        let mut artefacts: Vec<Artifact> = Vec::new();

        // Walk the workdir top-level. minimap2 writes `out.sam` (and
        // optionally a `.paf` / `.bam` if the user passed extras).
        let entries = match fs::read_dir(&job.workdir) {
            Ok(e) => e,
            Err(e) => {
                tracing::warn!(target: "valenx-minimap2", ?e, "workdir read failed");
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
                // `out.sam` — the canonical aligned-reads file.
                Some("sam") => (
                    ArtifactKind::Native,
                    "minimap2 aligned reads (SAM)".to_string(),
                ),
                // BAM appears when users pipe through samtools.
                Some("bam") => (
                    ArtifactKind::Native,
                    "minimap2 aligned reads (BAM)".to_string(),
                ),
                // PAF — minimap2's native long-read alignment format,
                // emitted when users skip the `-a` flag.
                Some("paf") => (
                    ArtifactKind::Native,
                    "minimap2 aligned reads (PAF)".to_string(),
                ),
                Some("log") => (ArtifactKind::Log, "minimap2 log".to_string()),
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
        // The bio-specific Capability variants land in a follow-up
        // task; ribbon contributions are already enough for the
        // registry to surface the adapter.
        Capabilities {
            capabilities: Vec::new(),
            ribbon_contributions: vec!["bio.minimap2.align"],
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn info_is_bio_domain() {
        let info = Minimap2Adapter::new().info();
        assert_eq!(info.id, "minimap2");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "MIT");
        assert_eq!(info.display_name, "minimap2");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = Minimap2Adapter::new().info();
        // 2.24 is the floor we test against; 3.0 reserves room for
        // an eventual major bump.
        assert_eq!(info.version_range.min_inclusive, Version::new(2, 24, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(3, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = Minimap2Adapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.minimap2.align"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = Minimap2Adapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }
}
