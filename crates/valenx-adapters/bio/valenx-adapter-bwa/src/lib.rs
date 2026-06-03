//! # valenx-adapter-bwa
//!
//! Adapter for [BWA](https://github.com/lh3/bwa) — Heng Li's
//! Burrows-Wheeler Aligner. The de-facto standard for mapping
//! Illumina-style short reads (50–250 bp) against a reference genome,
//! shipped as the `bwa` binary on every bioinformatics workstation.
//! BWA-MEM is the workhorse algorithm: seed-and-extend with maximal
//! exact matches, fast enough to drive whole-genome resequencing
//! pipelines.
//!
//! **Phase 18 — subprocess wrapper around `bwa mem`.** The user
//! supplies a reference FASTA plus 1 (single-end) or 2 (paired-end)
//! FASTQ files via `[bio.bwa]` in `case.toml`. `prepare()` builds the
//! BWT index next to the reference (`bwa index`) unless
//! `skip_index = true`, then composes the `bwa mem` invocation.
//! `run()` streams the alignment via the shared subprocess runner;
//! BWA prints progress to stderr line-by-line, so the line handler
//! can lift "Processed N reads" markers to progress hints.
//!
//! On `collect()` we look for the canonical `out.sam` aligned-reads
//! file and any auxiliary BWA outputs in the workdir.

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

use crate::case_input::BwaInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(BwaAdapter::new())
}

pub struct BwaAdapter;

impl BwaAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for BwaAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "bwa";
/// BWA's binary candidates. `bwa` is the canonical Linux / macOS
/// install name from Bioconda, Homebrew, and source builds.
const BINARIES: &[&str] = &["bwa"];

/// The aligned-reads filename written by `bwa mem -o`. Pinned so the
/// `prepare()` invocation, the `collect()` walk, and the artifact
/// label all agree on what to look for.
const OUT_SAM: &str = "out.sam";

impl Adapter for BwaAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "BWA",
            // BWA 0.7.x is the long-running stable line that every
            // distro ships; `bwa mem -o <file>` (used to keep the
            // subprocess runner streaming stderr cleanly without
            // having to redirect stdout to a file) landed in 0.7.16.
            // The upper bound 1.0 reserves room for an eventual major
            // bump.
            version_range: VersionRange {
                min_inclusive: Version::new(0, 7, 0),
                max_exclusive: Version::new(1, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            tool_license: "GPL-3.0",
            docs_url: "https://bio-bwa.sourceforge.net/bwa.shtml",
            homepage_url: "https://github.com/lh3/bwa",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(BINARIES) {
            Some(binary_path) => {
                // `bwa` (no args) prints a banner including the
                // version on stderr; our generic flag-based detector
                // reads stderr too, so a bare-name invocation works.
                // Try the conventional flag first regardless — older
                // BWAs ignore unknown flags and still print the
                // banner, while newer versions might learn `--version`
                // for free.
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
                hint: "BWA 0.7+ required; install via `apt install bwa`, \
                       `brew install bwa`, or `conda install -c bioconda bwa`"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = BwaInput::from_case_dir(&case.path)?;

        fs::create_dir_all(workdir)?;

        // Resolve the reference path against the case directory via
        // `confined_join` — rejects absolute paths and `..` traversal
        // out of the case sandbox. Same threat model as every other
        // Phase 17/18 bio adapter: a shared-case bundle should not be
        // able to point `reference` at `/etc/passwd`.
        let source_reference = confined_join(&case.path, &input.reference)?;
        if !source_reference.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.bwa].reference `{}` not found (resolved {})",
                    input.reference.display(),
                    source_reference.display()
                ),
            });
        }

        // Resolve each read file against the case directory via
        // `confined_join`. We *don't* copy them into the workdir —
        // short-read FASTQs routinely run to tens of GB, and BWA reads
        // them by path so duplicating them is pure waste. The sandbox
        // check still applies: a shared bundle can't point a read at
        // `/etc/passwd`.
        let mut resolved_reads: Vec<PathBuf> = Vec::with_capacity(input.reads.len());
        for read in &input.reads {
            let resolved = confined_join(&case.path, read)?;
            if !resolved.is_file() {
                return Err(AdapterError::InvalidCase {
                    case_path: case.path.join("case.toml"),
                    reason: format!(
                        "[bio.bwa].reads entry `{}` not found (resolved {})",
                        read.display(),
                        resolved.display()
                    ),
                });
            }
            resolved_reads.push(resolved);
        }

        let binary_path = find_on_path(BINARIES).ok_or_else(|| AdapterError::ToolNotInstalled {
            name: INFO_ID,
            hint: "BWA 0.7+ required; install via `apt install bwa`, \
                       `brew install bwa`, or `conda install -c bioconda bwa`"
                .into(),
        })?;

        // Build the BWT index next to the reference unless the user
        // opted out. Index files (`.bwt`, `.pac`, `.ann`, `.amb`,
        // `.sa`) are written next to `<reference>.fa` so successive
        // runs can opt out by setting `skip_index = true`. Running
        // the indexer here (in prepare) keeps the subsequent `run()`
        // call as a single `bwa mem` invocation, which lets the
        // subprocess runner stream stderr line-by-line without
        // chaining commands through a shell.
        if !input.skip_index {
            let index_status = std::process::Command::new(&binary_path)
                .arg("index")
                .arg(&source_reference)
                .stdin(std::process::Stdio::null())
                .stdout(std::process::Stdio::piped())
                .stderr(std::process::Stdio::piped())
                .output();
            match index_status {
                Ok(out) if out.status.success() => {}
                Ok(out) => {
                    let stderr = String::from_utf8_lossy(&out.stderr);
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "`bwa index {}` failed (exit {}): {}",
                        source_reference.display(),
                        out.status.code().unwrap_or(-1),
                        stderr.lines().next().unwrap_or("(no stderr)")
                    )));
                }
                Err(e) => {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "spawning `bwa index {}` failed: {e}",
                        source_reference.display()
                    )));
                }
            }
        }

        // Compose `bwa mem -t <N> -o out.sam <reference> <reads...>`.
        // `-o <file>` (BWA >= 0.7.16) avoids stdout redirection so
        // the shared subprocess runner can keep streaming stderr,
        // which is where `bwa mem` writes its progress chatter.
        let mut native_command: Vec<OsString> = vec![
            binary_path.into_os_string(),
            OsString::from("mem"),
            OsString::from("-t"),
            OsString::from(input.threads.to_string()),
            OsString::from("-o"),
            OsString::from(OUT_SAM),
        ];
        // Round-4 fix: extra_args go AFTER positionals — see
        // security/code-review.md. Pushing them before the reference
        // would let a hostile case.toml inject options that swallow
        // the positionals (e.g. `extra_args = ["-x", "pacbio"]`
        // would be fine, but `extra_args = ["-h"]` would shadow the
        // reference path with `-h`).
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
            // Whole-exome runs finish in minutes; whole-genome runs
            // run for hours on a single node. 4 hours is a generous
            // default that covers the long tail without being absurd.
            estimated_runtime: Some(Duration::from_secs(4 * 60 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting BWA", |line| {
            let mut hint = subprocess::Hint::default();
            // `bwa mem` emits "[M::process] read N sequences (X bp)..."
            // lines on stderr roughly every batch; lift to a 50%
            // progress tick so the UI shows forward motion. The
            // sentinel "Real time:" line marks the end-of-run timing
            // summary — pin it at 95%.
            if line.contains("Real time:") {
                hint.progress = Some((95.0, line.to_string()));
            } else if line.contains("[M::process]") || line.contains("read sequences") {
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
        // Provenance: hash the staged out.sam if present (the canonical
        // run output). Falls back to case.toml when the alignment hasn't
        // produced a SAM yet — keeps the provenance block well-formed
        // for partial / failed runs.
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
            "BWA",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);
        let mut artefacts: Vec<Artifact> = Vec::new();

        // Walk the workdir top-level. BWA only writes `out.sam`; the
        // BWT index sits next to the reference (outside the workdir),
        // so the only artifact to surface from a typical run is the
        // SAM. We also pick up any `.log` the user redirected stderr
        // to, in case future cases configure that.
        let entries = match fs::read_dir(&job.workdir) {
            Ok(e) => e,
            Err(e) => {
                tracing::warn!(target: "valenx-bwa", ?e, "workdir read failed");
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
                // `out.sam` — the aligned-reads file. SAM is plain
                // text; downstream tools (samtools view -bS) convert
                // to BAM, but the aligner output stays SAM.
                Some("sam") => (ArtifactKind::Native, "BWA aligned reads (SAM)".to_string()),
                // BWA can also emit BAM if the user piped through
                // samtools; surface as Native too.
                Some("bam") => (ArtifactKind::Native, "BWA aligned reads (BAM)".to_string()),
                Some("log") => (ArtifactKind::Log, "BWA log".to_string()),
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
            ribbon_contributions: vec!["bio.bwa.align"],
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn info_is_bio_domain() {
        let info = BwaAdapter::new().info();
        assert_eq!(info.id, "bwa");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "GPL-3.0");
        assert_eq!(info.display_name, "BWA");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = BwaAdapter::new().info();
        // BWA 0.7.x is the de facto stable line; 1.0 reserves room
        // for an eventual major bump.
        assert_eq!(info.version_range.min_inclusive, Version::new(0, 7, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(1, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = BwaAdapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.bwa.align"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = BwaAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }
}
