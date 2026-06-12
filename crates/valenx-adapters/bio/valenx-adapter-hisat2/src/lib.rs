//! # valenx-adapter-hisat2
//!
//! Adapter for [HISAT2](http://daehwankimlab.github.io/hisat2/) — Daehwan
//! Kim's graph-based splice-aware short-read aligner. HISAT2 is the
//! standard RNA-seq aligner for the modern Tuxedo / "new Tuxedo"
//! pipeline (HISAT2 → StringTie → Ballgown), built on a hierarchical
//! graph FM-index that handles both genome-wide alignment and the
//! per-transcript splice graph.
//!
//! **Phase 18.6 — subprocess wrapper around `hisat2`.** The user
//! supplies a reference FASTA plus 1 (single-end) or 2 (paired-end)
//! FASTQ files via `[bio.hisat2]` in `case.toml`. `prepare()` builds
//! the graph FM-index in the workdir (`hisat2-build`) unless
//! `skip_index = true`, then composes the `hisat2` invocation.
//! `run()` streams the alignment via the shared subprocess runner;
//! HISAT2 prints its summary stats on stderr at the end so the line
//! handler can lift the "overall alignment rate" line to a
//! near-completion progress hint.
//!
//! On `collect()` we report the canonical `out.sam` aligned-reads
//! file as a `Tabular` artifact.

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

use crate::case_input::Hisat2Input;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(Hisat2Adapter::new())
}

pub struct Hisat2Adapter;

impl Hisat2Adapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for Hisat2Adapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "hisat2";
/// HISAT2's binary candidates. `hisat2` is the canonical Linux /
/// macOS install name from Bioconda and source builds.
const BINARIES: &[&str] = &["hisat2"];

/// The aligned-reads filename we tell `hisat2 -S` to write. Pinned so
/// the `prepare()` invocation, the `collect()` walk, and the artifact
/// label all agree on what to look for.
const OUT_SAM: &str = "out.sam";

impl Adapter for Hisat2Adapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "HISAT2",
            // HISAT2 2.2.x is the long-running stable line that every
            // distro ships; 3.0 reserves room for an eventual major
            // bump.
            version_range: VersionRange {
                min_inclusive: Version::new(2, 2, 0),
                max_exclusive: Version::new(3, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            tool_license: "GPL-3.0",
            docs_url: "http://daehwankimlab.github.io/hisat2/manual/",
            homepage_url: "http://daehwankimlab.github.io/hisat2/",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(BINARIES) {
            Some(binary_path) => {
                // `hisat2 --version` prints the version on stdout;
                // the combined scanner picks it up cleanly.
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
                hint: "HISAT2 2.2+ required; install via `apt install hisat2`, \
                       `brew install hisat2`, or `conda install -c bioconda hisat2`"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = Hisat2Input::from_case_dir(&case.path)?;

        fs::create_dir_all(workdir)?;

        // Resolve the reference path against the case directory if
        // relative.
        let source_reference = if input.reference.is_absolute() {
            input.reference.clone()
        } else {
            valenx_core::adapter_helpers::confined_join(&case.path, &input.reference)?
        };
        if !source_reference.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.hisat2].reference `{}` not found (resolved {})",
                    input.reference.display(),
                    source_reference.display()
                ),
            });
        }

        // Resolve each read file against the case directory via
        // `confined_join` — rejects absolute paths and `..` traversal
        // out of the case sandbox (round-6 hardening). Same policy
        // as BWA / Bowtie2 — we *don't* copy them into the workdir
        // because RNA-seq FASTQs routinely hit tens of GB.
        let mut resolved_reads: Vec<PathBuf> = Vec::with_capacity(input.reads.len());
        for read in &input.reads {
            let resolved = confined_join(&case.path, read)?;
            if !resolved.is_file() {
                return Err(AdapterError::InvalidCase {
                    case_path: case.path.join("case.toml"),
                    reason: format!(
                        "[bio.hisat2].reads entry `{}` not found (resolved {})",
                        read.display(),
                        resolved.display()
                    ),
                });
            }
            resolved_reads.push(resolved);
        }

        let binary_path = find_on_path(BINARIES).ok_or_else(|| AdapterError::ToolNotInstalled {
            name: INFO_ID,
            hint: "HISAT2 2.2+ required; install via `apt install hisat2`, \
                       `brew install hisat2`, or `conda install -c bioconda hisat2`"
                .into(),
        })?;

        // The HISAT2 graph index basename is conventionally the
        // reference filename without its extension — `ref.fa` ->
        // `ref`. The index files themselves (`<base>.1.ht2` …
        // `<base>.8.ht2`) sit in the workdir (when we build) or
        // wherever the user pre-built one (when skip_index = true).
        let index_basename: String = source_reference
            .file_stem()
            .and_then(|s| s.to_str())
            .map(|s| s.to_string())
            .unwrap_or_else(|| "ref".to_string());

        // Build the graph FM-index in the workdir unless the user
        // opted out. We invoke `hisat2-build` synchronously here so
        // the subsequent `run()` is a single `hisat2` call, which
        // lets the shared subprocess runner stream stderr
        // line-by-line without chaining commands through a shell.
        if !input.skip_index {
            let index_status = std::process::Command::new("hisat2-build")
                .arg(&source_reference)
                .arg(&index_basename)
                .current_dir(workdir)
                .stdin(std::process::Stdio::null())
                .stdout(std::process::Stdio::piped())
                .stderr(std::process::Stdio::piped())
                .output();
            match index_status {
                Ok(out) if out.status.success() => {}
                Ok(out) => {
                    let stderr = String::from_utf8_lossy(&out.stderr);
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "`hisat2-build {} {}` failed (exit {}): {}",
                        source_reference.display(),
                        index_basename,
                        out.status.code().unwrap_or(-1),
                        stderr.lines().next().unwrap_or("(no stderr)")
                    )));
                }
                Err(e) => {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "spawning `hisat2-build {}` failed: {e}",
                        source_reference.display()
                    )));
                }
            }
        }

        // Compose the alignment invocation:
        //   hisat2 -x <index_base> -p <threads> -S out.sam
        //          [--rna-strandness <s> if not unstranded]
        //          (-U <reads[0]>) | (-1 <reads[0]> -2 <reads[1]>)
        //          [extras...]
        //
        // HISAT2 dispatches between single-end (`-U`) and paired-end
        // (`-1` / `-2`) based on which flag the reads come in on.
        let mut native_command: Vec<OsString> = vec![
            binary_path.into_os_string(),
            OsString::from("-x"),
            OsString::from(&index_basename),
            OsString::from("-p"),
            OsString::from(input.threads.to_string()),
            OsString::from("-S"),
            OsString::from(OUT_SAM),
        ];
        if input.strandness != "unstranded" {
            native_command.push(OsString::from("--rna-strandness"));
            native_command.push(OsString::from(&input.strandness));
        }
        if resolved_reads.len() == 1 {
            native_command.push(OsString::from("-U"));
            native_command.push(resolved_reads.remove(0).into_os_string());
        } else {
            native_command.push(OsString::from("-1"));
            native_command.push(resolved_reads[0].clone().into_os_string());
            native_command.push(OsString::from("-2"));
            native_command.push(resolved_reads[1].clone().into_os_string());
        }
        for arg in &input.extra_args {
            native_command.push(OsString::from(arg));
        }

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            // RNA-seq runs span minutes (small library) to many hours
            // (full mammalian transcriptome on a single node); 4
            // hours mirrors BWA / Bowtie2.
            estimated_runtime: Some(Duration::from_secs(4 * 60 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting HISAT2", |line| {
            let mut hint = subprocess::Hint::default();
            // HISAT2 emits its summary stats at the end of the run
            // on stderr — "N reads; of these:", "% overall alignment
            // rate", etc. The "overall alignment rate" line is the
            // very last useful marker so we pin it at 95%.
            if line.contains("overall alignment rate") {
                hint.progress = Some((95.0, line.to_string()));
            } else if line.contains("reads; of these:") {
                hint.progress = Some((50.0, line.to_string()));
            } else if line.contains("Error") || line.contains("(ERR):") {
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
        // to case.toml when the alignment hasn't produced a SAM yet
        // so the provenance block stays well-formed for partial /
        // failed runs.
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
            "HISAT2",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);
        let mut artefacts: Vec<Artifact> = Vec::new();

        // Walk the workdir top-level. HISAT2 writes `out.sam` in the
        // workdir alongside the graph FM-index files (`<base>.*.ht2`).
        // We only surface the SAM and any log file.
        let entries = match fs::read_dir(&job.workdir) {
            Ok(e) => e,
            Err(e) => {
                tracing::warn!(target: "valenx-hisat2", ?e, "workdir read failed");
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
                // `out.sam` — the aligned-reads file. Tabular per
                // output spec.
                Some("sam") => (ArtifactKind::Tabular, "HISAT2 aligned reads".to_string()),
                Some("log") => (ArtifactKind::Log, "HISAT2 log".to_string()),
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
            ribbon_contributions: vec!["bio.hisat2.align"],
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn info_is_bio_domain() {
        let info = Hisat2Adapter::new().info();
        assert_eq!(info.id, "hisat2");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "GPL-3.0");
        assert_eq!(info.display_name, "HISAT2");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = Hisat2Adapter::new().info();
        // 2.2.x is the de facto stable line; 3.0 reserves room for
        // an eventual major bump.
        assert_eq!(info.version_range.min_inclusive, Version::new(2, 2, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(3, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = Hisat2Adapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.hisat2.align"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = Hisat2Adapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }
}
