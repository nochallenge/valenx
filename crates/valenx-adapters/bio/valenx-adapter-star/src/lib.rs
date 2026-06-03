//! # valenx-adapter-star
//!
//! Adapter for [STAR](https://github.com/alexdobin/STAR) — Alex Dobin's
//! Spliced Transcripts Alignment to a Reference. STAR is the most
//! widely cited spliced RNA-seq aligner: a suffix-array seed-and-extend
//! design that handles both the genome alignment and the per-junction
//! splice graph in one pass, with throughput high enough to drive
//! GTEx-scale population RNA-seq.
//!
//! **Phase 18.6 — subprocess wrapper around `STAR`.** The user
//! supplies a genome index directory (or a FASTA + the directory
//! that should hold the index) plus 1 (single-end) or 2 (paired-end)
//! FASTQ files via `[bio.star]` in `case.toml`. `prepare()` runs
//! `STAR --runMode genomeGenerate` synchronously unless `skip_index
//! = true`, then composes the `STAR --runMode alignReads` invocation.
//! `run()` streams the alignment via the shared subprocess runner;
//! STAR prints periodic progress lines on stderr that the line
//! handler can lift to UI hints.
//!
//! On `collect()` we report the canonical alignment outputs:
//! `star_Aligned.out.bam` / `star_Aligned.out.sortedByCoord.out.bam`
//! / `star_Aligned.out.sam` as `Tabular` artifacts, plus
//! `star_Log.final.out` as a `Log`. STAR also leaves a flotilla of
//! intermediate `*_Log.out`, `*_SJ.out.tab`, and progress files in
//! the workdir; these are intentionally not surfaced — only the
//! three meaningful artifacts make it into the results.

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

use crate::case_input::StarInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(StarAdapter::new())
}

pub struct StarAdapter;

impl StarAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for StarAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "star";
/// STAR's binary candidate. Note the binary is **uppercase** `STAR`
/// — this is not a typo. The README, every release, and every
/// distro package use `STAR` as the on-PATH name; lowercase `star`
/// would collide with `glibc`'s `star` archiver on some systems and
/// the project chose the uppercase form to disambiguate.
const BINARIES: &[&str] = &["STAR"];

/// Filename prefix we tell `STAR --outFileNamePrefix` to use. STAR
/// emits a deterministic family of files all sharing this prefix
/// (`star_Aligned.*`, `star_Log.*`, `star_SJ.out.tab`, …); pinning
/// the prefix means `prepare()`, `run()`, and `collect()` all agree
/// on which files to look at.
const OUT_PREFIX: &str = "star_";

/// Translate the case-input `output_type` enum string into the one or
/// two argv tokens that follow `--outSAMtype` on the STAR command
/// line. STAR's flag takes "BAM" or "SAM" plus an optional sort
/// modifier — three combinations matter, and we surface them as a
/// single ergonomic enum so users don't have to know the exact
/// CLI shape.
fn samtype_args(s: &str) -> Vec<&'static str> {
    match s {
        "BAM_Unsorted" => vec!["BAM", "Unsorted"],
        "BAM_SortedByCoordinate" => vec!["BAM", "SortedByCoordinate"],
        "SAM" => vec!["SAM"],
        // Defensive default — case_input validation rejects
        // unknown values up front so this branch is unreachable in
        // practice, but returning a safe default keeps the function
        // total.
        _ => vec!["BAM", "SortedByCoordinate"],
    }
}

impl Adapter for StarAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "STAR",
            // STAR 2.7.x is the long-running stable line that every
            // distro ships; 3.0 reserves room for an eventual major
            // bump. `--runMode genomeGenerate` and the
            // `--outSAMtype BAM SortedByCoordinate` shape we depend
            // on have been stable since 2.7.0.
            version_range: VersionRange {
                min_inclusive: Version::new(2, 7, 0),
                max_exclusive: Version::new(3, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            tool_license: "MIT",
            docs_url: "https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf",
            homepage_url: "https://github.com/alexdobin/STAR",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(BINARIES) {
            Some(binary_path) => {
                // `STAR --version` prints a single line containing
                // the version on stdout (e.g. "2.7.11a"); the
                // combined scanner picks it up cleanly.
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
                hint: "STAR 2.7+ required; install via `apt install rna-star`, \
                       `brew install star`, or `conda install -c bioconda star`"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = StarInput::from_case_dir(&case.path)?;

        fs::create_dir_all(workdir)?;

        // Resolve genome_dir against the case directory. STAR
        // demands an existing directory at index-generation time
        // (the dir holds the index files); when we're building it
        // ourselves we create it first.
        // Round-9 hardening: relative `genome_dir` flows into STAR's
        // `--genomeDir` flag; wrap with `confined_join`.
        let resolved_genome_dir = if input.genome_dir.is_absolute() {
            input.genome_dir.clone()
        } else {
            confined_join(&case.path, &input.genome_dir)?
        };

        // Resolve the reference path against the case directory if
        // present and relative.
        // Round-9 hardening: wrap relative reference with `confined_join`.
        let resolved_reference: Option<PathBuf> = match &input.reference {
            Some(r) => {
                let p = if r.is_absolute() {
                    r.clone()
                } else {
                    confined_join(&case.path, r)?
                };
                if !p.is_file() {
                    return Err(AdapterError::InvalidCase {
                        case_path: case.path.join("case.toml"),
                        reason: format!(
                            "[bio.star].reference `{}` not found (resolved {})",
                            r.display(),
                            p.display()
                        ),
                    });
                }
                Some(p)
            }
            None => None,
        };

        // Resolve each read file against the case directory via
        // `confined_join` — rejects absolute paths and `..` traversal
        // out of the case sandbox (round-6 hardening). RNA-seq FASTQs
        // routinely run to tens of GB so we read them by path.
        let mut resolved_reads: Vec<PathBuf> = Vec::with_capacity(input.reads.len());
        for read in &input.reads {
            let resolved = confined_join(&case.path, read)?;
            if !resolved.is_file() {
                return Err(AdapterError::InvalidCase {
                    case_path: case.path.join("case.toml"),
                    reason: format!(
                        "[bio.star].reads entry `{}` not found (resolved {})",
                        read.display(),
                        resolved.display()
                    ),
                });
            }
            resolved_reads.push(resolved);
        }

        // Resolve the splice-junction GTF against the case directory.
        // Round-9 hardening: wrap relative sjdb_gtf with `confined_join`.
        let resolved_sjdb_gtf: Option<PathBuf> = match &input.sjdb_gtf {
            Some(g) => {
                let p = if g.is_absolute() {
                    g.clone()
                } else {
                    confined_join(&case.path, g)?
                };
                if !p.is_file() {
                    return Err(AdapterError::InvalidCase {
                        case_path: case.path.join("case.toml"),
                        reason: format!(
                            "[bio.star].sjdb_gtf `{}` not found (resolved {})",
                            g.display(),
                            p.display()
                        ),
                    });
                }
                Some(p)
            }
            None => None,
        };

        let binary_path = find_on_path(BINARIES).ok_or_else(|| AdapterError::ToolNotInstalled {
            name: INFO_ID,
            hint: "STAR 2.7+ required; install via `apt install rna-star`, \
                       `brew install star`, or `conda install -c bioconda star`"
                .into(),
        })?;

        // Build the genome index synchronously unless the user opted
        // out. STAR demands the genome dir already exist, so create
        // it first.
        if !input.skip_index {
            // case_input validation guarantees this is Some when
            // skip_index = false — defensively unwrap with a clear
            // error if something slipped through.
            let reference_for_index =
                resolved_reference
                    .as_ref()
                    .ok_or_else(|| AdapterError::InvalidCase {
                        case_path: case.path.join("case.toml"),
                        reason: "[bio.star].reference is required when skip_index = false".into(),
                    })?;

            fs::create_dir_all(&resolved_genome_dir)?;

            let mut index_cmd = std::process::Command::new(&binary_path);
            index_cmd
                .arg("--runMode")
                .arg("genomeGenerate")
                .arg("--genomeDir")
                .arg(&resolved_genome_dir)
                .arg("--genomeFastaFiles")
                .arg(reference_for_index)
                .arg("--runThreadN")
                .arg(input.threads.to_string());
            if let Some(gtf) = &resolved_sjdb_gtf {
                index_cmd
                    .arg("--sjdbGTFfile")
                    .arg(gtf)
                    .arg("--sjdbOverhang")
                    .arg("100");
            }
            let index_status = index_cmd
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
                        "`STAR --runMode genomeGenerate --genomeDir {}` failed (exit {}): {}",
                        resolved_genome_dir.display(),
                        out.status.code().unwrap_or(-1),
                        stderr.lines().next().unwrap_or("(no stderr)")
                    )));
                }
                Err(e) => {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "spawning `STAR --runMode genomeGenerate --genomeDir {}` failed: {e}",
                        resolved_genome_dir.display()
                    )));
                }
            }
        }

        // Compose the alignment invocation:
        //   STAR --runMode alignReads --genomeDir <genome_dir>
        //        --readFilesIn <reads...>
        //        --runThreadN N
        //        --outSAMtype <samtype...>
        //        --outFileNamePrefix star_
        //        [extras...]
        let mut native_command: Vec<OsString> = vec![
            binary_path.into_os_string(),
            OsString::from("--runMode"),
            OsString::from("alignReads"),
            OsString::from("--genomeDir"),
            resolved_genome_dir.into_os_string(),
            OsString::from("--readFilesIn"),
        ];
        // STAR's --readFilesIn takes the read files as separate argv
        // tokens (R1 [R2]), unlike most aligners that use distinct
        // -1 / -2 flags.
        for read in &resolved_reads {
            native_command.push(read.clone().into_os_string());
        }
        native_command.push(OsString::from("--runThreadN"));
        native_command.push(OsString::from(input.threads.to_string()));
        native_command.push(OsString::from("--outSAMtype"));
        for tok in samtype_args(&input.output_type) {
            native_command.push(OsString::from(tok));
        }
        native_command.push(OsString::from("--outFileNamePrefix"));
        native_command.push(OsString::from(OUT_PREFIX));
        for arg in &input.extra_args {
            native_command.push(OsString::from(arg));
        }

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            // RNA-seq runs span minutes (small library) to many hours
            // (full mammalian transcriptome on a single node); 4
            // hours mirrors BWA / Bowtie2 / HISAT2.
            estimated_runtime: Some(Duration::from_secs(4 * 60 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting STAR", |line| {
            let mut hint = subprocess::Hint::default();
            // STAR's progress chatter is on stderr — periodic
            // "Started job on" / "Started mapping on" / "Finished on"
            // markers wrap the run, with sparse "Number of input
            // reads" / "Uniquely mapped reads %" lines in between.
            // The "Finished on" line is the very last marker so we
            // pin it at 95%.
            if line.contains("Finished on") || line.contains("ALL DONE") {
                hint.progress = Some((95.0, line.to_string()));
            } else if line.contains("Started mapping") || line.contains("Number of input reads") {
                hint.progress = Some((50.0, line.to_string()));
            } else if line.contains("EXITING")
                || line.contains("FATAL ERROR")
                || line.contains("ERROR:")
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
        // Provenance: hash whichever alignment-output file is
        // present, falling back to case.toml so the provenance block
        // stays well-formed for partial / failed runs.
        let case_hash_input = {
            let bam_sorted = job
                .workdir
                .join(format!("{OUT_PREFIX}Aligned.sortedByCoord.out.bam"));
            let bam_unsorted = job.workdir.join(format!("{OUT_PREFIX}Aligned.out.bam"));
            let sam = job.workdir.join(format!("{OUT_PREFIX}Aligned.out.sam"));
            if bam_sorted.is_file() {
                bam_sorted
            } else if bam_unsorted.is_file() {
                bam_unsorted
            } else if sam.is_file() {
                sam
            } else {
                job.workdir.join("case.toml")
            }
        };
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "STAR",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);
        let mut artefacts: Vec<Artifact> = Vec::new();

        // Walk the workdir. STAR generates many intermediate files
        // (`star__STARtmp/`, `star_Log.out`, `star_Log.progress.out`,
        // `star_SJ.out.tab`, …); we only surface the three
        // meaningful ones per the output spec.
        let entries = match fs::read_dir(&job.workdir) {
            Ok(e) => e,
            Err(e) => {
                tracing::warn!(target: "valenx-star", ?e, "workdir read failed");
                return Ok(results);
            }
        };
        for entry in entries.flatten() {
            let path = entry.path();
            if !path.is_file() {
                continue;
            }
            let name = match path.file_name().and_then(|s| s.to_str()) {
                Some(n) => n.to_string(),
                None => continue,
            };

            // STAR's `--outSAMtype BAM` modes write
            // `<prefix>Aligned.out.bam` (Unsorted) or
            // `<prefix>Aligned.sortedByCoord.out.bam` (Sorted). The
            // SAM mode writes `<prefix>Aligned.out.sam`. We surface
            // all three as Tabular per the output spec — STAR's BAM
            // is BAM-format SAM, so the generic alignment-table
            // category fits even though the bytes are binary.
            let (kind, label) = if name == format!("{OUT_PREFIX}Aligned.sortedByCoord.out.bam") {
                (
                    ArtifactKind::Tabular,
                    "STAR aligned reads (BAM, sorted)".to_string(),
                )
            } else if name == format!("{OUT_PREFIX}Aligned.out.bam") {
                (
                    ArtifactKind::Tabular,
                    "STAR aligned reads (BAM)".to_string(),
                )
            } else if name == format!("{OUT_PREFIX}Aligned.out.sam") {
                (
                    ArtifactKind::Tabular,
                    "STAR aligned reads (SAM)".to_string(),
                )
            } else if name == format!("{OUT_PREFIX}Log.final.out") {
                (ArtifactKind::Log, "STAR alignment summary".to_string())
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
        Capabilities {
            capabilities: Vec::new(),
            ribbon_contributions: vec!["bio.star.align"],
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn info_is_bio_domain() {
        let info = StarAdapter::new().info();
        assert_eq!(info.id, "star");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "MIT");
        assert_eq!(info.display_name, "STAR");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = StarAdapter::new().info();
        // STAR 2.7.x is the de facto stable line; 3.0 reserves
        // room for an eventual major bump.
        assert_eq!(info.version_range.min_inclusive, Version::new(2, 7, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(3, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = StarAdapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.star.align"]);
    }

    #[test]
    fn samtype_args_maps_correctly() {
        // The three supported output_type values translate to the
        // expected `--outSAMtype` argv tokens. BAM modes split
        // into two args (format + sort), SAM is a single token.
        assert_eq!(samtype_args("BAM_Unsorted"), vec!["BAM", "Unsorted"]);
        assert_eq!(
            samtype_args("BAM_SortedByCoordinate"),
            vec!["BAM", "SortedByCoordinate"]
        );
        assert_eq!(samtype_args("SAM"), vec!["SAM"]);
        // Defensive default for unrecognised values — same as the
        // case_input default. case_input validation rejects unknown
        // values before they reach this function in practice.
        assert_eq!(samtype_args("nonsense"), vec!["BAM", "SortedByCoordinate"]);
    }

    /// Round-9 RED→GREEN: `[bio.star].genome_dir` used to be joined
    /// with bare `case.path.join`. Wrap with `confined_join`.
    #[test]
    fn prepare_rejects_genome_dir_traversing_outside_case_dir() {
        use valenx_test_utils::tempdir;
        let d = tempdir("star-genome-trav");
        std::fs::write(d.join("ref.fa"), b">chr1\nACGT\n").unwrap();
        std::fs::write(d.join("r1.fq"), b"@x\nACGT\n+\nIIII\n").unwrap();
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "star.align"

[bio.star]
genome_dir  = "../../etc"
reference   = "ref.fa"
reads       = ["r1.fq"]
threads     = 1
output_type = "BAM_SortedByCoordinate"
"#,
        )
        .unwrap();
        let case = Case {
            id: "star-genome-trav".into(),
            path: d.clone(),
        };
        let workdir = d.join("workdir");
        let err = StarAdapter::new().prepare(&case, &workdir).unwrap_err();
        let msg = format!("{err}");
        assert!(
            msg.contains("..") || msg.contains("stay within") || msg.contains("escape"),
            "expected confined_join rejection, got: {msg}"
        );
        let _ = std::fs::remove_dir_all(&d);
    }
}
