//! # valenx-adapter-gatk
//!
//! Adapter for the [Broad Institute GATK](https://gatk.broadinstitute.org/)
//! HaplotypeCaller — the canonical Java-based germline variant caller
//! and the second half (with [DeepVariant]) of the modern short-read
//! variant-calling toolkit. GATK is the reference implementation for
//! the GATK Best Practices pipelines used in clinical and research
//! resequencing.
//!
//! [DeepVariant]: ../valenx_adapter_deepvariant/index.html
//!
//! **Phase 19 — subprocess wrapper around `gatk HaplotypeCaller`.**
//! The user supplies a reference FASTA, a sorted/indexed BAM, an
//! output VCF path, and a JVM heap size via `[bio.gatk]` in
//! `case.toml`. `prepare()` composes the
//! `gatk --java-options -Xmx<heap> HaplotypeCaller …` invocation.
//! `run()` streams via the shared subprocess runner; GATK prints
//! `INFO/WARN/ERROR` markers to stderr line-by-line, so the line
//! handler can lift "ProgressMeter" markers to progress hints.
//!
//! ## How `--java-options` is passed
//!
//! GATK accepts two equivalent shapes:
//!
//! - `gatk --java-options -Xmx8g HaplotypeCaller …` (two args)
//! - `gatk --java-options=-Xmx8g HaplotypeCaller …` (long-flag-equals
//!   form, single arg)
//!
//! Shell quotes around `-Xmx8g` are commonly seen in GATK docs but
//! aren't load-bearing for argv splitting — `-Xmx8g` has no
//! whitespace, so the shell wouldn't split it regardless. The quotes
//! exist mostly to inoculate the value against any future glob /
//! history-expansion surprises.
//!
//! When we drive GATK through `Command::new` we don't go through a
//! shell at all, so argv splitting isn't a concern: passing
//! `--java-options` and `-Xmx<heap>` as two separate argv entries is
//! the cleanest form. We use that two-arg pattern (no `=` sign) for
//! symmetry with the rest of the adapter set, where `--flag value` is
//! the canonical shape.

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

use crate::case_input::GatkInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(GatkAdapter::new())
}

pub struct GatkAdapter;

impl GatkAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for GatkAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "gatk";
/// GATK's binary candidates. The Broad's release ships a Bash
/// launcher named `gatk` that wraps the JAR — Bioconda, the official
/// tarball, and source builds all install under that name.
const BINARIES: &[&str] = &["gatk"];

impl Adapter for GatkAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "GATK",
            // GATK 4.x is the current major line (the "GATK4 rewrite"
            // landed in 2018 and replaced the older 3.x Java codebase
            // entirely). We cap at 5.0 to reserve room for a future
            // major bump.
            version_range: VersionRange {
                min_inclusive: Version::new(4, 0, 0),
                max_exclusive: Version::new(5, 0, 0),
            },
            physics: &[Physics::Bio],
            // GATK4 is BSD-3-Clause licensed (the older 3.x line was
            // mixed). HaplotypeCaller specifically is in the
            // BSD-licensed core.
            license_mode: LicenseMode::Subprocess,
            tool_license: "BSD-3-Clause",
            docs_url: "https://gatk.broadinstitute.org/",
            homepage_url: "https://github.com/broadinstitute/gatk",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(BINARIES) {
            Some(binary_path) => {
                // `gatk --version` prints "The Genome Analysis Toolkit
                // (GATK) v4.4.0.0\n…" on stdout. The combined scanner
                // picks up the SemVer prefix.
                let found_version = detect_tool_version_semver(&binary_path, &["--version"]);
                let mut warnings = Vec::new();
                if found_version.is_none() {
                    // Binary on PATH but `--version` failed to produce
                    // a parseable string. The most common cause is a
                    // missing JDK — `gatk` is a thin wrapper script
                    // that delegates to `java -jar gatk-package*.jar`
                    // and bails noisily when java is absent. Surface
                    // the install hint instead of letting the user
                    // discover the Java requirement at run time.
                    warnings.push(
                        "GATK is on PATH but `gatk --version` failed — Java may be \
                         missing or misconfigured. Install JDK 17+ for HaplotypeCaller \
                         to run."
                            .into(),
                    );
                }
                Ok(ProbeReport {
                    ok: true,
                    found_version,
                    binary_path: Some(binary_path),
                    warnings,
                    required_env: Vec::new(),
                })
            }
            None => Err(AdapterError::ToolNotInstalled {
                name: INFO_ID,
                hint: "GATK 4.x required; install via \
                       `conda install -c bioconda gatk4`, the official \
                       release tarball from \
                       https://github.com/broadinstitute/gatk/releases, \
                       or the Docker image broadinstitute/gatk:latest"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = GatkInput::from_case_dir(&case.path)?;

        fs::create_dir_all(workdir)?;

        // Resolve every input path against the case directory via
        // `confined_join` — rejects absolute paths and `..` traversal
        // out of the case sandbox. A shared case bundle should not be
        // able to point `reference`, `input_bam`, or `intervals` at
        // `/etc/passwd` or similar.
        let source_reference = confined_join(&case.path, &input.reference)?;
        if !source_reference.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.gatk].reference `{}` not found (resolved {})",
                    input.reference.display(),
                    source_reference.display()
                ),
            });
        }

        let source_input_bam = confined_join(&case.path, &input.input_bam)?;
        if !source_input_bam.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.gatk].input_bam `{}` not found (resolved {})",
                    input.input_bam.display(),
                    source_input_bam.display()
                ),
            });
        }

        let source_intervals = match &input.intervals {
            Some(p) => {
                let resolved = confined_join(&case.path, p)?;
                if !resolved.is_file() {
                    return Err(AdapterError::InvalidCase {
                        case_path: case.path.join("case.toml"),
                        reason: format!(
                            "[bio.gatk].intervals `{}` not found (resolved {})",
                            p.display(),
                            resolved.display()
                        ),
                    });
                }
                Some(resolved)
            }
            None => None,
        };

        // Validate `output_vcf` is sandboxed too. We don't pre-create
        // it (HaplotypeCaller writes it), but the path needs to stay
        // inside the case dir so a hostile bundle can't write to
        // `/etc/cron.d/owned`. The string itself is what gets passed
        // to `gatk -O` below — keep the sandbox check on the resolved
        // form so we fail prepare rather than silently writing
        // outside the workdir.
        let _output_vcf_sandboxed = confined_join(&case.path, &input.output_vcf)?;

        let binary_path = find_on_path(BINARIES).ok_or_else(|| AdapterError::ToolNotInstalled {
            name: INFO_ID,
            hint: "GATK 4.x required; install via \
                       `conda install -c bioconda gatk4`, the official \
                       release tarball from \
                       https://github.com/broadinstitute/gatk/releases, \
                       or the Docker image broadinstitute/gatk:latest"
                .into(),
        })?;

        // Compose:
        //   gatk --java-options -Xmx<heap> HaplotypeCaller
        //     -R <reference> -I <input_bam> -O <output_vcf>
        //     [-L <intervals>] [extras...]
        //
        // `--java-options` and the `-Xmx<heap>` payload travel as two
        // separate argv entries — see the module-level note on why
        // that's the right shape when bypassing a shell.
        let mut native_command: Vec<OsString> = vec![
            binary_path.into_os_string(),
            OsString::from("--java-options"),
            OsString::from(format!("-Xmx{}", input.java_heap)),
            OsString::from("HaplotypeCaller"),
            OsString::from("-R"),
            source_reference.into_os_string(),
            OsString::from("-I"),
            source_input_bam.into_os_string(),
            OsString::from("-O"),
            OsString::from(&input.output_vcf),
        ];
        if let Some(intervals) = source_intervals {
            native_command.push(OsString::from("-L"));
            native_command.push(intervals.into_os_string());
        }
        for arg in &input.extra_args {
            native_command.push(OsString::from(arg));
        }

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            // Whole-genome HaplotypeCaller can run for many hours on a
            // single node; whole-exome runs an hour or so. 8 hours is
            // a generous default that covers the long tail.
            estimated_runtime: Some(Duration::from_secs(8 * 60 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting GATK HaplotypeCaller", |line| {
            let mut hint = subprocess::Hint::default();
            // GATK's stderr carries `INFO ProgressMeter` lines roughly
            // every 10 seconds — lift those to a 50% progress tick so
            // the UI shows forward motion. The `Done.` sentinel marks
            // the end-of-run summary; pin it at 95%.
            if line.contains("Done.") || line.contains("HaplotypeCaller done") {
                hint.progress = Some((95.0, line.to_string()));
            } else if line.contains("ProgressMeter") {
                hint.progress = Some((50.0, line.to_string()));
            } else if line.contains(" ERROR ") || line.contains("Exception") {
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
        // Recover the `-O` output path so we can surface it as an
        // artifact and use it for provenance hashing.
        let output_path = output_after_flag(job, "-O");

        let case_hash_input = output_path
            .clone()
            .filter(|p| p.is_file())
            .unwrap_or_else(|| job.workdir.join("case.toml"));
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "GATK",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);
        let mut artefacts: Vec<Artifact> = Vec::new();

        if let Some(out) = output_path {
            if out.is_file() {
                artefacts.push(Artifact {
                    path: out,
                    kind: ArtifactKind::Tabular,
                    checksum: None,
                    label: "GATK HaplotypeCaller VCF".to_string(),
                });
            }
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
            ribbon_contributions: vec!["bio.gatk.haplotype-caller"],
        }
    }
}

/// Walk the prepared command for the value following `flag`. Used
/// from `collect()` to recover the `-O` output path so we can surface
/// it as an artifact.
fn output_after_flag(job: &PreparedJob, flag: &str) -> Option<PathBuf> {
    let mut iter = job.native_command.iter().peekable();
    while let Some(arg) = iter.next() {
        if arg.to_str() == Some(flag) {
            if let Some(val) = iter.next() {
                return Some(PathBuf::from(val));
            }
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn info_is_bio_domain() {
        let info = GatkAdapter::new().info();
        assert_eq!(info.id, "gatk");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "BSD-3-Clause");
        assert_eq!(info.display_name, "GATK");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = GatkAdapter::new().info();
        // GATK4 is the current major line; 5.0 reserves room for an
        // eventual major bump.
        assert_eq!(info.version_range.min_inclusive, Version::new(4, 0, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(5, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = GatkAdapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.gatk.haplotype-caller"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = GatkAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }
}
