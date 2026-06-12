//! # valenx-adapter-deepvariant
//!
//! Adapter for [Google DeepVariant](https://github.com/google/deepvariant)
//! — the deep-learning germline variant caller from Google Health that
//! reframes variant calling as an image-classification problem. A
//! pileup of reads is rendered into an RGB-encoded "tensor" image and
//! a CNN trained on Genome in a Bottle truth sets calls the genotype
//! at each candidate site. DeepVariant typically posts the highest F1
//! scores in the precisionFDA Truth Challenge benchmarks, especially
//! on long-read data (PacBio HiFi and ONT R10.4) where it
//! significantly out-performs allele-frequency-based callers.
//!
//! **Phase 19 — subprocess wrapper around `run_deepvariant`.** The
//! `run_deepvariant` Python entry point is the official one-shot
//! driver that chains DeepVariant's three internal stages
//! (`make_examples` → `call_variants` → `postprocess_variants`) into
//! a single command. The user supplies a reference FASTA, a sorted
//! BAM, an output VCF path, the model type, and a shard count via
//! `[bio.deepvariant]` in `case.toml`. `prepare()` composes the
//! invocation; `run()` streams via the shared subprocess runner.
//!
//! ## Probe / install hints
//!
//! Most production users run DeepVariant via Docker
//! (`docker run google/deepvariant:latest`) or Singularity, since the
//! native Python install pulls in a full TensorFlow / numpy /
//! pysam stack. The adapter probes for the bare `run_deepvariant`
//! binary on PATH; if missing, the install hint surfaces both the
//! direct binary and the container-wrapper paths so users on
//! Docker-only setups know what to set up.

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

use crate::case_input::DeepVariantInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(DeepVariantAdapter::new())
}

pub struct DeepVariantAdapter;

impl DeepVariantAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for DeepVariantAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "deepvariant";
/// DeepVariant's binary candidates. The Bash launcher
/// `run_deepvariant` is the canonical name in both the native install
/// and the Docker / Singularity images that wrap it.
const BINARIES: &[&str] = &["run_deepvariant"];

/// Install hint surfaced when the binary isn't on PATH. Mentions both
/// the direct binary install (rare in practice) and the Docker /
/// Singularity wrapper paths most users actually run.
const INSTALL_HINT: &str = "DeepVariant 1.6+ required; either install the `run_deepvariant` \
     binary on PATH, or use the official container — \
     `docker run google/deepvariant:latest /opt/deepvariant/bin/run_deepvariant …` \
     (or the matching Singularity image). \
     See https://github.com/google/deepvariant for the up-to-date image tags.";

impl Adapter for DeepVariantAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "DeepVariant",
            // DeepVariant 1.6 (2023) is the floor we test against —
            // the model bundles for ONT_R104 landed in 1.6 and the
            // CLI shape we drive is stable from that release on. The
            // upper bound 2.0 reserves room for an eventual major
            // bump.
            version_range: VersionRange {
                min_inclusive: Version::new(1, 6, 0),
                max_exclusive: Version::new(2, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            tool_license: "BSD-3-Clause",
            docs_url: "https://github.com/google/deepvariant",
            homepage_url: "https://github.com/google/deepvariant",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(BINARIES) {
            Some(binary_path) => {
                // `run_deepvariant --version` prints the version on
                // stdout. The combined scanner picks up the SemVer
                // prefix.
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
                hint: INSTALL_HINT.into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = DeepVariantInput::from_case_dir(&case.path)?;

        fs::create_dir_all(workdir)?;

        // Resolve every input path against the case directory if
        // relative. Same convention as every other Phase 17/18 bio
        // adapter — `reference = "ref.fa"` next to `case.toml`.
        let source_reference = if input.reference.is_absolute() {
            input.reference.clone()
        } else {
            valenx_core::adapter_helpers::confined_join(&case.path, &input.reference)?
        };
        if !source_reference.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.deepvariant].reference `{}` not found (resolved {})",
                    input.reference.display(),
                    source_reference.display()
                ),
            });
        }

        let source_input_bam = if input.input_bam.is_absolute() {
            input.input_bam.clone()
        } else {
            valenx_core::adapter_helpers::confined_join(&case.path, &input.input_bam)?
        };
        if !source_input_bam.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.deepvariant].input_bam `{}` not found (resolved {})",
                    input.input_bam.display(),
                    source_input_bam.display()
                ),
            });
        }

        let binary_path = find_on_path(BINARIES).ok_or_else(|| AdapterError::ToolNotInstalled {
            name: INFO_ID,
            hint: INSTALL_HINT.into(),
        })?;

        // Compose:
        //   run_deepvariant --model_type=<type> --ref=<reference>
        //     --reads=<input_bam> --output_vcf=<output_vcf>
        //     --num_shards=<N> [extras...]
        //
        // DeepVariant uses the `--flag=value` shape for every
        // argument; we follow it verbatim so users who paste a stock
        // DeepVariant invocation into `extra_args` end up with a
        // consistent command line.
        let mut native_command: Vec<OsString> = vec![
            binary_path.into_os_string(),
            OsString::from(format!("--model_type={}", input.model_type)),
            OsString::from(format!("--ref={}", source_reference.display())),
            OsString::from(format!("--reads={}", source_input_bam.display())),
            OsString::from(format!("--output_vcf={}", input.output_vcf.display())),
            OsString::from(format!("--num_shards={}", input.num_shards)),
        ];
        for arg in &input.extra_args {
            native_command.push(OsString::from(arg));
        }

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            // DeepVariant on a 30x WGS sample runs many hours on CPU,
            // an hour or so on GPU. 12 hours covers both the long tail
            // of CPU-only runs and a typical multi-sample batch.
            estimated_runtime: Some(Duration::from_secs(12 * 60 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting DeepVariant", |line| {
            let mut hint = subprocess::Hint::default();
            // DeepVariant's three stages each print sentinel lines we
            // can lift to coarse progress ticks. `make_examples` is
            // the longest stage (typically 60–80% of wall time);
            // `call_variants` is the inference pass; `postprocess`
            // writes the final VCF.
            if line.contains("Done.") || line.contains("postprocess_variants") {
                hint.progress = Some((95.0, line.to_string()));
            } else if line.contains("call_variants") {
                hint.progress = Some((80.0, line.to_string()));
            } else if line.contains("make_examples") {
                hint.progress = Some((40.0, line.to_string()));
            } else if line.contains(" ERROR ") || line.contains("Traceback") {
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
        // Recover the `--output_vcf=<path>` value from the prepared
        // command. Unlike GATK / bcftools (which use `-O <file>`),
        // DeepVariant uses the `--flag=value` shape, so we scan for
        // the prefix instead of the standalone flag.
        let output_path = output_after_eq_flag(job, "--output_vcf=");

        let case_hash_input = output_path
            .clone()
            .filter(|p| p.is_file())
            .unwrap_or_else(|| job.workdir.join("case.toml"));
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "DeepVariant",
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
                    label: "DeepVariant VCF".to_string(),
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
            ribbon_contributions: vec!["bio.deepvariant.call"],
        }
    }
}

/// Walk the prepared command for the value following `--flag=` (the
/// shape DeepVariant uses for every argument). Used from `collect()`
/// to recover the `--output_vcf=` path so we can surface it as an
/// artifact.
fn output_after_eq_flag(job: &PreparedJob, prefix: &str) -> Option<PathBuf> {
    for arg in &job.native_command {
        if let Some(s) = arg.to_str() {
            if let Some(value) = s.strip_prefix(prefix) {
                return Some(PathBuf::from(value));
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
        let info = DeepVariantAdapter::new().info();
        assert_eq!(info.id, "deepvariant");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "BSD-3-Clause");
        assert_eq!(info.display_name, "DeepVariant");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = DeepVariantAdapter::new().info();
        // 1.6 is the floor we test against; 2.0 reserves room for an
        // eventual major bump.
        assert_eq!(info.version_range.min_inclusive, Version::new(1, 6, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(2, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = DeepVariantAdapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.deepvariant.call"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = DeepVariantAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }
}
