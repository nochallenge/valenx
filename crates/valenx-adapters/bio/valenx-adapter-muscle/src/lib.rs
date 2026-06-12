//! # valenx-adapter-muscle
//!
//! Adapter for [MUSCLE](https://drive5.com/muscle/) — Robert Edgar's
//! multiple-sequence-alignment package. MUSCLE 5 is a complete
//! rewrite around a probabilistic-consistency aligner with a
//! divide-and-conquer "super5" mode for very large inputs (10k+
//! sequences). Faster and more accurate than the long-running
//! MUSCLE 3.x line; the canonical MSA tool alongside MAFFT.
//!
//! **Phase 18 — subprocess wrapper around `muscle`.** The user
//! supplies a multi-FASTA via `[bio.muscle]` in `case.toml`.
//! `prepare()` resolves the input against the case directory, picks
//! the mode (`align` by default), and composes the `muscle`
//! invocation. MUSCLE 5 writes its output to the path given via
//! `-output <file>`, so the standard subprocess runner shape works
//! cleanly — `run()` streams stderr through the line handler for
//! progress, and `collect()` surfaces the canonical `aligned.fa`
//! plus any `.log` files.

#![forbid(unsafe_code)]
#![allow(missing_docs)]

pub mod case_input;
pub mod native;

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

use crate::case_input::MuscleInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(MuscleAdapter::new())
}

pub struct MuscleAdapter;

impl MuscleAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for MuscleAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "muscle";
/// MUSCLE's binary candidates. The MUSCLE 5 packaged builds install
/// under the canonical `muscle` name on every platform.
const BINARIES: &[&str] = &["muscle"];

/// The aligned-FASTA filename written by `muscle -output`. Pinned so
/// `prepare()`, `collect()`, and the artifact label all agree on
/// what to look for.
const OUT_FA: &str = "aligned.fa";

impl Adapter for MuscleAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "MUSCLE",
            // MUSCLE 5.x is a complete rewrite from the long-running
            // 3.x line; floor at 5.1 (the first widely-distributed
            // 5.x release on Bioconda) and reserve room for an
            // eventual MUSCLE 6.
            version_range: VersionRange {
                min_inclusive: Version::new(5, 1, 0),
                max_exclusive: Version::new(6, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            tool_license: "GPL-3.0",
            docs_url: "https://drive5.com/muscle5/manual/",
            homepage_url: "https://drive5.com/muscle/",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(BINARIES) {
            Some(binary_path) => {
                let found_version =
                    detect_tool_version_semver(&binary_path, &["-version", "--version", ""]);
                Ok(ProbeReport {
                    ok: true,
                    found_version,
                    binary_path: Some(binary_path),
                    warnings: Vec::new(),
                    required_env: Vec::new(),
                })
            }
            // Native Rust fallback: progressive + iterative MSA via valenx-align.
            None => Ok(ProbeReport {
                ok: true,
                found_version: None,
                binary_path: None,
                warnings: vec![
                    "muscle binary not found; using native Rust progressive+iterative MSA \
                     (valenx-align). Install MUSCLE 5 via conda/bioconda for the full \
                     MUSCLE 5 probabilistic-consistency aligner."
                        .to_string(),
                ],
                required_env: Vec::new(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = MuscleInput::from_case_dir(&case.path)?;

        fs::create_dir_all(workdir)?;

        // Resolve the input FASTA against the case directory if
        // relative. Same convention as every other Phase 17/18 bio
        // adapter — `input = "seqs.fa"` next to `case.toml`.
        let source_input = if input.input.is_absolute() {
            input.input.clone()
        } else {
            valenx_core::adapter_helpers::confined_join(&case.path, &input.input)?
        };
        if !source_input.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.muscle].input `{}` not found (resolved {})",
                    input.input.display(),
                    source_input.display()
                ),
            });
        }

        // Write native_params.toml — run() reads this for the native path.
        let native_params = native::NativeMsaParams {
            input_path: source_input
                .to_str()
                .ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "input path is not valid UTF-8: {}",
                        source_input.display()
                    ))
                })?
                .to_string(),
            output_name: OUT_FA.to_string(),
            // In "align" mode, use iterative refinement; "super5" is
            // the large-dataset mode — same progressive alignment, no
            // extra refinement passes in native mode.
            refine: input.mode == "align",
            max_iterations: 8,
        };
        native::write_params(workdir, &native_params)?;

        let native_command: Vec<OsString> = match find_on_path(BINARIES) {
            Some(binary_path) => {
                let mut cmd: Vec<OsString> = vec![
                    binary_path.into_os_string(),
                    OsString::from(format!("-{}", input.mode)),
                    source_input.into_os_string(),
                    OsString::from("-output"),
                    OsString::from(OUT_FA),
                ];
                if let Some(t) = input.threads {
                    cmd.push(OsString::from("-threads"));
                    cmd.push(OsString::from(t.to_string()));
                }
                for arg in &input.extra_args {
                    cmd.push(OsString::from(arg));
                }
                cmd
            }
            None => vec![OsString::from(native::NATIVE_SENTINEL)],
        };

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            estimated_runtime: Some(Duration::from_secs(4 * 60 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        // Native Rust path: progressive + iterative MSA, no subprocess.
        if job.native_command.first().map(|s| s.as_os_str())
            == Some(native::NATIVE_SENTINEL.as_ref())
        {
            return native::run_native(&job.workdir, ctx);
        }

        let report = subprocess::run(job, ctx, "starting MUSCLE", |line| {
            let mut hint = subprocess::Hint::default();
            // MUSCLE 5 logs progress markers like "00:00:14   123Mb
            // CPU 0% Refining" and "00:00:25   456Mb CPU 0% Done."
            // on stdout. Lift them to coarse UI ticks; pin "Done."
            // at 95% so the progress bar visibly approaches
            // completion before the run wraps.
            if line.contains("Done.") {
                hint.progress = Some((95.0, line.to_string()));
            } else if line.contains("Refining") || line.contains("Aligning") {
                hint.progress = Some((70.0, line.to_string()));
            } else if line.contains("Iter") || line.contains("guide tree") {
                hint.progress = Some((40.0, line.to_string()));
            } else if line.contains("ERROR") || line.contains("FATAL") {
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
        // Provenance: hash the staged aligned.fa if present, falling
        // back to case.toml when the run hasn't produced one yet.
        let case_hash_input = {
            let fa = job.workdir.join(OUT_FA);
            if fa.is_file() {
                fa
            } else {
                job.workdir.join("case.toml")
            }
        };
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "MUSCLE",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);
        let mut artefacts: Vec<Artifact> = Vec::new();

        // Walk the workdir top-level. MUSCLE writes `aligned.fa`
        // (our `-output` target); a `.log` file may appear if future
        // cases configure stderr redirection.
        let entries = match fs::read_dir(&job.workdir) {
            Ok(e) => e,
            Err(e) => {
                tracing::warn!(target: "valenx-muscle", ?e, "workdir read failed");
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
                // The aligned-FASTA we requested via `-output`. `fa`
                // is the canonical extension; some pipelines use
                // `fasta` / `aln` instead.
                Some("fa") | Some("fasta") | Some("aln") => {
                    (ArtifactKind::Native, "MUSCLE alignment (FASTA)".to_string())
                }
                // MUSCLE 5 can emit `.efa` (ensemble FASTA) when
                // running with `-perm` to produce alignment
                // ensembles.
                Some("efa") => (
                    ArtifactKind::Native,
                    "MUSCLE ensemble alignment".to_string(),
                ),
                Some("log") => (ArtifactKind::Log, "MUSCLE log".to_string()),
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
            ribbon_contributions: vec!["bio.muscle.msa"],
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn info_is_bio_domain() {
        let info = MuscleAdapter::new().info();
        assert_eq!(info.id, "muscle");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "GPL-3.0");
        assert_eq!(info.display_name, "MUSCLE");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = MuscleAdapter::new().info();
        // MUSCLE 5.x is the modern rewrite; 6.0 reserves room for an
        // eventual major bump.
        assert_eq!(info.version_range.min_inclusive, Version::new(5, 1, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(6, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = MuscleAdapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.muscle.msa"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = MuscleAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }

    #[test]
    fn probe_always_succeeds() {
        let report = MuscleAdapter::new().probe().unwrap();
        assert!(report.ok, "probe() must return ok=true in all modes");
    }

    #[test]
    fn native_sentinel_is_stable() {
        assert_eq!(native::NATIVE_SENTINEL, "valenx:native:msa");
    }
}
