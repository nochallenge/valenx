//! # valenx-adapter-hmmer
//!
//! Adapter for [HMMER](http://hmmer.org/) — Sean Eddy's profile-HMM
//! search suite. HMMER builds and queries profile hidden Markov models
//! over protein and nucleotide sequences; the canonical entry points
//! are `hmmsearch` (search a profile against a sequence database) and
//! `hmmscan` (scan a query against a profile database). The pair
//! powers Pfam, SMART, and most of the structured protein-family
//! annotation pipelines in modern bioinformatics.
//!
//! **Phase 18 — subprocess wrapper around `hmmsearch` / `hmmscan`.**
//! The user picks the subcommand via `tool` in `[bio.hmmer]` and
//! supplies a profile + a sequences FASTA. `prepare()` resolves both
//! against the case directory and composes the invocation. Both tools
//! emit a tabular summary via `--tblout <file>` and a verbose
//! human-readable report via `-o <file>`; we pin both to fixed names
//! in the workdir so `collect()` can find them deterministically.
//!
//! Because both HMMER outputs go to files (not stdout), `run()` uses
//! the shared [`valenx_core::subprocess::run`] runner — same shape as
//! BWA / minimap2 / MUSCLE.

#![forbid(unsafe_code)]
#![allow(missing_docs)]

pub mod case_input;
pub mod native;

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

use crate::case_input::HmmerInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(HmmerAdapter::new())
}

pub struct HmmerAdapter;

impl HmmerAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for HmmerAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "hmmer";
/// HMMER's probe binary. We look for `hmmsearch` because it ships
/// with every HMMER install — if it's missing, the suite isn't
/// installed. `hmmscan` would do equally well; pick one as the
/// canonical signal.
const PROBE_BINARIES: &[&str] = &["hmmsearch"];

/// The tabular-output filename written by `--tblout`. Pinned so
/// `prepare()` (which records the command), and `collect()` (which
/// labels the artifact) all agree on what to look for.
const OUT_TBLOUT: &str = "tblout.txt";

/// The verbose human-readable report written by `-o`. Pinned for the
/// same reason as `OUT_TBLOUT`.
const OUT_REPORT: &str = "hmmer.out";

impl Adapter for HmmerAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "HMMER",
            // HMMER 3.3 is the long-running stable line every distro
            // ships (3.3 in 2019, 3.3.2 in 2020, 3.4 in 2023). Floor
            // at 3.3.0 covers every reasonably modern install; upper
            // bound 4.0 reserves room for the next major.
            version_range: VersionRange {
                min_inclusive: Version::new(3, 3, 0),
                max_exclusive: Version::new(4, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            tool_license: "BSD-3-Clause",
            docs_url: "http://eddylab.org/software/hmmer/Userguide.pdf",
            homepage_url: "http://hmmer.org/",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(PROBE_BINARIES) {
            Some(binary_path) => {
                // `hmmsearch -h` prints "# HMMER 3.3.2 (Nov 2020)" on
                // stdout; the combined stdout+stderr scanner picks
                // the version up either way.
                let found_version = detect_tool_version_semver(&binary_path, &["-h", "--version"]);
                Ok(ProbeReport {
                    ok: true,
                    found_version,
                    binary_path: Some(binary_path),
                    warnings: Vec::new(),
                    required_env: Vec::new(),
                })
            }
            // Native Rust fallback: profile-HMM via valenx-align.
            // Profile must be a FASTA alignment in native mode.
            None => Ok(ProbeReport {
                ok: true,
                found_version: None,
                binary_path: None,
                warnings: vec![
                    "hmmsearch binary not found; using native Rust profile-HMM Viterbi search \
                     (valenx-align). Install HMMER 3.3+ via apt/brew/conda for calibrated \
                     E-values and pre-built .hmm profile support. Native mode requires the \
                     `profile` field to point to a FASTA multiple alignment (.fa/.fasta/.aln), \
                     not a prebuilt .hmm database."
                        .to_string(),
                ],
                required_env: Vec::new(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = HmmerInput::from_case_dir(&case.path)?;

        fs::create_dir_all(workdir)?;

        // Resolve the profile path against the case directory if
        // relative. Same convention as every other Phase 17/18 bio
        // adapter — `profile = "Pfam-A.hmm"` next to `case.toml`.
        let source_profile = resolve_input(&case.path, &input.profile, "profile")?;
        let source_sequences = resolve_input(&case.path, &input.sequences, "sequences")?;

        // Write native_params.toml for the native path.
        let native_hmmer_params = native::NativeHmmerParams {
            profile_fasta_path: source_profile
                .to_str()
                .ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "profile path is not valid UTF-8: {}",
                        source_profile.display()
                    ))
                })?
                .to_string(),
            sequences_path: source_sequences
                .to_str()
                .ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "sequences path is not valid UTF-8: {}",
                        source_sequences.display()
                    ))
                })?
                .to_string(),
            min_score: -100.0,
            tblout_name: OUT_TBLOUT.to_string(),
            report_name: OUT_REPORT.to_string(),
        };
        native::write_params(workdir, &native_hmmer_params)?;

        // Build the command: real binary or native sentinel.
        let native_command: Vec<OsString> = match find_on_path(&[input.tool.as_str()]) {
            Some(binary_path) => {
                let mut cmd: Vec<OsString> = vec![
                    binary_path.into_os_string(),
                    OsString::from("--cpu"),
                    OsString::from(input.cpus.to_string()),
                    OsString::from("-E"),
                    OsString::from(format_evalue(input.evalue)),
                    OsString::from("--tblout"),
                    OsString::from(OUT_TBLOUT),
                    OsString::from("-o"),
                    OsString::from(OUT_REPORT),
                ];
                // Round-4 fix: extra_args after positionals.
                cmd.push(source_profile.into_os_string());
                cmd.push(source_sequences.into_os_string());
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
            // Small profile vs small DB runs in seconds; full
            // Pfam-vs-proteome runs in tens of minutes; running the
            // full UniRef as a sequence DB can be hours. 4 hours is
            // a generous default that covers the long tail without
            // being absurd.
            estimated_runtime: Some(Duration::from_secs(4 * 60 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        // Native Rust path: profile-HMM Viterbi search, no subprocess.
        if job.native_command.first().map(|s| s.as_os_str())
            == Some(native::NATIVE_SENTINEL.as_ref())
        {
            return native::run_native(&job.workdir, ctx);
        }

        let report = subprocess::run(job, ctx, "starting HMMER", |line| {
            let mut hint = subprocess::Hint::default();
            // HMMER's stdout chatter is sparse: a banner at startup,
            // optional `# ` progress markers when `--cpu > 1`, and an
            // "Internal pipeline statistics summary" near the end.
            // Lift those to coarse UI ticks.
            if line.contains("Internal pipeline statistics summary")
                || line.starts_with("# CPU time:")
            {
                hint.progress = Some((95.0, line.to_string()));
            } else if line.contains("Scoring profiles") || line.contains("Searching") {
                hint.progress = Some((50.0, line.to_string()));
            } else if line.contains("Error") || line.contains("FATAL") {
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
        // Provenance: hash the staged tblout if present (the canonical
        // run output). Falls back to case.toml when the search hasn't
        // produced a tblout yet — keeps the provenance block well-formed
        // for partial / failed runs.
        let case_hash_input = {
            let tblout = job.workdir.join(OUT_TBLOUT);
            if tblout.is_file() {
                tblout
            } else {
                job.workdir.join("case.toml")
            }
        };
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "HMMER",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);
        let mut artefacts: Vec<Artifact> = Vec::new();

        // Walk the workdir top-level. HMMER writes the canonical
        // `tblout.txt` (tabular) + `hmmer.out` (verbose report). We
        // also pick up any `.log` the user redirected stderr to in
        // case future cases configure that.
        let entries = match fs::read_dir(&job.workdir) {
            Ok(e) => e,
            Err(e) => {
                tracing::warn!(target: "valenx-hmmer", ?e, "workdir read failed");
                return Ok(results);
            }
        };
        for entry in entries.flatten() {
            let path = entry.path();
            if !path.is_file() {
                continue;
            }
            let file_name = path
                .file_name()
                .and_then(|s| s.to_str())
                .map(|s| s.to_string());
            let ext = path
                .extension()
                .and_then(|s| s.to_str())
                .map(|s| s.to_ascii_lowercase());

            // Pin the canonical outputs by exact filename so they get
            // the right artifact kind. Tabular for the tblout (parser-
            // friendly), Log for the verbose report (human-readable).
            let (kind, label) = match (file_name.as_deref(), ext.as_deref()) {
                (Some(name), _) if name == OUT_TBLOUT => {
                    (ArtifactKind::Tabular, "HMMER hits (tblout)".to_string())
                }
                (Some(name), _) if name == OUT_REPORT => {
                    (ArtifactKind::Log, "HMMER report".to_string())
                }
                (_, Some("log")) => (ArtifactKind::Log, "HMMER log".to_string()),
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
            ribbon_contributions: vec!["bio.hmmer.search"],
        }
    }
}

/// Resolve a case-relative path against the case directory and
/// confirm it points at a regular file. Adapter-internal — the same
/// shape every other Phase 18 adapter inlines, factored out here so
/// the two HMMER inputs (profile + sequences) share one validator.
///
/// Round-11 fix (R11-5): pre-fix the relative-path branch did a bare
/// `case_dir.join(raw)`, which silently accepted `..` traversal out
/// of the case sandbox (e.g. `profile = "../../etc/passwd"`). Route
/// through `confined_join` so the same threat model every other
/// Phase 18 adapter (BWA / minimap2 / MAFFT / muscle / samtools)
/// applies to the HMMER inputs too.
fn resolve_input(case_dir: &Path, raw: &Path, label: &str) -> Result<PathBuf, AdapterError> {
    let resolved = if raw.is_absolute() {
        raw.to_path_buf()
    } else {
        confined_join(case_dir, raw).map_err(|e| match e {
            // Re-wrap so the adapter's error carries the same
            // `[bio.hmmer].{label}` framing every other diagnostic
            // here uses; the inner message still names the
            // offending path.
            AdapterError::InvalidCase { case_path, reason } => AdapterError::InvalidCase {
                case_path,
                reason: format!("[bio.hmmer].{label}: {reason}"),
            },
            other => other,
        })?
    };
    if !resolved.is_file() {
        return Err(AdapterError::InvalidCase {
            case_path: case_dir.join("case.toml"),
            reason: format!(
                "[bio.hmmer].{label} `{}` not found (resolved {})",
                raw.display(),
                resolved.display()
            ),
        });
    }
    Ok(resolved)
}

/// Format an E-value for the HMMER CLI. HMMER's `-E` flag accepts
/// both decimal (`0.01`) and scientific (`1e-10`) notation; Rust's
/// default float formatting renders very small or very large numbers
/// in scientific automatically and very moderate ones in decimal,
/// which is exactly the right shape.
fn format_evalue(e: f64) -> String {
    // Use `{}` (Display) — drops trailing zeros and switches to
    // scientific for extreme values automatically.
    format!("{e}")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn info_is_bio_domain() {
        let info = HmmerAdapter::new().info();
        assert_eq!(info.id, "hmmer");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "BSD-3-Clause");
        assert_eq!(info.display_name, "HMMER");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = HmmerAdapter::new().info();
        // HMMER 3.3 is the floor we test against; 4.0 reserves room
        // for an eventual major bump.
        assert_eq!(info.version_range.min_inclusive, Version::new(3, 3, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(4, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = HmmerAdapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.hmmer.search"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = HmmerAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }

    #[test]
    fn probe_always_succeeds() {
        let report = HmmerAdapter::new().probe().unwrap();
        assert!(report.ok, "probe() must return ok=true in all modes");
    }

    #[test]
    fn native_sentinel_is_stable() {
        assert_eq!(native::NATIVE_SENTINEL, "valenx:native:hmmer");
    }

    /// Round-11 RED→GREEN — R11-5. Pre-fix `resolve_input` did a
    /// bare `case_dir.join(raw)` for relative paths, which silently
    /// accepted `..` traversal (e.g. `profile = "../../etc/passwd"`).
    /// Every other Phase 18 adapter (BWA / minimap2 / MAFFT / muscle
    /// / samtools) routes its case-relative inputs through
    /// `confined_join`; HMMER was the missing one. The fix wires it
    /// through too so a shared-case bundle cannot point HMMER at
    /// arbitrary host files.
    #[test]
    fn resolve_input_rejects_parent_dir_traversal() {
        // Build a real case directory so the `is_file()` check at
        // the end of resolve_input never short-circuits the
        // confined_join error path.
        let mut case_dir = std::env::temp_dir();
        let nanos = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map(|d| d.as_nanos())
            .unwrap_or(0);
        case_dir.push(format!("valenx-r11-hmmer-{nanos}"));
        std::fs::create_dir_all(&case_dir).expect("create case dir");

        // A relative path that escapes the case sandbox via `..`.
        let evil = PathBuf::from("../../etc/passwd");
        let err = resolve_input(&case_dir, &evil, "profile")
            .expect_err("`..` traversal must be rejected by confined_join");
        match err {
            AdapterError::InvalidCase { reason, .. } => {
                assert!(
                    reason.contains("[bio.hmmer].profile"),
                    "reason must carry adapter framing; got: {reason}"
                );
            }
            other => panic!("expected AdapterError::InvalidCase, got {other:?}"),
        }

        let _ = std::fs::remove_dir_all(&case_dir);
    }
}
