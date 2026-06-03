//! # valenx-adapter-mmseqs2
//!
//! Adapter for [MMseqs2](https://mmseqs.com/) — Söding lab's "many vs.
//! many sequence searching" toolkit. MMseqs2 is the modern fast
//! alternative to BLAST for protein-vs-protein search and large-scale
//! clustering. Its prefilter / vectorised k-mer step is what makes
//! ColabFold's MSA generation tractable, so this adapter is also a
//! prerequisite for any downstream protein-prediction workflow.
//!
//! **Phase 18.5 — subprocess wrapper around `mmseqs <action>`.** The
//! user picks one of three high-level "easy-" workflows via `action` in
//! `[bio.mmseqs2]`:
//!
//! - `easy-search`    — exhaustive iterative search of `query` against
//!   `target`. Output is a BLAST-format-8 hit table.
//! - `easy-linsearch` — linear-time prefilter variant of `easy-search`.
//! - `easy-cluster`   — clustering of `query` alone; no `target`.
//!
//! Per-action command construction is factored into [`build_command`]
//! so the dispatch is one self-contained match without polluting the
//! Adapter::prepare body. Unknown actions return `InvalidCase` so a
//! schema drift never panics.

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

use crate::case_input::Mmseqs2Input;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(Mmseqs2Adapter::new())
}

pub struct Mmseqs2Adapter;

impl Mmseqs2Adapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for Mmseqs2Adapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "mmseqs2";
/// MMseqs2's binary candidates. The on-disk binary is just `mmseqs`
/// (no version suffix); Bioconda, Homebrew, and source builds all
/// install under that name.
const BINARIES: &[&str] = &["mmseqs"];

/// Scratch / temp directory MMseqs2 needs as its last positional
/// argument for every easy-* workflow. Pinned per-run inside the
/// workdir so the runner is fully self-contained.
const TMP_DIR: &str = "tmp";

impl Adapter for Mmseqs2Adapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "MMseqs2",
            // MMseqs2 ships with sequential integer "release" version
            // numbers (e.g. 14-7e284, 15-6f452); semver-wise we treat
            // the first integer as the major and gate at >= 14.
            version_range: VersionRange {
                min_inclusive: Version::new(14, 0, 0),
                max_exclusive: Version::new(17, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            tool_license: "MIT",
            docs_url: "https://github.com/soedinglab/mmseqs2/wiki",
            homepage_url: "https://mmseqs.com/",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(BINARIES) {
            Some(binary_path) => {
                // `mmseqs version` (positional, not `--version`) prints
                // the version banner on stdout; the combined scanner
                // picks it up.
                let found_version =
                    detect_tool_version_semver(&binary_path, &["version", "--version"]);
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
                hint: "MMseqs2 14+ required; install via `brew install mmseqs2`, \
                       or `conda install -c bioconda mmseqs2`"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = Mmseqs2Input::from_case_dir(&case.path)?;

        fs::create_dir_all(workdir)?;

        // Resolve every input path against the case directory if
        // relative, but build per-action so the validation messages
        // can name the right field.
        let resolved_query = if input.query.is_absolute() {
            input.query.clone()
        } else {
            valenx_core::adapter_helpers::confined_join(
            &case.path,
            &input.query,
        )?
        };
        if !resolved_query.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.mmseqs2].query `{}` not found (resolved {})",
                    input.query.display(),
                    resolved_query.display()
                ),
            });
        }

        let resolved_target: Option<PathBuf> = match &input.target {
            Some(p) => {
                // Round-9 hardening: relative `target` flows into the
                // mmseqs2 command line; wrap with `confined_join`.
                let resolved = if p.is_absolute() {
                    p.clone()
                } else {
                    valenx_core::adapter_helpers::confined_join(&case.path, p)?
                };
                if !resolved.is_file() {
                    return Err(AdapterError::InvalidCase {
                        case_path: case.path.join("case.toml"),
                        reason: format!(
                            "[bio.mmseqs2].target `{}` not found (resolved {})",
                            p.display(),
                            resolved.display()
                        ),
                    });
                }
                Some(resolved)
            }
            None => None,
        };

        let binary_path = find_on_path(BINARIES).ok_or_else(|| AdapterError::ToolNotInstalled {
            name: INFO_ID,
            hint: "MMseqs2 14+ required; install via `brew install mmseqs2`, \
                       or `conda install -c bioconda mmseqs2`"
                .into(),
        })?;

        let native_command = build_command(
            &binary_path,
            &resolved_query,
            resolved_target.as_deref(),
            &input,
            &case.path,
        )?;

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            // Searches against UniRef90 take hours on a single node;
            // 4 hours mirrors the rest of the bio adapters' generous
            // default for long-tail workloads.
            estimated_runtime: Some(Duration::from_secs(4 * 60 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting MMseqs2", |line| {
            let mut hint = subprocess::Hint::default();
            // MMseqs2 prints "Step N of M" markers and per-stage
            // timing summaries; the "Time for processing" line marks
            // the end of the run.
            if line.contains("Time for processing") {
                hint.progress = Some((95.0, line.to_string()));
            } else if line.contains("Step ") && line.contains(" of ") {
                hint.progress = Some((50.0, line.to_string()));
            } else if line.contains("Error") || line.contains("ERROR") {
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
        // Re-derive the action and output path from the prepared
        // command. The first arg after the binary is always the
        // subcommand; the output path slot is action-dependent.
        let action = job
            .native_command
            .get(1)
            .and_then(|s| s.to_str())
            .unwrap_or("");
        let output_path = mmseqs_output_path(job, action);

        let case_hash_input = output_path
            .clone()
            .filter(|p| p.is_file())
            .unwrap_or_else(|| job.workdir.join("case.toml"));
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "MMseqs2",
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
                let label = match action {
                    "easy-search" => "MMseqs2 easy-search hits".to_string(),
                    "easy-linsearch" => "MMseqs2 easy-linsearch hits".to_string(),
                    "easy-cluster" => "MMseqs2 easy-cluster output".to_string(),
                    _ => "MMseqs2 output".to_string(),
                };
                artefacts.push(Artifact {
                    path: out.clone(),
                    kind: ArtifactKind::Tabular,
                    checksum: None,
                    label,
                });
            }
        }

        artefacts.sort_by(|a, b| a.path.cmp(&b.path));
        results.artifacts = artefacts;
        Ok(results)
    }

    fn capabilities(&self) -> Capabilities {
        Capabilities {
            capabilities: Vec::new(),
            ribbon_contributions: vec!["bio.mmseqs2.search"],
        }
    }
}

/// Compose the `mmseqs` invocation for the given action.
///
/// Each easy-* workflow takes a fixed positional shape ending in a
/// scratch tmp dir. Sensitivity is meaningless for `easy-linsearch`
/// (the linear-time prefilter has no `-s` knob) so we pass it only
/// where MMseqs2 expects it.
///
/// Shapes:
///
/// - `easy-search`:    `mmseqs easy-search <query> <target> <output> tmp -s <s> --threads N [extras...]`
/// - `easy-linsearch`: `mmseqs easy-linsearch <query> <target> <output> tmp --threads N [extras...]`
/// - `easy-cluster`:   `mmseqs easy-cluster <query> <output> tmp -s <s> --threads N [extras...]`
pub fn build_command(
    binary_path: &Path,
    resolved_query: &Path,
    resolved_target: Option<&Path>,
    case: &Mmseqs2Input,
    case_path: &Path,
) -> Result<Vec<OsString>, AdapterError> {
    let mut cmd: Vec<OsString> = vec![
        binary_path.as_os_str().to_owned(),
        OsString::from(&case.action),
    ];

    match case.action.as_str() {
        "easy-search" => {
            cmd.push(resolved_query.as_os_str().to_owned());
            cmd.push(
                resolved_target
                    .expect("case_input enforces target for easy-search")
                    .as_os_str()
                    .to_owned(),
            );
            cmd.push(case.output.as_os_str().to_owned());
            cmd.push(OsString::from(TMP_DIR));
            cmd.push(OsString::from("-s"));
            cmd.push(OsString::from(format_sensitivity(case.sensitivity)));
            cmd.push(OsString::from("--threads"));
            cmd.push(OsString::from(case.threads.to_string()));
            for arg in &case.extra_args {
                cmd.push(OsString::from(arg));
            }
        }
        "easy-linsearch" => {
            cmd.push(resolved_query.as_os_str().to_owned());
            cmd.push(
                resolved_target
                    .expect("case_input enforces target for easy-linsearch")
                    .as_os_str()
                    .to_owned(),
            );
            cmd.push(case.output.as_os_str().to_owned());
            cmd.push(OsString::from(TMP_DIR));
            cmd.push(OsString::from("--threads"));
            cmd.push(OsString::from(case.threads.to_string()));
            for arg in &case.extra_args {
                cmd.push(OsString::from(arg));
            }
        }
        "easy-cluster" => {
            cmd.push(resolved_query.as_os_str().to_owned());
            // The output slot for easy-cluster is the *basename* —
            // MMseqs2 writes <basename>_cluster.tsv, _rep_seq.fasta,
            // and _all_seqs.fasta beside it.
            cmd.push(case.output.as_os_str().to_owned());
            cmd.push(OsString::from(TMP_DIR));
            cmd.push(OsString::from("-s"));
            cmd.push(OsString::from(format_sensitivity(case.sensitivity)));
            cmd.push(OsString::from("--threads"));
            cmd.push(OsString::from(case.threads.to_string()));
            for arg in &case.extra_args {
                cmd.push(OsString::from(arg));
            }
        }
        // Defensive — case_input rejects unknown actions, but a
        // future schema-only change shouldn't reach a panic. Surface
        // a soft InvalidCase so the user gets a helpful error if the
        // schema ever drifts.
        other => {
            return Err(AdapterError::InvalidCase {
                case_path: case_path.join("case.toml"),
                reason: format!(
                    "internal: unknown mmseqs2 action `{other}` slipped past schema validation"
                ),
            });
        }
    }
    Ok(cmd)
}

/// Format the sensitivity float for the CLI, dropping a trailing
/// `.0` on whole-number values so the command line stays compact and
/// matches MMseqs2's own examples (`-s 7.5`, `-s 4`).
fn format_sensitivity(s: f64) -> String {
    if s.fract() == 0.0 {
        format!("{}", s as i64)
    } else {
        format!("{s}")
    }
}

/// Recover the output path from the prepared command. For both
/// search variants, the output is the third positional after the
/// subcommand (`mmseqs <action> <query> <target> <output> ...`); for
/// `easy-cluster` it's the second (`mmseqs easy-cluster <query>
/// <output> ...`).
fn mmseqs_output_path(job: &PreparedJob, action: &str) -> Option<PathBuf> {
    let position = match action {
        "easy-search" | "easy-linsearch" => 4, // [bin, action, query, target, output]
        "easy-cluster" => 3,                   // [bin, action, query, output]
        _ => return None,
    };
    job.native_command.get(position).map(PathBuf::from)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn info_is_bio_domain() {
        let info = Mmseqs2Adapter::new().info();
        assert_eq!(info.id, "mmseqs2");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "MIT");
        assert_eq!(info.display_name, "MMseqs2");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = Mmseqs2Adapter::new().info();
        assert_eq!(info.version_range.min_inclusive, Version::new(14, 0, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(17, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = Mmseqs2Adapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.mmseqs2.search"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = Mmseqs2Adapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }

    /// Round-9 RED→GREEN: `[bio.mmseqs2].target` (optional second
    /// FASTA for the `easy-search`-style actions) used to be joined
    /// with bare `case.path.join`. Wrap with `confined_join`.
    #[test]
    fn prepare_rejects_target_traversing_outside_case_dir() {
        use valenx_test_utils::tempdir;
        let d = tempdir("mmseqs2-target-trav");
        std::fs::write(d.join("q.fa"), b">x\nACGT\n").unwrap();
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "mmseqs2.search"

[bio.mmseqs2]
action = "easy-search"
query  = "q.fa"
target = "../../etc/passwd"
output = "results.m8"
"#,
        )
        .unwrap();
        let case = Case {
            id: "mmseqs2-target-trav".into(),
            path: d.clone(),
        };
        let workdir = d.join("workdir");
        let err = Mmseqs2Adapter::new().prepare(&case, &workdir).unwrap_err();
        let msg = format!("{err}");
        assert!(
            msg.contains("..") || msg.contains("stay within") || msg.contains("escape"),
            "expected confined_join rejection, got: {msg}"
        );
        let _ = std::fs::remove_dir_all(&d);
    }
}
