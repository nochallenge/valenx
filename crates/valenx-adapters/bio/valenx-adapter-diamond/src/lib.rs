//! # valenx-adapter-diamond
//!
//! Adapter for [DIAMOND](https://github.com/bbuchfink/diamond) —
//! Buchfink, Reuter & Drost's ultra-fast protein aligner. DIAMOND
//! covers the same ground as BLASTP / BLASTX but is two to three
//! orders of magnitude faster, which is what makes it the go-to tool
//! for whole-metagenome searches and for indexing UniRef-scale
//! databases up front.
//!
//! **Phase 18.5 — subprocess wrapper around `diamond <action>`.** The
//! user picks the mode via `action` in `[bio.diamond]`:
//!
//! - `blastp`  — protein-vs-protein search.
//! - `blastx`  — translated nucleotide-vs-protein search.
//! - `makedb`  — build a `.dmnd` database from a FASTA.
//!
//! Per-action command construction is factored into [`build_command`].
//! Two DIAMOND quirks the adapter handles:
//!
//! 1. The `--default` sensitivity flag does not exist — DIAMOND's
//!    out-of-the-box default has no flag. So when the user picks
//!    `sensitivity = "default"` the adapter omits the flag entirely.
//! 2. In `makedb` mode the schema field roles flip: `query` is the
//!    *input* FASTA and `database` is the *output* DB basename
//!    (DIAMOND appends `.dmnd`). That's how the upstream CLI is
//!    actually shaped, and the adapter mirrors it directly so the
//!    schema names stay stable across actions.

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

use crate::case_input::DiamondInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(DiamondAdapter::new())
}

pub struct DiamondAdapter;

impl DiamondAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for DiamondAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "diamond";
/// DIAMOND's binary candidates. `diamond` is the canonical Linux /
/// macOS install name from Bioconda and source builds.
const BINARIES: &[&str] = &["diamond"];

impl Adapter for DiamondAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "DIAMOND",
            // DIAMOND 2.1.x is the long-running stable line; 3.0
            // reserves room for an eventual major bump.
            version_range: VersionRange {
                min_inclusive: Version::new(2, 1, 0),
                max_exclusive: Version::new(3, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            tool_license: "GPL-3.0",
            docs_url: "https://github.com/bbuchfink/diamond/wiki",
            homepage_url: "https://github.com/bbuchfink/diamond",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(BINARIES) {
            Some(binary_path) => {
                // `diamond version` (positional) prints the version
                // banner on stdout. `--version` also works on >= 2.x.
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
                hint: "DIAMOND 2.1+ required; install via `apt install diamond-aligner`, \
                       `brew install diamond`, or `conda install -c bioconda diamond`"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = DiamondInput::from_case_dir(&case.path)?;

        fs::create_dir_all(workdir)?;

        // Resolve the query path (always an existing file: input
        // FASTA for blastp/blastx, source FASTA for makedb).
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
                    "[bio.diamond].query `{}` not found (resolved {})",
                    input.query.display(),
                    resolved_query.display()
                ),
            });
        }

        // The `database` field's role flips per action:
        // - blastp/blastx: an existing `.dmnd` DB to search against.
        // - makedb:        the *output* basename DIAMOND will write
        //                  `<basename>.dmnd` to. Don't probe for
        //                  existence here.
        // Round-9 hardening: relative `database` values flow into the
        // `-d` flag of DIAMOND; wrap with `confined_join` so a hostile
        // case can't aim makedb writes or blastp reads at arbitrary
        // paths. Absolute paths stay supported (admin-managed shared
        // DBs).
        let resolved_database: PathBuf = if input.database.is_absolute() {
            input.database.clone()
        } else {
            valenx_core::adapter_helpers::confined_join(&case.path, &input.database)?
        };
        if input.action == "blastp" || input.action == "blastx" {
            // For search modes DIAMOND accepts the DB path either with
            // or without the `.dmnd` suffix. Be lenient: accept either
            // a literal match or the suffixed form.
            let with_suffix = {
                let mut p = resolved_database.clone();
                let new_ext = p
                    .extension()
                    .and_then(|s| s.to_str())
                    .map(|s| {
                        if s == "dmnd" {
                            "dmnd".to_string()
                        } else {
                            format!("{s}.dmnd")
                        }
                    })
                    .unwrap_or_else(|| "dmnd".to_string());
                p.set_extension(new_ext);
                p
            };
            if !resolved_database.is_file() && !with_suffix.is_file() {
                return Err(AdapterError::InvalidCase {
                    case_path: case.path.join("case.toml"),
                    reason: format!(
                        "[bio.diamond].database `{}` not found (resolved {})",
                        input.database.display(),
                        resolved_database.display()
                    ),
                });
            }
        }

        let binary_path = find_on_path(BINARIES).ok_or_else(|| AdapterError::ToolNotInstalled {
            name: INFO_ID,
            hint: "DIAMOND 2.1+ required; install via `apt install diamond-aligner`, \
                       `brew install diamond`, or `conda install -c bioconda diamond`"
                .into(),
        })?;

        let native_command = build_command(
            &binary_path,
            &resolved_query,
            &resolved_database,
            &input,
            &case.path,
        )?;

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            // makedb on a UniRef-scale FASTA takes ~30 minutes;
            // searches scale with query / DB size. 4 hours is the
            // generous default shared with the rest of the bio
            // adapters.
            estimated_runtime: Some(Duration::from_secs(4 * 60 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting DIAMOND", |line| {
            let mut hint = subprocess::Hint::default();
            // DIAMOND prints "Total time = ..." at the end of every
            // run; "Reported N pairwise alignments" arrives just
            // before that for searches. Lift either marker as a
            // near-completion progress hint.
            if line.contains("Total time") {
                hint.progress = Some((95.0, line.to_string()));
            } else if line.contains("Reported ") || line.contains("Processing query block") {
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
        // Re-derive the action from the prepared command. The first
        // arg after the binary is always the subcommand.
        let action = job
            .native_command
            .get(1)
            .and_then(|s| s.to_str())
            .unwrap_or("")
            .to_string();

        // For blastp/blastx the artifact is the `-o <output>` hit
        // table; for makedb it's the `<database>.dmnd` we've just
        // written. Use the prepared command to recover both.
        let case_hash_input;
        let mut artefacts: Vec<Artifact> = Vec::new();
        match action.as_str() {
            "blastp" | "blastx" => {
                let output_path = output_after_flag(job, "-o");
                case_hash_input = output_path
                    .clone()
                    .filter(|p| p.is_file())
                    .unwrap_or_else(|| job.workdir.join("case.toml"));
                if let Some(out) = output_path {
                    if out.is_file() {
                        artefacts.push(Artifact {
                            path: out,
                            kind: ArtifactKind::Tabular,
                            checksum: None,
                            label: format!("DIAMOND {action} hits"),
                        });
                    }
                }
            }
            "makedb" => {
                // After `-d <basename>` DIAMOND writes
                // `<basename>.dmnd`. Recover the basename, append the
                // suffix, surface as a Native artifact.
                let dmnd_path = output_after_flag(job, "-d").map(|base| {
                    let mut p = base.clone();
                    let new_ext = p
                        .extension()
                        .and_then(|s| s.to_str())
                        .map(|s| {
                            if s == "dmnd" {
                                "dmnd".to_string()
                            } else {
                                format!("{s}.dmnd")
                            }
                        })
                        .unwrap_or_else(|| "dmnd".to_string());
                    p.set_extension(new_ext);
                    p
                });
                case_hash_input = dmnd_path
                    .clone()
                    .filter(|p| p.is_file())
                    .unwrap_or_else(|| job.workdir.join("case.toml"));
                if let Some(out) = dmnd_path {
                    if out.is_file() {
                        artefacts.push(Artifact {
                            path: out,
                            kind: ArtifactKind::Native,
                            checksum: None,
                            label: "DIAMOND .dmnd database".to_string(),
                        });
                    }
                }
            }
            _ => {
                case_hash_input = job.workdir.join("case.toml");
            }
        }

        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "DIAMOND",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);
        artefacts.sort_by(|a, b| a.path.cmp(&b.path));
        results.artifacts = artefacts;
        Ok(results)
    }

    fn capabilities(&self) -> Capabilities {
        Capabilities {
            capabilities: Vec::new(),
            ribbon_contributions: vec!["bio.diamond.search"],
        }
    }
}

/// Compose the `diamond` invocation for the given action.
///
/// Each action gets its own command shape:
///
/// - `blastp`/`blastx`: `diamond <action> -q <query> -d <database>
///                      -o <output> [--<sensitivity>] -p N [extras...]`
///   The `--<sensitivity>` flag is dropped when the user picked
///   `default` (DIAMOND has no `--default` flag).
/// - `makedb`: `diamond makedb --in <query> -d <database> -p N
///             [extras...]`. Here `query` is the *input* FASTA and
///   `database` is the *output* DB basename.
pub fn build_command(
    binary_path: &Path,
    resolved_query: &Path,
    resolved_database: &Path,
    case: &DiamondInput,
    case_path: &Path,
) -> Result<Vec<OsString>, AdapterError> {
    let mut cmd: Vec<OsString> = vec![
        binary_path.as_os_str().to_owned(),
        OsString::from(&case.action),
    ];

    match case.action.as_str() {
        "blastp" | "blastx" => {
            cmd.push(OsString::from("-q"));
            cmd.push(resolved_query.as_os_str().to_owned());
            cmd.push(OsString::from("-d"));
            cmd.push(resolved_database.as_os_str().to_owned());
            cmd.push(OsString::from("-o"));
            cmd.push(case.output.as_os_str().to_owned());
            // DIAMOND has no `--default` flag — the out-of-the-box
            // default is "no flag at all". Skip emitting the
            // sensitivity flag in that case.
            if case.sensitivity != "default" {
                cmd.push(OsString::from(format!("--{}", case.sensitivity)));
            }
            cmd.push(OsString::from("-p"));
            cmd.push(OsString::from(case.threads.to_string()));
            for arg in &case.extra_args {
                cmd.push(OsString::from(arg));
            }
        }
        "makedb" => {
            cmd.push(OsString::from("--in"));
            cmd.push(resolved_query.as_os_str().to_owned());
            cmd.push(OsString::from("-d"));
            cmd.push(resolved_database.as_os_str().to_owned());
            cmd.push(OsString::from("-p"));
            cmd.push(OsString::from(case.threads.to_string()));
            for arg in &case.extra_args {
                cmd.push(OsString::from(arg));
            }
        }
        // Defensive — case_input rejects unknown actions, but a
        // future schema-only change shouldn't reach a panic.
        other => {
            return Err(AdapterError::InvalidCase {
                case_path: case_path.join("case.toml"),
                reason: format!(
                    "internal: unknown diamond action `{other}` slipped past schema validation"
                ),
            });
        }
    }
    Ok(cmd)
}

/// Walk the prepared command for the value following `flag`. Used
/// from `collect()` to recover the `-o` / `-d` path so we can surface
/// the resulting artifact.
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
        let info = DiamondAdapter::new().info();
        assert_eq!(info.id, "diamond");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "GPL-3.0");
        assert_eq!(info.display_name, "DIAMOND");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = DiamondAdapter::new().info();
        assert_eq!(info.version_range.min_inclusive, Version::new(2, 1, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(3, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = DiamondAdapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.diamond.search"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = DiamondAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }

    /// Round-9 RED→GREEN: `[bio.diamond].database` used to be joined
    /// with bare `case.path.join`. Wrap with `confined_join` so a
    /// hostile case can't route makedb writes to `../../../usr/bin/diamond.dmnd`
    /// or read from `/etc/passwd.dmnd` siblings.
    #[test]
    fn prepare_rejects_database_traversing_outside_case_dir() {
        use valenx_test_utils::tempdir;
        let d = tempdir("diamond-db-trav");
        std::fs::write(d.join("query.fa"), b">x\nACGT\n").unwrap();
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "diamond.search"

[bio.diamond]
action   = "blastp"
query    = "query.fa"
database = "../../etc/passwd"
output   = "results.tsv"
"#,
        )
        .unwrap();
        let case = Case {
            id: "diamond-db-trav".into(),
            path: d.clone(),
        };
        let workdir = d.join("workdir");
        let err = DiamondAdapter::new().prepare(&case, &workdir).unwrap_err();
        let msg = format!("{err}");
        assert!(
            msg.contains("..") || msg.contains("stay within") || msg.contains("escape"),
            "expected confined_join rejection, got: {msg}"
        );
        let _ = std::fs::remove_dir_all(&d);
    }
}
