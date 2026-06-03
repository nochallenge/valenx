//! # valenx-adapter-bcftools
//!
//! Adapter for [bcftools](https://samtools.github.io/bcftools/) — the
//! htslib companion tool for VCF / BCF data, the standard variant-call
//! format multitool that ships alongside samtools. This adapter wraps
//! the four most common subcommands: `view` (VCF↔BCF and subsetting),
//! `call` (multiallelic-caller variant calling from a BAM), `filter`
//! (soft / hard filters on existing variants), and `concat`
//! (concatenate per-region or per-chunk VCFs back into one).
//!
//! **Phase 19 — subprocess wrapper around `bcftools <action>`.** The
//! user picks the subcommand via `action` in `[bio.bcftools]`. Every
//! action writes its output to a file via `-o <output>`, so the shared
//! [`valenx_core::subprocess::run`] runner is enough — there's no
//! MAFFT-style stdout-capture detour like the samtools `flagstat`
//! action needs.
//!
//! Each-action command construction is factored into [`build_command`]
//! so the dispatch is one self-contained match without polluting the
//! Adapter::prepare body.

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

use crate::case_input::BcftoolsInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(BcftoolsAdapter::new())
}

pub struct BcftoolsAdapter;

impl BcftoolsAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for BcftoolsAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "bcftools";
/// bcftools' binary candidates. Bioconda, Homebrew, and source builds
/// all install under the canonical name.
const BINARIES: &[&str] = &["bcftools"];

impl Adapter for BcftoolsAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "bcftools",
            // bcftools 1.x has tracked samtools' release line; 1.17
            // (2023) is the floor we test against — it carries the
            // modernised CLI we recommend (consistent `--threads`,
            // `-O <type>` for output format, `-o <file>` everywhere).
            // The upper bound 2.0 reserves room for an eventual major
            // bump.
            version_range: VersionRange {
                min_inclusive: Version::new(1, 17, 0),
                max_exclusive: Version::new(2, 0, 0),
            },
            physics: &[Physics::Bio],
            // htslib (and bcftools with it) is MIT/Expat. The htslib
            // README tracks the exact terms.
            license_mode: LicenseMode::Subprocess,
            tool_license: "MIT",
            docs_url: "http://www.htslib.org/doc/bcftools.html",
            homepage_url: "https://www.htslib.org/",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(BINARIES) {
            Some(binary_path) => {
                // `bcftools --version` prints "bcftools 1.17 ..." on
                // stdout; the combined scanner picks it up cleanly.
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
                hint: "bcftools 1.17+ required; install via `apt install bcftools`, \
                       `brew install bcftools`, or `conda install -c bioconda bcftools`"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = BcftoolsInput::from_case_dir(&case.path)?;

        fs::create_dir_all(workdir)?;

        // Resolve every input path against the case directory via
        // `confined_join` — rejects absolute paths and `..` traversal
        // out of the case sandbox. Same threat model as every other
        // Phase 17/18 bio adapter: a shared-case bundle should not be
        // able to point any path-shaped field (input, inputs[],
        // reference) at `/etc/passwd`.
        //
        // Round-8 sibling-field sweep: pre-fix, only the *first* path
        // field happened to use confined_join in some adapters; the
        // others fell through to plain `case.path.join`. Every
        // user-supplied path field now goes through confined_join.
        let resolved_input: Option<PathBuf> = match &input.input {
            Some(p) => {
                let resolved = confined_join(&case.path, p)?;
                if !resolved.is_file() {
                    return Err(AdapterError::InvalidCase {
                        case_path: case.path.join("case.toml"),
                        reason: format!(
                            "[bio.bcftools].input `{}` not found (resolved {})",
                            p.display(),
                            resolved.display()
                        ),
                    });
                }
                Some(resolved)
            }
            None => None,
        };

        let mut resolved_inputs: Vec<PathBuf> = Vec::with_capacity(input.inputs.len());
        for entry in &input.inputs {
            let resolved = confined_join(&case.path, entry)?;
            if !resolved.is_file() {
                return Err(AdapterError::InvalidCase {
                    case_path: case.path.join("case.toml"),
                    reason: format!(
                        "[bio.bcftools].inputs entry `{}` not found (resolved {})",
                        entry.display(),
                        resolved.display()
                    ),
                });
            }
            resolved_inputs.push(resolved);
        }

        let resolved_reference: Option<PathBuf> = match &input.reference {
            Some(p) => {
                let resolved = confined_join(&case.path, p)?;
                if !resolved.is_file() {
                    return Err(AdapterError::InvalidCase {
                        case_path: case.path.join("case.toml"),
                        reason: format!(
                            "[bio.bcftools].reference `{}` not found (resolved {})",
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
            hint: "bcftools 1.17+ required; install via `apt install bcftools`, \
                       `brew install bcftools`, or `conda install -c bioconda bcftools`"
                .into(),
        })?;

        let native_command = build_command(
            &binary_path,
            resolved_input.as_deref(),
            &resolved_inputs,
            resolved_reference.as_deref(),
            &input,
            &case.path,
        )?;

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            // Variant calling on a single chromosome can run for
            // hours; whole-genome calls run longer still. 4 hours is
            // the same generous default as samtools/bwa.
            estimated_runtime: Some(Duration::from_secs(4 * 60 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting bcftools", |line| {
            let mut hint = subprocess::Hint::default();
            // bcftools prints sparse progress on stderr: "[E::...]"
            // markers for hard errors, "Lines total/split-ok/skipped"
            // summaries for `call` and `concat`. The shared runner
            // routes stderr through Warn-level logging on its own;
            // this stdout handler just lifts a couple of well-known
            // markers if they appear.
            if line.contains("Lines total") || line.contains("Real time:") {
                hint.progress = Some((95.0, line.to_string()));
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
        // Re-derive the action and output path from the prepared
        // command. The first arg after the binary is always the
        // subcommand name; the value after `-o` is always the output
        // path.
        let action = job
            .native_command
            .get(1)
            .and_then(|s| s.to_str())
            .unwrap_or("");
        let output_path = output_after_flag(job, "-o");

        // Provenance: hash the produced output if present. Falls back
        // to case.toml when the run hasn't produced anything yet —
        // keeps the provenance block well-formed for partial / failed
        // runs.
        let case_hash_input = output_path
            .clone()
            .filter(|p| p.is_file())
            .unwrap_or_else(|| job.workdir.join("case.toml"));
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "bcftools",
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
                    path: out.clone(),
                    kind: ArtifactKind::Tabular,
                    checksum: None,
                    label: format!("bcftools {action} output"),
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
            ribbon_contributions: vec!["bio.bcftools.call"],
        }
    }
}

/// Compose the `bcftools` invocation for the given action.
///
/// Each action gets its own command shape:
///
/// - `view`:   `bcftools view -O v --threads N -o <output> <input> [extras...]`
/// - `call`:   `bcftools call --threads N -O v -m -v -f <reference>
///             -o <output> <input> [extras...]`
/// - `filter`: `bcftools filter --threads N -O v -o <output> <input> [extras...]`
/// - `concat`: `bcftools concat --threads N -O v -o <output> <inputs...> [extras...]`
pub fn build_command(
    binary_path: &Path,
    resolved_input: Option<&Path>,
    resolved_inputs: &[PathBuf],
    resolved_reference: Option<&Path>,
    case: &BcftoolsInput,
    case_path: &Path,
) -> Result<Vec<OsString>, AdapterError> {
    let mut cmd: Vec<OsString> = vec![
        binary_path.as_os_str().to_owned(),
        OsString::from(&case.action),
    ];

    match case.action.as_str() {
        "view" => {
            cmd.push(OsString::from("-O"));
            cmd.push(OsString::from("v"));
            cmd.push(OsString::from("--threads"));
            cmd.push(OsString::from(case.threads.to_string()));
            cmd.push(OsString::from("-o"));
            cmd.push(case.output.as_os_str().to_owned());
            // case_input enforces input.is_some() for view.
            cmd.push(
                resolved_input
                    .expect("case_input enforces input for view")
                    .as_os_str()
                    .to_owned(),
            );
            for arg in &case.extra_args {
                cmd.push(OsString::from(arg));
            }
        }
        "call" => {
            cmd.push(OsString::from("--threads"));
            cmd.push(OsString::from(case.threads.to_string()));
            cmd.push(OsString::from("-O"));
            cmd.push(OsString::from("v"));
            // `-m` selects the multiallelic caller (the modern
            // default), `-v` keeps only variant sites — the standard
            // pair for "give me a clean variant-only VCF from this
            // BAM."
            cmd.push(OsString::from("-m"));
            cmd.push(OsString::from("-v"));
            cmd.push(OsString::from("-f"));
            cmd.push(
                resolved_reference
                    .expect("case_input enforces reference for call")
                    .as_os_str()
                    .to_owned(),
            );
            cmd.push(OsString::from("-o"));
            cmd.push(case.output.as_os_str().to_owned());
            cmd.push(
                resolved_input
                    .expect("case_input enforces input for call")
                    .as_os_str()
                    .to_owned(),
            );
            for arg in &case.extra_args {
                cmd.push(OsString::from(arg));
            }
        }
        "filter" => {
            cmd.push(OsString::from("--threads"));
            cmd.push(OsString::from(case.threads.to_string()));
            cmd.push(OsString::from("-O"));
            cmd.push(OsString::from("v"));
            cmd.push(OsString::from("-o"));
            cmd.push(case.output.as_os_str().to_owned());
            cmd.push(
                resolved_input
                    .expect("case_input enforces input for filter")
                    .as_os_str()
                    .to_owned(),
            );
            for arg in &case.extra_args {
                cmd.push(OsString::from(arg));
            }
        }
        "concat" => {
            cmd.push(OsString::from("--threads"));
            cmd.push(OsString::from(case.threads.to_string()));
            cmd.push(OsString::from("-O"));
            cmd.push(OsString::from("v"));
            cmd.push(OsString::from("-o"));
            cmd.push(case.output.as_os_str().to_owned());
            for entry in resolved_inputs {
                cmd.push(entry.as_os_str().to_owned());
            }
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
                    "internal: unknown bcftools action `{other}` slipped past schema validation"
                ),
            });
        }
    }
    Ok(cmd)
}

/// Walk the prepared command for the value following `flag`. Used
/// from `collect()` to recover the `-o` output path so we can surface
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
        let info = BcftoolsAdapter::new().info();
        assert_eq!(info.id, "bcftools");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "MIT");
        assert_eq!(info.display_name, "bcftools");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = BcftoolsAdapter::new().info();
        // 1.17 is the floor we test against; 2.0 reserves room for an
        // eventual major bump.
        assert_eq!(info.version_range.min_inclusive, Version::new(1, 17, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(2, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = BcftoolsAdapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.bcftools.call"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = BcftoolsAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }

    #[test]
    fn prepare_rejects_traversal_reference_path() {
        // Round-8 RED→GREEN: every path-shaped field in [bio.bcftools]
        // (input, inputs[], reference) now routes through
        // `confined_join`, not the bare `case.path.join(absolute)`
        // pattern that pre-fix happily returned `/etc/passwd` verbatim.
        // This test exercises `reference` — one representative sibling
        // path covers all three.
        //
        // We use `../etc/passwd` (relative traversal) rather than an
        // absolute `/etc/passwd` so the test works the same on Windows,
        // where Rust's `Path::is_absolute()` returns `false` for
        // `/etc/passwd`. confined_join rejects `..` components on
        // every platform.
        use valenx_test_utils::tempdir;
        let d = tempdir("bcftools-traversal");
        std::fs::write(d.join("aligned.bam"), b"placeholder").unwrap();
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "bcftools.call"

[bio.bcftools]
action    = "call"
input     = "aligned.bam"
output    = "calls.vcf"
reference = "../etc/passwd"
"#,
        )
        .unwrap();
        let case = Case {
            id: "bcftools-traversal".into(),
            path: d.clone(),
        };
        let workdir = d.join("workdir");
        let err = BcftoolsAdapter::new()
            .prepare(&case, &workdir)
            .unwrap_err();
        let msg = format!("{err}");
        assert!(
            msg.contains("absolute") || msg.contains("escape") || msg.contains("`..`") || msg.contains("traversal"),
            "expected confined_join rejection, got: {msg}"
        );
        let _ = std::fs::remove_dir_all(&d);
    }
}
