//! # valenx-adapter-igv
//!
//! Adapter for [IGV](https://igv.org/) — the Integrative Genomics
//! Viewer from the Broad Institute. The interactive desktop viewer
//! is a separate animal; this adapter wraps `igvtools`, IGV's
//! headless CLI multitool, the way `samtools` wraps htslib.
//! `igvtools` indexes BAM / VCF files, generates TDF density
//! tracks, sorts BAM / VCF, and produces tile-server inputs — the
//! file-prep layer that feeds an IGV web or desktop session.
//!
//! **Phase 23 — multi-action subprocess wrapper around
//! `igvtools <action>`.** Sister adapter to samtools / bcftools:
//! per-action dispatch on `case.action` via [`build_command`], shared
//! [`valenx_core::subprocess::run`] runner. Like the bcftools
//! post-fix shape, `build_command` returns
//! `Result<Vec<OsString>, AdapterError>` and surfaces an
//! `InvalidCase` for any unknown action — never panics — so a
//! schema drift can't bubble up as a process abort.
//!
//! Each-action command:
//!
//! - `index`: `igvtools index <input>` — writes `<input>.bai`
//!   (BAM) or `<input>.idx` (VCF) sidecar next to the input. No
//!   `output` flag; `collect()` walks the input's directory for
//!   the produced sidecar.
//! - `count`: `igvtools count -w <window_size> <input> <output>
//!   [extras...]` — TDF density file. `collect()` reports the
//!   `<output>` path.
//! - `sort` : `igvtools sort <input> <output> [extras...]` — sorted
//!   SAM / BAM / VCF. `collect()` reports `<output>`.
//! - `tile` : `igvtools tile <input> <output> [extras...]` — TDF
//!   tile from a coverage track. `collect()` reports `<output>`.

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

use crate::case_input::IgvInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(IgvAdapter::new())
}

pub struct IgvAdapter;

impl IgvAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for IgvAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "igv";
/// igvtools binary candidates. Bioconda, the IGV distribution, and
/// most package managers install the launcher under the canonical
/// lowercase name.
const BINARIES: &[&str] = &["igvtools"];

impl Adapter for IgvAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "IGV",
            // IGV / igvtools 2.16 (early 2023) is the floor we test
            // against — it carries the modernised `igvtools count` /
            // `tile` flag set we lean on. The upper bound 3.0
            // reserves room for the eventual major bump.
            version_range: VersionRange {
                min_inclusive: Version::new(2, 16, 0),
                max_exclusive: Version::new(3, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            tool_license: "MIT",
            docs_url: "https://software.broadinstitute.org/software/igv/",
            homepage_url: "https://igv.org/",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(BINARIES) {
            Some(binary_path) => {
                // `igvtools version` prints the version string;
                // `--version` is also accepted on recent builds.
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
                hint: "igvtools 2.16+ required; install via \
                       `conda install -c bioconda igvtools` or download from \
                       https://igv.org/doc/desktop/#DownloadPage/"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = IgvInput::from_case_dir(&case.path)?;

        fs::create_dir_all(workdir)?;

        // Resolve the input path against the case directory if
        // relative. Same convention as every other Phase 17/18 bio
        // adapter — `input = "aligned.bam"` next to `case.toml`.
        let source_input = if input.input.is_absolute() {
            input.input.clone()
        } else {
            valenx_core::adapter_helpers::confined_join(
            &case.path,
            &input.input,
        )?
        };
        if !source_input.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.igv].input `{}` not found (resolved {})",
                    input.input.display(),
                    source_input.display()
                ),
            });
        }

        let binary_path = find_on_path(BINARIES).ok_or_else(|| AdapterError::ToolNotInstalled {
            name: INFO_ID,
            hint: "igvtools 2.16+ required; install via \
                       `conda install -c bioconda igvtools` or download from \
                       https://igv.org/doc/desktop/#DownloadPage/"
                .into(),
        })?;

        let native_command = build_command(&binary_path, &source_input, &input, &case.path)?;

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            // igvtools indexing is fast (seconds for typical BAMs);
            // count / tile on whole-genome BAMs can run for tens of
            // minutes. 1 hour is a generous default.
            estimated_runtime: Some(Duration::from_secs(60 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting igvtools", |line| {
            let mut hint = subprocess::Hint::default();
            // igvtools prints sparse progress on stdout — most
            // useful markers are "Done." for completion and
            // "Indexing" / "Counting" for action progress.
            if line.contains("Done.") || line.contains("Wrote") {
                hint.progress = Some((95.0, line.to_string()));
            } else if line.contains("Indexing")
                || line.contains("Counting")
                || line.contains("Sorting")
                || line.contains("Tiling")
            {
                hint.progress = Some((50.0, line.to_string()));
            } else if line.contains("ERROR") || line.contains("Exception") {
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
        // Re-derive the action from the prepared command (same
        // technique as samtools). The first arg after the binary is
        // always the subcommand name.
        let action = job
            .native_command
            .get(1)
            .and_then(|s| s.to_str())
            .unwrap_or("");

        // Provenance: hash whichever output the action produces.
        // Falls back to case.toml when the run hasn't produced
        // anything yet — keeps the provenance block well-formed for
        // partial / failed runs.
        let case_hash_input = primary_output_path(job, action)
            .filter(|p| p.is_file())
            .unwrap_or_else(|| job.workdir.join("case.toml"));
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "igvtools",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);
        let mut artefacts: Vec<Artifact> = Vec::new();

        match action {
            "index" => {
                // `igvtools index <input>` writes a sidecar next to
                // the input, NOT in the workdir:
                //   - BAM → `<input>.bai`
                //   - VCF → `<input>.idx`
                //   - some inputs → `<basename>.idx` next to input
                // Look for both extensions so the user gets the
                // correct one back regardless of input format.
                if let Some(input_path) = positional_input(job) {
                    for ext in &["bai", "idx"] {
                        let sidecar = with_extra_extension(&input_path, ext);
                        if sidecar.is_file() {
                            artefacts.push(Artifact {
                                path: sidecar,
                                kind: ArtifactKind::Native,
                                checksum: None,
                                label: format!("igvtools index ({ext})"),
                            });
                        }
                    }
                }
            }
            "count" | "tile" => {
                // count / tile produce a TDF density file at the
                // positional output path (last positional after the
                // input). Surface as Native — TDF is IGV's binary
                // format.
                if let Some(out) = positional_output(job) {
                    if out.is_file() {
                        artefacts.push(Artifact {
                            path: out.clone(),
                            kind: ArtifactKind::Native,
                            checksum: None,
                            label: format!("igvtools {action} output"),
                        });
                    }
                }
            }
            "sort" => {
                // sort writes a sorted file at the positional output
                // path. Output format mirrors the input — SAM stays
                // SAM, BAM stays BAM, VCF stays VCF. Surface as
                // Tabular when the extension is `.sam` / `.vcf` and
                // Native otherwise (BAM / BCF binary formats).
                if let Some(out) = positional_output(job) {
                    if out.is_file() {
                        let kind = match out
                            .extension()
                            .and_then(|s| s.to_str())
                            .map(|s| s.to_ascii_lowercase())
                            .as_deref()
                        {
                            Some("sam") | Some("vcf") => ArtifactKind::Tabular,
                            _ => ArtifactKind::Native,
                        };
                        artefacts.push(Artifact {
                            path: out.clone(),
                            kind,
                            checksum: None,
                            label: format!("igvtools sort output ({})", out.display()),
                        });
                    }
                }
            }
            _ => {}
        }

        artefacts.sort_by(|a, b| a.path.cmp(&b.path));
        results.artifacts = artefacts;
        Ok(results)
    }

    fn capabilities(&self) -> Capabilities {
        // Bio-specific Capability variants land in a follow-up task;
        // ribbon contributions are already enough for the registry to
        // surface the adapter without crashing the UI's
        // capability-index builder.
        Capabilities {
            capabilities: Vec::new(),
            ribbon_contributions: vec!["bio.igv.index"],
        }
    }
}

/// Compose the `igvtools` invocation for the given action.
///
/// Returns the command vector or an `InvalidCase` error if the action
/// slipped past the schema-level validation in
/// [`IgvInput::from_case_dir`]. Mirrors the bcftools post-fix shape:
/// soft-error on unknown action, never panic.
///
/// Each action gets its own command shape:
///
/// - `index`: `igvtools index <input>`
/// - `count`: `igvtools count -w <window_size> <input> <output> [extras...]`
/// - `sort` : `igvtools sort <input> <output> [extras...]`
/// - `tile` : `igvtools tile <input> <output> [extras...]`
pub fn build_command(
    binary_path: &Path,
    source_input: &Path,
    case: &IgvInput,
    case_path: &Path,
) -> Result<Vec<OsString>, AdapterError> {
    let mut cmd: Vec<OsString> = vec![
        binary_path.as_os_str().to_owned(),
        OsString::from(&case.action),
    ];

    match case.action.as_str() {
        "index" => {
            // `igvtools index <input>` — sidecar lands next to the
            // input automatically. No output flag.
            cmd.push(source_input.as_os_str().to_owned());
            // `index` has no extras path in the typical flow, but
            // pass through anyway for forward-compat.
            for arg in &case.extra_args {
                cmd.push(OsString::from(arg));
            }
        }
        "count" => {
            // case_input enforces output.is_some() for count.
            let output = case
                .output
                .as_ref()
                .expect("case_input enforces output for count");
            cmd.push(OsString::from("-w"));
            cmd.push(OsString::from(case.window_size.to_string()));
            cmd.push(source_input.as_os_str().to_owned());
            cmd.push(output.as_os_str().to_owned());
            for arg in &case.extra_args {
                cmd.push(OsString::from(arg));
            }
        }
        "sort" => {
            let output = case
                .output
                .as_ref()
                .expect("case_input enforces output for sort");
            cmd.push(source_input.as_os_str().to_owned());
            cmd.push(output.as_os_str().to_owned());
            for arg in &case.extra_args {
                cmd.push(OsString::from(arg));
            }
        }
        "tile" => {
            let output = case
                .output
                .as_ref()
                .expect("case_input enforces output for tile");
            cmd.push(source_input.as_os_str().to_owned());
            cmd.push(output.as_os_str().to_owned());
            for arg in &case.extra_args {
                cmd.push(OsString::from(arg));
            }
        }
        // Soft-error on unknown action. case_input rejects unknown
        // actions during parse, but if a future schema-only change
        // adds an entry without a corresponding command shape,
        // surface an InvalidCase rather than panic.
        other => {
            return Err(AdapterError::InvalidCase {
                case_path: case_path.join("case.toml"),
                reason: format!(
                    "internal: unknown igvtools action `{other}` slipped past schema validation"
                ),
            });
        }
    }
    Ok(cmd)
}

/// Locate the positional input in the prepared command. For every
/// action we wrap, the input is the first non-flag positional after
/// the subcommand. `count` puts `-w <n>` before the input, so we
/// peek/skip flag-value pairs.
fn positional_input(job: &PreparedJob) -> Option<PathBuf> {
    let mut iter = job.native_command.iter().skip(2).peekable();
    while let Some(arg) = iter.next() {
        let s = arg.to_str().unwrap_or("");
        if s.starts_with('-') {
            // Skip the value of any known flag-with-value pair.
            // For our shapes, `-w` is the only flag that takes a
            // value before the input.
            if matches!(s, "-w") {
                let _ = iter.next();
            }
            continue;
        }
        return Some(PathBuf::from(arg));
    }
    None
}

/// Locate the positional output in the prepared command. For
/// `count` / `sort` / `tile` it's the second non-flag positional
/// after the subcommand (immediately after the input).
fn positional_output(job: &PreparedJob) -> Option<PathBuf> {
    let mut iter = job.native_command.iter().skip(2).peekable();
    let mut positional_seen = 0;
    while let Some(arg) = iter.next() {
        let s = arg.to_str().unwrap_or("");
        if s.starts_with('-') {
            if matches!(s, "-w") {
                let _ = iter.next();
            }
            continue;
        }
        positional_seen += 1;
        if positional_seen == 2 {
            return Some(PathBuf::from(arg));
        }
    }
    None
}

/// Pick the primary output path to hash for provenance, given the
/// action.
///
/// - `index`: the produced `.bai` / `.idx` sidecar next to the input.
/// - `count` / `sort` / `tile`: the explicit positional output.
fn primary_output_path(job: &PreparedJob, action: &str) -> Option<PathBuf> {
    match action {
        "index" => {
            // Pick whichever sidecar exists next to the input. `.bai`
            // first (BAM is the more common path); fall back to
            // `.idx` for VCF.
            let input = positional_input(job)?;
            let bai = with_extra_extension(&input, "bai");
            if bai.is_file() {
                return Some(bai);
            }
            let idx = with_extra_extension(&input, "idx");
            if idx.is_file() {
                return Some(idx);
            }
            // Return `.bai` as the canonical guess so live_provenance
            // still gets a deterministic path even if the run didn't
            // produce anything.
            Some(bai)
        }
        "count" | "sort" | "tile" => positional_output(job),
        _ => None,
    }
}

/// Append `.<ext>` to a path. Used to derive the `.bai` / `.idx`
/// sidecar paths for the `index` action.
fn with_extra_extension(path: &Path, ext: &str) -> PathBuf {
    let mut s = path.as_os_str().to_owned();
    s.push(".");
    s.push(ext);
    PathBuf::from(s)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn info_is_bio_domain() {
        let info = IgvAdapter::new().info();
        assert_eq!(info.id, "igv");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "MIT");
        assert_eq!(info.display_name, "IGV");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = IgvAdapter::new().info();
        // igvtools >= 2.16 (modernised count / tile flags); upper
        // bound 3.0 reserves room for the next major.
        assert_eq!(info.version_range.min_inclusive, Version::new(2, 16, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(3, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = IgvAdapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.igv.index"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = IgvAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }
}
