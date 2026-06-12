//! # valenx-adapter-lineardesign
//!
//! Adapter for [LinearDesign](https://github.com/LinearDesignSoftware/LinearDesign)
//! — Baidu Research's joint codon + secondary-structure mRNA design
//! tool. Given a target protein, LinearDesign co-optimises codon
//! adaptation index (CAI) against minimum free energy (MFE) of the
//! resulting mRNA's secondary structure, exposing a single tunable
//! `--lambda` knob that trades the two objectives off. The design
//! workhorse behind the modern mRNA-vaccine era — Apache-2.0,
//! actively maintained, sister to BWA / Smoldyn / NAMD in the bio
//! ecosystem.
//!
//! **Phase 43 — subprocess wrapper around `lineardesign`.** The user
//! supplies a protein FASTA via `[bio.lineardesign].protein` in
//! `case.toml`; `prepare()` composes
//! `lineardesign --aa <protein> --lambda <λ> --codon_usage <table>
//!  --output_basename <stem> [extras...]`, `run()` streams progress
//! via the shared subprocess runner, and `collect()` walks the workdir
//! for the `<basename>*.fasta` (designed mRNA), `<basename>*.txt`
//! (CAI/MFE report), and `*.log` outputs.

#![forbid(unsafe_code)]
#![allow(missing_docs)]

pub mod case_input;

use std::ffi::OsString;
use std::fs;
use std::path::Path;
use std::time::Duration;

use semver::Version;

use valenx_core::{
    adapter_helpers::{find_on_path, live_provenance},
    error::RunPhase,
    subprocess, Adapter, AdapterError, AdapterInfo, Capabilities, Case, LicenseMode, Physics,
    PreparedJob, ProbeReport, RunContext, RunReport, VersionRange,
};
use valenx_fields::{
    artifact::{Artifact, ArtifactKind},
    Results,
};

use crate::case_input::LinearDesignInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(LinearDesignAdapter::new())
}

pub struct LinearDesignAdapter;

impl LinearDesignAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for LinearDesignAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "lineardesign";
/// LinearDesign's binary candidate. The upstream repo ships a single
/// `lineardesign` driver; conda-forge / source builds expose the same
/// canonical lowercase name.
const BINARIES: &[&str] = &["lineardesign"];
/// Python-interpreter candidates probed only when `lineardesign` is
/// missing — surfaces a more useful "you have Python but not the
/// LinearDesign repo" hint for the common case where the user has
/// half-installed the tool.
const PYTHON_BINARIES: &[&str] = &["python3", "python"];

impl Adapter for LinearDesignAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "LinearDesign",
            // LinearDesign 1.0 is the canonical release line tagged in
            // the upstream repo (it tracks the published Nature paper).
            // Upper bound 2.0 reserves room for the next major bump.
            version_range: VersionRange {
                min_inclusive: Version::new(1, 0, 0),
                max_exclusive: Version::new(2, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            tool_license: "Apache-2.0",
            docs_url: "https://github.com/LinearDesignSoftware/LinearDesign",
            homepage_url: "https://github.com/LinearDesignSoftware/LinearDesign",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(BINARIES) {
            Some(binary_path) => Ok(ProbeReport {
                ok: true,
                // Upstream `lineardesign` doesn't expose a `--version`
                // banner; skip detection rather than emit a wrong
                // value. The version_range still gates compatibility
                // when the user supplies one explicitly.
                found_version: None,
                binary_path: Some(binary_path),
                warnings: Vec::new(),
                required_env: Vec::new(),
            }),
            None => {
                // If Python is on PATH but the LinearDesign driver
                // isn't, surface the half-installed-tool hint as a
                // warning rather than a hard ToolNotInstalled — the
                // user still benefits from validation flowing through
                // even when the binary's missing.
                if find_on_path(PYTHON_BINARIES).is_some() {
                    Ok(ProbeReport {
                        ok: false,
                        found_version: None,
                        binary_path: None,
                        warnings: vec!["LinearDesign not found on PATH; clone \
                             https://github.com/LinearDesignSoftware/LinearDesign \
                             and add the bin directory to PATH"
                            .into()],
                        required_env: Vec::new(),
                    })
                } else {
                    Err(AdapterError::ToolNotInstalled {
                        name: INFO_ID,
                        hint: "LinearDesign 1.0+ required; clone \
                               https://github.com/LinearDesignSoftware/LinearDesign \
                               and add the bin directory to PATH"
                            .into(),
                    })
                }
            }
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = LinearDesignInput::from_case_dir(&case.path)?;

        // Round-4 security: reject `output_basename = "../etc/passwd"`
        // and friends before the value flows into any path join.
        // Same pattern as the round-3 fix in bionetgen/iqtree/art/fasttree.
        valenx_core::adapter_helpers::validate_output_basename(
            &input.output_basename,
            "[bio.lineardesign].output_basename",
        )
        .map_err(|e| AdapterError::InvalidCase {
            case_path: case.path.join("case.toml"),
            reason: format!("{e}"),
        })?;

        fs::create_dir_all(workdir)?;

        // Resolve the protein FASTA against the case directory if
        // relative. LinearDesign reads it in place via `--aa <path>`;
        // we don't stage it into the workdir, just validate it exists
        // so the failure is fast and obvious.
        let source_protein = if input.protein.is_absolute() {
            input.protein.clone()
        } else {
            valenx_core::adapter_helpers::confined_join(&case.path, &input.protein)?
        };
        if !source_protein.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.lineardesign].protein `{}` not found (resolved {})",
                    input.protein.display(),
                    source_protein.display()
                ),
            });
        }

        let binary_path = find_on_path(BINARIES).ok_or_else(|| AdapterError::ToolNotInstalled {
            name: INFO_ID,
            hint: "LinearDesign 1.0+ required; clone \
                       https://github.com/LinearDesignSoftware/LinearDesign \
                       and add the bin directory to PATH"
                .into(),
        })?;

        // Compose `lineardesign --aa <protein> --lambda <λ>
        //          --codon_usage <table> --output_basename <stem>
        //          [extras...]`. Each `--flag value` is two separate
        // OsString args, not `--flag=value` — LinearDesign's argparse
        // accepts both but the split form survives quoting better
        // when arguments contain spaces or shell metachars.
        let mut native_command: Vec<OsString> = vec![
            binary_path.into_os_string(),
            OsString::from("--aa"),
            source_protein.into_os_string(),
            OsString::from("--lambda"),
            OsString::from(format!("{}", input.lambda_param)),
            OsString::from("--codon_usage"),
            OsString::from(&input.codon_usage),
            OsString::from("--output_basename"),
            OsString::from(&input.output_basename),
        ];
        for arg in &input.extra_args {
            native_command.push(OsString::from(arg));
        }

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            // A single short ORF (~300 aa) finishes in under a minute;
            // long full-spike-protein designs with restrictive lambda
            // settings can run for many minutes. 30 minutes covers
            // the long tail.
            estimated_runtime: Some(Duration::from_secs(30 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting LinearDesign", |line| {
            let mut hint = subprocess::Hint::default();
            // LinearDesign emits a startup banner ("LinearDesign"),
            // per-stage progress ("CAI: ...", "MFE: ..."), and a
            // sentinel near the end of the run.
            if line.contains("Design complete") || line.contains("design complete") {
                hint.progress = Some((95.0, line.to_string()));
            } else if line.contains("LinearDesign") || line.contains("lineardesign") {
                hint.progress = Some((5.0, line.to_string()));
            } else if line.contains("CAI") || line.contains("MFE") {
                hint.progress = Some((50.0, line.to_string()));
            } else if line.contains("ERROR") || line.contains("Error:") {
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
        // Provenance: hash the case.toml as the canonical input
        // descriptor. LinearDesign's outputs are named off
        // `output_basename`, so we re-parse the case to filter
        // typed outputs to that stem.
        let case_hash_input = job.workdir.join("case.toml");
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "lineardesign",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);
        let mut artefacts: Vec<Artifact> = Vec::new();

        // The basename comes from the staged case.toml in the workdir
        // when present (the orchestrator copies it there), or we just
        // accept everything matching the typed extensions otherwise.
        let basename = read_output_basename(&job.workdir);

        let entries = match fs::read_dir(&job.workdir) {
            Ok(e) => e,
            Err(e) => {
                tracing::warn!(target: "valenx-lineardesign", ?e, "workdir read failed");
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
            let stem = path
                .file_stem()
                .and_then(|s| s.to_str())
                .unwrap_or("")
                .to_string();
            let stem_matches_basename = match basename.as_deref() {
                Some(b) => stem.starts_with(b),
                None => true,
            };
            match ext.as_deref() {
                Some("fasta") => {
                    if !stem_matches_basename {
                        continue;
                    }
                    artefacts.push(Artifact {
                        path,
                        kind: ArtifactKind::Native,
                        checksum: None,
                        label: "LinearDesign optimized mRNA".to_string(),
                    });
                }
                Some("txt") => {
                    if !stem_matches_basename {
                        continue;
                    }
                    artefacts.push(Artifact {
                        path,
                        kind: ArtifactKind::Tabular,
                        checksum: None,
                        label: "LinearDesign report".to_string(),
                    });
                }
                Some("log") => {
                    artefacts.push(Artifact {
                        path,
                        kind: ArtifactKind::Log,
                        checksum: None,
                        label: "LinearDesign log".to_string(),
                    });
                }
                _ => continue,
            }
        }
        artefacts.sort_by(|a, b| a.path.cmp(&b.path));
        results.artifacts = artefacts;
        Ok(results)
    }

    fn capabilities(&self) -> Capabilities {
        Capabilities {
            capabilities: Vec::new(),
            ribbon_contributions: vec!["bio.lineardesign.design"],
        }
    }
}

/// Pull `output_basename` out of the staged `case.toml` for
/// `collect()`-time output filtering. Returns `None` if the case isn't
/// staged or the field is missing — `collect()` then accepts every
/// matching extension regardless of stem.
fn read_output_basename(workdir: &Path) -> Option<String> {
    // Round-23 sweep: bound staged case.toml at MAX_PROJECT_FILE_BYTES.
    let text = valenx_core::io_caps::read_capped_to_string(
        &workdir.join("case.toml"),
        valenx_core::project::loader::MAX_PROJECT_FILE_BYTES as usize,
    )
    .ok()?;
    let parsed: toml::Value = toml::from_str(&text).ok()?;
    parsed
        .get("bio")?
        .get("lineardesign")?
        .get("output_basename")?
        .as_str()
        .map(|s| s.to_string())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn info_is_bio_domain() {
        let info = LinearDesignAdapter::new().info();
        assert_eq!(info.id, "lineardesign");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "Apache-2.0");
        assert_eq!(info.display_name, "LinearDesign");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = LinearDesignAdapter::new().info();
        // LinearDesign 1.0 is the canonical release line tagged in
        // the upstream repo; 2.0 reserves room for the next major.
        assert_eq!(info.version_range.min_inclusive, Version::new(1, 0, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(2, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = LinearDesignAdapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.lineardesign.design"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = LinearDesignAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }
}
