//! # valenx-adapter-dssr
//!
//! Adapter for [DSSR](https://x3dna.org/dssr/) (Dissecting the
//! Spatial Structure of RNA) — Xiang-Jun Lu's reference toolkit for
//! DNA / RNA structural-feature annotation. DSSR reads a
//! nucleic-acid PDB and emits a single JSON file enumerating every
//! detected feature: base pairs (Watson-Crick, Hoogsteen, sugar-edge,
//! ...), multiplets, double helices, stems, hairpin / internal /
//! junction loops, kissing loops, A-minor motifs, ribose zippers,
//! pseudoknots, splayed-apart conformations, and more. It is the
//! standard machine-readable feature-extraction step in modern
//! RNA-structure pipelines.
//!
//! **Phase 39 — subprocess wrapper around `x3dna-dssr`.** The user
//! supplies an input `.pdb` and an `output_json` path via
//! `[bio.dssr]` in `case.toml`. `prepare()` composes a
//! `x3dna-dssr -i=<input_pdb> -o=<output_json> --json [extras...]`
//! invocation; `run()` streams the run via the shared subprocess
//! runner.
//!
//! On `collect()` we report the JSON output as a Tabular artefact
//! (DSSR's JSON is the canonical machine-readable summary; we tag
//! it Tabular rather than Native so downstream serdes can key off a
//! consistent kind).
//!
//! ## License flag
//!
//! DSSR ships under the X3DNA family non-OSS license — non-commercial
//! / academic use only. We surface this via
//! `tool_license = "DSSR-License"` and emit a probe warning when
//! the binary is found, with the literal string `"academic"` as a
//! stable anchor for tests and downstream license-aware filters.

#![forbid(unsafe_code)]
#![allow(missing_docs)]

pub mod case_input;

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

use crate::case_input::DssrInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(DssrAdapter::new())
}

pub struct DssrAdapter;

impl DssrAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for DssrAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "dssr";
/// DSSR's binary candidate. The canonical entry point shipped by
/// the X3DNA family is `x3dna-dssr`.
const BINARIES: &[&str] = &["x3dna-dssr"];

/// The probe-warning surfaced whenever DSSR is detected. Anchors a
/// stable "academic / non-commercial only" reminder for downstream
/// tooling and tests; the literal string `"academic"` is part of
/// the asserted contract.
const LICENSE_WARNING: &str = "DSSR is licensed for non-commercial / academic use only. \
     Confirm your use case complies with the DSSR license before \
     redistributing analyses or derived data.";

impl Adapter for DssrAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "DSSR",
            // DSSR 2.x is the modern stable line that ships with
            // X3DNA 2.4+. Floor 2.0; upper bound 3.0 reserves room
            // for an eventual major bump.
            version_range: VersionRange {
                min_inclusive: Version::new(2, 0, 0),
                max_exclusive: Version::new(3, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            // DSSR's terms aren't a recognised SPDX identifier; the
            // closest accurate label is the project's own
            // "DSSR-License" name.
            tool_license: "DSSR-License",
            docs_url: "https://x3dna.org/dssr/",
            homepage_url: "https://x3dna.org/dssr/",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(BINARIES) {
            Some(binary_path) => {
                // `x3dna-dssr --version` prints a banner with the
                // DSSR release on stdout. The generic detector tries
                // common version flags.
                let found_version = detect_tool_version_semver(&binary_path, &["--version", ""]);
                Ok(ProbeReport {
                    ok: true,
                    found_version,
                    binary_path: Some(binary_path),
                    // Always surface the license reminder when DSSR
                    // is detected — non-OSS academic use only.
                    warnings: vec![LICENSE_WARNING.to_string()],
                    required_env: Vec::new(),
                })
            }
            None => Err(AdapterError::ToolNotInstalled {
                name: INFO_ID,
                hint: "DSSR 2.0+ required; download from \
                       https://x3dna.org/dssr/ (registration required, \
                       academic-use license) and ensure `x3dna-dssr` \
                       is on PATH"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = DssrInput::from_case_dir(&case.path)?;

        // Round-10 H3: `output_json` is `PathBuf` and pre-fix flowed
        // into `workdir.join(&input.output_json)`. Validate as a
        // basename — DSSR writes a single JSON file, not a directory.
        if let Some(s) = input.output_json.to_str() {
            valenx_core::adapter_helpers::validate_output_basename(
                s,
                "[bio.dssr].output_json",
            )
            .map_err(|e| AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!("{e}"),
            })?;
        } else {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: "[bio.dssr].output_json: non-UTF-8 path rejected".into(),
            });
        }

        fs::create_dir_all(workdir)?;

        // Resolve the input PDB against the case directory if
        // relative.
        let source_pdb = if input.input_pdb.is_absolute() {
            input.input_pdb.clone()
        } else {
            valenx_core::adapter_helpers::confined_join(
            &case.path,
            &input.input_pdb,
        )?
        };
        if !source_pdb.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.dssr].input_pdb `{}` not found (resolved {})",
                    input.input_pdb.display(),
                    source_pdb.display()
                ),
            });
        }

        // Output paths are scoped to the workdir if relative.
        let output_json = if input.output_json.is_absolute() {
            input.output_json.clone()
        } else {
            workdir.join(&input.output_json)
        };

        let binary_path = find_on_path(BINARIES).ok_or_else(|| AdapterError::ToolNotInstalled {
            name: INFO_ID,
            hint: "DSSR 2.0+ required; download from \
                       https://x3dna.org/dssr/ (registration required, \
                       academic-use license) and ensure `x3dna-dssr` \
                       is on PATH"
                .into(),
        })?;

        // Compose `x3dna-dssr -i=<input_pdb> -o=<output_json> --json [extras...]`.
        // DSSR uses `key=value` flag form (no space between flag and
        // value) on its short-form options; we emit a single argv
        // entry per flag to match.
        let mut native_command: Vec<OsString> = vec![
            binary_path.into_os_string(),
            OsString::from(format!("-i={}", source_pdb.display())),
            OsString::from(format!("-o={}", output_json.display())),
            OsString::from("--json"),
        ];
        for arg in &input.extra_args {
            native_command.push(OsString::from(arg));
        }

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            // DSSR on a single structure runs in seconds; long-tail
            // is large RNAs / multi-model NMR ensembles, which still
            // finish well inside an hour.
            estimated_runtime: Some(Duration::from_secs(60 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting DSSR", |line| {
            let mut hint = subprocess::Hint::default();
            // DSSR's progress chatter is loose; surface the rare
            // structural anchors so the UI shows forward motion.
            if line.contains("Number of base pairs") || line.contains("Total number of nucleotides")
            {
                hint.progress = Some((50.0, line.to_string()));
            } else if line.contains("Time used") || line.contains("Done") {
                hint.progress = Some((95.0, line.to_string()));
            }
            if line.contains("ERROR") || line.contains("error:") || line.contains("FATAL") {
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
        // Provenance: hash the JSON output if produced, else
        // case.toml so the provenance block stays well-formed for
        // partial / failed runs.
        let case_hash_input = {
            let json = read_output_json_path(&job.workdir);
            match json {
                Some(p) if p.is_file() => p,
                _ => job.workdir.join("case.toml"),
            }
        };
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "DSSR",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);
        let mut artefacts: Vec<Artifact> = Vec::new();

        // Re-derive the JSON output path from the staged case.toml
        // so we report the exact file DSSR was told to write.
        // Failure to read the case.toml is non-fatal — we just
        // skip the artefact entry.
        if let Some(output_json) = read_output_json_path(&job.workdir) {
            if output_json.is_file() {
                artefacts.push(Artifact {
                    path: output_json,
                    kind: ArtifactKind::Tabular,
                    checksum: None,
                    label: "DSSR analysis (JSON)".to_string(),
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
            ribbon_contributions: vec!["bio.dssr.analyze"],
        }
    }
}

/// Pull the JSON output path out of the staged `case.toml` so
/// `collect()` can report it. Returns the resolved (workdir-rooted
/// when relative) path, or None if the file doesn't exist or can't
/// be parsed.
fn read_output_json_path(workdir: &Path) -> Option<std::path::PathBuf> {
    let case_toml = workdir.join("case.toml");
    let text = valenx_core::io_caps::read_capped_to_string(&case_toml, valenx_core::project::loader::MAX_PROJECT_FILE_BYTES as usize).ok()?;
    let parsed: toml::Value = toml::from_str(&text).ok()?;
    let raw = parsed
        .get("bio")?
        .get("dssr")?
        .get("output_json")?
        .as_str()?
        .to_string();
    let p = std::path::PathBuf::from(raw);
    if p.is_absolute() {
        Some(p)
    } else {
        Some(workdir.join(p))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn info_is_bio_domain() {
        let info = DssrAdapter::new().info();
        assert_eq!(info.id, "dssr");
        assert_eq!(info.physics, &[Physics::Bio]);
        // DSSR's custom non-OSS license, not a recognised SPDX
        // identifier — pin the project's own label.
        assert_eq!(info.tool_license, "DSSR-License");
        assert_eq!(info.display_name, "DSSR");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = DssrAdapter::new().info();
        // DSSR 2.0 is the floor; upper bound 3.0 reserves room for
        // the next major.
        assert_eq!(info.version_range.min_inclusive, Version::new(2, 0, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(3, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = DssrAdapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.dssr.analyze"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = DssrAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }

    #[test]
    fn probe_warning_mentions_academic() {
        // The license-flag warning is mandatory: DSSR is non-OSS
        // academic-use, and we surface that on every successful
        // probe. The literal "academic" anchor is what downstream
        // tooling and license-aware filters key off — pin it.
        assert!(
            LICENSE_WARNING.contains("academic"),
            "probe warning must contain `academic` anchor; got: {LICENSE_WARNING}"
        );
    }

    /// Round-10 H3 RED→GREEN: `output_json` flowed into
    /// `workdir.join(...)` with no validation. Hostile
    /// `output_json = "../etc/passwd"` is now rejected.
    #[test]
    fn prepare_rejects_output_json_path_traversal() {
        use valenx_test_utils::tempdir;
        let d = tempdir("dssr-output-trav");
        std::fs::write(d.join("rna.pdb"), b"FAKE\n").unwrap();
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "dssr.analyze"

[bio.dssr]
input_pdb   = "rna.pdb"
output_json = "../etc/passwd"
"#,
        )
        .unwrap();
        let case = Case {
            id: "trav".into(),
            path: d.clone(),
        };
        let workdir = d.join("workdir");
        let err = DssrAdapter::new()
            .prepare(&case, &workdir)
            .unwrap_err();
        let msg = format!("{err}");
        assert!(
            msg.contains("[bio.dssr].output_json"),
            "expected [bio.dssr].output_json in error, got: {msg}"
        );
        let _ = std::fs::remove_dir_all(&d);
    }
}
