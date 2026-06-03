//! # valenx-adapter-x3dna
//!
//! Adapter for [X3DNA](https://x3dna.org/) — Wilma Olson and Xiang-Jun
//! Lu's reference toolkit for DNA / RNA structural-geometry
//! analysis. X3DNA reads a nucleic-acid PDB, identifies base pairs,
//! and computes the canonical helical-step parameters (twist, roll,
//! tilt, slide, shift, rise) plus per-base intra-pair parameters
//! (buckle, propeller, opening, shear, stretch, stagger). It is the
//! workhorse behind structural-bioinformatics pipelines that need
//! quantitative DNA geometry — bending studies, drug-DNA / protein-DNA
//! complex analysis, RNA tertiary-structure annotation.
//!
//! **Phase 39 — subprocess wrapper around the `analyze` binary.** The
//! user supplies an input `.pdb` and an `output_basename` via
//! `[bio.x3dna]` in `case.toml`. `prepare()` composes a
//! `analyze <input_pdb> [extras...]` invocation; `run()` streams the
//! run via the shared subprocess runner.
//!
//! On `collect()` we walk the workdir for `<output_basename>*.par`
//! (the base-step parameter table) and `*.out` (the per-run log
//! `analyze` writes alongside).
//!
//! ## License flag
//!
//! X3DNA ships under a custom non-OSS license that restricts use to
//! non-commercial / academic contexts. We surface this accurately
//! via `tool_license = "X3DNA-License"` and emit a probe warning
//! when the binary is found, with the literal string `"academic"`
//! as a stable anchor for tests and downstream license-aware
//! filters.

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

use crate::case_input::X3dnaInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(X3dnaAdapter::new())
}

pub struct X3dnaAdapter;

impl X3dnaAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for X3dnaAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "x3dna";
/// X3DNA's binary candidate. The canonical entry point is the
/// `analyze` driver shipped with the X3DNA distribution; it must be
/// on PATH after the user runs `x3dna_setup` from the install dir.
const BINARIES: &[&str] = &["analyze"];

/// The probe-warning surfaced whenever X3DNA is detected. Anchors a
/// stable "academic / non-commercial only" reminder for downstream
/// tooling and tests; the literal string `"academic"` is part of the
/// asserted contract.
const LICENSE_WARNING: &str = "X3DNA is licensed for non-commercial / academic use only. \
     Confirm your use case complies with the X3DNA license before \
     redistributing analyses or derived data.";

impl Adapter for X3dnaAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "X3DNA",
            // X3DNA 2.4 (2020) is the modern stable release and the
            // floor we test against; upper bound 3.0 reserves room
            // for an eventual major bump.
            version_range: VersionRange {
                min_inclusive: Version::new(2, 4, 0),
                max_exclusive: Version::new(3, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            // X3DNA's terms aren't a recognised SPDX identifier;
            // the closest accurate label is the project's own
            // "X3DNA-License" name.
            tool_license: "X3DNA-License",
            docs_url: "https://x3dna.org/",
            homepage_url: "https://x3dna.org/",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(BINARIES) {
            Some(binary_path) => {
                // `analyze` (no args) prints a usage banner with
                // the X3DNA release on stdout / stderr; the generic
                // detector tries both `--version` and a bare-name
                // scan to lift the version where possible.
                let found_version = detect_tool_version_semver(&binary_path, &["--version", ""]);
                Ok(ProbeReport {
                    ok: true,
                    found_version,
                    binary_path: Some(binary_path),
                    // Always surface the license reminder when
                    // X3DNA is detected — non-OSS academic use only.
                    warnings: vec![LICENSE_WARNING.to_string()],
                    required_env: Vec::new(),
                })
            }
            None => Err(AdapterError::ToolNotInstalled {
                name: INFO_ID,
                hint: "X3DNA 2.4+ required; download from \
                       https://x3dna.org/ (registration required, \
                       academic-use license) and run `x3dna_setup` \
                       so `analyze` is on PATH"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = X3dnaInput::from_case_dir(&case.path)?;

        // Round-4 security: reject `output_basename = "../etc/passwd"`
        // and friends before the value flows into any path join.
        // Same pattern as the round-3 fix in bionetgen/iqtree/art/fasttree.
        valenx_core::adapter_helpers::validate_output_basename(
            &input.output_basename,
            "[bio.x3dna].output_basename",
        )
        .map_err(|e| AdapterError::InvalidCase {
            case_path: case.path.join("case.toml"),
            reason: format!("{e}"),
        })?;

        fs::create_dir_all(workdir)?;

        // Resolve the input PDB against the case directory if
        // relative. Same convention as every other Phase 17/18 bio
        // adapter — `input_pdb = "structure.pdb"` next to `case.toml`.
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
                    "[bio.x3dna].input_pdb `{}` not found (resolved {})",
                    input.input_pdb.display(),
                    source_pdb.display()
                ),
            });
        }

        let binary_path = find_on_path(BINARIES).ok_or_else(|| AdapterError::ToolNotInstalled {
            name: INFO_ID,
            hint: "X3DNA 2.4+ required; download from \
                       https://x3dna.org/ (registration required, \
                       academic-use license) and run `x3dna_setup` \
                       so `analyze` is on PATH"
                .into(),
        })?;

        // Compose `analyze <input_pdb> [extras...]`.
        // `analyze` is positional-only — it derives every output
        // filename from the input basename, so we just hand it the
        // PDB and any user-supplied extras.
        let mut native_command: Vec<OsString> =
            vec![binary_path.into_os_string(), source_pdb.into_os_string()];
        for arg in &input.extra_args {
            native_command.push(OsString::from(arg));
        }

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            // Single-structure analysis runs in seconds; long-tail
            // is large RNAs / multi-model NMR ensembles, which still
            // finish well inside an hour.
            estimated_runtime: Some(Duration::from_secs(60 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting X3DNA analyze", |line| {
            let mut hint = subprocess::Hint::default();
            // analyze's progress chatter is loose; surface the rare
            // structural anchors so the UI shows forward motion.
            if line.contains("Number of base-pairs") || line.contains("Total number") {
                hint.progress = Some((50.0, line.to_string()));
            } else if line.contains("complete") || line.contains("Done") {
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
        // Provenance: hash the case.toml as the canonical input
        // descriptor. analyze derives output filenames from the
        // input PDB basename, so we can't pin a single fixed-name
        // artifact for the prov hash.
        let case_hash_input = job.workdir.join("case.toml");
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "X3DNA",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);
        let mut artefacts: Vec<Artifact> = Vec::new();

        // Read the staged case.toml back out so we can restrict the
        // collected `.par` files to those whose stem starts with the
        // configured `output_basename`. Failure to read the case is
        // non-fatal — we then accept every `.par` as a candidate.
        let basename = read_output_basename(&job.workdir);

        let entries = match fs::read_dir(&job.workdir) {
            Ok(e) => e,
            Err(e) => {
                tracing::warn!(target: "valenx-x3dna", ?e, "workdir read failed");
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
            match ext.as_deref() {
                Some("par") => {
                    let stem = path
                        .file_stem()
                        .and_then(|s| s.to_str())
                        .unwrap_or("")
                        .to_string();
                    let stem_ok = match basename.as_deref() {
                        Some(b) => stem.starts_with(b),
                        None => true,
                    };
                    if !stem_ok {
                        continue;
                    }
                    artefacts.push(Artifact {
                        path,
                        kind: ArtifactKind::Tabular,
                        checksum: None,
                        label: "X3DNA base-step parameters".to_string(),
                    });
                }
                Some("out") => {
                    artefacts.push(Artifact {
                        path,
                        kind: ArtifactKind::Log,
                        checksum: None,
                        label: "X3DNA log".to_string(),
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
            ribbon_contributions: vec!["bio.x3dna.analyze"],
        }
    }
}

/// Pull `output_basename` out of the staged `case.toml` for
/// `collect()`-time `.par`-stem filtering. Returns None if the file
/// doesn't exist or can't be parsed — collect falls back to
/// accepting every `.par` in that case.
fn read_output_basename(workdir: &Path) -> Option<String> {
    let case_toml = workdir.join("case.toml");
    let text = valenx_core::io_caps::read_capped_to_string(&case_toml, valenx_core::project::loader::MAX_PROJECT_FILE_BYTES as usize).ok()?;
    let parsed: toml::Value = toml::from_str(&text).ok()?;
    parsed
        .get("bio")?
        .get("x3dna")?
        .get("output_basename")?
        .as_str()
        .map(|s| s.to_string())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn info_is_bio_domain() {
        let info = X3dnaAdapter::new().info();
        assert_eq!(info.id, "x3dna");
        assert_eq!(info.physics, &[Physics::Bio]);
        // X3DNA's custom non-OSS license, not a recognised SPDX
        // identifier — pin the project's own label.
        assert_eq!(info.tool_license, "X3DNA-License");
        assert_eq!(info.display_name, "X3DNA");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = X3dnaAdapter::new().info();
        // X3DNA 2.4 (2020) is the floor; upper bound 3.0 reserves
        // room for the next major.
        assert_eq!(info.version_range.min_inclusive, Version::new(2, 4, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(3, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = X3dnaAdapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.x3dna.analyze"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = X3dnaAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }

    #[test]
    fn probe_warning_mentions_academic() {
        // The license-flag warning is mandatory: X3DNA is non-OSS
        // academic-use, and we surface that on every successful
        // probe. The literal "academic" anchor is what downstream
        // tooling and license-aware filters key off — pin it.
        assert!(
            LICENSE_WARNING.contains("academic"),
            "probe warning must contain `academic` anchor; got: {LICENSE_WARNING}"
        );
    }
}
