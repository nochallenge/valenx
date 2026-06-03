//! # valenx-adapter-bionetgen
//!
//! Adapter for [BioNetGen](https://bionetgen.org/) — the rule-based
//! modeling language and tool suite for combinatorially-complex
//! signaling networks. The user writes BNGL (BioNetGen Language)
//! files describing molecular species, sites, and reaction *rules*,
//! and BioNetGen expands the rules into the underlying reaction
//! network and (optionally) integrates it deterministically (ODE) or
//! stochastically (SSA).
//!
//! **Phase 32 — subprocess wrapper around `BNG2.pl`.** `BNG2.pl` is
//! the canonical Perl driver: it reads a `.bngl` model, executes the
//! actions inside (`generate_network`, `simulate`, `simulate_ssa`,
//! `parameter_scan`), and writes a stack of fixed-suffix outputs
//! prefixed by `-o <output_basename>`. The user supplies the model,
//! the basename, and an optional `generate_only = true` flag (which
//! adds `--no-execute`, skipping simulate actions to emit just the
//! expanded reaction network).
//!
//! On `collect()` we walk the workdir for the canonical
//! `<basename>*.net` (reaction network), `<basename>*.gdat` (species
//! trajectories), and `<basename>*.cdat` (concentrations) outputs.

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

use crate::case_input::BioNetGenInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(BioNetGenAdapter::new())
}

pub struct BioNetGenAdapter;

impl BioNetGenAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for BioNetGenAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "bionetgen";
/// BioNetGen's canonical entry point. Perl-based; `BNG2.pl` lives on
/// PATH after a stock conda or source install.
const BINARIES: &[&str] = &["BNG2.pl"];

impl Adapter for BioNetGenAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "BioNetGen",
            // BioNetGen 2.8+ is the recent stable line; 3.0 reserves
            // room for an eventual major bump (NFsim integration is
            // tracking 2.x).
            version_range: VersionRange {
                min_inclusive: Version::new(2, 8, 0),
                max_exclusive: Version::new(3, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            tool_license: "MIT",
            docs_url: "https://bionetgen.org/index.php/Documentation",
            homepage_url: "https://bionetgen.org/",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(BINARIES) {
            Some(binary_path) => {
                // `BNG2.pl --version` prints something like
                // "BioNetGen version 2.8.5" on stdout / stderr.
                let found_version = detect_tool_version_semver(&binary_path, &["--version", ""]);
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
                hint: "BioNetGen 2.8+ required; install via \
                       `conda install -c bioconda bionetgen` or download \
                       a release from https://github.com/RuleWorld/bionetgen/releases"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = BioNetGenInput::from_case_dir(&case.path)?;

        // Round-3 security fix: `output_basename` flows into the
        // workdir path that BNG2.pl prefixes onto every output file. A
        // hostile case.toml setting it to `"../../etc/cron.d/x"` would
        // otherwise let it write outside the workdir.
        valenx_core::adapter_helpers::validate_output_basename(
            &input.output_basename,
            "[bio.bionetgen].output_basename",
        )
        .map_err(|e| AdapterError::InvalidCase {
            case_path: case.path.join("case.toml"),
            reason: format!("{e}"),
        })?;

        fs::create_dir_all(workdir)?;

        // Stage `case.toml` into the workdir so collect() can recover
        // the configured `output_basename` for prefix-filtering output
        // artifacts. Without this stage, the basename filter silently
        // degrades to "match everything".
        let staged_case_toml = workdir.join("case.toml");
        let source_case_toml = case.path.join("case.toml");
        if source_case_toml.is_file() {
            fs::copy(&source_case_toml, &staged_case_toml)
                .map_err(|e| AdapterError::Other(anyhow::anyhow!("stage case.toml: {e}")))?;
        }

        // Resolve the model path against the case directory if
        // relative. Same convention as every other Phase 17/18 bio
        // adapter — `model = "egfr.bngl"` next to `case.toml`.
        let source_model = if input.model.is_absolute() {
            input.model.clone()
        } else {
            valenx_core::adapter_helpers::confined_join(
            &case.path,
            &input.model,
        )?
        };
        if !source_model.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.bionetgen].model `{}` not found (resolved {})",
                    input.model.display(),
                    source_model.display()
                ),
            });
        }

        let binary_path = find_on_path(BINARIES).ok_or_else(|| AdapterError::ToolNotInstalled {
            name: INFO_ID,
            hint: "BioNetGen 2.8+ required; install via \
                       `conda install -c bioconda bionetgen` or download \
                       a release from https://github.com/RuleWorld/bionetgen/releases"
                .into(),
        })?;

        // Compose `BNG2.pl [--no-execute] -o <basename> <model> [extras...]`.
        // `--no-execute` skips simulate / scan / fitting actions and
        // emits just the expanded reaction network. `-o` pins every
        // output file's prefix so collect() walks deterministically.
        //
        // Round-3 fix: extras intentionally come AFTER the positional
        // `<model>` so a hostile case.toml can't slip a phantom
        // positional in via `extra_args` and shift `<model>` onto a
        // different argument slot.
        let mut native_command: Vec<OsString> = vec![binary_path.into_os_string()];
        if input.generate_only {
            native_command.push(OsString::from("--no-execute"));
        }
        native_command.push(OsString::from("-o"));
        native_command.push(OsString::from(&input.output_basename));
        native_command.push(source_model.into_os_string());
        for arg in &input.extra_args {
            native_command.push(OsString::from(arg));
        }

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            // Network generation on a small EGFR-style model finishes
            // in seconds; a full `parameter_scan` over a stiff signaling
            // network can run for hours. 4 hours covers the long tail.
            estimated_runtime: Some(Duration::from_secs(4 * 60 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting BioNetGen", |line| {
            let mut hint = subprocess::Hint::default();
            // BNG2.pl emits one-liners per action: "Read file <path>",
            // "ACTION: generate_network(...)", "Wrote network in
            // <T> CPU s.", "ACTION: simulate(...)" etc. Lift the
            // obvious milestones to coarse UI ticks.
            if line.contains("Wrote network") || line.contains("simulation completed") {
                hint.progress = Some((75.0, line.to_string()));
            } else if line.starts_with("ACTION:") {
                hint.progress = Some((25.0, line.to_string()));
            } else if line.contains("ABORT") || line.contains("ERROR") || line.contains("error:") {
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
        // Re-read the staged input so we can match output files by
        // basename. Failure is non-fatal — collect still surfaces
        // every recognised extension under a generic label.
        let basename: Option<String> = BioNetGenInput::from_case_dir(&job.workdir)
            .ok()
            .map(|i| i.output_basename);

        // Provenance: hash the canonical `<basename>.net` if present,
        // else case.toml so the provenance block stays well-formed
        // even for partial / failed runs.
        let case_hash_input: PathBuf = match &basename {
            Some(b) => {
                let net = job.workdir.join(format!("{b}.net"));
                if net.is_file() {
                    net
                } else {
                    job.workdir.join("case.toml")
                }
            }
            None => job.workdir.join("case.toml"),
        };
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "BioNetGen",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);
        let mut artefacts: Vec<Artifact> = Vec::new();

        // Walk the workdir top-level. BNG2.pl writes a stack of
        // `<basename>...` files: `<basename>.net` (reaction network),
        // `<basename>.gdat` (species trajectories), `<basename>.cdat`
        // (concentrations); `parameter_scan` emits per-trial variants
        // sharing the basename prefix (e.g. `<basename>_001.gdat`).
        let entries = match fs::read_dir(&job.workdir) {
            Ok(e) => e,
            Err(e) => {
                tracing::warn!(target: "valenx-bionetgen", ?e, "workdir read failed");
                return Ok(results);
            }
        };
        for entry in entries.flatten() {
            let path = entry.path();
            if !path.is_file() {
                continue;
            }
            let name = path
                .file_name()
                .and_then(|s| s.to_str())
                .map(|s| s.to_string());
            let Some(name) = name else { continue };

            // Restrict to files starting with the basename when we
            // know it; otherwise accept by extension alone.
            if let Some(b) = &basename {
                if !name.starts_with(b) {
                    continue;
                }
            }
            let ext = path
                .extension()
                .and_then(|s| s.to_str())
                .map(|s| s.to_ascii_lowercase());
            let (kind, label) = match ext.as_deref() {
                Some("net") => (
                    ArtifactKind::Native,
                    "BioNetGen reaction network".to_string(),
                ),
                Some("gdat") => (
                    ArtifactKind::Tabular,
                    "BioNetGen species trajectories".to_string(),
                ),
                Some("cdat") => (
                    ArtifactKind::Tabular,
                    "BioNetGen concentrations".to_string(),
                ),
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
        Capabilities {
            capabilities: Vec::new(),
            ribbon_contributions: vec!["bio.bionetgen.simulate"],
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn info_is_bio_domain() {
        let info = BioNetGenAdapter::new().info();
        assert_eq!(info.id, "bionetgen");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "MIT");
        assert_eq!(info.display_name, "BioNetGen");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = BioNetGenAdapter::new().info();
        assert_eq!(info.version_range.min_inclusive, Version::new(2, 8, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(3, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = BioNetGenAdapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.bionetgen.simulate"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = BioNetGenAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }

    /// Round-3 security fix: `output_basename` must be validated before
    /// being appended to the workdir path. Without the check, a hostile
    /// `output_basename = "../../etc/cron.d/x"` would let BNG2.pl write
    /// outside the workdir.
    #[test]
    fn prepare_rejects_traversal_in_output_basename() {
        let case_dir = std::env::temp_dir().join(format!(
            "valenx-bng-trav-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(&case_dir).unwrap();
        // Stage the model file so we get past the model-exists check.
        fs::write(case_dir.join("model.bngl"), b"# placeholder\n").unwrap();
        fs::write(
            case_dir.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "bionetgen.simulate"

[bio.bionetgen]
model = "model.bngl"
output_basename = "../../etc/cron.d/x"
"#,
        )
        .unwrap();
        let workdir = case_dir.join("work");
        let case = Case {
            id: "bng-trav".into(),
            path: case_dir.clone(),
        };
        let err = BioNetGenAdapter::new()
            .prepare(&case, &workdir)
            .expect_err("traversal output_basename must be rejected");
        match err {
            AdapterError::InvalidCase { reason, .. } => {
                let msg = reason;
                assert!(
                    msg.contains("output_basename") && (msg.contains("..") || msg.contains("separators")),
                    "expected traversal-rejection message, got: {msg}"
                );
            }
            other => panic!("unexpected error: {other:?}"),
        }
        let _ = fs::remove_dir_all(&case_dir);
    }
}
