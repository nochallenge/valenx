//! # valenx-adapter-physicell
//!
//! Adapter for [PhysiCell](https://physicell.org/) — Paul Macklin's
//! agent-based, off-lattice multicellular simulator. PhysiCell models
//! tens to hundreds of thousands of individual cells (each an agent
//! with state, mechanics, secretion, and phenotype) coupled to a
//! reaction-diffusion microenvironment for substrates like oxygen and
//! drugs. The canonical use case is tumour growth and immunology.
//!
//! **Phase 32 — subprocess wrapper around the per-project binary.**
//! Unlike a typical CLI tool, PhysiCell models compile to a project-
//! specific C++ executable: the user clones the framework, edits the
//! project's `custom_modules/` source, runs `make`, and ends up with
//! e.g. `./project` next to the project directory. The adapter
//! therefore takes both a `binary` path and the run-time XML
//! configuration via `[bio.physicell]` in `case.toml`.
//!
//! `probe()` checks `find_on_path(&["physicell"])` (which most installs
//! won't have — that's expected) and returns `ok: true` either way,
//! attaching a warning that PhysiCell models compile per-project. The
//! real validation happens in `prepare()`, which checks the user's
//! `binary` exists on disk before composing the command.
//!
//! On `collect()` we walk the `output/` subdirectory PhysiCell writes
//! its run snapshots to: `*.xml` and `*.mat` are SVG/MAT-format tissue
//! snapshots, `*.csv` are scalar summaries.

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

use crate::case_input::PhysiCellInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(PhysiCellAdapter::new())
}

pub struct PhysiCellAdapter;

impl PhysiCellAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for PhysiCellAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "physicell";
/// Optional PATH lookup. Most installs won't ship a generic
/// `physicell` binary — projects compile per-model — but if the user
/// has wrapped one we'll find it.
const BINARIES: &[&str] = &["physicell"];

impl Adapter for PhysiCellAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "PhysiCell",
            // PhysiCell 1.13+ is the recent stable line (1.14 ships as
            // of 2024 with the `cell-rules` DSL); 2.0 reserves room
            // for an eventual major bump.
            version_range: VersionRange {
                min_inclusive: Version::new(1, 13, 0),
                max_exclusive: Version::new(2, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            tool_license: "BSD-3-Clause",
            docs_url: "https://physicell.org/documentation/",
            homepage_url: "https://physicell.org/",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        // PhysiCell models compile per-project — there's no canonical
        // binary on PATH for most installs. Try the lookup but report
        // ok: true either way, attaching a warning so the user knows
        // we're going to validate `binary` at prepare time instead.
        match find_on_path(BINARIES) {
            Some(binary_path) => {
                let found_version = detect_tool_version_semver(&binary_path, &["--version", ""]);
                Ok(ProbeReport {
                    ok: true,
                    found_version,
                    binary_path: Some(binary_path),
                    warnings: Vec::new(),
                    required_env: Vec::new(),
                })
            }
            None => Ok(ProbeReport {
                ok: true,
                found_version: None,
                binary_path: None,
                warnings: vec![
                    "PhysiCell models compile per-project; no generic `physicell` \
                     binary on PATH is expected. The adapter validates the \
                     [bio.physicell].binary field at prepare time."
                        .to_string(),
                ],
                required_env: Vec::new(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = PhysiCellInput::from_case_dir(&case.path)?;

        fs::create_dir_all(workdir)?;

        // Round-3 security fix: `binary` flows straight into
        // `Command::new`, so a hostile case.toml could otherwise point
        // it at e.g. `/usr/bin/curl` and turn "Run case" into arbitrary
        // exec. PhysiCell models compile to a per-project binary, so
        // we can't allow-list a name — but we CAN demand the binary
        // live inside the case directory, which means the attacker
        // needs write access to a project's case folder before they
        // can take this path (at which point they already control the
        // project).
        let source_binary = if input.binary.is_absolute() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.physicell].binary `{}` must be a relative path \
                     inside the case directory, not an absolute path. \
                     PhysiCell models compile per-project — copy or symlink \
                     the compiled `./project` binary next to `case.toml`.",
                    input.binary.display()
                ),
            });
        } else if input
            .binary
            .components()
            .any(|c| matches!(c, std::path::Component::ParentDir))
        {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.physicell].binary `{}` must not traverse out of the \
                     case directory via `..` — keep the compiled binary \
                     inside the case folder.",
                    input.binary.display()
                ),
            });
        } else {
            // Round-9 classification: KEEP `case.path.join` here —
            // the two branches above already reject absolute paths
            // and `..` traversal explicitly with adapter-specific
            // messaging, so this is equivalent to `confined_join` for
            // this field.
            case.path.join(&input.binary)
        };
        if !source_binary.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.physicell].binary `{}` not found (resolved {}). \
                     PhysiCell models compile per-project — clone the \
                     framework, edit the project's `custom_modules/` source, \
                     run `make`, and point this field at the resulting \
                     executable.",
                    input.binary.display(),
                    source_binary.display()
                ),
            });
        }

        let source_config = if input.config.is_absolute() {
            input.config.clone()
        } else {
            valenx_core::adapter_helpers::confined_join(
            &case.path,
            &input.config,
        )?
        };
        if !source_config.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.physicell].config `{}` not found (resolved {})",
                    input.config.display(),
                    source_config.display()
                ),
            });
        }

        // Compose `<binary> <config> [extras...]`. PhysiCell binaries
        // accept the XML settings file as a positional argument.
        let mut native_command: Vec<OsString> = vec![
            source_binary.into_os_string(),
            source_config.into_os_string(),
        ];
        for arg in &input.extra_args {
            native_command.push(OsString::from(arg));
        }

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            // 2D demos finish in minutes; 3D tumour simulations with
            // tens of thousands of agents and reaction-diffusion
            // microenvironments routinely run for hours. 8 hours
            // covers the long tail without being absurd.
            estimated_runtime: Some(Duration::from_secs(8 * 60 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting PhysiCell", |line| {
            let mut hint = subprocess::Hint::default();
            // PhysiCell prints periodic "current simulated time:" /
            // "Total runtime:" lines on stdout. Lift the obvious
            // milestones to coarse UI ticks.
            if line.contains("Total runtime") || line.contains("simulation finished") {
                hint.progress = Some((95.0, line.to_string()));
            } else if line.contains("current simulated time") || line.contains("saving SVG output")
            {
                hint.progress = Some((50.0, line.to_string()));
            } else if line.to_ascii_lowercase().contains("error") || line.contains("ABORT") {
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
        let output_dir = job.workdir.join("output");

        // Provenance: hash the run's XML snapshot manifest if present
        // (PhysiCell writes `output/initial.xml` early); else
        // case.toml so the provenance block stays well-formed for
        // partial / failed runs.
        let case_hash_input = {
            let initial = output_dir.join("initial.xml");
            if initial.is_file() {
                initial
            } else {
                job.workdir.join("case.toml")
            }
        };
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "PhysiCell",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);
        let mut artefacts: Vec<Artifact> = Vec::new();

        // Walk `output/`. PhysiCell drops a stack of per-snapshot
        // files there: `output<N>.xml` (manifest), `output<N>_*.mat`
        // (cell + microenvironment state in MATLAB v4 binary), and
        // optional `*.csv` scalar summaries.
        let entries = match fs::read_dir(&output_dir) {
            Ok(e) => e,
            Err(e) => {
                tracing::warn!(
                    target: "valenx-physicell",
                    ?e,
                    "output dir read failed (`{}` missing?)",
                    output_dir.display()
                );
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
                Some("xml") | Some("mat") => (
                    ArtifactKind::Native,
                    "PhysiCell tissue snapshot".to_string(),
                ),
                Some("csv") => (ArtifactKind::Tabular, "PhysiCell scalar table".to_string()),
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
            ribbon_contributions: vec!["bio.physicell.simulate"],
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn info_is_bio_domain() {
        let info = PhysiCellAdapter::new().info();
        assert_eq!(info.id, "physicell");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "BSD-3-Clause");
        assert_eq!(info.display_name, "PhysiCell");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = PhysiCellAdapter::new().info();
        assert_eq!(info.version_range.min_inclusive, Version::new(1, 13, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(2, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = PhysiCellAdapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.physicell.simulate"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = PhysiCellAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }

    /// Round-3 security fix: an absolute binary path in case.toml
    /// would otherwise let a hostile case point `binary = "/usr/bin/curl"`
    /// and turn "Run case" into arbitrary code execution.
    #[test]
    fn prepare_rejects_absolute_binary_path() {
        let case_dir = std::env::temp_dir().join(format!(
            "valenx-physicell-absbin-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&case_dir).unwrap();
        // Cross-platform "absolute path that exists": on Unix use
        // /bin/sh, on Windows use C:\Windows\System32\cmd.exe — but
        // we only need the path to PARSE as absolute, so we pick a
        // platform-appropriate string. The validation rejects before
        // it checks file existence.
        let abs_binary = if cfg!(windows) {
            r#"binary = "C:/Windows/System32/cmd.exe""#
        } else {
            r#"binary = "/usr/bin/curl""#
        };
        std::fs::write(
            case_dir.join("case.toml"),
            format!(
                r#"[case]
physics = "bio"
solver  = "physicell.simulate"

[bio.physicell]
{abs_binary}
config = "config/PhysiCell_settings.xml"
"#
            ),
        )
        .unwrap();
        let workdir = case_dir.join("workdir");
        let case = Case {
            id: "absbin".into(),
            path: case_dir.clone(),
        };
        let result = PhysiCellAdapter::new().prepare(&case, &workdir);
        let err = result.expect_err("absolute binary path should be rejected");
        let msg = format!("{err}");
        assert!(
            msg.contains("absolute") || msg.contains("relative"),
            "expected absolute-path rejection message, got: {msg}"
        );
        let _ = std::fs::remove_dir_all(&case_dir);
    }

    /// Round-3 fix: `binary = "../../etc/passwd"` (or any `..` traversal)
    /// must also be rejected — case-scoping is meaningless if we just
    /// resolve relative against case.path.
    #[test]
    fn prepare_rejects_parent_dir_traversal_in_binary() {
        let case_dir = std::env::temp_dir().join(format!(
            "valenx-physicell-traversal-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&case_dir).unwrap();
        std::fs::write(
            case_dir.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "physicell.simulate"

[bio.physicell]
binary = "../../usr/bin/curl"
config = "config/PhysiCell_settings.xml"
"#,
        )
        .unwrap();
        let workdir = case_dir.join("workdir");
        let case = Case {
            id: "traversal".into(),
            path: case_dir.clone(),
        };
        let result = PhysiCellAdapter::new().prepare(&case, &workdir);
        let err = result.expect_err("parent-dir traversal in binary must be rejected");
        let msg = format!("{err}");
        assert!(
            msg.contains("..") || msg.contains("traverse"),
            "expected traversal rejection message, got: {msg}"
        );
        let _ = std::fs::remove_dir_all(&case_dir);
    }
}
