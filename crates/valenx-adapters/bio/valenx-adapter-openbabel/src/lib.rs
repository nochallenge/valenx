//! # valenx-adapter-openbabel
//!
//! Adapter for [Open Babel](https://openbabel.org/) — the chemical
//! toolbox for converting between ~120 chemistry file formats. The
//! everyday workhorse for "I have a ligand in MOL2, my docking pipeline
//! wants SDF" type problems: SMILES <-> MOL <-> SDF <-> PDB <-> XYZ
//! <-> CIF and many more, plus 3D-coordinate generation, hydrogen
//! addition, and canonicalisation.
//!
//! **Phase 24 — single-binary CLI wrapper.** The user supplies an
//! `input` file plus an `output` filename via `[bio.openbabel]` in
//! `case.toml`. `prepare()` resolves the paths against the case
//! directory and composes the `obabel` invocation; `run()` streams the
//! conversion via the shared subprocess runner.
//!
//! Open Babel sniffs format from extension by default — the explicit
//! `input_format` / `output_format` knobs are escape hatches for
//! stripped filenames or non-standard extensions. `gen_3d = true`
//! lifts 2D structures into 3D; `add_hydrogens = true` adds explicit
//! hydrogens before write.

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

use crate::case_input::OpenBabelInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(OpenBabelAdapter::new())
}

pub struct OpenBabelAdapter;

impl OpenBabelAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for OpenBabelAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "openbabel";
/// Open Babel's CLI binary. `obabel` is the canonical install name on
/// every distro (`apt install openbabel`, `brew install open-babel`,
/// `conda install -c conda-forge openbabel`).
const BINARIES: &[&str] = &["obabel"];

impl Adapter for OpenBabelAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "Open Babel",
            // Open Babel 3.1.x is the long-running stable line every
            // distro ships; 3.0 introduced the modern `obabel` CLI
            // surface and 3.1 polished it. Upper bound 4.0 reserves
            // room for the next major bump.
            version_range: VersionRange {
                min_inclusive: Version::new(3, 1, 0),
                max_exclusive: Version::new(4, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            tool_license: "GPL-2.0",
            docs_url: "https://openbabel.org/docs/dev/",
            homepage_url: "https://openbabel.org/",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(BINARIES) {
            Some(binary_path) => {
                // `obabel -V` prints `Open Babel <version> -- ...`.
                let found_version = detect_tool_version_semver(&binary_path, &["-V", "--version"]);
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
                hint: "Open Babel 3.1+ required; install via `apt install \
                       openbabel`, `brew install open-babel`, or `conda \
                       install -c conda-forge openbabel`"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = OpenBabelInput::from_case_dir(&case.path)?;

        // Round-10 H3: `output` is `PathBuf` and pre-fix flowed into
        // `workdir.join(&input.output)`. Validate as a basename
        // before the join so `output = "../etc/passwd"` is rejected.
        if let Some(s) = input.output.to_str() {
            valenx_core::adapter_helpers::validate_output_basename(s, "[bio.openbabel].output")
                .map_err(|e| AdapterError::InvalidCase {
                    case_path: case.path.join("case.toml"),
                    reason: format!("{e}"),
                })?;
        } else {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: "[bio.openbabel].output: non-UTF-8 path rejected".into(),
            });
        }

        fs::create_dir_all(workdir)?;

        // Resolve the input path against the case directory if relative.
        // The user authors `input = "ligand.mol2"` and expects it to
        // live alongside `case.toml`.
        let source_input = if input.input.is_absolute() {
            input.input.clone()
        } else {
            valenx_core::adapter_helpers::confined_join(&case.path, &input.input)?
        };
        if !source_input.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.openbabel].input `{}` not found (resolved {})",
                    input.input.display(),
                    source_input.display()
                ),
            });
        }

        // The output path lives inside the workdir by convention. We
        // accept absolute paths verbatim (the user owns the path).
        let dest_output = if input.output.is_absolute() {
            input.output.clone()
        } else {
            workdir.join(&input.output)
        };

        let binary_path = find_on_path(BINARIES).ok_or_else(|| AdapterError::ToolNotInstalled {
            name: INFO_ID,
            hint: "Open Babel 3.1+ required; install via `apt install \
                       openbabel`, `brew install open-babel`, or `conda \
                       install -c conda-forge openbabel`"
                .into(),
        })?;

        // Compose the `obabel` invocation:
        //   obabel <input> -O <output> [-i <fmt>] [-o <fmt>]
        //          [--gen3D] [-h] [extras...]
        let mut native_command: Vec<OsString> = vec![
            binary_path.into_os_string(),
            source_input.into_os_string(),
            OsString::from("-O"),
            dest_output.clone().into_os_string(),
        ];
        if let Some(ref fmt) = input.input_format {
            native_command.push(OsString::from("-i"));
            native_command.push(OsString::from(fmt));
        }
        if let Some(ref fmt) = input.output_format {
            native_command.push(OsString::from("-o"));
            native_command.push(OsString::from(fmt));
        }
        if input.gen_3d {
            native_command.push(OsString::from("--gen3D"));
        }
        if input.add_hydrogens {
            native_command.push(OsString::from("-h"));
        }
        for arg in &input.extra_args {
            native_command.push(OsString::from(arg));
        }

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            // Format conversion is fast; even multi-thousand-molecule
            // SDF batches finish in seconds. 5 minutes covers
            // `--gen3D` runs (which call out to a forcefield optimiser)
            // without being absurd.
            estimated_runtime: Some(Duration::from_secs(5 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting Open Babel", |line| {
            let mut hint = subprocess::Hint::default();
            // Open Babel emits "N molecules converted" on stderr at
            // the end. Lift that to a 95% progress tick.
            if line.contains("molecules converted") {
                hint.progress = Some((95.0, line.to_string()));
            } else if line.contains("Error") || line.contains("error") {
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
        // Re-parse the case to learn which output the user asked for
        // (collect runs against a workdir, not a Case, so we hop back
        // via the staged case.toml — the same trick PyMOL / VMD use).
        // If the case.toml isn't accessible from the workdir we fall
        // back to the workdir-walk path used by the script-driven
        // siblings.
        let case_toml = job.workdir.join("case.toml");
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "Open Babel",
            "unknown",
            &case_toml,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);

        // Walk the workdir and surface any file at the top level as the
        // converted output. Open Babel writes exactly the file the user
        // asked for, so a single-file walk suffices — but we list every
        // top-level file so cases with absolute output paths still
        // surface auxiliary artifacts (logs, cached fingerprints) the
        // user dropped in the workdir.
        let entries = match fs::read_dir(&job.workdir) {
            Ok(e) => e,
            Err(e) => {
                tracing::warn!(target: "valenx-openbabel", ?e, "workdir read failed");
                return Ok(results);
            }
        };
        let mut artefacts: Vec<Artifact> = Vec::new();
        for entry in entries.flatten() {
            let path = entry.path();
            if !path.is_file() {
                continue;
            }
            // Skip the staged case.toml — internal infrastructure, not
            // a deliverable.
            if path
                .file_name()
                .and_then(|s| s.to_str())
                .map(|s| s == "case.toml")
                .unwrap_or(false)
            {
                continue;
            }
            artefacts.push(Artifact {
                path,
                kind: ArtifactKind::Native,
                checksum: None,
                label: "Open Babel converted file".to_string(),
            });
        }
        artefacts.sort_by(|a, b| a.path.cmp(&b.path));
        results.artifacts = artefacts;
        Ok(results)
    }

    fn capabilities(&self) -> Capabilities {
        Capabilities {
            capabilities: Vec::new(),
            ribbon_contributions: vec!["bio.openbabel.convert"],
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn info_is_bio_domain() {
        let info = OpenBabelAdapter::new().info();
        assert_eq!(info.id, "openbabel");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "GPL-2.0");
        assert_eq!(info.display_name, "Open Babel");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = OpenBabelAdapter::new().info();
        assert_eq!(info.version_range.min_inclusive, Version::new(3, 1, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(4, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = OpenBabelAdapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.openbabel.convert"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = OpenBabelAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }

    /// Round-10 H3 RED→GREEN: `output` flowed into
    /// `workdir.join(...)` with no validation. Hostile
    /// `output = "../etc/passwd"` is now rejected.
    #[test]
    fn prepare_rejects_output_path_traversal() {
        use valenx_test_utils::tempdir;
        let d = tempdir("openbabel-output-trav");
        std::fs::write(d.join("lig.mol2"), b"# fake\n").unwrap();
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "openbabel.convert"

[bio.openbabel]
input  = "lig.mol2"
output = "../etc/passwd"
"#,
        )
        .unwrap();
        let case = Case {
            id: "trav".into(),
            path: d.clone(),
        };
        let workdir = d.join("workdir");
        let err = OpenBabelAdapter::new()
            .prepare(&case, &workdir)
            .unwrap_err();
        let msg = format!("{err}");
        assert!(
            msg.contains("[bio.openbabel].output"),
            "expected [bio.openbabel].output in error, got: {msg}"
        );
        let _ = std::fs::remove_dir_all(&d);
    }
}
