//! # valenx-adapter-psi4
//!
//! Adapter for [Psi4](https://psicode.org/) — an open-source ab
//! initio / DFT / post-Hartree-Fock quantum-chemistry package. Psi4
//! drives Hartree-Fock, density-functional theory, MP2/CCSD/CCSD(T),
//! multireference methods, and a stack of analysis tools, all through
//! a Python-driven "Psithon" input file. It is the de-facto open
//! quantum-chemistry workhorse in academic computational chemistry
//! pipelines.
//!
//! **Phase 25 — subprocess wrapper around the `psi4` binary.** The
//! user supplies a Psithon input file and an output path via
//! `[bio.psi4]` in `case.toml`. `prepare()` resolves the input,
//! composes a `psi4 -i <input> -o <output> -n <threads> [-m <memory>]
//! [extras...]` invocation, and stages everything in the workdir.
//! `run()` streams the run via the shared subprocess runner — Psi4
//! does most of its chatter on stdout, line by line, so the standard
//! handler picks up "==> Iterations" / "@DF-RHF iter" / "Job ended"
//! markers for free.
//!
//! On `collect()` we surface the canonical text output as the run's
//! `Log` artifact, and walk the workdir for the conventional binary
//! companions Psi4 emits when the user asks for them — `*.fchk`
//! (formatted checkpoint files for orbital exchange) and `*.molden`
//! (MO data for visualisation in Avogadro / VMD / Jmol).

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

use crate::case_input::Psi4Input;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(Psi4Adapter::new())
}

pub struct Psi4Adapter;

impl Psi4Adapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for Psi4Adapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "psi4";
/// Psi4's binary candidate. Conda / source / Bioconda / Homebrew all
/// install under the canonical `psi4` name.
const BINARIES: &[&str] = &["psi4"];

/// The default memory string the adapter pads onto the CLI. When the
/// case opts into the default we omit `-m` entirely so Psi4 uses its
/// own internal default (currently "500 mb"); when the case asks for
/// something different we pass `-m <value>` through.
const DEFAULT_MEMORY: &str = "1 gb";

impl Adapter for Psi4Adapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "Psi4",
            // Psi4's 1.x line has been the stable series for years;
            // 1.8 is the floor we test against (current: 1.9+).
            // Upper bound 2.0 reserves room for a major bump.
            version_range: VersionRange {
                min_inclusive: Version::new(1, 8, 0),
                max_exclusive: Version::new(2, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            tool_license: "LGPL-3.0",
            docs_url: "https://psicode.org/psi4manual/master/",
            homepage_url: "https://psicode.org/",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(BINARIES) {
            Some(binary_path) => {
                // `psi4 --version` prints something like "1.9.1" on
                // stdout; the helper handles the detection.
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
                hint: "Psi4 1.8+ required; install via \
                       `conda install -c psi4 psi4` or build from source"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = Psi4Input::from_case_dir(&case.path)?;

        // Round-10 H3: `output` is `PathBuf` and pre-fix flowed into
        // `workdir.join(&input.output)`. Validate as a basename
        // before the join so `output = "../etc/passwd"` is rejected.
        if let Some(s) = input.output.to_str() {
            valenx_core::adapter_helpers::validate_output_basename(s, "[bio.psi4].output")
                .map_err(|e| AdapterError::InvalidCase {
                    case_path: case.path.join("case.toml"),
                    reason: format!("{e}"),
                })?;
        } else {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: "[bio.psi4].output: non-UTF-8 path rejected".into(),
            });
        }

        fs::create_dir_all(workdir)?;

        // Resolve the input Psithon file against the case directory if
        // relative. Same convention as every other Phase 17/18 bio
        // adapter — `input = "input.dat"` next to `case.toml`.
        let source_input = if input.input.is_absolute() {
            input.input.clone()
        } else {
            valenx_core::adapter_helpers::confined_join(&case.path, &input.input)?
        };
        if !source_input.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.psi4].input `{}` not found (resolved {})",
                    input.input.display(),
                    source_input.display()
                ),
            });
        }

        // The output path is *generated*, not consumed. Resolve
        // relative paths into the workdir so the artifact lands next
        // to whatever else Psi4 writes (scratch files, fchk, molden).
        let output_path = if input.output.is_absolute() {
            input.output.clone()
        } else {
            workdir.join(&input.output)
        };

        let binary_path = find_on_path(BINARIES).ok_or_else(|| AdapterError::ToolNotInstalled {
            name: INFO_ID,
            hint: "Psi4 1.8+ required; install via \
                       `conda install -c psi4 psi4` or build from source"
                .into(),
        })?;

        // Compose `psi4 -i <input> -o <output> -n <threads>
        //              [-m <memory>] [extras...]`.
        // `-m` is only emitted when the user asked for something
        // other than the documented default — Psi4's own default is
        // "500 mb"; passing "1 gb" (our adapter default) every time
        // would override Psi4's internal default with a fixed value
        // even when the user didn't ask for one. Skipping `-m` when
        // the case left the field at the default keeps Psi4's own
        // behaviour intact.
        let mut native_command: Vec<OsString> = vec![
            binary_path.into_os_string(),
            OsString::from("-i"),
            source_input.into_os_string(),
            OsString::from("-o"),
            output_path.clone().into_os_string(),
            OsString::from("-n"),
            OsString::from(input.threads.to_string()),
        ];
        if input.memory != DEFAULT_MEMORY {
            native_command.push(OsString::from("-m"));
            native_command.push(OsString::from(&input.memory));
        }
        for arg in &input.extra_args {
            native_command.push(OsString::from(arg));
        }

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            // Single-point energies finish in seconds; CCSD(T) on
            // mid-sized basis sets runs for hours, and full
            // multireference jobs can run for days. 12 hours is a
            // generous default that covers the typical long tail.
            estimated_runtime: Some(Duration::from_secs(12 * 60 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting Psi4", |line| {
            let mut hint = subprocess::Hint::default();
            // Psi4's progress markers on stdout: a "Memory set to" /
            // "Threads set to" banner at startup, "==> Iterations"
            // when SCF kicks in, "@DF-RHF iter" lines for each SCF
            // cycle, and a "*** Psi4 exiting successfully" sentinel
            // at the end. Lift the obvious milestones into UI ticks.
            if line.contains("Psi4 exiting successfully") {
                hint.progress = Some((95.0, line.to_string()));
            } else if line.contains("==> Iterations") || line.contains("Final Energy:") {
                hint.progress = Some((70.0, line.to_string()));
            } else if line.contains("==> Geometry") || line.contains("Threads set to") {
                hint.progress = Some((10.0, line.to_string()));
            }
            if line.contains("Error:") || line.contains("RuntimeError") {
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
        // Provenance: hash the staged case.toml as the canonical input
        // descriptor. Psi4's actual output isn't fixed-name (the user
        // chooses it via `output = "..."`) so we walk the workdir for
        // it on the artifact side rather than try to pre-compute its
        // path here.
        let case_hash_input = job.workdir.join("case.toml");
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "Psi4",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);
        let mut artefacts: Vec<Artifact> = Vec::new();

        // Walk the workdir top-level. Psi4 writes its output to the
        // path passed via `-o`; auxiliary outputs (`.fchk` formatted
        // checkpoints, `.molden` orbital data) appear when the
        // Psithon script asks for them.
        let entries = match fs::read_dir(&job.workdir) {
            Ok(e) => e,
            Err(e) => {
                tracing::warn!(target: "valenx-psi4", ?e, "workdir read failed");
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
                // The text output we asked Psi4 to write via `-o`.
                // We assume the conventional `.dat` / `.out` suffix;
                // either way we mark it as Log so the UI knows it's
                // textual rather than binary.
                Some("dat") | Some("out") => (ArtifactKind::Log, "Psi4 output".to_string()),
                Some("log") => (ArtifactKind::Log, "Psi4 log".to_string()),
                // Formatted checkpoint files — orbital data in
                // Gaussian's text format. Native (binary-equivalent
                // text format consumed by Multiwfn / cubegen).
                Some("fchk") => (
                    ArtifactKind::Native,
                    "Psi4 formatted checkpoint".to_string(),
                ),
                // Molden orbital files for visualisation in Avogadro
                // / VMD / Jmol.
                Some("molden") => (ArtifactKind::Native, "Psi4 Molden orbitals".to_string()),
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
            ribbon_contributions: vec!["bio.psi4.compute"],
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn info_is_bio_domain() {
        let info = Psi4Adapter::new().info();
        assert_eq!(info.id, "psi4");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "LGPL-3.0");
        assert_eq!(info.display_name, "Psi4");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = Psi4Adapter::new().info();
        assert_eq!(info.version_range.min_inclusive, Version::new(1, 8, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(2, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = Psi4Adapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.psi4.compute"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = Psi4Adapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }

    /// Round-10 H3 RED→GREEN: `output` flowed into
    /// `workdir.join(...)` with no validation. Hostile
    /// `output = "../etc/passwd"` is now rejected.
    #[test]
    fn prepare_rejects_output_path_traversal() {
        use valenx_test_utils::tempdir;
        let d = tempdir("psi4-output-trav");
        std::fs::write(d.join("inp.dat"), b"# test\n").unwrap();
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "psi4.compute"

[bio.psi4]
input   = "inp.dat"
output  = "../etc/passwd"
threads = 1
"#,
        )
        .unwrap();
        let case = Case {
            id: "trav".into(),
            path: d.clone(),
        };
        let workdir = d.join("workdir");
        let err = Psi4Adapter::new().prepare(&case, &workdir).unwrap_err();
        let msg = format!("{err}");
        assert!(
            msg.contains("[bio.psi4].output"),
            "expected [bio.psi4].output in error, got: {msg}"
        );
        let _ = std::fs::remove_dir_all(&d);
    }
}
