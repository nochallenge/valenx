//! # valenx-adapter-rosetta
//!
//! Adapter for [Rosetta](https://www.rosettacommons.org/) — the
//! flagship modeling suite from RosettaCommons. Rosetta drives
//! protein design, structure prediction, docking, ligand binding,
//! and a long tail of related modeling tasks through its
//! `rosetta_scripts` binary, which reads an XML protocol describing
//! the modeling pipeline (filters, movers, scorefunctions) and
//! applies it to an input `.pdb`.
//!
//! **Phase 38 — subprocess wrapper around `rosetta_scripts`.** The
//! user supplies an XML protocol, an input `.pdb`, an output
//! basename, the number of decoys (`nstruct`), and the path to the
//! Rosetta `database/` data directory via `[bio.rosetta]` in
//! `case.toml`. `prepare()` composes a
//! `rosetta_scripts -database <path> -parser:protocol <xml>
//! -in:file:s <pdb> -out:prefix <basename> -nstruct <N> [extras...]`
//! invocation and stages everything in the workdir. `run()` streams
//! the run via the shared subprocess runner; `rosetta_scripts`
//! prints progress chatter to stdout (`apply` / `Finished`
//! markers) line-by-line.
//!
//! On `collect()` we walk the workdir for `<output_basename>*.pdb`
//! decoys and the canonical `score.sc` scorefile.
//!
//! ## License flag
//!
//! Rosetta ships under the RosettaCommons license — non-commercial
//! / academic use only without a separate commercial agreement. We
//! surface this accurately via `tool_license = "Rosetta-License"`
//! and emit a probe warning whenever the binary is found, with the
//! literal string `"academic"` as a stable anchor for tests and
//! downstream license-aware filters.

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

use crate::case_input::RosettaInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(RosettaAdapter::new())
}

pub struct RosettaAdapter;

impl RosettaAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for RosettaAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "rosetta";
/// Rosetta's binary candidates. Rosetta source builds emit
/// platform-suffixed names (`rosetta_scripts.linuxgccrelease`,
/// `.macosclangrelease`) by default; conda / packaged distributions
/// install a bare `rosetta_scripts` shim. Probe all three so the
/// most common install layouts are covered.
const BINARIES: &[&str] = &[
    "rosetta_scripts",
    "rosetta_scripts.linuxgccrelease",
    "rosetta_scripts.macosclangrelease",
];

/// The probe-warning surfaced whenever Rosetta is detected. Anchors
/// a stable "academic / non-commercial only" reminder for downstream
/// tooling and tests; the literal string `"academic"` is part of
/// the asserted contract.
const LICENSE_WARNING: &str = "Rosetta is licensed for non-commercial / academic use only. \
     Confirm your use case complies with the RosettaCommons license \
     before redistributing designs or derived data.";

impl Adapter for RosettaAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "Rosetta",
            // Rosetta's stable line is the 3.x series; 3.13 (2021)
            // is the floor we test against and covers every recent
            // release through 3.14. Upper bound 4.0 reserves room
            // for an eventual major bump.
            version_range: VersionRange {
                min_inclusive: Version::new(3, 13, 0),
                max_exclusive: Version::new(4, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            // Rosetta's terms aren't a recognised SPDX identifier;
            // the closest accurate label is the project's own
            // "Rosetta-License" name. Surfacing it here keeps
            // license-aware tooling honest.
            tool_license: "Rosetta-License",
            docs_url: "https://www.rosettacommons.org/docs/latest/",
            homepage_url: "https://www.rosettacommons.org/",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(BINARIES) {
            Some(binary_path) => {
                // `rosetta_scripts -version` (one dash) prints a
                // banner with the Rosetta release on stdout. The
                // generic detector tries common version flags.
                let found_version =
                    detect_tool_version_semver(&binary_path, &["-version", "--version"]);
                Ok(ProbeReport {
                    ok: true,
                    found_version,
                    binary_path: Some(binary_path),
                    // Always surface the license reminder when
                    // Rosetta is detected — it's a custom non-OSS
                    // license and we'd rather over-warn than have a
                    // user ship commercial output without checking.
                    warnings: vec![LICENSE_WARNING.to_string()],
                    required_env: Vec::new(),
                })
            }
            None => Err(AdapterError::ToolNotInstalled {
                name: INFO_ID,
                hint: "Rosetta 3.13+ required; obtain from \
                       https://www.rosettacommons.org/software/license-and-download \
                       (registration required, academic-use license) and \
                       ensure `rosetta_scripts` is on PATH"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = RosettaInput::from_case_dir(&case.path)?;

        // Round-4 security: reject `output_basename = "../etc/passwd"`
        // and friends before the value flows into any path join.
        // Same pattern as the round-3 fix in bionetgen/iqtree/art/fasttree.
        valenx_core::adapter_helpers::validate_output_basename(
            &input.output_basename,
            "[bio.rosetta].output_basename",
        )
        .map_err(|e| AdapterError::InvalidCase {
            case_path: case.path.join("case.toml"),
            reason: format!("{e}"),
        })?;

        fs::create_dir_all(workdir)?;

        // Resolve the protocol XML against the case directory if
        // relative. Same convention as every other Phase 17/18 bio
        // adapter — `protocol = "design.xml"` next to `case.toml`.
        let source_protocol = if input.protocol.is_absolute() {
            input.protocol.clone()
        } else {
            valenx_core::adapter_helpers::confined_join(&case.path, &input.protocol)?
        };
        if !source_protocol.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.rosetta].protocol `{}` not found (resolved {})",
                    input.protocol.display(),
                    source_protocol.display()
                ),
            });
        }

        // Resolve the input PDB against the case directory if
        // relative.
        let source_pdb = if input.input_pdb.is_absolute() {
            input.input_pdb.clone()
        } else {
            valenx_core::adapter_helpers::confined_join(&case.path, &input.input_pdb)?
        };
        if !source_pdb.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.rosetta].input_pdb `{}` not found (resolved {})",
                    input.input_pdb.display(),
                    source_pdb.display()
                ),
            });
        }

        // The database is a directory the user surfaces verbatim;
        // it's not staged into the workdir (Rosetta data dirs run
        // to gigabytes and are read-only at runtime). Round-8
        // sibling-field sweep: any user-supplied relative path is
        // still funneled through `confined_join` to refuse `..`
        // traversal out of the case sandbox. Absolute paths (the
        // common case for shared system installs) continue to be
        // forwarded verbatim — `confined_join` returns an error for
        // absolutes, so the absolute branch stays here.
        let database = if input.database.is_absolute() {
            input.database.clone()
        } else {
            valenx_core::adapter_helpers::confined_join(&case.path, &input.database)?
        };

        let binary_path = find_on_path(BINARIES).ok_or_else(|| AdapterError::ToolNotInstalled {
            name: INFO_ID,
            hint: "Rosetta 3.13+ required; obtain from \
                       https://www.rosettacommons.org/software/license-and-download \
                       (registration required, academic-use license) and \
                       ensure `rosetta_scripts` is on PATH"
                .into(),
        })?;

        // Compose `rosetta_scripts -database <path>
        //                          -parser:protocol <xml>
        //                          -in:file:s <pdb>
        //                          -out:prefix <basename>
        //                          -nstruct <N>
        //                          [extras...]`.
        let mut native_command: Vec<OsString> = vec![
            binary_path.into_os_string(),
            OsString::from("-database"),
            database.into_os_string(),
            OsString::from("-parser:protocol"),
            source_protocol.into_os_string(),
            OsString::from("-in:file:s"),
            source_pdb.into_os_string(),
            OsString::from("-out:prefix"),
            OsString::from(&input.output_basename),
            OsString::from("-nstruct"),
            OsString::from(input.nstruct.to_string()),
        ];
        for arg in &input.extra_args {
            native_command.push(OsString::from(arg));
        }

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            // Single-decoy designs finish in seconds-to-minutes;
            // a 1000-decoy FastDesign run takes hours. 12 hours is
            // a generous default that covers the typical long tail.
            estimated_runtime: Some(Duration::from_secs(12 * 60 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting Rosetta", |line| {
            let mut hint = subprocess::Hint::default();
            // Rosetta's progress markers on stdout: `protocols.jd2`
            // banners at startup, `apply` lines when a mover runs,
            // `Finished` / `successfully completed` at end-of-run.
            if line.contains("successfully completed") || line.contains("Finished") {
                hint.progress = Some((95.0, line.to_string()));
            } else if line.contains("apply()") || line.contains("protocols.rosetta_scripts") {
                hint.progress = Some((50.0, line.to_string()));
            } else if line.contains("protocols.jd2") {
                hint.progress = Some((10.0, line.to_string()));
            }
            if line.contains("ERROR") || line.contains("FATAL") {
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
        // Provenance: hash the canonical scorefile if one was
        // produced, else case.toml so the provenance block stays
        // well-formed for partial / failed runs.
        let case_hash_input = {
            let score = job.workdir.join("score.sc");
            if score.is_file() {
                score
            } else {
                job.workdir.join("case.toml")
            }
        };
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "Rosetta",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);
        let mut artefacts: Vec<Artifact> = Vec::new();

        // Read the staged case.toml back out so we can restrict the
        // collected `.pdb` outputs to those whose stem starts with
        // the configured `output_basename`. Failure to read the
        // case is non-fatal — we then accept every `.pdb` as a
        // potential decoy.
        let basename = read_output_basename(&job.workdir);

        let entries = match fs::read_dir(&job.workdir) {
            Ok(e) => e,
            Err(e) => {
                tracing::warn!(target: "valenx-rosetta", ?e, "workdir read failed");
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
            // The canonical scorefile — always pick it up.
            if name.as_deref() == Some("score.sc") {
                artefacts.push(Artifact {
                    path,
                    kind: ArtifactKind::Tabular,
                    checksum: None,
                    label: "Rosetta scores".to_string(),
                });
                continue;
            }
            let ext = path
                .extension()
                .and_then(|s| s.to_str())
                .map(|s| s.to_ascii_lowercase());
            if ext.as_deref() == Some("pdb") {
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
                    kind: ArtifactKind::Native,
                    checksum: None,
                    label: "Rosetta designed structure".to_string(),
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
            ribbon_contributions: vec!["bio.rosetta.protocol"],
        }
    }
}

/// Pull `output_basename` out of the staged `case.toml` for
/// `collect()`-time decoy filtering. Returns None if the file
/// doesn't exist or can't be parsed — collect falls back to
/// accepting every `.pdb` in that case.
fn read_output_basename(workdir: &Path) -> Option<String> {
    let case_toml = workdir.join("case.toml");
    let text = valenx_core::io_caps::read_capped_to_string(
        &case_toml,
        valenx_core::project::loader::MAX_PROJECT_FILE_BYTES as usize,
    )
    .ok()?;
    let parsed: toml::Value = toml::from_str(&text).ok()?;
    parsed
        .get("bio")?
        .get("rosetta")?
        .get("output_basename")?
        .as_str()
        .map(|s| s.to_string())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn info_is_bio_domain() {
        let info = RosettaAdapter::new().info();
        assert_eq!(info.id, "rosetta");
        assert_eq!(info.physics, &[Physics::Bio]);
        // Rosetta's custom non-OSS license, not a recognised SPDX
        // identifier — pin the project's own label.
        assert_eq!(info.tool_license, "Rosetta-License");
        assert_eq!(info.display_name, "Rosetta");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = RosettaAdapter::new().info();
        // Rosetta 3.13 (2021) is the floor; upper bound 4.0
        // reserves room for the next major.
        assert_eq!(info.version_range.min_inclusive, Version::new(3, 13, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(4, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = RosettaAdapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.rosetta.protocol"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = RosettaAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }

    #[test]
    fn probe_warning_mentions_academic() {
        // The license-flag warning is mandatory: Rosetta is
        // non-OSS academic-use, and we surface that on every
        // successful probe. The literal "academic" anchor is what
        // downstream tooling and license-aware filters key off —
        // pin it.
        assert!(
            LICENSE_WARNING.contains("academic"),
            "probe warning must contain `academic` anchor; got: {LICENSE_WARNING}"
        );
    }

    #[test]
    fn prepare_rejects_relative_traversal_database_path() {
        // Round-8 RED→GREEN: a relative `database` entry now routes
        // through `confined_join` (absolute paths to system installs
        // still pass through verbatim — confined_join refuses
        // absolutes but the absolute branch explicitly forwards them).
        // This test confirms `database = "../etc/passwd"` (relative
        // traversal escape) is rejected.
        use valenx_test_utils::tempdir;
        let d = tempdir("rosetta-traversal");
        std::fs::write(d.join("protocol.xml"), b"<ROSETTASCRIPTS/>").unwrap();
        std::fs::write(d.join("input.pdb"), b"HEADER").unwrap();
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "rosetta.protocol"

[bio.rosetta]
protocol        = "protocol.xml"
input_pdb       = "input.pdb"
database        = "../etc/passwd"
output_basename = "out"
nstruct         = 1
"#,
        )
        .unwrap();
        let case = Case {
            id: "rosetta-traversal".into(),
            path: d.clone(),
        };
        let workdir = d.join("workdir");
        let err = RosettaAdapter::new().prepare(&case, &workdir).unwrap_err();
        let msg = format!("{err}");
        assert!(
            msg.contains("absolute")
                || msg.contains("escape")
                || msg.contains("traversal")
                || msg.contains(".."),
            "expected confined_join rejection on database traversal, got: {msg}"
        );
        let _ = std::fs::remove_dir_all(&d);
    }
}
