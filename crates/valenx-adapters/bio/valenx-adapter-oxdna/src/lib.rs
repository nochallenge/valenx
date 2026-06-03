//! # valenx-adapter-oxdna
//!
//! Adapter for [oxDNA](https://lorenzo-rovigatti.github.io/oxDNA/) —
//! the coarse-grained molecular-dynamics engine for DNA / RNA. Where
//! the all-atom MD adapters (GROMACS / LAMMPS / OpenMM) treat every
//! heavy atom as a node, oxDNA collapses each nucleotide into a
//! handful of beads and runs orders of magnitude faster, making it
//! the de-facto choice for DNA-origami self-assembly studies and
//! large-scale aptamer thermodynamics.
//!
//! **Phase 17 — subprocess wrapper for user-provided control files.**
//! oxDNA reads everything (initial conformation, topology,
//! integration parameters, force-field selection) from a single
//! `input.dat` control file. The user references that path via
//! `[bio.oxdna].input` in `case.toml`; `prepare()` stages the file
//! (and an optional explicit topology) into the workdir and `run()`
//! invokes `oxDNA input.dat` via the shared subprocess runner.
//!
//! On `collect()` we list the three customary outputs as `Native`
//! artifacts: `last_conf.dat` (final configuration),
//! `trajectory.dat` (per-step trajectory dump), and `energy.dat`
//! (per-step thermodynamics). Structured trajectory parsing — the
//! step where `energy.dat` graduates into typed `ScalarRecord`
//! entries with a `TimeKey::Iteration` axis — lands in a follow-up
//! phase alongside the MDAnalysis adapter.

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

use crate::case_input::OxDnaInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(OxDnaAdapter::new())
}

pub struct OxDnaAdapter;

impl OxDnaAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for OxDnaAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "oxdna";
/// oxDNA's executable. The repo distributes the binary under the
/// CamelCase `oxDNA` name on every platform; lowercase isn't a
/// canonical install path so we keep the candidate list narrow.
const BINARIES: &[&str] = &["oxDNA"];

impl Adapter for OxDnaAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "oxDNA",
            // oxDNA 3.5 is the first release with the consolidated
            // CMake build + CUDA backend we lean on; the upper bound
            // bumps when a 4.x line lands.
            version_range: VersionRange {
                min_inclusive: Version::new(3, 5, 0),
                max_exclusive: Version::new(4, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            tool_license: "GPL-3.0-or-later",
            docs_url: "https://lorenzo-rovigatti.github.io/oxDNA/",
            homepage_url: "https://github.com/lorenzo-rovigatti/oxDNA",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(BINARIES) {
            Some(binary_path) => {
                // oxDNA prints a version banner in its startup
                // output; `--version` flags vary by build but the
                // detector tries both forms.
                let found_version = detect_tool_version_semver(&binary_path, &["--version", "-v"]);
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
                hint: "oxDNA 3.5+ required; build from \
                       https://github.com/lorenzo-rovigatti/oxDNA"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = OxDnaInput::from_case_dir(&case.path)?;

        fs::create_dir_all(workdir)?;

        // Stage the input.dat. Resolved against the case directory —
        // same convention as every other Phase 17 bio adapter.
        // `confined_join` rejects absolute paths and `..` traversal so
        // a malicious case bundle can't smuggle arbitrary host files
        // into the workdir.
        let source_input = confined_join(&case.path, &input.input)?;
        if !source_input.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.oxdna].input `{}` not found (resolved {})",
                    input.input.display(),
                    source_input.display()
                ),
            });
        }
        let input_filename = input
            .input
            .file_name()
            .ok_or_else(|| AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.oxdna].input path `{}` has no filename",
                    input.input.display()
                ),
            })?;
        let dest_input = workdir.join(input_filename);
        if source_input != dest_input {
            fs::copy(&source_input, &dest_input)?;
        }

        // Optionally stage the .top file. oxDNA reads the topology
        // path from the `topology = ...` line inside input.dat; if
        // the user passes it explicitly via `[bio.oxdna].topology`
        // we honour that and copy it next to the staged input so
        // the relative path resolves cleanly.
        if let Some(topology) = &input.topology {
            let source_top = confined_join(&case.path, topology)?;
            if !source_top.is_file() {
                return Err(AdapterError::InvalidCase {
                    case_path: case.path.join("case.toml"),
                    reason: format!(
                        "[bio.oxdna].topology `{}` not found (resolved {})",
                        topology.display(),
                        source_top.display()
                    ),
                });
            }
            let top_filename = topology
                .file_name()
                .ok_or_else(|| AdapterError::InvalidCase {
                    case_path: case.path.join("case.toml"),
                    reason: format!(
                        "[bio.oxdna].topology path `{}` has no filename",
                        topology.display()
                    ),
                })?;
            let dest_top = workdir.join(top_filename);
            if source_top != dest_top {
                fs::copy(&source_top, &dest_top)?;
            }
        }

        let binary_path = find_on_path(BINARIES).ok_or_else(|| AdapterError::ToolNotInstalled {
            name: INFO_ID,
            hint: "oxDNA 3.5+ required; build from \
                       https://github.com/lorenzo-rovigatti/oxDNA"
                .into(),
        })?;

        let native_command: Vec<OsString> =
            vec![binary_path.into_os_string(), OsString::from(input_filename)];

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            // Coarse-grained MD is fast — small duplex equilibration
            // finishes in seconds, large origami runs span hours
            // (still much faster than all-atom). 4 hours covers the
            // long tail; longer runs override via their own progress
            // reporting.
            estimated_runtime: Some(Duration::from_secs(4 * 60 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting oxDNA", |line| {
            let mut hint = subprocess::Hint::default();
            // oxDNA writes periodic step markers as it integrates
            // — lines starting with "Step " or containing "MC step"
            // (depending on integrator) are the canonical signal
            // that the simulation is progressing. We can't compute
            // a real percentage without the total step count, so
            // pin the spinner at 50% mid-run.
            if line.contains("MC step") || line.starts_with("Step ") {
                hint.progress = Some((50.0, line.trim().to_string()));
            } else if line.contains("ERROR") || line.contains("FATAL") {
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
        // Provenance: hash the staged input.dat (the canonical
        // "this case is configured this way" input). We don't know
        // the user's mesh / lock files so leave those empty.
        let input_path = first_input_in_workdir(&job.workdir);
        let case_hash_input = input_path
            .clone()
            .unwrap_or_else(|| job.workdir.join("case.toml"));
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "oxDNA",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);

        // Walk the workdir top level and pick up oxDNA's three
        // customary outputs (last_conf.dat, trajectory.dat,
        // energy.dat) plus any .top topology file we staged. Other
        // .dat files are surfaced as Native too — oxDNA users
        // sometimes configure custom output filenames via input.dat.
        let entries = match fs::read_dir(&job.workdir) {
            Ok(e) => e,
            Err(e) => {
                tracing::warn!(target: "valenx-oxdna", ?e, "workdir read failed");
                return Ok(results);
            }
        };
        let mut artefacts: Vec<Artifact> = Vec::new();
        for entry in entries.flatten() {
            let path = entry.path();
            if !path.is_file() {
                continue;
            }
            let name = path
                .file_name()
                .and_then(|s| s.to_str())
                .unwrap_or("")
                .to_string();
            let ext = path
                .extension()
                .and_then(|s| s.to_str())
                .map(|s| s.to_ascii_lowercase());
            let (kind, label) = match (name.as_str(), ext.as_deref()) {
                // Three canonical oxDNA outputs — recognised by
                // exact filename so we can label them precisely.
                ("last_conf.dat", _) => (
                    ArtifactKind::Native,
                    "oxDNA final configuration".to_string(),
                ),
                ("trajectory.dat", _) => (ArtifactKind::Native, "oxDNA trajectory".to_string()),
                ("energy.dat", _) => (ArtifactKind::Native, "oxDNA energy log".to_string()),
                // input.dat — the staged control file. Surface so
                // the user can re-inspect parameters from the
                // results pane.
                ("input.dat", _) => (ArtifactKind::Other, "oxDNA input control file".to_string()),
                // .top — the topology file we may have staged
                // alongside input.dat.
                (_, Some("top")) => (ArtifactKind::Native, "oxDNA topology".to_string()),
                // Any other .dat — oxDNA users sometimes rename
                // the output files in input.dat. List as Native
                // with a generic label so they don't disappear.
                (_, Some("dat")) => (ArtifactKind::Native, "oxDNA output (.dat)".to_string()),
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
        // The bio-specific Capability variants land in a follow-up
        // task; for now we publish an empty capability vector and
        // a single ribbon contribution so the registry can wire the
        // adapter in without crashing the UI's capability-index
        // builder.
        Capabilities {
            capabilities: Vec::new(),
            ribbon_contributions: vec!["bio.oxdna.batch"],
        }
    }
}

/// Lift the staged input.dat out of a workdir for provenance
/// hashing. Returns the file named `input.dat` at the top level
/// when present; falls back to the lexicographically-first `.dat`
/// file otherwise.
fn first_input_in_workdir(workdir: &Path) -> Option<PathBuf> {
    let canonical = workdir.join("input.dat");
    if canonical.is_file() {
        return Some(canonical);
    }
    let entries = fs::read_dir(workdir).ok()?;
    let mut hits: Vec<PathBuf> = entries
        .flatten()
        .map(|e| e.path())
        .filter(|p| {
            p.extension()
                .and_then(|s| s.to_str())
                .map(|s| s.eq_ignore_ascii_case("dat"))
                .unwrap_or(false)
        })
        .collect();
    hits.sort();
    hits.into_iter().next()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn info_is_bio_domain() {
        let info = OxDnaAdapter::new().info();
        assert_eq!(info.id, "oxdna");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "GPL-3.0-or-later");
        assert_eq!(info.display_name, "oxDNA");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = OxDnaAdapter::new().info();
        // We support oxDNA >= 3.5 (consolidated CMake + CUDA);
        // expect to revisit upper bound when 4.x lands.
        assert_eq!(info.version_range.min_inclusive, Version::new(3, 5, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(4, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = OxDnaAdapter::new().capabilities();
        // Capability variants land in a future task; ribbon
        // contributions are already enough for the registry to
        // surface the adapter.
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.oxdna.batch"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = OxDnaAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }

    #[test]
    fn first_input_in_workdir_prefers_canonical_name() {
        let tmp = std::env::temp_dir().join(format!(
            "valenx-oxdna-input-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(&tmp).unwrap();
        // Both files exist; the canonical name must win.
        fs::write(tmp.join("a_alphabetically_first.dat"), b"junk").unwrap();
        fs::write(tmp.join("input.dat"), b"# control file").unwrap();
        let f = first_input_in_workdir(&tmp).expect("found");
        assert!(f.ends_with("input.dat"));
        let _ = fs::remove_dir_all(&tmp);
    }

    #[test]
    fn first_input_in_workdir_falls_back_to_dat() {
        let tmp = std::env::temp_dir().join(format!(
            "valenx-oxdna-fallback-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(&tmp).unwrap();
        // No input.dat — fall back to the first .dat file.
        fs::write(tmp.join("custom_run.dat"), b"# user-named control").unwrap();
        fs::write(tmp.join("notes.txt"), b"placeholder").unwrap();
        let f = first_input_in_workdir(&tmp).expect("found");
        assert!(f.ends_with("custom_run.dat"));
        let _ = fs::remove_dir_all(&tmp);
    }

    /// `collect()` must classify the three canonical outputs by
    /// their exact filenames; topology and other .dat files surface
    /// with appropriate generic labels.
    #[test]
    fn collect_classifies_oxdna_outputs() {
        let tmp = std::env::temp_dir().join(format!(
            "valenx-oxdna-collect-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(&tmp).unwrap();
        fs::write(tmp.join("input.dat"), b"# control file").unwrap();
        fs::write(tmp.join("duplex.top"), b"# topology").unwrap();
        fs::write(tmp.join("last_conf.dat"), b"# final conf").unwrap();
        fs::write(tmp.join("trajectory.dat"), b"# traj frames").unwrap();
        fs::write(tmp.join("energy.dat"), b"# step energy").unwrap();
        fs::write(tmp.join("custom_dump.dat"), b"# user-renamed output").unwrap();
        fs::write(tmp.join("ignore.bin"), b"...").unwrap();

        let job = PreparedJob {
            workdir: tmp.clone(),
            native_command: vec![],
            environment: Vec::new(),
            estimated_runtime: None,
            kill_on_drop: true,
        };
        let results = OxDnaAdapter::new().collect(&job).unwrap();
        let labels: Vec<&str> = results.artifacts.iter().map(|a| a.label.as_str()).collect();
        assert!(labels.contains(&"oxDNA final configuration"));
        assert!(labels.contains(&"oxDNA trajectory"));
        assert!(labels.contains(&"oxDNA energy log"));
        assert!(labels.contains(&"oxDNA topology"));
        assert!(labels.contains(&"oxDNA input control file"));
        assert!(labels.contains(&"oxDNA output (.dat)"));
        // .bin must not surface — guards the deny-by-default path.
        assert!(!results
            .artifacts
            .iter()
            .any(|a| a.path.extension().is_some_and(|e| e == "bin")));
        let _ = fs::remove_dir_all(&tmp);
    }

    /// All three canonical outputs classify as `Native` (binary-ish
    /// solver-format output); the topology classifies as `Native`
    /// too. Pin the contract.
    #[test]
    fn collect_artifact_kinds_pin_native() {
        let tmp = std::env::temp_dir().join(format!(
            "valenx-oxdna-kinds-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(&tmp).unwrap();
        fs::write(tmp.join("last_conf.dat"), b"").unwrap();
        fs::write(tmp.join("trajectory.dat"), b"").unwrap();
        fs::write(tmp.join("energy.dat"), b"").unwrap();
        fs::write(tmp.join("system.top"), b"").unwrap();

        let job = PreparedJob {
            workdir: tmp.clone(),
            native_command: vec![],
            environment: Vec::new(),
            estimated_runtime: None,
            kill_on_drop: true,
        };
        let results = OxDnaAdapter::new().collect(&job).unwrap();
        for art in &results.artifacts {
            assert_eq!(art.kind, ArtifactKind::Native, "label was {}", art.label);
        }
        let _ = fs::remove_dir_all(&tmp);
    }
}
