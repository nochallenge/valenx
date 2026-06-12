//! # valenx-adapter-precice
//!
//! Meta-adapter for preCICE — partitioned multi-physics coupling.
//!
//! **Phase 9 — staging + concurrent orchestration scaffold.** The
//! adapter's `prepare()` validates the config and stages every
//! participant's directory; the [`orchestrator`] module provides
//! the `submit_all` / poll / cancel / `wait_until_terminal` flow
//! that fans the participants out across a
//! `valenx_core::Executor` concurrently. Today's `run()` still validates via
//! `precice-tools check`; wiring the orchestrator into the run
//! pipeline so a single Adapter::run drives the whole multi-
//! participant coupling is the next chunk.
//!
//! 1. Parses `[coupling]` + `[[coupling.participant]]` from the
//!    case's `case.toml`.
//! 2. Stages `precice-config.xml` from the case dir into the
//!    workdir.
//! 3. Copies every participant's case directory into a
//!    subdirectory of the workdir so a future orchestrator can
//!    launch them in place.
//! 4. Runs `precice-tools check` on the staged config to catch
//!    malformed XML before the user queues the run.
//! 5. Emits a JSON manifest listing the participants, adapter IDs,
//!    and exchanged fields — the [`orchestrator`] reads this to
//!    dispatch the individual solvers.
//!
//! See [`orchestrator`] for the concurrent participant model.

#![forbid(unsafe_code)]
#![allow(missing_docs)]

pub mod case_input;
pub mod orchestrator;

use std::ffi::OsString;
use std::fs;
use std::path::Path;
use std::time::Duration;

use semver::Version;
use serde::{Deserialize, Serialize};

use valenx_core::{
    adapter_helpers::{find_on_path, first_workdir_match},
    error::RunPhase,
    subprocess, Adapter, AdapterError, AdapterInfo, Capabilities, Capability, Case, LicenseMode,
    Physics, PreparedJob, ProbeReport, RunContext, RunReport, VersionRange,
};
use valenx_fields::{
    artifact::{Artifact, ArtifactKind},
    Results,
};

use crate::case_input::{CouplingInput, Participant};

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(PreciceAdapter::new())
}

pub struct PreciceAdapter;

impl PreciceAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for PreciceAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "precice";
const BINARIES: &[&str] = &["precice-tools"];
pub const MANIFEST_FILENAME: &str = "valenx_coupling.json";

/// Manifest the meta-adapter writes for the run orchestrator.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct CouplingManifest {
    pub valenx_adapter: String,
    pub config: String,
    pub max_coupling_iterations: u32,
    pub participants: Vec<ParticipantManifest>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ParticipantManifest {
    pub name: String,
    pub adapter_id: String,
    pub staged_dir: String,
    pub writes: Vec<String>,
    pub reads: Vec<String>,
}

impl Adapter for PreciceAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "preCICE",
            version_range: VersionRange {
                min_inclusive: Version::new(3, 0, 0),
                max_exclusive: Version::new(4, 0, 0),
            },
            physics: &[Physics::MultiPhysics],
            license_mode: LicenseMode::DynamicLinked,
            tool_license: "LGPL-3.0-or-later",
            docs_url: "https://precice.org/docs.html",
            homepage_url: "https://precice.org/",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        // Preferred probe: the `precice-tools` CLI that ships with
        // preCICE 3.x. If it's absent, accept the libprecice.so
        // being discoverable in a future iteration; for now the
        // absence of `precice-tools` is ToolNotInstalled.
        match find_on_path(BINARIES) {
            Some(binary_path) => {
                let found_version = valenx_core::adapter_helpers::detect_tool_version_semver(
                    &binary_path,
                    &["--version", "-v"],
                );
                Ok(ProbeReport {
                    ok: true,
                    found_version,
                    binary_path: Some(binary_path),
                    warnings: vec!["preCICE meta-adapter: full participant orchestration \
                     lands in the Phase 9 tail per RFC 0007; today the \
                     adapter stages the config + manifest and validates \
                     the XML"
                        .into()],
                    required_env: Vec::new(),
                })
            }
            None => Err(AdapterError::ToolNotInstalled {
                name: INFO_ID,
                hint: "preCICE 3.0+ required; install from precice.org \
                       (precice-tools on PATH)"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let (_header, input) = CouplingInput::from_case_dir(&case.path)?;

        fs::create_dir_all(workdir)?;

        // 1. Stage the preCICE XML config.
        let config_src = if input.config_path.is_absolute() {
            input.config_path.clone()
        } else {
            valenx_core::adapter_helpers::confined_join(&case.path, &input.config_path)?
        };
        if !config_src.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[coupling] config {} not found (resolved {})",
                    input.config_path.display(),
                    config_src.display()
                ),
            });
        }
        let config_dst = workdir.join(
            input
                .config_path
                .file_name()
                .map(|s| s.to_string_lossy().to_string())
                .unwrap_or_else(|| "precice-config.xml".to_string()),
        );
        if config_src != config_dst {
            fs::copy(&config_src, &config_dst)?;
        }

        // 2. Stage each participant's case directory.
        let mut participant_manifests: Vec<ParticipantManifest> = Vec::new();
        for p in &input.participants {
            let participant_dir = stage_participant(case, workdir, p)?;
            participant_manifests.push(ParticipantManifest {
                name: p.name.clone(),
                adapter_id: p.adapter_id.clone(),
                staged_dir: participant_dir
                    .file_name()
                    .map(|s| s.to_string_lossy().to_string())
                    .unwrap_or_default(),
                writes: p.writes.clone(),
                reads: p.reads.clone(),
            });
        }

        // 3. Write the manifest the Phase 9-tail run orchestrator
        //    consumes.
        let manifest = CouplingManifest {
            valenx_adapter: INFO_ID.into(),
            config: config_dst
                .file_name()
                .map(|s| s.to_string_lossy().to_string())
                .unwrap_or_default(),
            max_coupling_iterations: input.max_coupling_iterations,
            participants: participant_manifests,
        };
        let manifest_path = workdir.join(MANIFEST_FILENAME);
        let manifest_bytes = serde_json::to_vec_pretty(&manifest).map_err(|e| {
            AdapterError::Other(anyhow::anyhow!("serialise coupling manifest: {e}"))
        })?;
        valenx_core::io_caps::atomic_write_bytes(&manifest_path, &manifest_bytes)?;

        // 4. Run `precice-tools check config.xml` if the tool is on
        //    PATH — that's the current run action. The real
        //    orchestrator slot opens when Phase 9 tail ships.
        let binary_path = find_on_path(BINARIES).ok_or_else(|| AdapterError::ToolNotInstalled {
            name: INFO_ID,
            hint: "no precice-tools on PATH; install preCICE 3.0+".into(),
        })?;
        let config_arg = config_dst
            .file_name()
            .map(|s| s.to_os_string())
            .unwrap_or_else(|| OsString::from("precice-config.xml"));
        let native_command: Vec<OsString> = vec![
            binary_path.into_os_string(),
            OsString::from("check"),
            config_arg,
        ];

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            estimated_runtime: Some(Duration::from_secs(30)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "validating preCICE config", |line| {
            let mut hint = subprocess::Hint::default();
            if line.contains("CONFIGURATION") {
                hint.progress = Some((50.0, line.to_string()));
            }
            if line.contains("ERROR") || line.contains("WARNING") {
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
        // case_path: the staged precice-config.xml — that's the
        // canonical hashable input that defines the coupling. The
        // Valenx manifest (MANIFEST_FILENAME) is a derived artifact,
        // not a user input, so it sits in the mesh slot instead —
        // captures the participant list / field exchanges that the
        // orchestrator dispatches against.
        let case_path = first_workdir_match(&job.workdir, &["xml"])
            .unwrap_or_else(|| job.workdir.join("(no-config-xml-found)"));
        let manifest_path = job.workdir.join(MANIFEST_FILENAME);
        let mesh_path = if manifest_path.is_file() {
            Some(manifest_path)
        } else {
            None
        };
        let prov = valenx_core::adapter_helpers::live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "preCICE",
            "unknown",
            &case_path,
            mesh_path.as_deref(),
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);

        let manifest = job.workdir.join(MANIFEST_FILENAME);
        if manifest.is_file() {
            results.artifacts.push(Artifact {
                path: manifest,
                kind: ArtifactKind::Other,
                checksum: None,
                label: "preCICE coupling manifest (Valenx)".into(),
            });
        }
        if let Ok(entries) = fs::read_dir(&job.workdir) {
            for entry in entries.flatten() {
                let path = entry.path();
                let Some(ext) = path
                    .extension()
                    .and_then(|s| s.to_str())
                    .map(|s| s.to_ascii_lowercase())
                else {
                    continue;
                };
                if ext == "xml" {
                    results.artifacts.push(Artifact {
                        path: path.clone(),
                        kind: ArtifactKind::Other,
                        checksum: None,
                        label: "preCICE config XML".into(),
                    });
                }
            }
        }

        results.artifacts.sort_by(|a, b| a.path.cmp(&b.path));
        Ok(results)
    }

    fn capabilities(&self) -> Capabilities {
        Capabilities {
            capabilities: vec![
                Capability::CouplingFluidStructure,
                Capability::CouplingConjugateHeat,
                Capability::CouplingReactiveFlow,
            ],
            ribbon_contributions: vec![
                "coupling.precice.fsi",
                "coupling.precice.cht",
                "coupling.precice.reactive",
            ],
        }
    }
}

fn stage_participant(
    case: &Case,
    workdir: &Path,
    p: &Participant,
) -> Result<std::path::PathBuf, AdapterError> {
    // Round-9 hardening: relative `case_dir` flows into a recursive
    // copy into the workdir; wrap with `confined_join` so a hostile
    // case can't aim it at `../../etc`.
    let src = if p.case_dir.is_absolute() {
        p.case_dir.clone()
    } else {
        valenx_core::adapter_helpers::confined_join(&case.path, &p.case_dir)?
    };
    let target_name = p
        .case_dir
        .file_name()
        .or_else(|| Some(std::ffi::OsStr::new(&p.name)))
        .map(|s| s.to_string_lossy().to_string())
        .unwrap_or_else(|| p.name.clone());
    let dst = workdir.join(target_name);

    if !src.is_dir() {
        // Create an empty placeholder rather than failing — some
        // participants' cases are generated at run time. The
        // orchestrator will see an empty dir and surface a sensible
        // error when the participant adapter opens it.
        fs::create_dir_all(&dst)?;
        return Ok(dst);
    }
    // Round-9 M25: drop the local recursive copy in favour of the
    // shared `valenx_core::adapter_helpers::copy_dir_recursive`, which
    // refuses symlinks (so a poisoned participant case can't link
    // `mesh/` to `/etc`), caps recursion depth (kills pathological
    // cycles), and caps per-file size (kills disk-fill attacks).
    valenx_core::adapter_helpers::copy_dir_recursive(&src, &dst)?;
    Ok(dst)
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_core::adapter_helpers::sha256_hex_file;
    use valenx_test_utils::tempdir;

    #[test]
    fn info_is_multi_physics() {
        let info = PreciceAdapter::new().info();
        assert_eq!(info.id, "precice");
        assert_eq!(info.physics, &[Physics::MultiPhysics]);
    }

    #[test]
    fn collect_uses_live_provenance_with_real_case_hash() {
        let workdir = tempdir("precice-collect-prov");
        let xml_path = workdir.join("precice-config.xml");
        let xml_bytes = b"<?xml version=\"1.0\"?>\n<precice-configuration/>\n";
        std::fs::write(&xml_path, xml_bytes).expect("write xml");

        let job = PreparedJob {
            workdir: workdir.clone(),
            native_command: Vec::new(),
            environment: Vec::new(),
            estimated_runtime: None,
            kill_on_drop: false,
        };
        let results = PreciceAdapter::new().collect(&job).expect("collect");
        let prov = &results.provenance;

        assert_eq!(prov.adapter, INFO_ID);
        assert!(!prov.adapter_version.is_empty());
        assert_eq!(prov.tool, "preCICE");
        assert!(!prov.run_id.is_empty(), "run_id empty — stub still wired?");
        assert_eq!(prov.case_hash, sha256_hex_file(&xml_path));

        cleanup_lp(&workdir);
    }

    fn cleanup_lp(d: &std::path::Path) {
        let _ = std::fs::remove_dir_all(d);
    }

    #[test]
    fn copy_dir_recursive_replicates_tree() {
        // Round-9 M25: now exercises the shared
        // `valenx_core::adapter_helpers::copy_dir_recursive` (the
        // local copy was deleted); behaviour for plain trees is the
        // same as the round-6 implementation.
        let tmp_src = std::env::temp_dir().join(format!(
            "precice-src-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let tmp_dst = std::env::temp_dir().join(format!(
            "precice-dst-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(tmp_src.join("sub")).unwrap();
        std::fs::write(tmp_src.join("a.txt"), b"one\n").unwrap();
        std::fs::write(tmp_src.join("sub").join("b.txt"), b"two\n").unwrap();
        valenx_core::adapter_helpers::copy_dir_recursive(&tmp_src, &tmp_dst).unwrap();
        assert!(tmp_dst.join("a.txt").is_file());
        assert!(tmp_dst.join("sub").join("b.txt").is_file());
        let _ = std::fs::remove_dir_all(&tmp_src);
        let _ = std::fs::remove_dir_all(&tmp_dst);
    }

    /// Round-9 RED→GREEN: participant `case_dir` used to be joined
    /// with bare `case.path.join`. Wrap with `confined_join`.
    #[test]
    fn prepare_rejects_participant_case_dir_traversing_outside_case_dir() {
        let d = tempdir("precice-participant-trav");
        std::fs::write(d.join("precice-config.xml"), b"<?xml version=\"1.0\"?>\n").unwrap();
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
format  = "1.0"
name    = "trav"
physics = "multi-physics"
solver  = "precice"
mesh    = "(none)"

[coupling]
config = "precice-config.xml"

[[coupling.participant]]
name       = "Fluid"
adapter_id = "openfoam"
case_dir   = "../../etc"
writes     = ["Force"]
reads      = ["Displacement"]

[[coupling.participant]]
name       = "Solid"
adapter_id = "calculix"
case_dir   = "solid"
writes     = ["Displacement"]
reads      = ["Force"]
"#,
        )
        .unwrap();
        let case = Case {
            id: "precice-participant-trav".into(),
            path: d.clone(),
        };
        let workdir = d.join("workdir");
        let err = PreciceAdapter::new().prepare(&case, &workdir).unwrap_err();
        let msg = format!("{err}");
        assert!(
            msg.contains("..") || msg.contains("stay within") || msg.contains("escape"),
            "expected confined_join rejection, got: {msg}"
        );
        let _ = std::fs::remove_dir_all(&d);
    }
}
