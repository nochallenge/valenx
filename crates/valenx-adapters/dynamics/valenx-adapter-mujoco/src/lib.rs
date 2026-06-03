//! # valenx-adapter-mujoco
//!
//! Adapter for MuJoCo multi-body dynamics, driven via its Python
//! bindings. **Phase 8 — live for MJCF / URDF playback with
//! constant control signals and trajectory capture.** Inverse
//! dynamics, controller libraries, and RL-style rollouts extend
//! the case shape later.

#![forbid(unsafe_code)]
#![allow(missing_docs)]

pub mod case_input;
pub mod python_script;
pub mod summary_parser;

use std::ffi::OsString;
use std::fs;
use std::path::Path;
use std::time::Duration;

use semver::Version;

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

use crate::case_input::DynamicsInput;
use crate::python_script::{SCRIPT_FILENAME, SUMMARY_FILENAME, TIMESERIES_FILENAME};

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(MuJoCoAdapter::new())
}

pub struct MuJoCoAdapter;

impl MuJoCoAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for MuJoCoAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "mujoco";
const PYTHON_BINARIES: &[&str] = &["python3", "python"];

impl Adapter for MuJoCoAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "MuJoCo",
            version_range: VersionRange {
                min_inclusive: Version::new(3, 0, 0),
                max_exclusive: Version::new(4, 0, 0),
            },
            physics: &[Physics::Robotics],
            license_mode: LicenseMode::Subprocess,
            tool_license: "Apache-2.0",
            docs_url: "https://mujoco.readthedocs.io/",
            homepage_url: "https://mujoco.org/",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(PYTHON_BINARIES) {
            Some(binary_path) => {
                let found_version = valenx_core::adapter_helpers::detect_tool_version_semver(
                    &binary_path,
                    &["--version", "-V"],
                );
                Ok(ProbeReport {
                    ok: true,
                    found_version,
                    binary_path: Some(binary_path),
                    warnings: vec!["probe checks for Python; MuJoCo bindings (`pip install \
                     mujoco`) must also be importable at run time"
                        .into()],
                    required_env: Vec::new(),
                })
            }
            None => Err(AdapterError::ToolNotInstalled {
                name: INFO_ID,
                hint: "MuJoCo bindings require Python 3.9+; install via \
                       `pip install mujoco`"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let (_header, input) = DynamicsInput::from_case_dir(&case.path)?;

        fs::create_dir_all(workdir)?;

        // Stage model file.
        // Round-9 hardening: `model_source` (extracted from
        // `input.model`) is user-supplied data and gets copied into
        // the workdir; wrap relative paths with `confined_join`.
        let model_source = input.model.path();
        let source_abs = if model_source.is_absolute() {
            model_source.to_path_buf()
        } else {
            valenx_core::adapter_helpers::confined_join(&case.path, model_source)?
        };
        if !source_abs.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[dynamics.model] file {} not found (resolved {})",
                    model_source.display(),
                    source_abs.display()
                ),
            });
        }
        let file_name = model_source
            .file_name()
            .ok_or_else(|| AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!("model path `{}` has no filename", model_source.display()),
            })?;
        let staged = workdir.join(file_name);
        if source_abs != staged {
            fs::copy(&source_abs, &staged)?;
        }

        // Rewrite the input so the generated Python uses the staged
        // filename, not the original absolute / case-relative path.
        let write_input = rewrite_model_to_basename(input);

        let script_path = workdir.join(SCRIPT_FILENAME);
        python_script::write_to_file(&write_input, &script_path)?;

        let binary_path =
            find_on_path(PYTHON_BINARIES).ok_or_else(|| AdapterError::ToolNotInstalled {
                name: INFO_ID,
                hint: "no python / python3 on PATH".into(),
            })?;

        let native_command: Vec<OsString> = vec![
            binary_path.into_os_string(),
            OsString::from(SCRIPT_FILENAME),
        ];

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            estimated_runtime: Some(Duration::from_secs(120)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting MuJoCo", |line| {
            let mut hint = subprocess::Hint::default();
            if line.contains("[valenx] mujoco done") {
                hint.progress = Some((95.0, line.to_string()));
            } else if line.contains("Traceback") || line.contains("Error") {
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
        // case_path: the generated Python script (canonical hashable
        // input — captures the staged-MJCF reference + duration +
        // timestep). mesh_path: the staged MJCF (.xml) or URDF
        // (.urdf) model that the script loads.
        let case_path = job.workdir.join(SCRIPT_FILENAME);
        let mesh_path = first_workdir_match(&job.workdir, &["xml", "urdf"]);
        let prov = valenx_core::adapter_helpers::live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "MuJoCo",
            "unknown",
            &case_path,
            mesh_path.as_deref(),
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);

        let summary_path = job.workdir.join(SUMMARY_FILENAME);
        if summary_path.is_file() {
            if let Ok(summary) = summary_parser::parse_file(&summary_path) {
                results.meta.description = Some(format!(
                    "MuJoCo · nq={} · {} steps over {:.2} s",
                    summary.nq, summary.step_count, summary.duration_s
                ));
            }
            results.artifacts.push(Artifact {
                path: summary_path,
                kind: ArtifactKind::Tabular,
                checksum: None,
                label: "MuJoCo summary".into(),
            });
        }
        let ts_path = job.workdir.join(TIMESERIES_FILENAME);
        if ts_path.is_file() {
            // Parse the JSONL trajectory into per-DOF ScalarRecords
            // keyed by physical time. Lets the report layer chart
            // joint positions / velocities / controls without
            // re-reading the file. Failures are skipped silently;
            // the artifact stays listed.
            //
            // Round-21 L2: bound the read at
            // MAX_MUJOCO_TIMESTEP_BYTES (64 MiB). Pre-fix a long
            // rollout with many DOFs (or a hostile / corrupted
            // workdir) could produce multi-GB JSONL that would
            // slurp into memory before the line parser ran.
            if let Ok(text) = valenx_core::io_caps::read_capped_to_string(
                &ts_path,
                valenx_core::io_caps::MAX_MUJOCO_TIMESTEP_BYTES as usize,
            ) {
                load_mujoco_trajectory_into_results(&mut results, &text);
            }
            results.artifacts.push(Artifact {
                path: ts_path,
                kind: ArtifactKind::Tabular,
                checksum: None,
                label: "MuJoCo trajectory JSONL".into(),
            });
        }
        let script_path = job.workdir.join(SCRIPT_FILENAME);
        if script_path.is_file() {
            results.artifacts.push(Artifact {
                path: script_path,
                kind: ArtifactKind::Other,
                checksum: None,
                label: "MuJoCo script (generated)".into(),
            });
        }

        results.artifacts.sort_by(|a, b| a.path.cmp(&b.path));
        Ok(results)
    }

    fn capabilities(&self) -> Capabilities {
        Capabilities {
            capabilities: vec![Capability::Custom("dynamics.multibody.contact")],
            ribbon_contributions: vec!["dynamics.mujoco.simulate"],
        }
    }
}

/// Parse MuJoCo's `trajectory.jsonl` (one JSON object per line, each
/// `{"t": <s>, "qpos": [...], "qvel": [...], "ctrl": [...]}`) and
/// emit one [`valenx_fields::ScalarRecord`] per (component, timestep)
/// into the catalog.
///
/// Component names are `qpos[i]` / `qvel[i]` / `ctrl[i]`. For a
/// 7-DOF arm logging at every step that's ~21 records per step
/// (which can grow large for long runs); the catalog handles that
/// fine and the report layer can filter by name prefix.
fn load_mujoco_trajectory_into_results(results: &mut valenx_fields::Results, jsonl_text: &str) {
    use valenx_fields::units::SECOND;
    use valenx_fields::ScalarRecord;

    for line in jsonl_text.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        let row: serde_json::Value = match serde_json::from_str(trimmed) {
            Ok(v) => v,
            Err(_) => continue,
        };
        let Some(obj) = row.as_object() else {
            continue;
        };
        let t = match obj.get("t").and_then(|v| v.as_f64()) {
            Some(v) => v,
            None => continue,
        };
        let timekey = valenx_fields::TimeKey::Time {
            value: t,
            units: SECOND,
        };
        // Three known array fields. Each becomes one record per index.
        for (key, _units_label) in [("qpos", "rad"), ("qvel", "rad/s"), ("ctrl", "")] {
            if let Some(arr) = obj.get(key).and_then(|v| v.as_array()) {
                for (i, v) in arr.iter().enumerate() {
                    let Some(value) = v.as_f64() else { continue };
                    results.scalars.insert(ScalarRecord {
                        name: format!("{key}[{i}]"),
                        value,
                        units: valenx_fields::units::DIMENSIONLESS,
                        time: timekey,
                        source: valenx_fields::scalar::ScalarSource::Extracted,
                        description: None,
                    });
                }
            }
        }
    }
}

fn rewrite_model_to_basename(input: DynamicsInput) -> DynamicsInput {
    use case_input::ModelSource;
    let rewritten = match input.model {
        ModelSource::Mjcf { path } => ModelSource::Mjcf {
            path: std::path::PathBuf::from(
                path.file_name()
                    .map(|s| s.to_string_lossy().into_owned())
                    .unwrap_or_default(),
            ),
        },
        ModelSource::Urdf { path } => ModelSource::Urdf {
            path: std::path::PathBuf::from(
                path.file_name()
                    .map(|s| s.to_string_lossy().into_owned())
                    .unwrap_or_default(),
            ),
        },
    };
    DynamicsInput {
        model: rewritten,
        ..input
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_core::adapter_helpers::sha256_hex_file;
    use valenx_test_utils::tempdir;

    #[test]
    fn info_is_robotics_domain() {
        let info = MuJoCoAdapter::new().info();
        assert_eq!(info.id, "mujoco");
        assert_eq!(info.physics, &[Physics::Robotics]);
    }

    #[test]
    fn loads_mujoco_trajectory_into_scalar_catalog() {
        use valenx_fields::provenance::Sha256Hex;
        let prov = valenx_fields::Provenance {
            adapter: "mujoco".into(),
            adapter_version: "0".into(),
            tool: "MuJoCo".into(),
            tool_version: "3".into(),
            case_hash: Sha256Hex::new(""),
            mesh_hash: Sha256Hex::new(""),
            input_hash: Sha256Hex::new(""),
            tools_lock_hash: Sha256Hex::new(""),
            run_id: "00000000-0000-0000-0000-000000000000".into(),
            wall_time_seconds: 0.0,
            completed_at: "1970-01-01T00:00:00Z".into(),
            ancestors: Vec::new(),
        };
        let mut results = valenx_fields::Results::empty("mujoco-test", prov);
        // Two timesteps, 2-DOF system.
        let jsonl = "\
{\"t\":0.0,\"qpos\":[0.0,0.0],\"qvel\":[0.0,0.0],\"ctrl\":[0.5]}
{\"t\":0.01,\"qpos\":[0.001,-0.002],\"qvel\":[0.1,-0.2],\"ctrl\":[0.5]}
";
        super::load_mujoco_trajectory_into_results(&mut results, jsonl);
        // 2 rows × (2 qpos + 2 qvel + 1 ctrl) = 10 records.
        assert_eq!(results.scalars.len(), 10);
        // Component-by-index naming.
        let qpos0_at_t1 = results
            .scalars
            .all("qpos[0]")
            .iter()
            .find(|r| {
                matches!(
                    r.time,
                    valenx_fields::TimeKey::Time { value, .. } if (value - 0.01).abs() < 1e-12
                )
            })
            .expect("qpos[0] at t=0.01");
        assert!((qpos0_at_t1.value - 0.001).abs() < 1e-12);
    }

    #[test]
    fn mujoco_trajectory_skips_malformed_lines() {
        use valenx_fields::provenance::Sha256Hex;
        let prov = valenx_fields::Provenance {
            adapter: "mujoco".into(),
            adapter_version: "0".into(),
            tool: "MuJoCo".into(),
            tool_version: "3".into(),
            case_hash: Sha256Hex::new(""),
            mesh_hash: Sha256Hex::new(""),
            input_hash: Sha256Hex::new(""),
            tools_lock_hash: Sha256Hex::new(""),
            run_id: "00000000-0000-0000-0000-000000000000".into(),
            wall_time_seconds: 0.0,
            completed_at: "1970-01-01T00:00:00Z".into(),
            ancestors: Vec::new(),
        };
        let mut results = valenx_fields::Results::empty("mujoco-test", prov);
        let jsonl = "\
{\"t\":0.0,\"qpos\":[1.0],\"qvel\":[0.0],\"ctrl\":[]}
not-json
{\"t\":0.01,\"qpos\":[2.0],\"qvel\":[0.5],\"ctrl\":[]}
";
        super::load_mujoco_trajectory_into_results(&mut results, jsonl);
        // 2 valid rows × 2 components (qpos + qvel; ctrl is empty) = 4.
        assert_eq!(results.scalars.len(), 4);
    }

    #[test]
    fn rewrite_model_strips_directory() {
        use crate::case_input::ModelSource;
        let input = DynamicsInput {
            model: ModelSource::Mjcf {
                path: std::path::PathBuf::from("models/robot/arm.xml"),
            },
            duration_s: 1.0,
            timestep_s: None,
            ctrl: Default::default(),
            record_every_s: 0.01,
            initial_qpos: Vec::new(),
            initial_qvel: Vec::new(),
        };
        let rewritten = rewrite_model_to_basename(input);
        match rewritten.model {
            ModelSource::Mjcf { path } => {
                assert_eq!(path, std::path::PathBuf::from("arm.xml"));
            }
            other => panic!("wrong variant: {other:?}"),
        }
    }

    #[test]
    fn collect_uses_live_provenance_with_real_case_hash() {
        let workdir = tempdir("mujoco-collect-prov");
        let script_path = workdir.join(SCRIPT_FILENAME);
        let script_bytes = b"import mujoco\n# trivial\n";
        std::fs::write(&script_path, script_bytes).expect("write script");

        let job = PreparedJob {
            workdir: workdir.clone(),
            native_command: Vec::new(),
            environment: Vec::new(),
            estimated_runtime: None,
            kill_on_drop: false,
        };
        let results = MuJoCoAdapter::new().collect(&job).expect("collect");
        let prov = &results.provenance;

        assert_eq!(prov.adapter, INFO_ID);
        assert!(!prov.adapter_version.is_empty());
        assert_eq!(prov.tool, "MuJoCo");
        assert!(!prov.run_id.is_empty(), "run_id empty — stub still wired?");
        assert_eq!(prov.case_hash, sha256_hex_file(&script_path));

        cleanup_lp(&workdir);
    }

    fn cleanup_lp(d: &std::path::Path) {
        let _ = std::fs::remove_dir_all(d);
    }

    /// Round-9 RED→GREEN: `[dynamics.model.path]` (Mjcf/Urdf) used to
    /// be joined with bare `case.path.join`. Wrap with `confined_join`.
    #[test]
    fn prepare_rejects_model_path_traversing_outside_case_dir() {
        use valenx_test_utils::tempdir;
        let d = tempdir("mujoco-model-trav");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
format  = "1.0"
name    = "trav"
physics = "robotics"
solver  = "mujoco"
mesh    = "(none)"

[dynamics]
duration_s     = 1.0
record_every_s = 0.1

[dynamics.model]
kind = "mjcf"
path = "../../etc/passwd"
"#,
        )
        .unwrap();
        let case = Case {
            id: "mujoco-model-trav".into(),
            path: d.clone(),
        };
        let workdir = d.join("workdir");
        let err = MuJoCoAdapter::new().prepare(&case, &workdir).unwrap_err();
        let msg = format!("{err}");
        assert!(
            msg.contains("..") || msg.contains("stay within") || msg.contains("escape"),
            "expected confined_join rejection, got: {msg}"
        );
        let _ = std::fs::remove_dir_all(&d);
    }
}
