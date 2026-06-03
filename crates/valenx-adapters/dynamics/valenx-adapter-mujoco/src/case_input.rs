//! Parse the `[dynamics]` section of a MuJoCo case.
//!
//! Phase 8 MVP: load an MJCF / URDF file, optionally inject a
//! constant control signal per actuator, step for `duration_s`
//! seconds, and record state / qpos / qvel / ctrl time series.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};
use valenx_core::{AdapterError, CaseDef, CaseHeader};

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct DynamicsInput {
    pub model: ModelSource,
    pub duration_s: f64,
    pub timestep_s: Option<f64>,
    pub ctrl: BTreeMap<String, f64>,
    pub record_every_s: f64,
    pub initial_qpos: Vec<f64>,
    pub initial_qvel: Vec<f64>,
}

/// Where the MuJoCo model comes from.
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "kebab-case")]
pub enum ModelSource {
    /// Native MJCF `.xml`.
    Mjcf { path: PathBuf },
    /// URDF auto-converted by MuJoCo's loader.
    Urdf { path: PathBuf },
}

impl ModelSource {
    pub fn path(&self) -> &Path {
        match self {
            Self::Mjcf { path } | Self::Urdf { path } => path,
        }
    }
}

impl DynamicsInput {
    pub fn from_case_dir(case_dir: &Path) -> Result<(CaseHeader, Self), AdapterError> {
        let case_toml = case_dir.join("case.toml");
        let text = valenx_core::io_caps::read_capped_to_string(&case_toml, valenx_core::project::loader::MAX_PROJECT_FILE_BYTES as usize)?;
        let case_def: CaseDef = toml::from_str(&text).map_err(|e| AdapterError::InvalidCase {
            case_path: case_toml.clone(),
            reason: format!("parse: {e}"),
        })?;
        let input = Self::from_case_def(&case_def).map_err(|e| with_case_path(e, &case_toml))?;
        Ok((case_def.case, input))
    }

    pub fn from_case_def(case_def: &CaseDef) -> Result<Self, AdapterError> {
        if case_def.case.physics != "robotics" && case_def.case.physics != "dynamics" {
            return Err(invalid(format!(
                "mujoco adapter handles physics=\"robotics\" or \"dynamics\"; \
                 got physics=\"{}\"",
                case_def.case.physics
            )));
        }

        let dyn_tbl = case_def
            .section("dynamics")
            .and_then(|v| v.as_table())
            .ok_or_else(|| invalid("missing [dynamics] section"))?;

        let model = parse_model_source(dyn_tbl)?;

        let duration_s = dyn_tbl.get("duration_s").and_then(as_f64).unwrap_or(5.0);
        let timestep_s = dyn_tbl.get("timestep_s").and_then(as_f64);
        let record_every_s = dyn_tbl
            .get("record_every_s")
            .and_then(as_f64)
            .unwrap_or(0.01);

        let ctrl: BTreeMap<String, f64> = dyn_tbl
            .get("ctrl")
            .and_then(|v| v.as_table())
            .map(|tbl| {
                tbl.iter()
                    .filter_map(|(k, v)| {
                        v.as_float()
                            .or_else(|| v.as_integer().map(|i| i as f64))
                            .map(|f| (k.clone(), f))
                    })
                    .collect()
            })
            .unwrap_or_default();

        let initial_qpos: Vec<f64> = dyn_tbl
            .get("initial_qpos")
            .and_then(|v| v.as_array())
            .map(|arr| arr.iter().filter_map(as_f64).collect())
            .unwrap_or_default();
        let initial_qvel: Vec<f64> = dyn_tbl
            .get("initial_qvel")
            .and_then(|v| v.as_array())
            .map(|arr| arr.iter().filter_map(as_f64).collect())
            .unwrap_or_default();

        Ok(Self {
            model,
            duration_s,
            timestep_s,
            ctrl,
            record_every_s,
            initial_qpos,
            initial_qvel,
        })
    }
}

fn parse_model_source(tbl: &toml::value::Table) -> Result<ModelSource, AdapterError> {
    let model_tbl = tbl
        .get("model")
        .and_then(|v| v.as_table())
        .ok_or_else(|| invalid("missing [dynamics.model]"))?;
    let kind = model_tbl
        .get("kind")
        .and_then(|v| v.as_str())
        .unwrap_or_else(|| {
            // Auto-detect from extension if available.
            model_tbl
                .get("path")
                .and_then(|v| v.as_str())
                .map(|p| if p.ends_with(".urdf") { "urdf" } else { "mjcf" })
                .unwrap_or("mjcf")
        });
    let path = model_tbl
        .get("path")
        .and_then(|v| v.as_str())
        .map(PathBuf::from)
        .ok_or_else(|| invalid("[dynamics.model] needs `path`"))?;
    match kind {
        "mjcf" | "xml" => Ok(ModelSource::Mjcf { path }),
        "urdf" => Ok(ModelSource::Urdf { path }),
        other => Err(invalid(format!(
            "[dynamics.model] unknown kind \"{other}\" (supported: mjcf, urdf)"
        ))),
    }
}

fn as_f64(v: &toml::Value) -> Option<f64> {
    v.as_float().or_else(|| v.as_integer().map(|i| i as f64))
}

fn invalid(reason: impl Into<String>) -> AdapterError {
    AdapterError::InvalidCase {
        case_path: PathBuf::new(),
        reason: reason.into(),
    }
}

fn with_case_path(err: AdapterError, path: &Path) -> AdapterError {
    if let AdapterError::InvalidCase { case_path, reason } = err {
        if case_path.as_os_str().is_empty() {
            AdapterError::InvalidCase {
                case_path: path.to_path_buf(),
                reason,
            }
        } else {
            AdapterError::InvalidCase { case_path, reason }
        }
    } else {
        err
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_case() -> CaseDef {
        toml::from_str(
            r#"
[case]
format  = "1.0"
name    = "pendulum"
physics = "robotics"
solver  = "mujoco.mjcf"
mesh    = "(none)"

[dynamics]
duration_s     = 10.0
record_every_s = 0.05
initial_qpos   = [0.5]
initial_qvel   = [0.0]

[dynamics.model]
kind = "mjcf"
path = "pendulum.xml"

[dynamics.ctrl]
motor = 0.1
"#,
        )
        .unwrap()
    }

    #[test]
    fn parses_mjcf_case() {
        let cd = sample_case();
        let input = DynamicsInput::from_case_def(&cd).expect("parse");
        match &input.model {
            ModelSource::Mjcf { path } => assert_eq!(path, &PathBuf::from("pendulum.xml")),
            other => panic!("wrong variant: {other:?}"),
        }
        assert_eq!(input.initial_qpos, vec![0.5]);
        assert_eq!(input.ctrl.len(), 1);
        assert!((input.ctrl["motor"] - 0.1).abs() < 1e-9);
    }

    #[test]
    fn auto_detects_urdf_extension() {
        let text = r#"
[case]
format = "1.0"
name = "franka"
physics = "robotics"
solver = "mujoco.urdf"
mesh = "(none)"

[dynamics]
[dynamics.model]
path = "arm.urdf"
"#;
        let cd: CaseDef = toml::from_str(text).unwrap();
        let input = DynamicsInput::from_case_def(&cd).unwrap();
        assert!(matches!(input.model, ModelSource::Urdf { .. }));
    }

    #[test]
    fn rejects_non_robotics_physics() {
        let mut cd = sample_case();
        cd.case.physics = "cfd".into();
        assert!(matches!(
            DynamicsInput::from_case_def(&cd),
            Err(AdapterError::InvalidCase { .. })
        ));
    }

    #[test]
    fn dynamics_physics_tag_is_accepted() {
        let mut cd = sample_case();
        cd.case.physics = "dynamics".into();
        DynamicsInput::from_case_def(&cd).expect("accept \"dynamics\"");
    }
}
