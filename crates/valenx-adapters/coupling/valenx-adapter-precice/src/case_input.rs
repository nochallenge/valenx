//! Parse the `[coupling]` section of a preCICE meta-case.
//!
//! Phase 9 MVP scope: the user provides an existing
//! `precice-config.xml` alongside the case, plus declares the
//! participants so the adapter can validate that each participant
//! adapter is registered. The real concurrent orchestration
//! (spawning each participant's solver in parallel against the
//! shared coupling interface) is tracked as the Phase 9 tail and
//! lives behind [`RFC 0007`](../../../rfcs/0007-coupling-adapters.md).

use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};
use valenx_core::{AdapterError, CaseDef, CaseHeader};

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct CouplingInput {
    pub config_path: PathBuf,
    pub participants: Vec<Participant>,
    pub max_coupling_iterations: u32,
}

/// One participating solver + the Valenx adapter that owns it.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Participant {
    /// Name the participant uses inside `precice-config.xml`.
    pub name: String,
    /// Valenx adapter ID (must match a registered `AdapterInfo::id`).
    pub adapter_id: String,
    /// Case directory for this participant, relative to the
    /// coupling case or absolute.
    pub case_dir: PathBuf,
    /// Fields this participant *writes* to the shared interface.
    pub writes: Vec<String>,
    /// Fields this participant *reads* from the shared interface.
    pub reads: Vec<String>,
}

impl CouplingInput {
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
        if case_def.case.physics != "multi-physics" && case_def.case.physics != "coupling" {
            return Err(invalid(format!(
                "precice adapter handles physics=\"multi-physics\" or \"coupling\"; \
                 got physics=\"{}\"",
                case_def.case.physics
            )));
        }

        let coupling = case_def
            .section("coupling")
            .and_then(|v| v.as_table())
            .ok_or_else(|| invalid("missing [coupling] section"))?;

        let config_path = coupling
            .get("config")
            .and_then(|v| v.as_str())
            .map(PathBuf::from)
            .unwrap_or_else(|| PathBuf::from("precice-config.xml"));

        let max_coupling_iterations = coupling
            .get("max_coupling_iterations")
            .and_then(|v| v.as_integer())
            .map(|i| i.max(1) as u32)
            .unwrap_or(100);

        let participants: Vec<Participant> = coupling
            .get("participant")
            .and_then(|v| v.as_array())
            .map(|arr| {
                arr.iter()
                    .filter_map(|v| v.as_table())
                    .filter_map(|tbl| {
                        Some(Participant {
                            name: tbl.get("name").and_then(|v| v.as_str())?.to_string(),
                            adapter_id: tbl.get("adapter_id").and_then(|v| v.as_str())?.to_string(),
                            case_dir: tbl
                                .get("case_dir")
                                .and_then(|v| v.as_str())
                                .map(PathBuf::from)?,
                            writes: tbl
                                .get("writes")
                                .and_then(|v| v.as_array())
                                .map(|arr| {
                                    arr.iter()
                                        .filter_map(|v| v.as_str())
                                        .map(String::from)
                                        .collect()
                                })
                                .unwrap_or_default(),
                            reads: tbl
                                .get("reads")
                                .and_then(|v| v.as_array())
                                .map(|arr| {
                                    arr.iter()
                                        .filter_map(|v| v.as_str())
                                        .map(String::from)
                                        .collect()
                                })
                                .unwrap_or_default(),
                        })
                    })
                    .collect()
            })
            .unwrap_or_default();

        if participants.len() < 2 {
            return Err(invalid(
                "preCICE coupling needs at least two [[coupling.participant]] \
                 entries — the meta-adapter can't orchestrate a single-solver \
                 case",
            ));
        }

        Ok(Self {
            config_path,
            participants,
            max_coupling_iterations,
        })
    }
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
name    = "fsi-flap"
physics = "multi-physics"
solver  = "precice.coupling"
mesh    = "(none)"

[coupling]
config = "precice-config.xml"
max_coupling_iterations = 50

[[coupling.participant]]
name       = "Fluid"
adapter_id = "openfoam"
case_dir   = "./fluid"
writes     = ["Forces"]
reads      = ["Displacement"]

[[coupling.participant]]
name       = "Solid"
adapter_id = "calculix"
case_dir   = "./solid"
writes     = ["Displacement"]
reads      = ["Forces"]
"#,
        )
        .unwrap()
    }

    #[test]
    fn parses_fsi_case() {
        let cd = sample_case();
        let input = CouplingInput::from_case_def(&cd).expect("parse");
        assert_eq!(input.participants.len(), 2);
        assert_eq!(input.participants[0].adapter_id, "openfoam");
        assert_eq!(input.participants[1].writes, vec!["Displacement"]);
        assert_eq!(input.max_coupling_iterations, 50);
    }

    #[test]
    fn rejects_single_participant() {
        let text = r#"
[case]
format = "1.0"
name = "x"
physics = "multi-physics"
solver = "precice.coupling"
mesh = "(none)"

[coupling]

[[coupling.participant]]
name       = "Only"
adapter_id = "openfoam"
case_dir   = "./only"
"#;
        let cd: CaseDef = toml::from_str(text).unwrap();
        let err = CouplingInput::from_case_def(&cd).unwrap_err();
        match err {
            AdapterError::InvalidCase { reason, .. } => {
                assert!(reason.contains("two"));
            }
            other => panic!("wrong error: {other:?}"),
        }
    }

    #[test]
    fn rejects_non_multi_physics_physics() {
        let mut cd = sample_case();
        cd.case.physics = "cfd".into();
        assert!(matches!(
            CouplingInput::from_case_def(&cd),
            Err(AdapterError::InvalidCase { .. })
        ));
    }

    #[test]
    fn coupling_physics_tag_is_accepted() {
        let mut cd = sample_case();
        cd.case.physics = "coupling".into();
        CouplingInput::from_case_def(&cd).expect("accept \"coupling\"");
    }
}
