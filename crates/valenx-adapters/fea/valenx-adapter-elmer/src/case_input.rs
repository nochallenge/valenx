//! Parse the `[heat]` (or future `[em]`, `[flow]`, `[coupled]`)
//! section of an Elmer case into a typed spec the SIF writer
//! consumes.
//!
//! **Phase 3 MVP scope:** steady-state heat equation with isotropic
//! material and Dirichlet / Neumann boundary conditions. Coupled
//! problems, time-dependent runs, and richer equation solvers
//! (Navier-Stokes, Maxwell, Linear Elasticity) land in follow-ups
//! and will grow the `Equation` enum rather than replace it.

use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};
use valenx_core::{AdapterError, CaseDef, CaseHeader};

/// Every knob the SIF writer needs.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ElmerInput {
    pub equation: Equation,
    pub mesh_dir: PathBuf,
    pub material: Material,
    pub boundaries: Vec<BoundaryCondition>,
    pub simulation: Simulation,
    pub output_basename: String,
    /// Steady vs. transient time-stepping. Steady runs iterate to
    /// convergence; transient runs march `t = 0 → end_time` in
    /// fixed-size `delta_t` steps with the BDF integrator.
    #[serde(default)]
    pub time: TimeMode,
    /// Initial temperature applied to the body at `t = 0` for
    /// transient runs. `None` lets Elmer default to zero (which is
    /// almost never physically right; users typically want 293 K).
    /// Ignored for steady runs.
    #[serde(default)]
    pub initial_temperature: Option<f64>,
}

/// Steady vs. transient. Mirrors the OpenFOAM / CalculiX shape so
/// the canonical `case.toml` schema reads consistently across
/// adapters: a `[<equation>.transient]` block opts in to a transient
/// run with `end_time` + `delta_t`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "kebab-case")]
pub enum TimeMode {
    #[default]
    Steady,
    Transient {
        end_time: f64,
        delta_t: f64,
    },
}

impl TimeMode {
    /// `Simulation Type` value for the SIF Simulation block.
    pub fn simulation_type(self) -> &'static str {
        match self {
            Self::Steady => "Steady State",
            Self::Transient { .. } => "Transient",
        }
    }

    /// Whether the writer should emit `Timestep Sizes` /
    /// `Timestep Intervals` / a real BDF order.
    pub fn is_transient(self) -> bool {
        matches!(self, Self::Transient { .. })
    }
}

/// Which physics the case asks Elmer to solve.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum Equation {
    #[default]
    HeatEquation,
}

/// Isotropic material properties. Heat-equation-relevant fields
/// only today; future equations add their own.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Material {
    pub name: String,
    pub density: f64,
    pub heat_capacity: f64,
    pub heat_conductivity: f64,
}

/// One `Boundary Condition` block — either a fixed temperature or a
/// heat flux.
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "kebab-case")]
pub enum BoundaryCondition {
    /// Dirichlet — fix temperature to a value on the named patch.
    Temperature {
        name: String,
        target: u32,
        value: f64,
    },
    /// Neumann — fix inbound heat flux on the named patch.
    HeatFlux {
        name: String,
        target: u32,
        value: f64,
    },
}

/// Solver-control knobs.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Simulation {
    pub max_output_level: u32,
    pub steady_state_max_iterations: u32,
    pub convergence_tolerance: f64,
}

impl Default for Simulation {
    fn default() -> Self {
        Self {
            max_output_level: 5,
            steady_state_max_iterations: 20,
            convergence_tolerance: 1e-5,
        }
    }
}

impl ElmerInput {
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
        if case_def.case.physics != "fea" && case_def.case.physics != "multi-physics" {
            return Err(invalid(format!(
                "elmer adapter handles physics=\"fea\" or \"multi-physics\" cases; \
                 got physics=\"{}\"",
                case_def.case.physics
            )));
        }

        // Pick the equation block. MVP: heat equation only.
        let heat = case_def.section("heat").and_then(|v| v.as_table());
        let (equation, equation_tbl) = match heat {
            Some(tbl) => (Equation::HeatEquation, tbl),
            None => {
                return Err(invalid(
                    "missing [heat] section — elmer adapter MVP only solves the \
                     heat equation today; richer equations land in follow-ups",
                ));
            }
        };

        let mesh_dir = equation_tbl
            .get("mesh_dir")
            .and_then(|v| v.as_str())
            .map(PathBuf::from)
            .unwrap_or_else(|| PathBuf::from("mesh"));

        let material_tbl = equation_tbl
            .get("material")
            .and_then(|v| v.as_table())
            .ok_or_else(|| invalid("missing [heat.material]"))?;
        let material = Material {
            name: material_tbl
                .get("name")
                .and_then(|v| v.as_str())
                .unwrap_or("material")
                .to_string(),
            density: material_tbl
                .get("density")
                .and_then(as_f64)
                .or_else(|| material_tbl.get("rho").and_then(as_f64))
                .unwrap_or(1.0),
            heat_capacity: material_tbl
                .get("heat_capacity")
                .and_then(as_f64)
                .or_else(|| material_tbl.get("Cp").and_then(as_f64))
                .unwrap_or(1.0),
            heat_conductivity: material_tbl
                .get("heat_conductivity")
                .and_then(as_f64)
                .or_else(|| material_tbl.get("k").and_then(as_f64))
                .ok_or_else(|| invalid("[heat.material] needs `heat_conductivity` (or `k`)"))?,
        };

        let boundaries: Vec<BoundaryCondition> = equation_tbl
            .get("boundaries")
            .and_then(|v| v.as_array())
            .map(|arr| {
                arr.iter()
                    .filter_map(|v| v.as_table())
                    .filter_map(|tbl| parse_boundary(tbl).ok())
                    .collect()
            })
            .unwrap_or_default();

        let sim_tbl = equation_tbl.get("simulation").and_then(|v| v.as_table());
        let simulation = Simulation {
            max_output_level: sim_tbl
                .and_then(|t| t.get("max_output_level"))
                .and_then(|v| v.as_integer())
                .map(|i| i as u32)
                .unwrap_or(5),
            steady_state_max_iterations: sim_tbl
                .and_then(|t| t.get("steady_state_max_iterations"))
                .and_then(|v| v.as_integer())
                .map(|i| i as u32)
                .unwrap_or(20),
            convergence_tolerance: sim_tbl
                .and_then(|t| t.get("convergence_tolerance"))
                .and_then(as_f64)
                .unwrap_or(1e-5),
        };

        let output_basename = equation_tbl
            .get("output_basename")
            .and_then(|v| v.as_str())
            .unwrap_or("case")
            .to_string();

        // Optional `[<equation>.transient]` sub-block opts the run
        // into transient mode. Both fields are required when the
        // block is present — partial config produces an InvalidCase
        // pointing at the offending line rather than silently
        // running steady.
        let time = match equation_tbl.get("transient").and_then(|v| v.as_table()) {
            None => TimeMode::Steady,
            Some(t) => {
                let end_time = t
                    .get("end_time")
                    .and_then(as_f64)
                    .ok_or_else(|| invalid("[heat.transient] needs `end_time`"))?;
                let delta_t = t
                    .get("delta_t")
                    .and_then(as_f64)
                    .ok_or_else(|| invalid("[heat.transient] needs `delta_t`"))?;
                if !(end_time > 0.0 && delta_t > 0.0) {
                    return Err(invalid(format!(
                        "[heat.transient] end_time and delta_t must both be > 0; \
                         got end_time={end_time}, delta_t={delta_t}"
                    )));
                }
                if delta_t > end_time {
                    return Err(invalid(format!(
                        "[heat.transient] delta_t={delta_t} is larger than \
                         end_time={end_time} — that's a single step, almost \
                         certainly a mistake"
                    )));
                }
                TimeMode::Transient { end_time, delta_t }
            }
        };

        // Initial temperature lives at the equation level (not in
        // [heat.transient]) so users can specify it for steady runs
        // too if they want it written into the SIF as a hint, even
        // though Elmer ignores Initial Condition for steady solves.
        let initial_temperature = equation_tbl
            .get("initial_temperature")
            .and_then(as_f64)
            .or_else(|| equation_tbl.get("T0").and_then(as_f64));

        Ok(Self {
            equation,
            mesh_dir,
            material,
            boundaries,
            simulation,
            output_basename,
            time,
            initial_temperature,
        })
    }
}

fn parse_boundary(tbl: &toml::value::Table) -> Result<BoundaryCondition, AdapterError> {
    let name = tbl
        .get("name")
        .and_then(|v| v.as_str())
        .ok_or_else(|| invalid("boundary needs `name`"))?
        .to_string();
    let target = tbl
        .get("target")
        .and_then(|v| v.as_integer())
        .ok_or_else(|| invalid("boundary needs `target` (Elmer boundary id)"))?;
    let target = target.max(1) as u32;

    if let Some(value) = tbl.get("temperature").and_then(as_f64) {
        Ok(BoundaryCondition::Temperature {
            name,
            target,
            value,
        })
    } else if let Some(value) = tbl.get("heat_flux").and_then(as_f64) {
        Ok(BoundaryCondition::HeatFlux {
            name,
            target,
            value,
        })
    } else {
        Err(invalid(
            "boundary needs one of `temperature` or `heat_flux`",
        ))
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
name    = "heat-sink"
physics = "fea"
solver  = "elmer.heat"
mesh    = "primary"

[heat]
mesh_dir = "mesh"
output_basename = "heatsink"

[heat.material]
name               = "aluminium"
density            = 2700.0
heat_capacity      = 897.0
heat_conductivity  = 237.0

[[heat.boundaries]]
name        = "inlet"
target      = 1
temperature = 350.0

[[heat.boundaries]]
name      = "outlet"
target    = 2
heat_flux = 0.0

[heat.simulation]
max_output_level             = 4
steady_state_max_iterations  = 30
convergence_tolerance        = 1e-6
"#,
        )
        .unwrap()
    }

    #[test]
    fn parses_heat_case() {
        let cd = sample_case();
        let input = ElmerInput::from_case_def(&cd).expect("parse");
        assert_eq!(input.equation, Equation::HeatEquation);
        assert_eq!(input.material.name, "aluminium");
        assert!((input.material.heat_conductivity - 237.0).abs() < 1e-6);
        assert_eq!(input.boundaries.len(), 2);
        match &input.boundaries[0] {
            BoundaryCondition::Temperature {
                name,
                target,
                value,
            } => {
                assert_eq!(name, "inlet");
                assert_eq!(*target, 1);
                assert!((*value - 350.0).abs() < 1e-6);
            }
            other => panic!("expected Temperature BC, got {other:?}"),
        }
        match &input.boundaries[1] {
            BoundaryCondition::HeatFlux { value, .. } => {
                assert!(value.abs() < 1e-9);
            }
            other => panic!("expected HeatFlux BC, got {other:?}"),
        }
        assert_eq!(input.simulation.steady_state_max_iterations, 30);
        assert_eq!(input.output_basename, "heatsink");
        // Sample case has no [heat.transient] block, so time defaults
        // to Steady and initial_temperature stays None.
        assert_eq!(input.time, TimeMode::Steady);
        assert!(input.initial_temperature.is_none());
    }

    #[test]
    fn parses_transient_heat_case() {
        let text = r#"
[case]
format  = "1.0"
name    = "cooldown"
physics = "fea"
solver  = "elmer.heat"
mesh    = "primary"

[heat]
mesh_dir = "mesh"
output_basename = "cooldown"
initial_temperature = 600.0

[heat.material]
name              = "steel"
density           = 7850.0
heat_capacity     = 460.0
heat_conductivity = 50.0

[[heat.boundaries]]
name        = "outer"
target      = 1
temperature = 293.15

[heat.transient]
end_time = 60.0
delta_t  = 0.1
"#;
        let cd: CaseDef = toml::from_str(text).unwrap();
        let input = ElmerInput::from_case_def(&cd).expect("parse");
        match input.time {
            TimeMode::Transient { end_time, delta_t } => {
                assert!((end_time - 60.0).abs() < 1e-12);
                assert!((delta_t - 0.1).abs() < 1e-12);
            }
            other => panic!("expected Transient, got {other:?}"),
        }
        assert!(input.time.is_transient());
        assert_eq!(input.time.simulation_type(), "Transient");
        assert_eq!(input.initial_temperature, Some(600.0));
    }

    #[test]
    fn transient_block_rejects_partial_config() {
        let text = r#"
[case]
format = "1.0"
name = "broken"
physics = "fea"
solver = "elmer.heat"
mesh = "primary"

[heat]
[heat.material]
name = "x"
heat_conductivity = 1.0

[heat.transient]
delta_t = 0.1
"#;
        let cd: CaseDef = toml::from_str(text).unwrap();
        let err = ElmerInput::from_case_def(&cd).unwrap_err();
        match err {
            AdapterError::InvalidCase { reason, .. } => {
                assert!(
                    reason.contains("end_time"),
                    "expected end_time in error: {reason}"
                );
            }
            other => panic!("wrong error: {other:?}"),
        }
    }

    #[test]
    fn transient_block_rejects_negative_or_inverted_steps() {
        // delta_t > end_time → almost-certainly-a-mistake
        let text = r#"
[case]
format = "1.0"
name = "broken"
physics = "fea"
solver = "elmer.heat"
mesh = "primary"

[heat]
[heat.material]
name = "x"
heat_conductivity = 1.0

[heat.transient]
end_time = 0.1
delta_t  = 1.0
"#;
        let cd: CaseDef = toml::from_str(text).unwrap();
        let err = ElmerInput::from_case_def(&cd).unwrap_err();
        match err {
            AdapterError::InvalidCase { reason, .. } => {
                assert!(reason.contains("delta_t"), "got: {reason}");
            }
            other => panic!("wrong error: {other:?}"),
        }
    }

    #[test]
    fn t0_alias_works_for_initial_temperature() {
        let text = r#"
[case]
format = "1.0"
name = "warm-start"
physics = "fea"
solver = "elmer.heat"
mesh = "primary"

[heat]
T0 = 293.15

[heat.material]
name = "x"
heat_conductivity = 1.0
"#;
        let cd: CaseDef = toml::from_str(text).unwrap();
        let input = ElmerInput::from_case_def(&cd).expect("parse");
        assert_eq!(input.initial_temperature, Some(293.15));
    }

    #[test]
    fn rejects_non_fea_physics() {
        let mut cd = sample_case();
        cd.case.physics = "cfd".into();
        assert!(matches!(
            ElmerInput::from_case_def(&cd),
            Err(AdapterError::InvalidCase { .. })
        ));
    }

    #[test]
    fn requires_heat_conductivity() {
        let text = r#"
[case]
format = "1.0"
name = "x"
physics = "fea"
solver = "elmer.heat"
mesh = "primary"

[heat]
[heat.material]
name = "mystery"
"#;
        let cd: CaseDef = toml::from_str(text).unwrap();
        let err = ElmerInput::from_case_def(&cd).unwrap_err();
        match err {
            AdapterError::InvalidCase { reason, .. } => {
                assert!(reason.contains("heat_conductivity"));
            }
            other => panic!("expected InvalidCase, got {other:?}"),
        }
    }

    #[test]
    fn multi_physics_tag_is_also_accepted() {
        let mut cd = sample_case();
        cd.case.physics = "multi-physics".into();
        ElmerInput::from_case_def(&cd).expect("accept multi-physics");
    }
}
