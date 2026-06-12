//! Parse the `[structural]` section of a CalculiX case.toml into a
//! typed [`LinearStaticInput`] that the `.inp` writer consumes.
//!
//! Scope (Phase 3 MVP):
//!
//! - Linear static analysis (`*STEP` without `NLGEOM`).
//! - Isotropic elastic material (`E`, `nu`, `rho` optional).
//! - Tet4 + Hex8 element types.
//! - Dirichlet displacement BCs (fix a DOF to a value).
//! - Concentrated loads via `*CLOAD` and pressure loads via
//!   `*DLOAD`.
//!
//! Nonlinear steps, contact, thermal coupling, and frequency
//! analyses land in follow-ups.

use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};
use valenx_core::{AdapterError, CaseDef, CaseHeader};

/// Everything CCX needs for one linear-static step.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct LinearStaticInput {
    pub analysis: AnalysisKind,
    pub material: Material,
    pub mesh_source: PathBuf,
    pub boundaries: Vec<Boundary>,
    pub loads: Vec<Load>,
    pub step: Step,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum AnalysisKind {
    #[default]
    LinearStatic,
    /// Implicit linear dynamic (`*DYNAMIC`). Marches structural
    /// equations through time with a Hilber-Hughes-Taylor integrator;
    /// good for impact, drop-tests, harmonic excitation.
    LinearDynamic,
    Modal,
    /// Steady-state heat conduction (`*HEAT TRANSFER, STEADY STATE`).
    Thermal,
    /// Transient heat conduction (`*HEAT TRANSFER`). Time-marches
    /// the temperature field; good for cool-down / heat-up studies.
    ThermalTransient,
}

impl AnalysisKind {
    pub fn from_str_lenient(s: &str) -> Self {
        match s {
            "linear-static" | "static" | "linear_static" => Self::LinearStatic,
            "linear-dynamic" | "dynamic" | "transient-dynamic" | "linear_dynamic" => {
                Self::LinearDynamic
            }
            "modal" | "frequency" | "eigen" => Self::Modal,
            "thermal" | "heat" | "thermal-steady" | "heat-steady" => Self::Thermal,
            "thermal-transient" | "heat-transient" | "transient-thermal" => Self::ThermalTransient,
            _ => Self::LinearStatic,
        }
    }

    /// CalculiX `*STEP` child that matches this analysis.
    pub fn ccx_card(self) -> &'static str {
        match self {
            Self::LinearStatic => "*STATIC",
            Self::LinearDynamic => "*DYNAMIC",
            Self::Modal => "*FREQUENCY",
            Self::Thermal => "*HEAT TRANSFER, STEADY STATE",
            Self::ThermalTransient => "*HEAT TRANSFER",
        }
    }

    /// Whether the data row after this analysis card needs a
    /// `<time_increment>, <time_total>` line.
    ///
    /// True for the analyses that step through time:
    /// - `*STATIC` uses pseudo-time for load ramping (0.1, 1.0 →
    ///   ten increments to full load).
    /// - `*DYNAMIC` uses real time for transient structural marching.
    /// - `*HEAT TRANSFER` (no STEADY STATE qualifier) uses real time
    ///   for transient conduction.
    ///
    /// False for `*FREQUENCY` (eigen-extraction is single-shot) and
    /// `*HEAT TRANSFER, STEADY STATE` (converges, no time stepping).
    pub fn needs_increment_line(self) -> bool {
        matches!(
            self,
            Self::LinearStatic | Self::LinearDynamic | Self::ThermalTransient
        )
    }

    /// Whether this analysis is real-time-marching (as opposed to
    /// pseudo-time or eigen-extraction). True for `*DYNAMIC` and
    /// transient `*HEAT TRANSFER`. The run loop uses this to decide
    /// progress reporting + convergence semantics.
    pub fn is_time_marching(self) -> bool {
        matches!(self, Self::LinearDynamic | Self::ThermalTransient)
    }
}

/// Isotropic linear-elastic material. Units are SI (Pa for modulus,
/// kg/m^3 for density) consistent with canonical `Units`.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Material {
    pub name: String,
    pub youngs_modulus: f64,
    pub poissons_ratio: f64,
    pub density: Option<f64>,
}

/// Where the nodes + elements come from. Today we support the
/// canonical JSON form that `valenx-adapter-gmsh` writes via
/// `collect()`; a raw `.inp` mesh-only include lands next.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum MeshSourceKind {
    /// `mesh.canonical.json` — a serde-serialised `valenx_mesh::Mesh`.
    Canonical,
    /// Raw CalculiX include — user-provided `mesh.inp` with NODE +
    /// ELEMENT blocks.
    CcxInclude,
}

/// One displacement constraint on a node set.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Boundary {
    /// Node set (must exist as a `*NSET` or be a physical group in
    /// the mesh).
    pub nset: String,
    /// Degrees of freedom to fix. `1..=3` are translation, `4..=6`
    /// are rotation (for shell elements).
    pub dof_start: u8,
    pub dof_end: u8,
    /// Value (displacement magnitude). Usually `0.0` for fixed
    /// supports.
    pub value: f64,
}

/// Concentrated force on a node set.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Load {
    pub nset: String,
    pub dof: u8,
    pub force: f64,
}

/// Step control — increment sizing + output requests.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Step {
    pub time_total: f64,
    pub time_increment: f64,
    pub output_fields: Vec<OutputField>,
    /// Geometric-nonlinearity (large displacement / large rotation)
    /// flag. CalculiX writes this as `*STEP, NLGEOM`. Default false
    /// = small-displacement formulation (matches the linear-elastic
    /// MVP).
    ///
    /// Only meaningful for the structural analyses
    /// (LinearStatic / LinearDynamic). Thermal + Modal ignore the
    /// flag — the inp_writer omits the keyword in those cases.
    #[serde(default)]
    pub nlgeom: bool,
    /// Optional minimum auto-increment size. CalculiX uses this with
    /// `*STATIC, NLGEOM` to bound the iteration step from below;
    /// when None we omit the field and CalculiX picks its default.
    #[serde(default)]
    pub inc_min: Option<f64>,
    /// Optional maximum auto-increment size. Pairs with `inc_min`.
    #[serde(default)]
    pub inc_max: Option<f64>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum OutputField {
    /// Nodal displacements.
    U,
    /// Element stresses.
    S,
    /// Reaction forces.
    Rf,
    /// Temperatures (thermal steps).
    Nt,
}

impl OutputField {
    pub fn ccx_code(self) -> &'static str {
        match self {
            Self::U => "U",
            Self::S => "S",
            Self::Rf => "RF",
            Self::Nt => "NT",
        }
    }
}

impl LinearStaticInput {
    /// Load + parse a case directory.
    pub fn from_case_dir(case_dir: &Path) -> Result<(CaseHeader, Self), AdapterError> {
        let case_toml = case_dir.join("case.toml");
        let text = valenx_core::io_caps::read_capped_to_string(
            &case_toml,
            valenx_core::project::loader::MAX_PROJECT_FILE_BYTES as usize,
        )?;
        let case_def: CaseDef = toml::from_str(&text).map_err(|e| AdapterError::InvalidCase {
            case_path: case_toml.clone(),
            reason: format!("parse: {e}"),
        })?;
        let input = Self::from_case_def(&case_def).map_err(|e| with_case_path(e, &case_toml))?;
        Ok((case_def.case, input))
    }

    pub fn from_case_def(case_def: &CaseDef) -> Result<Self, AdapterError> {
        if case_def.case.physics != "fea" {
            return Err(invalid(format!(
                "calculix adapter only handles physics=\"fea\" cases; got physics=\"{}\"",
                case_def.case.physics
            )));
        }

        let structural = case_def
            .section("structural")
            .and_then(|v| v.as_table())
            .ok_or_else(|| invalid("missing [structural] section"))?;

        let analysis = structural
            .get("analysis")
            .and_then(|v| v.as_str())
            .map(AnalysisKind::from_str_lenient)
            .unwrap_or_default();

        let material_tbl = structural
            .get("material")
            .and_then(|v| v.as_table())
            .ok_or_else(|| invalid("missing [structural.material]"))?;
        let material = Material {
            name: material_tbl
                .get("name")
                .and_then(|v| v.as_str())
                .unwrap_or("material")
                .to_string(),
            youngs_modulus: material_tbl
                .get("youngs_modulus")
                .and_then(as_f64)
                .or_else(|| material_tbl.get("E").and_then(as_f64))
                .ok_or_else(|| invalid("[structural.material] needs `youngs_modulus` (or `E`)"))?,
            poissons_ratio: material_tbl
                .get("poissons_ratio")
                .and_then(as_f64)
                .or_else(|| material_tbl.get("nu").and_then(as_f64))
                .unwrap_or(0.3),
            density: material_tbl
                .get("density")
                .and_then(as_f64)
                .or_else(|| material_tbl.get("rho").and_then(as_f64)),
        };

        let mesh_source = structural
            .get("mesh_source")
            .and_then(|v| v.as_str())
            .map(PathBuf::from)
            .unwrap_or_else(|| PathBuf::from("mesh.canonical.json"));

        let boundaries: Vec<Boundary> = structural
            .get("boundaries")
            .and_then(|v| v.as_array())
            .map(|arr| {
                arr.iter()
                    .filter_map(|v| v.as_table())
                    .map(|tbl| Boundary {
                        nset: tbl
                            .get("nset")
                            .and_then(|v| v.as_str())
                            .unwrap_or("ALL")
                            .to_string(),
                        dof_start: tbl
                            .get("dof_start")
                            .and_then(|v| v.as_integer())
                            .unwrap_or(1)
                            .max(1) as u8,
                        dof_end: tbl.get("dof_end").and_then(|v| v.as_integer()).unwrap_or(3) as u8,
                        value: tbl.get("value").and_then(as_f64).unwrap_or(0.0),
                    })
                    .collect()
            })
            .unwrap_or_default();

        let loads: Vec<Load> = structural
            .get("loads")
            .and_then(|v| v.as_array())
            .map(|arr| {
                arr.iter()
                    .filter_map(|v| v.as_table())
                    .filter_map(|tbl| {
                        Some(Load {
                            nset: tbl.get("nset").and_then(|v| v.as_str())?.to_string(),
                            dof: tbl.get("dof").and_then(|v| v.as_integer())? as u8,
                            force: tbl.get("force").and_then(as_f64)?,
                        })
                    })
                    .collect()
            })
            .unwrap_or_default();

        let step_tbl = structural.get("step").and_then(|v| v.as_table());
        let step = Step {
            time_total: step_tbl
                .and_then(|t| t.get("time_total"))
                .and_then(as_f64)
                .unwrap_or(1.0),
            time_increment: step_tbl
                .and_then(|t| t.get("time_increment"))
                .and_then(as_f64)
                .unwrap_or(1.0),
            output_fields: step_tbl
                .and_then(|t| t.get("output_fields"))
                .and_then(|v| v.as_array())
                .map(|arr| {
                    arr.iter()
                        .filter_map(|v| v.as_str())
                        .map(|s| match s.to_ascii_lowercase().as_str() {
                            "u" | "displacement" => OutputField::U,
                            "s" | "stress" => OutputField::S,
                            "rf" | "reaction" => OutputField::Rf,
                            "nt" | "temperature" => OutputField::Nt,
                            _ => OutputField::U,
                        })
                        .collect()
                })
                .unwrap_or_else(|| vec![OutputField::U, OutputField::S, OutputField::Rf]),
            nlgeom: step_tbl
                .and_then(|t| t.get("nlgeom"))
                .and_then(|v| v.as_bool())
                .unwrap_or(false),
            inc_min: step_tbl.and_then(|t| t.get("inc_min")).and_then(as_f64),
            inc_max: step_tbl.and_then(|t| t.get("inc_max")).and_then(as_f64),
        };

        Ok(Self {
            analysis,
            material,
            mesh_source,
            boundaries,
            loads,
            step,
        })
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
name    = "cantilever"
physics = "fea"
solver  = "calculix.static"
mesh    = "primary"

[structural]
analysis = "linear-static"
mesh_source = "mesh.canonical.json"

[structural.material]
name = "steel"
E    = 210e9
nu   = 0.3
density = 7850.0

[[structural.boundaries]]
nset      = "fixed"
dof_start = 1
dof_end   = 3
value     = 0.0

[[structural.loads]]
nset  = "loaded-face"
dof   = 2
force = -1000.0

[structural.step]
time_total     = 1.0
time_increment = 1.0
output_fields  = ["U", "S", "RF"]
"#,
        )
        .unwrap()
    }

    #[test]
    fn parses_linear_static_case() {
        let cd = sample_case();
        let input = LinearStaticInput::from_case_def(&cd).expect("parse");
        assert_eq!(input.analysis, AnalysisKind::LinearStatic);
        assert!((input.material.youngs_modulus - 210e9).abs() < 1.0);
        assert_eq!(input.material.name, "steel");
        assert_eq!(input.boundaries.len(), 1);
        assert_eq!(input.boundaries[0].nset, "fixed");
        assert_eq!(input.loads.len(), 1);
        assert!((input.loads[0].force - -1000.0).abs() < 1e-9);
        assert_eq!(input.step.output_fields.len(), 3);
    }

    #[test]
    fn rejects_non_fea_physics() {
        let mut cd = sample_case();
        cd.case.physics = "cfd".into();
        assert!(matches!(
            LinearStaticInput::from_case_def(&cd),
            Err(AdapterError::InvalidCase { .. })
        ));
    }

    #[test]
    fn requires_youngs_modulus() {
        let text = r#"
[case]
format = "1.0"
name = "bad"
physics = "fea"
solver = "calculix.static"
mesh = "primary"

[structural]
[structural.material]
name = "mystery"
nu   = 0.3
"#;
        let cd: CaseDef = toml::from_str(text).unwrap();
        let err = LinearStaticInput::from_case_def(&cd).unwrap_err();
        match err {
            AdapterError::InvalidCase { reason, .. } => {
                assert!(reason.contains("youngs_modulus"));
            }
            other => panic!("expected InvalidCase, got {other:?}"),
        }
    }

    #[test]
    fn analysis_kind_is_lenient() {
        assert_eq!(
            AnalysisKind::from_str_lenient("static"),
            AnalysisKind::LinearStatic
        );
        assert_eq!(AnalysisKind::from_str_lenient("modal"), AnalysisKind::Modal);
        assert_eq!(
            AnalysisKind::from_str_lenient("heat"),
            AnalysisKind::Thermal
        );
        // Transient aliases land on the new variants.
        assert_eq!(
            AnalysisKind::from_str_lenient("dynamic"),
            AnalysisKind::LinearDynamic
        );
        assert_eq!(
            AnalysisKind::from_str_lenient("linear-dynamic"),
            AnalysisKind::LinearDynamic
        );
        assert_eq!(
            AnalysisKind::from_str_lenient("thermal-transient"),
            AnalysisKind::ThermalTransient
        );
        assert_eq!(
            AnalysisKind::from_str_lenient("heat-transient"),
            AnalysisKind::ThermalTransient
        );
    }

    #[test]
    fn ccx_card_picks_right_analysis() {
        assert_eq!(AnalysisKind::LinearStatic.ccx_card(), "*STATIC");
        assert_eq!(AnalysisKind::LinearDynamic.ccx_card(), "*DYNAMIC");
        assert_eq!(AnalysisKind::Modal.ccx_card(), "*FREQUENCY");
        assert_eq!(
            AnalysisKind::Thermal.ccx_card(),
            "*HEAT TRANSFER, STEADY STATE"
        );
        assert_eq!(AnalysisKind::ThermalTransient.ccx_card(), "*HEAT TRANSFER");
    }

    #[test]
    fn needs_increment_line_excludes_eigen_and_steady_thermal() {
        assert!(AnalysisKind::LinearStatic.needs_increment_line());
        assert!(AnalysisKind::LinearDynamic.needs_increment_line());
        assert!(AnalysisKind::ThermalTransient.needs_increment_line());
        // No time stepping for these:
        assert!(!AnalysisKind::Modal.needs_increment_line());
        assert!(!AnalysisKind::Thermal.needs_increment_line());
    }

    #[test]
    fn is_time_marching_only_for_real_time_analyses() {
        // Pseudo-time and one-shot analyses are not time-marching.
        assert!(!AnalysisKind::LinearStatic.is_time_marching());
        assert!(!AnalysisKind::Modal.is_time_marching());
        assert!(!AnalysisKind::Thermal.is_time_marching());
        // Real-time-marching analyses:
        assert!(AnalysisKind::LinearDynamic.is_time_marching());
        assert!(AnalysisKind::ThermalTransient.is_time_marching());
    }

    #[test]
    fn parses_dynamic_case_with_transient_step() {
        let text = r#"
[case]
format  = "1.0"
name    = "drop-test"
physics = "fea"
solver  = "calculix.dynamic"
mesh    = "primary"

[structural]
analysis    = "linear-dynamic"
mesh_source = "mesh.canonical.json"

[structural.material]
name    = "steel"
E       = 210e9
nu      = 0.3
density = 7850.0

[[structural.boundaries]]
nset      = "fixed"
dof_start = 1
dof_end   = 3
value     = 0.0

[[structural.loads]]
nset  = "tip"
dof   = 2
force = -1000.0

[structural.step]
time_total     = 0.01
time_increment = 1e-4
output_fields  = ["U", "S"]
"#;
        let cd: CaseDef = toml::from_str(text).unwrap();
        let input = LinearStaticInput::from_case_def(&cd).expect("parse");
        assert_eq!(input.analysis, AnalysisKind::LinearDynamic);
        assert!(input.analysis.needs_increment_line());
        assert!(input.analysis.is_time_marching());
        assert!((input.step.time_total - 0.01).abs() < 1e-12);
        assert!((input.step.time_increment - 1e-4).abs() < 1e-12);
    }

    #[test]
    fn step_nlgeom_defaults_to_false() {
        // Existing linear-static cases without an explicit `nlgeom`
        // key must continue parsing with NLGEOM disabled.
        let text = r#"
[case]
format = "1.0"
name = "x"
physics = "fea"
solver = "calculix.static"
mesh = "(none)"

[structural]
mesh_source = "canonical"

[structural.material]
name = "steel"
youngs_modulus = 2.1e11
poissons_ratio = 0.3

[structural.step]
time_total = 1.0
time_increment = 1.0
"#;
        let cd: CaseDef = toml::from_str(text).unwrap();
        let input = LinearStaticInput::from_case_def(&cd).expect("parse");
        assert!(!input.step.nlgeom);
        assert!(input.step.inc_min.is_none());
        assert!(input.step.inc_max.is_none());
    }

    #[test]
    fn step_picks_up_nlgeom_and_inc_bounds_when_set() {
        let text = r#"
[case]
format = "1.0"
name = "x"
physics = "fea"
solver = "calculix.static"
mesh = "(none)"

[structural]
mesh_source = "canonical"

[structural.material]
name = "steel"
youngs_modulus = 2.1e11
poissons_ratio = 0.3

[structural.step]
time_total = 2.0
time_increment = 0.1
nlgeom = true
inc_min = 1e-6
inc_max = 0.5
"#;
        let cd: CaseDef = toml::from_str(text).unwrap();
        let input = LinearStaticInput::from_case_def(&cd).expect("parse");
        assert!(input.step.nlgeom);
        assert!((input.step.inc_min.unwrap() - 1e-6).abs() < 1e-18);
        assert!((input.step.inc_max.unwrap() - 0.5).abs() < 1e-12);
    }
}
