//! Parse the `[md]` section of a LAMMPS case into a typed
//! [`LammpsInput`]. The input-deck writer consumes it.
//!
//! **Phase 5 MVP scope:** Lennard-Jones fluid in a box with NVE
//! integration. EAM / ReaxFF / external topologies extend the
//! `Potential` / `Initialization` enums without replacing them.

use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};
use valenx_core::{AdapterError, CaseDef, CaseHeader};

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct LammpsInput {
    pub units: Units,
    pub boundary: [BoundaryCondition; 3],
    pub atom_style: String,
    pub initialization: Initialization,
    pub potential: Potential,
    pub ensemble: Ensemble,
    pub run_steps: u64,
    pub timestep: f64,
    pub initial_temperature: Option<f64>,
    pub thermo_every: u64,
    pub dump_every: u64,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Default, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum Units {
    #[default]
    /// Lennard-Jones reduced units (default). Length σ, energy ε,
    /// mass m, time τ = σ√(m/ε).
    Lj,
    /// SI-ish metals: Å, eV, amu, ps.
    Metal,
    /// CGS-ish: g, cm, s, kcal/mol.
    Real,
    /// Full SI.
    Si,
}

impl Units {
    pub fn from_str_lenient(s: &str) -> Self {
        match s.to_ascii_lowercase().as_str() {
            "metal" => Self::Metal,
            "real" => Self::Real,
            "si" => Self::Si,
            _ => Self::Lj,
        }
    }

    pub fn lammps_keyword(self) -> &'static str {
        match self {
            Self::Lj => "lj",
            Self::Metal => "metal",
            Self::Real => "real",
            Self::Si => "si",
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Default, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum BoundaryCondition {
    #[default]
    /// Periodic.
    P,
    /// Fixed (non-periodic, with wall).
    F,
    /// Shrink-wrapped (non-periodic, box tracks atoms).
    S,
}

impl BoundaryCondition {
    pub fn from_str_lenient(s: &str) -> Self {
        match s.to_ascii_lowercase().as_str() {
            "f" | "fixed" => Self::F,
            "s" | "shrink" | "shrink-wrapped" => Self::S,
            _ => Self::P,
        }
    }

    pub fn lammps_char(self) -> char {
        match self {
            Self::P => 'p',
            Self::F => 'f',
            Self::S => 's',
        }
    }
}

/// How atoms come into existence.
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "kebab-case")]
pub enum Initialization {
    /// `lattice fcc <density>` + `region box block … ` + `create_atoms`.
    LjFccBox {
        density: f64,
        /// Size in units of `σ` (or lattice units for non-LJ).
        size: [f64; 3],
        num_types: u32,
    },
    /// User-provided LAMMPS data file (atoms + velocities).
    ReadData { path: PathBuf },
}

/// Interatomic potential.
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "kebab-case")]
pub enum Potential {
    /// 12-6 Lennard-Jones with standard pair_coeff row.
    LjCut {
        epsilon: f64,
        sigma: f64,
        cutoff: f64,
    },
    /// Embedded-atom method — user supplies the `.eam` / `.eam.alloy`
    /// file inside the case directory.
    Eam {
        path: PathBuf,
        elements: Vec<String>,
    },
}

/// Integration ensemble.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum Ensemble {
    #[default]
    /// Microcanonical — energy, volume, N conserved.
    Nve,
    /// Canonical — T + V + N conserved via Nose-Hoover.
    NvtNose,
    /// Isothermal-isobaric.
    NptParrinelloRahman,
}

impl LammpsInput {
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
        if case_def.case.physics != "molecular-dynamics" {
            return Err(invalid(format!(
                "lammps adapter only handles physics=\"molecular-dynamics\" cases; \
                 got physics=\"{}\"",
                case_def.case.physics
            )));
        }

        let md = case_def
            .section("md")
            .and_then(|v| v.as_table())
            .ok_or_else(|| invalid("missing [md] section"))?;

        let units = md
            .get("units")
            .and_then(|v| v.as_str())
            .map(Units::from_str_lenient)
            .unwrap_or_default();

        let atom_style = md
            .get("atom_style")
            .and_then(|v| v.as_str())
            .unwrap_or("atomic")
            .to_string();

        let boundary = match md.get("boundary").and_then(|v| v.as_array()) {
            Some(arr) if arr.len() == 3 => [0, 1, 2].map(|i| {
                arr[i]
                    .as_str()
                    .map(BoundaryCondition::from_str_lenient)
                    .unwrap_or_default()
            }),
            _ => [
                BoundaryCondition::P,
                BoundaryCondition::P,
                BoundaryCondition::P,
            ],
        };

        let initialization = parse_initialization(md)?;
        let potential = parse_potential(md)?;

        let ensemble = md
            .get("ensemble")
            .and_then(|v| v.as_str())
            .map(|s| match s.to_ascii_lowercase().as_str() {
                "nvt" | "nvt-nose" | "nose" => Ensemble::NvtNose,
                "npt" | "npt-parrinello" => Ensemble::NptParrinelloRahman,
                _ => Ensemble::Nve,
            })
            .unwrap_or_default();

        let timestep = md.get("timestep").and_then(as_f64).unwrap_or(0.005);
        let run_steps = md
            .get("run_steps")
            .and_then(|v| v.as_integer())
            .map(|i| i.max(1) as u64)
            .unwrap_or(1000);
        let initial_temperature = md.get("initial_temperature").and_then(as_f64);
        let thermo_every = md
            .get("thermo_every")
            .and_then(|v| v.as_integer())
            .map(|i| i as u64)
            .unwrap_or(100);
        let dump_every = md
            .get("dump_every")
            .and_then(|v| v.as_integer())
            .map(|i| i as u64)
            .unwrap_or(100);

        Ok(Self {
            units,
            boundary,
            atom_style,
            initialization,
            potential,
            ensemble,
            run_steps,
            timestep,
            initial_temperature,
            thermo_every,
            dump_every,
        })
    }
}

fn parse_initialization(md: &toml::value::Table) -> Result<Initialization, AdapterError> {
    let init_tbl = md
        .get("init")
        .and_then(|v| v.as_table())
        .ok_or_else(|| invalid("missing [md.init] section"))?;
    let kind = init_tbl
        .get("kind")
        .and_then(|v| v.as_str())
        .unwrap_or("lj-fcc-box");
    match kind {
        "lj-fcc-box" | "fcc" => {
            let density = init_tbl.get("density").and_then(as_f64).unwrap_or(0.8442);
            let size = init_tbl
                .get("size")
                .and_then(|v| v.as_array())
                .and_then(|arr| {
                    if arr.len() == 3 {
                        Some([
                            arr[0]
                                .as_float()
                                .or_else(|| arr[0].as_integer().map(|i| i as f64))?,
                            arr[1]
                                .as_float()
                                .or_else(|| arr[1].as_integer().map(|i| i as f64))?,
                            arr[2]
                                .as_float()
                                .or_else(|| arr[2].as_integer().map(|i| i as f64))?,
                        ])
                    } else {
                        None
                    }
                })
                .unwrap_or([10.0, 10.0, 10.0]);
            let num_types = init_tbl
                .get("num_types")
                .and_then(|v| v.as_integer())
                .map(|i| i.max(1) as u32)
                .unwrap_or(1);
            Ok(Initialization::LjFccBox {
                density,
                size,
                num_types,
            })
        }
        "read-data" | "data" => {
            let path = init_tbl
                .get("path")
                .and_then(|v| v.as_str())
                .ok_or_else(|| invalid("[md.init] kind=\"read-data\" needs `path`"))?;
            Ok(Initialization::ReadData {
                path: PathBuf::from(path),
            })
        }
        other => Err(invalid(format!(
            "[md.init] unknown kind \"{other}\" (supported: lj-fcc-box, read-data)"
        ))),
    }
}

fn parse_potential(md: &toml::value::Table) -> Result<Potential, AdapterError> {
    let pot_tbl = md
        .get("potential")
        .and_then(|v| v.as_table())
        .ok_or_else(|| invalid("missing [md.potential] section"))?;
    let kind = pot_tbl
        .get("kind")
        .and_then(|v| v.as_str())
        .unwrap_or("lj-cut");
    match kind {
        "lj-cut" | "lj" => Ok(Potential::LjCut {
            epsilon: pot_tbl.get("epsilon").and_then(as_f64).unwrap_or(1.0),
            sigma: pot_tbl.get("sigma").and_then(as_f64).unwrap_or(1.0),
            cutoff: pot_tbl.get("cutoff").and_then(as_f64).unwrap_or(2.5),
        }),
        "eam" | "eam-alloy" => {
            let path = pot_tbl
                .get("path")
                .and_then(|v| v.as_str())
                .ok_or_else(|| invalid("[md.potential] kind=\"eam\" needs `path`"))?;
            let elements: Vec<String> = pot_tbl
                .get("elements")
                .and_then(|v| v.as_array())
                .map(|arr| {
                    arr.iter()
                        .filter_map(|v| v.as_str())
                        .map(String::from)
                        .collect()
                })
                .unwrap_or_default();
            Ok(Potential::Eam {
                path: PathBuf::from(path),
                elements,
            })
        }
        other => Err(invalid(format!(
            "[md.potential] unknown kind \"{other}\" (supported: lj-cut, eam)"
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
name    = "lj-fluid"
physics = "molecular-dynamics"
solver  = "lammps.nve"
mesh    = "(none)"

[md]
units        = "lj"
atom_style   = "atomic"
boundary     = ["p", "p", "p"]
timestep     = 0.005
run_steps    = 2000
initial_temperature = 1.0
thermo_every = 50
dump_every   = 100
ensemble     = "nve"

[md.init]
kind      = "lj-fcc-box"
density   = 0.8442
size      = [10, 10, 10]
num_types = 1

[md.potential]
kind    = "lj-cut"
epsilon = 1.0
sigma   = 1.0
cutoff  = 2.5
"#,
        )
        .unwrap()
    }

    #[test]
    fn parses_lj_case() {
        let cd = sample_case();
        let input = LammpsInput::from_case_def(&cd).expect("parse");
        assert_eq!(input.units, Units::Lj);
        assert_eq!(input.boundary[0], BoundaryCondition::P);
        assert_eq!(input.ensemble, Ensemble::Nve);
        assert_eq!(input.run_steps, 2000);
        match &input.initialization {
            Initialization::LjFccBox {
                density,
                size,
                num_types,
            } => {
                assert!((density - 0.8442).abs() < 1e-6);
                assert_eq!(size, &[10.0, 10.0, 10.0]);
                assert_eq!(*num_types, 1);
            }
            other => panic!("wrong init: {other:?}"),
        }
        match &input.potential {
            Potential::LjCut { cutoff, .. } => assert!((cutoff - 2.5).abs() < 1e-6),
            other => panic!("wrong potential: {other:?}"),
        }
    }

    #[test]
    fn rejects_non_md_physics() {
        let mut cd = sample_case();
        cd.case.physics = "fea".into();
        assert!(matches!(
            LammpsInput::from_case_def(&cd),
            Err(AdapterError::InvalidCase { .. })
        ));
    }

    #[test]
    fn boundary_parser_is_lenient() {
        assert_eq!(
            BoundaryCondition::from_str_lenient("shrink"),
            BoundaryCondition::S
        );
        assert_eq!(
            BoundaryCondition::from_str_lenient("FIXED"),
            BoundaryCondition::F
        );
        assert_eq!(
            BoundaryCondition::from_str_lenient("whatever"),
            BoundaryCondition::P
        );
    }

    #[test]
    fn read_data_init_requires_path() {
        let text = r#"
[case]
format = "1.0"
name = "x"
physics = "molecular-dynamics"
solver = "lammps.nve"
mesh = "(none)"

[md]
[md.init]
kind = "read-data"

[md.potential]
kind = "lj-cut"
"#;
        let cd: CaseDef = toml::from_str(text).unwrap();
        let err = LammpsInput::from_case_def(&cd).unwrap_err();
        match err {
            AdapterError::InvalidCase { reason, .. } => {
                assert!(reason.contains("path"));
            }
            other => panic!("wrong error: {other:?}"),
        }
    }
}
