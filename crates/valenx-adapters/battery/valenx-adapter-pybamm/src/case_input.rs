//! Parse the `[battery]` section of a PyBaMM case into a typed
//! [`BatteryInput`]. Phase 7 MVP: single-protocol discharge
//! (CC / CCCV) on one cell, reporting voltage + SOC over time.

use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};
use valenx_core::{AdapterError, CaseDef, CaseHeader};

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct BatteryInput {
    pub model: ModelKind,
    pub parameter_set: String,
    pub protocol: Protocol,
    pub initial_soc: f64,
    pub time_horizon_s: f64,
    pub sample_every_s: f64,
}

/// Which electrochemical model to solve.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum ModelKind {
    /// Single-Particle Model — fastest, crude.
    Spm,
    /// SPM with electrolyte — middle ground.
    Spme,
    #[default]
    /// Doyle-Fuller-Newman pseudo-2D — the gold standard.
    Dfn,
}

impl ModelKind {
    pub fn pybamm_class(self) -> &'static str {
        match self {
            Self::Spm => "SPM",
            Self::Spme => "SPMe",
            Self::Dfn => "DFN",
        }
    }

    pub fn from_str_lenient(s: &str) -> Self {
        match s.to_ascii_uppercase().as_str() {
            "SPM" => Self::Spm,
            "SPME" => Self::Spme,
            _ => Self::Dfn,
        }
    }
}

/// A cycling protocol. MVP shapes. Real drive-cycle files land
/// later via a `DriveCycleFile { path }` variant.
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "kebab-case")]
pub enum Protocol {
    /// Constant-current discharge at `current_a` Amperes until the
    /// cell hits `cutoff_voltage_v`.
    CcDischarge {
        current_a: f64,
        cutoff_voltage_v: f64,
    },
    /// Constant-current, constant-voltage charge. `current_a` is the
    /// initial CC current; `voltage_v` is the CV plateau; `end_current_a`
    /// is the termination condition.
    Cccv {
        current_a: f64,
        voltage_v: f64,
        end_current_a: f64,
    },
}

impl BatteryInput {
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
        if case_def.case.physics != "battery" {
            return Err(invalid(format!(
                "pybamm adapter only handles physics=\"battery\" cases; \
                 got physics=\"{}\"",
                case_def.case.physics
            )));
        }

        let battery = case_def
            .section("battery")
            .and_then(|v| v.as_table())
            .ok_or_else(|| invalid("missing [battery] section"))?;

        let model = battery
            .get("model")
            .and_then(|v| v.as_str())
            .map(ModelKind::from_str_lenient)
            .unwrap_or_default();

        let parameter_set = battery
            .get("parameter_set")
            .and_then(|v| v.as_str())
            .unwrap_or("Chen2020")
            .to_string();

        let protocol = parse_protocol(battery)?;

        let initial_soc = battery.get("initial_soc").and_then(as_f64).unwrap_or(1.0);
        let time_horizon_s = battery
            .get("time_horizon_s")
            .and_then(as_f64)
            .or_else(|| {
                battery
                    .get("time_horizon_h")
                    .and_then(as_f64)
                    .map(|h| h * 3600.0)
            })
            .unwrap_or(3600.0);
        let sample_every_s = battery
            .get("sample_every_s")
            .and_then(as_f64)
            .unwrap_or(10.0);

        Ok(Self {
            model,
            parameter_set,
            protocol,
            initial_soc,
            time_horizon_s,
            sample_every_s,
        })
    }
}

fn parse_protocol(battery: &toml::value::Table) -> Result<Protocol, AdapterError> {
    let tbl = battery
        .get("protocol")
        .and_then(|v| v.as_table())
        .ok_or_else(|| invalid("missing [battery.protocol]"))?;
    let kind = tbl
        .get("kind")
        .and_then(|v| v.as_str())
        .unwrap_or("cc-discharge");
    match kind {
        "cc-discharge" | "cc" | "discharge" => Ok(Protocol::CcDischarge {
            current_a: tbl
                .get("current_a")
                .and_then(as_f64)
                .ok_or_else(|| invalid("[battery.protocol] CC needs `current_a`"))?,
            cutoff_voltage_v: tbl.get("cutoff_voltage_v").and_then(as_f64).unwrap_or(2.5),
        }),
        "cccv" | "ccc-v" => Ok(Protocol::Cccv {
            current_a: tbl
                .get("current_a")
                .and_then(as_f64)
                .ok_or_else(|| invalid("[battery.protocol] CCCV needs `current_a`"))?,
            voltage_v: tbl
                .get("voltage_v")
                .and_then(as_f64)
                .ok_or_else(|| invalid("[battery.protocol] CCCV needs `voltage_v`"))?,
            end_current_a: tbl.get("end_current_a").and_then(as_f64).unwrap_or(0.05),
        }),
        other => Err(invalid(format!(
            "[battery.protocol] unknown kind \"{other}\" (supported: cc-discharge, cccv)"
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
name    = "18650-discharge"
physics = "battery"
solver  = "pybamm.dfn"
mesh    = "(none)"

[battery]
model           = "DFN"
parameter_set   = "Chen2020"
initial_soc     = 1.0
time_horizon_h  = 3
sample_every_s  = 30

[battery.protocol]
kind              = "cc-discharge"
current_a         = 5.0
cutoff_voltage_v  = 2.5
"#,
        )
        .unwrap()
    }

    #[test]
    fn parses_cc_discharge() {
        let cd = sample_case();
        let input = BatteryInput::from_case_def(&cd).expect("parse");
        assert_eq!(input.model, ModelKind::Dfn);
        assert_eq!(input.parameter_set, "Chen2020");
        assert!((input.time_horizon_s - 3.0 * 3600.0).abs() < 1e-6);
        match input.protocol {
            Protocol::CcDischarge {
                current_a,
                cutoff_voltage_v,
            } => {
                assert!((current_a - 5.0).abs() < 1e-6);
                assert!((cutoff_voltage_v - 2.5).abs() < 1e-6);
            }
            other => panic!("wrong protocol: {other:?}"),
        }
    }

    #[test]
    fn rejects_non_battery_physics() {
        let mut cd = sample_case();
        cd.case.physics = "fea".into();
        assert!(matches!(
            BatteryInput::from_case_def(&cd),
            Err(AdapterError::InvalidCase { .. })
        ));
    }

    #[test]
    fn model_parsing_is_lenient() {
        assert_eq!(ModelKind::from_str_lenient("spm"), ModelKind::Spm);
        assert_eq!(ModelKind::from_str_lenient("SPMe"), ModelKind::Spme);
        assert_eq!(ModelKind::from_str_lenient("DFN"), ModelKind::Dfn);
        assert_eq!(ModelKind::from_str_lenient("whatever"), ModelKind::Dfn);
    }

    #[test]
    fn cccv_protocol_parses() {
        let text = r#"
[case]
format = "1.0"
name = "charge"
physics = "battery"
solver = "pybamm.dfn"
mesh = "(none)"

[battery]
[battery.protocol]
kind           = "cccv"
current_a      = -2.5
voltage_v      = 4.2
end_current_a  = 0.1
"#;
        let cd: CaseDef = toml::from_str(text).unwrap();
        let input = BatteryInput::from_case_def(&cd).unwrap();
        match input.protocol {
            Protocol::Cccv {
                current_a,
                voltage_v,
                end_current_a,
            } => {
                assert!((current_a + 2.5).abs() < 1e-6);
                assert!((voltage_v - 4.2).abs() < 1e-6);
                assert!((end_current_a - 0.1).abs() < 1e-6);
            }
            other => panic!("wrong protocol: {other:?}"),
        }
    }
}
