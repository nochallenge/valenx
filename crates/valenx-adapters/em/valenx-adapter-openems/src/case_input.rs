//! Parse the `[em]` section of an openEMS case into a typed
//! [`EmInput`]. MVP scope: rectangular FDTD volume with a
//! Gaussian-pulse excitation, a single dielectric block, and
//! absorbing boundary conditions. Antennas, dispersive materials,
//! and S-parameter sweeps extend the enums rather than replace
//! them.

use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};
use valenx_core::{AdapterError, CaseDef, CaseHeader};

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct EmInput {
    pub domain: Domain,
    pub excitation: Excitation,
    pub simulation: Simulation,
    pub boundary: BoundaryCondition,
    pub probes: Vec<Probe>,
}

/// Computational domain. `Box` is a rectangular FDTD grid; future
/// variants carry shape primitives for real antennas.
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "kebab-case")]
pub enum Domain {
    Box {
        /// Bounds in metres.
        min: [f64; 3],
        max: [f64; 3],
        /// Grid cell size in metres.
        cell_size: f64,
    },
}

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "kebab-case")]
pub enum Excitation {
    /// Differentiated Gaussian pulse centred at `center_freq`.
    Gauss {
        center_freq_hz: f64,
        bandwidth_hz: f64,
    },
    /// Sinusoidal continuous-wave source.
    Sine { freq_hz: f64 },
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Simulation {
    /// Total simulated time in seconds.
    pub sim_time_s: f64,
    /// End-of-run residual threshold (FDTD energy decay fraction).
    pub end_criterion: f64,
    /// Number of time-harmonic frequencies to post-process.
    pub post_freqs: Vec<f64>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Default, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum BoundaryCondition {
    #[default]
    /// Mur first-order absorbing boundary — cheap but lossy at
    /// glancing angles; fine for most antenna / EMC smoke cases.
    Mur,
    /// Perfectly Matched Layer — richer absorption; slower.
    Pml,
    /// Perfect Electric Conductor (metal wall).
    Pec,
}

impl BoundaryCondition {
    pub fn from_str_lenient(s: &str) -> Self {
        match s.to_ascii_lowercase().as_str() {
            "pml" => Self::Pml,
            "pec" | "metal" => Self::Pec,
            _ => Self::Mur,
        }
    }

    pub fn octave_keyword(self) -> &'static str {
        match self {
            Self::Mur => "MUR",
            Self::Pml => "PML_8",
            Self::Pec => "PEC",
        }
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Probe {
    pub name: String,
    /// Probe position in metres (x, y, z).
    pub position: [f64; 3],
    /// One of E_x / E_y / E_z / H_x / H_y / H_z.
    pub component: String,
}

impl EmInput {
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
        if case_def.case.physics != "electromagnetics" && case_def.case.physics != "em" {
            return Err(invalid(format!(
                "openems adapter only handles physics=\"electromagnetics\" / \"em\" cases; \
                 got physics=\"{}\"",
                case_def.case.physics
            )));
        }

        let em = case_def
            .section("em")
            .and_then(|v| v.as_table())
            .ok_or_else(|| invalid("missing [em] section"))?;

        let domain = parse_domain(em)?;
        let excitation = parse_excitation(em)?;

        let sim_tbl = em.get("simulation").and_then(|v| v.as_table());
        let simulation = Simulation {
            sim_time_s: sim_tbl
                .and_then(|t| t.get("sim_time_s"))
                .and_then(as_f64)
                .unwrap_or(10e-9),
            end_criterion: sim_tbl
                .and_then(|t| t.get("end_criterion"))
                .and_then(as_f64)
                .unwrap_or(1e-3),
            post_freqs: sim_tbl
                .and_then(|t| t.get("post_freqs"))
                .and_then(|v| v.as_array())
                .map(|arr| {
                    arr.iter()
                        .filter_map(|v| v.as_float().or_else(|| v.as_integer().map(|i| i as f64)))
                        .collect()
                })
                .unwrap_or_default(),
        };

        let boundary = em
            .get("boundary")
            .and_then(|v| v.as_str())
            .map(BoundaryCondition::from_str_lenient)
            .unwrap_or_default();

        let probes: Vec<Probe> = em
            .get("probes")
            .and_then(|v| v.as_array())
            .map(|arr| {
                arr.iter()
                    .filter_map(|v| v.as_table())
                    .filter_map(|tbl| {
                        Some(Probe {
                            name: tbl.get("name").and_then(|v| v.as_str())?.to_string(),
                            position: tbl.get("position").and_then(|v| v.as_array()).and_then(
                                |arr| {
                                    if arr.len() == 3 {
                                        Some([
                                            arr[0].as_float().or_else(|| {
                                                arr[0].as_integer().map(|i| i as f64)
                                            })?,
                                            arr[1].as_float().or_else(|| {
                                                arr[1].as_integer().map(|i| i as f64)
                                            })?,
                                            arr[2].as_float().or_else(|| {
                                                arr[2].as_integer().map(|i| i as f64)
                                            })?,
                                        ])
                                    } else {
                                        None
                                    }
                                },
                            )?,
                            component: tbl
                                .get("component")
                                .and_then(|v| v.as_str())
                                .unwrap_or("E_z")
                                .to_string(),
                        })
                    })
                    .collect()
            })
            .unwrap_or_default();

        Ok(Self {
            domain,
            excitation,
            simulation,
            boundary,
            probes,
        })
    }
}

fn parse_domain(em: &toml::value::Table) -> Result<Domain, AdapterError> {
    let tbl = em
        .get("domain")
        .and_then(|v| v.as_table())
        .ok_or_else(|| invalid("missing [em.domain]"))?;
    let kind = tbl.get("kind").and_then(|v| v.as_str()).unwrap_or("box");
    if kind != "box" {
        return Err(invalid(format!(
            "[em.domain] kind \"{kind}\" not supported yet; Phase 6 MVP is \
             rectangular box only"
        )));
    }
    let min = tbl
        .get("min")
        .and_then(|v| v.as_array())
        .and_then(|a| parse_vec3(a.as_slice()))
        .ok_or_else(|| invalid("[em.domain] needs `min = [x, y, z]`"))?;
    let max = tbl
        .get("max")
        .and_then(|v| v.as_array())
        .and_then(|a| parse_vec3(a.as_slice()))
        .ok_or_else(|| invalid("[em.domain] needs `max = [x, y, z]`"))?;
    let cell_size = tbl.get("cell_size").and_then(as_f64).unwrap_or(0.005);
    Ok(Domain::Box {
        min,
        max,
        cell_size,
    })
}

fn parse_excitation(em: &toml::value::Table) -> Result<Excitation, AdapterError> {
    let tbl = em
        .get("excitation")
        .and_then(|v| v.as_table())
        .ok_or_else(|| invalid("missing [em.excitation]"))?;
    let kind = tbl.get("kind").and_then(|v| v.as_str()).unwrap_or("gauss");
    match kind {
        "gauss" | "gaussian" => Ok(Excitation::Gauss {
            center_freq_hz: tbl
                .get("center_freq_hz")
                .and_then(as_f64)
                .or_else(|| tbl.get("center_freq").and_then(as_f64))
                .ok_or_else(|| invalid("[em.excitation] Gauss needs `center_freq_hz`"))?,
            bandwidth_hz: tbl
                .get("bandwidth_hz")
                .and_then(as_f64)
                .or_else(|| tbl.get("bandwidth").and_then(as_f64))
                .ok_or_else(|| invalid("[em.excitation] Gauss needs `bandwidth_hz`"))?,
        }),
        "sine" | "cw" => Ok(Excitation::Sine {
            freq_hz: tbl
                .get("freq_hz")
                .and_then(as_f64)
                .ok_or_else(|| invalid("[em.excitation] Sine needs `freq_hz`"))?,
        }),
        other => Err(invalid(format!(
            "[em.excitation] unknown kind \"{other}\" (supported: gauss, sine)"
        ))),
    }
}

fn parse_vec3(arr: &[toml::Value]) -> Option<[f64; 3]> {
    if arr.len() != 3 {
        return None;
    }
    let mut out = [0.0; 3];
    for (i, v) in arr.iter().enumerate() {
        out[i] = v.as_float().or_else(|| v.as_integer().map(|i| i as f64))?;
    }
    Some(out)
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
name    = "fdtd-box"
physics = "electromagnetics"
solver  = "openems.fdtd"
mesh    = "(none)"

[em]
boundary = "pml"

[em.domain]
kind      = "box"
min       = [-0.05, -0.05, -0.05]
max       = [ 0.05,  0.05,  0.05]
cell_size = 0.002

[em.excitation]
kind           = "gauss"
center_freq_hz = 5e9
bandwidth_hz   = 5e9

[em.simulation]
sim_time_s    = 10e-9
end_criterion = 1e-3
post_freqs    = [1e9, 5e9, 10e9]

[[em.probes]]
name      = "center"
position  = [0, 0, 0]
component = "E_z"
"#,
        )
        .unwrap()
    }

    #[test]
    fn parses_box_fdtd_case() {
        let cd = sample_case();
        let input = EmInput::from_case_def(&cd).expect("parse");
        match input.domain {
            Domain::Box {
                min,
                max,
                cell_size,
            } => {
                assert_eq!(min, [-0.05, -0.05, -0.05]);
                assert_eq!(max, [0.05, 0.05, 0.05]);
                assert!((cell_size - 0.002).abs() < 1e-9);
            }
        }
        match input.excitation {
            Excitation::Gauss {
                center_freq_hz,
                bandwidth_hz,
            } => {
                assert!((center_freq_hz - 5e9).abs() < 1.0);
                assert!((bandwidth_hz - 5e9).abs() < 1.0);
            }
            _ => panic!(),
        }
        assert_eq!(input.boundary, BoundaryCondition::Pml);
        assert_eq!(input.probes.len(), 1);
        assert_eq!(input.simulation.post_freqs.len(), 3);
    }

    #[test]
    fn rejects_non_em_physics() {
        let mut cd = sample_case();
        cd.case.physics = "cfd".into();
        assert!(matches!(
            EmInput::from_case_def(&cd),
            Err(AdapterError::InvalidCase { .. })
        ));
    }

    #[test]
    fn em_short_physics_tag_accepted() {
        let mut cd = sample_case();
        cd.case.physics = "em".into();
        EmInput::from_case_def(&cd).expect("accept \"em\"");
    }

    #[test]
    fn gauss_requires_both_freqs() {
        let text = r#"
[case]
format = "1.0"
name = "x"
physics = "electromagnetics"
solver = "openems.fdtd"
mesh = "(none)"

[em]
[em.domain]
kind = "box"
min = [0, 0, 0]
max = [1, 1, 1]

[em.excitation]
kind = "gauss"
center_freq_hz = 1e9
"#;
        let cd: CaseDef = toml::from_str(text).unwrap();
        let err = EmInput::from_case_def(&cd).unwrap_err();
        match err {
            AdapterError::InvalidCase { reason, .. } => {
                assert!(reason.contains("bandwidth"));
            }
            other => panic!("wrong error: {other:?}"),
        }
    }
}
