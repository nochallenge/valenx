//! Parse the `[chemistry]` section of a Cantera case into a typed
//! [`ChemistryInput`] the Python script generator consumes.
//!
//! **Phase 4 MVP scope:** equilibrium calculator. Given a mixture
//! and thermodynamic state, solve for equilibrium composition via
//! `gas.equilibrate(basis)`. Zero-D reactor networks and 1-D flames
//! extend the enum without replacing it.

use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};
use valenx_core::{AdapterError, CaseDef, CaseHeader};

/// Everything the Python script generator needs.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ChemistryInput {
    pub mechanism: Mechanism,
    pub analysis: Analysis,
    pub initial: ThermoState,
    /// Time-integration parameters for the reactor-network analyses.
    /// Required when [`Analysis::is_reactor`] returns true; ignored
    /// otherwise.
    #[serde(default)]
    pub reactor: Option<ReactorNetwork>,
}

/// Time-integration knobs for the 0-D batch-reactor analyses.
///
/// Cantera's reactor-network solver wraps a stiff ODE integrator that
/// auto-picks step sizes; the user only has to say "how long" and
/// "how many points to record."
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
pub struct ReactorNetwork {
    /// Total integration time in seconds. Auto-ignition delays for
    /// hydrocarbons live in the milliseconds-to-seconds range; pick
    /// `1e-3` for fast premixed cases, `10.0` for slow oxidation.
    pub end_time_s: f64,
    /// Number of evenly-spaced samples to record. Includes both the
    /// initial state (t=0) and the final state (t=end_time_s), so a
    /// minimum of 2 is required.
    pub n_samples: usize,
}

/// Mechanism file reference. Bundled mechanisms (gri30, air, h2o2)
/// resolve by name; user-supplied YAML / CTI / CHEMKIN files are
/// loaded by path.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "kebab-case")]
pub enum Mechanism {
    /// One of Cantera's built-in mechanisms.
    Bundled(String),
    /// External YAML / CTI file.
    External(PathBuf),
}

impl Mechanism {
    /// The Python-side argument to `cantera.Solution(...)`.
    ///
    /// Round-19 H2: both branches now route the user-controlled string
    /// through [`valenx_core::adapter_helpers::python_str_repr`] before
    /// embedding in the literal. Pre-fix, a `mechanism = "x\");\nimport os\nx=\""`
    /// would close the literal and inject statement-level Python the
    /// user controls — the helper escapes `"`, `\`, `\n`, `\r`, `\t`,
    /// and any control byte so the payload stays inside the quoted
    /// region. Path-shaped names also forward-slash-normalise after
    /// escaping (the legacy gmsh convention) so manifests authored on
    /// Windows still resolve on POSIX hosts.
    pub fn as_python_arg(&self) -> String {
        use valenx_core::adapter_helpers::python_str_repr;
        match self {
            Self::Bundled(name) => format!("\"{}\"", python_str_repr(name)),
            Self::External(path) => {
                // Use forward-slash-normalised path, same rule as our
                // gmsh generator. Slash-normalise BEFORE escaping so
                // the `\` → `/` transform doesn't run after the
                // backslash-escape has already doubled them.
                let slashed = path.to_string_lossy().replace('\\', "/");
                format!("\"{}\"", python_str_repr(&slashed))
            }
        }
    }
}

/// Which class of calculation. Each variant maps to a fundamentally
/// different Python code path, so making the enum explicit keeps
/// each branch focused.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum Analysis {
    /// Constant-property equilibrium. Default basis is TP.
    #[default]
    EquilibriumTP,
    /// Constant-enthalpy, constant-pressure equilibrium (adiabatic
    /// flame temperature).
    EquilibriumHP,
    /// Constant-internal-energy, constant-volume equilibrium.
    EquilibriumUV,
    /// 0-D batch reactor at constant pressure. Mass / volume vary,
    /// pressure is held; the integrator marches species
    /// concentrations and temperature through the supplied
    /// `end_time_s` window. Useful for ignition-delay /
    /// flame-speed precursor calculations.
    BatchReactorConstP,
    /// 0-D batch reactor at constant volume. Pressure / temperature
    /// vary, volume is held. Closer to a closed-vessel detonation
    /// experiment; ignition delay times here are nearly identical to
    /// const-P at moderate pressures.
    BatchReactorConstV,
}

impl Analysis {
    pub fn from_str_lenient(s: &str) -> Self {
        match s.to_ascii_lowercase().as_str() {
            "equilibrium-hp" | "hp" | "adiabatic-flame" => Self::EquilibriumHP,
            "equilibrium-uv" | "uv" => Self::EquilibriumUV,
            "batch-reactor-const-p"
            | "batch-reactor-constp"
            | "reactor-const-p"
            | "ignition-delay" => Self::BatchReactorConstP,
            "batch-reactor-const-v" | "batch-reactor-constv" | "reactor-const-v" => {
                Self::BatchReactorConstV
            }
            _ => Self::EquilibriumTP,
        }
    }

    /// Cantera basis string for the equilibrium analyses. Returns a
    /// sentinel `"TP"` for reactor variants — callers should branch
    /// on [`Self::is_reactor`] before using the result for control
    /// flow.
    pub fn equilibrate_basis(self) -> &'static str {
        match self {
            Self::EquilibriumTP => "TP",
            Self::EquilibriumHP => "HP",
            Self::EquilibriumUV => "UV",
            Self::BatchReactorConstP | Self::BatchReactorConstV => "TP",
        }
    }

    /// True for the time-integrated reactor-network analyses. The
    /// Python script generator emits a different code path for these
    /// (no `gas.equilibrate(...)` call; instead a Cantera ReactorNet
    /// + advance loop).
    pub fn is_reactor(self) -> bool {
        matches!(self, Self::BatchReactorConstP | Self::BatchReactorConstV)
    }

    /// Cantera reactor class name for the reactor variants. Returns
    /// `None` for equilibrium analyses where the concept doesn't
    /// apply.
    pub fn reactor_class(self) -> Option<&'static str> {
        match self {
            Self::BatchReactorConstP => Some("IdealGasConstPressureReactor"),
            Self::BatchReactorConstV => Some("IdealGasReactor"),
            _ => None,
        }
    }
}

/// Initial thermodynamic state. Composition is a Cantera-flavoured
/// string (`"CH4:1, O2:2, N2:7.52"`). We don't parse or validate it
/// — Cantera does that at run time and produces much better error
/// messages than we could.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ThermoState {
    pub temperature_k: f64,
    pub pressure_pa: f64,
    /// Composition string per Cantera's X (mole-fraction) syntax.
    pub composition: String,
}

impl ChemistryInput {
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
        if case_def.case.physics != "chemistry" {
            return Err(invalid(format!(
                "cantera adapter only handles physics=\"chemistry\" cases; \
                 got physics=\"{}\"",
                case_def.case.physics
            )));
        }

        let chem = case_def
            .section("chemistry")
            .and_then(|v| v.as_table())
            .ok_or_else(|| invalid("missing [chemistry] section"))?;

        let mechanism = match chem.get("mechanism") {
            Some(v) => {
                let raw = v.as_str().ok_or_else(|| {
                    invalid("[chemistry] mechanism must be a string (name or path)")
                })?;
                // Cantera looks up single-filename mechanisms (e.g.
                // `gri30.yaml`) on its built-in data path, so we
                // only treat a mechanism as External when it carries
                // a directory separator or a relative `./` prefix.
                let looks_like_path =
                    raw.contains('/') || raw.contains('\\') || raw.starts_with('.');
                if looks_like_path {
                    Mechanism::External(PathBuf::from(raw))
                } else {
                    Mechanism::Bundled(raw.to_string())
                }
            }
            None => Mechanism::Bundled("gri30.yaml".to_string()),
        };

        let analysis = chem
            .get("analysis")
            .and_then(|v| v.as_str())
            .map(Analysis::from_str_lenient)
            .unwrap_or_default();

        let initial_tbl = chem
            .get("initial")
            .and_then(|v| v.as_table())
            .ok_or_else(|| invalid("missing [chemistry.initial]"))?;
        let temperature_k = initial_tbl
            .get("T")
            .and_then(as_f64)
            .or_else(|| initial_tbl.get("temperature_k").and_then(as_f64))
            .or_else(|| initial_tbl.get("temperature").and_then(as_f64))
            .ok_or_else(|| invalid("[chemistry.initial] needs T (Kelvin)"))?;
        let pressure_pa = initial_tbl
            .get("P")
            .and_then(as_f64)
            .or_else(|| initial_tbl.get("pressure_pa").and_then(as_f64))
            .or_else(|| initial_tbl.get("pressure").and_then(as_f64))
            .ok_or_else(|| invalid("[chemistry.initial] needs P (Pa)"))?;
        let composition = initial_tbl
            .get("composition")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                invalid(
                    "[chemistry.initial] needs `composition` — \
                     Cantera mole-fraction syntax, e.g. \"CH4:1, O2:2, N2:7.52\"",
                )
            })?
            .to_string();

        // Reactor-network parameters: required when analysis is a
        // reactor variant; ignored otherwise so equilibrium cases
        // don't have to declare a [chemistry.reactor] block.
        let reactor = if analysis.is_reactor() {
            let rt = chem
                .get("reactor")
                .and_then(|v| v.as_table())
                .ok_or_else(|| {
                    invalid(
                        "[chemistry.reactor] block required for reactor analyses — \
                         add `end_time_s = ...` and `n_samples = ...`",
                    )
                })?;
            let end_time_s = rt
                .get("end_time_s")
                .and_then(as_f64)
                .or_else(|| rt.get("end_time").and_then(as_f64))
                .ok_or_else(|| invalid("[chemistry.reactor] needs `end_time_s` (seconds, > 0)"))?;
            // Reject NaN, negatives, and zero — wrap the comparison
            // explicitly so clippy doesn't read `!(x > 0.0)` as a
            // negated-comparison-on-PartialOrd lint.
            if !(end_time_s.is_finite() && end_time_s > 0.0) {
                return Err(invalid(format!(
                    "[chemistry.reactor].end_time_s must be > 0; got {end_time_s}"
                )));
            }
            let n_samples = rt
                .get("n_samples")
                .and_then(|v| v.as_integer())
                .map(|i| i as usize)
                .unwrap_or(100);
            if n_samples < 2 {
                return Err(invalid(format!(
                    "[chemistry.reactor].n_samples must be >= 2; got {n_samples}"
                )));
            }
            Some(ReactorNetwork {
                end_time_s,
                n_samples,
            })
        } else {
            None
        };

        Ok(Self {
            mechanism,
            analysis,
            initial: ThermoState {
                temperature_k,
                pressure_pa,
                composition,
            },
            reactor,
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
name    = "ch4-air"
physics = "chemistry"
solver  = "cantera.equilibrium"
mesh    = "(none)"

[chemistry]
mechanism = "gri30.yaml"
analysis  = "equilibrium-hp"

[chemistry.initial]
T           = 300.0
P           = 101325.0
composition = "CH4:1, O2:2, N2:7.52"
"#,
        )
        .unwrap()
    }

    #[test]
    fn parses_equilibrium_hp_case() {
        let cd = sample_case();
        let input = ChemistryInput::from_case_def(&cd).expect("parse");
        assert_eq!(input.mechanism, Mechanism::Bundled("gri30.yaml".into()));
        assert_eq!(input.analysis, Analysis::EquilibriumHP);
        assert!((input.initial.temperature_k - 300.0).abs() < 1e-6);
        assert!((input.initial.pressure_pa - 101325.0).abs() < 1e-3);
        assert_eq!(input.initial.composition, "CH4:1, O2:2, N2:7.52");
    }

    #[test]
    fn detects_external_mechanism_path() {
        let text = r#"
[case]
format = "1.0"
name = "x"
physics = "chemistry"
solver = "cantera.equilibrium"
mesh = "(none)"

[chemistry]
mechanism = "path/to/custom.yaml"

[chemistry.initial]
T = 300
P = 101325
composition = "N2:1"
"#;
        let cd: CaseDef = toml::from_str(text).unwrap();
        let input = ChemistryInput::from_case_def(&cd).expect("parse");
        match input.mechanism {
            Mechanism::External(path) => assert_eq!(path, PathBuf::from("path/to/custom.yaml")),
            other => panic!("wrong mechanism variant: {other:?}"),
        }
    }

    #[test]
    fn rejects_non_chemistry_physics() {
        let mut cd = sample_case();
        cd.case.physics = "cfd".into();
        assert!(matches!(
            ChemistryInput::from_case_def(&cd),
            Err(AdapterError::InvalidCase { .. })
        ));
    }

    #[test]
    fn missing_temperature_errors() {
        let text = r#"
[case]
format = "1.0"
name = "x"
physics = "chemistry"
solver = "cantera.equilibrium"
mesh = "(none)"

[chemistry]
mechanism = "gri30.yaml"

[chemistry.initial]
P = 101325
composition = "N2:1"
"#;
        let cd: CaseDef = toml::from_str(text).unwrap();
        let err = ChemistryInput::from_case_def(&cd).unwrap_err();
        match err {
            AdapterError::InvalidCase { reason, .. } => {
                assert!(reason.contains("T") || reason.contains("temperature"));
            }
            other => panic!("wrong error: {other:?}"),
        }
    }

    #[test]
    fn analysis_string_parsing_is_lenient() {
        assert_eq!(
            Analysis::from_str_lenient("equilibrium-tp"),
            Analysis::EquilibriumTP
        );
        assert_eq!(
            Analysis::from_str_lenient("adiabatic-flame"),
            Analysis::EquilibriumHP
        );
        assert_eq!(Analysis::from_str_lenient("UV"), Analysis::EquilibriumUV);
        assert_eq!(
            Analysis::from_str_lenient("batch-reactor-const-p"),
            Analysis::BatchReactorConstP
        );
        assert_eq!(
            Analysis::from_str_lenient("ignition-delay"),
            Analysis::BatchReactorConstP
        );
        assert_eq!(
            Analysis::from_str_lenient("batch-reactor-const-v"),
            Analysis::BatchReactorConstV
        );
        assert_eq!(
            Analysis::from_str_lenient("whatever"),
            Analysis::EquilibriumTP
        );
    }

    #[test]
    fn analysis_is_reactor_predicate_only_true_for_reactor_variants() {
        assert!(!Analysis::EquilibriumTP.is_reactor());
        assert!(!Analysis::EquilibriumHP.is_reactor());
        assert!(!Analysis::EquilibriumUV.is_reactor());
        assert!(Analysis::BatchReactorConstP.is_reactor());
        assert!(Analysis::BatchReactorConstV.is_reactor());
    }

    #[test]
    fn analysis_reactor_class_picks_right_cantera_class() {
        assert_eq!(
            Analysis::BatchReactorConstP.reactor_class(),
            Some("IdealGasConstPressureReactor")
        );
        assert_eq!(
            Analysis::BatchReactorConstV.reactor_class(),
            Some("IdealGasReactor")
        );
        // Equilibrium variants don't have a reactor class.
        assert_eq!(Analysis::EquilibriumTP.reactor_class(), None);
        assert_eq!(Analysis::EquilibriumHP.reactor_class(), None);
    }

    #[test]
    fn parses_batch_reactor_case_with_explicit_block() {
        let text = r#"
[case]
format  = "1.0"
name    = "ch4-ignition"
physics = "chemistry"
solver  = "cantera.equilibrium"
mesh    = "(none)"

[chemistry]
mechanism = "gri30.yaml"
analysis  = "batch-reactor-const-p"

[chemistry.initial]
T           = 1200.0
P           = 101325.0
composition = "CH4:1, O2:2, N2:7.52"

[chemistry.reactor]
end_time_s = 0.05
n_samples  = 200
"#;
        let cd: CaseDef = toml::from_str(text).unwrap();
        let input = ChemistryInput::from_case_def(&cd).expect("parse");
        assert_eq!(input.analysis, Analysis::BatchReactorConstP);
        let reactor = input.reactor.expect("reactor block must populate");
        assert!((reactor.end_time_s - 0.05).abs() < 1e-12);
        assert_eq!(reactor.n_samples, 200);
    }

    #[test]
    fn batch_reactor_without_reactor_block_errors() {
        let text = r#"
[case]
format = "1.0"
name = "x"
physics = "chemistry"
solver = "cantera.equilibrium"
mesh = "(none)"

[chemistry]
mechanism = "gri30.yaml"
analysis  = "batch-reactor-const-v"

[chemistry.initial]
T = 300
P = 101325
composition = "N2:1"
"#;
        let cd: CaseDef = toml::from_str(text).unwrap();
        let err = ChemistryInput::from_case_def(&cd).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("[chemistry.reactor]"), "got: {msg}");
    }

    #[test]
    fn batch_reactor_rejects_zero_or_negative_end_time() {
        let text = r#"
[case]
format = "1.0"
name = "x"
physics = "chemistry"
solver = "cantera.equilibrium"
mesh = "(none)"

[chemistry]
mechanism = "gri30.yaml"
analysis  = "batch-reactor-const-p"

[chemistry.initial]
T = 1200
P = 101325
composition = "CH4:1, O2:2"

[chemistry.reactor]
end_time_s = -0.1
"#;
        let cd: CaseDef = toml::from_str(text).unwrap();
        let err = ChemistryInput::from_case_def(&cd).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("end_time_s must be > 0"), "got: {msg}");
    }

    #[test]
    fn batch_reactor_defaults_n_samples_to_100_when_omitted() {
        let text = r#"
[case]
format = "1.0"
name = "x"
physics = "chemistry"
solver = "cantera.equilibrium"
mesh = "(none)"

[chemistry]
mechanism = "gri30.yaml"
analysis  = "batch-reactor-const-p"

[chemistry.initial]
T = 1200
P = 101325
composition = "CH4:1, O2:2"

[chemistry.reactor]
end_time_s = 0.05
"#;
        let cd: CaseDef = toml::from_str(text).unwrap();
        let input = ChemistryInput::from_case_def(&cd).expect("parse");
        assert_eq!(input.reactor.unwrap().n_samples, 100);
    }

    #[test]
    fn equilibrium_case_does_not_require_a_reactor_block() {
        // Equilibrium analyses ignore [chemistry.reactor] — the
        // existing equilibrium cases must keep parsing without one.
        let cd = sample_case();
        let input = ChemistryInput::from_case_def(&cd).expect("parse");
        assert!(input.reactor.is_none());
    }
}
