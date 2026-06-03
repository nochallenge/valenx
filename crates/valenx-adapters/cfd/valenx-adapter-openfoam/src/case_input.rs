//! Parse the canonical `case.toml` (per RFC 0001) into a typed
//! `SimpleFoamInput` that the dict generator consumes.
//!
//! Keeping this layer separate from the generator means the unit
//! tests can assert the parsing logic independently of any file I/O,
//! and an adapter author can inspect a typed case before committing
//! to running it.

use std::collections::BTreeMap;
use std::path::Path;

use serde::{Deserialize, Serialize};
use valenx_core::adapter_helpers::validate_structured_identifier;
use valenx_core::{AdapterError, CaseDef, CaseHeader};

/// Everything an OpenFOAM incompressible solver needs to write a case.
///
/// Despite the name (kept for compatibility) this struct now drives all
/// three live incompressible solvers: `simpleFoam` (steady RANS),
/// `pimpleFoam` (transient PIMPLE, with or without RANS), and `icoFoam`
/// (transient laminar PISO). The `solver` + `time` fields select the
/// flavour; the dict writer in `simple_foam.rs` dispatches on them.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SimpleFoamInput {
    /// Which OpenFOAM binary the case is targeting. Resolved from the
    /// `case.solver` string in `case.toml`.
    pub solver: SolverKind,
    /// Steady vs. transient time-stepping. Constrained by `solver`:
    /// `simpleFoam` requires Steady, the others require Transient.
    pub time: TimeMode,
    /// Steady-state iteration cap (`controlDict.endTime` for simpleFoam).
    /// Ignored when `time` is Transient — the transient block carries
    /// its own end_time / delta_t.
    pub iterations: u64,
    pub residual_target: f64,
    pub turbulence: TurbulenceModel,
    pub schemes: SchemePreset,
    pub fluid: Fluid,
    pub boundaries: BTreeMap<String, Boundary>,
    /// Compressible-fluid properties. Required when `solver` is
    /// `RhoSimpleFoam` (or any future compressible solver); ignored
    /// for the incompressible solvers. Defaults to [`Thermo::air`]
    /// when the `[flow.thermo]` block is omitted.
    pub thermo: Thermo,
    /// Inlet / initial temperature in K. Used to populate `0/T` for
    /// compressible runs and as the default for any `temperature`
    /// boundary that omits an explicit value. 293.15 K (20 °C) is
    /// the standard default.
    pub t_inlet: f64,
}

/// Which OpenFOAM solver this case is for. The string `openfoam.<name>`
/// from `case.toml` parses into one of these variants; unknown names
/// are rejected up-front in `prepare()` so the dict writer never has to
/// handle the empty case.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub enum SolverKind {
    /// Steady-state incompressible RANS (default for legacy cases).
    SimpleFoam,
    /// Transient incompressible PIMPLE (merges PISO and SIMPLE) — handles
    /// laminar or RANS with arbitrary time-step size.
    PimpleFoam,
    /// Transient incompressible laminar PISO. Strictly laminar; the
    /// case parser refuses RANS turbulence with a clear error rather
    /// than silently dropping it.
    IcoFoam,
    /// Steady-state compressible RANS for transonic / supersonic
    /// external aero. Uses `thermophysicalProperties` instead of
    /// `transportProperties`, needs a temperature field at `0/T`,
    /// and emits a perfect-gas + Sutherland transport block.
    RhoSimpleFoam,
}

impl SolverKind {
    /// Parse the `case.solver` string. Accepts `"openfoam.<solver>"`
    /// (the canonical adapter-prefixed form) and bare `"<solver>"`
    /// (convenience for adapter-internal callers that already know
    /// they're talking to OpenFOAM).
    pub fn from_solver_str(solver: &str) -> Option<Self> {
        let bare = solver.strip_prefix("openfoam.").unwrap_or(solver);
        match bare {
            "simpleFoam" => Some(Self::SimpleFoam),
            "pimpleFoam" => Some(Self::PimpleFoam),
            "icoFoam" => Some(Self::IcoFoam),
            "rhoSimpleFoam" => Some(Self::RhoSimpleFoam),
            _ => None,
        }
    }

    /// The OpenFOAM binary this kind invokes. Used both as the
    /// `application` entry in `controlDict` and as the executable
    /// name `prepare()` looks for on PATH.
    pub fn binary(&self) -> &'static str {
        match self {
            Self::SimpleFoam => "simpleFoam",
            Self::PimpleFoam => "pimpleFoam",
            Self::IcoFoam => "icoFoam",
            Self::RhoSimpleFoam => "rhoSimpleFoam",
        }
    }

    /// Whether this solver is steady-state.
    pub fn is_steady(&self) -> bool {
        matches!(self, Self::SimpleFoam | Self::RhoSimpleFoam)
    }

    /// Whether this solver supports RANS turbulence models. icoFoam is
    /// strictly laminar; the others accept any RANS model.
    pub fn supports_rans(&self) -> bool {
        !matches!(self, Self::IcoFoam)
    }

    /// Whether this solver is compressible (needs `thermophysicalProperties`,
    /// a temperature field, and density-based fvSolution blocks).
    /// Today only `rhoSimpleFoam`; future `sonicFoam` / `rhoPimpleFoam` /
    /// `chtMultiRegionFoam` would join.
    pub fn is_compressible(&self) -> bool {
        matches!(self, Self::RhoSimpleFoam)
    }
}

/// Compressible-fluid thermo properties for `thermophysicalProperties`.
/// Defaults match air at room conditions so users can omit the
/// `[flow.thermo]` block and still get a sensible run.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct Thermo {
    /// Molecular weight in g/mol (28.97 for air).
    pub molar_weight: f64,
    /// Specific heat at constant pressure in J/(kg·K) (1005 for air).
    pub cp: f64,
    /// Heat of formation in J/kg (≈ 0 for non-reacting gas).
    pub hf: f64,
    /// Sutherland's reference dynamic viscosity in Pa·s (1.4584e-6 for air).
    pub mu_ref: f64,
    /// Sutherland's reference temperature in K (110.4 for air).
    pub t_ref: f64,
    /// Prandtl number (0.7 for air).
    pub prandtl: f64,
}

impl Thermo {
    /// Air at standard conditions — used as the default when
    /// `[flow.thermo]` isn't supplied.
    pub fn air() -> Self {
        Self {
            molar_weight: 28.97,
            cp: 1005.0,
            hf: 0.0,
            mu_ref: 1.4584e-6,
            t_ref: 110.4,
            prandtl: 0.7,
        }
    }
}

/// Steady vs. transient time-stepping. Steady cases iterate to
/// convergence with no real time scale; transient cases march from
/// `t=0` to `end_time` in fixed-size `delta_t` steps and write fields
/// at `write_interval` intervals.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "kebab-case")]
pub enum TimeMode {
    /// Steady-state — endTime/deltaT are nominal iteration counters.
    Steady,
    /// Transient — fields are written every `write_interval` seconds
    /// from `t=0` to `t=end_time`.
    Transient {
        end_time: f64,
        delta_t: f64,
        write_interval: f64,
    },
}

impl TimeMode {
    /// Default transient parameters when `[solve.transient]` is omitted.
    /// 1-second sim, 1 ms steps, 100 ms snapshots — short enough to
    /// keep test iteration tight, long enough to be interesting.
    pub fn default_transient() -> Self {
        Self::Transient {
            end_time: 1.0,
            delta_t: 1e-3,
            write_interval: 0.1,
        }
    }
}

/// RANS turbulence model. The adapter currently supports a curated
/// subset; anything else falls back to laminar with a warning.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub enum TurbulenceModel {
    Laminar,
    KEpsilon,
    KOmega,
    KOmegaSST,
    SpalartAllmaras,
}

impl TurbulenceModel {
    pub fn from_str_lenient(s: &str) -> Self {
        match s {
            "laminar" => Self::Laminar,
            "kEpsilon" | "k-epsilon" => Self::KEpsilon,
            "kOmega" | "k-omega" => Self::KOmega,
            "kOmegaSST" | "k-omega-sst" | "sst" => Self::KOmegaSST,
            "SpalartAllmaras" | "spalart-allmaras" | "sa" => Self::SpalartAllmaras,
            _ => Self::Laminar,
        }
    }

    pub fn is_rans(&self) -> bool {
        !matches!(self, Self::Laminar)
    }

    /// The OpenFOAM `simulationType` value.
    pub fn simulation_type(&self) -> &'static str {
        if matches!(self, Self::Laminar) {
            "laminar"
        } else {
            "RAS"
        }
    }

    /// The OpenFOAM `RASModel` value used inside
    /// `constant/turbulenceProperties` when `simulationType` is `RAS`.
    pub fn ras_model(&self) -> &'static str {
        match self {
            Self::Laminar => "laminar", // unused
            Self::KEpsilon => "kEpsilon",
            Self::KOmega => "kOmega",
            Self::KOmegaSST => "kOmegaSST",
            Self::SpalartAllmaras => "SpalartAllmaras",
        }
    }
}

/// Discretization scheme preset. The real scheme registry is richer
/// (dozens of combinations); these are the two that cover 80% of
/// steady-state RANS cases.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum SchemePreset {
    /// Robust first-order upwind — converges easily, diffusive.
    UpwindFirstOrder,
    /// Second-order linearUpwind — accurate, needs better mesh.
    LinearSecondOrder,
}

impl SchemePreset {
    pub fn from_str_lenient(s: &str) -> Self {
        match s {
            "upwind-first-order" | "upwind" | "first-order" => Self::UpwindFirstOrder,
            "linear-second-order" | "linearUpwind" | "second-order" => Self::LinearSecondOrder,
            _ => Self::UpwindFirstOrder,
        }
    }
}

/// Fluid properties used by `constant/transportProperties`.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Fluid {
    pub name: String,
    /// Reference density (kg/m^3). Incompressible solvers still need
    /// it for post-processing force coefficients.
    pub rho: f64,
    /// Kinematic viscosity (m^2/s).
    pub nu: f64,
}

/// A boundary patch condition. The full list of OpenFOAM BCs is
/// huge; these five cover typical external-aero and internal-flow
/// setups.
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(tag = "type", rename_all = "kebab-case")]
pub enum Boundary {
    /// Dirichlet velocity, zero-gradient pressure.
    VelocityInlet {
        velocity: [f64; 3],
        #[serde(default)]
        turbulence_intensity: Option<f64>,
    },
    /// Zero-gradient velocity, Dirichlet pressure.
    PressureOutlet { pressure: f64 },
    /// No-slip wall: velocity zero, pressure zero-gradient, wall
    /// function on turbulence.
    NoSlip,
    /// Symmetry plane.
    Symmetry,
    /// 2D empty patch — the mandatory closeout on 2D cases.
    Empty,
}

impl SimpleFoamInput {
    /// Load and parse a case directory's `case.toml`.
    pub fn from_case_dir(case_dir: &Path) -> Result<(CaseHeader, Self), AdapterError> {
        let case_toml = case_dir.join("case.toml");
        let text = valenx_core::io_caps::read_capped_to_string(&case_toml, valenx_core::project::loader::MAX_PROJECT_FILE_BYTES as usize)?;
        let case_def: CaseDef = toml::from_str(&text).map_err(|e| AdapterError::InvalidCase {
            case_path: case_toml.clone(),
            reason: format!("parse: {e}"),
        })?;
        let input = Self::from_case_def(&case_def).map_err(|e| set_case_path(e, &case_toml))?;
        Ok((case_def.case, input))
    }

    /// Pull the simpleFoam-relevant fields out of a parsed `CaseDef`.
    /// Unknown values are coerced to sensible defaults rather than
    /// hard-failing — the generator still runs, and warnings surface
    /// via `AdapterError::InvalidCase` only when an essential field
    /// is missing.
    ///
    /// The returned `InvalidCase` errors carry an empty `case_path`;
    /// call sites that know the source path should rewrite it via
    /// [`set_case_path`] so the UI can click through to the offending
    /// file.
    pub fn from_case_def(case_def: &CaseDef) -> Result<Self, AdapterError> {
        if case_def.case.physics != "cfd" {
            return Err(invalid(format!(
                "openfoam adapter only handles physics=\"cfd\" cases; \
                 got physics=\"{}\"",
                case_def.case.physics
            )));
        }

        // Solver kind determines downstream constraints (steady-only,
        // laminar-only, etc.) so resolve it first.
        let solver = SolverKind::from_solver_str(&case_def.case.solver).ok_or_else(|| {
            invalid(format!(
                "unknown OpenFOAM solver \"{}\" — supported: simpleFoam, \
                 pimpleFoam, icoFoam",
                case_def.case.solver
            ))
        })?;

        let flow = case_def
            .section("flow")
            .and_then(|v| v.as_table())
            .ok_or_else(|| invalid("missing [flow] section"))?;
        let turbulence = flow
            .get("turbulence")
            .and_then(|v| v.as_str())
            .map(TurbulenceModel::from_str_lenient)
            .unwrap_or(TurbulenceModel::Laminar);
        if !solver.supports_rans() && turbulence.is_rans() {
            return Err(invalid(format!(
                "solver \"{}\" is laminar-only; case requested turbulence \
                 model \"{}\". Use pimpleFoam (transient) or simpleFoam \
                 (steady) for RANS, or set turbulence = \"laminar\".",
                solver.binary(),
                flow.get("turbulence")
                    .and_then(|v| v.as_str())
                    .unwrap_or("(unset)")
            )));
        }
        let schemes = flow
            .get("schemes")
            .and_then(|v| v.as_str())
            .map(SchemePreset::from_str_lenient)
            .unwrap_or(SchemePreset::UpwindFirstOrder);
        let fluid_table = flow
            .get("fluid")
            .and_then(|v| v.as_table())
            .ok_or_else(|| invalid("missing [flow.fluid] section"))?;
        let fluid = Fluid {
            name: fluid_table
                .get("name")
                .and_then(|v| v.as_str())
                .unwrap_or("fluid")
                .to_string(),
            rho: fluid_table
                .get("rho")
                .and_then(|v| v.as_float())
                .or_else(|| {
                    fluid_table
                        .get("rho")
                        .and_then(|v| v.as_integer())
                        .map(|i| i as f64)
                })
                .unwrap_or(1.225),
            nu: fluid_table
                .get("nu")
                .and_then(|v| v.as_float())
                .or_else(|| {
                    fluid_table
                        .get("nu")
                        .and_then(|v| v.as_integer())
                        .map(|i| i as f64)
                })
                .unwrap_or(1.5e-5),
        };

        let solve = case_def
            .section("solve")
            .and_then(|v| v.as_table())
            .ok_or_else(|| invalid("missing [solve] section"))?;
        let iterations = solve
            .get("iterations")
            .and_then(|v| v.as_integer())
            .map(|i| i.max(0) as u64)
            .unwrap_or(1000);
        let residual_target = solve
            .get("residual_target")
            .and_then(|v| v.as_float())
            .unwrap_or(1e-5);

        // Time mode. Steady solvers ignore [solve.transient] (with a
        // courteous error if the user combined steady + transient
        // params, which is almost certainly a mistake). Transient
        // solvers default to a sensible 1-second / 1 ms run if the
        // user omits the block, which keeps quick test cases trivial.
        let transient_table = solve.get("transient").and_then(|v| v.as_table());
        let time = if solver.is_steady() {
            if transient_table.is_some() {
                return Err(invalid(format!(
                    "solver \"{}\" is steady-state but case provided a \
                     [solve.transient] block — drop the block or switch \
                     to pimpleFoam/icoFoam.",
                    solver.binary()
                )));
            }
            TimeMode::Steady
        } else {
            match transient_table {
                None => TimeMode::default_transient(),
                Some(t) => {
                    let end_time = t
                        .get("end_time")
                        .and_then(|v| v.as_float())
                        .or_else(|| {
                            t.get("end_time")
                                .and_then(|v| v.as_integer())
                                .map(|i| i as f64)
                        })
                        .unwrap_or(1.0);
                    let delta_t = t
                        .get("delta_t")
                        .and_then(|v| v.as_float())
                        .or_else(|| {
                            t.get("delta_t")
                                .and_then(|v| v.as_integer())
                                .map(|i| i as f64)
                        })
                        .unwrap_or(1e-3);
                    let write_interval = t
                        .get("write_interval")
                        .and_then(|v| v.as_float())
                        .or_else(|| {
                            t.get("write_interval")
                                .and_then(|v| v.as_integer())
                                .map(|i| i as f64)
                        })
                        .unwrap_or(0.1);
                    if !(end_time > 0.0 && delta_t > 0.0 && write_interval > 0.0) {
                        return Err(invalid(format!(
                            "[solve.transient] requires end_time, delta_t, \
                             and write_interval all > 0; got \
                             end_time={end_time}, delta_t={delta_t}, \
                             write_interval={write_interval}"
                        )));
                    }
                    if delta_t > end_time {
                        return Err(invalid(format!(
                            "[solve.transient] delta_t={delta_t} is larger \
                             than end_time={end_time} — that's a single \
                             time step, almost certainly a mistake"
                        )));
                    }
                    TimeMode::Transient {
                        end_time,
                        delta_t,
                        write_interval,
                    }
                }
            }
        };

        let boundaries_table = case_def
            .section("boundaries")
            .and_then(|v| v.as_table())
            .ok_or_else(|| invalid("missing [boundaries] section"))?;
        let mut boundaries: BTreeMap<String, Boundary> = BTreeMap::new();
        for (name, v) in boundaries_table {
            // Round-15 M3: validate the user-supplied boundary key
            // before it flows into ~10 dict files via the
            // `for (name, boundary) in &input.boundaries` loops in
            // simple_foam.rs. Pre-fix a hostile key like
            // `inlet\n}\n\ninjectedPatch\n{\n  type wall;` would
            // silently inject a sibling patch into the boundaryField
            // block of every 0/<field> file. Catching it at case-load
            // time gives the user a clean InvalidCase error mentioning
            // the offending field rather than a corrupted dict that
            // simpleFoam refuses to parse mid-run.
            validate_structured_identifier(name, &format!("boundaries.{name}"))
                .map_err(|e| {
                    if let AdapterError::InvalidCase { reason, .. } = e {
                        invalid(reason)
                    } else {
                        e
                    }
                })?;
            let table = v
                .as_table()
                .ok_or_else(|| invalid(format!("[boundaries.{name}] must be a table")))?;
            let kind = table
                .get("type")
                .and_then(|v| v.as_str())
                .ok_or_else(|| invalid(format!("[boundaries.{name}] missing `type`")))?;
            let boundary = match kind {
                "velocity-inlet" => {
                    let velocity = parse_vec3(table.get("velocity")).ok_or_else(|| {
                        invalid(format!(
                            "[boundaries.{name}] velocity must be a 3-element array"
                        ))
                    })?;
                    let turbulence_intensity =
                        table.get("turbulence_intensity").and_then(|v| v.as_float());
                    Boundary::VelocityInlet {
                        velocity,
                        turbulence_intensity,
                    }
                }
                "pressure-outlet" => {
                    let pressure = table
                        .get("pressure")
                        .and_then(|v| v.as_float())
                        .unwrap_or(0.0);
                    Boundary::PressureOutlet { pressure }
                }
                "no-slip" | "wall" => Boundary::NoSlip,
                "symmetry" => Boundary::Symmetry,
                "empty" => Boundary::Empty,
                other => {
                    return Err(invalid(format!(
                        "[boundaries.{name}] unknown type \"{other}\" \
                         (supported: velocity-inlet, pressure-outlet, no-slip, \
                         symmetry, empty)"
                    )));
                }
            };
            boundaries.insert(name.clone(), boundary);
        }

        // Compressible thermo. Optional `[flow.thermo]` lets users
        // override; default = air. We always parse it (even for
        // incompressible runs) so downstream code can read `input.thermo`
        // without branching.
        let thermo = match flow.get("thermo").and_then(|v| v.as_table()) {
            None => Thermo::air(),
            Some(t) => Thermo {
                molar_weight: t
                    .get("molar_weight")
                    .and_then(as_f64_lenient)
                    .unwrap_or(28.97),
                cp: t.get("cp").and_then(as_f64_lenient).unwrap_or(1005.0),
                hf: t.get("hf").and_then(as_f64_lenient).unwrap_or(0.0),
                mu_ref: t
                    .get("mu_ref")
                    .and_then(as_f64_lenient)
                    .unwrap_or(1.4584e-6),
                t_ref: t.get("t_ref").and_then(as_f64_lenient).unwrap_or(110.4),
                prandtl: t.get("prandtl").and_then(as_f64_lenient).unwrap_or(0.7),
            },
        };
        let t_inlet = flow
            .get("t_inlet")
            .and_then(as_f64_lenient)
            .unwrap_or(293.15);

        Ok(Self {
            solver,
            time,
            iterations,
            residual_target,
            turbulence,
            schemes,
            fluid,
            boundaries,
            thermo,
            t_inlet,
        })
    }
}

fn as_f64_lenient(v: &toml::Value) -> Option<f64> {
    v.as_float().or_else(|| v.as_integer().map(|i| i as f64))
}

/// Build an `AdapterError::InvalidCase` with an empty `case_path` —
/// upstream call sites that know the path rewrite it via
/// [`set_case_path`].
fn invalid(reason: impl Into<String>) -> AdapterError {
    AdapterError::InvalidCase {
        case_path: std::path::PathBuf::new(),
        reason: reason.into(),
    }
}

/// If `err` is an `InvalidCase` with an empty `case_path`, rewrite it
/// to point at `path` so the UI can link users straight to the
/// offending file.
pub fn set_case_path(err: AdapterError, path: &Path) -> AdapterError {
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

fn parse_vec3(v: Option<&toml::Value>) -> Option<[f64; 3]> {
    let arr = v?.as_array()?;
    if arr.len() != 3 {
        return None;
    }
    let mut out = [0.0f64; 3];
    for (i, el) in arr.iter().enumerate() {
        out[i] = el
            .as_float()
            .or_else(|| el.as_integer().map(|x| x as f64))?;
    }
    Some(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_case() -> CaseDef {
        let text = r#"
[case]
format  = "1.0"
name    = "demo"
physics = "cfd"
solver  = "openfoam.simpleFoam"
mesh    = "default"

[flow]
turbulence = "kOmegaSST"
schemes    = "upwind-first-order"

[flow.fluid]
name = "air"
rho  = 1.225
nu   = 1.5e-5

[boundaries.inlet]
type     = "velocity-inlet"
velocity = [50.0, 0.0, 0.0]
turbulence_intensity = 0.05

[boundaries.outlet]
type     = "pressure-outlet"
pressure = 0.0

[boundaries.walls]
type = "no-slip"

[solve]
iterations      = 2000
residual_target = 1e-5
"#;
        toml::from_str(text).expect("parse test case.toml")
    }

    #[test]
    fn parses_canonical_case() {
        let cd = sample_case();
        let input = SimpleFoamInput::from_case_def(&cd).expect("parse");
        assert_eq!(input.solver, SolverKind::SimpleFoam);
        assert_eq!(input.time, TimeMode::Steady);
        assert_eq!(input.iterations, 2000);
        assert!((input.residual_target - 1e-5).abs() < 1e-12);
        assert_eq!(input.turbulence, TurbulenceModel::KOmegaSST);
        assert_eq!(input.schemes, SchemePreset::UpwindFirstOrder);
        assert_eq!(input.fluid.name, "air");
        assert!((input.fluid.nu - 1.5e-5).abs() < 1e-12);
        assert_eq!(input.boundaries.len(), 3);
        match input.boundaries.get("inlet").unwrap() {
            Boundary::VelocityInlet { velocity, .. } => {
                assert_eq!(velocity, &[50.0, 0.0, 0.0]);
            }
            _ => panic!("inlet should be velocity-inlet"),
        }
        assert!(matches!(
            input.boundaries.get("outlet").unwrap(),
            Boundary::PressureOutlet { .. }
        ));
        assert!(matches!(
            input.boundaries.get("walls").unwrap(),
            Boundary::NoSlip
        ));
    }

    #[test]
    fn solver_kind_parses_with_or_without_prefix() {
        assert_eq!(
            SolverKind::from_solver_str("openfoam.simpleFoam"),
            Some(SolverKind::SimpleFoam)
        );
        assert_eq!(
            SolverKind::from_solver_str("openfoam.pimpleFoam"),
            Some(SolverKind::PimpleFoam)
        );
        assert_eq!(
            SolverKind::from_solver_str("openfoam.icoFoam"),
            Some(SolverKind::IcoFoam)
        );
        // Bare names (adapter-internal use) also resolve.
        assert_eq!(
            SolverKind::from_solver_str("pimpleFoam"),
            Some(SolverKind::PimpleFoam)
        );
        assert_eq!(
            SolverKind::from_solver_str("rhoSimpleFoam"),
            Some(SolverKind::RhoSimpleFoam)
        );
        assert_eq!(SolverKind::from_solver_str("openfoam.foamRun"), None);
    }

    #[test]
    fn parses_transient_pimple_foam_case() {
        let text = r#"
[case]
format  = "1.0"
name    = "demo-transient"
physics = "cfd"
solver  = "openfoam.pimpleFoam"
mesh    = "default"

[flow]
turbulence = "kOmegaSST"
schemes    = "linear-second-order"

[flow.fluid]
name = "air"
rho  = 1.225
nu   = 1.5e-5

[boundaries.inlet]
type     = "velocity-inlet"
velocity = [10.0, 0.0, 0.0]

[boundaries.outlet]
type     = "pressure-outlet"
pressure = 0.0

[boundaries.walls]
type = "no-slip"

[solve]
residual_target = 1e-5

[solve.transient]
end_time       = 2.0
delta_t        = 5e-4
write_interval = 0.05
"#;
        let cd: CaseDef = toml::from_str(text).unwrap();
        let input = SimpleFoamInput::from_case_def(&cd).expect("parse");
        assert_eq!(input.solver, SolverKind::PimpleFoam);
        match input.time {
            TimeMode::Transient {
                end_time,
                delta_t,
                write_interval,
            } => {
                assert!((end_time - 2.0).abs() < 1e-12);
                assert!((delta_t - 5e-4).abs() < 1e-12);
                assert!((write_interval - 0.05).abs() < 1e-12);
            }
            other => panic!("expected Transient, got {other:?}"),
        }
        // RANS still allowed under pimpleFoam.
        assert_eq!(input.turbulence, TurbulenceModel::KOmegaSST);
    }

    #[test]
    fn pimple_foam_defaults_transient_block_when_omitted() {
        let text = r#"
[case]
format  = "1.0"
name    = "demo"
physics = "cfd"
solver  = "openfoam.pimpleFoam"
mesh    = "default"

[flow]
turbulence = "laminar"
schemes    = "upwind-first-order"

[flow.fluid]
name = "water"
rho  = 1000.0
nu   = 1e-6

[boundaries.inlet]
type     = "velocity-inlet"
velocity = [1.0, 0.0, 0.0]

[boundaries.outlet]
type     = "pressure-outlet"
pressure = 0.0

[boundaries.walls]
type = "no-slip"

[solve]
residual_target = 1e-5
"#;
        let cd: CaseDef = toml::from_str(text).unwrap();
        let input = SimpleFoamInput::from_case_def(&cd).expect("parse");
        // Defaults: 1 s / 1 ms / 100 ms snapshots.
        assert!(matches!(input.time, TimeMode::Transient { .. }));
    }

    #[test]
    fn ico_foam_rejects_rans_turbulence() {
        let text = r#"
[case]
format  = "1.0"
name    = "demo-laminar"
physics = "cfd"
solver  = "openfoam.icoFoam"
mesh    = "default"

[flow]
turbulence = "kEpsilon"
schemes    = "upwind-first-order"

[flow.fluid]
name = "water"
rho  = 1000.0
nu   = 1e-6

[boundaries.inlet]
type     = "velocity-inlet"
velocity = [1.0, 0.0, 0.0]

[boundaries.walls]
type = "no-slip"

[solve]
residual_target = 1e-5

[solve.transient]
end_time       = 0.1
delta_t        = 1e-4
write_interval = 0.01
"#;
        let cd: CaseDef = toml::from_str(text).unwrap();
        let err = SimpleFoamInput::from_case_def(&cd).unwrap_err();
        match err {
            AdapterError::InvalidCase { reason, .. } => {
                assert!(
                    reason.contains("icoFoam") && reason.contains("laminar"),
                    "got: {reason}"
                );
            }
            other => panic!("wrong error: {other:?}"),
        }
    }

    #[test]
    fn simple_foam_with_transient_block_is_rejected() {
        let text = r#"
[case]
format  = "1.0"
name    = "wrong"
physics = "cfd"
solver  = "openfoam.simpleFoam"
mesh    = "default"

[flow]
turbulence = "laminar"
schemes    = "upwind-first-order"

[flow.fluid]
name = "air"
rho  = 1.225
nu   = 1.5e-5

[boundaries.x]
type = "empty"

[solve]
iterations = 100
residual_target = 1e-5

[solve.transient]
end_time       = 1.0
delta_t        = 1e-3
write_interval = 0.1
"#;
        let cd: CaseDef = toml::from_str(text).unwrap();
        let err = SimpleFoamInput::from_case_def(&cd).unwrap_err();
        match err {
            AdapterError::InvalidCase { reason, .. } => {
                assert!(reason.contains("steady-state"), "got: {reason}");
            }
            other => panic!("wrong error: {other:?}"),
        }
    }

    #[test]
    fn unknown_openfoam_solver_is_rejected() {
        let mut cd = sample_case();
        cd.case.solver = "openfoam.shockFoam".into();
        let err = SimpleFoamInput::from_case_def(&cd).unwrap_err();
        match err {
            AdapterError::InvalidCase { reason, .. } => {
                assert!(reason.contains("shockFoam"), "got: {reason}");
            }
            other => panic!("wrong error: {other:?}"),
        }
    }

    #[test]
    fn rejects_non_cfd_physics() {
        let mut cd = sample_case();
        cd.case.physics = "fea".into();
        assert!(matches!(
            SimpleFoamInput::from_case_def(&cd),
            Err(AdapterError::InvalidCase { .. })
        ));
    }

    // -----------------------------------------------------------------
    // Round-15 M3 RED→GREEN: OpenFOAM dict boundary-name injection
    // via the boundary key. The simpleFoam writer interpolates
    // `name` (from `for (name, boundary) in &input.boundaries`) into
    // ~10 dict files; pre-fix a hostile `inlet\n}\n\ninjectedPatch\n
    // {\n  type wall;` key would silently inject a sibling patch.
    // The validation lives in `from_case_def` so the failure surfaces
    // at case-load time with a clear field name, not at dict-write
    // time with a corrupted OpenFOAM run.
    // -----------------------------------------------------------------

    #[test]
    fn boundary_name_with_newline_is_rejected_at_load_time() {
        let text = r#"
[case]
format = "1.0"
name = "x"
physics = "cfd"
solver = "openfoam.simpleFoam"
mesh = "default"

[flow]
turbulence = "laminar"
schemes    = "upwind-first-order"

[flow.fluid]
name = "water"
rho  = 1000.0
nu   = 1e-6

[solve]
iterations      = 10
residual_target = 1e-3
"#;
        // Build the case manually — we need a raw key with embedded
        // newlines / braces, which is hard to express in TOML literal
        // syntax. Parse the base case, then inject the hostile key
        // via the BTreeMap, then call from_case_def.
        let mut cd: CaseDef = toml::from_str(text).unwrap();
        // Use a string-quoted TOML key to express the hostile bytes.
        let extra = "[boundaries.\"inlet\\n}\\n\\ninjectedPatch\\n{\\n  type wall;\"]\ntype = \"no-slip\"\n";
        let combined = format!("{text}\n{extra}");
        let parse_result: Result<CaseDef, _> = toml::from_str(&combined);
        // Pull the bad key directly from the combined parse if it
        // worked; otherwise fall back to manual injection.
        if let Ok(cd2) = parse_result {
            cd = cd2;
        } else {
            // TOML rejected the embedded newline in the key. Inject
            // directly into the BTreeMap representation by building
            // the section by hand.
            let mut sections = std::collections::BTreeMap::new();
            let mut bdy_table = toml::value::Table::new();
            let mut inlet_val = toml::value::Table::new();
            inlet_val.insert("type".into(), toml::Value::String("no-slip".into()));
            bdy_table.insert(
                "inlet\n}\n\ninjectedPatch\n{\n  type wall;".into(),
                toml::Value::Table(inlet_val),
            );
            sections.insert("boundaries".into(), toml::Value::Table(bdy_table));
            cd.sections = sections;
        }
        let err = SimpleFoamInput::from_case_def(&cd).expect_err("must reject hostile key");
        match err {
            AdapterError::InvalidCase { reason, .. } => {
                assert!(
                    reason.contains("boundaries.") || reason.contains("ASCII alphanumeric"),
                    "expected validation error mentioning the field; got: {reason}"
                );
            }
            other => panic!("wrong error: {other:?}"),
        }
    }

    #[test]
    fn boundary_name_with_brace_is_rejected_at_load_time() {
        // OpenFOAM dict uses `{` / `}` for block delimiters; a key
        // containing `}` would close the boundaryField block early.
        let mut cd = sample_case();
        let mut sections = cd.sections.clone();
        let bdy_table = sections.get_mut("boundaries").unwrap().as_table_mut().unwrap();
        let mut hostile = toml::value::Table::new();
        hostile.insert("type".into(), toml::Value::String("no-slip".into()));
        bdy_table.insert("evil}injected".into(), toml::Value::Table(hostile));
        cd.sections = sections;
        let err = SimpleFoamInput::from_case_def(&cd).expect_err("must reject");
        assert!(matches!(err, AdapterError::InvalidCase { .. }));
    }

    #[test]
    fn boundary_name_alphanumeric_is_accepted() {
        // Round-15 M3 must not regress legitimate boundary names:
        // alphanumeric + `.`, `-`, `_` are all valid OpenFOAM patch
        // names.
        let text = r#"
[case]
format = "1.0"
name = "demo"
physics = "cfd"
solver = "openfoam.simpleFoam"
mesh = "default"

[flow]
turbulence = "laminar"
schemes    = "upwind-first-order"

[flow.fluid]
name = "air"
rho  = 1.225
nu   = 1.5e-5

[boundaries.inlet-1]
type = "velocity-inlet"
velocity = [1.0, 0.0, 0.0]

[boundaries."outlet.upper"]
type = "pressure-outlet"
pressure = 0.0

[boundaries.walls_no_slip]
type = "no-slip"

[solve]
iterations = 100
residual_target = 1e-5
"#;
        let cd: CaseDef = toml::from_str(text).unwrap();
        let input = SimpleFoamInput::from_case_def(&cd).expect("alphanumeric+.-_ must pass");
        assert!(input.boundaries.contains_key("inlet-1"));
        assert!(input.boundaries.contains_key("outlet.upper"));
        assert!(input.boundaries.contains_key("walls_no_slip"));
    }

    #[test]
    fn rejects_unknown_boundary_type() {
        let text = r#"
[case]
format = "1.0"
name = "x"
physics = "cfd"
solver = "openfoam.simpleFoam"
mesh = "default"

[flow]
turbulence = "laminar"
schemes    = "upwind-first-order"

[flow.fluid]
name = "water"
rho  = 1000.0
nu   = 1e-6

[boundaries.mystery]
type = "teleporter"

[solve]
iterations      = 10
residual_target = 1e-3
"#;
        let cd: CaseDef = toml::from_str(text).unwrap();
        let err = SimpleFoamInput::from_case_def(&cd).unwrap_err();
        match err {
            AdapterError::InvalidCase { reason, .. } => {
                assert!(reason.contains("teleporter"), "got: {reason}");
            }
            other => panic!("wrong error: {other:?}"),
        }
    }
}
