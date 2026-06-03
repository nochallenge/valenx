//! # valenx-adapter-cantera
//!
//! Adapter for Cantera — reaction kinetics + thermodynamics +
//! transport.
//!
//! **Phase 4 — live for equilibrium.** `prepare()` writes a
//! deterministic Python script that drives Cantera's Python API
//! for the selected analysis. `run()` spawns `python
//! valenx_cantera.py` via the shared subprocess runner. `collect()`
//! parses `summary.json` (initial + final thermodynamic states +
//! filtered mole-fraction map) and attaches it alongside the script.
//!
//! Today: TP / HP / UV equilibrium. Tomorrow: 0-D reactor networks
//! and 1-D freely-propagating flames (each extends the
//! `Analysis` enum + adds a Python code path; the surrounding
//! prepare / run / collect plumbing stays as is).

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
    adapter_helpers::find_on_path, error::RunPhase, subprocess, Adapter, AdapterError, AdapterInfo,
    Capabilities, Capability, Case, LicenseMode, Physics, PreparedJob, ProbeReport, RunContext,
    RunReport, VersionRange,
};
use valenx_fields::{
    artifact::{Artifact, ArtifactKind},
    Results,
};

use crate::case_input::ChemistryInput;
use crate::python_script::{SCRIPT_FILENAME, SUMMARY_FILENAME};

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(CanteraAdapter::new())
}

pub struct CanteraAdapter;

impl CanteraAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for CanteraAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "cantera";
/// The Python interpreter we invoke. Order matters: `python` on
/// Windows usually resolves to the latest 3.x; on Linux it's often
/// legacy-2. `python3` is the safer bet where both exist.
const PYTHON_BINARIES: &[&str] = &["python3", "python"];

impl Adapter for CanteraAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "Cantera",
            version_range: VersionRange {
                min_inclusive: Version::new(3, 0, 0),
                max_exclusive: Version::new(4, 0, 0),
            },
            physics: &[Physics::Chemistry],
            license_mode: LicenseMode::Subprocess,
            tool_license: "BSD-3-Clause",
            docs_url: "https://cantera.org/documentation/",
            homepage_url: "https://cantera.org/",
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
                    warnings: vec!["probe checks for `python` on PATH — Cantera itself \
                     (`pip install cantera`) must also be importable for \
                     runs to succeed"
                        .into()],
                    required_env: Vec::new(),
                })
            }
            None => Err(AdapterError::ToolNotInstalled {
                name: INFO_ID,
                hint: "Python 3.9+ with Cantera installed; \
                       `pip install cantera` after ensuring python3 is on PATH"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let (_header, input) = ChemistryInput::from_case_dir(&case.path)?;

        fs::create_dir_all(workdir)?;

        // If the mechanism is an external file, stage it alongside
        // the script so the relative path in the generated Python
        // resolves regardless of cwd.
        // Round-9 hardening: `Mechanism::External(path)` is user data
        // and gets copied/read; wrap relative paths with `confined_join`
        // so a hostile case can't aim it at `../../etc/passwd`.
        if let case_input::Mechanism::External(path) = &input.mechanism {
            let source = if path.is_absolute() {
                path.clone()
            } else {
                valenx_core::adapter_helpers::confined_join(&case.path, path)?
            };
            if !source.is_file() {
                return Err(AdapterError::InvalidCase {
                    case_path: case.path.join("case.toml"),
                    reason: format!(
                        "[chemistry] mechanism file {} not found (resolved {})",
                        path.display(),
                        source.display()
                    ),
                });
            }
            let file_name = path.file_name().ok_or_else(|| AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!("mechanism path `{}` has no filename", path.display()),
            })?;
            let dest = workdir.join(file_name);
            if source != dest {
                fs::copy(&source, &dest)?;
            }
        }

        // Rewrite the mechanism path so the emitted Python refers to
        // the workdir-local copy.
        let write_input = maybe_rewrite_external_mechanism(input);

        let script_path = workdir.join(SCRIPT_FILENAME);
        python_script::write_to_file(&write_input, &script_path)?;

        let binary_path =
            find_on_path(PYTHON_BINARIES).ok_or_else(|| AdapterError::ToolNotInstalled {
                name: INFO_ID,
                hint: "no python / python3 on PATH — install Python 3.9+ first".into(),
            })?;

        let native_command: Vec<OsString> = vec![
            binary_path.into_os_string(),
            OsString::from(SCRIPT_FILENAME),
        ];

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            // Equilibrium on GRI-3.0 is sub-second; give ourselves a
            // generous ceiling for large mechanisms.
            estimated_runtime: Some(Duration::from_secs(30)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting Cantera", |line| {
            let mut hint = subprocess::Hint::default();
            if line.contains("[valenx] equilibrium reached") {
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
        // Real provenance: hash the generated Python script (the
        // canonical "this case is configured this way" input —
        // captures mechanism + initial state + analysis kind) and
        // any external mechanism YAML the user supplied.
        let case_path = job.workdir.join(SCRIPT_FILENAME);
        // The mechanism YAML lives next to the case input; if it
        // got staged into the workdir it'll be the only .yaml file
        // there (the script itself is .py).
        let mesh_path = job.workdir.join("mechanism.yaml");
        let prov = valenx_core::adapter_helpers::live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "Cantera",
            "unknown",
            &case_path,
            if mesh_path.exists() {
                Some(mesh_path.as_path())
            } else {
                None
            },
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);

        let summary_path = job.workdir.join(SUMMARY_FILENAME);
        if summary_path.is_file() {
            match summary_parser::parse_file(&summary_path) {
                Ok(summary) => {
                    let desc = match &summary.final_ {
                        Some(frame) => format!(
                            "Cantera {} · T_final={:.1} K · {} species kept",
                            summary.analysis,
                            frame.temperature_k,
                            summary.species_count_kept.unwrap_or(0),
                        ),
                        None => format!("Cantera {}", summary.analysis),
                    };
                    results.meta.description = Some(desc);
                    // Populate the scalar catalog with the summary
                    // values so the report layer can read them
                    // without re-parsing the JSON. Cantera is 0D
                    // chemistry so there's no mesh and no Field
                    // catalog to populate; ScalarRecord is the
                    // right type for thermo + species data.
                    populate_scalars_from_summary(&mut results, &summary);
                    results.artifacts.push(Artifact {
                        path: summary_path.clone(),
                        kind: ArtifactKind::Tabular,
                        checksum: None,
                        label: format!(
                            "Cantera summary ({} species)",
                            summary.mole_fractions.len()
                        ),
                    });
                }
                Err(e) => {
                    tracing::warn!(target: "valenx-cantera", ?e, "summary parse failed");
                    results.artifacts.push(Artifact {
                        path: summary_path,
                        kind: ArtifactKind::Tabular,
                        checksum: None,
                        label: format!("Cantera summary (parse error: {e})"),
                    });
                }
            }
        }

        let script_path = job.workdir.join(SCRIPT_FILENAME);
        if script_path.is_file() {
            results.artifacts.push(Artifact {
                path: script_path,
                kind: ArtifactKind::Other,
                checksum: None,
                label: "Cantera script (generated)".into(),
            });
        }

        results.artifacts.sort_by(|a, b| a.path.cmp(&b.path));
        Ok(results)
    }

    fn capabilities(&self) -> Capabilities {
        Capabilities {
            capabilities: vec![
                Capability::ChemKinetics,
                Capability::ChemEquilibrium,
                Capability::ChemTransport,
                Capability::ChemCombustion,
            ],
            ribbon_contributions: vec![
                "chem.cantera.equilibrium",
                "chem.cantera.reactor",
                "chem.cantera.flame",
            ],
        }
    }
}

/// If the mechanism is `External(path)`, rewrite it to reference the
/// filename only — the adapter copies the external file alongside
/// the script during prepare(). The workdir-local copy lets the
/// Python script resolve regardless of the user's cwd.
/// Materialise the parsed [`summary_parser::CanteraSummary`] into the
/// canonical `Results.scalars` catalog. Cantera is zero-dimensional
/// chemistry — there's no mesh and no Field catalog to populate, so
/// the scalar-record path is the right home for thermo + species data.
///
/// Emits one record per:
/// - initial-state thermo (T, P, H, S, ρ) tagged with `_initial` suffix
/// - final-state thermo same shape, `_final` suffix
/// - mean molecular weight (no suffix — single value)
/// - every species mole fraction (`X_<species>`)
fn populate_scalars_from_summary(results: &mut Results, summary: &summary_parser::CanteraSummary) {
    use valenx_fields::units::*;
    use valenx_fields::ScalarRecord;

    // Define units inline here — kelvin / pascal / J/kg / kg/m^3.
    // valenx_fields' built-in constants cover the SI base units; for
    // "joules per kilogram" etc. we synthesise via Units::new with
    // the dimensional vector [L, M, T, I, Θ, N, J].
    let kelvin = KELVIN;
    let pascal = Units::new([-1, 1, -2, 0, 0, 0, 0], 1.0, Some("Pa"));
    let j_per_kg = Units::new([2, 0, -2, 0, 0, 0, 0], 1.0, Some("J/kg"));
    let j_per_kg_k = Units::new([2, 0, -2, 0, -1, 0, 0], 1.0, Some("J/(kg·K)"));
    let kg_per_m3 = Units::new([-3, 1, 0, 0, 0, 0, 0], 1.0, Some("kg/m³"));
    let g_per_mol = Units::new([0, 1, 0, 0, 0, -1, 0], 1e-3, Some("g/mol"));

    if let Some(initial) = summary.initial.as_ref() {
        results.scalars.insert(ScalarRecord::extracted(
            "T_initial",
            initial.temperature_k,
            kelvin,
        ));
        results.scalars.insert(ScalarRecord::extracted(
            "P_initial",
            initial.pressure_pa,
            pascal,
        ));
        results.scalars.insert(ScalarRecord::extracted(
            "H_initial",
            initial.enthalpy_mass,
            j_per_kg,
        ));
        results.scalars.insert(ScalarRecord::extracted(
            "S_initial",
            initial.entropy_mass,
            j_per_kg_k,
        ));
        results.scalars.insert(ScalarRecord::extracted(
            "rho_initial",
            initial.density,
            kg_per_m3,
        ));
    }
    if let Some(final_) = summary.final_.as_ref() {
        results.scalars.insert(ScalarRecord::extracted(
            "T_final",
            final_.temperature_k,
            kelvin,
        ));
        results.scalars.insert(ScalarRecord::extracted(
            "P_final",
            final_.pressure_pa,
            pascal,
        ));
        results.scalars.insert(ScalarRecord::extracted(
            "H_final",
            final_.enthalpy_mass,
            j_per_kg,
        ));
        results.scalars.insert(ScalarRecord::extracted(
            "S_final",
            final_.entropy_mass,
            j_per_kg_k,
        ));
        results.scalars.insert(ScalarRecord::extracted(
            "rho_final",
            final_.density,
            kg_per_m3,
        ));
    }
    if let Some(mw) = summary.mean_molecular_weight {
        results.scalars.insert(ScalarRecord::extracted(
            "mean_molecular_weight",
            mw,
            g_per_mol,
        ));
    }
    // One scalar per species mole fraction. Names prefixed `X_` so
    // they don't collide with thermo names if a future schema adds
    // an `H` or `S` species.
    for (species, fraction) in &summary.mole_fractions {
        results.scalars.insert(ScalarRecord::extracted(
            &format!("X_{species}"),
            *fraction,
            DIMENSIONLESS,
        ));
    }
}

fn maybe_rewrite_external_mechanism(input: ChemistryInput) -> ChemistryInput {
    use case_input::Mechanism;
    match input.mechanism {
        Mechanism::External(path) => ChemistryInput {
            mechanism: Mechanism::External(std::path::PathBuf::from(
                path.file_name()
                    .map(|s| s.to_string_lossy().into_owned())
                    .unwrap_or_default(),
            )),
            ..input
        },
        Mechanism::Bundled(_) => input,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn info_is_chemistry_domain() {
        let info = CanteraAdapter::new().info();
        assert_eq!(info.id, "cantera");
        assert_eq!(info.physics, &[Physics::Chemistry]);
        assert_eq!(info.tool_license, "BSD-3-Clause");
    }

    #[test]
    fn populate_scalars_extracts_thermo_and_species() {
        use std::collections::BTreeMap;
        use summary_parser::{CanteraSummary, ThermoFrame};
        use valenx_fields::provenance::Sha256Hex;
        let mut species = BTreeMap::new();
        species.insert("N2".to_string(), 0.78);
        species.insert("O2".to_string(), 0.21);
        species.insert("Ar".to_string(), 0.01);
        let summary = CanteraSummary {
            valenx_adapter: "cantera".into(),
            analysis: "equilibrium-tp".into(),
            mechanism: "gri30.yaml".into(),
            initial: Some(ThermoFrame {
                temperature_k: 300.0,
                pressure_pa: 101325.0,
                enthalpy_mass: 1000.0,
                entropy_mass: 6890.0,
                density: 1.225,
            }),
            final_: Some(ThermoFrame {
                temperature_k: 2400.0,
                pressure_pa: 101325.0,
                enthalpy_mass: 1500000.0,
                entropy_mass: 8200.0,
                density: 0.151,
            }),
            mole_fractions: species,
            mean_molecular_weight: Some(28.97),
            species_count_kept: Some(3),
            species_count_total: Some(53),
        };
        let prov = valenx_fields::Provenance {
            adapter: "cantera".into(),
            adapter_version: "0".into(),
            tool: "Cantera".into(),
            tool_version: "3.0".into(),
            case_hash: Sha256Hex::new(""),
            mesh_hash: Sha256Hex::new(""),
            input_hash: Sha256Hex::new(""),
            tools_lock_hash: Sha256Hex::new(""),
            run_id: "00000000-0000-0000-0000-000000000000".into(),
            wall_time_seconds: 0.0,
            completed_at: "1970-01-01T00:00:00Z".into(),
            ancestors: Vec::new(),
        };
        let mut results = Results::empty("cantera-test", prov);
        super::populate_scalars_from_summary(&mut results, &summary);

        // 5 initial + 5 final + 1 mean MW + 3 species = 14 records.
        assert_eq!(results.scalars.len(), 14);
        // Spot-check key entries.
        let t_final = results.scalars.get("T_final").expect("T_final");
        assert!((t_final.value - 2400.0).abs() < 1e-9);
        let x_n2 = results.scalars.get("X_N2").expect("X_N2");
        assert!((x_n2.value - 0.78).abs() < 1e-9);
    }

    #[test]
    fn rewrite_external_strips_directory() {
        let input = ChemistryInput {
            mechanism: case_input::Mechanism::External(std::path::PathBuf::from(
                "sub/dir/custom.yaml",
            )),
            analysis: case_input::Analysis::EquilibriumTP,
            initial: case_input::ThermoState {
                temperature_k: 300.0,
                pressure_pa: 101325.0,
                composition: "N2:1".into(),
            },
            reactor: None,
        };
        let rewritten = maybe_rewrite_external_mechanism(input);
        match rewritten.mechanism {
            case_input::Mechanism::External(path) => {
                assert_eq!(path, std::path::PathBuf::from("custom.yaml"));
            }
            other => panic!("wrong variant: {other:?}"),
        }
    }

    #[test]
    fn rewrite_bundled_is_noop() {
        let input = ChemistryInput {
            mechanism: case_input::Mechanism::Bundled("gri30.yaml".into()),
            analysis: case_input::Analysis::EquilibriumTP,
            initial: case_input::ThermoState {
                temperature_k: 300.0,
                pressure_pa: 101325.0,
                composition: "N2:1".into(),
            },
            reactor: None,
        };
        let rewritten = maybe_rewrite_external_mechanism(input);
        assert!(matches!(
            rewritten.mechanism,
            case_input::Mechanism::Bundled(_)
        ));
    }

    /// Round-9 RED→GREEN: `Mechanism::External(path)` used to be
    /// joined with bare `case.path.join`. Wrap relative paths with
    /// `confined_join` so a hostile case can't ask Cantera to read
    /// `../../etc/passwd`.
    #[test]
    fn prepare_rejects_external_mechanism_traversing_outside_case_dir() {
        let case_dir = std::env::temp_dir().join(format!(
            "valenx-cantera-mech-trav-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&case_dir).unwrap();
        std::fs::write(
            case_dir.join("case.toml"),
            r#"[case]
format  = "1.0"
name    = "trav"
physics = "chemistry"
solver  = "cantera.equilibrium"
mesh    = "(none)"

[chemistry]
mechanism = "../../etc/passwd"
analysis  = "equilibrium-tp"

[chemistry.initial]
T           = 300.0
P           = 101325.0
composition = "N2:1"
"#,
        )
        .unwrap();
        let workdir = case_dir.join("workdir");
        let case = Case {
            id: "cantera-mech-trav".into(),
            path: case_dir.clone(),
        };
        let err = CanteraAdapter::new()
            .prepare(&case, &workdir)
            .expect_err("must reject ../../etc/passwd mechanism");
        let msg = format!("{err}");
        assert!(
            msg.contains("..") || msg.contains("stay within") || msg.contains("escape"),
            "expected confined_join rejection, got: {msg}"
        );
        let _ = std::fs::remove_dir_all(&case_dir);
    }
}
