//! # valenx-adapter-pybamm
//!
//! Adapter for PyBaMM — Python-native battery modelling. **Phase 7
//! — live for single-protocol discharge / charge on DFN / SPM /
//! SPMe with built-in parameter sets.** Multi-cycle aging studies
//! and drive-cycle files extend the `Protocol` enum.

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

use crate::case_input::BatteryInput;
use crate::python_script::{SCRIPT_FILENAME, SUMMARY_FILENAME, TIMESERIES_FILENAME};

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(PyBammAdapter::new())
}

pub struct PyBammAdapter;

impl PyBammAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for PyBammAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "pybamm";
const PYTHON_BINARIES: &[&str] = &["python3", "python"];

impl Adapter for PyBammAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "PyBaMM",
            version_range: VersionRange {
                min_inclusive: Version::new(23, 0, 0),
                max_exclusive: Version::new(26, 0, 0),
            },
            physics: &[Physics::Battery],
            license_mode: LicenseMode::Subprocess,
            tool_license: "BSD-3-Clause",
            docs_url: "https://docs.pybamm.org/",
            homepage_url: "https://www.pybamm.org/",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(PYTHON_BINARIES) {
            Some(binary_path) => {
                // We probe for the Python interpreter (PyBaMM is a
                // Python lib loaded at run time). Reporting Python's
                // version is still useful — operators can confirm
                // PyBaMM's `requires-python` is satisfied without
                // dropping into a shell.
                let found_version = valenx_core::adapter_helpers::detect_tool_version_semver(
                    &binary_path,
                    &["--version", "-V"],
                );
                Ok(ProbeReport {
                    ok: true,
                    found_version,
                    binary_path: Some(binary_path),
                    warnings: vec!["probe checks for Python; PyBaMM itself (`pip install \
                     pybamm`) must also be importable for runs to succeed"
                        .into()],
                    required_env: Vec::new(),
                })
            }
            None => Err(AdapterError::ToolNotInstalled {
                name: INFO_ID,
                hint: "PyBaMM requires Python 3.9+; install via \
                       `pip install pybamm`"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let (_header, input) = BatteryInput::from_case_dir(&case.path)?;

        fs::create_dir_all(workdir)?;

        let script_path = workdir.join(SCRIPT_FILENAME);
        python_script::write_to_file(&input, &script_path)?;

        let binary_path =
            find_on_path(PYTHON_BINARIES).ok_or_else(|| AdapterError::ToolNotInstalled {
                name: INFO_ID,
                hint: "no python / python3 on PATH".into(),
            })?;

        let native_command: Vec<OsString> = vec![
            binary_path.into_os_string(),
            OsString::from(SCRIPT_FILENAME),
        ];

        // DFN single-discharge on Chen2020 takes a few seconds on a
        // decent laptop; SPM ~1 s. Budget 5 min for aging / long
        // horizons.
        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            estimated_runtime: Some(Duration::from_secs(300)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting PyBaMM", |line| {
            let mut hint = subprocess::Hint::default();
            if line.contains("[valenx] pybamm done") {
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
        // case_path: the generated Python script (canonical hashable
        // input — captures model + parameters + protocol). PyBaMM
        // doesn't ingest a separate mesh/topology file (the cell
        // geometry is parameterised inside the script), so
        // mesh_path stays None.
        let case_path = job.workdir.join(SCRIPT_FILENAME);
        let prov = valenx_core::adapter_helpers::live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "PyBaMM",
            "unknown",
            &case_path,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);

        let summary_path = job.workdir.join(SUMMARY_FILENAME);
        if summary_path.is_file() {
            if let Ok(summary) = summary_parser::parse_file(&summary_path) {
                results.meta.description = Some(format!(
                    "PyBaMM {} · {} samples · {:.3} V → {:.3} V over {:.1} s",
                    summary.model,
                    summary.samples,
                    summary.voltage_start_v,
                    summary.voltage_end_v,
                    summary.duration_s,
                ));
            }
            results.artifacts.push(Artifact {
                path: summary_path,
                kind: ArtifactKind::Tabular,
                checksum: None,
                label: "PyBaMM summary".into(),
            });
        }
        let ts_path = job.workdir.join(TIMESERIES_FILENAME);
        if ts_path.is_file() {
            // Parse the time-series CSV into ScalarRecord entries
            // keyed by physical time (TimeKey::Time). Lets the report
            // layer chart voltage(t) + current(t) without re-reading
            // the file. Failures are skipped silently.
            //
            // Round-23 named finding: bound the read at
            // MAX_PYBAMM_TIMESERIES_BYTES (256 MiB) — pre-fix a
            // poisoned or runaway discharge.csv would slurp into
            // memory before the line iterator parsed the header.
            if let Ok(text) = valenx_core::io_caps::read_capped_to_string(
                &ts_path,
                valenx_core::io_caps::MAX_PYBAMM_TIMESERIES_BYTES as usize,
            ) {
                load_pybamm_timeseries_into_results(&mut results, &text);
            }
            results.artifacts.push(Artifact {
                path: ts_path,
                kind: ArtifactKind::Tabular,
                checksum: None,
                label: "PyBaMM discharge time series".into(),
            });
        }
        let script_path = job.workdir.join(SCRIPT_FILENAME);
        if script_path.is_file() {
            results.artifacts.push(Artifact {
                path: script_path,
                kind: ArtifactKind::Other,
                checksum: None,
                label: "PyBaMM script (generated)".into(),
            });
        }

        results.artifacts.sort_by(|a, b| a.path.cmp(&b.path));
        Ok(results)
    }

    fn capabilities(&self) -> Capabilities {
        Capabilities {
            capabilities: vec![Capability::BatteryDfn, Capability::BatterySpm],
            ribbon_contributions: vec![
                "battery.pybamm.dfn",
                "battery.pybamm.spm",
                "battery.pybamm.spme",
            ],
        }
    }
}

/// Parse PyBaMM's `timeseries.csv` (`t_s,voltage_v,current_a` shape)
/// and insert one `ScalarRecord` per (column, timestep) into the
/// `Results.scalars` catalog. The first column is treated as time
/// (in seconds); subsequent columns become per-row records keyed
/// to that time.
///
/// Multi-column CSVs from future PyBaMM workflows (capacity, SOC,
/// temperature) work the same way without code changes.
fn load_pybamm_timeseries_into_results(results: &mut valenx_fields::Results, csv_text: &str) {
    use valenx_fields::units::{Units, SECOND};
    use valenx_fields::ScalarRecord;

    let mut lines = csv_text.lines();
    let header_line = match lines.find(|l| !l.trim().is_empty()) {
        Some(h) => h,
        None => return,
    };
    let columns: Vec<String> = header_line
        .split(',')
        .map(|s| s.trim().to_string())
        .collect();
    if columns.len() < 2 {
        return;
    }
    // Lookup table for a few canonical PyBaMM column names. Anything
    // unknown lands as DIMENSIONLESS — better than guessing wrong.
    let units_for = |name: &str| -> Units {
        match name {
            "t_s" | "time_s" => SECOND,
            "voltage_v" | "v" => Units::new([2, 1, -3, -1, 0, 0, 0], 1.0, Some("V")),
            "current_a" | "i" => Units::new([0, 0, 0, 1, 0, 0, 0], 1.0, Some("A")),
            "capacity_ah" => Units::new([0, 0, 1, 1, 0, 0, 0], 3600.0, Some("A·h")),
            "temperature_k" | "t_k" => Units::new([0, 0, 0, 0, 1, 0, 0], 1.0, Some("K")),
            _ => valenx_fields::units::DIMENSIONLESS,
        }
    };

    for line in lines {
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        let toks: Vec<&str> = trimmed.split(',').map(|s| s.trim()).collect();
        if toks.len() != columns.len() {
            continue;
        }
        // First column is time. Skip rows where it doesn't parse.
        let time_s: f64 = match toks[0].parse() {
            Ok(v) => v,
            Err(_) => continue,
        };
        let timekey = valenx_fields::TimeKey::Time {
            value: time_s,
            units: SECOND,
        };
        for (col, tok) in columns.iter().zip(toks.iter()).skip(1) {
            let value: f64 = match tok.parse() {
                Ok(v) => v,
                Err(_) => continue,
            };
            results.scalars.insert(ScalarRecord {
                name: col.clone(),
                value,
                units: units_for(col),
                time: timekey,
                source: valenx_fields::scalar::ScalarSource::Extracted,
                description: None,
            });
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_test_utils::tempdir;

    #[test]
    fn info_is_battery_domain() {
        let info = PyBammAdapter::new().info();
        assert_eq!(info.id, "pybamm");
        assert_eq!(info.physics, &[Physics::Battery]);
        assert_eq!(info.tool_license, "BSD-3-Clause");
    }

    #[test]
    fn collect_uses_live_provenance_with_real_case_hash() {
        let workdir = tempdir("pybamm-collect-prov");
        let script_path = workdir.join(SCRIPT_FILENAME);
        let script_bytes = b"import pybamm\n# trivial\n";
        std::fs::write(&script_path, script_bytes).expect("write script");

        let job = PreparedJob {
            workdir: workdir.clone(),
            native_command: Vec::new(),
            environment: Vec::new(),
            estimated_runtime: None,
            kill_on_drop: false,
        };
        let results = PyBammAdapter::new().collect(&job).expect("collect");
        let prov = &results.provenance;

        assert_eq!(prov.adapter, INFO_ID);
        assert!(!prov.adapter_version.is_empty());
        assert_eq!(prov.tool, "PyBaMM");
        assert!(!prov.run_id.is_empty(), "run_id empty — stub still wired?");
        assert_eq!(
            prov.case_hash,
            valenx_core::adapter_helpers::sha256_hex_file(&script_path)
        );

        cleanup_lp(&workdir);
    }

    fn cleanup_lp(d: &std::path::Path) {
        let _ = std::fs::remove_dir_all(d);
    }

    #[test]
    fn loads_pybamm_csv_into_scalar_catalog() {
        use valenx_fields::provenance::Sha256Hex;
        let prov = valenx_fields::Provenance {
            adapter: "pybamm".into(),
            adapter_version: "0".into(),
            tool: "PyBaMM".into(),
            tool_version: "24".into(),
            case_hash: Sha256Hex::new(""),
            mesh_hash: Sha256Hex::new(""),
            input_hash: Sha256Hex::new(""),
            tools_lock_hash: Sha256Hex::new(""),
            run_id: "00000000-0000-0000-0000-000000000000".into(),
            wall_time_seconds: 0.0,
            completed_at: "1970-01-01T00:00:00Z".into(),
            ancestors: Vec::new(),
        };
        let mut results = Results::empty("pybamm-test", prov);

        let csv = "t_s,voltage_v,current_a\n0,4.2,-5\n10,4.18,-5\n20,4.15,-5\n";
        super::load_pybamm_timeseries_into_results(&mut results, csv);

        // 3 rows × 2 non-time columns = 6 scalar records.
        assert_eq!(results.scalars.len(), 6);
        // Each voltage/current entry is keyed by its physical time.
        let voltage_at_t0 = results
            .scalars
            .all("voltage_v")
            .iter()
            .find(|r| {
                matches!(
                    r.time,
                    valenx_fields::TimeKey::Time { value, .. } if (value - 0.0).abs() < 1e-12
                )
            })
            .expect("voltage at t=0");
        assert!((voltage_at_t0.value - 4.2).abs() < 1e-12);
    }

    #[test]
    fn pybamm_csv_skips_malformed_rows() {
        use valenx_fields::provenance::Sha256Hex;
        let prov = valenx_fields::Provenance {
            adapter: "pybamm".into(),
            adapter_version: "0".into(),
            tool: "PyBaMM".into(),
            tool_version: "24".into(),
            case_hash: Sha256Hex::new(""),
            mesh_hash: Sha256Hex::new(""),
            input_hash: Sha256Hex::new(""),
            tools_lock_hash: Sha256Hex::new(""),
            run_id: "00000000-0000-0000-0000-000000000000".into(),
            wall_time_seconds: 0.0,
            completed_at: "1970-01-01T00:00:00Z".into(),
            ancestors: Vec::new(),
        };
        let mut results = Results::empty("pybamm-test", prov);
        // Mix valid + malformed rows. Parser should keep only the
        // valid ones rather than fail outright.
        let csv = "t_s,voltage_v\n0,4.2\nbad-time,9.9\n10,4.1\n";
        super::load_pybamm_timeseries_into_results(&mut results, csv);
        assert_eq!(results.scalars.len(), 2);
    }
}
