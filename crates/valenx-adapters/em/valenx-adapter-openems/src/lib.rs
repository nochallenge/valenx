//! # valenx-adapter-openems
//!
//! Subprocess adapter for openEMS FDTD electromagnetics, driven via
//! a generated Octave script. **Phase 6 — live for rectangular FDTD
//! cases with Gauss / Sine excitation and Mur / PML / PEC
//! boundaries.**
//!
//! Antennas with explicit geometry, dispersive materials, and
//! S-parameter sweeps join by extending the `Domain` / `Material` /
//! `Probe` enums. The surrounding prepare / run / collect plumbing
//! stays as-is.

#![forbid(unsafe_code)]
#![allow(missing_docs)]

pub mod case_input;
pub mod octave_script;

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

use crate::case_input::EmInput;
use crate::octave_script::{SCRIPT_FILENAME, SIM_DIR};

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(OpenEmsAdapter::new())
}

pub struct OpenEmsAdapter;

impl OpenEmsAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for OpenEmsAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "openems";
/// Octave first; fall back to MATLAB if the user keeps a commercial
/// licence around. openEMS is happy with either.
const BINARIES: &[&str] = &["octave-cli", "octave", "matlab"];

impl Adapter for OpenEmsAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "openEMS",
            version_range: VersionRange {
                min_inclusive: Version::new(0, 0, 35),
                max_exclusive: Version::new(1, 0, 0),
            },
            physics: &[Physics::Em],
            license_mode: LicenseMode::Subprocess,
            tool_license: "GPL-3.0-or-later",
            docs_url: "https://docs.openems.de/",
            homepage_url: "https://openems.de/",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(BINARIES) {
            Some(binary_path) => {
                let found_version = valenx_core::adapter_helpers::detect_tool_version_semver(
                    &binary_path,
                    &["--version", "-v"],
                );
                Ok(ProbeReport {
                    ok: true,
                    found_version,
                    binary_path: Some(binary_path),
                    warnings: vec!["probe checks for Octave / MATLAB on PATH — openEMS + \
                     CSXCAD packages must also be installed in the same \
                     environment"
                        .into()],
                    required_env: Vec::new(),
                })
            }
            None => Err(AdapterError::ToolNotInstalled {
                name: INFO_ID,
                hint: "Octave + openEMS required; `apt install octave openems` \
                       on Debian-likes, or build from openems.de"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let (_header, input) = EmInput::from_case_dir(&case.path)?;

        fs::create_dir_all(workdir)?;

        let script_path = workdir.join(SCRIPT_FILENAME);
        octave_script::write_to_file(&input, &script_path)?;

        let binary_path = find_on_path(BINARIES).ok_or_else(|| AdapterError::ToolNotInstalled {
            name: INFO_ID,
            hint: "no Octave / MATLAB on PATH".into(),
        })?;
        let is_matlab = binary_path
            .file_name()
            .and_then(|s| s.to_str())
            .map(|s| s.starts_with("matlab"))
            .unwrap_or(false);

        let native_command: Vec<OsString> = if is_matlab {
            // MATLAB: `matlab -nodisplay -r "run('valenx_openems.m'); exit;"`.
            vec![
                binary_path.into_os_string(),
                OsString::from("-batch"),
                OsString::from(format!("run('{SCRIPT_FILENAME}')")),
            ]
        } else {
            // Octave: --no-gui --eval "valenx_openems"
            vec![
                binary_path.into_os_string(),
                OsString::from("--no-gui"),
                OsString::from("--quiet"),
                OsString::from(SCRIPT_FILENAME),
            ]
        };

        // 10 ns of 1 GHz excitation at 2 mm resolution on a 10-cm box
        // is ~30 s on a modern laptop; use 2 min as a generous UI
        // ceiling.
        let estimated_runtime = Some(Duration::from_secs(120));

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            estimated_runtime,
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting Octave/openEMS", |line| {
            let mut hint = subprocess::Hint::default();
            if let Some(pct) = openems_progress_hint(line) {
                hint.progress = Some((pct, line.to_string()));
            }
            if line.contains("error:") || line.contains("ERROR") {
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
        // case_path: the generated Octave script (canonical hashable
        // input — captures domain + materials + probes + excitation).
        // openEMS doesn't ingest a separate geometry/mesh file; the
        // FDTD grid + structures are defined inside the script via
        // CSXCAD calls, so mesh_path stays None.
        let case_path = job.workdir.join(SCRIPT_FILENAME);
        let prov = valenx_core::adapter_helpers::live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "openEMS",
            "unknown",
            &case_path,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);

        // The script writes everything into `workdir/sim/`.
        let sim_root = job.workdir.join(SIM_DIR);
        if sim_root.is_dir() {
            collect_into(&sim_root, &mut results);
            // For every `<probe>.csv` the post-processing block
            // produced, parse it into per-row ScalarRecords so the
            // report layer can chart probe waveforms without
            // re-reading the file. Failures (missing file, bad row)
            // are silent — the .csv stays listed as an artifact.
            load_probe_csvs_into_results(&sim_root, &mut results);
        }

        // Always attach the generated script for transparency.
        let script_path = job.workdir.join(SCRIPT_FILENAME);
        if script_path.is_file() {
            results.artifacts.push(Artifact {
                path: script_path,
                kind: ArtifactKind::Other,
                checksum: None,
                label: "openEMS Octave script (generated)".into(),
            });
        }

        results.artifacts.sort_by(|a, b| a.path.cmp(&b.path));
        Ok(results)
    }

    fn capabilities(&self) -> Capabilities {
        Capabilities {
            capabilities: vec![Capability::EmFdtdTimeDomain],
            ribbon_contributions: vec!["em.openems.fdtd"],
        }
    }
}

fn collect_into(dir: &Path, results: &mut Results) {
    let Ok(entries) = fs::read_dir(dir) else {
        return;
    };
    for entry in entries.flatten() {
        let path = entry.path();
        let Ok(ftype) = entry.file_type() else {
            continue;
        };
        if ftype.is_dir() {
            collect_into(&path, results);
            continue;
        }
        if !ftype.is_file() {
            continue;
        }
        let ext = path
            .extension()
            .and_then(|s| s.to_str())
            .map(|s| s.to_ascii_lowercase())
            .unwrap_or_default();
        let (kind, label) = match ext.as_str() {
            "h5" | "hdf5" => (ArtifactKind::Native, "openEMS HDF5 field dump".to_string()),
            "xml" => (ArtifactKind::Other, "openEMS FDTD XML".to_string()),
            "vtr" | "vts" | "vtu" => (
                ArtifactKind::VizData,
                "openEMS VTK field snapshot".to_string(),
            ),
            "dat" | "csv" => (
                ArtifactKind::Tabular,
                "openEMS probe time series".to_string(),
            ),
            _ => continue,
        };
        results.artifacts.push(Artifact {
            path,
            kind,
            checksum: None,
            label,
        });
    }
}

fn openems_progress_hint(line: &str) -> Option<f32> {
    if line.contains("Processing Properties") {
        Some(10.0)
    } else if line.contains("Initialization of FDTD") {
        Some(25.0)
    } else if line.contains("Timestep:") {
        // Keep mid-run ticks at 50 % so the bar doesn't oscillate.
        Some(50.0)
    } else if line.contains("Estimated runtime") {
        Some(60.0)
    } else if line.contains("Total simulation time") {
        Some(92.0)
    } else if line.contains("[valenx] openEMS finished") {
        Some(98.0)
    } else {
        None
    }
}

/// For every `<probe>.csv` file the post-processing block produced
/// inside the simulation directory, parse the `t_s,value` rows and
/// emit one [`valenx_fields::ScalarRecord`] per timestep into the
/// catalog. Probe name (the file stem) becomes the record name;
/// time goes into `TimeKey::Time { value, units: SECOND }`.
///
/// Failures (missing file, bad header, unparseable row) are skipped
/// silently — the .csv stays listed as an artifact for the user to
/// inspect by hand.
fn load_probe_csvs_into_results(sim_root: &Path, results: &mut Results) {
    use valenx_fields::units::SECOND;
    use valenx_fields::ScalarRecord;

    let Ok(entries) = fs::read_dir(sim_root) else {
        return;
    };
    for entry in entries.flatten() {
        let path = entry.path();
        if path.extension().and_then(|s| s.to_str()) != Some("csv") {
            continue;
        }
        let probe_name = match path
            .file_stem()
            .and_then(|s| s.to_str())
            .map(|s| s.to_string())
        {
            Some(n) => n,
            None => continue,
        };
        // Round-21 L1: bound the per-probe CSV read at
        // MAX_OPENEMS_CSV_BYTES (64 MiB). Pre-fix a poisoned workdir
        // with a multi-GB CSV would slurp before the line iterator
        // even parsed the header.
        let Ok(text) = valenx_core::io_caps::read_capped_to_string(
            &path,
            valenx_core::io_caps::MAX_OPENEMS_CSV_BYTES as usize,
        ) else {
            continue;
        };
        let mut lines = text.lines();
        // Skip header — should be `t_s,value`. If it's missing or
        // shaped differently we still try to parse the body but
        // silently bail on rows that don't fit `<float>,<float>`.
        let header = lines.next().unwrap_or("").trim();
        if !header.starts_with("t_s") {
            // Not our format; skip.
            continue;
        }
        for line in lines {
            let trimmed = line.trim();
            if trimmed.is_empty() {
                continue;
            }
            let mut iter = trimmed.split(',');
            let (Some(t_str), Some(v_str)) = (iter.next(), iter.next()) else {
                continue;
            };
            let (Ok(t_s), Ok(value)) = (t_str.trim().parse::<f64>(), v_str.trim().parse::<f64>())
            else {
                continue;
            };
            results.scalars.insert(ScalarRecord {
                name: probe_name.clone(),
                value,
                units: valenx_fields::units::DIMENSIONLESS,
                time: valenx_fields::TimeKey::Time {
                    value: t_s,
                    units: SECOND,
                },
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
    fn info_is_em_domain() {
        let info = OpenEmsAdapter::new().info();
        assert_eq!(info.id, "openems");
        assert_eq!(info.physics, &[Physics::Em]);
    }

    #[test]
    fn collect_uses_live_provenance_with_real_case_hash() {
        let workdir = tempdir("openems-collect-prov");
        let script_path = workdir.join(SCRIPT_FILENAME);
        let script_bytes = b"% openEMS placeholder script\n";
        std::fs::write(&script_path, script_bytes).expect("write script");

        let job = PreparedJob {
            workdir: workdir.clone(),
            native_command: Vec::new(),
            environment: Vec::new(),
            estimated_runtime: None,
            kill_on_drop: false,
        };
        let results = OpenEmsAdapter::new().collect(&job).expect("collect");
        let prov = &results.provenance;

        assert_eq!(prov.adapter, INFO_ID);
        assert!(!prov.adapter_version.is_empty());
        assert_eq!(prov.tool, "openEMS");
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
    fn loads_probe_csvs_into_scalar_catalog() {
        use valenx_fields::provenance::Sha256Hex;
        let prov = valenx_fields::Provenance {
            adapter: "openems".into(),
            adapter_version: "0".into(),
            tool: "openEMS".into(),
            tool_version: "0".into(),
            case_hash: Sha256Hex::new(""),
            mesh_hash: Sha256Hex::new(""),
            input_hash: Sha256Hex::new(""),
            tools_lock_hash: Sha256Hex::new(""),
            run_id: "00000000-0000-0000-0000-000000000000".into(),
            wall_time_seconds: 0.0,
            completed_at: "1970-01-01T00:00:00Z".into(),
            ancestors: Vec::new(),
        };
        let mut results = Results::empty("openems-test", prov);

        // Set up a fake sim/ dir with a probe CSV.
        let tmp = std::env::temp_dir().join(format!(
            "valenx-openems-probe-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&tmp).unwrap();
        std::fs::write(
            tmp.join("center.csv"),
            "t_s,value\n0,0\n1e-9,0.5\n2e-9,1.0\n",
        )
        .unwrap();
        // A non-probe CSV should be ignored (header doesn't start
        // with t_s) — we leave a malformed file in to verify.
        std::fs::write(tmp.join("notes.csv"), "x,y,z\n1,2,3\n").unwrap();

        super::load_probe_csvs_into_results(&tmp, &mut results);

        // Only the center probe contributes. 3 rows × 1 record each.
        assert_eq!(results.scalars.len(), 3);
        let at_t1 = results
            .scalars
            .all("center")
            .iter()
            .find(|r| {
                matches!(
                    r.time,
                    valenx_fields::TimeKey::Time { value, .. } if (value - 1e-9).abs() < 1e-18
                )
            })
            .expect("center at t=1 ns");
        assert!((at_t1.value - 0.5).abs() < 1e-12);

        let _ = std::fs::remove_dir_all(&tmp);
    }

    #[test]
    fn progress_hints_monotonic() {
        let pts = [
            openems_progress_hint("Processing Properties"),
            openems_progress_hint("Initialization of FDTD engine..."),
            openems_progress_hint("Timestep: 500"),
            openems_progress_hint("Estimated runtime: 12 seconds"),
            openems_progress_hint("Total simulation time: 5.000000e+01 s"),
            openems_progress_hint("[valenx] openEMS finished; sim path = sim"),
        ];
        let mut last = 0.0f32;
        for (i, p) in pts.iter().enumerate() {
            let v = p.expect("banner");
            assert!(v >= last, "step {i}: {last} -> {v}");
            last = v;
        }
    }
}
