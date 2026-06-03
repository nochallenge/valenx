//! # valenx-adapter-cello
//!
//! Adapter for [Cello v2](https://github.com/CIDARLAB/Cello-v2) —
//! the canonical genetic-circuit DNA compiler from CIDAR. Cello
//! consumes a Verilog netlist describing the desired logic
//! function plus a triplet of JSON constraint files (a user
//! constraint file pinning the chassis / library, an input sensor
//! file pinning the input promoters, an output device file pinning
//! the reporter), and emits a fully assembled DNA construct that
//! implements the logic in a living cell. The compiler runs a
//! simulated-annealing optimization over the gate-assignment
//! problem and outputs a Graphviz `.dot` netlist, a circuit
//! diagram PNG, and a human-readable report.
//!
//! **Phase 33 — subprocess wrapper around `java -jar cello.jar`.**
//! Cello is JAR-distributed (no `cello` launcher binary on PATH);
//! the user supplies the absolute path to the jar via
//! `[bio.cello].jar` in `case.toml`. We probe `java` itself but
//! not the jar — different sites pin different Cello releases.
//!
//! `prepare()` composes a
//! `java -jar <jar> -inputNetlist <verilog> -targetDataFile <ucf>
//!  -inputSensorFile <in> -outputDeviceFile <out>
//!  -outputDir <output_basename> [extras...]`
//! invocation; `run()` streams the run via the shared subprocess
//! runner. `collect()` walks the workdir for the canonical
//! `<output_basename>*.txt` (report), `<output_basename>*.png`
//! (circuit diagram), and `<output_basename>*.dot` (netlist).

#![forbid(unsafe_code)]
#![allow(missing_docs)]

pub mod case_input;

use std::ffi::OsString;
use std::fs;
use std::path::Path;
use std::time::Duration;

use semver::Version;

use valenx_core::{
    adapter_helpers::{confined_join, find_on_path, live_provenance},
    error::RunPhase,
    subprocess, Adapter, AdapterError, AdapterInfo, Capabilities, Case, LicenseMode, Physics,
    PreparedJob, ProbeReport, RunContext, RunReport, VersionRange,
};
use valenx_fields::{
    artifact::{Artifact, ArtifactKind},
    Results,
};

use crate::case_input::CelloInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(CelloAdapter::new())
}

pub struct CelloAdapter;

impl CelloAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for CelloAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "cello";
/// Cello v2 is JAR-distributed — we probe `java` itself, not a
/// `cello` launcher. The user supplies the jar path via case
/// input.
const BINARIES: &[&str] = &["java"];

impl Adapter for CelloAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "Cello",
            // Cello v2 is the modern Java rewrite (2020+); the v1
            // line was Python and is deprecated. Upper bound 3.0
            // reserves room for an eventual major bump.
            version_range: VersionRange {
                min_inclusive: Version::new(2, 0, 0),
                max_exclusive: Version::new(3, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            tool_license: "BSD-3-Clause",
            docs_url: "https://github.com/CIDARLAB/Cello-v2",
            homepage_url: "https://www.cellocad.org/",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(BINARIES) {
            Some(binary_path) => Ok(ProbeReport {
                ok: true,
                // Cello's version comes from the jar itself; we
                // surface no version here. The user pins the
                // release implicitly by the jar they point at.
                found_version: None,
                binary_path: Some(binary_path),
                warnings: vec![
                    "probe found `java` on PATH but cannot verify the cello.jar \
                     release without invoking it; ensure `[bio.cello].jar` \
                     points at a valid Cello v2 distribution"
                        .into(),
                ],
                required_env: Vec::new(),
            }),
            None => Err(AdapterError::ToolNotInstalled {
                name: INFO_ID,
                hint: "Java 8+ JRE required to run Cello v2; install via \
                       your package manager (`apt install default-jre`, \
                       `brew install openjdk`, etc.) and ensure `java` is \
                       on PATH"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = CelloInput::from_case_dir(&case.path)?;

        // Round-4 security: reject `output_basename = "../etc/passwd"`
        // and friends before the value flows into any path join.
        // Same pattern as the round-3 fix in bionetgen/iqtree/art/fasttree.
        valenx_core::adapter_helpers::validate_output_basename(
            &input.output_basename,
            "[bio.cello].output_basename",
        )
        .map_err(|e| AdapterError::InvalidCase {
            case_path: case.path.join("case.toml"),
            reason: format!("{e}"),
        })?;

        fs::create_dir_all(workdir)?;

        // Resolve every case-supplied input path against the case
        // directory. The verilog / UCF / input-sensor / output-device
        // JSONs are user-authored data — sandbox them via
        // `confined_join` so a hostile shared bundle can't point at
        // `/etc/passwd`. The `jar` field is a system-install path
        // (the user points it at their `/opt/cello/cello.jar`), so
        // we still accept absolute paths there as documented; that
        // path is privileged-by-convention and not user data.
        let source_jar = if input.jar.is_absolute() {
            input.jar.clone()
        } else {
            case.path.join(&input.jar)
        };
        if !source_jar.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.cello].jar `{}` not found (resolved {})",
                    input.jar.display(),
                    source_jar.display()
                ),
            });
        }

        let source_verilog = confined_join(&case.path, &input.verilog)?;
        if !source_verilog.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.cello].verilog `{}` not found (resolved {})",
                    input.verilog.display(),
                    source_verilog.display()
                ),
            });
        }

        let source_ucf = confined_join(&case.path, &input.user_constraints)?;
        if !source_ucf.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.cello].user_constraints `{}` not found (resolved {})",
                    input.user_constraints.display(),
                    source_ucf.display()
                ),
            });
        }

        let source_in = confined_join(&case.path, &input.input_sensors)?;
        if !source_in.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.cello].input_sensors `{}` not found (resolved {})",
                    input.input_sensors.display(),
                    source_in.display()
                ),
            });
        }

        let source_out = confined_join(&case.path, &input.output_devices)?;
        if !source_out.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.cello].output_devices `{}` not found (resolved {})",
                    input.output_devices.display(),
                    source_out.display()
                ),
            });
        }

        let java_path = find_on_path(BINARIES).ok_or_else(|| AdapterError::ToolNotInstalled {
            name: INFO_ID,
            hint: "Java 8+ JRE required to run Cello v2; install via your \
                   package manager and ensure `java` is on PATH"
                .into(),
        })?;

        // Compose the canonical Cello v2 invocation. Cello dumps
        // results into the directory named by `-outputDir`, which
        // we anchor on the user's `output_basename` so collect()'s
        // prefix filter has a stable target.
        let mut native_command: Vec<OsString> = vec![
            java_path.into_os_string(),
            OsString::from("-jar"),
            source_jar.into_os_string(),
            OsString::from("-inputNetlist"),
            source_verilog.into_os_string(),
            OsString::from("-targetDataFile"),
            source_ucf.into_os_string(),
            OsString::from("-inputSensorFile"),
            source_in.into_os_string(),
            OsString::from("-outputDeviceFile"),
            source_out.into_os_string(),
            OsString::from("-outputDir"),
            OsString::from(&input.output_basename),
        ];
        for arg in &input.extra_args {
            native_command.push(OsString::from(arg));
        }

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            // Cello's simulated-annealing assignment search is
            // CPU-bound: small circuits compile in seconds, larger
            // ones (10+ gates over a deep gate library) can run for
            // an hour or more. 4 hours is a generous default.
            estimated_runtime: Some(Duration::from_secs(4 * 60 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting Cello", |line| {
            let mut hint = subprocess::Hint::default();
            // Cello's stdout chatter: a startup banner, per-stage
            // progress ("Stage: ASSIGN", "Stage: PLACE"), and a
            // final "Cello finished" / "Output written to..."
            // sentinel.
            if line.contains("Cello finished") || line.contains("Output written") {
                hint.progress = Some((95.0, line.to_string()));
            } else if line.contains("Stage: ") {
                hint.progress = Some((50.0, line.to_string()));
            } else if line.contains("Exception")
                || line.contains("Error")
                || line.contains("java.lang.")
            {
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
        // Provenance: hash the case.toml as the canonical input
        // descriptor. Cello writes results into a directory it
        // creates under the workdir, so we can't anchor on a
        // single fixed-name artifact.
        let case_hash_input = job.workdir.join("case.toml");
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "Cello",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);
        let mut artefacts: Vec<Artifact> = Vec::new();

        // Walk the workdir top-level + the basename-named
        // subdirectory Cello creates. Cello's three canonical
        // output families are `.txt` (human-readable report),
        // `.png` (rendered circuit diagram), and `.dot` (Graphviz
        // netlist). Restrict to outputs whose stem starts with the
        // configured `output_basename` so unrelated files don't
        // leak in.
        let basename = read_output_basename(&job.workdir);

        let mut roots: Vec<std::path::PathBuf> = vec![job.workdir.to_path_buf()];
        if let Some(b) = basename.as_deref() {
            let nested = job.workdir.join(b);
            if nested.is_dir() {
                roots.push(nested);
            }
        }
        for root in roots {
            let entries = match fs::read_dir(&root) {
                Ok(e) => e,
                Err(e) => {
                    tracing::warn!(target: "valenx-cello", ?e, "workdir read failed");
                    continue;
                }
            };
            for entry in entries.flatten() {
                let path = entry.path();
                if !path.is_file() {
                    continue;
                }
                let ext = path
                    .extension()
                    .and_then(|s| s.to_str())
                    .map(|s| s.to_ascii_lowercase());
                let stem = path
                    .file_stem()
                    .and_then(|s| s.to_str())
                    .unwrap_or("")
                    .to_string();
                // Files inside the basename-named subdirectory
                // are accepted unconditionally (Cello generates
                // many short-named files there); top-level
                // matches still need the prefix filter.
                let inside_basename_dir = root != job.workdir;
                let stem_ok = if inside_basename_dir {
                    true
                } else {
                    match basename.as_deref() {
                        Some(b) => stem.starts_with(b),
                        None => true,
                    }
                };
                if !stem_ok {
                    continue;
                }
                let (kind, label) = match ext.as_deref() {
                    Some("txt") => (ArtifactKind::Log, "Cello report".to_string()),
                    Some("png") => (ArtifactKind::Native, "Cello circuit diagram".to_string()),
                    Some("dot") => (ArtifactKind::Native, "Cello netlist".to_string()),
                    _ => continue,
                };
                artefacts.push(Artifact {
                    path,
                    kind,
                    checksum: None,
                    label,
                });
            }
        }
        artefacts.sort_by(|a, b| a.path.cmp(&b.path));
        results.artifacts = artefacts;
        Ok(results)
    }

    fn capabilities(&self) -> Capabilities {
        Capabilities {
            capabilities: Vec::new(),
            ribbon_contributions: vec!["bio.cello.compile"],
        }
    }
}

/// Re-read the `[bio.cello].output_basename` from a staged
/// `case.toml` for collect()-time output filtering. Returns `None`
/// when the case.toml is missing or unparseable — collect() then
/// accepts every matching file at the workdir top level.
fn read_output_basename(workdir: &Path) -> Option<String> {
    // Round-23 sweep: bound staged case.toml at MAX_PROJECT_FILE_BYTES.
    let text = valenx_core::io_caps::read_capped_to_string(
        &workdir.join("case.toml"),
        valenx_core::project::loader::MAX_PROJECT_FILE_BYTES as usize,
    )
    .ok()?;
    let parsed: toml::Value = toml::from_str(&text).ok()?;
    parsed
        .get("bio")?
        .get("cello")?
        .get("output_basename")?
        .as_str()
        .map(str::to_string)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn info_is_bio_domain() {
        let info = CelloAdapter::new().info();
        assert_eq!(info.id, "cello");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "BSD-3-Clause");
        assert_eq!(info.display_name, "Cello");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = CelloAdapter::new().info();
        // Cello v2 is the modern Java rewrite (2020+); upper
        // bound 3.0 reserves room for the next major.
        assert_eq!(info.version_range.min_inclusive, Version::new(2, 0, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(3, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = CelloAdapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.cello.compile"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = CelloAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }
}
