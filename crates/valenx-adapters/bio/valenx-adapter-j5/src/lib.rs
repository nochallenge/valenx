//! # valenx-adapter-j5
//!
//! Adapter for [j5](https://j5.jbei.org/) — JBEI's canonical
//! DNA-assembly automation tool. j5 consumes a target circuit
//! design (CSV row per cassette) plus a parts library
//! (CSV row per part / oligo), then plans the optimal
//! Gibson / Golden-Gate / SLIC / SLIM assembly strategy and writes
//! the per-step protocol + GenBank construct files.
//!
//! **Phase 33 — subprocess wrapper around `java -jar j5.jar`.**
//! j5 is JAR-distributed (no `j5` launcher binary on PATH); the
//! user supplies the absolute path to `j5.jar` via
//! `[bio.j5].jar` in `case.toml`. We probe that `java` is on
//! PATH but not the jar itself — different sites pin different
//! j5 releases under different paths.
//!
//! `prepare()` composes a
//! `java -jar <jar> -d <design_csv> -p <parts_csv> -o <output_basename> [extras...]`
//! invocation; `run()` streams the run via the shared subprocess
//! runner. `collect()` walks the workdir for the canonical
//! `<output_basename>*.csv` (assembly plan) and
//! `<output_basename>*.gb` (GenBank constructs).

#![forbid(unsafe_code)]
#![allow(missing_docs)]

pub mod case_input;

use std::ffi::OsString;
use std::fs;
use std::path::Path;
use std::time::Duration;

use semver::Version;

use valenx_core::{
    adapter_helpers::{find_on_path, live_provenance},
    error::RunPhase,
    subprocess, Adapter, AdapterError, AdapterInfo, Capabilities, Case, LicenseMode, Physics,
    PreparedJob, ProbeReport, RunContext, RunReport, VersionRange,
};
use valenx_fields::{
    artifact::{Artifact, ArtifactKind},
    Results,
};

use crate::case_input::J5Input;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(J5Adapter::new())
}

pub struct J5Adapter;

impl J5Adapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for J5Adapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "j5";
/// j5 is JAR-distributed — we probe `java` itself, not a
/// `j5` launcher. The user supplies the jar path via case input.
const BINARIES: &[&str] = &["java"];

impl Adapter for J5Adapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "j5",
            // j5 has been on a 1.x line for over a decade; the
            // upper bound 2.0 reserves room for an eventual major
            // bump.
            version_range: VersionRange {
                min_inclusive: Version::new(1, 0, 0),
                max_exclusive: Version::new(2, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            tool_license: "BSD-3-Clause",
            docs_url: "https://j5.jbei.org/index.php/Main_Page",
            homepage_url: "https://j5.jbei.org/",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(BINARIES) {
            Some(binary_path) => Ok(ProbeReport {
                ok: true,
                // j5's version comes from the jar itself, not from
                // `java`; we surface no version here. The user pins
                // the j5 release implicitly by the jar they point at.
                found_version: None,
                binary_path: Some(binary_path),
                warnings: vec!["probe found `java` on PATH but cannot verify the j5.jar \
                     release without invoking it; ensure `[bio.j5].jar` points \
                     at a valid j5 distribution"
                    .into()],
                required_env: Vec::new(),
            }),
            None => Err(AdapterError::ToolNotInstalled {
                name: INFO_ID,
                hint: "Java 8+ JRE required to run j5; install via your \
                       package manager (`apt install default-jre`, \
                       `brew install openjdk`, etc.) and ensure `java` is \
                       on PATH"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = J5Input::from_case_dir(&case.path)?;

        // Round-4 security: reject `output_basename = "../etc/passwd"`
        // and friends before the value flows into any path join.
        // Same pattern as the round-3 fix in bionetgen/iqtree/art/fasttree.
        valenx_core::adapter_helpers::validate_output_basename(
            &input.output_basename,
            "[bio.j5].output_basename",
        )
        .map_err(|e| AdapterError::InvalidCase {
            case_path: case.path.join("case.toml"),
            reason: format!("{e}"),
        })?;

        fs::create_dir_all(workdir)?;

        // Resolve the jar path against the case directory if relative.
        // Almost always absolute (jars live under /opt or similar),
        // but support the relative form too.
        let source_jar = if input.jar.is_absolute() {
            input.jar.clone()
        } else {
            case.path.join(&input.jar)
        };
        if !source_jar.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.j5].jar `{}` not found (resolved {})",
                    input.jar.display(),
                    source_jar.display()
                ),
            });
        }

        // Resolve the two CSV inputs against the case directory.
        let source_design = if input.design_csv.is_absolute() {
            input.design_csv.clone()
        } else {
            valenx_core::adapter_helpers::confined_join(
            &case.path,
            &input.design_csv,
        )?
        };
        if !source_design.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.j5].design_csv `{}` not found (resolved {})",
                    input.design_csv.display(),
                    source_design.display()
                ),
            });
        }
        let source_parts = if input.parts_csv.is_absolute() {
            input.parts_csv.clone()
        } else {
            valenx_core::adapter_helpers::confined_join(
            &case.path,
            &input.parts_csv,
        )?
        };
        if !source_parts.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.j5].parts_csv `{}` not found (resolved {})",
                    input.parts_csv.display(),
                    source_parts.display()
                ),
            });
        }

        let java_path = find_on_path(BINARIES).ok_or_else(|| AdapterError::ToolNotInstalled {
            name: INFO_ID,
            hint: "Java 8+ JRE required to run j5; install via your \
                   package manager and ensure `java` is on PATH"
                .into(),
        })?;

        // Compose `java -jar <jar> -d <design> -p <parts>
        //          -o <output_basename> [extras...]`. Outputs
        // land in the cwd (the workdir), tagged with the
        // basename stem.
        let mut native_command: Vec<OsString> = vec![
            java_path.into_os_string(),
            OsString::from("-jar"),
            source_jar.into_os_string(),
            OsString::from("-d"),
            source_design.into_os_string(),
            OsString::from("-p"),
            source_parts.into_os_string(),
            OsString::from("-o"),
            OsString::from(&input.output_basename),
        ];
        for arg in &input.extra_args {
            native_command.push(OsString::from(arg));
        }

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            // j5 assembly planning is interactive-fast for typical
            // designs (seconds) and bound by JVM startup; large
            // multi-cassette libraries can run for several minutes.
            // 30 minutes is a generous default.
            estimated_runtime: Some(Duration::from_secs(30 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting j5", |line| {
            let mut hint = subprocess::Hint::default();
            // j5 writes a typical Java application banner on
            // startup, then per-cassette progress lines, then a
            // "j5 finished" / "Wrote ..." sentinel at end-of-run.
            if line.contains("j5 finished") || line.contains("Wrote ") {
                hint.progress = Some((95.0, line.to_string()));
            } else if line.contains("Processing") || line.contains("Cassette") {
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
        // descriptor. j5's output filenames are derived from the
        // basename plus the cassette name, so we can't anchor on
        // a single fixed-name artifact.
        let case_hash_input = job.workdir.join("case.toml");
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "j5",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);
        let mut artefacts: Vec<Artifact> = Vec::new();

        // Restrict CSV / GenBank outputs to those whose stem starts
        // with the configured `output_basename` so the user's input
        // CSVs don't pollute the artefact list.
        let basename = read_output_basename(&job.workdir);

        let entries = match fs::read_dir(&job.workdir) {
            Ok(e) => e,
            Err(e) => {
                tracing::warn!(target: "valenx-j5", ?e, "workdir read failed");
                return Ok(results);
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
            let stem_ok = match basename.as_deref() {
                Some(b) => stem.starts_with(b),
                None => true,
            };
            if !stem_ok {
                continue;
            }
            let (kind, label) = match ext.as_deref() {
                Some("csv") => (ArtifactKind::Tabular, "j5 assembly plan".to_string()),
                Some("gb") => (ArtifactKind::Native, "j5 GenBank output".to_string()),
                _ => continue,
            };
            artefacts.push(Artifact {
                path,
                kind,
                checksum: None,
                label,
            });
        }
        artefacts.sort_by(|a, b| a.path.cmp(&b.path));
        results.artifacts = artefacts;
        Ok(results)
    }

    fn capabilities(&self) -> Capabilities {
        Capabilities {
            capabilities: Vec::new(),
            ribbon_contributions: vec!["bio.j5.assemble"],
        }
    }
}

/// Re-read the `[bio.j5].output_basename` from a staged `case.toml`
/// for collect()-time output filtering. Returns `None` when the
/// case.toml is missing or unparseable — collect() then accepts
/// every CSV / GenBank file in the workdir (best-effort).
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
        .get("j5")?
        .get("output_basename")?
        .as_str()
        .map(str::to_string)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn info_is_bio_domain() {
        let info = J5Adapter::new().info();
        assert_eq!(info.id, "j5");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "BSD-3-Clause");
        assert_eq!(info.display_name, "j5");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = J5Adapter::new().info();
        // j5 has been on a 1.x line for over a decade; upper
        // bound 2.0 reserves room for the next major.
        assert_eq!(info.version_range.min_inclusive, Version::new(1, 0, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(2, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = J5Adapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.j5.assemble"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = J5Adapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }
}
