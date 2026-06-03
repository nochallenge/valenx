//! # valenx-adapter-jalview
//!
//! Adapter for [Jalview](https://www.jalview.org/) — the
//! Barton group's Java-based multiple sequence alignment viewer.
//! Jalview consumes a multiple sequence alignment (FASTA /
//! Clustal / Stockholm / etc.) and renders it as an image
//! (PNG / SVG), an HTML page, or re-emits the alignment in any
//! of its supported formats.
//!
//! **Phase 41 — subprocess wrapper around `java -jar jalview.jar`
//! in headless `-nodisplay` mode.** Jalview is JAR-distributed
//! (no `jalview` launcher binary on PATH for headless work);
//! the user supplies the absolute path to `jalview.jar` via
//! `[bio.jalview].jar` in `case.toml`. We probe that `java` is
//! on PATH but not the jar itself — different sites pin different
//! Jalview releases under different paths.
//!
//! `prepare()` composes a
//! `java -jar <jar> -nodisplay -open <input>
//!  -<output_format> <output_basename>.<ext> [extras...]`
//! invocation; `run()` streams the run via the shared subprocess
//! runner. `collect()` walks the workdir for the canonical
//! `<output_basename>*` artefacts (image, HTML, re-formatted
//! alignment) plus any `*.log` Jalview emits.

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

use crate::case_input::JalviewInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(JalviewAdapter::new())
}

pub struct JalviewAdapter;

impl JalviewAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for JalviewAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "jalview";
/// Jalview is JAR-distributed — we probe `java` itself, not a
/// `jalview` launcher. The user supplies the jar path via case
/// input.
const BINARIES: &[&str] = &["java"];

impl Adapter for JalviewAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "Jalview",
            // Jalview has been on the 2.11.x line since 2018;
            // the upper bound 3.0 reserves room for an eventual
            // major bump.
            version_range: VersionRange {
                min_inclusive: Version::new(2, 11, 0),
                max_exclusive: Version::new(3, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            tool_license: "GPL-3.0",
            docs_url: "https://www.jalview.org/help.html",
            homepage_url: "https://www.jalview.org/",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(BINARIES) {
            Some(binary_path) => Ok(ProbeReport {
                ok: true,
                // Jalview's version comes from the jar itself, not
                // from `java`; we surface no version here. The user
                // pins the Jalview release implicitly by the jar
                // they point at.
                found_version: None,
                binary_path: Some(binary_path),
                warnings: vec![
                    "probe found `java` on PATH but cannot verify the jalview.jar \
                     release without invoking it; ensure `[bio.jalview].jar` points \
                     at a valid Jalview distribution"
                        .into(),
                ],
                required_env: Vec::new(),
            }),
            None => Err(AdapterError::ToolNotInstalled {
                name: INFO_ID,
                hint: "Java 11+ JRE required to run Jalview; install via your \
                       package manager (`apt install default-jre`, \
                       `brew install openjdk`, etc.) and ensure `java` is \
                       on PATH"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = JalviewInput::from_case_dir(&case.path)?;

        // Round-4 security: reject `output_basename = "../etc/passwd"`
        // and friends before the value flows into any path join.
        // Same pattern as the round-3 fix in bionetgen/iqtree/art/fasttree.
        valenx_core::adapter_helpers::validate_output_basename(
            &input.output_basename,
            "[bio.jalview].output_basename",
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
                    "[bio.jalview].jar `{}` not found (resolved {})",
                    input.jar.display(),
                    source_jar.display()
                ),
            });
        }

        // Resolve the alignment input against the case directory.
        let source_input = if input.input.is_absolute() {
            input.input.clone()
        } else {
            valenx_core::adapter_helpers::confined_join(
            &case.path,
            &input.input,
        )?
        };
        if !source_input.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.jalview].input `{}` not found (resolved {})",
                    input.input.display(),
                    source_input.display()
                ),
            });
        }

        let java_path = find_on_path(BINARIES).ok_or_else(|| AdapterError::ToolNotInstalled {
            name: INFO_ID,
            hint: "Java 11+ JRE required to run Jalview; install via your \
                   package manager and ensure `java` is on PATH"
                .into(),
        })?;

        // Map output_format → file extension. Jalview's headless
        // flag and the natural file suffix differ for `clustal`
        // (the alignment file convention is `.aln`).
        let ext = output_extension(&input.output_format);
        let output_filename = format!("{}.{}", input.output_basename, ext);
        let output_flag = format!("-{}", input.output_format);

        // Compose `java -jar <jar> -nodisplay -open <input>
        //          -<format> <basename>.<ext> [extras...]`. Outputs
        // land in the cwd (the workdir).
        let mut native_command: Vec<OsString> = vec![
            java_path.into_os_string(),
            OsString::from("-jar"),
            source_jar.into_os_string(),
            OsString::from("-nodisplay"),
            OsString::from("-open"),
            source_input.into_os_string(),
            OsString::from(&output_flag),
            OsString::from(&output_filename),
        ];
        for arg in &input.extra_args {
            native_command.push(OsString::from(arg));
        }

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            // Headless Jalview rendering is interactive-fast for
            // typical alignments (seconds) and bound by JVM
            // startup; very large multi-thousand-sequence
            // alignments can run for several minutes. 30 minutes
            // is a generous default.
            estimated_runtime: Some(Duration::from_secs(30 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting Jalview", |line| {
            let mut hint = subprocess::Hint::default();
            // Jalview writes a typical Java application banner on
            // startup, then alignment-loading progress lines, then
            // a "Wrote" / image-saved sentinel at end-of-run.
            if line.contains("Wrote ") || line.contains("Saved ") {
                hint.progress = Some((95.0, line.to_string()));
            } else if line.contains("Loading") || line.contains("Opening") {
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
        // descriptor. Jalview's output filename is derived from
        // the basename plus the format-specific extension, so we
        // anchor on the case.toml.
        let case_hash_input = job.workdir.join("case.toml");
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "jalview",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);
        let mut artefacts: Vec<Artifact> = Vec::new();

        // Restrict format-specific outputs to those whose stem
        // starts with the configured `output_basename` so the
        // user's input alignment doesn't pollute the artefact
        // list. Logs (`*.log`) are accepted regardless of stem
        // since Jalview names them on its own.
        let basename = read_output_basename(&job.workdir);

        let entries = match fs::read_dir(&job.workdir) {
            Ok(e) => e,
            Err(e) => {
                tracing::warn!(target: "valenx-jalview", ?e, "workdir read failed");
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
            let stem_matches_basename = match basename.as_deref() {
                Some(b) => stem.starts_with(b),
                None => true,
            };
            let (kind, label) = match ext.as_deref() {
                Some("png") if stem_matches_basename => {
                    (ArtifactKind::Native, "Jalview alignment image".to_string())
                }
                Some("svg") if stem_matches_basename => {
                    (ArtifactKind::Native, "Jalview SVG".to_string())
                }
                Some("html") if stem_matches_basename => {
                    (ArtifactKind::Native, "Jalview HTML".to_string())
                }
                Some("fasta") if stem_matches_basename => {
                    (ArtifactKind::Native, "Jalview FASTA".to_string())
                }
                Some("aln") if stem_matches_basename => {
                    (ArtifactKind::Tabular, "Jalview alignment".to_string())
                }
                Some("log") => (ArtifactKind::Log, "Jalview log".to_string()),
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
            ribbon_contributions: vec!["bio.jalview.view"],
        }
    }
}

/// Map a Jalview output format flag (`png`, `html`, `svg`,
/// `fasta`, `clustal`) to the conventional file extension.
/// `clustal` is the format flag but `.aln` is the file
/// convention. Anything we don't recognise falls back to using
/// the format string itself as the extension (defensive default).
fn output_extension(format: &str) -> String {
    match format {
        "png" => "png".to_string(),
        "html" => "html".to_string(),
        "svg" => "svg".to_string(),
        "fasta" => "fasta".to_string(),
        "clustal" => "aln".to_string(),
        other => other.to_string(),
    }
}

/// Re-read the `[bio.jalview].output_basename` from a staged
/// `case.toml` for collect()-time output filtering. Returns
/// `None` when the case.toml is missing or unparseable —
/// collect() then accepts every recognised file in the workdir
/// (best-effort).
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
        .get("jalview")?
        .get("output_basename")?
        .as_str()
        .map(str::to_string)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn info_is_bio_domain() {
        let info = JalviewAdapter::new().info();
        assert_eq!(info.id, "jalview");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "GPL-3.0");
        assert_eq!(info.display_name, "Jalview");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = JalviewAdapter::new().info();
        // Jalview has been on the 2.11.x line since 2018; the
        // upper bound 3.0 reserves room for the next major.
        assert_eq!(info.version_range.min_inclusive, Version::new(2, 11, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(3, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = JalviewAdapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.jalview.view"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = JalviewAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }
}
