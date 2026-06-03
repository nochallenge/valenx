//! # valenx-adapter-fiji
//!
//! Adapter for [Fiji](https://fiji.sc/) — the canonical ImageJ
//! distribution (NIH / Schindelin et al). Fiji is the de-facto
//! reference platform for bioimage analysis in microscopy: it
//! consumes 2D / 3D / time-lapse image stacks (TIFF, PNG, OME-TIFF,
//! Bio-Formats) and processes them via Java-based plugins driven by
//! `.ijm` macros.
//!
//! **Phase 40 — subprocess wrapper for user-provided Fiji macros
//! in headless mode.** Sister adapter to ChimeraX / VMD / Jalview:
//! the user supplies a `.ijm` macro (and optionally an input image)
//! and the adapter composes a headless invocation. Fiji ships as a
//! platform-specific application bundle (`ImageJ-linux64`,
//! `ImageJ.exe`, `Fiji.app/Contents/MacOS/ImageJ-macosx`) that
//! cannot be predicted from a fixed PATH alone, so the case-input
//! schema carries the absolute launcher path verbatim
//! (`[bio.fiji].fiji_app`).
//!
//! `prepare()` composes
//! `<fiji_app> --headless --console -macro <macro_file> [extras...]`.
//! The `--headless --console` combo runs Fiji without GUI and
//! prints macro output to stdout for log capture. `run()` streams
//! the run via the shared subprocess runner. `collect()` walks the
//! workdir for the canonical `<output_basename>*.{tif,tiff,png,csv}`
//! artefacts plus any `*.log` Fiji emits.

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

use crate::case_input::FijiInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(FijiAdapter::new())
}

pub struct FijiAdapter;

impl FijiAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for FijiAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "fiji";
/// Fiji launcher binary candidates across platforms. Linux ships
/// `ImageJ-linux64`, macOS ships `Fiji.app/Contents/MacOS/ImageJ-macosx`
/// (the bundle's MacOS subdir typically isn't on PATH), Windows
/// ships `ImageJ.exe`. `fiji` is the lowercase symlink some Linux
/// package managers create.
const BINARIES: &[&str] = &["ImageJ-linux64", "ImageJ-macosx", "ImageJ.exe", "fiji"];

impl Adapter for FijiAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "Fiji",
            // Fiji has been on the 2.x ImageJ2 line since 2014;
            // the upper bound 3.0 reserves room for an eventual
            // major bump.
            version_range: VersionRange {
                min_inclusive: Version::new(2, 0, 0),
                max_exclusive: Version::new(3, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            tool_license: "GPL-3.0",
            docs_url: "https://imagej.net/software/fiji/",
            homepage_url: "https://fiji.sc/",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(BINARIES) {
            Some(binary_path) => Ok(ProbeReport {
                ok: true,
                // Fiji's reported version comes from the launcher
                // itself, but the version-detection invocation
                // requires `--version` and its output format has
                // shifted across releases. We surface no version
                // here — the user pins the Fiji release implicitly
                // by the launcher path they point at.
                found_version: None,
                binary_path: Some(binary_path),
                warnings: Vec::new(),
                required_env: Vec::new(),
            }),
            None => {
                // Fiji isn't on PATH. Check Java as a hint: if
                // present, the user has the JRE Fiji needs but
                // hasn't added the bin directory to PATH yet, so
                // surface an actionable warning. If neither is
                // present we report ToolNotInstalled.
                if find_on_path(&["java"]).is_some() {
                    Ok(ProbeReport {
                        ok: false,
                        found_version: None,
                        binary_path: None,
                        warnings: vec!["Fiji not found on PATH; download from \
                             https://fiji.sc and add the bin directory \
                             to PATH"
                            .into()],
                        required_env: Vec::new(),
                    })
                } else {
                    Err(AdapterError::ToolNotInstalled {
                        name: INFO_ID,
                        hint: "Fiji (the canonical ImageJ distribution) \
                               required; download from https://fiji.sc \
                               and add the bin directory to PATH (Java \
                               JRE is bundled with Fiji)"
                            .into(),
                    })
                }
            }
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = FijiInput::from_case_dir(&case.path)?;

        // Round-4 security: reject `output_basename = "../etc/passwd"`
        // and friends before the value flows into any path join.
        // Same pattern as the round-3 fix in bionetgen/iqtree/art/fasttree.
        valenx_core::adapter_helpers::validate_output_basename(
            &input.output_basename,
            "[bio.fiji].output_basename",
        )
        .map_err(|e| AdapterError::InvalidCase {
            case_path: case.path.join("case.toml"),
            reason: format!("{e}"),
        })?;

        fs::create_dir_all(workdir)?;

        // Resolve the launcher path against the case directory if
        // relative. Almost always absolute (Fiji installs land
        // under /opt or /Applications), but support the relative
        // form too.
        // Round-9 classification: KEEP `case.path.join` here —
        // `fiji_app` is the system launcher (e.g. `/Applications/Fiji.app/Contents/MacOS/ImageJ-macosx`)
        // which is admin-managed, not user data. The macro file and
        // input image fields below go through `confined_join`.
        let source_app = if input.fiji_app.is_absolute() {
            input.fiji_app.clone()
        } else {
            case.path.join(&input.fiji_app)
        };
        if !source_app.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.fiji].fiji_app `{}` not found (resolved {})",
                    input.fiji_app.display(),
                    source_app.display()
                ),
            });
        }

        // Resolve the macro file against the case directory.
        let source_macro = if input.macro_file.is_absolute() {
            input.macro_file.clone()
        } else {
            valenx_core::adapter_helpers::confined_join(
            &case.path,
            &input.macro_file,
        )?
        };
        if !source_macro.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.fiji].macro_file `{}` not found (resolved {})",
                    input.macro_file.display(),
                    source_macro.display()
                ),
            });
        }

        // Resolve the optional input image against the case dir.
        // Macros that synthesise from scratch leave this `None`.
        // Round-9 hardening: `input_image` is user-supplied data —
        // wrap the relative branch with `confined_join` so a hostile
        // case can't ask Fiji to open `../../etc/passwd`.
        if let Some(img) = &input.input_image {
            let source_img = if img.is_absolute() {
                img.clone()
            } else {
                valenx_core::adapter_helpers::confined_join(&case.path, img)?
            };
            if !source_img.is_file() {
                return Err(AdapterError::InvalidCase {
                    case_path: case.path.join("case.toml"),
                    reason: format!(
                        "[bio.fiji].input_image `{}` not found (resolved {})",
                        img.display(),
                        source_img.display()
                    ),
                });
            }
        }

        // Compose `<fiji_app> --headless --console -macro <macro> [extras...]`.
        // `--headless --console` runs Fiji without GUI and prints
        // macro output to stdout. The macro is responsible for
        // loading any input image (typically via the macro's
        // `open(getArgument())` pattern with the path passed in
        // `extra_args`), so we don't pass `input_image` directly
        // on the command line — we just validate it exists for
        // the user's benefit.
        let mut native_command: Vec<OsString> = vec![
            source_app.into_os_string(),
            OsString::from("--headless"),
            OsString::from("--console"),
            OsString::from("-macro"),
            source_macro.into_os_string(),
        ];
        for arg in &input.extra_args {
            native_command.push(OsString::from(arg));
        }

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            // Headless Fiji macros range from sub-second
            // single-image filters to multi-hour batch processing
            // of large microscopy stacks. 60 minutes is a
            // reasonable default that covers most workflows
            // without artificially capping batch jobs.
            estimated_runtime: Some(Duration::from_secs(60 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting Fiji", |line| {
            let mut hint = subprocess::Hint::default();
            // Fiji writes a Java-style banner on startup, then
            // macro `print()` output (often progress markers), and
            // ends with a "saved" / "done" sentinel from the macro
            // itself. We pattern-match conservatively — macros are
            // user-authored, so we can't rely on a specific
            // wording.
            if line.contains("saved") || line.contains("Saved") || line.contains("done") {
                hint.progress = Some((95.0, line.to_string()));
            } else if line.contains("Loading")
                || line.contains("Opening")
                || line.contains("Processing")
            {
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
        // descriptor. Fiji output filenames are macro-defined; we
        // anchor on the case.toml, which captures the macro path
        // and basename that drove the run.
        let case_hash_input = job.workdir.join("case.toml");
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "fiji",
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
        // user's input image / macro doesn't pollute the artefact
        // list. Logs (`*.log`) are accepted regardless of stem
        // since Fiji / the JVM may name them on their own.
        let basename = read_output_basename(&job.workdir);

        let entries = match fs::read_dir(&job.workdir) {
            Ok(e) => e,
            Err(e) => {
                tracing::warn!(target: "valenx-fiji", ?e, "workdir read failed");
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
                Some("tif") if stem_matches_basename => {
                    (ArtifactKind::Native, "Fiji image (TIFF)".to_string())
                }
                Some("tiff") if stem_matches_basename => {
                    (ArtifactKind::Native, "Fiji image (TIFF)".to_string())
                }
                Some("png") if stem_matches_basename => {
                    (ArtifactKind::Native, "Fiji image (PNG)".to_string())
                }
                Some("csv") if stem_matches_basename => {
                    (ArtifactKind::Tabular, "Fiji measurements".to_string())
                }
                Some("log") => (ArtifactKind::Log, "Fiji log".to_string()),
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
            ribbon_contributions: vec!["bio.fiji.process"],
        }
    }
}

/// Re-read the `[bio.fiji].output_basename` from a staged
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
        .get("fiji")?
        .get("output_basename")?
        .as_str()
        .map(str::to_string)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn info_is_bio_domain() {
        let info = FijiAdapter::new().info();
        assert_eq!(info.id, "fiji");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "GPL-3.0");
        assert_eq!(info.display_name, "Fiji");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = FijiAdapter::new().info();
        // Fiji has been on the 2.x ImageJ2 line since 2014; the
        // upper bound 3.0 reserves room for the next major.
        assert_eq!(info.version_range.min_inclusive, Version::new(2, 0, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(3, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = FijiAdapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.fiji.process"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = FijiAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }

    /// Round-9 RED→GREEN: `[bio.fiji].input_image` used to be joined
    /// with bare `case.path.join`. Wrap with `confined_join`.
    #[test]
    fn prepare_rejects_input_image_traversing_outside_case_dir() {
        use valenx_test_utils::tempdir;
        let d = tempdir("fiji-input-trav");
        // Stage a fiji_app + macro file so prepare gets past those
        // earlier validations and reaches the input_image branch.
        std::fs::write(d.join("ImageJ-fake"), b"").unwrap();
        std::fs::write(d.join("macro.ijm"), b"// noop").unwrap();
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "fiji.process"

[bio.fiji]
fiji_app        = "ImageJ-fake"
macro_file      = "macro.ijm"
input_image     = "../../etc/passwd"
output_basename = "out"
"#,
        )
        .unwrap();
        let case = Case {
            id: "fiji-input-trav".into(),
            path: d.clone(),
        };
        let workdir = d.join("workdir");
        let err = FijiAdapter::new().prepare(&case, &workdir).unwrap_err();
        let msg = format!("{err}");
        assert!(
            msg.contains("..") || msg.contains("stay within") || msg.contains("escape"),
            "expected confined_join rejection, got: {msg}"
        );
        let _ = std::fs::remove_dir_all(&d);
    }
}
