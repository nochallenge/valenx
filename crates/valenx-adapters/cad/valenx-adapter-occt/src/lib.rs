//! # valenx-adapter-occt
//!
//! Phase 2 scaffold + Python subprocess wrapper around
//! [pythonocc-core](https://github.com/tpaviot/pythonocc-core), the
//! canonical Python binding for OpenCASCADE Technology.
//!
//! The original LGPL-linked design (load OCCT's BRep kernel directly
//! into the Valenx binary via an `occt-sys` shim) lives on as a
//! future Phase 2 enhancement; today's adapter spawns a Python
//! interpreter that imports `OCC.Core` to do the work, keeping the
//! LGPL kernel out of the Valenx binary.
//!
//! `prepare()` resolves a user-supplied `.py` script against the
//! case directory, optionally stages an input geometry file
//! (`.step` / `.iges` / `.brep`) alongside it, writes a
//! `valenx_params.json` parameters file, and builds a
//! `python <script>` command. `run()` invokes that via the shared
//! [`valenx_core::subprocess`] runner. `collect()` walks the workdir
//! for output geometry / mesh files matching the user-declared
//! `output_basename` and surfaces them as typed artifacts.

#![forbid(unsafe_code)]
#![allow(missing_docs)]

pub mod case_input;

use std::ffi::OsString;
use std::fs;
use std::path::{Path, PathBuf};
use std::time::Duration;

use semver::Version;

use valenx_core::{
    adapter_helpers::{confined_join, find_on_path, live_provenance},
    error::RunPhase,
    subprocess, Adapter, AdapterError, AdapterInfo, Capabilities, Capability, Case, LicenseMode,
    Physics, PreparedJob, ProbeReport, RunContext, RunReport, VersionRange,
};
use valenx_fields::{
    artifact::{Artifact, ArtifactKind},
    Results,
};

use crate::case_input::OcctInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(OcctAdapter::new())
}

pub struct OcctAdapter;

impl OcctAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for OcctAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "occt";
/// Python interpreter candidates. `python3` first because on Linux
/// `python` may still be Python 2 on legacy distros; on Windows
/// `python` typically resolves to the Windows Store / 3.x install.
const PYTHON_BINARIES: &[&str] = &["python3", "python"];

impl Adapter for OcctAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "OpenCASCADE (OCCT)",
            // OCCT 7.6 onward is what current pythonocc-core wheels
            // are built against; we revisit the upper bound when
            // 8.x lands and the binding catches up.
            version_range: VersionRange {
                min_inclusive: Version::new(7, 6, 0),
                max_exclusive: Version::new(8, 0, 0),
            },
            physics: &[Physics::Geometry],
            license_mode: LicenseMode::Subprocess,
            tool_license: "LGPL-2.1-or-later WITH OCCT-exception-1.0",
            docs_url: "https://dev.opencascade.org/doc/refman/html/",
            homepage_url: "https://dev.opencascade.org/",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(PYTHON_BINARIES) {
            Some(binary_path) => {
                // Try `import OCC.Core` to confirm pythonocc-core
                // is actually available from this interpreter.
                // Failure degrades to a warning so the probe still
                // surfaces a useful state — the adapter is fine if
                // the user installs pythonocc-core later.
                let importable = pythonocc_importable(&binary_path);
                let mut warnings = Vec::new();
                if !importable {
                    warnings.push(
                        "probe found `python` on PATH but could not import \
                         `OCC.Core` — install pythonocc-core with \
                         `pip install pythonocc-core` or `conda install -c \
                         conda-forge pythonocc-core` for runs to succeed"
                            .into(),
                    );
                }
                Ok(ProbeReport {
                    ok: true,
                    // pythonocc-core's package version is decoupled
                    // from the OCCT kernel version; we surface
                    // `None` here rather than report something
                    // misleading. Detailed kernel version probing
                    // can land in a follow-up.
                    found_version: None,
                    binary_path: Some(binary_path),
                    warnings,
                    required_env: Vec::new(),
                })
            }
            None => Err(AdapterError::ToolNotInstalled {
                name: INFO_ID,
                hint: "Python 3.9+ with pythonocc-core installed; \
                       `pip install pythonocc-core` (or `conda install -c \
                       conda-forge pythonocc-core`) after ensuring python3 \
                       is on PATH"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = OcctInput::from_case_dir(&case.path)?;

        // Round-4 security: reject `output_basename = "../etc/passwd"`
        // and friends before the value flows into any path join.
        // Same pattern as the round-3 fix in bionetgen/iqtree/art/fasttree.
        valenx_core::adapter_helpers::validate_output_basename(
            &input.output_basename,
            "[cad.occt].output_basename",
        )
        .map_err(|e| AdapterError::InvalidCase {
            case_path: case.path.join("case.toml"),
            reason: format!("{e}"),
        })?;

        fs::create_dir_all(workdir)?;

        // Enforce a `.py` extension on the script. Keeps the
        // contract obvious (we feed the file to a Python interpreter)
        // and gives a structured error pointing at case.toml when
        // the user typos something else.
        let ext_ok = input
            .script
            .extension()
            .and_then(|e| e.to_str())
            .map(|e| e.eq_ignore_ascii_case("py"))
            .unwrap_or(false);
        if !ext_ok {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[cad.occt].script `{}` must end in `.py` — the OCCT \
                     adapter feeds the file to a Python interpreter that \
                     imports OCC.Core",
                    input.script.display()
                ),
            });
        }

        // Resolve + stage the script. `confined_join` rejects
        // absolute paths and `..` traversal so the staged copy
        // stays confined to the case directory.
        let source_script = confined_join(&case.path, &input.script)?;
        if !source_script.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[cad.occt].script `{}` not found (resolved {})",
                    input.script.display(),
                    source_script.display()
                ),
            });
        }
        let script_filename =
            input
                .script
                .file_name()
                .ok_or_else(|| AdapterError::InvalidCase {
                    case_path: case.path.join("case.toml"),
                    reason: format!(
                        "[cad.occt].script path `{}` has no filename",
                        input.script.display()
                    ),
                })?;
        let dest_script = workdir.join(script_filename);
        if source_script != dest_script {
            fs::copy(&source_script, &dest_script)?;
        }

        // Optional input geometry — same confined-join treatment.
        let staged_input_geometry: Option<PathBuf> = if let Some(g) = &input.input_geometry {
            let source_geom = confined_join(&case.path, g)?;
            if !source_geom.is_file() {
                return Err(AdapterError::InvalidCase {
                    case_path: case.path.join("case.toml"),
                    reason: format!(
                        "[cad.occt].input_geometry `{}` not found (resolved {})",
                        g.display(),
                        source_geom.display()
                    ),
                });
            }
            let geom_filename = g.file_name().ok_or_else(|| AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[cad.occt].input_geometry `{}` has no filename",
                    g.display()
                ),
            })?;
            let dest_geom = workdir.join(geom_filename);
            if source_geom != dest_geom {
                fs::copy(&source_geom, &dest_geom)?;
            }
            Some(PathBuf::from(geom_filename))
        } else {
            None
        };

        // Hand-rolled JSON for the parameters file. Keeps the
        // serde_json dep out of this crate and the format is
        // trivial — two scalar fields and one optional string.
        let params_path = workdir.join("valenx_params.json");
        let params_json = build_params_json(&input.output_basename, staged_input_geometry.as_ref());
        valenx_core::io_caps::atomic_write_str(&params_path, &params_json)?;

        // Resolve the Python binary. If the user pinned a specific
        // interpreter, honor it; otherwise walk PATH.
        // Round-4 security: validate python interpreter spec
        // against the allow-list AND resolve to a real binary
        // in one step. Closes the arbitrary-binary-exec class
        // that round-3 only patched in 8 of the 48 affected
        // adapters.
        let binary_path = valenx_core::adapter_helpers::resolve_python_binary(
            &input.python,
            PYTHON_BINARIES,
        )
        // Round-5: do NOT rewrap as ToolNotInstalled — the resolver
        // returns InvalidCase for allow-list rejections (which a hint
        // string would have hidden) and a clear Other for PATH lookup
        // failures. Pass the error through unchanged.
        ?;

        let native_command: Vec<OsString> = vec![
            binary_path.into_os_string(),
            OsString::from(script_filename),
        ];

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            // OCCT scripts span seconds (a single boolean op) to
            // many minutes (assemblies with mesh export). Pick a
            // 10-minute default; long runs can self-report progress.
            estimated_runtime: Some(Duration::from_secs(10 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting OCCT", |line| {
            let mut hint = subprocess::Hint::default();
            // Convention: the user-supplied script can emit a
            // sentinel line `[valenx] occt done` to signal
            // completion before exit; we lift that to a 95%
            // progress tick so the UI doesn't sit at
            // "indeterminate" until the process actually exits.
            if line.contains("[valenx] occt done") {
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
        // Re-derive the expected output_basename from the staged
        // params file. The case directory isn't on PreparedJob, so
        // we read what `prepare()` wrote into the workdir.
        let output_basename = read_output_basename(&job.workdir).unwrap_or_else(|| "model".into());

        // Provenance: prefer the staged script as the canonical
        // hashable input; fall back to the params file.
        let case_hash_input = first_script_in_workdir(&job.workdir)
            .or_else(|| {
                let p = job.workdir.join("valenx_params.json");
                p.is_file().then_some(p)
            })
            .unwrap_or_else(|| job.workdir.join("(no-input-found)"));
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "OCCT",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);
        let mut artefacts: Vec<Artifact> = Vec::new();

        let entries = match fs::read_dir(&job.workdir) {
            Ok(e) => e,
            Err(e) => {
                tracing::warn!(target: "valenx-occt", ?e, "workdir read failed");
                return Ok(results);
            }
        };
        for entry in entries.flatten() {
            let path = entry.path();
            if !path.is_file() {
                continue;
            }
            let stem = path
                .file_stem()
                .and_then(|s| s.to_str())
                .unwrap_or("")
                .to_string();
            let ext = path
                .extension()
                .and_then(|s| s.to_str())
                .map(|s| s.to_ascii_lowercase());

            // Files matching `<output_basename>*.<ext>`. Accept
            // stems that equal the basename or start with it
            // followed by any suffix (e.g. "model_part1.step").
            let matches_base =
                stem == output_basename || stem.starts_with(&format!("{output_basename}_"));

            let (kind, label) = match (matches_base, ext.as_deref()) {
                (true, Some("step")) | (true, Some("stp")) => {
                    (ArtifactKind::Native, "OCCT STEP geometry".to_string())
                }
                (true, Some("iges")) | (true, Some("igs")) => {
                    (ArtifactKind::Native, "OCCT IGES geometry".to_string())
                }
                (true, Some("brep")) => (ArtifactKind::Native, "OCCT BRep geometry".to_string()),
                (true, Some("stl")) => (ArtifactKind::Native, "OCCT STL mesh".to_string()),
                (true, Some("csv")) => {
                    (ArtifactKind::Tabular, "OCCT measurement table".to_string())
                }
                // `.log` files surface regardless of stem — the
                // subprocess runner / user script may pick its own
                // log filename and we still want the operator to
                // see it.
                (_, Some("log")) => (ArtifactKind::Log, "OCCT log".to_string()),
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
            capabilities: vec![
                Capability::GeoStep,
                Capability::GeoIges,
                Capability::GeoBRep,
                Capability::GeoStl,
            ],
            ribbon_contributions: vec!["cad.occt.run"],
        }
    }
}

/// Minimal hand-rolled JSON encoder for the params file. The schema
/// is `{"output_basename": "<stem>"[, "input_geometry": "<filename>"]}`
/// — `input_geometry` is omitted when `None` (we don't write
/// `null`).
fn build_params_json(output_basename: &str, input_geometry: Option<&PathBuf>) -> String {
    let mut s = String::new();
    s.push('{');
    s.push_str("\"output_basename\": ");
    push_json_string(&mut s, output_basename);
    if let Some(g) = input_geometry {
        s.push_str(", \"input_geometry\": ");
        push_json_string(&mut s, &g.to_string_lossy());
    }
    s.push('}');
    s
}

/// JSON-encode a string. Handles the small set of characters that
/// matter for our payload (filenames + a short identifier) — we
/// don't carry arbitrary user prose through here so a minimal
/// escape table keeps the helper tiny.
fn push_json_string(s: &mut String, value: &str) {
    s.push('"');
    for ch in value.chars() {
        match ch {
            '"' => s.push_str("\\\""),
            '\\' => s.push_str("\\\\"),
            '\n' => s.push_str("\\n"),
            '\r' => s.push_str("\\r"),
            '\t' => s.push_str("\\t"),
            c if (c as u32) < 0x20 => {
                s.push_str(&format!("\\u{:04x}", c as u32));
            }
            c => s.push(c),
        }
    }
    s.push('"');
}

/// Read `output_basename` back out of the staged params file in a
/// workdir. Returns `None` if the file is missing or malformed —
/// callers fall back to a sensible default rather than failing
/// `collect()`.
fn read_output_basename(workdir: &Path) -> Option<String> {
    let text = valenx_core::io_caps::read_capped_to_string(
        &workdir.join("valenx_params.json"),
        valenx_core::io_caps::MAX_ADAPTER_PARAMS_BYTES as usize,
    )
    .ok()?;
    // Tiny ad-hoc extractor; keeps the serde_json dep out of the
    // crate. Looks for `"output_basename": "<value>"`.
    let key = "\"output_basename\"";
    let key_pos = text.find(key)?;
    let after_key = &text[key_pos + key.len()..];
    let colon_pos = after_key.find(':')?;
    let after_colon = &after_key[colon_pos + 1..];
    let q1 = after_colon.find('"')?;
    let after_q1 = &after_colon[q1 + 1..];
    let q2 = after_q1.find('"')?;
    Some(after_q1[..q2].to_string())
}

/// Lift the staged Python script out of a workdir for provenance
/// hashing. Returns the lexicographically-first `.py` file at the
/// top level, or `None` if none exists yet.
fn first_script_in_workdir(workdir: &Path) -> Option<PathBuf> {
    let entries = fs::read_dir(workdir).ok()?;
    let mut hits: Vec<PathBuf> = entries
        .flatten()
        .map(|e| e.path())
        .filter(|p| {
            p.extension()
                .and_then(|s| s.to_str())
                .map(|s| s.eq_ignore_ascii_case("py"))
                .unwrap_or(false)
        })
        .collect();
    hits.sort();
    hits.into_iter().next()
}

/// Run `python -c "import OCC.Core"` and return whether it
/// succeeded. We don't extract a version because pythonocc-core's
/// release line is decoupled from the underlying OCCT kernel
/// version; reporting a misleading number is worse than reporting
/// `None`.
fn pythonocc_importable(python_binary: &Path) -> bool {
    std::process::Command::new(python_binary)
        .arg("-c")
        .arg("import OCC.Core")
        .stdin(std::process::Stdio::null())
        .stdout(std::process::Stdio::null())
        .stderr(std::process::Stdio::null())
        .status()
        .map(|s| s.success())
        .unwrap_or(false)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn info_is_geometry_domain() {
        let info = OcctAdapter::new().info();
        assert_eq!(info.id, "occt");
        assert_eq!(info.physics, &[Physics::Geometry]);
        assert_eq!(
            info.tool_license,
            "LGPL-2.1-or-later WITH OCCT-exception-1.0"
        );
        assert_eq!(info.display_name, "OpenCASCADE (OCCT)");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = OcctAdapter::new().info();
        // pythonocc-core wheels track OCCT 7.6+; revisit when 8.x lands.
        assert_eq!(info.version_range.min_inclusive, Version::new(7, 6, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(8, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = OcctAdapter::new().capabilities();
        // Geometry-format capabilities still apply: the adapter
        // emits STEP / IGES / BRep / STL as outputs.
        assert!(caps.capabilities.contains(&Capability::GeoStep));
        assert!(caps.capabilities.contains(&Capability::GeoIges));
        assert!(caps.capabilities.contains(&Capability::GeoBRep));
        assert!(caps.capabilities.contains(&Capability::GeoStl));
        assert_eq!(caps.ribbon_contributions, vec!["cad.occt.run"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = OcctAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }
}
