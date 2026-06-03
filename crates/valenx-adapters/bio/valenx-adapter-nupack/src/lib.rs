//! # valenx-adapter-nupack
//!
//! Adapter for [NUPACK](https://www.nupack.org/) — Niles Pierce's
//! lab at Caltech ships the canonical software stack for nucleic-
//! acid thermodynamics, complex-equilibrium analysis, and
//! sequence design. Where ViennaRNA and RNAstructure focus on a
//! single sequence's MFE secondary structure, NUPACK extends to
//! multi-strand systems (test tubes, complexes, target structures)
//! and inverts the relationship: given a target structure, NUPACK's
//! design routines invert the thermodynamic model to propose
//! sequences that fold into it.
//!
//! **Phase 28 — subprocess wrapper for user-provided Python scripts.**
//! NUPACK's modern API surface is a Python package; there's no
//! canonical CLI. The user authors `design.py` (or whatever
//! filename) referenced from `[bio.nupack].script`, which calls
//! into `nupack.tubes`, `nupack.complex`, `nupack.thermodynamics`,
//! etc. `prepare()` stages the script (and optional input file) into
//! the workdir, drops a flat `valenx_params.json` with the parsed
//! knobs, and `run()` invokes `python <script>` via the shared
//! subprocess runner.
//!
//! ## `valenx_params.json`
//!
//! Same convention as RFdiffusion / DeepChem / ESMFold — a flat
//! JSON object the user's script reads with
//! `json.load(open("valenx_params.json"))`:
//!
//! ```json
//! {
//!   "input":           "complex.fa",
//!   "output_basename": "design",
//!   "temperature":     37.0,
//!   "sodium":          1.0
//! }
//! ```
//!
//! Scripts pass the values through to NUPACK themselves. Keeping
//! the adapter free of upstream API churn means `case.toml` knobs
//! actually reach the model.
//!
//! ## License flag
//!
//! NUPACK ships under a custom non-OSS license that restricts
//! commercial redistribution to academic / non-commercial contexts.
//! We surface this via `tool_license = "NUPACK-License"` and emit a
//! probe warning whenever `nupack` is importable. The probe-warning
//! text contains the literal string `"academic"` as a stable anchor
//! for tests and downstream license-aware filters.
//!
//! On `collect()` we walk the workdir for files starting with
//! `<output_basename>*` and classify them by extension: `.json`
//! (analysis log), `.npc` (NUPACK proprietary container), `.csv` /
//! `.tsv` (tabular) — anything else with the basename prefix is
//! surfaced as Native with an "NUPACK output" label.

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
    subprocess, Adapter, AdapterError, AdapterInfo, Capabilities, Case, LicenseMode, Physics,
    PreparedJob, ProbeReport, RunContext, RunReport, VersionRange,
};
use valenx_fields::{
    artifact::{Artifact, ArtifactKind},
    Results,
};

use crate::case_input::NupackInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(NupackAdapter::new())
}

pub struct NupackAdapter;

impl NupackAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for NupackAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "nupack";
/// Python interpreter candidates. `python3` first because on Linux
/// `python` may still be Python 2 on legacy distros; on Windows
/// `python` typically resolves to the Windows Store / 3.x install.
const PYTHON_BINARIES: &[&str] = &["python3", "python"];

/// The probe-warning surfaced whenever NUPACK is detected. The
/// literal string `"academic"` is part of the asserted contract — it
/// anchors the license reminder so downstream license-aware filters
/// and tests can key off a stable substring.
const LICENSE_WARNING: &str = "NUPACK is licensed for non-commercial / academic use only. \
     Caltech's NUPACK license restricts redistribution + commercial \
     use; confirm your use case complies before publishing analyses.";

impl Adapter for NupackAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "NUPACK",
            // NUPACK 4.0 (2020) is the modern Python-first rewrite;
            // every supported feature in this adapter targets the
            // 4.x API surface. Upper-bound 5.0 reserves room for an
            // eventual major bump.
            version_range: VersionRange {
                min_inclusive: Version::new(4, 0, 0),
                max_exclusive: Version::new(5, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            // NUPACK's terms aren't a recognised SPDX identifier; the
            // closest accurate label is the project's own custom
            // license. Mislabeling as MIT / BSD would be misleading.
            tool_license: "NUPACK-License",
            docs_url: "https://docs.nupack.org/",
            homepage_url: "https://www.nupack.org/",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(PYTHON_BINARIES) {
            Some(binary_path) => {
                // Try `import nupack; print(nupack.__version__)`
                // first — that's the only string that confirms NUPACK
                // itself is installed (vs. just having Python on
                // PATH). If the import fails, surface a successful
                // probe with an install-hint warning so the user
                // gets actionable guidance.
                let found_version = detect_nupack_version(&binary_path);
                let mut warnings = vec![LICENSE_WARNING.to_string()];
                if found_version.is_none() {
                    warnings.push(
                        "probe found `python` on PATH but could not import \
                         `nupack` — install NUPACK from \
                         https://www.nupack.org/ (registration required, \
                         academic-use license) and `pip install nupack-*.whl` \
                         for runs to succeed"
                            .into(),
                    );
                }
                Ok(ProbeReport {
                    ok: true,
                    found_version,
                    binary_path: Some(binary_path),
                    warnings,
                    required_env: Vec::new(),
                })
            }
            None => Err(AdapterError::ToolNotInstalled {
                name: INFO_ID,
                hint: "Python 3.9+ with NUPACK installed; download NUPACK from \
                       https://www.nupack.org/ (registration required, \
                       academic-use license) and `pip install nupack-*.whl` \
                       after ensuring python3 is on PATH"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = NupackInput::from_case_dir(&case.path)?;

        // Round-4 security: reject `output_basename = "../etc/passwd"`
        // and friends before the value flows into any path join.
        // Same pattern as the round-3 fix in bionetgen/iqtree/art/fasttree.
        valenx_core::adapter_helpers::validate_output_basename(
            &input.output_basename,
            "[bio.nupack].output_basename",
        )
        .map_err(|e| AdapterError::InvalidCase {
            case_path: case.path.join("case.toml"),
            reason: format!("{e}"),
        })?;

        fs::create_dir_all(workdir)?;

        // Stage the user-supplied Python script. Resolves against the
        // case directory — same convention as every other Python-script
        // bio adapter. `confined_join` rejects absolute paths and `..`
        // traversal so a malicious case bundle can't smuggle arbitrary
        // host files into the workdir.
        let source_script = confined_join(&case.path, &input.script)?;
        if !source_script.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.nupack].script `{}` not found (resolved {})",
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
                        "[bio.nupack].script path `{}` has no filename",
                        input.script.display()
                    ),
                })?;
        let dest_script = workdir.join(script_filename);
        if source_script != dest_script {
            fs::copy(&source_script, &dest_script)?;
        }

        // Stage the optional input file (typically a multi-strand
        // FASTA / target-structure spec). Resolves against the case
        // directory; staging filename mirrors source.
        let staged_input: Option<OsString> = if let Some(ref input_path) = input.input {
            let source_input = confined_join(&case.path, input_path)?;
            if !source_input.is_file() {
                return Err(AdapterError::InvalidCase {
                    case_path: case.path.join("case.toml"),
                    reason: format!(
                        "[bio.nupack].input `{}` not found (resolved {})",
                        input_path.display(),
                        source_input.display()
                    ),
                });
            }
            let fname = input_path
                .file_name()
                .ok_or_else(|| AdapterError::InvalidCase {
                    case_path: case.path.join("case.toml"),
                    reason: format!(
                        "[bio.nupack].input path `{}` has no filename",
                        input_path.display()
                    ),
                })?
                .to_os_string();
            let dest = workdir.join(&fname);
            if source_input != dest {
                fs::copy(&source_input, &dest)?;
            }
            Some(fname)
        } else {
            None
        };

        // Drop a flat `valenx_params.json` so the user's script can
        // read the parsed `[bio.nupack]` knobs without reparsing
        // case.toml. Hand-rolled to keep the dep list aligned with
        // the RDKit / RFdiffusion / DeepChem siblings — no serde_json
        // dependency just for a 4-key flat object.
        let input_field = match &staged_input {
            Some(f) => json_string(&f.to_string_lossy()),
            None => "null".to_string(),
        };
        let params_json = format!(
            "{{\n  \"input\": {input_field},\n  \"output_basename\": {basename},\n  \"temperature\": {temperature},\n  \"sodium\": {sodium}\n}}\n",
            basename = json_string(&input.output_basename),
            temperature = format_number(input.temperature),
            sodium = format_number(input.sodium),
        );
        valenx_core::io_caps::atomic_write_str(&workdir.join("valenx_params.json"), &params_json)?;

        // Resolve the Python binary. Same logic as every other
        // Python-script adapter: bare `python` / `python3` walks
        // PATH; absolute paths or pinned interpreters are honored
        // verbatim.
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

        // Stash the output_basename so collect() can filter for
        // `<basename>*` files without reparsing case.toml.
        let environment: Vec<(OsString, OsString)> = vec![(
            OsString::from("VALENX_NUPACK_OUTPUT_BASENAME"),
            OsString::from(&input.output_basename),
        )];

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment,
            // NUPACK analysis on a small tube finishes in seconds;
            // multi-target sequence design with hundreds of strands
            // can run for hours. 4 hours is a generous default that
            // covers the long tail without being absurd.
            estimated_runtime: Some(Duration::from_secs(4 * 60 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting NUPACK", |line| {
            let mut hint = subprocess::Hint::default();
            // Convention: the user-supplied script can emit a
            // sentinel line `[valenx] nupack done` to signal
            // completion before exit; lift to a 95% progress tick.
            if line.contains("[valenx] nupack done") {
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
        // Recover the basename so we can filter for `<basename>*`
        // files. Without it we'd surface every `.json` / `.csv` in
        // the workdir, including the staged `valenx_params.json`.
        let basename = job
            .environment
            .iter()
            .find(|(k, _)| k == "VALENX_NUPACK_OUTPUT_BASENAME")
            .map(|(_, v)| v.to_string_lossy().to_string());

        // Provenance: hash the staged script (the canonical
        // "this case is configured this way" input).
        let script_path = first_script_in_workdir(&job.workdir);
        let case_hash_input = script_path
            .clone()
            .unwrap_or_else(|| job.workdir.join("case.toml"));
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "NUPACK",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);
        let mut artefacts: Vec<Artifact> = Vec::new();

        // Walk the workdir top level. Filter to files starting with
        // `<output_basename>` so we don't surface the staged
        // `valenx_params.json` / staged input file as if they were
        // NUPACK results.
        let entries = match fs::read_dir(&job.workdir) {
            Ok(e) => e,
            Err(e) => {
                tracing::warn!(target: "valenx-nupack", ?e, "workdir read failed");
                return Ok(results);
            }
        };
        for entry in entries.flatten() {
            let path = entry.path();
            if !path.is_file() {
                continue;
            }
            let stem_or_filename = path.file_name().and_then(|s| s.to_str()).unwrap_or("");
            if let Some(b) = &basename {
                if !stem_or_filename.starts_with(b.as_str()) {
                    continue;
                }
            } else {
                // No basename recorded — collect() falling back to
                // empty results is preferable to surfacing every
                // file in the workdir.
                continue;
            }
            let ext = path
                .extension()
                .and_then(|s| s.to_str())
                .map(|s| s.to_ascii_lowercase());
            let (kind, label) = match ext.as_deref() {
                Some("json") => (ArtifactKind::Log, "NUPACK analysis log".to_string()),
                Some("npc") => (
                    ArtifactKind::Native,
                    "NUPACK proprietary output".to_string(),
                ),
                Some("csv") | Some("tsv") => {
                    (ArtifactKind::Tabular, "NUPACK tabular output".to_string())
                }
                _ => (ArtifactKind::Native, "NUPACK output".to_string()),
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
            ribbon_contributions: vec!["bio.nupack.analyze"],
        }
    }
}

/// Escape a string for embedding inside a JSON string literal.
/// Mirrors the helper used by RFdiffusion / DeepChem / ESMFold so we
/// don't have to pull in a serde_json dep just to emit a 4-key flat
/// object.
fn json_string(s: &str) -> String {
    let mut out = String::with_capacity(s.len() + 2);
    out.push('"');
    for c in s.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            '\u{0008}' => out.push_str("\\b"),
            '\u{000C}' => out.push_str("\\f"),
            c if (c as u32) < 0x20 => {
                out.push_str(&format!("\\u{:04x}", c as u32));
            }
            c => out.push(c),
        }
    }
    out.push('"');
    out
}

/// Render a finite f64 as a JSON number literal. Whole-number values
/// render with an explicit `.0` so the file always parses as a JSON
/// number (rather than an integer that some downstream tooling might
/// coerce back into Python's `int`).
fn format_number(value: f64) -> String {
    if value.fract() == 0.0 {
        format!("{value:.1}")
    } else {
        format!("{value}")
    }
}

/// Lift the staged Python script out of the workdir for provenance
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

/// Run `python -c "import nupack; print(nupack.__version__)"` and
/// parse a `semver::Version` out of stdout. Returns `None` on any
/// failure (interpreter unusable, NUPACK not importable, version
/// string malformed); `probe()` falls back to a warning in that case.
fn detect_nupack_version(python_binary: &Path) -> Option<Version> {
    let output = std::process::Command::new(python_binary)
        .arg("-c")
        .arg("import nupack; print(nupack.__version__)")
        .stdin(std::process::Stdio::null())
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::piped())
        .output()
        .ok()?;
    if !output.status.success() {
        return None;
    }
    let stdout = String::from_utf8(output.stdout).ok()?;
    let raw = valenx_core::adapter_helpers::extract_semver(&stdout)?;
    let dots = raw.chars().filter(|c| *c == '.').count();
    let normalised: String = match dots {
        0 => format!("{raw}.0.0"),
        1 => format!("{raw}.0"),
        _ => raw,
    };
    Version::parse(&normalised).ok()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn info_is_bio_domain() {
        let info = NupackAdapter::new().info();
        assert_eq!(info.id, "nupack");
        assert_eq!(info.physics, &[Physics::Bio]);
        // The license identifier must surface NUPACK's custom non-OSS
        // license rather than mislabel as MIT / BSD.
        assert_eq!(info.tool_license, "NUPACK-License");
        assert_eq!(info.display_name, "NUPACK");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = NupackAdapter::new().info();
        // NUPACK 4.0 is the modern Python-first rewrite; 5.0 reserves
        // room for an eventual major bump.
        assert_eq!(info.version_range.min_inclusive, Version::new(4, 0, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(5, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = NupackAdapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.nupack.analyze"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = NupackAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }

    #[test]
    fn probe_warning_mentions_academic() {
        // The license-flag warning is mandatory: NUPACK is non-OSS
        // academic-use, and we surface that on every successful
        // probe. The literal "academic" anchor is what downstream
        // tooling and license-aware filters key off — pin it.
        assert!(
            LICENSE_WARNING.contains("academic"),
            "probe warning must contain `academic` anchor; got: {LICENSE_WARNING}"
        );
    }
}
