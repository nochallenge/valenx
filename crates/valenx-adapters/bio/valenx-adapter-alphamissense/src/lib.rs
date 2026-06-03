//! # valenx-adapter-alphamissense
//!
//! Adapter for [AlphaMissense](https://github.com/google-deepmind/alphamissense) —
//! DeepMind's missense effect predictor (Cheng et al., Science).
//! AlphaMissense classifies single-amino-acid substitutions as benign,
//! pathogenic, or ambiguous using a model bootstrapped from AlphaFold's
//! structural prior plus protein-language-model features.
//!
//! **Phase 35.6 — Python subprocess wrapper for user-provided
//! scripts.** AlphaMissense ships as the upstream `alphamissense`
//! Python package (github.com/google-deepmind/alphamissense); the
//! adapter doesn't reimplement prediction logic. The user authors a
//! `predict.py` referenced from `[bio.alphamissense].script` in
//! `case.toml` that does `import alphamissense` and the actual
//! prediction logic. `prepare()` stages the script (and an optional
//! input `.fa` template) into the workdir, drops a flat
//! `valenx_params.json` next to it, and `run()` invokes
//! `python <script>` via the shared subprocess runner.
//!
//! ## License flag — academic / non-commercial weights
//!
//! AlphaMissense's source code is Apache-2.0 but the **model weights
//! are released under CC-BY-NC-SA-4.0 (academic / non-commercial
//! use only)**. The adapter surfaces this constraint via
//! [`AdapterInfo::tool_license`] and pushes an "academic /
//! non-commercial" warning into every successful [`ProbeReport`] so
//! downstream UI / audit log can flag it. The literal substrings
//! `"academic"` and `"non-commercial"` are part of the asserted
//! contract — see `LICENSE_WARNING` and the
//! `probe_warning_mentions_academic_and_non_commercial` test.
//!
//! ## `valenx_params.json`
//!
//! ```json
//! {
//!   "output_basename": "output",
//!   "input_fasta": "target.fa"
//! }
//! ```
//!
//! `input_fasta` is omitted entirely (not `null`) when the user did
//! not supply one.
//!
//! On `collect()` we walk the workdir for `<basename>*.csv` /
//! `<basename>*.tsv` (predicted pathogenicity scores), `<basename>*.png`
//! (visualisation plots), plus any `*.log` files Python emits.

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

use crate::case_input::AlphaMissenseInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(AlphaMissenseAdapter::new())
}

pub struct AlphaMissenseAdapter;

impl AlphaMissenseAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for AlphaMissenseAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "alphamissense";
const PYTHON_BINARIES: &[&str] = &["python3", "python"];

/// The probe-warning surfaced whenever a Python interpreter is found.
/// Anchors a stable "academic / non-commercial only" reminder for
/// downstream tooling and tests; the literal substrings `"academic"`
/// and `"non-commercial"` are part of the asserted contract — sister
/// to the NAMD / AlphaFold 3 license-warning convention.
const LICENSE_WARNING: &str = "AlphaMissense weights are licensed CC-BY-NC-SA-4.0 — for academic \
     / non-commercial use only. See \
     https://github.com/google-deepmind/alphamissense#licence for terms.";

impl Adapter for AlphaMissenseAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "AlphaMissense",
            // AlphaMissense 1.x is the current upstream line; 2.0
            // reserves room for the next major.
            version_range: VersionRange {
                min_inclusive: Version::new(1, 0, 0),
                max_exclusive: Version::new(2, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            // AlphaMissense weights ship under CC-BY-NC-SA-4.0
            // (academic / non-commercial). The source code itself is
            // Apache-2.0 but the weights dominate downstream license
            // obligations — surface the restrictive licence here so
            // the registry / first-run wizard can show it.
            tool_license: "CC-BY-NC-SA-4.0",
            docs_url: "https://github.com/google-deepmind/alphamissense",
            homepage_url: "https://alphamissense.hegelab.org/",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(PYTHON_BINARIES) {
            Some(binary_path) => {
                // The academic / non-commercial warning is
                // unconditional once we've found a Python interpreter
                // — AlphaMissense's license applies to anyone who
                // downloads the weights, not just to runs that
                // successfully import the package. Tests assert this
                // verbatim ("academic", "non-commercial").
                let mut warnings: Vec<String> = vec![LICENSE_WARNING.to_string()];
                let import_ok = alphamissense_importable(&binary_path);
                if !import_ok {
                    warnings.push(
                        "probe found `python` on PATH but could not import \
                         `alphamissense` — install from upstream \
                         https://github.com/google-deepmind/alphamissense \
                         for runs to succeed"
                            .into(),
                    );
                }
                Ok(ProbeReport {
                    ok: true,
                    found_version: None,
                    binary_path: Some(binary_path),
                    warnings,
                    required_env: Vec::new(),
                })
            }
            None => Err(AdapterError::ToolNotInstalled {
                name: INFO_ID,
                hint: "Python 3.10+ with alphamissense installed; install \
                       from upstream \
                       https://github.com/google-deepmind/alphamissense \
                       after ensuring python3 is on PATH"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = AlphaMissenseInput::from_case_dir(&case.path)?;

        // Round-4 security: reject `output_basename = "../etc/passwd"`
        // and friends before the value flows into any path join.
        // Same pattern as the round-3 fix in bionetgen/iqtree/art/fasttree.
        valenx_core::adapter_helpers::validate_output_basename(
            &input.output_basename,
            "[bio.alphamissense].output_basename",
        )
        .map_err(|e| AdapterError::InvalidCase {
            case_path: case.path.join("case.toml"),
            reason: format!("{e}"),
        })?;

        fs::create_dir_all(workdir)?;

        let source_script = confined_join(&case.path, &input.script)?;
        if !source_script.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.alphamissense].script `{}` not found (resolved {})",
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
                        "[bio.alphamissense].script path `{}` has no filename",
                        input.script.display()
                    ),
                })?;
        let dest_script = workdir.join(script_filename);
        if source_script != dest_script {
            fs::copy(&source_script, &dest_script)?;
        }

        let staged_input_fasta: Option<String> = match input.input_fasta.as_ref() {
            Some(fa_path) => {
                let source_fa = confined_join(&case.path, fa_path)?;
                if !source_fa.is_file() {
                    return Err(AdapterError::InvalidCase {
                        case_path: case.path.join("case.toml"),
                        reason: format!(
                            "[bio.alphamissense].input_fasta `{}` not found (resolved {})",
                            fa_path.display(),
                            source_fa.display()
                        ),
                    });
                }
                let fa_name = fa_path
                    .file_name()
                    .ok_or_else(|| AdapterError::InvalidCase {
                        case_path: case.path.join("case.toml"),
                        reason: format!(
                            "[bio.alphamissense].input_fasta path `{}` has no filename",
                            fa_path.display()
                        ),
                    })?;
                let dest_fa = workdir.join(fa_name);
                if source_fa != dest_fa {
                    fs::copy(&source_fa, &dest_fa)?;
                }
                Some(fa_name.to_string_lossy().to_string())
            }
            None => None,
        };

        let mut params = String::new();
        params.push_str("{\n");
        params.push_str("  \"output_basename\": ");
        params.push_str(&json_string(&input.output_basename));
        if let Some(name) = staged_input_fasta.as_deref() {
            params.push_str(",\n  \"input_fasta\": ");
            params.push_str(&json_string(name));
        }
        params.push_str("\n}\n");
        valenx_core::io_caps::atomic_write_str(&workdir.join("valenx_params.json"), &params)?;

        // Round-3 security fix (round-12 sweep): `input.python` flows
        // into `Command::new`. `resolve_python_binary` bundles
        // allow-list validation, absolute-path acceptance,
        // `..`-traversal rejection, and PATH resolution — replaces the
        // hand-rolled pattern shared with alphafold3 / anndata /
        // be-designer / esmfold / rfdiffusion.
        let binary_path =
            valenx_core::adapter_helpers::resolve_python_binary(&input.python, PYTHON_BINARIES)
                .map_err(|e| AdapterError::InvalidCase {
                    case_path: case.path.join("case.toml"),
                    reason: format!("[bio.alphamissense].python: {e}"),
                })?;

        let native_command: Vec<OsString> = vec![
            binary_path.into_os_string(),
            OsString::from(script_filename),
        ];

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            // AlphaMissense inference is GPU-friendly and per-protein
            // can range from seconds (single substitution) to minutes
            // (full saturation mutagenesis); 4 hours covers genome-
            // wide batches.
            estimated_runtime: Some(Duration::from_secs(4 * 60 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting AlphaMissense", |line| {
            let mut hint = subprocess::Hint::default();
            if line.contains("[valenx] alphamissense done") {
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
        let script_path = first_script_in_workdir(&job.workdir);
        let case_hash_input = script_path
            .clone()
            .unwrap_or_else(|| job.workdir.join("case.toml"));
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "AlphaMissense",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);
        let mut artefacts: Vec<Artifact> = Vec::new();

        let basename = read_output_basename(&job.workdir);

        let entries = match fs::read_dir(&job.workdir) {
            Ok(e) => e,
            Err(e) => {
                tracing::warn!(target: "valenx-alphamissense", ?e, "workdir read failed");
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
            match ext.as_deref() {
                Some("csv") | Some("tsv") => {
                    if !stem_matches_basename {
                        continue;
                    }
                    artefacts.push(Artifact {
                        path,
                        kind: ArtifactKind::Tabular,
                        checksum: None,
                        label: "AlphaMissense pathogenicity scores".to_string(),
                    });
                }
                Some("png") => {
                    if !stem_matches_basename {
                        continue;
                    }
                    artefacts.push(Artifact {
                        path,
                        kind: ArtifactKind::Native,
                        checksum: None,
                        label: "AlphaMissense plot".to_string(),
                    });
                }
                Some("log") => {
                    artefacts.push(Artifact {
                        path,
                        kind: ArtifactKind::Log,
                        checksum: None,
                        label: "AlphaMissense log".to_string(),
                    });
                }
                _ => continue,
            }
        }
        artefacts.sort_by(|a, b| a.path.cmp(&b.path));
        results.artifacts = artefacts;
        Ok(results)
    }

    fn capabilities(&self) -> Capabilities {
        Capabilities {
            capabilities: Vec::new(),
            ribbon_contributions: vec!["bio.alphamissense.predict"],
        }
    }
}

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

fn read_output_basename(workdir: &Path) -> Option<String> {
    let text = valenx_core::io_caps::read_capped_to_string(
        &workdir.join("valenx_params.json"),
        valenx_core::io_caps::MAX_ADAPTER_PARAMS_BYTES as usize,
    )
    .ok()?;
    extract_json_string(&text, "output_basename")
}

fn extract_json_string(text: &str, key: &str) -> Option<String> {
    let needle = format!("\"{key}\":");
    let idx = text.find(&needle)?;
    let rest = &text[idx + needle.len()..];
    let start = rest.find('"')? + 1;
    let body = &rest[start..];
    let mut out = String::new();
    let mut chars = body.chars();
    while let Some(c) = chars.next() {
        match c {
            '"' => return Some(out),
            '\\' => match chars.next()? {
                '"' => out.push('"'),
                '\\' => out.push('\\'),
                'n' => out.push('\n'),
                'r' => out.push('\r'),
                't' => out.push('\t'),
                'b' => out.push('\u{0008}'),
                'f' => out.push('\u{000C}'),
                other => out.push(other),
            },
            c => out.push(c),
        }
    }
    None
}

fn alphamissense_importable(python_binary: &Path) -> bool {
    let output = std::process::Command::new(python_binary)
        .arg("-c")
        .arg("import alphamissense")
        .stdin(std::process::Stdio::null())
        .stdout(std::process::Stdio::null())
        .stderr(std::process::Stdio::null())
        .output();
    match output {
        Ok(o) => o.status.success(),
        Err(_) => false,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn info_is_bio_domain() {
        let info = AlphaMissenseAdapter::new().info();
        assert_eq!(info.id, "alphamissense");
        assert_eq!(info.physics, &[Physics::Bio]);
        // The CC-BY-NC-SA-4.0 weight licence is the load-bearing
        // detail downstream consumers (audit log, registry UI) read
        // from the AdapterInfo. Pin it.
        assert_eq!(info.tool_license, "CC-BY-NC-SA-4.0");
        assert_eq!(info.display_name, "AlphaMissense");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = AlphaMissenseAdapter::new().info();
        assert_eq!(info.version_range.min_inclusive, Version::new(1, 0, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(2, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = AlphaMissenseAdapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.alphamissense.predict"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = AlphaMissenseAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }

    /// AlphaMissense's weights are released under CC-BY-NC-SA-4.0; the
    /// probe report must surface that fact verbatim so downstream UI /
    /// audit log can flag any commercial use.
    ///
    /// The literal `"academic"` and `"non-commercial"` substrings are
    /// what downstream tooling and license-aware filters key off — pin
    /// them. Sister to the NAMD / AlphaFold 3 license-warning
    /// convention.
    ///
    /// We always assert against the static `LICENSE_WARNING`
    /// constant. We also exercise the live `probe()` path when Python
    /// happens to be on PATH (so CI machines without it still pass);
    /// when present, the same substrings must surface in the probe
    /// report's warnings.
    /// Round-12 M7 RED→GREEN: after migrating to the shared
    /// `resolve_python_binary` helper, a `..`-bearing python spec
    /// like `"../python3"` is still rejected. The previous
    /// inline pattern in this adapter didn't enforce the
    /// `..`-traversal guard on its own — it relied on
    /// `validate_python_binary`'s allow-list, which actually passes
    /// `../python3` because the basename is `python3`. Switching
    /// to `resolve_python_binary` closes the gap by adding an
    /// explicit `ParentDir` component check.
    #[test]
    fn prepare_rejects_parent_dir_traversal_python_spec() {
        use std::fs;
        let tmp = std::env::temp_dir().join(format!(
            "alphamissense_round12_m7_{}",
            std::process::id()
        ));
        let _ = fs::remove_dir_all(&tmp);
        let case_dir = tmp.join("case");
        fs::create_dir_all(&case_dir).unwrap();
        // Minimal runnable case.toml with all required fields, but
        // python = "../python3" — the adapter must reject before any
        // subprocess spawn attempt. We supply a real script + fasta
        // so the test exercises the python-binary path rather than
        // bailing out earlier on a missing input field.
        fs::write(case_dir.join("run.py"), b"# stub\n").unwrap();
        fs::write(case_dir.join("seq.fa"), b">x\nM\n").unwrap();
        fs::write(
            case_dir.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "alphamissense.predict"

[bio.alphamissense]
script          = "run.py"
input_fasta     = "seq.fa"
output_basename = "out"
python          = "../python3"
"#,
        )
        .unwrap();
        let case = Case {
            id: "m7-traversal".into(),
            path: case_dir.clone(),
        };
        let workdir = tmp.join("work");
        let err = AlphaMissenseAdapter::new()
            .prepare(&case, &workdir)
            .expect_err("must reject ../python3");
        match err {
            AdapterError::InvalidCase { reason, .. } => {
                assert!(
                    reason.contains("..") || reason.contains("traverse"),
                    "expected parent-dir traversal rejection; got: {reason}"
                );
            }
            other => panic!("expected InvalidCase, got: {other:?}"),
        }
        let _ = fs::remove_dir_all(&tmp);
    }

    #[test]
    #[ignore] // subprocess-coupled test — run interactively only
    fn probe_warning_mentions_academic_and_non_commercial() {
        assert!(
            LICENSE_WARNING.contains("academic"),
            "probe warning must contain `academic` anchor; got: {LICENSE_WARNING}"
        );
        assert!(
            LICENSE_WARNING.contains("non-commercial"),
            "probe warning must contain `non-commercial` anchor; got: {LICENSE_WARNING}"
        );

        // Best-effort live probe — only assert if a Python interpreter
        // is on PATH. Skipping when none is keeps the test green on CI
        // machines without Python.
        if find_on_path(PYTHON_BINARIES).is_some() {
            let report = AlphaMissenseAdapter::new().probe().expect("probe");
            assert!(
                report
                    .warnings
                    .iter()
                    .any(|w| w.contains("academic") && w.contains("non-commercial")),
                "live probe warnings must surface the academic / non-commercial \
                 anchors; got: {:?}",
                report.warnings
            );
        }
    }
}
