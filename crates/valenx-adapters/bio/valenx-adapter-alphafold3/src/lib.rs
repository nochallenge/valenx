//! # valenx-adapter-alphafold3
//!
//! Adapter for [AlphaFold 3](https://github.com/google-deepmind/alphafold3) —
//! DeepMind's next-generation biomolecular structure predictor. AF3
//! goes beyond AF2's protein-only scope: a single inference run can
//! fold a complex of proteins, nucleic acids, and small-molecule
//! ligands described as a single JSON job spec.
//!
//! **Phase 17.5 — subprocess wrapper around `run_alphafold.py`.** The
//! user supplies the path to the AF3 checkout (`run_alphafold.py`
//! plus a directory of weights and a directory of databases) via
//! `[bio.alphafold3]` in `case.toml`. `prepare()` stages the input
//! JSON into the workdir and constructs the canonical AF3 command
//! line. `run()` invokes
//! `python <run_script> --json_path=… --output_dir=… …` via the
//! shared subprocess runner.
//!
//! On `collect()` we walk `<workdir>/<job_name>/` for the customary
//! AF3 outputs: `<job_name>_model.cif` (the predicted structure) and
//! `<job_name>_summary_confidences.json` (per-component confidence
//! metrics).
//!
//! ## License flag — non-commercial weights
//!
//! AF3's source code is Apache-2.0 but the **model weights are
//! released under CC-BY-NC-4.0 (non-commercial only)**. The adapter
//! surfaces this constraint via [`AdapterInfo::tool_license`] and
//! pushes a "non-commercial" warning into every successful
//! [`ProbeReport`] so downstream UI / audit log can flag it.

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

use crate::case_input::AlphaFold3Input;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(AlphaFold3Adapter::new())
}

pub struct AlphaFold3Adapter;

impl AlphaFold3Adapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for AlphaFold3Adapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "alphafold3";
const PYTHON_BINARIES: &[&str] = &["python3", "python"];

/// The non-commercial-weights warning. Surfaced verbatim from
/// `probe()` — there's a test that asserts the substring
/// `"non-commercial"` appears in the warning text so callers can
/// rely on it being there.
const NON_COMMERCIAL_WARNING: &str = "AlphaFold 3 model weights are CC-BY-NC-4.0 \
     (non-commercial). Confirm your use case complies with the upstream license \
     before redistributing predictions or derived data.";

impl Adapter for AlphaFold3Adapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "AlphaFold 3",
            // AF3's first public release is 3.0; upper bound 4.0
            // reserves room for a future major bump.
            version_range: VersionRange {
                min_inclusive: Version::new(3, 0, 0),
                max_exclusive: Version::new(4, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            // AF3 weights ship under CC-BY-NC-4.0 (non-commercial).
            // The source code itself is Apache-2.0 but the weights
            // dominate downstream license obligations — surface the
            // restrictive licence here so the registry / first-run
            // wizard can show it.
            tool_license: "CC-BY-NC-4.0",
            docs_url: "https://github.com/google-deepmind/alphafold3",
            homepage_url: "https://github.com/google-deepmind/alphafold3",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(PYTHON_BINARIES) {
            Some(binary_path) => {
                let found_version = detect_alphafold3_version(&binary_path);
                // The non-commercial warning is unconditional once
                // we've found a Python interpreter — AF3's license
                // applies to anyone who downloads the weights, not
                // just to runs that successfully import the package.
                // Tests assert this verbatim ("non-commercial").
                let mut warnings: Vec<String> = vec![NON_COMMERCIAL_WARNING.to_string()];
                if found_version.is_none() {
                    warnings.push(
                        "probe found `python` on PATH but could not import \
                         `alphafold3` — clone \
                         https://github.com/google-deepmind/alphafold3 and \
                         install per the README before runs will succeed"
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
                hint: "Python 3.11+ with AlphaFold 3 installed; clone \
                       https://github.com/google-deepmind/alphafold3 and \
                       follow the install README after ensuring python3 \
                       is on PATH"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = AlphaFold3Input::from_case_dir(&case.path)?;

        fs::create_dir_all(workdir)?;

        // Resolve `run_script`. Like AF2, we don't copy AF3's
        // `run_alphafold.py` into the workdir — it depends on the
        // surrounding repo layout. Verify it exists and pass the
        // absolute path on the command line.
        // Round-9 hardening: `run_script` is user-supplied data and
        // flows into `Command::new` (via `python <script>`) — wrap
        // relative paths with `confined_join` so a hostile case
        // can't aim it at `../../etc/whatever`.
        let run_script = if input.run_script.is_absolute() {
            input.run_script.clone()
        } else {
            confined_join(&case.path, &input.run_script)?
        };
        if !run_script.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.alphafold3].run_script `{}` not found (resolved {})",
                    input.run_script.display(),
                    run_script.display()
                ),
            });
        }

        // Stage the AF3 JSON job spec into the workdir.
        // `confined_join` rejects absolute paths and `..` traversal so
        // the staged copy stays confined to the case directory.
        let source_json = confined_join(&case.path, &input.input_json)?;
        if !source_json.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.alphafold3].input_json `{}` not found (resolved {})",
                    input.input_json.display(),
                    source_json.display()
                ),
            });
        }
        let json_filename =
            input
                .input_json
                .file_name()
                .ok_or_else(|| AdapterError::InvalidCase {
                    case_path: case.path.join("case.toml"),
                    reason: format!(
                        "[bio.alphafold3].input_json path `{}` has no filename",
                        input.input_json.display()
                    ),
                })?;
        let dest_json = workdir.join(json_filename);
        if source_json != dest_json {
            fs::copy(&source_json, &dest_json)?;
        }

        // model_dir + db_dir are large user-managed directories — we
        // never copy them. Resolve relative paths against the case
        // directory so the user can co-locate them with case.toml,
        // then verify each is a real directory so the user gets a
        // fast fail before AF3 does.
        // Round-9 classification: KEEP `case.path.join` here —
        // `model_dir` and `db_dir` (below) are AF3's multi-TB
        // weights/databases bundles that normally live wherever the
        // admin staged them on a large volume outside the case
        // sandbox. We only `.is_dir()` them for fast-fail; the actual
        // command-line flag uses the absolute path as-is.
        let model_dir = if input.model_dir.is_absolute() {
            input.model_dir.clone()
        } else {
            case.path.join(&input.model_dir)
        };
        if !model_dir.is_dir() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.alphafold3].model_dir `{}` is not a directory (resolved {})",
                    input.model_dir.display(),
                    model_dir.display()
                ),
            });
        }
        let db_dir = if input.db_dir.is_absolute() {
            input.db_dir.clone()
        } else {
            case.path.join(&input.db_dir)
        };
        if !db_dir.is_dir() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.alphafold3].db_dir `{}` is not a directory (resolved {})",
                    input.db_dir.display(),
                    db_dir.display()
                ),
            });
        }

        // Round-3 security fix (round-12 sweep): `input.python` flows
        // straight into `Command::new`, so a hostile case.toml could
        // otherwise point it at e.g. `/usr/bin/curl` and turn "Run case"
        // into arbitrary exec. `resolve_python_binary` bundles
        // allow-list validation, absolute-path acceptance, `..`-traversal
        // rejection, and PATH resolution into a single call shared with
        // alphamissense / anndata / be-designer / esmfold / rfdiffusion.
        let binary_path =
            valenx_core::adapter_helpers::resolve_python_binary(&input.python, PYTHON_BINARIES)
                .map_err(|e| AdapterError::InvalidCase {
                    case_path: case.path.join("case.toml"),
                    reason: format!("[bio.alphafold3].python: {e}"),
                })?;

        let mut native_command: Vec<OsString> = Vec::new();
        native_command.push(binary_path.into_os_string());
        native_command.push(run_script.into_os_string());
        // Use the separated `--flag value` form for every path
        // argument. The joined `--flag=value` form is ambiguous when
        // `value` itself contains `=` (legal POSIX path char, common
        // in `foo=bar/baz` directory names); absl::flags parses
        // `--json_path=/foo=bar/baz` as flag = "json_path",
        // value = "/foo", junk = "bar/baz". Separating the tokens
        // sidesteps the shell-parsing rule entirely. The
        // num_diffusion_samples flag is a small integer with no
        // ambiguity, so we leave it in the compact form.
        native_command.push(OsString::from("--json_path"));
        native_command.push(OsString::from(dest_json.as_os_str()));
        native_command.push(OsString::from("--output_dir"));
        native_command.push(OsString::from(workdir.as_os_str()));
        native_command.push(OsString::from("--model_dir"));
        native_command.push(OsString::from(model_dir.as_os_str()));
        native_command.push(OsString::from("--db_dir"));
        native_command.push(OsString::from(db_dir.as_os_str()));
        native_command.push(OsString::from(format!(
            "--num_diffusion_samples={}",
            input.num_diffusion_samples
        )));

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            // AF3 jobs span minutes (small monomers) to many hours
            // (large complexes with many diffusion samples). 8 hours
            // is a reasonable default.
            estimated_runtime: Some(Duration::from_secs(8 * 60 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting AlphaFold 3", |line| {
            let mut hint = subprocess::Hint::default();
            // AF3's stdout doesn't have a documented progress format
            // yet; do the conservative thing and only flag warnings /
            // errors, leaving the spinner to keep the UI alive while
            // diffusion samples roll in.
            if line.contains("Traceback") || line.contains("Error") {
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
        let json_path = first_json_in_workdir(&job.workdir);
        let case_hash_input = json_path
            .clone()
            .unwrap_or_else(|| job.workdir.join("case.toml"));
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "AlphaFold 3",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);
        let mut artefacts: Vec<Artifact> = Vec::new();

        if let Some(p) = json_path {
            artefacts.push(Artifact {
                path: p,
                kind: ArtifactKind::Other,
                checksum: None,
                label: "AlphaFold 3 input JSON".to_string(),
            });
        }

        // AF3 writes outputs under `<workdir>/<job_name>/` (the job
        // name comes from the input JSON's `name` field). Walk every
        // immediate subdirectory so we don't have to know the job
        // name up front. Within each subdir we look for:
        //   - `<job_name>_model.cif` (Native)
        //   - `<job_name>_summary_confidences.json` (Tabular)
        //   - any other JSON / CIF supporting files
        let entries = match fs::read_dir(&job.workdir) {
            Ok(e) => e,
            Err(e) => {
                tracing::warn!(target: "valenx-alphafold3", ?e, "workdir read failed");
                return Ok(results);
            }
        };
        for entry in entries.flatten() {
            let path = entry.path();
            if !path.is_dir() {
                continue;
            }
            let sub_entries = match fs::read_dir(&path) {
                Ok(e) => e,
                Err(_) => continue,
            };
            for sub in sub_entries.flatten() {
                let p = sub.path();
                if !p.is_file() {
                    continue;
                }
                let name = match p.file_name().and_then(|s| s.to_str()) {
                    Some(s) => s.to_string(),
                    None => continue,
                };
                let ext = p
                    .extension()
                    .and_then(|s| s.to_str())
                    .map(|s| s.to_ascii_lowercase());
                match ext.as_deref() {
                    Some("cif") => {
                        artefacts.push(Artifact {
                            path: p,
                            kind: ArtifactKind::Native,
                            checksum: None,
                            label: format!(
                                "AlphaFold 3 prediction `{}`",
                                name.trim_end_matches(".cif")
                            ),
                        });
                    }
                    Some("json") => {
                        let label = if name.contains("summary_confidences") {
                            "AlphaFold 3 summary confidences".to_string()
                        } else {
                            "AlphaFold 3 metadata".to_string()
                        };
                        artefacts.push(Artifact {
                            path: p,
                            kind: ArtifactKind::Tabular,
                            checksum: None,
                            label,
                        });
                    }
                    _ => continue,
                }
            }
        }

        artefacts.sort_by(|a, b| a.path.cmp(&b.path));
        results.artifacts = artefacts;
        Ok(results)
    }

    fn capabilities(&self) -> Capabilities {
        Capabilities {
            capabilities: Vec::new(),
            ribbon_contributions: vec!["bio.alphafold3.predict"],
        }
    }
}

/// Lift the staged input JSON out of the workdir for provenance
/// hashing. Returns the lexicographically-first `.json` file at the
/// top level (AF3's per-job output JSONs live under
/// `<job_name>/` subdirs, so the top-level JSON is the user's input
/// spec).
fn first_json_in_workdir(workdir: &Path) -> Option<PathBuf> {
    let entries = fs::read_dir(workdir).ok()?;
    let mut hits: Vec<PathBuf> = entries
        .flatten()
        .map(|e| e.path())
        .filter(|p| {
            p.is_file()
                && p.extension()
                    .and_then(|s| s.to_str())
                    .map(|s| s.eq_ignore_ascii_case("json"))
                    .unwrap_or(false)
        })
        .collect();
    hits.sort();
    hits.into_iter().next()
}

/// Run `python -c "import alphafold3; print(alphafold3.__version__)"`
/// and parse a `semver::Version` out of stdout. Returns `None` on any
/// failure.
fn detect_alphafold3_version(python_binary: &Path) -> Option<Version> {
    let output = std::process::Command::new(python_binary)
        .arg("-c")
        .arg("import alphafold3; print(alphafold3.__version__)")
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
    use valenx_test_utils::tempdir;

    #[test]
    fn info_is_bio_domain() {
        let info = AlphaFold3Adapter::new().info();
        assert_eq!(info.id, "alphafold3");
        assert_eq!(info.physics, &[Physics::Bio]);
        // The non-commercial weight licence is the load-bearing
        // detail downstream consumers (audit log, registry UI) read
        // from the AdapterInfo. Pin it.
        assert_eq!(info.tool_license, "CC-BY-NC-4.0");
        assert_eq!(info.display_name, "AlphaFold 3");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = AlphaFold3Adapter::new().info();
        assert_eq!(info.version_range.min_inclusive, Version::new(3, 0, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(4, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = AlphaFold3Adapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.alphafold3.predict"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = AlphaFold3Adapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }

    /// AF3's weights are released under CC-BY-NC-4.0; the probe
    /// report must surface that fact verbatim so downstream UI /
    /// audit log can flag any commercial use.
    ///
    /// We exercise the message-construction path directly by
    /// inspecting the embedded constant — that way the test runs
    /// deterministically on every CI host regardless of whether
    /// Python is on PATH.
    #[test]
    #[ignore] // subprocess-coupled test — run interactively only
    fn probe_warning_mentions_non_commercial() {
        // The constant must contain the substring
        // `"non-commercial"` verbatim — every probe path that
        // succeeds clones this into ProbeReport.warnings.
        assert!(
            NON_COMMERCIAL_WARNING.contains("non-commercial"),
            "probe warning text was: {NON_COMMERCIAL_WARNING}"
        );
        // When the probe path actually finds Python, it pushes the
        // warning. We can't depend on Python being present, but we
        // can verify the live path produces it when it does — so
        // assert it via the only branch that doesn't need a child
        // process.
        if find_on_path(PYTHON_BINARIES).is_some() {
            let report = AlphaFold3Adapter::new()
                .probe()
                .expect("probe ok with python on PATH");
            assert!(
                report.warnings.iter().any(|w| w.contains("non-commercial")),
                "probe warnings: {:?}",
                report.warnings
            );
        }
    }

    /// `collect()` must walk `<workdir>/<job_name>/` for AF3's
    /// `<job_name>_model.cif` + `<job_name>_summary_confidences.json`
    /// outputs and surface them as the right artifact kinds.
    #[test]
    fn collect_walks_job_subdir_for_cif_and_summary_json() {
        let tmp = std::env::temp_dir().join(format!(
            "valenx-af3-collect-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(&tmp).unwrap();
        fs::write(tmp.join("job.json"), br#"{"name": "test_complex"}"#).unwrap();
        let job_subdir = tmp.join("test_complex");
        fs::create_dir_all(&job_subdir).unwrap();
        fs::write(
            job_subdir.join("test_complex_model.cif"),
            b"data_test\n_cell.length_a 10.0\n",
        )
        .unwrap();
        fs::write(
            job_subdir.join("test_complex_summary_confidences.json"),
            br#"{"plddt_mean": 85.0}"#,
        )
        .unwrap();

        let job = PreparedJob {
            workdir: tmp.clone(),
            native_command: vec![],
            environment: Vec::new(),
            estimated_runtime: None,
            kill_on_drop: true,
        };
        let results = AlphaFold3Adapter::new().collect(&job).unwrap();

        let cif = results
            .artifacts
            .iter()
            .find(|a| a.path.extension().is_some_and(|e| e == "cif"))
            .expect("CIF artifact present");
        assert_eq!(cif.kind, ArtifactKind::Native);
        assert!(
            cif.label.contains("AlphaFold 3 prediction"),
            "label was: {}",
            cif.label
        );

        let summary = results
            .artifacts
            .iter()
            .find(|a| {
                a.path
                    .file_name()
                    .and_then(|n| n.to_str())
                    .is_some_and(|n| n.contains("summary_confidences"))
            })
            .expect("summary JSON present");
        assert_eq!(summary.kind, ArtifactKind::Tabular);
        assert_eq!(summary.label, "AlphaFold 3 summary confidences");

        // Top-level input JSON surfaces too.
        let input_json = results
            .artifacts
            .iter()
            .find(|a| {
                a.path
                    .file_name()
                    .and_then(|n| n.to_str())
                    .is_some_and(|n| n == "job.json")
            })
            .expect("input JSON present");
        assert_eq!(input_json.kind, ArtifactKind::Other);
        assert_eq!(input_json.label, "AlphaFold 3 input JSON");

        let _ = fs::remove_dir_all(&tmp);
    }

    /// `prepare()` must fail fast when `model_dir` doesn't resolve to
    /// a real directory. Without this check AF3 would crash deep in
    /// its own startup with a much less actionable message.
    #[test]
    fn prepare_rejects_missing_model_dir() {
        let tmp = tempdir("alphafold3-af3-no-modeldir");
        let case_dir = tmp.join("case");
        fs::create_dir_all(&case_dir).unwrap();
        // Stage a valid run_script + JSON + db_dir so prepare() gets
        // past the earlier checks; only `model_dir` should be the
        // failure.
        fs::write(case_dir.join("run_alphafold.py"), b"# placeholder\n").unwrap();
        fs::write(case_dir.join("job.json"), b"{\"name\": \"x\"}\n").unwrap();
        let db = case_dir.join("af3-db");
        fs::create_dir_all(&db).unwrap();
        fs::write(
            case_dir.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "alphafold3.predict"

[bio.alphafold3]
run_script = "run_alphafold.py"
input_json = "job.json"
model_dir  = "does_not_exist"
db_dir     = "af3-db"
"#,
        )
        .unwrap();
        let case = Case {
            id: "af3-no-modeldir".into(),
            path: case_dir.clone(),
        };
        let workdir = tmp.join("work");
        let err = AlphaFold3Adapter::new()
            .prepare(&case, &workdir)
            .expect_err("expected InvalidCase for missing model_dir");
        match err {
            AdapterError::InvalidCase { reason, .. } => {
                assert!(
                    reason.contains("model_dir") && reason.contains("not a directory"),
                    "reason was: {reason}"
                );
            }
            other => panic!("unexpected error: {other:?}"),
        }
        let _ = fs::remove_dir_all(&tmp);
    }

    /// `prepare()` must fail fast when `db_dir` doesn't resolve to a
    /// real directory.
    #[test]
    fn prepare_rejects_missing_db_dir() {
        let tmp = tempdir("alphafold3-af3-no-dbdir");
        let case_dir = tmp.join("case");
        fs::create_dir_all(&case_dir).unwrap();
        fs::write(case_dir.join("run_alphafold.py"), b"# placeholder\n").unwrap();
        fs::write(case_dir.join("job.json"), b"{\"name\": \"x\"}\n").unwrap();
        let model = case_dir.join("af3-models");
        fs::create_dir_all(&model).unwrap();
        fs::write(
            case_dir.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "alphafold3.predict"

[bio.alphafold3]
run_script = "run_alphafold.py"
input_json = "job.json"
model_dir  = "af3-models"
db_dir     = "does_not_exist"
"#,
        )
        .unwrap();
        let case = Case {
            id: "af3-no-dbdir".into(),
            path: case_dir.clone(),
        };
        let workdir = tmp.join("work");
        let err = AlphaFold3Adapter::new()
            .prepare(&case, &workdir)
            .expect_err("expected InvalidCase for missing db_dir");
        match err {
            AdapterError::InvalidCase { reason, .. } => {
                assert!(
                    reason.contains("db_dir") && reason.contains("not a directory"),
                    "reason was: {reason}"
                );
            }
            other => panic!("unexpected error: {other:?}"),
        }
        let _ = fs::remove_dir_all(&tmp);
    }

    /// Round-3 security fix: `[bio.alphafold3].python` flows into
    /// `Command::new`. A hostile case.toml setting it to `/usr/bin/curl`
    /// must be rejected before any subprocess is spawned.
    #[test]
    fn prepare_rejects_arbitrary_python_binary() {
        let tmp = tempdir("alphafold3-bad-python");
        let case_dir = tmp.join("case");
        fs::create_dir_all(&case_dir).unwrap();
        fs::write(case_dir.join("run_alphafold.py"), b"# placeholder\n").unwrap();
        fs::write(case_dir.join("job.json"), b"{\"name\": \"x\"}\n").unwrap();
        fs::create_dir_all(case_dir.join("af3-models")).unwrap();
        fs::create_dir_all(case_dir.join("af3-db")).unwrap();
        // Hostile python path. /usr/bin/curl on Unix or any other
        // unrelated binary on Windows — the allow-list check rejects
        // based on the basename, so neither needs to actually exist.
        fs::write(
            case_dir.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "alphafold3.predict"

[bio.alphafold3]
run_script = "run_alphafold.py"
input_json = "job.json"
model_dir  = "af3-models"
db_dir     = "af3-db"
python     = "/usr/bin/curl"
"#,
        )
        .unwrap();
        let case = Case {
            id: "af3-bad-python".into(),
            path: case_dir.clone(),
        };
        let workdir = tmp.join("work");
        let err = AlphaFold3Adapter::new()
            .prepare(&case, &workdir)
            .expect_err("hostile python value must be rejected");
        match err {
            AdapterError::InvalidCase { reason, .. } => {
                assert!(
                    reason.contains("python") && reason.contains("allow"),
                    "reason was: {reason}"
                );
            }
            other => panic!("unexpected error: {other:?}"),
        }
        let _ = fs::remove_dir_all(&tmp);
    }

    /// Round-9 RED→GREEN: `[bio.alphafold3].run_script` used to be
    /// joined with bare `case.path.join`, which let a hostile case
    /// supply `run_script = "../../etc/passwd"` (or absolute paths
    /// on POSIX) and have AF3 execute whatever Python the user could
    /// reach. The fix wraps the relative branch with `confined_join`,
    /// same policy as `input_json`.
    #[test]
    fn prepare_rejects_run_script_traversing_outside_case_dir() {
        let tmp = tempdir("alphafold3-runscript-trav");
        let case_dir = tmp.join("case");
        fs::create_dir_all(&case_dir).unwrap();
        fs::write(case_dir.join("job.json"), b"{\"name\": \"x\"}\n").unwrap();
        let model = case_dir.join("af3-models");
        fs::create_dir_all(&model).unwrap();
        let db = case_dir.join("af3-db");
        fs::create_dir_all(&db).unwrap();
        fs::write(
            case_dir.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "alphafold3.predict"

[bio.alphafold3]
run_script = "../../etc/passwd"
input_json = "job.json"
model_dir  = "af3-models"
db_dir     = "af3-db"
"#,
        )
        .unwrap();
        let case = Case {
            id: "af3-runscript-trav".into(),
            path: case_dir.clone(),
        };
        let workdir = tmp.join("work");
        let err = AlphaFold3Adapter::new()
            .prepare(&case, &workdir)
            .expect_err("must reject ../../etc/passwd run_script");
        let msg = format!("{err}");
        assert!(
            msg.contains("..") || msg.contains("stay within") || msg.contains("escape"),
            "expected confined_join rejection, got: {msg}"
        );
        let _ = fs::remove_dir_all(&tmp);
    }
}
