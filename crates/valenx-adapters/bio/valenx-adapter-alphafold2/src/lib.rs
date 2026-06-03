//! # valenx-adapter-alphafold2
//!
//! Adapter for [AlphaFold 2](https://github.com/google-deepmind/alphafold) —
//! DeepMind's flagship protein structure predictor. AF2 is heavier
//! than ESMFold / OpenFold: it depends on a large reference database
//! (~2 TB) of MSA + template structures the user provides via the
//! `data_dir` flag.
//!
//! **Phase 17.5 — subprocess wrapper around `run_alphafold.py`.** The
//! user supplies the path to AF2's checkout (`run_alphafold.py` plus
//! the matching `data_dir`) via `[bio.alphafold2]` in `case.toml`.
//! `prepare()` stages the FASTA query into the workdir and
//! constructs the canonical AF2 command line. `run()` invokes
//! `python <run_script> --fasta_paths=… --output_dir=… …` via the
//! shared subprocess runner.
//!
//! On `collect()` we walk `<workdir>/<query_name>/` for the
//! customary AF2 outputs: `ranked_*.pdb` (the predicted models in
//! confidence order) and `ranking_debug.json` (the per-model
//! ranking metadata). PDBs are parsed via
//! [`valenx_bio::format::pdb::read`] and surfaced as typed
//! [`ArtifactKind::Native`] artifacts.

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

use crate::case_input::AlphaFold2Input;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(AlphaFold2Adapter::new())
}

pub struct AlphaFold2Adapter;

impl AlphaFold2Adapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for AlphaFold2Adapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "alphafold2";
const PYTHON_BINARIES: &[&str] = &["python3", "python"];

impl Adapter for AlphaFold2Adapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "AlphaFold 2",
            // AF2 2.3 is the first release with the multimer v3
            // weights surface we lean on. Upper bound 3.0 reserves
            // room for AF3, which lives in its own adapter.
            version_range: VersionRange {
                min_inclusive: Version::new(2, 3, 0),
                max_exclusive: Version::new(3, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            tool_license: "Apache-2.0",
            docs_url: "https://github.com/google-deepmind/alphafold",
            homepage_url: "https://deepmind.google/",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(PYTHON_BINARIES) {
            Some(binary_path) => {
                let found_version = detect_alphafold_version(&binary_path);
                let mut warnings = Vec::new();
                if found_version.is_none() {
                    warnings.push(
                        "probe found `python` on PATH but could not import \
                         `alphafold` — clone \
                         https://github.com/google-deepmind/alphafold and \
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
                hint: "Python 3.9+ with AlphaFold 2 installed; clone \
                       https://github.com/google-deepmind/alphafold and \
                       follow the install README after ensuring python3 \
                       is on PATH"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = AlphaFold2Input::from_case_dir(&case.path)?;

        fs::create_dir_all(workdir)?;

        // Resolve `run_script`. Unlike the script-stage adapters we
        // do NOT copy AF2's `run_alphafold.py` into the workdir —
        // it depends on the surrounding repo layout. Just verify it
        // exists and pass an absolute path on the command line.
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
                    "[bio.alphafold2].run_script `{}` not found (resolved {})",
                    input.run_script.display(),
                    run_script.display()
                ),
            });
        }

        // Stage the FASTA query into the workdir so it sits next to
        // the prediction outputs. `confined_join` rejects absolute
        // paths and `..` traversal so the staged copy stays confined
        // to the case directory.
        let source_fasta = confined_join(&case.path, &input.query_fasta)?;
        if !source_fasta.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.alphafold2].query_fasta `{}` not found (resolved {})",
                    input.query_fasta.display(),
                    source_fasta.display()
                ),
            });
        }
        let fasta_filename =
            input
                .query_fasta
                .file_name()
                .ok_or_else(|| AdapterError::InvalidCase {
                    case_path: case.path.join("case.toml"),
                    reason: format!(
                        "[bio.alphafold2].query_fasta path `{}` has no filename",
                        input.query_fasta.display()
                    ),
                })?;
        let dest_fasta = workdir.join(fasta_filename);
        if source_fasta != dest_fasta {
            fs::copy(&source_fasta, &dest_fasta)?;
        }

        // `data_dir` is large (~2 TB) and we never copy it; just
        // verify it exists. Surface a probe-style InvalidCase if it
        // doesn't so the user gets a fast fail before AF2 does.
        // Round-9 classification: KEEP `case.path.join` here —
        // `data_dir` is the AF2 reference databases bundle (UniRef,
        // PDB70, BFD, etc.) that lives wherever the admin staged it
        // on a multi-TB volume; expecting it to sit inside the case
        // sandbox would be wrong, and we never feed it to a subprocess
        // as an argument that could traverse — it's only `.is_dir()`'d
        // for early validation. The actual command-line flag uses the
        // absolute path as-is.
        let data_dir = if input.data_dir.is_absolute() {
            input.data_dir.clone()
        } else {
            case.path.join(&input.data_dir)
        };
        if !data_dir.is_dir() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.alphafold2].data_dir `{}` is not a directory (resolved {})",
                    input.data_dir.display(),
                    data_dir.display()
                ),
            });
        }

        // Resolve the Python binary.
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

        // Build the `run_alphafold.py` flag set. Pass arguments as
        // OsStrings so paths with non-UTF8 characters survive
        // round-tripping through the executor.
        let mut native_command: Vec<OsString> = Vec::new();
        native_command.push(binary_path.into_os_string());
        native_command.push(run_script.into_os_string());
        // Round-14 L11 (mirrors AF3's round-12 fix): use the
        // separated `--flag value` form for every path argument. The
        // joined `--flag=value` form is ambiguous when `value` itself
        // contains `=` (legal POSIX path char, common in `foo=bar/baz`
        // directory names); absl::flags parses
        // `--fasta_paths=/foo=bar/baz` as flag = "fasta_paths",
        // value = "/foo", junk = "bar/baz". Separating the tokens
        // sidesteps the shell-parsing rule entirely. The non-path
        // flags (max_template_date, model_preset) stay in the
        // compact form — they're short literal values where the
        // `=` ambiguity doesn't apply.
        native_command.push(OsString::from("--fasta_paths"));
        native_command.push(OsString::from(dest_fasta.as_os_str()));
        native_command.push(OsString::from("--output_dir"));
        native_command.push(OsString::from(workdir.as_os_str()));
        native_command.push(OsString::from("--data_dir"));
        native_command.push(OsString::from(data_dir.as_os_str()));
        native_command.push(OsString::from(format!(
            "--max_template_date={}",
            input.max_template_date
        )));
        native_command.push(OsString::from(format!(
            "--model_preset={}",
            input.model_preset
        )));

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            // AF2 runtime is dominated by MSA generation + folding;
            // 8 hours is a generous default for a single monomer on
            // consumer hardware.
            estimated_runtime: Some(Duration::from_secs(8 * 60 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting AlphaFold 2", |line| {
            let mut hint = subprocess::Hint::default();
            // AF2's stdout marches through coarse phases (jackhmmer,
            // hhsearch, model inference, relax). Map a few well-known
            // log lines to progress ticks so the UI shows movement.
            if line.contains("Searching") || line.contains("jackhmmer") {
                hint.progress = Some((30.0, line.to_string()));
            } else if line.contains("Predicting") || line.contains("Running model") {
                hint.progress = Some((60.0, line.to_string()));
            } else if line.contains("Relaxation") || line.contains("relaxed") {
                hint.progress = Some((90.0, line.to_string()));
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
        let fasta_path = first_fasta_in_workdir(&job.workdir);
        let case_hash_input = fasta_path
            .clone()
            .unwrap_or_else(|| job.workdir.join("case.toml"));
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "AlphaFold 2",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);
        let mut artefacts: Vec<Artifact> = Vec::new();

        if let Some(p) = fasta_path {
            artefacts.push(Artifact {
                path: p,
                kind: ArtifactKind::Other,
                checksum: None,
                label: "AlphaFold 2 input FASTA".to_string(),
            });
        }

        // AF2 writes outputs into `<workdir>/<query_name>/` —
        // `ranked_*.pdb` (per-model ordered by confidence) and
        // `ranking_debug.json`. Walk every immediate subdirectory so
        // we don't have to know the query name up front.
        let entries = match fs::read_dir(&job.workdir) {
            Ok(e) => e,
            Err(e) => {
                tracing::warn!(target: "valenx-alphafold2", ?e, "workdir read failed");
                return Ok(results);
            }
        };
        let mut ranked_pdbs: Vec<PathBuf> = Vec::new();
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
                    Some("pdb") if name.starts_with("ranked_") => {
                        ranked_pdbs.push(p);
                    }
                    Some("json") if name == "ranking_debug.json" => {
                        artefacts.push(Artifact {
                            path: p,
                            kind: ArtifactKind::Tabular,
                            checksum: None,
                            label: "AlphaFold 2 ranking metadata".to_string(),
                        });
                    }
                    _ => continue,
                }
            }
        }
        ranked_pdbs.sort();
        for (idx, path) in ranked_pdbs.into_iter().enumerate() {
            let model_n = idx;
            let stem = path
                .file_stem()
                .and_then(|s| s.to_str())
                .unwrap_or("ranked")
                .to_string();
            // Round-22 M2: cap the per-PDB read at MAX_PDB_FILE_BYTES
            // (256 MiB) so a poisoned workdir with a multi-GB `.pdb`
            // can't OOM `collect()` before the parser runs.
            let label = match valenx_core::io_caps::read_capped_to_string(
                &path,
                valenx_core::io_caps::MAX_PDB_FILE_BYTES as usize,
            ) {
                Ok(text) => match valenx_bio::format::pdb::read(&stem, &text) {
                    Ok(structure) => format!(
                        "AlphaFold 2 ranked model {} (predicted, {} atoms, {} residues)",
                        model_n,
                        structure.atom_count(),
                        structure.residue_count()
                    ),
                    Err(e) => format!(
                        "AlphaFold 2 ranked model {} (predicted, parse warning: {})",
                        model_n,
                        e.to_string().lines().next().unwrap_or("invalid")
                    ),
                },
                Err(_) => format!("AlphaFold 2 ranked model {model_n} (predicted)"),
            };
            artefacts.push(Artifact {
                path,
                kind: ArtifactKind::Native,
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
            ribbon_contributions: vec!["bio.alphafold2.predict"],
        }
    }
}

/// Lift the staged FASTA out of the workdir for provenance hashing.
fn first_fasta_in_workdir(workdir: &Path) -> Option<PathBuf> {
    let entries = fs::read_dir(workdir).ok()?;
    let mut hits: Vec<PathBuf> = entries
        .flatten()
        .map(|e| e.path())
        .filter(|p| {
            p.extension()
                .and_then(|s| s.to_str())
                .map(|s| {
                    let s = s.to_ascii_lowercase();
                    matches!(s.as_str(), "fasta" | "fa" | "faa" | "fna")
                })
                .unwrap_or(false)
        })
        .collect();
    hits.sort();
    hits.into_iter().next()
}

/// Run `python -c "import alphafold; print(alphafold.__version__)"` —
/// AF2's package surfaces `__version__` as a top-level attribute
/// since 2.0. Returns `None` on any failure.
fn detect_alphafold_version(python_binary: &Path) -> Option<Version> {
    let output = std::process::Command::new(python_binary)
        .arg("-c")
        .arg("import alphafold; print(alphafold.__version__)")
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

    /// Minimal valid PDB record covering one residue.
    const SAMPLE_PDB: &str = "\
ATOM      1  N   ALA A   1      11.104  13.207   2.063  1.00  0.00           N
ATOM      2  CA  ALA A   1      11.804  13.793   3.215  1.00  0.00           C
ATOM      3  C   ALA A   1      11.072  15.058   3.668  1.00  0.00           C
ATOM      4  O   ALA A   1       9.835  15.117   3.586  1.00  0.00           O
ATOM      5  CB  ALA A   1      11.916  12.789   4.357  1.00  0.00           C
END
";

    #[test]
    fn info_is_bio_domain() {
        let info = AlphaFold2Adapter::new().info();
        assert_eq!(info.id, "alphafold2");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "Apache-2.0");
        assert_eq!(info.display_name, "AlphaFold 2");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = AlphaFold2Adapter::new().info();
        assert_eq!(info.version_range.min_inclusive, Version::new(2, 3, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(3, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = AlphaFold2Adapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.alphafold2.predict"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = AlphaFold2Adapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }

    /// `collect()` must walk `<workdir>/<query_name>/` for ranked PDBs
    /// and `ranking_debug.json` and surface them as the right artifact
    /// kinds.
    #[test]
    fn collect_walks_query_subdir_for_ranked_models() {
        let tmp = std::env::temp_dir().join(format!(
            "valenx-af2-collect-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(&tmp).unwrap();
        fs::write(tmp.join("query.fasta"), b">x\nACDEF\n").unwrap();
        let q_subdir = tmp.join("query");
        fs::create_dir_all(&q_subdir).unwrap();
        fs::write(q_subdir.join("ranked_0.pdb"), SAMPLE_PDB).unwrap();
        fs::write(q_subdir.join("ranked_1.pdb"), SAMPLE_PDB).unwrap();
        fs::write(
            q_subdir.join("ranking_debug.json"),
            br#"{"plddts": {"model_1": 85.0}}"#,
        )
        .unwrap();

        let job = PreparedJob {
            workdir: tmp.clone(),
            native_command: vec![],
            environment: Vec::new(),
            estimated_runtime: None,
            kill_on_drop: true,
        };
        let results = AlphaFold2Adapter::new().collect(&job).unwrap();

        let pdbs: Vec<&Artifact> = results
            .artifacts
            .iter()
            .filter(|a| a.path.extension().is_some_and(|e| e == "pdb"))
            .collect();
        assert_eq!(pdbs.len(), 2);
        for art in &pdbs {
            assert_eq!(art.kind, ArtifactKind::Native);
            assert!(
                art.label.contains("AlphaFold 2 ranked model"),
                "label was: {}",
                art.label
            );
            assert!(
                art.label.contains("5 atoms") && art.label.contains("1 residues"),
                "label was: {}",
                art.label
            );
        }

        let json = results
            .artifacts
            .iter()
            .find(|a| a.path.extension().is_some_and(|e| e == "json"))
            .expect("ranking JSON present");
        assert_eq!(json.kind, ArtifactKind::Tabular);
        assert_eq!(json.label, "AlphaFold 2 ranking metadata");

        let fasta = results
            .artifacts
            .iter()
            .find(|a| a.path.extension().is_some_and(|e| e == "fasta"))
            .expect("FASTA present");
        assert_eq!(fasta.kind, ArtifactKind::Other);
        assert_eq!(fasta.label, "AlphaFold 2 input FASTA");

        let _ = fs::remove_dir_all(&tmp);
    }

    /// A malformed PDB shouldn't crash collect — it should degrade
    /// to a parse-warning label so the UI still surfaces the raw
    /// file.
    #[test]
    fn collect_pdb_parse_failure_degrades_gracefully() {
        let tmp = std::env::temp_dir().join(format!(
            "valenx-af2-bad-pdb-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(&tmp).unwrap();
        let q_subdir = tmp.join("query");
        fs::create_dir_all(&q_subdir).unwrap();
        // ATOM lines must be >= 78 cols; this one is far too short.
        fs::write(
            q_subdir.join("ranked_0.pdb"),
            b"ATOM      1  N   ALA A   1\n",
        )
        .unwrap();

        let job = PreparedJob {
            workdir: tmp.clone(),
            native_command: vec![],
            environment: Vec::new(),
            estimated_runtime: None,
            kill_on_drop: true,
        };
        let results = AlphaFold2Adapter::new().collect(&job).unwrap();
        let pdb_art = results
            .artifacts
            .iter()
            .find(|a| a.path.extension().is_some_and(|e| e == "pdb"))
            .expect("artifact still surfaced");
        assert_eq!(pdb_art.kind, ArtifactKind::Native);
        assert!(
            pdb_art.label.contains("parse warning"),
            "label was: {}",
            pdb_art.label
        );
        let _ = fs::remove_dir_all(&tmp);
    }

    /// `prepare()` must fail fast when `data_dir` doesn't resolve to
    /// a real directory, mirroring the FASTA-not-found pattern.
    /// Without this check AF2 would launch and crash deep in its own
    /// startup with a much less actionable message.
    #[test]
    fn prepare_rejects_missing_data_dir() {
        let tmp = tempdir("alphafold2-af2-no-datadir");
        let case_dir = tmp.join("case");
        fs::create_dir_all(&case_dir).unwrap();
        // Stage a valid run_script + FASTA so prepare() gets past the
        // earlier checks; only `data_dir` should be the failure.
        fs::write(case_dir.join("run_alphafold.py"), b"# placeholder\n").unwrap();
        fs::write(case_dir.join("query.fasta"), b">x\nACDEF\n").unwrap();
        fs::write(
            case_dir.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "alphafold2.predict"

[bio.alphafold2]
run_script        = "run_alphafold.py"
query_fasta       = "query.fasta"
data_dir          = "does_not_exist"
max_template_date = "2022-01-01"
model_preset      = "monomer"
"#,
        )
        .unwrap();
        let case = Case {
            id: "af2-no-datadir".into(),
            path: case_dir.clone(),
        };
        let workdir = tmp.join("work");
        let err = AlphaFold2Adapter::new()
            .prepare(&case, &workdir)
            .expect_err("expected InvalidCase for missing data_dir");
        match err {
            AdapterError::InvalidCase { reason, .. } => {
                assert!(
                    reason.contains("data_dir") && reason.contains("not a directory"),
                    "reason was: {reason}"
                );
            }
            other => panic!("unexpected error: {other:?}"),
        }
        let _ = fs::remove_dir_all(&tmp);
    }

    /// Round-22 M2 RED→GREEN: a `.pdb` file in the query subdir that
    /// is larger than `MAX_PDB_FILE_BYTES` (256 MiB) must produce a
    /// generic fallback label without slurping the file into memory.
    /// Pre-fix `collect()` did a bare `fs::read_to_string(&path)` for
    /// each ranked PDB and would have allocated the full file size
    /// before the PDB parser saw the first `ATOM` line.
    ///
    /// Uses `set_len` to create a sparse over-cap file without
    /// writing 500 MiB of zeros on every CI run.
    #[test]
    fn collect_skips_oversize_pdb_file() {
        let tmp = std::env::temp_dir().join(format!(
            "valenx-af2-r22m2-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(&tmp).unwrap();
        // Mandatory FASTA so collect() walks the query subdir.
        fs::write(tmp.join("query.fasta"), b">x\nACDEF\n").unwrap();
        let q_subdir = tmp.join("query");
        fs::create_dir_all(&q_subdir).unwrap();
        let pdb_path = q_subdir.join("ranked_0.pdb");
        // Past the 256 MiB MAX_PDB_FILE_BYTES cap.
        let f = fs::File::create(&pdb_path).unwrap();
        f.set_len(valenx_core::io_caps::MAX_PDB_FILE_BYTES + 1)
            .unwrap();
        drop(f);

        let job = PreparedJob {
            workdir: tmp.clone(),
            native_command: vec![],
            environment: Vec::new(),
            estimated_runtime: None,
            kill_on_drop: true,
        };
        let results = AlphaFold2Adapter::new().collect(&job).unwrap();

        let pdb_art = results
            .artifacts
            .iter()
            .find(|a| a.path.extension().is_some_and(|e| e == "pdb"))
            .expect("artifact still surfaced (discovery walk runs before read)");
        assert_eq!(pdb_art.kind, ArtifactKind::Native);
        // Pre-fix: label would have shown "X atoms, Y residues" after
        // slurping 256+ MiB into memory. Post-fix: the read errors out
        // and we fall through to the generic predicted-model label
        // (`(predicted)` — see the `Err(_) => format!(...)` branch in
        // `collect()`).
        assert!(
            pdb_art.label.ends_with("(predicted)"),
            "expected generic fallback label, got: {}",
            pdb_art.label
        );
        let _ = fs::remove_dir_all(&tmp);
    }

    /// Round-9 RED→GREEN: `[bio.alphafold2].run_script` used to be
    /// joined with bare `case.path.join`, which on POSIX silently
    /// accepted `/etc/passwd` (absolute paths replace the prefix) and
    /// on Windows accepted `..\..\..` traversal out of the case
    /// sandbox. The fix wraps the relative branch with
    /// `confined_join`, which rejects both shapes with a clear
    /// `InvalidCase` message — same policy as `query_fasta`.
    #[test]
    fn prepare_rejects_run_script_traversing_outside_case_dir() {
        let tmp = tempdir("alphafold2-runscript-trav");
        let case_dir = tmp.join("case");
        fs::create_dir_all(&case_dir).unwrap();
        fs::write(case_dir.join("query.fasta"), b">x\nACDEF\n").unwrap();
        fs::write(
            case_dir.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "alphafold2.predict"

[bio.alphafold2]
run_script        = "../../etc/passwd"
query_fasta       = "query.fasta"
data_dir          = "/tmp/af2-data"
max_template_date = "2022-01-01"
model_preset      = "monomer"
"#,
        )
        .unwrap();
        let case = Case {
            id: "af2-runscript-trav".into(),
            path: case_dir.clone(),
        };
        let workdir = tmp.join("work");
        let err = AlphaFold2Adapter::new()
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
