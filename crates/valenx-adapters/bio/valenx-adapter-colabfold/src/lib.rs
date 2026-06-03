//! # valenx-adapter-colabfold
//!
//! Adapter for [ColabFold](https://github.com/sokrypton/ColabFold) —
//! the open-source protein structure prediction pipeline that wraps
//! AlphaFold2 / OmegaFold / ESMFold with a faster MSA generation step
//! (MMseqs2 in place of the original AlphaFold's heavy
//! database-search). The bread-and-butter use case: hand it a FASTA
//! query, get back ranked PDB models with per-residue pLDDT
//! confidence scores.
//!
//! **Phase 17 — subprocess wrapper for the public CLI.** The user
//! references their input FASTA via `[bio.colabfold].input_fasta`
//! in `case.toml`; `prepare()` stages the FASTA into the workdir;
//! `run()` invokes `colabfold_batch <input.fasta> result/` via the
//! shared subprocess runner. Optional `num_recycles` / `num_models`
//! knobs surface as `--num-recycle` / `--num-models` flags only when
//! the user overrides the typical defaults, keeping the common case's
//! command line clean.
//!
//! On `collect()` we walk `result/` for `*.pdb` outputs (one per
//! predicted model). Each PDB is parsed via
//! [`valenx_bio::format::pdb::read`] and surfaced as a typed
//! [`ArtifactKind::Native`] artifact with a "ColabFold model N
//! (predicted)" label. Per-model pLDDT extraction from the customary
//! `<id>_relaxed_rank_001_<plddt>_..pdb` filenames is a Phase 17.5
//! enhancement.

#![forbid(unsafe_code)]
#![allow(missing_docs)]

pub mod case_input;

use std::ffi::OsString;
use std::fs;
use std::path::{Path, PathBuf};
use std::time::Duration;

use semver::Version;

use valenx_core::{
    adapter_helpers::{confined_join, detect_tool_version_semver, find_on_path, live_provenance},
    error::RunPhase,
    subprocess, Adapter, AdapterError, AdapterInfo, Capabilities, Case, LicenseMode, Physics,
    PreparedJob, ProbeReport, RunContext, RunReport, VersionRange,
};
use valenx_fields::{
    artifact::{Artifact, ArtifactKind},
    Results,
};

use crate::case_input::ColabFoldInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(ColabFoldAdapter::new())
}

pub struct ColabFoldAdapter;

impl ColabFoldAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for ColabFoldAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "colabfold";
/// ColabFold ships its public CLI as `colabfold_batch` (snake-case)
/// in the upstream pip package — no alternate spellings.
const BINARIES: &[&str] = &["colabfold_batch"];
/// ColabFold's CLI defaults — match these to elide the `--num-recycle`
/// / `--num-models` flags from the run command in the typical case.
const DEFAULT_NUM_RECYCLES: u32 = 3;
const DEFAULT_NUM_MODELS: u32 = 5;
/// Subdirectory ColabFold writes outputs into. Hard-wired to
/// `result/` because the run command always passes that name as the
/// second positional argument.
const RESULT_DIR: &str = "result";

impl Adapter for ColabFoldAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "ColabFold",
            // ColabFold 1.5 is the first release with the stable
            // MMseqs2 API surface and the `colabfold_batch` CLI we
            // lean on; upper bound 2.0 reserves room for an upcoming
            // major bump.
            version_range: VersionRange {
                min_inclusive: Version::new(1, 5, 0),
                max_exclusive: Version::new(2, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            tool_license: "MIT",
            docs_url: "https://github.com/sokrypton/ColabFold",
            homepage_url: "https://github.com/sokrypton/ColabFold",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(BINARIES) {
            Some(binary_path) => {
                // ColabFold's CLI prints version on `--version`.
                let found_version = detect_tool_version_semver(&binary_path, &["--version", "-v"]);
                Ok(ProbeReport {
                    ok: true,
                    found_version,
                    binary_path: Some(binary_path),
                    warnings: Vec::new(),
                    required_env: Vec::new(),
                })
            }
            None => Err(AdapterError::ToolNotInstalled {
                name: INFO_ID,
                hint: "ColabFold 1.5+ required; install via \
                       `pip install colabfold` and ensure \
                       `colabfold_batch` is on PATH"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = ColabFoldInput::from_case_dir(&case.path)?;

        fs::create_dir_all(workdir)?;

        // Stage the FASTA. `confined_join` rejects absolute paths and
        // `..` traversal so the staged copy stays confined to the case
        // directory.
        let source_fasta = confined_join(&case.path, &input.input_fasta)?;
        if !source_fasta.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.colabfold].input_fasta `{}` not found (resolved {})",
                    input.input_fasta.display(),
                    source_fasta.display()
                ),
            });
        }
        let fasta_filename =
            input
                .input_fasta
                .file_name()
                .ok_or_else(|| AdapterError::InvalidCase {
                    case_path: case.path.join("case.toml"),
                    reason: format!(
                        "[bio.colabfold].input_fasta path `{}` has no filename",
                        input.input_fasta.display()
                    ),
                })?;
        let dest_fasta = workdir.join(fasta_filename);
        if source_fasta != dest_fasta {
            fs::copy(&source_fasta, &dest_fasta)?;
        }

        let binary_path = find_on_path(BINARIES).ok_or_else(|| AdapterError::ToolNotInstalled {
            name: INFO_ID,
            hint: "could not locate `colabfold_batch` on PATH".into(),
        })?;

        // Build the command. ColabFold's positional contract is
        // `colabfold_batch <input> <result_dir>`; only the override
        // flags surface when the user diverges from the defaults so
        // the typical command line stays clean.
        let mut native_command: Vec<OsString> = vec![
            binary_path.into_os_string(),
            OsString::from(fasta_filename),
            OsString::from(RESULT_DIR),
        ];
        if input.num_recycles != DEFAULT_NUM_RECYCLES {
            native_command.push(OsString::from("--num-recycle"));
            native_command.push(OsString::from(input.num_recycles.to_string()));
        }
        if input.num_models != DEFAULT_NUM_MODELS {
            native_command.push(OsString::from("--num-models"));
            native_command.push(OsString::from(input.num_models.to_string()));
        }

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            // Structure prediction varies enormously by sequence
            // length and GPU. Short peptides finish in minutes;
            // 1000+ residue targets on consumer GPUs run for hours.
            // 4 hours is a reasonable default; long runs override
            // through their own progress reporting.
            estimated_runtime: Some(Duration::from_secs(4 * 60 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting ColabFold", |line| {
            let mut hint = subprocess::Hint::default();
            // ColabFold's stdout marches through three coarse phases:
            // MSA generation -> folding -> Amber relax. Map each to a
            // progress tick so the UI shows movement. The exact log
            // strings are emitted by colabfold_batch itself; matching
            // them substring-wise keeps us robust to formatting
            // tweaks across point releases.
            if line.contains("Generating MSAs") {
                hint.progress = Some((30.0, line.to_string()));
            } else if line.contains("Folding") {
                hint.progress = Some((60.0, line.to_string()));
            } else if line.contains("relaxed") {
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
        // Provenance: hash the staged FASTA (the canonical
        // "this case is configured this way" input). Falls back to
        // case.toml when the FASTA isn't present yet.
        let fasta_path = first_fasta_in_workdir(&job.workdir);
        let case_hash_input = fasta_path
            .clone()
            .unwrap_or_else(|| job.workdir.join("case.toml"));
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "ColabFold",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);
        let mut artefacts: Vec<Artifact> = Vec::new();

        // Surface the staged FASTA so the user can find their query
        // sequence next to the predictions.
        if let Some(p) = fasta_path {
            artefacts.push(Artifact {
                path: p,
                kind: ArtifactKind::Other,
                checksum: None,
                label: "ColabFold input FASTA".to_string(),
            });
        }

        // Walk result/ for the customary outputs. ColabFold writes
        // a stable mix of files: `*.pdb` (one per ranked model),
        // `*.json` (per-residue confidence scores), `*.png` (PAE
        // plots), `log.txt` (run log).
        let result_dir = job.workdir.join(RESULT_DIR);
        if let Ok(entries) = fs::read_dir(&result_dir) {
            // Track model index for PDBs ranked by filename; collect
            // up front so we can sort lexicographically (which lines
            // up with rank because ColabFold encodes rank in the
            // filename as `rank_001`, `rank_002`, …).
            let mut pdb_paths: Vec<PathBuf> = Vec::new();
            for entry in entries.flatten() {
                let path = entry.path();
                if !path.is_file() {
                    continue;
                }
                let ext = path
                    .extension()
                    .and_then(|s| s.to_str())
                    .map(|s| s.to_ascii_lowercase());
                match ext.as_deref() {
                    Some("pdb") => pdb_paths.push(path),
                    Some("json") => artefacts.push(Artifact {
                        path,
                        kind: ArtifactKind::Other,
                        checksum: None,
                        label: "ColabFold confidence JSON".to_string(),
                    }),
                    Some("png") => artefacts.push(Artifact {
                        path,
                        kind: ArtifactKind::Image,
                        checksum: None,
                        label: "ColabFold PAE plot".to_string(),
                    }),
                    Some("txt") | Some("log") => artefacts.push(Artifact {
                        path,
                        kind: ArtifactKind::Log,
                        checksum: None,
                        label: "ColabFold log".to_string(),
                    }),
                    _ => continue,
                }
            }
            pdb_paths.sort();
            for (idx, path) in pdb_paths.into_iter().enumerate() {
                // Soft-validate via the canonical PDB reader. Failed
                // parses degrade to a parse-warning label so the user
                // can still surface the raw file. Phase 17.5 will
                // extract the pLDDT score from the filename.
                let model_n = idx + 1;
                let stem = path
                    .file_stem()
                    .and_then(|s| s.to_str())
                    .unwrap_or("model")
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
                            "ColabFold model {} (predicted, {} atoms, {} residues)",
                            model_n,
                            structure.atom_count(),
                            structure.residue_count()
                        ),
                        Err(e) => format!(
                            "ColabFold model {} (predicted, parse warning: {})",
                            model_n,
                            e.to_string().lines().next().unwrap_or("invalid")
                        ),
                    },
                    Err(_) => format!("ColabFold model {model_n} (predicted)"),
                };
                artefacts.push(Artifact {
                    path,
                    kind: ArtifactKind::Native,
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
        // The bio-specific Capability variants land in a follow-up
        // task; ribbon contributions are already enough for the
        // registry to surface the adapter.
        Capabilities {
            capabilities: Vec::new(),
            ribbon_contributions: vec!["bio.colabfold.predict"],
        }
    }
}

/// Lift the staged FASTA out of the workdir for provenance hashing.
/// Returns the lexicographically-first `.fasta` / `.fa` / `.faa` /
/// `.fna` file at the top level, or `None` if none exists yet.
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

#[cfg(test)]
mod tests {
    use super::*;

    /// Minimal valid PDB record covering one residue. Mirrors the
    /// fixture used by the OpenMM adapter so collect()'s PDB-parse
    /// path exercises against the same canonical reader.
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
        let info = ColabFoldAdapter::new().info();
        assert_eq!(info.id, "colabfold");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "MIT");
        assert_eq!(info.display_name, "ColabFold");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = ColabFoldAdapter::new().info();
        assert_eq!(info.version_range.min_inclusive, Version::new(1, 5, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(2, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = ColabFoldAdapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.colabfold.predict"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = ColabFoldAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }

    #[test]
    fn first_fasta_in_workdir_picks_lexicographic_first() {
        let tmp = std::env::temp_dir().join(format!(
            "valenx-colabfold-fasta-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(&tmp).unwrap();
        fs::write(tmp.join("z_late.fasta"), b">x\nACGT\n").unwrap();
        fs::write(tmp.join("a_first.fa"), b">x\nACGT\n").unwrap();
        fs::write(tmp.join("not_fasta.txt"), b"placeholder").unwrap();
        let f = first_fasta_in_workdir(&tmp).expect("found");
        assert!(f.ends_with("a_first.fa"));
        let _ = fs::remove_dir_all(&tmp);
    }

    #[test]
    fn first_fasta_in_workdir_returns_none_when_empty() {
        let tmp = std::env::temp_dir().join(format!(
            "valenx-colabfold-empty-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(&tmp).unwrap();
        fs::write(tmp.join("notes.md"), b"placeholder").unwrap();
        assert!(first_fasta_in_workdir(&tmp).is_none());
        let _ = fs::remove_dir_all(&tmp);
    }

    /// `collect()` must walk `result/` for PDB / JSON / PNG / log
    /// outputs and surface them with the right artifact kind and
    /// labels. The PDB content here parses cleanly so the model label
    /// must include atom + residue counts from the parser.
    #[test]
    fn collect_walks_result_dir_and_classifies_outputs() {
        let tmp = std::env::temp_dir().join(format!(
            "valenx-colabfold-collect-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(&tmp).unwrap();
        fs::write(tmp.join("query.fasta"), b">x\nACDEF\n").unwrap();
        let result_dir = tmp.join(RESULT_DIR);
        fs::create_dir_all(&result_dir).unwrap();
        // Two ranked PDBs — ColabFold's customary filename pattern.
        fs::write(
            result_dir.join("query_relaxed_rank_001_alphafold2_model_3.pdb"),
            SAMPLE_PDB,
        )
        .unwrap();
        fs::write(
            result_dir.join("query_relaxed_rank_002_alphafold2_model_1.pdb"),
            SAMPLE_PDB,
        )
        .unwrap();
        fs::write(
            result_dir.join("query_scores.json"),
            br#"{"plddt": [85.0, 90.0]}"#,
        )
        .unwrap();
        fs::write(result_dir.join("query_pae.png"), b"\x89PNG\r\n").unwrap();
        fs::write(result_dir.join("log.txt"), b"line one\nline two\n").unwrap();

        let job = PreparedJob {
            workdir: tmp.clone(),
            native_command: vec![],
            environment: Vec::new(),
            estimated_runtime: None,
            kill_on_drop: true,
        };
        let results = ColabFoldAdapter::new().collect(&job).unwrap();

        // Two PDB models, both Native, each with the model-N label.
        let pdbs: Vec<&Artifact> = results
            .artifacts
            .iter()
            .filter(|a| a.path.extension().is_some_and(|e| e == "pdb"))
            .collect();
        assert_eq!(pdbs.len(), 2);
        for art in &pdbs {
            assert_eq!(art.kind, ArtifactKind::Native);
            assert!(
                art.label.contains("ColabFold model"),
                "label was: {}",
                art.label
            );
            assert!(art.label.contains("predicted"), "label was: {}", art.label);
            assert!(
                art.label.contains("5 atoms") && art.label.contains("1 residues"),
                "label was: {}",
                art.label
            );
        }
        // Lexicographic sort + 1-indexed model numbering: rank_001
        // is model 1, rank_002 is model 2.
        let model1 = pdbs
            .iter()
            .find(|a| {
                a.path
                    .file_name()
                    .is_some_and(|n| n.to_string_lossy().contains("rank_001"))
            })
            .expect("rank_001 PDB present");
        assert!(model1.label.contains("model 1"));
        let model2 = pdbs
            .iter()
            .find(|a| {
                a.path
                    .file_name()
                    .is_some_and(|n| n.to_string_lossy().contains("rank_002"))
            })
            .expect("rank_002 PDB present");
        assert!(model2.label.contains("model 2"));

        // JSON, PNG, and log classify into the expected artifact kinds.
        let json = results
            .artifacts
            .iter()
            .find(|a| a.path.extension().is_some_and(|e| e == "json"))
            .expect("JSON artifact");
        assert_eq!(json.kind, ArtifactKind::Other);
        assert_eq!(json.label, "ColabFold confidence JSON");

        let png = results
            .artifacts
            .iter()
            .find(|a| a.path.extension().is_some_and(|e| e == "png"))
            .expect("PNG artifact");
        assert_eq!(png.kind, ArtifactKind::Image);

        let log = results
            .artifacts
            .iter()
            .find(|a| a.path.extension().is_some_and(|e| e == "txt"))
            .expect("log artifact");
        assert_eq!(log.kind, ArtifactKind::Log);

        // Input FASTA is also surfaced (Other kind).
        let fasta = results
            .artifacts
            .iter()
            .find(|a| a.path.extension().is_some_and(|e| e == "fasta"))
            .expect("FASTA artifact");
        assert_eq!(fasta.kind, ArtifactKind::Other);
        assert_eq!(fasta.label, "ColabFold input FASTA");

        let _ = fs::remove_dir_all(&tmp);
    }

    /// A malformed PDB shouldn't crash collect — it should degrade
    /// to a parse-warning label so the UI still surfaces the raw
    /// file. Mirrors the same edge-case handling in the OpenMM
    /// adapter.
    #[test]
    fn collect_pdb_parse_failure_degrades_gracefully() {
        let tmp = std::env::temp_dir().join(format!(
            "valenx-colabfold-bad-pdb-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(&tmp).unwrap();
        let result_dir = tmp.join(RESULT_DIR);
        fs::create_dir_all(&result_dir).unwrap();
        // ATOM lines must be >= 78 cols; this one is far too short.
        fs::write(
            result_dir.join("broken_rank_001.pdb"),
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
        let results = ColabFoldAdapter::new().collect(&job).unwrap();
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

    /// `collect()` against a workdir with no `result/` subdirectory
    /// must not crash — it should yield only the input FASTA (when
    /// staged) plus any other top-level artifacts. Common during a
    /// failed prepare() + retry flow.
    #[test]
    fn collect_with_missing_result_dir_returns_input_only() {
        let tmp = std::env::temp_dir().join(format!(
            "valenx-colabfold-no-result-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(&tmp).unwrap();
        fs::write(tmp.join("query.fasta"), b">x\nM\n").unwrap();
        let job = PreparedJob {
            workdir: tmp.clone(),
            native_command: vec![],
            environment: Vec::new(),
            estimated_runtime: None,
            kill_on_drop: true,
        };
        let results = ColabFoldAdapter::new().collect(&job).unwrap();
        // Only the input FASTA surfaces — no PDBs / JSON / PNG.
        assert_eq!(results.artifacts.len(), 1);
        assert_eq!(results.artifacts[0].label, "ColabFold input FASTA");
        let _ = fs::remove_dir_all(&tmp);
    }
}
