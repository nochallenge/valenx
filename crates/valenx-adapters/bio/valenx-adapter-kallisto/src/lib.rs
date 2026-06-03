//! # valenx-adapter-kallisto
//!
//! Adapter for [Kallisto](https://pachterlab.github.io/kallisto/) —
//! Lior Pachter's pseudoalignment-based RNA-seq quantification tool.
//! Kallisto skips the explicit alignment step entirely: it builds a
//! transcript-level de Bruijn graph index and asks "which set of
//! transcripts is this read compatible with?", then runs an EM
//! model on those compatibility classes to estimate per-transcript
//! abundance. The result is a quantification an order of magnitude
//! faster than alignment-based pipelines, with comparable accuracy.
//!
//! **Phase 20 — subprocess wrapper around `kallisto`.** The user
//! supplies a transcriptome FASTA plus 1 (single-end) or 2
//! (paired-end) FASTQ files via `[bio.kallisto]` in `case.toml`.
//! `prepare()` builds the kallisto index file (`kallisto index`)
//! unless `skip_index = true`, then composes the `kallisto quant`
//! invocation against `output_dir`. `run()` streams the
//! quantification via the shared subprocess runner.
//!
//! On `collect()` we surface kallisto's canonical outputs from
//! `output_dir`: `abundance.tsv` (per-transcript abundance table),
//! `abundance.h5` (HDF5 with bootstrap samples), and `run_info.json`
//! (recorded run statistics).

#![forbid(unsafe_code)]
#![allow(missing_docs)]

pub mod case_input;

use std::ffi::OsString;
use std::fs;
use std::path::{Path, PathBuf};
use std::time::Duration;

use semver::Version;

use valenx_core::{
    adapter_helpers::{
        confined_join, detect_tool_version_semver, find_on_path, live_provenance,
        validate_output_dir,
    },
    error::RunPhase,
    subprocess, Adapter, AdapterError, AdapterInfo, Capabilities, Case, LicenseMode, Physics,
    PreparedJob, ProbeReport, RunContext, RunReport, VersionRange,
};
use valenx_fields::{
    artifact::{Artifact, ArtifactKind},
    Results,
};

use crate::case_input::KallistoInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(KallistoAdapter::new())
}

pub struct KallistoAdapter;

impl KallistoAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for KallistoAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "kallisto";
/// Kallisto's binary candidates. `kallisto` is the canonical Linux /
/// macOS install name from Bioconda and source builds.
const BINARIES: &[&str] = &["kallisto"];

impl Adapter for KallistoAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "Kallisto",
            // Kallisto 0.50.x is the long-running stable line (the
            // pseudoalignment rewrite landed in 0.50.0); 1.0
            // reserves room for an eventual major bump.
            version_range: VersionRange {
                min_inclusive: Version::new(0, 50, 0),
                max_exclusive: Version::new(1, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            tool_license: "BSD-2-Clause",
            docs_url: "https://pachterlab.github.io/kallisto/manual",
            homepage_url: "https://pachterlab.github.io/kallisto/",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(BINARIES) {
            Some(binary_path) => {
                // `kallisto version` (no `--`) prints the version on
                // stdout. Try the conventional `--version` flag too
                // for forward compatibility.
                let found_version =
                    detect_tool_version_semver(&binary_path, &["version", "--version"]);
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
                hint: "Kallisto 0.50+ required; install via \
                       `conda install -c bioconda kallisto` or \
                       `brew install kallisto`"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = KallistoInput::from_case_dir(&case.path)?;

        // Round-10 H3: `index` is `PathBuf` and pre-fix flowed into
        // `workdir.join(&input.index)`. Validate as a basename —
        // kallisto's index is a single `.idx` file (unlike salmon's
        // directory-based index), so basename-only is correct.
        if let Some(s) = input.index.to_str() {
            valenx_core::adapter_helpers::validate_output_basename(
                s,
                "[bio.kallisto].index",
            )
            .map_err(|e| AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!("{e}"),
            })?;
        } else {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: "[bio.kallisto].index: non-UTF-8 path rejected".into(),
            });
        }

        fs::create_dir_all(workdir)?;

        // Stage `case.toml` into the workdir so collect() can recover
        // the configured `output_basename` for prefix-filtering output
        // artifacts. Without this stage, the basename filter silently
        // degrades to "match everything".
        let staged_case_toml = workdir.join("case.toml");
        let source_case_toml = case.path.join("case.toml");
        if source_case_toml.is_file() {
            fs::copy(&source_case_toml, &staged_case_toml)
                .map_err(|e| AdapterError::Other(anyhow::anyhow!("stage case.toml: {e}")))?;
        }

        // Resolve the transcriptome path against the case directory if
        // relative.
        // Round-9 hardening: relative `transcriptome` flows into the
        // kallisto command line as `-i <path>` (index build) or the
        // index input; wrap with `confined_join`. Absolute paths
        // stay supported (shared reference bundles).
        let source_transcriptome = if input.transcriptome.is_absolute() {
            input.transcriptome.clone()
        } else {
            valenx_core::adapter_helpers::confined_join(&case.path, &input.transcriptome)?
        };
        if !source_transcriptome.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.kallisto].transcriptome `{}` not found (resolved {})",
                    input.transcriptome.display(),
                    source_transcriptome.display()
                ),
            });
        }

        // Resolve the index path against the workdir if relative.
        // Unlike Salmon's directory-based index, kallisto's index is
        // a single `.idx` file we either create here or read.
        let resolved_index = if input.index.is_absolute() {
            input.index.clone()
        } else {
            workdir.join(&input.index)
        };

        // Resolve each read file against the case directory via
        // `confined_join` — rejects absolute paths and `..` traversal
        // out of the case sandbox (round-6 hardening). We *don't*
        // copy them into the workdir — RNA-seq FASTQs routinely run
        // to tens of GB and kallisto reads them by path.
        let mut resolved_reads: Vec<PathBuf> = Vec::with_capacity(input.reads.len());
        for read in &input.reads {
            let resolved = confined_join(&case.path, read)?;
            if !resolved.is_file() {
                return Err(AdapterError::InvalidCase {
                    case_path: case.path.join("case.toml"),
                    reason: format!(
                        "[bio.kallisto].reads entry `{}` not found (resolved {})",
                        read.display(),
                        resolved.display()
                    ),
                });
            }
            resolved_reads.push(resolved);
        }

        // Round-5: validate the user-supplied `output_dir` is a safe
        // workdir-relative subdirectory before joining it in. Rejects
        // absolute paths AND `..` traversal even though we resolve
        // against the workdir (not the case dir).
        validate_output_dir(&input.output_dir, "[bio.kallisto].output_dir")?;

        // Resolve the quant output directory against the workdir if
        // relative. We don't require it to pre-exist — kallisto
        // creates it.
        let resolved_output_dir = if input.output_dir.is_absolute() {
            input.output_dir.clone()
        } else {
            workdir.join(&input.output_dir)
        };

        let binary_path = find_on_path(BINARIES).ok_or_else(|| AdapterError::ToolNotInstalled {
            name: INFO_ID,
            hint: "Kallisto 0.50+ required; install via \
                       `conda install -c bioconda kallisto` or \
                       `brew install kallisto`"
                .into(),
        })?;

        // Build the kallisto index unless the user opted out. Running
        // it here (in prepare) keeps the subsequent `run()` call as a
        // single `kallisto quant` invocation, which lets the
        // subprocess runner stream stderr line-by-line without
        // chaining commands through a shell.
        if !input.skip_index {
            let index_status = std::process::Command::new(&binary_path)
                .arg("index")
                .arg("-i")
                .arg(&resolved_index)
                .arg(&source_transcriptome)
                .current_dir(workdir)
                .stdin(std::process::Stdio::null())
                .stdout(std::process::Stdio::piped())
                .stderr(std::process::Stdio::piped())
                .output();
            match index_status {
                Ok(out) if out.status.success() => {}
                Ok(out) => {
                    let stderr = String::from_utf8_lossy(&out.stderr);
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "`kallisto index -i {} {}` failed (exit {}): {}",
                        resolved_index.display(),
                        source_transcriptome.display(),
                        out.status.code().unwrap_or(-1),
                        stderr.lines().next().unwrap_or("(no stderr)")
                    )));
                }
                Err(e) => {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "spawning `kallisto index -i {}` failed: {e}",
                        resolved_index.display()
                    )));
                }
            }
        }

        // Compose the quantification invocation. Kallisto's flag
        // surface differs slightly between paired-end (no `--single`,
        // reads as positional argv) and single-end (needs `--single`
        // plus `-l` mean and `-s` stdev for the fragment-length
        // distribution).
        let mut native_command: Vec<OsString> = vec![
            binary_path.into_os_string(),
            OsString::from("quant"),
            OsString::from("-i"),
            resolved_index.into_os_string(),
            OsString::from("-o"),
            resolved_output_dir.into_os_string(),
            OsString::from("-t"),
            OsString::from(input.threads.to_string()),
        ];
        if resolved_reads.len() == 2 {
            // Paired-end: <reads[0]> <reads[1]> as positional args.
            native_command.push(resolved_reads[0].clone().into_os_string());
            native_command.push(resolved_reads[1].clone().into_os_string());
        } else {
            // Single-end: --single -l <mean> -s <sd> <reads[0]>. The
            // case-input parser already verified the fragment stats
            // are present and positive when reads.len() == 1.
            let l = input
                .fragment_length
                .expect("single-end fragment_length validated by case_input");
            let s = input
                .fragment_sd
                .expect("single-end fragment_sd validated by case_input");
            native_command.push(OsString::from("--single"));
            native_command.push(OsString::from("-l"));
            native_command.push(OsString::from(format!("{l}")));
            native_command.push(OsString::from("-s"));
            native_command.push(OsString::from(format!("{s}")));
            native_command.push(resolved_reads[0].clone().into_os_string());
        }
        for arg in &input.extra_args {
            native_command.push(OsString::from(arg));
        }

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            // Kallisto pseudoalignment is famously fast — even a
            // full-mammalian RNA-seq library finishes in minutes —
            // but 4 hours mirrors the rest of the RNA-seq adapter
            // family and absorbs the long tail (giant indices,
            // contended I/O).
            estimated_runtime: Some(Duration::from_secs(4 * 60 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting Kallisto", |line| {
            let mut hint = subprocess::Hint::default();
            // Kallisto emits `[index]`, `[quant]`, and
            // `[bstrp]` (bootstrap) bracketed phase markers on
            // stderr; the final `[  bstrp]` / "running EM"
            // sequence wraps the quantification.
            if line.contains("[bstrp]") || line.contains("running EM") {
                hint.progress = Some((95.0, line.to_string()));
            } else if line.contains("[quant]") || line.contains("processed") {
                hint.progress = Some((50.0, line.to_string()));
            } else if line.contains("Error") || line.contains("[error]") || line.contains("ERROR:")
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
        // Kallisto's outputs land inside the user-supplied
        // `output_dir`. We resolve it the same way `prepare()` did so
        // the `collect()` walk hits the same directory regardless of
        // whether the user gave us an absolute path or one relative
        // to the workdir.
        let input_for_paths = KallistoInput::from_case_dir(&job.workdir).ok();
        let output_dir = match input_for_paths {
            Some(i) => {
                if i.output_dir.is_absolute() {
                    i.output_dir
                } else {
                    job.workdir.join(&i.output_dir)
                }
            }
            None => job.workdir.clone(),
        };

        // Provenance: hash the staged `abundance.tsv` if present (the
        // canonical run output). Falls back to case.toml when the
        // quantification hasn't produced one yet.
        let case_hash_input = {
            let tsv = output_dir.join("abundance.tsv");
            if tsv.is_file() {
                tsv
            } else {
                job.workdir.join("case.toml")
            }
        };
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "Kallisto",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);
        let mut artefacts: Vec<Artifact> = Vec::new();

        // Walk the output directory top-level. Kallisto writes three
        // canonical files: `abundance.tsv` (the per-transcript
        // abundance table), `abundance.h5` (HDF5 with bootstrap
        // samples for downstream Sleuth), and `run_info.json` (run
        // statistics). We surface all three per the output spec.
        let entries = match fs::read_dir(&output_dir) {
            Ok(e) => e,
            Err(e) => {
                tracing::warn!(target: "valenx-kallisto", ?e, "output_dir read failed");
                return Ok(results);
            }
        };
        for entry in entries.flatten() {
            let path = entry.path();
            if !path.is_file() {
                continue;
            }
            let name = match path.file_name().and_then(|s| s.to_str()) {
                Some(n) => n.to_string(),
                None => continue,
            };
            let (kind, label) = match name.as_str() {
                "abundance.tsv" => (
                    ArtifactKind::Tabular,
                    "Kallisto transcript abundance".to_string(),
                ),
                "abundance.h5" => (ArtifactKind::Native, "Kallisto HDF5 abundance".to_string()),
                "run_info.json" => (ArtifactKind::Log, "Kallisto run info".to_string()),
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
            ribbon_contributions: vec!["bio.kallisto.quant"],
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn info_is_bio_domain() {
        let info = KallistoAdapter::new().info();
        assert_eq!(info.id, "kallisto");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "BSD-2-Clause");
        assert_eq!(info.display_name, "Kallisto");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = KallistoAdapter::new().info();
        // Kallisto 0.50.x is the de facto stable line; 1.0 reserves
        // room for an eventual major bump.
        assert_eq!(info.version_range.min_inclusive, Version::new(0, 50, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(1, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = KallistoAdapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.kallisto.quant"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = KallistoAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }

    /// Round-9 RED→GREEN: `[bio.kallisto].transcriptome` used to be
    /// joined with bare `case.path.join`. Wrap with `confined_join`.
    #[test]
    fn prepare_rejects_transcriptome_traversing_outside_case_dir() {
        use valenx_test_utils::tempdir;
        let d = tempdir("kallisto-trans-trav");
        std::fs::write(d.join("r1.fq"), b"@x\nACGT\n+\nIIII\n").unwrap();
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "kallisto.quant"

[bio.kallisto]
transcriptome = "../../etc/passwd"
index         = "idx.kdx"
reads         = ["r1.fq"]
output_dir    = "out"
threads       = 1
fragment_length = 200.0
fragment_sd     = 20.0
"#,
        )
        .unwrap();
        let case = Case {
            id: "kallisto-trans-trav".into(),
            path: d.clone(),
        };
        let workdir = d.join("workdir");
        let err = KallistoAdapter::new()
            .prepare(&case, &workdir)
            .unwrap_err();
        let msg = format!("{err}");
        assert!(
            msg.contains("..") || msg.contains("stay within") || msg.contains("escape"),
            "expected confined_join rejection, got: {msg}"
        );
        let _ = std::fs::remove_dir_all(&d);
    }

    /// Round-10 H3 RED→GREEN: `index` is `PathBuf` and pre-fix
    /// flowed into `workdir.join(&input.index)`. Hostile
    /// `index = "../etc/passwd"` is now rejected.
    #[test]
    fn prepare_rejects_index_path_traversal() {
        use valenx_test_utils::tempdir;
        let d = tempdir("kallisto-index-trav");
        std::fs::write(d.join("trans.fa"), b">x\nACGT\n").unwrap();
        std::fs::write(d.join("r1.fq"), b"@x\nACGT\n+\nIIII\n").unwrap();
        // Single-end runs need fragment_length + fragment_sd; supply
        // them so the case-input parser doesn't short-circuit before
        // the `index` field is reached.
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "kallisto.quant"

[bio.kallisto]
transcriptome    = "trans.fa"
index            = "../etc/passwd"
reads            = ["r1.fq"]
output_dir       = "quant"
fragment_length  = 200.0
fragment_sd      = 20.0
"#,
        )
        .unwrap();
        let case = Case {
            id: "trav".into(),
            path: d.clone(),
        };
        let workdir = d.join("workdir");
        let err = KallistoAdapter::new()
            .prepare(&case, &workdir)
            .unwrap_err();
        let msg = format!("{err}");
        assert!(
            msg.contains("[bio.kallisto].index"),
            "expected [bio.kallisto].index in error, got: {msg}"
        );
        let _ = std::fs::remove_dir_all(&d);
    }
}
