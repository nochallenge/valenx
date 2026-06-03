//! # valenx-adapter-salmon
//!
//! Adapter for [Salmon](https://combine-lab.github.io/salmon/) — Rob
//! Patro's transcript-level quantification tool. Salmon estimates
//! per-transcript abundance directly from RNA-seq reads using
//! quasi-mapping plus a two-phase EM model, and is the de-facto
//! standard for fast RNA-seq quantification (the alternative being
//! Kallisto's pseudoalignment approach).
//!
//! **Phase 20 — subprocess wrapper around `salmon`.** The user
//! supplies a transcriptome FASTA plus 1 (single-end) or 2
//! (paired-end) FASTQ files via `[bio.salmon]` in `case.toml`.
//! `prepare()` builds the salmon index in `index_dir` (`salmon
//! index`) unless `skip_index = true`, then composes the `salmon
//! quant` invocation against `output_dir`. `run()` streams the
//! quantification via the shared subprocess runner; salmon prints
//! progress to stderr.
//!
//! On `collect()` we surface salmon's canonical outputs from
//! `output_dir`: `quant.sf` (the per-transcript abundance table)
//! and `cmd_info.json` (the recorded invocation).

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

use crate::case_input::SalmonInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(SalmonAdapter::new())
}

pub struct SalmonAdapter;

impl SalmonAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for SalmonAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "salmon";
/// Salmon's binary candidates. `salmon` is the canonical Linux /
/// macOS install name from Bioconda and source builds.
const BINARIES: &[&str] = &["salmon"];

impl Adapter for SalmonAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "Salmon",
            // Salmon 1.10.x is the long-running stable line from
            // late 2022 onward; 2.0 reserves room for an eventual
            // major bump.
            version_range: VersionRange {
                min_inclusive: Version::new(1, 10, 0),
                max_exclusive: Version::new(2, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            tool_license: "GPL-3.0",
            docs_url: "https://salmon.readthedocs.io/",
            homepage_url: "https://combine-lab.github.io/salmon/",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(BINARIES) {
            Some(binary_path) => {
                // `salmon --version` prints the version on stdout.
                let found_version = detect_tool_version_semver(&binary_path, &["--version"]);
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
                hint: "Salmon 1.10+ required; install via \
                       `conda install -c bioconda salmon` or \
                       `brew install salmon`"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = SalmonInput::from_case_dir(&case.path)?;

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
        // salmon command line as `-t <path>` (index) or the index
        // input; wrap with `confined_join`. Absolute paths stay
        // supported (admin-managed shared references).
        let source_transcriptome = if input.transcriptome.is_absolute() {
            input.transcriptome.clone()
        } else {
            valenx_core::adapter_helpers::confined_join(&case.path, &input.transcriptome)?
        };
        if !source_transcriptome.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.salmon].transcriptome `{}` not found (resolved {})",
                    input.transcriptome.display(),
                    source_transcriptome.display()
                ),
            });
        }

        // Round-10 H3: validate `index_dir` before joining. The index
        // is a directory salmon writes into / reads from; multi-
        // component relative subpaths (`indices/v45`) are legitimate,
        // so `validate_output_dir` (allows multi-component, rejects
        // `..` traversal + absolute paths) is the correct guard —
        // not the basename helper.
        valenx_core::adapter_helpers::validate_output_dir(
            &input.index_dir,
            "[bio.salmon].index_dir",
        )?;

        // Resolve the index directory against the workdir if relative.
        // Same convention as STAR — the index is a directory we either
        // build into or read pre-built from.
        let resolved_index_dir = if input.index_dir.is_absolute() {
            input.index_dir.clone()
        } else {
            workdir.join(&input.index_dir)
        };

        // Resolve each read file against the case directory via
        // `confined_join` — rejects absolute paths and `..` traversal
        // out of the case sandbox (round-6 hardening). We *don't*
        // copy them into the workdir — RNA-seq FASTQs routinely run
        // to tens of GB and salmon reads them by path.
        let mut resolved_reads: Vec<PathBuf> = Vec::with_capacity(input.reads.len());
        for read in &input.reads {
            let resolved = confined_join(&case.path, read)?;
            if !resolved.is_file() {
                return Err(AdapterError::InvalidCase {
                    case_path: case.path.join("case.toml"),
                    reason: format!(
                        "[bio.salmon].reads entry `{}` not found (resolved {})",
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
        validate_output_dir(&input.output_dir, "[bio.salmon].output_dir")?;

        // Resolve the quant output directory against the workdir if
        // relative. We don't require it to pre-exist — salmon creates
        // it.
        let resolved_output_dir = if input.output_dir.is_absolute() {
            input.output_dir.clone()
        } else {
            workdir.join(&input.output_dir)
        };

        let binary_path = find_on_path(BINARIES).ok_or_else(|| AdapterError::ToolNotInstalled {
            name: INFO_ID,
            hint: "Salmon 1.10+ required; install via \
                       `conda install -c bioconda salmon` or \
                       `brew install salmon`"
                .into(),
        })?;

        // Build the salmon index in `index_dir` unless the user opted
        // out. Running the indexer here (in prepare) keeps the
        // subsequent `run()` call as a single `salmon quant`
        // invocation, which lets the subprocess runner stream stderr
        // line-by-line without chaining commands through a shell.
        if !input.skip_index {
            let index_status = std::process::Command::new(&binary_path)
                .arg("index")
                .arg("-t")
                .arg(&source_transcriptome)
                .arg("-i")
                .arg(&resolved_index_dir)
                .arg("-p")
                .arg(input.threads.to_string())
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
                        "`salmon index -t {} -i {}` failed (exit {}): {}",
                        source_transcriptome.display(),
                        resolved_index_dir.display(),
                        out.status.code().unwrap_or(-1),
                        stderr.lines().next().unwrap_or("(no stderr)")
                    )));
                }
                Err(e) => {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "spawning `salmon index -t {}` failed: {e}",
                        source_transcriptome.display()
                    )));
                }
            }
        }

        // Compose the quantification invocation:
        //   salmon quant -i <index_dir> -l <libtype> -p <threads>
        //                -o <output_dir>
        //                (-r <reads[0]>) | (-1 <reads[0]> -2 <reads[1]>)
        //                [extras...]
        let mut native_command: Vec<OsString> = vec![
            binary_path.into_os_string(),
            OsString::from("quant"),
            OsString::from("-i"),
            resolved_index_dir.into_os_string(),
            OsString::from("-l"),
            OsString::from(&input.libtype),
            OsString::from("-p"),
            OsString::from(input.threads.to_string()),
            OsString::from("-o"),
            resolved_output_dir.into_os_string(),
        ];
        if resolved_reads.len() == 1 {
            native_command.push(OsString::from("-r"));
            native_command.push(resolved_reads.remove(0).into_os_string());
        } else {
            native_command.push(OsString::from("-1"));
            native_command.push(resolved_reads[0].clone().into_os_string());
            native_command.push(OsString::from("-2"));
            native_command.push(resolved_reads[1].clone().into_os_string());
        }
        for arg in &input.extra_args {
            native_command.push(OsString::from(arg));
        }

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            // Salmon quantification typically completes in tens of
            // minutes on a single bulk RNA-seq library; 4 hours
            // mirrors the rest of the RNA-seq adapter family for the
            // long tail.
            estimated_runtime: Some(Duration::from_secs(4 * 60 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting Salmon", |line| {
            let mut hint = subprocess::Hint::default();
            // Salmon emits its progress via tqdm-like markers on
            // stderr — "processed N reads" mid-run, then a final
            // summary block. We pin the "Mapping rate" / "writing
            // output" markers near completion.
            if line.contains("writing output") || line.contains("Mapping rate") {
                hint.progress = Some((95.0, line.to_string()));
            } else if line.contains("processed") && line.contains("reads") {
                hint.progress = Some((50.0, line.to_string()));
            } else if line.contains("Error") || line.contains("[error]") {
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
        // Salmon's outputs land inside the user-supplied
        // `output_dir`. We resolve it the same way `prepare()` did so
        // the `collect()` walk hits the same directory regardless of
        // whether the user gave us an absolute path or one relative to
        // the workdir.
        let input_for_paths = SalmonInput::from_case_dir(&job.workdir).ok();
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

        // Provenance: hash the staged `quant.sf` if present (the
        // canonical run output). Falls back to case.toml when the
        // quantification hasn't produced one yet.
        let case_hash_input = {
            let qsf = output_dir.join("quant.sf");
            if qsf.is_file() {
                qsf
            } else {
                job.workdir.join("case.toml")
            }
        };
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "Salmon",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);
        let mut artefacts: Vec<Artifact> = Vec::new();

        // Walk the output directory top-level. Salmon writes a
        // handful of files (`quant.sf`, `cmd_info.json`,
        // `lib_format_counts.json`, `logs/`, plus `aux_info/`); we
        // surface the two canonical artifacts per the output spec.
        let entries = match fs::read_dir(&output_dir) {
            Ok(e) => e,
            Err(e) => {
                tracing::warn!(target: "valenx-salmon", ?e, "output_dir read failed");
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
                "quant.sf" => (
                    ArtifactKind::Tabular,
                    "Salmon transcript quantification".to_string(),
                ),
                "cmd_info.json" => (ArtifactKind::Log, "Salmon command info".to_string()),
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
            ribbon_contributions: vec!["bio.salmon.quant"],
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn info_is_bio_domain() {
        let info = SalmonAdapter::new().info();
        assert_eq!(info.id, "salmon");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "GPL-3.0");
        assert_eq!(info.display_name, "Salmon");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = SalmonAdapter::new().info();
        // Salmon 1.10.x is the de facto stable line; 2.0 reserves
        // room for an eventual major bump.
        assert_eq!(info.version_range.min_inclusive, Version::new(1, 10, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(2, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = SalmonAdapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.salmon.quant"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = SalmonAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }

    /// Round-9 RED→GREEN: `[bio.salmon].transcriptome` used to be
    /// joined with bare `case.path.join`. Wrap with `confined_join`.
    #[test]
    fn prepare_rejects_transcriptome_traversing_outside_case_dir() {
        use valenx_test_utils::tempdir;
        let d = tempdir("salmon-trans-trav");
        std::fs::write(d.join("r1.fq"), b"@x\nACGT\n+\nIIII\n").unwrap();
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "salmon.quant"

[bio.salmon]
transcriptome = "../../etc/passwd"
index_dir     = "idx"
reads         = ["r1.fq"]
output_dir    = "out"
threads       = 1
libtype       = "A"
"#,
        )
        .unwrap();
        let case = Case {
            id: "salmon-trans-trav".into(),
            path: d.clone(),
        };
        let workdir = d.join("workdir");
        let err = SalmonAdapter::new()
            .prepare(&case, &workdir)
            .unwrap_err();
        let msg = format!("{err}");
        assert!(
            msg.contains("..") || msg.contains("stay within") || msg.contains("escape"),
            "expected confined_join rejection, got: {msg}"
        );
        let _ = std::fs::remove_dir_all(&d);
    }

    /// Round-10 H3 RED→GREEN: `index_dir` flowed into
    /// `workdir.join(...)` with no validation. Hostile
    /// `index_dir = "../etc"` is now rejected by
    /// `validate_output_dir` (which allows multi-component subpaths
    /// but rejects `..` traversal + absolute paths).
    #[test]
    fn prepare_rejects_index_dir_path_traversal() {
        use valenx_test_utils::tempdir;
        let d = tempdir("salmon-idx-trav");
        std::fs::write(d.join("trans.fa"), b">x\nACGT\n").unwrap();
        std::fs::write(d.join("r1.fq"), b"@x\nACGT\n+\nIIII\n").unwrap();
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "salmon.quant"

[bio.salmon]
transcriptome = "trans.fa"
index_dir     = "../etc"
reads         = ["r1.fq"]
output_dir    = "out"
threads       = 1
libtype       = "A"
"#,
        )
        .unwrap();
        let case = Case {
            id: "trav".into(),
            path: d.clone(),
        };
        let workdir = d.join("workdir");
        let err = SalmonAdapter::new()
            .prepare(&case, &workdir)
            .unwrap_err();
        let msg = format!("{err}");
        assert!(
            msg.contains("[bio.salmon].index_dir"),
            "expected [bio.salmon].index_dir in error, got: {msg}"
        );
        let _ = std::fs::remove_dir_all(&d);
    }
}
