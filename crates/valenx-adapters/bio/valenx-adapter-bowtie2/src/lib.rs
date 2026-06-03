//! # valenx-adapter-bowtie2
//!
//! Adapter for [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/)
//! — Langmead & Salzberg's gapped short-read aligner. Bowtie2 is the
//! second of the "big two" Illumina short-read aligners alongside BWA;
//! its FM-index design and end-to-end / local alignment modes make it
//! the workhorse for RNA-seq, ChIP-seq, and bisulfite-sequencing
//! pipelines.
//!
//! **Phase 18.5 — subprocess wrapper around `bowtie2`.** The user
//! supplies a reference FASTA plus 1 (single-end) or 2 (paired-end)
//! FASTQ files via `[bio.bowtie2]` in `case.toml`. `prepare()` builds
//! the FM-index next to the reference (`bowtie2-build`) unless
//! `skip_index = true`, then composes the `bowtie2` invocation.
//! `run()` streams the alignment via the shared subprocess runner;
//! Bowtie2 prints its summary stats on stderr at the end so the line
//! handler can lift the "overall alignment rate" line to a
//! near-completion progress hint.
//!
//! On `collect()` we report the canonical `out.sam` aligned-reads file
//! and any auxiliary log the user may have configured.

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

use crate::case_input::Bowtie2Input;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(Bowtie2Adapter::new())
}

pub struct Bowtie2Adapter;

impl Bowtie2Adapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for Bowtie2Adapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "bowtie2";
/// Bowtie2's binary candidates. `bowtie2` is the canonical Linux /
/// macOS install name from Bioconda, Homebrew, and source builds.
const BINARIES: &[&str] = &["bowtie2"];

/// The aligned-reads filename we tell `bowtie2 -S` to write. Pinned so
/// the `prepare()` invocation, the `collect()` walk, and the artifact
/// label all agree on what to look for.
const OUT_SAM: &str = "out.sam";

impl Adapter for Bowtie2Adapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "Bowtie2",
            // Bowtie2 2.5.x is the current stable line that every
            // distro ships; 3.0 reserves room for an eventual major
            // bump.
            version_range: VersionRange {
                min_inclusive: Version::new(2, 5, 0),
                max_exclusive: Version::new(3, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            tool_license: "GPL-3.0",
            docs_url: "https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml",
            homepage_url: "https://bowtie-bio.sourceforge.net/bowtie2/",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(BINARIES) {
            Some(binary_path) => {
                // `bowtie2 --version` prints the version on stdout;
                // the combined scanner picks it up cleanly.
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
                hint: "Bowtie2 2.5+ required; install via `apt install bowtie2`, \
                       `brew install bowtie2`, or `conda install -c bioconda bowtie2`"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = Bowtie2Input::from_case_dir(&case.path)?;

        fs::create_dir_all(workdir)?;

        // Resolve the reference path against the case directory if
        // relative.
        let source_reference = if input.reference.is_absolute() {
            input.reference.clone()
        } else {
            valenx_core::adapter_helpers::confined_join(
            &case.path,
            &input.reference,
        )?
        };
        if !source_reference.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.bowtie2].reference `{}` not found (resolved {})",
                    input.reference.display(),
                    source_reference.display()
                ),
            });
        }

        // Resolve each read file against the case directory via
        // `confined_join` — rejects absolute paths and `..`
        // traversal out of the case sandbox (round-6 hardening; the
        // pre-round-6 `case.path.join(read)` accepted both). Same
        // policy as BWA — we *don't* copy them into the workdir
        // because short-read FASTQs routinely hit tens of GB.
        let mut resolved_reads: Vec<PathBuf> = Vec::with_capacity(input.reads.len());
        for read in &input.reads {
            let resolved = confined_join(&case.path, read)?;
            if !resolved.is_file() {
                return Err(AdapterError::InvalidCase {
                    case_path: case.path.join("case.toml"),
                    reason: format!(
                        "[bio.bowtie2].reads entry `{}` not found (resolved {})",
                        read.display(),
                        resolved.display()
                    ),
                });
            }
            resolved_reads.push(resolved);
        }

        let binary_path = find_on_path(BINARIES).ok_or_else(|| AdapterError::ToolNotInstalled {
            name: INFO_ID,
            hint: "Bowtie2 2.5+ required; install via `apt install bowtie2`, \
                       `brew install bowtie2`, or `conda install -c bioconda bowtie2`"
                .into(),
        })?;

        // The Bowtie2 index basename is conventionally the reference
        // filename without its extension — `ref.fa` -> `ref`. The
        // index files themselves (`<base>.1.bt2`, `.2.bt2`, `.3.bt2`,
        // `.4.bt2`, `.rev.1.bt2`, `.rev.2.bt2`) sit in the workdir
        // (when we build) or next to the reference (when the user
        // pre-built one and set `skip_index = true`).
        let index_basename: String = source_reference
            .file_stem()
            .and_then(|s| s.to_str())
            .map(|s| s.to_string())
            .unwrap_or_else(|| "ref".to_string());

        // Build the FM-index in the workdir unless the user opted
        // out. We invoke `bowtie2-build` synchronously here so the
        // subsequent `run()` is a single `bowtie2` call, which lets
        // the shared subprocess runner stream stderr line-by-line
        // without chaining commands through a shell.
        if !input.skip_index {
            let index_status = std::process::Command::new("bowtie2-build")
                .arg(&source_reference)
                .arg(&index_basename)
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
                        "`bowtie2-build {} {}` failed (exit {}): {}",
                        source_reference.display(),
                        index_basename,
                        out.status.code().unwrap_or(-1),
                        stderr.lines().next().unwrap_or("(no stderr)")
                    )));
                }
                Err(e) => {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "spawning `bowtie2-build {}` failed: {e}",
                        source_reference.display()
                    )));
                }
            }
        }

        // Compose the alignment invocation:
        //   bowtie2 -x <index_base> --<preset> -p <threads> -S out.sam
        //           (-U <reads[0]>) | (-1 <reads[0]> -2 <reads[1]>)
        //           [extras...]
        //
        // Bowtie2 dispatches between single-end (`-U`) and paired-end
        // (`-1` / `-2`) based on which flag the reads come in on.
        let preset_flag = format!("--{}", input.preset);
        let mut native_command: Vec<OsString> = vec![
            binary_path.into_os_string(),
            OsString::from("-x"),
            OsString::from(&index_basename),
            OsString::from(preset_flag),
            OsString::from("-p"),
            OsString::from(input.threads.to_string()),
            OsString::from("-S"),
            OsString::from(OUT_SAM),
        ];
        if resolved_reads.len() == 1 {
            native_command.push(OsString::from("-U"));
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
            // Whole-exome runs finish in minutes; whole-genome runs
            // run for hours on a single node. 4 hours mirrors BWA's
            // generous default.
            estimated_runtime: Some(Duration::from_secs(4 * 60 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting Bowtie2", |line| {
            let mut hint = subprocess::Hint::default();
            // Bowtie2 emits its summary stats at the end of the run on
            // stderr — "N reads; of these:", "% overall alignment
            // rate", etc. The "overall alignment rate" line is the
            // very last useful marker so we pin it at 95%.
            if line.contains("overall alignment rate") {
                hint.progress = Some((95.0, line.to_string()));
            } else if line.contains("reads; of these:") {
                hint.progress = Some((50.0, line.to_string()));
            } else if line.contains("Error") || line.contains("ERR:") {
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
        // Provenance: hash the staged out.sam if present; fall back to
        // case.toml when the alignment hasn't produced a SAM yet so
        // the provenance block stays well-formed for partial / failed
        // runs.
        let case_hash_input = {
            let sam = job.workdir.join(OUT_SAM);
            if sam.is_file() {
                sam
            } else {
                job.workdir.join("case.toml")
            }
        };
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "Bowtie2",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);
        let mut artefacts: Vec<Artifact> = Vec::new();

        // Walk the workdir top-level. Bowtie2 only writes `out.sam`
        // here; the FM-index (`<base>.*.bt2`) and any logs the user
        // redirected stderr to also live here, so we surface them too.
        let entries = match fs::read_dir(&job.workdir) {
            Ok(e) => e,
            Err(e) => {
                tracing::warn!(target: "valenx-bowtie2", ?e, "workdir read failed");
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
            let (kind, label) = match ext.as_deref() {
                // `out.sam` — the aligned-reads file. Tabular per
                // output spec.
                Some("sam") => (ArtifactKind::Tabular, "Bowtie2 aligned reads".to_string()),
                Some("log") => (ArtifactKind::Log, "Bowtie2 log".to_string()),
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
            ribbon_contributions: vec!["bio.bowtie2.align"],
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn info_is_bio_domain() {
        let info = Bowtie2Adapter::new().info();
        assert_eq!(info.id, "bowtie2");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "GPL-3.0");
        assert_eq!(info.display_name, "Bowtie2");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = Bowtie2Adapter::new().info();
        // 2.5.x is the de facto stable line; 3.0 reserves room for an
        // eventual major bump.
        assert_eq!(info.version_range.min_inclusive, Version::new(2, 5, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(3, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = Bowtie2Adapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.bowtie2.align"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = Bowtie2Adapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }

    #[test]
    fn prepare_rejects_absolute_reads_entry() {
        // Round-6 RED→GREEN: the per-read resolution loop in `prepare()`
        // now goes through `confined_join`, which refuses any
        // `[bio.bowtie2].reads` entry that's absolute or escapes the
        // case directory. Pre-fix, `case.path.join("/etc/passwd")`
        // returned `/etc/passwd` on POSIX and let a hostile case bundle
        // hand the bowtie2 subprocess a path well outside the
        // case sandbox.
        use valenx_test_utils::tempdir;
        let d = tempdir("bowtie2-traversal");
        // A plausible-looking reference file in the case dir so the
        // reference resolution passes; the test targets the reads
        // loop specifically.
        std::fs::write(d.join("ref.fa"), ">chr1\nACGT\n").unwrap();
        // The poisoned reads entry: POSIX-absolute path that pre-fix
        // would be accepted verbatim by `case.path.join`.
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "bowtie2.align"

[bio.bowtie2]
reference = "ref.fa"
reads     = ["/etc/passwd"]
"#,
        )
        .unwrap();
        let case = Case {
            id: "bowtie2-traversal".into(),
            path: d.clone(),
        };
        let workdir = d.join("workdir");
        let err = Bowtie2Adapter::new().prepare(&case, &workdir).unwrap_err();
        let msg = format!("{err}");
        // `confined_join` rejects absolute paths with a message
        // that mentions the absolute path. Either error variant is
        // acceptable here; we just need to see the absolute path
        // refused.
        assert!(
            msg.contains("absolute") || msg.contains("/etc/passwd") || msg.contains("escape"),
            "expected confined_join rejection, got: {msg}"
        );
        let _ = std::fs::remove_dir_all(&d);
    }
}
