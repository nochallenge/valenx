//! # valenx-adapter-fasttree
//!
//! Adapter for [FastTree](http://www.microbesonline.org/fasttree/) —
//! Morgan Price's approximate-maximum-likelihood phylogenetics tool,
//! the canonical "fast tree on a big alignment" workhorse. Where
//! IQ-TREE and RAxML-NG search for the global ML tree (slow but
//! rigorous), FastTree builds an initial tree by neighbour-joining
//! and locally improves it with NNI / SPR moves until convergence.
//! On a 50k-sequence alignment FastTree finishes in minutes versus
//! hours for the full-ML tools.
//!
//! **Phase 30 — subprocess wrapper around `FastTree` / `fasttree`.**
//! The user supplies a multi-FASTA / PHYLIP alignment plus an
//! `output` path via `[bio.fasttree]` in `case.toml`. `prepare()`
//! resolves the input against the case directory and composes the
//! invocation:
//!
//! - nucleotide: `<binary> -nt [-gtr] [-gamma] <alignment>` → stdout
//! - amino-acid: `<binary> [-gamma] <alignment>` → stdout
//!
//! FastTree writes the Newick tree to **stdout** with no `-o` flag.
//! `run()` follows the same stdout-redirect pattern MAFFT uses —
//! spawn the child directly, hand its stdout to a `File`, and let
//! stderr carry the progress chatter line-by-line. `collect()`
//! surfaces the redirected output as the canonical Newick artifact.

#![forbid(unsafe_code)]
#![allow(missing_docs)]

pub mod case_input;
pub mod native;

use std::ffi::OsString;
use std::fs;
use std::io::BufReader;
use std::path::Path;
use std::process::{Command, Stdio};
use std::sync::mpsc;
use std::thread;
use std::time::{Duration, Instant};

use semver::Version;

use valenx_core::{
    adapter::LogLevel,
    adapter_helpers::{detect_tool_version_semver, find_on_path, live_provenance},
    error::RunPhase,
    io_caps::read_capped_lines_bounded,
    subprocess::{KillOnDropChild, MAX_LINE_BYTES, SUBPROCESS_CHANNEL_CAPACITY},
    Adapter, AdapterError, AdapterInfo, Capabilities, Case, LicenseMode, Physics, PreparedJob,
    ProbeReport, RunContext, RunReport, VersionRange,
};
use valenx_fields::{
    artifact::{Artifact, ArtifactKind},
    Results,
};

use crate::case_input::FastTreeInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(FastTreeAdapter::new())
}

pub struct FastTreeAdapter;

impl FastTreeAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for FastTreeAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "fasttree";
/// FastTree's binary name varies by distro: most ship the upstream
/// `FastTree` (capital F, no hyphen) name, but some Debian/Ubuntu
/// packages and Bioconda's older recipes install lowercase
/// `fasttree`. Probe both, prefer the upstream casing.
const BINARIES: &[&str] = &["FastTree", "fasttree"];

impl Adapter for FastTreeAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "FastTree",
            // FastTree's 2.1.x line has been the stable series for
            // over a decade (2.1.11 is the most recent point release
            // as of 2018). Floor at 2.1 covers every modern install;
            // upper bound 3.0 reserves room for a hypothetical major.
            version_range: VersionRange {
                min_inclusive: Version::new(2, 1, 0),
                max_exclusive: Version::new(3, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            tool_license: "GPL-2.0",
            docs_url: "http://www.microbesonline.org/fasttree/",
            homepage_url: "http://www.microbesonline.org/fasttree/",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(BINARIES) {
            Some(binary_path) => {
                // FastTree with no args prints "FastTree version X.Y.Z"
                // on stderr and a usage banner; the generic detector's
                // combined-stream scan picks up the leading SemVer.
                let found_version = detect_tool_version_semver(&binary_path, &[""]);
                Ok(ProbeReport {
                    ok: true,
                    found_version,
                    binary_path: Some(binary_path),
                    warnings: Vec::new(),
                    required_env: Vec::new(),
                })
            }
            // Native Rust fallback: BIONJ + NNI/SPR ML via valenx-phylo.
            None => Ok(ProbeReport {
                ok: true,
                found_version: None,
                binary_path: None,
                warnings: vec![
                    "FastTree binary not found; using native Rust BIONJ + ML topology \
                     optimisation (valenx-phylo). Install FastTree 2.1+ via apt/brew/conda \
                     for the full FastTree approximate-ML algorithm. Native mode uses JC69 \
                     for nucleotide (GTR not yet implemented); amino-acid inputs use \
                     p-distance + BIONJ only."
                        .to_string(),
                ],
                required_env: Vec::new(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = FastTreeInput::from_case_dir(&case.path)?;

        fs::create_dir_all(workdir)?;

        // Resolve the alignment against the case directory if relative.
        let source_alignment = if input.alignment.is_absolute() {
            input.alignment.clone()
        } else {
            valenx_core::adapter_helpers::confined_join(&case.path, &input.alignment)?
        };
        if !source_alignment.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.fasttree].alignment `{}` not found (resolved {})",
                    input.alignment.display(),
                    source_alignment.display()
                ),
            });
        }

        // Round-3 security fix: validate `output` is a single
        // path-component basename. It later becomes
        // `workdir.join(output)` in run() — without validation, a
        // hostile case.toml setting `output = "../../etc/cron.d/x"`
        // would let FastTree's stdout redirect write outside workdir.
        if let Some(name_str) = input.output.to_str() {
            valenx_core::adapter_helpers::validate_output_basename(
                name_str,
                "[bio.fasttree].output",
            )
            .map_err(|e| AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!("{e}"),
            })?;
        }

        // Write native_params.toml for the native path.
        let output_name = input
            .output
            .to_str()
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "output path is not valid UTF-8: {}",
                    input.output.display()
                ))
            })?
            .to_string();

        let native_params = native::NativeFasttreeParams {
            alignment_path: source_alignment
                .to_str()
                .ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "alignment path is not valid UTF-8: {}",
                        source_alignment.display()
                    ))
                })?
                .to_string(),
            output_name: output_name.clone(),
            seq_type: input.seq_type.clone(),
            ml_refine: true,
        };
        native::write_params(workdir, &native_params)?;

        // Stash the output path in the environment so run() can
        // resolve it without re-parsing the case. We use a custom env
        // var name (not exported to the child) — Valenx's env vector
        // is just adapter-internal scratch space here.
        let environment: Vec<(OsString, OsString)> = vec![(
            OsString::from("VALENX_FASTTREE_OUTPUT"),
            OsString::from(input.output.as_os_str()),
        )];

        // Compose the FastTree invocation. NO output flag — FastTree
        // writes Newick to stdout; the run() below redirects stdout
        // to `<workdir>/<output>`. Order:
        //   nucleotide: `-nt [-gtr] [-gamma] <alignment>`
        //   amino-acid: `[-gamma] <alignment>`
        //
        // Round-3 fix: extras MUST come after the positional
        // `<alignment>` so a hostile `extra_args = ["phantom"]` can't
        // shift the alignment onto a different argument slot.
        let native_command: Vec<OsString> = match find_on_path(BINARIES) {
            Some(binary_path) => {
                let mut cmd: Vec<OsString> = vec![binary_path.into_os_string()];
                if input.seq_type == "nt" {
                    cmd.push(OsString::from("-nt"));
                    if input.use_gtr {
                        cmd.push(OsString::from("-gtr"));
                    }
                }
                if input.gamma {
                    cmd.push(OsString::from("-gamma"));
                }
                cmd.push(source_alignment.into_os_string());
                for arg in &input.extra_args {
                    cmd.push(OsString::from(arg));
                }
                cmd
            }
            None => vec![OsString::from(native::NATIVE_SENTINEL)],
        };

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment,
            // FastTree on a 5k-sequence alignment finishes in seconds;
            // a 50k-sequence run can take an hour or two. 4 hours is
            // a generous default that still fails fast on stuck runs.
            estimated_runtime: Some(Duration::from_secs(4 * 60 * 60)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        if job.native_command.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "PreparedJob.native_command is empty — prepare() should \
                 have populated it"
            )));
        }

        // Native Rust path: BIONJ + ML topology inference, no subprocess.
        if job.native_command[0] == native::NATIVE_SENTINEL {
            return native::run_native(&job.workdir, ctx);
        }

        // Recover the output filename that prepare() stashed in the
        // environment vector; resolve against the workdir.
        let output_rel = job
            .environment
            .iter()
            .find(|(k, _)| k == "VALENX_FASTTREE_OUTPUT")
            .map(|(_, v)| v.clone())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "FastTree output filename missing from PreparedJob.environment — \
                     prepare() should have populated VALENX_FASTTREE_OUTPUT"
                ))
            })?;
        let out_path = job.workdir.join(std::path::PathBuf::from(&output_rel));
        let out_file = std::fs::File::create(&out_path).map_err(|e| {
            AdapterError::Other(anyhow::anyhow!(
                "open {} for write: {e}",
                out_path.display()
            ))
        })?;

        let program = &job.native_command[0];
        let args: Vec<&OsString> = job.native_command.iter().skip(1).collect();

        ctx.report_progress(0.0, "starting FastTree");
        ctx.log(
            LogLevel::Info,
            &format!(
                "spawning {} with {} arg(s) in {}; stdout -> {}",
                program.to_string_lossy(),
                args.len(),
                job.workdir.display(),
                out_path.display()
            ),
        );

        let mut cmd = Command::new(program);
        for a in &args {
            cmd.arg(a);
        }
        // Don't propagate the VALENX_FASTTREE_OUTPUT scratch var to
        // the child — FastTree doesn't read it, and it would only
        // pollute the child's environment.
        for (k, v) in &job.environment {
            if k == "VALENX_FASTTREE_OUTPUT" {
                continue;
            }
            cmd.env(k, v);
        }
        cmd.current_dir(&job.workdir)
            .stdin(Stdio::null())
            .stdout(Stdio::from(out_file))
            .stderr(Stdio::piped());

        let raw_child = cmd.spawn().map_err(|e| AdapterError::Run {
            exit_code: -1,
            stderr: format!("failed to spawn {}: {e}", program.to_string_lossy()),
            phase: RunPhase::Startup,
        })?;
        // Round-24 H2: KillOnDropChild guard.
        let mut kill_guard = KillOnDropChild::new(raw_child, true);

        // stdout is going to the file; only stderr needs a reader.
        let stderr = kill_guard
            .inner_mut()
            .stderr
            .take()
            .ok_or_else(|| AdapterError::Other(anyhow::anyhow!("child stderr not captured")))?;
        // Round-24 H2: bounded sync_channel + capped lines.
        let (tx, rx) = mpsc::sync_channel::<String>(SUBPROCESS_CHANNEL_CAPACITY);
        let se_thread = thread::spawn(move || {
            let reader = BufReader::new(stderr);
            for line in read_capped_lines_bounded(reader, MAX_LINE_BYTES) {
                let bytes = match line {
                    Ok(b) => b,
                    Err(_) => break,
                };
                let mut s = String::from_utf8_lossy(&bytes).into_owned();
                if s.ends_with('\n') {
                    s.pop();
                    if s.ends_with('\r') {
                        s.pop();
                    }
                }
                if tx.send(s).is_err() {
                    break;
                }
            }
        });

        let start = Instant::now();
        let mut warnings: Vec<String> = Vec::new();
        let mut stderr_tail: Vec<String> = Vec::new();
        const STDERR_TAIL_MAX: usize = 64;

        loop {
            if ctx.check_cancel().is_err() {
                let _ = kill_guard.inner_mut().kill();
                let _ = kill_guard.inner_mut().wait();
                return Err(AdapterError::Cancelled);
            }
            match rx.recv_timeout(Duration::from_millis(200)) {
                Ok(line) => {
                    process_stderr(&line, ctx, &mut warnings);
                    if stderr_tail.len() >= STDERR_TAIL_MAX {
                        stderr_tail.remove(0);
                    }
                    stderr_tail.push(line);
                }
                Err(mpsc::RecvTimeoutError::Timeout) => {
                    if let Some(status) = kill_guard
                        .inner_mut()
                        .try_wait()
                        .map_err(AdapterError::Io)?
                    {
                        // Drain remaining lines so nothing's lost.
                        for line in rx.try_iter() {
                            process_stderr(&line, ctx, &mut warnings);
                            if stderr_tail.len() >= STDERR_TAIL_MAX {
                                stderr_tail.remove(0);
                            }
                            stderr_tail.push(line);
                        }
                        let _ = se_thread.join();
                        return finalize(status, start.elapsed(), warnings, stderr_tail);
                    }
                }
                Err(mpsc::RecvTimeoutError::Disconnected) => {
                    let status = kill_guard.inner_mut().wait().map_err(AdapterError::Io)?;
                    let _ = se_thread.join();
                    return finalize(status, start.elapsed(), warnings, stderr_tail);
                }
            }
        }
    }

    fn collect(&self, job: &PreparedJob) -> Result<Results, AdapterError> {
        // Look up the output filename from the environment vector
        // (same channel prepare() used to hand it to run()).
        let output_rel = job
            .environment
            .iter()
            .find(|(k, _)| k == "VALENX_FASTTREE_OUTPUT")
            .map(|(_, v)| v.clone());
        let out_path = output_rel.map(|rel| job.workdir.join(std::path::PathBuf::from(&rel)));

        // Provenance: hash the staged Newick tree if present, falling
        // back to case.toml when the run hasn't produced one yet.
        let case_hash_input = match &out_path {
            Some(p) if p.is_file() => p.clone(),
            _ => job.workdir.join("case.toml"),
        };
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "FastTree",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);
        let mut artefacts: Vec<Artifact> = Vec::new();

        // Surface the redirected output as the canonical Newick
        // artifact. We only emit it if the file actually exists —
        // a partial / failed run may not have produced one.
        if let Some(p) = out_path {
            if p.is_file() {
                artefacts.push(Artifact {
                    path: p,
                    kind: ArtifactKind::Native,
                    checksum: None,
                    label: "FastTree Newick tree".to_string(),
                });
            }
        }
        artefacts.sort_by(|a, b| a.path.cmp(&b.path));
        results.artifacts = artefacts;
        Ok(results)
    }

    fn capabilities(&self) -> Capabilities {
        // Bio-specific Capability variants land in a follow-up task;
        // ribbon contributions are already enough for the registry
        // to surface the adapter.
        Capabilities {
            capabilities: Vec::new(),
            ribbon_contributions: vec!["bio.fasttree.tree"],
        }
    }
}

/// Mirror of MAFFT's stderr-line handler. Logs the line and lifts
/// FastTree's progress markers to coarse UI ticks.
fn process_stderr(line: &str, ctx: &mut RunContext<'_>, warnings: &mut Vec<String>) {
    ctx.log(LogLevel::Info, line);
    // FastTree writes progress on stderr: an "Initial topology"
    // marker after NJ, "Refining topology" during NNI / SPR rounds,
    // and a "Total time:" line when it finishes.
    if line.contains("Read alignment")
        || line.contains("Initial topology")
        || line.contains("Joined") && line.contains("of")
    {
        ctx.report_progress(20.0, line);
    } else if line.contains("Refining topology") || line.contains("ME NNIs") {
        ctx.report_progress(60.0, line);
    } else if line.contains("ML NNIs") || line.contains("Branch lengths") {
        ctx.report_progress(80.0, line);
    } else if line.contains("Total time:") || line.contains("TreeLogLk") {
        ctx.report_progress(95.0, line);
    }
    if line.contains("ERROR") || line.contains("Error:") || line.contains("WARNING") {
        warnings.push(line.trim().to_string());
    }
}

/// Bridge a child's exit status into RunReport / Run-error, mirroring
/// MAFFT's finalize so FastTree's failure mode matches every other
/// Subprocess-mode adapter.
fn finalize(
    status: std::process::ExitStatus,
    wall_time: Duration,
    warnings: Vec<String>,
    stderr_tail: Vec<String>,
) -> Result<RunReport, AdapterError> {
    let exit_code = status.code().unwrap_or(-1);
    if !status.success() {
        let stderr = if stderr_tail.is_empty() {
            format!("FastTree exited {exit_code} with no stderr output")
        } else {
            stderr_tail.join("\n")
        };
        return Err(AdapterError::Run {
            exit_code,
            stderr,
            phase: RunPhase::Solve,
        });
    }
    Ok(RunReport {
        exit_code,
        wall_time,
        converged: Some(true),
        residual_history: Vec::new(),
        warnings,
        final_phase: Some(RunPhase::Shutdown),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn info_is_bio_domain() {
        let info = FastTreeAdapter::new().info();
        assert_eq!(info.id, "fasttree");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "GPL-2.0");
        assert_eq!(info.display_name, "FastTree");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = FastTreeAdapter::new().info();
        // 2.1.x is the long-standing stable line; 3.0 reserves room
        // for a hypothetical future major bump.
        assert_eq!(info.version_range.min_inclusive, Version::new(2, 1, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(3, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = FastTreeAdapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.fasttree.tree"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = FastTreeAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }

    #[test]
    fn probe_always_succeeds() {
        let report = FastTreeAdapter::new().probe().unwrap();
        assert!(report.ok, "probe() must return ok=true in all modes");
    }

    #[test]
    fn native_sentinel_is_stable() {
        assert_eq!(native::NATIVE_SENTINEL, "valenx:native:fasttree");
    }
}
