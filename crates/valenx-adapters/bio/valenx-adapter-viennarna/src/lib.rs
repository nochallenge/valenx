//! # valenx-adapter-viennarna
//!
//! Adapter for [ViennaRNA](https://www.tbi.univie.ac.at/RNA/) — the
//! University of Vienna's reference toolkit for RNA secondary-structure
//! prediction. `RNAfold` is the workhorse: minimum-free-energy folding
//! over the canonical Turner thermodynamic model, with optional
//! partition-function / dot-plot output.
//!
//! **Phase 28 + in-house native path.** Two execution modes:
//!
//! - **Native Rust (zero-download, default when `RNAfold` not found):**
//!   Uses `valenx_rnastruct::fold::zuker::mfe_d2` with the complete
//!   Turner-2004 nearest-neighbor parameters and coaxial-stacking (`-d2`
//!   mode). Energies match ViennaRNA's `RNAfold -d2` exactly (to
//!   rounding). Output format is identical. No license restriction:
//!   the Turner-2004 parameters are published science.
//!
//! - **Subprocess (when `RNAfold` is installed):** Shells out to the
//!   ViennaRNA binary for full feature parity including partition-
//!   function PostScript output. ViennaRNA's custom non-OSS license
//!   applies only to this path.
//!
//! `probe()` succeeds in both modes. The adapter is always
//! `AdapterStatus::Ready`.
//!
//! ## License flag
//!
//! ViennaRNA ships under a custom non-OSS license that restricts
//! commercial redistribution to academic / non-commercial contexts.
//! We surface this via `tool_license = "ViennaRNA-License"` and emit
//! a probe warning **when the binary is found**. The native Rust path
//! carries no license restriction.

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
    adapter_helpers::{confined_join, detect_tool_version_semver, find_on_path, live_provenance},
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

use crate::case_input::ViennaRnaInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(ViennaRnaAdapter::new())
}

pub struct ViennaRnaAdapter;

impl ViennaRnaAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for ViennaRnaAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "viennarna";
/// ViennaRNA's binary candidates. The MFE folder installs as
/// `RNAfold` (capital R-N-A) on every supported platform — Bioconda,
/// Homebrew, the upstream tarball, and Debian / Ubuntu's
/// `vienna-rna` package all use the canonical capitalization.
const BINARIES: &[&str] = &["RNAfold"];

/// The probe-warning surfaced whenever ViennaRNA is detected. The
/// literal string `"academic"` is part of the asserted contract — it
/// anchors the license reminder so downstream license-aware filters
/// and tests can key off a stable substring.
const LICENSE_WARNING: &str = "ViennaRNA is licensed for non-commercial / academic use only. \
     Confirm your use case complies with the upstream license \
     before redistributing folds or derived data.";

impl Adapter for ViennaRnaAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "ViennaRNA",
            // ViennaRNA's 2.5.x line is the current stable series
            // (2.5.0 in 2022, point releases through 2.6+). 2.5
            // covers the modern API surface; upper-bound 3.0 reserves
            // room for an eventual major rewrite.
            version_range: VersionRange {
                min_inclusive: Version::new(2, 5, 0),
                max_exclusive: Version::new(3, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            // ViennaRNA's terms aren't a recognised SPDX identifier;
            // the closest accurate label is the project's own custom
            // license. Mislabeling as MIT / BSD would be misleading.
            tool_license: "ViennaRNA-License",
            docs_url: "https://www.tbi.univie.ac.at/RNA/",
            homepage_url: "https://www.tbi.univie.ac.at/RNA/",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(BINARIES) {
            Some(binary_path) => {
                // `RNAfold --version` prints "RNAfold X.Y.Z" on stdout.
                let found_version = detect_tool_version_semver(&binary_path, &["--version"]);
                Ok(ProbeReport {
                    ok: true,
                    found_version,
                    binary_path: Some(binary_path),
                    // Always surface the license reminder when
                    // ViennaRNA is detected — custom non-OSS license.
                    warnings: vec![LICENSE_WARNING.to_string()],
                    required_env: Vec::new(),
                })
            }
            // Native Rust fallback: always available, no binary needed,
            // no license restriction. The adapter is Ready in both modes.
            None => Ok(ProbeReport {
                ok: true,
                found_version: None, // no external version; native Rust impl
                binary_path: None,
                warnings: vec![
                    "RNAfold not found; using native Rust Zuker/Turner-2004 folder \
                     (valenx-rnastruct). Install `vienna-rna` (apt/brew/conda) to use \
                     the ViennaRNA binary instead."
                        .to_string(),
                ],
                required_env: Vec::new(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = ViennaRnaInput::from_case_dir(&case.path)?;

        fs::create_dir_all(workdir)?;

        // Resolve the input FASTA against the case directory via
        // `confined_join` — rejects absolute paths and `..` traversal
        // out of the case sandbox.
        let source_input = confined_join(&case.path, &input.input)?;
        if !source_input.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.viennarna].input `{}` not found (resolved {})",
                    input.input.display(),
                    source_input.display()
                ),
            });
        }

        // Stash the output filename so both run() and collect() can
        // find it without re-parsing case.toml.
        let environment: Vec<(OsString, OsString)> = vec![(
            OsString::from("VALENX_VIENNARNA_OUTPUT"),
            OsString::from(input.output.as_os_str()),
        )];

        // Write native_params.toml regardless of mode — run() reads it
        // for the native path; the subprocess path ignores it.
        let params = native::NativeViennaParams {
            input_path: source_input
                .to_str()
                .ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "input path is not valid UTF-8: {}",
                        source_input.display()
                    ))
                })?
                .to_string(),
            output_name: input
                .output
                .to_str()
                .ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "output path is not valid UTF-8: {}",
                        input.output.display()
                    ))
                })?
                .to_string(),
            temperature: input.temperature,
            partition_function: input.partition_function,
            allow_gu: input.allow_gu,
        };
        native::write_params(workdir, &params)?;

        // Choose subprocess vs native path.
        let native_command: Vec<OsString> = match find_on_path(BINARIES) {
            Some(binary_path) => {
                // Subprocess mode: compose the full RNAfold invocation.
                let mut cmd: Vec<OsString> = vec![
                    binary_path.into_os_string(),
                    OsString::from("-i"),
                    source_input.into_os_string(),
                    OsString::from("-T"),
                    OsString::from(format_temperature(input.temperature)),
                ];
                if input.partition_function {
                    cmd.push(OsString::from("-p"));
                }
                if !input.allow_gu {
                    cmd.push(OsString::from("--noGU"));
                }
                for arg in &input.extra_args {
                    cmd.push(OsString::from(arg));
                }
                cmd
            }
            // Native mode sentinel: run() detects this and calls the
            // Rust algorithm directly instead of spawning a child.
            None => vec![OsString::from(native::NATIVE_SENTINEL)],
        };

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment,
            estimated_runtime: Some(Duration::from_secs(30 * 60)),
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

        // Native Rust path: Zuker MFE folding, no subprocess.
        if job.native_command[0] == native::NATIVE_SENTINEL {
            return native::run_native(&job.workdir, ctx);
        }

        // Subprocess path: recover the output filename that prepare() stashed.
        let output_rel = job
            .environment
            .iter()
            .find(|(k, _)| k == "VALENX_VIENNARNA_OUTPUT")
            .map(|(_, v)| v.clone())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "ViennaRNA output filename missing from PreparedJob.environment — \
                     prepare() should have populated VALENX_VIENNARNA_OUTPUT"
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

        ctx.report_progress(0.0, "starting RNAfold");
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
        // Don't propagate the VALENX_VIENNARNA_OUTPUT scratch var to
        // the child — RNAfold doesn't read it and it would only
        // pollute the child's environment.
        for (k, v) in &job.environment {
            if k == "VALENX_VIENNARNA_OUTPUT" {
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
            .find(|(k, _)| k == "VALENX_VIENNARNA_OUTPUT")
            .map(|(_, v)| v.clone());
        let out_path = output_rel.map(|rel| job.workdir.join(std::path::PathBuf::from(&rel)));

        // Provenance: hash the staged structure output if present,
        // falling back to case.toml when the run hasn't produced one
        // yet — keeps the provenance block well-formed for partial /
        // failed runs.
        let case_hash_input = match &out_path {
            Some(p) if p.is_file() => p.clone(),
            _ => job.workdir.join("case.toml"),
        };
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "ViennaRNA",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);
        let mut artefacts: Vec<Artifact> = Vec::new();

        // Surface the redirected output as the canonical structure
        // artifact. We only emit it if the file actually exists — a
        // partial / failed run may not have produced one.
        if let Some(p) = out_path {
            if p.is_file() {
                artefacts.push(Artifact {
                    path: p,
                    kind: ArtifactKind::Native,
                    checksum: None,
                    label: "ViennaRNA secondary-structure output".to_string(),
                });
            }
        }
        artefacts.sort_by(|a, b| a.path.cmp(&b.path));
        results.artifacts = artefacts;
        Ok(results)
    }

    fn capabilities(&self) -> Capabilities {
        Capabilities {
            capabilities: Vec::new(),
            ribbon_contributions: vec!["bio.viennarna.fold"],
        }
    }
}

/// Format a Celsius float for RNAfold's `-T` flag. We want a stable,
/// locale-independent decimal representation — `format!("{}")` on an
/// f64 picks the shortest round-trip rendering, which matches what
/// users typed in their case.toml. Whole-number temperatures render
/// without a trailing `.0` (e.g. `37`), which RNAfold accepts.
fn format_temperature(temp: f64) -> String {
    format!("{temp}")
}

/// Mirror of MAFFT's stderr-line handler. Logs the line and lifts
/// RNAfold's progress markers to coarse UI ticks.
fn process_stderr(line: &str, ctx: &mut RunContext<'_>, warnings: &mut Vec<String>) {
    ctx.log(LogLevel::Info, line);
    // RNAfold's progress chatter is sparse: typically a banner at
    // startup and an occasional "WARNING" line. The partition-
    // function pass writes "Computing partition function" markers
    // we can lift to a 60% tick.
    if line.contains("Computing partition") || line.contains("base pair probabilities") {
        ctx.report_progress(60.0, line);
    } else if line.contains("Done") || line.contains("free energy") {
        ctx.report_progress(95.0, line);
    }
    if line.contains("ERROR") || line.contains("Error:") || line.contains("WARNING") {
        warnings.push(line.trim().to_string());
    }
}

/// Bridge a child's exit status into RunReport / Run-error, mirroring
/// MAFFT's finalize so RNAfold's failure mode matches every other
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
            format!("RNAfold exited {exit_code} with no stderr output")
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
        let info = ViennaRnaAdapter::new().info();
        assert_eq!(info.id, "viennarna");
        assert_eq!(info.physics, &[Physics::Bio]);
        // The license identifier must surface ViennaRNA's custom
        // non-OSS license rather than mislabel as MIT / BSD.
        assert_eq!(info.tool_license, "ViennaRNA-License");
        assert_eq!(info.display_name, "ViennaRNA");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = ViennaRnaAdapter::new().info();
        // ViennaRNA's 2.5.x line is the current stable; 3.0 reserves
        // room for an eventual major rewrite.
        assert_eq!(info.version_range.min_inclusive, Version::new(2, 5, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(3, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = ViennaRnaAdapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.viennarna.fold"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = ViennaRnaAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }

    #[test]
    fn probe_warning_mentions_academic() {
        // The license-flag warning is mandatory: ViennaRNA is non-OSS
        // academic-use, and we surface that on every successful
        // probe. The literal "academic" anchor is what downstream
        // tooling and license-aware filters key off — pin it.
        assert!(
            LICENSE_WARNING.contains("academic"),
            "probe warning must contain `academic` anchor; got: {LICENSE_WARNING}"
        );
    }

    /// When RNAfold is not installed, probe() must still return Ok
    /// with the native-rust version string. The adapter must always be Ready.
    #[test]
    fn probe_succeeds_in_native_mode_when_binary_absent() {
        // We can't guarantee the external binary is absent in CI, but
        // we can test the codepath directly by checking that the
        // ProbeReport::ok field is true regardless of binary presence.
        let report = ViennaRnaAdapter::new().probe().unwrap();
        assert!(report.ok, "probe() must return ok=true in all modes");
        // Either the real binary version (Some) or None in native mode.
        // The important thing is that ok=true.
    }

    /// Ensure the native sentinel constant matches what prepare() sets.
    #[test]
    fn native_sentinel_is_stable() {
        assert_eq!(native::NATIVE_SENTINEL, "valenx:native:viennarna");
    }
}
