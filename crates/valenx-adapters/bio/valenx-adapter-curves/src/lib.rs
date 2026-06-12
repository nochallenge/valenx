//! # valenx-adapter-curves
//!
//! Adapter for [Curves+](https://bisi.ibcp.fr/tools/curves_plus/) —
//! Richard Lavery's reference toolkit for DNA helical-axis analysis.
//! Curves+ fits a curvilinear helical axis through a nucleic-acid
//! structure and reports per-base axis-curvature, base-pair
//! parameters relative to that axis, and a `.cda` file describing
//! the axis itself for downstream visualisation. It is the canonical
//! tool for "is this DNA bent, and if so, how" questions in
//! protein-DNA / drug-DNA structural studies.
//!
//! **Phase 39 — subprocess wrapper around `Cur+` with stdin-piped
//! parameters.** Curves+ takes its parameters as a Fortran-style
//! `&inp ... &end` namelist block on stdin followed by strand /
//! axis residue cards. We script those via the **CTFFIND-style
//! stdin-feed pattern**: `prepare()` writes the parameter body to
//! `curves_params.txt` in the workdir, stashes the filename under
//! `VALENX_CURVES_PARAMS_FILE`, and `run()` opens the file and
//! pipes its contents into `Cur+`'s stdin via `Stdio::from(file)`.
//!
//! ## Why a custom run() instead of `subprocess::run`
//!
//! The shared `subprocess::run` helper closes the child's stdin
//! (`Stdio::null()`). Curves+ on a closed stdin reads EOF before
//! parsing its first parameter and exits with an error. The custom
//! run() opens the parameters file with `File::open()` and hands
//! its FD to the child via `Stdio::from(file)` — Curves+ sees a
//! pipe pre-loaded with the namelist body and parses it as if a
//! human had typed it.
//!
//! ## License flag
//!
//! Curves+ ships under a custom non-OSS license that restricts use
//! to non-commercial / academic contexts. We surface this via
//! `tool_license = "Curves-License"` and emit a probe warning when
//! the binary is found, with the literal string `"academic"` as a
//! stable anchor for tests and downstream license-aware filters.
//!
//! On `collect()` we walk the workdir for `<output_basename>*.lis`
//! (the helical-analysis log) and `<output_basename>*.cda` (the
//! axis-curve data file).

#![forbid(unsafe_code)]
#![allow(missing_docs)]

pub mod case_input;

use std::ffi::OsString;
use std::fs::{self, File};
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

use crate::case_input::CurvesInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(CurvesAdapter::new())
}

pub struct CurvesAdapter;

impl CurvesAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for CurvesAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "curves";
/// Curves+'s binary candidate. The canonical entry point is the
/// short `Cur+` name shipped with the Curves+ distribution.
const BINARIES: &[&str] = &["Cur+"];

/// Filename of the parameters file the adapter writes into the
/// workdir during `prepare()`. Pinned so the env var, the file
/// write in `prepare()`, and the stdin redirect in `run()` all
/// agree.
const PARAMS_FILENAME: &str = "curves_params.txt";

/// Environment-variable name `prepare()` stashes the parameters
/// filename under so `run()` can recover it. We strip the var
/// before spawning the child so Curves+ doesn't see it.
const PARAMS_ENV_VAR: &str = "VALENX_CURVES_PARAMS_FILE";

/// The probe-warning surfaced whenever Curves+ is detected. Anchors
/// a stable "academic / non-commercial only" reminder for downstream
/// tooling and tests; the literal string `"academic"` is part of
/// the asserted contract.
const LICENSE_WARNING: &str = "Curves+ is licensed for non-commercial / academic use only. \
     Confirm your use case complies with the Curves+ license before \
     redistributing analyses or derived data.";

impl Adapter for CurvesAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "Curves+",
            // Curves+ 2.x is the modern stable line; 2.0 is the
            // floor we test against. Upper bound 3.0 reserves room
            // for an eventual major bump.
            version_range: VersionRange {
                min_inclusive: Version::new(2, 0, 0),
                max_exclusive: Version::new(3, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            // Curves+'s terms aren't a recognised SPDX identifier;
            // the closest accurate label is the project's own
            // "Curves-License" name.
            tool_license: "Curves-License",
            docs_url: "https://bisi.ibcp.fr/tools/curves_plus/",
            homepage_url: "https://bisi.ibcp.fr/tools/curves_plus/",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(BINARIES) {
            Some(binary_path) => {
                // `Cur+` (no args / closed stdin) prints a banner
                // before exiting; the generic detector tries the
                // conventional flags then falls back to a bare-name
                // scan.
                let found_version = detect_tool_version_semver(&binary_path, &["--version", ""]);
                Ok(ProbeReport {
                    ok: true,
                    found_version,
                    binary_path: Some(binary_path),
                    // Always surface the license reminder when
                    // Curves+ is detected — non-OSS academic use
                    // only.
                    warnings: vec![LICENSE_WARNING.to_string()],
                    required_env: Vec::new(),
                })
            }
            None => Err(AdapterError::ToolNotInstalled {
                name: INFO_ID,
                hint: "Curves+ 2.0+ required; download from \
                       https://bisi.ibcp.fr/tools/curves_plus/ \
                       (registration required, academic-use license) \
                       and ensure `Cur+` is on PATH"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = CurvesInput::from_case_dir(&case.path)?;

        // Round-4 security: reject `output_basename = "../etc/passwd"`
        // and friends before the value flows into any path join.
        // Same pattern as the round-3 fix in bionetgen/iqtree/art/fasttree.
        valenx_core::adapter_helpers::validate_output_basename(
            &input.output_basename,
            "[bio.curves].output_basename",
        )
        .map_err(|e| AdapterError::InvalidCase {
            case_path: case.path.join("case.toml"),
            reason: format!("{e}"),
        })?;

        fs::create_dir_all(workdir)?;

        // Resolve the input PDB against the case directory if
        // relative.
        let source_pdb = if input.input_pdb.is_absolute() {
            input.input_pdb.clone()
        } else {
            valenx_core::adapter_helpers::confined_join(&case.path, &input.input_pdb)?
        };
        if !source_pdb.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.curves].input_pdb `{}` not found (resolved {})",
                    input.input_pdb.display(),
                    source_pdb.display()
                ),
            });
        }

        let binary_path = find_on_path(BINARIES).ok_or_else(|| AdapterError::ToolNotInstalled {
            name: INFO_ID,
            hint: "Curves+ 2.0+ required; download from \
                       https://bisi.ibcp.fr/tools/curves_plus/ \
                       (registration required, academic-use license) \
                       and ensure `Cur+` is on PATH"
                .into(),
        })?;

        // Write the namelist body Curves+ expects on stdin. The
        // shape mirrors the canonical Curves+ tutorial input — an
        // `&inp` block naming the PDB and output basename, then
        // strand cards giving Curves+ the residue range to analyse.
        let params_path = workdir.join(PARAMS_FILENAME);
        let strand_length = input
            .last_residue
            .saturating_sub(input.first_residue)
            .saturating_add(1);
        let params_body = build_params_body(
            &source_pdb,
            &input.output_basename,
            strand_length,
            input.first_residue,
            input.last_residue,
        );
        valenx_core::io_caps::atomic_write_bytes(&params_path, params_body.as_bytes()).map_err(
            |e| AdapterError::Other(anyhow::anyhow!("write {} body: {e}", params_path.display())),
        )?;

        // Stash the params filename under a sentinel env var so
        // run() can recover it. The custom run() strips this var
        // before spawning so Curves+ never sees it — the env table
        // here is purely an adapter-internal channel.
        let mut native_command: Vec<OsString> = vec![binary_path.into_os_string()];
        for arg in &input.extra_args {
            native_command.push(OsString::from(arg));
        }

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: vec![(
                OsString::from(PARAMS_ENV_VAR),
                OsString::from(PARAMS_FILENAME),
            )],
            // Curves+ on a single structure runs in seconds; long
            // tail is multi-model NMR ensembles, which still finish
            // well inside an hour.
            estimated_runtime: Some(Duration::from_secs(60 * 60)),
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

        // Recover the params filename from the sentinel env var, then
        // strip it from the env table so Curves+ doesn't see it.
        let mut params_filename: Option<OsString> = None;
        let mut filtered_env: Vec<(OsString, OsString)> = Vec::with_capacity(job.environment.len());
        for (k, v) in &job.environment {
            if k == OsString::from(PARAMS_ENV_VAR).as_os_str() {
                params_filename = Some(v.clone());
            } else {
                filtered_env.push((k.clone(), v.clone()));
            }
        }
        let params_filename = params_filename.ok_or_else(|| {
            AdapterError::Other(anyhow::anyhow!(
                "missing {PARAMS_ENV_VAR} in PreparedJob.environment — \
                 prepare() should have populated it"
            ))
        })?;
        let params_path = job.workdir.join(&params_filename);
        let params_file = File::open(&params_path).map_err(|e| {
            AdapterError::Other(anyhow::anyhow!(
                "open {} for read: {e}",
                params_path.display()
            ))
        })?;

        let program = &job.native_command[0];
        let args: Vec<&OsString> = job.native_command.iter().skip(1).collect();

        ctx.report_progress(0.0, "starting Curves+");
        ctx.log(
            LogLevel::Info,
            &format!(
                "spawning {} with {} arg(s) in {}",
                program.to_string_lossy(),
                args.len(),
                job.workdir.display()
            ),
        );

        let mut cmd = Command::new(program);
        for a in &args {
            cmd.arg(a);
        }
        for (k, v) in &filtered_env {
            cmd.env(k, v);
        }
        cmd.current_dir(&job.workdir)
            .stdin(Stdio::from(params_file))
            .stdout(Stdio::piped())
            .stderr(Stdio::piped());

        let raw_child = cmd.spawn().map_err(|e| AdapterError::Run {
            exit_code: -1,
            stderr: format!("failed to spawn {}: {e}", program.to_string_lossy()),
            phase: RunPhase::Startup,
        })?;
        // Round-24 H2: KillOnDropChild guard.
        let mut kill_guard = KillOnDropChild::new(raw_child, true);

        // Curves+ chats on stdout (per-step "Reading PDB" / "Helical
        // axis" lines). Stderr carries the rare hard errors. Drain
        // both via channel threads so a slow consumer doesn't
        // deadlock the child.
        let stdout = kill_guard
            .inner_mut()
            .stdout
            .take()
            .ok_or_else(|| AdapterError::Other(anyhow::anyhow!("child stdout not captured")))?;
        let stderr = kill_guard
            .inner_mut()
            .stderr
            .take()
            .ok_or_else(|| AdapterError::Other(anyhow::anyhow!("child stderr not captured")))?;
        // Round-24 H2: bounded sync_channel + capped lines.
        let (tx, rx) = mpsc::sync_channel::<(LogLevel, String)>(SUBPROCESS_CHANNEL_CAPACITY);
        let tx_out = tx.clone();
        let so_thread = thread::spawn(move || {
            let reader = BufReader::new(stdout);
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
                if tx_out.send((LogLevel::Info, s)).is_err() {
                    break;
                }
            }
        });
        let tx_err = tx;
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
                if tx_err.send((LogLevel::Warn, s)).is_err() {
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
                Ok((level, line)) => {
                    process_line(level, &line, ctx, &mut warnings);
                    if level == LogLevel::Warn {
                        if stderr_tail.len() >= STDERR_TAIL_MAX {
                            stderr_tail.remove(0);
                        }
                        stderr_tail.push(line);
                    }
                }
                Err(mpsc::RecvTimeoutError::Timeout) => {
                    if let Some(status) = kill_guard
                        .inner_mut()
                        .try_wait()
                        .map_err(AdapterError::Io)?
                    {
                        // Drain any remaining lines so nothing's lost.
                        for (level, line) in rx.try_iter() {
                            process_line(level, &line, ctx, &mut warnings);
                            if level == LogLevel::Warn {
                                if stderr_tail.len() >= STDERR_TAIL_MAX {
                                    stderr_tail.remove(0);
                                }
                                stderr_tail.push(line);
                            }
                        }
                        let _ = so_thread.join();
                        let _ = se_thread.join();
                        return finalize(status, start.elapsed(), warnings, stderr_tail);
                    }
                }
                Err(mpsc::RecvTimeoutError::Disconnected) => {
                    let status = kill_guard.inner_mut().wait().map_err(AdapterError::Io)?;
                    let _ = so_thread.join();
                    let _ = se_thread.join();
                    return finalize(status, start.elapsed(), warnings, stderr_tail);
                }
            }
        }
    }

    fn collect(&self, job: &PreparedJob) -> Result<Results, AdapterError> {
        // Provenance: hash the staged parameters file as the
        // canonical input descriptor. Falls back to case.toml when
        // the params file isn't present (e.g. partial / failed run).
        let case_hash_input = {
            let p = job.workdir.join(PARAMS_FILENAME);
            if p.is_file() {
                p
            } else {
                job.workdir.join("case.toml")
            }
        };
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "Curves+",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);
        let mut artefacts: Vec<Artifact> = Vec::new();

        // Read the staged case.toml back out so we can restrict the
        // collected `.lis` / `.cda` outputs to those whose stem
        // starts with the configured `output_basename`. Failure to
        // read the case is non-fatal — we then accept every
        // candidate.
        let basename = read_output_basename(&job.workdir);

        let entries = match fs::read_dir(&job.workdir) {
            Ok(e) => e,
            Err(e) => {
                tracing::warn!(target: "valenx-curves", ?e, "workdir read failed");
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
                Some("lis") => (ArtifactKind::Log, "Curves+ helical analysis".to_string()),
                Some("cda") => (ArtifactKind::Tabular, "Curves+ axis curve data".to_string()),
                _ => continue,
            };
            let stem = path
                .file_stem()
                .and_then(|s| s.to_str())
                .unwrap_or("")
                .to_string();
            let stem_ok = match basename.as_deref() {
                Some(b) => stem.starts_with(b),
                None => true,
            };
            if !stem_ok {
                continue;
            }
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
            ribbon_contributions: vec!["bio.curves.analyze"],
        }
    }
}

/// Build the body of the parameters file we feed to Curves+ on
/// stdin. The shape mirrors the canonical Curves+ tutorial input —
/// an `&inp` namelist block naming the PDB and output basename,
/// then strand / axis residue cards giving Curves+ the residue
/// range to analyse.
fn build_params_body(
    input_pdb: &Path,
    output_basename: &str,
    strand_length: u32,
    first_residue: u32,
    last_residue: u32,
) -> String {
    let pdb_disp = input_pdb.display();
    let mut body = String::new();
    body.push_str(&format!(
        "&inp file={pdb_disp},\n     lis={output_basename},\n     ions=.t.,\n     test=.f.,\n     line=.f.,\n     fit=.t.,\n     axfrm=.t.,\n     ends=.f., &end\n"
    ));
    body.push_str("2 1 -1 0 0\n");
    body.push_str(&format!(
        "1 {strand_length} 0 {first_residue} {last_residue}\n"
    ));
    body
}

/// Per-line forwarder used by run()'s drain loop. Logs the line at
/// the right level and lifts Curves+'s progress markers to coarse
/// UI ticks; collects stderr lines into the warnings vector for
/// the RunReport.
fn process_line(level: LogLevel, line: &str, ctx: &mut RunContext<'_>, warnings: &mut Vec<String>) {
    ctx.log(level, line);
    if line.contains("Reading PDB") || line.contains("Reading file") {
        ctx.report_progress(20.0, line);
    } else if line.contains("Helical axis") || line.contains("Curve fit") {
        ctx.report_progress(60.0, line);
    } else if line.contains("Curves+ ended") || line.contains("complete") {
        ctx.report_progress(95.0, line);
    }
    if level == LogLevel::Warn || line.contains("ERROR") || line.contains("error:") {
        warnings.push(line.trim().to_string());
    }
}

/// Bridge a child's exit status into the adapter's RunReport / Run
/// error, mirroring `subprocess::finalize` so Curves+'s failure
/// mode matches every other Subprocess-mode adapter.
fn finalize(
    status: std::process::ExitStatus,
    wall_time: Duration,
    warnings: Vec<String>,
    stderr_tail: Vec<String>,
) -> Result<RunReport, AdapterError> {
    let exit_code = status.code().unwrap_or(-1);
    if !status.success() {
        let stderr = if stderr_tail.is_empty() {
            format!("Curves+ exited {exit_code} with no stderr output")
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

/// Pull `output_basename` out of the staged `case.toml` for
/// `collect()`-time stem filtering. Returns None if the file
/// doesn't exist or can't be parsed — collect falls back to
/// accepting every `.lis` / `.cda` in that case.
fn read_output_basename(workdir: &Path) -> Option<String> {
    let case_toml = workdir.join("case.toml");
    let text = valenx_core::io_caps::read_capped_to_string(
        &case_toml,
        valenx_core::project::loader::MAX_PROJECT_FILE_BYTES as usize,
    )
    .ok()?;
    let parsed: toml::Value = toml::from_str(&text).ok()?;
    parsed
        .get("bio")?
        .get("curves")?
        .get("output_basename")?
        .as_str()
        .map(|s| s.to_string())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn info_is_bio_domain() {
        let info = CurvesAdapter::new().info();
        assert_eq!(info.id, "curves");
        assert_eq!(info.physics, &[Physics::Bio]);
        // Curves+'s custom non-OSS license, not a recognised SPDX
        // identifier — pin the project's own label.
        assert_eq!(info.tool_license, "Curves-License");
        assert_eq!(info.display_name, "Curves+");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = CurvesAdapter::new().info();
        // Curves+ 2.0 is the floor; upper bound 3.0 reserves room
        // for the next major.
        assert_eq!(info.version_range.min_inclusive, Version::new(2, 0, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(3, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = CurvesAdapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.curves.analyze"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = CurvesAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }

    #[test]
    fn probe_warning_mentions_academic() {
        // The license-flag warning is mandatory: Curves+ is
        // non-OSS academic-use, and we surface that on every
        // successful probe. The literal "academic" anchor is what
        // downstream tooling and license-aware filters key off —
        // pin it.
        assert!(
            LICENSE_WARNING.contains("academic"),
            "probe warning must contain `academic` anchor; got: {LICENSE_WARNING}"
        );
    }
}
