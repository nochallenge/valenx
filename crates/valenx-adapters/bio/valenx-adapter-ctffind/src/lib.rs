//! # valenx-adapter-ctffind
//!
//! Adapter for [CTFFIND](https://grigoriefflab.umassmed.edu/ctffind4)
//! — Niko Grigorieff's contrast transfer function (CTF) estimation
//! tool. CTFFIND is the gold standard for fitting per-micrograph CTF
//! parameters (defocus, astigmatism, phase shift) in single-particle
//! cryo-EM workflows. RELION, cryoSPARC, EMAN2, and most automated
//! pipelines all wrap CTFFIND under the hood.
//!
//! **Phase 36 — subprocess wrapper around `ctffind` with stdin-piped
//! parameters.** CTFFIND's CLI is interactive: it prompts the user
//! line-by-line for each microscope parameter on startup. We script
//! the responses by writing them to a parameters text file in the
//! workdir and piping that file into the child's stdin. This mirrors
//! the MAFFT custom-run pattern, but with stdin redirect instead of
//! stdout redirect.
//!
//! ## Why a custom run() instead of `subprocess::run`
//!
//! The shared `subprocess::run` helper closes the child's stdin
//! (`Stdio::null()`). CTFFIND running on a closed stdin reads EOF
//! before its first prompt and exits with an error. Custom run()
//! opens the parameters file with `File::open()` and hands its FD
//! to the child via `Stdio::from(file)` — the child sees a pipe
//! pre-loaded with one parameter per prompt and responds to each
//! line as if a human had typed it.
//!
//! ## License flag
//!
//! CTFFIND is distributed under the Janelia Research Campus
//! non-commercial / academic-only license. We surface this accurately
//! via a `tool_license` value of `Janelia-License` and emit a probe
//! warning when the binary is found. The probe-warning text contains
//! the literal string `"academic"` as a stable anchor for tests and
//! downstream filters.
//!
//! On `collect()` we surface the diagnostic image (`output_diagnostic`)
//! as `Native` and the per-micrograph parameter file (`output_txt`)
//! as `Tabular`.

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

use crate::case_input::CtffindInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(CtffindAdapter::new())
}

pub struct CtffindAdapter;

impl CtffindAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for CtffindAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "ctffind";
/// CTFFIND's binary candidate. Both source and Bioconda installs use
/// the canonical `ctffind` name.
const BINARIES: &[&str] = &["ctffind"];

/// Filename of the parameters file the adapter writes into the workdir
/// during `prepare()`. Pinned so the env var, the file write in
/// `prepare()`, and the stdin redirect in `run()` all agree.
const PARAMS_FILENAME: &str = "ctffind_params.txt";

/// Environment-variable name `prepare()` stashes the parameters
/// filename under so `run()` can recover it. We strip the var before
/// spawning the child so CTFFIND doesn't see it.
const PARAMS_ENV_VAR: &str = "VALENX_CTFFIND_PARAMS_FILE";

/// The probe-warning surfaced whenever CTFFIND is detected. Anchors a
/// stable "academic / non-commercial only" reminder for downstream
/// tooling and tests; the literal string `"academic"` is part of the
/// asserted contract.
const LICENSE_WARNING: &str = "CTFFIND is licensed for non-commercial / academic use only. \
     Confirm your use case complies with the Janelia license before \
     redistributing CTF estimates or derived data.";

impl Adapter for CtffindAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "CTFFIND",
            // CTFFIND4 is the long-running stable line. 4.1 is the
            // floor we test against (the modern point-group); upper
            // bound 5.0 reserves room for the announced CTFFIND5
            // line.
            version_range: VersionRange {
                min_inclusive: Version::new(4, 1, 0),
                max_exclusive: Version::new(5, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            // CTFFIND's terms aren't a recognised SPDX identifier;
            // the closest accurate label is the project's own
            // Janelia name.
            tool_license: "Janelia-License",
            docs_url: "https://grigoriefflab.umassmed.edu/ctffind4",
            homepage_url: "https://grigoriefflab.umassmed.edu/ctffind4",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(BINARIES) {
            Some(binary_path) => {
                // CTFFIND prints its version on `--version` for
                // recent builds; older builds print it as part of
                // the interactive banner. The combined scanner
                // handles both.
                let found_version = detect_tool_version_semver(&binary_path, &["--version", ""]);
                Ok(ProbeReport {
                    ok: true,
                    found_version,
                    binary_path: Some(binary_path),
                    // Always surface the license reminder when
                    // CTFFIND is detected — non-OSS academic use
                    // only.
                    warnings: vec![LICENSE_WARNING.to_string()],
                    required_env: Vec::new(),
                })
            }
            None => Err(AdapterError::ToolNotInstalled {
                name: INFO_ID,
                hint: "CTFFIND 4.1+ required; download from \
                       https://grigoriefflab.umassmed.edu/ctffind4 \
                       (registration required, academic-use license)"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = CtffindInput::from_case_dir(&case.path)?;

        // Round-10 H3: `output_diagnostic` and `output_txt` are
        // `PathBuf` from case.toml and flowed into `workdir.join(…)`
        // with no validation. CTFFIND writes single files (not
        // directories) for each, so `validate_output_basename` is the
        // correct guard — rejects `..` traversal, absolute paths,
        // path separators.
        for (field, value) in [
            ("[bio.ctffind].output_diagnostic", &input.output_diagnostic),
            ("[bio.ctffind].output_txt", &input.output_txt),
        ] {
            if let Some(s) = value.to_str() {
                valenx_core::adapter_helpers::validate_output_basename(s, field).map_err(|e| {
                    AdapterError::InvalidCase {
                        case_path: case.path.join("case.toml"),
                        reason: format!("{e}"),
                    }
                })?;
            } else {
                return Err(AdapterError::InvalidCase {
                    case_path: case.path.join("case.toml"),
                    reason: format!("{field}: non-UTF-8 path rejected"),
                });
            }
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

        // Resolve the input micrograph against the case directory.
        let source_micrograph = if input.micrograph.is_absolute() {
            input.micrograph.clone()
        } else {
            valenx_core::adapter_helpers::confined_join(&case.path, &input.micrograph)?
        };
        if !source_micrograph.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.ctffind].micrograph `{}` not found (resolved {})",
                    input.micrograph.display(),
                    source_micrograph.display()
                ),
            });
        }

        // Output paths are scoped to the workdir if relative.
        let output_diagnostic = if input.output_diagnostic.is_absolute() {
            input.output_diagnostic.clone()
        } else {
            workdir.join(&input.output_diagnostic)
        };
        let output_txt = if input.output_txt.is_absolute() {
            input.output_txt.clone()
        } else {
            workdir.join(&input.output_txt)
        };

        let binary_path = find_on_path(BINARIES).ok_or_else(|| AdapterError::ToolNotInstalled {
            name: INFO_ID,
            hint: "CTFFIND 4.1+ required; download from \
                       https://grigoriefflab.umassmed.edu/ctffind4 \
                       (registration required, academic-use license)"
                .into(),
        })?;

        // Write the interactive-prompt script. CTFFIND's prompt
        // sequence (v4.1) is, in order: input image, output
        // diagnostic image, pixel size, voltage, Cs, amplitude
        // contrast, box size, min res, max res, min defocus, max
        // defocus, defocus step, find additional phase shift,
        // expert options, then a series of expert sub-prompts.
        // The output_txt path is derived by CTFFIND from the
        // diagnostic-image basename in v4.1, so we don't list it
        // separately here.
        let params_path = workdir.join(PARAMS_FILENAME);
        let params_body = build_params_body(
            &source_micrograph,
            &output_diagnostic,
            input.pixel_size,
            input.voltage,
            input.cs,
            input.amplitude_contrast,
        );
        valenx_core::io_caps::atomic_write_bytes(&params_path, params_body.as_bytes()).map_err(
            |e| AdapterError::Other(anyhow::anyhow!("write {} body: {e}", params_path.display())),
        )?;

        // Stash the params filename under a sentinel env var so
        // run() can recover it. The custom run() strips this var
        // before spawning so CTFFIND never sees it — the env table
        // here is purely an adapter-internal channel.
        let mut native_command: Vec<OsString> = vec![binary_path.into_os_string()];
        for arg in &input.extra_args {
            native_command.push(OsString::from(arg));
        }

        // Touch the output_txt path so collect() can find it even on
        // failed / partial runs. CTFFIND will overwrite it on success.
        let _ = File::create(&output_txt);

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: vec![(
                OsString::from(PARAMS_ENV_VAR),
                OsString::from(PARAMS_FILENAME),
            )],
            // CTFFIND on a single micrograph runs in seconds; on a
            // batched directory it scales with micrograph count.
            // 1 hour is generous for the long tail.
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
        // strip it from the env table so CTFFIND doesn't see it.
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

        ctx.report_progress(0.0, "starting CTFFIND");
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

        // CTFFIND chats on stdout (interactive prompt echoes plus the
        // per-stage "Estimating defocus..." progress lines). Stderr
        // carries the rare hard errors. Drain both via channel
        // threads so a slow consumer doesn't deadlock the child.
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
        // Provenance: hash the staged parameters file as the canonical
        // input descriptor. Falls back to case.toml when the params
        // file isn't present (e.g. partial / failed run).
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
            "CTFFIND",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);
        let mut artefacts: Vec<Artifact> = Vec::new();

        // Re-derive the output paths from the case so we report the
        // exact files CTFFIND wrote. Same resolution logic as
        // prepare(): absolute paths land where they say; relative
        // paths sit under the workdir.
        let case_path = job.workdir.join("case.toml");
        if case_path.is_file() {
            // case.toml is staged into the workdir at prepare time
            // by the executor; if it isn't here we fall back to a
            // workdir walk for the canonical filenames.
            if let Ok(input) = CtffindInput::from_case_dir(&job.workdir) {
                let output_diagnostic = if input.output_diagnostic.is_absolute() {
                    input.output_diagnostic.clone()
                } else {
                    job.workdir.join(&input.output_diagnostic)
                };
                let output_txt = if input.output_txt.is_absolute() {
                    input.output_txt.clone()
                } else {
                    job.workdir.join(&input.output_txt)
                };
                if output_diagnostic.is_file() {
                    artefacts.push(Artifact {
                        path: output_diagnostic,
                        kind: ArtifactKind::Native,
                        checksum: None,
                        label: "CTFFIND diagnostic image".to_string(),
                    });
                }
                if output_txt.is_file() {
                    artefacts.push(Artifact {
                        path: output_txt,
                        kind: ArtifactKind::Tabular,
                        checksum: None,
                        label: "CTFFIND parameters".to_string(),
                    });
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
            ribbon_contributions: vec!["bio.ctffind.estimate"],
        }
    }
}

/// Build the body of the parameters file we feed to CTFFIND on stdin.
/// Each line answers one of CTFFIND v4.1's interactive prompts in
/// order. The precise prompt sequence varies slightly between
/// builds; the values we supply are the conventional defaults from
/// single-particle cryo-EM workflows.
fn build_params_body(
    micrograph: &Path,
    output_diagnostic: &Path,
    pixel_size: f64,
    voltage: f64,
    cs: f64,
    amplitude_contrast: f64,
) -> String {
    let mut body = String::new();
    let micrograph_disp = micrograph.display();
    let output_diagnostic_disp = output_diagnostic.display();
    // Input image
    body.push_str(&format!("{micrograph_disp}\n"));
    // Output diagnostic image
    body.push_str(&format!("{output_diagnostic_disp}\n"));
    // Pixel size (Angstrom / px)
    body.push_str(&format!("{pixel_size}\n"));
    // Acceleration voltage (kV)
    body.push_str(&format!("{voltage}\n"));
    // Spherical aberration Cs (mm)
    body.push_str(&format!("{cs}\n"));
    // Amplitude contrast (fraction)
    body.push_str(&format!("{amplitude_contrast}\n"));
    // Box size for spectrum (px)
    body.push_str("512\n");
    // Minimum resolution (A)
    body.push_str("30.0\n");
    // Maximum resolution (A)
    body.push_str("5.0\n");
    // Minimum defocus to search (A)
    body.push_str("5000.0\n");
    // Maximum defocus to search (A)
    body.push_str("50000.0\n");
    // Defocus search step (A)
    body.push_str("500.0\n");
    // Find additional phase shift?
    body.push_str("no\n");
    // Expert options?
    body.push_str("yes\n");
    // Resample if pixel size > 1.5 A?
    body.push_str("no\n");
    // Known defocus range?
    body.push_str("no\n");
    // Use slower exhaustive search?
    body.push_str("yes\n");
    // Restrain astigmatism?
    body.push_str("no\n");
    // Find additional phase shift? (expert pass)
    body.push_str("no\n");
    body
}

/// Per-line forwarder used by run()'s drain loop. Logs the line at
/// the right level and lifts CTFFIND's progress markers to coarse
/// UI ticks; collects stderr lines into the warnings vector for
/// the RunReport.
fn process_line(level: LogLevel, line: &str, ctx: &mut RunContext<'_>, warnings: &mut Vec<String>) {
    ctx.log(level, line);
    if line.contains("Estimating defocus") || line.contains("Searching") {
        ctx.report_progress(40.0, line);
    } else if line.contains("Refining") {
        ctx.report_progress(70.0, line);
    } else if line.contains("Estimated values:") || line.contains("CTFFIND finished") {
        ctx.report_progress(95.0, line);
    }
    if level == LogLevel::Warn || line.contains("ERROR") || line.contains("error:") {
        warnings.push(line.trim().to_string());
    }
}

/// Bridge a child's exit status into the adapter's RunReport / Run
/// error, mirroring `subprocess::finalize` so CTFFIND's failure mode
/// matches every other Subprocess-mode adapter.
fn finalize(
    status: std::process::ExitStatus,
    wall_time: Duration,
    warnings: Vec<String>,
    stderr_tail: Vec<String>,
) -> Result<RunReport, AdapterError> {
    let exit_code = status.code().unwrap_or(-1);
    if !status.success() {
        let stderr = if stderr_tail.is_empty() {
            format!("CTFFIND exited {exit_code} with no stderr output")
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
        let info = CtffindAdapter::new().info();
        assert_eq!(info.id, "ctffind");
        assert_eq!(info.physics, &[Physics::Bio]);
        // The license identifier surfaces the Janelia non-OSS
        // license rather than mislabeling as MIT / BSD.
        assert_eq!(info.tool_license, "Janelia-License");
        assert_eq!(info.display_name, "CTFFIND");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = CtffindAdapter::new().info();
        // CTFFIND 4.1+ is the modern stable line; upper bound 5.0
        // reserves room for the announced CTFFIND5 line.
        assert_eq!(info.version_range.min_inclusive, Version::new(4, 1, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(5, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = CtffindAdapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.ctffind.estimate"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = CtffindAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }

    #[test]
    fn probe_warning_mentions_academic() {
        // The license-flag warning is mandatory: CTFFIND is non-OSS
        // academic-use, and we surface that on every successful
        // probe. The literal "academic" anchor is what downstream
        // tooling and license-aware filters key off — pin it.
        assert!(
            LICENSE_WARNING.contains("academic"),
            "probe warning must contain `academic` anchor; got: {LICENSE_WARNING}"
        );
    }

    /// Round-10 H3 RED→GREEN: `output_diagnostic` flowed into
    /// `workdir.join(...)` with no validation. Hostile
    /// `output_diagnostic = "../etc/passwd"` is now rejected.
    #[test]
    fn prepare_rejects_output_diagnostic_path_traversal() {
        use valenx_test_utils::tempdir;
        let d = tempdir("ctffind-output-trav");
        std::fs::write(d.join("mic.mrc"), b"FAKE").unwrap();
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "ctffind.estimate"

[bio.ctffind]
micrograph         = "mic.mrc"
output_diagnostic  = "../etc/passwd"
output_txt         = "ctf.txt"
pixel_size         = 1.0
voltage            = 300.0
cs                 = 2.7
amplitude_contrast = 0.07
"#,
        )
        .unwrap();
        let case = Case {
            id: "trav".into(),
            path: d.clone(),
        };
        let workdir = d.join("workdir");
        let err = CtffindAdapter::new().prepare(&case, &workdir).unwrap_err();
        let msg = format!("{err}");
        assert!(
            msg.contains("[bio.ctffind].output_diagnostic"),
            "expected [bio.ctffind].output_diagnostic in error, got: {msg}"
        );
        let _ = std::fs::remove_dir_all(&d);
    }

    /// RED→GREEN (round-24 H2): the drain loops in `run()` now use
    /// `read_capped_lines_bounded(reader, MAX_LINE_BYTES)` so a
    /// hostile child that floods stdout with bytes and never emits
    /// a `\n` can no longer OOM the drain thread. Test exercises the
    /// helper directly via the exact same shape the adapter uses
    /// (BufReader wrapping a Read source, MAX_LINE_BYTES cap).
    ///
    /// Pre-fix `BufReader::lines().map_while(Result::ok)` would
    /// allocate one giant `String` for the entire stream — at the
    /// adapter's 1 MiB cap this test would now allocate up to ~1 MiB
    /// before erroring; pre-fix it would have happily allocated the
    /// whole 4 MiB. Drop the cap to 64 KiB so the test stays small
    /// while still demonstrating the bound.
    #[test]
    fn drain_helper_caps_unbounded_line() {
        use std::io::Cursor;
        const TEST_CAP: usize = 64 * 1024;
        let payload = vec![b'x'; 256 * 1024]; // 4x cap, no newline
        let reader = BufReader::new(Cursor::new(payload));
        let mut errs = 0;
        let mut oks = 0;
        for line in read_capped_lines_bounded(reader, TEST_CAP) {
            match line {
                Ok(_) => oks += 1,
                Err(e) => {
                    assert_eq!(e.kind(), std::io::ErrorKind::InvalidData);
                    errs += 1;
                }
            }
        }
        assert_eq!(oks, 0, "no successful lines on flooded stream");
        assert_eq!(errs, 1, "exactly one cap error, then iteration ends");
    }

    /// RED→GREEN (round-24 H2): the sync_channel bound is shared
    /// across all 10 adapters. We import `SUBPROCESS_CHANNEL_CAPACITY`
    /// to compile-anchor the symbol so a future workspace-wide
    /// refactor that drops the constant breaks this test. Pre-fix
    /// `mpsc::channel()` had no bound, so the import didn't exist
    /// in this crate.
    #[test]
    fn sync_channel_bound_is_imported() {
        // Compile-time anchor: the const must be a positive integer
        // and the bounded sync_channel must accept it as a capacity.
        let cap = SUBPROCESS_CHANNEL_CAPACITY;
        assert!(cap > 0);
        let (tx, _rx) = std::sync::mpsc::sync_channel::<String>(cap);
        // sync_channel with bound=N accepts N sends without a recv.
        for _ in 0..cap {
            tx.try_send("ok".into()).expect("under bound");
        }
        // Next try_send fills the buffer.
        assert!(matches!(
            tx.try_send("overflow".into()),
            Err(std::sync::mpsc::TrySendError::Full(_))
        ));
    }
}
