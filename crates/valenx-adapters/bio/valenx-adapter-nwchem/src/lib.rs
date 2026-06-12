//! # valenx-adapter-nwchem
//!
//! Adapter for [NWChem](https://nwchemgit.github.io/) — the
//! Pacific-Northwest National Laboratory's massively-parallel ab
//! initio quantum chemistry suite. NWChem covers Hartree-Fock, DFT,
//! coupled-cluster, plane-wave DFT, classical molecular dynamics, and
//! QM/MM, with a domain-specific input language (`geometry`, `basis`,
//! `task` blocks in a `.nw` file). It scales to thousands of MPI
//! ranks on national HPC clusters.
//!
//! **Phase 25 — subprocess wrapper around the `nwchem` binary.** The
//! user supplies a `.nw` input file and an output path via
//! `[bio.nwchem]` in `case.toml`, plus an optional MPI process count.
//! `prepare()` resolves the input, picks `nwchem <input>` for serial
//! runs or `mpirun -n <N> nwchem <input>` for parallel ones, and
//! checks for `mpirun` on PATH when `mpi_procs > 1`.
//!
//! NWChem (like MAFFT) writes its run output to **stdout** rather
//! than to an `-o`-style flag, so we mirror MAFFT's pattern: a custom
//! `run()` that opens the output file, hands its FD to the child as
//! stdout, and lets the child write line-by-line straight to disk.
//!
//! On `collect()` we surface the staged output file as the run's
//! `Log` artifact.

#![forbid(unsafe_code)]
#![allow(missing_docs)]

pub mod case_input;

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

use crate::case_input::NwchemInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(NwchemAdapter::new())
}

pub struct NwchemAdapter;

impl NwchemAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for NwchemAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "nwchem";
/// NWChem's binary candidate.
const BINARIES: &[&str] = &["nwchem"];
/// `mpirun` is the de-facto cross-MPI-implementation launcher (OpenMPI,
/// MPICH, Intel MPI all alias it). We refuse parallel runs when it's
/// missing from PATH rather than silently dropping back to serial.
const MPI_LAUNCHER: &[&str] = &["mpirun"];

impl Adapter for NwchemAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "NWChem",
            // NWChem's 7.x line landed in 2020 and 7.2 (Dec 2022) is
            // the floor we test against; upper bound 8.0 reserves
            // room for an eventual major bump.
            version_range: VersionRange {
                min_inclusive: Version::new(7, 2, 0),
                max_exclusive: Version::new(8, 0, 0),
            },
            physics: &[Physics::Bio],
            license_mode: LicenseMode::Subprocess,
            tool_license: "ECL-2.0",
            docs_url: "https://nwchemgit.github.io/",
            homepage_url: "https://nwchemgit.github.io/",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(BINARIES) {
            Some(binary_path) => {
                // `nwchem -V` prints something like
                // "Northwest Computational Chemistry Package (NWChem) 7.2.0"
                // on stdout. The helper handles the rest.
                let found_version = detect_tool_version_semver(&binary_path, &["-V", ""]);
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
                hint: "NWChem 7.2+ required; install via \
                       `apt install nwchem`, `conda install -c conda-forge nwchem`, \
                       or build from source"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let input = NwchemInput::from_case_dir(&case.path)?;

        // Round-10 H3: `output` is `PathBuf` and pre-fix flowed into
        // `workdir.join(&input.output)`. Validate as a basename
        // before the join so `output = "../etc/passwd"` is rejected.
        if let Some(s) = input.output.to_str() {
            valenx_core::adapter_helpers::validate_output_basename(s, "[bio.nwchem].output")
                .map_err(|e| AdapterError::InvalidCase {
                    case_path: case.path.join("case.toml"),
                    reason: format!("{e}"),
                })?;
        } else {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: "[bio.nwchem].output: non-UTF-8 path rejected".into(),
            });
        }

        fs::create_dir_all(workdir)?;

        // Resolve the input `.nw` file against the case directory if
        // relative.
        let source_input = if input.input.is_absolute() {
            input.input.clone()
        } else {
            valenx_core::adapter_helpers::confined_join(&case.path, &input.input)?
        };
        if !source_input.is_file() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.nwchem].input `{}` not found (resolved {})",
                    input.input.display(),
                    source_input.display()
                ),
            });
        }

        let nwchem_path = find_on_path(BINARIES).ok_or_else(|| AdapterError::ToolNotInstalled {
            name: INFO_ID,
            hint: "NWChem 7.2+ required; install via \
                       `apt install nwchem`, `conda install -c conda-forge nwchem`, \
                       or build from source"
                .into(),
        })?;

        // Compose the invocation. Two shapes:
        //
        // 1. Serial:   `nwchem <input>`
        // 2. Parallel: `mpirun -n <N> nwchem <input>` — `mpirun` must
        //    be on PATH; if it isn't we fail at prepare() with a
        //    helpful install hint rather than letting the child fail
        //    later with a less obvious "command not found".
        let mut native_command: Vec<OsString> = Vec::new();
        if input.mpi_procs > 1 {
            let mpi = find_on_path(MPI_LAUNCHER).ok_or_else(|| AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "[bio.nwchem].mpi_procs = {} requires `mpirun` on PATH; \
                     install OpenMPI (`apt install openmpi-bin`) or MPICH \
                     (`apt install mpich`)",
                    input.mpi_procs
                ),
            })?;
            native_command.push(mpi.into_os_string());
            native_command.push(OsString::from("-n"));
            native_command.push(OsString::from(input.mpi_procs.to_string()));
            native_command.push(nwchem_path.into_os_string());
        } else {
            native_command.push(nwchem_path.into_os_string());
        }
        // Round-4 fix: extra_args after positionals (source_input) — see
        // security/code-review.md. Pre-fix the loop sat between the
        // nwchem binary and source_input, letting a hostile case.toml
        // supply e.g. extra_args = ["--help"] that swallows source_input.
        native_command.push(source_input.into_os_string());
        for arg in &input.extra_args {
            native_command.push(OsString::from(arg));
        }

        // Stash the resolved output path in the environment so run()
        // can find it without re-parsing the case TOML. We use a
        // workdir-relative resolution: relative paths land next to
        // case.toml in the workdir, absolute ones go where asked.
        let output_path = if input.output.is_absolute() {
            input.output.clone()
        } else {
            workdir.join(&input.output)
        };

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            // We pass the resolved output path through the environment
            // so run() doesn't have to re-derive it. Using an env var
            // keeps the PreparedJob struct stable.
            environment: vec![(
                OsString::from("VALENX_NWCHEM_OUTPUT"),
                output_path.into_os_string(),
            )],
            // SCF on small basis sets finishes in seconds; CCSD(T) on
            // a few-hundred-atom system runs for hours; plane-wave
            // DFT on full unit cells can run for days. 24 hours
            // covers the long tail.
            estimated_runtime: Some(Duration::from_secs(24 * 60 * 60)),
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

        // Recover the staged output path. prepare() stashed it in the
        // environment so run() doesn't have to re-derive it.
        let output_path = job
            .environment
            .iter()
            .find(|(k, _)| k == "VALENX_NWCHEM_OUTPUT")
            .map(|(_, v)| std::path::PathBuf::from(v))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "VALENX_NWCHEM_OUTPUT missing from PreparedJob.environment"
                ))
            })?;

        // Open the output file. NWChem writes its run report to
        // stdout, so we hand the file's FD to the child as stdout.
        let out_file = std::fs::File::create(&output_path).map_err(|e| {
            AdapterError::Other(anyhow::anyhow!(
                "open {} for write: {e}",
                output_path.display()
            ))
        })?;

        let program = &job.native_command[0];
        let args: Vec<&OsString> = job.native_command.iter().skip(1).collect();

        ctx.report_progress(0.0, "starting NWChem");
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
        // Forward only "real" environment vars to the child. The
        // VALENX_NWCHEM_OUTPUT shim is private to this adapter and
        // shouldn't leak into NWChem's environment.
        for (k, v) in &job.environment {
            if k != "VALENX_NWCHEM_OUTPUT" {
                cmd.env(k, v);
            }
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
        // Recover the output path from the same shim used in run().
        let output_path = job
            .environment
            .iter()
            .find(|(k, _)| k == "VALENX_NWCHEM_OUTPUT")
            .map(|(_, v)| std::path::PathBuf::from(v));

        // Provenance: hash the staged output if present, falling back
        // to case.toml when the run hasn't produced one yet.
        let case_hash_input = match &output_path {
            Some(p) if p.is_file() => p.clone(),
            _ => job.workdir.join("case.toml"),
        };
        let prov = live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "NWChem",
            "unknown",
            &case_hash_input,
            None,
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);
        let mut artefacts: Vec<Artifact> = Vec::new();

        if let Some(p) = output_path {
            if p.is_file() {
                artefacts.push(Artifact {
                    path: p,
                    kind: ArtifactKind::Log,
                    checksum: None,
                    label: "NWChem output".to_string(),
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
            ribbon_contributions: vec!["bio.nwchem.compute"],
        }
    }
}

/// Mirror of MAFFT's stderr handler — log every line and lift the
/// few obvious progress markers NWChem emits to coarse UI ticks.
fn process_stderr(line: &str, ctx: &mut RunContext<'_>, warnings: &mut Vec<String>) {
    ctx.log(LogLevel::Info, line);
    // NWChem's progress chatter on stderr is sparse (most output goes
    // to stdout). The visible markers are MPI startup banners and
    // any failure messages. We lift the obvious milestones to ticks.
    if line.contains("CITATION") || line.contains("Total times") {
        ctx.report_progress(95.0, line);
    }
    if line.contains("ERROR") || line.contains("Error:") || line.contains("FATAL") {
        warnings.push(line.trim().to_string());
    }
}

/// Bridge a child's exit status into the adapter's RunReport / Run
/// error, mirroring MAFFT's `finalize` so NWChem's failure mode
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
            format!("NWChem exited {exit_code} with no stderr output")
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
        let info = NwchemAdapter::new().info();
        assert_eq!(info.id, "nwchem");
        assert_eq!(info.physics, &[Physics::Bio]);
        assert_eq!(info.tool_license, "ECL-2.0");
        assert_eq!(info.display_name, "NWChem");
    }

    #[test]
    fn info_version_range_matches_supported_band() {
        let info = NwchemAdapter::new().info();
        assert_eq!(info.version_range.min_inclusive, Version::new(7, 2, 0));
        assert_eq!(info.version_range.max_exclusive, Version::new(8, 0, 0));
    }

    #[test]
    fn capabilities_publishes_ribbon_contribution() {
        let caps = NwchemAdapter::new().capabilities();
        assert!(caps.capabilities.is_empty());
        assert_eq!(caps.ribbon_contributions, vec!["bio.nwchem.compute"]);
    }

    #[test]
    fn license_mode_is_subprocess() {
        let info = NwchemAdapter::new().info();
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
    }

    /// Round-10 H3 RED→GREEN: `output` flowed into
    /// `workdir.join(...)` with no validation. Hostile
    /// `output = "../etc/passwd"` is now rejected.
    #[test]
    fn prepare_rejects_output_path_traversal() {
        use valenx_test_utils::tempdir;
        let d = tempdir("nwchem-output-trav");
        std::fs::write(d.join("inp.nw"), b"# test\n").unwrap();
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "nwchem.compute"

[bio.nwchem]
input     = "inp.nw"
output    = "../etc/passwd"
mpi_procs = 1
"#,
        )
        .unwrap();
        let case = Case {
            id: "trav".into(),
            path: d.clone(),
        };
        let workdir = d.join("workdir");
        let err = NwchemAdapter::new().prepare(&case, &workdir).unwrap_err();
        let msg = format!("{err}");
        assert!(
            msg.contains("[bio.nwchem].output"),
            "expected [bio.nwchem].output in error, got: {msg}"
        );
        let _ = std::fs::remove_dir_all(&d);
    }
}
