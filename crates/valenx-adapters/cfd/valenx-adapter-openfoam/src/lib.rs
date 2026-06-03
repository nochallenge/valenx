//! # valenx-adapter-openfoam
//!
//! Adapter for OpenFOAM, running under the Subprocess license mode
//! defined in [RFC 0002](../../../rfcs/0002-adapter-contract.md).
//!
//! Translates Valenx's canonical `Case` + `Mesh` into OpenFOAM
//! dictionaries (`system/`, `constant/`, `0/` inside the workdir),
//! invokes the appropriate solver binary as a child process, streams
//! progress and logs back to `valenx-core`, and collects VTK outputs
//! into canonical `Results`.
//!
//! The OpenFOAM binaries are never linked against — only executed as
//! child processes with argv and files as the only interface. This
//! is the mechanism that keeps the GPL tool isolated from Valenx's
//! Apache-licensed binary.

#![forbid(unsafe_code)]
#![allow(missing_docs)] // relaxed during pre-alpha; see valenx-fields for rationale

pub mod case_input;
pub mod dict;
pub mod log_parser;
pub mod simple_foam;

use std::ffi::OsString;
use std::path::{Path, PathBuf};
use std::time::Duration;

use semver::Version;

use valenx_core::{
    adapter_helpers::find_on_path,
    error::RunPhase,
    io_caps::{read_capped_to_bytes, MAX_VTK_FILE_BYTES},
    subprocess, Adapter, AdapterError, AdapterInfo, Capabilities, Capability, Case, LicenseMode,
    Physics, PreparedJob, ProbeReport, RunContext, RunReport, VersionRange,
};
use valenx_fields::Results;

use crate::case_input::{set_case_path, SimpleFoamInput};
use crate::log_parser::{parse_line, signal_to_sample, LogSignal};

/// Entry point the registry calls at startup.
pub fn adapter() -> Box<dyn Adapter> {
    Box::new(OpenFoamAdapter::new())
}

/// The OpenFOAM adapter.
pub struct OpenFoamAdapter;

impl OpenFoamAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for OpenFoamAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "openfoam";
/// Binaries we hunt on PATH during `probe()`. The list is ordered by
/// the most likely "smoke test" name — finding any one of them is
/// enough to confirm OpenFOAM is installed and on PATH.
const SOLVER_BINARIES: &[&str] = &[
    "simpleFoam",
    "pimpleFoam",
    "icoFoam",
    "rhoSimpleFoam",
    "foamRun",
];

fn version_range() -> VersionRange {
    // Supported OpenFOAM.com releases: v2306 .. v2506 (exclusive).
    // OpenFOAM's `vYYMM` scheme is mapped into SemVer by encoding the
    // year-month as `<year>.<month>.0`.
    VersionRange {
        min_inclusive: Version::new(23, 6, 0),
        max_exclusive: Version::new(25, 6, 0),
    }
}

impl Adapter for OpenFoamAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "OpenFOAM",
            version_range: version_range(),
            physics: &[Physics::Cfd],
            license_mode: LicenseMode::Subprocess,
            tool_license: "GPL-3.0-only",
            docs_url: "https://www.openfoam.com/documentation/",
            homepage_url: "https://www.openfoam.com/",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        // Look for any of the known solver binaries on PATH and try
        // to extract a version string from their `-help` output. The
        // adapter_helpers::detect_tool_version helper handles the
        // spawn + parse, returning None on any failure (binary
        // missing, exit non-zero, no semver in stdout/stderr).
        match find_on_path(SOLVER_BINARIES) {
            Some(binary_path) => {
                // OpenFOAM solvers print the version as part of the
                // `-help` banner: `Version: 11`. Try a few common
                // flags so the v0 detection covers ESI / Foundation
                // forks + older / newer versions without per-fork
                // branching.
                // Try several common --version-style flags. Returns
                // a String we parse into a semver::Version for the
                // ProbeReport's typed slot. Two-part versions like
                // "11" / "3.0" need padding to "11.0.0" / "3.0.0"
                // because semver::Version requires major.minor.patch.
                let raw = ["-help", "-version", "--version"].iter().find_map(|flag| {
                    valenx_core::adapter_helpers::detect_tool_version(&binary_path, flag)
                });
                let found_version = raw.as_deref().and_then(|s| {
                    // Pad shorter forms to major.minor.patch.
                    let dots = s.chars().filter(|c| *c == '.').count();
                    let normalised: String = match dots {
                        0 => format!("{s}.0.0"),
                        1 => format!("{s}.0"),
                        _ => s.to_string(),
                    };
                    semver::Version::parse(&normalised).ok()
                });
                let mut warnings = Vec::new();
                if found_version.is_none() {
                    warnings.push(
                        "version extraction from solver `-help` returned no semver — \
                         binary may be a wrapper or non-standard fork"
                            .into(),
                    );
                }
                Ok(ProbeReport {
                    ok: true,
                    found_version,
                    binary_path: Some(binary_path),
                    warnings,
                    required_env: vec![
                        ("WM_PROJECT_DIR", "set by OpenFOAM's etc/bashrc".into()),
                        ("FOAM_USER_LIBBIN", "set by OpenFOAM's etc/bashrc".into()),
                    ],
                })
            }
            None => Err(AdapterError::ToolNotInstalled {
                name: INFO_ID,
                hint: "no simpleFoam or foamRun binary on PATH; \
                       source OpenFOAM's etc/bashrc or install via first-run wizard"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        // 1. Load canonical case.toml. The parser resolves the solver
        //    string into a `SolverKind`, so prepare() doesn't need its
        //    own solver-string check.
        let case_toml = case.path.join("case.toml");
        let (_header, input) = SimpleFoamInput::from_case_dir(&case.path)?;

        // 2. Generate every dict file under `workdir/`.
        std::fs::create_dir_all(workdir)?;
        simple_foam::write_case(&input, workdir).map_err(AdapterError::Io)?;

        // 2a. If the user supplied a gmsh `.msh` file (either already
        //     in the workdir or alongside `case.toml`), convert it
        //     into `constant/polyMesh/` via `gmshToFoam`. Skipped
        //     silently when there's already a polyMesh or no .msh
        //     to be found — the solver will complain at run time
        //     about a missing mesh, which is the correct UX for
        //     "you forgot to generate a mesh".
        ensure_poly_mesh(&case.path, workdir)?;

        // 3. Resolve the specific solver binary (simpleFoam / pimpleFoam
        //    / icoFoam) — fall back to any OpenFOAM binary on PATH so
        //    the error message lands on a meaningful step rather than
        //    the generic ToolNotInstalled. If we can't even find a
        //    fallback, surface the canonical "tool not installed" so
        //    the first-run wizard can pick it up.
        let target_binary = input.solver.binary();
        let binary_path = find_on_path(&[target_binary])
            .or_else(|| find_on_path(SOLVER_BINARIES))
            .ok_or_else(|| AdapterError::ToolNotInstalled {
                name: INFO_ID,
                hint: format!(
                    "no OpenFOAM solver (looked for {target_binary} and friends) \
                     on PATH — source OpenFOAM's etc/bashrc before running, \
                     or install OpenFOAM via the first-run wizard"
                ),
            })?;
        // Loud-but-recoverable warning: the case asked for solver X
        // but we found a different binary. We still ran the dict-
        // emission step so users can fix their PATH and retry without
        // regenerating the case.
        let resolved_binary = binary_path
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or(target_binary);
        if resolved_binary != target_binary {
            return Err(set_case_path(
                AdapterError::InvalidCase {
                    case_path: std::path::PathBuf::new(),
                    reason: format!(
                        "case requested OpenFOAM solver `{target_binary}`, but \
                         it is not on PATH (found `{resolved_binary}` instead). \
                         Install the solver or change `case.solver` to match \
                         what's available."
                    ),
                },
                &case_toml,
            ));
        }

        // 4. Build a runnable invocation. OpenFOAM solvers interpret
        //    the cwd as the case directory; we pass `-case` explicitly
        //    so the process is location-agnostic.
        let native_command: Vec<OsString> = vec![
            binary_path.into_os_string(),
            OsString::from("-case"),
            workdir.as_os_str().to_os_string(),
        ];

        // 5. Rough runtime estimate so the UI can size its progress
        //    indicators. Steady solvers: 15 ms per iteration. Transient
        //    solvers: 15 ms per time step (end_time / delta_t). Both
        //    are mid-range-laptop guesses; the live progress bar takes
        //    over after the first few residuals.
        let estimated_steps = match input.time {
            crate::case_input::TimeMode::Steady => input.iterations,
            crate::case_input::TimeMode::Transient {
                end_time, delta_t, ..
            } => {
                if delta_t > 0.0 {
                    (end_time / delta_t).ceil() as u64
                } else {
                    1
                }
            }
        };
        let estimated_runtime = Some(Duration::from_millis(estimated_steps.saturating_mul(15)));

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            estimated_runtime,
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        run_prepared_job(job, ctx)
    }

    fn collect(&self, job: &PreparedJob) -> Result<Results, AdapterError> {
        // Walk the workdir for artifacts OpenFOAM left behind, then
        // parse any ASCII `.vtu` files into the canonical Field
        // catalog so the viewport / report layer can render them.
        // Older runs that haven't had `foamToVTK` applied still get
        // their raw artifacts (logs, residual histories, polyMesh
        // dirs) registered for the user to reach manually.
        // Hash the OpenFOAM case dictionaries that the adapter
        // wrote during prepare. controlDict + system + constant are
        // the canonical "this case is configured this way" inputs;
        // we hash controlDict because it's the smallest stable
        // identifier (varies per case but doesn't shift with mesh
        // refinement).
        let case_path = job.workdir.join("system").join("controlDict");
        let mesh_path = job.workdir.join("constant").join("polyMesh").join("points");
        let prov = valenx_core::adapter_helpers::live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "OpenFOAM",
            "unknown",
            &case_path,
            if mesh_path.exists() {
                Some(mesh_path.as_path())
            } else {
                None
            },
            None,
            0.0,
        );
        let mut results = Results::empty("openfoam", prov);
        if job.workdir.is_dir() {
            results.artifacts = discover_artifacts(&job.workdir);
            // Parse every .vtu artifact into the canonical Field
            // catalog. Failures on individual files are logged-and-
            // skipped rather than fatal — a broken VTK shouldn't
            // wipe out residual history + the rest of the artifacts.
            load_vtk_fields_into_results(&mut results);
        }
        Ok(results)
    }

    fn capabilities(&self) -> Capabilities {
        Capabilities {
            capabilities: vec![
                Capability::CfdSteady,
                Capability::CfdTransient,
                Capability::CfdIncompressible,
                Capability::CfdCompressible,
                Capability::CfdMultiphase,
                Capability::CfdTurbulenceRans,
                Capability::CfdTurbulenceLes,
            ],
            ribbon_contributions: vec![
                "cfd.openfoam.simple",
                "cfd.openfoam.piso",
                "cfd.openfoam.pimple",
                "cfd.openfoam.rhoSimple",
            ],
        }
    }
}

// (uuid_placeholder removed — collect() now uses
// adapter_helpers::live_provenance which generates a fresh
// UUIDv4-shaped run_id per call.)

// ---------------------------------------------------------------------------
// polyMesh materialisation via OpenFOAM's `gmshToFoam`.
// ---------------------------------------------------------------------------

/// Ensure `workdir/constant/polyMesh/` exists. Three cases:
///
/// 1. Already present → Ok(()).
/// 2. A `mesh.msh` (gmsh v4.1) file is sitting in the workdir or the
///    case directory → copy it into the workdir if needed, then
///    invoke `gmshToFoam mesh.msh` with the workdir as cwd.
/// 3. No polyMesh and no `.msh` → also Ok(()), on the theory that
///    the user may have their own mesh-generation pipeline we
///    don't know about. simpleFoam will emit a clear error at run
///    time if the mesh isn't there.
///
/// Case (2) fails loudly when `gmshToFoam` isn't on `PATH` (it
/// ships with OpenFOAM, so that's the same install). Structured
/// `AdapterError::Run` surfaces a non-zero gmshToFoam exit with
/// its stderr tail so the UI can show it.
fn ensure_poly_mesh(case_path: &Path, workdir: &Path) -> Result<(), AdapterError> {
    // Case (1): polyMesh already there.
    if workdir
        .join("constant")
        .join("polyMesh")
        .join("points")
        .is_file()
    {
        return Ok(());
    }

    // Case (2): find a .msh to convert.
    let candidates = [workdir.join("mesh.msh"), case_path.join("mesh.msh")];
    let msh_source = candidates.iter().find(|p| p.is_file());
    let Some(msh_source) = msh_source else {
        // Case (3): nothing to do.
        return Ok(());
    };

    // Stage the mesh in the workdir if it isn't already.
    let msh_in_workdir = workdir.join("mesh.msh");
    if msh_source != &msh_in_workdir {
        std::fs::copy(msh_source, &msh_in_workdir)?;
    }

    let binary_path =
        find_on_path(&["gmshToFoam"]).ok_or_else(|| AdapterError::ToolNotInstalled {
            name: INFO_ID,
            hint: "`gmshToFoam` not found on PATH — install OpenFOAM or \
                   pre-generate `constant/polyMesh/` manually (blockMesh, \
                   cfMesh, snappyHexMesh, …)"
                .into(),
        })?;

    let output = std::process::Command::new(&binary_path)
        .current_dir(workdir)
        .arg("mesh.msh")
        .output()
        .map_err(|e| AdapterError::Run {
            exit_code: -1,
            stderr: format!("gmshToFoam failed to spawn: {e}"),
            phase: RunPhase::Startup,
        })?;

    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
        return Err(AdapterError::Run {
            exit_code: output.status.code().unwrap_or(-1),
            stderr: if stderr.is_empty() {
                "gmshToFoam failed with no stderr output".to_string()
            } else {
                stderr
            },
            phase: RunPhase::Startup,
        });
    }

    // Sanity check — if gmshToFoam succeeded we should now have a
    // points file in place. If we don't, surface that as a parse
    // error rather than letting simpleFoam fail a frame later.
    if !workdir
        .join("constant")
        .join("polyMesh")
        .join("points")
        .is_file()
    {
        return Err(AdapterError::ParseOutput {
            file: workdir.join("constant").join("polyMesh").join("points"),
            reason: "gmshToFoam exited 0 but no constant/polyMesh/points \
                     was produced"
                .to_string(),
        });
    }

    Ok(())
}

// ---------------------------------------------------------------------------
// Artifact discovery
// ---------------------------------------------------------------------------

/// Recursively walk the workdir and classify every file as either
/// viz data, a log, or other. Hashing and richer metadata are left
/// for a future pass; the UI mostly needs paths, kinds, and labels.
fn discover_artifacts(workdir: &Path) -> Vec<valenx_fields::artifact::Artifact> {
    use valenx_fields::artifact::Artifact;

    let mut out: Vec<Artifact> = Vec::new();
    let mut stack: Vec<PathBuf> = vec![workdir.to_path_buf()];
    while let Some(dir) = stack.pop() {
        let entries = match std::fs::read_dir(&dir) {
            Ok(e) => e,
            Err(_) => continue,
        };
        for entry in entries.flatten() {
            let path = entry.path();
            let file_type = match entry.file_type() {
                Ok(t) => t,
                Err(_) => continue,
            };
            if file_type.is_dir() {
                stack.push(path);
                continue;
            }
            if !file_type.is_file() {
                continue;
            }
            if let Some(artifact) = classify_file(&path, workdir) {
                out.push(artifact);
            }
        }
    }
    // Deterministic ordering so UI lists don't flicker across
    // collect() calls.
    out.sort_by(|a, b| a.path.cmp(&b.path));
    out
}

fn classify_file(path: &Path, workdir: &Path) -> Option<valenx_fields::artifact::Artifact> {
    use valenx_fields::artifact::{Artifact, ArtifactKind};

    let lower = path
        .extension()
        .and_then(|s| s.to_str())
        .map(|s| s.to_ascii_lowercase());
    let rel = path.strip_prefix(workdir).unwrap_or(path).to_path_buf();
    let rel_display = rel.display().to_string();

    let (kind, label) = match lower.as_deref() {
        Some("vtk") | Some("vtu") | Some("vtp") | Some("pvd") => {
            (ArtifactKind::VizData, format!("VTK: {rel_display}"))
        }
        Some("log") => (ArtifactKind::Log, format!("log: {rel_display}")),
        Some("csv") | Some("tsv") | Some("dat") => {
            (ArtifactKind::Tabular, format!("table: {rel_display}"))
        }
        Some("png") | Some("jpg") | Some("jpeg") | Some("svg") => {
            (ArtifactKind::Image, format!("image: {rel_display}"))
        }
        _ => {
            // Nameless "log" files OpenFOAM leaves (e.g. `log.simpleFoam`)
            let name = path.file_name()?.to_string_lossy();
            if name.starts_with("log.") {
                (ArtifactKind::Log, format!("log: {rel_display}"))
            } else {
                return None;
            }
        }
    };

    Some(Artifact {
        path: path.to_path_buf(),
        kind,
        checksum: None,
        label,
    })
}

// ---------------------------------------------------------------------------
// VTK ingestion — parse every .vtu artifact into canonical Fields.
// ---------------------------------------------------------------------------

/// For every `.vtu` artifact already collected, parse + convert into
/// canonical [`valenx_fields::Field`]s and insert them into the
/// `Results::fields` catalog. Per-file failures are logged-and-skipped
/// (a single malformed VTK shouldn't lose us everything else).
///
/// Time keys are derived from the filename when possible: an OpenFOAM
/// `cavity_500.vtu` becomes `TimeKey::Iteration(500)`, distinguishing
/// time-series snapshots in the catalog. Unparseable filenames fall
/// back to `TimeKey::Steady` and overwrite each other in the catalog
/// (the user only ever looks up the latest one).
fn load_vtk_fields_into_results(results: &mut Results) {
    use valenx_fields::artifact::ArtifactKind;

    // Gather VTK artifacts up-front so we can sort by filename — that
    // gives chronological order for OpenFOAM's `<case>_<N>.vtu` naming
    // and predictable order for any other writer. Both .vtu (XML) and
    // .vtk (legacy binary) extensions are accepted; vtk_dispatch routes
    // each file to the right reader by sniffing the magic prefix.
    let mut vtk_paths: Vec<PathBuf> = results
        .artifacts
        .iter()
        .filter(|a| a.kind == ArtifactKind::VizData)
        .filter(|a| {
            a.path
                .extension()
                .and_then(|s| s.to_str())
                .map(|s| {
                    let s = s.to_ascii_lowercase();
                    s == "vtu" || s == "vtk"
                })
                .unwrap_or(false)
        })
        .map(|a| a.path.clone())
        .collect();
    vtk_paths.sort();

    for path in vtk_paths {
        // Round-20 L1: cap the per-VTK read so a corrupted (or
        // adversarial) workdir with a multi-GB `.vtu` / `.vtk` file
        // can't OOM the renderer before `vtk_dispatch` validates the
        // magic bytes. Same cap as the SU2 sister site.
        let Ok(bytes) = read_capped_to_bytes(&path, MAX_VTK_FILE_BYTES) else {
            continue;
        };
        // Treat the file stem as the Mesh id so users can tell which
        // VTU/VTK produced which fields when there are multiple time
        // steps in the catalog.
        let stem = path.file_stem().and_then(|s| s.to_str()).unwrap_or("vtk");
        let (_mesh, fields) = match valenx_fields::vtk_dispatch::load_canonical(&bytes, stem) {
            Ok(pair) => pair,
            Err(_e) => {
                // Skip silently — adapters can't propagate per-file
                // warnings into Results today, and the artifact is
                // still listed for users to inspect by hand.
                continue;
            }
        };
        let timekey = time_key_from_filename(stem);
        for mut f in fields {
            f.time = timekey;
            results.fields.insert(f);
        }
    }
}

/// Derive a `TimeKey` from an OpenFOAM-style VTK filename.
///
/// Convention: `<case>_<N>` where `<N>` is the time-step index
/// (integer for steady simpleFoam, monotonic for transient). We
/// strip everything up to and including the LAST underscore and
/// try to parse the rest as `u64`; success → `Iteration(N)`,
/// failure → `Steady`.
fn time_key_from_filename(stem: &str) -> valenx_fields::TimeKey {
    let after_last_underscore = stem.rsplit('_').next().unwrap_or("");
    if let Ok(n) = after_last_underscore.parse::<u64>() {
        valenx_fields::TimeKey::Iteration(n)
    } else {
        valenx_fields::TimeKey::Steady
    }
}

// ---------------------------------------------------------------------------
// Subprocess driver — delegates to `valenx_core::subprocess::run`
// and layers OpenFOAM-specific log parsing on top via a LineHandler.
// ---------------------------------------------------------------------------

/// Drive a prepared OpenFOAM job through to completion, streaming
/// residuals and progress back through `ctx`. All the generic
/// spawn / pipe / cancel plumbing lives in `valenx_core::subprocess`
/// so this function only owns the OpenFOAM-flavoured parsing.
///
/// Steady solvers (simpleFoam) emit `Time = 1, 2, 3, ...` — pseudo-
/// time iteration counts. Transient solvers (pimpleFoam, icoFoam)
/// emit `Time = 0.0005, 0.001, 0.0015, ...` — real seconds. The run
/// loop treats both uniformly: a step counter advances every time
/// we see a Time line (regardless of whether the value jumped from
/// 1→2 or 0→0.0005), and the progress percentage divides by the
/// total estimated steps from `estimated_runtime`. The progress
/// label echoes the real-valued time so the UI shows what the solver
/// reported, not just an opaque step number.
fn run_prepared_job(job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
    let total_steps_estimate = job
        .estimated_runtime
        .map(|d| ((d.as_millis() / 15).max(1)) as u64)
        .unwrap_or(1000);

    // Step counter — increments on each Time line. Used for the
    // `iteration` field of every ResidualSample we emit, so the
    // residual chart is monotonic and well-spaced regardless of
    // whether the underlying time is integer or sub-second.
    let mut step_count: u64 = 0;
    // Transient runs are detected by observing a fractional or
    // sub-second `Time = ...` value in the log. Steady simpleFoam
    // increments by integer pseudo-time so this never trips. Once
    // set we report `converged: None` instead of `Some(false)` —
    // transient runs march to a fixed end_time and "did residuals
    // drop below threshold" is the wrong question to ask.
    let mut is_transient = false;
    let mut residual_history: Vec<valenx_core::ResidualSample> = Vec::new();

    let report = subprocess::run(job, ctx, "starting solver", |line| {
        let signal = parse_line(line);
        let mut hint = subprocess::Hint::default();
        match signal {
            LogSignal::Time(t) => {
                step_count = step_count.saturating_add(1);
                if is_transient_time(t) {
                    is_transient = true;
                }
                let pct = if total_steps_estimate > 0 {
                    (step_count as f32 / total_steps_estimate as f32) * 100.0
                } else {
                    0.0
                };
                // Format the time so steady cases stay readable
                // ("Time = 250") and transient cases preserve enough
                // precision to see motion ("Time = 0.0005").
                let label = format_time_label(t);
                hint.progress = Some((pct, label));
            }
            LogSignal::Residual { .. } => {
                if let Some(sample) = signal_to_sample(&signal, step_count) {
                    residual_history.push(sample);
                }
            }
            LogSignal::Execution { .. } => { /* purely informational */ }
            LogSignal::Other => {
                if line.contains("--> FOAM Warning") || line.contains("--> FOAM FATAL") {
                    hint.warning = Some(line.trim().to_string());
                }
            }
        }
        hint
    })?;

    // Convergence semantic differs by time mode:
    // - Steady: did the last reported initial residual drop below the
    //   target for every tracked field? `last_residual_below` tells us.
    // - Transient: not applicable. The run reaches `end_time` whether
    //   residuals are large or small; saying `Some(false)` after a
    //   successful pimpleFoam run is misleading. The UI renders `None`
    //   as "convergence unknown" which is the honest answer.
    let converged = if is_transient {
        None
    } else {
        Some(last_residual_below(&residual_history, 1e-5))
    };
    Ok(RunReport {
        exit_code: report.exit_code,
        wall_time: report.wall_time,
        converged,
        residual_history,
        warnings: report.warnings,
        final_phase: Some(RunPhase::Shutdown),
    })
}

/// Heuristic: does this `Time = …` value look transient?
///
/// Steady simpleFoam emits whole-number pseudo-time (`1, 2, 3, …`).
/// Transient pimpleFoam / icoFoam emit either:
///
/// - sub-1.0 values (e.g. `0.0005`), or
/// - whole-second values produced by accumulating fractional `delta_t`
///   that happen to land on an integer (e.g. `delta_t = 0.0005` after
///   2000 steps lands on `Time = 1`).
///
/// The first shape is unambiguous — anything fractional is transient.
/// The second shape collides with steady cases; we can't disambiguate
/// from a single value. We accept the false negative on the boundary
/// case (a transient run that lands exactly on integer-second marks
/// for many consecutive steps) — by the next sub-step the heuristic
/// fires again and `is_transient` flips before the run completes.
fn is_transient_time(t: f64) -> bool {
    // Fractional value (e.g. 0.0005) → definitely transient.
    if t > 0.0 && (t - t.trunc()).abs() > 1e-9 {
        return true;
    }
    // Sub-1 non-zero value (e.g. 0.5) — also fractional but the trunc
    // check above already caught it. Whole-second integer values
    // remain ambiguous and stay false here.
    false
}

/// Format an OpenFOAM `Time = ...` value for the UI progress label.
/// Whole numbers stay integer-shaped so simpleFoam reads as
/// `"Time = 250"`; sub-second values use up to 6 significant digits
/// so transient runs read as `"Time = 0.000500"` (six s.f. is enough
/// to distinguish adjacent millisecond steps without spamming
/// 12-digit decimals at the user).
fn format_time_label(t: f64) -> String {
    if t == 0.0 {
        return "Time = 0".to_string();
    }
    let abs = t.abs();
    if abs >= 1.0 && (t - t.trunc()).abs() < 1e-9 {
        // Integer-valued: simpleFoam-flavoured.
        format!("Time = {}", t as u64)
    } else if abs >= 1e-4 {
        // Plain decimal with trimmed trailing zeros.
        let s = format!("{t:.6}");
        let s = s.trim_end_matches('0').trim_end_matches('.');
        format!("Time = {s}")
    } else {
        // Very small times — fall back to scientific notation.
        format!("Time = {t:.3e}")
    }
}

fn last_residual_below(history: &[valenx_core::ResidualSample], threshold: f64) -> bool {
    if history.is_empty() {
        return false;
    }
    // A steady solve has "converged" when the last-reported initial
    // residual for every tracked field is below the threshold.
    use std::collections::BTreeMap;
    let mut last_per_field: BTreeMap<&'static str, f64> = BTreeMap::new();
    for s in history {
        last_per_field.insert(s.field, s.value);
    }
    last_per_field.values().all(|v| *v < threshold)
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_test_utils::tempdir;

    #[test]
    fn info_declares_subprocess_mode() {
        let a = OpenFoamAdapter::new();
        let info = a.info();
        assert_eq!(info.id, "openfoam");
        assert_eq!(info.license_mode, LicenseMode::Subprocess);
        assert_eq!(info.physics, &[Physics::Cfd]);
    }

    #[test]
    fn version_range_contains_2406() {
        let r = version_range();
        // v2406 → 24.6.0
        assert!(r.contains(&Version::new(24, 6, 0)));
        // v2506 is the exclusive upper bound
        assert!(!r.contains(&Version::new(25, 6, 0)));
        // v2212 is below the floor
        assert!(!r.contains(&Version::new(22, 12, 0)));
    }

    #[test]
    fn capabilities_listed() {
        let a = OpenFoamAdapter::new();
        let caps = a.capabilities();
        assert!(caps.capabilities.contains(&Capability::CfdSteady));
        assert!(caps.ribbon_contributions.contains(&"cfd.openfoam.simple"));
    }

    #[test]
    fn classify_file_sorts_into_kinds() {
        use valenx_fields::artifact::ArtifactKind;
        let workdir = std::path::PathBuf::from("/tmp/workdir");
        let vtk =
            super::classify_file(&workdir.join("VTK/case_500.vtu"), &workdir).expect("classified");
        assert_eq!(vtk.kind, ArtifactKind::VizData);
        assert!(vtk.label.contains("VTK"));

        let log =
            super::classify_file(&workdir.join("log.simpleFoam"), &workdir).expect("classified");
        assert_eq!(log.kind, ArtifactKind::Log);

        let csv = super::classify_file(&workdir.join("forces.csv"), &workdir).expect("classified");
        assert_eq!(csv.kind, ArtifactKind::Tabular);

        // Unclassified noise is silently skipped — we don't want the
        // artifact list to scream at users with every temp file.
        let skip = super::classify_file(&workdir.join("system/fvSchemes"), &workdir);
        assert!(skip.is_none());
    }

    #[test]
    fn collect_loads_vtk_fields_into_catalog() {
        // Smoke test for the VTU → Field pipeline. Drop a minimal
        // ASCII .vtu in a fake workdir, run collect(), assert the
        // Field catalog has the expected scalar + vector entries.
        let tmp = tempdir("openfoam-valenx-of-collect-vtk");
        let vtk_dir = tmp.join("VTK");
        std::fs::create_dir_all(&vtk_dir).unwrap();
        let vtk_path = vtk_dir.join("cavity_500.vtu");
        std::fs::write(
            &vtk_path,
            r#"<?xml version="1.0"?>
<VTKFile type="UnstructuredGrid" version="2.0" byte_order="LittleEndian">
  <UnstructuredGrid>
    <Piece NumberOfPoints="3" NumberOfCells="1">
      <Points>
        <DataArray type="Float32" NumberOfComponents="3" format="ascii">
          0 0 0  1 0 0  0 1 0
        </DataArray>
      </Points>
      <Cells>
        <DataArray Name="connectivity" format="ascii">0 1 2</DataArray>
        <DataArray Name="offsets" format="ascii">3</DataArray>
        <DataArray Name="types" format="ascii">5</DataArray>
      </Cells>
      <PointData>
        <DataArray type="Float32" Name="p" NumberOfComponents="1" format="ascii">10 20 30</DataArray>
      </PointData>
    </Piece>
  </UnstructuredGrid>
</VTKFile>"#,
        )
        .unwrap();

        let job = PreparedJob {
            workdir: tmp.clone(),
            native_command: vec![],
            environment: vec![],
            estimated_runtime: None,
            kill_on_drop: false,
        };
        let adapter = OpenFoamAdapter::new();
        let results = adapter.collect(&job).expect("collect");

        // Artifact appears in the list (existing behaviour).
        assert!(
            results.artifacts.iter().any(|a| a.path == vtk_path),
            "VTK artifact missing from results"
        );
        // And — the new behaviour — its fields land in the catalog.
        assert!(
            !results.fields.is_empty(),
            "expected p field to be loaded into Results.fields catalog"
        );
        let p = results
            .fields
            .by_name("p")
            .next()
            .expect("p field in catalog");
        assert_eq!(p.kind, valenx_fields::FieldKind::Scalar);
        assert_eq!(p.location, valenx_fields::Location::OnNode);
        // Filename `cavity_500.vtu` → TimeKey::Iteration(500).
        assert_eq!(p.time, valenx_fields::TimeKey::Iteration(500));
        cleanup(&tmp);
    }

    #[test]
    fn time_key_from_filename_handles_canonical_shapes() {
        use valenx_fields::TimeKey;
        // OpenFOAM's foamToVTK convention: `<case>_<N>`.
        assert_eq!(
            super::time_key_from_filename("cavity_500"),
            TimeKey::Iteration(500)
        );
        assert_eq!(
            super::time_key_from_filename("cavity_0"),
            TimeKey::Iteration(0)
        );
        // Multi-underscore: only the last segment counts.
        assert_eq!(
            super::time_key_from_filename("flow_around_cylinder_120"),
            TimeKey::Iteration(120)
        );
        // No trailing number → Steady.
        assert_eq!(super::time_key_from_filename("snapshot"), TimeKey::Steady);
        assert_eq!(super::time_key_from_filename(""), TimeKey::Steady);
    }

    #[test]
    fn load_vtk_fields_into_results_picks_up_legacy_binary_artifact() {
        // Synthesise a minimal VTK legacy binary file (1 tet + 1
        // scalar field), register it as a VizData artifact, and
        // verify the collector lands the field in results.fields.
        // Catches regressions in the .vtu/.vtk extension filter and
        // the dispatcher routing.
        use valenx_fields::artifact::{Artifact, ArtifactKind};
        use valenx_fields::provenance::Sha256Hex;

        let tmp = std::env::temp_dir().join(format!(
            "valenx-of-vtk-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&tmp).unwrap();
        let vtk_path = tmp.join("cavity_0.vtk");
        let mut buf: Vec<u8> = Vec::new();
        buf.extend_from_slice(b"# vtk DataFile Version 3.0\n");
        buf.extend_from_slice(b"openfoam smoke\n");
        buf.extend_from_slice(b"BINARY\n");
        buf.extend_from_slice(b"DATASET UNSTRUCTURED_GRID\n");
        buf.extend_from_slice(b"POINTS 4 float\n");
        for v in [
            0.0_f32, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0,
        ] {
            buf.extend_from_slice(&v.to_be_bytes());
        }
        buf.push(b'\n');
        buf.extend_from_slice(b"CELLS 1 5\n");
        for v in [4u32, 0, 1, 2, 3] {
            buf.extend_from_slice(&v.to_be_bytes());
        }
        buf.push(b'\n');
        buf.extend_from_slice(b"CELL_TYPES 1\n");
        buf.extend_from_slice(&10u32.to_be_bytes());
        buf.push(b'\n');
        buf.extend_from_slice(b"POINT_DATA 4\n");
        buf.extend_from_slice(b"SCALARS p float 1\n");
        buf.extend_from_slice(b"LOOKUP_TABLE default\n");
        for v in [0.0_f32, 1.0, 2.0, 3.0] {
            buf.extend_from_slice(&v.to_be_bytes());
        }
        buf.push(b'\n');
        std::fs::write(&vtk_path, buf).unwrap();

        let prov = valenx_fields::Provenance {
            adapter: "openfoam".into(),
            adapter_version: "0".into(),
            tool: "OpenFOAM".into(),
            tool_version: "0".into(),
            case_hash: Sha256Hex::new(""),
            mesh_hash: Sha256Hex::new(""),
            input_hash: Sha256Hex::new(""),
            tools_lock_hash: Sha256Hex::new(""),
            run_id: "00000000-0000-0000-0000-000000000000".into(),
            wall_time_seconds: 0.0,
            completed_at: "1970-01-01T00:00:00Z".into(),
            ancestors: Vec::new(),
        };
        let mut results = Results::empty("cavity", prov);
        results.artifacts.push(Artifact {
            path: vtk_path.clone(),
            kind: ArtifactKind::VizData,
            checksum: None,
            label: "smoke vtk".into(),
        });

        super::load_vtk_fields_into_results(&mut results);
        // The "p" scalar field must land in the catalog with the
        // file-stem-derived TimeKey::Iteration(0).
        assert!(
            results.fields.names().any(|n| n == "p"),
            "expected p in fields; got names: {:?}",
            results.fields.names().collect::<Vec<_>>()
        );

        let _ = std::fs::remove_dir_all(&tmp);
    }

    #[test]
    fn discover_artifacts_walks_subdirs() {
        use valenx_fields::artifact::ArtifactKind;

        let tmp = std::env::temp_dir().join(format!(
            "valenx-of-collect-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(tmp.join("VTK")).unwrap();
        std::fs::write(tmp.join("log.simpleFoam"), b"Time = 1\n").unwrap();
        std::fs::write(tmp.join("VTK").join("case_500.vtu"), b"<vtk />").unwrap();
        std::fs::write(tmp.join("forces.csv"), b"a,b,c\n").unwrap();

        let artifacts = super::discover_artifacts(&tmp);
        // 3 classified files — log, vtu, csv. fvSchemes and friends
        // get skipped as uninteresting.
        assert_eq!(artifacts.len(), 3);
        let has_vtk = artifacts.iter().any(|a| a.kind == ArtifactKind::VizData);
        let has_log = artifacts.iter().any(|a| a.kind == ArtifactKind::Log);
        let has_csv = artifacts.iter().any(|a| a.kind == ArtifactKind::Tabular);
        assert!(has_vtk && has_log && has_csv);

        let _ = std::fs::remove_dir_all(tmp);
    }

    #[test]
    fn ensure_poly_mesh_noop_when_polymesh_exists() {
        // If `constant/polyMesh/points` is already in place, the
        // helper must not try to spawn gmshToFoam — so this test
        // passes even on machines where OpenFOAM isn't installed.
        let tmp = tempdir("openfoam-ensure_poly_mesh_present");
        std::fs::create_dir_all(tmp.join("constant").join("polyMesh")).unwrap();
        std::fs::write(
            tmp.join("constant").join("polyMesh").join("points"),
            "(0 0 0)\n",
        )
        .unwrap();
        let result = super::ensure_poly_mesh(std::path::Path::new("/nonexistent"), &tmp);
        assert!(result.is_ok(), "got {result:?}");
        cleanup(&tmp);
    }

    #[test]
    fn ensure_poly_mesh_noop_when_no_mesh() {
        // With no polyMesh and no .msh, the helper should return
        // Ok(()) and let simpleFoam complain at run time about the
        // missing mesh — it's not the adapter's job to second-guess
        // user workflows.
        let case_dir = tempdir("openfoam-ensure_poly_mesh_no_case");
        let workdir = tempdir("openfoam-ensure_poly_mesh_no_workdir");
        let result = super::ensure_poly_mesh(&case_dir, &workdir);
        assert!(result.is_ok(), "got {result:?}");
        // No polyMesh was created.
        assert!(!workdir
            .join("constant")
            .join("polyMesh")
            .join("points")
            .is_file());
        cleanup(&case_dir);
        cleanup(&workdir);
    }

    #[test]
    fn ensure_poly_mesh_errors_cleanly_without_gmshtofoam() {
        // With a mesh.msh staged but gmshToFoam missing, we expect
        // a structured ToolNotInstalled error — not a panic, not a
        // silent Ok.
        let case_dir = tempdir("openfoam-ensure_poly_mesh_with_msh");
        let workdir = tempdir("openfoam-ensure_poly_mesh_with_msh_wd");
        std::fs::write(
            case_dir.join("mesh.msh"),
            b"$MeshFormat\n4.1 0 8\n$EndMeshFormat\n",
        )
        .unwrap();

        let result = super::ensure_poly_mesh(&case_dir, &workdir);
        match result {
            Err(AdapterError::ToolNotInstalled { name, hint }) => {
                assert_eq!(name, "openfoam");
                assert!(
                    hint.contains("gmshToFoam"),
                    "hint should mention gmshToFoam, got: {hint}"
                );
                // Also — the staging step should still have copied
                // the .msh into the workdir so a future retry with
                // gmshToFoam installed works without re-staging.
                assert!(workdir.join("mesh.msh").is_file());
            }
            Err(AdapterError::Run { .. }) => {
                // Machine has gmshToFoam but the fake mesh.msh is
                // invalid — the helper will surface the non-zero
                // exit as a Run error. Also acceptable.
            }
            Ok(()) => {
                // Machine has gmshToFoam AND somehow accepted the
                // trivial mesh. Equally fine.
            }
            Err(other) => panic!("unexpected error: {other:?}"),
        }
        cleanup(&case_dir);
        cleanup(&workdir);
    }

    fn cleanup(d: &std::path::Path) {
        let _ = std::fs::remove_dir_all(d);
    }

    #[test]
    fn is_transient_time_detects_fractional_values() {
        // Whole-number pseudo-time (steady simpleFoam) → not transient.
        assert!(!is_transient_time(0.0));
        assert!(!is_transient_time(1.0));
        assert!(!is_transient_time(250.0));
        // Fractional time (transient pimpleFoam / icoFoam) → transient.
        assert!(is_transient_time(0.0005));
        assert!(is_transient_time(0.5));
        assert!(is_transient_time(1.25));
        assert!(is_transient_time(0.001));
        // Edge case: very small float noise around a whole second
        // shouldn't trip the heuristic — we're tolerant within 1e-9.
        assert!(!is_transient_time(1.0 + 1e-10));
        // But a clearly-fractional value at the same magnitude does.
        assert!(is_transient_time(1.0 + 1e-6));
    }

    #[test]
    fn format_time_label_handles_steady_and_transient_shapes() {
        // Steady: integer pseudo-time prints as a whole number.
        assert_eq!(format_time_label(0.0), "Time = 0");
        assert_eq!(format_time_label(1.0), "Time = 1");
        assert_eq!(format_time_label(250.0), "Time = 250");
        // Transient: sub-second values keep enough digits to be
        // distinguishable but trim trailing zeros.
        assert_eq!(format_time_label(0.5), "Time = 0.5");
        assert_eq!(format_time_label(0.05), "Time = 0.05");
        assert_eq!(format_time_label(0.0005), "Time = 0.0005");
        // Very small times (microseconds and below) read better in
        // scientific notation.
        assert_eq!(format_time_label(5e-6), "Time = 5.000e-6");
        // Mixed: fractional > 1 keeps decimals.
        assert_eq!(format_time_label(1.25), "Time = 1.25");
    }

    #[test]
    fn last_residual_below_convergence() {
        use valenx_core::ResidualSample;
        let hist = vec![
            ResidualSample {
                iteration: 1,
                field: "Ux",
                value: 0.1,
            },
            ResidualSample {
                iteration: 1,
                field: "p",
                value: 0.2,
            },
            ResidualSample {
                iteration: 500,
                field: "Ux",
                value: 1e-6,
            },
            ResidualSample {
                iteration: 500,
                field: "p",
                value: 2e-6,
            },
        ];
        assert!(super::last_residual_below(&hist, 1e-5));
        assert!(!super::last_residual_below(&hist, 1e-7));
    }

    #[test]
    fn last_residual_below_needs_all_fields() {
        use valenx_core::ResidualSample;
        // Ux is converged, p isn't → not converged.
        let hist = vec![
            ResidualSample {
                iteration: 500,
                field: "Ux",
                value: 1e-6,
            },
            ResidualSample {
                iteration: 500,
                field: "p",
                value: 1e-2,
            },
        ];
        assert!(!super::last_residual_below(&hist, 1e-5));
    }

    #[test]
    fn prepare_emits_transient_case_tree_against_fixture() {
        // Mirror of `prepare_emits_full_case_tree_against_fixture` but
        // for the new transient (pimpleFoam) fixture. Verifies the
        // controlDict picks the right solver name and the writer emits
        // a real-seconds endTime rather than an iteration count.
        let manifest = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        let workspace_root = manifest
            .ancestors()
            .nth(4)
            .expect("climb up to workspace root")
            .to_path_buf();
        let case_dir = workspace_root
            .join("tests")
            .join("fixtures")
            .join("minimal.valenx")
            .join("cases")
            .join("cfd-transient");
        if !case_dir.join("case.toml").is_file() {
            eprintln!("skipping: fixture not found at {}", case_dir.display());
            return;
        }

        let case = Case {
            id: "cfd-transient".into(),
            path: case_dir,
        };
        let workdir = std::env::temp_dir().join(format!(
            "valenx-of-prepare-transient-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));

        let adapter = OpenFoamAdapter::new();
        let _ = match adapter.prepare(&case, &workdir) {
            Ok(p) => p,
            Err(AdapterError::ToolNotInstalled { .. }) | Err(AdapterError::InvalidCase { .. }) => {
                // pimpleFoam isn't on PATH in CI, OR is on PATH but
                // some other binary won (unlikely but possible). Either
                // way the dict files should already be on disk.
                assert!(
                    workdir.join("system").join("controlDict").is_file(),
                    "controlDict missing — dict writer didn't run before \
                     the binary lookup"
                );
                let control = std::fs::read_to_string(workdir.join("system/controlDict"))
                    .expect("read controlDict");
                assert!(
                    control.contains("application     pimpleFoam;"),
                    "transient fixture should select pimpleFoam: {control}"
                );
                let _ = std::fs::remove_dir_all(&workdir);
                return;
            }
            Err(e) => panic!("prepare failed: {e:?}"),
        };

        // pimpleFoam is on PATH and won the lookup. Verify the
        // controlDict + fvSolution are transient-shaped.
        let control = std::fs::read_to_string(workdir.join("system/controlDict")).unwrap();
        assert!(control.contains("application     pimpleFoam;"));
        assert!(control.contains("writeControl    adjustableRunTime;"));
        let fv_solution = std::fs::read_to_string(workdir.join("system/fvSolution")).unwrap();
        assert!(fv_solution.contains("PIMPLE\n"));
        assert!(!fv_solution.contains("relaxationFactors\n"));

        let _ = std::fs::remove_dir_all(&workdir);
    }

    #[test]
    fn prepare_emits_full_case_tree_against_fixture() {
        // Resolve the workspace-level fixture. CARGO_MANIFEST_DIR =
        // crates/valenx-adapters/cfd/valenx-adapter-openfoam at test
        // time; the fixture lives four directories up at
        // `tests/fixtures/minimal.valenx/cases/cfd-steady`.
        let manifest = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        let workspace_root = manifest
            .ancestors()
            .nth(4)
            .expect("climb up to workspace root")
            .to_path_buf();
        let case_dir = workspace_root
            .join("tests")
            .join("fixtures")
            .join("minimal.valenx")
            .join("cases")
            .join("cfd-steady");
        if !case_dir.join("case.toml").is_file() {
            // Workspace layout changed; skip rather than fail the
            // crate-level test — the roundtrip integration test in
            // valenx-core will already have flagged the same drift.
            eprintln!("skipping: fixture not found at {}", case_dir.display());
            return;
        }

        let case = Case {
            id: "cfd-steady".into(),
            path: case_dir,
        };
        let workdir = std::env::temp_dir().join(format!(
            "valenx-of-prepare-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));

        let adapter = OpenFoamAdapter::new();
        let prepared = match adapter.prepare(&case, &workdir) {
            Ok(p) => p,
            Err(AdapterError::ToolNotInstalled { .. }) => {
                // OpenFOAM isn't on PATH in CI — the dict files should
                // still be written. Assert that and return.
                assert!(workdir.join("system").join("controlDict").is_file());
                let _ = std::fs::remove_dir_all(&workdir);
                return;
            }
            Err(e) => panic!("prepare failed: {e:?}"),
        };

        assert_eq!(prepared.workdir, workdir);
        assert!(prepared
            .native_command
            .iter()
            .any(|s| s.to_string_lossy().contains("Foam")));
        for rel in &[
            "system/controlDict",
            "system/fvSchemes",
            "system/fvSolution",
            "constant/transportProperties",
            "constant/turbulenceProperties",
            "0/U",
            "0/p",
        ] {
            assert!(
                workdir.join(rel).is_file(),
                "prepare should have emitted {rel}"
            );
        }

        let _ = std::fs::remove_dir_all(&workdir);
    }
}
