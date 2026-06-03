//! # valenx-adapter-calculix
//!
//! Adapter for CalculiX (CCX solver + CGX pre/post).
//!
//! **Phase 3 — live for linear-static FEA.** `prepare()` reads the
//! `[structural]` section of the case, loads the canonical mesh the
//! companion gmsh adapter wrote to disk, and emits a self-contained
//! `job.inp` Abaqus-flavoured deck. `run()` spawns `ccx` via the
//! shared [`valenx_core::subprocess`] runner. `collect()` walks the
//! workdir for `.frd`, `.dat`, and `.log` artifacts.
//!
//! Scope today: linear static, modal, steady-state thermal with
//! isotropic elastic materials, Dirichlet BCs, concentrated loads.
//! Nonlinear, contact, and transient lands in follow-ups.

#![forbid(unsafe_code)]
#![allow(missing_docs)]

pub mod case_input;
pub mod inp_writer;

use std::ffi::OsString;
use std::fs;
use std::path::{Path, PathBuf};
use std::time::Duration;

use semver::Version;

use valenx_core::{
    adapter_helpers::{confined_join, find_on_path, first_workdir_match},
    error::RunPhase,
    subprocess, Adapter, AdapterError, AdapterInfo, Capabilities, Capability, Case, LicenseMode,
    Physics, PreparedJob, ProbeReport, RunContext, RunReport, VersionRange,
};
use valenx_fields::{
    artifact::{Artifact, ArtifactKind},
    Results,
};
use valenx_mesh::Mesh;

use crate::case_input::LinearStaticInput;

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(CalculixAdapter::new())
}

pub struct CalculixAdapter;

impl CalculixAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for CalculixAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "calculix";
const BINARIES: &[&str] = &["ccx", "ccx_2.22", "ccx_2.21", "ccx_2.20", "ccx_2.19"];
/// The root input deck we emit. CalculiX takes the job name without
/// the `.inp` suffix on the command line.
const INP_JOBNAME: &str = "job";

impl Adapter for CalculixAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "CalculiX",
            version_range: VersionRange {
                min_inclusive: Version::new(2, 19, 0),
                max_exclusive: Version::new(3, 0, 0),
            },
            physics: &[Physics::Fea],
            license_mode: LicenseMode::Subprocess,
            tool_license: "GPL-2.0-or-later",
            docs_url: "http://www.dhondt.de/",
            homepage_url: "http://www.calculix.de/",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(BINARIES) {
            Some(binary_path) => {
                // CalculiX `ccx -v` prints "CalculiX Version 2.21" to
                // stderr. The combined-stream capture inside the
                // helper picks it up without per-adapter branching.
                let found_version = valenx_core::adapter_helpers::detect_tool_version_semver(
                    &binary_path,
                    &["-v", "--version"],
                );
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
                hint: "CalculiX (ccx) 2.19+ required; install from calculix.de or \
                       your distribution"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let (header, input) = LinearStaticInput::from_case_dir(&case.path)?;

        fs::create_dir_all(workdir)?;

        // 1. Resolve + load the canonical mesh the [structural]
        //    section points at. Look in both the case dir and the
        //    workdir so the output of a preceding gmsh run is
        //    visible without the user copying anything.
        let mesh = load_mesh(&input.mesh_source, &case.path, workdir)?;
        if mesh.nodes.is_empty() {
            return Err(AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!(
                    "canonical mesh at {} has zero nodes — was the mesher run?",
                    input.mesh_source.display()
                ),
            });
        }

        // 2. Emit the .inp file.
        let inp_path = workdir.join(format!("{INP_JOBNAME}.inp"));
        inp_writer::write_to_file(&mesh, &input, &header.name, &inp_path)?;

        // 3. Resolve ccx and build the command line. CalculiX takes
        //    the job name WITHOUT the `.inp` suffix.
        let binary_path =
            find_first_on_path(BINARIES).ok_or_else(|| AdapterError::ToolNotInstalled {
                name: INFO_ID,
                hint: "ccx not found on PATH; install CalculiX from calculix.de".into(),
            })?;

        let native_command: Vec<OsString> = vec![
            binary_path.into_os_string(),
            OsString::from("-i"),
            OsString::from(INP_JOBNAME),
        ];

        // 4. Estimate runtime — rough but helps the UI show progress.
        //    Linear static on a 100k-node mesh typically takes ~30 s
        //    on a mid-range laptop; scale roughly linearly.
        let n_nodes = mesh.nodes.len() as u64;
        let estimated_runtime = Some(Duration::from_millis((n_nodes / 3_000).max(5) * 1_000));

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            estimated_runtime,
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        let report = subprocess::run(job, ctx, "starting ccx", |line| {
            let mut hint = subprocess::Hint::default();
            if let Some(pct) = ccx_progress_hint(line) {
                hint.progress = Some((pct, line.to_string()));
            }
            // `*ERROR` / `*WARNING` banners go straight to the warnings
            // vector so the UI can render them separately from the
            // raw log.
            if line.contains("*ERROR") || line.contains("*WARNING") {
                hint.warning = Some(line.trim().to_string());
            }
            hint
        })?;
        Ok(RunReport {
            exit_code: report.exit_code,
            wall_time: report.wall_time,
            // ccx exits non-zero on non-convergence; a successful
            // shared::run return already implies the solve finished.
            converged: Some(true),
            residual_history: Vec::new(),
            warnings: report.warnings,
            final_phase: Some(RunPhase::Shutdown),
        })
    }

    fn collect(&self, job: &PreparedJob) -> Result<Results, AdapterError> {
        // Real provenance: hash the .inp input + the .frd output
        // (results) when both exist. Adapter writes the .inp on
        // prepare and CalculiX itself produces the .frd on run.
        let case_path = first_workdir_match(&job.workdir, &["inp"])
            .unwrap_or_else(|| job.workdir.join("(no-inp-found)"));
        let mesh_path = first_workdir_match(&job.workdir, &["msh", "vtu"]);
        let prov = valenx_core::adapter_helpers::live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "CalculiX",
            "unknown",
            &case_path,
            mesh_path.as_deref(),
            None,
            0.0,
        );
        let mut results = Results::empty(INFO_ID, prov);

        let extensions_to_keep = [
            ("frd", ArtifactKind::VizData, "CalculiX .frd result"),
            ("dat", ArtifactKind::Tabular, "CalculiX .dat print"),
            ("sta", ArtifactKind::Log, "CalculiX .sta step log"),
            ("cvg", ArtifactKind::Log, "CalculiX .cvg convergence"),
            ("12d", ArtifactKind::Native, "CalculiX .12d restart"),
            ("inp", ArtifactKind::Other, "CalculiX .inp input"),
        ];

        if let Ok(entries) = fs::read_dir(&job.workdir) {
            for entry in entries.flatten() {
                let path = entry.path();
                let ext = path
                    .extension()
                    .and_then(|s| s.to_str())
                    .map(|s| s.to_ascii_lowercase());
                let Some(ext) = ext else { continue };
                if let Some(&(_, kind, label)) =
                    extensions_to_keep.iter().find(|(e, ..)| *e == ext.as_str())
                {
                    results.artifacts.push(Artifact {
                        path: path.clone(),
                        kind,
                        checksum: None,
                        label: label.to_string(),
                    });
                }
            }
        }
        results.artifacts.sort_by(|a, b| a.path.cmp(&b.path));
        // Parse every ASCII .frd artifact into the canonical Field
        // catalog. Per-file failures are skipped silently (binary
        // .frd, malformed mid-write file, etc.) — the artifact is
        // still listed for users to inspect by hand.
        load_frd_fields_into_results(&mut results);
        // Auto-derive the von-Mises stress field from any 6-component
        // stress tensor that landed. Users get the "what's the worst
        // stress" answer in the Results pane without bouncing into
        // ParaView. The original tensor stays in the catalog;
        // <name>_vms is appended.
        derive_vms_fields(&mut results);
        Ok(results)
    }

    fn capabilities(&self) -> Capabilities {
        Capabilities {
            capabilities: vec![
                Capability::FeaLinearStatic,
                Capability::FeaNonlinearStatic,
                Capability::FeaModal,
                Capability::FeaHarmonic,
                Capability::FeaTransient,
                Capability::FeaThermal,
                Capability::FeaContact,
            ],
            ribbon_contributions: vec![
                "fea.calculix.static",
                "fea.calculix.modal",
                "fea.calculix.thermal",
            ],
        }
    }
}

/// Wrap `find_on_path` so our test helpers can stub it in the future
/// without touching every call site.
fn find_first_on_path(names: &[&str]) -> Option<PathBuf> {
    find_on_path(names)
}

/// Map a CalculiX stdout line to a coarse progress percentage. Based
/// on the messages ccx 2.19+ emits in the default verbosity.
/// For every ASCII `.frd` artifact already collected, parse + convert
/// into canonical [`valenx_fields::Field`]s and insert into the
/// `Results::fields` catalog. Per-file failures are skipped silently
/// (binary .frd / malformed mid-write file / unknown block type) —
/// the artifact stays listed for the user to inspect by hand.
fn load_frd_fields_into_results(results: &mut valenx_fields::Results) {
    use valenx_fields::artifact::ArtifactKind;

    let mut frd_paths: Vec<PathBuf> = results
        .artifacts
        .iter()
        .filter(|a| a.kind == ArtifactKind::VizData)
        .filter(|a| {
            a.path
                .extension()
                .and_then(|s| s.to_str())
                .map(|s| s.eq_ignore_ascii_case("frd"))
                .unwrap_or(false)
        })
        .map(|a| a.path.clone())
        .collect();
    frd_paths.sort();

    for path in frd_paths {
        // Round-18 M1: cap the .frd read at `MAX_FRD_FILE_BYTES` so
        // a hostile multi-GB artifact in the workdir doesn't OOM the
        // post-processor. Silently skip oversize files (mirrors how
        // the loop already silently skips unreadable / binary / mid-
        // write `.frd` files — the user still sees the path listed
        // in `results.artifacts`).
        let Ok(text) = valenx_core::io_caps::read_capped_to_string(
            &path,
            valenx_core::io_caps::MAX_FRD_FILE_BYTES,
        ) else {
            continue;
        };
        let parsed = match valenx_fields::frd::parse_ascii(&text) {
            Ok(p) => p,
            Err(_) => continue,
        };
        let fields = valenx_fields::frd::to_canonical_fields(&parsed);
        for f in fields {
            results.fields.insert(f);
        }
    }
}

/// Walk every loaded field; for each 6-component stress tensor,
/// compute the von-Mises scalar via [`valenx_fields::stress::von_mises_from_components`]
/// and append `<name>_vms` to the catalog. Idempotent —
/// `_vms`-suffixed fields are skipped so calling twice doesn't
/// produce `S_vms_vms`.
fn derive_vms_fields(results: &mut valenx_fields::Results) {
    use valenx_fields::FieldKind;
    let names: Vec<String> = results.fields.names().map(|s| s.to_string()).collect();
    let mut to_insert: Vec<valenx_fields::Field> = Vec::new();
    for name in &names {
        if name.ends_with("_vms") {
            continue;
        }
        // Each name may have multiple time-step entries; derive
        // VMS for every one.
        for f in results.fields.by_name(name) {
            let is_voigt = matches!(f.kind, FieldKind::Vector { dim: 6 })
                || matches!(f.kind, FieldKind::Tensor { rows: 3, cols: 3 });
            if !is_voigt {
                continue;
            }
            if let Some(vms) = valenx_fields::stress::von_mises_from_components(f) {
                to_insert.push(vms);
            }
        }
    }
    for f in to_insert {
        results.fields.insert(f);
    }
}

fn ccx_progress_hint(line: &str) -> Option<f32> {
    if line.contains("*INFO reading input") {
        Some(5.0)
    } else if line.contains("filling the matrix") {
        Some(20.0)
    } else if line.contains("factoring the system") {
        Some(40.0)
    } else if line.contains("solving the system") {
        Some(60.0)
    } else if line.contains("creating frd") || line.contains("writing results") {
        Some(90.0)
    } else if line.contains("Job finished") {
        Some(99.0)
    } else {
        None
    }
}

/// Load the canonical mesh referenced by the case.
///
/// `source` is sandboxed via [`confined_join`] before being resolved
/// against either the workdir or the case directory — a hostile case
/// bundle cannot point `mesh_source` at `/etc/passwd` or traverse out
/// of the case sandbox.
fn load_mesh(source: &Path, case_dir: &Path, workdir: &Path) -> Result<Mesh, AdapterError> {
    // Sandbox the case-relative path. The workdir-relative variant
    // shares the same constraint (workdir is also confined).
    let case_candidate = confined_join(case_dir, source)?;
    let workdir_candidate = confined_join(workdir, source)?;
    let candidates = [workdir_candidate, case_candidate];
    let Some(path) = candidates.iter().find(|p| p.is_file()) else {
        return Err(AdapterError::InvalidCase {
            case_path: case_dir.join("case.toml"),
            reason: format!(
                "[structural] mesh_source {} not found in workdir or case dir — \
                 run the mesher first, or copy the mesh alongside case.toml",
                source.display()
            ),
        });
    };
    let ext = path
        .extension()
        .and_then(|s| s.to_str())
        .map(|s| s.to_ascii_lowercase())
        .unwrap_or_default();

    match ext.as_str() {
        "json" => {
            // Round-18 M3 (R12 M5 sister): cap the JSON mesh read at
            // the shared `MAX_MESH_JSON_BYTES` so a hostile
            // mesh.canonical.json doesn't OOM the loader before serde
            // sees it.
            let text = valenx_core::io_caps::read_capped_to_string(
                path,
                valenx_core::io_caps::MAX_MESH_JSON_BYTES,
            )?;
            serde_json::from_str::<Mesh>(&text).map_err(|e| AdapterError::ParseOutput {
                file: path.clone(),
                reason: format!("not a canonical Mesh JSON: {e}"),
            })
        }
        other => Err(AdapterError::InvalidCase {
            case_path: case_dir.join("case.toml"),
            reason: format!(
                "unsupported mesh_source extension `.{other}` — today the \
                 calculix adapter only reads the canonical JSON \
                 (`mesh.canonical.json`) produced by gmsh's collect()"
            ),
        }),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn info_declares_fea_domain() {
        let info = CalculixAdapter::new().info();
        assert_eq!(info.id, "calculix");
        assert_eq!(info.physics, &[Physics::Fea]);
    }

    #[test]
    fn collect_loads_frd_fields_into_catalog() {
        // Smoke test for the .frd → Field pipeline. Drop a minimal
        // ASCII .frd in a fake workdir, run collect(), assert the
        // DISP field lands in the Field catalog as a vector
        // OnNode field at TimeKey::Iteration(1).
        let tmp = std::env::temp_dir().join(format!(
            "valenx-ccx-collect-frd-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&tmp).unwrap();
        let frd_path = tmp.join("job.frd");
        std::fs::write(
            &frd_path,
            "    1C  UTILITY_PROGRAMS\n    \
                2C\n \
                -1   1  0 0 0\n \
                -1   2  1 0 0\n \
                -3\n    \
                1PSTEP    1    1    1\n  \
                100CL  102  0  4  1  1\n \
                -4  DISP    4    1\n \
                -5  D1    1    2    1    0\n \
                -5  D2    1    2    2    0\n \
                -5  D3    1    2    3    0\n \
                -1   1  0 0 0\n \
                -1   2  1.5e-4 0 0\n \
                -3\n",
        )
        .unwrap();

        let job = PreparedJob {
            workdir: tmp.clone(),
            native_command: vec![],
            environment: vec![],
            estimated_runtime: None,
            kill_on_drop: false,
        };
        let adapter = CalculixAdapter::new();
        let results = adapter.collect(&job).expect("collect");

        // .frd registered as artifact (existing behaviour).
        assert!(results.artifacts.iter().any(|a| a.path == frd_path));
        // DISP field landed in catalog (new behaviour).
        assert!(
            !results.fields.is_empty(),
            "expected DISP field in Results.fields"
        );
        let disp = results.fields.by_name("DISP").next().expect("DISP field");
        assert_eq!(disp.kind, valenx_fields::FieldKind::Vector { dim: 3 });
        assert_eq!(disp.location, valenx_fields::Location::OnNode);
        assert_eq!(disp.time, valenx_fields::TimeKey::Iteration(1));
        let _ = std::fs::remove_dir_all(&tmp);
    }

    #[test]
    fn ccx_progress_hints_are_monotonic() {
        let pts = [
            ccx_progress_hint("*INFO reading input"),
            ccx_progress_hint("filling the matrix"),
            ccx_progress_hint("factoring the system"),
            ccx_progress_hint("solving the system"),
            ccx_progress_hint("creating frd"),
            ccx_progress_hint("Job finished"),
        ];
        let mut last = 0.0f32;
        for (i, p) in pts.iter().enumerate() {
            let v = p.expect("known ccx token");
            assert!(v >= last, "step {i}: {last} -> {v}");
            last = v;
        }
    }

    #[test]
    fn load_mesh_missing_file_is_invalid_case() {
        let case_dir = std::env::temp_dir().join(format!(
            "ccx-load-test-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&case_dir).unwrap();
        let workdir = std::env::temp_dir().join(format!(
            "ccx-load-wd-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&workdir).unwrap();
        let err = load_mesh(
            std::path::Path::new("mesh.canonical.json"),
            &case_dir,
            &workdir,
        )
        .unwrap_err();
        assert!(matches!(err, AdapterError::InvalidCase { .. }));
        let _ = std::fs::remove_dir_all(&case_dir);
        let _ = std::fs::remove_dir_all(&workdir);
    }

    /// Round-18 M3 RED→GREEN: `load_mesh` reading a `.json` mesh
    /// source must refuse a file above `MAX_MESH_JSON_BYTES` (sister
    /// to R12 M5).
    #[test]
    fn load_mesh_rejects_oversize_json() {
        use std::io::{Seek, SeekFrom, Write};
        let case_dir = std::env::temp_dir().join(format!(
            "ccx-mesh-oversize-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&case_dir).unwrap();
        let workdir = std::env::temp_dir().join(format!(
            "ccx-mesh-oversize-wd-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&workdir).unwrap();
        // Sparse-file trick: stat reports past-the-cap but disk
        // footprint stays minimal.
        let mesh_path = workdir.join("mesh.canonical.json");
        {
            let mut f = std::fs::OpenOptions::new()
                .write(true)
                .create_new(true)
                .open(&mesh_path)
                .unwrap();
            f.seek(SeekFrom::Start(
                valenx_core::io_caps::MAX_MESH_JSON_BYTES as u64 + 1,
            ))
            .unwrap();
            f.write_all(b"x").unwrap();
        }
        let err = load_mesh(
            std::path::Path::new("mesh.canonical.json"),
            &case_dir,
            &workdir,
        )
        .expect_err("must reject oversize mesh json");
        match err {
            AdapterError::Io(_) => {} // expected: cap rejection surfaces as Io
            other => panic!("expected Io (size cap), got: {other:?}"),
        }
        let _ = std::fs::remove_dir_all(&case_dir);
        let _ = std::fs::remove_dir_all(&workdir);
    }

    #[test]
    fn derive_vms_fields_appends_vms_for_voigt_stress() {
        // Synthesize Results with a 6-component stress tensor at
        // a single node (uniaxial 100 MPa Sxx). derive_vms_fields
        // should append "S_vms" with the textbook value 100 MPa.
        use valenx_fields::{
            provenance::Sha256Hex, FieldKind, Location, RegionRef, Results, TimeKey,
        };
        let prov = valenx_fields::Provenance {
            adapter: "calculix".into(),
            adapter_version: "0".into(),
            tool: "CalculiX".into(),
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
        let mut results = Results::empty("test", prov);
        let stress = valenx_fields::Field {
            name: "S".into(),
            kind: FieldKind::Vector { dim: 6 },
            location: Location::OnNode,
            region: RegionRef("default".into()),
            units: valenx_fields::units::Units::new([-1, 1, -2, 0, 0, 0, 0], 1.0, Some("Pa")),
            time: TimeKey::Steady,
            data: vec![100.0e6, 0.0, 0.0, 0.0, 0.0, 0.0],
            range: None,
        };
        results.fields.insert(stress);
        super::derive_vms_fields(&mut results);
        let vms_present = results.fields.names().any(|n| n == "S_vms");
        assert!(vms_present, "S_vms not appended");
        let vms = results
            .fields
            .by_name("S_vms")
            .next()
            .expect("at least one S_vms entry");
        assert_eq!(vms.data.len(), 1);
        assert!(
            (vms.data[0] - 100.0e6).abs() < 1e-3,
            "uniaxial VMS should equal Sxx; got {}",
            vms.data[0]
        );
    }

    #[test]
    fn derive_vms_fields_is_idempotent() {
        // Calling twice must not produce S_vms_vms.
        use valenx_fields::{
            provenance::Sha256Hex, FieldKind, Location, RegionRef, Results, TimeKey,
        };
        let prov = valenx_fields::Provenance {
            adapter: "calculix".into(),
            adapter_version: "0".into(),
            tool: "CalculiX".into(),
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
        let mut results = Results::empty("test", prov);
        results.fields.insert(valenx_fields::Field {
            name: "S".into(),
            kind: FieldKind::Vector { dim: 6 },
            location: Location::OnNode,
            region: RegionRef("default".into()),
            units: valenx_fields::units::DIMENSIONLESS,
            time: TimeKey::Steady,
            data: vec![1.0; 6],
            range: None,
        });
        super::derive_vms_fields(&mut results);
        super::derive_vms_fields(&mut results);
        let names: Vec<&str> = results.fields.names().collect();
        assert!(names.contains(&"S_vms"));
        assert!(!names.contains(&"S_vms_vms"));
    }

    #[test]
    fn derive_vms_fields_skips_non_voigt_inputs() {
        // Scalar / 3-vector fields shouldn't trigger VMS derivation.
        use valenx_fields::{
            provenance::Sha256Hex, FieldKind, Location, RegionRef, Results, TimeKey,
        };
        let prov = valenx_fields::Provenance {
            adapter: "calculix".into(),
            adapter_version: "0".into(),
            tool: "CalculiX".into(),
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
        let mut results = Results::empty("test", prov);
        results.fields.insert(valenx_fields::Field {
            name: "U".into(),
            kind: FieldKind::Vector { dim: 3 },
            location: Location::OnNode,
            region: RegionRef("default".into()),
            units: valenx_fields::units::DIMENSIONLESS,
            time: TimeKey::Steady,
            data: vec![1.0; 9], // 3 nodes x 3 components
            range: None,
        });
        super::derive_vms_fields(&mut results);
        let names: Vec<&str> = results.fields.names().collect();
        assert!(!names.contains(&"U_vms"));
    }
}
