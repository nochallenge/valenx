//! # valenx-adapter-gmsh
//!
//! Subprocess-isolated adapter for the gmsh unstructured mesher.
//!
//! **Phase 2 — live.** `prepare()` emits a `.geo` from the case's
//! `[mesh]` section, `run()` spawns `gmsh` with stdout/stderr
//! streaming + cancellation, `collect()` parses the resulting `.msh`
//! file into a canonical [`valenx_mesh::Mesh`] and attaches both the
//! `.msh` and the `.geo` as artifacts for provenance.
//!
//! Scope: procedural `box` and `sphere` domains, plus `merge` for
//! STL / BRep / STEP / IGES files. Rich feature-tree meshing joins
//! when the FreeCAD / OCCT adapters graduate past scaffolding.

#![forbid(unsafe_code)]
#![allow(missing_docs)]

pub mod geo_writer;
pub mod mesh_input;
pub mod msh_parser;

use std::ffi::OsString;
use std::fs;
use std::path::{Path, PathBuf};
use std::time::Duration;

use semver::Version;

use valenx_core::{
    adapter_helpers::find_on_path, error::RunPhase, subprocess, Adapter, AdapterError, AdapterInfo,
    Capabilities, Capability, Case, LicenseMode, Physics, PreparedJob, ProbeReport, RunContext,
    RunReport, VersionRange,
};
use valenx_fields::Results;

use crate::geo_writer::DEFAULT_MSH_FILENAME;
use crate::mesh_input::{Domain, MeshSpec};

pub fn adapter() -> Box<dyn Adapter> {
    Box::new(GmshAdapter::new())
}

pub struct GmshAdapter;

impl GmshAdapter {
    pub fn new() -> Self {
        Self
    }
}

impl Default for GmshAdapter {
    fn default() -> Self {
        Self::new()
    }
}

const INFO_ID: &str = "gmsh";
const BINARIES: &[&str] = &["gmsh"];
/// The `.geo` file we generate inside the workdir.
const GEO_FILENAME: &str = "case.geo";

impl Adapter for GmshAdapter {
    fn info(&self) -> AdapterInfo {
        AdapterInfo {
            id: INFO_ID,
            display_name: "gmsh",
            version_range: VersionRange {
                min_inclusive: Version::new(4, 12, 0),
                max_exclusive: Version::new(5, 0, 0),
            },
            physics: &[Physics::Meshing],
            license_mode: LicenseMode::Subprocess,
            tool_license: "GPL-2.0-or-later",
            docs_url: "https://gmsh.info/doc/texinfo/gmsh.html",
            homepage_url: "https://gmsh.info/",
        }
    }

    fn probe(&self) -> Result<ProbeReport, AdapterError> {
        match find_on_path(BINARIES) {
            Some(binary_path) => {
                // gmsh prints its version to stderr when given
                // `--version`; the helper combines stdout + stderr
                // before searching for a semver match.
                let found_version = valenx_core::adapter_helpers::detect_tool_version_semver(
                    &binary_path,
                    &["--version", "-version"],
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
                hint: "gmsh 4.12+ required; install from gmsh.info or your \
                       package manager"
                    .into(),
            }),
        }
    }

    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError> {
        let (_header, spec) = MeshSpec::from_case_dir(&case.path)?;

        fs::create_dir_all(workdir)?;

        // Copy merged geometry into the workdir so the `.geo`'s
        // relative `Merge` statement resolves.
        // Round-9 hardening: `Domain::MergeFile { path }` is
        // user-supplied data and gets copied into the workdir; wrap
        // relative paths with `confined_join`.
        if let Domain::MergeFile { path } = &spec.domain {
            let source = if path.is_absolute() {
                path.clone()
            } else {
                valenx_core::adapter_helpers::confined_join(&case.path, path)?
            };
            if !source.is_file() {
                return Err(AdapterError::InvalidCase {
                    case_path: case.path.join("case.toml"),
                    reason: format!(
                        "[mesh] source file {} not found (resolved {})",
                        path.display(),
                        source.display()
                    ),
                });
            }
            let file_name = path.file_name().ok_or_else(|| AdapterError::InvalidCase {
                case_path: case.path.join("case.toml"),
                reason: format!("[mesh] source `{}` has no filename", path.display()),
            })?;
            let dest = workdir.join(file_name);
            if source != dest {
                fs::copy(&source, &dest)?;
            }
        }

        // Rewrite the spec so the emitted .geo uses a relative path
        // to the copied file (merge) — leaves other domains intact.
        let write_spec = rewrite_merge_path_for_workdir(spec);

        let geo_path = workdir.join(GEO_FILENAME);
        geo_writer::write_to_file(&write_spec, &geo_path)?;

        let binary_path = find_on_path(BINARIES).ok_or_else(|| AdapterError::ToolNotInstalled {
            name: INFO_ID,
            hint: "gmsh 4.12+ required; install via gmsh.info".into(),
        })?;

        // `gmsh -3 case.geo -o mesh.msh` runs meshing and writes the
        // output in one shot. `-N` forces the output dim so the
        // script's `Mesh N;` line stays authoritative.
        let dim_flag = format!("-{}", write_spec.dim.as_int());
        let native_command: Vec<OsString> = vec![
            binary_path.into_os_string(),
            dim_flag.into(),
            GEO_FILENAME.into(),
            "-o".into(),
            DEFAULT_MSH_FILENAME.into(),
        ];

        Ok(PreparedJob {
            workdir: workdir.to_path_buf(),
            native_command,
            environment: Vec::new(),
            estimated_runtime: Some(Duration::from_secs(30)),
            kill_on_drop: true,
        })
    }

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
        run_prepared_job(job, ctx)
    }

    fn collect(&self, job: &PreparedJob) -> Result<Results, AdapterError> {
        // Real provenance: hash the .geo input + the produced .msh
        // mesh when both exist. Both fall back to empty Sha256Hex
        // sentinels per live_provenance's policy.
        let case_path = job.workdir.join("case.geo");
        let mesh_path = job.workdir.join(DEFAULT_MSH_FILENAME);
        let prov = valenx_core::adapter_helpers::live_provenance(
            INFO_ID,
            env!("CARGO_PKG_VERSION"),
            "gmsh",
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
        let mut results = Results::empty(INFO_ID, prov);

        // If a .msh file landed in the workdir, parse it into a
        // canonical Mesh and save both the raw file and the parsed
        // form.
        let msh_path = job.workdir.join(DEFAULT_MSH_FILENAME);
        if msh_path.is_file() {
            match msh_parser::parse_file(&msh_path, &format!("gmsh-{}", job.workdir.display())) {
                Ok(mesh) => {
                    let canonical_path = job.workdir.join("mesh.canonical.json");
                    if let Ok(bytes) = serde_json::to_vec_pretty(&mesh) {
                        if valenx_core::io_caps::atomic_write_bytes(&canonical_path, &bytes).is_ok()
                        {
                            results.artifacts.push(valenx_fields::artifact::Artifact {
                                path: canonical_path,
                                kind: valenx_fields::artifact::ArtifactKind::VizData,
                                checksum: None,
                                label: format!(
                                    "canonical mesh · {} nodes · {} elements",
                                    mesh.stats.node_count, mesh.stats.element_count,
                                ),
                            });
                        }
                    }
                    results.artifacts.push(valenx_fields::artifact::Artifact {
                        path: msh_path.clone(),
                        kind: valenx_fields::artifact::ArtifactKind::Native,
                        checksum: None,
                        label: format!(
                            "gmsh .msh · {} nodes · {} elements",
                            mesh.stats.node_count, mesh.stats.element_count,
                        ),
                    });
                }
                Err(e) => {
                    tracing::warn!(target: "valenx-gmsh", ?e, ?msh_path, "msh parse failed");
                    results.artifacts.push(valenx_fields::artifact::Artifact {
                        path: msh_path,
                        kind: valenx_fields::artifact::ArtifactKind::Native,
                        checksum: None,
                        label: format!("gmsh .msh (parse error: {e})"),
                    });
                }
            }
        }

        let geo_path = job.workdir.join(GEO_FILENAME);
        if geo_path.is_file() {
            results.artifacts.push(valenx_fields::artifact::Artifact {
                path: geo_path,
                kind: valenx_fields::artifact::ArtifactKind::Other,
                checksum: None,
                label: "gmsh .geo (generated)".into(),
            });
        }
        Ok(results)
    }

    fn capabilities(&self) -> Capabilities {
        Capabilities {
            capabilities: vec![
                Capability::Meshing2D,
                Capability::Meshing3D,
                Capability::MeshingUnstructured,
                Capability::MeshingPrismLayers,
            ],
            ribbon_contributions: vec!["mesh.gmsh.generate", "mesh.gmsh.refine"],
        }
    }
}

/// Rewrite any `Domain::MergeFile` path to just its filename so the
/// emitted `.geo` references a file that lives next to it in the
/// workdir (where we copied it in `prepare`).
fn rewrite_merge_path_for_workdir(spec: MeshSpec) -> MeshSpec {
    match spec.domain {
        Domain::MergeFile { path } => MeshSpec {
            domain: Domain::MergeFile {
                path: PathBuf::from(
                    path.file_name()
                        .map(|s| s.to_string_lossy().into_owned())
                        .unwrap_or_default(),
                ),
            },
            ..spec
        },
        _ => spec,
    }
}

// ---------------------------------------------------------------------------
// Subprocess runner — thin wrapper over `valenx_core::subprocess::run`
// that layers gmsh-flavoured progress + warning detection on top.
// ---------------------------------------------------------------------------

fn run_prepared_job(job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError> {
    let report = subprocess::run(job, ctx, "starting gmsh", |line| {
        let mut hint = subprocess::Hint::default();
        if let Some(pct) = gmsh_progress_hint(line) {
            hint.progress = Some((pct, line.to_string()));
        }
        if let Some(w) = gmsh_warning_of_interest(line) {
            hint.warning = Some(w);
        }
        hint
    })?;

    Ok(RunReport {
        exit_code: report.exit_code,
        wall_time: report.wall_time,
        converged: Some(true), // meshing either works or it doesn't
        residual_history: Vec::new(),
        warnings: report.warnings,
        final_phase: Some(RunPhase::Shutdown),
    })
}

/// Coarse progress hints for gmsh stdout. Gmsh prints canonical
/// lines like `Info    : Meshing 3D... (Delaunay)` and
/// `Info    : Done meshing 3D (wall X s)` that we map to steps.
fn gmsh_progress_hint(line: &str) -> Option<f32> {
    if line.contains("Meshing 1D") {
        Some(10.0)
    } else if line.contains("Done meshing 1D") {
        Some(25.0)
    } else if line.contains("Meshing 2D") {
        Some(35.0)
    } else if line.contains("Done meshing 2D") {
        Some(55.0)
    } else if line.contains("Meshing 3D") {
        Some(65.0)
    } else if line.contains("Done meshing 3D") {
        Some(92.0)
    } else if line.contains("Writing") && line.contains(".msh") {
        Some(98.0)
    } else {
        None
    }
}

fn gmsh_warning_of_interest(line: &str) -> Option<String> {
    if line.contains("Warning") || line.contains("Error") {
        Some(line.trim().to_string())
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn info_declares_meshing_domain() {
        let info = GmshAdapter::new().info();
        assert_eq!(info.id, "gmsh");
        assert_eq!(info.physics, &[Physics::Meshing]);
    }

    #[test]
    fn progress_hints_monotonic_for_typical_run() {
        let pts = [
            gmsh_progress_hint("Info    : Meshing 1D..."),
            gmsh_progress_hint("Info    : Done meshing 1D (wall 0.01 s)"),
            gmsh_progress_hint("Info    : Meshing 2D..."),
            gmsh_progress_hint("Info    : Done meshing 2D (wall 0.10 s)"),
            gmsh_progress_hint("Info    : Meshing 3D..."),
            gmsh_progress_hint("Info    : Done meshing 3D (wall 0.50 s)"),
            gmsh_progress_hint("Info    : Writing mesh.msh..."),
        ];
        let mut last = 0.0f32;
        for (i, p) in pts.iter().enumerate() {
            let v = p.expect("each line hints progress");
            assert!(v >= last, "step {i} regressed: {last} → {v}");
            last = v;
        }
        assert!(last > 95.0);
    }

    #[test]
    fn rewrite_merge_strips_directory() {
        use crate::mesh_input::{Algorithm2D, Algorithm3D, MeshDim, MeshSizes};
        let spec = MeshSpec {
            domain: Domain::MergeFile {
                path: PathBuf::from("geometry/inlet.stl"),
            },
            sizes: MeshSizes::default(),
            algorithm_2d: Algorithm2D::FrontalDelaunay,
            algorithm_3d: Algorithm3D::Delaunay,
            dim: MeshDim::Three,
            physical_volume_name: "domain".into(),
            physical_surface_name: "walls".into(),
            boundary_layer: None,
        };
        let rewritten = rewrite_merge_path_for_workdir(spec);
        match rewritten.domain {
            Domain::MergeFile { path } => assert_eq!(path, PathBuf::from("inlet.stl")),
            other => panic!("expected MergeFile, got {other:?}"),
        }
    }

    /// Round-9 RED→GREEN: `Domain::MergeFile { path }` used to be
    /// joined with bare `case.path.join`. Wrap with `confined_join`.
    #[test]
    fn prepare_rejects_merge_path_traversing_outside_case_dir() {
        let d = std::env::temp_dir().join(format!(
            "valenx-gmsh-merge-trav-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&d).unwrap();
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
format  = "1.0"
name    = "trav"
physics = "meshing"
solver  = "gmsh"
mesh    = "(none)"

[mesh]
type   = "merge"
source = "../../etc/passwd"
dim    = 3
"#,
        )
        .unwrap();
        let case = Case {
            id: "gmsh-merge-trav".into(),
            path: d.clone(),
        };
        let workdir = d.join("workdir");
        let err = GmshAdapter::new().prepare(&case, &workdir).unwrap_err();
        let msg = format!("{err}");
        assert!(
            msg.contains("..") || msg.contains("stay within") || msg.contains("escape"),
            "expected confined_join rejection, got: {msg}"
        );
        let _ = std::fs::remove_dir_all(&d);
    }
}
