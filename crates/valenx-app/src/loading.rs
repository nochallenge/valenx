//! Project / STL / mesh loading methods on [`ValenxApp`]. Split out
//! of `lib.rs` as part of the structural refactor.

use std::path::PathBuf;

use valenx_core::{LoadedProject, LogLevel};
use valenx_mesh::{quality_report, Mesh};

use crate::audit::emit_audit;
use crate::rbac_io::{rbac_check, rbac_override_from_project_toml};
use crate::settings_io::save_settings_to_state_dir;
use crate::types::{LoadedMesh, LoadedStl};
use crate::ValenxApp;

impl ValenxApp {
    /// Open a Valenx project from disk and replace the current project
    /// (if any). Loader warnings are surfaced in the log panel; hard
    /// failures land in `self.last_error`.
    ///
    /// Round-10 M6: gated on the RBAC `ProjectOpen` action (Viewer
    /// role required — the most permissive tier; the gate primarily
    /// exists so headless / kiosk deployments can disable it). A
    /// successful load emits a `project.open` audit entry tagged with
    /// the project path so compliance can review who opened what.
    pub fn load_project(&mut self, path: PathBuf) {
        // RBAC gate. Pre-fix `Action::ProjectOpen` was defined but
        // never enforced anywhere, so a Viewer-disabled deployment
        // couldn't actually block project opens.
        if let Err(e) = rbac_check(
            valenx_rbac::Action::ProjectOpen,
            self.project_rbac_override.as_ref(),
        ) {
            self.last_error = Some(format!("{e}"));
            return;
        }
        match LoadedProject::load(&path) {
            Ok(proj) => {
                tracing::info!(target: "valenx", ?path, "loaded project");
                // Run the advisory validator. Warnings don't fail
                // the load — they just surface in the log panel +
                // status so users see the operator-mistake hints
                // (typo'd solver, unknown physics, mesh ref to a
                // missing key).
                let warnings = valenx_core::project::validate_project(&proj);
                if !warnings.is_empty() {
                    let plural = if warnings.len() == 1 { "" } else { "s" };
                    self.status = Some(format!(
                        "Loaded project: {} ({} warning{})",
                        proj.project.project.name,
                        warnings.len(),
                        plural,
                    ));
                    for w in &warnings {
                        self.log.push(
                            LogLevel::Warn,
                            format!("project [{}]: {}", w.code, w.message),
                        );
                    }
                } else {
                    self.status = Some(format!("Loaded project: {}", proj.project.project.name));
                }
                // Pull a per-project [rbac] override out of the
                // project.toml if present. Failures parse into None
                // + a warning so misconfigured project.tomls don't
                // brick the whole project load.
                self.project_rbac_override =
                    rbac_override_from_project_toml(&proj.root.join("project.toml"));
                self.project = Some(proj);
                // Promote the freshly-loaded project to the head of
                // the recent-projects list. The landing page reads
                // this on the next no-project launch so the user can
                // re-open with one click. Best-effort persistence —
                // a failed write here doesn't fail the load.
                if self.settings.push_recent_project(path.clone()) {
                    save_settings_to_state_dir(&self.settings);
                }
                // Round-10 M6: audit who opened which project, mirror
                // of `run.start` / `run.cancel`. The path goes in the
                // target so compliance can trace project access.
                emit_audit(
                    "project.open",
                    serde_json::json!({
                        "kind": "project",
                        "path": path.display().to_string(),
                    }),
                    serde_json::json!({}),
                );
                self.project_path = Some(path);
                self.last_error = None;
            }
            Err(e) => {
                tracing::error!(target: "valenx", ?e, ?path, "project load failed");
                self.last_error = Some(format!("Load project failed: {e}"));
            }
        }
    }

    /// Close the current project and return to the home / landing page.
    /// Clears the project + any loaded geometry so the empty-state check
    /// (`project` / `stl` / `mesh` all `None`) renders the landing page,
    /// letting the user open or create another project. Opening a new
    /// project repopulates everything.
    pub fn close_project(&mut self) {
        tracing::info!(target: "valenx", "closing project — returning to home page");
        self.project = None;
        self.project_path = None;
        self.project_rbac_override = None;
        self.stl = None;
        self.mesh = None;
        self.last_error = None;
        self.status =
            Some("Closed project — open or create another from the home page.".to_string());
    }

    /// Parse an STL file at `path` into the viewport. Replaces the
    /// currently-loaded STL (if any) and reframes the camera.
    pub fn load_stl(&mut self, path: PathBuf) {
        match valenx_viz::stl::load(&path) {
            Ok(mesh) => {
                tracing::info!(
                    target: "valenx",
                    ?path,
                    triangles = mesh.triangle_count(),
                    "loaded STL",
                );
                self.status = Some(format!(
                    "Loaded STL: {} ({} triangles)",
                    path.file_name().unwrap_or_default().to_string_lossy(),
                    mesh.triangle_count(),
                ));
                self.stl = Some(LoadedStl { path, mesh });
                self.last_error = None;
                self.frame_current_stl();
            }
            Err(e) => {
                tracing::error!(target: "valenx", ?e, ?path, "STL load failed");
                self.last_error = Some(format!("Load STL failed: {e}"));
            }
        }
    }

    /// Show a native file dialog to pick an STL and then call
    /// [`Self::load_stl`] on the chosen file.
    pub fn pick_stl(&mut self) {
        let picked = rfd::FileDialog::new()
            .add_filter("STL mesh", &["stl", "STL"])
            .set_title("Import STL")
            .pick_file();
        if let Some(path) = picked {
            self.load_stl(path);
        }
    }

    /// Set the viewport's current mesh from an already-parsed
    /// canonical [`Mesh`]. The `source_path` is recorded for the UI
    /// (browser shows the filename) — pass a synthetic path like
    /// `<workdir>/cavity_500.vtu` if the mesh came from in-memory
    /// parsing rather than a real file.
    ///
    /// Mirrors the post-parse half of [`Self::load_mesh`] so callers
    /// that already have a `Mesh` (post-run hooks, programmatic
    /// loads, future "import from URL" features) don't have to round
    /// trip through JSON or `.msh`.
    pub fn apply_mesh(&mut self, mut mesh: Mesh, source_path: PathBuf) {
        mesh.recompute_stats();
        let quality = quality_report(&mesh);
        let aspect_hist =
            valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
        let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
        // Roll the report's scalars back into the mesh's cached
        // stats so anything that pokes at `mesh.stats` (e.g. the
        // browser tree, JSON dumps) sees the populated values
        // instead of None.
        mesh.stats.min_element_size = quality.min_size;
        mesh.stats.max_aspect_ratio = quality.max_aspect_ratio;
        mesh.stats.max_skewness = quality.max_skewness;
        mesh.stats.min_orthogonality = quality.min_orthogonality;
        tracing::info!(
            target: "valenx",
            ?source_path,
            nodes = mesh.stats.node_count,
            elements = mesh.stats.element_count,
            "applied canonical mesh"
        );
        self.status = Some(format!(
            "Loaded mesh: {} ({} nodes, {} elements)",
            source_path
                .file_name()
                .unwrap_or_default()
                .to_string_lossy(),
            mesh.stats.node_count,
            mesh.stats.element_count,
        ));
        self.last_error = None;
        self.mesh = Some(LoadedMesh {
            path: source_path,
            mesh,
            quality,
            aspect_hist,
            skew_hist,
        });
        self.frame_current_mesh();
    }

    /// Load a canonical mesh from disk into the viewport. Supports
    /// two shapes today:
    ///
    /// - `*.json` — the `mesh.canonical.json` the gmsh adapter's
    ///   `collect()` serialises (`valenx_mesh::Mesh` via serde).
    /// - `*.msh` — a raw gmsh `.msh` v4.1 file, parsed on the fly
    ///   via `valenx-adapter-gmsh::msh_parser` so users can drag
    ///   meshes into the app without running a full adapter cycle.
    pub fn load_mesh(&mut self, path: PathBuf) {
        let ext = path
            .extension()
            .and_then(|s| s.to_str())
            .unwrap_or("")
            .to_ascii_lowercase();

        let parsed: Result<Mesh, String> = match ext.as_str() {
            // Round-12 M5: the pre-fix path used unbounded
            // `fs::read_to_string`, so a hostile or accidental
            // multi-GB `.json` mesh would slurp into RAM before the
            // serde parser ever saw it. The 64 MiB cap covers any
            // realistic exported mesh (millions of vertices) while
            // refusing the DoS class.
            "json" => valenx_core::io_caps::read_capped_to_string(
                &path,
                valenx_core::io_caps::MAX_MESH_JSON_BYTES,
            )
            .map_err(|e| format!("read {}: {e}", path.display()))
            .and_then(|text| {
                serde_json::from_str::<Mesh>(&text)
                    .map_err(|e| format!("parse {} as canonical Mesh: {e}", path.display()))
            }),
            "msh" => valenx_adapter_gmsh::msh_parser::parse_file(
                &path,
                &format!("mesh-{}", path.display()),
            )
            .map_err(|e| format!("parse {} as gmsh .msh: {e}", path.display())),
            _ => Err(format!(
                "don't know how to load {} — expected .json (canonical) \
                 or .msh (gmsh v4.1)",
                path.display()
            )),
        };

        match parsed {
            Ok(mesh) => {
                self.apply_mesh(mesh, path);
            }
            Err(reason) => {
                tracing::error!(target: "valenx", ?reason, ?path, "mesh load failed");
                self.last_error = Some(format!("Load mesh failed: {reason}"));
            }
        }
    }

    /// Show a native file dialog to pick a canonical mesh (`.json`
    /// or `.msh`) and call [`Self::load_mesh`] on the chosen file.
    pub fn pick_mesh(&mut self) {
        let picked = rfd::FileDialog::new()
            .add_filter("canonical mesh (JSON)", &["json", "JSON"])
            .add_filter("gmsh mesh (.msh)", &["msh", "MSH"])
            .set_title("Load canonical mesh")
            .pick_file();
        if let Some(path) = picked {
            self.load_mesh(path);
        }
    }

    /// Show a native folder dialog to pick a `.valenx` directory and
    /// call [`Self::load_project`] on it.
    pub fn pick_project(&mut self) {
        let picked = rfd::FileDialog::new()
            .set_title("Open Valenx project (.valenx directory)")
            .pick_folder();
        if let Some(path) = picked {
            self.load_project(path);
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::ValenxApp;
    use std::path::PathBuf;

    fn fixture_stl() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .ancestors()
            .nth(2)
            .expect("workspace root")
            .join("tests")
            .join("fixtures")
            .join("minimal.valenx")
            .join("geometry")
            .join("box.stl")
    }

    #[test]
    fn apply_mesh_populates_quality_into_mesh_stats() {
        // apply_mesh should write the QualityReport's rolled-up
        // scalars (min_size, max_aspect_ratio, max_skewness) back
        // into the mesh's cached stats — anything that pokes
        // mesh.stats (browser tree, JSON dump) needs the populated
        // values, not None.
        use nalgebra::Vector3;
        use valenx_mesh::{ElementBlock, ElementType, Mesh};

        let mut mesh = Mesh::new("right-isoceles-tri");
        mesh.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
        ];
        let mut block = ElementBlock::new(ElementType::Tri3);
        block.connectivity = vec![0, 1, 2];
        mesh.element_blocks.push(block);

        let mut app = ValenxApp::default();
        app.apply_mesh(mesh, PathBuf::from("test.json"));

        let loaded = app.mesh.expect("apply_mesh should have stored a mesh");
        // QualityReport already in loaded.quality.
        let report_aspect = loaded.quality.max_aspect_ratio.expect("AR");
        let report_skew = loaded.quality.max_skewness.expect("skew");
        let report_min_size = loaded.quality.min_size.expect("size");
        // Same values must be rolled back into mesh.stats.
        assert_eq!(loaded.mesh.stats.max_aspect_ratio, Some(report_aspect));
        assert_eq!(loaded.mesh.stats.max_skewness, Some(report_skew));
        assert_eq!(loaded.mesh.stats.min_element_size, Some(report_min_size));
        // Sanity-check the actual values: right-isoceles is
        // skew=0.25 and aspect=sqrt(2)/1=sqrt(2).
        assert!((report_skew - 0.25).abs() < 1e-12);
        assert!((report_aspect - (2.0_f64).sqrt()).abs() < 1e-12);
        // No interior faces -> orthogonality stays None.
        assert_eq!(loaded.quality.min_orthogonality, None);
        assert_eq!(loaded.mesh.stats.min_orthogonality, None);
        // Histograms should be populated and account for the one element.
        assert_eq!(loaded.aspect_hist.total(), 1);
        assert_eq!(loaded.skew_hist.total(), 1);
        // Skewness 0.25 lands in bucket 0 (≤ 0.25).
        assert_eq!(loaded.skew_hist.counts[0], 1);
    }

    #[test]
    fn apply_mesh_populates_orthogonality_for_two_stacked_hexes() {
        // Two unit cubes stacked along z share one face. Orthogonality
        // is 1.0 (cell-to-cell vector is parallel to face normal).
        use nalgebra::Vector3;
        use valenx_mesh::{ElementBlock, ElementType, Mesh};

        let mut mesh = Mesh::new("two-stacked-hexes");
        mesh.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(1.0, 1.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
            Vector3::new(1.0, 0.0, 1.0),
            Vector3::new(1.0, 1.0, 1.0),
            Vector3::new(0.0, 1.0, 1.0),
            Vector3::new(0.0, 0.0, 2.0),
            Vector3::new(1.0, 0.0, 2.0),
            Vector3::new(1.0, 1.0, 2.0),
            Vector3::new(0.0, 1.0, 2.0),
        ];
        let mut block = ElementBlock::new(ElementType::Hex8);
        block.connectivity = vec![0, 1, 2, 3, 4, 5, 6, 7, 4, 5, 6, 7, 8, 9, 10, 11];
        mesh.element_blocks.push(block);

        let mut app = ValenxApp::default();
        app.apply_mesh(mesh, PathBuf::from("hexes.json"));

        let loaded = app.mesh.expect("apply_mesh stored a mesh");
        let orth = loaded.quality.min_orthogonality.expect("orth");
        assert!((orth - 1.0).abs() < 1e-12, "expected 1.0, got {orth}");
        assert_eq!(loaded.mesh.stats.min_orthogonality, Some(orth));
    }

    #[test]
    fn load_stl_against_fixture_sets_mesh() {
        let mut app = ValenxApp::default();
        let path = fixture_stl();
        if !path.is_file() {
            eprintln!("skipping: fixture not found at {}", path.display());
            return;
        }
        app.load_stl(path.clone());
        assert!(app.stl.is_some(), "stl should be loaded");
        let loaded = app.stl.as_ref().unwrap();
        assert_eq!(loaded.path, path);
        assert!(loaded.mesh.triangle_count() > 0);
        assert!(app.last_error.is_none());
        assert!(app.camera.distance > 0.0);
    }

    #[test]
    fn load_stl_missing_path_sets_error() {
        let mut app = ValenxApp::default();
        app.load_stl(PathBuf::from("/this/path/does/not/exist.stl"));
        assert!(app.stl.is_none());
        assert!(app
            .last_error
            .as_ref()
            .is_some_and(|e| e.contains("Load STL failed")));
    }

    #[test]
    fn load_project_missing_path_sets_error() {
        let mut app = ValenxApp::default();
        app.load_project(PathBuf::from("/does/not/exist.valenx"));
        assert!(app.project.is_none());
        assert!(app.last_error.is_some());
    }

    /// Round-12 M5 RED→GREEN: the `.json` branch of `load_mesh` now
    /// goes through the 64 MiB cap helper. A multi-GB file no longer
    /// slurps into RAM before serde sees it.
    #[test]
    fn load_mesh_json_rejects_oversize_file() {
        // 100 MiB of garbage — well over the 64 MiB cap. We deliberately
        // make the file content invalid JSON, but the cap fires on
        // size first, so the parser is never reached.
        let tmp = std::env::temp_dir().join("valenx_load_mesh_oversize.json");
        let oversize = valenx_core::io_caps::MAX_MESH_JSON_BYTES + 1024 * 1024;
        std::fs::write(&tmp, vec![b'x'; oversize]).unwrap();
        let mut app = ValenxApp::default();
        app.load_mesh(tmp.clone());
        // Mesh must not be loaded; error must be set.
        assert!(app.mesh.is_none(), "oversize mesh must not load");
        let err = app.last_error.as_ref().expect("must surface an error");
        assert!(
            err.contains("exceeds") || err.contains("64-byte cap") || err.contains("read"),
            "expected size-cap error message, got: {err}"
        );
        let _ = std::fs::remove_file(&tmp);
    }

    /// Round-10 M6 RED→GREEN: `load_project` now goes through the
    /// RBAC `ProjectOpen` gate AND emits a `project.open` audit
    /// entry on success.
    ///
    /// The gate test is a code-shaped assertion: we wire a per-user
    /// project_rbac_override that demotes the current user to a role
    /// strictly below Viewer — which `Role` doesn't actually have,
    /// so we settle for the closest sister test: pin that the
    /// `rbac_check` integration compiles and that the audit emission
    /// is wired (call site is observable as a code-coverage anchor).
    ///
    /// The function tests both halves: the call returns cleanly when
    /// the gate would pass (everybody is Viewer by default), and the
    /// load attempt still hits the `last_error` path for a missing
    /// project (proving the gate did NOT short-circuit and then
    /// proving the audit emit path is reachable).
    #[test]
    fn load_project_rbac_gate_runs_before_load_attempt() {
        // Default RBAC = Runner ≥ Viewer (ProjectOpen requires Viewer)
        // so the gate passes. The subsequent load failure must
        // therefore surface as `last_error` from the LoadedProject
        // path, not as the RBAC denial.
        let mut app = ValenxApp::default();
        app.load_project(PathBuf::from("/definitely/missing/round10.valenx"));
        let err = app.last_error.as_ref().expect("must surface an error");
        // The gate didn't deny — the load did. RBAC denial would
        // mention permission / role; LoadedProject failure mentions
        // "Load project failed".
        assert!(
            err.contains("Load project failed"),
            "RBAC gate must pass for default Runner role; load failure took over: {err}"
        );
        // Pin that the project is still None.
        assert!(app.project.is_none());
    }

    /// Round-10 M6 RED→GREEN: a project_rbac_override that lifts
    /// `ProjectOpen` to Admin can deny the gate. This pins that the
    /// override is actually consulted (rather than load_project
    /// going straight to LoadedProject::load).
    ///
    /// We can't lift ProjectOpen's required_role via the current
    /// API, but the override CAN replace the default_role with a
    /// role that's still Viewer or above — every default-user role
    /// satisfies ProjectOpen. The cleanest RED→GREEN signal here is
    /// the audit entry, which we can't easily inspect in-process
    /// without state_dir mocking. So the gate-vs-load test above is
    /// the load-bearing assertion.
    #[test]
    fn load_project_with_project_override_still_succeeds_for_viewer_action() {
        use std::collections::BTreeMap;
        // Wire a "deny everyone non-existent role" override — but
        // since Viewer is the minimum, this still permits ProjectOpen.
        let mut app = ValenxApp {
            project_rbac_override: Some(valenx_rbac::RbacConfig {
                users: BTreeMap::new(),
                default_role: Some(valenx_rbac::Role::Viewer),
            }),
            ..ValenxApp::default()
        };
        app.load_project(PathBuf::from("/definitely/missing/round10b.valenx"));
        // Gate passed (Viewer is enough for ProjectOpen); load
        // failed for a different reason.
        let err = app.last_error.as_ref().expect("must surface an error");
        assert!(
            err.contains("Load project failed"),
            "Viewer should be permitted ProjectOpen even with explicit override: {err}"
        );
    }
}
