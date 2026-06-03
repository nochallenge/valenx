//! Single-case run orchestration methods on [`ValenxApp`]. Split out
//! of `lib.rs` as part of the structural refactor.

use std::sync::mpsc::TryRecvError;

use valenx_core::{Case, LogLevel};

use crate::audit::emit_audit;
use crate::file_browser::{open_path_or_copy, POPUP_DISABLED_PREFIX};
use crate::history::save_run_history_to_state_dir;
use crate::mesh_loader::{latest_snapshot_in_workdir, load_mesh_from_vtk};
use crate::rbac_io::rbac_check;
#[allow(unused_imports)] // `RunHandle` is referenced from a doc link below
use crate::run::RunHandle;
use crate::run::{self, RunEvent};
use crate::solver_parse::adapter_id_from_solver;
use crate::types::RunHistoryEntry;
use crate::ValenxApp;

impl ValenxApp {
    /// Back-compat: runs the first case in the project. Kept because
    /// the command palette + older UI paths still reference it.
    pub fn run_first_case(&mut self) {
        if let Some(project) = &self.project {
            if let Some(first) = project.case_names().first() {
                self.selected_case = Some(first.clone());
                self.run_selected_case();
                return;
            }
        }
        self.last_error = Some("Load a .valenx project with at least one case first.".into());
    }

    /// Run the currently-selected case forced through a specific
    /// adapter (chosen by the command palette's per-adapter `Run:`
    /// entry). Validates that the case's solver actually maps to
    /// `adapter_id`; surfaces a clear `last_error` on mismatch
    /// rather than running a case the user didn't expect.
    pub fn run_selected_case_with_adapter(&mut self, adapter_id: &str) {
        let project = match &self.project {
            Some(p) => p,
            None => {
                self.last_error = Some("Load a .valenx project first.".into());
                return;
            }
        };
        let case_name = match self
            .selected_case
            .clone()
            .or_else(|| project.case_names().first().cloned())
        {
            Some(n) => n,
            None => {
                self.last_error = Some("Project has no cases.".into());
                return;
            }
        };
        let case_def = match project.cases.get(&case_name) {
            Some(cd) => cd,
            None => {
                self.last_error = Some(format!(
                    "Case `{case_name}` referenced but not loaded in the project."
                ));
                return;
            }
        };
        let solver_adapter = adapter_id_from_solver(&case_def.case.solver);
        if solver_adapter != adapter_id {
            self.last_error = Some(format!(
                "Case `{case_name}` has solver `{}` (adapter `{solver_adapter}`), \
                 not `{adapter_id}`. Edit the case's `solver` field to switch \
                 adapters, or pick a case that already targets `{adapter_id}`.",
                case_def.case.solver
            ));
            return;
        }
        // Solver matches — delegate to the regular run path which
        // handles registry lookup, RBAC, audit, workdir creation,
        // and run-handle bookkeeping.
        self.run_selected_case();
    }

    /// Open a folder picker and scaffold a new `.valenx` project
    /// rooted at the chosen directory whose default case targets
    /// `adapter_id`. Drives the same template library
    /// (`valenx-init` uses) so palette + ribbon + CLI converge on
    /// one source of truth. Sets `status` on success (and auto-loads
    /// the project so the user can run it immediately) or
    /// `last_error` on any IO / unknown-adapter failure.
    pub fn new_case_for_adapter(&mut self, adapter_id: &str) {
        let case_dir = match valenx_core::init_templates::case_dir_for(adapter_id) {
            Some(d) => d,
            None => {
                self.last_error = Some(format!(
                    "No starter template for adapter `{adapter_id}` — bring \
                     your own case.toml and load it via File → Open project."
                ));
                return;
            }
        };
        let case_body = match valenx_core::init_templates::case_toml_body(adapter_id) {
            Some(b) => b,
            None => {
                self.last_error = Some(format!(
                    "Internal: case_dir resolved but case.toml body missing for `{adapter_id}`."
                ));
                return;
            }
        };

        let picked = rfd::FileDialog::new()
            .set_title(format!(
                "New {adapter_id} project — pick destination directory"
            ))
            .pick_folder();
        let dir = match picked {
            Some(d) => d,
            None => return, // user cancelled — silent
        };

        let project_toml_path = dir.join("project.toml");
        if project_toml_path.exists() {
            self.last_error = Some(format!(
                "{} already exists — pick an empty directory or load it via File → Open project.",
                project_toml_path.display()
            ));
            return;
        }
        let project_name = dir
            .file_name()
            .and_then(|s| s.to_str())
            .unwrap_or("untitled");
        if let Err(e) = std::fs::create_dir_all(&dir) {
            self.last_error = Some(format!("create {}: {e}", dir.display()));
            return;
        }
        let rendered = match valenx_core::init_templates::project_toml(project_name, case_dir) {
            Ok(s) => s,
            Err(e) => {
                self.last_error = Some(format!(
                    "Invalid project / case name: {e}. Pick a destination directory whose name \
                     uses only ASCII letters, digits, `_`, `.`, or `-`."
                ));
                return;
            }
        };
        if let Err(e) = valenx_core::io_caps::atomic_write_str(&project_toml_path, &rendered) {
            self.last_error = Some(format!("write {}: {e}", project_toml_path.display()));
            return;
        }
        let case_dir_path = dir.join("cases").join(case_dir);
        if let Err(e) = std::fs::create_dir_all(&case_dir_path) {
            self.last_error = Some(format!("create {}: {e}", case_dir_path.display()));
            return;
        }
        let case_toml_path = case_dir_path.join("case.toml");
        if let Err(e) = valenx_core::io_caps::atomic_write_str(&case_toml_path, &case_body) {
            self.last_error = Some(format!("write {}: {e}", case_toml_path.display()));
            return;
        }
        self.status = Some(format!(
            "Scaffolded {adapter_id} project at {}",
            dir.display()
        ));
        self.last_error = None;
        // Auto-load so the user lands inside the freshly-scaffolded
        // project without an extra File → Open click.
        self.load_project(dir);
    }

    /// Kick off a run of the currently-selected case. Resolves the
    /// adapter by parsing the case's `solver` field (prefix before
    /// the first dot — e.g. `"openfoam.simpleFoam"` → `"openfoam"`).
    /// Non-blocking: the work happens on a background thread, the
    /// UI polls events from the returned [`RunHandle`].
    pub fn run_selected_case(&mut self) {
        if let Err(e) = rbac_check(
            valenx_rbac::Action::RunCase,
            self.project_rbac_override.as_ref(),
        ) {
            self.last_error = Some(format!("{e}"));
            return;
        }
        if self.run_handle.is_some() {
            self.last_error = Some("A run is already in progress.".into());
            return;
        }
        let project = match &self.project {
            Some(p) => p,
            None => {
                self.last_error = Some("Load a .valenx project first.".into());
                return;
            }
        };

        // Resolve the case name: explicit selection, else first.
        let case_name = match self
            .selected_case
            .clone()
            .or_else(|| project.case_names().first().cloned())
        {
            Some(n) => n,
            None => {
                self.last_error = Some("Project has no cases.".into());
                return;
            }
        };
        let case_def = match project.cases.get(&case_name) {
            Some(cd) => cd,
            None => {
                self.last_error = Some(format!(
                    "Case `{case_name}` referenced but not loaded in the project."
                ));
                return;
            }
        };

        let adapter_id = adapter_id_from_solver(&case_def.case.solver);
        let adapter = match self.registry.get(adapter_id) {
            Some(e) if e.status.is_ready() => e.adapter.clone(),
            Some(_) => {
                self.last_error = Some(format!(
                    "`{adapter_id}` adapter is registered but not Ready — \
                     re-probe from Settings after installing the tool."
                ));
                return;
            }
            None => {
                self.last_error = Some(format!(
                    "No adapter registered for solver `{}` (looked up id `{adapter_id}`).",
                    case_def.case.solver
                ));
                return;
            }
        };

        let case_path = project.root.join("cases").join(&case_name);
        let workdir = std::env::temp_dir().join(format!(
            "valenx-run-{}-{}",
            case_name,
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .map(|d| d.as_secs())
                .unwrap_or(0)
        ));
        if let Err(e) = std::fs::create_dir_all(&workdir) {
            self.last_error = Some(format!("create workdir: {e}"));
            return;
        }

        self.residuals.clear();
        self.last_run_report = None;
        self.last_run_error = None;
        self.run_progress = 0.0;
        self.run_message = format!("starting {adapter_id}…");

        let case = Case {
            id: case_name.clone(),
            path: case_path,
        };
        self.running_case_name = Some(case_name.clone());
        // Audit: who started which case via which adapter.
        emit_audit(
            "run.start",
            serde_json::json!({"kind": "case", "case": case_name}),
            serde_json::json!({"adapter": adapter_id}),
        );
        self.run_handle = Some(run::spawn(adapter, case, workdir));
        self.status = Some(format!("Running case `{case_name}` via `{adapter_id}`"));
    }

    /// Cancel the currently-running adapter (subject to the RBAC
    /// `CancelRun` action). Records an audit-log entry and surfaces
    /// errors via `last_error` rather than panicking.
    pub fn cancel_run(&mut self) {
        if let Err(e) = rbac_check(
            valenx_rbac::Action::CancelRun,
            self.project_rbac_override.as_ref(),
        ) {
            self.last_error = Some(format!("{e}"));
            return;
        }
        if let Some(h) = &self.run_handle {
            h.cancel();
            self.run_message = "cancelling…".into();
            // Audit: who cancelled which running case.
            if let Some(name) = &self.running_case_name {
                emit_audit(
                    "run.cancel",
                    serde_json::json!({"kind": "case", "case": name}),
                    serde_json::json!({}),
                );
            }
        }
    }

    /// Open the last prepared workdir in the host's file browser
    /// (Explorer / Finder / xdg-open). No-ops with a structured
    /// `last_error` when nothing's been prepared yet — the menu
    /// item is gated on the same condition, but the command palette
    /// can still fire this without a workdir present.
    pub fn open_prepare_workdir(&mut self) {
        let Some(path) = self.last_prepare_workdir.clone() else {
            self.last_error =
                Some("Nothing prepared yet — run `Run → Prepare selected case` first.".into());
            return;
        };
        let disable = self.settings.disable_file_browser_popups;
        match open_path_or_copy(&path, disable) {
            Ok(()) => {}
            Err(reason) if reason.starts_with(POPUP_DISABLED_PREFIX) => {
                // Kill-switch path: not an error — surface the path
                // as a neutral status line so the user gets the
                // workdir without an Explorer popup.
                self.status = Some(reason);
            }
            Err(reason) => {
                // Surface the workdir alongside the error so the user has
                // a fallback path even when the OS launcher is broken
                // (headless Linux without xdg-open, etc.).
                self.last_error = Some(format!(
                    "Couldn't open file browser: {reason}. The workdir is at {}",
                    path.display()
                ));
            }
        }
    }

    /// Open the last run's workdir in the host's file browser. The
    /// run workdir holds the solver's actual output — `.vtu` / `.frd`
    /// / `log.simpleFoam` / etc. — so this is the post-run "where
    /// did the results go?" answer.
    pub fn open_run_workdir(&mut self) {
        let Some(path) = self.last_run_workdir.clone() else {
            self.last_error =
                Some("No completed run yet — finish a run before asking for its workdir.".into());
            return;
        };
        let disable = self.settings.disable_file_browser_popups;
        match open_path_or_copy(&path, disable) {
            Ok(()) => {}
            Err(reason) if reason.starts_with(POPUP_DISABLED_PREFIX) => {
                // Kill-switch path: not an error — surface the path
                // as a neutral status line.
                self.status = Some(reason);
            }
            Err(reason) => {
                self.last_error = Some(format!(
                    "Couldn't open file browser: {reason}. The workdir is at {}",
                    path.display()
                ));
            }
        }
    }

    /// Prepare the currently-selected case without running it. Mirrors
    /// [`Self::run_selected_case`] up to adapter resolution, then calls
    /// `adapter.prepare()` synchronously and stops. The emitted dict
    /// tree (or whatever artifacts the adapter writes during prepare)
    /// lives in a temp workdir and the path is stashed in
    /// `last_prepare_workdir` so the UI can surface it.
    ///
    /// Useful when the user has a `case.toml` they want to inspect or
    /// edit before launching the solver — or when they don't have the
    /// solver installed at all and just want to see what the adapter
    /// would write.
    pub fn prepare_selected_case(&mut self) {
        if let Err(e) = rbac_check(
            valenx_rbac::Action::PrepareCase,
            self.project_rbac_override.as_ref(),
        ) {
            self.last_error = Some(format!("{e}"));
            return;
        }
        if self.run_handle.is_some() {
            self.last_error = Some(
                "A run is already in progress — cancel it before preparing a new case.".into(),
            );
            return;
        }
        let project = match &self.project {
            Some(p) => p,
            None => {
                self.last_error = Some("Load a .valenx project first.".into());
                return;
            }
        };

        // Resolve the case name: explicit selection, else first.
        let case_name = match self
            .selected_case
            .clone()
            .or_else(|| project.case_names().first().cloned())
        {
            Some(n) => n,
            None => {
                self.last_error = Some("Project has no cases.".into());
                return;
            }
        };
        let case_def = match project.cases.get(&case_name) {
            Some(cd) => cd,
            None => {
                self.last_error = Some(format!(
                    "Case `{case_name}` referenced but not loaded in the project."
                ));
                return;
            }
        };

        let adapter_id = adapter_id_from_solver(&case_def.case.solver);
        let adapter = match self.registry.get(adapter_id) {
            Some(e) if e.status.is_ready() => e.adapter.clone(),
            Some(_) => {
                self.last_error = Some(format!(
                    "`{adapter_id}` adapter is registered but not Ready — \
                     re-probe from Settings after installing the tool."
                ));
                return;
            }
            None => {
                self.last_error = Some(format!(
                    "No adapter registered for solver `{}` (looked up id `{adapter_id}`).",
                    case_def.case.solver
                ));
                return;
            }
        };

        let case_path = project.root.join("cases").join(&case_name);
        let workdir = std::env::temp_dir().join(format!(
            "valenx-prepare-{}-{}",
            case_name,
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .map(|d| d.as_secs())
                .unwrap_or(0)
        ));
        if let Err(e) = std::fs::create_dir_all(&workdir) {
            self.last_error = Some(format!("create workdir: {e}"));
            return;
        }

        let case = Case {
            id: case_name.clone(),
            path: case_path,
        };
        match adapter.prepare(&case, &workdir) {
            Ok(prepared) => {
                self.last_prepare_workdir = Some(workdir.clone());
                // Stash the prepared job so the user can later click
                // "Run from prepared workdir" — that path skips the
                // prepare step entirely so any hand-edits to the dicts
                // survive.
                self.last_prepared_job = Some((adapter_id.to_string(), prepared));
                self.last_error = None;
                self.status = Some(format!(
                    "Prepared `{case_name}` via `{adapter_id}` → {}",
                    workdir.display()
                ));
                emit_audit(
                    "prepare.start",
                    serde_json::json!({"kind": "case", "case": case_name}),
                    serde_json::json!({
                        "adapter": adapter_id,
                        "workdir": workdir.display().to_string(),
                        "result": "ok",
                    }),
                );
            }
            Err(e) => {
                // Many adapters (notably OpenFOAM) write the dict tree
                // BEFORE the binary lookup, so a ToolNotInstalled at
                // this stage is meaningful: the dicts are on disk for
                // inspection. Surface the workdir alongside the error
                // so users can still find what was generated.
                self.last_prepare_workdir = Some(workdir.clone());
                // No PreparedJob to stash — the caller should re-run
                // prepare once the underlying issue is fixed.
                self.last_prepared_job = None;
                self.last_error = Some(format!(
                    "Prepare for `{case_name}` failed: {e} (partial output may be in {})",
                    workdir.display()
                ));
                emit_audit(
                    "prepare.start",
                    serde_json::json!({"kind": "case", "case": case_name}),
                    serde_json::json!({
                        "adapter": adapter_id,
                        "workdir": workdir.display().to_string(),
                        "result": "failed",
                        "error": e.to_string(),
                    }),
                );
            }
        }
    }

    /// Spawn the solver against the workdir produced by the most
    /// recent successful [`Self::prepare_selected_case`], skipping
    /// the prepare step. The user's edits to the dict files inside
    /// the workdir survive — that's the whole point of this
    /// "prepare → edit → run" workflow.
    pub fn run_from_prepared_workdir(&mut self) {
        if let Err(e) = rbac_check(
            valenx_rbac::Action::RunCase,
            self.project_rbac_override.as_ref(),
        ) {
            self.last_error = Some(format!("{e}"));
            return;
        }
        if self.run_handle.is_some() {
            self.last_error = Some("A run is already in progress.".into());
            return;
        }
        let Some((adapter_id, prepared)) = self.last_prepared_job.clone() else {
            self.last_error = Some(
                "Nothing prepared yet — click `Run → Prepare selected case` first, \
                 then optionally edit the generated files in the workdir."
                    .into(),
            );
            return;
        };
        let adapter = match self.registry.get(&adapter_id) {
            Some(e) if e.status.is_ready() => e.adapter.clone(),
            Some(_) => {
                self.last_error = Some(format!(
                    "`{adapter_id}` adapter is no longer Ready — re-probe from \
                     Settings, then re-prepare the case."
                ));
                return;
            }
            None => {
                self.last_error = Some(format!(
                    "Adapter `{adapter_id}` was unregistered since the last prepare \
                     — re-prepare to get a fresh job."
                ));
                return;
            }
        };

        self.residuals.clear();
        self.last_run_report = None;
        self.last_run_error = None;
        self.run_progress = 0.0;
        self.run_message = format!("re-running {adapter_id} from prepared workdir…");

        let workdir_display = prepared.workdir.display().to_string();
        // Re-runs use the same selected_case key for history tracking.
        // If the user re-prepares against a different case the history
        // entry will be overwritten — last-run-wins is the intended
        // semantic.
        self.running_case_name = self.selected_case.clone();
        self.run_handle = Some(run::spawn_prepared(adapter, prepared));
        self.status = Some(format!(
            "Re-running prepared workdir via `{adapter_id}` → {workdir_display}"
        ));
    }

    /// Drain the run channel. Runs at most `MAX_EVENTS_PER_FRAME`
    /// events so a chatty log doesn't starve the UI thread.
    pub(crate) fn pump_run_events(&mut self) {
        const MAX_EVENTS_PER_FRAME: usize = 256;
        let mut drained = 0;
        let mut finished = false;

        if let Some(h) = &self.run_handle {
            while drained < MAX_EVENTS_PER_FRAME {
                match h.rx.try_recv() {
                    Ok(ev) => {
                        drained += 1;
                        match ev {
                            RunEvent::Starting => {
                                self.run_message = "starting…".into();
                            }
                            RunEvent::Progress { pct, message } => {
                                self.run_progress = pct;
                                self.run_message = message;
                            }
                            RunEvent::LogLine { level, line } => {
                                self.residuals.ingest_log_line(&line);
                                self.log.push(level, line);
                            }
                            RunEvent::Finished(report) => {
                                self.run_progress = 100.0;
                                self.run_message = format!(
                                    "done — exit {} in {:?}",
                                    report.exit_code, report.wall_time
                                );
                                // Stash a history entry under the
                                // running case name so the browser
                                // can show ✓/✗ next to the case row.
                                // Persist the map to disk so the badge
                                // survives an app restart.
                                if let Some(name) = self.running_case_name.clone() {
                                    self.run_history.insert(
                                        name,
                                        RunHistoryEntry {
                                            succeeded: report.exit_code == 0,
                                            wall_time: report.wall_time,
                                            converged: report.converged,
                                        },
                                    );
                                    save_run_history_to_state_dir(&self.run_history);
                                }
                                self.last_run_report = Some(report);
                                // Don't set finished=true yet —
                                // Collected may still be on its way.
                                // pump_run_events catches the channel
                                // disconnect when the worker thread
                                // exits, which is what actually
                                // triggers the join+cleanup below.
                            }
                            RunEvent::Collected(results) => {
                                let n_fields = results.fields.len();
                                let n_artifacts = results.artifacts.len();
                                let n_scalars = results.scalars.len();
                                // Persist a sidecar `results.json` next
                                // to the run's workdir so the sweep
                                // exporter (and anything else that
                                // wants to consume Results across
                                // app restarts) can pick it up
                                // without re-running the adapter.
                                if let Some(handle) = &self.run_handle {
                                    let target = handle.workdir.join("results.json");
                                    match serde_json::to_string_pretty(&*results) {
                                        Ok(json) => {
                                            if let Err(e) =
                                                valenx_core::io_caps::atomic_write_str(&target, &json)
                                            {
                                                self.log.push(
                                                    LogLevel::Warn,
                                                    format!("persist {}: {e}", target.display()),
                                                );
                                            }
                                        }
                                        Err(e) => {
                                            self.log.push(
                                                LogLevel::Warn,
                                                format!("serialise results: {e}"),
                                            );
                                        }
                                    }
                                }
                                self.last_run_results = Some(results);
                                if n_fields > 0 || n_scalars > 0 {
                                    self.run_message = format!(
                                        "done — {n_fields} fields, {n_scalars} scalars, {n_artifacts} artifacts collected"
                                    );
                                }
                            }
                            RunEvent::Failed(msg) => {
                                // Failure: record the case as
                                // unsucessful so the badge in the
                                // browser shows ✗. Persist the map
                                // so the badge survives an app
                                // restart.
                                if let Some(name) = self.running_case_name.clone() {
                                    self.run_history.insert(
                                        name,
                                        RunHistoryEntry {
                                            succeeded: false,
                                            wall_time: std::time::Duration::ZERO,
                                            converged: None,
                                        },
                                    );
                                    save_run_history_to_state_dir(&self.run_history);
                                }
                                self.last_run_error = Some(msg.clone());
                                self.run_message = format!("failed: {msg}");
                                finished = true;
                            }
                        }
                    }
                    Err(TryRecvError::Empty) => break,
                    Err(TryRecvError::Disconnected) => {
                        finished = true;
                        break;
                    }
                }
            }
        }

        if finished {
            // Cleanly join the thread so the handle doesn't leak.
            // Capture the adapter-id + workdir so post-run hooks can
            // dispatch on what just ran before the handle drops.
            let completed_adapter_id;
            let completed_workdir;
            if let Some(mut handle) = self.run_handle.take() {
                completed_adapter_id = Some(handle.adapter_id);
                completed_workdir = Some(handle.workdir.clone());
                if let Some(j) = handle.thread.take() {
                    let _ = j.join();
                }
            } else {
                completed_adapter_id = None;
                completed_workdir = None;
            }

            if let (Some(id), Some(workdir)) = (completed_adapter_id, completed_workdir) {
                // Stash the run workdir so the user can "Open in file
                // browser" it from the Results pane. We set this even
                // when the run failed — failed runs typically still
                // leave a partial dict tree + log on disk, and the
                // user often wants to dig in to figure out why.
                self.last_run_workdir = Some(workdir.clone());
                self.on_run_finished(id, &workdir);
            }
            // Clear running_case_name now that the run is fully
            // wound up. Subsequent runs set it again at spawn time.
            self.running_case_name = None;
        }
    }

    /// Post-run hook. Called once per completed run (success or
    /// failure) so the app can auto-open produced artifacts.
    ///
    /// Branches:
    /// - `gmsh` / `netgen`: load `mesh.canonical.json` (gmsh falls
    ///   back to a raw `mesh.msh` if the canonical serialisation
    ///   somehow didn't land) into the viewport.
    /// - `freecad`: load the exported `output.stl` so STEP imports
    ///   are immediately visible.
    /// - `openfoam`: parse the latest `.vtu` in the workdir and load
    ///   its mesh into the viewport — gives the user instant visual
    ///   confirmation that the solver wrote real geometry. Field-
    ///   colored rendering is the next layer up.
    fn on_run_finished(&mut self, adapter_id: &str, workdir: &std::path::Path) {
        match adapter_id {
            "gmsh" | "netgen" => {
                let canonical = workdir.join("mesh.canonical.json");
                if canonical.is_file() {
                    tracing::info!(
                        target: "valenx",
                        ?canonical,
                        adapter = adapter_id,
                        "auto-loading canonical mesh from finished mesh-adapter run"
                    );
                    self.load_mesh(canonical);
                } else if adapter_id == "gmsh" {
                    // Fall back to the raw .msh if the canonical
                    // serialisation somehow didn't land. Netgen
                    // doesn't have a canonical loader for `.vol`
                    // outside the adapter's collect() yet, so this
                    // fallback is gmsh-only.
                    let raw = workdir.join("mesh.msh");
                    if raw.is_file() {
                        tracing::info!(
                            target: "valenx",
                            ?raw,
                            "canonical mesh not found — auto-loading raw .msh"
                        );
                        self.load_mesh(raw);
                    }
                }
            }
            "freecad" => {
                let exported_stl = workdir.join("output.stl");
                if exported_stl.is_file() {
                    tracing::info!(
                        target: "valenx",
                        ?exported_stl,
                        "auto-loading FreeCAD export into viewport"
                    );
                    self.load_stl(exported_stl);
                }
            }
            "openfoam" | "elmer" | "su2" | "calculix" | "code-aster" => {
                // Snapshot auto-load: prefer a `.pvd` time-series
                // manifest if the solver wrote one (it's the curated
                // multi-step view), falling back to the lexicographic
                // latest `.vtu`/`.vtk` walk for solvers that emit raw
                // snapshots without a manifest.
                if let Some(snap_path) = latest_snapshot_in_workdir(workdir) {
                    match load_mesh_from_vtk(&snap_path) {
                        Ok(mesh) => {
                            tracing::info!(
                                target: "valenx",
                                ?snap_path,
                                "auto-loading mesh from latest snapshot"
                            );
                            self.apply_mesh(mesh, snap_path);
                        }
                        Err(reason) => {
                            tracing::warn!(
                                target: "valenx",
                                ?snap_path,
                                ?reason,
                                "could not auto-load snapshot into viewport"
                            );
                        }
                    }
                }
            }
            _ => {}
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::ValenxApp;

    #[test]
    fn run_selected_without_project_errors_cleanly() {
        let mut app = ValenxApp::default();
        app.run_selected_case();
        assert!(app.run_handle.is_none());
        assert!(
            app.last_error
                .as_ref()
                .is_some_and(|e| e.contains("project")),
            "got {:?}",
            app.last_error
        );
    }

    #[test]
    fn prepare_selected_without_project_errors_cleanly() {
        let mut app = ValenxApp::default();
        app.prepare_selected_case();
        // No workdir created and no run started.
        assert!(app.last_prepare_workdir.is_none());
        assert!(app.run_handle.is_none());
        assert!(
            app.last_error
                .as_ref()
                .is_some_and(|e| e.contains("project")),
            "got {:?}",
            app.last_error
        );
    }

    #[test]
    fn open_prepare_workdir_without_prepare_errors_cleanly() {
        let mut app = ValenxApp::default();
        app.open_prepare_workdir();
        assert!(
            app.last_error
                .as_ref()
                .is_some_and(|e| e.contains("prepared")),
            "got {:?}",
            app.last_error
        );
    }

    #[test]
    fn open_run_workdir_without_run_errors_cleanly() {
        let mut app = ValenxApp::default();
        app.open_run_workdir();
        assert!(
            app.last_error
                .as_ref()
                .is_some_and(|e| e.contains("completed run")),
            "got {:?}",
            app.last_error
        );
    }

    #[test]
    fn run_from_prepared_without_prepare_errors_cleanly() {
        let mut app = ValenxApp::default();
        app.run_from_prepared_workdir();
        assert!(app.run_handle.is_none());
        assert!(
            app.last_error
                .as_ref()
                .is_some_and(|e| e.contains("Nothing prepared")),
            "got {:?}",
            app.last_error
        );
    }

    #[test]
    fn run_without_project_surfaces_error() {
        let mut app = ValenxApp::default();
        app.run_first_case();
        assert!(app.run_handle.is_none());
        assert!(app
            .last_error
            .as_ref()
            .is_some_and(|e| e.contains("project")));
    }
}
