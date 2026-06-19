//! # valenx-app
//!
//! The Valenx desktop application. Owns the main window, the ribbon,
//! the browser tree, the viewport host, the residual panel, the
//! command palette, and the run orchestration that ties them
//! together.
//!
//! Design specified in [DESIGN.md § 6 Screens](../DESIGN.md#6-screens).
//! Architecture positioned in [ARCHITECTURE.md § 2](../ARCHITECTURE.md).
//!
//! ## Year-1 shell (Phase 1)
//!
//! Four panels + one overlay:
//!
//! 1. **Ribbon** on top with File / View / Run / Help menus.
//! 2. **Browser tree** on the left listing project, cases, geometry,
//!    results, and — powered by the live [`AdapterRegistry`] — every
//!    registered adapter with a status colour.
//! 3. **Viewport** in the centre, an interactive **shaded** 3D render
//!    with back-face culling and a dot-product light model (or
//!    wireframe, toggleable). `wgpu` swap is a rasteriser change,
//!    not an API rewrite — the projection math already ships in
//!    `valenx-viz::projection`.
//! 4. **Residual chart** at the bottom using `egui_plot`, live-fed
//!    from the OpenFOAM solver's stdout through the channel bridge
//!    in the `run` module.
//! 5. **Command palette** overlay (Ctrl+P) with fuzzy-search over
//!    every action the app knows how to perform.
//!
//! End-to-end CFD thread works: Open a `.valenx` project → Run →
//! the first case prepares an OpenFOAM simpleFoam deck, spawns the
//! solver on a background thread, streams residuals into the chart,
//! and reports the final convergence state.

#![forbid(unsafe_code)]
#![allow(missing_docs)] // relaxed during pre-alpha; see valenx-fields for rationale

pub mod aero;
pub mod aero_workbench;
pub mod animate_workbench;
pub(crate) mod background;
pub mod cad_workbench;
pub mod car_workbench;
pub mod cfd_workbench;
pub mod draft2d_workbench;
pub mod fem_workbench;
pub mod headless;
pub mod hvac_workbench;
pub mod interior_workbench;
pub mod neuro_workbench;
pub mod reinforcement_workbench;
pub mod render_workbench;
pub mod reverse_workbench;
pub mod variant_effect_workbench;
pub mod windturbine_workbench;

pub mod assistant_workbench;
pub mod astro;
pub mod astro_workbench;
pub mod cam_overlay;
pub mod collision_workbench;
pub mod commands;
#[cfg(test)]
mod coverage_ui_tests;
pub mod docking;
pub mod draft_overlay;
pub mod drone_workbench;
pub mod engine_workbench;
pub mod fasteners_workbench;
pub mod fields_workbench;
pub mod first_run;
pub mod fixedwing_workbench;
pub mod fourbar_workbench;
pub mod frames_workbench;
pub mod gasdynamics_workbench;
pub mod gears_workbench;
pub mod genetics;
pub mod genetics_workbench;
pub mod geomatics_workbench;
pub mod heatexchanger_workbench;
pub mod keyboard_help;
pub mod landing_page;
pub mod log_panel;
pub mod marine_workbench;
pub mod mesh_toolbox;
pub mod new_project_dialog;
pub mod panel_help;
pub mod pbr_forward_pass;
pub mod piping_workbench;
pub mod project_tabs;
pub mod rail_workbench;
pub mod reactdyn_workbench;
pub mod residuals;
pub mod rocket_mesh;
pub mod rocket_workbench;
pub mod run;
pub mod scene_overlay;
pub mod settings;
pub mod setup;
pub mod sheetmetal_workbench;
pub mod shortcuts;
pub mod sketch_overlay;
pub mod solarpv_workbench;
pub mod springs_workbench;
pub mod theme;
pub mod tooltips;
pub mod types;
pub mod undo;
pub mod viewport;
pub mod viewport_2d;
pub mod viewport_kind;
pub mod welcome_tour;
pub mod wgpu_renderer;
pub mod workbench_ui;

// Concern-focused helper modules — what used to be a single
// `helpers.rs` bag-of-everything (1422 LOC, 36 fns spanning 8+
// unrelated concerns). Sibling modules let callers `use
// crate::history::save_run_history_to_state_dir` and have the
// import name actually tell them which concern they're reaching
// into.
pub(crate) mod adapter_status;
pub mod audit;
pub mod file_browser;
pub(crate) mod histograms;
pub mod history;
pub(crate) mod mesh_loader;
pub mod rbac_io;
pub mod settings_io;
pub mod solver_parse;
pub mod state_paths;
pub mod time_format;

// Concern-focused impl ValenxApp blocks split out of this file to
// keep the root module readable. Each module holds one `impl
// ValenxApp` (and, for `update`, one `impl eframe::App for
// ValenxApp`). The split is structural only — methods keep their
// existing signatures + visibility — so callers compose the same
// regardless of which module a method actually lives in.
mod audit_history;
mod loading;
mod run_actions;
mod sweep;
mod update;
mod view_actions;

use std::path::PathBuf;
use std::sync::Arc;

use valenx_core::{AdapterRegistry, AdapterStatus, LoadedProject, RunReport};
use valenx_viz::OrbitCamera;

use crate::commands::CommandPalette;
use crate::log_panel::LogPanel;
use crate::residuals::ResidualHistory;
use crate::run::{RunHandle, SweepHandle};
use crate::settings::Settings;
use crate::viewport::ShadingMode;
use crate::wgpu_renderer::WgpuRenderer;

// Re-export the moved items so the public API (`use
// valenx_app::adapter_id_from_solver;`, `use valenx_app::run;`,
// `use valenx_app::BottomTab;`, …) stays stable across the
// extraction. The implementations now live in concern-focused
// sibling modules.
pub use crate::audit::emit_audit;
pub use crate::file_browser::open_path_in_file_browser;
pub use crate::history::{
    load_run_history_from_state_dir, load_sweep_history_from_state_dir,
    save_run_history_to_state_dir, save_sweep_history_to_state_dir,
};
pub use crate::rbac_io::{load_rbac_config, load_rbac_outcome, RbacLoadOutcome};
pub use crate::settings_io::{load_settings_from_state_dir, save_settings_to_state_dir};
pub use crate::setup::{crashes_dir, run};
pub use crate::solver_parse::{adapter_id_from_solver, derived_inputs_from_case_toml};
pub use crate::state_paths::state_dir;
pub use crate::time_format::format_time_key;
pub use crate::types::{BottomTab, LoadedMesh, LoadedStl, RunHistoryEntry, SweepHistoryEntry};

/// Root application state.
#[derive(Default)]
pub struct ValenxApp {
    /// Opt-in dockable / tiling central-panel layout (View → Docked
    /// layout). Default-built tile tree; only painted when
    /// [`ValenxApp::docked_layout`] is on. See [`docking`].
    pub(crate) docking: docking::DockingState,
    /// When true, the central panel renders the [`docking`] tile tree
    /// (resizable splits / tabs / close / drag) instead of the classic
    /// single-viewport layout. Default `false` (classic layout).
    pub(crate) docked_layout: bool,
    pub(crate) project: Option<LoadedProject>,
    pub(crate) project_path: Option<PathBuf>,
    /// RBAC override block parsed from the loaded project's
    /// `project.toml`. Merged on top of the global `<state_dir>/rbac.json`
    /// at every permission check, so a sensitive project can promote
    /// or demote per-user roles without rewriting the global config.
    /// `None` when no project is loaded or when project.toml has no
    /// `[rbac]` block.
    pub(crate) project_rbac_override: Option<valenx_rbac::RbacConfig>,
    pub(crate) stl: Option<LoadedStl>,
    pub(crate) mesh: Option<LoadedMesh>,
    pub(crate) camera: OrbitCamera,
    pub(crate) shading: ShadingMode,
    pub(crate) last_error: Option<String>,
    pub(crate) status: Option<String>,
    pub(crate) about_open: bool,

    pub(crate) registry: AdapterRegistry,
    pub(crate) residuals: ResidualHistory,
    pub(crate) log: LogPanel,
    pub(crate) bottom_tab: BottomTab,

    /// Which case the user clicked on in the browser, if any. `None`
    /// falls back to the first case in the project when a run is
    /// started.
    pub(crate) selected_case: Option<String>,

    pub(crate) run_handle: Option<RunHandle>,
    /// Live threaded sweep runner. `Some(_)` while a sweep is
    /// executing; cleared when the worker thread finishes / fails.
    pub(crate) sweep_handle: Option<SweepHandle>,
    /// Per-sweep progress: (succeeded, failed, total). Updated as
    /// `SweepEvent::JobFinished` events come in. The numbers persist
    /// across the sweep_handle's lifetime so the status pane can
    /// keep showing the last result after the worker exits.
    pub(crate) sweep_progress: (usize, usize, usize),
    /// Status text for the active sweep — surfaced near the run
    /// progress in the UI.
    pub(crate) sweep_message: String,
    pub(crate) run_progress: f32,
    pub(crate) run_message: String,
    pub(crate) last_run_report: Option<Box<RunReport>>,
    pub(crate) last_run_error: Option<String>,

    /// Last successful prepare-only workdir, if any. Set by
    /// [`Self::prepare_selected_case`] so the UI can show the path
    /// and the "Open in file browser" action can act on it. `None`
    /// until the user clicks "Prepare".
    pub(crate) last_prepare_workdir: Option<PathBuf>,

    /// PreparedJob from the most recent successful prepare, kept so
    /// [`Self::run_from_prepared_workdir`] can run the solver against
    /// the user's hand-edited dicts without re-emitting them. The
    /// adapter id that produced this job lives alongside it because
    /// `spawn_prepared` needs to look the adapter back up in the
    /// registry.
    pub(crate) last_prepared_job: Option<(String, valenx_core::PreparedJob)>,

    /// Last completed run's workdir, captured when the run handle
    /// drops at the end of [`Self::pump_run_events`]. Mirrors
    /// `last_prepare_workdir` for the run pipeline so users can
    /// "Open in file browser" the dir holding their .vtu / .frd /
    /// .log artifacts after the solver finishes. `None` until the
    /// first run completes.
    pub(crate) last_run_workdir: Option<PathBuf>,

    /// Results bundle from the most recent successful run, populated
    /// when the worker thread sends `RunEvent::Collected`. Carries
    /// the parsed Field catalog (e.g. OpenFOAM's VTU fields), scalar
    /// records, artifact list, and provenance. `None` until the
    /// first run completes successfully.
    pub(crate) last_run_results: Option<Box<valenx_fields::Results>>,

    /// Which field the viewport's colour overlay is showing. Set by
    /// clicking a field name in the Results pane. `None` falls back
    /// to "first scalar OnNode field that matches the mesh" — same
    /// auto-pick used before the field selector landed.
    pub(crate) selected_field_name: Option<String>,

    /// Index into the selected field's time series — `0` = first
    /// snapshot, `1` = second, etc. Driven by the slider in the
    /// Results pane. Clamped every frame so the index can't outrun
    /// the time-series length when the user switches fields with
    /// different snapshot counts.
    pub(crate) selected_time_index: usize,

    /// Per-case run history — last outcome + wall time. Keyed by
    /// case name (project-local). Populated when a run finishes;
    /// surfaces in the case browser as a small ✓/✗ badge so users
    /// can see at a glance which cases they've already exercised
    /// without scrolling logs. Persisted to
    /// `<state_dir>/run-history.json` after every run so it
    /// survives app restarts.
    pub(crate) run_history: std::collections::BTreeMap<String, RunHistoryEntry>,
    /// Per-case sweep history. Mirrors `run_history` but for the
    /// sweep pipeline — recorded when a sweep finishes (sync or
    /// async) so the case browser can show "you swept this with N
    /// derived cases at `<ts>`". Persisted to
    /// `<state_dir>/sweep-history.json` so it survives an app
    /// restart.
    pub(crate) sweep_history: std::collections::BTreeMap<String, SweepHistoryEntry>,

    /// Case name currently being run, captured at spawn time so the
    /// Finished/Failed handlers can record the outcome under the
    /// right key even if the user has moved their cursor / changed
    /// `selected_case` while the solver was running.
    pub(crate) running_case_name: Option<String>,

    pub(crate) palette: CommandPalette,
    pub(crate) settings: Settings,
    pub(crate) settings_open: bool,
    pub(crate) theme_applied: bool,

    pub(crate) wgpu_renderer: Option<WgpuRenderer>,

    /// Whether the right-side Mesh Toolbox panel is visible. Defaults
    /// to `true` so it surfaces automatically as soon as a mesh /
    /// STL is loaded; the View menu and the command palette can hide
    /// it for users who want a clean viewport.
    pub(crate) show_mesh_toolbox: bool,
    /// Whether the left-side Browser panel is visible. Defaults to
    /// `true`; the ribbon toggle, View menu, and command palette can
    /// hide it to give the viewport the full width.
    pub(crate) show_browser: bool,
    /// Whether the viewport cursor snaps to the ground grid (Fusion-style):
    /// the live cursor coordinate snaps to the nearest grid node, with a
    /// marker drawn there. Defaults to `true`; toggled from the View menu.
    pub(crate) snap_to_grid: bool,
    /// Receiver for background adapter-probe results (see
    /// [`valenx_core::AdapterRegistry::spawn_probe_all`]). `Some` while the
    /// background probe is in flight; drained each frame in `update` and
    /// cleared to `None` when the probe thread finishes. Probing off the
    /// main thread keeps startup instant — it fixed a ~30s cold-start
    /// freeze (141 external tools probed synchronously in `new`).
    pub(crate) adapter_probe_rx:
        Option<std::sync::mpsc::Receiver<(&'static str, valenx_core::AdapterStatus)>>,
    /// Form-input state for the toolbox panel (translate deltas,
    /// scale factors, rotation axis + angle, mirror plane, cut-
    /// plane point + normal, repair tolerance). Cleared back to
    /// defaults on app construction; persisted across panel toggles.
    pub(crate) mesh_toolbox: crate::mesh_toolbox::MeshToolboxState,

    /// First CAD operand (operand "A" for boolean ops). Set when the
    /// user creates a primitive through the Part workbench section
    /// with the "Create as second" toggle off, and rewritten every
    /// time a boolean op runs (the result replaces operand A).
    pub(crate) current_solid: Option<valenx_cad::Solid>,
    /// Second CAD operand (operand "B"). Set when the user creates
    /// a primitive with the "Create as second" toggle on. Cleared
    /// whenever a boolean op consumes it so the toolbox is honest
    /// about needing a new B for the next op.
    pub(crate) second_solid: Option<valenx_cad::Solid>,

    /// First-launch wizard state. Loaded from
    /// `<state_dir>/first-run.json` on startup; defaults to a
    /// fresh, never-completed decision when the file doesn't exist.
    pub(crate) first_run_decision: valenx_first_run::FirstRunDecision,
    /// Whether the wizard's egui modal is open right now. Always
    /// initialised to `false` — the wizard never auto-opens because
    /// Valenx ships native Rust engines for every major simulation
    /// domain (external adapters are an optional power-user surface,
    /// so pushing first-time users to install OpenFOAM / GROMACS /
    /// Python contradicts the value proposition). Re-openable from
    /// the Settings menu's "Re-probe external tools" entry and the
    /// command palette.
    pub(crate) first_run_open: bool,
    /// Cached environment report. Built lazily on the frame the
    /// wizard opens, so the registry's probe results survive across
    /// frames without re-probing every redraw.
    pub(crate) first_run_report: Option<valenx_first_run::EnvironmentReport>,

    /// Loaded locale catalogue. Populated in `new()` from the
    /// embedded en-US baseline; future versions will pick the
    /// locale matching the user's OS preference and fall back to
    /// en-US when a translation is missing. Wrap in
    /// `Option<Arc<…>>`-style sharing if hot-swap becomes a
    /// requirement (it isn't yet — we set the locale once at
    /// startup).
    pub(crate) catalogue: valenx_i18n::LocaleCatalogue,

    /// Phase 21 — Macro recorder. UI panels append actions via
    /// [`Self::record_macro_action`] when the user clicks a
    /// recordable button. `start_recording` / `stop_recording`
    /// flip the recorder state.
    pub(crate) macro_recorder: valenx_macro::MacroRecorder,

    /// Phase 22 — Add-on registry. Owns the in-memory list of
    /// installed add-ons + dispatches install/update/uninstall via
    /// the manual install-by-directory flow.
    pub(crate) addons: valenx_addons::AddonRegistry,
    /// Whether the Add-on Manager panel is visible.
    pub(crate) show_addon_manager: bool,

    /// Whether the right-side Genetics Workbench panel is
    /// visible. Defaults to `false` (the CAD-side Mesh Toolbox is the
    /// default right panel); flipped on from the View menu / command
    /// palette. The two right-side workbenches are independent — both
    /// can be open at once, egui docks them side by side.
    pub(crate) show_genetics_workbench: bool,
    /// Form + result state for the thirteen Genetics-workbench panels
    /// (one per computational-biology crate). See
    /// [`crate::genetics_workbench`].
    pub(crate) genetics: crate::genetics_workbench::GeneticsWorkbenchState,

    /// Whether the right-side Aerodynamics / Wind
    /// Tunnel workbench panel is visible. Defaults to `false`; flipped
    /// on from the View menu. Independent of the Mesh Toolbox and the
    /// Genetics workbench — egui docks them side by side.
    pub(crate) show_aero_workbench: bool,
    /// Form + result state for the Wind Tunnel workbench — the eight
    /// workflow sections wrapping the `valenx-aero` CFD engine. See
    /// [`crate::aero_workbench`].
    pub(crate) aero: crate::aero_workbench::AeroWorkbenchState,
    /// The aero flow-visualization field overlay, if one is active.
    /// When `Some`, the viewport colours the loaded mesh by this scalar
    /// field through the per-vertex colour ramp — it takes priority
    /// over the post-run results overlay. Pushed by the Wind Tunnel
    /// workbench's "Show field in 3-D viewport" button.
    pub(crate) aero_field_overlay: Option<valenx_fields::Field>,

    /// Whether the right-side FEM Workbench panel is visible. Defaults
    /// to `false`; flipped on from the View menu. Independent of the
    /// other workbenches — egui docks them side by side.
    pub(crate) show_fem_workbench: bool,
    /// Form + result state for the FEM Workbench — native linear-static
    /// and modal finite-element analysis wrapping the `valenx-fem`
    /// in-process solvers (no external solver, no input deck). See
    /// [`crate::fem_workbench`].
    pub(crate) fem: crate::fem_workbench::FemWorkbenchState,

    /// Whether the right-side CFD Workbench panel is visible. Defaults
    /// to `false`; flipped on from the View menu. Independent of the
    /// other workbenches — egui docks them side by side.
    pub(crate) show_cfd_workbench: bool,
    /// Form + result state for the CFD Workbench — native 2-D
    /// incompressible laminar CFD (SIMPLE) wrapping `valenx-cfd-native`.
    /// See [`crate::cfd_workbench`].
    pub(crate) cfd: crate::cfd_workbench::CfdWorkbenchState,

    /// Whether the right-side Reaction Dynamics workbench is visible.
    /// Defaults to `false`; flipped on from the View menu. Independent of
    /// the other workbenches — egui docks them side by side.
    pub(crate) show_reactdyn_workbench: bool,
    /// Form + result state for the Reaction Dynamics workbench — native
    /// ab-initio MD (AIMD) wrapping `valenx-reactdyn`. See
    /// [`crate::reactdyn_workbench`].
    pub(crate) reactdyn: crate::reactdyn_workbench::ReactdynWorkbenchState,

    /// Whether the right-side Springs Workbench is visible. Defaults to
    /// `false`; flipped on from the View menu. Independent of the other
    /// workbenches — egui docks them side by side.
    pub(crate) show_springs_workbench: bool,
    /// Form + result state for the Springs Workbench — native helical-spring
    /// design wrapping `valenx-springs`. See [`crate::springs_workbench`].
    pub(crate) springs: crate::springs_workbench::SpringsWorkbenchState,

    /// Whether the right-side Marine / Hull Workbench is visible. Off by
    /// default; toggled from the View menu.
    pub(crate) show_marine_workbench: bool,
    /// Form + result state for the Marine / Hull Workbench — native
    /// box-form hull hydrostatics wrapping `valenx-marine`. See
    /// [`crate::marine_workbench`].
    pub(crate) marine: crate::marine_workbench::MarineWorkbenchState,

    /// Whether the right-side Gears Workbench is visible. Defaults to
    /// `false`; flipped on from the View menu. Independent of the other
    /// workbenches — egui docks them side by side.
    pub(crate) show_gears_workbench: bool,
    /// Form + result state for the Gears Workbench — native involute-gear
    /// design wrapping `valenx-gears`. See [`crate::gears_workbench`].
    pub(crate) gears: crate::gears_workbench::GearsWorkbenchState,

    /// Whether the right-side Drone Workbench is visible. Defaults to
    /// `false`; flipped on from the View menu. Independent of the other
    /// workbenches — egui docks them side by side.
    pub(crate) show_drone_workbench: bool,
    /// Form + result state for the Drone Workbench — native multirotor
    /// hover performance wrapping `valenx-drone`. See [`crate::drone_workbench`].
    pub(crate) drone: crate::drone_workbench::DroneWorkbenchState,

    /// Whether the right-side Geomatics Workbench is visible. Defaults to
    /// `false`; flipped on from the View menu. Independent of the other
    /// workbenches — egui docks them side by side.
    pub(crate) show_geomatics_workbench: bool,
    /// Form + result state for the Geomatics Workbench — native geodesic
    /// calculations wrapping `valenx-geomatics`. See
    /// [`crate::geomatics_workbench`].
    pub(crate) geomatics: crate::geomatics_workbench::GeomaticsWorkbenchState,

    /// Whether the right-side Four-Bar Linkage Workbench is visible.
    /// Defaults to `false`; flipped on from the View menu. Independent of the
    /// other workbenches — egui docks them side by side.
    pub(crate) show_fourbar_workbench: bool,
    /// Form + result state for the Four-Bar Linkage Workbench — native planar
    /// four-bar mechanism kinematics wrapping `valenx-kinematics`. See
    /// [`crate::fourbar_workbench`].
    pub(crate) fourbar: crate::fourbar_workbench::FourBarWorkbenchState,

    /// Whether the right-side Piping Workbench is visible. Defaults to
    /// `false`; flipped on from the View menu. Independent of the other
    /// workbenches — egui docks them side by side.
    pub(crate) show_piping_workbench: bool,
    /// Form + result state for the Piping Workbench — native pipe-section
    /// sizing wrapping `valenx-piping`. See [`crate::piping_workbench`].
    pub(crate) piping: crate::piping_workbench::PipingWorkbenchState,

    /// Whether the right-side Rail / Train Workbench is visible. Defaults to
    /// `false`; flipped on from the View menu. Independent of the other
    /// workbenches — egui docks them side by side.
    pub(crate) show_rail_workbench: bool,
    /// Form + result state for the Rail / Train Workbench — native train
    /// resistance + tractive effort wrapping `valenx-rail`. See
    /// [`crate::rail_workbench`].
    pub(crate) rail: crate::rail_workbench::RailWorkbenchState,

    /// Whether the right-side Collision Workbench is visible. Defaults to
    /// `false`; flipped on from the View menu. Independent of the other
    /// workbenches — egui docks them side by side.
    pub(crate) show_collision_workbench: bool,
    /// Form + result state for the Collision Workbench — native AABB
    /// geometry + overlap tests wrapping `valenx-collision`. See
    /// [`crate::collision_workbench`].
    pub(crate) collision: crate::collision_workbench::CollisionWorkbenchState,

    /// Whether the right-side Solar PV Workbench is visible. Defaults to
    /// `false`; flipped on from the View menu. Independent of the other
    /// workbenches — egui docks them side by side.
    pub(crate) show_solarpv_workbench: bool,
    /// Form + result state for the Solar PV Workbench — native single-diode
    /// photovoltaic cell performance wrapping `valenx-solarpv`. See
    /// [`crate::solarpv_workbench`].
    pub(crate) solarpv: crate::solarpv_workbench::SolarPvWorkbenchState,

    /// Whether the right-side Sheet Metal Workbench is visible. Defaults to
    /// `false`; flipped on from the View menu. Independent of the other
    /// workbenches — egui docks them side by side.
    pub(crate) show_sheetmetal_workbench: bool,
    /// Form + result state for the Sheet Metal Workbench — native bend
    /// allowance / deduction wrapping `valenx-sheet-metal`. See
    /// [`crate::sheetmetal_workbench`].
    pub(crate) sheetmetal: crate::sheetmetal_workbench::SheetmetalWorkbenchState,

    /// Whether the right-side Field Statistics Workbench is visible. Defaults
    /// to `false`; flipped on from the View menu. Independent of the other
    /// workbenches — egui docks them side by side.
    pub(crate) show_fields_workbench: bool,
    /// Form + result state for the Field Statistics Workbench — descriptive
    /// statistics over a pasted number list, via `valenx-fields`. See
    /// [`crate::fields_workbench`].
    pub(crate) fields: crate::fields_workbench::FieldsWorkbenchState,

    /// Whether the right-side Fasteners Workbench is visible. Defaults to
    /// `false`; flipped on from the View menu. Independent of the other
    /// workbenches — egui docks them side by side.
    pub(crate) show_fasteners_workbench: bool,
    /// Form + result state for the Fasteners Workbench — ISO 4017 hex-bolt
    /// dimensions wrapping `valenx-fasteners`. See
    /// [`crate::fasteners_workbench`].
    pub(crate) fasteners: crate::fasteners_workbench::FastenersWorkbenchState,

    /// Whether the right-side Fixed-Wing / Aircraft Workbench is visible.
    /// Defaults to `false`; flipped on from the View menu. Independent of the
    /// other workbenches — egui docks them side by side.
    pub(crate) show_fixedwing_workbench: bool,
    /// Form + result state for the Fixed-Wing / Aircraft Workbench — native
    /// preliminary aircraft point-performance wrapping `valenx-fixedwing`.
    /// See [`crate::fixedwing_workbench`].
    pub(crate) fixedwing: crate::fixedwing_workbench::FixedWingWorkbenchState,

    /// Whether the right-side Frames Workbench is visible. Defaults to
    /// `false`; flipped on from the View menu. Independent of the other
    /// workbenches — egui docks them side by side.
    pub(crate) show_frames_workbench: bool,
    /// Form + result state for the Frames Workbench — structural
    /// cross-section properties wrapping `valenx-frames`. See
    /// [`crate::frames_workbench`].
    pub(crate) frames: crate::frames_workbench::FramesWorkbenchState,

    /// Whether the right-side Gas Dynamics workbench is visible. Defaults to
    /// `false`; flipped on from the View menu. Independent of the other
    /// workbenches — egui docks them side by side.
    pub(crate) show_gasdynamics_workbench: bool,
    /// Form + result state for the Gas Dynamics workbench — 1-D
    /// compressible-flow relations wrapping `valenx-gasdynamics`. See
    /// [`crate::gasdynamics_workbench`].
    pub(crate) gasdynamics: crate::gasdynamics_workbench::GasDynamicsWorkbenchState,

    /// Whether the right-side Neural-Interface (BCI stimulation) workbench is
    /// visible. Defaults to `false`; flipped on from the View menu.
    pub(crate) show_neuro_workbench: bool,
    /// Form + result state for the Neural-Interface workbench, wrapping
    /// `valenx-neuro`. See [`crate::neuro_workbench`].
    pub(crate) neuro: crate::neuro_workbench::NeuroWorkbenchState,

    /// Whether the right-side Wind Turbine workbench is visible. Defaults to
    /// `false`; flipped on from the View menu. Independent of the other
    /// workbenches — egui docks them side by side.
    pub(crate) show_windturbine_workbench: bool,
    /// Form + result state for the Wind Turbine workbench — native
    /// actuator-disc wind-turbine power wrapping `valenx-windturbine`. See
    /// [`crate::windturbine_workbench`].
    pub(crate) windturbine: crate::windturbine_workbench::WindTurbineWorkbenchState,

    /// Whether the right-side Parametric-CAD workbench is visible. Defaults to
    /// `false`; flipped on from the View menu.
    pub(crate) show_cad_workbench: bool,
    /// Form + result state for the Parametric-CAD workbench, wrapping
    /// `valenx-solvespace-3d`. See [`crate::cad_workbench`].
    pub(crate) cad: crate::cad_workbench::CadWorkbenchState,

    /// Whether the right-side 2D Drafting workbench is visible. Defaults to
    /// `false`; flipped on from the View menu.
    pub(crate) show_draft2d_workbench: bool,
    /// State for the 2D Drafting workbench, wrapping `valenx-librecad-2d`. See
    /// [`crate::draft2d_workbench`].
    pub(crate) draft2d: crate::draft2d_workbench::Draft2dWorkbenchState,

    /// Whether the right-side Reinforcement workbench is visible. Defaults to
    /// `false`; flipped on from the View menu.
    pub(crate) show_reinforcement_workbench: bool,
    /// State for the Reinforcement workbench, wrapping `valenx-reinforcement`.
    /// See [`crate::reinforcement_workbench`].
    pub(crate) reinforcement: crate::reinforcement_workbench::ReinforcementWorkbenchState,

    /// Whether the right-side Path-Traced Render workbench is visible. Defaults
    /// to `false`; flipped on from the View menu.
    pub(crate) show_render_workbench: bool,
    /// State for the Render workbench, wrapping `valenx-pathtrace`. See
    /// [`crate::render_workbench`].
    pub(crate) render: crate::render_workbench::RenderWorkbenchState,

    /// Whether the right-side HVAC workbench is visible. Defaults to `false`;
    /// flipped on from the View menu.
    pub(crate) show_hvac_workbench: bool,
    /// State for the HVAC workbench, wrapping `valenx-hvac`. See
    /// [`crate::hvac_workbench`].
    pub(crate) hvac: crate::hvac_workbench::HvacWorkbenchState,

    /// Whether the right-side Reverse-Engineering workbench is visible.
    /// Defaults to `false`; flipped on from the View menu.
    pub(crate) show_reverse_workbench: bool,
    /// State for the Reverse-Engineering workbench, wrapping `valenx-reverse`.
    /// See [`crate::reverse_workbench`].
    pub(crate) reverse: crate::reverse_workbench::ReverseWorkbenchState,

    /// Whether the right-side Interior-Design workbench is visible. Defaults to
    /// `false`; flipped on from the View menu.
    pub(crate) show_interior_workbench: bool,
    /// State for the Interior-Design workbench, wrapping `valenx-interior`. See
    /// [`crate::interior_workbench`].
    pub(crate) interior: crate::interior_workbench::InteriorWorkbenchState,

    /// Whether the right-side Animation workbench is visible. Defaults to
    /// `false`; flipped on from the View menu.
    pub(crate) show_animate_workbench: bool,
    /// State for the Animation workbench, wrapping `valenx-animate`. See
    /// [`crate::animate_workbench`].
    pub(crate) animate: crate::animate_workbench::AnimateWorkbenchState,

    /// Whether the right-side Variant-Effect workbench is visible. Defaults to
    /// `false`; flipped on from the View menu.
    pub(crate) show_variant_effect_workbench: bool,
    /// State for the Variant-Effect workbench, wrapping `valenx-variant-effect`.
    /// See [`crate::variant_effect_workbench`].
    pub(crate) variant_effect: crate::variant_effect_workbench::VariantEffectWorkbenchState,

    /// Whether the right-side Astro / Launch workbench panel is visible.
    /// Defaults to `false`; flipped on from the View menu (Ctrl+4).
    /// Independent of the Mesh Toolbox, Genetics and Wind Tunnel
    /// workbenches — egui docks them side by side.
    pub(crate) show_astro_workbench: bool,
    /// Form + result state for the Astro / Launch workbench — the launch
    /// ascent simulator + the closed-form mission planners wrapping the
    /// `valenx-astro` crate. See [`crate::astro_workbench`].
    pub(crate) astro: crate::astro_workbench::AstroWorkbenchState,

    /// Whether the right-side Rocket workbench panel is visible. Defaults
    /// to `false`; flipped on from the View menu. Surfaces the
    /// `valenx-rocket-demo` coupled design→simulate pipeline.
    pub(crate) show_rocket_workbench: bool,
    /// Form + result state for the Rocket workbench — the reactive
    /// design→simulate panel wrapping `valenx-rocket-demo`. See
    /// [`crate::rocket_workbench`].
    pub(crate) rocket: crate::rocket_workbench::RocketWorkbenchState,

    /// Whether the right-side Engine workbench panel is visible — the
    /// reactive engine design → analyze → optimize → export loop. On by
    /// default (set in [`ValenxApp::new`]).
    pub(crate) show_engine_workbench: bool,
    /// Form + result state for the Engine workbench. See
    /// [`crate::engine_workbench`].
    pub(crate) engine: crate::engine_workbench::EngineWorkbenchState,

    /// Whether the right-side Heat Exchanger workbench is visible. Defaults to
    /// `false`; flipped on from the View menu.
    pub(crate) show_heatexchanger_workbench: bool,
    /// State for the Heat Exchanger workbench, wrapping `valenx-heatexchanger`.
    /// See [`crate::heatexchanger_workbench`].
    pub(crate) heatexchanger: crate::heatexchanger_workbench::HeatExchangerWorkbenchState,

    /// Whether the right-side Car workbench panel is visible. Defaults to
    /// `false`; toggled from the View menu. Wraps `valenx-vehicle`'s
    /// performance model. See [`crate::car_workbench`].
    pub(crate) show_car_workbench: bool,
    /// Form + result state for the Car workbench (design → simulate over
    /// `valenx-vehicle`).
    pub(crate) car: crate::car_workbench::CarWorkbenchState,

    /// Whether the right-side Assistant activity sidebar is visible. On by
    /// default (set in [`ValenxApp::new`]) so the app narrates its own work
    /// via the live feed.
    pub(crate) show_assistant_panel: bool,
    /// State for the Assistant activity sidebar (the live `.jsonl` feed
    /// path). See [`crate::assistant_workbench`].
    pub(crate) assistant: crate::assistant_workbench::AssistantWorkbenchState,

    /// Whether the keyboard-shortcut cheat-sheet overlay is open.
    /// Toggled by the `?` key + by Help → Keyboard shortcuts.
    pub(crate) keyboard_help_open: bool,

    /// Whether the per-panel contextual help popup is open. Mapped
    /// to F1 + Help → Panel help.
    pub(crate) panel_help_open: bool,

    /// Which panel's help text the F1 popup shows. Resolved at the
    /// moment F1 is pressed from "what workbench is active right
    /// now"; defaults to "Sequence" when nothing else is up.
    pub(crate) panel_help_target: String,

    /// Whether the first-launch welcome tour is currently open. Auto-set
    /// to `true` on a fresh install (gated by `settings.welcome_tour_completed`);
    /// re-openable from the Help menu.
    pub(crate) welcome_tour_open: bool,

    /// Tour navigation state — which step the user is on, and
    /// whether they've finished. See [`crate::welcome_tour::TourState`].
    pub(crate) welcome_tour_state: crate::welcome_tour::TourState,

    /// "File → New Project…" modal state. `Some(_)` while the dialog
    /// is open; `None` once the user clicks Create / Cancel / closes
    /// the window. Triggered by the File menu, the command palette,
    /// and the Ctrl+N shortcut. See [`crate::new_project_dialog`].
    pub(crate) new_project_dialog: Option<crate::new_project_dialog::NewProjectDialog>,

    /// One-line notice rendered inline on the welcome / landing page
    /// (next to the recent-projects list). Set by the host when a
    /// landing-page action produces a result that doesn't belong in
    /// the top status bar — currently the "removed missing project
    /// from recents" confirmation. Cleared as soon as the user takes
    /// another action on the landing page.
    pub(crate) landing_inline_message: Option<String>,

    /// Memoised command-palette entry list, keyed by
    /// `(registry.len(), show_non_oss_adapters)`. `build_visible_commands`
    /// allocates ~360 `String`s per call and used to run every frame;
    /// the cache invalidates only when the registry grows (rare —
    /// re-probe / load) or the OSS-only toggle flips in Settings.
    /// `None` until the first palette render fills it.
    pub(crate) palette_cache: Option<(usize, bool, Vec<crate::commands::CommandKind>)>,

    // ── Swappable viewport system (cloud/viewport) ────────────────────────
    /// Which viewport implementation is rendered in the central panel.
    ///
    /// Defaults to `Viewport3D`; switches to `Viewport2dDna` when the
    /// user first enables the Genetics Workbench (and can be overridden
    /// at any time from **View → Central viewport**).
    pub(crate) active_viewport: crate::viewport_kind::ViewportKind,

    /// Open project tabs (Chrome-style) plus the active index. Drives
    /// which workbench the tab strip shows. See [`crate::project_tabs`].
    pub(crate) tab_bar: crate::project_tabs::TabBar,

    /// Persistent state for the 2D DNA / plasmid viewport. Survives
    /// viewport-kind switches so pan, zoom, and sub-view selection are
    /// remembered when the user returns to the 2D view.
    pub(crate) viewport_2d: crate::viewport_2d::Viewport2dState,
}

impl ValenxApp {
    /// Build a fresh `ValenxApp`. Restores persisted user settings
    /// from the OS state directory and, when `initial_stl` is `Some`,
    /// queues that file for loading on the next frame.
    pub fn new(cc: &eframe::CreationContext<'_>, initial_stl: Option<PathBuf>) -> Self {
        let mut app = Self::default();
        // Restore persisted user preferences from disk before
        // applying any defaults — the user's last theme / shading
        // choice is what they expect to see when the app reopens.
        if let Some(loaded) = load_settings_from_state_dir() {
            app.settings = loaded;
        }
        // Bridge the `force_external_vina` toggle to the Vina adapter
        // via the process-wide atomic in valenx-core. The adapter
        // reads it in `run()` to decide whether to skip the native
        // engine even when the case picked `engine = "native"`.
        // Atomic rather than env var because `std::env::set_var` is
        // `unsafe` to call once threads exist (Linux glibc) — and
        // egui has already spawned its render thread by here.
        valenx_core::set_force_external_vina(app.settings.force_external_vina);
        app.shading = app.settings.default_shading;
        // Mesh Toolbox is on by default — surfaces automatically as
        // soon as the user drops an STL or loads a canonical mesh.
        // Hidden via View menu / palette for a clean viewport.
        app.show_mesh_toolbox = true;
        app.show_browser = true;
        // Assistant activity sidebar on by default — the desktop app
        // narrates its own work via a live feed (empty until appended to).
        app.show_assistant_panel = true;
        // Rocket workbench on by default too, so the Valenx LV-1 ascent
        // plot is visible at launch without hunting the View menu.
        app.show_rocket_workbench = true;
        // Engine-design workbench on by default too — design → analyze →
        // optimize → export an engine, visible at launch.
        app.show_engine_workbench = true;
        app.snap_to_grid = true;
        app.init_registry();
        // Restore the per-case run-history map from disk so the
        // browser's ✓/✗ badges survive app restarts. A missing /
        // unparseable file is silently treated as "no history yet";
        // we'll start fresh and over-write on the next run.
        if let Some(history) = load_run_history_from_state_dir() {
            app.run_history = history;
        }
        if let Some(history) = load_sweep_history_from_state_dir() {
            app.sweep_history = history;
        }
        // First-launch wizard. Load any persisted decision so a user
        // who already dismissed it stays dismissed. Auto-show used to
        // gate on `should_auto_show` (open if not yet completed) — that
        // pushed users to install OpenFOAM / GROMACS / Python toolchains
        // on first launch, which contradicts Valenx's value prop (we
        // ship native Rust engines for every major simulation domain;
        // external adapters are an optional power-user surface). The
        // wizard now stays accessible from the Settings menu's
        // "Re-probe external tools" entry + the command palette, but
        // doesn't auto-open on first launch.
        if let Some(decision) = first_run::load_first_run_from_state_dir() {
            app.first_run_decision = decision;
        }
        app.first_run_open = false;

        // Welcome tour — auto-open on the next launch after install
        // (or after a settings reset) until the user dismisses it.
        // Independent from the wizard above: the wizard probes the
        // adapter environment, the tour orients the new user in the
        // four workbenches.
        if !app.settings.welcome_tour_completed {
            app.welcome_tour_open = true;
        }
        // Restore the cheat-sheet overlay's open state — a user who
        // pinned it from the palette had it remember the choice.
        app.keyboard_help_open = app.settings.keyboard_shortcuts_overlay_open;
        // Default the F1-help target so first-frame F1 lands on
        // something — the value gets recomputed on every press of F1
        // from the active workbench.
        app.panel_help_target = "Sequence".to_string();

        // Load the en-US locale catalogue. The baseline is baked
        // into the binary via `include_str!` at compile time so we
        // don't need to ship the .ftl file alongside. Future
        // versions (v0.2.0+) will let the user pick a locale from
        // Settings; for v0.1.0 every install renders en-US.
        app.catalogue = valenx_i18n::embedded_en_us();
        if let Some(render_state) = cc.wgpu_render_state.as_ref() {
            app.wgpu_renderer = Some(WgpuRenderer::new(render_state));
        }
        // Surface RBAC parse errors at startup. A malformed
        // `<state_dir>/rbac.json` would otherwise silently fall back
        // to a default config; we instead surface the parse error so
        // the operator notices, and the loader fails closed to the
        // Viewer role (see `RbacLoadOutcome::into_active_config`).
        if let RbacLoadOutcome::ParseError(msg) = load_rbac_outcome() {
            app.last_error = Some(format!(
                "RBAC config failed to parse — running with read-only Viewer role until fixed: {msg}"
            ));
        }
        if let Some(path) = initial_stl {
            app.load_stl(path);
        }
        app
    }

    /// Build a fresh [`valenx_first_run::EnvironmentReport`] from
    /// the registry's current probe results. Caches the report on
    /// `self.first_run_report` so subsequent frames skip the
    /// rebuild — call [`Self::invalidate_first_run_report`] after
    /// a re-probe to refresh.
    ///
    /// Honours `settings.show_non_oss_adapters`: when `false` (the
    /// default), adapters whose `tool_license` doesn't satisfy
    /// `is_oss_license()` are omitted from the report, so the wizard's
    /// "Detected N of M" headline + per-row table reflect the
    /// OSS-only filtered totals. The setting toggle invalidates the
    /// cached report so a flip in Settings reflects on the next open.
    pub(crate) fn ensure_first_run_report(&mut self) -> &valenx_first_run::EnvironmentReport {
        if self.first_run_report.is_none() {
            let show_non_oss = self.settings.show_non_oss_adapters;
            let probes: Vec<(
                String,
                String,
                valenx_first_run::AdapterAvailability,
                Option<String>,
            )> = self
                .registry
                .iter()
                .filter(|(_, entry)| {
                    show_non_oss
                        || valenx_core::adapter_helpers::is_oss_license(
                            entry.adapter.info().tool_license,
                        )
                })
                .map(|(_, entry)| {
                    let info = entry.adapter.info();
                    let availability = match &entry.status {
                        AdapterStatus::Ready { .. } => {
                            valenx_first_run::AdapterAvailability::Installed
                        }
                        AdapterStatus::Missing { .. } => {
                            valenx_first_run::AdapterAvailability::Missing
                        }
                        AdapterStatus::Outdated { .. } => {
                            valenx_first_run::AdapterAvailability::Outdated
                        }
                        AdapterStatus::Broken { .. } => {
                            valenx_first_run::AdapterAvailability::Broken
                        }
                        AdapterStatus::Disabled => valenx_first_run::AdapterAvailability::Disabled,
                    };
                    let detected_version = match &entry.status {
                        AdapterStatus::Ready { report } => {
                            report.found_version.as_ref().map(|v| v.to_string())
                        }
                        AdapterStatus::Outdated { found, .. } => Some(found.clone()),
                        _ => None,
                    };
                    (
                        info.id.to_string(),
                        info.display_name.to_string(),
                        availability,
                        detected_version,
                    )
                })
                .collect();
            self.first_run_report = Some(first_run::build_report(&probes));
        }
        self.first_run_report.as_ref().expect("just populated")
    }

    /// Drop the cached first-run report so the next render rebuilds
    /// it from the current registry status. Called after a re-probe
    /// — without this the wizard would show stale "Missing" rows.
    pub fn invalidate_first_run_report(&mut self) {
        self.first_run_report = None;
    }

    /// Rerun probes on every registered adapter. Called from the
    /// Re-probe button in the browser and from Settings when the
    /// "re-probe on close" switch is on.
    pub fn reprobe(&mut self) {
        let ready = self.registry.probe_all();
        self.status = Some(format!(
            "Adapter probe complete · {ready} of {} ready",
            self.registry.len()
        ));
    }

    /// Walk every registered adapter and return a JSON-friendly catalog.
    /// LLM clients call this FIRST to learn what the tool can do.
    pub fn list_capabilities(&self) -> serde_json::Value {
        let adapters: Vec<valenx_core::AdapterDescriptor> = self
            .registry
            .iter()
            .map(|(_, entry)| valenx_core::AdapterDescriptor::from_info(&entry.adapter.info()))
            .collect();
        serde_json::json!({
            "valenx_version": env!("CARGO_PKG_VERSION"),
            "adapter_count": adapters.len(),
            "adapters": adapters,
        })
    }

    /// Pop up the Settings window.
    pub fn open_settings(&mut self) {
        self.settings_open = true;
    }

    /// Pop up the About dialog.
    pub fn open_about(&mut self) {
        self.about_open = true;
    }

    /// Open the "File → New Project…" modal. Idempotent — re-opening
    /// while already open is a no-op (keeps the user's in-progress
    /// edits intact).
    pub fn open_new_project_dialog(&mut self) {
        if self.new_project_dialog.is_none() {
            self.new_project_dialog = Some(crate::new_project_dialog::NewProjectDialog::default());
        }
    }

    /// Make the right-side Mesh Toolbox SidePanel visible. Mirrors the
    /// `toggle_mesh_toolbox` action in semantics ("show this panel"),
    /// but with a non-toggling signature so callers (notably the
    /// headless screenshot harness) get a deterministic post-condition
    /// without needing to read the prior state.
    pub fn enable_mesh_toolbox(&mut self) {
        self.show_mesh_toolbox = true;
    }

    /// Make the right-side Genetics Workbench SidePanel visible and
    /// select `panel` as its active tab. The screenshot harness drives
    /// this once per genetics panel — one PNG per active selection.
    pub fn enable_genetics_workbench(&mut self, panel: crate::genetics_workbench::GeneticsPanel) {
        self.show_genetics_workbench = true;
        self.genetics.active = panel;
    }

    /// Make the right-side Aerodynamics / Wind Tunnel workbench
    /// SidePanel visible. The eight workflow sections all render
    /// inside one scrollable column when the panel is shown.
    pub fn enable_aero_workbench(&mut self) {
        self.show_aero_workbench = true;
    }

    /// Make the right-side FEM Workbench SidePanel visible. The
    /// linear-static + modal analyses render inside one scrollable
    /// column when the panel is shown.
    pub fn enable_fem_workbench(&mut self) {
        self.show_fem_workbench = true;
    }

    /// Make the right-side CFD Workbench SidePanel visible. The
    /// lid-driven-cavity + channel-flow cases render inside one
    /// scrollable column when the panel is shown.
    pub fn enable_cfd_workbench(&mut self) {
        self.show_cfd_workbench = true;
    }

    /// Make the right-side Reaction Dynamics workbench SidePanel visible.
    pub fn enable_reactdyn_workbench(&mut self) {
        self.show_reactdyn_workbench = true;
    }

    /// Make the right-side Springs Workbench SidePanel visible.
    pub fn enable_springs_workbench(&mut self) {
        self.show_springs_workbench = true;
    }

    /// Make the right-side Gears Workbench SidePanel visible.
    pub fn enable_gears_workbench(&mut self) {
        self.show_gears_workbench = true;
    }

    /// Make the right-side Geomatics Workbench SidePanel visible.
    pub fn enable_geomatics_workbench(&mut self) {
        self.show_geomatics_workbench = true;
    }

    /// Make the right-side Piping Workbench SidePanel visible.
    pub fn enable_piping_workbench(&mut self) {
        self.show_piping_workbench = true;
    }

    /// Make the right-side Collision Workbench SidePanel visible.
    pub fn enable_collision_workbench(&mut self) {
        self.show_collision_workbench = true;
    }

    /// Make the right-side Sheet Metal Workbench SidePanel visible.
    pub fn enable_sheetmetal_workbench(&mut self) {
        self.show_sheetmetal_workbench = true;
    }

    /// Make the right-side Field Statistics Workbench SidePanel visible.
    pub fn enable_fields_workbench(&mut self) {
        self.show_fields_workbench = true;
    }

    /// Make the right-side Fasteners Workbench SidePanel visible.
    pub fn enable_fasteners_workbench(&mut self) {
        self.show_fasteners_workbench = true;
    }

    /// Make the right-side Frames Workbench SidePanel visible.
    pub fn enable_frames_workbench(&mut self) {
        self.show_frames_workbench = true;
    }

    /// Make the right-side Astro / Launch workbench SidePanel visible.
    /// The ascent simulator + the mission planners render inside one
    /// scrollable column (across two tabs) when the panel is shown.
    pub fn enable_astro_workbench(&mut self) {
        self.show_astro_workbench = true;
    }

    /// Register every compiled-in adapter and probe them all.
    /// Tomorrow this list comes from the plugin loader; for now it's
    /// the hard-coded set of adapters that have graduated from
    /// scaffold to real implementations.
    fn init_registry(&mut self) {
        // CFD
        self.registry
            .register(Arc::new(valenx_adapter_openfoam::OpenFoamAdapter::new()));
        self.registry
            .register(Arc::new(valenx_adapter_su2::Su2Adapter::new()));
        // Mesh / CAD
        self.registry
            .register(Arc::new(valenx_adapter_gmsh::GmshAdapter::new()));
        self.registry
            .register(Arc::new(valenx_adapter_netgen::NetgenAdapter::new()));
        self.registry
            .register(Arc::new(valenx_adapter_freecad::FreeCadAdapter::new()));
        self.registry
            .register(Arc::new(valenx_adapter_occt::OcctAdapter::new()));
        // FEA
        self.registry
            .register(Arc::new(valenx_adapter_calculix::CalculixAdapter::new()));
        self.registry
            .register(Arc::new(valenx_adapter_elmer::ElmerAdapter::new()));
        self.registry
            .register(Arc::new(valenx_adapter_code_aster::CodeAsterAdapter::new()));
        self.registry.register(Arc::new(
            valenx_adapter_openradioss::OpenRadiossAdapter::new(),
        ));
        // Chemistry
        self.registry
            .register(Arc::new(valenx_adapter_cantera::CanteraAdapter::new()));
        // MD
        self.registry
            .register(Arc::new(valenx_adapter_lammps::LammpsAdapter::new()));
        self.registry
            .register(Arc::new(valenx_adapter_gromacs::GromacsAdapter::new()));
        // EM
        self.registry
            .register(Arc::new(valenx_adapter_openems::OpenEmsAdapter::new()));
        self.registry
            .register(Arc::new(valenx_adapter_meep::MeepAdapter::new()));
        // Battery
        self.registry
            .register(Arc::new(valenx_adapter_pybamm::PyBammAdapter::new()));
        // Dynamics / coupling
        self.registry
            .register(Arc::new(valenx_adapter_mujoco::MuJoCoAdapter::new()));
        self.registry
            .register(Arc::new(valenx_adapter_precice::PreciceAdapter::new()));
        // Biology (Phase 17)
        self.registry
            .register(Arc::new(valenx_adapter_biopython::BiopythonAdapter::new()));
        self.registry
            .register(Arc::new(valenx_adapter_rdkit::RdkitAdapter::new()));
        self.registry
            .register(Arc::new(valenx_adapter_openmm::OpenMmAdapter::new()));
        self.registry
            .register(Arc::new(valenx_adapter_chimerax::ChimeraXAdapter::new()));
        self.registry
            .register(Arc::new(valenx_adapter_oxdna::OxDnaAdapter::new()));
        self.registry
            .register(Arc::new(valenx_adapter_mdanalysis::MdAnalysisAdapter::new()));
        self.registry
            .register(Arc::new(valenx_adapter_colabfold::ColabFoldAdapter::new()));
        // Biology (Phase 18) — sequence alignment toolkit
        self.registry
            .register(Arc::new(valenx_adapter_bwa::BwaAdapter::new()));
        self.registry
            .register(Arc::new(valenx_adapter_minimap2::Minimap2Adapter::new()));
        self.registry
            .register(Arc::new(valenx_adapter_mafft::MafftAdapter::new()));
        self.registry
            .register(Arc::new(valenx_adapter_muscle::MuscleAdapter::new()));
        self.registry
            .register(Arc::new(valenx_adapter_hmmer::HmmerAdapter::new()));
        self.registry
            .register(Arc::new(valenx_adapter_samtools::SamtoolsAdapter::new()));
        // Biology (Phase 17.5) — structure prediction expansion
        self.registry
            .register(Arc::new(valenx_adapter_esmfold::EsmFoldAdapter::new()));
        self.registry
            .register(Arc::new(valenx_adapter_openfold::OpenFoldAdapter::new()));
        self.registry
            .register(Arc::new(valenx_adapter_alphafold2::AlphaFold2Adapter::new()));
        self.registry
            .register(Arc::new(valenx_adapter_alphafold3::AlphaFold3Adapter::new()));
        // Biology (Phase 19) — variant calling toolkit
        self.registry
            .register(Arc::new(valenx_adapter_bcftools::BcftoolsAdapter::new()));
        self.registry
            .register(Arc::new(valenx_adapter_gatk::GatkAdapter::new()));
        self.registry.register(Arc::new(
            valenx_adapter_deepvariant::DeepVariantAdapter::new(),
        ));
        // Biology (Phase 23) — molecular viewers
        self.registry
            .register(Arc::new(valenx_adapter_pymol::PymolAdapter::new()));
        self.registry
            .register(Arc::new(valenx_adapter_vmd::VmdAdapter::new()));
        self.registry
            .register(Arc::new(valenx_adapter_igv::IgvAdapter::new()));
        // Biology (Phase 27) — protein design
        self.registry.register(Arc::new(
            valenx_adapter_rfdiffusion::RfDiffusionAdapter::new(),
        ));
        self.registry.register(Arc::new(
            valenx_adapter_proteinmpnn::ProteinMpnnAdapter::new(),
        ));
        // Biology (Phase 34) — molecular docking
        self.registry
            .register(Arc::new(valenx_adapter_vina::VinaAdapter::new()));
        self.registry
            .register(Arc::new(valenx_adapter_autodock4::AutoDock4Adapter::new()));
        // Biology (Phase 24) — cheminformatics expansion
        self.registry
            .register(Arc::new(valenx_adapter_deepchem::DeepChemAdapter::new()));
        self.registry
            .register(Arc::new(valenx_adapter_openbabel::OpenBabelAdapter::new()));
        self.registry
            .register(Arc::new(valenx_adapter_avogadro::AvogadroAdapter::new()));
        // Biology (Phase 22) — workflow managers
        self.registry
            .register(Arc::new(valenx_adapter_nextflow::NextflowAdapter::new()));
        self.registry
            .register(Arc::new(valenx_adapter_snakemake::SnakemakeAdapter::new()));
        // Biology (Phase 19.5) — single-cell genomics
        self.registry
            .register(Arc::new(valenx_adapter_scanpy::ScanpyAdapter::new()));
        self.registry
            .register(Arc::new(valenx_adapter_scvi::ScviAdapter::new()));
        // Biology (Phase 27.5) — protein design expansion
        self.registry
            .register(Arc::new(valenx_adapter_chroma::ChromaAdapter::new()));
        self.registry
            .register(Arc::new(valenx_adapter_esm_if::EsmIfAdapter::new()));
        self.registry
            .register(Arc::new(valenx_adapter_rfantibody::RfAntibodyAdapter::new()));
        // Biology (Phase 18.5) — aligners expansion
        self.registry
            .register(Arc::new(valenx_adapter_bowtie2::Bowtie2Adapter::new()));
        self.registry
            .register(Arc::new(valenx_adapter_mmseqs2::Mmseqs2Adapter::new()));
        self.registry
            .register(Arc::new(valenx_adapter_diamond::DiamondAdapter::new()));
        // Biology (Phase 18.6) — RNA-seq alignment
        self.registry
            .register(Arc::new(valenx_adapter_hisat2::Hisat2Adapter::new()));
        self.registry
            .register(Arc::new(valenx_adapter_star::StarAdapter::new()));
        // Biology (Phase 20) — transcript quantification
        self.registry
            .register(Arc::new(valenx_adapter_salmon::SalmonAdapter::new()));
        self.registry
            .register(Arc::new(valenx_adapter_kallisto::KallistoAdapter::new()));
        // Biology (Phase 30) — phylogenetics
        self.registry
            .register(Arc::new(valenx_adapter_iqtree::IqTreeAdapter::new()));
        self.registry
            .register(Arc::new(valenx_adapter_raxml_ng::RaxmlNgAdapter::new()));
        self.registry
            .register(Arc::new(valenx_adapter_fasttree::FastTreeAdapter::new()));
        // Biology (Phase 28) — RNA secondary structure
        self.registry
            .register(Arc::new(valenx_adapter_viennarna::ViennaRnaAdapter::new()));
        self.registry.register(Arc::new(
            valenx_adapter_rnastructure::RnaStructureAdapter::new(),
        ));
        self.registry
            .register(Arc::new(valenx_adapter_nupack::NupackAdapter::new()));
        // Biology (Phase 25) — quantum chemistry
        self.registry
            .register(Arc::new(valenx_adapter_psi4::Psi4Adapter::new()));
        self.registry
            .register(Arc::new(valenx_adapter_nwchem::NwchemAdapter::new()));
        self.registry
            .register(Arc::new(valenx_adapter_xtb::XtbAdapter::new()));
        // Biology (Phase 27.6) — EvolutionaryScale models
        self.registry
            .register(Arc::new(valenx_adapter_esm3::Esm3Adapter::new()));
        self.registry
            .register(Arc::new(valenx_adapter_esmc::EsmcAdapter::new()));
        // Biology (Phase 32) — systems biology
        self.registry
            .register(Arc::new(valenx_adapter_copasi::CopasiAdapter::new()));
        self.registry
            .register(Arc::new(valenx_adapter_bionetgen::BioNetGenAdapter::new()));
        self.registry
            .register(Arc::new(valenx_adapter_physicell::PhysiCellAdapter::new()));
        // Biology (Phase 36) — cryo-EM
        self.registry
            .register(Arc::new(valenx_adapter_relion::RelionAdapter::new()));
        self.registry
            .register(Arc::new(valenx_adapter_eman2::Eman2Adapter::new()));
        self.registry
            .register(Arc::new(valenx_adapter_ctffind::CtffindAdapter::new()));
        // Biology (Phase 31) — sequencing read simulators
        self.registry
            .register(Arc::new(valenx_adapter_art::ArtAdapter::new()));
        self.registry
            .register(Arc::new(valenx_adapter_wgsim::WgsimAdapter::new()));
        self.registry
            .register(Arc::new(valenx_adapter_badread::BadreadAdapter::new()));
        // Biology (Phase 35) — CRISPR design
        self.registry
            .register(Arc::new(valenx_adapter_chopchop::ChopchopAdapter::new()));
        self.registry
            .register(Arc::new(valenx_adapter_crispor::CrisporAdapter::new()));
        self.registry.register(Arc::new(
            valenx_adapter_cas_offinder::CasOffinderAdapter::new(),
        ));
        // Biology (Phase 38) — Rosetta family
        self.registry
            .register(Arc::new(valenx_adapter_rosetta::RosettaAdapter::new()));
        self.registry
            .register(Arc::new(valenx_adapter_pyrosetta::PyRosettaAdapter::new()));
        // Biology (Phase 29) — population genetics
        self.registry
            .register(Arc::new(valenx_adapter_slim::SlimAdapter::new()));
        self.registry
            .register(Arc::new(valenx_adapter_msprime::MsprimeAdapter::new()));
        self.registry
            .register(Arc::new(valenx_adapter_tskit::TskitAdapter::new()));
        // Biology (Phase 30.5) — Bayesian phylogenetics
        self.registry
            .register(Arc::new(valenx_adapter_beast2::Beast2Adapter::new()));
        self.registry
            .register(Arc::new(valenx_adapter_mrbayes::MrBayesAdapter::new()));
        // Biology (Phase 39) — DNA structural geometry
        self.registry
            .register(Arc::new(valenx_adapter_x3dna::X3dnaAdapter::new()));
        self.registry
            .register(Arc::new(valenx_adapter_curves::CurvesAdapter::new()));
        self.registry
            .register(Arc::new(valenx_adapter_dssr::DssrAdapter::new()));
        // Biology (Phase 5.5) — MD analysis expansion
        self.registry
            .register(Arc::new(valenx_adapter_plumed::PlumedAdapter::new()));
        self.registry
            .register(Arc::new(valenx_adapter_prody::ProdyAdapter::new()));
        self.registry
            .register(Arc::new(valenx_adapter_cpptraj::CpptrajAdapter::new()));
        // Biology (Phase 33) — synthetic biology
        self.registry
            .register(Arc::new(valenx_adapter_pysbol::PySbolAdapter::new()));
        self.registry
            .register(Arc::new(valenx_adapter_j5::J5Adapter::new()));
        self.registry
            .register(Arc::new(valenx_adapter_cello::CelloAdapter::new()));
        // Biology (Phase 18.7) — alignment toolkit expansion
        self.registry
            .register(Arc::new(valenx_adapter_blast::BlastAdapter::new()));
        self.registry
            .register(Arc::new(valenx_adapter_clustalo::ClustaloAdapter::new()));
        self.registry
            .register(Arc::new(valenx_adapter_tcoffee::TCoffeeAdapter::new()));
        // Biology (Phase 19.6) — single-cell genomics expansion
        self.registry
            .register(Arc::new(valenx_adapter_seurat::SeuratAdapter::new()));
        self.registry
            .register(Arc::new(valenx_adapter_anndata::AnnDataAdapter::new()));
        // Biology (Phase 5.6) — bio MD engines
        self.registry
            .register(Arc::new(valenx_adapter_namd::NamdAdapter::new()));
        self.registry
            .register(Arc::new(valenx_adapter_amber_sander::SanderAdapter::new()));
        self.registry
            .register(Arc::new(valenx_adapter_hoomd::HoomdAdapter::new()));
        // Biology (Phase 5.7) — MD analysis sister
        self.registry
            .register(Arc::new(valenx_adapter_mdtraj::MdtrajAdapter::new()));
        // Biology (Phase 17.7) — structure prediction + search expansion
        self.registry.register(Arc::new(
            valenx_adapter_rosettafold::RoseTTAFoldAdapter::new(),
        ));
        self.registry
            .register(Arc::new(valenx_adapter_omegafold::OmegaFoldAdapter::new()));
        self.registry
            .register(Arc::new(valenx_adapter_foldseek::FoldseekAdapter::new()));
        // Biology (Phase 32.5) — spatial stochastic reaction-diffusion
        self.registry
            .register(Arc::new(valenx_adapter_smoldyn::SmoldynAdapter::new()));
        self.registry
            .register(Arc::new(valenx_adapter_mcell::McellAdapter::new()));
        // Biology (Phase 41) — sequence editors / plasmid design
        self.registry
            .register(Arc::new(valenx_adapter_pydna::PydnaAdapter::new()));
        self.registry
            .register(Arc::new(valenx_adapter_jalview::JalviewAdapter::new()));
        // Biology (Phase 40) — microscopy / bioimage analysis
        self.registry
            .register(Arc::new(valenx_adapter_fiji::FijiAdapter::new()));
        self.registry.register(Arc::new(
            valenx_adapter_cellprofiler::CellProfilerAdapter::new(),
        ));
        self.registry
            .register(Arc::new(valenx_adapter_ilastik::IlastikAdapter::new()));
        // Biology (Phase 22.5) — workflow expansion
        self.registry
            .register(Arc::new(valenx_adapter_planemo::PlanemoAdapter::new()));
        self.registry
            .register(Arc::new(valenx_adapter_cromwell::CromwellAdapter::new()));
        self.registry
            .register(Arc::new(valenx_adapter_cwltool::CwltoolAdapter::new()));
        // Biology (Phase 42) — modern web 3D molecular visualization
        self.registry
            .register(Arc::new(valenx_adapter_molstar::MolstarAdapter::new()));
        self.registry
            .register(Arc::new(valenx_adapter_ngl::NglAdapter::new()));
        // Biology (Phase 43) — mRNA design / codon optimization
        self.registry
            .register(Arc::new(valenx_adapter_dnachisel::DnaChiselAdapter::new()));
        self.registry.register(Arc::new(
            valenx_adapter_lineardesign::LinearDesignAdapter::new(),
        ));
        self.registry
            .register(Arc::new(valenx_adapter_icodon::IcodonAdapter::new()));
        // Biology (Phase 44.5) — RNA folding expansion
        self.registry
            .register(Arc::new(valenx_adapter_mfold::MfoldAdapter::new()));
        self.registry
            .register(Arc::new(valenx_adapter_eternafold::EternaFoldAdapter::new()));
        self.registry
            .register(Arc::new(valenx_adapter_linearfold::LinearFoldAdapter::new()));
        // Biology (Phase 35.5) — base + prime editing design
        self.registry.register(Arc::new(
            valenx_adapter_be_designer::BeDesignerAdapter::new(),
        ));
        self.registry
            .register(Arc::new(valenx_adapter_be_hive::BeHiveAdapter::new()));
        self.registry.register(Arc::new(
            valenx_adapter_primedesign::PrimeDesignAdapter::new(),
        ));
        self.registry
            .register(Arc::new(valenx_adapter_pegfinder::PegFinderAdapter::new()));
        // Biology (Phase 35.6) — edit-outcome prediction
        self.registry
            .register(Arc::new(valenx_adapter_indelphi::IndelphiAdapter::new()));
        self.registry
            .register(Arc::new(valenx_adapter_forecast::ForecastAdapter::new()));
        self.registry.register(Arc::new(
            valenx_adapter_alphamissense::AlphaMissenseAdapter::new(),
        ));
        self.registry
            .register(Arc::new(valenx_adapter_crispritz::CrispritzAdapter::new()));
        // Biology (Phase 45) — pharmacokinetics + RNA tertiary
        self.registry
            .register(Arc::new(valenx_adapter_pksim::PkSimAdapter::new()));
        self.registry
            .register(Arc::new(valenx_adapter_simrna::SimRnaAdapter::new()));
        // Probe adapters on a BACKGROUND thread so the first frame is not
        // blocked. Probing all 141 external tools (a PATH search + a
        // version-spawn for each installed one) sequentially on the main
        // thread froze startup for ~5s warm / ~30s cold. Results stream in
        // and are applied per frame in `update`; the native adapter paths
        // work regardless of probe state, so nothing waits on this.
        self.adapter_probe_rx = Some(self.registry.spawn_probe_all());
        tracing::info!(
            target: "valenx",
            registered = self.registry.len(),
            "registry initialised — probing adapters in background"
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_state_is_empty() {
        let app = ValenxApp::default();
        assert!(app.project.is_none());
        assert!(app.stl.is_none());
        assert!(app.last_error.is_none());
        assert!(app.run_handle.is_none());
        assert!(app.residuals.is_empty());
        assert!(app.selected_case.is_none());
    }

    #[test]
    fn default_state_has_no_prepared_workdir() {
        // Regression: prepare_selected_case must not have side-
        // effects from app construction. The field starts None and
        // only gets populated after a successful prepare call.
        let app = ValenxApp::default();
        assert!(app.last_prepare_workdir.is_none());
        // Same for the run workdir companion field.
        assert!(app.last_run_workdir.is_none());
    }

    #[test]
    fn default_state_has_no_prepared_job() {
        let app = ValenxApp::default();
        assert!(app.last_prepared_job.is_none());
        // last_run_results is the new sibling of last_run_report,
        // populated only after a successful run's collect() completes.
        assert!(app.last_run_results.is_none());
        // selected_field_name starts None — the field overlay falls
        // back to "first scalar OnNode field that matches the mesh"
        // until the user clicks one in the Results pane.
        assert!(app.selected_field_name.is_none());
        // selected_time_index defaults to 0 — first snapshot in the
        // selected field's time series.
        assert_eq!(app.selected_time_index, 0);
        // run_history starts empty and running_case_name unset.
        assert!(app.run_history.is_empty());
        assert!(app.running_case_name.is_none());
    }

    #[test]
    fn run_history_entry_round_trips_through_map() {
        // Smoke test: the map stores exactly what we put in. The
        // event-handler logic in pump_run_events that actually
        // writes these is integration-tested elsewhere; this test
        // just locks in the types so the browser badge code can
        // rely on the shape.
        let mut app = ValenxApp::default();
        app.run_history.insert(
            "cfd-steady".into(),
            RunHistoryEntry {
                succeeded: true,
                wall_time: std::time::Duration::from_secs(42),
                converged: Some(true),
            },
        );
        let entry = app.run_history.get("cfd-steady").expect("inserted");
        assert!(entry.succeeded);
        assert_eq!(entry.converged, Some(true));
        assert_eq!(entry.wall_time, std::time::Duration::from_secs(42));
    }

    // No direct test for `open_path_in_file_browser` — calling it
    // with any path actually spawns the host's file browser, which
    // pops up an error dialog or a window in CI. The function's
    // public Result<(), String> signature is the contract; the
    // no-workdir error paths above test the callers, which is what
    // matters for app behaviour.
}
