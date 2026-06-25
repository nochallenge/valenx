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
pub mod agent_commands;
pub mod animate_workbench;
pub mod antenna_workbench;
pub(crate) mod background;
pub mod batterypack_workbench;
pub mod beam_workbench;
pub mod blackhole_workbench;
pub mod cad_workbench;
pub mod car_workbench;
pub mod cfd_workbench;
pub mod confidence;
pub mod draft2d_workbench;
pub mod fem_workbench;
pub mod headless;
pub mod heatpump_workbench;
pub mod hvac_workbench;
pub mod inductionmotor_workbench;
pub mod interior_workbench;
pub mod neuro_workbench;
pub mod reinforcement_workbench;
pub mod render_workbench;
pub mod reverse_workbench;
pub mod variant_effect_workbench;
pub mod windturbine_workbench;

pub mod acidbase_workbench;
pub mod acoustics_workbench;
pub mod assistant_workbench;
pub mod astro;
pub mod astro_workbench;
pub mod bearing_workbench;
pub mod beltdrive_workbench;
pub mod bjt_workbench;
pub mod bmr_workbench;
pub mod bolt_workbench;
pub mod bonemech_workbench;
pub mod bracket_product;
pub mod brake_workbench;
pub mod buckling_workbench;
pub mod cam_overlay;
pub mod capacitor_workbench;
pub mod chaindrive_workbench;
pub mod clutch_workbench;
pub mod coil_workbench;
pub mod collision_workbench;
pub mod columnsteel_workbench;
pub mod combustion_workbench;
pub mod commands;
pub mod conveyor_workbench;
pub mod cosim_workbench;
#[cfg(test)]
mod coverage_ui_tests;
pub mod creep_workbench;
pub mod dcmotor_workbench;
pub mod dna_product;
pub mod dock_layout;
pub mod docking;
pub mod draft_overlay;
pub mod drone_workbench;
pub mod electrochem_workbench;
pub mod engine_workbench;
pub mod enzymekinetics_workbench;
pub mod fanlaws_workbench;
pub mod fasteners_workbench;
pub mod fatigue_workbench;
pub mod fields_workbench;
pub mod first_run;
pub mod fixedwing_workbench;
pub mod fluidstatics_workbench;
pub mod flywheel_workbench;
pub mod fourbar_workbench;
#[cfg(test)]
mod widget_naming_tests;
// Mechanical + civil batch — surface valenx-shaftdesign, -screwthread,
// -pulley, -spring-design, -springcombination, -vibration, -rivet,
// -soilbearing as reactive right-side workbenches.
pub mod fracture_workbench;
pub mod frames_workbench;
pub mod gasdynamics_workbench;
pub mod gearbox_workbench;
// Science batch 5 — surface valenx-camdynamics, -battery-ecm, -diffusion,
// -dimensional, -fft as reactive right-side workbenches.
pub mod batteryecm_workbench;
pub mod camdynamics_workbench;
pub mod diffusion_workbench;
pub mod dimensional_workbench;
pub mod fft_workbench;
pub mod gears_workbench;
pub mod geartooth_workbench;
pub mod genetics;
pub mod genetics_workbench;
pub mod geomatics_workbench;
// EE / DSP workbenches (electronics batch) — surface valenx-opamp, -led,
// -thermocouple, -transmissionline, -powerfactor, -resistor-network,
// -rectifier, -filter as reactive right-side workbenches.
pub mod filter_workbench;
pub mod heatexchanger_workbench;
pub mod heattransfer_workbench;
pub mod hydraulics_workbench;
pub mod inclinedplane_workbench;
pub mod insulation_workbench;
pub mod keyboard_help;
pub mod landing_page;
pub mod leadscrew_workbench;
pub mod led_workbench;
pub mod leverage_workbench;
pub mod log_panel;
pub mod marine_workbench;
pub mod materials;
pub mod mbd_workbench;
pub mod mesh_prims;
pub mod mesh_toolbox;
pub mod mission_planner_workbench;
pub mod missionsim_workbench;
pub mod mohr_workbench;
/// Richer molecular-viewer representations (sticks / cartoon / marching-cubes
/// surface) extending [`genetics::molecule_view`]; pure mesh generators wired
/// into the Macromolecular-Structure panel's representation picker.
pub mod molviz;
pub mod mosfet_workbench;
pub mod new_project_dialog;
pub mod opamp_workbench;
pub mod optics_workbench;
pub mod orifice_workbench;
pub mod param_sketch_panel;
pub mod pbr_forward_pass;
pub mod pharmacokinetics_workbench;
pub mod photogrammetry_workbench;
pub mod pipeflow_workbench;
pub mod pipenetwork_workbench;
pub mod piping_workbench;
// Science / bio / civil batch — surface valenx-retainingwall, -openchannel,
// -weir, -thermocycle, -queueing, -radioactivity, -osmosis, -thermoreg,
// -hemodynamics, -popdynamics as reactive right-side workbenches.
pub mod autonomy_workbench;
pub mod fluids_workbench;
pub mod hemodynamics_workbench;
pub mod ocean_workbench;
pub mod openchannel_workbench;
pub mod osmosis_workbench;
pub mod plate_workbench;
pub mod pneumatics_workbench;
pub mod popdynamics_workbench;
pub mod powerfactor_workbench;
pub mod ppi_workbench;
pub mod pressurevessel_workbench;
/// Per-file registry of agent-bridge `show_3d` mesh producers (replaces the old
/// per-kind reducer arms; new 3-D tools register from their own module).
pub mod products_registry;
pub mod project_library;
pub mod project_navigator;
pub mod project_tabs;
pub mod projectile_workbench;
pub mod psychrometrics_workbench;
pub mod pulley_workbench;
pub mod pump_workbench;
pub mod queueing_workbench;
pub mod radioactivity_workbench;
pub mod rail_workbench;
pub mod rcbeam_workbench;
pub mod reactdyn_workbench;
pub mod rectifier_workbench;
pub mod refrigeration_workbench;
pub mod residuals;
pub mod resistornetwork_workbench;
pub mod retainingwall_workbench;
pub mod rivet_workbench;
pub mod rocket_mesh;
pub mod rocket_workbench;
pub mod rom_workbench;
pub mod rotor_workbench;
pub mod run;
pub mod scene_overlay;
pub mod screwthread_workbench;
pub mod sensors_workbench;
pub mod settings;
pub mod setup;
pub mod shaftdesign_workbench;
pub mod sheetmetal_workbench;
pub mod shortcuts;
pub mod sketch_overlay;
pub mod soilbearing_workbench;
pub mod solarpv_workbench;
pub mod springcombination_workbench;
pub mod springdesign_workbench;
pub mod springs_workbench;
pub mod statics_workbench;
pub mod straingauge_workbench;
pub mod strainrosette_workbench;
pub mod survivability_workbench;
pub mod thermalexpansion_workbench;
pub mod thermistor_workbench;
pub mod thermocouple_workbench;
pub mod thermocycle_workbench;
pub mod thermoreg_workbench;
pub mod threephase_workbench;
pub mod torsion_workbench;
pub mod transformer_workbench;
pub mod transmissionline_workbench;
pub mod truss_workbench;
pub mod types;
pub mod uas_workbench;
pub mod undo;
pub mod uq_workbench;
pub mod vibration_workbench;
pub mod viewport;
pub mod viewport_2d;
pub mod viewport_kind;
pub mod weir_workbench;
pub mod welcome_tour;
pub mod wgpu_renderer;
pub mod workbench_chrome;
pub mod workbench_focus;

// Concern-focused helper modules — what used to be a single
// `helpers.rs` bag-of-everything (1422 LOC, 36 fns spanning 8+
// unrelated concerns). Sibling modules let callers `use
// crate::history::save_run_history_to_state_dir` and have the
// import name actually tell them which concern they're reaching
// into.
pub(crate) mod adapter_status;
pub mod audit;
pub mod file_browser;
pub mod history;
pub(crate) mod mesh_loader;
pub mod rbac_io;
pub mod settings_io;
pub mod state_paths;

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
pub use crate::state_paths::state_dir;
pub use crate::types::{BottomTab, LoadedMesh, LoadedStl, RunHistoryEntry, SweepHistoryEntry};

// Stage A1 of the app split (docs/refactor/2026-06-20-valenx-app-split.md):
// these leaf modules moved to `valenx-app-core`. Re-export them here so
// the existing public API (`valenx_app::theme`, `valenx_app::tooltips`,
// `valenx_app::panel_help`, `valenx_app::workbench_ui`,
// `valenx_app::format_time_key`, `valenx_app::adapter_id_from_solver`, …)
// stays stable and in-crate `crate::theme::…` paths keep resolving.
pub use valenx_app_core::solver_parse::{adapter_id_from_solver, derived_inputs_from_case_toml};
pub use valenx_app_core::time_format::format_time_key;
pub use valenx_app_core::{
    histograms, menu_ui, panel_help, solver_parse, theme, time_format, tooltips, workbench_ui,
};

/// The **finished build result** an external agent posts into a "Workbench +
/// Agent" unit's workspace tile, so the agent's output (e.g. a sized rocket, a
/// gear train) shows up in *this* pane instead of staying a placeholder.
///
/// Set per unit `n` by the [`crate::agent_commands::AgentCommand::ShowProduct`]
/// (text card) or [`crate::agent_commands::AgentCommand::Show3d`] (live 3-D
/// view) bridge commands and rendered by [`crate::dock_layout`]'s
/// `render_workspace_body`:
///
/// - a text result is a bold `title` heading over a list of plain-text `lines`
///   (one row each);
/// - a **3-D** result additionally carries a [`LoadedMesh`] in [`Self::mesh`]
///   and a fixed [`OrbitCamera`] in [`Self::camera`], which the pane renders as
///   an actual lit 3-D view (same look as the central viewport) at a fixed 3/4
///   angle.
///
/// Not `Clone`/`Default`/`Debug`: [`LoadedMesh`] owns a `valenx_mesh::Mesh` +
/// quality reports and implements none of them. The only writer is the
/// agent-command reducer's `insert` and the only reader is
/// `render_workspace_body`'s `get` (both move/borrow, never clone or format),
/// so none of those bounds is needed.
pub struct WorkspaceProduct {
    /// Card heading (rendered bold), e.g. the product name.
    pub title: String,
    /// Result rows shown under the heading, one `ui.label` per entry. Empty
    /// for a pure 3-D product.
    pub lines: Vec<String>,
    /// When `Some`, the pane renders this mesh as a live lit 3-D view (using
    /// [`Self::camera`]) instead of a text card. Built by the `show_3d`
    /// command (e.g. the LV-1 rocket via
    /// `crate::rocket_workbench::lv1_loaded_mesh`).
    pub mesh: Option<LoadedMesh>,
    /// Optional per-vertex base colours for [`Self::mesh`], one `[r, g, b]`
    /// in `[0, 1]` per surface vertex of the mesh's renderable triangle skin
    /// (the order [`crate::wgpu_renderer::triangles_to_vertices`] emits:
    /// triangle-major, then the three vertices of each triangle). `None` for
    /// plain meshes (rocket / gear / bracket / rcbeam) — those render in the
    /// neutral brushed-metal colour. The FEM cantilever product sets this to a
    /// von-Mises stress colormap so the deformed shape reads as a stress map
    /// rather than flat grey. When the length matches the mesh's surface vertex
    /// count the tile renders with [`crate::wgpu_renderer::triangles_to_vertices_colored`];
    /// otherwise it falls back to the plain metal path. Ignored when
    /// [`Self::mesh`] is `None`.
    pub vertex_colors: Option<Vec<[f32; 3]>>,
    /// Fixed camera the 3-D view is rendered from (a pleasant 3/4 angle for
    /// Stage 1 — per-tile orbit is a later stage). Ignored when
    /// [`Self::mesh`] is `None`.
    pub camera: OrbitCamera,
    /// When `Some`, the pane renders a **2-D engineering drawing** painted with
    /// egui (no wgpu) — an RC-beam section + rebar, or a DNA construct map —
    /// instead of the 3-D viewport or the text card. Set by the
    /// [`crate::agent_commands::AgentCommand::Show2d`] bridge command and painted
    /// by [`crate::dock_layout`]'s `render_workspace_body` (a branch *between*
    /// the 3-D viewport and the text card). `None` for plain / mesh / text
    /// products. See [`Workspace2dKind`].
    pub kind2d: Option<Workspace2dKind>,
    /// Transient status line for the tile's "Export STL" action, shown next to
    /// the button so the user sees where the mesh was written (e.g.
    /// `"saved → C:\\…\\Downloads\\valenx_rocket.stl"`) or why it failed
    /// (`"export failed: …"`). `None` until the first export. Only meaningful
    /// for mesh products ([`Self::mesh`] is `Some`); other product kinds never
    /// show the button and so never set it. Set by
    /// [`crate::dock_layout`]'s `render_workspace_body` on button click.
    pub last_export: Option<String>,
    /// When `Some`, the pane renders this raster **image** scaled to fit the
    /// tile — the path-traced `render` view (a small Cornell-box framebuffer)
    /// is the canonical producer. The pixels are carried as a CPU
    /// [`egui::ColorImage`]; the first frame uploads them to a GPU texture
    /// cached in [`Self::image_texture`] and then draws with `ui.image`.
    /// Rendered by [`crate::dock_layout`]'s `render_workspace_body` in a branch
    /// *before* the text-card fall-through (a peer to the 3-D and 2-D
    /// branches). `None` for mesh / 2-D / text products.
    pub image: Option<egui::ColorImage>,
    /// Lazily-uploaded GPU texture for [`Self::image`] — `None` until the first
    /// frame that renders the image, which uploads the [`egui::ColorImage`]
    /// once (keyed by the tile id) and caches the handle here so subsequent
    /// frames reuse it instead of re-allocating a texture every repaint. The
    /// handle frees the GPU texture when the product is dropped (RAII). Never
    /// set for non-image products.
    pub image_texture: Option<egui::TextureHandle>,
    /// When `Some`, this mesh product is **animated**: the toolbar shows
    /// Play/Pause + speed + reset controls and, while playing, the per-tile 3-D
    /// view re-poses the mesh nodes each frame from [`ProductAnimation::t`]
    /// (see [`crate::dock_layout`]'s render path). `None` for static products
    /// (the existing byte-identical static render). Only meaningful when
    /// [`Self::mesh`] is `Some` — a 2-D / text / image product never animates.
    /// The 2-stage spur reducer (`gears_workbench::gear_product`) sets this to a
    /// [`ProductMotion::RigidParts`] so each gear counter-rotates and the teeth
    /// visibly mesh; other mesh products leave it `None` for now.
    pub animation: Option<ProductAnimation>,
}

/// Real-time **animation state** for an animated [`WorkspaceProduct`].
///
/// The toolbar advances [`Self::t`] (seconds of *animation* time) each frame
/// while [`Self::playing`], scaled by [`Self::speed`]; the render path then
/// poses the product's mesh from `t` via [`Self::motion`]. Kept a plain
/// `Clone + PartialEq` value (no GPU / mesh handles) so products stay cheap to
/// snapshot.
#[derive(Clone, Debug, PartialEq)]
pub struct ProductAnimation {
    /// Whether the clock is currently running (Play vs Pause).
    pub playing: bool,
    /// Playback-speed multiplier applied to wall-clock dt, in `0.0..=4.0`
    /// (the toolbar slider's range). `1.0` is real-time relative to the motion's
    /// own `rad_per_s`.
    pub speed: f32,
    /// Elapsed **animation** time in seconds. Drives every motion's angle
    /// (`rad_per_s * t`). Reset to `0.0` by the toolbar's reset button.
    pub t: f32,
    /// What moves and how.
    pub motion: ProductMotion,
}

/// The kind of rigid motion an animated product performs.
#[derive(Clone, Debug, PartialEq)]
pub enum ProductMotion {
    /// Spin the **whole mesh** about a fixed world axis through `pivot` — a
    /// generic turntable for any single-body product. Angle = `rad_per_s * t`.
    Turntable {
        /// World-space rotation axis (need not be unit length; normalised at
        /// pose time).
        axis: [f32; 3],
        /// World-space point the axis passes through.
        pivot: [f32; 3],
        /// Angular rate (radians per second of animation time).
        rad_per_s: f32,
    },
    /// Per-part **rigid rotation**: each [`RigidPart`] spins its own contiguous
    /// node range about its own axis/pivot. Used by the gear train so meshing
    /// gears counter-rotate independently. Nodes outside every part's range stay
    /// fixed.
    RigidParts(Vec<RigidPart>),
}

/// One independently-rotating body inside a [`ProductMotion::RigidParts`]
/// product: the half-open `node_range` of mesh nodes to rotate, about `axis`
/// through `pivot`, at `rad_per_s`. The ranges are recorded at mesh-fusion time
/// (the order parts are concatenated into the combined node array), so they tile
/// the mesh's node count without overlap.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidPart {
    /// Half-open `[start, end)` index range into the mesh's `nodes` array.
    pub node_range: std::ops::Range<usize>,
    /// World-space rotation axis (normalised at pose time).
    pub axis: [f32; 3],
    /// World-space point the axis passes through (this body's shaft centre).
    pub pivot: [f32; 3],
    /// Signed angular rate (rad/s); the sign encodes spin direction so meshing
    /// pairs counter-rotate.
    pub rad_per_s: f32,
}

impl WorkspaceProduct {
    /// Give a mesh product a **default inspect-spin** animation if it has none,
    /// so every bridge-rendered 3-D product carries the Play/Pause + speed
    /// toolbar (an idle turntable the user can start to inspect the part from
    /// all sides).
    ///
    /// Sets [`Self::animation`] to a paused [`ProductMotion::Turntable`] about
    /// `+Z` through the mesh's axis-aligned bounding-box centre at
    /// `0.4 rad/s` (~1 revolution / 15 s — a gentle inspect spin), but **only**
    /// when both:
    ///
    /// - [`Self::animation`] is currently `None` (so a product that already has
    ///   real motion — e.g. the gear train's [`ProductMotion::RigidParts`] — is
    ///   left untouched), and
    /// - [`Self::mesh`] is `Some` (a 2-D / card / image product, which never
    ///   animates, is left untouched).
    ///
    /// `playing` starts `false`: the control is present but the product does
    /// **not** auto-spin until the user (or an `animate` command) presses Play.
    /// A no-op in every other case, so it is safe to call unconditionally right
    /// after a product is built. The pivot is the `(min + max) * 0.5` midpoint
    /// of the mesh's node AABB (via `mesh_loader::mesh_bounding_box`); an empty
    /// mesh falls back to the origin.
    pub fn ensure_default_animation(&mut self) {
        if self.animation.is_some() {
            return;
        }
        let Some(loaded) = self.mesh.as_ref() else {
            return;
        };
        let pivot = match crate::mesh_loader::mesh_bounding_box(&loaded.mesh) {
            Some((min, max)) => [
                (min[0] + max[0]) * 0.5,
                (min[1] + max[1]) * 0.5,
                (min[2] + max[2]) * 0.5,
            ],
            None => [0.0, 0.0, 0.0],
        };
        self.animation = Some(ProductAnimation {
            playing: false,
            speed: 1.0,
            t: 0.0,
            motion: ProductMotion::Turntable {
                axis: [0.0, 0.0, 1.0],
                pivot,
                rad_per_s: 0.4,
            },
        });
    }
}

/// Which **2-D engineering drawing** a [`WorkspaceProduct`] carries, with the
/// small plain-data view the egui painter needs (no wgpu / mesh types, so the
/// whole thing is cheaply `Clone`). Painted by [`crate::dock_layout`]'s
/// `render_workspace_body`.
#[derive(Clone, Debug, PartialEq)]
pub enum Workspace2dKind {
    /// A reinforced-concrete **beam section + rebar** drawing — a filled
    /// concrete rectangle with the tension bars near the bottom, dimension
    /// lines, and the flexural numbers.
    RcSection(RcSectionView),
    /// A DNA **construct map** — a horizontal baseline = the construct with each
    /// feature (ATG / ORF / His6 / stop) drawn as a coloured block proportional
    /// to its nt span, plus a nt ruler.
    DnaMap(DnaMapView),
    /// A **2-D CAD drawing** — the LibreCAD-style entity soup (lines / circles /
    /// arcs / polylines) of the 2-D drafting workbench, painted as a
    /// fit-to-tile drawing.
    Draft2d(Draft2dView),
    /// An **interior floor-plan** — one or more room wall polygons plus the
    /// placed-furniture rectangles, painted as a fit-to-tile plan.
    FloorPlan(FloorPlanView),
    /// A **2-D line / bar chart** — one or more `(x, y)` series drawn as
    /// polylines (or vertical bars) inside a framed, auto-scaled plot area with
    /// gridlines, tick labels, axis labels and a small legend. The natural
    /// presentation for products that are really *plots* (an FFT magnitude
    /// spectrum, a concentration-vs-position profile, a population-vs-time
    /// trajectory) rather than a 3-D blob or a bare text card. See
    /// [`ChartData`].
    Chart(ChartData),
}

/// Plain-data backing for a [`Workspace2dKind::Chart`] — everything the egui
/// chart painter needs and nothing it doesn't (no GPU / mesh types, so it stays
/// cheaply `Clone`). The painter auto-scales the axes from the union of every
/// series' point range, so the producer only supplies the raw `(x, y)` data and
/// the labels.
#[derive(Clone, Debug, PartialEq)]
pub struct ChartData {
    /// Chart heading, drawn at the top of the plot (e.g. `"FFT magnitude
    /// spectrum"`).
    pub title: String,
    /// Label for the horizontal axis (e.g. `"frequency (Hz)"`).
    pub x_label: String,
    /// Label for the vertical axis (e.g. `"|X[k]|"`).
    pub y_label: String,
    /// The plotted series, in draw order. Each is a polyline or a bar group; the
    /// painter colours them from a small fixed palette and lists them in the
    /// legend. An empty list paints just the framed, label-only plot.
    pub series: Vec<ChartSeries>,
}

/// One labelled data series on a [`ChartData`] — a sequence of `[x, y]` points
/// drawn either as a connected polyline (`bars == false`) or as vertical bars
/// rising from `y = 0` (`bars == true`).
#[derive(Clone, Debug, PartialEq)]
pub struct ChartSeries {
    /// Legend label for the series (e.g. `"infectious"`).
    pub label: String,
    /// The `[x, y]` data points, in x order. The painter maps them through the
    /// auto-scaled axis transform; fewer than two points still draws the
    /// markers / bars it can.
    pub points: Vec<[f64; 2]>,
    /// When `true`, draw each point as a vertical bar from the x-axis baseline
    /// up to its `y` value (a spectrum / histogram look); when `false`, connect
    /// the points with a polyline (a curve-vs-time look).
    pub bars: bool,
}

/// Plain-data view for the RC-beam **section drawing** ([`Workspace2dKind::RcSection`]).
/// Geometry in millimetres plus the already-formatted flexural readout rows
/// (the same `lines` the 3-D / text product carries) so the painter can show
/// the key numbers without re-deriving them.
#[derive(Clone, Debug, PartialEq)]
pub struct RcSectionView {
    /// Section width `b` (mm) — the horizontal extent of the drawn rectangle.
    pub width_mm: f64,
    /// Section depth `h` (mm) — the vertical extent of the drawn rectangle.
    pub depth_mm: f64,
    /// Clear cover to the bar centres (mm) — how far the rebar sits in from the
    /// faces.
    pub cover_mm: f64,
    /// Diameter of each tension bar (mm) — sets the drawn circle size.
    pub bar_dia_mm: f64,
    /// Number of tension bars drawn across the bottom.
    pub n_bars: usize,
    /// The flexural readout rows (Mn, φMn, ρ vs ρ_bal, utilisation, …), shown as
    /// a small text block beside the section.
    pub lines: Vec<String>,
}

/// One labelled feature span on a [`DnaMapView`] — a half-open `[start, end)`
/// nucleotide interval with a display label and an RGB block colour.
#[derive(Clone, Debug, PartialEq)]
pub struct DnaFeature {
    /// Short label drawn by the block (e.g. `"ATG"`, `"ORF"`, `"His6"`, `"stop"`).
    pub label: String,
    /// Inclusive start nucleotide index (0-based).
    pub start: usize,
    /// Exclusive end nucleotide index.
    pub end: usize,
    /// Block fill colour as `[r, g, b]` (0–255).
    pub color: [u8; 3],
}

/// Plain-data view for the DNA **construct map** ([`Workspace2dKind::DnaMap`]).
/// The total construct length plus its feature spans; the painter lays the
/// baseline across the tile and draws each feature proportional to its span.
#[derive(Clone, Debug, PartialEq)]
pub struct DnaMapView {
    /// Total construct length in nucleotides — the baseline's full extent.
    pub total_nt: usize,
    /// The labelled feature spans (ATG / ORF / His6 / stop), each a
    /// `[start, end)` nt interval with a colour. Ordered for drawing.
    pub features: Vec<DnaFeature>,
}

/// One drawing primitive on a [`Draft2dView`] — the plain-data (`f64`,
/// drawing-unit) form of a `valenx_librecad_2d::Entity2D`, carrying only the
/// geometry the tile painter needs (no layer / DXF metadata). The painter maps
/// these to screen pixels with a fit-to-tile transform.
#[derive(Clone, Debug, PartialEq)]
pub enum Draft2dEntity {
    /// A straight segment from `a` to `b` (drawing units).
    Line {
        /// Start point `[x, y]`.
        a: [f64; 2],
        /// End point `[x, y]`.
        b: [f64; 2],
    },
    /// A full circle of `radius` about `centre` (drawing units).
    Circle {
        /// Centre `[x, y]`.
        centre: [f64; 2],
        /// Radius in drawing units.
        radius: f64,
    },
    /// A circular arc about `centre` from `start_angle_deg` to `end_angle_deg`
    /// (degrees, CCW), tessellated by the painter.
    Arc {
        /// Centre `[x, y]`.
        centre: [f64; 2],
        /// Radius in drawing units.
        radius: f64,
        /// Start angle (degrees, CCW from +x).
        start_angle_deg: f64,
        /// End angle (degrees, CCW from +x).
        end_angle_deg: f64,
    },
    /// A polyline through `vertices`; `closed` joins the last vertex back to
    /// the first.
    Polyline {
        /// Ordered `[x, y]` vertices (drawing units).
        vertices: Vec<[f64; 2]>,
        /// Whether the last vertex connects back to the first.
        closed: bool,
    },
}

/// Plain-data view for the **2-D CAD drawing** ([`Workspace2dKind::Draft2d`]).
/// The drawing's entities plus the model's overall extent (so the painter can
/// fit it to the tile), and the already-formatted readout rows (entity count,
/// extent) shown beside the drawing.
#[derive(Clone, Debug, PartialEq)]
pub struct Draft2dView {
    /// The drawing primitives to paint, in draw order.
    pub entities: Vec<Draft2dEntity>,
    /// Axis-aligned drawing-unit bounds `((min_x, min_y), (max_x, max_y))` of
    /// all entities — the box the painter scales to fit. A degenerate (empty)
    /// drawing leaves this at the default unit box.
    pub bounds: ([f64; 2], [f64; 2]),
    /// Readout rows shown beside the drawing (entity count, extent).
    pub lines: Vec<String>,
}

/// One furniture rectangle on a [`FloorPlanView`] — an axis-aligned footprint
/// centred at `centre` with `size` `[w, d]` (metres) and a short `label`.
#[derive(Clone, Debug, PartialEq)]
pub struct FloorPlanItem {
    /// Footprint centre `[x, y]` (metres, plan coordinates).
    pub centre: [f64; 2],
    /// Footprint size `[width, depth]` (metres).
    pub size: [f64; 2],
    /// Short label drawn in the rectangle (the furniture kind).
    pub label: String,
}

/// Plain-data view for the **interior floor-plan**
/// ([`Workspace2dKind::FloorPlan`]). The room wall polygons plus the placed
/// furniture footprints, the overall plan extent (so the painter can fit it to
/// the tile), and the readout rows (room / piece counts).
#[derive(Clone, Debug, PartialEq)]
pub struct FloorPlanView {
    /// Each room's wall polygon as ordered `[x, y]` vertices (metres). Drawn as
    /// a closed loop of wall segments.
    pub rooms: Vec<Vec<[f64; 2]>>,
    /// The placed-furniture footprints.
    pub furniture: Vec<FloorPlanItem>,
    /// Axis-aligned plan-unit bounds `((min_x, min_y), (max_x, max_y))` of the
    /// rooms — the box the painter scales to fit.
    pub bounds: ([f64; 2], [f64; 2]),
    /// Readout rows shown beside the plan (room + furniture counts).
    pub lines: Vec<String>,
}

/// Memoised command-palette entry list with its cache key. The key is
/// `(registry.len(), library.content_rev(), show_non_oss_adapters,
/// focus_category)`; the value is the built [`crate::commands::CommandKind`]
/// list. Aliased so [`ValenxApp::palette_cache`] stays under clippy's
/// type-complexity bar. See the cache-build site in `update.rs`.
type PaletteCache = Option<(
    usize,
    u64,
    bool,
    Option<String>,
    Vec<crate::commands::CommandKind>,
)>;

/// Root application state.
#[derive(Default)]
pub struct ValenxApp {
    /// Opt-in dockable / tiling central-panel layout (View → Docked
    /// layout). Default-built tile tree; only painted when
    /// [`ValenxApp::docked_layout`] is on. See [`docking`].
    pub docking: docking::DockingState,
    /// When true, the central panel renders the [`docking`] tile tree
    /// (resizable splits / tabs / close / drag) instead of the classic
    /// single-viewport layout. Default `false` (classic layout).
    pub docked_layout: bool,
    /// Per-panel chrome state for the right-side workbenches — keyed by the
    /// panel's stable `SidePanel` id, it records whether each is collapsed
    /// and whether it is docked / floating / popped out into its own OS
    /// window. Driven by [`workbench_chrome::workbench_shell`]; empty until a
    /// panel's header is first interacted with. See [`workbench_chrome`].
    pub workbench_chrome:
        std::collections::HashMap<String, crate::workbench_chrome::PanelChromeState>,

    /// Opt-in **dockable / tileable layout for the right-side workbench
    /// panels** (View → "Dockable panel layout (beta)"). When `true`, the
    /// run of `draw_<x>_workbench` dispatch in `update.rs` is replaced by a
    /// single [`egui_tiles`] tree (`dock_layout::draw_dock_layout`) that
    /// hosts every open workbench as a draggable / reorderable / splittable
    /// tile. Default `false` — the classic stacked right-side `SidePanel`
    /// layout is unchanged and stays the default. Distinct from
    /// [`ValenxApp::docked_layout`], which tiles the *central* viewport.
    ///
    /// **Per-tab.** This holds the *active* tab's value; switching tabs swaps
    /// it (with [`Self::dock_tree`] / [`Self::viewport_hidden`] /
    /// [`Self::viewport_collapsed`]) via
    /// [`project_tabs::WorkspaceDoc`], so a tab carrying a "Workbench +
    /// Agent" grid keeps it while a freshly-opened tab starts with the dock
    /// off. A newly-opened tab installs a default document → this is `false`.
    pub dock_enabled: bool,
    /// The lazily-built [`egui_tiles`] tree backing the dockable workbench
    /// layout. `None` until the first frame [`ValenxApp::dock_enabled`] is on
    /// (built from whichever workbenches are open then); thereafter
    /// `dock_layout::draw_dock_layout` syncs it each frame — adding a tile
    /// when a workbench is opened and dropping the tile when one is closed.
    /// Panes are panel-id `String`s (e.g. `"valenx_engine_workbench"`).
    ///
    /// **Per-tab.** This is the *active* tab's tree; switching tabs `take`s it
    /// into the outgoing tab's [`project_tabs::WorkspaceDoc`] and installs the
    /// incoming tab's (each tab owns its own dock layout). A newly-opened tab
    /// starts with `None` (a clean workspace — no other tab's agent grid).
    pub dock_tree: Option<egui_tiles::Tree<String>>,
    /// Monotonic counter for **"Workbench + Agent"** units launched into the
    /// dock (View → "New Workbench + Agent" / "Open 6 …"). Each unit is a
    /// paired `"workspace:<n>"` (empty build canvas) + `"agent:<n>"` (Claude
    /// chat) tile; this is the highest `n` handed out so far (default `0` =
    /// none). See [`dock_layout`]. Unlike the `DOCKABLE_PANELS` tiles, these
    /// are not gated on a `show_*` flag — they persist in [`Self::dock_tree`]
    /// until the user closes them.
    ///
    /// **Global (deliberately not per-tab).** Although [`Self::dock_tree`] is
    /// per-tab, this counter stays on the app so the `<n>` minted for each
    /// unit is unique **across all tabs** — `agent:<n>` is the unit's chat
    /// channel id (`valenx_chat_*_u<n>`), and per-tab counters would collide
    /// two different tabs' "Agent 1" onto the same channel. So each tab's
    /// [`Self::dock_tree`] may contain different `agent:<n>` panes, but every
    /// `<n>` is globally distinct.
    pub wb_agent_counter: usize,
    /// When `true`, the active tab's dock is a **clean agent product tab**: its
    /// [`Self::dock_tree`] holds **only** that unit's `[workspace:n | agent:n]`
    /// pair, and `dock_layout::draw_dock_layout` must NOT sync the flag-gated
    /// `dock_layout::DOCKABLE_PANELS` (notably the global
    /// `"valenx_assistant_panel"`) into it — so the agent-built product tab
    /// shows exactly one chat (its own `agent:n`) beside its workspace, never
    /// the global Assistant pane too. Set by
    /// `dock_layout::set_clean_workbench_agent_dock` (called from the `new_unit`
    /// bridge command) and `false` for every other tab, so the landing tab and
    /// manually-opened "Workbench + Agent" units keep the existing behaviour
    /// (the assistant tile may share their grid).
    ///
    /// **Per-tab.** Swapped in/out with [`Self::dock_tree`] /
    /// [`Self::dock_enabled`] via [`project_tabs::WorkspaceDoc`]; a
    /// newly-opened tab installs a default document → this is `false`.
    pub dock_agent_only: bool,
    /// Per-channel **cursor** for the agent-drives-valenx command bridge: how
    /// many lines of channel `n`'s command file
    /// ([`crate::agent_commands::cmd_path`]) have already been applied. On the
    /// first poll for a channel the cursor is seeded to the file's current line
    /// count so pre-existing history is **not** replayed on launch; thereafter
    /// only genuinely-new appended lines run. Defaults empty (derive). See
    /// [`crate::agent_commands::poll_and_apply_agent_commands`].
    pub agent_cmd_cursor: std::collections::HashMap<usize, usize>,
    /// Cursor for the **global** (no-`_u`-suffix) agent-command channel
    /// (`<base-dir>/valenx_chat_cmd.jsonl`): how many lines of that file have
    /// already been applied. Unlike [`Self::agent_cmd_cursor`] this channel is
    /// **not** keyed per unit and is polled on every poll regardless of
    /// [`Self::wb_agent_counter`], so an external agent can `new_unit` to
    /// bootstrap its own Workbench+Agent unit before any unit exists. `None`
    /// until the first poll sees the file → starts at line 0 (stale history is
    /// wiped at launch by [`crate::agent_commands::clear_command_files`]). See
    /// [`crate::agent_commands::poll_and_apply_agent_commands`].
    pub agent_global_cmd_cursor: Option<usize>,
    /// Last time the agent-command files were polled, used to throttle the disk
    /// reads to ~1/sec. `None` until the first poll (derive default). See
    /// [`crate::agent_commands`].
    pub last_agent_poll: Option<std::time::Instant>,
    /// Per-unit chat **input buffers** for the "Workbench + Agent" `agent:<n>`
    /// tiles, keyed by unit number `n`. Each unit's chat `TextEdit` binds to its
    /// own entry here (via `unit_chat_inputs.entry(n).or_default()`) so the six
    /// agent chats don't share one input box and mirror each other's typing.
    /// The classic base Assistant panel keeps using
    /// `crate::assistant_workbench::AssistantWorkbenchState::input` instead.
    /// Defaults empty (derive); an entry is created lazily the first time a unit
    /// chat is drawn. See [`dock_layout`] and
    /// `crate::assistant_workbench::assistant_chat_ui`.
    pub unit_chat_inputs: std::collections::HashMap<usize, String>,
    /// Per-unit **finished build result** for the "Workbench + Agent"
    /// `workspace:<n>` tiles, keyed by unit number `n`. An external agent posts
    /// one via [`crate::agent_commands::AgentCommand::ShowProduct`]; the matching
    /// `workspace:<n>` tile then renders it as a result card (replacing the empty
    /// "the agent's output will appear here" placeholder). Defaults empty
    /// (derive); the same `n` the bridge uses to post Notes to `feed_u<n>`. See
    /// [`WorkspaceProduct`] and [`crate::dock_layout`].
    pub workspace_products: std::collections::HashMap<usize, WorkspaceProduct>,

    /// **Lazy-build queue** for `new_unit`: unit number `n` → the product `kind`
    /// string its tab should render, deferred until that `workspace:<n>` pane is
    /// first viewed (or the unit is animated). `new_unit` inserts here and opens
    /// the tab INSTANTLY instead of building the 3-D product up front, so an
    /// agent fleet can open 130+ tabs in a burst without building every product
    /// at once (which briefly hung the app). The deferred build runs in
    /// `agent_commands::materialize_pending` (crate-private), called from
    /// `render_workspace_body` (first render of the pane) and the `animate`
    /// reducer; it moves the entry into [`Self::workspace_products`] and removes
    /// it here. Defaults empty (derive). A `kind`-less `new_unit` inserts
    /// nothing, so there is nothing to materialize.
    pub pending_products: std::collections::HashMap<usize, String>,

    pub project: Option<LoadedProject>,
    pub project_path: Option<PathBuf>,
    /// RBAC override block parsed from the loaded project's
    /// `project.toml`. Merged on top of the global `<state_dir>/rbac.json`
    /// at every permission check, so a sensitive project can promote
    /// or demote per-user roles without rewriting the global config.
    /// `None` when no project is loaded or when project.toml has no
    /// `[rbac]` block.
    pub project_rbac_override: Option<valenx_rbac::RbacConfig>,
    pub stl: Option<LoadedStl>,
    pub mesh: Option<LoadedMesh>,
    pub camera: OrbitCamera,
    pub shading: ShadingMode,
    pub last_error: Option<String>,
    pub status: Option<String>,
    pub about_open: bool,

    pub registry: AdapterRegistry,
    pub residuals: ResidualHistory,
    pub log: LogPanel,
    pub bottom_tab: BottomTab,
    /// When `true`, the bottom Residuals / Log dock collapses to just
    /// its thin header strip (the tab selectors + the collapse/expand
    /// toggle); the content body — residual plot, log text, or the
    /// empty-state placeholder — is skipped and the panel stops
    /// reserving vertical space. Toggled by the AI-drivable
    /// "Collapse panel" / "Expand panel" button in the header row.
    /// Defaults to `false` (expanded) via `#[derive(Default)]`.
    pub bottom_panel_collapsed: bool,

    /// When `true`, the left Browser panel collapses to a thin vertical
    /// bar holding only the AI-drivable "Expand panel" button; the heavy
    /// browser body (open-tabs list, navigator, Cases / Geometry / Mesh
    /// / Results sections) is skipped and the panel stops reserving its
    /// wide default width. Mirrors `bottom_panel_collapsed` for the
    /// bottom dock. Toggled by the "Collapse panel" / "Expand panel"
    /// button; separate from `show_browser` (the show/hide toggle).
    /// Defaults to `false` (expanded) via `#[derive(Default)]`.
    pub browser_collapsed: bool,

    /// Which case the user clicked on in the browser, if any. `None`
    /// falls back to the first case in the project when a run is
    /// started.
    pub selected_case: Option<String>,

    pub run_handle: Option<RunHandle>,
    /// Live threaded sweep runner. `Some(_)` while a sweep is
    /// executing; cleared when the worker thread finishes / fails.
    pub sweep_handle: Option<SweepHandle>,
    /// Per-sweep progress: (succeeded, failed, total). Updated as
    /// `SweepEvent::JobFinished` events come in. The numbers persist
    /// across the sweep_handle's lifetime so the status pane can
    /// keep showing the last result after the worker exits.
    pub sweep_progress: (usize, usize, usize),
    /// Status text for the active sweep — surfaced near the run
    /// progress in the UI.
    pub sweep_message: String,
    pub run_progress: f32,
    pub run_message: String,
    pub last_run_report: Option<Box<RunReport>>,
    pub last_run_error: Option<String>,

    /// Last successful prepare-only workdir, if any. Set by
    /// [`Self::prepare_selected_case`] so the UI can show the path
    /// and the "Open in file browser" action can act on it. `None`
    /// until the user clicks "Prepare".
    pub last_prepare_workdir: Option<PathBuf>,

    /// PreparedJob from the most recent successful prepare, kept so
    /// [`Self::run_from_prepared_workdir`] can run the solver against
    /// the user's hand-edited dicts without re-emitting them. The
    /// adapter id that produced this job lives alongside it because
    /// `spawn_prepared` needs to look the adapter back up in the
    /// registry.
    pub last_prepared_job: Option<(String, valenx_core::PreparedJob)>,

    /// Last completed run's workdir, captured when the run handle
    /// drops at the end of `Self::pump_run_events`. Mirrors
    /// `last_prepare_workdir` for the run pipeline so users can
    /// "Open in file browser" the dir holding their .vtu / .frd /
    /// .log artifacts after the solver finishes. `None` until the
    /// first run completes.
    pub last_run_workdir: Option<PathBuf>,

    /// Results bundle from the most recent successful run, populated
    /// when the worker thread sends `RunEvent::Collected`. Carries
    /// the parsed Field catalog (e.g. OpenFOAM's VTU fields), scalar
    /// records, artifact list, and provenance. `None` until the
    /// first run completes successfully.
    pub last_run_results: Option<Box<valenx_fields::Results>>,

    /// Which field the viewport's colour overlay is showing. Set by
    /// clicking a field name in the Results pane. `None` falls back
    /// to "first scalar OnNode field that matches the mesh" — same
    /// auto-pick used before the field selector landed.
    pub selected_field_name: Option<String>,

    /// Index into the selected field's time series — `0` = first
    /// snapshot, `1` = second, etc. Driven by the slider in the
    /// Results pane. Clamped every frame so the index can't outrun
    /// the time-series length when the user switches fields with
    /// different snapshot counts.
    pub selected_time_index: usize,

    /// Per-case run history — last outcome + wall time. Keyed by
    /// case name (project-local). Populated when a run finishes;
    /// surfaces in the case browser as a small ✓/✗ badge so users
    /// can see at a glance which cases they've already exercised
    /// without scrolling logs. Persisted to
    /// `<state_dir>/run-history.json` after every run so it
    /// survives app restarts.
    pub run_history: std::collections::BTreeMap<String, RunHistoryEntry>,
    /// Per-case sweep history. Mirrors `run_history` but for the
    /// sweep pipeline — recorded when a sweep finishes (sync or
    /// async) so the case browser can show "you swept this with N
    /// derived cases at `<ts>`". Persisted to
    /// `<state_dir>/sweep-history.json` so it survives an app
    /// restart.
    pub sweep_history: std::collections::BTreeMap<String, SweepHistoryEntry>,

    /// Case name currently being run, captured at spawn time so the
    /// Finished/Failed handlers can record the outcome under the
    /// right key even if the user has moved their cursor / changed
    /// `selected_case` while the solver was running.
    pub running_case_name: Option<String>,

    pub palette: CommandPalette,
    pub settings: Settings,
    pub settings_open: bool,
    pub theme_applied: bool,

    pub wgpu_renderer: Option<WgpuRenderer>,

    /// Whether the right-side Mesh Toolbox panel is visible. Defaults
    /// to `true` so it surfaces automatically as soon as a mesh /
    /// STL is loaded; the View menu and the command palette can hide
    /// it for users who want a clean viewport.
    pub show_mesh_toolbox: bool,
    /// Whether the left-side Browser panel is visible. Defaults to
    /// `true`; the ribbon toggle, View menu, and command palette can
    /// hide it to give the viewport the full width.
    pub show_browser: bool,
    /// Whether the viewport cursor snaps to the ground grid (Fusion-style):
    /// the live cursor coordinate snaps to the nearest grid node, with a
    /// marker drawn there. Defaults to `true`; toggled from the View menu.
    pub snap_to_grid: bool,
    /// When `true`, the central 2D/3D viewport body is hidden entirely — the
    /// central area shows a "viewport hidden" placeholder and the wgpu / 2D
    /// render is skipped. Driven by the viewport header's ✕ and its `⋯` menu,
    /// and by View → "Hide 3D viewport". **Per-tab** (swapped with the dock
    /// state via [`project_tabs::WorkspaceDoc`]); a newly-opened tab starts
    /// `false` (viewport shown). Defaults to `false`. See `update`.
    pub viewport_hidden: bool,
    /// When `true` (and not [`Self::viewport_hidden`]), only the central
    /// viewport's slim chrome header is drawn and its 2D/3D body is skipped —
    /// the viewport is "rolled up" to just its title + controls. Toggled by
    /// the header's `−` (minimize) icon. **Per-tab** (swapped with the dock
    /// state via [`project_tabs::WorkspaceDoc`]); a newly-opened tab starts
    /// `false` (body shown). Defaults to `false`. See `update`.
    pub viewport_collapsed: bool,
    /// Receiver for background adapter-probe results (see
    /// [`valenx_core::AdapterRegistry::spawn_probe_all`]). `Some` while the
    /// background probe is in flight; drained each frame in `update` and
    /// cleared to `None` when the probe thread finishes. Probing off the
    /// main thread keeps startup instant — it fixed a ~30s cold-start
    /// freeze (141 external tools probed synchronously in `new`).
    pub adapter_probe_rx:
        Option<std::sync::mpsc::Receiver<(&'static str, valenx_core::AdapterStatus)>>,
    /// Form-input state for the toolbox panel (translate deltas,
    /// scale factors, rotation axis + angle, mirror plane, cut-
    /// plane point + normal, repair tolerance). Cleared back to
    /// defaults on app construction; persisted across panel toggles.
    pub mesh_toolbox: crate::mesh_toolbox::MeshToolboxState,

    /// First CAD operand (operand "A" for boolean ops). Set when the
    /// user creates a primitive through the Part workbench section
    /// with the "Create as second" toggle off, and rewritten every
    /// time a boolean op runs (the result replaces operand A).
    pub current_solid: Option<valenx_cad::Solid>,
    /// Second CAD operand (operand "B"). Set when the user creates
    /// a primitive with the "Create as second" toggle on. Cleared
    /// whenever a boolean op consumes it so the toolbox is honest
    /// about needing a new B for the next op.
    pub second_solid: Option<valenx_cad::Solid>,

    /// First-launch wizard state. Loaded from
    /// `<state_dir>/first-run.json` on startup; defaults to a
    /// fresh, never-completed decision when the file doesn't exist.
    pub first_run_decision: valenx_first_run::FirstRunDecision,
    /// Whether the wizard's egui modal is open right now. Always
    /// initialised to `false` — the wizard never auto-opens because
    /// Valenx ships native Rust engines for every major simulation
    /// domain (external adapters are an optional power-user surface,
    /// so pushing first-time users to install OpenFOAM / GROMACS /
    /// Python contradicts the value proposition). Re-openable from
    /// the Settings menu's "Re-probe external tools" entry and the
    /// command palette.
    pub first_run_open: bool,
    /// Cached environment report. Built lazily on the frame the
    /// wizard opens, so the registry's probe results survive across
    /// frames without re-probing every redraw.
    pub first_run_report: Option<valenx_first_run::EnvironmentReport>,

    /// Loaded locale catalogue. Populated in `new()` from the
    /// embedded en-US baseline; future versions will pick the
    /// locale matching the user's OS preference and fall back to
    /// en-US when a translation is missing. Wrap in
    /// `Option<Arc<…>>`-style sharing if hot-swap becomes a
    /// requirement (it isn't yet — we set the locale once at
    /// startup).
    pub catalogue: valenx_i18n::LocaleCatalogue,

    /// Phase 21 — Macro recorder. UI panels append actions via
    /// `macro_recorder.record` when the user clicks a
    /// recordable button. `start_recording` / `stop_recording`
    /// flip the recorder state.
    pub macro_recorder: valenx_macro::MacroRecorder,

    /// Phase 22 — Add-on registry. Owns the in-memory list of
    /// installed add-ons + dispatches install/update/uninstall via
    /// the manual install-by-directory flow.
    pub addons: valenx_addons::AddonRegistry,
    /// Whether the Add-on Manager panel is visible.
    pub show_addon_manager: bool,

    /// Whether the right-side Genetics Workbench panel is
    /// visible. Defaults to `false` (the CAD-side Mesh Toolbox is the
    /// default right panel); flipped on from the View menu / command
    /// palette. The two right-side workbenches are independent — both
    /// can be open at once, egui docks them side by side.
    pub show_genetics_workbench: bool,
    /// Form + result state for the thirteen Genetics-workbench panels
    /// (one per computational-biology crate). See
    /// [`crate::genetics_workbench`].
    pub genetics: crate::genetics_workbench::GeneticsWorkbenchState,

    /// Whether the right-side Aerodynamics / Wind
    /// Tunnel workbench panel is visible. Defaults to `false`; flipped
    /// on from the View menu. Independent of the Mesh Toolbox and the
    /// Genetics workbench — egui docks them side by side.
    pub show_aero_workbench: bool,
    /// Form + result state for the Wind Tunnel workbench — the eight
    /// workflow sections wrapping the `valenx-aero` CFD engine. See
    /// [`crate::aero_workbench`].
    pub aero: crate::aero_workbench::AeroWorkbenchState,
    /// The aero flow-visualization field overlay, if one is active.
    /// When `Some`, the viewport colours the loaded mesh by this scalar
    /// field through the per-vertex colour ramp — it takes priority
    /// over the post-run results overlay. Pushed by the Wind Tunnel
    /// workbench's "Show field in 3-D viewport" button.
    pub aero_field_overlay: Option<valenx_fields::Field>,

    /// Whether the right-side FEM Workbench panel is visible. Defaults
    /// to `false`; flipped on from the View menu. Independent of the
    /// other workbenches — egui docks them side by side.
    pub show_fem_workbench: bool,
    /// Form + result state for the FEM Workbench — native linear-static
    /// and modal finite-element analysis wrapping the `valenx-fem`
    /// in-process solvers (no external solver, no input deck). See
    /// [`crate::fem_workbench`].
    pub fem: crate::fem_workbench::FemWorkbenchState,

    /// Whether the right-side Black-Hole / Relativity workbench is visible.
    /// Defaults to `false`; flipped on from the View menu. Independent of the
    /// other workbenches — egui docks them side by side.
    pub show_blackhole_workbench: bool,
    /// Form + result state for the Black-Hole / Relativity workbench — native
    /// general relativity (Kerr–Newman observables, thermodynamics, shadow
    /// ray-tracer) over the in-process `valenx-relativity` engine. See
    /// [`crate::blackhole_workbench`].
    pub blackhole: crate::blackhole_workbench::BlackHoleWorkbenchState,

    /// Whether the right-side Rotor / Drone (BEMT) workbench is visible.
    /// Defaults to `false`; flipped on from the View menu or opened by the
    /// agent bridge under the id `"rotor"`. Independent of the other
    /// workbenches — egui docks them side by side.
    pub show_rotor_workbench: bool,
    /// Form + result state for the Rotor / Drone (BEMT) workbench — native
    /// propeller / rotor blade-element-momentum-theory performance over the
    /// in-process `valenx-rotor` engine. See [`crate::rotor_workbench`].
    pub rotor: crate::rotor_workbench::RotorWorkbenchState,

    /// Whether the right-side Induction Motor workbench is visible. Defaults
    /// to `false`; flipped on from the View menu.
    pub show_inductionmotor_workbench: bool,
    /// State for the Induction Motor workbench, wrapping
    /// `valenx-inductionmotor`. See [`crate::inductionmotor_workbench`].
    pub inductionmotor: crate::inductionmotor_workbench::InductionMotorWorkbenchState,

    /// Whether the right-side CFD Workbench panel is visible. Defaults
    /// to `false`; flipped on from the View menu. Independent of the
    /// other workbenches — egui docks them side by side.
    pub show_cfd_workbench: bool,
    /// Form + result state for the CFD Workbench — native 2-D
    /// incompressible laminar CFD (SIMPLE) wrapping `valenx-cfd-native`.
    /// See [`crate::cfd_workbench`].
    pub cfd: crate::cfd_workbench::CfdWorkbenchState,

    /// Whether the right-side Reaction Dynamics workbench is visible.
    /// Defaults to `false`; flipped on from the View menu. Independent of
    /// the other workbenches — egui docks them side by side.
    pub show_reactdyn_workbench: bool,
    /// Form + result state for the Reaction Dynamics workbench — native
    /// ab-initio MD (AIMD) wrapping `valenx-reactdyn`. See
    /// [`crate::reactdyn_workbench`].
    pub reactdyn: crate::reactdyn_workbench::ReactdynWorkbenchState,

    /// Whether the right-side Springs Workbench is visible. Defaults to
    /// `false`; flipped on from the View menu. Independent of the other
    /// workbenches — egui docks them side by side.
    pub show_springs_workbench: bool,
    /// Form + result state for the Springs Workbench — native helical-spring
    /// design wrapping `valenx-springs`. See [`crate::springs_workbench`].
    pub springs: crate::springs_workbench::SpringsWorkbenchState,

    /// Whether the right-side Bearing workbench is visible. Defaults to
    /// `false`; flipped on from the View menu.
    pub show_bearing_workbench: bool,
    /// State for the Bearing workbench, wrapping `valenx-bearing`. See
    /// [`crate::bearing_workbench`].
    pub bearing: crate::bearing_workbench::BearingWorkbenchState,

    /// Whether the right-side Belt Drive workbench is visible. Defaults to
    /// `false`; flipped on from the View menu.
    pub show_beltdrive_workbench: bool,
    /// State for the Belt Drive workbench, wrapping `valenx-beltdrive`. See
    /// [`crate::beltdrive_workbench`].
    pub beltdrive: crate::beltdrive_workbench::BeltDriveWorkbenchState,

    /// Whether the right-side Buckling workbench is visible. Defaults to
    /// `false`; flipped on from the View menu.
    pub show_buckling_workbench: bool,
    /// State for the Buckling workbench, wrapping `valenx-buckling`. See
    /// [`crate::buckling_workbench`].
    pub buckling: crate::buckling_workbench::BucklingWorkbenchState,

    /// Whether the right-side Brake workbench is visible. Defaults to
    /// `false`; flipped on from the View menu.
    pub show_brake_workbench: bool,
    /// State for the Brake workbench, wrapping `valenx-brake`. See
    /// [`crate::brake_workbench`].
    pub brake: crate::brake_workbench::BrakeWorkbenchState,

    /// Whether the right-side Fatigue workbench is visible. Defaults to
    /// `false`; flipped on from the View menu.
    pub show_fatigue_workbench: bool,
    /// State for the Fatigue workbench, wrapping `valenx-fatigue`. See
    /// [`crate::fatigue_workbench`].
    pub fatigue: crate::fatigue_workbench::FatigueWorkbenchState,

    /// Whether the right-side Gear Tooth workbench is visible. Defaults to
    /// `false`; flipped on from the View menu.
    pub show_geartooth_workbench: bool,
    /// State for the Gear Tooth workbench, wrapping `valenx-geartooth`. See
    /// [`crate::geartooth_workbench`].
    pub geartooth: crate::geartooth_workbench::GeartoothWorkbenchState,

    /// Whether the right-side Pharmacokinetics workbench is visible. Defaults
    /// to `false`; flipped on from the View menu.
    pub show_pharmacokinetics_workbench: bool,
    /// State for the Pharmacokinetics workbench, wrapping
    /// `valenx-pharmacokinetics`. See [`crate::pharmacokinetics_workbench`].
    pub pharmacokinetics: crate::pharmacokinetics_workbench::PharmacokineticsWorkbenchState,

    /// Whether the right-side Pipe Network workbench is visible. Defaults to
    /// `false`; flipped on from the View menu.
    pub show_pipenetwork_workbench: bool,
    /// State for the Pipe Network workbench, wrapping `valenx-pipenetwork`.
    /// See [`crate::pipenetwork_workbench`].
    pub pipenetwork: crate::pipenetwork_workbench::PipeNetworkWorkbenchState,

    /// Whether the right-side RC Beam workbench is visible. Defaults to
    /// `false`; flipped on from the View menu.
    pub show_rcbeam_workbench: bool,
    /// State for the RC Beam workbench, wrapping `valenx-rcbeam`. See
    /// [`crate::rcbeam_workbench`].
    pub rcbeam: crate::rcbeam_workbench::RcBeamWorkbenchState,

    /// Whether the right-side Marine / Hull Workbench is visible. Off by
    /// default; toggled from the View menu.
    pub show_marine_workbench: bool,
    /// Form + result state for the Marine / Hull Workbench — native
    /// box-form hull hydrostatics wrapping `valenx-marine`. See
    /// [`crate::marine_workbench`].
    pub marine: crate::marine_workbench::MarineWorkbenchState,

    /// Whether the right-side Capacitor workbench is visible. Defaults to
    /// `false`; flipped on from the View menu.
    pub show_capacitor_workbench: bool,
    /// State for the Capacitor workbench, wrapping `valenx-capacitor`. See
    /// [`crate::capacitor_workbench`].
    pub capacitor: crate::capacitor_workbench::CapacitorWorkbenchState,

    /// Whether the right-side Fan Laws workbench is visible. Defaults to
    /// `false`; flipped on from the View menu.
    pub show_fanlaws_workbench: bool,
    /// State for the Fan Laws workbench, wrapping `valenx-fanlaws`. See
    /// [`crate::fanlaws_workbench`].
    pub fanlaws: crate::fanlaws_workbench::FanLawsWorkbenchState,

    /// Whether the right-side Creep workbench is visible. Defaults to
    /// `false`; flipped on from the View menu.
    pub show_creep_workbench: bool,
    /// State for the Creep workbench, wrapping `valenx-creep`. See
    /// [`crate::creep_workbench`].
    pub creep: crate::creep_workbench::CreepWorkbenchState,

    /// Whether the right-side Electrochemistry workbench is visible. Defaults
    /// to `false`; flipped on from the View menu.
    pub show_electrochem_workbench: bool,
    /// State for the Electrochemistry workbench, wrapping `valenx-electrochem`.
    /// See [`crate::electrochem_workbench`].
    pub electrochem: crate::electrochem_workbench::ElectrochemWorkbenchState,

    /// Whether the right-side Enzyme Kinetics workbench is visible. Defaults
    /// to `false`; flipped on from the View menu.
    pub show_enzymekinetics_workbench: bool,
    /// State for the Enzyme Kinetics workbench, wrapping
    /// `valenx-enzymekinetics`. See [`crate::enzymekinetics_workbench`].
    pub enzymekinetics: crate::enzymekinetics_workbench::EnzymeKineticsWorkbenchState,

    /// Whether the right-side Gears Workbench is visible. Defaults to
    /// `false`; flipped on from the View menu. Independent of the other
    /// workbenches — egui docks them side by side.
    pub show_gears_workbench: bool,
    /// Form + result state for the Gears Workbench — native involute-gear
    /// design wrapping `valenx-gears`. See [`crate::gears_workbench`].
    pub gears: crate::gears_workbench::GearsWorkbenchState,

    /// Whether the right-side Pneumatics workbench is visible (View menu).
    pub show_pneumatics_workbench: bool,
    /// State for the Pneumatics workbench. See [`crate::pneumatics_workbench`].
    pub pneumatics: crate::pneumatics_workbench::PneumaticsWorkbenchState,

    /// Whether the right-side Psychrometrics workbench is visible (View menu).
    pub show_psychrometrics_workbench: bool,
    /// State for the Psychrometrics workbench. See [`crate::psychrometrics_workbench`].
    pub psychrometrics: crate::psychrometrics_workbench::PsychrometricsWorkbenchState,

    /// Whether the right-side Thermistor workbench is visible (View menu).
    pub show_thermistor_workbench: bool,
    /// State for the Thermistor workbench. See [`crate::thermistor_workbench`].
    pub thermistor: crate::thermistor_workbench::ThermistorWorkbenchState,

    /// Whether the right-side Strain Gauge workbench is visible (View menu).
    pub show_straingauge_workbench: bool,
    /// State for the Strain Gauge workbench. See [`crate::straingauge_workbench`].
    pub straingauge: crate::straingauge_workbench::StrainGaugeWorkbenchState,

    /// Whether the right-side Drone Workbench is visible. Defaults to
    /// `false`; flipped on from the View menu. Independent of the other
    /// workbenches — egui docks them side by side.
    pub show_drone_workbench: bool,
    /// Form + result state for the Drone Workbench — native multirotor
    /// hover performance wrapping `valenx-drone`. See [`crate::drone_workbench`].
    pub drone: crate::drone_workbench::DroneWorkbenchState,

    /// Whether the right-side Acoustics workbench is visible. Defaults to
    /// `false`; flipped on from the View menu.
    pub show_acoustics_workbench: bool,
    /// State for the Acoustics workbench, wrapping `valenx-acoustics`. See
    /// [`crate::acoustics_workbench`].
    pub acoustics: crate::acoustics_workbench::AcousticsWorkbenchState,

    /// Whether the right-side Acid-Base workbench is visible. Defaults to
    /// `false`; flipped on from the View menu.
    pub show_acidbase_workbench: bool,
    /// State for the Acid-Base workbench, wrapping `valenx-acidbase`. See
    /// [`crate::acidbase_workbench`].
    pub acidbase: crate::acidbase_workbench::AcidBaseWorkbenchState,

    /// Whether the right-side BJT workbench is visible. Defaults to `false`;
    /// flipped on from the View menu.
    pub show_bjt_workbench: bool,
    /// State for the BJT workbench, wrapping `valenx-bjt`. See
    /// [`crate::bjt_workbench`].
    pub bjt: crate::bjt_workbench::BjtWorkbenchState,

    /// Whether the right-side BMR / TDEE workbench is visible. Defaults to
    /// `false`; flipped on from the View menu.
    pub show_bmr_workbench: bool,
    /// State for the BMR / TDEE workbench, wrapping `valenx-bmr`. See
    /// [`crate::bmr_workbench`].
    pub bmr: crate::bmr_workbench::BmrWorkbenchState,

    /// Whether the right-side Bolted Joint workbench is visible. Defaults to
    /// `false`; flipped on from the View menu.
    pub show_bolt_workbench: bool,
    /// State for the Bolted Joint workbench, wrapping `valenx-bolt`. See
    /// [`crate::bolt_workbench`].
    pub bolt: crate::bolt_workbench::BoltWorkbenchState,

    /// Whether the right-side **Parametric Sketch (constraints)** panel is
    /// visible. Defaults to `false`; flipped on from
    /// **Part Design → "Parametric Sketch (constraints)"**. This panel is a
    /// first-class, discoverable host for the in-house `valenx-sketch`
    /// constraint sketcher — it shares its sketch state with the Mesh
    /// Toolbox's Sketcher section (`mesh_toolbox.sketcher`), so there is no
    /// separate state struct. See [`crate::param_sketch_panel`].
    pub show_param_sketch: bool,

    /// Whether the right-side Geomatics Workbench is visible. Defaults to
    /// `false`; flipped on from the View menu. Independent of the other
    /// workbenches — egui docks them side by side.
    pub show_geomatics_workbench: bool,
    /// Form + result state for the Geomatics Workbench — native geodesic
    /// calculations wrapping `valenx-geomatics`. See
    /// [`crate::geomatics_workbench`].
    pub geomatics: crate::geomatics_workbench::GeomaticsWorkbenchState,

    /// Whether the right-side Op-Amp workbench is visible (View menu). Off by default.
    pub show_opamp_workbench: bool,
    /// Form + result state for the Op-Amp workbench — ideal closed-loop gain /
    /// bandwidth on `valenx-opamp`. See [`crate::opamp_workbench`].
    pub opamp: crate::opamp_workbench::OpAmpWorkbenchState,
    /// Whether the right-side LED workbench is visible (View menu). Off by default.
    pub show_led_workbench: bool,
    /// Form + result state for the LED workbench — series resistor sizing on
    /// `valenx-led`. See [`crate::led_workbench`].
    pub led: crate::led_workbench::LedWorkbenchState,
    /// Whether the right-side Thermocouple workbench is visible (View menu). Off by default.
    pub show_thermocouple_workbench: bool,
    /// Form + result state for the Thermocouple workbench — Seebeck EMF on
    /// `valenx-thermocouple`. See [`crate::thermocouple_workbench`].
    pub thermocouple: crate::thermocouple_workbench::ThermocoupleWorkbenchState,
    /// Whether the right-side Transmission Line workbench is visible (View menu). Off by default.
    pub show_transmissionline_workbench: bool,
    /// Form + result state for the Transmission Line workbench — reflection / VSWR
    /// on `valenx-transmissionline`. See [`crate::transmissionline_workbench`].
    pub transmissionline: crate::transmissionline_workbench::TransmissionLineWorkbenchState,
    /// Whether the right-side Power Factor workbench is visible (View menu). Off by default.
    pub show_powerfactor_workbench: bool,
    /// Form + result state for the Power Factor workbench — AC power triangle +
    /// correction on `valenx-powerfactor`. See [`crate::powerfactor_workbench`].
    pub powerfactor: crate::powerfactor_workbench::PowerFactorWorkbenchState,
    /// Whether the right-side Resistor Network workbench is visible (View menu). Off by default.
    pub show_resistornetwork_workbench: bool,
    /// Form + result state for the Resistor Network workbench — series / parallel /
    /// divider on `valenx-resistor-network`. See [`crate::resistornetwork_workbench`].
    pub resistornetwork: crate::resistornetwork_workbench::ResistorNetworkWorkbenchState,
    /// Whether the right-side Rectifier workbench is visible (View menu). Off by default.
    pub show_rectifier_workbench: bool,
    /// Form + result state for the Rectifier workbench — rectifier figures + ripple
    /// on `valenx-rectifier`. See [`crate::rectifier_workbench`].
    pub rectifier: crate::rectifier_workbench::RectifierWorkbenchState,
    /// Whether the right-side Filter workbench is visible (View menu). Off by default.
    pub show_filter_workbench: bool,
    /// Form + result state for the Filter workbench — RC / RLC response on
    /// `valenx-filter`. See [`crate::filter_workbench`].
    pub filter: crate::filter_workbench::FilterWorkbenchState,

    /// Whether the right-side Heat Transfer workbench is visible. Defaults
    /// to `false`; flipped on from the View menu.
    pub show_heattransfer_workbench: bool,
    /// State for the Heat Transfer workbench, wrapping `valenx-heat-transfer`.
    /// See [`crate::heattransfer_workbench`].
    pub heattransfer: crate::heattransfer_workbench::HeatTransferWorkbenchState,

    /// Whether the right-side Four-Bar Linkage Workbench is visible.
    /// Defaults to `false`; flipped on from the View menu. Independent of the
    /// other workbenches — egui docks them side by side.
    pub show_fourbar_workbench: bool,
    /// Form + result state for the Four-Bar Linkage Workbench — native planar
    /// four-bar mechanism kinematics wrapping `valenx-kinematics`. See
    /// [`crate::fourbar_workbench`].
    pub fourbar: crate::fourbar_workbench::FourBarWorkbenchState,

    /// Whether the right-side Shaft Design workbench is visible (View menu). Off by default.
    pub show_shaftdesign_workbench: bool,
    /// State for the Shaft Design workbench — combined bending + torsion shaft
    /// sizing on `valenx-shaftdesign`. See [`crate::shaftdesign_workbench`].
    pub shaftdesign: crate::shaftdesign_workbench::ShaftDesignWorkbenchState,
    /// Whether the right-side Power Screw workbench is visible (View menu). Off by default.
    pub show_screwthread_workbench: bool,
    /// State for the Power Screw workbench — square-thread lead-screw torque on
    /// `valenx-screwthread`. See [`crate::screwthread_workbench`].
    pub screwthread: crate::screwthread_workbench::ScrewThreadWorkbenchState,
    /// Whether the right-side Pulley System workbench is visible (View menu). Off by default.
    pub show_pulley_workbench: bool,
    /// State for the Pulley System workbench — block-and-tackle mechanical
    /// advantage on `valenx-pulley`. See [`crate::pulley_workbench`].
    pub pulley: crate::pulley_workbench::PulleyWorkbenchState,
    /// Whether the right-side Spring Design workbench is visible (View menu). Off by default.
    pub show_springdesign_workbench: bool,
    /// State for the Spring Design workbench — helical compression spring on
    /// `valenx-spring-design`. See [`crate::springdesign_workbench`].
    pub springdesign: crate::springdesign_workbench::SpringDesignWorkbenchState,
    /// Whether the right-side Spring Combination workbench is visible (View menu). Off by default.
    pub show_springcombination_workbench: bool,
    /// State for the Spring Combination workbench — series / parallel spring
    /// networks on `valenx-springcombination`. See [`crate::springcombination_workbench`].
    pub springcombination: crate::springcombination_workbench::SpringCombinationWorkbenchState,
    /// Whether the right-side Vibration workbench is visible (View menu). Off by default.
    pub show_vibration_workbench: bool,
    /// State for the Vibration workbench — single-DOF forced-vibration response
    /// on `valenx-vibration`. See [`crate::vibration_workbench`].
    pub vibration: crate::vibration_workbench::VibrationWorkbenchState,
    /// Whether the right-side Riveted Joint workbench is visible (View menu). Off by default.
    pub show_rivet_workbench: bool,
    /// State for the Riveted Joint workbench — rivet-joint strength + failure
    /// mode on `valenx-rivet`. See [`crate::rivet_workbench`].
    pub rivet: crate::rivet_workbench::RivetWorkbenchState,
    /// Whether the right-side Soil Bearing workbench is visible (View menu). Off by default.
    pub show_soilbearing_workbench: bool,
    /// State for the Soil Bearing workbench — Terzaghi strip-footing bearing
    /// capacity on `valenx-soilbearing`. See [`crate::soilbearing_workbench`].
    pub soilbearing: crate::soilbearing_workbench::SoilBearingWorkbenchState,

    /// Whether the right-side Piping Workbench is visible. Defaults to
    /// `false`; flipped on from the View menu. Independent of the other
    /// workbenches — egui docks them side by side.
    pub show_piping_workbench: bool,
    /// Form + result state for the Piping Workbench — native pipe-section
    /// sizing wrapping `valenx-piping`. See [`crate::piping_workbench`].
    pub piping: crate::piping_workbench::PipingWorkbenchState,

    /// Whether the right-side Retaining Wall workbench is visible (View menu). Off by default.
    pub show_retainingwall_workbench: bool,
    /// State for the Retaining Wall workbench — Rankine earth pressure on
    /// `valenx-retainingwall`. See [`crate::retainingwall_workbench`].
    pub retainingwall: crate::retainingwall_workbench::RetainingWallWorkbenchState,
    /// Whether the right-side Open Channel workbench is visible (View menu). Off by default.
    pub show_openchannel_workbench: bool,
    /// State for the Open Channel workbench — Manning flow + Froude on
    /// `valenx-openchannel`. See [`crate::openchannel_workbench`].
    pub openchannel: crate::openchannel_workbench::OpenChannelWorkbenchState,
    /// Whether the right-side Weir Flow workbench is visible (View menu). Off by default.
    pub show_weir_workbench: bool,
    /// State for the Weir Flow workbench — sharp-crested weir discharge on
    /// `valenx-weir`. See [`crate::weir_workbench`].
    pub weir: crate::weir_workbench::WeirWorkbenchState,
    /// Whether the right-side Thermodynamic Cycle workbench is visible (View menu). Off by default.
    pub show_thermocycle_workbench: bool,
    /// State for the Thermodynamic Cycle workbench — ideal-cycle efficiency on
    /// `valenx-thermocycle`. See [`crate::thermocycle_workbench`].
    pub thermocycle: crate::thermocycle_workbench::ThermoCycleWorkbenchState,
    /// Whether the right-side Queueing workbench is visible (View menu). Off by default.
    pub show_queueing_workbench: bool,
    /// State for the Queueing workbench — M/M/1 steady-state metrics on
    /// `valenx-queueing`. See [`crate::queueing_workbench`].
    pub queueing: crate::queueing_workbench::QueueingWorkbenchState,
    /// Whether the right-side Radioactive Decay workbench is visible (View menu). Off by default.
    pub show_radioactivity_workbench: bool,
    /// State for the Radioactive Decay workbench — single-nuclide decay on
    /// `valenx-radioactivity`. See [`crate::radioactivity_workbench`].
    pub radioactivity: crate::radioactivity_workbench::RadioactivityWorkbenchState,
    /// Whether the right-side Osmosis workbench is visible (View menu). Off by default.
    pub show_osmosis_workbench: bool,
    /// State for the Osmosis workbench — van't Hoff + Starling on
    /// `valenx-osmosis`. See [`crate::osmosis_workbench`].
    pub osmosis: crate::osmosis_workbench::OsmosisWorkbenchState,
    /// Whether the right-side Thermoregulation workbench is visible (View menu). Off by default.
    pub show_thermoreg_workbench: bool,
    /// State for the Thermoregulation workbench — single-node heat balance on
    /// `valenx-thermoreg`. See [`crate::thermoreg_workbench`].
    pub thermoreg: crate::thermoreg_workbench::ThermoRegWorkbenchState,
    /// Whether the right-side Hemodynamics workbench is visible (View menu). Off by default.
    pub show_hemodynamics_workbench: bool,
    /// State for the Hemodynamics workbench — cardiac output / Poiseuille /
    /// Windkessel on `valenx-hemodynamics`. See [`crate::hemodynamics_workbench`].
    pub hemodynamics: crate::hemodynamics_workbench::HemodynamicsWorkbenchState,
    /// Whether the right-side Population Dynamics workbench is visible (View menu). Off by default.
    pub show_popdynamics_workbench: bool,
    /// State for the Population Dynamics workbench — SIR / logistic /
    /// Lotka-Volterra on `valenx-popdynamics`. See [`crate::popdynamics_workbench`].
    pub popdynamics: crate::popdynamics_workbench::PopDynamicsWorkbenchState,

    /// Whether the right-side Rail / Train Workbench is visible. Defaults to
    /// `false`; flipped on from the View menu. Independent of the other
    /// workbenches — egui docks them side by side.
    pub show_rail_workbench: bool,
    /// Form + result state for the Rail / Train Workbench — native train
    /// resistance + tractive effort wrapping `valenx-rail`. See
    /// [`crate::rail_workbench`].
    pub rail: crate::rail_workbench::RailWorkbenchState,

    /// Whether the right-side Bone Mechanics workbench is visible. Defaults to
    /// `false`; flipped on from the View menu.
    pub show_bonemech_workbench: bool,
    /// State for the Bone Mechanics workbench, wrapping `valenx-bonemech`. See
    /// [`crate::bonemech_workbench`].
    pub bonemech: crate::bonemech_workbench::BonemechWorkbenchState,

    /// Whether the right-side Chain Drive workbench is visible. Defaults to
    /// `false`; flipped on from the View menu.
    pub show_chaindrive_workbench: bool,
    /// State for the Chain Drive workbench, wrapping `valenx-chaindrive`. See
    /// [`crate::chaindrive_workbench`].
    pub chaindrive: crate::chaindrive_workbench::ChainDriveWorkbenchState,

    /// Whether the right-side Clutch workbench is visible. Defaults to `false`;
    /// flipped on from the View menu.
    pub show_clutch_workbench: bool,
    /// State for the Clutch workbench, wrapping `valenx-clutch`. See
    /// [`crate::clutch_workbench`].
    pub clutch: crate::clutch_workbench::ClutchWorkbenchState,

    /// Whether the right-side Solenoid Coil workbench is visible. Defaults to
    /// `false`; flipped on from the View menu.
    pub show_coil_workbench: bool,
    /// State for the Solenoid Coil workbench, wrapping `valenx-coil`. See
    /// [`crate::coil_workbench`].
    pub coil: crate::coil_workbench::CoilWorkbenchState,

    /// Whether the right-side Steel Column workbench is visible. Defaults to
    /// `false`; flipped on from the View menu.
    pub show_columnsteel_workbench: bool,
    /// State for the Steel Column workbench, wrapping `valenx-columnsteel`. See
    /// [`crate::columnsteel_workbench`].
    pub columnsteel: crate::columnsteel_workbench::ColumnSteelWorkbenchState,

    /// Whether the right-side Collision Workbench is visible. Defaults to
    /// `false`; flipped on from the View menu. Independent of the other
    /// workbenches — egui docks them side by side.
    pub show_collision_workbench: bool,
    /// Form + result state for the Collision Workbench — native AABB
    /// geometry + overlap tests wrapping `valenx-collision`. See
    /// [`crate::collision_workbench`].
    pub collision: crate::collision_workbench::CollisionWorkbenchState,

    /// Whether the right-side Statics workbench is visible (View menu).
    pub show_statics_workbench: bool,
    /// State for the Statics workbench. See [`crate::statics_workbench`].
    pub statics: crate::statics_workbench::StaticsWorkbenchState,

    /// Whether the right-side Projectile workbench is visible (View menu).
    pub show_projectile_workbench: bool,
    /// State for the Projectile workbench. See [`crate::projectile_workbench`].
    pub projectile: crate::projectile_workbench::ProjectileWorkbenchState,

    /// Whether the right-side Conveyor workbench is visible (View menu).
    pub show_conveyor_workbench: bool,
    /// State for the Conveyor workbench. See [`crate::conveyor_workbench`].
    pub conveyor: crate::conveyor_workbench::ConveyorWorkbenchState,

    /// Whether the right-side Fluid Statics workbench is visible (View menu).
    pub show_fluidstatics_workbench: bool,
    /// State for the Fluid Statics workbench. See [`crate::fluidstatics_workbench`].
    pub fluidstatics: crate::fluidstatics_workbench::FluidStaticsWorkbenchState,

    /// Whether the right-side Plate Bending workbench is visible (View menu).
    pub show_plate_workbench: bool,
    /// State for the Plate Bending workbench. See [`crate::plate_workbench`].
    pub plate: crate::plate_workbench::PlateWorkbenchState,

    /// Whether the right-side Strain Rosette workbench is visible (View menu).
    pub show_strainrosette_workbench: bool,
    /// State for the Strain Rosette workbench. See [`crate::strainrosette_workbench`].
    pub strainrosette: crate::strainrosette_workbench::StrainRosetteWorkbenchState,

    /// Whether the right-side Transformer workbench is visible (View menu).
    pub show_transformer_workbench: bool,
    /// State for the Transformer workbench. See [`crate::transformer_workbench`].
    pub transformer: crate::transformer_workbench::TransformerWorkbenchState,

    /// Whether the right-side Three-Phase workbench is visible (View menu).
    pub show_threephase_workbench: bool,
    /// State for the Three-Phase workbench. See [`crate::threephase_workbench`].
    pub threephase: crate::threephase_workbench::ThreePhaseWorkbenchState,

    /// Whether the right-side Solar PV Workbench is visible. Defaults to
    /// `false`; flipped on from the View menu. Independent of the other
    /// workbenches — egui docks them side by side.
    pub show_solarpv_workbench: bool,
    /// Form + result state for the Solar PV Workbench — native single-diode
    /// photovoltaic cell performance wrapping `valenx-solarpv`. See
    /// [`crate::solarpv_workbench`].
    pub solarpv: crate::solarpv_workbench::SolarPvWorkbenchState,

    /// Whether the right-side Sheet Metal Workbench is visible. Defaults to
    /// `false`; flipped on from the View menu. Independent of the other
    /// workbenches — egui docks them side by side.
    pub show_sheetmetal_workbench: bool,
    /// Form + result state for the Sheet Metal Workbench — native bend
    /// allowance / deduction wrapping `valenx-sheet-metal`. See
    /// [`crate::sheetmetal_workbench`].
    pub sheetmetal: crate::sheetmetal_workbench::SheetmetalWorkbenchState,

    /// Whether the right-side Truss Workbench is visible. Defaults to
    /// `false`; flipped on from the View menu. Independent of the other
    /// workbenches — egui docks them side by side.
    pub show_truss_workbench: bool,
    /// Form + result state for the Truss Workbench — native planar
    /// pin-jointed truss analysis wrapping `valenx-truss`. See
    /// [`crate::truss_workbench`].
    pub truss: crate::truss_workbench::TrussWorkbenchState,

    /// Whether the right-side Field Statistics Workbench is visible. Defaults
    /// to `false`; flipped on from the View menu. Independent of the other
    /// workbenches — egui docks them side by side.
    pub show_fields_workbench: bool,
    /// Form + result state for the Field Statistics Workbench — descriptive
    /// statistics over a pasted number list, via `valenx-fields`. See
    /// [`crate::fields_workbench`].
    pub fields: crate::fields_workbench::FieldsWorkbenchState,

    /// Whether the right-side Gearbox workbench is visible. Defaults to
    /// `false`; flipped on from the View menu.
    pub show_gearbox_workbench: bool,
    /// State for the Gearbox workbench, wrapping `valenx-gearbox`. See
    /// [`crate::gearbox_workbench`].
    pub gearbox: crate::gearbox_workbench::GearboxWorkbenchState,

    /// Whether the right-side Cam Dynamics workbench is visible (View menu). Off by default.
    pub show_camdynamics_workbench: bool,
    /// State for the Cam Dynamics workbench — cam-follower rise kinematics on
    /// `valenx-camdynamics`. See [`crate::camdynamics_workbench`].
    pub camdynamics: crate::camdynamics_workbench::CamDynamicsWorkbenchState,
    /// Whether the right-side Battery ECM workbench is visible (View menu). Off by default.
    pub show_batteryecm_workbench: bool,
    /// State for the Battery ECM workbench — first-order Thevenin terminal voltage
    /// on `valenx-battery-ecm`. See [`crate::batteryecm_workbench`].
    pub batteryecm: crate::batteryecm_workbench::BatteryEcmWorkbenchState,
    /// Whether the right-side Diffusion workbench is visible (View menu). Off by default.
    pub show_diffusion_workbench: bool,
    /// State for the Diffusion workbench — Fickian flux + Gaussian spread on
    /// `valenx-diffusion`. See [`crate::diffusion_workbench`].
    pub diffusion: crate::diffusion_workbench::DiffusionWorkbenchState,
    /// Whether the right-side Dimensionless Numbers workbench is visible (View menu). Off by default.
    pub show_dimensional_workbench: bool,
    /// State for the Dimensionless Numbers workbench — similitude groups +
    /// regime classifiers on `valenx-dimensional`. See [`crate::dimensional_workbench`].
    pub dimensional: crate::dimensional_workbench::DimensionalWorkbenchState,
    /// Whether the right-side FFT / Spectrum workbench is visible (View menu). Off by default.
    pub show_fft_workbench: bool,
    /// State for the FFT / Spectrum workbench — DFT of a synthesized tone on
    /// `valenx-fft`. See [`crate::fft_workbench`].
    pub fft: crate::fft_workbench::FftWorkbenchState,

    /// Whether the right-side Fasteners Workbench is visible. Defaults to
    /// `false`; flipped on from the View menu. Independent of the other
    /// workbenches — egui docks them side by side.
    pub show_fasteners_workbench: bool,
    /// Form + result state for the Fasteners Workbench — ISO 4017 hex-bolt
    /// dimensions wrapping `valenx-fasteners`. See
    /// [`crate::fasteners_workbench`].
    pub fasteners: crate::fasteners_workbench::FastenersWorkbenchState,

    /// Whether the right-side Fixed-Wing / Aircraft Workbench is visible.
    /// Defaults to `false`; flipped on from the View menu. Independent of the
    /// other workbenches — egui docks them side by side.
    pub show_fixedwing_workbench: bool,
    /// Form + result state for the Fixed-Wing / Aircraft Workbench — native
    /// preliminary aircraft point-performance wrapping `valenx-fixedwing`.
    /// See [`crate::fixedwing_workbench`].
    pub fixedwing: crate::fixedwing_workbench::FixedWingWorkbenchState,

    /// Whether the right-side Combustion workbench is visible (View menu).
    pub show_combustion_workbench: bool,
    /// State for the Combustion workbench. See [`crate::combustion_workbench`].
    pub combustion: crate::combustion_workbench::CombustionWorkbenchState,

    /// Whether the right-side Flywheel workbench is visible (View menu).
    pub show_flywheel_workbench: bool,
    /// State for the Flywheel workbench. See [`crate::flywheel_workbench`].
    pub flywheel: crate::flywheel_workbench::FlywheelWorkbenchState,

    /// Whether the right-side Fracture Mechanics workbench is visible (View menu).
    pub show_fracture_workbench: bool,
    /// State for the Fracture workbench. See [`crate::fracture_workbench`].
    pub fracture: crate::fracture_workbench::FractureWorkbenchState,

    /// Whether the right-side Hydraulics workbench is visible (View menu).
    pub show_hydraulics_workbench: bool,
    /// State for the Hydraulics workbench. See [`crate::hydraulics_workbench`].
    pub hydraulics: crate::hydraulics_workbench::HydraulicsWorkbenchState,

    /// Whether the right-side Inclined Plane workbench is visible (View menu).
    pub show_inclinedplane_workbench: bool,
    /// State for the Inclined Plane workbench. See [`crate::inclinedplane_workbench`].
    pub inclinedplane: crate::inclinedplane_workbench::InclinedPlaneWorkbenchState,

    /// Whether the right-side Insulation workbench is visible (View menu).
    pub show_insulation_workbench: bool,
    /// State for the Insulation workbench. See [`crate::insulation_workbench`].
    pub insulation: crate::insulation_workbench::InsulationWorkbenchState,

    /// Whether the right-side Lead Screw workbench is visible (View menu).
    pub show_leadscrew_workbench: bool,
    /// State for the Lead Screw workbench. See [`crate::leadscrew_workbench`].
    pub leadscrew: crate::leadscrew_workbench::LeadscrewWorkbenchState,

    /// Whether the right-side Lever workbench is visible (View menu).
    pub show_leverage_workbench: bool,
    /// State for the Lever workbench. See [`crate::leverage_workbench`].
    pub leverage: crate::leverage_workbench::LeverageWorkbenchState,

    /// Whether the right-side Mohr's Circle workbench is visible (View menu).
    pub show_mohr_workbench: bool,
    /// State for the Mohr's Circle workbench. See [`crate::mohr_workbench`].
    pub mohr: crate::mohr_workbench::MohrWorkbenchState,

    /// Whether the right-side MOSFET workbench is visible (View menu).
    pub show_mosfet_workbench: bool,
    /// State for the MOSFET workbench. See [`crate::mosfet_workbench`].
    pub mosfet: crate::mosfet_workbench::MosfetWorkbenchState,

    /// Whether the right-side Optics workbench is visible (View menu).
    pub show_optics_workbench: bool,
    /// State for the Optics workbench. See [`crate::optics_workbench`].
    pub optics: crate::optics_workbench::OpticsWorkbenchState,

    /// Whether the right-side Orifice Meter workbench is visible (View menu).
    pub show_orifice_workbench: bool,
    /// State for the Orifice Meter workbench. See [`crate::orifice_workbench`].
    pub orifice: crate::orifice_workbench::OrificeWorkbenchState,

    /// Whether the right-side Pressure Vessel workbench is visible (View menu).
    pub show_pressurevessel_workbench: bool,
    /// State for the Pressure Vessel workbench. See [`crate::pressurevessel_workbench`].
    pub pressurevessel: crate::pressurevessel_workbench::PressureVesselWorkbenchState,

    /// Whether the right-side Torsion workbench is visible (View menu).
    pub show_torsion_workbench: bool,
    /// State for the Torsion workbench. See [`crate::torsion_workbench`].
    pub torsion: crate::torsion_workbench::TorsionWorkbenchState,

    /// Whether the right-side Refrigeration workbench is visible (View menu).
    pub show_refrigeration_workbench: bool,
    /// State for the Refrigeration workbench. See [`crate::refrigeration_workbench`].
    pub refrigeration: crate::refrigeration_workbench::RefrigerationWorkbenchState,

    /// Whether the right-side Frames Workbench is visible. Defaults to
    /// `false`; flipped on from the View menu. Independent of the other
    /// workbenches — egui docks them side by side.
    pub show_frames_workbench: bool,
    /// Form + result state for the Frames Workbench — structural
    /// cross-section properties wrapping `valenx-frames`. See
    /// [`crate::frames_workbench`].
    pub frames: crate::frames_workbench::FramesWorkbenchState,

    /// Whether the right-side DC Motor Workbench is visible. Defaults to
    /// `false`; flipped on from the View menu. Independent of the other
    /// workbenches — egui docks them side by side.
    pub show_dcmotor_workbench: bool,
    /// Form + result state for the DC Motor Workbench — native brushed-DC-
    /// motor performance wrapping `valenx-dcmotor`. See
    /// [`crate::dcmotor_workbench`].
    pub dcmotor: crate::dcmotor_workbench::DcMotorWorkbenchState,

    /// Whether the right-side Gas Dynamics workbench is visible. Defaults to
    /// `false`; flipped on from the View menu. Independent of the other
    /// workbenches — egui docks them side by side.
    pub show_gasdynamics_workbench: bool,
    /// Form + result state for the Gas Dynamics workbench — 1-D
    /// compressible-flow relations wrapping `valenx-gasdynamics`. See
    /// [`crate::gasdynamics_workbench`].
    pub gasdynamics: crate::gasdynamics_workbench::GasDynamicsWorkbenchState,

    /// Whether the right-side Thermal Expansion workbench is visible. Defaults
    /// to `false`; flipped on from the View menu. Independent of the other
    /// workbenches — egui docks them side by side.
    pub show_thermalexpansion_workbench: bool,
    /// Form + result state for the Thermal Expansion workbench — linear
    /// expansion + constrained stress wrapping `valenx-thermalexpansion`. See
    /// [`crate::thermalexpansion_workbench`].
    pub thermalexpansion: crate::thermalexpansion_workbench::ThermalExpansionWorkbenchState,

    /// Whether the right-side Neural-Interface (BCI stimulation) workbench is
    /// visible. Defaults to `false`; flipped on from the View menu.
    pub show_neuro_workbench: bool,
    /// Form + result state for the Neural-Interface workbench, wrapping
    /// `valenx-neuro`. See [`crate::neuro_workbench`].
    pub neuro: crate::neuro_workbench::NeuroWorkbenchState,

    /// Whether the right-side Wind Turbine workbench is visible. Defaults to
    /// `false`; flipped on from the View menu. Independent of the other
    /// workbenches — egui docks them side by side.
    pub show_windturbine_workbench: bool,
    /// Form + result state for the Wind Turbine workbench — native
    /// actuator-disc wind-turbine power wrapping `valenx-windturbine`. See
    /// [`crate::windturbine_workbench`].
    pub windturbine: crate::windturbine_workbench::WindTurbineWorkbenchState,

    /// Whether the right-side Parametric-CAD workbench is visible. Defaults to
    /// `false`; flipped on from the View menu.
    pub show_cad_workbench: bool,
    /// Form + result state for the Parametric-CAD workbench, wrapping
    /// `valenx-solvespace-3d`. See [`crate::cad_workbench`].
    pub cad: crate::cad_workbench::CadWorkbenchState,

    /// Whether the right-side Antenna workbench is visible. Defaults to
    /// `false`; flipped on from the View menu.
    pub show_antenna_workbench: bool,
    /// State for the Antenna workbench, wrapping `valenx-antenna`. See
    /// [`crate::antenna_workbench`].
    pub antenna: crate::antenna_workbench::AntennaWorkbenchState,

    /// Whether the right-side 2D Drafting workbench is visible. Defaults to
    /// `false`; flipped on from the View menu.
    pub show_draft2d_workbench: bool,
    /// State for the 2D Drafting workbench, wrapping `valenx-librecad-2d`. See
    /// [`crate::draft2d_workbench`].
    pub draft2d: crate::draft2d_workbench::Draft2dWorkbenchState,

    /// Whether the right-side Reinforcement workbench is visible. Defaults to
    /// `false`; flipped on from the View menu.
    pub show_reinforcement_workbench: bool,
    /// State for the Reinforcement workbench, wrapping `valenx-reinforcement`.
    /// See [`crate::reinforcement_workbench`].
    pub reinforcement: crate::reinforcement_workbench::ReinforcementWorkbenchState,

    /// Whether the right-side Path-Traced Render workbench is visible. Defaults
    /// to `false`; flipped on from the View menu.
    pub show_render_workbench: bool,
    /// State for the Render workbench, wrapping `valenx-pathtrace`. See
    /// [`crate::render_workbench`].
    pub render: crate::render_workbench::RenderWorkbenchState,

    /// Whether the right-side HVAC workbench is visible. Defaults to `false`;
    /// flipped on from the View menu.
    pub show_hvac_workbench: bool,
    /// State for the HVAC workbench, wrapping `valenx-hvac`. See
    /// [`crate::hvac_workbench`].
    pub hvac: crate::hvac_workbench::HvacWorkbenchState,

    /// Whether the right-side Beam Workbench is visible. Defaults to `false`;
    /// flipped on from the View menu.
    pub show_beam_workbench: bool,
    /// State for the Beam Workbench, wrapping `valenx-beam`. See
    /// [`crate::beam_workbench`].
    pub beam: crate::beam_workbench::BeamWorkbenchState,

    /// Whether the right-side Reverse-Engineering workbench is visible.
    /// Defaults to `false`; flipped on from the View menu.
    pub show_reverse_workbench: bool,
    /// State for the Reverse-Engineering workbench, wrapping `valenx-reverse`.
    /// See [`crate::reverse_workbench`].
    pub reverse: crate::reverse_workbench::ReverseWorkbenchState,

    /// Whether the right-side Pump workbench is visible. Defaults to `false`;
    /// flipped on from the View menu.
    pub show_pump_workbench: bool,
    /// State for the Pump workbench, wrapping `valenx-pump`. See
    /// [`crate::pump_workbench`].
    pub pump: crate::pump_workbench::PumpWorkbenchState,

    /// Whether the right-side Interior-Design workbench is visible. Defaults to
    /// `false`; flipped on from the View menu.
    pub show_interior_workbench: bool,
    /// State for the Interior-Design workbench, wrapping `valenx-interior`. See
    /// [`crate::interior_workbench`].
    pub interior: crate::interior_workbench::InteriorWorkbenchState,

    /// Whether the right-side Animation workbench is visible. Defaults to
    /// `false`; flipped on from the View menu.
    pub show_animate_workbench: bool,
    /// State for the Animation workbench, wrapping `valenx-animate`. See
    /// [`crate::animate_workbench`].
    pub animate: crate::animate_workbench::AnimateWorkbenchState,

    /// Whether the right-side Variant-Effect workbench is visible. Defaults to
    /// `false`; flipped on from the View menu.
    pub show_variant_effect_workbench: bool,
    /// State for the Variant-Effect workbench, wrapping `valenx-variant-effect`.
    /// See [`crate::variant_effect_workbench`].
    pub variant_effect: crate::variant_effect_workbench::VariantEffectWorkbenchState,

    /// Whether the right-side Heat Pump workbench is visible. Defaults to
    /// `false`; flipped on from the View menu.
    pub show_heatpump_workbench: bool,
    /// State for the Heat Pump workbench, wrapping `valenx-heatpump`. See
    /// [`crate::heatpump_workbench`].
    pub heatpump: crate::heatpump_workbench::HeatPumpWorkbenchState,

    /// Whether the right-side Astro / Launch workbench panel is visible.
    /// Defaults to `false`; flipped on from the View menu (Ctrl+4).
    /// Independent of the Mesh Toolbox, Genetics and Wind Tunnel
    /// workbenches — egui docks them side by side.
    pub show_astro_workbench: bool,
    /// Form + result state for the Astro / Launch workbench — the launch
    /// ascent simulator + the closed-form mission planners wrapping the
    /// `valenx-astro` crate. See [`crate::astro_workbench`].
    pub astro: crate::astro_workbench::AstroWorkbenchState,

    /// Whether the right-side Pipe Flow workbench is visible. Defaults to
    /// `false`; flipped on from the View menu.
    pub show_pipeflow_workbench: bool,
    /// State for the Pipe Flow workbench, wrapping `valenx-pipeflow`. See
    /// [`crate::pipeflow_workbench`].
    pub pipeflow: crate::pipeflow_workbench::PipeFlowWorkbenchState,

    /// Whether the right-side Rocket workbench panel is visible. Defaults
    /// to `false`; flipped on from the View menu. Surfaces the
    /// `valenx-rocket-demo` coupled design→simulate pipeline.
    pub show_rocket_workbench: bool,
    /// Form + result state for the Rocket workbench — the reactive
    /// design→simulate panel wrapping `valenx-rocket-demo`. See
    /// [`crate::rocket_workbench`].
    pub rocket: crate::rocket_workbench::RocketWorkbenchState,

    /// Whether the right-side Battery Pack workbench is visible. Defaults to
    /// `false`; flipped on from the View menu.
    pub show_batterypack_workbench: bool,
    /// State for the Battery Pack workbench, wrapping `valenx-batterypack`.
    /// See [`crate::batterypack_workbench`].
    pub batterypack: crate::batterypack_workbench::BatteryPackWorkbenchState,

    /// Whether the right-side Engine workbench panel is visible — the
    /// reactive engine design → analyze → optimize → export loop. On by
    /// default (set in [`ValenxApp::new`]).
    pub show_engine_workbench: bool,
    /// Form + result state for the Engine workbench. See
    /// [`crate::engine_workbench`].
    pub engine: crate::engine_workbench::EngineWorkbenchState,

    /// Whether the right-side Heat Exchanger workbench is visible. Defaults to
    /// `false`; flipped on from the View menu.
    pub show_heatexchanger_workbench: bool,
    /// State for the Heat Exchanger workbench, wrapping `valenx-heatexchanger`.
    /// See [`crate::heatexchanger_workbench`].
    pub heatexchanger: crate::heatexchanger_workbench::HeatExchangerWorkbenchState,

    /// Whether the right-side Car workbench panel is visible. Defaults to
    /// `false`; toggled from the View menu. Wraps `valenx-vehicle`'s
    /// performance model. See [`crate::car_workbench`].
    pub show_car_workbench: bool,
    /// Form + result state for the Car workbench (design → simulate over
    /// `valenx-vehicle`).
    pub car: crate::car_workbench::CarWorkbenchState,

    /// Whether the right-side Sensors workbench panel is visible. Defaults to
    /// `false`; toggled from the View menu or opened by the agent bridge under
    /// the id `"sensors"`. Wraps `valenx-sensors` (LiDAR + radar).
    /// See [`crate::sensors_workbench`].
    pub show_sensors_workbench: bool,
    /// Form + result state for the Sensors workbench (LiDAR scan / radar
    /// measurement over the in-process `valenx-sensors` engine).
    pub sensors: crate::sensors_workbench::SensorsWorkbenchState,

    /// Whether the right-side Fluids (SPH) workbench panel is visible. Defaults
    /// to `false`; toggled from the View menu or opened by the agent bridge
    /// under the id `"fluids"`. Wraps `valenx-fluids` (SPH particle simulation).
    /// See [`crate::fluids_workbench`].
    pub show_fluids_workbench: bool,
    /// Form + result state for the Fluids (SPH) workbench (particle simulation
    /// over the in-process `valenx-fluids` engine).
    pub fluids: crate::fluids_workbench::FluidsWorkbenchState,

    /// Whether the right-side Ocean workbench panel is visible. Defaults to
    /// `false`; toggled from the View menu or opened by the agent bridge under
    /// the id `"ocean"`. Wraps `valenx-ocean` (Gerstner wave field + quasi-static
    /// Archimedes buoyancy). See [`crate::ocean_workbench`].
    pub show_ocean_workbench: bool,
    /// Form + result state for the Ocean workbench (wave-height profile +
    /// floating-body settle over the in-process `valenx-ocean` engine).
    pub ocean: crate::ocean_workbench::OceanWorkbenchState,

    /// Whether the right-side ROM (reduced-order model) workbench panel is
    /// visible. Defaults to `false`; toggled from the View menu or opened by the
    /// agent bridge under the id `"rom"`. Wraps `valenx-rom` (POD / DMD / DEIM).
    /// See [`crate::rom_workbench`].
    pub show_rom_workbench: bool,
    /// Form + result state for the ROM workbench (POD energy spectrum +
    /// snapshot reconstruction over the in-process `valenx-rom` engine).
    pub rom: crate::rom_workbench::RomWorkbenchState,

    /// Whether the right-side UQ (uncertainty quantification) workbench panel is
    /// visible. Defaults to `false`; toggled from the View menu or opened by the
    /// agent bridge under the id `"uq"`. Wraps `valenx-uq` (Monte-Carlo
    /// propagation + Sobol sensitivity + FORM reliability). See
    /// [`crate::uq_workbench`].
    pub show_uq_workbench: bool,
    /// Form + result state for the UQ workbench (output histogram + Sobol bar
    /// chart + FORM Pf over the in-process `valenx-uq` engine).
    pub uq: crate::uq_workbench::UqWorkbenchState,

    /// Whether the right-side UAS (small-UAS design + defensive counter-UAS)
    /// workbench panel is visible. Defaults to `false`; toggled from the View
    /// menu or opened by the agent bridge under the id `"uas"` (aliases
    /// `"drone"` / `"counteruas"`). Wraps `valenx-uas` (multirotor / fixed-wing
    /// performance + trade study + defensive intercept GEOMETRY — no weapon
    /// employment). See [`crate::uas_workbench`].
    pub show_uas_workbench: bool,
    /// Form + result state for the UAS workbench (performance readout + trade
    /// Pareto scatter + counter-UAS intercept plan view over the in-process
    /// `valenx-uas` engine).
    pub uas: crate::uas_workbench::UasWorkbenchState,

    /// Whether the right-side Mission-Simulation (general discrete-event /
    /// agent constructive simulation) workbench panel is visible. Defaults to
    /// `false`; toggled from the View menu or opened by the agent bridge under
    /// the id `"missionsim"` (aliases `"mission"` / `"wargame"`). Wraps
    /// `valenx-mission-sim`: a discrete-event scheduler, analytic movers, and
    /// range-based detection, with ABSTRACT probabilistic engagement (a Pk input
    /// plus the Lanchester square-law ODE; no lethality, targeting, or
    /// kill-chain). See [`crate::missionsim_workbench`].
    pub show_missionsim_workbench: bool,
    /// Form + result state for the Mission-Simulation workbench (plan-view entity
    /// tracks + Lanchester force-vs-time plot + outcome metrics over the
    /// in-process `valenx-mission-sim` engine).
    pub missionsim: crate::missionsim_workbench::MissionSimWorkbenchState,

    /// Whether the right-side **Mission Planner** workbench panel is visible.
    /// Defaults to `false`; toggled from the View menu or opened by the agent
    /// bridge under the id `"missionplanner"`. A geographic (lat/lon) map with
    /// entities following waypoint routes, played back in real time
    /// (`valenx-mission-sim::planner`). Movement + routes only (Stage 1) — no
    /// engagement / sensors / orbits. See [`crate::mission_planner_workbench`].
    pub show_mission_planner_workbench: bool,
    /// Live scenario + playback state for the Mission Planner workbench.
    pub mission_planner: crate::mission_planner_workbench::MissionPlannerWorkbenchState,

    /// Whether the right-side Survivability / protection workbench panel is
    /// visible. Defaults to `false`; toggled from the View menu or opened by the
    /// agent bridge under the id `"survivability"` (aliases `"protection"` /
    /// `"blast"`). Wraps `valenx-survivability` — the DEFENSIVE / protective side
    /// of the shared blast/impact physics: free-field blast loading, SDOF
    /// protective response + the pressure-impulse iso-damage diagram, minimum
    /// armor sizing, and an occupant tolerance screen. Every output is "minimum
    /// protection to survive threat X"; no penetration / lethality is modeled.
    /// See [`crate::survivability_workbench`].
    pub show_survivability_workbench: bool,
    /// Form + result state for the Survivability / protection workbench
    /// (Friedlander pressure-time curve + P-I iso-damage diagram with the design
    /// point + armor / occupant readouts over the in-process
    /// `valenx-survivability` crate).
    pub survivability: crate::survivability_workbench::SurvivabilityWorkbenchState,

    /// Whether the right-side Photogrammetry / SfM scan workbench panel is
    /// visible. Defaults to `false`; toggled from the View menu or opened by the
    /// agent bridge under the id `"photogrammetry"` (aliases `"sfm"` / `"scan"`).
    /// Wraps `valenx-photogrammetry` (COLMAP-style structure-from-motion:
    /// features + matching + two-view geometry + incremental mapper + bundle
    /// adjustment). See [`crate::photogrammetry_workbench`].
    pub show_photogrammetry_workbench: bool,
    /// Form + result state for the Photogrammetry workbench (synthetic-scene SfM
    /// recovery: recovered sparse cloud + camera poses + reprojection error over
    /// the in-process `valenx-photogrammetry` mapper).
    pub photogrammetry: crate::photogrammetry_workbench::PhotogrammetryWorkbenchState,

    /// Whether the right-side Co-Simulation (FMI / HELICS) workbench panel is
    /// visible. Defaults to `false`; toggled from the View menu or opened by
    /// the agent bridge under the id `"cosim"` (aliases `"co-simulation"` /
    /// `"fmi"`). Wraps `valenx-adapter-fmi` (the in-house native co-sim master
    /// — Jacobi / Gauss-Seidel over a Subsystem graph — plus a strongly-coupled
    /// fixed-point implicit coupler). See [`crate::cosim_workbench`].
    pub show_cosim_workbench: bool,
    /// Form + result state for the Co-Simulation workbench (two coupled
    /// mass-spring-dampers exchanged through the in-process valenx-adapter-fmi
    /// coordinator: exchanged-signal history + coupling iterations + error vs a
    /// monolithic reference).
    pub cosim: crate::cosim_workbench::CosimWorkbenchState,

    /// Whether the right-side Protein-interaction (PPI / interactome) workbench
    /// panel is visible. Defaults to `false`; toggled from the View menu or
    /// opened by the agent bridge under the id `"ppi"` (aliases `"interactome"`
    /// / `"network"`). Wraps `valenx-ppi` (the in-house sequence-first
    /// coevolution PPI engine — APC-corrected mutual-information over a paired
    /// MSA folded into a fused [0,1] score, plus an all-vs-all interactome
    /// screen). See [`crate::ppi_workbench`].
    pub show_ppi_workbench: bool,
    /// Form + result state for the PPI / interactome workbench (a named demo
    /// host × pathogen interactome screened through the in-process valenx-ppi
    /// engine: the real scored interaction network + degree / betweenness
    /// centrality + BFS shortest path computed over it).
    pub ppi: crate::ppi_workbench::PpiWorkbenchState,

    /// Whether the right-side Multibody-dynamics (robot / contact) workbench
    /// panel is visible. Defaults to `false`; toggled from the View menu or
    /// opened by the agent bridge under the id `"mbd"` (aliases `"multibody"` /
    /// `"robot"`). Wraps `valenx-mbd` (the in-house planar constrained-DAE
    /// multibody solver, its Featherstone articulated-body algorithm, and its
    /// penalty-contact + Coulomb-friction model). See [`crate::mbd_workbench`].
    pub show_mbd_workbench: bool,
    /// Form + result state for the Multibody-dynamics workbench (an
    /// energy-conserving articulated rod pendulum advanced by the real
    /// constrained-DAE `System`, and a body dropped onto a plane through the
    /// real penalty-contact + Coulomb-friction path — trajectory + contact-force
    /// history + energy / penetration diagnostics).
    pub mbd: crate::mbd_workbench::MbdWorkbenchState,

    /// Whether the right-side Autonomy V&V workbench panel is visible. Defaults
    /// to `false`; toggled from the View menu or opened by the agent bridge
    /// under the id `"autonomy"`. Wraps `valenx-autonomy-vnv` (scenario-based
    /// verification of an autonomous vehicle with simulated sensors). See
    /// [`crate::autonomy_workbench`].
    pub show_autonomy_workbench: bool,
    /// Form + result state for the Autonomy V&V workbench (scenario → trace →
    /// requirement report over the in-process `valenx-autonomy-vnv` framework).
    pub autonomy: crate::autonomy_workbench::AutonomyWorkbenchState,

    /// Whether the right-side Assistant activity sidebar is visible. On by
    /// default (set in [`ValenxApp::new`]) so the app narrates its own work
    /// via the live feed.
    pub show_assistant_panel: bool,
    /// State for the Assistant activity sidebar (the live `.jsonl` feed
    /// path). See [`crate::assistant_workbench`].
    pub assistant: crate::assistant_workbench::AssistantWorkbenchState,

    /// Whether the keyboard-shortcut cheat-sheet overlay is open.
    /// Toggled by the `?` key + by Help → Keyboard shortcuts.
    pub keyboard_help_open: bool,

    /// Whether the per-panel contextual help popup is open. Mapped
    /// to F1 + Help → Panel help.
    pub panel_help_open: bool,

    /// Which panel's help text the F1 popup shows. Resolved at the
    /// moment F1 is pressed from "what workbench is active right
    /// now"; defaults to "Sequence" when nothing else is up.
    pub panel_help_target: String,

    /// Whether the first-launch welcome tour is currently open. Auto-set
    /// to `true` on a fresh install (gated by `settings.welcome_tour_completed`);
    /// re-openable from the Help menu.
    pub welcome_tour_open: bool,

    /// Tour navigation state — which step the user is on, and
    /// whether they've finished. See [`crate::welcome_tour::TourState`].
    pub welcome_tour_state: crate::welcome_tour::TourState,

    /// "File → New Project…" modal state. `Some(_)` while the dialog
    /// is open; `None` once the user clicks Create / Cancel / closes
    /// the window. Triggered by the File menu, the command palette,
    /// and the Ctrl+N shortcut. See [`crate::new_project_dialog`].
    pub new_project_dialog: Option<crate::new_project_dialog::NewProjectDialog>,

    /// One-line notice rendered inline on the welcome / landing page
    /// (next to the recent-projects list). Set by the host when a
    /// landing-page action produces a result that doesn't belong in
    /// the top status bar — currently the "removed missing project
    /// from recents" confirmation. Cleared as soon as the user takes
    /// another action on the landing page.
    pub landing_inline_message: Option<String>,

    /// Memoised command-palette entry list, keyed by
    /// `(registry.len(), library.content_rev(), show_non_oss_adapters,
    /// focus_category)`.
    /// `build_visible_commands` allocates ~360 `String`s per call and used
    /// to run every frame; the cache invalidates when the registry grows
    /// (rare — re-probe / load), when the saved-project list changes (the
    /// launcher lists one entry per saved project, so add/remove/**rename**
    /// must rebuild), when the OSS-only toggle flips in Settings, or when the
    /// **domain-focus filter** changes (it narrows the launcher's workbench-tab
    /// entries). `None` until the first palette render fills it.
    ///
    /// The second key is [`crate::project_library::ProjectLibrary::content_rev`],
    /// a content fingerprint over each project's `(id, name)` — so an in-place
    /// rename (which leaves `projects.len()` unchanged) now flips the key and
    /// the launcher shows the new name on the next frame. The fourth key is
    /// [`ValenxApp::focus_category`] so flipping the focus rebuilds immediately.
    /// See the cache-build site in `update.rs`. Type aliased as [`PaletteCache`].
    pub palette_cache: PaletteCache,

    // ── Swappable viewport system (cloud/viewport) ────────────────────────
    /// Which viewport implementation is rendered in the central panel.
    ///
    /// Defaults to `Viewport3D`; switches to `Viewport2dDna` when the
    /// user first enables the Genetics Workbench (and can be overridden
    /// at any time from **View → Central viewport**).
    pub active_viewport: crate::viewport_kind::ViewportKind,

    /// Open project tabs (Chrome-style) plus the active index. Drives
    /// which workbench the tab strip shows. See [`crate::project_tabs`].
    pub tab_bar: crate::project_tabs::TabBar,

    /// Pending tab-close confirmation. `Some(i)` while the "Close tab?"
    /// modal is open for tab index `i` (set by the strip's ✕ / right-click
    /// "Close"); cleared on Cancel, or on confirm right after the tab is
    /// actually closed. Closing a tab discards its (unsaved) workspace
    /// document, so the close is gated behind this explicit confirm. See
    /// [`crate::project_tabs::draw_tab_strip`].
    pub tab_close_confirm: Option<usize>,

    /// Pending **"Close all tabs?"** confirmation. `Some(())` while the
    /// "Close all N tabs?" modal is open (set by the strip toolbar's
    /// "Close all tabs" button); cleared on Cancel, or on confirm right after
    /// every tab is closed via [`crate::project_tabs::TabBar::close_all`].
    /// Closing all tabs discards each tab's unsaved workspace document, so the
    /// batch close is gated behind this explicit confirm. See
    /// [`crate::project_tabs::draw_tab_strip`].
    pub tab_close_all_confirm: Option<()>,

    /// Pending "Save as project…" prompt opened from a tab's right-click
    /// menu. `Some` while the modal is up; carries the source tab index, the
    /// in-progress project name, and the chosen destination folder id (None =
    /// unfiled). On confirm it calls `library.add_project` + persists. Drawn
    /// by `project_tabs`'s `draw_save_as_project`.
    pub tab_save_as_project: Option<crate::project_tabs::SaveAsProjectPrompt>,

    /// The foldered, persistent **project library** — saved projects the
    /// user manages from the Browser's "Projects" navigator (search /
    /// folders / pin / reorder / reopen-as-tab). Loaded from
    /// `<state_dir>/library.json` in [`ValenxApp::new`]; persisted on every
    /// mutation. See [`crate::project_library`] / [`crate::project_navigator`].
    pub library: crate::project_library::ProjectLibrary,

    /// Transient UI state for the project navigator (search box text, the
    /// inline-rename editor target + buffer, the "New folder" name prompt).
    /// Never persisted — rebuilt empty each launch. See
    /// [`crate::project_navigator::NavigatorState`].
    pub nav_state: crate::project_navigator::NavigatorState,

    /// Transient UI state for the Browser panel's **"Open Tabs"** list (the
    /// VS-Code-style "Open Editors" pane mirroring every open tab) — currently
    /// just the search-box text. Never persisted. See
    /// [`crate::project_tabs::OpenTabsState`] / `draw_open_tabs_list`.
    pub open_tabs_state: crate::project_tabs::OpenTabsState,

    /// Persistent state for the 2D DNA / plasmid viewport. Survives
    /// viewport-kind switches so pan, zoom, and sub-view selection are
    /// remembered when the user returns to the 2D view.
    pub viewport_2d: crate::viewport_2d::Viewport2dState,

    /// **Domain-focus filter** — the working domain the workbench surfaces are
    /// narrowed to, or [`None`] for "All" (show everything, the default).
    ///
    /// A pure-UI focus layer: when `Some(category)`, the Ctrl+P launcher's
    /// `OpenWorkbenchTab` entries, the tab strip's "From template" menu, and the
    /// View menu's primary-workbench toggles show only workbenches whose
    /// [`crate::project_tabs::TabKind::group`] category equals it; "All"
    /// ([`None`]) shows everything. Nothing is removed or feature-gated — see
    /// [`crate::workbench_focus`]. **In-session only** (transient view
    /// preference, not written to the settings file), so it resets to "All" on
    /// relaunch.
    pub focus_category: Option<String>,
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
        // Restore the foldered project library from `<state_dir>/library.json`
        // so the Browser's "Projects" navigator shows the user's saved
        // projects on launch. A missing / corrupt file yields an empty
        // library (see `ProjectLibrary::load`).
        app.library = crate::project_library::ProjectLibrary::load();
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
        // Wipe any leftover agent-command files from a previous run so the
        // agent-drives-valenx bridge starts each session with a clean slate
        // (the cursor logic then runs every newly-appended command from line 0
        // without replaying stale history). `app.assistant` is initialized in
        // `Self::default()` above, so its inbox path resolves here.
        crate::agent_commands::clear_command_files(&app);

        // Agent-bridge wake thread. The in-`update()` heartbeat
        // (`ctx.request_repaint_after(POLL_INTERVAL)` in `update.rs`) only
        // fires *while `update()` runs*, and egui is reactive: when valenx is
        // idle (occluded, in the background, or otherwise not receiving input —
        // the normal case while an external agent drives it) frames stop, so
        // the agent-command poll stops with them and appended commands sit
        // unread. A detached background thread holding a clone of the egui
        // `Context` (`egui::Context` is `Send + Sync`) pokes the event loop on
        // a fixed cadence regardless of window/focus state: cross-thread
        // `request_repaint()` routes through eframe's
        // `set_request_repaint_callback`, which posts a winit
        // `UserEvent::RequestRepaint` via the `EventLoopProxy` — delivered even
        // when no OS input is arriving. Spawned once per process (guarded
        // below) because `eframe::run_native` may invoke this constructor more
        // than once under `run_and_return`.
        //
        // CAVEAT (honest): this does NOT guarantee polling while the window is
        // *minimized* on eframe 0.28. The run-and-return event loop hard-gates
        // painting on `window.is_minimized()` (eframe `native/run.rs`): when
        // iconified it drops the scheduled repaint without calling
        // `request_redraw`, so `update()` (and the poll) do not run until the
        // window is restored. There is no `NativeOptions`/`ViewportBuilder`
        // option to disable that gate. What the thread DOES guarantee is an
        // immediate flush the instant the window is restored, and continuous
        // ~3 Hz polling whenever the window is merely idle/background/occluded
        // but not minimized. Needs runtime verification.
        {
            use std::sync::atomic::{AtomicBool, Ordering};
            static WAKE_THREAD_STARTED: AtomicBool = AtomicBool::new(false);
            if WAKE_THREAD_STARTED
                .compare_exchange(false, true, Ordering::SeqCst, Ordering::SeqCst)
                .is_ok()
            {
                let ctx = cc.egui_ctx.clone();
                std::thread::Builder::new()
                    .name("valenx-agent-wake".to_owned())
                    .spawn(move || loop {
                        std::thread::sleep(std::time::Duration::from_millis(300));
                        ctx.request_repaint();
                    })
                    // A failure to spawn the wake thread is non-fatal: the app
                    // still runs and the in-`update()` heartbeat covers the
                    // interactive case, so we only log rather than abort launch.
                    .unwrap_or_else(|e| {
                        tracing::warn!(
                            target: "valenx",
                            "failed to spawn agent-bridge wake thread: {e}"
                        );
                        // Return a dummy handle so the closure's type is
                        // `JoinHandle<()>`; the thread simply never started.
                        std::thread::spawn(|| {})
                    });
            }
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

    /// Build a `LoadedMesh` whose node cloud spans the given AABB, so a test
    /// can assert `ensure_default_animation`'s pivot is the box centre. The
    /// mesh carries only nodes (no elements) — the quality / histogram
    /// companions run fine on an element-less mesh (empty results), which is
    /// all `ensure_default_animation` (node-only AABB) needs.
    fn loaded_mesh_spanning(min: [f64; 3], max: [f64; 3]) -> LoadedMesh {
        let mut mesh = valenx_mesh::Mesh::new("test-aabb");
        // Two opposite corners are enough to fix the AABB; add the origin to
        // prove the centre is the midpoint of the extent, not of the points.
        mesh.nodes
            .push(nalgebra::Vector3::new(min[0], min[1], min[2]));
        mesh.nodes
            .push(nalgebra::Vector3::new(max[0], max[1], max[2]));
        mesh.nodes.push(nalgebra::Vector3::new(0.0, 0.0, 0.0));
        let quality = valenx_mesh::quality_report(&mesh);
        let aspect_hist =
            valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
        let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
        LoadedMesh {
            path: PathBuf::from("<test>/aabb"),
            mesh,
            quality,
            aspect_hist,
            skew_hist,
        }
    }

    /// A mesh `WorkspaceProduct` (carries `mesh: Some`) with no animation yet —
    /// the shape every bridge-rendered registry product has before
    /// `ensure_default_animation` runs.
    fn mesh_product_without_animation(min: [f64; 3], max: [f64; 3]) -> WorkspaceProduct {
        WorkspaceProduct {
            title: "Mesh part".into(),
            lines: Vec::new(),
            mesh: Some(loaded_mesh_spanning(min, max)),
            vertex_colors: None,
            camera: OrbitCamera::default(),
            kind2d: None,
            last_export: None,
            image: None,
            image_texture: None,
            animation: None,
        }
    }

    #[test]
    fn ensure_default_animation_adds_turntable_at_aabb_centre() {
        // A mesh product with no animation gets a paused +Z turntable whose
        // pivot is the AABB centre — here the box [-2,-4,-6]..[4,6,10] (which
        // already contains the helper's origin node, so it is the full extent)
        // centres on [1, 1, 2].
        let mut product = mesh_product_without_animation([-2.0, -4.0, -6.0], [4.0, 6.0, 10.0]);
        product.ensure_default_animation();
        let anim = product
            .animation
            .as_ref()
            .expect("a default animation was attached to the mesh product");
        assert!(!anim.playing, "the default inspect-spin starts paused");
        assert_eq!(anim.speed, 1.0, "default speed is 1.0×");
        assert_eq!(anim.t, 0.0, "the clock starts at zero");
        match &anim.motion {
            ProductMotion::Turntable {
                axis,
                pivot,
                rad_per_s,
            } => {
                assert_eq!(*axis, [0.0, 0.0, 1.0], "spins about +Z");
                assert_eq!(*rad_per_s, 0.4, "~1 rev / 15 s gentle inspect rate");
                for (got, want) in pivot.iter().zip([1.0_f32, 1.0, 2.0].iter()) {
                    assert!(
                        (got - want).abs() < 1e-5,
                        "pivot {pivot:?} is the AABB centre"
                    );
                }
            }
            other => panic!("expected a Turntable, got {other:?}"),
        }
    }

    #[test]
    fn ensure_default_animation_is_a_no_op_on_an_already_animated_product() {
        // A product that already animates (the gear's RigidParts) keeps its own
        // motion untouched — the default is NOT overwritten onto it.
        let parts = vec![RigidPart {
            node_range: 0..3,
            axis: [0.0, 1.0, 0.0],
            pivot: [5.0, 0.0, 5.0],
            rad_per_s: -2.0,
        }];
        let mut product = mesh_product_without_animation([-1.0, -1.0, -1.0], [1.0, 1.0, 1.0]);
        product.animation = Some(ProductAnimation {
            playing: true,
            speed: 3.0,
            t: 1.5,
            motion: ProductMotion::RigidParts(parts.clone()),
        });
        product.ensure_default_animation();
        let anim = product.animation.as_ref().expect("animation retained");
        assert!(anim.playing, "the existing playing flag is preserved");
        assert_eq!(anim.speed, 3.0, "existing speed untouched");
        assert_eq!(anim.t, 1.5, "existing clock untouched");
        match &anim.motion {
            ProductMotion::RigidParts(p) => assert_eq!(*p, parts, "RigidParts left intact"),
            other => panic!("RigidParts must not be replaced by a Turntable, got {other:?}"),
        }
    }

    #[test]
    fn ensure_default_animation_is_a_no_op_on_a_mesh_less_product() {
        // A 2-D / card / image product (mesh: None) never animates — no control
        // is conjured onto it.
        let mut product = WorkspaceProduct {
            title: "Card".into(),
            lines: vec!["a result".into()],
            mesh: None,
            vertex_colors: None,
            camera: OrbitCamera::default(),
            kind2d: None,
            last_export: None,
            image: None,
            image_texture: None,
            animation: None,
        };
        product.ensure_default_animation();
        assert!(
            product.animation.is_none(),
            "no animation is attached to a mesh-less product"
        );
    }
}
