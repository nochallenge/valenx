//! Chrome-style **project tabs** — an open-many, switch-between strip of
//! project workspaces.
//!
//! valenx's domain tools are independent right-dock workbench panels,
//! toggled from the View menu. The tab strip is a thin navigation layer
//! over them: each tab owns one [`TabKind`], **activating** a tab shows
//! that kind's workbench and hides the others (so the user works one
//! project at a time, like browser tabs), the `➕ New tab` button opens a
//! fresh **blank** project, and the painter-drawn `✕` closes one.
//!
//! ## New-tab behaviour
//!
//! `➕ New tab` creates a **blank, named** project immediately — no forced
//! workbench pick and no folder dialog. The user types a name (e.g.
//! "boat") and starts building. The categorised workbench launcher is kept
//! as a *secondary* affordance: the `＋ from template ▾` menu opens a tab
//! pre-bound to a specific workbench kind. That menu's body is wrapped in
//! [`crate::menu_ui::scrollable_menu`] so the long category list stays
//! on-screen and scrolls.
//!
//! ## Naming, renaming, and saving
//!
//! A tab can be renamed inline (double-click the title, or right-click →
//! *Rename*); the name persists on the [`ProjectTab`]. Tabs can be saved to
//! disk individually (`<state_dir>/tabs/<name>.json`) or as a whole named
//! *session* group (`<state_dir>/sessions/<name>.json`), and reopened later
//! from the `Open saved ▾` menu. Persistence is plain serde-JSON through the
//! crash-safe [`crate::state_paths::atomic_write`].
//!
//! The strip is **additive and non-breaking**: a fresh app starts with zero
//! tabs and the existing default layout untouched. Tab mode only engages
//! once the user opens the first tab.
//!
//! ## Per-tab workspace documents
//!
//! Each tab owns its own **scene / project document** — the loaded
//! geometry, mesh, camera, and selected case/field/time. These live in
//! [`WorkspaceDoc`]; the [`TabBar`] keeps one `WorkspaceDoc` per tab in
//! [`TabBar::docs`] (invariant: `docs.len() == tabs.len()`). The *active*
//! tab's document is checked **out** into the live [`ValenxApp`] fields
//! (`app.project`, `app.mesh`, `app.camera`, …) so the rest of the app
//! reads/writes one plain set of fields; `docs[active]` is a default
//! placeholder while that tab is checked out. Switching tabs swaps the
//! live fields back into the old tab's slot and installs the new tab's
//! document ([`switch_active_to`]), so opening a blank "+ New tab" gives a
//! genuinely empty scene and switching back restores the prior geometry.
//!
//! Per-*workbench* parameter state (e.g. the rocket-design inputs, the CFD
//! config) is **not** yet per-tab — only the scene/project document above
//! is isolated. App-global runtime (the adapter registry, residual/log
//! panels, run/sweep handles, settings, the `show_*` workbench flags) stays
//! shared across tabs by design. The `docs` vector is runtime-only and is
//! not serialised: a [`SavedSession`] still persists just `{name, tabs,
//! active}`, and a restored session rebuilds `docs` as fresh defaults.

use eframe::egui;
use serde::{Deserialize, Serialize};

use crate::state_paths::{atomic_write, state_dir};
use crate::types::{LoadedMesh, LoadedStl};
use crate::viewport_kind::ViewportKind;
use crate::ValenxApp;
use std::path::PathBuf;
use valenx_core::LoadedProject;
use valenx_viz::OrbitCamera;

/// One tab's **workspace document** — the per-tab scene / project state.
///
/// These are exactly the [`ValenxApp`] fields that make up "what's loaded
/// in this project": the loaded project + its path + RBAC override, the
/// dropped STL, the canonical mesh, the orbit camera, and the
/// currently-selected case / field / time-index, plus the last run's
/// results + workdir. The *active* tab keeps its document checked out in
/// the live `ValenxApp` fields; every other tab parks its document in
/// [`TabBar::docs`].
///
/// It **also** carries the per-tab *dockable layout / viewport view state*:
/// whether the dockable workbench layout is on (`dock_enabled`),
/// whether the central 3-D viewport is hidden / collapsed
/// (`viewport_hidden` / `viewport_collapsed`), and the
/// dock's own tile tree (`dock_tree`). This is what makes the
/// "Workbench + Agent" grid **per-tab**: a tab that has six agent units
/// keeps them, while a freshly-opened tab gets a clean view (dock off,
/// viewport shown, no tree) and so shows its workbench + the 3-D viewport
/// rather than another tab's agent grid. Note the per-unit chat-channel
/// counter [`ValenxApp::wb_agent_counter`] is intentionally **not** here —
/// it stays global so `agent:n` ids never collide across tabs.
///
/// Construction is via [`Default`] (a fresh, empty document — no project,
/// no mesh, the default camera, dock off, viewport shown, no dock tree).
/// `WorkspaceDoc::capture` moves the live fields *out* of an app (leaving
/// them empty / default), and `WorkspaceDoc::install` moves a document's
/// fields back *in*. Neither clones the meshes or the dock tree — they are
/// `move`d through `Option`/`Box`/`take`.
#[derive(Default)]
pub struct WorkspaceDoc {
    project: Option<LoadedProject>,
    project_path: Option<PathBuf>,
    project_rbac_override: Option<valenx_rbac::RbacConfig>,
    stl: Option<LoadedStl>,
    mesh: Option<LoadedMesh>,
    camera: OrbitCamera,
    selected_case: Option<String>,
    selected_field_name: Option<String>,
    selected_time_index: usize,
    last_run_results: Option<Box<valenx_fields::Results>>,
    last_run_workdir: Option<PathBuf>,
    /// Per-tab: is the dockable workbench layout (incl. any "Workbench +
    /// Agent" grid) on for this tab? Mirrors [`ValenxApp::dock_enabled`].
    dock_enabled: bool,
    /// Per-tab: is the central 3-D viewport hidden for this tab? Mirrors
    /// [`ValenxApp::viewport_hidden`].
    viewport_hidden: bool,
    /// Per-tab: is the central viewport rolled up to its header for this
    /// tab? Mirrors [`ValenxApp::viewport_collapsed`].
    viewport_collapsed: bool,
    /// Per-tab: this tab's dock tile tree (the layout of its docked
    /// workbenches / agent units). Mirrors [`ValenxApp::dock_tree`].
    /// `take`n in/out — never cloned.
    dock_tree: Option<egui_tiles::Tree<String>>,
}

impl WorkspaceDoc {
    /// Move this tab's scene/project + dock/view fields **out** of `app`,
    /// leaving the live fields empty / default (so the caller can immediately
    /// [`install`](Self::install) another document into the now-cleared
    /// app). No mesh and no dock tree is cloned: every `Option`/`Box` is
    /// `take`n and the camera is swapped for its default. The dock booleans
    /// reset to `false` (dock off / viewport shown) in the live app, and the
    /// dock tree is `take`n out wholesale.
    fn capture(app: &mut ValenxApp) -> WorkspaceDoc {
        WorkspaceDoc {
            project: app.project.take(),
            project_path: app.project_path.take(),
            project_rbac_override: app.project_rbac_override.take(),
            stl: app.stl.take(),
            mesh: app.mesh.take(),
            camera: std::mem::take(&mut app.camera),
            selected_case: app.selected_case.take(),
            selected_field_name: app.selected_field_name.take(),
            selected_time_index: std::mem::take(&mut app.selected_time_index),
            last_run_results: app.last_run_results.take(),
            last_run_workdir: app.last_run_workdir.take(),
            // Per-tab dock / viewport view state. The dock_tree MUST be
            // `take`n (not cloned) — egui_tiles trees are not Clone here and
            // the tab owns its layout outright.
            dock_enabled: std::mem::take(&mut app.dock_enabled),
            viewport_hidden: std::mem::take(&mut app.viewport_hidden),
            viewport_collapsed: std::mem::take(&mut app.viewport_collapsed),
            dock_tree: app.dock_tree.take(),
        }
    }

    /// Move this document's fields **into** `app`, replacing whatever the
    /// live scene/project + dock/view fields hold. Pair with
    /// [`capture`](Self::capture) (capture the outgoing tab first, then
    /// install the incoming one) so a scene is never lost. Installing a
    /// fresh [`Default`] document is exactly how a newly-opened tab gets a
    /// clean view: dock off, viewport shown, no dock tree. Consumes `self`.
    fn install(self, app: &mut ValenxApp) {
        app.project = self.project;
        app.project_path = self.project_path;
        app.project_rbac_override = self.project_rbac_override;
        app.stl = self.stl;
        app.mesh = self.mesh;
        app.camera = self.camera;
        app.selected_case = self.selected_case;
        app.selected_field_name = self.selected_field_name;
        app.selected_time_index = self.selected_time_index;
        app.last_run_results = self.last_run_results;
        app.last_run_workdir = self.last_run_workdir;
        // Per-tab dock / viewport view state (swapped in so the active tab's
        // dock grid shows and a clean tab's does not).
        app.dock_enabled = self.dock_enabled;
        app.viewport_hidden = self.viewport_hidden;
        app.viewport_collapsed = self.viewport_collapsed;
        app.dock_tree = self.dock_tree;
    }
}

/// A project kind a tab can hold. [`TabKind::Blank`] is an empty project
/// (the default `➕ New tab`); every other variant maps to exactly one
/// primary workbench panel (the `show_*` flag it drives on [`ValenxApp`]).
#[derive(Clone, Copy, PartialEq, Eq, Debug, Serialize, Deserialize)]
pub enum TabKind {
    /// A blank, empty project — no workbench is forced open. The user names
    /// it and chooses what to build from the menus.
    Blank,
    // -- Aerospace --
    Rocket,
    Engine,
    Astro,
    Aero,
    Gasdynamics,
    // -- Simulation --
    Cfd,
    Fem,
    Reactdyn,
    Fields,
    // -- CAD & mesh --
    Cad,
    MeshToolbox,
    Sheetmetal,
    Reverse,
    Draft2d,
    Render,
    Animate,
    // -- Machine design --
    Springs,
    Gears,
    Fasteners,
    Frames,
    Collision,
    // -- Civil & AEC --
    Piping,
    Hvac,
    Reinforcement,
    Interior,
    Geomatics,
    // -- Life sciences --
    Genetics,
    Neuro,
    VariantEffect,
}

impl TabKind {
    /// Every *template* kind (i.e. excluding [`TabKind::Blank`]), in
    /// `＋ from template`-menu order (grouped via [`Self::group`]).
    pub const TEMPLATES: [TabKind; 29] = [
        TabKind::Rocket,
        TabKind::Engine,
        TabKind::Astro,
        TabKind::Aero,
        TabKind::Gasdynamics,
        TabKind::Cfd,
        TabKind::Fem,
        TabKind::Reactdyn,
        TabKind::Fields,
        TabKind::Cad,
        TabKind::MeshToolbox,
        TabKind::Sheetmetal,
        TabKind::Reverse,
        TabKind::Draft2d,
        TabKind::Render,
        TabKind::Animate,
        TabKind::Springs,
        TabKind::Gears,
        TabKind::Fasteners,
        TabKind::Frames,
        TabKind::Collision,
        TabKind::Piping,
        TabKind::Hvac,
        TabKind::Reinforcement,
        TabKind::Interior,
        TabKind::Geomatics,
        TabKind::Genetics,
        TabKind::Neuro,
        TabKind::VariantEffect,
    ];

    /// Group header shown in the `＋ from template` new-tab menu. Blank is
    /// not menu-grouped (it has its own dedicated button).
    pub fn group(self) -> &'static str {
        match self {
            TabKind::Blank => "General",
            TabKind::Rocket
            | TabKind::Engine
            | TabKind::Astro
            | TabKind::Aero
            | TabKind::Gasdynamics => "Aerospace",
            TabKind::Cfd | TabKind::Fem | TabKind::Reactdyn | TabKind::Fields => "Simulation",
            TabKind::Cad
            | TabKind::MeshToolbox
            | TabKind::Sheetmetal
            | TabKind::Reverse
            | TabKind::Draft2d
            | TabKind::Render
            | TabKind::Animate => "CAD & mesh",
            TabKind::Springs
            | TabKind::Gears
            | TabKind::Fasteners
            | TabKind::Frames
            | TabKind::Collision => "Machine design",
            TabKind::Piping
            | TabKind::Hvac
            | TabKind::Reinforcement
            | TabKind::Interior
            | TabKind::Geomatics => "Civil & AEC",
            TabKind::Genetics | TabKind::Neuro | TabKind::VariantEffect => "Life sciences",
        }
    }

    /// Tab + menu label.
    pub fn label(self) -> &'static str {
        match self {
            TabKind::Blank => "Untitled",
            TabKind::Rocket => "Rocket",
            TabKind::Engine => "Engine",
            TabKind::Astro => "Astro / Launch",
            TabKind::Aero => "Aerodynamics",
            TabKind::Gasdynamics => "Gas dynamics",
            TabKind::Cfd => "CFD",
            TabKind::Fem => "FEM",
            TabKind::Reactdyn => "Reaction dynamics",
            TabKind::Fields => "Field statistics",
            TabKind::Cad => "Parametric CAD",
            TabKind::MeshToolbox => "Mesh toolbox",
            TabKind::Sheetmetal => "Sheet metal",
            TabKind::Reverse => "Reverse engineering",
            TabKind::Draft2d => "2D drafting",
            TabKind::Render => "Path-traced render",
            TabKind::Animate => "Animation",
            TabKind::Springs => "Springs",
            TabKind::Gears => "Gears",
            TabKind::Fasteners => "Fasteners",
            TabKind::Frames => "Frames / sections",
            TabKind::Collision => "Collision",
            TabKind::Piping => "Piping",
            TabKind::Hvac => "HVAC",
            TabKind::Reinforcement => "Reinforcement",
            TabKind::Interior => "Interior design",
            TabKind::Geomatics => "Geomatics",
            TabKind::Genetics => "Genetics",
            TabKind::Neuro => "Neural interface",
            TabKind::VariantEffect => "Variant effect",
        }
    }

    /// Turn this kind's workbench panel **on**. Callers clear every panel
    /// first via [`clear_all_workbenches`] so exactly one is left visible.
    /// [`TabKind::Blank`] opens no workbench (an empty project).
    fn show(self, app: &mut ValenxApp) {
        match self {
            TabKind::Blank => {} // empty project — nothing forced open
            TabKind::Rocket => app.show_rocket_workbench = true,
            TabKind::Engine => app.show_engine_workbench = true,
            TabKind::Astro => app.show_astro_workbench = true,
            TabKind::Aero => app.show_aero_workbench = true,
            TabKind::Gasdynamics => app.show_gasdynamics_workbench = true,
            TabKind::Cfd => app.show_cfd_workbench = true,
            TabKind::Fem => app.show_fem_workbench = true,
            TabKind::Reactdyn => app.show_reactdyn_workbench = true,
            TabKind::Fields => app.show_fields_workbench = true,
            TabKind::Cad => app.show_cad_workbench = true,
            TabKind::MeshToolbox => app.show_mesh_toolbox = true,
            TabKind::Sheetmetal => app.show_sheetmetal_workbench = true,
            TabKind::Reverse => app.show_reverse_workbench = true,
            TabKind::Draft2d => app.show_draft2d_workbench = true,
            TabKind::Render => app.show_render_workbench = true,
            TabKind::Animate => app.show_animate_workbench = true,
            TabKind::Springs => app.show_springs_workbench = true,
            TabKind::Gears => app.show_gears_workbench = true,
            TabKind::Fasteners => app.show_fasteners_workbench = true,
            TabKind::Frames => app.show_frames_workbench = true,
            TabKind::Collision => app.show_collision_workbench = true,
            TabKind::Piping => app.show_piping_workbench = true,
            TabKind::Hvac => app.show_hvac_workbench = true,
            TabKind::Reinforcement => app.show_reinforcement_workbench = true,
            TabKind::Interior => app.show_interior_workbench = true,
            TabKind::Geomatics => app.show_geomatics_workbench = true,
            TabKind::Genetics => app.show_genetics_workbench = true,
            TabKind::Neuro => app.show_neuro_workbench = true,
            TabKind::VariantEffect => app.show_variant_effect_workbench = true,
        }
    }

    /// Map a short, **case-insensitive** workbench id string to a [`TabKind`]
    /// (e.g. `"rocket"` → [`TabKind::Rocket`], `"varianteffect"` →
    /// [`TabKind::VariantEffect`]). Returns `None` for an unknown id — callers
    /// (the agent-drives-valenx bridge in [`crate::agent_commands`]) then fall
    /// back to a blank tab / skip rather than panicking.
    ///
    /// Accepts a couple of friendly aliases where the workbench has more than
    /// one common name: `mesh`/`meshtoolbox`, `variant`/`varianteffect`. This is
    /// the inverse of the ids an external agent is told to emit; keep it in sync
    /// with [`Self::TEMPLATES`].
    pub fn from_id(s: &str) -> Option<TabKind> {
        match s.trim().to_ascii_lowercase().as_str() {
            "rocket" => Some(TabKind::Rocket),
            "engine" => Some(TabKind::Engine),
            "astro" => Some(TabKind::Astro),
            "aero" => Some(TabKind::Aero),
            "gasdynamics" => Some(TabKind::Gasdynamics),
            "cfd" => Some(TabKind::Cfd),
            "fem" => Some(TabKind::Fem),
            "reactdyn" => Some(TabKind::Reactdyn),
            "fields" => Some(TabKind::Fields),
            "cad" => Some(TabKind::Cad),
            "mesh" | "meshtoolbox" => Some(TabKind::MeshToolbox),
            "sheetmetal" => Some(TabKind::Sheetmetal),
            "reverse" => Some(TabKind::Reverse),
            "draft2d" => Some(TabKind::Draft2d),
            "render" => Some(TabKind::Render),
            "animate" => Some(TabKind::Animate),
            "springs" => Some(TabKind::Springs),
            "gears" => Some(TabKind::Gears),
            "fasteners" => Some(TabKind::Fasteners),
            "frames" => Some(TabKind::Frames),
            "collision" => Some(TabKind::Collision),
            "piping" => Some(TabKind::Piping),
            "hvac" => Some(TabKind::Hvac),
            "reinforcement" => Some(TabKind::Reinforcement),
            "interior" => Some(TabKind::Interior),
            "geomatics" => Some(TabKind::Geomatics),
            "genetics" => Some(TabKind::Genetics),
            "neuro" => Some(TabKind::Neuro),
            "variant" | "varianteffect" => Some(TabKind::VariantEffect),
            _ => None,
        }
    }

    /// The central viewport this kind prefers (genetics is the 2D DNA
    /// view; everything else — including a blank project — is the 3D
    /// viewport).
    fn viewport(self) -> ViewportKind {
        match self {
            TabKind::Genetics => ViewportKind::Viewport2dDna,
            _ => ViewportKind::Viewport3D,
        }
    }
}

/// One open project tab: its kind plus a user-facing title. The two
/// `edit_*` fields drive inline rename and are transient (never persisted).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ProjectTab {
    /// The project kind this tab hosts.
    pub kind: TabKind,
    /// Title shown on the tab (defaults to the kind label, or "Untitled N"
    /// for a blank tab).
    pub title: String,
    /// `true` while the title is being edited inline. Transient.
    #[serde(skip)]
    pub editing: bool,
    /// Scratch buffer backing the inline rename [`egui::TextEdit`].
    /// Transient.
    #[serde(skip)]
    pub edit_buf: String,
}

impl ProjectTab {
    /// A tab of `kind` with the given title and no active edit.
    fn new(kind: TabKind, title: impl Into<String>) -> Self {
        ProjectTab {
            kind,
            title: title.into(),
            editing: false,
            edit_buf: String::new(),
        }
    }
}

/// A saved set of tabs plus the active index — the on-disk form of a tab
/// *group* (a named session). A single saved tab uses the same envelope
/// with one entry.
#[derive(Clone, Serialize, Deserialize)]
pub struct SavedSession {
    /// Display name of the session/group.
    pub name: String,
    /// Tabs, left to right.
    pub tabs: Vec<ProjectTab>,
    /// Active tab index within `tabs`, if any.
    pub active: Option<usize>,
}

/// The project-tab strip state, owned by [`ValenxApp`].
#[derive(Default)]
pub struct TabBar {
    /// Open tabs, left to right.
    pub tabs: Vec<ProjectTab>,
    /// Per-tab workspace documents, **index-aligned with [`Self::tabs`]**
    /// (invariant: `docs.len() == tabs.len()`). `docs[i]` holds tab `i`'s
    /// parked scene/project document — except for the *active* tab, whose
    /// document is checked out into the live [`ValenxApp`] fields while
    /// `docs[active]` is a default placeholder. Runtime-only: never
    /// serialised (a [`SavedSession`] stores just `{name, tabs, active}`),
    /// so a restored/appended session rebuilds these as fresh defaults.
    pub docs: Vec<WorkspaceDoc>,
    /// Index of the active tab in [`Self::tabs`], or `None` when empty.
    pub active: Option<usize>,
    /// Monotonic counter feeding the default "Untitled N" name for blank
    /// tabs, so successive blanks get distinct titles.
    pub blank_counter: usize,
}

impl TabBar {
    /// Open a blank, empty project tab with an auto-generated "Untitled N"
    /// name, make it active, and return its index. This is what the default
    /// `New tab` button does — no workbench is forced open. Pushes a fresh
    /// [`WorkspaceDoc`] alongside the tab to keep `docs.len() ==
    /// tabs.len()`. Does **not** swap the live document — the caller (see
    /// `apply_intent`) runs [`switch_active_to`] right after so the
    /// previous tab's scene is parked and this blank tab starts empty.
    pub fn open_blank(&mut self) -> usize {
        self.blank_counter += 1;
        let title = format!("Untitled {}", self.blank_counter);
        self.tabs.push(ProjectTab::new(TabKind::Blank, title));
        self.docs.push(WorkspaceDoc::default());
        let idx = self.tabs.len() - 1;
        self.active = Some(idx);
        idx
    }

    /// Open a new tab of `kind` (titled with the kind label), make it
    /// active, and return its index. Used by the `From template` menu.
    /// Pushes a fresh [`WorkspaceDoc`] to preserve the `docs.len() ==
    /// tabs.len()` invariant (the live-document swap is done by the caller
    /// via [`switch_active_to`]).
    pub fn open(&mut self, kind: TabKind) -> usize {
        self.tabs.push(ProjectTab::new(kind, kind.label()));
        self.docs.push(WorkspaceDoc::default());
        let idx = self.tabs.len() - 1;
        self.active = Some(idx);
        idx
    }

    /// Close the tab at `idx` (and its parked document) and pick a sensible
    /// new active tab (the previous neighbour, or `None` when the strip
    /// empties). Keeps `docs` index-aligned with `tabs`. The live document
    /// is reconciled by the caller (see the `close` branch of
    /// `apply_intent`), which installs the new active tab's document.
    pub fn close(&mut self, idx: usize) {
        if idx >= self.tabs.len() {
            return;
        }
        self.tabs.remove(idx);
        if idx < self.docs.len() {
            self.docs.remove(idx);
        }
        self.active = if self.tabs.is_empty() {
            None
        } else {
            Some(idx.min(self.tabs.len() - 1))
        };
    }

    /// The active tab's kind, if any.
    pub fn active_kind(&self) -> Option<TabKind> {
        self.active.and_then(|i| self.tabs.get(i)).map(|t| t.kind)
    }

    /// Snapshot the whole strip as a [`SavedSession`] named `name` — the
    /// in-memory form of a saved group.
    pub fn snapshot(&self, name: impl Into<String>) -> SavedSession {
        SavedSession {
            name: name.into(),
            tabs: self.tabs.clone(),
            active: self.active,
        }
    }

    /// Replace the whole strip with the tabs from `session`, clearing the
    /// transient edit state and clamping `active` into range. Used when the
    /// user reopens a saved group. Rebuilds [`Self::docs`] as one fresh
    /// default per restored tab (documents are not serialised), preserving
    /// the `docs.len() == tabs.len()` invariant.
    pub fn restore(&mut self, session: SavedSession) {
        self.tabs = session
            .tabs
            .into_iter()
            .map(|mut t| {
                t.editing = false;
                t.edit_buf.clear();
                t
            })
            .collect();
        self.docs = (0..self.tabs.len())
            .map(|_| WorkspaceDoc::default())
            .collect();
        self.active = match session.active {
            Some(i) if i < self.tabs.len() => Some(i),
            _ if self.tabs.is_empty() => None,
            _ => Some(0),
        };
    }

    /// Append the tabs from `session` after the current ones (used to
    /// reopen a *single* saved tab without discarding the open set), make
    /// the first appended tab active, and return its index if any. Pushes a
    /// fresh default [`WorkspaceDoc`] for each appended tab so `docs` stays
    /// index-aligned with `tabs`.
    pub fn append(&mut self, session: SavedSession) -> Option<usize> {
        if session.tabs.is_empty() {
            return None;
        }
        let first = self.tabs.len();
        for mut t in session.tabs {
            t.editing = false;
            t.edit_buf.clear();
            self.tabs.push(t);
            self.docs.push(WorkspaceDoc::default());
        }
        self.active = Some(first);
        Some(first)
    }
}

// ---------------------------------------------------------------------------
// Persistence — save / load single tabs and whole groups under the state dir.
// ---------------------------------------------------------------------------

/// Directory holding single-tab saves: `<state_dir>/tabs`.
fn tabs_dir() -> Option<std::path::PathBuf> {
    state_dir().map(|d| d.join("tabs"))
}

/// Directory holding tab-group / session saves: `<state_dir>/sessions`.
fn sessions_dir() -> Option<std::path::PathBuf> {
    state_dir().map(|d| d.join("sessions"))
}

/// Sanitise a user-supplied save name into a safe single-path-segment file
/// stem: keep alphanumerics, space, dash, underscore and dot; map every
/// other char to `_`; collapse to "untitled" if nothing usable remains.
/// This keeps the name from escaping the target directory (no `/`, `\`,
/// `..`, drive letters) regardless of what the user typed.
fn sanitize_name(name: &str) -> String {
    let cleaned: String = name
        .chars()
        .map(|c| {
            if c.is_alphanumeric() || matches!(c, ' ' | '-' | '_' | '.') {
                c
            } else {
                '_'
            }
        })
        .collect();
    let trimmed = cleaned.trim().trim_matches('.').trim();
    if trimmed.is_empty() {
        "untitled".to_string()
    } else {
        trimmed.to_string()
    }
}

/// Serialize `session` to pretty JSON. Separated out so the round-trip tests
/// can exercise (de)serialisation without touching the filesystem.
fn to_json(session: &SavedSession) -> Result<String, serde_json::Error> {
    serde_json::to_string_pretty(session)
}

/// Parse a [`SavedSession`] from JSON, clearing transient edit flags.
fn from_json(text: &str) -> Result<SavedSession, serde_json::Error> {
    let mut s: SavedSession = serde_json::from_str(text)?;
    for t in &mut s.tabs {
        t.editing = false;
        t.edit_buf.clear();
    }
    Ok(s)
}

// The save/load/list functions split into a thin public wrapper (which
// resolves the real per-OS state dir) and an `*_in(dir, …)` inner that takes
// the directory explicitly. The inner form keeps the I/O testable without
// mutating the process-global state-dir env var — the round-trip tests pass a
// throwaway temp dir directly, so they stay parallel-safe.

/// Persist a session into `dir`, named `<stem>.json`. Best-effort: a write
/// failure is swallowed (returns `false`).
fn save_session_in(dir: &std::path::Path, stem: &str, session: &SavedSession) -> bool {
    let Ok(text) = to_json(session) else {
        return false;
    };
    atomic_write(&dir.join(format!("{stem}.json")), &text).is_ok()
}

/// Persist a single tab to `<state_dir>/tabs/<name>.json`. Best-effort: a
/// missing state dir or write failure is swallowed (returns `false`),
/// mirroring the rest of valenx's state persistence.
pub fn save_single_tab(tab: &ProjectTab) -> bool {
    let Some(dir) = tabs_dir() else {
        return false;
    };
    let session = SavedSession {
        name: tab.title.clone(),
        tabs: vec![tab.clone()],
        active: Some(0),
    };
    save_session_in(&dir, &sanitize_name(&tab.title), &session)
}

/// Persist the whole tab strip as a named group to
/// `<state_dir>/sessions/<name>.json`. Best-effort (returns `false` on a
/// missing state dir or write failure).
pub fn save_group(bar: &TabBar, name: &str) -> bool {
    let Some(dir) = sessions_dir() else {
        return false;
    };
    save_session_in(&dir, &sanitize_name(name), &bar.snapshot(name))
}

/// Load a [`SavedSession`] from a JSON file, bounded to a sane size so a
/// corrupt/hostile file can't OOM the load. Returns `None` on any error.
fn load_session_file(path: &std::path::Path) -> Option<SavedSession> {
    let meta = std::fs::metadata(path).ok()?;
    // Sessions are tiny; cap reads well above any realistic size.
    if meta.len() > crate::settings_io::MAX_STATE_FILE_BYTES as u64 {
        return None;
    }
    let text = std::fs::read_to_string(path).ok()?;
    from_json(&text).ok()
}

/// Load `<dir>/<sanitized name>.json` as a session.
fn load_session_in(dir: &std::path::Path, name: &str) -> Option<SavedSession> {
    load_session_file(&dir.join(format!("{}.json", sanitize_name(name))))
}

/// List saved single-tab names (the file stems under `<state_dir>/tabs`),
/// sorted. Empty when the dir is absent.
pub fn list_saved_tabs() -> Vec<String> {
    list_json_stems(tabs_dir())
}

/// List saved group/session names (file stems under
/// `<state_dir>/sessions`), sorted. Empty when the dir is absent.
pub fn list_saved_groups() -> Vec<String> {
    list_json_stems(sessions_dir())
}

/// Collect the `*.json` file stems in `dir`, sorted alphabetically.
fn list_json_stems(dir: Option<std::path::PathBuf>) -> Vec<String> {
    let Some(dir) = dir else {
        return Vec::new();
    };
    let Ok(rd) = std::fs::read_dir(&dir) else {
        return Vec::new();
    };
    let mut out: Vec<String> = rd
        .flatten()
        .filter_map(|e| {
            let p = e.path();
            if p.extension().is_some_and(|x| x == "json") {
                p.file_stem().and_then(|s| s.to_str()).map(str::to_string)
            } else {
                None
            }
        })
        .collect();
    out.sort();
    out
}

/// Load a saved single tab by name (its sanitised stem) from
/// `<state_dir>/tabs/<name>.json`.
pub fn load_saved_tab(name: &str) -> Option<SavedSession> {
    load_session_in(&tabs_dir()?, name)
}

/// Load a saved group by name from `<state_dir>/sessions/<name>.json`.
pub fn load_saved_group(name: &str) -> Option<SavedSession> {
    load_session_in(&sessions_dir()?, name)
}

// ---------------------------------------------------------------------------
// Workbench reconciliation.
// ---------------------------------------------------------------------------

/// Hide every project workbench panel. The active tab (if any) then
/// re-shows exactly one via [`TabKind::show`] (or none, for a blank tab).
fn clear_all_workbenches(app: &mut ValenxApp) {
    app.show_rocket_workbench = false;
    app.show_engine_workbench = false;
    app.show_astro_workbench = false;
    app.show_aero_workbench = false;
    app.show_gasdynamics_workbench = false;
    app.show_cfd_workbench = false;
    app.show_fem_workbench = false;
    app.show_reactdyn_workbench = false;
    app.show_fields_workbench = false;
    app.show_cad_workbench = false;
    app.show_mesh_toolbox = false;
    app.show_sheetmetal_workbench = false;
    app.show_reverse_workbench = false;
    app.show_draft2d_workbench = false;
    app.show_render_workbench = false;
    app.show_animate_workbench = false;
    app.show_springs_workbench = false;
    app.show_gears_workbench = false;
    app.show_fasteners_workbench = false;
    app.show_frames_workbench = false;
    app.show_collision_workbench = false;
    app.show_piping_workbench = false;
    app.show_hvac_workbench = false;
    app.show_reinforcement_workbench = false;
    app.show_interior_workbench = false;
    app.show_geomatics_workbench = false;
    app.show_genetics_workbench = false;
    app.show_neuro_workbench = false;
    app.show_variant_effect_workbench = false;
}

/// Reconcile the visible workbench + central viewport with the active
/// tab: clear every project panel, then show the active tab's kind and
/// switch the viewport to match. A blank tab shows no workbench (just the
/// 3D viewport). With no active tab, everything stays hidden (the user
/// closed the last tab).
pub fn sync_active(app: &mut ValenxApp) {
    let kind = app.tab_bar.active_kind();
    clear_all_workbenches(app);
    if let Some(kind) = kind {
        kind.show(app);
        app.active_viewport = kind.viewport();
    }
}

/// Switch the active tab to `new_idx`, swapping the per-tab workspace
/// document so each tab keeps its own scene/project.
///
/// The currently-active tab's live scene (`app.project`, `app.mesh`,
/// `app.camera`, …) is captured (`WorkspaceDoc::capture`) back into its
/// `docs` slot, then `docs[new_idx]` is taken and
/// installed (`WorkspaceDoc::install`) into the live fields — so the
/// outgoing tab keeps its geometry and the incoming tab shows its own
/// (empty for a fresh blank tab). Finally `active` is set and the visible
/// workbench + viewport are reconciled via [`sync_active`].
///
/// `new_idx` out of range is ignored (no-op). When the previous `active`
/// was `None` (pre-tab mode — the user just opened the first tab), the old
/// live scene is intentionally discarded: the new tab installs its own
/// (default/empty) document so "+ New tab" starts fresh.
pub fn switch_active_to(app: &mut ValenxApp, new_idx: usize) {
    if new_idx >= app.tab_bar.docs.len() {
        return;
    }
    // Park the outgoing tab's live scene back into its slot (if any).
    if let Some(a) = app.tab_bar.active {
        if a < app.tab_bar.docs.len() {
            app.tab_bar.docs[a] = WorkspaceDoc::capture(app);
        } else {
            // Defensive: drop the live scene if the old index is stale.
            let _ = WorkspaceDoc::capture(app);
        }
    } else {
        // Pre-tab mode: discard the previous live scene so the first tab
        // starts from a clean document.
        let _ = WorkspaceDoc::capture(app);
    }
    // Install the incoming tab's document into the live fields.
    let doc = std::mem::take(&mut app.tab_bar.docs[new_idx]);
    doc.install(app);
    app.tab_bar.active = Some(new_idx);
    sync_active(app);
}

/// Park the **currently-active** tab's live scene back into its `docs` slot
/// (a no-op when there is no active tab or the index is stale). This is the
/// "capture the outgoing tab before we re-point `active`" step shared by every
/// open/restore/append path — extracted so the agent-drives-valenx bridge
/// ([`crate::agent_commands`]) can reuse the exact same reconcile the UI uses
/// rather than duplicating the snippet.
pub(crate) fn park_active_doc(app: &mut ValenxApp) {
    if let Some(a) = app.tab_bar.active {
        if a < app.tab_bar.docs.len() {
            app.tab_bar.docs[a] = WorkspaceDoc::capture(app);
        }
    }
}

/// Reconcile the live workspace document with whatever `app.tab_bar.active`
/// already points at, **discarding** the current live scene.
///
/// Unlike [`switch_active_to`] (which parks the outgoing scene), this drops
/// the live fields and installs `docs[active]` (or clears them to a default
/// empty document when there is no active tab). It's the right reconcile
/// after operations that rebuild / replace the tab set and have *already*
/// set `active` themselves — restoring a saved group, appending a saved
/// tab, or closing a tab — where the outgoing live scene either no longer
/// has a home or is being deliberately replaced. Always ends with
/// [`sync_active`]. Exposed `pub(crate)` so [`crate::agent_commands`]'s
/// `NewTab` reducer can finish an open exactly as the tab strip does.
pub(crate) fn install_active_doc(app: &mut ValenxApp) {
    // Drop the current live scene (its tab is gone / being replaced).
    let _ = WorkspaceDoc::capture(app);
    match app.tab_bar.active {
        Some(i) if i < app.tab_bar.docs.len() => {
            let doc = std::mem::take(&mut app.tab_bar.docs[i]);
            doc.install(app);
        }
        // No active tab (or stale index): leave the live fields empty.
        _ => WorkspaceDoc::default().install(app),
    }
    sync_active(app);
}

/// Actually close the tab at `idx` (and discard its workspace document),
/// reconciling the live scene afterwards. Called once the user confirms the
/// "Close tab?" modal. Preserves the active tab's live scene first so
/// closing a *non-active* tab doesn't lose the active tab's geometry;
/// closing the active tab discards its scene (its slot is removed) and the
/// neighbour's document is installed.
fn perform_close(app: &mut ValenxApp, idx: usize) {
    if idx >= app.tab_bar.tabs.len() {
        return;
    }
    park_active_doc(app);
    app.tab_bar.close(idx);
    install_active_doc(app);
}

/// Working state of the "Save as project…" modal (opened from a tab's
/// right-click menu). Held on [`crate::ValenxApp::tab_save_as_project`] while
/// the dialog is up. On confirm, the tab at `tab_idx` is cloned into the
/// foldered project library under `name` in `folder` (None = unfiled).
pub struct SaveAsProjectPrompt {
    /// Index of the source tab being saved.
    pub tab_idx: usize,
    /// In-progress project name (seeded from the tab title).
    pub name: String,
    /// Chosen destination folder id, or `None` for "(unfiled)".
    pub folder: Option<String>,
}

/// What a single frame of the tab strip wants to do, accumulated while the
/// read-only borrow of the tab vec is live and applied afterwards.
#[derive(Default)]
struct StripIntent {
    activate: Option<usize>,
    /// **Request** to close the tab at this index — opens the "Close tab?"
    /// confirmation modal rather than closing immediately. The real close
    /// only happens once the user confirms (see [`perform_close`]).
    request_close: Option<usize>,
    open_template: Option<TabKind>,
    open_blank: bool,
    save_tab: Option<usize>,
    /// Open the "Save as project…" prompt for the tab at this index — adds
    /// the tab to the foldered project library (see [`crate::project_library`]).
    save_as_project: Option<usize>,
    save_group: bool,
    open_saved_group: Option<String>,
    open_saved_tab: Option<String>,
    /// Commit an inline rename: (tab index, new title).
    commit_rename: Option<(usize, String)>,
    /// Begin an inline rename of the tab at this index.
    begin_rename: Option<usize>,
    /// Open a paired "Workbench + Agent" unit (an empty workspace tile + a
    /// Claude chat tile) in the dockable region via
    /// [`ValenxApp::add_workbench_agent_pair`]. Used by any caller that wants
    /// the simple "new bottom row" placement (e.g. the View menu).
    open_wb_agent: bool,
    /// Open a paired "Workbench + Agent" unit at a **chosen** grid position,
    /// picked from the tab-strip "+ Workbench+Agent" placement dropdown.
    /// Routed to [`ValenxApp::add_workbench_agent_pair_at`].
    add_wb_agent_at: Option<crate::dock_layout::UnitAddTarget>,
}

/// Draw the project-tab strip (a slim panel just below the ribbon) and
/// apply any click this frame (open blank / open template / activate /
/// request-close / rename / save / open-saved), then render the
/// "Close tab?" confirmation modal if a close is pending.
pub fn draw_tab_strip(app: &mut ValenxApp, ctx: &egui::Context) {
    let mut intent = StripIntent::default();

    egui::TopBottomPanel::top("valenx_project_tabs").show(ctx, |ui| {
        ui.horizontal(|ui| {
            // Primary: instant blank named project (no forced workbench, no
            // folder dialog). Plain ASCII label so no font-glyph "tofu" box.
            if ui
                .button("+ New tab")
                .on_hover_text("New blank project — name it and start building")
                .clicked()
            {
                intent.open_blank = true;
            }

            // Paired "Workbench + Agent" unit — an empty workspace tile + a
            // Claude chat tile dropped into the dockable region (turns the
            // dockable layout on). A dropdown lets the user PLACE the new unit
            // precisely: a brand-new row at top/bottom, or into an existing row
            // (left/right end) of the current grid. The row list is read live
            // from the dock tree (`dock_grid_rows`) — safe here because the
            // dock_tree is owned (`Some`) during the tab strip; the dock itself
            // renders later in `update.rs`. Plain ASCII labels (no glyph carets)
            // so nothing renders as a "tofu" box. Body wrapped in
            // `scrollable_menu` so a tall grid's row list stays on-screen.
            use crate::dock_layout::UnitAddTarget;
            ui.menu_button("+ Workbench+Agent", |ui| {
                crate::menu_ui::scrollable_menu(ui, |ui| {
                    if ui
                        .button("New row at top")
                        .on_hover_text("Add the unit as a new first row")
                        .clicked()
                    {
                        intent.add_wb_agent_at = Some(UnitAddTarget::NewRowTop);
                        ui.close_menu();
                    }
                    if ui
                        .button("New row at bottom")
                        .on_hover_text("Add the unit as a new last row")
                        .clicked()
                    {
                        intent.add_wb_agent_at = Some(UnitAddTarget::NewRowBottom);
                        ui.close_menu();
                    }
                    // Live grid shape: one entry per row, with its unit count.
                    let rows = app.dock_grid_rows();
                    if !rows.is_empty() {
                        ui.separator();
                        ui.label(egui::RichText::new("Add into a row:").weak().small());
                        for (i, units) in rows.iter().enumerate() {
                            ui.horizontal(|ui| {
                                ui.label(format!("Row {} ({} units)", i + 1, units));
                                if ui
                                    .small_button("left")
                                    .on_hover_text("Add at the left end of this row")
                                    .clicked()
                                {
                                    intent.add_wb_agent_at = Some(UnitAddTarget::RowStart(i));
                                    ui.close_menu();
                                }
                                if ui
                                    .small_button("right")
                                    .on_hover_text("Add at the right end of this row")
                                    .clicked()
                                {
                                    intent.add_wb_agent_at = Some(UnitAddTarget::RowEnd(i));
                                    ui.close_menu();
                                }
                            });
                        }
                    }
                });
            });

            // Secondary: start a tab pre-bound to a workbench template. The
            // body is wrapped in `scrollable_menu` so the long category list
            // stays on-screen and scrolls instead of running off the bottom.
            // ASCII label (no glyph caret) so it never renders as tofu.
            ui.menu_button("From template", |ui| {
                crate::menu_ui::scrollable_menu(ui, |ui| {
                    let mut last_group = "";
                    for kind in TabKind::TEMPLATES {
                        let group = kind.group();
                        if group != last_group {
                            if !last_group.is_empty() {
                                ui.separator();
                            }
                            ui.label(egui::RichText::new(group).small().weak());
                            last_group = group;
                        }
                        if ui.button(kind.label()).clicked() {
                            intent.open_template = Some(kind);
                            ui.close_menu();
                        }
                    }
                });
            });

            // Open a previously-saved tab or group. ASCII label (no glyph
            // caret) so it never renders as tofu.
            ui.menu_button("Open saved", |ui| {
                crate::menu_ui::scrollable_menu(ui, |ui| {
                    let groups = list_saved_groups();
                    let tabs = list_saved_tabs();
                    if groups.is_empty() && tabs.is_empty() {
                        ui.label(egui::RichText::new("(nothing saved yet)").weak().small());
                    }
                    if !groups.is_empty() {
                        ui.label(egui::RichText::new("Groups (sessions)").small().weak());
                        for name in groups {
                            // Plain text item (the emoji prefix tofu'd).
                            if ui.button(format!("Group: {name}")).clicked() {
                                intent.open_saved_group = Some(name);
                                ui.close_menu();
                            }
                        }
                        if !tabs.is_empty() {
                            ui.separator();
                        }
                    }
                    if !tabs.is_empty() {
                        ui.label(egui::RichText::new("Single tabs").small().weak());
                        for name in tabs {
                            // Plain text item (the emoji prefix tofu'd).
                            if ui.button(format!("Tab: {name}")).clicked() {
                                intent.open_saved_tab = Some(name);
                                ui.close_menu();
                            }
                        }
                    }
                });
            });

            // Save the whole open set as a named group/session.
            if !app.tab_bar.tabs.is_empty()
                && ui
                    .button("Save group…")
                    .on_hover_text("Save all open tabs as a named session")
                    .clicked()
            {
                intent.save_group = true;
            }

            ui.separator();

            if app.tab_bar.tabs.is_empty() {
                // Plain ASCII (the arrow glyph tofu'd).
                ui.label(egui::RichText::new("New tab to begin").weak().small());
            }

            let active = app.tab_bar.active;
            // Iterate by index so the inline-edit buffer can be mutated.
            for i in 0..app.tab_bar.tabs.len() {
                let selected = active == Some(i);
                let editing = app.tab_bar.tabs[i].editing;

                if editing {
                    // Inline rename: a single-line text field, committed on
                    // Enter or focus loss.
                    let resp = ui.add(
                        egui::TextEdit::singleline(&mut app.tab_bar.tabs[i].edit_buf)
                            .desired_width(120.0)
                            .id_source(("tab_rename", i)),
                    );
                    resp.request_focus();
                    let lost_focus = resp.lost_focus();
                    let enter = lost_focus && ui.input(|inp| inp.key_pressed(egui::Key::Enter));
                    if enter || lost_focus {
                        intent.commit_rename = Some((i, app.tab_bar.tabs[i].edit_buf.clone()));
                    }
                } else {
                    let label = app.tab_bar.tabs[i].title.clone();
                    let group = app.tab_bar.tabs[i].kind.group();
                    let resp = ui
                        .selectable_label(selected, label)
                        .on_hover_text(format!("{group} — double-click to rename"));
                    if resp.clicked() {
                        intent.activate = Some(i);
                    }
                    if resp.double_clicked() {
                        intent.begin_rename = Some(i);
                    }
                    // Right-click context menu: rename / save / close. "Save
                    // this tab" is the escape hatch before a discard-on-close.
                    resp.context_menu(|ui| {
                        if ui.button("Rename").clicked() {
                            intent.begin_rename = Some(i);
                            ui.close_menu();
                        }
                        if ui.button("Save this tab").clicked() {
                            intent.save_tab = Some(i);
                            ui.close_menu();
                        }
                        // Add to the foldered project library (the Browser
                        // "Projects" navigator) via a name + folder prompt.
                        if ui
                            .button("Save as project…")
                            .on_hover_text("Add this tab to the Projects library (Browser panel)")
                            .clicked()
                        {
                            intent.save_as_project = Some(i);
                            ui.close_menu();
                        }
                        ui.separator();
                        if ui.button("Close").clicked() {
                            // Request a close — opens the confirm modal.
                            intent.request_close = Some(i);
                            ui.close_menu();
                        }
                    });
                }

                // Painter-drawn ✕ (reused from the workbench chrome) — never
                // a font-glyph "tofu" box. Requests a close (the confirm
                // modal gates the actual discard).
                if crate::workbench_chrome::close_x_button(ui, "Close tab").clicked() {
                    intent.request_close = Some(i);
                }
                ui.separator();
            }
        });
    });

    apply_intent(app, intent);
    draw_close_confirm(app, ctx);
    draw_save_as_project(app, ctx);
}

/// Render the "Close tab?" confirmation modal while
/// [`ValenxApp::tab_close_confirm`] is `Some`. Closing a tab discards its
/// (unsaved) workspace document, so the destructive close is gated behind
/// an explicit confirm. [Cancel] clears the pending index; [Close tab]
/// performs the real close (+ document removal + live-scene reconcile) via
/// [`perform_close`]. The dialog points the user at "Save this tab" (the
/// right-click escape hatch) so work can be preserved first.
fn draw_close_confirm(app: &mut ValenxApp, ctx: &egui::Context) {
    let Some(idx) = app.tab_close_confirm else {
        return;
    };
    // The index may have gone stale (tab removed another way) — bail safely.
    let Some(title) = app.tab_bar.tabs.get(idx).map(|t| t.title.clone()) else {
        app.tab_close_confirm = None;
        return;
    };

    let mut do_close = false;
    let mut do_cancel = false;
    egui::Window::new("Close tab?")
        .collapsible(false)
        .resizable(false)
        .anchor(egui::Align2::CENTER_CENTER, egui::vec2(0.0, 0.0))
        .show(ctx, |ui| {
            ui.label(format!(
                "Close \"{title}\"? This tab and its unsaved work will be permanently discarded."
            ));
            ui.add_space(4.0);
            ui.label(
                egui::RichText::new("Tip: right-click a tab to Save it first.")
                    .small()
                    .weak(),
            );
            ui.add_space(8.0);
            ui.horizontal(|ui| {
                if ui.button("Cancel").clicked() {
                    do_cancel = true;
                }
                // Red-ish destructive action.
                let close_btn = egui::Button::new(
                    egui::RichText::new("Close tab").color(egui::Color32::from_rgb(220, 80, 80)),
                );
                if ui.add(close_btn).clicked() {
                    do_close = true;
                }
            });
        });

    if do_cancel {
        app.tab_close_confirm = None;
    } else if do_close {
        perform_close(app, idx);
        app.tab_close_confirm = None;
    }
}

/// Render the "Save as project…" modal while
/// [`crate::ValenxApp::tab_save_as_project`] is `Some`. Mirrors
/// [`draw_close_confirm`]: an anchored, non-collapsible window with a name
/// field, a folder picker (existing library folders + "(unfiled)"), and
/// Save / Cancel. On Save it clones the source tab into the project library
/// under the entered name + folder (the tab is `Clone` + `Serialize`) and
/// persists `library.json`. The project name overrides the cloned tab's
/// title so the library entry reads as the user typed.
fn draw_save_as_project(app: &mut ValenxApp, ctx: &egui::Context) {
    let Some(prompt) = &app.tab_save_as_project else {
        return;
    };
    let idx = prompt.tab_idx;
    // Bail safely if the source tab vanished (closed another way).
    let Some(tab_title) = app.tab_bar.tabs.get(idx).map(|t| t.title.clone()) else {
        app.tab_save_as_project = None;
        return;
    };

    let mut do_save = false;
    let mut do_cancel = false;
    // Snapshot the folder list for the picker (immutable read of the lib).
    let folders: Vec<(String, String)> = app
        .library
        .sorted_folders()
        .into_iter()
        .map(|f| (f.id.clone(), f.name.clone()))
        .collect();

    egui::Window::new("Save as project")
        .collapsible(false)
        .resizable(false)
        .anchor(egui::Align2::CENTER_CENTER, egui::vec2(0.0, 0.0))
        .show(ctx, |ui| {
            ui.label(format!("Save \"{tab_title}\" to the Projects library."));
            ui.add_space(6.0);
            ui.horizontal(|ui| {
                ui.label("Name:");
                if let Some(p) = &mut app.tab_save_as_project {
                    ui.add(
                        egui::TextEdit::singleline(&mut p.name)
                            .desired_width(220.0)
                            .hint_text("Project name"),
                    );
                }
            });
            ui.add_space(4.0);
            ui.horizontal(|ui| {
                ui.label("Folder:");
                // Current selection label.
                let current = app
                    .tab_save_as_project
                    .as_ref()
                    .and_then(|p| p.folder.clone());
                let current_label = current
                    .as_ref()
                    .and_then(|fid| {
                        folders
                            .iter()
                            .find(|(id, _)| id == fid)
                            .map(|(_, n)| n.clone())
                    })
                    .unwrap_or_else(|| "(unfiled)".to_string());
                egui::ComboBox::from_id_source("save_as_project_folder")
                    .selected_text(current_label)
                    .show_ui(ui, |ui| {
                        if let Some(p) = &mut app.tab_save_as_project {
                            ui.selectable_value(&mut p.folder, None, "(unfiled)");
                            for (fid, fname) in &folders {
                                ui.selectable_value(&mut p.folder, Some(fid.clone()), fname);
                            }
                        }
                    });
            });
            if folders.is_empty() {
                ui.label(
                    egui::RichText::new("Tip: create folders in the Browser → Projects navigator.")
                        .small()
                        .weak(),
                );
            }
            ui.add_space(8.0);
            ui.horizontal(|ui| {
                if ui.button("Cancel").clicked() {
                    do_cancel = true;
                }
                if ui.button("Save").clicked() {
                    do_save = true;
                }
            });
        });

    if do_cancel {
        app.tab_save_as_project = None;
    } else if do_save {
        // Pull the entered name + folder, clone the source tab, override its
        // title with the project name, add to the library, persist.
        if let Some(prompt) = app.tab_save_as_project.take() {
            if let Some(tab) = app.tab_bar.tabs.get(prompt.tab_idx) {
                let mut saved_tab = tab.clone();
                // Transient inline-edit state should never persist.
                saved_tab.editing = false;
                saved_tab.edit_buf.clear();
                let trimmed = prompt.name.trim();
                if !trimmed.is_empty() {
                    saved_tab.title = trimmed.to_string();
                }
                app.library.add_project(saved_tab, prompt.folder);
                let _ = app.library.save();
            }
        }
    }
}

/// Apply this frame's accumulated [`StripIntent`] after the read-only
/// borrows in [`draw_tab_strip`] end. At most a couple of these fire per
/// frame in practice; each leaves `active` consistent.
fn apply_intent(app: &mut ValenxApp, intent: StripIntent) {
    if let Some((i, new_title)) = intent.commit_rename {
        if let Some(tab) = app.tab_bar.tabs.get_mut(i) {
            let trimmed = new_title.trim();
            if !trimmed.is_empty() {
                tab.title = trimmed.to_string();
            }
            tab.editing = false;
            tab.edit_buf.clear();
        }
    }
    if let Some(i) = intent.begin_rename {
        if let Some(tab) = app.tab_bar.tabs.get_mut(i) {
            tab.edit_buf = tab.title.clone();
            tab.editing = true;
        }
    }
    if let Some(i) = intent.save_tab {
        if let Some(tab) = app.tab_bar.tabs.get(i) {
            let _ = save_single_tab(tab);
        }
    }
    if let Some(i) = intent.save_as_project {
        // Open the "Save as project…" prompt seeded with the tab's title.
        if let Some(tab) = app.tab_bar.tabs.get(i) {
            app.tab_save_as_project = Some(SaveAsProjectPrompt {
                tab_idx: i,
                name: tab.title.clone(),
                folder: None,
            });
        }
    }
    if intent.save_group {
        // Name the group after the active tab (or "session"); a future
        // dialog could prompt, but auto-naming keeps the flow one click.
        let name = app
            .tab_bar
            .active
            .and_then(|i| app.tab_bar.tabs.get(i))
            .map(|t| t.title.clone())
            .unwrap_or_else(|| "session".to_string());
        let _ = save_group(&app.tab_bar, &name);
    }
    if let Some(name) = intent.open_saved_group {
        if let Some(session) = load_saved_group(&name) {
            // `restore` rebuilds the strip + fresh default docs and sets
            // `active`; the old live scene is discarded with it.
            app.tab_bar.restore(session);
            install_active_doc(app);
        }
    }
    if let Some(name) = intent.open_saved_tab {
        if let Some(session) = load_saved_tab(&name) {
            // Park the currently-active scene before `append` re-points
            // `active` at the first appended (fresh) tab, so switching back
            // restores it.
            park_active_doc(app);
            app.tab_bar.append(session);
            install_active_doc(app);
        }
    }
    if let Some(i) = intent.request_close {
        // The ✕ / right-click "Close" only *requests* a close; the real
        // close happens once the user confirms the "Close tab?" modal (a
        // tab's unsaved workspace document is discarded on close).
        if i < app.tab_bar.tabs.len() {
            app.tab_close_confirm = Some(i);
        }
    }
    if let Some(i) = intent.activate {
        // Swap documents so each tab keeps its own scene.
        if i < app.tab_bar.tabs.len() && app.tab_bar.active != Some(i) {
            switch_active_to(app, i);
        }
    }
    if let Some(kind) = intent.open_template {
        // Park the outgoing tab's scene, open the new tab (pushes a fresh
        // default doc + makes it active), then install that empty doc so the
        // new tab starts blank and the prior tab keeps its geometry.
        park_active_doc(app);
        app.tab_bar.open(kind);
        install_active_doc(app);
    }
    if intent.open_blank {
        park_active_doc(app);
        app.tab_bar.open_blank();
        install_active_doc(app);
    }
    if intent.open_wb_agent {
        // Drop a paired Workspace + Agent unit into the dockable region
        // (this also turns the dockable layout on). Independent of the
        // project-tab document state, so no tab/doc reconcile is needed.
        app.add_workbench_agent_pair();
    }
    if let Some(target) = intent.add_wb_agent_at {
        // Same as above, but the dropdown chose a precise grid position for
        // the new unit (new top/bottom row, or into an existing row's
        // left/right end). Also turns the dockable layout on.
        app.add_workbench_agent_pair_at(target);
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod tests {
    use super::*;

    #[test]
    fn templates_are_unique_and_grouped() {
        // No duplicate kinds in TEMPLATES.
        for (i, a) in TabKind::TEMPLATES.iter().enumerate() {
            for b in &TabKind::TEMPLATES[i + 1..] {
                assert_ne!(a, b, "duplicate kind in TEMPLATES: {a:?}");
            }
        }
        // Blank is intentionally NOT a template.
        assert!(
            !TabKind::TEMPLATES.contains(&TabKind::Blank),
            "Blank must not appear in the template menu"
        );
        // Every kind (incl. Blank) has a non-empty label and group.
        for k in TabKind::TEMPLATES {
            assert!(!k.label().is_empty());
            assert!(!k.group().is_empty());
        }
        assert!(!TabKind::Blank.label().is_empty());
        assert!(!TabKind::Blank.group().is_empty());
    }

    #[test]
    fn from_id_maps_every_template_kind_case_insensitively() {
        // The canonical id for each kind maps back to it (the inverse the
        // agent-drives-valenx bridge relies on), and matching is
        // case-insensitive / whitespace-tolerant.
        let canonical = [
            ("rocket", TabKind::Rocket),
            ("engine", TabKind::Engine),
            ("astro", TabKind::Astro),
            ("aero", TabKind::Aero),
            ("gasdynamics", TabKind::Gasdynamics),
            ("cfd", TabKind::Cfd),
            ("fem", TabKind::Fem),
            ("reactdyn", TabKind::Reactdyn),
            ("fields", TabKind::Fields),
            ("cad", TabKind::Cad),
            ("meshtoolbox", TabKind::MeshToolbox),
            ("sheetmetal", TabKind::Sheetmetal),
            ("reverse", TabKind::Reverse),
            ("draft2d", TabKind::Draft2d),
            ("render", TabKind::Render),
            ("animate", TabKind::Animate),
            ("springs", TabKind::Springs),
            ("gears", TabKind::Gears),
            ("fasteners", TabKind::Fasteners),
            ("frames", TabKind::Frames),
            ("collision", TabKind::Collision),
            ("piping", TabKind::Piping),
            ("hvac", TabKind::Hvac),
            ("reinforcement", TabKind::Reinforcement),
            ("interior", TabKind::Interior),
            ("geomatics", TabKind::Geomatics),
            ("genetics", TabKind::Genetics),
            ("neuro", TabKind::Neuro),
            ("varianteffect", TabKind::VariantEffect),
        ];
        // Every TEMPLATES kind is covered by the canonical table above.
        assert_eq!(canonical.len(), TabKind::TEMPLATES.len());
        for (id, kind) in canonical {
            assert_eq!(TabKind::from_id(id), Some(kind), "id {id} should map");
            // Case-insensitive + whitespace-tolerant.
            assert_eq!(TabKind::from_id(&id.to_uppercase()), Some(kind));
            assert_eq!(TabKind::from_id(&format!("  {id}  ")), Some(kind));
        }
        // Friendly aliases.
        assert_eq!(TabKind::from_id("mesh"), Some(TabKind::MeshToolbox));
        assert_eq!(TabKind::from_id("variant"), Some(TabKind::VariantEffect));
        // Unknown ids (and Blank, which has no id) map to None.
        assert_eq!(TabKind::from_id("blank"), None);
        assert_eq!(TabKind::from_id("nope"), None);
        assert_eq!(TabKind::from_id(""), None);
    }

    /// The core per-tab-document invariant: there is exactly one
    /// [`WorkspaceDoc`] slot per tab at all times.
    fn assert_docs_aligned(bar: &TabBar) {
        assert_eq!(
            bar.docs.len(),
            bar.tabs.len(),
            "docs must stay index-aligned with tabs"
        );
    }

    #[test]
    fn open_blank_pushes_a_named_blank_and_activates() {
        let mut bar = TabBar::default();
        assert_eq!(bar.active, None);
        let i = bar.open_blank();
        assert_eq!(i, 0);
        assert_eq!(bar.active, Some(0));
        assert_eq!(bar.active_kind(), Some(TabKind::Blank));
        assert_eq!(bar.tabs[0].title, "Untitled 1");
        // A workspace document slot is pushed alongside the tab.
        assert_docs_aligned(&bar);
        // Successive blanks get distinct auto-names.
        bar.open_blank();
        assert_eq!(bar.tabs[1].title, "Untitled 2");
        assert_docs_aligned(&bar);
    }

    #[test]
    fn open_template_pushes_and_activates() {
        let mut bar = TabBar::default();
        let i = bar.open(TabKind::Rocket);
        assert_eq!(i, 0);
        assert_eq!(bar.active, Some(0));
        assert_eq!(bar.active_kind(), Some(TabKind::Rocket));
        assert_eq!(bar.tabs[0].title, "Rocket");
        bar.open(TabKind::Genetics);
        assert_eq!(bar.active, Some(1));
        assert_eq!(bar.tabs.len(), 2);
        assert_docs_aligned(&bar);
    }

    #[test]
    fn close_picks_a_neighbour_then_empties() {
        let mut bar = TabBar::default();
        bar.open(TabKind::Rocket);
        bar.open(TabKind::Cad);
        bar.open(TabKind::Genetics); // active = 2
        assert_docs_aligned(&bar);
        bar.close(2);
        assert_eq!(bar.tabs.len(), 2);
        assert_eq!(bar.active, Some(1)); // clamped to last
        assert_docs_aligned(&bar);
        bar.close(0);
        assert_eq!(bar.active, Some(0));
        assert_docs_aligned(&bar);
        bar.close(0);
        assert_eq!(bar.tabs.len(), 0);
        assert_eq!(bar.active, None);
        assert_docs_aligned(&bar);
    }

    #[test]
    fn sync_active_shows_exactly_one_workbench() {
        let mut app = ValenxApp::default();
        app.tab_bar.open(TabKind::Rocket);
        sync_active(&mut app);
        assert!(app.show_rocket_workbench);
        assert!(!app.show_genetics_workbench);
        assert!(!app.show_cad_workbench);
        assert_eq!(app.active_viewport, ViewportKind::Viewport3D);

        // Switching to a genetics tab hides the rocket and flips the view.
        app.tab_bar.open(TabKind::Genetics);
        sync_active(&mut app);
        assert!(app.show_genetics_workbench);
        assert!(!app.show_rocket_workbench);
        assert_eq!(app.active_viewport, ViewportKind::Viewport2dDna);
    }

    #[test]
    fn blank_tab_shows_no_workbench_but_uses_3d_viewport() {
        let mut app = ValenxApp::default();
        // Pre-set a workbench so we can prove the blank tab clears it.
        app.show_cad_workbench = true;
        app.tab_bar.open_blank();
        sync_active(&mut app);
        assert_eq!(count_shown(&app), 0, "a blank tab opens no workbench");
        assert_eq!(app.active_viewport, ViewportKind::Viewport3D);
    }

    #[test]
    fn sync_with_no_active_hides_everything() {
        let mut app = ValenxApp::default();
        app.show_rocket_workbench = true;
        app.show_cad_workbench = true;
        // No tabs → active is None → everything cleared.
        sync_active(&mut app);
        assert!(!app.show_rocket_workbench);
        assert!(!app.show_cad_workbench);
    }

    /// A minimal, valid [`LoadedMesh`] for exercising the per-tab scene
    /// swap. The mesh itself is empty — we only care that `app.mesh` is
    /// `Some` vs `None` across tab switches.
    fn test_loaded_mesh() -> crate::types::LoadedMesh {
        let mesh = valenx_mesh::Mesh::new("test-tab-scene");
        let quality = valenx_mesh::quality_report(&mesh);
        let aspect_hist =
            valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
        let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
        crate::types::LoadedMesh {
            path: std::path::PathBuf::from("<test>/scene"),
            mesh,
            quality,
            aspect_hist,
            skew_hist,
        }
    }

    #[test]
    fn opening_a_blank_tab_yields_a_fresh_scene_and_switching_back_restores_it() {
        // The headline per-tab-document promise: a tab's loaded geometry is
        // isolated. Open a Rocket tab, load a mesh into the live scene, open
        // a blank tab — the new tab is genuinely empty (mesh == None) — then
        // switch back to the rocket tab and its mesh returns.
        let mut app = ValenxApp::default();

        // Open the first (rocket) tab through the real apply_intent path.
        apply_intent(
            &mut app,
            StripIntent {
                open_template: Some(TabKind::Rocket),
                ..Default::default()
            },
        );
        assert_eq!(app.tab_bar.active, Some(0));
        assert!(app.mesh.is_none(), "a fresh rocket tab starts with no mesh");

        // Load a mesh into the live scene (as if the user imported geometry).
        app.mesh = Some(test_loaded_mesh());
        assert!(app.mesh.is_some());

        // Open a blank tab — its scene must be empty, and the rocket tab's
        // mesh must be parked into docs[0] (not lost, not leaked).
        apply_intent(
            &mut app,
            StripIntent {
                open_blank: true,
                ..Default::default()
            },
        );
        assert_eq!(app.tab_bar.active, Some(1));
        assert!(
            app.mesh.is_none(),
            "the blank tab's scene is fresh — no mesh leaks across tabs"
        );

        // Switch back to the rocket tab: its mesh comes back.
        apply_intent(
            &mut app,
            StripIntent {
                activate: Some(0),
                ..Default::default()
            },
        );
        assert_eq!(app.tab_bar.active, Some(0));
        assert!(
            app.mesh.is_some(),
            "switching back to the rocket tab restores its mesh"
        );
        // And the blank tab's (empty) doc is preserved for next time.
        assert_docs_aligned(&app.tab_bar);
    }

    #[test]
    fn the_workbench_agent_dock_is_per_tab_and_a_new_tab_is_clean() {
        // The headline per-tab-DOCK promise (the user-confirmed bug fix): the
        // "Workbench + Agent" grid belongs to its tab, NOT the whole app.
        // Open a Workbench+Agent unit on tab A (dock on + a dock_tree), open a
        // fresh tab B — B is a CLEAN workspace (dock off, no tree, viewport
        // shown) — then switch back to A and its dock + tree return.
        let mut app = ValenxApp::default();

        // Tab A: a real project tab, then launch a Workbench+Agent unit into
        // it through the same intent the View / tab-strip button uses.
        apply_intent(
            &mut app,
            StripIntent {
                open_template: Some(TabKind::Rocket),
                ..Default::default()
            },
        );
        apply_intent(
            &mut app,
            StripIntent {
                open_wb_agent: true,
                ..Default::default()
            },
        );
        assert_eq!(app.tab_bar.active, Some(0));
        assert!(app.dock_enabled, "tab A has the dock on");
        assert!(app.dock_tree.is_some(), "tab A has a dock tree (the grid)");
        // The agent unit hides the viewport to fill the workspace.
        assert!(app.viewport_hidden, "tab A hid the viewport for the grid");
        // One unit minted; the counter is GLOBAL so this `n` is unique.
        assert_eq!(app.wb_agent_counter, 1);

        // Open a fresh blank tab B — it must start as a CLEAN workspace: the
        // dock off, NO dock tree (not tab A's agent grid), and the 3-D
        // viewport shown. This is the bug the user hit (the grid used to show
        // on every tab).
        apply_intent(
            &mut app,
            StripIntent {
                open_blank: true,
                ..Default::default()
            },
        );
        assert_eq!(app.tab_bar.active, Some(1));
        assert!(!app.dock_enabled, "a new tab starts with the dock OFF");
        assert!(
            app.dock_tree.is_none(),
            "a new tab has NO dock tree — not tab A's agent grid"
        );
        assert!(
            !app.viewport_hidden,
            "a new tab shows the 3-D viewport, not the agent grid"
        );
        assert!(!app.viewport_collapsed, "a new tab's viewport is expanded");
        // The global counter is NOT reset by opening a tab.
        assert_eq!(app.wb_agent_counter, 1, "the unit counter stays global");

        // Switch back to tab A: its dock state + tree return intact.
        apply_intent(
            &mut app,
            StripIntent {
                activate: Some(0),
                ..Default::default()
            },
        );
        assert_eq!(app.tab_bar.active, Some(0));
        assert!(app.dock_enabled, "switching back restores tab A's dock");
        assert!(
            app.dock_tree.is_some(),
            "switching back restores tab A's dock tree (the agent grid)"
        );
        assert!(app.viewport_hidden, "tab A's hidden-viewport state returns");
        assert_docs_aligned(&app.tab_bar);
    }

    #[test]
    fn wb_agent_counter_stays_global_across_tabs_so_channels_never_collide() {
        // The per-unit chat-channel id is `agent:<n>` → counters must be
        // GLOBAL, not per-tab, or two tabs' "Agent 1" would map to one
        // channel. Mint a unit on tab A, open tab B, mint a unit there — the
        // second `n` must be 2 (continued), never reset to 1.
        let mut app = ValenxApp::default();
        apply_intent(
            &mut app,
            StripIntent {
                open_template: Some(TabKind::Rocket),
                ..Default::default()
            },
        );
        apply_intent(
            &mut app,
            StripIntent {
                open_wb_agent: true,
                ..Default::default()
            },
        );
        assert_eq!(app.wb_agent_counter, 1);

        // New tab (clean dock), then mint a unit on it.
        apply_intent(
            &mut app,
            StripIntent {
                open_blank: true,
                ..Default::default()
            },
        );
        assert_eq!(app.wb_agent_counter, 1, "opening a tab never resets it");
        apply_intent(
            &mut app,
            StripIntent {
                open_wb_agent: true,
                ..Default::default()
            },
        );
        assert_eq!(
            app.wb_agent_counter, 2,
            "the second unit gets a globally-unique n (2), not a per-tab reset to 1"
        );
    }

    #[test]
    fn closing_the_active_tab_discards_its_scene_and_restores_the_neighbour() {
        // Closing a tab discards that tab's scene; the neighbour's document
        // (and only it) is installed into the live fields.
        let mut app = ValenxApp::default();
        // Tab 0: rocket with a mesh.
        apply_intent(
            &mut app,
            StripIntent {
                open_template: Some(TabKind::Rocket),
                ..Default::default()
            },
        );
        app.mesh = Some(test_loaded_mesh());
        // Tab 1: blank (empty scene), now active.
        apply_intent(
            &mut app,
            StripIntent {
                open_blank: true,
                ..Default::default()
            },
        );
        assert!(app.mesh.is_none());

        // perform_close(1) (the confirm modal's action): discards tab 1's
        // empty scene and installs tab 0's — whose mesh returns.
        perform_close(&mut app, 1);
        assert_eq!(app.tab_bar.active, Some(0));
        assert_eq!(app.tab_bar.tabs.len(), 1);
        assert_docs_aligned(&app.tab_bar);
        assert!(
            app.mesh.is_some(),
            "closing the blank tab brings back the rocket tab's mesh"
        );

        // Closing the last tab clears the live scene entirely.
        perform_close(&mut app, 0);
        assert_eq!(app.tab_bar.active, None);
        assert!(app.tab_bar.tabs.is_empty());
        assert!(app.mesh.is_none(), "no tabs left → empty live scene");
        assert_docs_aligned(&app.tab_bar);
    }

    #[test]
    fn requesting_a_close_opens_the_confirm_and_does_not_close_yet() {
        // The ✕ / right-click "Close" only *requests* a close: it sets the
        // pending confirm index and leaves the tab (and its scene) intact.
        let mut app = ValenxApp::default();
        apply_intent(
            &mut app,
            StripIntent {
                open_template: Some(TabKind::Rocket),
                ..Default::default()
            },
        );
        app.mesh = Some(test_loaded_mesh());

        apply_intent(
            &mut app,
            StripIntent {
                request_close: Some(0),
                ..Default::default()
            },
        );
        assert_eq!(
            app.tab_close_confirm,
            Some(0),
            "a close request opens the confirm modal"
        );
        assert_eq!(app.tab_bar.tabs.len(), 1, "the tab is not closed yet");
        assert!(
            app.mesh.is_some(),
            "the scene survives an un-confirmed close"
        );
    }

    #[test]
    fn open_wb_agent_intent_launches_a_workbench_agent_pair() {
        // The tab strip's "+ Workbench+Agent" button routes through
        // `StripIntent::open_wb_agent`, which must enable the dock and add a
        // paired workspace/agent unit (counter bumps from 0 → 1).
        let mut app = ValenxApp::default();
        assert_eq!(app.wb_agent_counter, 0);
        assert!(!app.dock_enabled);
        apply_intent(
            &mut app,
            StripIntent {
                open_wb_agent: true,
                ..Default::default()
            },
        );
        assert!(app.dock_enabled, "the pair launcher turns the dock on");
        assert_eq!(app.wb_agent_counter, 1, "one Workbench+Agent unit added");
    }

    #[test]
    fn add_wb_agent_at_intent_places_a_unit_at_the_chosen_spot() {
        // The "+ Workbench+Agent" dropdown routes through
        // `StripIntent::add_wb_agent_at`. Build a 3x2 grid first, then place a
        // new unit into the right end of row 0 — the grid must become [4, 3]
        // and the dock must be enabled.
        let mut app = ValenxApp::default();
        app.open_six_workbench_agents();
        assert_eq!(app.dock_grid_rows(), vec![3, 3]);

        apply_intent(
            &mut app,
            StripIntent {
                add_wb_agent_at: Some(crate::dock_layout::UnitAddTarget::RowEnd(0)),
                ..Default::default()
            },
        );
        assert!(app.dock_enabled, "placing a unit turns the dock on");
        assert_eq!(app.wb_agent_counter, 7, "one more unit minted");
        assert_eq!(app.dock_grid_rows(), vec![4, 3], "row 0 grew");
    }

    #[test]
    fn sanitize_name_strips_path_separators_and_dots() {
        assert_eq!(sanitize_name("boat"), "boat");
        assert_eq!(sanitize_name("my boat-1.2"), "my boat-1.2");
        // Path-escape attempts are neutralised: separators become `_` and a
        // leading/trailing `..` is trimmed, so the result can never be a
        // traversal component or contain a path separator.
        let escaped = sanitize_name("../../etc/passwd");
        assert!(!escaped.contains('/') && !escaped.contains('\\'));
        assert!(!escaped.starts_with("..") && !escaped.ends_with(".."));
        assert_eq!(sanitize_name("a/b\\c"), "a_b_c");
        assert!(!sanitize_name("..").contains(".."));
        // Empty / all-junk collapses to a usable stem.
        assert_eq!(sanitize_name("   "), "untitled");
        assert_eq!(sanitize_name("..."), "untitled");
    }

    #[test]
    fn group_json_round_trip_preserves_tabs_and_active() {
        let mut bar = TabBar::default();
        bar.open(TabKind::Rocket);
        bar.open_blank();
        bar.tabs[1].title = "boat".to_string();
        bar.open(TabKind::Genetics);
        bar.active = Some(1);
        // Set transient edit state to prove it does NOT survive the trip.
        bar.tabs[0].editing = true;
        bar.tabs[0].edit_buf = "scratch".to_string();

        let session = bar.snapshot("my session");
        let json = to_json(&session).expect("serialize");
        let parsed = from_json(&json).expect("deserialize");

        assert_eq!(parsed.name, "my session");
        assert_eq!(parsed.active, Some(1));
        assert_eq!(parsed.tabs.len(), 3);
        assert_eq!(parsed.tabs[0].kind, TabKind::Rocket);
        assert_eq!(parsed.tabs[1].kind, TabKind::Blank);
        assert_eq!(parsed.tabs[1].title, "boat");
        assert_eq!(parsed.tabs[2].kind, TabKind::Genetics);
        // Transient edit state is reset on load.
        assert!(!parsed.tabs[0].editing);
        assert!(parsed.tabs[0].edit_buf.is_empty());
    }

    #[test]
    fn single_tab_json_round_trip() {
        let tab = ProjectTab::new(TabKind::Cad, "bracket");
        let session = SavedSession {
            name: tab.title.clone(),
            tabs: vec![tab.clone()],
            active: Some(0),
        };
        let json = to_json(&session).expect("serialize");
        let parsed = from_json(&json).expect("deserialize");
        assert_eq!(parsed.tabs.len(), 1);
        assert_eq!(parsed.tabs[0].kind, TabKind::Cad);
        assert_eq!(parsed.tabs[0].title, "bracket");
        assert_eq!(parsed.active, Some(0));
    }

    #[test]
    fn restore_clamps_out_of_range_active() {
        let mut bar = TabBar::default();
        bar.open(TabKind::Rocket);
        // Session whose active index points past the end.
        let session = SavedSession {
            name: "bad".to_string(),
            tabs: vec![ProjectTab::new(TabKind::Cad, "a")],
            active: Some(7),
        };
        bar.restore(session);
        assert_eq!(bar.tabs.len(), 1);
        assert_eq!(bar.active, Some(0), "out-of-range active clamps to 0");
        assert_docs_aligned(&bar);

        // Empty session → active None.
        bar.restore(SavedSession {
            name: "empty".to_string(),
            tabs: vec![],
            active: Some(3),
        });
        assert!(bar.tabs.is_empty());
        assert_eq!(bar.active, None);
        assert_docs_aligned(&bar);
    }

    #[test]
    fn append_adds_after_existing_and_activates_first_added() {
        let mut bar = TabBar::default();
        bar.open(TabKind::Rocket);
        let session = SavedSession {
            name: "extra".to_string(),
            tabs: vec![
                ProjectTab::new(TabKind::Cad, "x"),
                ProjectTab::new(TabKind::Fem, "y"),
            ],
            active: Some(1),
        };
        let idx = bar.append(session);
        assert_eq!(idx, Some(1));
        assert_eq!(bar.tabs.len(), 3);
        assert_eq!(bar.active, Some(1));
        assert_eq!(bar.tabs[1].title, "x");
        assert_docs_aligned(&bar);
    }

    /// A unique throwaway directory under the system temp dir, removed when
    /// the returned guard drops. Used to exercise the on-disk save/load path
    /// without touching the process-global state-dir env var (so the tests
    /// stay parallel-safe — Rust runs them on threads in one process).
    struct TempDir(std::path::PathBuf);
    impl TempDir {
        fn new(tag: &str) -> Self {
            let p = std::env::temp_dir().join(format!(
                "valenx-{tag}-{}",
                std::time::SystemTime::now()
                    .duration_since(std::time::UNIX_EPOCH)
                    .unwrap()
                    .as_nanos()
            ));
            std::fs::create_dir_all(&p).unwrap();
            TempDir(p)
        }
        fn path(&self) -> &std::path::Path {
            &self.0
        }
    }
    impl Drop for TempDir {
        fn drop(&mut self) {
            let _ = std::fs::remove_dir_all(&self.0);
        }
    }

    #[test]
    fn save_then_load_group_round_trips_through_disk() {
        let dir = TempDir::new("tabs-group");

        let mut bar = TabBar::default();
        bar.open(TabKind::Rocket);
        bar.open_blank();
        bar.tabs[1].title = "boat".to_string();
        bar.active = Some(1);

        assert!(
            save_session_in(dir.path(), &sanitize_name("trip"), &bar.snapshot("trip")),
            "save should succeed"
        );
        let names = list_json_stems(Some(dir.path().to_path_buf()));
        assert!(
            names.contains(&"trip".to_string()),
            "saved name listed: {names:?}"
        );

        let loaded = load_session_in(dir.path(), "trip").expect("load session");
        assert_eq!(loaded.name, "trip");
        assert_eq!(loaded.tabs.len(), 2);
        assert_eq!(loaded.tabs[0].kind, TabKind::Rocket);
        assert_eq!(loaded.tabs[1].title, "boat");
        assert_eq!(loaded.active, Some(1));
    }

    #[test]
    fn save_then_load_single_tab_round_trips_through_disk() {
        let dir = TempDir::new("tabs-single");

        let tab = ProjectTab::new(TabKind::Genetics, "dna-1");
        let session = SavedSession {
            name: tab.title.clone(),
            tabs: vec![tab.clone()],
            active: Some(0),
        };
        assert!(
            save_session_in(dir.path(), &sanitize_name(&tab.title), &session),
            "save should succeed"
        );
        assert!(list_json_stems(Some(dir.path().to_path_buf())).contains(&"dna-1".to_string()));

        let loaded = load_session_in(dir.path(), "dna-1").expect("load session");
        assert_eq!(loaded.tabs.len(), 1);
        assert_eq!(loaded.tabs[0].kind, TabKind::Genetics);
        assert_eq!(loaded.tabs[0].title, "dna-1");
    }

    #[test]
    fn load_session_in_returns_none_for_missing_file() {
        let dir = TempDir::new("tabs-missing");
        assert!(load_session_in(dir.path(), "nope").is_none());
    }

    /// How many project workbench panels are currently shown.
    fn count_shown(app: &ValenxApp) -> usize {
        [
            app.show_rocket_workbench,
            app.show_engine_workbench,
            app.show_astro_workbench,
            app.show_aero_workbench,
            app.show_gasdynamics_workbench,
            app.show_cfd_workbench,
            app.show_fem_workbench,
            app.show_reactdyn_workbench,
            app.show_fields_workbench,
            app.show_cad_workbench,
            app.show_mesh_toolbox,
            app.show_sheetmetal_workbench,
            app.show_reverse_workbench,
            app.show_draft2d_workbench,
            app.show_render_workbench,
            app.show_animate_workbench,
            app.show_springs_workbench,
            app.show_gears_workbench,
            app.show_fasteners_workbench,
            app.show_frames_workbench,
            app.show_collision_workbench,
            app.show_piping_workbench,
            app.show_hvac_workbench,
            app.show_reinforcement_workbench,
            app.show_interior_workbench,
            app.show_geomatics_workbench,
            app.show_genetics_workbench,
            app.show_neuro_workbench,
            app.show_variant_effect_workbench,
        ]
        .into_iter()
        .filter(|&b| b)
        .count()
    }

    /// Render, once in a headless context, the workbench panel that a tab
    /// of `kind` owns. A blank tab owns no panel, so it's a no-op.
    fn draw_kind(kind: TabKind, app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| match kind {
            TabKind::Blank => {}
            TabKind::Rocket => crate::rocket_workbench::draw_rocket_workbench(app, ctx),
            TabKind::Engine => crate::engine_workbench::draw_engine_workbench(app, ctx),
            TabKind::Astro => crate::astro_workbench::draw_astro_workbench(app, ctx),
            TabKind::Aero => crate::aero_workbench::draw_aero_workbench(app, ctx),
            TabKind::Gasdynamics => {
                crate::gasdynamics_workbench::draw_gasdynamics_workbench(app, ctx)
            }
            TabKind::Cfd => crate::cfd_workbench::draw_cfd_workbench(app, ctx),
            TabKind::Fem => crate::fem_workbench::draw_fem_workbench(app, ctx),
            TabKind::Reactdyn => crate::reactdyn_workbench::draw_reactdyn_workbench(app, ctx),
            TabKind::Fields => crate::fields_workbench::draw_fields_workbench(app, ctx),
            TabKind::Cad => crate::cad_workbench::draw_cad_workbench(app, ctx),
            TabKind::MeshToolbox => crate::mesh_toolbox::draw_mesh_toolbox(app, ctx),
            TabKind::Sheetmetal => crate::sheetmetal_workbench::draw_sheetmetal_workbench(app, ctx),
            TabKind::Reverse => crate::reverse_workbench::draw_reverse_workbench(app, ctx),
            TabKind::Draft2d => crate::draft2d_workbench::draw_draft2d_workbench(app, ctx),
            TabKind::Render => crate::render_workbench::draw_render_workbench(app, ctx),
            TabKind::Animate => crate::animate_workbench::draw_animate_workbench(app, ctx),
            TabKind::Springs => crate::springs_workbench::draw_springs_workbench(app, ctx),
            TabKind::Gears => crate::gears_workbench::draw_gears_workbench(app, ctx),
            TabKind::Fasteners => crate::fasteners_workbench::draw_fasteners_workbench(app, ctx),
            TabKind::Frames => crate::frames_workbench::draw_frames_workbench(app, ctx),
            TabKind::Collision => crate::collision_workbench::draw_collision_workbench(app, ctx),
            TabKind::Piping => crate::piping_workbench::draw_piping_workbench(app, ctx),
            TabKind::Hvac => crate::hvac_workbench::draw_hvac_workbench(app, ctx),
            TabKind::Reinforcement => {
                crate::reinforcement_workbench::draw_reinforcement_workbench(app, ctx)
            }
            TabKind::Interior => crate::interior_workbench::draw_interior_workbench(app, ctx),
            TabKind::Geomatics => crate::geomatics_workbench::draw_geomatics_workbench(app, ctx),
            TabKind::Genetics => crate::genetics_workbench::draw_genetics_workbench(app, ctx),
            TabKind::Neuro => crate::neuro_workbench::draw_neuro_workbench(app, ctx),
            TabKind::VariantEffect => {
                crate::variant_effect_workbench::draw_variant_effect_workbench(app, ctx)
            }
        });
    }

    #[test]
    fn every_template_kind_activates_exactly_one_workbench_and_renders() {
        // The tab system's core promise: opening a tab of any template kind
        // activates exactly that kind's workbench (no leaks, no flags left
        // set) and the workbench renders without panicking.
        for kind in TabKind::TEMPLATES {
            let mut app = ValenxApp::default();
            app.tab_bar.open(kind);
            sync_active(&mut app);
            assert_eq!(
                count_shown(&app),
                1,
                "{kind:?} should show exactly one workbench"
            );
            draw_kind(kind, &mut app);
        }
    }
}

/// Headless render tests for the tab *strip itself* (its menus, the inline
/// rename text field, and the painter-drawn ✕). These mount the real
/// [`draw_tab_strip`] in a windowless [`egui::Context`]; nothing opens an OS
/// window and nothing reaches `rfd::FileDialog` (the strip has no file
/// dialog at all — saving goes straight to the state dir).
#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    /// Run the tab strip once in a headless context.
    fn draw_strip(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_tab_strip(app, ctx);
        });
    }

    #[test]
    fn empty_strip_draws_without_panic() {
        // A fresh app has no tabs; the strip shows the "New tab to begin"
        // hint plus the New-tab / template / Open-saved controls.
        let mut app = ValenxApp::default();
        assert!(app.tab_bar.tabs.is_empty());
        draw_strip(&mut app);
    }

    #[test]
    fn strip_with_mixed_tabs_draws_without_panic() {
        // A blank tab + several template tabs, with one active — exercises
        // the per-tab selectable label, the painter ✕, and the menus.
        let mut app = ValenxApp::default();
        app.tab_bar.open_blank();
        app.tab_bar.open(TabKind::Rocket);
        app.tab_bar.open(TabKind::Genetics);
        sync_active(&mut app);
        draw_strip(&mut app);
    }

    #[test]
    fn strip_renders_inline_rename_field_without_panic() {
        // With a tab in edit mode the strip mounts the rename TextEdit and
        // calls request_focus on it; that path must render headlessly.
        let mut app = ValenxApp::default();
        app.tab_bar.open_blank();
        app.tab_bar.tabs[0].editing = true;
        app.tab_bar.tabs[0].edit_buf = "boat".to_string();
        draw_strip(&mut app);
        // The field renders; with no synthesised key/focus events the edit
        // flag is left as-is (commit happens on Enter / focus loss).
        assert!(app.tab_bar.tabs[0].editing);
    }

    #[test]
    fn strip_renders_close_confirm_modal_without_panic() {
        // With a close pending, the strip mounts the "Close tab?" modal on
        // top of the strip; that path must render headlessly. With no
        // synthesised click on Cancel/Close, the pending index is left set.
        let mut app = ValenxApp::default();
        app.tab_bar.open(TabKind::Rocket);
        sync_active(&mut app);
        app.tab_close_confirm = Some(0);
        draw_strip(&mut app);
        assert_eq!(
            app.tab_close_confirm,
            Some(0),
            "no button was clicked, so the confirm stays open"
        );
        assert_eq!(app.tab_bar.tabs.len(), 1, "nothing closed without confirm");
    }

    #[test]
    fn strip_with_wb_agent_button_draws_without_launching() {
        // The "+ Workbench+Agent" menu button renders on the strip; with no
        // synthesised click on a menu item it must NOT fire the launcher (dock
        // stays off, counter stays 0). A menu_button's body only runs when the
        // popup is open, so a plain frame leaves everything untouched.
        let mut app = ValenxApp::default();
        draw_strip(&mut app);
        assert!(!app.dock_enabled, "no click → the dock stays off");
        assert_eq!(app.wb_agent_counter, 0, "no click → no unit launched");
    }

    #[test]
    fn strip_draws_with_a_populated_grid_for_the_placement_menu() {
        // With a live Workbench+Agent grid present, the tab strip still draws
        // headlessly — the "+ Workbench+Agent" dropdown reads `dock_grid_rows`
        // to build its "Add into a row:" list while the dock_tree is owned.
        // (The menu body itself only runs when the popup is open, but this
        // proves the strip + read path are panic-free with a grid present and
        // that the dock tree is left intact by merely drawing the strip.)
        let mut app = ValenxApp::default();
        app.open_six_workbench_agents();
        assert_eq!(app.dock_grid_rows(), vec![3, 3]);
        draw_strip(&mut app);
        // Drawing the strip neither launched another unit nor disturbed the
        // grid (the dock renders elsewhere, in update.rs).
        assert_eq!(app.wb_agent_counter, 6);
        assert_eq!(app.dock_grid_rows(), vec![3, 3]);
    }

    #[test]
    fn stale_close_confirm_index_clears_safely() {
        // If the pending index points past the end (the tab vanished another
        // way), the modal renderer clears it instead of indexing OOB.
        let mut app = ValenxApp::default();
        app.tab_close_confirm = Some(5);
        draw_strip(&mut app);
        assert_eq!(app.tab_close_confirm, None);
    }
}
