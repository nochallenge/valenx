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
//! ## Scope (v1)
//!
//! v1 is a *view* layer: the heavy per-domain state stays shared (one
//! rocket design, one CAD document, …), so two tabs of the same kind show
//! the same underlying project. Fully independent per-tab documents are a
//! later stage; the model here ([`TabBar`] owning a `Vec` of [`ProjectTab`]
//! plus an `active` index) is already shaped for it.

use eframe::egui;
use serde::{Deserialize, Serialize};

use crate::state_paths::{atomic_write, state_dir};
use crate::viewport_kind::ViewportKind;
use crate::ValenxApp;

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
#[derive(Clone, Serialize, Deserialize)]
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
    /// Index of the active tab in [`Self::tabs`], or `None` when empty.
    pub active: Option<usize>,
    /// Monotonic counter feeding the default "Untitled N" name for blank
    /// tabs, so successive blanks get distinct titles.
    pub blank_counter: usize,
}

impl TabBar {
    /// Open a blank, empty project tab with an auto-generated "Untitled N"
    /// name, make it active, and return its index. This is what the default
    /// `➕ New tab` button does — no workbench is forced open.
    pub fn open_blank(&mut self) -> usize {
        self.blank_counter += 1;
        let title = format!("Untitled {}", self.blank_counter);
        self.tabs.push(ProjectTab::new(TabKind::Blank, title));
        let idx = self.tabs.len() - 1;
        self.active = Some(idx);
        idx
    }

    /// Open a new tab of `kind` (titled with the kind label), make it
    /// active, and return its index. Used by the `＋ from template` menu.
    pub fn open(&mut self, kind: TabKind) -> usize {
        self.tabs.push(ProjectTab::new(kind, kind.label()));
        let idx = self.tabs.len() - 1;
        self.active = Some(idx);
        idx
    }

    /// Close the tab at `idx` and pick a sensible new active tab (the
    /// previous neighbour, or `None` when the strip empties).
    pub fn close(&mut self, idx: usize) {
        if idx >= self.tabs.len() {
            return;
        }
        self.tabs.remove(idx);
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
    /// user reopens a saved group.
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
        self.active = match session.active {
            Some(i) if i < self.tabs.len() => Some(i),
            _ if self.tabs.is_empty() => None,
            _ => Some(0),
        };
    }

    /// Append the tabs from `session` after the current ones (used to
    /// reopen a *single* saved tab without discarding the open set), make
    /// the first appended tab active, and return its index if any.
    pub fn append(&mut self, session: SavedSession) -> Option<usize> {
        if session.tabs.is_empty() {
            return None;
        }
        let first = self.tabs.len();
        for mut t in session.tabs {
            t.editing = false;
            t.edit_buf.clear();
            self.tabs.push(t);
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

/// What a single frame of the tab strip wants to do, accumulated while the
/// read-only borrow of the tab vec is live and applied afterwards.
#[derive(Default)]
struct StripIntent {
    activate: Option<usize>,
    close: Option<usize>,
    open_template: Option<TabKind>,
    open_blank: bool,
    save_tab: Option<usize>,
    save_group: bool,
    open_saved_group: Option<String>,
    open_saved_tab: Option<String>,
    /// Commit an inline rename: (tab index, new title).
    commit_rename: Option<(usize, String)>,
    /// Begin an inline rename of the tab at this index.
    begin_rename: Option<usize>,
}

/// Draw the project-tab strip (a slim panel just below the ribbon) and
/// apply any click this frame (open blank / open template / activate /
/// close / rename / save / open-saved).
pub fn draw_tab_strip(app: &mut ValenxApp, ctx: &egui::Context) {
    let mut intent = StripIntent::default();

    egui::TopBottomPanel::top("valenx_project_tabs").show(ctx, |ui| {
        ui.horizontal(|ui| {
            // Primary: instant blank named project (no forced workbench, no
            // folder dialog).
            if ui
                .button("➕ New tab")
                .on_hover_text("New blank project — name it and start building")
                .clicked()
            {
                intent.open_blank = true;
            }

            // Secondary: start a tab pre-bound to a workbench template. The
            // body is wrapped in `scrollable_menu` so the long category list
            // stays on-screen and scrolls instead of running off the bottom.
            ui.menu_button("＋ from template ▾", |ui| {
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

            // Open a previously-saved tab or group.
            ui.menu_button("Open saved ▾", |ui| {
                crate::menu_ui::scrollable_menu(ui, |ui| {
                    let groups = list_saved_groups();
                    let tabs = list_saved_tabs();
                    if groups.is_empty() && tabs.is_empty() {
                        ui.label(egui::RichText::new("(nothing saved yet)").weak().small());
                    }
                    if !groups.is_empty() {
                        ui.label(egui::RichText::new("Groups (sessions)").small().weak());
                        for name in groups {
                            if ui.button(format!("🗂 {name}")).clicked() {
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
                            if ui.button(format!("📄 {name}")).clicked() {
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
                ui.label(egui::RichText::new("← New tab to begin").weak().small());
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
                    // Right-click context menu: rename / save / close.
                    resp.context_menu(|ui| {
                        if ui.button("Rename").clicked() {
                            intent.begin_rename = Some(i);
                            ui.close_menu();
                        }
                        if ui.button("Save this tab").clicked() {
                            intent.save_tab = Some(i);
                            ui.close_menu();
                        }
                        ui.separator();
                        if ui.button("Close").clicked() {
                            intent.close = Some(i);
                            ui.close_menu();
                        }
                    });
                }

                // Painter-drawn ✕ (reused from the workbench chrome) — never
                // a font-glyph "tofu" box.
                if crate::workbench_chrome::close_x_button(ui, "Close tab").clicked() {
                    intent.close = Some(i);
                }
                ui.separator();
            }
        });
    });

    apply_intent(app, intent);
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
            app.tab_bar.restore(session);
            sync_active(app);
        }
    }
    if let Some(name) = intent.open_saved_tab {
        if let Some(session) = load_saved_tab(&name) {
            app.tab_bar.append(session);
            sync_active(app);
        }
    }
    if let Some(i) = intent.close {
        app.tab_bar.close(i);
        sync_active(app);
    }
    if let Some(i) = intent.activate {
        if i < app.tab_bar.tabs.len() {
            app.tab_bar.active = Some(i);
            sync_active(app);
        }
    }
    if let Some(kind) = intent.open_template {
        app.tab_bar.open(kind);
        sync_active(app);
    }
    if intent.open_blank {
        app.tab_bar.open_blank();
        sync_active(app);
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
    fn open_blank_pushes_a_named_blank_and_activates() {
        let mut bar = TabBar::default();
        assert_eq!(bar.active, None);
        let i = bar.open_blank();
        assert_eq!(i, 0);
        assert_eq!(bar.active, Some(0));
        assert_eq!(bar.active_kind(), Some(TabKind::Blank));
        assert_eq!(bar.tabs[0].title, "Untitled 1");
        // Successive blanks get distinct auto-names.
        bar.open_blank();
        assert_eq!(bar.tabs[1].title, "Untitled 2");
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
    }

    #[test]
    fn close_picks_a_neighbour_then_empties() {
        let mut bar = TabBar::default();
        bar.open(TabKind::Rocket);
        bar.open(TabKind::Cad);
        bar.open(TabKind::Genetics); // active = 2
        bar.close(2);
        assert_eq!(bar.tabs.len(), 2);
        assert_eq!(bar.active, Some(1)); // clamped to last
        bar.close(0);
        assert_eq!(bar.active, Some(0));
        bar.close(0);
        assert_eq!(bar.tabs.len(), 0);
        assert_eq!(bar.active, None);
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

        // Empty session → active None.
        bar.restore(SavedSession {
            name: "empty".to_string(),
            tabs: vec![],
            active: Some(3),
        });
        assert!(bar.tabs.is_empty());
        assert_eq!(bar.active, None);
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
}
