//! Chrome-style **project tabs** — an open-many, switch-between strip of
//! project workspaces, each bound to a project *kind*.
//!
//! valenx's domain tools are independent right-dock workbench panels,
//! toggled from the View menu. The tab strip is a thin navigation layer
//! over them: each tab owns one [`TabKind`], **activating** a tab shows
//! that kind's workbench and hides the others (so the user works one
//! project at a time, like browser tabs), the `➕` button opens a new
//! tab of any kind, and `✕` closes one.
//!
//! The strip is **additive and non-breaking**: a fresh app starts with
//! zero tabs and the existing default layout untouched. Tab mode only
//! engages once the user opens the first tab.
//!
//! ## Scope (v1)
//!
//! v1 is a *view* layer: the heavy per-domain state stays shared (one
//! rocket design, one CAD document, …), so two tabs of the same kind
//! show the same underlying project. Fully independent per-tab documents
//! are a later stage; the model here ([`TabBar`] owning a `Vec` of
//! [`ProjectTab`] plus an `active` index) is already shaped for it.

use eframe::egui;

use crate::viewport_kind::ViewportKind;
use crate::ValenxApp;

/// A project kind a tab can hold. Each maps to exactly one primary
/// workbench panel (the `show_*` flag it drives on [`ValenxApp`]).
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum TabKind {
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
    /// Every kind, in `➕`-menu order (grouped via [`Self::group`]).
    pub const ALL: [TabKind; 29] = [
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

    /// Group header shown in the `➕` new-tab menu.
    pub fn group(self) -> &'static str {
        match self {
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
    fn show(self, app: &mut ValenxApp) {
        match self {
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
    /// view; everything else is the 3D viewport).
    fn viewport(self) -> ViewportKind {
        match self {
            TabKind::Genetics => ViewportKind::Viewport2dDna,
            _ => ViewportKind::Viewport3D,
        }
    }
}

/// One open project tab: its kind plus a user-facing title.
pub struct ProjectTab {
    /// The project kind this tab hosts.
    pub kind: TabKind,
    /// Title shown on the tab (defaults to the kind label).
    pub title: String,
}

/// The project-tab strip state, owned by [`ValenxApp`].
#[derive(Default)]
pub struct TabBar {
    /// Open tabs, left to right.
    pub tabs: Vec<ProjectTab>,
    /// Index of the active tab in [`Self::tabs`], or `None` when empty.
    pub active: Option<usize>,
}

impl TabBar {
    /// Open a new tab of `kind`, make it active, and return its index.
    pub fn open(&mut self, kind: TabKind) -> usize {
        self.tabs.push(ProjectTab {
            kind,
            title: kind.label().to_string(),
        });
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
}

/// Hide every project workbench panel. The active tab (if any) then
/// re-shows exactly one via [`TabKind::show`].
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
/// switch the viewport to match. With no active tab, everything stays
/// hidden (the user closed the last tab).
pub fn sync_active(app: &mut ValenxApp) {
    let kind = app.tab_bar.active_kind();
    clear_all_workbenches(app);
    if let Some(kind) = kind {
        kind.show(app);
        app.active_viewport = kind.viewport();
    }
}

/// Draw the project-tab strip (a slim panel just below the ribbon) and
/// apply any click this frame (activate / close / open-new).
pub fn draw_tab_strip(app: &mut ValenxApp, ctx: &egui::Context) {
    let mut to_activate: Option<usize> = None;
    let mut to_close: Option<usize> = None;
    let mut to_open: Option<TabKind> = None;

    egui::TopBottomPanel::top("valenx_project_tabs").show(ctx, |ui| {
        ui.horizontal(|ui| {
            // ➕ new-tab menu, grouped by domain.
            ui.menu_button("➕ New tab", |ui| {
                let mut last_group = "";
                for kind in TabKind::ALL {
                    let group = kind.group();
                    if group != last_group {
                        if !last_group.is_empty() {
                            ui.separator();
                        }
                        ui.label(egui::RichText::new(group).small().weak());
                        last_group = group;
                    }
                    if ui.button(kind.label()).clicked() {
                        to_open = Some(kind);
                        ui.close_menu();
                    }
                }
            });
            ui.separator();

            if app.tab_bar.tabs.is_empty() {
                ui.label(
                    egui::RichText::new("← open a project tab to begin")
                        .weak()
                        .small(),
                );
            }

            let active = app.tab_bar.active;
            for (i, tab) in app.tab_bar.tabs.iter().enumerate() {
                let selected = active == Some(i);
                if ui
                    .selectable_label(selected, &tab.title)
                    .on_hover_text(tab.kind.group())
                    .clicked()
                {
                    to_activate = Some(i);
                }
                if ui.small_button("✕").on_hover_text("Close tab").clicked() {
                    to_close = Some(i);
                }
                ui.separator();
            }
        });
    });

    // Apply this frame's intent after the read-only borrows above end.
    // Order: close, then activate, then open — at most one fires per
    // frame in practice, but each leaves `active` consistent.
    if let Some(i) = to_close {
        app.tab_bar.close(i);
        sync_active(app);
    }
    if let Some(i) = to_activate {
        if i < app.tab_bar.tabs.len() {
            app.tab_bar.active = Some(i);
            sync_active(app);
        }
    }
    if let Some(kind) = to_open {
        app.tab_bar.open(kind);
        sync_active(app);
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod tests {
    use super::*;

    #[test]
    fn all_kinds_are_unique_and_grouped() {
        // No duplicate kinds in ALL.
        for (i, a) in TabKind::ALL.iter().enumerate() {
            for b in &TabKind::ALL[i + 1..] {
                assert_ne!(a, b, "duplicate kind in ALL: {a:?}");
            }
        }
        // Every kind has a non-empty label and group.
        for k in TabKind::ALL {
            assert!(!k.label().is_empty());
            assert!(!k.group().is_empty());
        }
    }

    #[test]
    fn open_pushes_and_activates() {
        let mut bar = TabBar::default();
        assert_eq!(bar.active, None);
        let i = bar.open(TabKind::Rocket);
        assert_eq!(i, 0);
        assert_eq!(bar.active, Some(0));
        assert_eq!(bar.active_kind(), Some(TabKind::Rocket));
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
    fn sync_with_no_active_hides_everything() {
        let mut app = ValenxApp::default();
        app.show_rocket_workbench = true;
        app.show_cad_workbench = true;
        // No tabs → active is None → everything cleared.
        sync_active(&mut app);
        assert!(!app.show_rocket_workbench);
        assert!(!app.show_cad_workbench);
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
    /// of `kind` owns.
    fn draw_kind(kind: TabKind, app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| match kind {
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
    fn every_tab_kind_activates_exactly_one_workbench_and_renders() {
        // The tab system's core promise: opening a tab of any kind activates
        // exactly that kind's workbench (no leaks, no flags left set) and the
        // workbench renders without panicking (no unreachable stub).
        for kind in TabKind::ALL {
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
