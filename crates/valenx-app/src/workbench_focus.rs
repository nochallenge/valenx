//! Frontend **domain-focus** filter — a pure-UI focus layer that narrows the
//! workbench menus/launcher to a single working domain.
//!
//! valenx ships 44 primary workbench [`TabKind`](crate::project_tabs::TabKind)s
//! spread across a handful of domains (Aerospace, Simulation, CAD & mesh,
//! Machine design, Civil & AEC, Life sciences, …). A user who is, say, doing
//! only bio work sees every domain's tools at once. The domain-focus filter
//! lets them pick a domain so the workbench surfaces show **only** that
//! domain's tools.
//!
//! ## What this is (and is not)
//!
//! This is a **pure-UI focus layer**, nothing more:
//! - Nothing is removed or feature-gated — every workbench stays compiled in,
//!   every solver stays reachable. The filter only hides menu *entries* that
//!   don't match the chosen domain.
//! - There is always an **"All"** escape (the default), which shows everything
//!   exactly as before. "All" is the [`None`] value of
//!   [`ValenxApp::focus_category`](crate::ValenxApp::focus_category).
//! - No compute / registry / category logic is reinvented: the domain of a
//!   workbench is read from the **existing**
//!   [`TabKind::group`](crate::project_tabs::TabKind::group) category fn.
//!
//! ## The category source
//!
//! [`TabKind::group`](crate::project_tabs::TabKind::group) already maps every
//! template workbench to a category string ("Aerospace", "Simulation", …).
//! That string **is** the domain. [`focus_categories`] collects the distinct
//! set (in `TEMPLATES` order, deduplicated) for the selector; [`in_focus`] is
//! the one predicate every filtered surface calls.
//!
//! ## Filtered surfaces
//!
//! Three workbench surfaces honour the focus (all via [`in_focus`]):
//! 1. the **Ctrl+P universal launcher**'s `OpenWorkbenchTab` entries
//!    (`crate::commands::build_visible_commands`);
//! 2. the tab strip's **"From template"** workbench menu
//!    (`crate::project_tabs::draw_tab_strip`);
//! 3. the **View menu**'s primary-workbench toggles (those that map to a
//!    `TabKind`) (`crate::update`).
//!
//! ## Persistence
//!
//! The focus is **in-session only** — it is a transient view preference held on
//! [`ValenxApp`](crate::ValenxApp) and is *not* written to the settings file,
//! so it resets to "All" on relaunch (matching how the other per-session view
//! toggles like the Open-Tabs search behave). Persisting it would mean adding a
//! field to the on-disk `Settings`; that was intentionally left out to keep
//! this a non-breaking, compute-free UI layer.

use crate::project_tabs::TabKind;
use crate::ValenxApp;
use eframe::egui;

/// The accessible Name (and visible prefix) of the focus selector combo box —
/// the string an AI driver / screen reader addresses it by. Kept as a single
/// const so the UI and the accessibility test agree on one literal.
pub const FOCUS_SELECTOR_LABEL: &str = "Focus";

/// The label shown for the "no filter" (show-everything) state — the always-
/// reachable escape. Selecting it sets the focus back to [`None`].
pub const ALL_LABEL: &str = "All";

/// Does a workbench `kind` belong to the focused domain?
///
/// `focus == None` is the "All" escape — every kind is in focus (today's
/// behaviour). `focus == Some(cat)` matches a kind iff its
/// [`TabKind::group`] category equals `cat`. The comparison is on the exact
/// category string the selector offered (which itself came from `group()`), so
/// it can never drift from the category source.
pub fn in_focus(kind: TabKind, focus: Option<&str>) -> bool {
    match focus {
        None => true,
        Some(cat) => kind.group() == cat,
    }
}

/// The distinct workbench-domain categories, in first-seen `TEMPLATES` order
/// and deduplicated — the choices the focus selector offers below "All".
///
/// Reads straight from [`TabKind::group`]; no category list is hand-maintained
/// here, so adding a `TabKind` (with a `group()`) automatically surfaces its
/// domain in the picker.
pub fn focus_categories() -> Vec<&'static str> {
    let mut out: Vec<&'static str> = Vec::new();
    for kind in TabKind::TEMPLATES {
        let g = kind.group();
        if !out.contains(&g) {
            out.push(g);
        }
    }
    out
}

/// The set of template [`TabKind`]s currently visible under `focus` — i.e.
/// every kind for which [`in_focus`] holds. `None` ⇒ all templates.
///
/// This is the canonical "what the filter shows" set the filtered surfaces are
/// expected to agree with; exposed mostly so tests can assert the predicate and
/// the surfaces line up.
pub fn visible_templates(focus: Option<&str>) -> Vec<TabKind> {
    TabKind::TEMPLATES
        .into_iter()
        .filter(|k| in_focus(*k, focus))
        .collect()
}

/// Draw the compact **`Focus: [All v]`** domain picker for the top menu bar.
///
/// A small [`egui::ComboBox`] listing "All" + every [`focus_categories`] entry.
/// Picking one sets [`ValenxApp::focus_category`](crate::ValenxApp::focus_category)
/// (`None` for "All"); the change takes effect immediately on the next frame's
/// menu/launcher draw.
///
/// **AI-drivability / accessibility:** the combo is built with
/// [`egui::ComboBox::from_label`] using [`FOCUS_SELECTOR_LABEL`], so the
/// accessibility (UIA) tree exposes it by that Name — an AI driver / screen
/// reader can find and set it by "Focus" (the standing AI-drivable-first rule).
/// Each option is a uniquely-named [`egui::Ui::selectable_value`].
pub fn focus_selector(app: &mut ValenxApp, ui: &mut egui::Ui) {
    // Current visible caption: the chosen category, or "All" when unfocused.
    let current = app.focus_category.clone();
    let selected_text = current.as_deref().unwrap_or(ALL_LABEL).to_string();

    egui::ComboBox::from_label(FOCUS_SELECTOR_LABEL)
        .selected_text(selected_text)
        .show_ui(ui, |ui| {
            // The always-reachable "All" escape first (sets focus = None).
            if ui
                .selectable_label(current.is_none(), ALL_LABEL)
                .on_hover_text("Show every workbench (no domain filter)")
                .clicked()
            {
                app.focus_category = None;
            }
            ui.separator();
            for cat in focus_categories() {
                let is_sel = current.as_deref() == Some(cat);
                if ui
                    .selectable_label(is_sel, cat)
                    .on_hover_text(format!("Show only {cat} workbenches"))
                    .clicked()
                {
                    app.focus_category = Some(cat.to_string());
                }
            }
        })
        .response
        .on_hover_text(
            "Domain focus: narrow the Tools / template menus and the Ctrl+P \
             launcher to one domain's workbenches. \"All\" shows everything. \
             Pure view filter — nothing is removed.",
        );
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `None` (= "All") is the escape: every template is in focus.
    #[test]
    fn none_focus_shows_all_templates() {
        let visible = visible_templates(None);
        assert_eq!(
            visible.len(),
            TabKind::TEMPLATES.len(),
            "All (None) must show every template kind"
        );
        for k in TabKind::TEMPLATES {
            assert!(in_focus(k, None), "{k:?} must be in focus under All");
        }
    }

    /// `Some("Life sciences")` ⇒ exactly the Bio-category TabKinds — no more,
    /// no fewer. (This is the spec's Bio filter check; the category label for
    /// the bio domain is `TabKind::group()`'s "Life sciences".)
    #[test]
    fn bio_focus_shows_exactly_the_bio_kinds() {
        let bio_cat = TabKind::Genetics.group();
        assert_eq!(bio_cat, "Life sciences", "bio domain label is from group()");

        let visible = visible_templates(Some(bio_cat));
        // The expected set, computed independently straight from group().
        let expected: Vec<TabKind> = TabKind::TEMPLATES
            .into_iter()
            .filter(|k| k.group() == bio_cat)
            .collect();

        assert_eq!(visible, expected, "Bio focus must be exactly the bio kinds");
        // And it really is the known bio workbenches (Genetics/Neuro/Variant/PPI).
        assert!(visible.contains(&TabKind::Genetics));
        assert!(visible.contains(&TabKind::Neuro));
        assert!(visible.contains(&TabKind::VariantEffect));
        assert!(visible.contains(&TabKind::Ppi));
        // A non-bio kind is excluded.
        assert!(!in_focus(TabKind::Rocket, Some(bio_cat)));
        assert!(!visible.contains(&TabKind::Cfd));
        // Every visible kind genuinely carries the bio category.
        for k in &visible {
            assert_eq!(k.group(), bio_cat);
        }
    }

    /// The Aerospace focus likewise yields exactly the Aerospace kinds, and the
    /// predicate agrees with the visible set element-by-element.
    #[test]
    fn aerospace_focus_matches_predicate() {
        let cat = "Aerospace";
        let visible = visible_templates(Some(cat));
        assert!(!visible.is_empty(), "Aerospace must have members");
        assert!(visible.contains(&TabKind::Rocket));
        assert!(visible.contains(&TabKind::Aero));
        for k in TabKind::TEMPLATES {
            assert_eq!(
                in_focus(k, Some(cat)),
                visible.contains(&k),
                "predicate and visible-set must agree for {k:?}"
            );
        }
    }

    /// `focus_categories` is the deduplicated set of every `group()` over the
    /// templates, in first-seen order, and excludes "All".
    #[test]
    fn focus_categories_are_distinct_and_cover_every_group() {
        let cats = focus_categories();
        // No duplicates.
        for (i, a) in cats.iter().enumerate() {
            for b in &cats[i + 1..] {
                assert_ne!(a, b, "duplicate category {a} in focus_categories");
            }
        }
        // "All" is the selector's separate escape, not a category.
        assert!(!cats.contains(&ALL_LABEL));
        // Every template's group is offered.
        for k in TabKind::TEMPLATES {
            assert!(
                cats.contains(&k.group()),
                "category for {k:?} ({}) missing from selector",
                k.group()
            );
        }
        // Every offered category is reachable by at least one template
        // (no phantom categories).
        for cat in &cats {
            assert!(
                TabKind::TEMPLATES.into_iter().any(|k| k.group() == *cat),
                "offered category {cat} has no template"
            );
        }
    }

    /// An unknown focus string matches nothing (defensive: a stale persisted /
    /// hand-set value never accidentally shows a random subset).
    #[test]
    fn unknown_focus_matches_no_template() {
        let visible = visible_templates(Some("Nonexistent domain"));
        assert!(visible.is_empty());
        for k in TabKind::TEMPLATES {
            assert!(!in_focus(k, Some("Nonexistent domain")));
        }
    }

    /// The selector's accessible Name is the documented literal — the string an
    /// AI driver / screen reader uses to find the combo. `focus_selector` builds
    /// it via `ComboBox::from_label(FOCUS_SELECTOR_LABEL)`, so asserting the
    /// const guards the AI-drivable contract without needing a live egui frame.
    #[test]
    fn focus_selector_labelled_by_name_is_stable() {
        assert_eq!(FOCUS_SELECTOR_LABEL, "Focus");
        assert!(!FOCUS_SELECTOR_LABEL.is_empty());
    }

    /// Render the focus selector headlessly in a real egui pass (accesskit
    /// enabled) and confirm the accessibility tree exposes a node named "Focus"
    /// — the on-frame proof of the AI-drivable / `labelled_by` contract.
    /// Mirrors the headless idiom in `widget_naming_tests` /
    /// `uq_workbench::headless_ui_tests`: `enable_accesskit()` → one
    /// `run(RawInput::default())` frame → sweep `accesskit_update.nodes`.
    #[test]
    fn focus_selector_exposes_accessible_name_on_frame() {
        use egui::accesskit::{Node, NodeId};

        let ctx = egui::Context::default();
        ctx.enable_accesskit();
        let out = ctx.run(egui::RawInput::default(), |ctx| {
            egui::CentralPanel::default().show(ctx, |ui| {
                // A default app is enough — the selector only reads
                // `focus_category` (None) and would write it on click.
                let mut app = ValenxApp::default();
                focus_selector(&mut app, ui);
            });
        });
        let nodes: Vec<(NodeId, Node)> = out
            .platform_output
            .accesskit_update
            .expect("accesskit tree is produced when enabled")
            .nodes;

        let found = nodes
            .iter()
            .any(|(_, n)| n.name().is_some_and(|name| name == FOCUS_SELECTOR_LABEL));
        assert!(
            found,
            "focus selector must expose an accessible node named '{FOCUS_SELECTOR_LABEL}'"
        );
    }
}
