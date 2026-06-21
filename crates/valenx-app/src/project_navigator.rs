//! The **project navigator** — an IDE-style "Projects" tree mounted at the
//! top of the Browser panel. It renders the [`crate::project_library`]
//! catalogue (saved projects, organised into folders, pinned, recent) and
//! lets the user search, open, and manage a hundred-plus projects without
//! drilling through the flat `Open saved ▾` menu.
//!
//! ## Layout
//!
//! ```text
//! Projects                         [search box]
//!   ★ Pinned        ← projects with `pinned = true`
//!   🕘 Recent       ← N most-recently-opened
//!   📁 All
//!       <folder 1>  ← CollapsingHeader per folder
//!         project rows…
//!       <folder 2>
//!         project rows…
//!       (unfiled)   ← projects with no folder
//!         project rows…
//! ```
//!
//! Each row is a [`egui::Ui::selectable_label`] showing the project name +
//! a small dim type label (the tab kind). Clicking a row **opens the
//! project as a tab** (reusing the exact tab-reconcile the tab strip uses:
//! `project_tabs::park_active_doc` → `tab_bar.open(kind)` → set the
//! saved title → `project_tabs::install_active_doc`), bumps the
//! project's `last_opened`, and persists the library.
//!
//! ## Intent pattern
//!
//! Like the tab strip ([`crate::project_tabs`]'s `StripIntent`), all row
//! interactions accumulate into a [`NavIntent`] **while the read borrow of
//! the library is live**, and `apply_nav_intent` mutates the library
//! afterwards — so a click that renames/moves/deletes never aliases the
//! borrow it was drawn from. Drag-and-drop reordering is deferred (no
//! `egui_dnd` dep in v1); the row context menu's **Move up / Move down**
//! cover reordering for now.
//!
//! The ★ / 🕘 / 📁 prefix glyphs render in egui 0.28's default font; if a
//! future font swap regresses them to "tofu" boxes, swap the section labels
//! for the plain-text fallbacks ("Pinned" / "Recent" / "All").

use eframe::egui;

use crate::commands::fuzzy_score;
use crate::project_library::ProjectLibrary;

/// Transient navigator UI state owned by [`crate::ValenxApp`]. None of this
/// is persisted — it's scratch state for the search box, the inline rename
/// editor, and the "New folder…" / "Move to ▸ New folder…" name prompt.
#[derive(Default)]
pub struct NavigatorState {
    /// Live contents of the search box; filters every section via
    /// `commands::fuzzy_score`.
    pub search: String,
    /// While `Some(id)`, the project with that id is being renamed inline;
    /// [`Self::rename_buf`] backs the text field.
    pub renaming: Option<String>,
    /// Scratch buffer for the inline project rename.
    pub rename_buf: String,
    /// While `Some(id)`, the folder with that id is being renamed inline.
    pub renaming_folder: Option<String>,
    /// Scratch buffer for the inline folder rename.
    pub folder_rename_buf: String,
    /// Drives the "New folder" name modal. `Some(target)` while the prompt
    /// is open; the payload says what to do once a name is entered.
    pub new_folder_prompt: Option<NewFolderTarget>,
    /// Scratch buffer for the "New folder" name prompt.
    pub new_folder_buf: String,
    /// `true` while the "Delete all projects and folders?" confirmation modal
    /// is open (the destructive whole-library wipe). Cleared on Cancel, or on
    /// confirm right after [`ProjectLibrary::clear_all`] runs.
    pub clear_all_prompt: bool,
}

/// What a freshly-created folder should do once the user names it in the
/// "New folder" prompt.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum NewFolderTarget {
    /// Just create the folder (the "📁 All ▸ New folder…" affordance).
    Create,
    /// Create the folder, then move project `id` into it (the row context
    /// menu's "Move to folder ▸ New folder…" affordance).
    CreateAndMove(String),
}

/// One frame's accumulated navigator actions, applied after the read borrow
/// of the library ends (see the module-level "Intent pattern" note).
#[derive(Default)]
pub struct NavIntent {
    /// Open the project with this id as a tab.
    pub open: Option<String>,
    /// Begin an inline rename of the project with this id.
    pub begin_rename: Option<String>,
    /// Commit an inline project rename: (project id, new name).
    pub commit_rename: Option<(String, String)>,
    /// Cancel any in-progress inline rename.
    pub cancel_rename: bool,
    /// Duplicate the project with this id.
    pub duplicate: Option<String>,
    /// Delete the project with this id.
    pub delete: Option<String>,
    /// Toggle the pinned flag: (project id, new value).
    pub set_pinned: Option<(String, bool)>,
    /// Move project to a folder (or unfiled): (project id, folder id option).
    pub move_to_folder: Option<(String, Option<String>)>,
    /// Reorder a project within its folder: (project id, `up`).
    pub reorder: Option<(String, bool)>,
    /// Open the "New folder" name prompt for the given purpose.
    pub open_new_folder_prompt: Option<NewFolderTarget>,
    /// Open the "Clear all projects" confirmation modal (the destructive
    /// batch-delete of the whole library).
    pub open_clear_all_prompt: bool,
    /// Confirmed: clear the entire library (all projects + folders).
    pub clear_all: bool,
    /// Begin an inline rename of the folder with this id.
    pub begin_folder_rename: Option<String>,
    /// Commit an inline folder rename: (folder id, new name).
    pub commit_folder_rename: Option<(String, String)>,
    /// Delete the folder with this id (its projects fall back to unfiled).
    pub delete_folder: Option<String>,
}

impl NavIntent {
    /// `true` when nothing happened this frame (lets the caller skip the
    /// apply + persist work entirely on idle frames).
    fn is_empty(&self) -> bool {
        self.open.is_none()
            && self.begin_rename.is_none()
            && self.commit_rename.is_none()
            && !self.cancel_rename
            && self.duplicate.is_none()
            && self.delete.is_none()
            && self.set_pinned.is_none()
            && self.move_to_folder.is_none()
            && self.reorder.is_none()
            && self.open_new_folder_prompt.is_none()
            && !self.open_clear_all_prompt
            && !self.clear_all
            && self.begin_folder_rename.is_none()
            && self.commit_folder_rename.is_none()
            && self.delete_folder.is_none()
    }
}

/// Draw the "Projects" navigator at the top of the Browser panel. Renders
/// the search box + the ★ Pinned / 🕘 Recent / 📁 All sections, then applies
/// this frame's [`NavIntent`] (opening/renaming/moving/etc.) and persists
/// the library if anything mutated it.
pub(crate) fn draw_navigator(app: &mut crate::ValenxApp, ui: &mut egui::Ui) {
    let mut intent = NavIntent::default();

    egui::CollapsingHeader::new("Projects")
        .default_open(true)
        .show(ui, |ui| {
            // ── Search box + "new folder" affordance ──────────────────────
            ui.horizontal(|ui| {
                ui.add(
                    egui::TextEdit::singleline(&mut app.nav_state.search)
                        .hint_text("Search projects…")
                        .desired_width(150.0),
                );
                if !app.nav_state.search.is_empty() && ui.small_button("Clear").clicked() {
                    app.nav_state.search.clear();
                }
            });

            // Pinned management buttons (kept ABOVE the scrolling tree, next to
            // the Projects header + search box): create a folder, or wipe the
            // whole library. Plain-ASCII, uniquely-named so the accessibility
            // tree exposes each by Name. "Clear all projects" opens a confirm
            // modal (it is a destructive, can't-be-undone batch delete) and is
            // disabled when the library is already empty.
            ui.horizontal(|ui| {
                if ui
                    .small_button("+ New folder…")
                    .on_hover_text("Create a folder to organise projects")
                    .clicked()
                {
                    intent.open_new_folder_prompt = Some(NewFolderTarget::Create);
                }
                let any = !app.library.projects.is_empty() || !app.library.folders.is_empty();
                if ui
                    .add_enabled(any, egui::Button::new("Clear all projects").small())
                    .on_hover_text("Delete every saved project and folder (cannot be undone)")
                    .clicked()
                {
                    intent.open_clear_all_prompt = true;
                }
            });

            if app.library.projects.is_empty() {
                ui.add_space(2.0);
                ui.label(
                    egui::RichText::new(
                        "(no saved projects yet — right-click a tab → \"Save as project…\")",
                    )
                    .weak()
                    .small(),
                );
            }

            let query = app.nav_state.search.clone();
            let matches = |name: &str| fuzzy_score(&query, name).is_some();

            // Borrow the library immutably for the whole draw; every row
            // interaction records into `intent` instead of mutating here.
            let lib = &app.library;
            let renaming = app.nav_state.renaming.clone();
            let renaming_folder = app.nav_state.renaming_folder.clone();

            // Only the PROJECT TREE (the ★ Pinned / 🕘 Recent / 📁 All
            // sections) scrolls vertically — the Projects header, search box,
            // and the pinned management buttons above stay fixed. With 130+
            // projects the tree would otherwise overrun the Browser panel.
            // `auto_shrink([false, false])` lets the scroll area claim the full
            // remaining panel height/width so the inner CollapsingHeaders keep
            // their normal layout and the overflow scrolls instead of clipping.
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    // ── ★ Pinned ──────────────────────────────────────────────────
                    let pinned: Vec<_> = lib
                        .pinned()
                        .into_iter()
                        .filter(|p| matches(&p.name))
                        .collect();
                    if !pinned.is_empty() {
                        egui::CollapsingHeader::new(format!("\u{2605} Pinned ({})", pinned.len()))
                            .default_open(true)
                            .id_source("nav_pinned")
                            .show(ui, |ui| {
                                for p in &pinned {
                                    project_row(
                                        ui,
                                        p,
                                        lib,
                                        renaming.as_deref(),
                                        &mut app.nav_state.rename_buf,
                                        &mut intent,
                                    );
                                }
                            });
                    }

                    // ── 🕘 Recent ─────────────────────────────────────────────────
                    let recent: Vec<_> = lib
                        .recent(RECENT_LEN)
                        .into_iter()
                        .filter(|p| matches(&p.name))
                        .collect();
                    if !recent.is_empty() {
                        egui::CollapsingHeader::new(format!("\u{1F558} Recent ({})", recent.len()))
                            .default_open(true)
                            .id_source("nav_recent")
                            .show(ui, |ui| {
                                for p in &recent {
                                    project_row(
                                        ui,
                                        p,
                                        lib,
                                        renaming.as_deref(),
                                        &mut app.nav_state.rename_buf,
                                        &mut intent,
                                    );
                                }
                            });
                    }

                    // ── 📁 All (folders + unfiled) ────────────────────────────────
                    egui::CollapsingHeader::new("\u{1F4C1} All")
                        .default_open(true)
                        .id_source("nav_all")
                        .show(ui, |ui| {
                            // (The "+ New folder…" affordance now lives in the pinned
                            // management row above the scrolling tree.)
                            for folder in lib.sorted_folders() {
                                let children: Vec<_> = lib
                                    .projects_in(Some(&folder.id))
                                    .into_iter()
                                    .filter(|p| matches(&p.name))
                                    .collect();
                                // While searching, hide folders with no matching child
                                // so the tree collapses to just the hits.
                                if children.is_empty() && !query.is_empty() {
                                    continue;
                                }
                                // Inline folder-rename editor takes over the header row.
                                if renaming_folder.as_deref() == Some(folder.id.as_str()) {
                                    let resp = ui.add(
                                        egui::TextEdit::singleline(
                                            &mut app.nav_state.folder_rename_buf,
                                        )
                                        .desired_width(150.0)
                                        .id_source(("nav_folder_rename", &folder.id)),
                                    );
                                    resp.request_focus();
                                    if resp.lost_focus() {
                                        intent.commit_folder_rename = Some((
                                            folder.id.clone(),
                                            app.nav_state.folder_rename_buf.clone(),
                                        ));
                                    }
                                    continue;
                                }
                                let header = egui::CollapsingHeader::new(format!(
                                    "{} ({})",
                                    folder.name,
                                    children.len()
                                ))
                                .id_source(("nav_folder", &folder.id))
                                .default_open(query.is_empty().not_default(children.len()));
                                let resp = header.show(ui, |ui| {
                                    for p in &children {
                                        project_row(
                                            ui,
                                            p,
                                            lib,
                                            renaming.as_deref(),
                                            &mut app.nav_state.rename_buf,
                                            &mut intent,
                                        );
                                    }
                                    if children.is_empty() {
                                        ui.label(egui::RichText::new("(empty)").weak().small());
                                    }
                                });
                                // Right-click the folder header → folder actions.
                                resp.header_response.context_menu(|ui| {
                                    if ui.button("Rename folder").clicked() {
                                        intent.begin_folder_rename = Some(folder.id.clone());
                                        ui.close_menu();
                                    }
                                    if ui
                                        .button("Delete folder")
                                        .on_hover_text("Projects inside move to (unfiled)")
                                        .clicked()
                                    {
                                        intent.delete_folder = Some(folder.id.clone());
                                        ui.close_menu();
                                    }
                                });
                            }

                            // Unfiled group (only shown if it has matching projects).
                            let unfiled: Vec<_> = lib
                                .projects_in(None)
                                .into_iter()
                                .filter(|p| matches(&p.name))
                                .collect();
                            if !unfiled.is_empty() {
                                egui::CollapsingHeader::new(format!(
                                    "(unfiled) ({})",
                                    unfiled.len()
                                ))
                                .id_source("nav_unfiled")
                                .default_open(true)
                                .show(ui, |ui| {
                                    for p in &unfiled {
                                        project_row(
                                            ui,
                                            p,
                                            lib,
                                            renaming.as_deref(),
                                            &mut app.nav_state.rename_buf,
                                            &mut intent,
                                        );
                                    }
                                });
                            }
                        }); // 📁 All CollapsingHeader
                }); // ScrollArea::vertical (the scrolling project tree)
        });

    // Apply this frame's actions after the read borrow ends, and persist
    // only when something actually changed.
    if !intent.is_empty() {
        apply_nav_intent(app, intent);
    }

    // The "New folder" name prompt modal (mirrors `draw_close_confirm`).
    draw_new_folder_prompt(app, ui.ctx());
    // The "Clear all projects" destructive-confirm modal.
    draw_clear_all_confirm(app, ui.ctx());
}

/// Number of entries shown in the 🕘 Recent section.
const RECENT_LEN: usize = 6;

/// Tiny extension so a folder defaults open when not searching but only if
/// it has children — keeps the All tree tidy. (A free fn would do; this
/// reads better at the call site.)
trait DefaultOpenExt {
    fn not_default(self, child_count: usize) -> bool;
}
impl DefaultOpenExt for bool {
    fn not_default(self, child_count: usize) -> bool {
        // `self` is `query.is_empty()`: when searching (`false`) always
        // expand so the user sees the hits; when idle, only auto-expand
        // non-empty folders.
        if self {
            child_count > 0
        } else {
            true
        }
    }
}

/// Render a single project row: a [`selectable_label`] with the project name
/// and a dim type tag, plus the right-click context menu. Interactions are
/// recorded into `intent` (never applied here). When the row is the one
/// being renamed, a single-line text editor replaces the label.
fn project_row(
    ui: &mut egui::Ui,
    p: &crate::project_library::SavedProject,
    lib: &ProjectLibrary,
    renaming: Option<&str>,
    rename_buf: &mut String,
    intent: &mut NavIntent,
) {
    if renaming == Some(p.id.as_str()) {
        let resp = ui.add(
            egui::TextEdit::singleline(rename_buf)
                .desired_width(160.0)
                .id_source(("nav_rename", &p.id)),
        );
        resp.request_focus();
        let enter = resp.lost_focus() && ui.input(|i| i.key_pressed(egui::Key::Enter));
        if enter {
            intent.commit_rename = Some((p.id.clone(), rename_buf.clone()));
        } else if resp.lost_focus() {
            // Focus lost without Enter (clicked elsewhere / Esc): commit the
            // edited buffer too so the rename isn't silently dropped.
            intent.commit_rename = Some((p.id.clone(), rename_buf.clone()));
        }
        return;
    }

    let kind_label = p.tab.kind.label();
    let pin_mark = if p.pinned { "\u{2605} " } else { "" };
    let resp = ui
        .horizontal(|ui| {
            let r = ui.selectable_label(false, format!("{pin_mark}{}", p.name));
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                ui.label(
                    egui::RichText::new(kind_label)
                        .small()
                        .color(egui::Color32::from_gray(120)),
                );
            });
            r
        })
        .inner;
    let resp = resp.on_hover_text(format!("{kind_label} — click to open as a tab"));

    if resp.clicked() {
        intent.open = Some(p.id.clone());
    }
    if resp.double_clicked() {
        intent.begin_rename = Some(p.id.clone());
    }

    resp.context_menu(|ui| {
        if ui.button("Open").clicked() {
            intent.open = Some(p.id.clone());
            ui.close_menu();
        }
        if ui.button("Rename").clicked() {
            intent.begin_rename = Some(p.id.clone());
            ui.close_menu();
        }
        if ui.button("Duplicate").clicked() {
            intent.duplicate = Some(p.id.clone());
            ui.close_menu();
        }
        ui.separator();
        // Move to folder ▸ (folders + New folder…).
        ui.menu_button("Move to folder", |ui| {
            crate::menu_ui::scrollable_menu(ui, |ui| {
                let is_unfiled = p.folder.is_none();
                if ui
                    .add_enabled(!is_unfiled, egui::Button::new("(unfiled)"))
                    .clicked()
                {
                    intent.move_to_folder = Some((p.id.clone(), None));
                    ui.close_menu();
                }
                for folder in lib.sorted_folders() {
                    let here = p.folder.as_deref() == Some(folder.id.as_str());
                    if ui
                        .add_enabled(!here, egui::Button::new(&folder.name))
                        .clicked()
                    {
                        intent.move_to_folder = Some((p.id.clone(), Some(folder.id.clone())));
                        ui.close_menu();
                    }
                }
                ui.separator();
                if ui.button("New folder…").clicked() {
                    intent.open_new_folder_prompt =
                        Some(NewFolderTarget::CreateAndMove(p.id.clone()));
                    ui.close_menu();
                }
            });
        });
        if p.pinned {
            if ui.button("Unpin").clicked() {
                intent.set_pinned = Some((p.id.clone(), false));
                ui.close_menu();
            }
        } else if ui.button("Pin").clicked() {
            intent.set_pinned = Some((p.id.clone(), true));
            ui.close_menu();
        }
        ui.separator();
        if ui.button("Move up").clicked() {
            intent.reorder = Some((p.id.clone(), true));
            ui.close_menu();
        }
        if ui.button("Move down").clicked() {
            intent.reorder = Some((p.id.clone(), false));
            ui.close_menu();
        }
        ui.separator();
        // Destructive — tint it like the close-tab action.
        let del = egui::Button::new(
            egui::RichText::new("Delete").color(egui::Color32::from_rgb(220, 80, 80)),
        );
        if ui.add(del).clicked() {
            intent.delete = Some(p.id.clone());
            ui.close_menu();
        }
    });
}

/// Apply one frame's [`NavIntent`] to the library + the navigator's
/// transient state, then persist the library if it changed. Mirrors
/// [`crate::project_tabs`]'s `apply_intent`.
fn apply_nav_intent(app: &mut crate::ValenxApp, intent: NavIntent) {
    let mut dirty = false;

    if let Some(id) = intent.open {
        if open_project_as_tab(app, &id) {
            // last_opened bumped + persisted inside open_project_as_tab.
        }
    }

    if let Some(id) = intent.begin_rename {
        if let Some(p) = app.library.get(&id) {
            app.nav_state.rename_buf = p.name.clone();
            app.nav_state.renaming = Some(id);
        }
    }
    if let Some((id, name)) = intent.commit_rename {
        app.library.rename(&id, &name);
        app.nav_state.renaming = None;
        app.nav_state.rename_buf.clear();
        dirty = true;
    }
    if intent.cancel_rename {
        app.nav_state.renaming = None;
        app.nav_state.rename_buf.clear();
    }

    if let Some(id) = intent.duplicate {
        app.library.duplicate(&id);
        dirty = true;
    }
    if let Some(id) = intent.delete {
        app.library.remove(&id);
        // If we were renaming the now-deleted row, drop the editor.
        if app.nav_state.renaming.as_deref() == Some(id.as_str()) {
            app.nav_state.renaming = None;
            app.nav_state.rename_buf.clear();
        }
        dirty = true;
    }
    if let Some((id, pinned)) = intent.set_pinned {
        app.library.set_pinned(&id, pinned);
        dirty = true;
    }
    if let Some((id, folder)) = intent.move_to_folder {
        app.library.move_to_folder(&id, folder);
        dirty = true;
    }
    if let Some((id, up)) = intent.reorder {
        app.library.reorder(&id, up);
        dirty = true;
    }

    if let Some(target) = intent.open_new_folder_prompt {
        app.nav_state.new_folder_prompt = Some(target);
        app.nav_state.new_folder_buf.clear();
    }
    if intent.open_clear_all_prompt {
        app.nav_state.clear_all_prompt = true;
    }
    if intent.clear_all {
        // Confirmed whole-library wipe: drop every project + folder, then
        // persist the now-empty catalogue. Also clear any in-progress inline
        // rename editors, since their target rows no longer exist.
        app.library.clear_all();
        app.nav_state.renaming = None;
        app.nav_state.rename_buf.clear();
        app.nav_state.renaming_folder = None;
        app.nav_state.folder_rename_buf.clear();
        dirty = true;
    }
    if let Some(id) = intent.begin_folder_rename {
        if let Some(f) = app.library.folders.iter().find(|f| f.id == id) {
            app.nav_state.folder_rename_buf = f.name.clone();
            app.nav_state.renaming_folder = Some(id);
        }
    }
    if let Some((id, name)) = intent.commit_folder_rename {
        app.library.rename_folder(&id, &name);
        app.nav_state.renaming_folder = None;
        app.nav_state.folder_rename_buf.clear();
        dirty = true;
    }
    if let Some(id) = intent.delete_folder {
        app.library.remove_folder(&id);
        if app.nav_state.renaming_folder.as_deref() == Some(id.as_str()) {
            app.nav_state.renaming_folder = None;
            app.nav_state.folder_rename_buf.clear();
        }
        dirty = true;
    }

    if dirty {
        let _ = app.library.save();
    }
}

/// Open the library project `id` as a project tab, reusing the **exact**
/// tab-reconcile the tab strip uses: park the active doc, open a fresh tab of
/// the saved kind, set its title to the saved project name, then install the
/// (empty) doc so the new tab starts clean while the previous tab keeps its
/// geometry. Bumps `last_opened` + persists the library. Returns `true` if a
/// project with that id existed.
///
/// Exposed `pub(crate)` so the command palette's universal launcher
/// ([`crate::commands::dispatch`] of [`crate::commands::CommandKind::OpenSavedProject`])
/// can open a saved project through the **exact** same path the navigator
/// row-click uses, rather than duplicating the reconcile + persist snippet.
pub(crate) fn open_project_as_tab(app: &mut crate::ValenxApp, id: &str) -> bool {
    // Pull the kind + title out before we touch the tab bar (ends the
    // library borrow cleanly).
    let Some((kind, title)) = app.library.get(id).map(|p| (p.tab.kind, p.name.clone())) else {
        return false;
    };

    crate::project_tabs::park_active_doc(app);
    let idx = app.tab_bar.open(kind);
    if let Some(t) = app.tab_bar.tabs.get_mut(idx) {
        t.title = title;
    }
    crate::project_tabs::install_active_doc(app);

    app.library.mark_opened(id);
    let _ = app.library.save();
    true
}

/// Render the "New folder" name prompt modal while
/// [`NavigatorState::new_folder_prompt`] is `Some`. Mirrors
/// [`crate::project_tabs`]'s `draw_close_confirm`: an anchored,
/// non-collapsible window with a name field + Create / Cancel. On Create it
/// adds the folder (and, for the [`NewFolderTarget::CreateAndMove`] variant,
/// moves the originating project into the new folder), then persists.
fn draw_new_folder_prompt(app: &mut crate::ValenxApp, ctx: &egui::Context) {
    let Some(target) = app.nav_state.new_folder_prompt.clone() else {
        return;
    };

    let mut do_create = false;
    let mut do_cancel = false;
    egui::Window::new("New folder")
        .collapsible(false)
        .resizable(false)
        .anchor(egui::Align2::CENTER_CENTER, egui::vec2(0.0, 0.0))
        .show(ctx, |ui| {
            ui.label("Folder name:");
            let resp = ui.add(
                egui::TextEdit::singleline(&mut app.nav_state.new_folder_buf)
                    .desired_width(220.0)
                    .hint_text("e.g. Aerospace"),
            );
            resp.request_focus();
            let enter = resp.lost_focus() && ui.input(|i| i.key_pressed(egui::Key::Enter));
            ui.add_space(8.0);
            ui.horizontal(|ui| {
                if ui.button("Cancel").clicked() {
                    do_cancel = true;
                }
                if ui.button("Create").clicked() || enter {
                    do_create = true;
                }
            });
        });

    if do_cancel {
        app.nav_state.new_folder_prompt = None;
        app.nav_state.new_folder_buf.clear();
    } else if do_create {
        let name = app.nav_state.new_folder_buf.clone();
        let fid = app.library.add_folder(&name);
        if let NewFolderTarget::CreateAndMove(pid) = &target {
            app.library.move_to_folder(pid, Some(fid));
        }
        let _ = app.library.save();
        app.nav_state.new_folder_prompt = None;
        app.nav_state.new_folder_buf.clear();
    }
}

/// Render the **"Delete all projects and folders?"** confirmation modal while
/// [`NavigatorState::clear_all_prompt`] is `true`. Mirrors
/// [`draw_new_folder_prompt`]: an anchored, non-collapsible window. Clearing
/// the library is a destructive, can't-be-undone batch delete, so it is gated
/// behind this explicit confirm. [Cancel] clears the pending flag; [Delete
/// all] routes a [`NavIntent::clear_all`] through [`apply_nav_intent`] (which
/// runs [`ProjectLibrary::clear_all`] + persists). The confirm button is named
/// distinctly from the toolbar's "Clear all projects" so the accessibility
/// tree never shows two same-named controls.
fn draw_clear_all_confirm(app: &mut crate::ValenxApp, ctx: &egui::Context) {
    if !app.nav_state.clear_all_prompt {
        return;
    }
    let n_projects = app.library.projects.len();
    let n_folders = app.library.folders.len();

    let mut do_clear = false;
    let mut do_cancel = false;
    egui::Window::new("Clear all projects?")
        .collapsible(false)
        .resizable(false)
        .anchor(egui::Align2::CENTER_CENTER, egui::vec2(0.0, 0.0))
        .show(ctx, |ui| {
            ui.label(format!(
                "Delete all {n_projects} projects and {n_folders} folders? This cannot be undone."
            ));
            ui.add_space(8.0);
            ui.horizontal(|ui| {
                if ui.button("Cancel").clicked() {
                    do_cancel = true;
                }
                // Red-ish destructive action, uniquely named.
                let del = egui::Button::new(
                    egui::RichText::new("Delete all").color(egui::Color32::from_rgb(220, 80, 80)),
                );
                if ui.add(del).clicked() {
                    do_clear = true;
                }
            });
        });

    if do_cancel {
        app.nav_state.clear_all_prompt = false;
    } else if do_clear {
        apply_nav_intent(
            app,
            NavIntent {
                clear_all: true,
                ..Default::default()
            },
        );
        app.nav_state.clear_all_prompt = false;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::project_library::ProjectLibrary;
    use crate::project_tabs::{ProjectTab, TabKind};

    fn tab(kind: TabKind, title: &str) -> ProjectTab {
        ProjectTab {
            kind,
            title: title.to_string(),
            group: None,
            workbench_kind: None,
            editing: false,
            edit_buf: String::new(),
        }
    }

    /// Apply just the library-mutating parts of a `NavIntent` against a bare
    /// library (no `ValenxApp`), so the intent→library wiring is testable
    /// without standing up egui / the whole app. This mirrors the
    /// library-side branches of `apply_nav_intent` exactly.
    fn apply_to_library(lib: &mut ProjectLibrary, intent: &NavIntent) {
        if let Some((id, name)) = &intent.commit_rename {
            lib.rename(id, name);
        }
        if let Some(id) = &intent.duplicate {
            lib.duplicate(id);
        }
        if let Some(id) = &intent.delete {
            lib.remove(id);
        }
        if let Some((id, pinned)) = &intent.set_pinned {
            lib.set_pinned(id, *pinned);
        }
        if let Some((id, folder)) = &intent.move_to_folder {
            lib.move_to_folder(id, folder.clone());
        }
        if let Some((id, up)) = &intent.reorder {
            lib.reorder(id, *up);
        }
        if intent.clear_all {
            lib.clear_all();
        }
    }

    #[test]
    fn nav_intent_empty_by_default() {
        assert!(NavIntent::default().is_empty());
    }

    #[test]
    fn nav_intent_set_pinned_applies_to_library() {
        let mut lib = ProjectLibrary::default();
        let id = lib.add_project(tab(TabKind::Rocket, "P"), None);
        assert!(!lib.get(&id).unwrap().pinned);

        let intent = NavIntent {
            set_pinned: Some((id.clone(), true)),
            ..Default::default()
        };
        assert!(!intent.is_empty());
        apply_to_library(&mut lib, &intent);
        assert!(lib.get(&id).unwrap().pinned);
    }

    #[test]
    fn nav_intent_rename_and_delete_apply() {
        let mut lib = ProjectLibrary::default();
        let id = lib.add_project(tab(TabKind::Cad, "Old"), None);

        apply_to_library(
            &mut lib,
            &NavIntent {
                commit_rename: Some((id.clone(), "Renamed".into())),
                ..Default::default()
            },
        );
        assert_eq!(lib.get(&id).unwrap().name, "Renamed");

        apply_to_library(
            &mut lib,
            &NavIntent {
                delete: Some(id.clone()),
                ..Default::default()
            },
        );
        assert!(lib.get(&id).is_none());
    }

    #[test]
    fn nav_intent_move_to_folder_applies() {
        let mut lib = ProjectLibrary::default();
        let fid = lib.add_folder("Group");
        let pid = lib.add_project(tab(TabKind::Fem, "Beam"), None);
        apply_to_library(
            &mut lib,
            &NavIntent {
                move_to_folder: Some((pid.clone(), Some(fid.clone()))),
                ..Default::default()
            },
        );
        assert_eq!(lib.get(&pid).unwrap().folder.as_deref(), Some(fid.as_str()));
    }

    #[test]
    fn nav_intent_reorder_applies() {
        let mut lib = ProjectLibrary::default();
        let a = lib.add_project(tab(TabKind::Rocket, "A"), None);
        let b = lib.add_project(tab(TabKind::Cad, "B"), None);
        apply_to_library(
            &mut lib,
            &NavIntent {
                reorder: Some((b.clone(), true)),
                ..Default::default()
            },
        );
        let names: Vec<String> = lib
            .projects_in(None)
            .iter()
            .map(|p| p.name.clone())
            .collect();
        assert_eq!(names, vec!["B", "A"]);
        let _ = a;
    }

    #[test]
    fn nav_intent_clear_all_is_not_empty_and_empties_library() {
        // The clear-all and open-prompt intents must register as "not empty"
        // so `draw_navigator` actually runs the apply (and persists).
        assert!(!NavIntent {
            clear_all: true,
            ..Default::default()
        }
        .is_empty());
        assert!(!NavIntent {
            open_clear_all_prompt: true,
            ..Default::default()
        }
        .is_empty());

        // Routing a `clear_all` intent through the library-side apply wipes
        // every project + folder.
        let mut lib = ProjectLibrary::default();
        let f = lib.add_folder("Group");
        let _a = lib.add_project(tab(TabKind::Rocket, "A"), Some(f));
        let _b = lib.add_project(tab(TabKind::Cad, "B"), None);
        apply_to_library(
            &mut lib,
            &NavIntent {
                clear_all: true,
                ..Default::default()
            },
        );
        assert!(lib.projects.is_empty());
        assert!(lib.folders.is_empty());
    }

    #[test]
    fn default_open_ext_expands_when_searching() {
        // query.is_empty() == false (searching): always expand.
        assert!(false.not_default(0));
        assert!(false.not_default(3));
        // query.is_empty() == true (idle): only expand non-empty folders.
        assert!(!true.not_default(0));
        assert!(true.not_default(2));
    }
}
