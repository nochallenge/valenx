//! Welcome / start page rendered in the central viewport area when
//! no project is loaded.
//!
//! Design pattern: VS Code's welcome tab + Fusion 360's start page.
//! The placeholder "No geometry loaded — drop a .stl here" surface
//! was the user's first impression for v0.1.0-alpha.1, which (a)
//! didn't lead with the value prop (native engines, no external
//! installs) and (b) hid the Create-New-Project / Open-Project
//! affordances behind the File menu and the Ctrl+N shortcut. This
//! module replaces that placeholder with a proper welcome screen.
//!
//! ## Layout
//!
//! Vertically-centred in the central viewport:
//!
//! - **Title row** — wordmark + version + tagline
//! - **Action row** — two big cards: New Project, Open Project
//! - **Recent projects** — most-recent-first list, click to re-open
//! - **Native-engines tagline** — the value-prop bumper
//! - **Quick links** — README / CHANGELOG / docs in the host file
//!   browser, with the repository URL as a fallback
//!
//! The render function is pure with respect to disk: it returns a
//! [`LandingAction`] enum the host dispatches against. Side-effecting
//! actions (file-browser launches, project loads) all funnel through
//! the host so this module stays testable.
//!
//! ## Why not a full panel?
//!
//! The landing page only shows when `self.project.is_none()`. The
//! moment a project loads, the central viewport switches back to its
//! normal STL/mesh rendering path. There's no toggle / menu item —
//! the page exists strictly as an empty-state for the central panel.

use std::path::{Path, PathBuf};

use eframe::egui;

/// Side-effect-bearing actions the welcome page asks its host to
/// perform. Returned from [`render`] each frame; `None` means "the
/// user just hovered / scrolled, do nothing".
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum LandingAction {
    /// User clicked the New Project card — open the same modal Ctrl+N
    /// opens.
    NewProject,
    /// User clicked the Open Project card — open the same folder
    /// picker Ctrl+O opens.
    OpenProject,
    /// User clicked a recent-project row — load that path. Mirrors
    /// the existing `app.load_project(path)` flow.
    OpenRecent(PathBuf),
    /// User clicked a recent-project row whose path no longer exists
    /// on disk. The host should drop the entry from
    /// `settings.recent_projects` + persist; the landing page renders
    /// a one-line "(removed missing project)" hint in the next frame.
    DropMissingRecent(PathBuf),
    /// User clicked one of the quick-action links — open `path_or_url`
    /// in the host's file browser (for filesystem paths) or default
    /// browser (for URLs). The host routes both cases through
    /// [`open_link`] which falls back to a platform-appropriate
    /// launcher.
    OpenLink(String),
}

impl LandingAction {
    /// The three quick-link targets the welcome page surfaces by
    /// default (Learn more / What's new / Documentation). Resolves
    /// each via `local_doc_or_url` so a source checkout opens the
    /// in-tree markdown while an installed build falls back to the
    /// repo URL. Exposed for unit testing — the host wires the same
    /// resolution inside [`render`].
    pub fn default_links_resolve(repo_url: &str) -> [(&'static str, String); 3] {
        let learn_more = local_doc_or_url("README.md", repo_url);
        let whats_new = local_doc_or_url(
            "CHANGELOG.md",
            &format!("{repo_url}/blob/master/CHANGELOG.md"),
        );
        let docs = local_doc_or_url(
            "docs/INSTALLER.md",
            &format!("{repo_url}/blob/master/docs/INSTALLER.md"),
        );
        [
            ("Learn more", learn_more),
            ("What's new", whats_new),
            ("Documentation", docs),
        ]
    }
}

/// Render the welcome page into the current `ui`'s available rect.
/// Returns `Some(action)` when the user clicked something this frame,
/// `None` otherwise.
///
/// `recent_projects` is borrowed from `settings.recent_projects` so
/// the host doesn't have to clone — the page only reads the list.
/// `version` is the app's `CARGO_PKG_VERSION` (passed in so this
/// module doesn't need its own `env!`).
/// `repo_url` is the workspace `repository` field, used as the
/// fallback "Learn more" target when README isn't on disk.
/// `inline_message` is an optional one-line notice rendered next to
/// the recent-projects list — used to confirm "(removed missing
/// project from recents)" after a `DropMissingRecent` action lands.
pub fn render(
    ui: &mut egui::Ui,
    version: &str,
    recent_projects: &[PathBuf],
    repo_url: &str,
    inline_message: Option<&str>,
) -> Option<LandingAction> {
    let mut action: Option<LandingAction> = None;

    // Centred column with a sensible max width so wide windows
    // don't stretch the cards across the screen. We pad the top
    // so the content sits in the upper-middle of the viewport
    // (where the eye lands first) rather than dead-centre.
    let available = ui.available_size();
    let top_pad = (available.y * 0.08).clamp(16.0, 96.0);
    ui.add_space(top_pad);

    let max_content_width = 720.0_f32.min(available.x - 32.0).max(320.0);

    ui.vertical_centered(|ui| {
        ui.set_max_width(max_content_width);

        // --- Title row -----------------------------------------------
        ui.heading(
            egui::RichText::new(format!("Valenx {version}"))
                .size(28.0)
                .strong(),
        );
        ui.add_space(2.0);
        ui.label(
            egui::RichText::new("Native open-source simulation suite")
                .size(14.0)
                .color(egui::Color32::from_gray(170)),
        );
        ui.add_space(28.0);

        // --- Action row (two big cards) ------------------------------
        // Horizontal layout, each card a fixed-width button. We use
        // `egui::Frame` to draw the card outline so the click target
        // is the whole rectangle, not just the button label.
        ui.horizontal(|ui| {
            // egui's horizontal layout doesn't centre its children;
            // a leading spacer balances the right-side natural gap
            // so the two cards land symmetrically inside the
            // centred column.
            let card_width = 220.0_f32;
            let total_cards_width = card_width * 2.0 + 16.0;
            let leading = ((max_content_width - total_cards_width) / 2.0).max(0.0);
            ui.add_space(leading);

            if action_card(ui, card_width, "+  New Project", "Ctrl+N").clicked() {
                action = Some(LandingAction::NewProject);
            }
            ui.add_space(16.0);
            if action_card(ui, card_width, "Open Project", "Ctrl+O").clicked() {
                action = Some(LandingAction::OpenProject);
            }
        });
        ui.add_space(28.0);

        // --- Recent projects -----------------------------------------
        ui.horizontal(|ui| {
            ui.label(
                egui::RichText::new("Recent projects")
                    .size(13.0)
                    .strong()
                    .color(egui::Color32::from_gray(200)),
            );
        });
        ui.add_space(4.0);
        ui.separator();
        ui.add_space(2.0);

        if recent_projects.is_empty() {
            ui.label(
                egui::RichText::new("(no recent projects yet — start by creating one)")
                    .size(12.0)
                    .italics()
                    .color(egui::Color32::from_gray(140)),
            );
        } else {
            // Vertical list, one row per project. We render each row
            // as a selectable label so hover + click feel right; the
            // path tail (last folder) is the headline, with the full
            // path as a weaker subtitle on the same row. Missing-on-
            // disk entries dispatch a `DropMissingRecent` instead of
            // `OpenRecent` so the host can prune the entry + show the
            // inline confirmation message rather than failing the
            // load with a top-status-bar error.
            for path in recent_projects {
                if let Some((picked, exists)) = recent_project_row(ui, path) {
                    action = Some(if exists {
                        LandingAction::OpenRecent(picked)
                    } else {
                        LandingAction::DropMissingRecent(picked)
                    });
                }
            }
        }
        if let Some(msg) = inline_message {
            ui.add_space(4.0);
            ui.label(
                egui::RichText::new(msg)
                    .size(11.0)
                    .italics()
                    .color(egui::Color32::from_rgb(220, 170, 80)),
            );
        }

        ui.add_space(28.0);

        // --- Native-engines tagline ----------------------------------
        // This is the BIG messaging shift. Replaces the previous
        // "drop a .stl" empty-state hint with the value prop: native
        // implementations, no external installs.
        ui.label(
            egui::RichText::new("Native engines included — no external installs required.")
                .size(13.0)
                .strong()
                .color(egui::Color32::from_rgb(140, 200, 230)),
        );
        ui.add_space(4.0);
        ui.label(
            egui::RichText::new(
                "Valenx ships native Rust implementations of CFD, FEA, molecular \
                 dynamics, quantum chemistry, RNA folding, protein structure \
                 prediction, and more.",
            )
            .size(12.0)
            .color(egui::Color32::from_gray(170)),
        );
        ui.add_space(20.0);

        // --- Quick links ---------------------------------------------
        // Three text links. Each routes through `LandingAction::OpenLink`
        // so the host's `open_link` helper handles the URL-vs-path
        // dispatch.
        ui.horizontal_wrapped(|ui| {
            // Prefer the local README / CHANGELOG / docs when present —
            // they ship inside the install — and fall back to the
            // public repo URL when this is a portable / source-only
            // checkout where the markdown isn't laid out next to the
            // binary. Resolved through `LandingAction::default_links_resolve`
            // so the unit tests exercise the same code path the render
            // does.
            let [(learn_label, learn_more), (whats_label, whats_new), (docs_label, docs)] =
                LandingAction::default_links_resolve(repo_url);

            if link_row(ui, learn_label).clicked() {
                action = Some(LandingAction::OpenLink(learn_more));
            }
            ui.label(egui::RichText::new("  ·  ").color(egui::Color32::from_gray(110)));
            if link_row(ui, whats_label).clicked() {
                action = Some(LandingAction::OpenLink(whats_new));
            }
            ui.label(egui::RichText::new("  ·  ").color(egui::Color32::from_gray(110)));
            if link_row(ui, docs_label).clicked() {
                action = Some(LandingAction::OpenLink(docs));
            }
        });
    });

    action
}

/// Render a single big action card. Returns the egui response so the
/// caller can check `.clicked()`.
///
/// The card is a hand-drawn rectangle (filled background + outline
/// stroke) plus a centred title + hint label drawn through a child
/// UI. Hovering brightens the fill + recolours the stroke so users
/// get the "this is a button" affordance.
fn action_card(ui: &mut egui::Ui, width: f32, title: &str, hint: &str) -> egui::Response {
    let height = 96.0;
    let (rect, response) = ui.allocate_exact_size(egui::vec2(width, height), egui::Sense::click());
    let (fill, stroke) = if response.hovered() {
        (
            egui::Color32::from_gray(52),
            egui::Stroke::new(1.5, egui::Color32::from_rgb(120, 170, 220)),
        )
    } else {
        (
            egui::Color32::from_gray(36),
            egui::Stroke::new(1.0, egui::Color32::from_gray(60)),
        )
    };
    let rounding = egui::Rounding::same(6.0);
    let painter = ui.painter_at(rect);
    painter.rect_filled(rect, rounding, fill);
    painter.rect_stroke(rect, rounding, stroke);
    // Draw the label + hint as two centred lines. Avoids needing a
    // child Ui (which would require non-deprecated `child_ui` plumbing
    // that varies across egui versions); the painter API is stable.
    let center = rect.center();
    painter.text(
        center - egui::vec2(0.0, 12.0),
        egui::Align2::CENTER_CENTER,
        title,
        egui::FontId::proportional(18.0),
        egui::Color32::from_gray(230),
    );
    painter.text(
        center + egui::vec2(0.0, 18.0),
        egui::Align2::CENTER_CENTER,
        hint,
        egui::FontId::proportional(12.0),
        egui::Color32::from_gray(160),
    );
    response
}

/// Render one recent-project row. Returns `Some((path, exists))` if
/// the user clicked it this frame — `exists == false` signals the
/// host should prune the entry from settings rather than load it.
fn recent_project_row(ui: &mut egui::Ui, path: &Path) -> Option<(PathBuf, bool)> {
    let label = path
        .file_name()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_else(|| path.display().to_string());
    let full = path.display().to_string();
    // Missing-on-disk projects render weaker so the user can spot a
    // moved / deleted entry. Clicking a missing entry tells the host
    // to drop it from the recents list (and the landing page shows
    // an inline confirmation on the next frame) — clicking an
    // existing entry routes through the normal load path.
    let exists = path.exists();
    let title_color = if exists {
        egui::Color32::from_gray(220)
    } else {
        egui::Color32::from_gray(140)
    };
    let response = ui.horizontal(|ui| {
        let title = egui::RichText::new(format!("> {label}")).color(title_color);
        let inner = ui.selectable_label(false, title);
        // Right-align the full path so the row reads "name . . . path".
        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
            let mut path_text = egui::RichText::new(full)
                .size(11.0)
                .color(egui::Color32::from_gray(120));
            if !exists {
                path_text = path_text.italics();
            }
            ui.label(path_text);
        });
        inner
    });
    let inner = response.inner;
    if !exists {
        // The hover-text confirms why the entry is greyed out before
        // the user commits to the click.
        inner
            .on_hover_text("Project path no longer exists on disk — click to remove from recents")
            .clicked()
            .then(|| (path.to_path_buf(), false))
    } else if inner.clicked() {
        Some((path.to_path_buf(), true))
    } else {
        None
    }
}

/// Render a clickable text link in the welcome page. Uses egui's
/// link styling so the hover / click affordance reads correctly.
fn link_row(ui: &mut egui::Ui, label: &str) -> egui::Response {
    ui.link(
        egui::RichText::new(label)
            .size(12.0)
            .color(egui::Color32::from_rgb(140, 180, 220)),
    )
}

/// If a path relative to the binary's working directory exists,
/// return its absolute string; otherwise return the supplied URL.
/// Used by the Learn-more / What's-new / Documentation links so a
/// source checkout opens the in-tree markdown while an installed
/// build falls back to the public repo URL.
fn local_doc_or_url(rel: &str, fallback_url: &str) -> String {
    if let Ok(abs) = std::env::current_dir().map(|d| d.join(rel)) {
        if abs.is_file() {
            return abs.display().to_string();
        }
    }
    fallback_url.to_string()
}

/// Cross-platform "open this URL or path in the user's preferred
/// app". The host wires the `LandingAction::OpenLink` variant to this
/// — it routes URLs through the default web browser and filesystem
/// paths through the host's file browser. Returns `Err(reason)` only
/// when the launcher fails to spawn.
///
/// Uses the `opener` crate which dispatches to:
///   - `ShellExecuteW` (Win32) on Windows — passes the URL/path as a
///     wide-char native arg, NOT through cmd.exe re-parsing. The old
///     `cmd /C start "" <target>` shell-out was vulnerable to cmd's
///     metachar handling (`&`, `^`, `|`, `>`, `<`, `%var%`) so paths
///     like `C:\Apps & Tools\Valenx\` would break or worse.
///   - `open(1)` on macOS.
///   - The xdg-open preference order on Linux / BSDs.
pub fn open_link(target: &str) -> Result<(), String> {
    opener::open(target).map_err(|e| format!("open link launcher: {e}"))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn local_doc_or_url_returns_url_when_file_missing() {
        // A clearly-impossible relative path falls through to the URL.
        let out = local_doc_or_url("__definitely_not_a_real_file__.md", "https://example/");
        assert_eq!(out, "https://example/");
    }

    #[test]
    fn default_links_resolve_returns_three_well_formed_entries() {
        // The host wires three quick-link rows (Learn more / What's new
        // / Documentation). They should resolve into three labelled
        // entries with non-empty targets — when no local README /
        // CHANGELOG / docs/INSTALLER.md is present, the targets fall
        // through to the repo URL (which is what an installed build
        // sees). When the workspace files ARE present (a dev / source
        // checkout), the targets resolve to local paths.
        let repo = "https://github.com/nochallenge/valenx";
        let links = LandingAction::default_links_resolve(repo);
        assert_eq!(links.len(), 3);
        // Labels are stable enough to assert against — they're the
        // strings shown in the UI.
        assert_eq!(links[0].0, "Learn more");
        assert_eq!(links[1].0, "What's new");
        assert_eq!(links[2].0, "Documentation");
        for (label, target) in &links {
            assert!(!target.is_empty(), "{label}: target is empty");
            // Either a URL (starts with https://) OR a filesystem path —
            // both are valid LandingAction::OpenLink payloads that
            // `open_link` accepts. We just check we're not handing back
            // garbage like a stray template placeholder.
            assert!(
                target.starts_with("https://") || std::path::Path::new(target).is_absolute(),
                "{label}: target should be URL or absolute path, got {target}"
            );
        }
    }

    #[test]
    fn default_links_use_supplied_repo_url() {
        // When the local docs are absent the fallback URLs must include
        // the host-supplied repo URL — proves the function isn't
        // silently hardcoding a different repo.
        // The Documentation link points at `<repo>/blob/master/docs/INSTALLER.md`
        // which always resolves to a URL during this test because we
        // change into an unrelated tempdir first (so the relative
        // "docs/INSTALLER.md" probe in `local_doc_or_url` fails).
        let tmp = tempfile::tempdir().expect("create tempdir");
        let prev_cwd = std::env::current_dir().expect("read cwd");
        std::env::set_current_dir(tmp.path()).expect("chdir to tempdir");
        let result = std::panic::catch_unwind(|| {
            let repo = "https://example.test/owner/repo";
            let links = LandingAction::default_links_resolve(repo);
            for (label, target) in &links {
                assert!(
                    target.starts_with(repo),
                    "{label}: expected URL prefixed with {repo}, got {target}"
                );
            }
        });
        // Always restore cwd before propagating any test failure so
        // the rest of the suite sees the right working directory.
        std::env::set_current_dir(prev_cwd).expect("restore cwd");
        result.expect("inner assertions");
    }

    #[test]
    fn landing_action_open_recent_dispatch_construction() {
        // Sanity-check that the host can match against the action
        // variants for OpenRecent + DropMissingRecent without surprise
        // — these are the two outcomes of clicking a recent-project
        // row (present vs missing on disk).
        let path = std::path::PathBuf::from("/some/path/to/project.valenx");
        let action = LandingAction::OpenRecent(path.clone());
        match action {
            LandingAction::OpenRecent(p) => assert_eq!(p, path),
            other => panic!("expected OpenRecent, got {other:?}"),
        }
        let action = LandingAction::DropMissingRecent(path.clone());
        match action {
            LandingAction::DropMissingRecent(p) => assert_eq!(p, path),
            other => panic!("expected DropMissingRecent, got {other:?}"),
        }
    }
}
