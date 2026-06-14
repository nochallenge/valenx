//! The **Assistant** activity sidebar — valenx narrates its own work.
//!
//! A right-side panel that displays a *live* feed of what the AI assistant
//! is doing, so the desktop app is self-describing: anyone watching can
//! see, in-app, what is being designed / simulated / built right now.
//!
//! The feed is a newline-delimited JSON (`.jsonl`) file — one
//! [`FeedEntry`] per line — that an external agent (or the app itself)
//! appends to. The panel re-reads it about once a second (via
//! [`egui::Context::request_repaint_after`]) so new entries appear live
//! without any user interaction. A missing or unreadable feed file is
//! treated as "no activity yet", never an error, and a half-written final
//! line (mid-append) is skipped rather than breaking the whole feed.
//!
//! The feed path is `$VALENX_ASSISTANT_FEED` when set, otherwise
//! `<state_dir>/assistant_feed.jsonl`.

use std::path::{Path, PathBuf};
use std::time::Duration;

use eframe::egui;
use serde::Deserialize;

use crate::ValenxApp;

/// One entry in the assistant activity feed.
#[derive(Clone, Debug, Default, Deserialize)]
pub struct FeedEntry {
    /// Short timestamp label (e.g. `"20:03"`). Free-form; shown verbatim.
    #[serde(default)]
    pub time: String,
    /// One-line headline of the activity.
    #[serde(default)]
    pub title: String,
    /// Optional longer detail shown under the title.
    #[serde(default)]
    pub detail: String,
    /// Category tag driving the accent colour: `build` / `result` / `ship`
    /// / `warn` (anything else → neutral).
    #[serde(default)]
    pub kind: String,
}

/// Persistent state for the Assistant panel — just the resolved feed path
/// (entries are re-read from disk each frame so the feed stays live).
pub struct AssistantWorkbenchState {
    feed_path: PathBuf,
}

impl Default for AssistantWorkbenchState {
    fn default() -> Self {
        Self {
            feed_path: assistant_feed_path(),
        }
    }
}

/// Resolve the assistant feed file: `$VALENX_ASSISTANT_FEED` if set,
/// otherwise `<state_dir>/assistant_feed.jsonl` (falling back to the
/// system temp dir when no per-user state dir resolves).
pub fn assistant_feed_path() -> PathBuf {
    if let Ok(p) = std::env::var("VALENX_ASSISTANT_FEED") {
        if !p.is_empty() {
            return PathBuf::from(p);
        }
    }
    crate::state_paths::state_dir()
        .map(|d| d.join("assistant_feed.jsonl"))
        .unwrap_or_else(|| std::env::temp_dir().join("valenx_assistant_feed.jsonl"))
}

/// Parse a `.jsonl` feed body into entries, skipping blank or malformed
/// lines (a half-written final line while the file is being appended must
/// not break the whole feed).
fn parse_feed(body: &str) -> Vec<FeedEntry> {
    body.lines()
        .map(str::trim)
        .filter(|l| !l.is_empty())
        .filter_map(|l| serde_json::from_str::<FeedEntry>(l).ok())
        .collect()
}

/// Read + parse the feed file. A missing / unreadable file yields an empty
/// feed (not an error) — the panel then shows "waiting".
fn load_feed(path: &Path) -> Vec<FeedEntry> {
    std::fs::read_to_string(path)
        .map(|s| parse_feed(&s))
        .unwrap_or_default()
}

/// Accent colour for an entry's `kind` tag.
fn accent(kind: &str) -> egui::Color32 {
    match kind {
        "build" => egui::Color32::from_rgb(55, 138, 221), // blue
        "result" => egui::Color32::from_rgb(80, 200, 140), // teal/green
        "ship" => egui::Color32::from_rgb(127, 119, 221), // purple
        "warn" => egui::Color32::from_rgb(220, 160, 60),  // amber
        _ => egui::Color32::from_rgb(150, 150, 145),      // neutral
    }
}

/// Draw the Assistant activity sidebar. A no-op when the
/// `show_assistant_panel` toggle is off.
pub fn draw_assistant_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_assistant_panel {
        return;
    }
    // Live: poll the feed about once a second even when otherwise idle, so
    // appended entries appear without the user touching anything.
    ctx.request_repaint_after(Duration::from_millis(1000));
    let entries = load_feed(&app.assistant.feed_path);

    egui::SidePanel::right("valenx_assistant_panel")
        .resizable(true)
        .default_width(330.0)
        .width_range(260.0..=560.0)
        .show(ctx, |ui| {
            ui.horizontal(|ui| {
                ui.heading("Assistant");
                ui.label(
                    egui::RichText::new("● live")
                        .small()
                        .color(egui::Color32::from_rgb(80, 200, 140)),
                );
            });
            ui.label(
                egui::RichText::new("what Claude is building, live in-app")
                    .weak()
                    .small(),
            );
            ui.separator();
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .stick_to_bottom(true)
                .show(ui, |ui| {
                    if entries.is_empty() {
                        ui.add_space(8.0);
                        ui.label(
                            egui::RichText::new("Waiting for activity…")
                                .weak()
                                .italics(),
                        );
                        return;
                    }
                    for e in &entries {
                        ui.add_space(8.0);
                        ui.horizontal(|ui| {
                            if !e.time.is_empty() {
                                ui.label(egui::RichText::new(&e.time).monospace().small().weak());
                            }
                            ui.label(
                                egui::RichText::new(&e.title)
                                    .strong()
                                    .color(accent(&e.kind)),
                            );
                        });
                        if !e.detail.is_empty() {
                            ui.label(egui::RichText::new(&e.detail).small());
                        }
                        ui.separator();
                    }
                });
        });
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_feed_reads_entries_and_skips_garbage() {
        let body = "\
{\"time\":\"20:01\",\"title\":\"Designed the LV-1\",\"detail\":\"two-stage kerolox\",\"kind\":\"build\"}

not json at all
{\"title\":\"Reached orbit\",\"kind\":\"result\"}
{\"time\":\"20:05\",\"title\":\"half-written";
        let entries = parse_feed(body);
        // Two well-formed lines parse; the blank, the garbage, and the
        // truncated final line are skipped.
        assert_eq!(entries.len(), 2);
        assert_eq!(entries[0].title, "Designed the LV-1");
        assert_eq!(entries[0].kind, "build");
        assert_eq!(entries[0].detail, "two-stage kerolox");
        // Missing fields default to empty, not an error.
        assert_eq!(entries[1].title, "Reached orbit");
        assert_eq!(entries[1].time, "");
        assert_eq!(entries[1].detail, "");
    }

    #[test]
    fn load_feed_missing_file_is_empty_not_an_error() {
        let path = std::env::temp_dir().join("valenx_assistant_feed_does_not_exist_xyz.jsonl");
        let _ = std::fs::remove_file(&path);
        assert!(load_feed(&path).is_empty());
    }

    #[test]
    fn accent_maps_known_kinds() {
        assert_ne!(accent("build"), accent("result"));
        assert_ne!(accent("ship"), accent("warn"));
        // Unknown kind falls back to the neutral colour.
        assert_eq!(accent("whatever"), accent(""));
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;
    use std::io::Write;

    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_assistant_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        app.show_assistant_panel = false;
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_waiting_state_without_panic() {
        // No feed file → the panel renders the "waiting" placeholder.
        let mut app = ValenxApp::default();
        app.show_assistant_panel = true;
        app.assistant.feed_path =
            std::env::temp_dir().join("valenx_assistant_feed_absent_for_test.jsonl");
        let _ = std::fs::remove_file(&app.assistant.feed_path);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_a_populated_feed_without_panic() {
        let mut f = tempfile::NamedTempFile::new().unwrap();
        let line = r#"{"time":"20:03","title":"Designed Valenx LV-1","detail":"reached 146 x 8267 km orbit","kind":"result"}"#;
        writeln!(f, "{line}").unwrap();
        let mut app = ValenxApp::default();
        app.show_assistant_panel = true;
        app.assistant.feed_path = f.path().to_path_buf();
        draw_workbench(&mut app);
    }
}
