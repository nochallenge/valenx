//! The **Assistant** chat sidebar — a two-way channel with Claude.
//!
//! A right-side panel that shows a *live* feed of what the AI assistant is
//! doing AND lets the user type messages back. The feed is a
//! newline-delimited JSON (`.jsonl`) file — one [`FeedEntry`] per line — that
//! an external agent (or the app itself) appends to; the panel re-reads it
//! about once a second so new entries appear live. Messages the user types are
//! appended to an *inbox* `.jsonl` (for an external agent to read) and echoed
//! into the feed so they appear in the panel immediately.
//!
//! Feed path:  `$VALENX_ASSISTANT_FEED`  else `<state_dir>/assistant_feed.jsonl`.
//! Inbox path: `$VALENX_ASSISTANT_INBOX` else `<state_dir>/assistant_inbox.jsonl`.

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
    /// Category tag driving the accent colour: `build` / `result` / `ship` /
    /// `warn` / `user` (anything else → neutral).
    #[serde(default)]
    pub kind: String,
}

/// Persistent state for the Assistant panel.
pub struct AssistantWorkbenchState {
    feed_path: PathBuf,
    inbox_path: PathBuf,
    /// The in-progress message in the chat input box.
    input: String,
}

impl Default for AssistantWorkbenchState {
    fn default() -> Self {
        Self {
            feed_path: assistant_feed_path(),
            inbox_path: assistant_inbox_path(),
            input: String::new(),
        }
    }
}

/// Resolve the assistant feed file: `$VALENX_ASSISTANT_FEED` if set, otherwise
/// `<state_dir>/assistant_feed.jsonl` (system-temp fallback).
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

/// Resolve the assistant **inbox** file (messages the user types in-app, for
/// an external agent to read): `$VALENX_ASSISTANT_INBOX` if set, otherwise
/// `<state_dir>/assistant_inbox.jsonl` (system-temp fallback).
pub fn assistant_inbox_path() -> PathBuf {
    if let Ok(p) = std::env::var("VALENX_ASSISTANT_INBOX") {
        if !p.is_empty() {
            return PathBuf::from(p);
        }
    }
    crate::state_paths::state_dir()
        .map(|d| d.join("assistant_inbox.jsonl"))
        .unwrap_or_else(|| std::env::temp_dir().join("valenx_assistant_inbox.jsonl"))
}

/// Append a single line to a `.jsonl` file (create + append; best-effort).
fn append_line(path: &Path, line: &str) {
    use std::io::Write;
    if let Ok(mut f) = std::fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(path)
    {
        let _ = writeln!(f, "{line}");
    }
}

/// Send a user message: append it to the inbox (for the external agent) and
/// echo it into the feed so it appears in the panel immediately.
fn send_to_assistant(inbox: &Path, feed: &Path, msg: &str) {
    append_line(inbox, &serde_json::json!({ "text": msg }).to_string());
    append_line(
        feed,
        &serde_json::json!({ "title": "You", "detail": msg, "kind": "user" }).to_string(),
    );
}

/// Parse a `.jsonl` feed body into entries, skipping blank or malformed lines
/// (a half-written final line while the file is being appended must not break
/// the whole feed).
fn parse_feed(body: &str) -> Vec<FeedEntry> {
    body.lines()
        .map(str::trim)
        .filter(|l| !l.is_empty())
        .filter_map(|l| serde_json::from_str::<FeedEntry>(l).ok())
        .collect()
}

/// Read + parse the feed file. A missing / unreadable file yields an empty
/// feed (not an error).
fn load_feed(path: &Path) -> Vec<FeedEntry> {
    std::fs::read_to_string(path)
        .map(|s| parse_feed(&s))
        .unwrap_or_default()
}

/// The most recent **build / result / ship** headline in the assistant feed,
/// as a short `"title — detail"` (or just `title`) string — `None` if the feed
/// has no such entry yet. Used by the dock's `"workspace:<n>"` placeholder
/// tile to surface the latest thing the agent reported building, without that
/// tile needing to know the feed's on-disk format. Skips the user's own
/// `"user"` echoes so it reflects the agent's progress, not the prompt.
pub(crate) fn latest_build_status(app: &ValenxApp) -> Option<String> {
    load_feed(&app.assistant.feed_path)
        .into_iter()
        .rev()
        .find(|e| matches!(e.kind.as_str(), "build" | "result" | "ship") && !e.title.is_empty())
        .map(|e| {
            if e.detail.is_empty() {
                e.title
            } else {
                format!("{} — {}", e.title, e.detail)
            }
        })
}

/// Accent colour for an entry's `kind` tag.
fn accent(kind: &str) -> egui::Color32 {
    match kind {
        "build" => egui::Color32::from_rgb(55, 138, 221), // blue
        "result" => egui::Color32::from_rgb(80, 200, 140), // teal/green
        "ship" => egui::Color32::from_rgb(127, 119, 221), // purple
        "warn" => egui::Color32::from_rgb(220, 160, 60),  // amber
        "user" => egui::Color32::from_rgb(120, 170, 255), // light blue (your messages)
        _ => egui::Color32::from_rgb(150, 150, 145),      // neutral
    }
}

/// Draw the Assistant chat sidebar. A no-op when the `show_assistant_panel`
/// toggle is off.
pub fn draw_assistant_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_assistant_panel {
        return;
    }
    // Live: poll the feed about once a second even when otherwise idle, so
    // appended entries (incl. replies) appear without the user touching
    // anything.
    ctx.request_repaint_after(Duration::from_millis(1000));

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_assistant_panel",
        "Assistant",
        assistant_workbench_body,
    );
    if close {
        app.show_assistant_panel = false;
    }
}

/// The Assistant activity-sidebar body — the bottom-pinned chat input plus
/// the scrolling live feed. Extracted from [`draw_assistant_workbench`] so
/// it can be hosted by the classic
/// [`crate::workbench_chrome::workbench_shell`] *or* the opt-in dockable
/// tile layout ([`crate::dock_layout`]). Reloads the feed + reschedules the
/// ~1 s idle repaint up front so the dock path stays live.
pub(crate) fn assistant_workbench_body(app: &mut ValenxApp, ui: &mut egui::Ui) {
    ui.ctx().request_repaint_after(Duration::from_millis(1000));
    let entries = load_feed(&app.assistant.feed_path);
    ui.label(
        egui::RichText::new("● live  ·  chat with Claude — type below")
            .weak()
            .small()
            .color(egui::Color32::from_rgb(80, 200, 140)),
    );
    ui.separator();
    // Chat input pinned to the bottom; the feed scrolls above it. The panel
    // id is scoped to the host `ui` so multiple Assistant bodies (e.g. six
    // "Agent N" dock tiles) each get a unique TopBottomPanel id instead of
    // colliding on one shared string.
    egui::TopBottomPanel::bottom(ui.id().with("valenx_assistant_input")).show_inside(ui, |ui| {
        ui.add_space(4.0);
        ui.horizontal(|ui| {
            let resp = ui.add(
                egui::TextEdit::singleline(&mut app.assistant.input)
                    .hint_text("Message Claude…")
                    .desired_width(f32::INFINITY),
            );
            let enter = resp.lost_focus() && ui.input(|i| i.key_pressed(egui::Key::Enter));
            if enter || ui.button("Send").clicked() {
                let msg = app.assistant.input.trim().to_string();
                if !msg.is_empty() {
                    send_to_assistant(&app.assistant.inbox_path, &app.assistant.feed_path, &msg);
                    app.assistant.input.clear();
                    ui.ctx().request_repaint();
                    resp.request_focus();
                }
            }
        });
        ui.add_space(2.0);
    });
    egui::ScrollArea::vertical()
        .auto_shrink([false, false])
        .stick_to_bottom(true)
        .show(ui, |ui| {
            if entries.is_empty() {
                ui.add_space(8.0);
                ui.label(
                    egui::RichText::new("No messages yet — say hi below.")
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
            // Interactive cue: if the newest message is the user's, show
            // a "responding" indicator until a reply lands in the feed.
            if entries.last().map(|e| e.kind.as_str()) == Some("user") {
                ui.add_space(6.0);
                ui.label(
                    egui::RichText::new("Claude is responding...")
                        .weak()
                        .italics(),
                );
            }
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
        assert_eq!(entries.len(), 2);
        assert_eq!(entries[0].title, "Designed the LV-1");
        assert_eq!(entries[0].kind, "build");
        assert_eq!(entries[0].detail, "two-stage kerolox");
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
        assert_ne!(accent("user"), accent("build"));
        assert_eq!(accent("whatever"), accent(""));
    }

    #[test]
    fn send_to_assistant_writes_inbox_and_echoes_feed() {
        let dir = std::env::temp_dir();
        let inbox = dir.join("valenx_test_assistant_inbox_xyz.jsonl");
        let feed = dir.join("valenx_test_assistant_feed_xyz.jsonl");
        let _ = std::fs::remove_file(&inbox);
        let _ = std::fs::remove_file(&feed);
        send_to_assistant(&inbox, &feed, "hello claude");
        let ib = std::fs::read_to_string(&inbox).unwrap();
        let fd = std::fs::read_to_string(&feed).unwrap();
        assert!(ib.contains("hello claude"));
        // The feed echo is a parseable FeedEntry with the "You" title.
        let echoed = parse_feed(&fd);
        assert_eq!(echoed.len(), 1);
        assert_eq!(echoed[0].title, "You");
        assert_eq!(echoed[0].detail, "hello claude");
        assert_eq!(echoed[0].kind, "user");
        let _ = std::fs::remove_file(&inbox);
        let _ = std::fs::remove_file(&feed);
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
    fn workbench_draws_empty_state_without_panic() {
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
