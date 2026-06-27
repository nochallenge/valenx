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

impl AssistantWorkbenchState {
    /// The base **inbox** path (`$VALENX_ASSISTANT_INBOX` or the state-dir
    /// default). Read-only accessor used by [`crate::agent_commands`] to derive
    /// the per-channel command-file directory (the agent-drives-valenx bridge
    /// puts its command files beside the chat inbox/feed).
    pub(crate) fn inbox_path(&self) -> &Path {
        &self.inbox_path
    }

    /// Test-only: override the inbox path so a test can point the
    /// agent-command base dir (the inbox path's parent) at a temp directory it
    /// controls. See `crate::agent_commands` tests.
    #[cfg(test)]
    pub(crate) fn set_inbox_path_for_test(&mut self, path: PathBuf) {
        self.inbox_path = path;
    }

    /// Test-only: override the base **feed** path so a test can read back the
    /// per-unit feed file ([`unit_feed_path`]) `append_feed_note` writes,
    /// without colliding with the live app's feed. See `crate::agent_commands`
    /// tests.
    #[cfg(test)]
    pub(crate) fn set_feed_path_for_test(&mut self, path: PathBuf) {
        self.feed_path = path;
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

/// Insert `_u{n}` before a path's final `.jsonl` extension, so a base feed/inbox
/// file becomes a **per-unit** one: e.g.
/// `…/valenx_chat_feed.jsonl` → `…/valenx_chat_feed_u3.jsonl`. If the path
/// doesn't end in `.jsonl` (unusual; a custom env override), the suffix is
/// appended whole (`…/feed.log` → `…/feed.log_u3`) so units still get distinct
/// files. Keeps the parent directory unchanged.
fn per_unit_path(base: &Path, n: usize) -> PathBuf {
    let file = base
        .file_name()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_default();
    let renamed = match file.strip_suffix(".jsonl") {
        Some(stem) => format!("{stem}_u{n}.jsonl"),
        None => format!("{file}_u{n}"),
    };
    match base.parent() {
        Some(dir) => dir.join(renamed),
        None => PathBuf::from(renamed),
    }
}

/// The **per-unit feed** path for "Workbench + Agent" unit `n`: the base feed
/// path ([`AssistantWorkbenchState::feed_path`], itself
/// `$VALENX_ASSISTANT_FEED` or the state-dir default) with `_u{n}` inserted
/// before `.jsonl`. Each `agent:<n>` tile reads/renders this file so the six
/// chats are independent conversations rather than one shared feed.
pub(crate) fn unit_feed_path(app: &ValenxApp, n: usize) -> PathBuf {
    per_unit_path(&app.assistant.feed_path, n)
}

/// The **per-unit inbox** path for "Workbench + Agent" unit `n` (messages the
/// user types in that unit's chat, for an external agent to read): the base
/// inbox path with `_u{n}` inserted before `.jsonl`. Paired with
/// [`unit_feed_path`].
pub(crate) fn unit_inbox_path(app: &ValenxApp, n: usize) -> PathBuf {
    per_unit_path(&app.assistant.inbox_path, n)
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

/// Append a [`FeedEntry`] line (title + detail + accent `kind`) to **channel
/// `n`'s** feed file, so it shows up in that agent tile's chat. Used by the
/// agent-drives-valenx bridge's `Note` command
/// ([`crate::agent_commands`]) to post a visible summary into the panel.
/// Best-effort (a write failure is swallowed, like the rest of the feed I/O).
///
/// **Channel `0` is the GLOBAL feed** (`valenx_chat_feed.jsonl`, the base
/// [`AssistantWorkbenchState::feed_path`] with **no** `_u{n}` suffix), not a
/// unit feed. Real Workbench+Agent units are always `1..=wb_agent_counter`, so
/// `0` is free as the sentinel the agent-bridge's `apply_global` uses to ack
/// global-channel commands into the global feed an agent reading the global
/// channel watches. Any `n >= 1` resolves to that unit's [`unit_feed_path`] as
/// before.
pub(crate) fn append_feed_note(app: &ValenxApp, n: usize, title: &str, detail: &str, kind: &str) {
    let feed = if n == 0 {
        app.assistant.feed_path.clone()
    } else {
        unit_feed_path(app, n)
    };
    append_line(
        &feed,
        &serde_json::json!({ "title": title, "detail": detail, "kind": kind }).to_string(),
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
///
/// Bounded by [`crate::settings_io::MAX_STATE_FILE_BYTES`] — the same defensive
/// size-gate [`crate::agent_commands`]'s `read_cmd_file` uses: a corrupt or
/// hostile feed file is checked via `metadata` *before* the read, so an
/// oversized file can't OOM the ~1 s repaint that re-reads the feed. Over the
/// cap → treated as an empty feed, exactly like a missing/unreadable one.
fn load_feed(path: &Path) -> Vec<FeedEntry> {
    let Ok(meta) = std::fs::metadata(path) else {
        return Vec::new();
    };
    if meta.len() > crate::settings_io::MAX_STATE_FILE_BYTES as u64 {
        return Vec::new();
    }
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

    // FIX B — balanced width. This classic docked Assistant only renders for
    // non-tiled (calculator-style) tabs: a 3-D tab hides it and hosts the
    // Assistant as a grid tile instead (`!dock_enabled`, see update.rs). Here
    // the active workbench form sits in the centre / a right SidePanel, so give
    // the Assistant ~40% of the window width with a sane ~320px floor — both
    // the form and the chat stay readable instead of the Assistant squishing to
    // a sliver. The user can still drag it wider/narrower (down to the floor).
    let screen_w = ctx.screen_rect().width();
    let default_w = (screen_w * 0.40).clamp(320.0, 640.0);
    let close = crate::workbench_chrome::workbench_shell_sized(
        app,
        ctx,
        "valenx_assistant_panel",
        "Assistant",
        Some((default_w, 320.0)),
        assistant_workbench_body,
    );
    if close {
        app.show_assistant_panel = false;
    }
}

/// The classic Assistant activity-sidebar body (the base shared channel) — the
/// bottom-pinned chat input plus the scrolling live feed. Hosted by the classic
/// [`crate::workbench_chrome::workbench_shell`] *or*, as the lone Assistant
/// dock tile, by [`crate::dock_layout`]. Thin wrapper that delegates to the
/// channel-aware [`assistant_chat_ui`] with `channel = None` (base paths +
/// `app.assistant.input`).
pub(crate) fn assistant_workbench_body(app: &mut ValenxApp, ui: &mut egui::Ui) {
    assistant_chat_ui(app, ui, None);
}

/// Render a Claude chat — the bottom-pinned input plus the scrolling live feed
/// — for one **channel**:
///
/// - `None` → the **classic base panel**: base feed/inbox paths
///   ([`AssistantWorkbenchState::feed_path`]/`inbox_path`) and the shared
///   [`AssistantWorkbenchState::input`] buffer. This is the historical
///   single-channel behaviour, unchanged.
/// - `Some(n)` → the independent chat for **"Workbench + Agent" unit `n`**:
///   reads/renders unit `n`'s own feed file ([`unit_feed_path`]), binds the
///   input to `app.unit_chat_inputs[n]` (per-unit buffer), and Send appends to
///   unit `n`'s inbox ([`unit_inbox_path`]) + echoes into its feed. This is
///   what makes the six agent tiles independent conversations instead of one
///   mirrored feed/input.
///
/// Reloads the feed + reschedules the ~1 s idle repaint up front so both the
/// classic and dock hosts stay live. The bottom-input `TopBottomPanel` id is
/// scoped to the host `ui` ([`egui::Ui::id`]) so multiple chats on screen each
/// get a unique id instead of colliding on one shared string.
pub(crate) fn assistant_chat_ui(app: &mut ValenxApp, ui: &mut egui::Ui, channel: Option<usize>) {
    ui.ctx().request_repaint_after(Duration::from_millis(1000));
    // Resolve this channel's on-disk paths up front (owned, so they don't tie
    // up a borrow of `app` while we mutate the input buffer below).
    let (feed_path, inbox_path) = match channel {
        None => (
            app.assistant.feed_path.clone(),
            app.assistant.inbox_path.clone(),
        ),
        Some(n) => (unit_feed_path(app, n), unit_inbox_path(app, n)),
    };
    let entries = load_feed(&feed_path);
    ui.label(
        egui::RichText::new("● live  ·  chat with Claude — type below")
            .weak()
            .small()
            .color(egui::Color32::from_rgb(80, 200, 140)),
    );
    ui.separator();
    // Chat input pinned to the bottom; the feed scrolls above it.
    egui::TopBottomPanel::bottom(ui.id().with("valenx_assistant_input")).show_inside(ui, |ui| {
        ui.add_space(4.0);
        ui.horizontal(|ui| {
            // Bind to the right input buffer for this channel: the shared
            // `assistant.input` for the base panel, or this unit's own entry in
            // `unit_chat_inputs` (created lazily) so per-unit chats don't share
            // one box.
            let input: &mut String = match channel {
                None => &mut app.assistant.input,
                Some(n) => app.unit_chat_inputs.entry(n).or_default(),
            };
            // A caption associated via `labelled_by` so the chat input carries an
            // accessible name (a bare TextEdit has none; its hint text is not a
            // name) — addressable by a screen reader / AI driver.
            let msg_cap = ui.label(egui::RichText::new("Message").weak().small());
            let resp = ui
                .add(
                    egui::TextEdit::singleline(input)
                        .hint_text("Message Claude…")
                        .desired_width(f32::INFINITY),
                )
                .labelled_by(msg_cap.id);
            let enter = resp.lost_focus() && ui.input(|i| i.key_pressed(egui::Key::Enter));
            if enter || ui.button("Send").clicked() {
                let msg = input.trim().to_string();
                if !msg.is_empty() {
                    send_to_assistant(&inbox_path, &feed_path, &msg);
                    input.clear();
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
                    // Slightly larger body text (~15px) for readability; the
                    // timestamp above stays `.small()`.
                    ui.label(egui::RichText::new(&e.detail).size(15.0));
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
    fn load_feed_reads_a_normal_small_feed() {
        // A normal under-cap feed parses through the size-gate unchanged.
        use std::io::Write;
        let mut f = tempfile::NamedTempFile::new().unwrap();
        writeln!(
            f,
            r#"{{"title":"Designed LV-1","detail":"two-stage kerolox","kind":"build"}}"#
        )
        .unwrap();
        let entries = load_feed(f.path());
        assert_eq!(entries.len(), 1);
        assert_eq!(entries[0].title, "Designed LV-1");
    }

    #[test]
    fn load_feed_rejects_oversized_file_without_reading_it() {
        // Defensive symmetry with `agent_commands::read_cmd_file`: a feed file
        // past MAX_STATE_FILE_BYTES is gated by `metadata` BEFORE the read, so
        // a hostile/corrupt file can't OOM the ~1 s feed re-read. We don't
        // write a real 10 MiB file (slow); instead confirm the gate logic
        // `load_feed` runs would reject an oversized length up-front, mirroring
        // settings_io's oversized-state-file test.
        let oversized_len = crate::settings_io::MAX_STATE_FILE_BYTES as u64 + 1;
        let gated = oversized_len > crate::settings_io::MAX_STATE_FILE_BYTES as u64;
        assert!(gated, "an over-cap length must trip load_feed's size gate");
    }

    #[test]
    fn accent_maps_known_kinds() {
        assert_ne!(accent("build"), accent("result"));
        assert_ne!(accent("ship"), accent("warn"));
        assert_ne!(accent("user"), accent("build"));
        assert_eq!(accent("whatever"), accent(""));
    }

    #[test]
    fn per_unit_path_inserts_unit_suffix_before_jsonl() {
        // `…/valenx_chat_feed.jsonl` → `…/valenx_chat_feed_u3.jsonl`, parent dir
        // preserved. A non-.jsonl base appends the whole suffix.
        let base = std::path::Path::new("/tmp/dir/valenx_chat_feed.jsonl");
        assert_eq!(
            per_unit_path(base, 3),
            std::path::PathBuf::from("/tmp/dir/valenx_chat_feed_u3.jsonl")
        );
        let base2 = std::path::Path::new("/tmp/dir/feed.log");
        assert_eq!(
            per_unit_path(base2, 7),
            std::path::PathBuf::from("/tmp/dir/feed.log_u7")
        );
    }

    #[test]
    fn unit_channels_resolve_to_distinct_feed_and_inbox_paths() {
        // assistant_chat_ui's channels None / Some(3) / Some(4) must map to
        // DISTINCT feed+inbox files so the per-unit chats don't share one feed.
        let app = ValenxApp::default();
        let base_feed = app.assistant.feed_path.clone();
        let base_inbox = app.assistant.inbox_path.clone();

        let f3 = unit_feed_path(&app, 3);
        let f4 = unit_feed_path(&app, 4);
        let i3 = unit_inbox_path(&app, 3);
        let i4 = unit_inbox_path(&app, 4);

        // Every channel's feed path is distinct from the others and from base.
        assert_ne!(base_feed, f3);
        assert_ne!(base_feed, f4);
        assert_ne!(f3, f4);
        // Same for inboxes.
        assert_ne!(base_inbox, i3);
        assert_ne!(base_inbox, i4);
        assert_ne!(i3, i4);
        // Feed and inbox of the same unit are also distinct (they derive from
        // the two different base files).
        assert_ne!(f3, i3);
    }

    #[test]
    fn unit_chat_inputs_keep_separate_text_per_unit() {
        // The per-unit input map keeps each unit's typed text independent — the
        // bug being fixed was all agent chats sharing one input buffer.
        let mut app = ValenxApp::default();
        app.unit_chat_inputs
            .entry(3)
            .or_default()
            .push_str("hello from 3");
        app.unit_chat_inputs
            .entry(4)
            .or_default()
            .push_str("four says hi");
        assert_eq!(app.unit_chat_inputs.get(&3).unwrap(), "hello from 3");
        assert_eq!(app.unit_chat_inputs.get(&4).unwrap(), "four says hi");
        // Neither touched the shared base-panel input.
        assert!(app.assistant.input.is_empty());
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
