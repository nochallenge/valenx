//! **Agent-drives-valenx bridge (v1).**
//!
//! Lets an *external* agent (the Claude session backing a "Workbench + Agent"
//! chat tile) actually **drive** valenx: open / focus / rename / close project
//! tabs, switch a tab's workbench, and post a visible summary into the chat —
//! all by appending newline-delimited JSON *commands* to a per-channel command
//! file that valenx polls each frame and executes through the **existing,
//! vetted** tab / dock methods. No raw app field is poked; every effect goes
//! through the same code path a user click would.
//!
//! ## Command file
//!
//! Each agent channel `n` (an `agent:<n>` dock tile) has its own command file,
//! a sibling of that channel's chat feed/inbox:
//! `<base-dir>/valenx_chat_cmd_u{n}.jsonl` (see [`cmd_path`]). The base dir is
//! the parent of [`crate::assistant_workbench`]'s inbox path (i.e.
//! `$VALENX_ASSISTANT_INBOX`'s directory, or the state-dir default), so the
//! command channel lives right next to the existing chat channels.
//!
//! ## Wire format
//!
//! One JSON object per line, **internally tagged** on `"cmd"` — e.g.
//! `{"cmd":"new_tab","name":"Rocket","workbench":"rocket"}`. See
//! [`AgentCommand`] for every variant. Unparseable lines are skipped (a
//! half-written final line while the agent is appending must not wedge the
//! channel), exactly like the feed parser.
//!
//! ## Replay safety
//!
//! [`poll_and_apply_agent_commands`] keeps a per-channel **cursor**
//! ([`ValenxApp::agent_cmd_cursor`]) of how many lines it has already applied.
//! On the **first** poll for a channel the cursor is set to the file's
//! *current* line count, so a command file left over from a previous run (or
//! the agent's own scrollback) is **not** re-executed on launch — only
//! genuinely new appended lines run. Disk reads are throttled to ~1/sec via
//! [`ValenxApp::last_agent_poll`].

use std::path::{Path, PathBuf};
use std::time::{Duration, Instant};

use serde::Deserialize;

use crate::project_tabs::{self, TabKind};
use crate::ValenxApp;

/// Throttle: poll the command files at most this often (matches the ~1 s the
/// chat feed itself re-reads at).
const POLL_INTERVAL: Duration = Duration::from_millis(1000);

/// Fixed stem for a per-channel command file (the `_u{n}.jsonl` suffix is
/// appended per channel by [`cmd_path`]).
const CMD_STEM: &str = "valenx_chat_cmd";

/// One command an external agent can append to its channel's command file to
/// drive valenx. **Internally tagged** on `"cmd"`, so each line is a flat
/// object like `{"cmd":"new_tab","name":"…"}`.
///
/// Every variant is honoured through an *existing* vetted tab/dock method (see
/// [`apply`]); none writes a raw [`ValenxApp`] field. Unknown command tags or
/// bad workbench ids are skipped without panicking.
#[derive(Debug, Clone, PartialEq, Eq, Deserialize)]
#[serde(tag = "cmd", rename_all = "snake_case")]
pub enum AgentCommand {
    /// Open a **new project tab** named `name`. If `workbench` is a known
    /// [`TabKind`] id (see [`TabKind::from_id`]) the tab is bound to that
    /// workbench; otherwise it opens blank. Mirrors the "From template" /
    /// "+ New tab" flow and makes the new tab active.
    NewTab {
        /// Title for the new tab.
        name: String,
        /// Optional workbench id to bind (case-insensitive; e.g. `"rocket"`).
        #[serde(default)]
        workbench: Option<String>,
    },
    /// Switch the **active** tab's workbench to `id` (a [`TabKind`] id),
    /// reconciling the visible panel + viewport. No-op (skipped) if no tab is
    /// active or `id` is unknown.
    OpenWorkbench {
        /// Workbench id to switch the active tab to (case-insensitive).
        id: String,
    },
    /// **Focus** (activate) the first tab whose title equals `name`. No-op if
    /// no tab matches.
    FocusTab {
        /// Title of the tab to focus.
        name: String,
    },
    /// **Rename** the active tab to `name`.
    RenameTab {
        /// New title for the active tab.
        name: String,
    },
    /// Post a visible **summary** line into this channel's chat feed (so the
    /// agent's narration shows up in the panel). `kind` is the feed accent tag
    /// (`build` / `result` / `ship` / `warn`); defaults to `ship`.
    Note {
        /// The message body.
        text: String,
        /// Optional feed accent tag.
        #[serde(default)]
        kind: Option<String>,
    },
    /// **Close** a tab — routed through the same "Close tab?" confirm modal a
    /// user ✕ opens (never a silent hard-delete). `name` picks the first tab
    /// with that title; omitted → the active tab.
    CloseTab {
        /// Title of the tab to close, or `None` for the active tab.
        #[serde(default)]
        name: Option<String>,
    },
}

/// The per-channel **command file** path for agent channel `n`:
/// `<base-dir>/valenx_chat_cmd_u{n}.jsonl`, where `<base-dir>` is the directory
/// holding the channel chat files (the parent of
/// [`crate::assistant_workbench`]'s inbox path). Sits beside that channel's
/// feed/inbox so all three channels share one directory.
pub fn cmd_path(app: &ValenxApp, n: usize) -> PathBuf {
    // Derive the base directory from the existing inbox path so the command
    // channel follows the same `$VALENX_ASSISTANT_INBOX` / state-dir override
    // the chat channels use.
    let base_dir: PathBuf = app
        .assistant
        .inbox_path()
        .parent()
        .map(Path::to_path_buf)
        .unwrap_or_else(std::env::temp_dir);
    base_dir.join(format!("{CMD_STEM}_u{n}.jsonl"))
}

/// Count the newline-delimited lines in a `.jsonl` body the same way the
/// applier iterates them — non-empty, trimmed lines — so the cursor and the
/// apply loop always agree on "how many lines are there".
fn line_count(body: &str) -> usize {
    body.lines()
        .map(str::trim)
        .filter(|l| !l.is_empty())
        .count()
}

/// Read a bounded command-file body (`None` if absent / unreadable / too big).
/// Bounded by [`crate::settings_io::MAX_STATE_FILE_BYTES`] so a corrupt or
/// hostile file can't OOM the poll.
fn read_cmd_file(path: &Path) -> Option<String> {
    let meta = std::fs::metadata(path).ok()?;
    if meta.len() > crate::settings_io::MAX_STATE_FILE_BYTES as u64 {
        return None;
    }
    std::fs::read_to_string(path).ok()
}

/// Poll every agent channel's command file and apply any **newly appended**
/// commands, advancing that channel's cursor. Throttled to ~1/sec via
/// [`ValenxApp::last_agent_poll`]; cheap to call every frame.
///
/// Call this from the `update()` loop just before the tab strip is drawn,
/// where a clean `&mut self` is available (no dock-tree borrow held).
pub fn poll_and_apply_agent_commands(app: &mut ValenxApp) {
    // Throttle disk reads.
    let now = Instant::now();
    if let Some(last) = app.last_agent_poll {
        if now.duration_since(last) < POLL_INTERVAL {
            return;
        }
    }
    app.last_agent_poll = Some(now);

    // Scan every channel handed out so far (1..=wb_agent_counter). Channels
    // with no command file are skipped cheaply.
    let highest = app.wb_agent_counter;
    for n in 1..=highest {
        let path = cmd_path(app, n);
        let Some(body) = read_cmd_file(&path) else {
            continue;
        };
        let total = line_count(&body);

        // First poll for this channel: adopt the current line count as the
        // cursor so pre-existing history is NOT replayed on launch.
        let start = match app.agent_cmd_cursor.get(&n) {
            Some(&c) => c.min(total), // clamp in case the file was truncated
            None => {
                app.agent_cmd_cursor.insert(n, total);
                total
            }
        };
        if start >= total {
            continue; // nothing new
        }

        // Apply only the genuinely new lines, skipping unparseable ones (like
        // the feed parser's filter_map(...ok)).
        for line in body
            .lines()
            .map(str::trim)
            .filter(|l| !l.is_empty())
            .skip(start)
        {
            if let Ok(cmd) = serde_json::from_str::<AgentCommand>(line) {
                apply(app, n, cmd);
            }
            // Unparseable line → skip safely.
        }
        // Advance the cursor past everything we just saw.
        app.agent_cmd_cursor.insert(n, total);
    }
}

/// Apply **one** [`AgentCommand`] for channel `n` through existing vetted
/// methods only. Every branch is a no-op-on-bad-input (unknown workbench id,
/// no active tab, no title match) rather than a panic.
fn apply(app: &mut ValenxApp, n: usize, cmd: AgentCommand) {
    match cmd {
        AgentCommand::NewTab { name, workbench } => {
            // Mirror draw_tab_strip's `open_template` intent exactly: park the
            // outgoing tab's scene, open the new tab (pushes a fresh doc + makes
            // it active), then install that empty doc so the new tab starts
            // clean and the prior tab keeps its geometry.
            let kind = workbench
                .as_deref()
                .and_then(TabKind::from_id)
                .unwrap_or(TabKind::Blank);
            project_tabs::park_active_doc(app);
            app.tab_bar.open(kind);
            // Title the freshly-opened (now active) tab.
            if let Some(idx) = app.tab_bar.active {
                if let Some(tab) = app.tab_bar.tabs.get_mut(idx) {
                    let trimmed = name.trim();
                    if !trimmed.is_empty() {
                        tab.title = trimmed.to_string();
                    }
                }
            }
            project_tabs::install_active_doc(app);
        }
        AgentCommand::OpenWorkbench { id } => {
            if let Some(kind) = TabKind::from_id(&id) {
                if let Some(idx) = app.tab_bar.active {
                    if let Some(tab) = app.tab_bar.tabs.get_mut(idx) {
                        tab.kind = kind;
                        project_tabs::sync_active(app);
                    }
                }
            }
        }
        AgentCommand::FocusTab { name } => {
            let target = app.tab_bar.tabs.iter().position(|t| t.title == name);
            if let Some(idx) = target {
                project_tabs::switch_active_to(app, idx);
            }
        }
        AgentCommand::RenameTab { name } => {
            if let Some(idx) = app.tab_bar.active {
                if let Some(tab) = app.tab_bar.tabs.get_mut(idx) {
                    let trimmed = name.trim();
                    if !trimmed.is_empty() {
                        tab.title = trimmed.to_string();
                    }
                }
            }
        }
        AgentCommand::Note { text, kind } => {
            crate::assistant_workbench::append_feed_note(
                app,
                n,
                "Claude",
                &text,
                kind.as_deref().unwrap_or("ship"),
            );
        }
        AgentCommand::CloseTab { name } => {
            // Resolve the target index (named tab, or the active one), then
            // open the confirm modal rather than hard-deleting — the user
            // confirms the destructive close.
            let idx = match name {
                Some(title) => app.tab_bar.tabs.iter().position(|t| t.title == title),
                None => app.tab_bar.active,
            };
            if let Some(idx) = idx {
                if idx < app.tab_bar.tabs.len() {
                    app.tab_close_confirm = Some(idx);
                }
            }
        }
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod tests {
    use super::*;

    #[test]
    fn parses_internally_tagged_commands() {
        // Each variant round-trips from its `{"cmd":...}` wire form.
        let nt: AgentCommand =
            serde_json::from_str(r#"{"cmd":"new_tab","name":"Rocket","workbench":"rocket"}"#)
                .unwrap();
        assert_eq!(
            nt,
            AgentCommand::NewTab {
                name: "Rocket".into(),
                workbench: Some("rocket".into()),
            }
        );
        // workbench is optional.
        let nt2: AgentCommand = serde_json::from_str(r#"{"cmd":"new_tab","name":"boat"}"#).unwrap();
        assert_eq!(
            nt2,
            AgentCommand::NewTab {
                name: "boat".into(),
                workbench: None
            }
        );
        let ow: AgentCommand =
            serde_json::from_str(r#"{"cmd":"open_workbench","id":"fem"}"#).unwrap();
        assert_eq!(ow, AgentCommand::OpenWorkbench { id: "fem".into() });
        let ft: AgentCommand =
            serde_json::from_str(r#"{"cmd":"focus_tab","name":"boat"}"#).unwrap();
        assert_eq!(
            ft,
            AgentCommand::FocusTab {
                name: "boat".into()
            }
        );
        let rt: AgentCommand =
            serde_json::from_str(r#"{"cmd":"rename_tab","name":"hull v2"}"#).unwrap();
        assert_eq!(
            rt,
            AgentCommand::RenameTab {
                name: "hull v2".into()
            }
        );
        let note: AgentCommand =
            serde_json::from_str(r#"{"cmd":"note","text":"reached orbit","kind":"result"}"#)
                .unwrap();
        assert_eq!(
            note,
            AgentCommand::Note {
                text: "reached orbit".into(),
                kind: Some("result".into())
            }
        );
        let ct: AgentCommand = serde_json::from_str(r#"{"cmd":"close_tab"}"#).unwrap();
        assert_eq!(ct, AgentCommand::CloseTab { name: None });
        let ctn: AgentCommand =
            serde_json::from_str(r#"{"cmd":"close_tab","name":"boat"}"#).unwrap();
        assert_eq!(
            ctn,
            AgentCommand::CloseTab {
                name: Some("boat".into())
            }
        );
    }

    #[test]
    fn unknown_command_tag_is_an_error_not_a_panic() {
        // An unknown `cmd` value just fails to parse (the poll loop skips it).
        assert!(serde_json::from_str::<AgentCommand>(r#"{"cmd":"nuke","name":"x"}"#).is_err());
        assert!(serde_json::from_str::<AgentCommand>("not json at all").is_err());
    }

    #[test]
    fn new_tab_opens_a_named_bound_tab() {
        // NewTab{name:"Rocket", workbench:"rocket"} applied to a default app
        // yields an active tab titled "Rocket" of kind Rocket.
        let mut app = ValenxApp::default();
        apply(
            &mut app,
            1,
            AgentCommand::NewTab {
                name: "Rocket".into(),
                workbench: Some("rocket".into()),
            },
        );
        let idx = app.tab_bar.active.expect("a tab is active after NewTab");
        assert_eq!(app.tab_bar.tabs[idx].title, "Rocket");
        assert_eq!(app.tab_bar.tabs[idx].kind, TabKind::Rocket);
        // The bound workbench is the one shown.
        assert!(app.show_rocket_workbench);
    }

    #[test]
    fn new_tab_with_unknown_workbench_opens_blank() {
        let mut app = ValenxApp::default();
        apply(
            &mut app,
            1,
            AgentCommand::NewTab {
                name: "Mystery".into(),
                workbench: Some("does-not-exist".into()),
            },
        );
        let idx = app.tab_bar.active.unwrap();
        assert_eq!(app.tab_bar.tabs[idx].title, "Mystery");
        assert_eq!(app.tab_bar.tabs[idx].kind, TabKind::Blank);
    }

    #[test]
    fn open_workbench_switches_the_active_tab_kind() {
        let mut app = ValenxApp::default();
        app.tab_bar.open(TabKind::Blank);
        project_tabs::sync_active(&mut app);
        apply(
            &mut app,
            1,
            AgentCommand::OpenWorkbench { id: "fem".into() },
        );
        let idx = app.tab_bar.active.unwrap();
        assert_eq!(app.tab_bar.tabs[idx].kind, TabKind::Fem);
        assert!(app.show_fem_workbench);
    }

    #[test]
    fn focus_tab_activates_the_named_tab() {
        let mut app = ValenxApp::default();
        // Two tabs via the reducer so titles are set.
        apply(
            &mut app,
            1,
            AgentCommand::NewTab {
                name: "first".into(),
                workbench: None,
            },
        );
        apply(
            &mut app,
            1,
            AgentCommand::NewTab {
                name: "second".into(),
                workbench: None,
            },
        );
        assert_eq!(app.tab_bar.active, Some(1)); // "second" active
        apply(
            &mut app,
            1,
            AgentCommand::FocusTab {
                name: "first".into(),
            },
        );
        assert_eq!(app.tab_bar.active, Some(0));
    }

    #[test]
    fn rename_tab_changes_the_active_title() {
        let mut app = ValenxApp::default();
        app.tab_bar.open(TabKind::Cad);
        apply(
            &mut app,
            1,
            AgentCommand::RenameTab {
                name: "my part".into(),
            },
        );
        let idx = app.tab_bar.active.unwrap();
        assert_eq!(app.tab_bar.tabs[idx].title, "my part");
    }

    #[test]
    fn close_tab_opens_the_confirm_modal_not_a_hard_delete() {
        let mut app = ValenxApp::default();
        app.tab_bar.open(TabKind::Rocket);
        app.tab_bar.open(TabKind::Cad); // active = 1
        assert_eq!(app.tab_bar.tabs.len(), 2);
        apply(&mut app, 1, AgentCommand::CloseTab { name: None });
        // Tab is NOT removed — the confirm modal is armed at the active index.
        assert_eq!(app.tab_bar.tabs.len(), 2);
        assert_eq!(app.tab_close_confirm, Some(1));
    }

    #[test]
    fn close_tab_by_name_targets_that_tab() {
        let mut app = ValenxApp::default();
        apply(
            &mut app,
            1,
            AgentCommand::NewTab {
                name: "keep".into(),
                workbench: None,
            },
        );
        apply(
            &mut app,
            1,
            AgentCommand::NewTab {
                name: "drop".into(),
                workbench: None,
            },
        );
        apply(
            &mut app,
            1,
            AgentCommand::CloseTab {
                name: Some("keep".into()),
            },
        );
        assert_eq!(app.tab_close_confirm, Some(0)); // "keep" is index 0
    }

    #[test]
    fn cmd_path_is_distinct_per_channel_and_a_jsonl() {
        let app = ValenxApp::default();
        let p1 = cmd_path(&app, 1);
        let p3 = cmd_path(&app, 3);
        assert_ne!(p1, p3);
        assert!(p1.to_string_lossy().ends_with("_u1.jsonl"));
        assert!(p3.to_string_lossy().ends_with("_u3.jsonl"));
        // It lives in the same directory as the chat inbox channel.
        assert_eq!(p1.parent(), app.assistant.inbox_path().parent());
    }

    #[test]
    fn line_count_ignores_blank_lines() {
        assert_eq!(line_count(""), 0);
        assert_eq!(line_count("a\n\nb\n"), 2);
        assert_eq!(line_count("  \n a \n"), 1);
    }

    #[test]
    fn first_poll_adopts_line_count_so_history_is_not_replayed() {
        // A command file that pre-exists at launch must NOT be executed: the
        // first poll sets the cursor to the current line count, runs nothing,
        // and only a later appended line takes effect.
        let mut app = ValenxApp::default();
        app.wb_agent_counter = 1;
        // Point the channel at a temp command file with two pre-existing cmds.
        let dir = std::env::temp_dir();
        let path = dir.join(format!("{CMD_STEM}_u1_replaytest.jsonl"));
        let _ = std::fs::remove_file(&path);
        std::fs::write(
            &path,
            "{\"cmd\":\"new_tab\",\"name\":\"old1\"}\n{\"cmd\":\"new_tab\",\"name\":\"old2\"}\n",
        )
        .unwrap();
        // Drive the cursor logic directly (cmd_path points elsewhere, so we
        // exercise the same read→count→cursor sequence by hand here).
        let body = std::fs::read_to_string(&path).unwrap();
        let total = line_count(&body);
        assert_eq!(total, 2);
        // First poll: no cursor yet → adopt the count, apply nothing.
        assert!(app.agent_cmd_cursor.get(&1).is_none());
        let start = match app.agent_cmd_cursor.get(&1) {
            Some(&c) => c.min(total),
            None => {
                app.agent_cmd_cursor.insert(1, total);
                total
            }
        };
        assert_eq!(start, 2, "first-poll cursor == current line count");
        assert!(start >= total, "nothing new to run on first poll");
        assert_eq!(app.tab_bar.tabs.len(), 0, "pre-existing history not run");
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn poll_runs_only_newly_appended_lines() {
        // End-to-end through poll_and_apply_agent_commands: seed a file, first
        // poll adopts (runs nothing), append one line, second poll runs just
        // that one. Uses the real cmd_path for channel 1.
        let mut app = ValenxApp::default();
        app.wb_agent_counter = 1;
        let path = cmd_path(&app, 1);
        let _ = std::fs::remove_file(&path);
        std::fs::write(&path, "{\"cmd\":\"new_tab\",\"name\":\"preexisting\"}\n").unwrap();

        // First poll: history adopted, nothing run.
        poll_and_apply_agent_commands(&mut app);
        assert_eq!(app.tab_bar.tabs.len(), 0, "history must not replay");
        assert_eq!(app.agent_cmd_cursor.get(&1), Some(&1));

        // Append a genuinely new command, force the throttle open, poll again.
        use std::io::Write as _;
        let mut f = std::fs::OpenOptions::new()
            .append(true)
            .open(&path)
            .unwrap();
        writeln!(
            f,
            "{{\"cmd\":\"new_tab\",\"name\":\"fresh\",\"workbench\":\"cad\"}}"
        )
        .unwrap();
        app.last_agent_poll = None; // bypass the 1s throttle for the test
        poll_and_apply_agent_commands(&mut app);

        assert_eq!(app.tab_bar.tabs.len(), 1, "only the new line ran");
        let idx = app.tab_bar.active.unwrap();
        assert_eq!(app.tab_bar.tabs[idx].title, "fresh");
        assert_eq!(app.tab_bar.tabs[idx].kind, TabKind::Cad);
        let _ = std::fs::remove_file(&path);
    }
}
