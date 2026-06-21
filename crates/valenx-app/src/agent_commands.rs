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
//! Stale command files are wiped **once at launch** by [`clear_command_files`]
//! (called from [`ValenxApp::new`]), so no command file left over from a
//! previous run (or the agent's own scrollback) survives into this session —
//! the only lines a channel ever holds are genuinely new commands the agent
//! appends *now*.
//!
//! Given that, [`poll_and_apply_agent_commands`] keeps a per-channel **cursor**
//! ([`ValenxApp::agent_cmd_cursor`]) of how many lines it has already applied
//! and, on the **first** poll for a channel, starts that cursor at **0** so
//! every appended command runs from the very first line. This is essential for
//! the live flow: the agent *creates* the command file by appending its
//! commands, so the first poll that sees the file is exactly the poll that must
//! run them — adopting the line count here would skip the agent's whole first
//! batch. Disk reads are throttled to ~1/sec via [`ValenxApp::last_agent_poll`].

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
    /// Render a **finished build result** into *this* unit's workspace tile
    /// (`workspace:<n>`), replacing its "the agent's output will appear here"
    /// placeholder with a result card: `title` as a bold heading over `lines`
    /// (one row each). Stored on [`crate::ValenxApp::workspace_products`] under
    /// the same channel `n` the bridge posts Notes to, so the agent's output
    /// lands in the workspace pane paired with its chat. This is the answer to
    /// "when the agent builds a rocket, show the rocket here, not a placeholder".
    ShowProduct {
        /// Card heading (rendered bold), e.g. the product name.
        title: String,
        /// Result rows shown under the heading, one per entry. Optional.
        #[serde(default)]
        lines: Vec<String>,
    },
    /// Render a **live 3-D model** into *this* unit's workspace tile
    /// (`workspace:<n>`): the tile shows the actual lit mesh (same look as the
    /// central viewport) at a fixed 3/4 camera, not a text card or a
    /// placeholder. `kind` selects the model — currently `"rocket"` (the LV-1
    /// launch vehicle); unknown kinds are skipped. Like every other command the
    /// effective channel is the command file's `n`; the optional wire `n` is
    /// accepted (and ignored) so an agent may include it for readability.
    ///
    /// Stored on [`crate::ValenxApp::workspace_products`] under the channel `n`,
    /// the same key [`ShowProduct`](AgentCommand::ShowProduct) and `Note` use,
    /// so the 3-D view lands in the workspace pane paired with its chat. This is
    /// the answer to "when the agent builds a rocket, show the *actual rocket*
    /// here — lit and in 3-D — not a card".
    ///
    /// Note the explicit `rename`: the enum's `rename_all = "snake_case"`
    /// would map `Show3d` to `"show3d"` (no underscore before the digit), but
    /// the wire tag is `"show_3d"`.
    #[serde(rename = "show_3d")]
    Show3d {
        /// Optional unit number for readability; ignored in favour of the
        /// command file's channel `n` (the bridge always routes by channel).
        #[serde(default)]
        n: Option<usize>,
        /// Which model to show. `"rocket"` → the LV-1; other values are
        /// skipped.
        kind: String,
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

/// **Wipe stale command files at launch.** Deletes every
/// `valenx_chat_cmd_u*.jsonl` in the command-channel base directory (the parent
/// of [`crate::assistant_workbench`]'s inbox path — the same dir [`cmd_path`]
/// uses) so that no command file left over from a previous run is replayed this
/// session. Best-effort: a missing dir or an un-deletable file is ignored.
///
/// This is what makes the "start the cursor at line 0" rule in
/// [`poll_and_apply_agent_commands`] safe — once the leftovers are gone, the
/// only lines a channel can contain are genuinely-new commands the agent
/// appends during *this* run, so running them all from line 0 replays nothing.
/// Call **once** at startup (see [`ValenxApp::new`]).
pub fn clear_command_files(app: &ValenxApp) {
    // Same base-dir derivation as `cmd_path` so we clean exactly the directory
    // the bridge writes into.
    let base_dir: PathBuf = app
        .assistant
        .inbox_path()
        .parent()
        .map(Path::to_path_buf)
        .unwrap_or_else(std::env::temp_dir);
    let Ok(entries) = std::fs::read_dir(&base_dir) else {
        return; // dir missing / unreadable → nothing to clean
    };
    let prefix = format!("{CMD_STEM}_u");
    for entry in entries.flatten() {
        let path = entry.path();
        // Match files named `valenx_chat_cmd_u*.jsonl`.
        let Some(name) = path.file_name().and_then(|n| n.to_str()) else {
            continue;
        };
        if name.starts_with(&prefix) && name.ends_with(".jsonl") {
            // Best-effort delete; a still-open handle elsewhere just leaves it.
            let _ = std::fs::remove_file(&path);
        }
    }
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

        // First poll for this channel: start at line 0 so EVERY appended command
        // runs. In the live flow the agent creates this file by appending its
        // commands, so the first poll that sees the file is the one that must run
        // them — there is no stale history to skip because `clear_command_files`
        // wiped any leftovers at launch.
        let start = match app.agent_cmd_cursor.get(&n) {
            Some(&c) => c.min(total), // clamp in case the file was truncated
            None => 0,
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
            // Rocket tab → load the 3-D LV-1 mesh into this (now active, empty)
            // tab's viewport immediately, so the central panel renders the
            // rocket this frame instead of the landing page. The per-tab mesh
            // starts `None` and the workbench's global first-open guard may
            // already be set from an earlier tab, so the workbench body itself
            // would never re-request the load — do it here in the reducer
            // (which runs before `show_landing` is computed). Rocket-specific
            // for now; other workbenches get their own product-load later.
            if kind == TabKind::Rocket {
                crate::rocket_workbench::ensure_lv1_3d_loaded(app);
            }
        }
        AgentCommand::OpenWorkbench { id } => {
            if let Some(kind) = TabKind::from_id(&id) {
                if let Some(idx) = app.tab_bar.active {
                    if let Some(tab) = app.tab_bar.tabs.get_mut(idx) {
                        tab.kind = kind;
                        project_tabs::sync_active(app);
                    }
                }
                // Same rocket-specific 3-D load as NewTab: switching the active
                // tab to the Rocket workbench should show the rocket model in
                // the centre, not the landing page.
                if kind == TabKind::Rocket {
                    crate::rocket_workbench::ensure_lv1_3d_loaded(app);
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
        AgentCommand::ShowProduct { title, lines } => {
            // Publish the finished result into THIS unit's workspace tile. The
            // reducer already knows the channel `n` (the same one Notes post to),
            // so the `workspace:<n>` pane and its paired `agent:<n>` chat agree.
            // A text card carries no mesh; the camera field is unused but must
            // be set (WorkspaceProduct is no longer Default).
            app.workspace_products.insert(
                n,
                crate::WorkspaceProduct {
                    title,
                    lines,
                    mesh: None,
                    camera: valenx_viz::OrbitCamera::default(),
                },
            );
        }
        AgentCommand::Show3d { n: _, kind } => {
            // Publish a LIVE 3-D model into THIS unit's workspace tile, keyed by
            // the channel `n` the reducer already knows (the wire `n` is ignored
            // — the bridge routes by file channel, same as every other command).
            // The `workspace:<n>` pane then renders the actual lit mesh at a
            // fixed 3/4 camera (see `dock_layout::render_workspace_body`).
            // Each `kind` builds its producer's LoadedMesh + a fixed 3/4 camera
            // and publishes a 3-D `WorkspaceProduct`; the render path
            // (`dock_layout::render_workspace_body`) is kind-agnostic. An
            // unknown kind is skipped safely (no panic, no placeholder churn),
            // consistent with the rest of the reducer's bad-input handling. Add
            // new model kinds as further `else if kind == "<x>"` arms.
            if kind == "rocket" {
                let mesh = crate::rocket_workbench::lv1_loaded_mesh();
                let camera = crate::rocket_workbench::lv1_camera(&mesh.mesh);
                app.workspace_products.insert(
                    n,
                    crate::WorkspaceProduct {
                        title: "Rocket".into(),
                        lines: vec![],
                        mesh: Some(mesh),
                        camera,
                    },
                );
            } else if kind == "gear" {
                let (mesh, lines) = crate::gears_workbench::gear_train_loaded_mesh();
                let camera = crate::gears_workbench::gear_train_camera(&mesh.mesh);
                app.workspace_products.insert(
                    n,
                    crate::WorkspaceProduct {
                        title: "2-stage spur reducer".into(),
                        lines,
                        mesh: Some(mesh),
                        camera,
                    },
                );
            } else if kind == "bracket" {
                let (mesh, lines) = crate::bracket_product::bracket_loaded_mesh();
                let camera = crate::bracket_product::bracket_camera(&mesh.mesh);
                app.workspace_products.insert(
                    n,
                    crate::WorkspaceProduct {
                        title: "L-bracket".into(),
                        lines,
                        mesh: Some(mesh),
                        camera,
                    },
                );
            } else if kind == "rcbeam" {
                let (mesh, lines) = crate::rcbeam_workbench::rcbeam_loaded_mesh();
                let camera = crate::rcbeam_workbench::rcbeam_camera(&mesh.mesh);
                app.workspace_products.insert(
                    n,
                    crate::WorkspaceProduct {
                        title: "RC beam (6 m, 25 kN/m)".into(),
                        lines,
                        mesh: Some(mesh),
                        camera,
                    },
                );
            } else if kind == "fem" {
                // Steel cantilever (1 m, 50×100 mm, 5 kN tip): the real
                // `valenx-fem` linear-static solve, shown as the GREY deformed
                // boundary skin (per-vertex stress colour is a deferred renderer
                // change) plus the FE-vs-analytical readout rows.
                let (mesh, lines) = crate::fem_workbench::fem_beam_loaded_mesh();
                let camera = crate::fem_workbench::fem_beam_camera(&mesh.mesh);
                app.workspace_products.insert(
                    n,
                    crate::WorkspaceProduct {
                        title: "FEM cantilever (steel, 5 kN tip)".into(),
                        lines,
                        mesh: Some(mesh),
                        camera,
                    },
                );
            } else if kind == "dna" {
                // Codon-optimised therapeutic-peptide construct + synthesis
                // screen. This is a TEXT product (mesh: None) — the workspace
                // card renders the sequence + CAI / GC / hairpin ΔG rows,
                // including the explicit "NOT a biosecurity screen" note.
                let (title, lines) = crate::dna_product::dna_construct_lines();
                app.workspace_products.insert(
                    n,
                    crate::WorkspaceProduct {
                        title,
                        lines,
                        mesh: None,
                        camera: valenx_viz::OrbitCamera::default(),
                    },
                );
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
        // The 3-D LV-1 mesh is loaded into the central viewport, so the central
        // panel's `show_landing` (project & stl & mesh all None) is false and
        // the rocket renders instead of the welcome page.
        let mesh = app
            .mesh
            .as_ref()
            .expect("rocket tab loads the 3-D mesh into the viewport");
        assert!(mesh.path.to_string_lossy().contains("valenx-lv1"));
        assert!(app.project.is_none() && app.stl.is_none()); // ⇒ show_landing == false
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

    /// Point `app`'s command base dir at a fresh, test-private temp directory
    /// (named by `tag`) so the command files this test reads/writes can't
    /// collide with another test's channel-1 file or a live app. Returns the
    /// base dir; the caller derives paths via `cmd_path`.
    fn isolate_cmd_dir(app: &mut ValenxApp, tag: &str) -> PathBuf {
        let dir = std::env::temp_dir().join(format!(
            "valenx_agentcmd_test_{}_{}",
            std::process::id(),
            tag
        ));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        app.assistant
            .set_inbox_path_for_test(dir.join("assistant_inbox.jsonl"));
        dir
    }

    #[test]
    fn first_poll_runs_all_appended_commands_from_line_zero() {
        // The live flow: a channel's command file is *created* by the agent
        // appending its commands, so the first poll that sees the file must run
        // ALL of them (cursor starts at 0). Stale-history replay is prevented by
        // `clear_command_files` at launch, not by adopting the line count here.
        let mut app = ValenxApp::default();
        app.wb_agent_counter = 1;
        let dir = isolate_cmd_dir(&mut app, "firstpoll");
        let path = cmd_path(&app, 1);
        // The agent appends two commands, creating the file.
        std::fs::write(
            &path,
            "{\"cmd\":\"new_tab\",\"name\":\"one\"}\n{\"cmd\":\"new_tab\",\"name\":\"two\"}\n",
        )
        .unwrap();

        // First poll for the channel: no cursor yet → start at 0 → both run.
        assert!(!app.agent_cmd_cursor.contains_key(&1));
        poll_and_apply_agent_commands(&mut app);

        assert_eq!(
            app.tab_bar.tabs.len(),
            2,
            "both freshly-appended commands ran on the first poll"
        );
        assert_eq!(app.tab_bar.tabs[0].title, "one");
        assert_eq!(app.tab_bar.tabs[1].title, "two");
        // Cursor advanced past everything just applied.
        assert_eq!(app.agent_cmd_cursor.get(&1), Some(&2));
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn poll_runs_only_newly_appended_lines() {
        // End-to-end through poll_and_apply_agent_commands: the first batch runs
        // in full (cursor from 0), then a later append runs ONLY the new line as
        // the cursor advances. Uses an isolated temp dir so the channel-1 file
        // is private to this test.
        let mut app = ValenxApp::default();
        app.wb_agent_counter = 1;
        let dir = isolate_cmd_dir(&mut app, "onlynew");
        let path = cmd_path(&app, 1);
        std::fs::write(&path, "{\"cmd\":\"new_tab\",\"name\":\"first\"}\n").unwrap();

        // First poll: the freshly-created file's one command runs.
        poll_and_apply_agent_commands(&mut app);
        assert_eq!(app.tab_bar.tabs.len(), 1, "first appended command ran");
        assert_eq!(app.tab_bar.tabs[0].title, "first");
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

        assert_eq!(
            app.tab_bar.tabs.len(),
            2,
            "only the new line ran the 2nd time"
        );
        let idx = app.tab_bar.active.unwrap();
        assert_eq!(app.tab_bar.tabs[idx].title, "fresh");
        assert_eq!(app.tab_bar.tabs[idx].kind, TabKind::Cad);
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn show_product_publishes_into_the_units_workspace() {
        // End-to-end through the REAL poll/reducer path: the agent appends a
        // `show_product` line on channel 1; the first poll runs it and the
        // finished result lands on `app.workspace_products[&1]` (so the
        // `workspace:1` tile renders a card instead of the placeholder).
        let mut app = ValenxApp::default();
        app.wb_agent_counter = 1;
        let dir = isolate_cmd_dir(&mut app, "showproduct");
        let path = cmd_path(&app, 1);
        std::fs::write(
            &path,
            "{\"cmd\":\"show_product\",\"title\":\"Rocket\",\"lines\":[\"thrust 1000 kN\",\"to orbit\"]}\n",
        )
        .unwrap();

        assert!(app.workspace_products.is_empty());
        poll_and_apply_agent_commands(&mut app);

        let product = app
            .workspace_products
            .get(&1)
            .expect("channel-1 product set by the show_product command");
        assert_eq!(product.title, "Rocket");
        assert_eq!(product.lines, vec!["thrust 1000 kN", "to orbit"]);
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn new_tab_rocket_loads_3d_mesh_through_real_poll_path() {
        // The bug: an agent-created "Rocket" tab showed the landing page because
        // its per-tab mesh was None while the workbench's global first-open
        // guard was already set. Fix: the NewTab reducer loads the LV-1 3-D mesh
        // for a rocket tab. Drive it end-to-end through the REAL
        // `poll_and_apply_agent_commands` path (a `new_tab` line with
        // workbench:"rocket") and assert the mesh is present and is the LV-1, so
        // the central panel's `show_landing` would be false.
        let mut app = ValenxApp::default();
        app.wb_agent_counter = 1;
        let dir = isolate_cmd_dir(&mut app, "rocket3d");
        let path = cmd_path(&app, 1);
        std::fs::write(
            &path,
            "{\"cmd\":\"new_tab\",\"name\":\"Rocket\",\"workbench\":\"rocket\"}\n",
        )
        .unwrap();

        assert!(app.mesh.is_none(), "no mesh before the command runs");
        poll_and_apply_agent_commands(&mut app);

        let idx = app.tab_bar.active.expect("rocket tab is active");
        assert_eq!(app.tab_bar.tabs[idx].kind, TabKind::Rocket);
        // The viewport now holds the LV-1 3-D mesh ⇒ show_landing is false.
        let mesh = app
            .mesh
            .as_ref()
            .expect("rocket tab loaded the 3-D mesh via the agent bridge");
        assert!(
            mesh.path.to_string_lossy().contains("valenx-lv1"),
            "loaded mesh is the LV-1 rocket (path = {:?})",
            mesh.path
        );
        assert!(
            app.project.is_none() && app.stl.is_none() && app.mesh.is_some(),
            "show_landing == false: the rocket renders, not the welcome page"
        );
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn show_3d_parses_with_and_without_n() {
        // The `show_3d` command parses from its wire form; `n` is optional.
        let s: AgentCommand = serde_json::from_str(r#"{"cmd":"show_3d","kind":"rocket"}"#).unwrap();
        assert_eq!(
            s,
            AgentCommand::Show3d {
                n: None,
                kind: "rocket".into(),
            }
        );
        let s2: AgentCommand =
            serde_json::from_str(r#"{"cmd":"show_3d","n":2,"kind":"rocket"}"#).unwrap();
        assert_eq!(
            s2,
            AgentCommand::Show3d {
                n: Some(2),
                kind: "rocket".into(),
            }
        );
    }

    #[test]
    fn show_3d_rocket_publishes_a_live_mesh_into_the_workspace() {
        // End-to-end through the REAL poll/reducer path: the agent appends a
        // `show_3d` rocket line on channel 1; the first poll runs it and the
        // unit's workspace product gains a live LV-1 mesh (so the
        // `workspace:1` tile renders the actual 3-D rocket, not a card).
        let mut app = ValenxApp::default();
        app.wb_agent_counter = 1;
        let dir = isolate_cmd_dir(&mut app, "show3d_rocket");
        let path = cmd_path(&app, 1);
        std::fs::write(&path, "{\"cmd\":\"show_3d\",\"kind\":\"rocket\"}\n").unwrap();

        assert!(app.workspace_products.is_empty());
        poll_and_apply_agent_commands(&mut app);

        let product = app
            .workspace_products
            .get(&1)
            .expect("channel-1 product set by the show_3d command");
        let mesh = product
            .mesh
            .as_ref()
            .expect("show_3d rocket attaches a live LoadedMesh");
        assert!(
            mesh.path.to_string_lossy().contains("valenx-lv1"),
            "the attached mesh is the LV-1 rocket (path = {:?})",
            mesh.path
        );
        // It's a 3-D product, so no text rows.
        assert_eq!(product.title, "Rocket");
        assert!(product.lines.is_empty());
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// Drive each newly-wired 3-D `kind` (`gear` / `bracket` / `rcbeam`)
    /// end-to-end through the REAL poll/reducer path and assert the unit's
    /// workspace product gains a live mesh with a non-empty `Tri3` surface
    /// plus its numeric readout rows — the same contract the rocket meets.
    fn assert_show_3d_kind_publishes_live_mesh(kind: &str, tag: &str, expect_line_substr: &str) {
        let mut app = ValenxApp::default();
        app.wb_agent_counter = 1;
        let dir = isolate_cmd_dir(&mut app, tag);
        let path = cmd_path(&app, 1);
        std::fs::write(
            &path,
            format!("{{\"cmd\":\"show_3d\",\"kind\":\"{kind}\"}}\n"),
        )
        .unwrap();

        assert!(app.workspace_products.is_empty());
        poll_and_apply_agent_commands(&mut app);

        let product = app
            .workspace_products
            .get(&1)
            .unwrap_or_else(|| panic!("channel-1 product set by show_3d {kind}"));
        let mesh = product
            .mesh
            .as_ref()
            .unwrap_or_else(|| panic!("show_3d {kind} attaches a live LoadedMesh"));
        // The mesh is non-empty and carries surface (Tri3) elements.
        assert!(!mesh.mesh.nodes.is_empty(), "{kind}: mesh has vertices");
        assert!(mesh.mesh.total_elements() > 0, "{kind}: mesh has triangles");
        // Tagged with its synthetic source path.
        assert!(
            mesh.path.to_string_lossy().contains(tag),
            "{kind}: path = {:?}",
            mesh.path
        );
        // 3-D products built from a producer carry numeric readout rows.
        assert!(
            product.lines.iter().any(|l| l.contains(expect_line_substr)),
            "{kind}: expected a line containing {expect_line_substr:?}, got {:?}",
            product.lines
        );
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn show_3d_gear_publishes_a_live_mesh_into_the_workspace() {
        assert_show_3d_kind_publishes_live_mesh("gear", "valenx-2stage", "ratio");
    }

    #[test]
    fn show_3d_bracket_publishes_a_live_mesh_into_the_workspace() {
        assert_show_3d_kind_publishes_live_mesh("bracket", "valenx-l-bracket", "M5");
    }

    #[test]
    fn show_3d_rcbeam_publishes_a_live_mesh_into_the_workspace() {
        assert_show_3d_kind_publishes_live_mesh("rcbeam", "valenx-rcbeam", "Mn");
    }

    #[test]
    fn show_3d_fem_publishes_a_live_mesh_into_the_workspace() {
        // The steel cantilever ships a non-empty Tri3 boundary skin (the deformed
        // shape, grey) plus a readout row reporting the FE max deflection.
        assert_show_3d_kind_publishes_live_mesh("fem", "valenx-fem-cantilever", "max deflection");
    }

    #[test]
    fn show_3d_dna_publishes_a_text_card_into_the_workspace() {
        // End-to-end through the REAL poll/reducer path: the agent appends a
        // `show_3d` dna line on channel 1; the first poll runs it and the unit's
        // workspace product is a TEXT card (mesh: None) with the codon-optimised
        // construct, a GC / CAI line, and the explicit no-biosecurity-screen note.
        let mut app = ValenxApp::default();
        app.wb_agent_counter = 1;
        let dir = isolate_cmd_dir(&mut app, "show3d_dna");
        let path = cmd_path(&app, 1);
        std::fs::write(&path, "{\"cmd\":\"show_3d\",\"kind\":\"dna\"}\n").unwrap();

        assert!(app.workspace_products.is_empty());
        poll_and_apply_agent_commands(&mut app);

        let product = app
            .workspace_products
            .get(&1)
            .expect("channel-1 product set by show_3d dna");
        // A text card: no mesh, non-empty rows.
        assert!(product.mesh.is_none(), "dna is a text card (no mesh)");
        assert!(!product.lines.is_empty(), "dna card has rows");
        // A GC / CAI line is present.
        assert!(
            product.lines.iter().any(|l| l.contains("GC content")),
            "a GC content line is present: {:?}",
            product.lines
        );
        assert!(
            product.lines.iter().any(|l| l.starts_with("CAI")),
            "a CAI line is present: {:?}",
            product.lines
        );
        // The explicit honesty note that this is NOT a biosecurity screen.
        assert!(
            product
                .lines
                .iter()
                .any(|l| l.contains("NOT a biosecurity screen")),
            "the no-biosecurity-screen note is present: {:?}",
            product.lines
        );
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn show_3d_unknown_kind_is_skipped() {
        // An unknown `kind` is a safe no-op (no product inserted, no panic),
        // consistent with the reducer's bad-input handling elsewhere.
        let mut app = ValenxApp::default();
        apply(
            &mut app,
            1,
            AgentCommand::Show3d {
                n: None,
                kind: "not-a-model".into(),
            },
        );
        assert!(app.workspace_products.get(&1).is_none());
    }

    #[test]
    fn show_product_lines_default_to_empty() {
        // `lines` is optional on the wire — a title-only card parses fine.
        let sp: AgentCommand =
            serde_json::from_str(r#"{"cmd":"show_product","title":"Gear train"}"#).unwrap();
        assert_eq!(
            sp,
            AgentCommand::ShowProduct {
                title: "Gear train".into(),
                lines: vec![],
            }
        );
    }

    #[test]
    fn clear_command_files_removes_matching_files() {
        // `clear_command_files` deletes every `valenx_chat_cmd_u*.jsonl` in the
        // command-channel base dir (the inbox path's parent) and leaves other
        // files alone. Build an app whose inbox path lives in a fresh temp dir
        // so the scan targets a directory we control.
        let dir =
            std::env::temp_dir().join(format!("valenx_clear_cmd_test_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();

        // Two command files (different channels) + one unrelated file.
        let cmd1 = dir.join(format!("{CMD_STEM}_u1.jsonl"));
        let cmd2 = dir.join(format!("{CMD_STEM}_u2.jsonl"));
        let keep = dir.join("assistant_feed.jsonl");
        std::fs::write(&cmd1, "{\"cmd\":\"note\",\"text\":\"x\"}\n").unwrap();
        std::fs::write(&cmd2, "{\"cmd\":\"note\",\"text\":\"y\"}\n").unwrap();
        std::fs::write(&keep, "unrelated\n").unwrap();

        // Point the app's assistant inbox at a path inside `dir` so its parent
        // (the base dir) is exactly `dir`.
        let mut app = ValenxApp::default();
        app.assistant
            .set_inbox_path_for_test(dir.join("assistant_inbox.jsonl"));
        assert_eq!(cmd_path(&app, 1), cmd1, "test base dir wired up correctly");

        clear_command_files(&app);

        assert!(!cmd1.exists(), "channel-1 command file deleted");
        assert!(!cmd2.exists(), "channel-2 command file deleted");
        assert!(keep.exists(), "unrelated file left intact");

        let _ = std::fs::remove_dir_all(&dir);
    }
}
