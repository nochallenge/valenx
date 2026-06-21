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
///
/// Public so the eframe `update()` loop can schedule an unconditional
/// heartbeat repaint on the *same* cadence (see `update.rs`). egui is
/// reactive — when valenx is idle/unfocused (the normal case while an
/// external agent drives it from the background) `update()` would otherwise
/// stop being called, so [`poll_and_apply_agent_commands`] would never run and
/// appended commands (incl. the **global** `new_unit`) would sit unprocessed.
/// Keeping the heartbeat in lockstep with this interval guarantees the poll
/// fires every ~1 s regardless of focus.
pub const POLL_INTERVAL: Duration = Duration::from_millis(1000);

/// Fixed stem for a per-channel command file (the `_u{n}.jsonl` suffix is
/// appended per channel by [`cmd_path`]). The **global** channel file (read by
/// [`apply_global`] for [`NewUnit`](AgentCommand::NewUnit)) is this stem with
/// **no** suffix: `<base-dir>/valenx_chat_cmd.jsonl` (see [`global_cmd_path`]).
const CMD_STEM: &str = "valenx_chat_cmd";

/// Upper bound on how many Workbench+Agent units the global
/// [`NewUnit`](AgentCommand::NewUnit) command will open, so a runaway / hostile
/// command file cannot spawn unbounded panes. A `new_unit` arriving once
/// [`crate::ValenxApp::wb_agent_counter`] has reached this is ignored.
const MAX_UNITS: usize = 200;

/// One command an external agent can append to its channel's command file to
/// drive valenx. **Internally tagged** on `"cmd"`, so each line is a flat
/// object like `{"cmd":"new_tab","name":"…"}`.
///
/// Every variant is honoured through an *existing* vetted tab/dock method (see
/// `apply`); none writes a raw [`ValenxApp`] field. Unknown command tags or
/// bad workbench ids are skipped without panicking.
// `Eq` is intentionally NOT derived: the `animate` command carries an
// `Option<f32>` (`speed`), and `f32` is not `Eq`. `PartialEq` is enough for the
// tests' `assert_eq!` round-trips, and nothing uses `AgentCommand` as a hash /
// set key, so dropping `Eq` is free.
#[derive(Debug, Clone, PartialEq, Deserialize)]
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
    /// Render a **2-D engineering drawing** into *this* unit's workspace tile
    /// (`workspace:<n>`): the tile paints a flat egui drawing (no wgpu) — these
    /// are the user's originally-spec'd outputs. `kind` selects the drawing:
    /// `"rcbeam"` → a reinforced-concrete **section + rebar** drawing (the
    /// canonical 300×550 section), `"dna"` → a DNA **construct map** (the
    /// 93-nt ATG·ORF·His6·stop construct as coloured feature blocks on a nt
    /// ruler). Unknown kinds are skipped. These are **distinct** from the
    /// `show_3d` views — `show_3d{kind:"rcbeam"}` is the lit 3-D solid and
    /// `show_3d{kind:"dna"}` is the text card; `show_2d` adds the 2-D drawings,
    /// so both coexist.
    ///
    /// As with [`Show3d`](AgentCommand::Show3d) the effective channel is the
    /// command file's `n`; the optional wire `n` is accepted (and ignored) for
    /// readability. The enum's `rename_all = "snake_case"` would map `Show2d`
    /// to `"show2d"` (no underscore before the digit), so the wire tag is
    /// pinned to `"show_2d"`.
    #[serde(rename = "show_2d")]
    Show2d {
        /// Optional unit number for readability; ignored in favour of the
        /// command file's channel `n` (the bridge always routes by channel).
        #[serde(default)]
        n: Option<usize>,
        /// Which 2-D drawing to show. `"rcbeam"` → the RC section + rebar,
        /// `"dna"` → the DNA construct map; other values are skipped.
        kind: String,
    },
    /// **Open a brand-new "Workbench + Agent" unit** and (optionally) build a
    /// product into it — entirely from the file, *no UI click*. This is the
    /// **bootstrap** command an external agent uses to mint its own tab before
    /// any unit exists, so it is honoured **only on the global channel**
    /// (`<base-dir>/valenx_chat_cmd.jsonl`, no `_u` suffix) by the `apply_global`
    /// handler; appearing on a per-unit file it is parsed but ignored (a
    /// per-unit channel already *is* a unit). Drives the very same path the
    /// "+ Workbench+Agent → New row at bottom" button uses
    /// (`ValenxApp::add_workbench_agent_pair_at` with `UnitAddTarget::NewRowBottom`),
    /// which bumps [`crate::ValenxApp::wb_agent_counter`] to the new unit `n`.
    ///
    /// Then, on that new unit `n`:
    /// - if `kind` is set, the named product is rendered into the unit's
    ///   `workspace:<n>` tile via the **same** `show_3d` / `show_2d` reducer
    ///   paths a running agent would use (registry mesh, the `dna` text card,
    ///   or a 2-D drawing — unknown kinds render nothing, no panic);
    /// - if `title` is set, it overrides the rendered product's heading (or, if
    ///   no product was rendered, a title-only card is shown so the workspace
    ///   names itself);
    /// - if `note` is set, it is posted to unit `n`'s chat feed (the same
    ///   `append_feed_note` path the [`Note`](AgentCommand::Note) command uses)
    ///   so the agent's narration shows up;
    /// - a `"Unit <n> ready"` confirmation note is **always** posted so the
    ///   agent can detect the unit opened.
    ///
    /// Bounded: refused once the `MAX_UNITS` cap of units exist, so a runaway
    /// command file can't spawn unbounded panes.
    NewUnit {
        /// Optional product to build into the new unit, by the same id
        /// `show_3d` / `show_2d` accept (e.g. `"rocket"`, `"gear"`, `"rcbeam"`,
        /// `"dna"`). Absent → an empty unit.
        #[serde(default)]
        kind: Option<String>,
        /// Optional heading for the unit's workspace product card.
        #[serde(default)]
        title: Option<String>,
        /// Optional narration line posted to the new unit's chat feed.
        #[serde(default)]
        note: Option<String>,
    },
    /// **Drive an animated product's playback clock** on a unit whose
    /// `workspace:<n>` product carries a [`crate::ProductAnimation`] (e.g. the
    /// meshing gear train) — the file-driven equivalent of clicking the tile's
    /// Play/Pause button and dragging its speed slider. `play` sets the
    /// playing/paused state and `speed` sets the playback multiplier (clamped to
    /// the toolbar's `0.0..=8.0`); either may be omitted to leave that field
    /// untouched. The effective unit is `n` when given, else the command file's
    /// channel — so it works on a per-unit channel (`n` absent) and on the
    /// global channel (where an explicit `n` is required to pick a unit). A unit
    /// with no animated product just gets a feed note saying so; nothing panics.
    #[serde(rename = "animate")]
    Animate {
        /// Optional target unit; defaults to the command file's channel.
        #[serde(default)]
        n: Option<usize>,
        /// When set, play (`true`) or pause (`false`) the product's clock.
        #[serde(default)]
        play: Option<bool>,
        /// When set, the playback-speed multiplier (clamped to `0.0..=8.0`).
        #[serde(default)]
        speed: Option<f32>,
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

/// The **global** command-channel file: `<base-dir>/valenx_chat_cmd.jsonl` (the
/// `valenx_chat_cmd` stem with **no** `_u{n}` suffix), in the same base
/// directory as the per-unit channels ([`cmd_path`]). Polled on **every** poll
/// regardless of
/// [`crate::ValenxApp::wb_agent_counter`] so an external agent can append a
/// [`NewUnit`](AgentCommand::NewUnit) to open its own unit before any unit
/// exists — the entry point for agent-per-tab product generation.
pub fn global_cmd_path(app: &ValenxApp) -> PathBuf {
    let base_dir: PathBuf = app
        .assistant
        .inbox_path()
        .parent()
        .map(Path::to_path_buf)
        .unwrap_or_else(std::env::temp_dir);
    base_dir.join(format!("{CMD_STEM}.jsonl"))
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
    // The global (no-`_u`-suffix) channel file is wiped too, so a stale
    // `new_unit` from a previous run is never replayed at launch.
    let global_name = format!("{CMD_STEM}.jsonl");
    for entry in entries.flatten() {
        let path = entry.path();
        // Match the per-unit `valenx_chat_cmd_u*.jsonl` files and the global
        // `valenx_chat_cmd.jsonl` file.
        let Some(name) = path.file_name().and_then(|n| n.to_str()) else {
            continue;
        };
        if (name.starts_with(&prefix) && name.ends_with(".jsonl")) || name == global_name {
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

    // GLOBAL channel first: read `<base-dir>/valenx_chat_cmd.jsonl` (no `_u`
    // suffix) on EVERY poll, NOT gated on `wb_agent_counter`, so an external
    // agent can `new_unit` to bootstrap its own Workbench+Agent unit before any
    // unit exists. Mirrors the per-unit append-only cursor read below, but with
    // its own persistent cursor `agent_global_cmd_cursor` and dispatch through
    // `apply_global`.
    {
        let path = global_cmd_path(app);
        if let Some(body) = read_cmd_file(&path) {
            let total = line_count(&body);
            // First poll that sees the file → start at line 0 so EVERY appended
            // command runs (stale history was wiped by `clear_command_files` at
            // launch). Thereafter only genuinely-new lines run; clamp in case
            // the file was truncated.
            let start = match app.agent_global_cmd_cursor {
                Some(c) => c.min(total),
                None => 0,
            };
            if start < total {
                for line in body
                    .lines()
                    .map(str::trim)
                    .filter(|l| !l.is_empty())
                    .skip(start)
                {
                    if let Ok(cmd) = serde_json::from_str::<AgentCommand>(line) {
                        apply_global(app, cmd);
                    }
                    // Unparseable line → skip safely.
                }
                app.agent_global_cmd_cursor = Some(total);
            }
        }
    }

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

/// Apply **one** command from the **global** channel
/// (`<base-dir>/valenx_chat_cmd.jsonl`). Only [`NewUnit`](AgentCommand::NewUnit)
/// is meaningful here — it opens a fresh Workbench+Agent unit (the bootstrap an
/// agent needs before any unit exists) and optionally builds a product into it.
/// Every other variant is a per-*unit* command and is ignored on the global
/// channel (it has no unit to act on); a malformed line never reaches here (the
/// poll loop skips unparseable lines). Like [`apply`], every branch is a
/// no-op-on-bad-input rather than a panic.
fn apply_global(app: &mut ValenxApp, cmd: AgentCommand) {
    let AgentCommand::NewUnit { kind, title, note } = cmd else {
        // Non-`new_unit` commands are per-unit; the global channel has no unit
        // to target, so they are ignored here.
        return;
    };

    // Bound the unit count so a runaway / hostile command file can't spawn
    // unbounded panes.
    if app.wb_agent_counter >= MAX_UNITS {
        return;
    }

    // Open the new unit in its OWN top-strip project tab so each agent unit
    // lands in an isolated workspace (its own dock tree), instead of stacking
    // another pane into whatever tab is current. Mirror the `NewTab` reducer's
    // park → open → title → install dance EXACTLY: park the outgoing tab's scene
    // FIRST, open a fresh Blank tab (pushes a doc + makes it active), title it,
    // then install that empty doc so the prior tab keeps its geometry and this
    // one starts clean. Only AFTER the install do we add the Workbench+Agent
    // pair, which builds this (now-active, empty) tab's dock tree — so routing
    // by the global unit `n` (cmd_u{n} / feed_u{n} / workspace_products[n]) is
    // unchanged. (Crucially NOT the rocket-specific `ensure_lv1_3d_loaded`
    // branch — a unit's product is rendered into its workspace tile, not the
    // central viewport.)
    project_tabs::park_active_doc(app);
    app.tab_bar.open(TabKind::Blank);
    // Resolve the new tab's title: the caller's `title`, else the product
    // `kind` (its label), else `None` → a placeholder rewritten with the unit
    // number once it's known below.
    let resolved_title: Option<String> = title
        .as_deref()
        .or(kind.as_deref())
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .map(str::to_string);
    let title_is_placeholder = resolved_title.is_none();
    if let Some(idx) = app.tab_bar.active {
        if let Some(tab) = app.tab_bar.tabs.get_mut(idx) {
            tab.title = resolved_title.unwrap_or_else(|| "Workbench".to_string());
        }
    }
    project_tabs::install_active_doc(app);

    // Add the unit through the SAME path the "+ Workbench+Agent → New row at
    // bottom" button uses; this bumps `wb_agent_counter` to the new unit `n` and
    // builds the pair into this tab's (now-installed, empty) dock tree.
    app.add_workbench_agent_pair_at(crate::dock_layout::UnitAddTarget::NewRowBottom);
    let n = app.wb_agent_counter;

    // If the tab title was a pure placeholder (no caller `title`, no product
    // `kind`), rewrite it now that the unit number is known so the tab names
    // itself rather than reading "Workbench".
    if title_is_placeholder {
        if let Some(idx) = app.tab_bar.active {
            if let Some(tab) = app.tab_bar.tabs.get_mut(idx) {
                tab.title = format!("Unit {n}");
            }
        }
    }

    // If a product kind was named, render it into the new unit's `workspace:<n>`
    // tile by delegating to the SAME reducer paths a running agent would use:
    // first the `show_3d` path (registry meshes + the `dna` text card), and if
    // that produced nothing, the `show_2d` path (the 2-D-only drawings). An
    // unknown kind renders nothing — a safe no-op, like the reducer elsewhere.
    if let Some(kind) = kind {
        apply(
            app,
            n,
            AgentCommand::Show3d {
                n: None,
                kind: kind.clone(),
            },
        );
        if !app.workspace_products.contains_key(&n) {
            apply(app, n, AgentCommand::Show2d { n: None, kind });
        }
        // Ensure the freshly-rendered product carries a default inspect-spin so
        // the new unit's tile shows the Play/Pause + speed controls. The
        // `Show3d` arm above already does this for registry meshes; this keeps
        // the guarantee at the `new_unit` render path too (idempotent — a no-op
        // once the product already animates or is mesh-less).
        if let Some(product) = app.workspace_products.get_mut(&n) {
            product.ensure_default_animation();
        }
    }

    // If a title was given, use it as the workspace card heading: override a
    // rendered product's title, or (when no product rendered) show a title-only
    // card so the workspace names itself.
    if let Some(title) = title {
        let trimmed = title.trim();
        if !trimmed.is_empty() {
            match app.workspace_products.get_mut(&n) {
                Some(product) => product.title = trimmed.to_string(),
                None => {
                    app.workspace_products.insert(
                        n,
                        crate::WorkspaceProduct {
                            title: trimmed.to_string(),
                            lines: Vec::new(),
                            mesh: None,
                            vertex_colors: None,
                            camera: valenx_viz::OrbitCamera::default(),
                            kind2d: None,
                            last_export: None,
                            image: None,
                            image_texture: None,
                            animation: None,
                        },
                    );
                }
            }
        }
    }

    // Post the agent's narration (if any) to the new unit's chat feed, through
    // the same path the `note` command uses.
    if let Some(note) = note {
        let trimmed = note.trim();
        if !trimmed.is_empty() {
            crate::assistant_workbench::append_feed_note(app, n, "Claude", trimmed, "build");
        }
    }

    // ALWAYS confirm the unit opened, so an agent polling the feed can detect
    // its new unit is live.
    crate::assistant_workbench::append_feed_note(
        app,
        n,
        "Claude",
        &format!("Unit {n} ready"),
        "ship",
    );
}

/// Apply **one** [`AgentCommand`] for channel `ch` through existing vetted
/// methods only. Every branch is a no-op-on-bad-input (unknown workbench id,
/// no active tab, no title match) rather than a panic.
fn apply(app: &mut ValenxApp, ch: usize, cmd: AgentCommand) {
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
                ch,
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
            // reducer already knows the channel `ch` (the same one Notes post to),
            // so the `workspace:<n>` pane and its paired `agent:<n>` chat agree.
            // A text card carries no mesh; the camera field is unused but must
            // be set (WorkspaceProduct is no longer Default).
            app.workspace_products.insert(
                ch,
                crate::WorkspaceProduct {
                    title,
                    lines,
                    mesh: None,
                    vertex_colors: None,
                    camera: valenx_viz::OrbitCamera::default(),
                    kind2d: None,
                    last_export: None,
                    image: None,
                    image_texture: None,
                    animation: None,
                },
            );
        }
        AgentCommand::Show3d { n: _, kind } => {
            // Publish a LIVE 3-D model into THIS unit's workspace tile, keyed by
            // the channel `ch` the reducer already knows (the wire `n` is ignored
            // — the bridge routes by file channel, same as every other command).
            // The `workspace:<n>` pane then renders the actual lit mesh at a
            // fixed 3/4 camera (see `dock_layout::render_workspace_body`).
            //
            // The mesh kinds (rocket / gear / bracket / rcbeam / fem) are looked
            // up in the per-file registry ([`crate::products_registry`]) instead
            // of an inline per-kind chain: each tool registers its own
            // pure builder in its own module, so a new 3-D tool is added there
            // (plus one table line) without touching this shared reducer. An
            // unknown kind resolves to `None` and is skipped safely (no panic,
            // no placeholder churn), consistent with the reducer's other
            // bad-input handling.
            if let Some(entry) = crate::products_registry::lookup(&kind) {
                let mut product = (entry.build)();
                // Give every bridge-rendered mesh product a default paused
                // inspect-spin (Turntable about +Z through the AABB centre) so
                // the tile shows the Play/Pause + speed controls. A no-op for a
                // product that already animates (the gear's RigidParts) or that
                // is mesh-less, so it is safe to call unconditionally here.
                product.ensure_default_animation();
                app.workspace_products.insert(ch, product);
            } else if kind == "dna" {
                // `dna` is NOT a registry 3-D mesh kind — it's a TEXT product
                // (mesh: None): the codon-optimised therapeutic-peptide
                // construct + synthesis screen, rendered as a card with the
                // sequence + CAI / GC / hairpin ΔG rows, including the explicit
                // "NOT a biosecurity screen" note. Kept inline here (migrate
                // later); the registry is 3-D-mesh-only for now.
                let (title, lines) = crate::dna_product::dna_construct_lines();
                app.workspace_products.insert(
                    ch,
                    crate::WorkspaceProduct {
                        title,
                        lines,
                        mesh: None,
                        vertex_colors: None,
                        camera: valenx_viz::OrbitCamera::default(),
                        kind2d: None,
                        last_export: None,
                        image: None,
                        image_texture: None,
                        animation: None,
                    },
                );
            }
        }
        AgentCommand::Show2d { n: _, kind } => {
            // Publish a 2-D engineering DRAWING into THIS unit's workspace tile,
            // keyed by the channel `ch` the reducer already knows (the wire `n`
            // is ignored — the bridge routes by file channel, like every other
            // command). The `workspace:<n>` pane then paints the flat egui
            // drawing (see `dock_layout::render_workspace_body`'s 2-D branch).
            // These are the user's originally-spec'd outputs and are distinct
            // from the `show_3d` views — both can coexist for the same kind. An
            // unknown kind is skipped safely (no panic, no placeholder churn).
            if kind == "rcbeam" {
                // Canonical 300×550 RC section + rebar, with the flexural rows.
                let (view, lines) = crate::rcbeam_workbench::rcbeam_section_view();
                app.workspace_products.insert(
                    ch,
                    crate::WorkspaceProduct {
                        title: "RC Beam — section".into(),
                        lines,
                        mesh: None,
                        vertex_colors: None,
                        camera: valenx_viz::OrbitCamera::default(),
                        kind2d: Some(crate::Workspace2dKind::RcSection(view)),
                        last_export: None,
                        image: None,
                        image_texture: None,
                        animation: None,
                    },
                );
            } else if kind == "dna" {
                // The 93-nt ATG·ORF·His6·stop construct as a feature map, with
                // the sequence / CAI / GC rows (incl. the no-biosecurity note).
                let (map, lines) = crate::dna_product::dna_construct_map();
                app.workspace_products.insert(
                    ch,
                    crate::WorkspaceProduct {
                        title: "DNA Construct — map".into(),
                        lines,
                        mesh: None,
                        vertex_colors: None,
                        camera: valenx_viz::OrbitCamera::default(),
                        kind2d: Some(crate::Workspace2dKind::DnaMap(map)),
                        last_export: None,
                        image: None,
                        image_texture: None,
                        animation: None,
                    },
                );
            }
        }
        AgentCommand::Animate { n, play, speed } => {
            // Drive the playback clock of THIS unit's animated product (the
            // file-driven Play/Pause + speed-slider). Route by `n` when the
            // command carries one (so it also works arriving on the global
            // channel), else the command file's channel `n`. Same vetted state a
            // user toolbar click would set; a missing animation is a no-op note,
            // never a panic.
            let target = n.unwrap_or(ch);
            if let Some(p) = app.workspace_products.get_mut(&target) {
                match p.animation.as_mut() {
                    Some(a) => {
                        if let Some(v) = play {
                            a.playing = v;
                        }
                        if let Some(s) = speed {
                            a.speed = s.clamp(0.0, 8.0);
                        }
                        let state = if a.playing { "playing" } else { "paused" };
                        crate::assistant_workbench::append_feed_note(
                            app,
                            target,
                            "Claude",
                            &format!("Animation {state}"),
                            "result",
                        );
                    }
                    None => {
                        crate::assistant_workbench::append_feed_note(
                            app,
                            target,
                            "Claude",
                            "(no animation on this product)",
                            "warn",
                        );
                    }
                }
            }
        }
        AgentCommand::NewUnit { .. } => {
            // `new_unit` is a **global-channel bootstrap** command (it opens a
            // brand-new unit), handled by `apply_global`. On a per-unit channel
            // there is no new unit to open, so it is a deliberate no-op here.
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
    fn show_3d_render_publishes_an_image_into_the_workspace() {
        // End-to-end through the REAL poll/reducer path: `show_3d` render on
        // channel 1 routes through the registry to the path-traced IMAGE
        // product — no mesh / no 2-D drawing, a non-empty `egui::ColorImage`.
        let mut app = ValenxApp::default();
        app.wb_agent_counter = 1;
        let dir = isolate_cmd_dir(&mut app, "show3d_render");
        let path = cmd_path(&app, 1);
        std::fs::write(&path, "{\"cmd\":\"show_3d\",\"kind\":\"render\"}\n").unwrap();

        assert!(app.workspace_products.is_empty());
        poll_and_apply_agent_commands(&mut app);

        let product = app
            .workspace_products
            .get(&1)
            .expect("channel-1 product set by show_3d render");
        assert!(product.mesh.is_none(), "render is an image (no mesh)");
        assert!(
            product.kind2d.is_none(),
            "render is an image (not a 2-D drawing)"
        );
        let image = product.image.as_ref().expect("render carries a ColorImage");
        let [w, h] = image.size;
        assert!(w > 0 && h > 0, "render image has non-zero size");
        assert_eq!(image.pixels.len(), w * h, "ColorImage pixel count = w·h");
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn show_3d_animate_publishes_a_timeline_card_into_the_workspace() {
        // `show_3d` animate routes through the registry to a DATA-ONLY text
        // card summarising the keyframe timeline (mesh None, kind2d None, rows
        // mentioning the keyframe count + duration).
        let mut app = ValenxApp::default();
        app.wb_agent_counter = 1;
        let dir = isolate_cmd_dir(&mut app, "show3d_animate");
        let path = cmd_path(&app, 1);
        std::fs::write(&path, "{\"cmd\":\"show_3d\",\"kind\":\"animate\"}\n").unwrap();

        poll_and_apply_agent_commands(&mut app);

        let product = app
            .workspace_products
            .get(&1)
            .expect("channel-1 product set by show_3d animate");
        assert!(product.mesh.is_none(), "animate is a text card (no mesh)");
        assert!(product.image.is_none(), "animate is a text card (no image)");
        assert!(
            product.kind2d.is_none(),
            "animate is a text card (not a 2-D drawing)"
        );
        assert!(
            product.lines.iter().any(|l| l.contains("keyframes")),
            "the card reports the keyframe count: {:?}",
            product.lines
        );
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn show_3d_draft2d_and_interior_publish_2d_drawings_into_the_workspace() {
        // `show_3d` draft2d / interior route through the registry to 2-D
        // DRAWING products (kind2d Some, mesh None, image None) — the same egui
        // 2-D branch as rcbeam / dna, with the matching plain-data view.
        for (kind, tag) in [
            ("draft2d", "show3d_draft2d"),
            ("interior", "show3d_interior"),
        ] {
            let mut app = ValenxApp::default();
            app.wb_agent_counter = 1;
            let dir = isolate_cmd_dir(&mut app, tag);
            let path = cmd_path(&app, 1);
            std::fs::write(
                &path,
                format!("{{\"cmd\":\"show_3d\",\"kind\":\"{kind}\"}}\n"),
            )
            .unwrap();

            poll_and_apply_agent_commands(&mut app);

            let product = app
                .workspace_products
                .get(&1)
                .unwrap_or_else(|| panic!("channel-1 product set by show_3d {kind}"));
            assert!(product.mesh.is_none(), "{kind} is a 2-D drawing (no mesh)");
            assert!(
                product.image.is_none(),
                "{kind} is a 2-D drawing (no image)"
            );
            match product.kind2d.as_ref() {
                Some(crate::Workspace2dKind::Draft2d(view)) if kind == "draft2d" => {
                    assert!(!view.entities.is_empty(), "draft2d drawing has entities");
                }
                Some(crate::Workspace2dKind::FloorPlan(plan)) if kind == "interior" => {
                    assert!(!plan.rooms.is_empty(), "interior plan has a room");
                }
                other => panic!("{kind}: unexpected kind2d {other:?}"),
            }
            let _ = std::fs::remove_dir_all(&dir);
        }
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
        assert!(!app.workspace_products.contains_key(&1));
    }

    #[test]
    fn show_3d_mesh_kinds_route_through_the_per_file_registry() {
        // The reducer's `show_3d` mesh path delegates to
        // `products_registry::lookup`. Assert the registry resolves every
        // migrated 3-D kind (and that an unknown kind returns None, the no-op
        // the reducer relies on), so the bridge and the registry stay in sync.
        for kind in ["rocket", "gear", "bracket", "rcbeam", "fem"] {
            assert!(
                crate::products_registry::lookup(kind).is_some(),
                "registry resolves the migrated 3-D kind {kind:?}"
            );
        }
        assert!(
            crate::products_registry::lookup("not-a-model").is_none(),
            "an unknown kind returns None (the reducer then skips it)"
        );
        // `dna` is a text card, not a registry 3-D mesh kind.
        assert!(crate::products_registry::lookup("dna").is_none());
    }

    #[test]
    fn show_2d_parses_with_and_without_n() {
        // `show_2d` parses from its wire form; `n` is optional.
        let s: AgentCommand = serde_json::from_str(r#"{"cmd":"show_2d","kind":"rcbeam"}"#).unwrap();
        assert_eq!(
            s,
            AgentCommand::Show2d {
                n: None,
                kind: "rcbeam".into(),
            }
        );
        let s2: AgentCommand =
            serde_json::from_str(r#"{"cmd":"show_2d","n":3,"kind":"dna"}"#).unwrap();
        assert_eq!(
            s2,
            AgentCommand::Show2d {
                n: Some(3),
                kind: "dna".into(),
            }
        );
    }

    #[test]
    fn show_2d_rcbeam_publishes_a_section_drawing_into_the_workspace() {
        // End-to-end through the REAL poll/reducer path: the agent appends a
        // `show_2d` rcbeam line on channel 1; the first poll runs it and the
        // unit's workspace product is a 2-D RC SECTION drawing (kind2d =
        // Some(RcSection), mesh None) carrying the flexural rows.
        let mut app = ValenxApp::default();
        app.wb_agent_counter = 1;
        let dir = isolate_cmd_dir(&mut app, "show2d_rcbeam");
        let path = cmd_path(&app, 1);
        std::fs::write(&path, "{\"cmd\":\"show_2d\",\"kind\":\"rcbeam\"}\n").unwrap();

        assert!(app.workspace_products.is_empty());
        poll_and_apply_agent_commands(&mut app);

        let product = app
            .workspace_products
            .get(&1)
            .expect("channel-1 product set by show_2d rcbeam");
        // It's a 2-D drawing: no mesh, a RcSection kind2d.
        assert!(product.mesh.is_none(), "rcbeam 2-D carries no 3-D mesh");
        match product.kind2d.as_ref() {
            Some(crate::Workspace2dKind::RcSection(view)) => {
                assert_eq!(view.width_mm, 300.0);
                assert_eq!(view.n_bars, 3);
                assert!(view.bar_dia_mm > 0.0);
                // The flexural rows came along.
                assert!(
                    product.lines.iter().any(|l| l.contains("nominal Mn")),
                    "section drawing carries the Mn row: {:?}",
                    product.lines
                );
            }
            other => panic!("expected RcSection kind2d, got {other:?}"),
        }
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn show_2d_dna_publishes_a_construct_map_into_the_workspace() {
        // End-to-end through the REAL poll/reducer path: the agent appends a
        // `show_2d` dna line on channel 1; the first poll runs it and the unit's
        // workspace product is a 2-D DNA CONSTRUCT MAP (kind2d = Some(DnaMap)
        // with features non-empty spanning total_nt, mesh None) carrying the
        // sequence rows.
        let mut app = ValenxApp::default();
        app.wb_agent_counter = 1;
        let dir = isolate_cmd_dir(&mut app, "show2d_dna");
        let path = cmd_path(&app, 1);
        std::fs::write(&path, "{\"cmd\":\"show_2d\",\"kind\":\"dna\"}\n").unwrap();

        assert!(app.workspace_products.is_empty());
        poll_and_apply_agent_commands(&mut app);

        let product = app
            .workspace_products
            .get(&1)
            .expect("channel-1 product set by show_2d dna");
        assert!(product.mesh.is_none(), "dna 2-D carries no 3-D mesh");
        match product.kind2d.as_ref() {
            Some(crate::Workspace2dKind::DnaMap(map)) => {
                assert_eq!(map.total_nt, 93, "the canonical construct is 93 nt");
                assert!(!map.features.is_empty(), "the map has feature spans");
                // Every feature is an in-bounds half-open interval and the spans
                // collectively cover the full construct (a 0-start feature and a
                // feature ending exactly at total_nt).
                for f in &map.features {
                    assert!(f.start < f.end, "{}: start < end", f.label);
                    assert!(f.end <= map.total_nt, "{}: end within total_nt", f.label);
                }
                assert!(
                    map.features.iter().any(|f| f.start == 0),
                    "a feature starts at 0 (the ATG/ORF)"
                );
                assert!(
                    map.features.iter().any(|f| f.end == map.total_nt),
                    "a feature reaches total_nt (the stop)"
                );
                // The codon-optimisation rows came along (incl. the honesty note).
                assert!(
                    product.lines.iter().any(|l| l.contains("GC content")),
                    "map carries the GC row: {:?}",
                    product.lines
                );
                assert!(
                    product
                        .lines
                        .iter()
                        .any(|l| l.contains("NOT a biosecurity screen")),
                    "the no-biosecurity-screen note is present: {:?}",
                    product.lines
                );
            }
            other => panic!("expected DnaMap kind2d, got {other:?}"),
        }
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn show_2d_unknown_kind_is_skipped() {
        // An unknown `kind` is a safe no-op (no product inserted, no panic).
        let mut app = ValenxApp::default();
        apply(
            &mut app,
            1,
            AgentCommand::Show2d {
                n: None,
                kind: "not-a-drawing".into(),
            },
        );
        assert!(!app.workspace_products.contains_key(&1));
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

    #[test]
    fn clear_command_files_also_removes_the_global_channel_file() {
        // The global (no-`_u`-suffix) channel file `valenx_chat_cmd.jsonl` is
        // wiped at launch too, so a stale `new_unit` is never replayed.
        let dir =
            std::env::temp_dir().join(format!("valenx_clear_global_test_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();

        let global = dir.join(format!("{CMD_STEM}.jsonl"));
        let per_unit = dir.join(format!("{CMD_STEM}_u1.jsonl"));
        std::fs::write(&global, "{\"cmd\":\"new_unit\"}\n").unwrap();
        std::fs::write(&per_unit, "{\"cmd\":\"note\",\"text\":\"x\"}\n").unwrap();

        let mut app = ValenxApp::default();
        app.assistant
            .set_inbox_path_for_test(dir.join("assistant_inbox.jsonl"));
        assert_eq!(
            global_cmd_path(&app),
            global,
            "test base dir wired up correctly"
        );

        clear_command_files(&app);

        assert!(!global.exists(), "global channel file deleted");
        assert!(!per_unit.exists(), "per-unit channel file deleted");

        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn new_unit_parses_bare_and_rich_wire_forms() {
        // The bootstrap command parses from its global-channel wire form: the
        // bare `{"cmd":"new_unit"}` (all fields default to None) and the rich
        // `{"cmd":"new_unit","kind":"rocket","title":"Rocket"}`.
        let bare: AgentCommand = serde_json::from_str(r#"{"cmd":"new_unit"}"#).unwrap();
        assert_eq!(
            bare,
            AgentCommand::NewUnit {
                kind: None,
                title: None,
                note: None,
            }
        );
        let rich: AgentCommand =
            serde_json::from_str(r#"{"cmd":"new_unit","kind":"rocket","title":"Rocket"}"#).unwrap();
        assert_eq!(
            rich,
            AgentCommand::NewUnit {
                kind: Some("rocket".into()),
                title: Some("Rocket".into()),
                note: None,
            }
        );
        // The full form with a narration note parses too.
        let full: AgentCommand = serde_json::from_str(
            r#"{"cmd":"new_unit","kind":"gear","title":"Gear train","note":"Designing the reducer…"}"#,
        )
        .unwrap();
        assert_eq!(
            full,
            AgentCommand::NewUnit {
                kind: Some("gear".into()),
                title: Some("Gear train".into()),
                note: Some("Designing the reducer…".into()),
            }
        );
    }

    #[test]
    fn new_unit_increments_the_unit_counter_by_one() {
        // Applying a bare `new_unit` through the global handler opens exactly one
        // unit (the counter goes up by 1) and turns the dock on, exactly as the
        // "+ Workbench+Agent → New row at bottom" button would.
        let mut app = ValenxApp::default();
        assert_eq!(app.wb_agent_counter, 0);
        apply_global(
            &mut app,
            AgentCommand::NewUnit {
                kind: None,
                title: None,
                note: None,
            },
        );
        assert_eq!(app.wb_agent_counter, 1, "exactly one unit opened");
        assert!(app.dock_enabled, "the dock is turned on for the unit grid");
        // A bare new_unit (no kind, no title) renders no workspace product.
        assert!(!app.workspace_products.contains_key(&1));
    }

    #[test]
    fn new_unit_with_kind_renders_a_product_into_the_new_unit() {
        // `new_unit` with `kind:"rocket"` opens unit 1 AND publishes the rocket
        // product into `workspace_products[&1]` via the same show_3d path a
        // running agent uses (a live LV-1 mesh).
        let mut app = ValenxApp::default();
        apply_global(
            &mut app,
            AgentCommand::NewUnit {
                kind: Some("rocket".into()),
                title: None,
                note: None,
            },
        );
        assert_eq!(app.wb_agent_counter, 1);
        let product = app
            .workspace_products
            .get(&1)
            .expect("new_unit{kind:rocket} renders a product into the new unit");
        let mesh = product
            .mesh
            .as_ref()
            .expect("the rocket kind attaches a live mesh");
        assert!(mesh.path.to_string_lossy().contains("valenx-lv1"));
    }

    #[test]
    fn new_unit_title_overrides_the_rendered_product_heading() {
        // When both `kind` and `title` are set, the product renders and its
        // heading is replaced by the caller's title.
        let mut app = ValenxApp::default();
        apply_global(
            &mut app,
            AgentCommand::NewUnit {
                kind: Some("rocket".into()),
                title: Some("LV-1 Heavy".into()),
                note: None,
            },
        );
        let product = app.workspace_products.get(&1).expect("product set");
        assert_eq!(product.title, "LV-1 Heavy", "title overrode the heading");
        assert!(product.mesh.is_some(), "the rocket mesh still rendered");
    }

    #[test]
    fn new_unit_title_only_shows_a_named_card() {
        // A `title` with no `kind` shows a title-only card so the workspace
        // names itself even with nothing built yet.
        let mut app = ValenxApp::default();
        apply_global(
            &mut app,
            AgentCommand::NewUnit {
                kind: None,
                title: Some("Empty bench".into()),
                note: None,
            },
        );
        let product = app
            .workspace_products
            .get(&1)
            .expect("a title-only card is shown");
        assert_eq!(product.title, "Empty bench");
        assert!(product.mesh.is_none() && product.lines.is_empty());
    }

    #[test]
    fn new_unit_is_bounded_at_max_units() {
        // Once `MAX_UNITS` units exist, a further `new_unit` is refused (the
        // counter does not advance) so a runaway file can't spawn unbounded
        // panes.
        let mut app = ValenxApp::default();
        app.wb_agent_counter = MAX_UNITS;
        apply_global(
            &mut app,
            AgentCommand::NewUnit {
                kind: None,
                title: None,
                note: None,
            },
        );
        assert_eq!(
            app.wb_agent_counter, MAX_UNITS,
            "new_unit refused at the cap"
        );
    }

    #[test]
    fn new_unit_on_a_per_unit_channel_is_a_no_op() {
        // `new_unit` is global-only: routed through the per-unit `apply` it does
        // nothing (no extra unit, no product), since a per-unit channel already
        // is a unit.
        let mut app = ValenxApp::default();
        app.wb_agent_counter = 3;
        apply(
            &mut app,
            3,
            AgentCommand::NewUnit {
                kind: Some("rocket".into()),
                title: Some("x".into()),
                note: None,
            },
        );
        assert_eq!(app.wb_agent_counter, 3, "per-unit new_unit opens no unit");
        assert!(app.workspace_products.is_empty());
    }

    #[test]
    fn global_poll_opens_a_unit_and_builds_a_product_from_the_file() {
        // END-TO-END through the REAL poll path: an agent appends a rich
        // `new_unit` line to the GLOBAL command file (no unit exists yet,
        // wb_agent_counter == 0). The first poll reads the global channel,
        // opens unit 1, renders the gear product into it, applies the title,
        // and advances the global cursor — all from the file, no UI click.
        let mut app = ValenxApp::default();
        let dir = isolate_cmd_dir(&mut app, "global_newunit");
        let path = global_cmd_path(&app);
        std::fs::write(
            &path,
            "{\"cmd\":\"new_unit\",\"kind\":\"gear\",\"title\":\"Reducer\",\"note\":\"building\"}\n",
        )
        .unwrap();

        assert_eq!(app.wb_agent_counter, 0, "no unit before the poll");
        assert!(app.agent_global_cmd_cursor.is_none());
        poll_and_apply_agent_commands(&mut app);

        assert_eq!(app.wb_agent_counter, 1, "the global new_unit opened unit 1");
        let product = app
            .workspace_products
            .get(&1)
            .expect("the gear product was built into the new unit");
        // The title overrode the gear product's heading, and the gear mesh is live.
        assert_eq!(product.title, "Reducer");
        assert!(product.mesh.is_some(), "gear kind attaches a live mesh");
        // The global cursor advanced past the one applied line.
        assert_eq!(app.agent_global_cmd_cursor, Some(1));
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn global_poll_runs_only_newly_appended_new_units() {
        // The global channel mirrors the per-unit append-only semantics: the
        // first poll runs the first line (cursor 0→1), then a later append runs
        // ONLY the new line (cursor 1→2), so two units open in total — never a
        // replay of the first.
        let mut app = ValenxApp::default();
        let dir = isolate_cmd_dir(&mut app, "global_onlynew");
        let path = global_cmd_path(&app);
        std::fs::write(&path, "{\"cmd\":\"new_unit\"}\n").unwrap();

        poll_and_apply_agent_commands(&mut app);
        assert_eq!(app.wb_agent_counter, 1, "first new_unit opened one unit");
        assert_eq!(app.agent_global_cmd_cursor, Some(1));

        use std::io::Write as _;
        let mut f = std::fs::OpenOptions::new()
            .append(true)
            .open(&path)
            .unwrap();
        writeln!(f, "{{\"cmd\":\"new_unit\"}}").unwrap();
        app.last_agent_poll = None; // bypass the 1s throttle for the test
        poll_and_apply_agent_commands(&mut app);

        assert_eq!(
            app.wb_agent_counter, 2,
            "only the newly-appended new_unit ran the 2nd time (no replay)"
        );
        assert_eq!(app.agent_global_cmd_cursor, Some(2));
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// Build a `WorkspaceProduct` carrying an animation that starts paused at
    /// `1.0×`, with a simple turntable motion — the fixture the `animate` tests
    /// drive Play/Pause + speed against.
    fn animated_product() -> crate::WorkspaceProduct {
        crate::WorkspaceProduct {
            title: "Spinner".into(),
            lines: Vec::new(),
            mesh: None,
            vertex_colors: None,
            camera: valenx_viz::OrbitCamera::default(),
            kind2d: None,
            last_export: None,
            image: None,
            image_texture: None,
            animation: Some(crate::ProductAnimation {
                playing: false,
                speed: 1.0,
                t: 0.0,
                motion: crate::ProductMotion::Turntable {
                    axis: [0.0, 1.0, 0.0],
                    pivot: [0.0, 0.0, 0.0],
                    rad_per_s: 1.0,
                },
            }),
        }
    }

    #[test]
    fn animate_round_trips_its_wire_form() {
        // The `animate` command parses from its `{"cmd":"animate",…}` wire form;
        // `n` / `play` / `speed` are each optional.
        let full: AgentCommand =
            serde_json::from_str(r#"{"cmd":"animate","n":1,"play":true,"speed":2.0}"#).unwrap();
        assert_eq!(
            full,
            AgentCommand::Animate {
                n: Some(1),
                play: Some(true),
                speed: Some(2.0),
            }
        );
        // A bare `animate` (drive the command file's own channel, toggle
        // nothing) parses with every field None.
        let bare: AgentCommand = serde_json::from_str(r#"{"cmd":"animate"}"#).unwrap();
        assert_eq!(
            bare,
            AgentCommand::Animate {
                n: None,
                play: None,
                speed: None,
            }
        );
    }

    #[test]
    fn animate_flips_playing_and_speed_on_an_animated_product() {
        // `animate{n:1,play:true,speed:2.0}` applied to a unit whose product has
        // an animation flips it to playing at 2.0×, through the same state a
        // user's Play button + speed slider would set.
        let mut app = ValenxApp::default();
        app.workspace_products.insert(1, animated_product());
        // Routed by the explicit `n` even though we hand `apply` a different
        // channel, so it works arriving on the global channel too.
        apply(
            &mut app,
            7,
            AgentCommand::Animate {
                n: Some(1),
                play: Some(true),
                speed: Some(2.0),
            },
        );
        let anim = app.workspace_products[&1]
            .animation
            .as_ref()
            .expect("the product keeps its animation");
        assert!(anim.playing, "play:true started the clock");
        assert_eq!(anim.speed, 2.0, "speed:2.0 was applied");
    }

    #[test]
    fn animate_clamps_speed_and_defaults_the_target_to_the_channel() {
        // `speed` is clamped to the toolbar's 0.0..=8.0, and a missing `n`
        // targets the command file's channel `ch` (here 1) — pausing it.
        let mut app = ValenxApp::default();
        app.workspace_products.insert(1, animated_product());
        // First start it, with an out-of-range speed that must clamp to 8.0.
        apply(
            &mut app,
            1,
            AgentCommand::Animate {
                n: None,
                play: Some(true),
                speed: Some(99.0),
            },
        );
        {
            let anim = app.workspace_products[&1].animation.as_ref().unwrap();
            assert!(anim.playing, "play:true via the channel default");
            assert_eq!(anim.speed, 8.0, "speed clamped to the 8.0 ceiling");
        }
        // Now pause it (speed left untouched → stays 8.0).
        apply(
            &mut app,
            1,
            AgentCommand::Animate {
                n: None,
                play: Some(false),
                speed: None,
            },
        );
        let anim = app.workspace_products[&1].animation.as_ref().unwrap();
        assert!(!anim.playing, "play:false paused it");
        assert_eq!(anim.speed, 8.0, "omitted speed left the value untouched");
    }

    #[test]
    fn animate_on_a_product_without_an_animation_is_a_safe_no_op() {
        // A unit whose product has no animation just gets a feed note (handled
        // inside `apply`); nothing panics and no animation appears.
        let mut app = ValenxApp::default();
        app.workspace_products.insert(
            2,
            crate::WorkspaceProduct {
                title: "Static".into(),
                lines: Vec::new(),
                mesh: None,
                vertex_colors: None,
                camera: valenx_viz::OrbitCamera::default(),
                kind2d: None,
                last_export: None,
                image: None,
                image_texture: None,
                animation: None,
            },
        );
        apply(
            &mut app,
            2,
            AgentCommand::Animate {
                n: None,
                play: Some(true),
                speed: Some(2.0),
            },
        );
        assert!(
            app.workspace_products[&2].animation.is_none(),
            "no animation was conjured onto a static product"
        );
        // An animate aimed at a unit with no product at all is also a no-op.
        apply(
            &mut app,
            99,
            AgentCommand::Animate {
                n: None,
                play: Some(true),
                speed: None,
            },
        );
        assert!(!app.workspace_products.contains_key(&99));
    }

    /// Collect the pane ids in a dock tile tree (mirrors `dock_layout`'s private
    /// `pane_ids` test helper) so a test can assert a unit's `workspace:<n>`
    /// pane is present.
    fn dock_pane_ids(tree: &egui_tiles::Tree<String>) -> std::collections::HashSet<String> {
        tree.tiles
            .tiles()
            .filter_map(|t| match t {
                egui_tiles::Tile::Pane(id) => Some(id.clone()),
                _ => None,
            })
            .collect()
    }

    #[test]
    fn global_new_unit_opens_its_own_tab_with_the_units_dock() {
        // One-tab-per-product: a global `new_unit` parks the current tab, opens a
        // fresh project tab, and builds the unit's `[workspace:n | agent:n]` pair
        // into THAT tab's (now live) dock tree. Drive it end-to-end through the
        // REAL poll path and assert the tab strip grew by one and the installed
        // (live) dock contains the new unit's `workspace:1` pane.
        let mut app = ValenxApp::default();
        let dir = isolate_cmd_dir(&mut app, "global_newunit_tab");
        let path = global_cmd_path(&app);
        std::fs::write(&path, "{\"cmd\":\"new_unit\",\"title\":\"Rocket bay\"}\n").unwrap();

        let tabs_before = app.tab_bar.tabs.len();
        assert_eq!(app.wb_agent_counter, 0, "no unit before the poll");
        poll_and_apply_agent_commands(&mut app);

        // Exactly one new tab, made active, titled from the command.
        assert_eq!(
            app.tab_bar.tabs.len(),
            tabs_before + 1,
            "new_unit grew the tab strip by one"
        );
        let active = app.tab_bar.active.expect("the new tab is active");
        assert_eq!(app.tab_bar.tabs[active].title, "Rocket bay");
        // The unit opened (counter bumped) and its dock is checked out into the
        // live `dock_tree` (this tab is active), carrying its workspace pane.
        assert_eq!(app.wb_agent_counter, 1);
        let tree = app
            .dock_tree
            .as_ref()
            .expect("the new tab's unit dock is installed live");
        assert!(
            dock_pane_ids(tree).contains("workspace:1"),
            "the installed dock holds the unit's workspace:1 pane"
        );
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn global_new_unit_without_title_or_kind_names_the_tab_by_unit_number() {
        // With neither `title` nor `kind`, the tab title is a placeholder that is
        // rewritten to "Unit <n>" once the unit number is known.
        let mut app = ValenxApp::default();
        let dir = isolate_cmd_dir(&mut app, "global_newunit_placeholder");
        let path = global_cmd_path(&app);
        std::fs::write(&path, "{\"cmd\":\"new_unit\"}\n").unwrap();

        poll_and_apply_agent_commands(&mut app);

        let active = app.tab_bar.active.expect("the new tab is active");
        assert_eq!(
            app.tab_bar.tabs[active].title, "Unit 1",
            "a placeholder tab names itself by unit number"
        );
        let _ = std::fs::remove_dir_all(&dir);
    }
}
