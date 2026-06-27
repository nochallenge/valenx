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

use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};
use std::path::{Path, PathBuf};
use std::time::{Duration, Instant};

use serde::Deserialize;

use valenx_viz::ViewDirection;

use crate::project_tabs::{self, TabKind};
use crate::ValenxApp;

/// A typed value an agent can assign to a labelled workbench control via
/// [`SetControl`](AgentCommand::SetControl). **Untagged** on the wire, so the
/// JSON literal is written directly — `42`, `4000`, `0.55`, `true`,
/// `"linear"` — and serde picks the first matching arm. Order matters: `bool`
/// is tried before the numbers (a JSON `true` is *only* a bool), then the
/// integer arm (so a whole number like `4000` arrives as `Int`, not a lossy
/// float), then the float arm, then the string fallback.
///
/// The receiving workbench's `agent_set` decides how to interpret the value for
/// a given control (e.g. a `usize` sample-count reads [`as_i64`](AgentValue::as_i64),
/// a `f64` coefficient reads [`as_f64`](AgentValue::as_f64), an enum-by-name
/// reads [`as_str`](AgentValue::as_str)); a value of the wrong shape for the
/// named control yields a fail-loud `Err(String)` (→ a `warn` feed note), never
/// a panic.
#[derive(Debug, Clone, PartialEq, Deserialize)]
#[serde(untagged)]
pub enum AgentValue {
    /// A boolean (`true` / `false`) — for checkbox / toggle controls.
    Bool(bool),
    /// A whole number — for integer controls (`usize` / `u32` / `u64` counts,
    /// sample sizes, seeds). Also coercible to `f64` for a float control.
    Int(i64),
    /// A real number — for floating-point controls.
    Float(f64),
    /// A string — for enum-by-name controls (a combo selection like
    /// `"linear"`), or a free-text field.
    Str(String),
}

impl AgentValue {
    /// Interpret this value as an `f64` for a floating-point control. NaN/Inf is
    /// rejected (fail-loud, mirroring `as_i64`) so a hostile value can never be
    /// written to a control and reported `Ok` — this single gate closes the
    /// write-then-Ok gap across every bare `as_f64()?` setter. An `Int` widens
    /// losslessly; a `Bool` / `Str` is a type error.
    pub fn as_f64(&self) -> Result<f64, String> {
        match self {
            AgentValue::Float(v) if v.is_finite() => Ok(*v),
            AgentValue::Float(v) => Err(format!("expected a finite number, got {v}")),
            AgentValue::Int(v) => Ok(*v as f64),
            other => Err(format!("expected a number, got {other:?}")),
        }
    }

    /// Interpret this value as an `i64` for an integer control. A whole-valued
    /// `Float` (e.g. `4000.0`) is accepted and truncated to the integer; a
    /// fractional `Float`, a `Bool`, or a `Str` is a type error (fail-loud).
    pub fn as_i64(&self) -> Result<i64, String> {
        match self {
            AgentValue::Int(v) => Ok(*v),
            // Accept a float that is exactly integral so `{"value": 4000.0}`
            // still sets a usize control; reject a fractional float loudly.
            AgentValue::Float(v) if v.fract() == 0.0 && v.is_finite() => Ok(*v as i64),
            other => Err(format!("expected an integer, got {other:?}")),
        }
    }

    /// Interpret this value as a `bool` for a toggle control. Only a JSON
    /// boolean qualifies; a number / string is a type error (fail-loud) so a
    /// stray `1` never silently flips a flag.
    pub fn as_bool(&self) -> Result<bool, String> {
        match self {
            AgentValue::Bool(v) => Ok(*v),
            other => Err(format!("expected a boolean, got {other:?}")),
        }
    }

    /// Interpret this value as a string slice for an enum-by-name / text
    /// control. Only a JSON string qualifies (a number is a type error) so an
    /// enum selection is always written explicitly as text.
    pub fn as_str(&self) -> Result<&str, String> {
        match self {
            AgentValue::Str(s) => Ok(s),
            other => Err(format!("expected a string, got {other:?}")),
        }
    }
}

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

/// Pick a **deterministic** header colour for a category group from its name,
/// used by [`apply_global`] when a [`NewUnit`](AgentCommand::NewUnit) `group`
/// mints a fresh band.
///
/// Same name → same colour (so a category reads identically every time it is
/// re-created in a session, and two units that mint the *same* group would get
/// the same colour even if they raced), while distinct names spread across the
/// hue wheel so categories are visually separable. Implementation: hash the
/// name with [`DefaultHasher`], take `hue = hash % 360`, and convert from HSV at
/// a fixed **mid** saturation `0.55` and value `0.85` — bright, pleasant,
/// mid-saturation tints that are never near-white (s ≠ 0) nor near-black
/// (v well above 0). The fixed S/V is what keeps every group legible against
/// the strip regardless of which hue the hash lands on.
fn color_for(name: &str) -> [u8; 3] {
    let mut hasher = DefaultHasher::new();
    name.hash(&mut hasher);
    let hue = (hasher.finish() % 360) as f32; // 0..360 degrees around the wheel
    hsv_to_rgb8(hue, 0.55, 0.85)
}

/// Convert an HSV colour (`hue` in degrees `0..360`, `s`/`v` in `0.0..=1.0`) to
/// an 8-bit RGB triple. A small dependency-free helper for [`color_for`]; the
/// standard piece-wise HSV→RGB conversion.
fn hsv_to_rgb8(hue: f32, s: f32, v: f32) -> [u8; 3] {
    let h = hue.rem_euclid(360.0) / 60.0; // sector 0..6
    let c = v * s; // chroma
    let x = c * (1.0 - (h.rem_euclid(2.0) - 1.0).abs());
    let m = v - c;
    let (r1, g1, b1) = match h as u32 {
        0 => (c, x, 0.0),
        1 => (x, c, 0.0),
        2 => (0.0, c, x),
        3 => (0.0, x, c),
        4 => (x, 0.0, c),
        _ => (c, 0.0, x), // 5 (and the h==6.0 wrap)
    };
    let to_u8 = |f: f32| ((f + m) * 255.0).round().clamp(0.0, 255.0) as u8;
    [to_u8(r1), to_u8(g1), to_u8(b1)]
}

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
        /// Optional **category group** to file the new unit's tab into — the
        /// named, coloured, collapsible Chrome-style band in the tab strip (see
        /// [`crate::project_tabs::TabGroup`]). When set, `apply_global` files
        /// the just-opened tab into an existing group of that name (so two
        /// `new_unit`s sharing a `group` land in ONE band), or mints a new band
        /// named `group` with a deterministic colour (`color_for`) when none
        /// exists yet. This is what lets an agent organise ~130 product tabs
        /// into a handful of manageable categories instead of one flat strip.
        /// Absent → the tab is left ungrouped (the prior behaviour), so older
        /// `new_unit` JSON without `group` still parses (serde default).
        #[serde(default)]
        group: Option<String>,
    },
    /// **Drive an animated product's playback clock** on a unit whose
    /// `workspace:<n>` product carries a [`crate::ProductAnimation`] (e.g. the
    /// meshing gear train) — the file-driven equivalent of clicking the tile's
    /// Play/Pause button and dragging its speed slider. `play` sets the
    /// playing/paused state and `speed` sets the playback multiplier (clamped to
    /// the toolbar's `0.0..=4.0`, with a non-finite value falling back to `1.0`);
    /// either may be omitted to leave that field
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
    /// **Run a command-palette action by its stable id.** Bridges the file
    /// channel into the *existing* command-palette registry
    /// ([`crate::commands::static_commands`]): `id` is matched against a
    /// [`crate::commands::Command::id`] (the stable `&'static str` like
    /// `"view.front"`, `"run.selected-case"`, `"settings.open"`) and the
    /// matching command is invoked through the **same** `(cmd.invoke)(app)`
    /// function pointer a user click / `Ctrl+P` selection runs — so the ~66
    /// palette actions become drivable over the robust polled file-bridge with
    /// **no duplicated action logic**. An unknown `id` posts an "unknown command
    /// id" feed note and is otherwise a no-op (never a panic); a successful run
    /// posts a `"ran <id>"` ack note.
    ///
    /// Honest scope: most ids are pure in-process state (the 8 camera views,
    /// shading, run/sweep/cancel, settings, audit-tail) and drive cleanly
    /// headless. A handful open a **native file dialog** (`file.new-project`,
    /// `file.open-project`, `file.import-stl`, `file.load-mesh`,
    /// `file.save-mesh-stl`, the HTML/CSV/audit-open file-browser actions) —
    /// they are still exposed (a user driving the GUI may want them) but are not
    /// usefully driven headless. The bridge does **not** try to suppress those
    /// dialogs.
    RunCommand {
        /// The stable command id to run (e.g. `"view.front"`).
        id: String,
    },
    /// **Enumerate the available command-palette ids** into this channel's chat
    /// feed, so an agent can *discover* what [`RunCommand`](AgentCommand::RunCommand)
    /// accepts without hard-coding the list. Posts a single feed note listing
    /// every [`crate::commands::static_commands`] id (the same registry
    /// `RunCommand` resolves against). No app state changes.
    ListCommands,
    /// **Set a labelled workbench parameter by its user-visible caption.** This
    /// closes the param-setting half of AI-drivability: where
    /// [`RunCommand`](AgentCommand::RunCommand) fires *actions*, `SetControl`
    /// writes a tool's *input values* over the robust polled file-bridge.
    ///
    /// `name` is the exact caption the user sees next to the control (e.g.
    /// `"Monte-Carlo samples N"`, `"a1 (coeff on x1)"`, `"sensor range (m)"`) —
    /// the **same** string the widget is `labelled_by`, so an agent that read a
    /// caption from the accessibility tree (or from
    /// [`ListControls`](AgentCommand::ListControls)) can set it by that name.
    /// `value` is the typed [`AgentValue`] to assign.
    ///
    /// `workbench` selects the target tool by id (a [`TabKind::from_id`] alias,
    /// e.g. `"uq"`); when omitted the **active tab's** workbench is used (so an
    /// agent driving the visible tool needn't name it). Resolution + the
    /// per-workbench validated assignment go through `set_control` → that
    /// workbench's `agent_set`; an unknown workbench, an unknown caption, or a
    /// value of the wrong type posts a `warn` feed note and changes nothing
    /// (never a panic). A successful set posts an ack note naming the control
    /// and its new value, so the change is visible in the unit's chat.
    SetControl {
        /// The control's user-visible caption (== its `labelled_by` text).
        name: String,
        /// The typed value to assign (untagged: `42`, `0.55`, `true`, `"sum"`).
        value: AgentValue,
        /// Optional target workbench id (default: the active tab's workbench).
        #[serde(default)]
        workbench: Option<String>,
    },
    /// **Enumerate the settable control captions** of a workbench into this
    /// channel's chat feed, so an agent can *discover* what
    /// [`SetControl`](AgentCommand::SetControl) accepts without hard-coding the
    /// names. `workbench` selects the tool by id (default: the active tab's
    /// workbench). Posts a single feed note listing every caption that
    /// workbench's `agent_set` recognises; a workbench with no `SetControl`
    /// support yet posts a note saying so. No app state changes.
    ListControls {
        /// Optional target workbench id (default: the active tab's workbench).
        #[serde(default)]
        workbench: Option<String>,
    },
    /// **Read a workbench's COMPUTED result back into this channel's chat feed.**
    /// This closes the live-driving loop: where [`SetControl`](AgentCommand::SetControl)
    /// writes a tool's inputs and [`RunCommand`](AgentCommand::RunCommand) fires
    /// the solve, `ReadReadout` reads the *answer* back so an agent can
    /// **self-verify** what it just drove — no screenshot, no GUI focus needed.
    ///
    /// `workbench` selects the target tool by id (a [`TabKind::from_id`] alias,
    /// e.g. `"uq"`); when omitted the **active tab's** workbench is used. The
    /// resolved workbench's read-only `agent_readout` returns its current result
    /// text (the same string the panel renders) — posted as a `result` feed note.
    /// A workbench that has not been run yet (or has no readout wired) posts a
    /// `warn` "not run yet / no readout" note instead; an unknown/unresolved
    /// workbench posts a `warn` and changes nothing (never a panic). Purely
    /// read-only: no app state is modified.
    ReadReadout {
        /// Optional target workbench id (default: the active tab's workbench).
        #[serde(default)]
        workbench: Option<String>,
    },
    /// **Click any button / widget in the active workbench by its accessible
    /// name** — the generic, no-per-workbench-wiring drive command. Where
    /// [`SetControl`](AgentCommand::SetControl) writes a tool's *inputs* and
    /// [`RunCommand`](AgentCommand::RunCommand) fires a *registered* action,
    /// `invoke_named` triggers the workbench's own in-panel buttons (its
    /// "▶ Compute" / "Run" / "Apply" / "Analyze") **by the exact caption a
    /// screen reader / UI-Automation client sees** — closing the AI-drivability
    /// gap for the ~40 workbenches that have named controls + buttons but no
    /// bespoke bridge run-id.
    ///
    /// Mechanism (identical to how an external UIA client clicks by name):
    /// `invoke_named` runs a **headless probe** of the active workbench's panel
    /// in a throwaway accesskit-enabled context, finds the node whose accessible
    /// `name` matches `name`, and queues an `accesskit` `Default` action on that
    /// node id ([`crate::ValenxApp::pending_accesskit_actions`]). On the **next**
    /// frame `crate::ValenxApp::raw_input_hook` injects it as an
    /// `AccessKitActionRequest`, and egui reports the matching button as
    /// `.clicked()` — the same path a real click takes (the bridge keeps frames
    /// alive via the unfocused `request_repaint`, so the queued action fires
    /// promptly). The node id is a deterministic hash of the widget's egui id, so
    /// the id resolved in the probe is the *same* id the live frame uses.
    ///
    /// Matching is exact, then a case-insensitive fallback. No match (or no
    /// active workbench) posts a `warn` feed note and changes nothing (never a
    /// panic); a match posts an ack naming the button it queued.
    InvokeNamed {
        /// The button/widget's accessible name (== the caption the user sees,
        /// the same string the UI-Automation tree exposes as `Name`).
        name: String,
    },
    /// **Enumerate the clickable button captions in the active workbench** into
    /// this channel's chat feed, so an agent can *discover* what
    /// [`InvokeNamed`](AgentCommand::InvokeNamed) accepts without hard-coding
    /// them. Runs the same headless accesskit probe and lists every node whose
    /// role is a button (i.e. invokable), by name. No app state changes; an
    /// empty/closed workbench posts a note saying so.
    ListButtons,
    /// **Read the active workbench panel's visible text back into this channel's
    /// chat feed**, generically — the accessibility-tree counterpart to
    /// [`ReadReadout`](AgentCommand::ReadReadout) that needs no per-workbench
    /// `agent_readout`. Runs the same headless accesskit probe and concatenates
    /// the readable text nodes (labels / values) of the active workbench, so an
    /// agent can self-verify a computed result by name even on a tool that has
    /// no bespoke readout wired. Posted as a `result` feed note (bounded in
    /// length). Purely read-only; a closed workbench posts a `warn`.
    ReadText,
    /// **Snap the central 3-D viewport camera to a canonical view.** `dir` is one
    /// of `front` / `back` / `left` / `right` / `top` / `bottom` / `iso`
    /// (case-insensitive), mapped to a [`valenx_viz::ViewDirection`] and applied
    /// via the **same** `OrbitCamera::set_view` the ViewCube buttons drive
    /// ([`crate::ValenxApp::camera_mut`]). This closes the camera half of the
    /// viewport AI-drivability gap over the robust polled file-bridge. An
    /// unrecognised `dir` posts a `warn` feed note and changes nothing (never a
    /// panic); a successful snap posts an ack note. The camera *target*/distance
    /// are untouched — only the orbit angles snap, exactly like the ViewCube.
    SetView {
        /// Canonical view name (case-insensitive): `front`, `back`, `left`,
        /// `right`, `top`, `bottom`, `iso`.
        dir: String,
    },
    /// **Orbit the central 3-D camera by a degree delta** — the file-driven
    /// equivalent of a middle-mouse drag. `dx_deg` changes azimuth, `dy_deg`
    /// changes elevation (clamped to `±89.9°` by `OrbitCamera::orbit`, the same
    /// vetted method the drag uses). Always succeeds (any finite delta is valid;
    /// a non-finite delta is ignored with a `warn` note); posts an ack note with
    /// the new azimuth/elevation.
    Orbit {
        /// Azimuth delta in degrees (horizontal orbit).
        dx_deg: f32,
        /// Elevation delta in degrees (vertical orbit, clamped to ±89.9°).
        dy_deg: f32,
    },
    /// **Dolly the central 3-D camera in/out** — the file-driven equivalent of a
    /// scroll-wheel zoom. `factor` is the fractional zoom `OrbitCamera::zoom`
    /// takes: positive zooms **in** (e.g. `0.1` = 10% closer), negative zooms
    /// out; the method clamps so the camera can't invert through the target. A
    /// non-finite `factor` is ignored with a `warn` note; a valid zoom posts an
    /// ack note with the new distance.
    Zoom {
        /// Fractional zoom amount (positive = in, negative = out).
        factor: f32,
    },
    /// **Frame the whole loaded model in the central 3-D viewport** — the
    /// file-driven equivalent of the "Fit / Frame all" action. Reframes the
    /// camera around the loaded mesh's (or STL's) axis-aligned bounding box via
    /// the **same** [`crate::ValenxApp::frame_current_mesh`] /
    /// [`crate::ValenxApp::frame_current_stl`] methods the menu uses. When
    /// **nothing** is loaded the camera is left unchanged and a `warn` feed note
    /// says so (never a panic); on success an ack note is posted.
    FrameAll,
    /// **Auto-tile / organize the open workspace panels into a balanced grid**
    /// — the file-driven equivalent of the tab-strip "Tile" button. Routes
    /// through the **same** `crate::dock_layout::ValenxApp::auto_tile_dock`
    /// the button calls: every open central-workspace panel (workbenches, the
    /// Assistant, any Workbench+Agent tiles) is reflowed into a
    /// `ceil(sqrt(N))`-column grid so all of them stay visible and legible
    /// instead of one being crushed to a sliver. A no-op (warn note) when no
    /// panel is open or on a clean agent product tab; on success an ack note
    /// reports how many panels were tiled.
    Tile,
    /// **Add a straight-line vertex to the in-house CAD sketch** by explicit
    /// model-space coordinates — the file-driven equivalent of a Line-tool click
    /// on the sketch canvas. Routes through the **same**
    /// [`crate::cad_workbench::CadWorkbenchState::sketch_add_point`] the mouse
    /// click uses (first point seeds the start anchor, each later point appends a
    /// `Line` segment; a no-op once the loop is closed). Always parses; posts an
    /// ack note with the new anchor count.
    AddSketchPoint {
        /// Sketch-plane X (model units).
        x: f64,
        /// Sketch-plane Y (model units).
        y: f64,
    },
    /// **Add a 3-point circular arc to the in-house CAD sketch** — the
    /// file-driven equivalent of the Arc tool's three clicks (`start`, then a
    /// point `via` ON the arc, then `end`). Routes through the **same**
    /// [`crate::cad_workbench::CadWorkbenchState::sketch_add_arc`] the canvas
    /// uses (if the sketch is empty `start` seeds the start anchor; a no-op once
    /// the loop is closed). Always parses; posts an ack note.
    AddSketchArc {
        /// Arc start point `[x, y]` (model units).
        start: [f64; 2],
        /// A point ON the arc between start and end `[x, y]` (model units).
        via: [f64; 2],
        /// Arc end point `[x, y]` (model units).
        end: [f64; 2],
    },
    /// **Extrude the current in-house CAD sketch profile into a solid** — the
    /// file-driven equivalent of the panel's "Extrude sketch" button. Sweeps the
    /// drawn profile (from
    /// [`crate::cad_workbench::CadWorkbenchState::sketch_points`]) along +Z by
    /// `height` through the **same**
    /// [`crate::cad_workbench::CadWorkbenchState::add_extrude_from_sketch`] +
    /// `request_rebuild` path the button uses. `height` must be `> 0` (and
    /// finite) — otherwise a `warn` feed note is posted and nothing changes
    /// (never a panic). A profile of fewer than 3 anchors is also reported. On a
    /// valid extrude an ack note is posted and the tree is flagged to rebuild
    /// into the viewport on the next frame.
    ExtrudeSketch {
        /// Extrude depth along +Z (model units); must be `> 0`.
        height: f64,
    },
    /// **Add a straight line to the in-house 2-D drafting drawing** by explicit
    /// endpoints — the file-driven equivalent of the 2-D form's line `+` button.
    /// Routes through the **same** `Drawing2D::add(Entity2D::Line { .. })` path
    /// (via [`crate::draft2d_workbench::Draft2dWorkbenchState::agent_add_line`],
    /// layer `"0"`). Always parses; posts an ack note with the new entity count.
    ///
    /// Note the explicit `rename`: the enum's `rename_all = "snake_case"` would
    /// map `Add2dLine` to `"add2d_line"` (no underscore before the digit), but
    /// the wire tag is pinned to `"add_2d_line"` (matching the `show_2d` style).
    #[serde(rename = "add_2d_line")]
    Add2dLine {
        /// First endpoint X (drawing units).
        x1: f64,
        /// First endpoint Y (drawing units).
        y1: f64,
        /// Second endpoint X (drawing units).
        x2: f64,
        /// Second endpoint Y (drawing units).
        y2: f64,
    },
    /// **Add a circle to the in-house 2-D drafting drawing** by centre + radius —
    /// the file-driven equivalent of the 2-D form's circle `+` button. Routes
    /// through the **same** `Drawing2D::add(Entity2D::Circle { .. })` path (via
    /// [`crate::draft2d_workbench::Draft2dWorkbenchState::agent_add_circle`],
    /// layer `"0"`). `r` must be `> 0` (and finite) — otherwise a `warn` feed
    /// note is posted and nothing changes (never a panic). On success an ack note
    /// with the new entity count is posted.
    ///
    /// Explicit `rename` for the same reason as [`Add2dLine`](AgentCommand::Add2dLine):
    /// `rename_all` would yield `"add2d_circle"`, but the wire tag is pinned to
    /// `"add_2d_circle"`.
    #[serde(rename = "add_2d_circle")]
    Add2dCircle {
        /// Circle centre X (drawing units).
        cx: f64,
        /// Circle centre Y (drawing units).
        cy: f64,
        /// Circle radius (drawing units); must be `> 0`.
        r: f64,
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

/// Is the agent file-bridge **active** this session? True when an external
/// agent is (or is about to be) driving valenx, used by the eframe `update()`
/// loop to decide whether to schedule the faster ~6 fps heartbeat that keeps
/// queued commands flowing while the window is unfocused/idle (see `update.rs`).
///
/// Two independent signals, either of which means "drive me":
/// - the chat-bridge **env channels** are configured
///   (`$VALENX_ASSISTANT_INBOX` / `$VALENX_ASSISTANT_FEED`) — set by the launcher
///   whenever the in-app Assistant bridge is wired up; or
/// - a **global command file** ([`global_cmd_path`]) already exists on disk —
///   the moment an agent appends its first command (even with no env override
///   and no unit yet), so a cold agent that just writes the global file is
///   served promptly.
///
/// Cheap: an env lookup plus one `Path::exists` (a single `stat`), called once
/// per frame. When neither holds the normal interactive build pays nothing
/// beyond that.
pub fn bridge_active(app: &ValenxApp) -> bool {
    std::env::var_os("VALENX_ASSISTANT_INBOX").is_some()
        || std::env::var_os("VALENX_ASSISTANT_FEED").is_some()
        || global_cmd_path(app).exists()
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
/// (`<base-dir>/valenx_chat_cmd.jsonl`).
///
/// [`NewUnit`](AgentCommand::NewUnit) is the global-channel **bootstrap** — it
/// opens a fresh Workbench+Agent unit (the only place that can, since a per-unit
/// channel already *is* a unit) and optionally builds a product into it. It is
/// handled inline below.
///
/// **Every other drive command is ALSO honoured here**, so an external agent can
/// drive the whole app from the one global file with no unit-bootstrap: the
/// app/tab-level commands (`new_tab`, `open_workbench`, `focus_tab`,
/// `rename_tab`, `close_tab`), the control/readout commands (`set_control`,
/// `run_command`, `read_readout`, `list_controls`, `list_commands`), the camera
/// /sketch/2-D commands, `note`, and `animate` (which carries its own target
/// `n`) are dispatched through the **same** [`apply`] reducer the per-unit
/// channel uses, against the **active tab / app** — with the special **channel
/// `0`** sentinel so their ack/`warn` feed notes land in the **GLOBAL** feed
/// (`valenx_chat_feed.jsonl`, see
/// [`crate::assistant_workbench::append_feed_note`]) that an agent reading the
/// global channel watches, not a per-unit feed.
///
/// The inherently per-unit *product* commands (`show_product`, `show_3d`,
/// `show_2d`) target a `workspace:<n>` pane that the global channel has no unit
/// for, so they are a deliberate no-op here (use `new_unit` to mint a unit
/// first, then drive them on that unit's channel). A malformed line never
/// reaches here (the poll loop skips unparseable lines). Like [`apply`], every
/// branch is a no-op-on-bad-input rather than a panic.
fn apply_global(app: &mut ValenxApp, cmd: AgentCommand) {
    let AgentCommand::NewUnit {
        kind,
        title,
        note,
        group,
    } = cmd
    else {
        match cmd {
            // Inherently per-unit product renders: no `workspace:<n>` pane on
            // the global channel, so skip (mint a unit with `new_unit` first).
            AgentCommand::ShowProduct { .. }
            | AgentCommand::Show3d { .. }
            | AgentCommand::Show2d { .. } => {}
            // `new_unit` is the early-bound arm above; it can't reach here.
            AgentCommand::NewUnit { .. } => {}
            // Everything else (tab ops, set/run/read/list, camera, sketch, 2-D,
            // note, animate) drives the active tab / app through the SAME
            // `apply` reducer, with channel `0` so acks go to the GLOBAL feed.
            other => apply(app, 0, other),
        }
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

    // CLEAN PRODUCT-TAB DOCK: `add_workbench_agent_pair_at` *appended* the
    // `[workspace:n | agent:n]` pair to this tab's dock tree, but the tab's
    // per-frame `sync_tree` would otherwise inject the global Assistant pane
    // (`valenx_assistant_panel`, on by default) beside the unit's own `agent:n`
    // chat — leaving the product tab with TWO chat panes. Replace the tree with
    // just this unit's pair and latch `dock_agent_only` so the sync never adds
    // the Assistant (or any other `DOCKABLE_PANELS`) into a product tab. The
    // landing tab and manually-opened Workbench+Agent units are untouched (their
    // `dock_agent_only` stays `false`).
    app.set_clean_workbench_agent_dock(n);

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

    // PER-TAB WORKBENCH LINK: record the product `kind` on this (now-active)
    // tab so [`project_tabs::sync_active`] re-shows exactly its
    // `show_<kind>_workbench` panel whenever the tab is active — the
    // inputs/calculations/readouts render on the right, alongside this unit's
    // dock (workspace render + agent chat) in the centre. Set BEFORE `kind` is
    // moved into `pending_products` below. A `kind`-less unit links nothing
    // (the dock fills the tab as before).
    if let Some(kind) = kind.as_deref() {
        let kind = kind.trim();
        if !kind.is_empty() {
            if let Some(idx) = app.tab_bar.active {
                if let Some(tab) = app.tab_bar.tabs.get_mut(idx) {
                    tab.workbench_kind = Some(kind.to_string());
                }
            }
            // Reconcile now so the linked panel is visible from this frame on
            // (not only after the next tab switch). `sync_active` clears every
            // panel then turns on just this tab's linked workbench.
            project_tabs::sync_active(app);
        }
    }

    // CATEGORY GROUP: file the just-opened unit tab into a named, coloured,
    // collapsible tab-strip band so an agent can organise ~130 product tabs into
    // a handful of categories. The unit is fully formed here (tab opened +
    // titled, dock cleaned, workbench linked), and the new tab is the **active**
    // one — `TabBar::open` set `active = Some(len-1)` and nothing since
    // (`set_clean_workbench_agent_dock`, the title rewrite, `sync_active`)
    // changes the active index, so `app.tab_bar.active` IS this unit's tab. All
    // grouping goes through the SAME vetted `TabBar` group methods a user's
    // strip context-menu drives (`new_group_with_tab` / `assign_to_group` /
    // `rename_group` / `set_group_color`) — no raw field poke. An empty `group`
    // name leaves the tab ungrouped (the prior behaviour).
    if let Some(group_name) = group.as_deref() {
        let group_name = group_name.trim();
        if !group_name.is_empty() {
            if let Some(idx) = app.tab_bar.active {
                // Reuse an existing band of this name so two `new_unit`s sharing
                // a `group` land in ONE band (no duplicate); else mint a fresh
                // band, then name + colour it deterministically so the same
                // category always reads the same colour across a session.
                let existing = app
                    .tab_bar
                    .groups
                    .iter()
                    .find(|g| g.name == group_name)
                    .map(|g| g.id.clone());
                match existing {
                    Some(gid) => app.tab_bar.assign_to_group(idx, &gid),
                    None => {
                        if let Some(gid) = app.tab_bar.new_group_with_tab(idx) {
                            app.tab_bar.rename_group(&gid, group_name);
                            app.tab_bar.set_group_color(&gid, color_for(group_name));
                        }
                    }
                }
            }
        }
    }

    // LAZY-BUILD: if a product kind was named, DON'T build it now — record it in
    // `pending_products` so the actual 3-D/2-D product is built only when this
    // unit's `workspace:<n>` pane is first VIEWED (see `materialize_pending`,
    // called from `render_workspace_body`) or when the unit is animated. Opening
    // the tab stays instant, so an agent fleet can `new_unit` 130+ tabs in a
    // burst without building every product at once (which briefly hung the app).
    // A `kind`-less `new_unit` records nothing — there is nothing to build.
    if let Some(kind) = kind {
        app.pending_products.insert(n, kind);
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
            // LAZY-BUILD: if this unit's product was only queued by `new_unit`
            // (its tab not yet viewed, so nothing in `workspace_products`),
            // build it now so the agent can animate it without first opening the
            // pane. A no-op once built / nothing pending; after it the product
            // (with its default inspect-spin) is present to toggle below.
            materialize_pending(app, target);
            if let Some(p) = app.workspace_products.get_mut(&target) {
                match p.animation.as_mut() {
                    Some(a) => {
                        if let Some(v) = play {
                            a.playing = v;
                        }
                        if let Some(s) = speed {
                            // Match the toolbar slider + `ProductAnimation::speed`
                            // range (`0.0..=4.0`), and guard a non-finite speed:
                            // `f32::NAN.clamp(..)` returns NaN, which would poison
                            // `anim.t` → NaN rotation coords on the next tick.
                            a.speed = if s.is_finite() {
                                s.clamp(0.0, 4.0)
                            } else {
                                1.0
                            };
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
        AgentCommand::RunCommand { id } => {
            // A few workbench RUN actions live only as in-panel buttons (no
            // palette command), but still need to be bridge-drivable. Handle
            // those ids here BEFORE the palette lookup; they resolve the active
            // workbench and fire the SAME `*_and_store` path the button calls.
            // Bridge-only on purpose — kept out of the user-facing palette.
            if matches!(id.as_str(), "missionsim.run" | "missionsim.run-monte-carlo") {
                run_missionsim_bridge(app, ch, &id);
            } else if matches!(
                id.as_str(),
                "missionplanner.play"
                    | "missionplanner.pause"
                    | "missionplanner.route"
                    | "missionplanner.los"
            ) {
                run_mission_planner_bridge(app, ch, &id);
            } else if matches!(id.as_str(), "morphogenesis.play" | "morphogenesis.pause") {
                run_morphogenesis_bridge(app, ch, &id);
            } else if id.as_str() == "topopt.run" {
                run_topopt_bridge(app, ch, &id);
            } else if id.as_str() == "brep.build" {
                run_brep_bridge(app, ch, &id);
            } else if id.as_str() == "thermo.compute" {
                run_thermo_bridge(app, ch, &id);
            } else if id.as_str() == "quantum.run" {
                run_quantum_bridge(app, ch, &id);
            } else if id.as_str() == "nodegraph.eval" {
                run_nodegraph_bridge(app, ch, &id);
            } else if id.as_str() == "bondgraph.solve" {
                run_bondgraph_bridge(app, ch, &id);
            } else if id.as_str() == "surrogate.train" {
                run_surrogate_bridge(app, ch, &id);
            } else if id.as_str() == "optics.compute" {
                run_optics_bridge(app, ch, &id);
            } else if id.as_str() == "acoustics.compute" {
                run_acoustics_bridge(app, ch, &id);
            } else if id.as_str() == "waveform.parse" {
                run_waveform_bridge(app, ch, &id);
            } else {
                // Resolve `id` against the EXISTING command-palette registry and
                // invoke the matching command through the SAME `(cmd.invoke)(app)`
                // function pointer a user click / Ctrl+P selection runs (see
                // `commands::dispatch`'s `Static` arm) — no action logic is
                // duplicated here. An unknown id is a feed note + no-op (no panic).
                match crate::commands::static_commands()
                    .iter()
                    .find(|c| c.id.0 == id)
                {
                    Some(cmd) => {
                        (cmd.invoke)(app);
                        crate::assistant_workbench::append_feed_note(
                            app,
                            ch,
                            "Claude",
                            &format!("ran {id}"),
                            "result",
                        );
                    }
                    None => {
                        crate::assistant_workbench::append_feed_note(
                            app,
                            ch,
                            "Claude",
                            &format!("unknown command id: {id}"),
                            "warn",
                        );
                    }
                }
            }
        }
        AgentCommand::ListCommands => {
            // Enumerate the SAME registry `RunCommand` resolves against so an
            // agent can discover the runnable ids. One feed note listing every
            // static command id; no app state changes.
            let ids: Vec<&str> = crate::commands::static_commands()
                .iter()
                .map(|c| c.id.0)
                .collect();
            let body = format!("commands ({}): {}", ids.len(), ids.join(", "));
            crate::assistant_workbench::append_feed_note(app, ch, "Claude", &body, "result");
        }
        AgentCommand::SetControl {
            name,
            value,
            workbench,
        } => {
            set_control(app, ch, &name, &value, workbench.as_deref());
        }
        AgentCommand::ListControls { workbench } => {
            list_controls(app, ch, workbench.as_deref());
        }
        AgentCommand::ReadReadout { workbench } => {
            read_readout(app, ch, workbench.as_deref());
        }
        AgentCommand::InvokeNamed { name } => {
            invoke_named(app, ch, &name);
        }
        AgentCommand::ListButtons => {
            list_buttons(app, ch);
        }
        AgentCommand::ReadText => {
            read_text(app, ch);
        }
        AgentCommand::SetView { dir } => {
            // Snap the central viewport camera to a canonical ViewCube
            // orientation through the SAME `OrbitCamera::set_view` the ViewCube
            // buttons drive. An unrecognised name is a fail-loud warn note.
            match view_direction_from_str(&dir) {
                Some(vd) => {
                    app.camera_mut().set_view(vd);
                    crate::assistant_workbench::append_feed_note(
                        app,
                        ch,
                        "Claude",
                        &format!("view → {dir}"),
                        "result",
                    );
                }
                None => crate::assistant_workbench::append_feed_note(
                    app,
                    ch,
                    "Claude",
                    &format!("set_view: unknown view direction: {dir:?} (use front/back/left/right/top/bottom/iso)"),
                    "warn",
                ),
            }
        }
        AgentCommand::Orbit { dx_deg, dy_deg } => {
            // Orbit by a degree delta through `OrbitCamera::orbit` (which clamps
            // elevation). Guard non-finite deltas so a stray NaN can't poison the
            // camera angles.
            if dx_deg.is_finite() && dy_deg.is_finite() {
                app.camera_mut().orbit(dx_deg, dy_deg);
                let cam = app.camera_mut();
                let (az, el) = (cam.azimuth_deg, cam.elevation_deg);
                crate::assistant_workbench::append_feed_note(
                    app,
                    ch,
                    "Claude",
                    &format!("orbit → az {az:.1}°, el {el:.1}°"),
                    "result",
                );
            } else {
                crate::assistant_workbench::append_feed_note(
                    app,
                    ch,
                    "Claude",
                    "orbit: non-finite delta ignored",
                    "warn",
                );
            }
        }
        AgentCommand::Zoom { factor } => {
            // Dolly in/out through `OrbitCamera::zoom` (which clamps so the
            // camera can't invert through the target). Guard a non-finite factor.
            if factor.is_finite() {
                app.camera_mut().zoom(factor);
                let dist = app.camera_mut().distance;
                crate::assistant_workbench::append_feed_note(
                    app,
                    ch,
                    "Claude",
                    &format!("zoom → distance {dist:.3}"),
                    "result",
                );
            } else {
                crate::assistant_workbench::append_feed_note(
                    app,
                    ch,
                    "Claude",
                    "zoom: non-finite factor ignored",
                    "warn",
                );
            }
        }
        AgentCommand::FrameAll => {
            // Reframe around the loaded model's AABB through the SAME
            // `frame_current_*` methods the menu uses. Prefer a loaded mesh,
            // fall back to a loaded STL; if neither is loaded leave the camera
            // unchanged and say so (fail-loud, no panic).
            if app.mesh.is_some() {
                app.frame_current_mesh();
                crate::assistant_workbench::append_feed_note(
                    app,
                    ch,
                    "Claude",
                    "framed loaded mesh",
                    "result",
                );
            } else if app.stl.is_some() {
                app.frame_current_stl();
                crate::assistant_workbench::append_feed_note(
                    app,
                    ch,
                    "Claude",
                    "framed loaded STL",
                    "result",
                );
            } else {
                crate::assistant_workbench::append_feed_note(
                    app,
                    ch,
                    "Claude",
                    "frame_all: no mesh or STL loaded (camera unchanged)",
                    "warn",
                );
            }
        }
        AgentCommand::Tile => {
            // Organize the open workspace panels into a balanced grid through the
            // SAME `auto_tile_dock` the tab-strip "Tile" button drives. A clean
            // agent product tab or an empty workspace is a deliberate no-op
            // (warn note, never a panic); otherwise report the panel count.
            if app.dock_agent_only {
                crate::assistant_workbench::append_feed_note(
                    app,
                    ch,
                    "Claude",
                    "tile: this is a single agent unit (nothing to re-grid)",
                    "warn",
                );
            } else {
                let count = app.open_tileable_count();
                if count == 0 {
                    crate::assistant_workbench::append_feed_note(
                        app,
                        ch,
                        "Claude",
                        "tile: no workspace panels open",
                        "warn",
                    );
                } else {
                    app.auto_tile_dock();
                    crate::assistant_workbench::append_feed_note(
                        app,
                        ch,
                        "Claude",
                        &format!("tiled {count} panel(s) into a balanced grid"),
                        "result",
                    );
                }
            }
        }
        AgentCommand::AddSketchPoint { x, y } => {
            // Append a Line-tool vertex to the in-house CAD sketch through the
            // SAME `sketch_add_point` the canvas click uses. Guard non-finite
            // coordinates (a real click can't produce them) so a hostile 1e400
            // can't seed a degenerate profile for a later extrude.
            if !(x.is_finite() && y.is_finite()) {
                crate::assistant_workbench::append_feed_note(
                    app,
                    ch,
                    "Claude",
                    &format!("add_sketch_point: non-finite coordinate ignored ({x}, {y})"),
                    "warn",
                );
            } else {
                app.cad.sketch_add_point(x, y);
                let n = app.cad.sketch_anchor_count();
                crate::assistant_workbench::append_feed_note(
                    app,
                    ch,
                    "Claude",
                    &format!("sketch point ({x}, {y}) — {n} anchor(s)"),
                    "result",
                );
            }
        }
        AgentCommand::AddSketchArc { start, via, end } => {
            // Append a 3-point arc to the in-house CAD sketch through the SAME
            // `sketch_add_arc` the Arc tool uses.
            app.cad.sketch_add_arc(start, via, end);
            let n = app.cad.sketch_anchor_count();
            crate::assistant_workbench::append_feed_note(
                app,
                ch,
                "Claude",
                &format!("sketch arc → {n} anchor(s)"),
                "result",
            );
        }
        AgentCommand::ExtrudeSketch { height } => {
            // Extrude the current sketch profile into a solid through the SAME
            // `add_extrude_from_sketch` + `request_rebuild` path the panel button
            // uses. Validate height > 0 (and finite) and a ≥3-anchor profile;
            // both failures are fail-loud warn notes (no panic, no state change).
            if !(height.is_finite() && height > 0.0) {
                crate::assistant_workbench::append_feed_note(
                    app,
                    ch,
                    "Claude",
                    &format!("extrude_sketch: height must be > 0 (got {height})"),
                    "warn",
                );
            } else {
                // `sketch_points()` borrows `app.cad`; clone the profile so the
                // following `&mut` call to `add_extrude_from_sketch` doesn't
                // overlap the borrow.
                let profile = app.cad.sketch_points();
                if profile.len() < 3 {
                    crate::assistant_workbench::append_feed_note(
                        app,
                        ch,
                        "Claude",
                        &format!(
                            "extrude_sketch: need ≥3 sketch anchors to extrude (have {})",
                            profile.len()
                        ),
                        "warn",
                    );
                } else {
                    app.cad.add_extrude_from_sketch(&profile, height);
                    app.cad.request_rebuild();
                    crate::assistant_workbench::append_feed_note(
                        app,
                        ch,
                        "Claude",
                        &format!("extruded sketch ({} anchors) by {height}", profile.len()),
                        "result",
                    );
                }
            }
        }
        AgentCommand::Add2dLine { x1, y1, x2, y2 } => {
            // Add a line to the in-house 2-D drawing through the SAME
            // `Drawing2D::add(Entity2D::Line)` path the form's `+` button uses.
            // Guard non-finite endpoints (matching add_2d_circle's r guard).
            if !(x1.is_finite() && y1.is_finite() && x2.is_finite() && y2.is_finite()) {
                crate::assistant_workbench::append_feed_note(
                    app,
                    ch,
                    "Claude",
                    "add_2d_line: non-finite coordinate ignored",
                    "warn",
                );
            } else {
                app.draft2d.agent_add_line([x1, y1], [x2, y2]);
                let n = app.draft2d.entity_count();
                crate::assistant_workbench::append_feed_note(
                    app,
                    ch,
                    "Claude",
                    &format!(
                        "2-D line added — {n} entit{}",
                        if n == 1 { "y" } else { "ies" }
                    ),
                    "result",
                );
            }
        }
        AgentCommand::Add2dCircle { cx, cy, r } => {
            // Add a circle to the in-house 2-D drawing through the SAME
            // `Drawing2D::add(Entity2D::Circle)` path the form's `+` button uses.
            // Validate r > 0 (the button also floors at 0.1) — fail-loud on a
            // non-positive radius rather than silently drawing a dot.
            if !(r.is_finite() && r > 0.0) {
                crate::assistant_workbench::append_feed_note(
                    app,
                    ch,
                    "Claude",
                    &format!("add_2d_circle: radius must be > 0 (got {r})"),
                    "warn",
                );
            } else {
                app.draft2d.agent_add_circle([cx, cy], r);
                let n = app.draft2d.entity_count();
                crate::assistant_workbench::append_feed_note(
                    app,
                    ch,
                    "Claude",
                    &format!(
                        "2-D circle added — {n} entit{}",
                        if n == 1 { "y" } else { "ies" }
                    ),
                    "result",
                );
            }
        }
        AgentCommand::NewUnit { .. } => {
            // `new_unit` is a **global-channel bootstrap** command (it opens a
            // brand-new unit), handled by `apply_global`. On a per-unit channel
            // there is no new unit to open, so it is a deliberate no-op here.
        }
    }
}

/// Map a case-insensitive view name to a [`ViewDirection`] for
/// [`SetView`](AgentCommand::SetView). Accepts the seven ViewCube faces
/// (`front` / `back` / `left` / `right` / `top` / `bottom`) plus `iso` (with
/// `isometric` / `home` as friendly aliases). Anything else → `None`, which the
/// caller turns into a fail-loud `warn` note.
fn view_direction_from_str(s: &str) -> Option<ViewDirection> {
    match s.trim().to_ascii_lowercase().as_str() {
        "front" => Some(ViewDirection::Front),
        "back" | "rear" => Some(ViewDirection::Back),
        "left" => Some(ViewDirection::Left),
        "right" => Some(ViewDirection::Right),
        "top" => Some(ViewDirection::Top),
        "bottom" => Some(ViewDirection::Bottom),
        "iso" | "isometric" | "home" => Some(ViewDirection::Iso),
        _ => None,
    }
}

/// Resolve the target workbench for a [`SetControl`](AgentCommand::SetControl) /
/// [`ListControls`](AgentCommand::ListControls): an explicit `workbench` id
/// ([`TabKind::from_id`], case-insensitive, alias-tolerant) when given, else the
/// **active tab's** workbench. `None` means neither resolved (no active tab, or
/// an unknown id) — the caller turns that into a fail-loud `warn` note.
fn resolve_target_kind(app: &ValenxApp, workbench: Option<&str>) -> Option<TabKind> {
    match workbench {
        Some(id) => TabKind::from_id(id),
        None => app.tab_bar.active_kind(),
    }
}

/// Fire a Mission-simulation RUN action from the bridge. The Run / Run
/// Monte-Carlo actions exist only as in-panel buttons; this routes the two
/// bridge ids (`missionsim.run`, `missionsim.run-monte-carlo`) to the SAME
/// `*_and_store` functions the buttons call, so a `RunCommand` drives them.
///
/// The active tab must be a [`TabKind::MissionSim`] (so the run lands on the
/// workbench the user is looking at); otherwise this posts a fail-loud `warn`
/// note and changes nothing. After a run it acks with the workbench status line
/// (which already carries the single-run summary, plus the Monte-Carlo headline
/// for the ensemble path); `read_readout` then returns the full MC summary.
fn run_missionsim_bridge(app: &mut ValenxApp, ch: usize, id: &str) {
    if resolve_target_kind(app, None) != Some(TabKind::MissionSim) {
        crate::assistant_workbench::append_feed_note(
            app,
            ch,
            "Claude",
            &format!("{id}: active tab is not the Mission-simulation workbench"),
            "warn",
        );
        return;
    }
    match id {
        "missionsim.run" => crate::missionsim_workbench::run_and_store(app),
        "missionsim.run-monte-carlo" => crate::missionsim_workbench::run_monte_carlo_and_store(app),
        _ => unreachable!("run_missionsim_bridge called with a non-missionsim id"),
    }
    let status = app.missionsim.status.clone();
    crate::assistant_workbench::append_feed_note(
        app,
        ch,
        "Claude",
        &format!("ran {id} \u{2014} {status}"),
        "result",
    );
}

/// Fire a Mission Planner action from the bridge. The Play / Pause toggles and
/// the `Compute route` / `Compute LoS` actions exist only as in-panel buttons;
/// this routes the bridge ids (`missionplanner.play`, `missionplanner.pause`,
/// `missionplanner.route`, `missionplanner.los`) to the SAME `play` / `pause` /
/// `route` / `los` functions the buttons call, so a `RunCommand` drives
/// real-time playback, computes the tactical A\* route, or computes line of
/// sight.
///
/// The active tab must be a [`TabKind::MissionPlanner`]; otherwise this posts a
/// fail-loud `warn` note and changes nothing. After acting it acks with the
/// planner readout (sim time + entity count + playing/paused + route status).
fn run_mission_planner_bridge(app: &mut ValenxApp, ch: usize, id: &str) {
    if resolve_target_kind(app, None) != Some(TabKind::MissionPlanner) {
        crate::assistant_workbench::append_feed_note(
            app,
            ch,
            "Claude",
            &format!("{id}: active tab is not the Mission Planner workbench"),
            "warn",
        );
        return;
    }
    match id {
        "missionplanner.play" => crate::mission_planner_workbench::play(app),
        "missionplanner.pause" => crate::mission_planner_workbench::pause(app),
        "missionplanner.route" => crate::mission_planner_workbench::route(app),
        "missionplanner.los" => crate::mission_planner_workbench::los(app),
        _ => unreachable!("run_mission_planner_bridge called with a non-missionplanner id"),
    }
    let readout = app
        .mission_planner
        .agent_readout()
        .unwrap_or_else(|| "(no readout)".to_string());
    crate::assistant_workbench::append_feed_note(
        app,
        ch,
        "Claude",
        &format!("ran {id} \u{2014} {readout}"),
        "result",
    );
}

/// Fire a Morphogenesis playback action from the bridge. The Play / Pause
/// toggles exist only as in-panel buttons; this routes the two bridge ids
/// (`morphogenesis.play`, `morphogenesis.pause`) to the SAME `play` / `pause`
/// functions the buttons call, so a `RunCommand` drives real-time growth.
///
/// The active tab must be a [`TabKind::Morphogenesis`]; otherwise this posts a
/// fail-loud `warn` note and changes nothing. After toggling it acks with the
/// morphogenesis readout (preset + step count + mean V + playing/paused).
fn run_morphogenesis_bridge(app: &mut ValenxApp, ch: usize, id: &str) {
    if resolve_target_kind(app, None) != Some(TabKind::Morphogenesis) {
        crate::assistant_workbench::append_feed_note(
            app,
            ch,
            "Claude",
            &format!("{id}: active tab is not the Morphogenesis workbench"),
            "warn",
        );
        return;
    }
    match id {
        "morphogenesis.play" => crate::morphogenesis_workbench::play(app),
        "morphogenesis.pause" => crate::morphogenesis_workbench::pause(app),
        _ => unreachable!("run_morphogenesis_bridge called with a non-morphogenesis id"),
    }
    let readout = app
        .morphogenesis
        .agent_readout()
        .unwrap_or_else(|| "(no readout)".to_string());
    crate::assistant_workbench::append_feed_note(
        app,
        ch,
        "Claude",
        &format!("ran {id} \u{2014} {readout}"),
        "result",
    );
}

/// Fire the Topology Optimization run from the bridge. The **Run optimization**
/// action exists only as an in-panel button; this routes the bridge id
/// `topopt.run` to the SAME `run` function the button calls, so a `RunCommand`
/// drives the full SIMP optimisation.
///
/// The active tab must be a [`TabKind::TopOpt`]; otherwise this posts a
/// fail-loud `warn` note and changes nothing. After running it acks with the
/// topopt readout (load case + grid + iterations + final compliance + volume).
fn run_topopt_bridge(app: &mut ValenxApp, ch: usize, id: &str) {
    if resolve_target_kind(app, None) != Some(TabKind::TopOpt) {
        crate::assistant_workbench::append_feed_note(
            app,
            ch,
            "Claude",
            &format!("{id}: active tab is not the Topology Optimization workbench"),
            "warn",
        );
        return;
    }
    crate::topopt_workbench::run(app);
    let readout = app
        .topopt
        .agent_readout()
        .unwrap_or_else(|| "(no readout)".to_string());
    crate::assistant_workbench::append_feed_note(
        app,
        ch,
        "Claude",
        &format!("ran {id} \u{2014} {readout}"),
        "result",
    );
}

/// Fire the Part B-Rep build from the bridge. The **Build** action exists
/// only as an in-panel button; this routes the bridge id `brep.build` to
/// the SAME `run` function the button calls, so a `RunCommand` drives the
/// full primitive-build → boolean → tessellation pipeline.
///
/// The active tab must be a [`TabKind::BrepCad`]; otherwise this posts a
/// fail-loud `warn` note and changes nothing. After building it acks with
/// the brep readout (primitives + boolean op + solid topology + mesh
/// counts + volume + bounding box).
fn run_brep_bridge(app: &mut ValenxApp, ch: usize, id: &str) {
    if resolve_target_kind(app, None) != Some(TabKind::BrepCad) {
        crate::assistant_workbench::append_feed_note(
            app,
            ch,
            "Claude",
            &format!("{id}: active tab is not the Part B-Rep workbench"),
            "warn",
        );
        return;
    }
    crate::brep_workbench::run(app);
    let readout = app
        .brep
        .agent_readout()
        .unwrap_or_else(|| "(no readout)".to_string());
    crate::assistant_workbench::append_feed_note(
        app,
        ch,
        "Claude",
        &format!("ran {id} \u{2014} {readout}"),
        "result",
    );
}

/// Fire the Thermodynamics compute from the bridge. The **Compute** action
/// exists only as an in-panel button; this routes the bridge id
/// `thermo.compute` to the SAME `run` function the button calls.
///
/// The active tab must be a [`TabKind::Thermo`]; otherwise this posts a
/// fail-loud `warn` note and changes nothing. After computing it acks with
/// the thermo readout (fluid + model + state + Z roots + Psat).
fn run_thermo_bridge(app: &mut ValenxApp, ch: usize, id: &str) {
    if resolve_target_kind(app, None) != Some(TabKind::Thermo) {
        crate::assistant_workbench::append_feed_note(
            app,
            ch,
            "Claude",
            &format!("{id}: active tab is not the Thermodynamics workbench"),
            "warn",
        );
        return;
    }
    crate::thermo_workbench::run(app);
    let readout = app
        .thermo
        .agent_readout()
        .unwrap_or_else(|| "(no readout)".to_string());
    crate::assistant_workbench::append_feed_note(
        app,
        ch,
        "Claude",
        &format!("ran {id} \u{2014} {readout}"),
        "result",
    );
}

/// Fire the Optics thin-lens compute from the bridge. The **Analyze** action
/// exists only as an in-panel button; this routes the bridge id
/// `optics.compute` to the SAME `run` function the button calls.
///
/// The active tab must be a [`TabKind::Optics`]; otherwise this posts a
/// fail-loud `warn` note and changes nothing. After computing it acks with
/// the optics readout (object distance + focal length + image distance +
/// magnification + real/virtual).
fn run_optics_bridge(app: &mut ValenxApp, ch: usize, id: &str) {
    if resolve_target_kind(app, None) != Some(TabKind::Optics) {
        crate::assistant_workbench::append_feed_note(
            app,
            ch,
            "Claude",
            &format!("{id}: active tab is not the Optics workbench"),
            "warn",
        );
        return;
    }
    crate::optics_workbench::run(app);
    let readout = app
        .optics
        .agent_readout()
        .unwrap_or_else(|| "(no readout)".to_string());
    crate::assistant_workbench::append_feed_note(
        app,
        ch,
        "Claude",
        &format!("ran {id} \u{2014} {readout}"),
        "result",
    );
}

/// Fire the Acoustics monopole-radiation compute from the bridge. The
/// **Radiate** action exists only as an in-panel button; this routes the
/// bridge id `acoustics.compute` to the SAME `run` function the button
/// calls.
///
/// The active tab must be a [`TabKind::Acoustics`]; otherwise this posts a
/// fail-loud `warn` note and changes nothing. After computing it acks with
/// the acoustics readout (source radius + surface velocity + frequency +
/// observer distance + radiated pressure + SPL).
fn run_acoustics_bridge(app: &mut ValenxApp, ch: usize, id: &str) {
    if resolve_target_kind(app, None) != Some(TabKind::Acoustics) {
        crate::assistant_workbench::append_feed_note(
            app,
            ch,
            "Claude",
            &format!("{id}: active tab is not the Acoustics workbench"),
            "warn",
        );
        return;
    }
    crate::acoustics_workbench::run(app);
    let readout = app
        .acoustics
        .agent_readout()
        .unwrap_or_else(|| "(no readout)".to_string());
    crate::assistant_workbench::append_feed_note(
        app,
        ch,
        "Claude",
        &format!("ran {id} \u{2014} {readout}"),
        "result",
    );
}

/// Fire the Waveform VCD parse from the bridge. The **Parse** action exists
/// only as an in-panel button; this routes the bridge id `waveform.parse` to
/// the SAME `run` function the button calls.
///
/// The active tab must be a [`TabKind::Waveform`]; otherwise this posts a
/// fail-loud `warn` note and changes nothing. After parsing it acks with the
/// waveform readout (signal count + per-signal name/width/transition count +
/// time range).
fn run_waveform_bridge(app: &mut ValenxApp, ch: usize, id: &str) {
    if resolve_target_kind(app, None) != Some(TabKind::Waveform) {
        crate::assistant_workbench::append_feed_note(
            app,
            ch,
            "Claude",
            &format!("{id}: active tab is not the Waveform workbench"),
            "warn",
        );
        return;
    }
    crate::waveform_workbench::run(app);
    let readout = app
        .waveform
        .agent_readout()
        .unwrap_or_else(|| "(no readout)".to_string());
    crate::assistant_workbench::append_feed_note(
        app,
        ch,
        "Claude",
        &format!("ran {id} \u{2014} {readout}"),
        "result",
    );
}

/// Fire the Quantum circuit run from the bridge. The **Run** action exists
/// only as an in-panel button; this routes the bridge id `quantum.run` to
/// the SAME `run` function the button calls.
///
/// The active tab must be a [`TabKind::Quantum`]; otherwise this posts a
/// fail-loud `warn` note and changes nothing. After running it acks with
/// the quantum readout (qubit count + gate count + basis-state probs).
fn run_quantum_bridge(app: &mut ValenxApp, ch: usize, id: &str) {
    if resolve_target_kind(app, None) != Some(TabKind::Quantum) {
        crate::assistant_workbench::append_feed_note(
            app,
            ch,
            "Claude",
            &format!("{id}: active tab is not the Quantum circuit workbench"),
            "warn",
        );
        return;
    }
    crate::quantum_workbench::run(app);
    let readout = app
        .quantum
        .agent_readout()
        .unwrap_or_else(|| "(no readout)".to_string());
    crate::assistant_workbench::append_feed_note(
        app,
        ch,
        "Claude",
        &format!("ran {id} \u{2014} {readout}"),
        "result",
    );
}

/// Fire the Node Graph topological evaluation from the bridge. The **Evaluate**
/// action exists only as an in-panel button; this routes the bridge id
/// `nodegraph.eval` to the SAME `run` function the button calls, so a
/// `RunCommand` drives the full evaluation pass.
///
/// The active tab must be a [`TabKind::NodeGraph`]; otherwise this posts a
/// fail-loud `warn` note and changes nothing. After evaluating it acks with the
/// node-graph readout (node count + edge count + Output value(s)).
fn run_nodegraph_bridge(app: &mut ValenxApp, ch: usize, id: &str) {
    if resolve_target_kind(app, None) != Some(TabKind::NodeGraph) {
        crate::assistant_workbench::append_feed_note(
            app,
            ch,
            "Claude",
            &format!("{id}: active tab is not the Node Graph workbench"),
            "warn",
        );
        return;
    }
    crate::nodegraph_workbench::run(app);
    let readout = app
        .nodegraph
        .agent_readout()
        .unwrap_or_else(|| "(no readout)".to_string());
    crate::assistant_workbench::append_feed_note(
        app,
        ch,
        "Claude",
        &format!("ran {id} \u{2014} {readout}"),
        "result",
    );
}

/// Fire the Bond Graph derive-then-integrate solve from the bridge. The
/// **Solve** action exists only as an in-panel button; this routes the bridge
/// id `bondgraph.solve` to the SAME `run` function the button calls, so a
/// `RunCommand` derives the bond-graph state equations and integrates them.
///
/// The active tab must be a [`TabKind::BondGraph`]; otherwise this posts a
/// fail-loud `warn` note and changes nothing. After solving it acks with the
/// bond-graph readout (preset + ODE order + natural frequency / damping +
/// final state).
fn run_bondgraph_bridge(app: &mut ValenxApp, ch: usize, id: &str) {
    if resolve_target_kind(app, None) != Some(TabKind::BondGraph) {
        crate::assistant_workbench::append_feed_note(
            app,
            ch,
            "Claude",
            &format!("{id}: active tab is not the Bond Graph workbench"),
            "warn",
        );
        return;
    }
    crate::bondgraph_workbench::run(app);
    let readout = app
        .bondgraph
        .agent_readout()
        .unwrap_or_else(|| "(no readout)".to_string());
    crate::assistant_workbench::append_feed_note(
        app,
        ch,
        "Claude",
        &format!("ran {id} \u{2014} {readout}"),
        "result",
    );
}

/// Fire the Surrogate Model sample-then-train from the bridge. The **Train**
/// action exists only as an in-panel button; this routes the bridge id
/// `surrogate.train` to the SAME `run` function the button calls, so a
/// `RunCommand` samples the truth and fits the MLP.
///
/// The active tab must be a [`TabKind::Surrogate`]; otherwise this posts a
/// fail-loud `warn` note and changes nothing. After training it acks with the
/// surrogate readout (train/test MSE + the live surrogate-vs-true prediction).
fn run_surrogate_bridge(app: &mut ValenxApp, ch: usize, id: &str) {
    if resolve_target_kind(app, None) != Some(TabKind::Surrogate) {
        crate::assistant_workbench::append_feed_note(
            app,
            ch,
            "Claude",
            &format!("{id}: active tab is not the Surrogate Model workbench"),
            "warn",
        );
        return;
    }
    crate::surrogate_workbench::run(app);
    let readout = app
        .surrogate
        .agent_readout()
        .unwrap_or_else(|| "(no readout)".to_string());
    crate::assistant_workbench::append_feed_note(
        app,
        ch,
        "Claude",
        &format!("ran {id} \u{2014} {readout}"),
        "result",
    );
}

/// Apply one [`SetControl`](AgentCommand::SetControl): resolve the target
/// workbench (explicit id or the active tab), then assign `value` to the
/// caption-named control through that workbench's **own** validated `agent_set`.
/// Every failure path posts a `warn` feed note to channel `ch` and changes
/// nothing — an unresolved/unsupported workbench, an unknown caption, or a
/// wrong-typed value — so a bad command is loud but never a panic. A successful
/// set posts an ack note naming the control and its new value.
///
/// Dispatch mirrors [`crate::project_tabs::set_workbench_flag`]: a `match` over
/// the resolved [`TabKind`] routes to the right workbench module. Only the
/// workbenches wired this round have an `agent_set`; the rest fall through to a
/// "no settable controls" note (the honest follow-up-sweep surface).
fn set_control(
    app: &mut ValenxApp,
    ch: usize,
    name: &str,
    value: &AgentValue,
    workbench: Option<&str>,
) {
    let warn = |app: &mut ValenxApp, msg: String| {
        crate::assistant_workbench::append_feed_note(app, ch, "Claude", &msg, "warn");
    };

    let Some(kind) = resolve_target_kind(app, workbench) else {
        warn(
            app,
            match workbench {
                Some(id) => format!("set_control: unknown workbench id: {id}"),
                None => "set_control: no active tab to target (pass a workbench id)".to_string(),
            },
        );
        return;
    };

    // Route to the resolved workbench's validated setter. Each arm calls a
    // `agent_set(name, value) -> Result<(), String>` that owns the caption ->
    // field mapping + range validation for that tool.
    let result: Result<(), String> = match kind {
        TabKind::Uq => app.uq.agent_set(name, value),
        TabKind::Uas => app.uas.agent_set(name, value),
        TabKind::MissionSim => app.missionsim.agent_set(name, value),
        TabKind::MissionPlanner => app.mission_planner.agent_set(name, value),
        TabKind::Morphogenesis => app.morphogenesis.agent_set(name, value),
        TabKind::Survivability => app.survivability.agent_set(name, value),
        TabKind::Draft2d => app.draft2d.agent_set(name, value),
        // ---- agent_set sweep, batch 1 ----
        TabKind::Cfd => app.cfd.agent_set(name, value),
        TabKind::Fem => app.fem.agent_set(name, value),
        TabKind::Gears => app.gears.agent_set(name, value),
        TabKind::Springs => app.springs.agent_set(name, value),
        TabKind::Gasdynamics => app.gasdynamics.agent_set(name, value),
        TabKind::Rotor => app.rotor.agent_set(name, value),
        TabKind::Fluids => app.fluids.agent_set(name, value),
        TabKind::Ocean => app.ocean.agent_set(name, value),
        TabKind::Rom => app.rom.agent_set(name, value),
        TabKind::Reinforcement => app.reinforcement.agent_set(name, value),
        TabKind::Sensors => app.sensors.agent_set(name, value),
        TabKind::Hvac => app.hvac.agent_set(name, value),
        TabKind::Frames => app.frames.agent_set(name, value),
        TabKind::Mbd => app.mbd.agent_set(name, value),
        TabKind::Piping => app.piping.agent_set(name, value),
        // ---- agent_set sweep, batch 2 ----
        TabKind::Aero => app.aero.agent_set(name, value),
        TabKind::Astro => app.astro.agent_set(name, value),
        TabKind::BlackHole => app.blackhole.agent_set(name, value),
        TabKind::Reactdyn => app.reactdyn.agent_set(name, value),
        TabKind::Render => app.render.agent_set(name, value),
        TabKind::Reverse => app.reverse.agent_set(name, value),
        TabKind::Sheetmetal => app.sheetmetal.agent_set(name, value),
        TabKind::Fields => app.fields.agent_set(name, value),
        TabKind::Animate => app.animate.agent_set(name, value),
        TabKind::Fasteners => app.fasteners.agent_set(name, value),
        TabKind::Collision => app.collision.agent_set(name, value),
        TabKind::Interior => app.interior.agent_set(name, value),
        TabKind::Geomatics => app.geomatics.agent_set(name, value),
        TabKind::Genetics => app.genetics.agent_set(name, value),
        TabKind::Neuro => app.neuro.agent_set(name, value),
        TabKind::Ppi => app.ppi.agent_set(name, value),
        TabKind::Autonomy => app.autonomy.agent_set(name, value),
        TabKind::Photogrammetry => app.photogrammetry.agent_set(name, value),
        TabKind::Cosim => app.cosim.agent_set(name, value),
        TabKind::VariantEffect => app.variant_effect.agent_set(name, value),
        TabKind::MeshToolbox => app.mesh_toolbox.agent_set(name, value),
        // ---- agent_set sweep, batch 3 (rocket / engine / cad) ----
        TabKind::Rocket => app.rocket.agent_set(name, value),
        TabKind::Engine => app.engine.agent_set(name, value),
        TabKind::Cad => app.cad.agent_set(name, value),
        TabKind::BrepCad => app.brep.agent_set(name, value),
        TabKind::Thermo => app.thermo.agent_set(name, value),
        TabKind::Quantum => app.quantum.agent_set(name, value),
        TabKind::TopOpt => app.topopt.agent_set(name, value),
        TabKind::NodeGraph => app.nodegraph.agent_set(name, value),
        TabKind::BondGraph => app.bondgraph.agent_set(name, value),
        TabKind::Surrogate => app.surrogate.agent_set(name, value),
        TabKind::Optics => app.optics.agent_set(name, value),
        TabKind::Acoustics => app.acoustics.agent_set(name, value),
        TabKind::Waveform => app.waveform.agent_set(name, value),
        other => Err(format!(
            "set_control: workbench {other:?} ({}) has no settable controls yet",
            kind.label()
        )),
    };

    match result {
        Ok(()) => crate::assistant_workbench::append_feed_note(
            app,
            ch,
            "Claude",
            &format!("set {name} = {value:?}"),
            "result",
        ),
        Err(e) => warn(app, format!("set_control: {e}")),
    }
}

/// Apply one [`ListControls`](AgentCommand::ListControls): resolve the target
/// workbench (explicit id or the active tab) and post a single feed note listing
/// every settable caption that workbench's `agent_set` recognises, so an agent
/// can discover the [`SetControl`](AgentCommand::SetControl) name space. A
/// workbench with no setter yet posts a "no settable controls" note; an
/// unresolved workbench posts a `warn`. No app state changes.
fn list_controls(app: &mut ValenxApp, ch: usize, workbench: Option<&str>) {
    let Some(kind) = resolve_target_kind(app, workbench) else {
        let msg = match workbench {
            Some(id) => format!("list_controls: unknown workbench id: {id}"),
            None => "list_controls: no active tab to target (pass a workbench id)".to_string(),
        };
        crate::assistant_workbench::append_feed_note(app, ch, "Claude", &msg, "warn");
        return;
    };

    let names: &[&str] = match kind {
        TabKind::Uq => crate::uq_workbench::UqWorkbenchState::agent_control_names(),
        TabKind::Uas => crate::uas_workbench::UasWorkbenchState::agent_control_names(),
        TabKind::MissionSim => {
            crate::missionsim_workbench::MissionSimWorkbenchState::agent_control_names()
        }
        TabKind::MissionPlanner => {
            crate::mission_planner_workbench::MissionPlannerWorkbenchState::agent_control_names()
        }
        TabKind::Morphogenesis => {
            crate::morphogenesis_workbench::MorphogenesisWorkbenchState::agent_control_names()
        }
        TabKind::Survivability => {
            crate::survivability_workbench::SurvivabilityWorkbenchState::agent_control_names()
        }
        TabKind::Draft2d => crate::draft2d_workbench::Draft2dWorkbenchState::agent_control_names(),
        // ---- agent_set sweep, batch 1 ----
        TabKind::Cfd => crate::cfd_workbench::CfdWorkbenchState::agent_control_names(),
        TabKind::Fem => crate::fem_workbench::FemWorkbenchState::agent_control_names(),
        TabKind::Gears => crate::gears_workbench::GearsWorkbenchState::agent_control_names(),
        TabKind::Springs => crate::springs_workbench::SpringsWorkbenchState::agent_control_names(),
        TabKind::Gasdynamics => {
            crate::gasdynamics_workbench::GasDynamicsWorkbenchState::agent_control_names()
        }
        TabKind::Rotor => crate::rotor_workbench::RotorWorkbenchState::agent_control_names(),
        TabKind::Fluids => crate::fluids_workbench::FluidsWorkbenchState::agent_control_names(),
        TabKind::Ocean => crate::ocean_workbench::OceanWorkbenchState::agent_control_names(),
        TabKind::Rom => crate::rom_workbench::RomWorkbenchState::agent_control_names(),
        TabKind::Reinforcement => {
            crate::reinforcement_workbench::ReinforcementWorkbenchState::agent_control_names()
        }
        TabKind::Sensors => crate::sensors_workbench::SensorsWorkbenchState::agent_control_names(),
        TabKind::Hvac => crate::hvac_workbench::HvacWorkbenchState::agent_control_names(),
        TabKind::Frames => crate::frames_workbench::FramesWorkbenchState::agent_control_names(),
        TabKind::Mbd => crate::mbd_workbench::MbdWorkbenchState::agent_control_names(),
        TabKind::Piping => crate::piping_workbench::PipingWorkbenchState::agent_control_names(),
        // ---- agent_set sweep, batch 2 ----
        TabKind::Aero => crate::aero_workbench::AeroWorkbenchState::agent_control_names(),
        TabKind::Astro => crate::astro_workbench::AstroWorkbenchState::agent_control_names(),
        TabKind::BlackHole => {
            crate::blackhole_workbench::BlackHoleWorkbenchState::agent_control_names()
        }
        TabKind::Reactdyn => {
            crate::reactdyn_workbench::ReactdynWorkbenchState::agent_control_names()
        }
        TabKind::Render => crate::render_workbench::RenderWorkbenchState::agent_control_names(),
        TabKind::Reverse => crate::reverse_workbench::ReverseWorkbenchState::agent_control_names(),
        TabKind::Sheetmetal => {
            crate::sheetmetal_workbench::SheetmetalWorkbenchState::agent_control_names()
        }
        TabKind::Fields => crate::fields_workbench::FieldsWorkbenchState::agent_control_names(),
        TabKind::Animate => crate::animate_workbench::AnimateWorkbenchState::agent_control_names(),
        TabKind::Fasteners => {
            crate::fasteners_workbench::FastenersWorkbenchState::agent_control_names()
        }
        TabKind::Collision => {
            crate::collision_workbench::CollisionWorkbenchState::agent_control_names()
        }
        TabKind::Interior => {
            crate::interior_workbench::InteriorWorkbenchState::agent_control_names()
        }
        TabKind::Geomatics => {
            crate::geomatics_workbench::GeomaticsWorkbenchState::agent_control_names()
        }
        TabKind::Genetics => {
            crate::genetics_workbench::GeneticsWorkbenchState::agent_control_names()
        }
        TabKind::Neuro => crate::neuro_workbench::NeuroWorkbenchState::agent_control_names(),
        TabKind::Ppi => crate::ppi_workbench::PpiWorkbenchState::agent_control_names(),
        TabKind::Autonomy => {
            crate::autonomy_workbench::AutonomyWorkbenchState::agent_control_names()
        }
        TabKind::Photogrammetry => {
            crate::photogrammetry_workbench::PhotogrammetryWorkbenchState::agent_control_names()
        }
        TabKind::Cosim => crate::cosim_workbench::CosimWorkbenchState::agent_control_names(),
        TabKind::VariantEffect => {
            crate::variant_effect_workbench::VariantEffectWorkbenchState::agent_control_names()
        }
        TabKind::MeshToolbox => crate::mesh_toolbox::MeshToolboxState::agent_control_names(),
        // ---- agent_set sweep, batch 3 (rocket / engine / cad) ----
        TabKind::Rocket => crate::rocket_workbench::RocketWorkbenchState::agent_control_names(),
        TabKind::Engine => crate::engine_workbench::EngineWorkbenchState::agent_control_names(),
        TabKind::Cad => crate::cad_workbench::CadWorkbenchState::agent_control_names(),
        TabKind::BrepCad => crate::brep_workbench::BrepWorkbenchState::agent_control_names(),
        TabKind::Thermo => crate::thermo_workbench::ThermoWorkbenchState::agent_control_names(),
        TabKind::Quantum => crate::quantum_workbench::QuantumWorkbenchState::agent_control_names(),
        TabKind::TopOpt => crate::topopt_workbench::TopOptWorkbenchState::agent_control_names(),
        TabKind::NodeGraph => {
            crate::nodegraph_workbench::NodeGraphWorkbenchState::agent_control_names()
        }
        TabKind::BondGraph => {
            crate::bondgraph_workbench::BondGraphWorkbenchState::agent_control_names()
        }
        TabKind::Surrogate => {
            crate::surrogate_workbench::SurrogateWorkbenchState::agent_control_names()
        }
        TabKind::Optics => crate::optics_workbench::OpticsWorkbenchState::agent_control_names(),
        TabKind::Acoustics => {
            crate::acoustics_workbench::AcousticsWorkbenchState::agent_control_names()
        }
        TabKind::Waveform => {
            crate::waveform_workbench::WaveformWorkbenchState::agent_control_names()
        }
        _ => &[],
    };

    let body = if names.is_empty() {
        format!(
            "no settable controls for workbench {} ({kind:?}) yet",
            kind.label()
        )
    } else {
        format!("controls ({}): {}", names.len(), names.join(", "))
    };
    crate::assistant_workbench::append_feed_note(app, ch, "Claude", &body, "result");
}

/// Apply one [`ReadReadout`](AgentCommand::ReadReadout): resolve the target
/// workbench (explicit id or the active tab), read its **computed result text**
/// through that workbench's read-only `agent_readout`, and post it to channel
/// `ch`'s chat feed as a `result` note — so an agent can read the answer back
/// and self-verify what it just drove (the read half of the live-driving loop).
///
/// Fail-loud, never a panic: an unresolved/unknown workbench posts a `warn`
/// note; a workbench that has not been run yet, or one with no readout wired,
/// posts a `warn` "not run yet / no readout" note. A successful read posts the
/// result text. No app state is modified (purely read-only).
///
/// Dispatch mirrors [`set_control`] / [`list_controls`]: a `match` over the
/// resolved [`TabKind`] routes to the right workbench's `agent_readout`. Only
/// the workbenches wired this round have one; the rest fall through to the
/// "no readout" note (the honest follow-up-sweep surface, exactly like
/// `agent_set` started).
fn read_readout(app: &mut ValenxApp, ch: usize, workbench: Option<&str>) {
    let Some(kind) = resolve_target_kind(app, workbench) else {
        let msg = match workbench {
            Some(id) => format!("read_readout: unknown workbench id: {id}"),
            None => "read_readout: no active tab to target (pass a workbench id)".to_string(),
        };
        crate::assistant_workbench::append_feed_note(app, ch, "Claude", &msg, "warn");
        return;
    };

    // Route to the resolved workbench's read-only readout. Each arm returns
    // `Option<String>` — `Some(text)` once the tool has a computed result (or its
    // error line), `None` when it has not been run yet. Only the workbenches
    // wired this round have an `agent_readout`; `other` means "no readout wired".
    let readout: Option<Option<String>> = match kind {
        TabKind::Uq => Some(app.uq.agent_readout()),
        TabKind::Cfd => Some(app.cfd.agent_readout()),
        TabKind::Fem => Some(app.fem.agent_readout()),
        TabKind::Gears => Some(app.gears.agent_readout()),
        TabKind::Springs => Some(app.springs.agent_readout()),
        TabKind::Gasdynamics => Some(app.gasdynamics.agent_readout()),
        TabKind::Fluids => Some(app.fluids.agent_readout()),
        TabKind::Uas => Some(app.uas.agent_readout()),
        TabKind::MissionSim => Some(app.missionsim.agent_readout()),
        TabKind::MissionPlanner => Some(app.mission_planner.agent_readout()),
        TabKind::Morphogenesis => Some(app.morphogenesis.agent_readout()),
        TabKind::Survivability => Some(app.survivability.agent_readout()),
        TabKind::Genetics => Some(app.genetics.agent_readout()),
        TabKind::TopOpt => Some(app.topopt.agent_readout()),
        TabKind::BrepCad => Some(app.brep.agent_readout()),
        TabKind::Thermo => Some(app.thermo.agent_readout()),
        TabKind::Quantum => Some(app.quantum.agent_readout()),
        TabKind::NodeGraph => Some(app.nodegraph.agent_readout()),
        TabKind::BondGraph => Some(app.bondgraph.agent_readout()),
        TabKind::Surrogate => Some(app.surrogate.agent_readout()),
        TabKind::Optics => Some(app.optics.agent_readout()),
        TabKind::Acoustics => Some(app.acoustics.agent_readout()),
        TabKind::Waveform => Some(app.waveform.agent_readout()),
        _ => None,
    };

    match readout {
        // The workbench is wired AND has a computed result → post it.
        Some(Some(text)) => {
            let body = format!("{} readout: {text}", kind.label());
            crate::assistant_workbench::append_feed_note(app, ch, "Claude", &body, "result");
        }
        // Wired but not run yet → fail-loud warn.
        Some(None) => {
            let body = format!(
                "read_readout: {} ({kind:?}) not run yet / no readout",
                kind.label()
            );
            crate::assistant_workbench::append_feed_note(app, ch, "Claude", &body, "warn");
        }
        // No readout wired for this workbench yet (the follow-up-sweep surface).
        None => {
            let body = format!(
                "read_readout: workbench {} ({kind:?}) has no readout wired yet",
                kind.label()
            );
            crate::assistant_workbench::append_feed_note(app, ch, "Claude", &body, "warn");
        }
    }
}

/// Draw the **active workbench's right-side panel** into a (probe) context, for
/// the headless accessibility probe behind `invoke_named` / [`list_buttons`] /
/// [`read_text`]. Dispatches the active [`TabKind`] to the SAME
/// `draw_<wb>_workbench(app, ctx)` entry point the live `update()` loop calls, so
/// the probe tree is byte-for-byte the panel the user sees (every `draw_*` fn
/// self-gates on its own `show_*` flag, which the caller force-sets for the probe
/// frame). A pure dispatch table — no per-workbench *logic* — mirroring the
/// `set_control` / `read_readout` maps in this file; a `TabKind` not yet mapped
/// (e.g. `Blank`) simply draws nothing, and the caller reports "no buttons".
fn draw_active_workbench_probe(app: &mut ValenxApp, ctx: &egui::Context, kind: TabKind) {
    match kind {
        TabKind::Rocket => crate::rocket_workbench::draw_rocket_workbench(app, ctx),
        TabKind::Engine => crate::engine_workbench::draw_engine_workbench(app, ctx),
        TabKind::Astro => crate::astro_workbench::draw_astro_workbench(app, ctx),
        TabKind::Aero => crate::aero_workbench::draw_aero_workbench(app, ctx),
        TabKind::Gasdynamics => crate::gasdynamics_workbench::draw_gasdynamics_workbench(app, ctx),
        TabKind::Rotor => crate::rotor_workbench::draw_rotor_workbench(app, ctx),
        TabKind::BlackHole => crate::blackhole_workbench::draw_blackhole_workbench(app, ctx),
        TabKind::Cfd => crate::cfd_workbench::draw_cfd_workbench(app, ctx),
        TabKind::Fem => crate::fem_workbench::draw_fem_workbench(app, ctx),
        TabKind::Reactdyn => crate::reactdyn_workbench::draw_reactdyn_workbench(app, ctx),
        TabKind::Fields => crate::fields_workbench::draw_fields_workbench(app, ctx),
        TabKind::Thermo => crate::thermo_workbench::draw_thermo_workbench(app, ctx),
        TabKind::Quantum => crate::quantum_workbench::draw_quantum_workbench(app, ctx),
        TabKind::Cad => crate::cad_workbench::draw_cad_workbench(app, ctx),
        TabKind::BrepCad => crate::brep_workbench::draw_brep_workbench(app, ctx),
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
        TabKind::Ppi => crate::ppi_workbench::draw_ppi_workbench(app, ctx),
        TabKind::Morphogenesis => {
            crate::morphogenesis_workbench::draw_morphogenesis_workbench(app, ctx)
        }
        TabKind::Sensors => crate::sensors_workbench::draw_sensors_workbench(app, ctx),
        TabKind::Autonomy => crate::autonomy_workbench::draw_autonomy_workbench(app, ctx),
        TabKind::Fluids => crate::fluids_workbench::draw_fluids_workbench(app, ctx),
        TabKind::Ocean => crate::ocean_workbench::draw_ocean_workbench(app, ctx),
        TabKind::Rom => crate::rom_workbench::draw_rom_workbench(app, ctx),
        TabKind::Uq => crate::uq_workbench::draw_uq_workbench(app, ctx),
        TabKind::Uas => crate::uas_workbench::draw_uas_workbench(app, ctx),
        TabKind::MissionSim => crate::missionsim_workbench::draw_missionsim_workbench(app, ctx),
        TabKind::MissionPlanner => {
            crate::mission_planner_workbench::draw_mission_planner_workbench(app, ctx)
        }
        TabKind::Survivability => {
            crate::survivability_workbench::draw_survivability_workbench(app, ctx)
        }
        TabKind::Photogrammetry => {
            crate::photogrammetry_workbench::draw_photogrammetry_workbench(app, ctx)
        }
        TabKind::Cosim => crate::cosim_workbench::draw_cosim_workbench(app, ctx),
        TabKind::Mbd => crate::mbd_workbench::draw_mbd_workbench(app, ctx),
        TabKind::TopOpt => crate::topopt_workbench::draw_topopt_workbench(app, ctx),
        TabKind::NodeGraph => crate::nodegraph_workbench::draw_nodegraph_workbench(app, ctx),
        TabKind::BondGraph => crate::bondgraph_workbench::draw_bondgraph_workbench(app, ctx),
        TabKind::Surrogate => crate::surrogate_workbench::draw_surrogate_workbench(app, ctx),
        TabKind::Optics => crate::optics_workbench::draw_optics_workbench(app, ctx),
        TabKind::Acoustics => crate::acoustics_workbench::draw_acoustics_workbench(app, ctx),
        TabKind::Waveform => crate::waveform_workbench::draw_waveform_workbench(app, ctx),
        // No standalone right-panel to probe (an empty project tab).
        TabKind::Blank => {}
    }
}

/// Run a **headless accessibility probe of the active workbench** and return its
/// emitted accesskit nodes — the SAME `(NodeId, Node)` tree a screen reader / a
/// UI-Automation client reads. The shared engine behind the three generic,
/// no-per-workbench-wiring commands (`invoke_named`, [`list_buttons`],
/// [`read_text`]).
///
/// Renders the active workbench's panel into a throwaway accesskit-enabled
/// [`egui::Context`] via [`draw_active_workbench_probe`] (force-setting the
/// panel's `show_*` flag for the probe frame, then restoring it so the probe
/// can't leave a panel toggled on). Returns `None` when there is no active
/// workbench to probe.
///
/// **Why the node ids are usable on the live frame:** an egui widget's
/// `accesskit::NodeId` is a deterministic hash of its egui `Id` (derived from
/// the widget's source location / label), so the id a button gets in this probe
/// context is identical to the id it gets in the real frame — the id resolved
/// here can be queued for `crate::ValenxApp::raw_input_hook` to inject.
fn probe_active_workbench(
    app: &mut ValenxApp,
) -> Option<Vec<(egui::accesskit::NodeId, egui::accesskit::Node)>> {
    let kind = app.tab_bar.active_kind()?;

    // Force the active workbench's `show_*` flag on for the probe (its `draw_*`
    // fn early-returns when hidden), remembering the prior value so we can leave
    // app state exactly as we found it — the probe must be side-effect-free.
    let was_shown = workbench_show_flag(app, kind);
    set_workbench_show_flag(app, kind, true);

    let ctx = egui::Context::default();
    ctx.enable_accesskit();
    let out = ctx.run(egui::RawInput::default(), |ctx| {
        draw_active_workbench_probe(app, ctx, kind);
    });

    // Restore the prior visibility so the probe leaves no trace.
    set_workbench_show_flag(app, kind, was_shown);

    out.platform_output.accesskit_update.map(|u| u.nodes)
}

/// Probe the active workbench panel headlessly and return its **readable text**
/// (static-label captions + the live value of any value-bearing control) as a
/// list of trimmed, non-empty strings. The crate-internal engine the
/// [`crate::self_test`] generic product check reads output from — the SAME
/// `probe_active_workbench` tree [`read_text`] posts to the feed, returned as
/// data instead of a feed note. `None` when there is no active workbench (or the
/// panel emitted no accesskit tree); `Some(vec)` (possibly empty) otherwise.
pub(crate) fn probe_active_workbench_text(app: &mut ValenxApp) -> Option<Vec<String>> {
    use egui::accesskit::Role;
    let nodes = probe_active_workbench(app)?;
    let mut parts: Vec<String> = Vec::new();
    for (_, n) in &nodes {
        let text: Option<String> = match n.role() {
            Role::StaticText => n.name().map(str::to_string),
            Role::TextInput | Role::SpinButton => n
                .value()
                .map(str::to_string)
                .or_else(|| n.name().map(str::to_string)),
            _ => None,
        };
        if let Some(t) = text {
            let t = t.trim();
            if !t.is_empty() {
                parts.push(t.to_string());
            }
        }
    }
    Some(parts)
}

/// Read the current value of the active workbench's `show_*` flag, so
/// [`probe_active_workbench`] can restore it after the probe. Mirrors
/// [`set_workbench_show_flag`].
fn workbench_show_flag(app: &ValenxApp, kind: TabKind) -> bool {
    match kind {
        TabKind::Rocket => app.show_rocket_workbench,
        TabKind::Engine => app.show_engine_workbench,
        TabKind::Astro => app.show_astro_workbench,
        TabKind::Aero => app.show_aero_workbench,
        TabKind::Gasdynamics => app.show_gasdynamics_workbench,
        TabKind::Rotor => app.show_rotor_workbench,
        TabKind::BlackHole => app.show_blackhole_workbench,
        TabKind::Cfd => app.show_cfd_workbench,
        TabKind::Fem => app.show_fem_workbench,
        TabKind::Reactdyn => app.show_reactdyn_workbench,
        TabKind::Fields => app.show_fields_workbench,
        TabKind::Thermo => app.show_thermo_workbench,
        TabKind::Quantum => app.show_quantum_workbench,
        TabKind::Cad => app.show_cad_workbench,
        TabKind::BrepCad => app.show_brep_workbench,
        TabKind::MeshToolbox => app.show_mesh_toolbox,
        TabKind::Sheetmetal => app.show_sheetmetal_workbench,
        TabKind::Reverse => app.show_reverse_workbench,
        TabKind::Draft2d => app.show_draft2d_workbench,
        TabKind::Render => app.show_render_workbench,
        TabKind::Animate => app.show_animate_workbench,
        TabKind::Springs => app.show_springs_workbench,
        TabKind::Gears => app.show_gears_workbench,
        TabKind::Fasteners => app.show_fasteners_workbench,
        TabKind::Frames => app.show_frames_workbench,
        TabKind::Collision => app.show_collision_workbench,
        TabKind::Piping => app.show_piping_workbench,
        TabKind::Hvac => app.show_hvac_workbench,
        TabKind::Reinforcement => app.show_reinforcement_workbench,
        TabKind::Interior => app.show_interior_workbench,
        TabKind::Geomatics => app.show_geomatics_workbench,
        TabKind::Genetics => app.show_genetics_workbench,
        TabKind::Neuro => app.show_neuro_workbench,
        TabKind::VariantEffect => app.show_variant_effect_workbench,
        TabKind::Ppi => app.show_ppi_workbench,
        TabKind::Morphogenesis => app.show_morphogenesis_workbench,
        TabKind::Sensors => app.show_sensors_workbench,
        TabKind::Autonomy => app.show_autonomy_workbench,
        TabKind::Fluids => app.show_fluids_workbench,
        TabKind::Ocean => app.show_ocean_workbench,
        TabKind::Rom => app.show_rom_workbench,
        TabKind::Uq => app.show_uq_workbench,
        TabKind::Uas => app.show_uas_workbench,
        TabKind::MissionSim => app.show_missionsim_workbench,
        TabKind::MissionPlanner => app.show_mission_planner_workbench,
        TabKind::Survivability => app.show_survivability_workbench,
        TabKind::Photogrammetry => app.show_photogrammetry_workbench,
        TabKind::Cosim => app.show_cosim_workbench,
        TabKind::Mbd => app.show_mbd_workbench,
        TabKind::TopOpt => app.show_topopt_workbench,
        TabKind::NodeGraph => app.show_nodegraph_workbench,
        TabKind::BondGraph => app.show_bondgraph_workbench,
        TabKind::Surrogate => app.show_surrogate_workbench,
        TabKind::Optics => app.show_optics_workbench,
        TabKind::Acoustics => app.show_acoustics_workbench,
        TabKind::Waveform => app.show_waveform_workbench,
        TabKind::Blank => false,
    }
}

/// Set the active workbench's `show_*` flag for the probe, then restore it.
/// Mirrors [`crate::project_tabs`]'s `TabKind::show` (which only ever sets
/// `true`); this variant takes the value so [`probe_active_workbench`] can both
/// force the panel open and put the flag back exactly as it was.
fn set_workbench_show_flag(app: &mut ValenxApp, kind: TabKind, v: bool) {
    match kind {
        TabKind::Rocket => app.show_rocket_workbench = v,
        TabKind::Engine => app.show_engine_workbench = v,
        TabKind::Astro => app.show_astro_workbench = v,
        TabKind::Aero => app.show_aero_workbench = v,
        TabKind::Gasdynamics => app.show_gasdynamics_workbench = v,
        TabKind::Rotor => app.show_rotor_workbench = v,
        TabKind::BlackHole => app.show_blackhole_workbench = v,
        TabKind::Cfd => app.show_cfd_workbench = v,
        TabKind::Fem => app.show_fem_workbench = v,
        TabKind::Reactdyn => app.show_reactdyn_workbench = v,
        TabKind::Fields => app.show_fields_workbench = v,
        TabKind::Thermo => app.show_thermo_workbench = v,
        TabKind::Quantum => app.show_quantum_workbench = v,
        TabKind::Cad => app.show_cad_workbench = v,
        TabKind::BrepCad => app.show_brep_workbench = v,
        TabKind::MeshToolbox => app.show_mesh_toolbox = v,
        TabKind::Sheetmetal => app.show_sheetmetal_workbench = v,
        TabKind::Reverse => app.show_reverse_workbench = v,
        TabKind::Draft2d => app.show_draft2d_workbench = v,
        TabKind::Render => app.show_render_workbench = v,
        TabKind::Animate => app.show_animate_workbench = v,
        TabKind::Springs => app.show_springs_workbench = v,
        TabKind::Gears => app.show_gears_workbench = v,
        TabKind::Fasteners => app.show_fasteners_workbench = v,
        TabKind::Frames => app.show_frames_workbench = v,
        TabKind::Collision => app.show_collision_workbench = v,
        TabKind::Piping => app.show_piping_workbench = v,
        TabKind::Hvac => app.show_hvac_workbench = v,
        TabKind::Reinforcement => app.show_reinforcement_workbench = v,
        TabKind::Interior => app.show_interior_workbench = v,
        TabKind::Geomatics => app.show_geomatics_workbench = v,
        TabKind::Genetics => app.show_genetics_workbench = v,
        TabKind::Neuro => app.show_neuro_workbench = v,
        TabKind::VariantEffect => app.show_variant_effect_workbench = v,
        TabKind::Ppi => app.show_ppi_workbench = v,
        TabKind::Morphogenesis => app.show_morphogenesis_workbench = v,
        TabKind::Sensors => app.show_sensors_workbench = v,
        TabKind::Autonomy => app.show_autonomy_workbench = v,
        TabKind::Fluids => app.show_fluids_workbench = v,
        TabKind::Ocean => app.show_ocean_workbench = v,
        TabKind::Rom => app.show_rom_workbench = v,
        TabKind::Uq => app.show_uq_workbench = v,
        TabKind::Uas => app.show_uas_workbench = v,
        TabKind::MissionSim => app.show_missionsim_workbench = v,
        TabKind::MissionPlanner => app.show_mission_planner_workbench = v,
        TabKind::Survivability => app.show_survivability_workbench = v,
        TabKind::Photogrammetry => app.show_photogrammetry_workbench = v,
        TabKind::Cosim => app.show_cosim_workbench = v,
        TabKind::Mbd => app.show_mbd_workbench = v,
        TabKind::TopOpt => app.show_topopt_workbench = v,
        TabKind::NodeGraph => app.show_nodegraph_workbench = v,
        TabKind::BondGraph => app.show_bondgraph_workbench = v,
        TabKind::Surrogate => app.show_surrogate_workbench = v,
        TabKind::Optics => app.show_optics_workbench = v,
        TabKind::Acoustics => app.show_acoustics_workbench = v,
        TabKind::Waveform => app.show_waveform_workbench = v,
        TabKind::Blank => {}
    }
}

/// Is this accesskit node a **clickable button** (so `invoke_named` /
/// [`list_buttons`] should consider it)? egui maps `ui.button(...)` to
/// [`Role::Button`] and a selectable label / toggle to
/// [`Role::ToggleButton`]; both are invokable via the `Default` action, so both
/// count. Other roles (labels, spin buttons, the window root) are excluded.
fn is_clickable(node: &egui::accesskit::Node) -> bool {
    use egui::accesskit::Role;
    matches!(node.role(), Role::Button | Role::ToggleButton)
}

/// Apply one [`InvokeNamed`](AgentCommand::InvokeNamed): probe the active
/// workbench's accessibility tree, resolve `name` to a clickable node, and queue
/// an `accesskit` `Default` action on it for `crate::ValenxApp::raw_input_hook`
/// to inject next frame (so the matching button reports `.clicked()` — the same
/// path a real click takes). Exact name match first, then case-insensitive. No
/// active workbench, or no matching button, posts a fail-loud `warn` note and
/// queues nothing (never a panic); a match posts an ack naming the queued button.
fn invoke_named(app: &mut ValenxApp, ch: usize, name: &str) {
    let warn = |app: &mut ValenxApp, msg: String| {
        crate::assistant_workbench::append_feed_note(app, ch, "Claude", &msg, "warn");
    };

    let Some(nodes) = probe_active_workbench(app) else {
        warn(
            app,
            "invoke_named: no active workbench to drive (open a workbench tab first)".to_string(),
        );
        return;
    };

    // Resolve the target node by accessible name: exact match preferred, then a
    // case-insensitive fallback so an agent that lower-cased a caption still
    // hits. Only clickable nodes (buttons / toggles) are eligible, so a matching
    // *label* never shadows the real button.
    let want = name.trim();
    let resolved = nodes
        .iter()
        .find(|(_, n)| is_clickable(n) && n.name() == Some(want))
        .or_else(|| {
            nodes.iter().find(|(_, n)| {
                is_clickable(n) && n.name().is_some_and(|s| s.eq_ignore_ascii_case(want))
            })
        })
        .map(|(id, _)| *id);

    match resolved {
        Some(id) => {
            // Queue the click for raw_input_hook to inject next frame. The bridge
            // keeps frames alive (unfocused request_repaint), so it fires promptly.
            app.pending_accesskit_actions
                .push((id, egui::accesskit::Action::Default));
            crate::assistant_workbench::append_feed_note(
                app,
                ch,
                "Claude",
                &format!("invoke_named: queued click on \u{201c}{want}\u{201d}"),
                "result",
            );
        }
        None => {
            // List what *is* clickable so the agent can correct the name.
            let avail: Vec<&str> = nodes
                .iter()
                .filter(|(_, n)| is_clickable(n))
                .filter_map(|(_, n)| n.name())
                .collect();
            warn(
                app,
                format!(
                    "invoke_named: no clickable button named \u{201c}{want}\u{201d} in the active workbench (available: {})",
                    if avail.is_empty() {
                        "none".to_string()
                    } else {
                        avail.join(", ")
                    }
                ),
            );
        }
    }
}

/// Apply one [`ListButtons`](AgentCommand::ListButtons): probe the active
/// workbench's accessibility tree and post a single feed note listing every
/// clickable button caption, so an agent can discover the
/// [`InvokeNamed`](AgentCommand::InvokeNamed) name space. No app state changes.
fn list_buttons(app: &mut ValenxApp, ch: usize) {
    let Some(nodes) = probe_active_workbench(app) else {
        crate::assistant_workbench::append_feed_note(
            app,
            ch,
            "Claude",
            "list_buttons: no active workbench (open a workbench tab first)",
            "warn",
        );
        return;
    };
    // Dedup while preserving order (egui can emit a caption more than once, e.g.
    // a label node duplicating a button's text).
    let mut names: Vec<&str> = Vec::new();
    for (_, n) in &nodes {
        if is_clickable(n) {
            if let Some(s) = n.name() {
                if !names.contains(&s) {
                    names.push(s);
                }
            }
        }
    }
    let body = if names.is_empty() {
        "buttons (0): (none — this workbench has no named buttons)".to_string()
    } else {
        format!("buttons ({}): {}", names.len(), names.join(", "))
    };
    crate::assistant_workbench::append_feed_note(app, ch, "Claude", &body, "result");
}

/// Apply one [`ReadText`](AgentCommand::ReadText): probe the active workbench's
/// accessibility tree and post its readable text (label / value nodes), so an
/// agent can self-verify a result generically without a per-workbench
/// `agent_readout`. Bounded in length to keep the feed note manageable. No app
/// state changes.
fn read_text(app: &mut ValenxApp, ch: usize) {
    /// Cap on the concatenated text so a huge panel can't post a giant note.
    const MAX_TEXT: usize = 2000;

    let Some(nodes) = probe_active_workbench(app) else {
        crate::assistant_workbench::append_feed_note(
            app,
            ch,
            "Claude",
            "read_text: no active workbench (open a workbench tab first)",
            "warn",
        );
        return;
    };

    use egui::accesskit::Role;
    // Collect readable text: static labels and the value of any value-bearing
    // node (spin buttons / text fields show their current value via `value()`).
    // Skip the window root and structural nodes. Dedup consecutive repeats.
    let mut parts: Vec<String> = Vec::new();
    for (_, n) in &nodes {
        let text: Option<String> = match n.role() {
            // egui maps `ui.label(...)` to `Role::StaticText` (there is no
            // `Role::Label` in accesskit 0.12).
            Role::StaticText => n.name().map(str::to_string),
            // Value-bearing controls expose their current value via `value()`;
            // fall back to the caption when a value isn't set.
            Role::TextInput | Role::SpinButton => n
                .value()
                .map(str::to_string)
                .or_else(|| n.name().map(str::to_string)),
            _ => None,
        };
        if let Some(t) = text {
            let t = t.trim();
            if !t.is_empty() && parts.last().map(String::as_str) != Some(t) {
                parts.push(t.to_string());
            }
        }
    }

    let mut body = parts.join(" \u{2022} ");
    if body.len() > MAX_TEXT {
        body.truncate(MAX_TEXT);
        body.push_str(" \u{2026}");
    }
    let body = if body.is_empty() {
        "read_text: (no readable text in the active workbench panel)".to_string()
    } else {
        format!("text: {body}")
    };
    crate::assistant_workbench::append_feed_note(app, ch, "Claude", &body, "result");
}

/// **LAZY-BUILD materialiser.** Build unit `n`'s deferred product the first time
/// it is needed — when its `workspace:<n>` pane is rendered, or when the unit is
/// animated. `new_unit` records only the product `kind` in
/// [`crate::ValenxApp::pending_products`] (so opening a tab is instant); this
/// turns that record into a live product on demand, through the EXACT same path
/// `new_unit` used to build inline: the `show_3d` reducer (registry meshes + the
/// `dna` text card), then the `show_2d` reducer as a fallback (2-D-only
/// drawings), then a default inspect-spin via `ensure_default_animation`.
///
/// Idempotent and cheap to call every frame: if the product already exists in
/// [`crate::ValenxApp::workspace_products`], or there is nothing pending for `n`,
/// it returns immediately without building. A title-only card the `new_unit`
/// title path may have inserted (when both `kind` and `title` were given) is
/// preserved — its heading is reapplied onto the freshly built product so the
/// caller's title still overrides the product's default heading.
pub(crate) fn materialize_pending(app: &mut ValenxApp, n: usize) {
    // Nothing queued → nothing to build. This is the sole early-out: a product
    // built directly (`show_3d`/`show_2d`) or already materialised has no pending
    // entry, so it is never rebuilt; only a unit `new_unit` queued is built here.
    // (We must NOT early-out merely because `workspace_products` has `n`: the
    // `new_unit` title path may have parked an instant *title-only placeholder*
    // card there while the real product is still pending.)
    let Some(kind) = app.pending_products.remove(&n) else {
        return;
    };

    // Keep the wire `kind` to stamp the confidence badge after the build: the
    // `kind` binding is consumed by the `Show2d` reducer path below.
    let kind_for_badge = kind.clone();

    // Take any title-only placeholder card the `new_unit` title path inserted
    // (when both `kind` and `title` were given): capture its heading and REMOVE
    // it so the build below sees an empty slot — exactly as the old inline
    // `new_unit` build did (so the `show_3d`→`show_2d` fallback keys off a truly
    // empty slot, e.g. for a 2-D-only `rcbeam`). The title is reapplied after the
    // build, preserving the title-overrides-heading contract.
    let title_override = app.workspace_products.remove(&n).map(|p| p.title);

    // Build through the SAME reducer paths `new_unit` used to: first `show_3d`
    // (registry meshes + the `dna` text card); if that produced nothing, the
    // `show_2d` 2-D drawings. An unknown kind builds nothing — a safe no-op.
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

    match app.workspace_products.get_mut(&n) {
        Some(product) => {
            // Default paused inspect-spin so the tile shows the Play/Pause +
            // speed controls (idempotent — a no-op once the product animates or
            // is mesh-less). Mirrors the old inline `new_unit` build.
            product.ensure_default_animation();
            // Reapply the caller's title (if any) so it still overrides the
            // product's default heading, exactly as the eager build did.
            if let Some(title) = title_override {
                product.title = title;
            }
        }
        None => {
            // The kind built nothing (unknown / not yet a registered product).
            // If the `new_unit` title path had parked a title-only card, restore
            // it so the workspace still names itself — matching the eager build,
            // where the title block ran after a no-op build and kept the card.
            if let Some(title) = title_override {
                app.workspace_products.insert(
                    n,
                    crate::WorkspaceProduct {
                        title,
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

    // BASE MATERIAL COLOUR: most catalogue products build a plain `Tri3` mesh
    // and leave `vertex_colors` `None`, so the tile renderer paints them in the
    // single neutral brushed-metal grey — every rocket / beam / motor reads the
    // same. Wired centrally here (the one place every product is built) so the
    // whole catalogue gains material variety without editing each producer: for
    // any product that has a `mesh` but authored NO colours of its own, fill a
    // uniform per-product base tone (steel / concrete / copper / … by category;
    // see `crate::materials::base_color_for`) over the renderer's coloured path.
    //
    // Untouched on purpose:
    //   • products that already set `vertex_colors` (the rebuilt per-part
    //     machinery / aero / marine / fasteners, the FEM von-Mises stress map,
    //     the molecule CPK / aero-Cp overlays) — they keep their richer shading;
    //   • card / 2-D / image products (`mesh: None`) — nothing to colour.
    //
    // The fill length is EXACTLY `3 × total_elements()`: every grey product is a
    // pure `Tri3` surface mesh, so `mesh_to_triangle_surface` emits one triangle
    // per element and the coloured path expects three vertices per triangle (see
    // `crate::wgpu_renderer::triangles_to_vertices_colored`). A length that did
    // not match would degrade safely to grey via `dock_layout`'s length guard.
    if let Some(product) = app.workspace_products.get_mut(&n) {
        if product.vertex_colors.is_none() {
            if let Some(loaded) = product.mesh.as_ref() {
                let tris = loaded.mesh.total_elements();
                let base = crate::materials::base_color_for(&kind_for_badge);
                product.vertex_colors = Some(vec![base; 3 * tris]);
            }
        }
    }

    // CONFIDENCE BADGE: stamp one honest validation line as the final readout
    // entry of every built product. Wired centrally here (not in each producer)
    // so the badge is uniform and conflict-free; it then flows into BOTH the
    // workspace readout tile and the agent-feed post below.
    if let Some(product) = app.workspace_products.get_mut(&n) {
        product
            .lines
            .push(crate::confidence::confidence_for(&kind_for_badge).badge_line());
    }

    // SHOW THE BUILD: now that the unit's real product exists, post its readout
    // to the unit's agent feed so opening the product tab shows the actual
    // numbers in the chat (not just the "Unit N ready" line). `materialize_pending`
    // runs at most once per unit (it `remove`d the pending entry up top), so this
    // posts exactly once. Only a product carrying a non-empty `lines` readout is
    // worth posting — a pure 3-D mesh (empty `lines`) or a kind that built
    // nothing posts nothing. The note body is `title` + each readout line; very
    // long bodies are truncated so a runaway readout can't bloat the feed.
    if let Some(product) = app.workspace_products.get(&n) {
        if !product.lines.is_empty() {
            let mut body = product.title.clone();
            for line in &product.lines {
                body.push('\n');
                body.push_str(line);
            }
            // Cap the feed note so an unexpectedly huge readout stays bounded.
            // Truncate at the nearest char boundary at or below the cap so a
            // multi-byte char straddling the limit never panics `truncate`.
            const MAX_FEED_NOTE_BYTES: usize = 4000;
            if body.len() > MAX_FEED_NOTE_BYTES {
                let mut cut = MAX_FEED_NOTE_BYTES;
                while cut > 0 && !body.is_char_boundary(cut) {
                    cut -= 1;
                }
                body.truncate(cut);
                body.push('…');
            }
            crate::assistant_workbench::append_feed_note(app, n, "Claude", &body, "result");
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
    fn run_command_parses_and_invokes_a_camera_view() {
        // The wire form round-trips.
        let rc: AgentCommand =
            serde_json::from_str(r#"{"cmd":"run_command","id":"view.front"}"#).unwrap();
        assert_eq!(
            rc,
            AgentCommand::RunCommand {
                id: "view.front".into()
            }
        );

        // Applying it actually invokes the palette command through the SAME
        // `(cmd.invoke)(app)` path a click uses: the default camera is
        // (az 45, el 25); `view.front`'s ViewDirection snaps it to (0, 0). The
        // changed azimuth/elevation prove the registry command really ran.
        let mut app = ValenxApp::default();
        assert_eq!(app.camera.azimuth_deg, 45.0);
        assert_eq!(app.camera.elevation_deg, 25.0);
        apply(
            &mut app,
            1,
            AgentCommand::RunCommand {
                id: "view.front".into(),
            },
        );
        assert_eq!(app.camera.azimuth_deg, 0.0, "view.front sets azimuth 0");
        assert_eq!(app.camera.elevation_deg, 0.0, "view.front sets elevation 0");
    }

    #[test]
    fn run_command_routed_through_poll_changes_app_state() {
        // End-to-end through the REAL poll/reducer path: an inbound
        // `{"cmd":"run_command","id":"view.top"}` on channel 1 runs the palette
        // command, snapping the camera to Top (az 0, el 90).
        let mut app = ValenxApp::default();
        app.wb_agent_counter = 1;
        let dir = isolate_cmd_dir(&mut app, "run_command_poll");
        let path = cmd_path(&app, 1);
        std::fs::write(&path, "{\"cmd\":\"run_command\",\"id\":\"view.top\"}\n").unwrap();

        poll_and_apply_agent_commands(&mut app);

        assert_eq!(app.camera.elevation_deg, 90.0, "view.top sets elevation 90");
        assert_eq!(app.camera.azimuth_deg, 0.0, "view.top sets azimuth 0");
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn run_command_unknown_id_posts_a_note_and_does_not_panic() {
        // An unknown id must NOT panic; it posts an "unknown command id" feed
        // note and leaves app state untouched. Point the per-unit feed at an
        // isolated dir so we can read the note back.
        let mut app = ValenxApp::default();
        let dir = isolate_cmd_dir(&mut app, "run_command_unknown");
        app.assistant
            .set_feed_path_for_test(dir.join("assistant_feed.jsonl"));

        let before_az = app.camera.azimuth_deg;
        apply(
            &mut app,
            1,
            AgentCommand::RunCommand {
                id: "no.such.command".into(),
            },
        );
        // No state change (nothing ran).
        assert_eq!(app.camera.azimuth_deg, before_az);

        // A `warn` note naming the bad id was posted to unit 1's feed.
        let feed_path = crate::assistant_workbench::unit_feed_path(&app, 1);
        let body = std::fs::read_to_string(&feed_path).expect("unit-1 feed file written");
        let posted = body
            .lines()
            .filter_map(|l| serde_json::from_str::<serde_json::Value>(l).ok())
            .any(|v| {
                v.get("kind").and_then(|k| k.as_str()) == Some("warn")
                    && v.get("detail").and_then(|d| d.as_str()).is_some_and(|d| {
                        d.contains("unknown command id") && d.contains("no.such.command")
                    })
            });
        assert!(
            posted,
            "unknown id posts a warn note naming the id; feed = {body:?}"
        );
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn list_commands_posts_a_non_empty_list_of_ids() {
        // `ListCommands` posts one feed note enumerating every static command
        // id (the same registry `RunCommand` resolves against), so an agent can
        // discover what's runnable. Assert the note is present, non-empty, and
        // mentions a known id.
        let mut app = ValenxApp::default();
        let dir = isolate_cmd_dir(&mut app, "list_commands");
        app.assistant
            .set_feed_path_for_test(dir.join("assistant_feed.jsonl"));

        // Sanity: the wire form parses.
        let lc: AgentCommand = serde_json::from_str(r#"{"cmd":"list_commands"}"#).unwrap();
        assert_eq!(lc, AgentCommand::ListCommands);

        apply(&mut app, 1, AgentCommand::ListCommands);

        let feed_path = crate::assistant_workbench::unit_feed_path(&app, 1);
        let body = std::fs::read_to_string(&feed_path).expect("unit-1 feed file written");
        let n_ids = crate::commands::static_commands().len();
        let posted = body
            .lines()
            .filter_map(|l| serde_json::from_str::<serde_json::Value>(l).ok())
            .any(|v| {
                v.get("detail").and_then(|d| d.as_str()).is_some_and(|d| {
                    d.contains(&format!("commands ({n_ids})"))
                        && d.contains("view.front")
                        && d.contains("run.selected-case")
                })
            });
        assert!(
            posted && n_ids > 0,
            "list_commands posts a non-empty id list; feed = {body:?}"
        );
        let _ = std::fs::remove_dir_all(&dir);
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
        assert_eq!(product.title, "Rocket");
        // The rocket now carries its flight readout (orbit / Δv / max-Q /
        // peak-g) into the unit's agent chat — not just "Unit N ready".
        assert!(
            !product.lines.is_empty(),
            "rocket product carries its flight readout"
        );
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
    fn materialize_pending_posts_the_products_readout_to_the_unit_feed() {
        // SHOW THE BUILD: when a `new_unit`-queued product is materialised on
        // first view, `materialize_pending` must post the product's own readout
        // (its `title` + `lines`) to that unit's agent feed — so opening the
        // product tab shows its real numbers in the chat, not just "Unit N
        // ready". Drive it directly: queue the `dna` kind (a pure text card with
        // a non-empty readout) on unit 1, materialise, then read the unit's feed
        // file back and assert a `result` note carrying the readout is present.
        let mut app = ValenxApp::default();
        app.wb_agent_counter = 1;
        let dir = isolate_cmd_dir(&mut app, "materialize_feed");
        // Point the per-unit FEED at the same isolated dir so we can read it
        // back without colliding with the live app's feed.
        app.assistant
            .set_feed_path_for_test(dir.join("assistant_feed.jsonl"));

        // What `new_unit` records: only the product kind is queued (lazy build).
        app.pending_products.insert(1, "dna".to_string());
        assert!(app.workspace_products.is_empty());

        crate::agent_commands::materialize_pending(&mut app, 1);

        // The product was built (a dna text card with a readout) …
        let product = app
            .workspace_products
            .get(&1)
            .expect("materialize built the dna product");
        assert!(!product.lines.is_empty(), "dna card has a readout");
        let first_line = product.lines[0].clone();

        // … and its readout was posted to unit 1's feed as a `result` note.
        let feed_path = crate::assistant_workbench::unit_feed_path(&app, 1);
        let body = std::fs::read_to_string(&feed_path).expect("unit-1 feed file written");
        // The feed line is a JSON object; assert the build-readout note is there
        // (kind "result", carrying the product title and its first readout row).
        let posted = body
            .lines()
            .filter_map(|l| serde_json::from_str::<serde_json::Value>(l).ok())
            .any(|v| {
                v.get("kind").and_then(|k| k.as_str()) == Some("result")
                    && v.get("detail")
                        .and_then(|d| d.as_str())
                        .is_some_and(|d| d.contains(&product.title) && d.contains(&first_line))
            });
        assert!(
            posted,
            "materialize posts a `result` note with the product readout; feed = {body:?}"
        );
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn materialize_pending_fills_a_base_colour_for_a_previously_grey_product() {
        // A catalogue product that builds a plain `Tri3` mesh but authors NO
        // colours of its own (e.g. `rocket`) used to render in the single
        // neutral grey. `materialize_pending` must now fill a uniform per-
        // product base tone over the renderer's coloured path — and the fill
        // length must be EXACTLY 3 × the mesh's triangle count, or the renderer
        // (dock_layout's length guard) would silently fall back to grey.
        let mut app = ValenxApp::default();
        let dir = isolate_cmd_dir(&mut app, "materialize_base_colour");

        app.pending_products.insert(1, "rocket".to_string());
        crate::agent_commands::materialize_pending(&mut app, 1);

        let product = app
            .workspace_products
            .get(&1)
            .expect("materialize built the rocket product");
        let mesh = product
            .mesh
            .as_ref()
            .expect("rocket product carries a mesh");
        let tris = mesh.mesh.total_elements();
        assert!(tris > 0, "rocket mesh has triangles");

        let colors = product
            .vertex_colors
            .as_ref()
            .expect("a previously-grey product now gets a base colour");
        assert_eq!(
            colors.len(),
            3 * tris,
            "base-colour fill is exactly 3 × triangle count so the renderer takes \
             the coloured path"
        );
        // The fill is the kind's category base colour, applied uniformly.
        let expected = crate::materials::base_color_for("rocket");
        assert!(
            colors.iter().all(|&c| c == expected),
            "every vertex carries the rocket (aerospace) base colour {expected:?}"
        );

        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn materialize_pending_does_not_override_a_per_part_coloured_product() {
        // The rebuilt per-part-coloured products (e.g. the `gearbox`, which sets
        // `vertex_colors: Some(..)` with a colour per gear/shaft/housing) must
        // keep their own richer shading — the central base-colour fill only
        // applies when a product authored NO colours, so a product that already
        // has them is left byte-for-byte untouched.
        let mut app = ValenxApp::default();
        let dir = isolate_cmd_dir(&mut app, "materialize_no_override");

        // The colours the gearbox producer authors for itself (the source of
        // truth the fill must NOT clobber).
        let authored = crate::gearbox_workbench::gearbox_product()
            .vertex_colors
            .expect("gearbox authors its own per-part colours");

        app.pending_products.insert(1, "gearbox".to_string());
        crate::agent_commands::materialize_pending(&mut app, 1);

        let product = app
            .workspace_products
            .get(&1)
            .expect("materialize built the gearbox product");
        let colors = product
            .vertex_colors
            .as_ref()
            .expect("gearbox keeps its per-part colours");
        assert_eq!(
            *colors, authored,
            "the per-part-coloured gearbox is NOT overridden by the base-colour fill"
        );
        // Sanity: it really is more than one flat colour (genuinely per-part),
        // so this isn't accidentally passing on a degenerate single-tone vec.
        assert!(
            colors.iter().any(|&c| c != colors[0]),
            "gearbox colours are genuinely per-part, not a single base tone"
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
                group: None,
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
                group: None,
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
                group: None,
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
                group: None,
            },
        );
        assert_eq!(app.wb_agent_counter, 1, "exactly one unit opened");
        assert!(app.dock_enabled, "the dock is turned on for the unit grid");
        // A bare new_unit (no kind, no title) renders no workspace product.
        assert!(!app.workspace_products.contains_key(&1));
    }

    #[test]
    fn new_unit_with_kind_queues_then_materialises_a_product_into_the_new_unit() {
        // LAZY-BUILD: `new_unit` with `kind:"rocket"` opens unit 1 INSTANTLY and
        // only QUEUES the rocket kind in `pending_products` — nothing is built
        // into `workspace_products` yet. `materialize_pending` then builds it
        // through the same show_3d path a running agent uses (a live LV-1 mesh)
        // and moves the entry across, clearing the pending queue.
        let mut app = ValenxApp::default();
        apply_global(
            &mut app,
            AgentCommand::NewUnit {
                kind: Some("rocket".into()),
                title: None,
                note: None,
                group: None,
            },
        );
        assert_eq!(app.wb_agent_counter, 1);
        // Lazy: queued, not built.
        assert_eq!(
            app.pending_products.get(&1).map(String::as_str),
            Some("rocket"),
            "new_unit queues the kind instead of building it"
        );
        assert!(
            !app.workspace_products.contains_key(&1),
            "no product is built until the pane is viewed"
        );

        // First view of the pane → materialise.
        materialize_pending(&mut app, 1);
        assert!(
            !app.pending_products.contains_key(&1),
            "the pending entry is consumed once materialised"
        );
        let product = app
            .workspace_products
            .get(&1)
            .expect("materialize_pending builds the queued product into the unit");
        let mesh = product
            .mesh
            .as_ref()
            .expect("the rocket kind attaches a live mesh");
        assert!(mesh.path.to_string_lossy().contains("valenx-lv1"));
        // It carries the default inspect-spin so the tile shows the controls.
        assert!(
            product.animation.is_some(),
            "the materialised mesh product gains a default animation"
        );
    }

    #[test]
    fn materialize_pending_preserves_a_title_override_onto_the_built_product() {
        // When both `kind` and `title` are set, `new_unit` queues the kind AND
        // shows an instant title-only card. `materialize_pending` then builds the
        // real product but re-applies the caller's title so it still overrides
        // the product's default heading (and the rocket mesh is live).
        let mut app = ValenxApp::default();
        apply_global(
            &mut app,
            AgentCommand::NewUnit {
                kind: Some("rocket".into()),
                title: Some("LV-1 Heavy".into()),
                note: None,
                group: None,
            },
        );
        // Instant title-only card; kind queued for lazy build.
        assert_eq!(app.workspace_products.get(&1).unwrap().title, "LV-1 Heavy");
        assert!(app.workspace_products.get(&1).unwrap().mesh.is_none());
        assert_eq!(
            app.pending_products.get(&1).map(String::as_str),
            Some("rocket")
        );

        materialize_pending(&mut app, 1);
        let product = app.workspace_products.get(&1).expect("product built");
        assert_eq!(
            product.title, "LV-1 Heavy",
            "title still overrode the heading"
        );
        assert!(
            product.mesh.is_some(),
            "the rocket mesh was built on materialise"
        );
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
                group: None,
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
                group: None,
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
                group: None,
            },
        );
        assert_eq!(app.wb_agent_counter, 3, "per-unit new_unit opens no unit");
        assert!(app.workspace_products.is_empty());
    }

    #[test]
    fn global_poll_opens_a_unit_and_queues_its_product_lazily() {
        // END-TO-END through the REAL poll path with LAZY-BUILD: an agent appends
        // a rich `new_unit` line to the GLOBAL command file (no unit exists yet,
        // wb_agent_counter == 0). The first poll reads the global channel, opens
        // unit 1 INSTANTLY, QUEUES the gear kind in `pending_products` (it does
        // NOT build the mesh up front — so a burst of new_units stays cheap),
        // shows the instant title-only card, and advances the global cursor —
        // all from the file, no UI click. The gear product is built only when
        // `materialize_pending` runs (first pane view / animate).
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
        // LAZY: the kind is queued, not built — no mesh in workspace_products yet.
        assert_eq!(
            app.pending_products.get(&1).map(String::as_str),
            Some("gear"),
            "the global new_unit queued the gear kind for lazy build"
        );
        let card = app
            .workspace_products
            .get(&1)
            .expect("the instant title-only card is shown");
        assert_eq!(card.title, "Reducer", "the title card names the unit");
        assert!(
            card.mesh.is_none(),
            "no mesh is built until the pane is viewed (lazy build)"
        );
        // The global cursor advanced past the one applied line.
        assert_eq!(app.agent_global_cmd_cursor, Some(1));

        // Materialise (as the first pane render would) → the gear mesh builds and
        // the title override is preserved.
        materialize_pending(&mut app, 1);
        assert!(!app.pending_products.contains_key(&1), "pending consumed");
        let product = app.workspace_products.get(&1).expect("gear product built");
        assert_eq!(product.title, "Reducer", "title still overrode the heading");
        assert!(product.mesh.is_some(), "gear kind attaches a live mesh");
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

    #[test]
    fn global_channel_drives_open_workbench_and_set_control_against_active_tab() {
        // FIX: the GLOBAL command file now honours the full drive command set
        // (not just `new_unit`), so an agent can drive EVERYTHING from the one
        // global channel with no unit bootstrap. END-TO-END through the REAL
        // poll path: with NO Workbench+Agent unit (wb_agent_counter == 0) and a
        // Blank active tab, the agent appends `open_workbench` (switch the
        // active tab to UQ) then `set_control` (write a UQ param) to the GLOBAL
        // file. The poll applies both against the active tab/app, and their ack
        // notes land in the GLOBAL feed (channel-0 sentinel), not a per-unit one.
        let mut app = ValenxApp::default();
        let dir = isolate_cmd_dir(&mut app, "global_full_cmds");
        // Route the GLOBAL feed at a file we can read back.
        let global_feed = dir.join("valenx_chat_feed.jsonl");
        app.assistant.set_feed_path_for_test(global_feed.clone());
        // A Blank active tab for `open_workbench` to switch.
        app.tab_bar.open(TabKind::Blank);
        project_tabs::sync_active(&mut app);
        assert_ne!(
            app.uq.params.n_samples, 256,
            "precondition: UQ sample count not already 256"
        );

        let path = global_cmd_path(&app);
        std::fs::write(
            &path,
            concat!(
                "{\"cmd\":\"open_workbench\",\"id\":\"uq\"}\n",
                "{\"cmd\":\"set_control\",\"name\":\"Monte-Carlo samples N\",\"value\":256,\"workbench\":\"uq\"}\n",
            ),
        )
        .unwrap();

        assert_eq!(
            app.wb_agent_counter, 0,
            "no unit exists — global drive only"
        );
        poll_and_apply_agent_commands(&mut app);

        // `open_workbench` switched the active tab to UQ...
        let idx = app.tab_bar.active.expect("a tab is active");
        assert_eq!(
            app.tab_bar.tabs[idx].kind,
            TabKind::Uq,
            "global open_workbench switched the active tab to UQ"
        );
        assert!(app.show_uq_workbench, "the UQ workbench panel is shown");
        // ...and `set_control` wrote the UQ param through the same reducer.
        assert_eq!(
            app.uq.params.n_samples, 256,
            "global set_control wrote the UQ sample count against the active tab"
        );

        // The ack notes were posted to the GLOBAL feed (channel 0), not a
        // per-unit feed — that is what an agent reading the global channel sees.
        let body = std::fs::read_to_string(&global_feed).expect("global feed written");
        let ack_count = body
            .lines()
            .filter_map(|l| serde_json::from_str::<serde_json::Value>(l).ok())
            .filter(|v| v.get("detail").and_then(|d| d.as_str()).is_some())
            .count();
        assert!(
            ack_count >= 1,
            "global-channel commands ack into the global feed; feed = {body:?}"
        );
        // No per-unit feed (u0 / u1) was created for these global acks.
        assert!(
            !crate::assistant_workbench::unit_feed_path(&app, 0).exists(),
            "global acks do not leak into a u0 unit feed"
        );
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn bridge_active_true_when_global_cmd_file_exists() {
        // `bridge_active` drives the faster ~6 fps heartbeat in update.rs. It is
        // true the moment a global command file exists on disk, so a cold agent
        // that just writes the global file is served promptly even with no env
        // override set. (The env-var signals are covered by their own presence;
        // here we exercise the on-disk-file signal in isolation.)
        let mut app = ValenxApp::default();
        let dir = isolate_cmd_dir(&mut app, "bridge_active");
        let path = global_cmd_path(&app);
        // Without the file (and absent env vars) the file-signal is off. We
        // can't assert the whole function is false here because the test
        // process may have the env vars set, so just assert the file flips it on.
        std::fs::write(&path, "{\"cmd\":\"note\",\"text\":\"hi\"}\n").unwrap();
        assert!(
            bridge_active(&app),
            "an existing global command file marks the bridge active"
        );
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
    fn animate_clamps_speed_to_slider_range_and_guards_non_finite() {
        // The stored speed must match the toolbar slider + `ProductAnimation`
        // range (`0.0..=4.0`), and a non-finite speed must not poison `anim.t`.
        let mut app = ValenxApp::default();
        app.workspace_products.insert(1, animated_product());

        // An over-range speed (8.0) clamps to the slider's max (4.0), not 8.0.
        apply(
            &mut app,
            1,
            AgentCommand::Animate {
                n: Some(1),
                play: None,
                speed: Some(8.0),
            },
        );
        assert_eq!(
            app.workspace_products[&1].animation.as_ref().unwrap().speed,
            4.0,
            "speed 8.0 clamps to the 0.0..=4.0 slider max"
        );

        // A non-finite speed (NaN) falls back to a finite default (1.0) rather
        // than storing NaN — `f32::NAN.clamp(..)` would otherwise return NaN.
        apply(
            &mut app,
            1,
            AgentCommand::Animate {
                n: Some(1),
                play: None,
                speed: Some(f32::NAN),
            },
        );
        let speed = app.workspace_products[&1].animation.as_ref().unwrap().speed;
        assert!(speed.is_finite(), "a non-finite speed must not be stored");
        assert_eq!(speed, 1.0, "non-finite speed falls back to 1.0");
    }

    #[test]
    fn animate_on_a_pending_unit_materialises_it_then_starts_playing() {
        // LAZY-BUILD: an agent can animate a unit whose tab was never viewed.
        // `new_unit{kind:rocket}` only queued the kind (nothing in
        // workspace_products). `animate{play:true}` must first build the product
        // via `materialize_pending` (which attaches the default inspect-spin),
        // then flip it to playing — all without the pane ever being rendered.
        let mut app = ValenxApp::default();
        apply_global(
            &mut app,
            AgentCommand::NewUnit {
                kind: Some("rocket".into()),
                title: None,
                note: None,
                group: None,
            },
        );
        assert!(
            !app.workspace_products.contains_key(&1),
            "lazy: the unit is pending, not built, before animate"
        );
        assert_eq!(
            app.pending_products.get(&1).map(String::as_str),
            Some("rocket")
        );

        apply(
            &mut app,
            1,
            AgentCommand::Animate {
                n: None,
                play: Some(true),
                speed: None,
            },
        );

        assert!(
            !app.pending_products.contains_key(&1),
            "animate materialised the pending unit"
        );
        let anim = app.workspace_products[&1]
            .animation
            .as_ref()
            .expect("materialise attached the default inspect-spin animation");
        assert!(
            anim.playing,
            "play:true started the freshly-built product's clock"
        );
    }

    #[test]
    fn animate_clamps_speed_and_defaults_the_target_to_the_channel() {
        // `speed` is clamped to the toolbar slider + `ProductAnimation` range
        // (`0.0..=4.0`), and a missing `n` targets the command file's channel
        // `ch` (here 1) — pausing it.
        let mut app = ValenxApp::default();
        app.workspace_products.insert(1, animated_product());
        // First start it, with an out-of-range speed that must clamp to 4.0.
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
            assert_eq!(anim.speed, 4.0, "speed clamped to the 4.0 ceiling");
        }
        // Now pause it (speed left untouched → stays 4.0).
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
        assert_eq!(anim.speed, 4.0, "omitted speed left the value untouched");
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

    #[test]
    fn new_unit_links_the_product_kind_to_the_tab_and_shows_its_workbench() {
        // WORKBENCH-TOOL-PER-TAB: a `new_unit{kind}` records the registry kind on
        // the freshly-opened product tab's `workbench_kind`, and `apply_global`
        // reconciles so exactly that one `show_*_workbench` flag is on — its
        // tool panel renders on the right alongside the unit's dock.
        let mut app = ValenxApp::default();
        apply_global(
            &mut app,
            AgentCommand::NewUnit {
                kind: Some("fem".into()),
                title: Some("Bracket FEA".into()),
                note: None,
                group: None,
            },
        );
        let active = app.tab_bar.active.expect("the product tab is active");
        assert_eq!(
            app.tab_bar.tabs[active].workbench_kind.as_deref(),
            Some("fem"),
            "the product kind is linked to the tab"
        );
        assert!(
            app.show_fem_workbench,
            "the linked FEM workbench panel is shown for the product tab"
        );
        // No other workbench leaked on (spot-check a spread).
        assert!(!app.show_rocket_workbench);
        assert!(!app.show_car_workbench);
        assert!(!app.show_dcmotor_workbench);
    }

    #[test]
    fn new_unit_group_files_tabs_into_named_collapsible_bands() {
        // Driving two new_units with the SAME group name lands both in ONE
        // band; a third with a different name mints a SECOND band. This is what
        // organises the ~130-product catalogue into manageable category groups.
        let mut app = ValenxApp::default();
        let unit = |kind: &str, group: &str| AgentCommand::NewUnit {
            kind: Some(kind.into()),
            title: Some(kind.into()),
            note: None,
            group: Some(group.into()),
        };
        apply_global(&mut app, unit("rocket", "Aerospace"));
        let g1 = app.tab_bar.tabs[app.tab_bar.active.unwrap()]
            .group
            .clone()
            .expect("first unit is filed into a band");
        apply_global(&mut app, unit("engine", "Aerospace"));
        let g2 = app.tab_bar.tabs[app.tab_bar.active.unwrap()]
            .group
            .clone()
            .expect("second unit is filed into a band");
        assert_eq!(g1, g2, "two units sharing a group name join ONE band");
        apply_global(&mut app, unit("gear", "Machine design"));
        let g3 = app.tab_bar.tabs[app.tab_bar.active.unwrap()]
            .group
            .clone()
            .expect("third unit is filed into a band");
        assert_ne!(g3, g1, "a different group name mints a SECOND band");
        assert_eq!(app.tab_bar.groups.len(), 2, "exactly two bands exist");
        let aero = app.tab_bar.groups.iter().find(|g| g.id == g1).unwrap();
        assert_eq!(aero.name, "Aerospace", "the band is named as requested");
        assert_ne!(aero.color, [0, 0, 0], "the band carries a colour");
    }

    #[test]
    fn new_unit_without_kind_links_no_workbench() {
        // A `kind`-less unit links nothing — `workbench_kind` stays `None` and no
        // workbench flag is forced on (the dock fills the tab as before).
        let mut app = ValenxApp::default();
        apply_global(
            &mut app,
            AgentCommand::NewUnit {
                kind: None,
                title: Some("Scratch".into()),
                note: None,
                group: None,
            },
        );
        let active = app.tab_bar.active.expect("the tab is active");
        assert!(
            app.tab_bar.tabs[active].workbench_kind.is_none(),
            "no kind → no per-tab workbench link"
        );
        assert!(!app.show_fem_workbench);
        assert!(!app.show_rocket_workbench);
    }

    #[test]
    fn new_unit_links_each_of_a_spread_of_kinds() {
        // Round-trip a representative spread through the product-tab path:
        // rocket / fem / dcmotor / pump / gears each links + shows its own panel.
        // The visible flag is read back through a small per-kind probe.
        fn flag_on(app: &ValenxApp, kind: &str) -> bool {
            match kind {
                "rocket" => app.show_rocket_workbench,
                "fem" => app.show_fem_workbench,
                "dcmotor" => app.show_dcmotor_workbench,
                "pump" => app.show_pump_workbench,
                "gears" => app.show_gears_workbench,
                _ => false,
            }
        }
        for kind in ["rocket", "fem", "dcmotor", "pump", "gears"] {
            let mut app = ValenxApp::default();
            apply_global(
                &mut app,
                AgentCommand::NewUnit {
                    kind: Some(kind.into()),
                    title: None,
                    note: None,
                    group: None,
                },
            );
            let active = app.tab_bar.active.expect("active");
            assert_eq!(
                app.tab_bar.tabs[active].workbench_kind.as_deref(),
                Some(kind),
                "{kind} linked to its product tab"
            );
            assert!(flag_on(&app, kind), "{kind}'s workbench panel is shown");
        }
    }

    // ---- SetControl / ListControls / AgentValue -----------------------------

    #[test]
    fn agent_value_parses_untagged_from_each_json_scalar() {
        // The untagged `AgentValue` reads a bare JSON scalar; the arm order means
        // a whole number is an Int (not a lossy Float) and `true` is a Bool.
        fn v(s: &str) -> AgentValue {
            serde_json::from_str(s).unwrap()
        }
        assert_eq!(v("true"), AgentValue::Bool(true));
        assert_eq!(v("4000"), AgentValue::Int(4000));
        assert_eq!(v("0.55"), AgentValue::Float(0.55));
        assert_eq!(v("\"linear\""), AgentValue::Str("linear".into()));

        // Coercions: Int widens to f64; an integral Float reads back as i64; a
        // fractional Float / bool / string is a typed error, never a panic.
        assert_eq!(AgentValue::Int(7).as_f64().unwrap(), 7.0);
        assert_eq!(AgentValue::Float(4000.0).as_i64().unwrap(), 4000);
        assert!(AgentValue::Float(1.5).as_i64().is_err());
        assert!(AgentValue::Bool(true).as_f64().is_err());
        assert!(AgentValue::Int(1).as_bool().is_err());
        assert!(AgentValue::Int(1).as_str().is_err());
    }

    #[test]
    fn set_control_command_round_trips_from_wire() {
        // The wire form `{"cmd":"set_control",...}` parses, with `workbench`
        // optional (defaulting to the active tab).
        let sc: AgentCommand = serde_json::from_str(
            r#"{"cmd":"set_control","name":"Monte-Carlo samples N","value":256,"workbench":"uq"}"#,
        )
        .unwrap();
        assert_eq!(
            sc,
            AgentCommand::SetControl {
                name: "Monte-Carlo samples N".into(),
                value: AgentValue::Int(256),
                workbench: Some("uq".into()),
            }
        );
        let sc2: AgentCommand = serde_json::from_str(
            r#"{"cmd":"set_control","name":"failure threshold t","value":1.5}"#,
        )
        .unwrap();
        assert_eq!(
            sc2,
            AgentCommand::SetControl {
                name: "failure threshold t".into(),
                value: AgentValue::Float(1.5),
                workbench: None,
            }
        );
        let lc: AgentCommand =
            serde_json::from_str(r#"{"cmd":"list_controls","workbench":"uq"}"#).unwrap();
        assert_eq!(
            lc,
            AgentCommand::ListControls {
                workbench: Some("uq".into())
            }
        );
    }

    #[test]
    fn set_control_sets_a_uq_param_and_a_run_uses_it() {
        // End-to-end through the REAL poll/reducer path: an inbound `set_control`
        // on channel 1 (workbench "uq") writes `Monte-Carlo samples N`; the new
        // value is visible in app state AND a subsequent run produces exactly that
        // many output samples (proving the set actually feeds the solver).
        let mut app = ValenxApp::default();
        app.wb_agent_counter = 1;
        let dir = isolate_cmd_dir(&mut app, "set_control_uq");
        let path = cmd_path(&app, 1);
        assert_ne!(
            app.uq.params.n_samples, 256,
            "precondition: not already 256"
        );
        std::fs::write(
            &path,
            "{\"cmd\":\"set_control\",\"name\":\"Monte-Carlo samples N\",\"value\":256,\"workbench\":\"uq\"}\n",
        )
        .unwrap();

        poll_and_apply_agent_commands(&mut app);

        assert_eq!(
            app.uq.params.n_samples, 256,
            "set_control wrote the UQ sample count"
        );
        // The set value drives the solver: a run yields 256 output samples.
        let res = app.uq.run().expect("uq run after set_control");
        assert_eq!(
            res.output_samples.len(),
            256,
            "the run used the agent-set sample count"
        );
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn set_control_targets_the_active_tab_when_workbench_omitted() {
        // With no `workbench` field, the active tab's kind selects the target.
        // Open a UQ tab, then set a coefficient without naming the workbench.
        let mut app = ValenxApp::default();
        app.tab_bar.open(TabKind::Uq);
        project_tabs::sync_active(&mut app);
        apply(
            &mut app,
            1,
            AgentCommand::SetControl {
                name: "a1 (coeff on x1)".into(),
                value: AgentValue::Float(2.5),
                workbench: None,
            },
        );
        assert_eq!(app.uq.params.a1, 2.5, "active-tab routing set the coeff");
    }

    #[test]
    fn set_control_enum_by_name_selects_the_model() {
        // The `response model g` combo is set by its menu word.
        let mut app = ValenxApp::default();
        apply(
            &mut app,
            1,
            AgentCommand::SetControl {
                name: "response model g".into(),
                value: AgentValue::Str("product".into()),
                workbench: Some("uq".into()),
            },
        );
        assert_eq!(
            app.uq.params.model,
            crate::uq_workbench::ModelPreset::Product
        );
    }

    #[test]
    fn set_control_unknown_name_posts_warn_note_and_does_not_panic() {
        // An unknown caption must NOT panic; it posts a `warn` feed note and
        // leaves the params untouched.
        let mut app = ValenxApp::default();
        let dir = isolate_cmd_dir(&mut app, "set_control_unknown");
        app.assistant
            .set_feed_path_for_test(dir.join("assistant_feed.jsonl"));
        let before = app.uq.params.n_samples;

        apply(
            &mut app,
            1,
            AgentCommand::SetControl {
                name: "no such control".into(),
                value: AgentValue::Int(5),
                workbench: Some("uq".into()),
            },
        );
        assert_eq!(app.uq.params.n_samples, before, "nothing changed");

        let feed_path = crate::assistant_workbench::unit_feed_path(&app, 1);
        let body = std::fs::read_to_string(&feed_path).expect("unit-1 feed written");
        let warned = body
            .lines()
            .filter_map(|l| serde_json::from_str::<serde_json::Value>(l).ok())
            .any(|v| {
                v.get("kind").and_then(|k| k.as_str()) == Some("warn")
                    && v.get("detail")
                        .and_then(|d| d.as_str())
                        .is_some_and(|d| d.contains("unknown UQ control"))
            });
        assert!(warned, "unknown name posts a warn note; feed = {body:?}");
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn set_control_type_mismatch_posts_warn_note_and_does_not_panic() {
        // A value of the wrong type (a string into the numeric sample count)
        // must NOT panic; it posts a `warn` note and leaves the field untouched.
        let mut app = ValenxApp::default();
        let dir = isolate_cmd_dir(&mut app, "set_control_typemismatch");
        app.assistant
            .set_feed_path_for_test(dir.join("assistant_feed.jsonl"));
        let before = app.uq.params.n_samples;

        apply(
            &mut app,
            1,
            AgentCommand::SetControl {
                name: "Monte-Carlo samples N".into(),
                value: AgentValue::Str("lots".into()),
                workbench: Some("uq".into()),
            },
        );
        assert_eq!(
            app.uq.params.n_samples, before,
            "bad-typed set changed nothing"
        );

        let feed_path = crate::assistant_workbench::unit_feed_path(&app, 1);
        let body = std::fs::read_to_string(&feed_path).expect("unit-1 feed written");
        let warned = body
            .lines()
            .filter_map(|l| serde_json::from_str::<serde_json::Value>(l).ok())
            .any(|v| {
                v.get("kind").and_then(|k| k.as_str()) == Some("warn")
                    && v.get("detail")
                        .and_then(|d| d.as_str())
                        .is_some_and(|d| d.contains("expected an integer"))
            });
        assert!(warned, "type mismatch posts a warn note; feed = {body:?}");
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn set_control_unknown_workbench_posts_warn_note() {
        // An unknown workbench id is a fail-loud `warn`, not a panic.
        let mut app = ValenxApp::default();
        let dir = isolate_cmd_dir(&mut app, "set_control_badwb");
        app.assistant
            .set_feed_path_for_test(dir.join("assistant_feed.jsonl"));
        apply(
            &mut app,
            1,
            AgentCommand::SetControl {
                name: "whatever".into(),
                value: AgentValue::Int(1),
                workbench: Some("no-such-workbench".into()),
            },
        );
        let feed_path = crate::assistant_workbench::unit_feed_path(&app, 1);
        let body = std::fs::read_to_string(&feed_path).expect("unit-1 feed written");
        assert!(
            body.contains("unknown workbench id"),
            "bad workbench id posts a warn note; feed = {body:?}"
        );
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn read_readout_parses_with_and_without_workbench() {
        // The wire form round-trips; `workbench` is optional.
        let r: AgentCommand =
            serde_json::from_str(r#"{"cmd":"read_readout","workbench":"uq"}"#).unwrap();
        assert_eq!(
            r,
            AgentCommand::ReadReadout {
                workbench: Some("uq".into())
            }
        );
        let r2: AgentCommand = serde_json::from_str(r#"{"cmd":"read_readout"}"#).unwrap();
        assert_eq!(r2, AgentCommand::ReadReadout { workbench: None });
    }

    #[test]
    fn read_readout_routed_through_poll_posts_the_uq_result() {
        // The FULL live-driving loop end-to-end through the REAL poll path on
        // channel 1: (1) `set_control` writes the UQ sample count, (2) the UQ
        // pipeline is run (the same `run_and_store` the Run button calls, which
        // folds the result into the workbench `status` summary), then (3)
        // `read_readout` reads that computed result back into the unit's chat
        // feed as a `result` note — so an agent can self-verify what it drove.
        let mut app = ValenxApp::default();
        app.wb_agent_counter = 1;
        // Make UQ the active tab so the no-`workbench` `read_readout` resolves to
        // it (exercises the active-tab routing too); the `set_control` line still
        // names "uq" explicitly.
        app.tab_bar.open(TabKind::Uq);
        project_tabs::sync_active(&mut app);
        let dir = isolate_cmd_dir(&mut app, "read_readout_uq");
        // Point the per-unit FEED at the isolated dir so we can read it back.
        app.assistant
            .set_feed_path_for_test(dir.join("assistant_feed.jsonl"));

        // Step 1: set the sample count via the poll path.
        let path = cmd_path(&app, 1);
        std::fs::write(
            &path,
            "{\"cmd\":\"set_control\",\"name\":\"Monte-Carlo samples N\",\"value\":256,\"workbench\":\"uq\"}\n",
        )
        .unwrap();
        poll_and_apply_agent_commands(&mut app);
        assert_eq!(app.uq.params.n_samples, 256, "set_control wrote the count");

        // Step 2: run the UQ pipeline (folds the result into `status`).
        crate::uq_workbench::run_and_store(&mut app);
        let status = app.uq.status.clone();
        assert!(
            !status.is_empty() && app.uq.result.is_some(),
            "uq produced a result + status summary"
        );

        // Step 3: read it back via the poll path (no `workbench` → active tab).
        use std::io::Write as _;
        let mut f = std::fs::OpenOptions::new()
            .append(true)
            .open(&path)
            .unwrap();
        writeln!(f, "{{\"cmd\":\"read_readout\"}}").unwrap();
        app.last_agent_poll = None; // bypass the 1s throttle for the test
        poll_and_apply_agent_commands(&mut app);

        // The feed received a `result` note carrying the UQ readout text.
        let feed_path = crate::assistant_workbench::unit_feed_path(&app, 1);
        let body = std::fs::read_to_string(&feed_path).expect("unit-1 feed written");
        let posted = body
            .lines()
            .filter_map(|l| serde_json::from_str::<serde_json::Value>(l).ok())
            .any(|v| {
                v.get("kind").and_then(|k| k.as_str()) == Some("result")
                    && v.get("detail")
                        .and_then(|d| d.as_str())
                        .is_some_and(|d| d.contains("readout") && d.contains(&status))
            });
        assert!(
            posted,
            "read_readout posts a `result` note with the UQ result text; feed = {body:?}"
        );
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn read_readout_not_run_yet_posts_a_warn_note() {
        // A workbench that has not been run yet (empty status / no result) posts a
        // fail-loud `warn` "not run yet / no readout" note rather than a panic or
        // a bogus empty result.
        let mut app = ValenxApp::default();
        let dir = isolate_cmd_dir(&mut app, "read_readout_notrun");
        app.assistant
            .set_feed_path_for_test(dir.join("assistant_feed.jsonl"));
        assert!(app.uq.status.is_empty() && app.uq.result.is_none());

        apply(
            &mut app,
            1,
            AgentCommand::ReadReadout {
                workbench: Some("uq".into()),
            },
        );

        let feed_path = crate::assistant_workbench::unit_feed_path(&app, 1);
        let body = std::fs::read_to_string(&feed_path).expect("unit-1 feed written");
        let warned = body
            .lines()
            .filter_map(|l| serde_json::from_str::<serde_json::Value>(l).ok())
            .any(|v| {
                v.get("kind").and_then(|k| k.as_str()) == Some("warn")
                    && v.get("detail")
                        .and_then(|d| d.as_str())
                        .is_some_and(|d| d.contains("not run yet"))
            });
        assert!(warned, "not-run posts a warn note; feed = {body:?}");
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn read_readout_unknown_workbench_posts_a_warn_note() {
        // An unknown workbench id is a fail-loud `warn`, not a panic.
        let mut app = ValenxApp::default();
        let dir = isolate_cmd_dir(&mut app, "read_readout_badwb");
        app.assistant
            .set_feed_path_for_test(dir.join("assistant_feed.jsonl"));
        apply(
            &mut app,
            1,
            AgentCommand::ReadReadout {
                workbench: Some("no-such-workbench".into()),
            },
        );
        let feed_path = crate::assistant_workbench::unit_feed_path(&app, 1);
        let body = std::fs::read_to_string(&feed_path).expect("unit-1 feed written");
        assert!(
            body.contains("unknown workbench id"),
            "bad workbench id posts a warn note; feed = {body:?}"
        );
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn set_control_routed_through_poll_sets_uas_and_survivability() {
        // Two more representative workbenches set via the REAL poll path on their
        // own channels, proving the active-tab-independent `workbench:` routing.
        let mut app = ValenxApp::default();
        app.wb_agent_counter = 1;
        let dir = isolate_cmd_dir(&mut app, "set_control_multi");
        let path = cmd_path(&app, 1);
        std::fs::write(
            &path,
            concat!(
                "{\"cmd\":\"set_control\",\"name\":\"sensor range (m)\",\"value\":1234.0,\"workbench\":\"uas\"}\n",
                "{\"cmd\":\"set_control\",\"name\":\"charge mass W (kg)\",\"value\":42,\"workbench\":\"survivability\"}\n",
            ),
        )
        .unwrap();

        poll_and_apply_agent_commands(&mut app);

        assert_eq!(app.uas.params.counter.sensor_range_m, 1234.0);
        assert_eq!(app.survivability.params.charge_kg, 42.0);
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn set_control_routed_through_poll_sets_batch1_workbenches() {
        // A representative slice of the batch-1 workbenches set via the REAL
        // poll path on one channel, proving each new `TabKind` arm routes to the
        // right state's `agent_set` (active-tab-independent `workbench:` routing).
        // These five expose their parameters publicly (`pub params` with `pub`
        // fields), so the landed value is asserted directly; the remaining
        // batch-1 workbenches with private fields have their own in-module
        // `agent_set` round-trip tests.
        let mut app = ValenxApp::default();
        app.wb_agent_counter = 1;
        let dir = isolate_cmd_dir(&mut app, "set_control_batch1");
        let path = cmd_path(&app, 1);
        std::fs::write(
            &path,
            concat!(
                "{\"cmd\":\"set_control\",\"name\":\"truncation rank k\",\"value\":5,\"workbench\":\"rom\"}\n",
                "{\"cmd\":\"set_control\",\"name\":\"number of steps\",\"value\":42,\"workbench\":\"fluids\"}\n",
                "{\"cmd\":\"set_control\",\"name\":\"number of waves N\",\"value\":8,\"workbench\":\"ocean\"}\n",
                "{\"cmd\":\"set_control\",\"name\":\"gravity g (m/s^2)\",\"value\":3.71,\"workbench\":\"mbd\"}\n",
                "{\"cmd\":\"set_control\",\"name\":\"blade count n\",\"value\":3,\"workbench\":\"rotor\"}\n",
            ),
        )
        .unwrap();

        poll_and_apply_agent_commands(&mut app);

        assert_eq!(app.rom.params.rank, 5);
        assert_eq!(app.fluids.params.num_steps, 42);
        assert_eq!(app.ocean.params.num_waves, 8);
        assert_eq!(app.mbd.params.gravity, 3.71);
        assert_eq!(app.rotor.blade_count, 3);
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn set_control_routed_through_poll_sets_batch2_workbenches() {
        // A representative slice of the batch-2 workbenches (this sweep) set via
        // the REAL `poll_and_apply_agent_commands` path on one channel, proving
        // each new `TabKind` arm routes to the right state's `agent_set`
        // (active-tab-independent `workbench:` routing). These three keep their
        // parameters private, so success is asserted through the publicly visible
        // ack note `set_control` posts on a routed-and-validated set: a MISSING
        // arm would instead post the "no settable controls yet" warn, so seeing
        // the `set <name> = <value>` ack proves the dispatch arm exists and ran.
        let mut app = ValenxApp::default();
        app.wb_agent_counter = 1;
        let dir = isolate_cmd_dir(&mut app, "set_control_batch2");
        app.assistant
            .set_feed_path_for_test(dir.join("assistant_feed.jsonl"));
        let path = cmd_path(&app, 1);
        std::fs::write(
            &path,
            concat!(
                "{\"cmd\":\"set_control\",\"name\":\"Furniture\",\"value\":\"sofa\",\"workbench\":\"interior\"}\n",
                "{\"cmd\":\"set_control\",\"name\":\"max bounces\",\"value\":4,\"workbench\":\"render\"}\n",
                "{\"cmd\":\"set_control\",\"name\":\"Variants (one per line)\",\"value\":\"p.R273H\",\"workbench\":\"varianteffect\"}\n",
            ),
        )
        .unwrap();

        poll_and_apply_agent_commands(&mut app);

        // Read channel-1's feed and confirm an ack note landed for each set (and
        // that none routed to the "no settable controls yet" fallthrough).
        let feed_path = crate::assistant_workbench::unit_feed_path(&app, 1);
        let body = std::fs::read_to_string(&feed_path).expect("unit-1 feed written");
        let details: Vec<String> = body
            .lines()
            .filter_map(|l| serde_json::from_str::<serde_json::Value>(l).ok())
            .filter_map(|v| v.get("detail").and_then(|d| d.as_str()).map(str::to_string))
            .collect();
        let has = |needle: &str| details.iter().any(|d| d.contains(needle));
        assert!(
            has("set Furniture = Str(\"sofa\")"),
            "interior arm routed + set; feed = {details:?}"
        );
        assert!(
            has("set max bounces = Int(4)"),
            "render arm routed + set; feed = {details:?}"
        );
        assert!(
            has("set Variants (one per line) = Str(\"p.R273H\")"),
            "variant-effect arm routed + set; feed = {details:?}"
        );
        assert!(
            !details.iter().any(|d| d.contains("no settable controls")),
            "no batch-2 set fell through to the unwired fallthrough; feed = {details:?}"
        );
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn set_control_routed_through_poll_sets_batch3_rocket_engine_cad() {
        // The final batch-3 workbenches (rocket / engine / cad) set via the REAL
        // `poll_and_apply_agent_commands` path on one channel, proving each new
        // `TabKind` arm routes to the right state's `agent_set`
        // (active-tab-independent `workbench:` routing). All three keep their
        // parameters private to their own module, so — exactly as batch-2 does —
        // success is asserted through the publicly visible ack note `set_control`
        // posts on a routed-and-validated set: a MISSING arm would instead post
        // the "no settable controls yet" warn, so seeing the `set <name> =
        // <value>` ack proves the dispatch arm exists and ran. Each chosen value
        // is a type the target's `agent_set` accepts without a validation error,
        // so it produces the ack (not a warn).
        let mut app = ValenxApp::default();
        app.wb_agent_counter = 1;
        let dir = isolate_cmd_dir(&mut app, "set_control_batch3");
        app.assistant
            .set_feed_path_for_test(dir.join("assistant_feed.jsonl"));
        let path = cmd_path(&app, 1);
        std::fs::write(
            &path,
            concat!(
                "{\"cmd\":\"set_control\",\"name\":\"strut count N\",\"value\":4,\"workbench\":\"rocket\"}\n",
                "{\"cmd\":\"set_control\",\"name\":\"chamber temp\",\"value\":3500.0,\"workbench\":\"engine\"}\n",
                "{\"cmd\":\"set_control\",\"name\":\"Material density\",\"value\":7850.0,\"workbench\":\"cad\"}\n",
            ),
        )
        .unwrap();

        poll_and_apply_agent_commands(&mut app);

        // Read channel-1's feed and confirm an ack note landed for each set (and
        // that none routed to the "no settable controls yet" fallthrough).
        let feed_path = crate::assistant_workbench::unit_feed_path(&app, 1);
        let body = std::fs::read_to_string(&feed_path).expect("unit-1 feed written");
        let details: Vec<String> = body
            .lines()
            .filter_map(|l| serde_json::from_str::<serde_json::Value>(l).ok())
            .filter_map(|v| v.get("detail").and_then(|d| d.as_str()).map(str::to_string))
            .collect();
        let has = |needle: &str| details.iter().any(|d| d.contains(needle));
        assert!(
            has("set strut count N = Int(4)"),
            "rocket arm routed + set; feed = {details:?}"
        );
        assert!(
            has("set chamber temp = Float(3500.0)"),
            "engine arm routed + set; feed = {details:?}"
        );
        assert!(
            has("set Material density = Float(7850.0)"),
            "cad arm routed + set; feed = {details:?}"
        );
        assert!(
            !details.iter().any(|d| d.contains("no settable controls")),
            "no batch-3 set fell through to the unwired fallthrough; feed = {details:?}"
        );
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn list_controls_posts_the_uq_caption_list() {
        // `list_controls` posts one feed note enumerating the workbench's settable
        // captions, so an agent can discover the SetControl name space.
        let mut app = ValenxApp::default();
        let dir = isolate_cmd_dir(&mut app, "list_controls_uq");
        app.assistant
            .set_feed_path_for_test(dir.join("assistant_feed.jsonl"));
        apply(
            &mut app,
            1,
            AgentCommand::ListControls {
                workbench: Some("uq".into()),
            },
        );
        let feed_path = crate::assistant_workbench::unit_feed_path(&app, 1);
        let body = std::fs::read_to_string(&feed_path).expect("unit-1 feed written");
        let n = crate::uq_workbench::UqWorkbenchState::agent_control_names().len();
        let posted = body
            .lines()
            .filter_map(|l| serde_json::from_str::<serde_json::Value>(l).ok())
            .any(|v| {
                v.get("detail").and_then(|d| d.as_str()).is_some_and(|d| {
                    d.contains(&format!("controls ({n})"))
                        && d.contains("Monte-Carlo samples N")
                        && d.contains("response model g")
                })
            });
        assert!(
            posted && n > 0,
            "list_controls posts the caption list; feed = {body:?}"
        );
        let _ = std::fs::remove_dir_all(&dir);
    }

    // ---- Viewport bridge: camera control + 2-D sketch placement (gap 4) ----

    #[test]
    fn viewport_commands_parse_from_wire_form() {
        // Each new command round-trips from its `{"cmd":...}` wire form.
        let sv: AgentCommand = serde_json::from_str(r#"{"cmd":"set_view","dir":"front"}"#).unwrap();
        assert_eq!(
            sv,
            AgentCommand::SetView {
                dir: "front".into()
            }
        );
        let orb: AgentCommand =
            serde_json::from_str(r#"{"cmd":"orbit","dx_deg":30.0,"dy_deg":-10.0}"#).unwrap();
        assert_eq!(
            orb,
            AgentCommand::Orbit {
                dx_deg: 30.0,
                dy_deg: -10.0
            }
        );
        let z: AgentCommand = serde_json::from_str(r#"{"cmd":"zoom","factor":0.25}"#).unwrap();
        assert_eq!(z, AgentCommand::Zoom { factor: 0.25 });
        let fa: AgentCommand = serde_json::from_str(r#"{"cmd":"frame_all"}"#).unwrap();
        assert_eq!(fa, AgentCommand::FrameAll);
        let sp: AgentCommand =
            serde_json::from_str(r#"{"cmd":"add_sketch_point","x":1.5,"y":2.5}"#).unwrap();
        assert_eq!(sp, AgentCommand::AddSketchPoint { x: 1.5, y: 2.5 });
        let sa: AgentCommand = serde_json::from_str(
            r#"{"cmd":"add_sketch_arc","start":[0.0,0.0],"via":[1.0,1.0],"end":[2.0,0.0]}"#,
        )
        .unwrap();
        assert_eq!(
            sa,
            AgentCommand::AddSketchArc {
                start: [0.0, 0.0],
                via: [1.0, 1.0],
                end: [2.0, 0.0]
            }
        );
        let ex: AgentCommand =
            serde_json::from_str(r#"{"cmd":"extrude_sketch","height":3.0}"#).unwrap();
        assert_eq!(ex, AgentCommand::ExtrudeSketch { height: 3.0 });
        let l2: AgentCommand =
            serde_json::from_str(r#"{"cmd":"add_2d_line","x1":0.0,"y1":0.0,"x2":5.0,"y2":5.0}"#)
                .unwrap();
        assert_eq!(
            l2,
            AgentCommand::Add2dLine {
                x1: 0.0,
                y1: 0.0,
                x2: 5.0,
                y2: 5.0
            }
        );
        let c2: AgentCommand =
            serde_json::from_str(r#"{"cmd":"add_2d_circle","cx":1.0,"cy":2.0,"r":4.0}"#).unwrap();
        assert_eq!(
            c2,
            AgentCommand::Add2dCircle {
                cx: 1.0,
                cy: 2.0,
                r: 4.0
            }
        );
    }

    #[test]
    fn set_view_changes_camera_angles_through_real_poll_path() {
        // End-to-end through `poll_and_apply_agent_commands`: a `set_view`
        // command snaps the central camera to the named canonical view. Start
        // from the default (az 45, el 25), set "front" (0, 0), then "top"
        // (el 90) and assert the azimuth/elevation actually changed.
        let mut app = ValenxApp::default();
        app.wb_agent_counter = 1;
        let dir = isolate_cmd_dir(&mut app, "setview");
        let path = cmd_path(&app, 1);
        std::fs::write(&path, "{\"cmd\":\"set_view\",\"dir\":\"front\"}\n").unwrap();

        // Sanity: default camera is NOT already on the Front view.
        assert!(app.camera_mut().azimuth_deg != 0.0 || app.camera_mut().elevation_deg != 0.0);
        poll_and_apply_agent_commands(&mut app);
        assert_eq!(app.camera_mut().azimuth_deg, 0.0, "front azimuth");
        assert_eq!(app.camera_mut().elevation_deg, 0.0, "front elevation");

        // A second view ("top") moves elevation to 90.
        use std::io::Write as _;
        let mut f = std::fs::OpenOptions::new()
            .append(true)
            .open(&path)
            .unwrap();
        writeln!(f, "{{\"cmd\":\"set_view\",\"dir\":\"top\"}}").unwrap();
        app.last_agent_poll = None; // bypass throttle
        poll_and_apply_agent_commands(&mut app);
        assert_eq!(app.camera_mut().elevation_deg, 90.0, "top elevation");
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn set_view_unknown_direction_posts_warn_and_leaves_camera() {
        // A bogus view name is a fail-loud `warn` note and changes nothing.
        let mut app = ValenxApp::default();
        let dir = isolate_cmd_dir(&mut app, "setview_bad");
        app.assistant
            .set_feed_path_for_test(dir.join("assistant_feed.jsonl"));
        let (az0, el0) = {
            let c = app.camera_mut();
            (c.azimuth_deg, c.elevation_deg)
        };
        apply(
            &mut app,
            1,
            AgentCommand::SetView {
                dir: "sideways".into(),
            },
        );
        // Camera untouched.
        assert_eq!(app.camera_mut().azimuth_deg, az0);
        assert_eq!(app.camera_mut().elevation_deg, el0);
        // A `warn` note was posted.
        let feed_path = crate::assistant_workbench::unit_feed_path(&app, 1);
        let body = std::fs::read_to_string(&feed_path).expect("unit-1 feed written");
        let warned = body
            .lines()
            .filter_map(|l| serde_json::from_str::<serde_json::Value>(l).ok())
            .any(|v| {
                v.get("kind").and_then(|k| k.as_str()) == Some("warn")
                    && v.get("detail")
                        .and_then(|d| d.as_str())
                        .is_some_and(|d| d.contains("unknown view direction"))
            });
        assert!(warned, "unknown view posts a warn note; feed = {body:?}");
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn orbit_and_zoom_drive_the_camera() {
        // Orbit changes azimuth/elevation; zoom changes distance — both through
        // the vetted `OrbitCamera` methods. A non-finite delta/factor is ignored.
        let mut app = ValenxApp::default();
        let az0 = app.camera_mut().azimuth_deg;
        let dist0 = app.camera_mut().distance;
        apply(
            &mut app,
            1,
            AgentCommand::Orbit {
                dx_deg: 20.0,
                dy_deg: 5.0,
            },
        );
        assert!(
            (app.camera_mut().azimuth_deg - (az0 + 20.0)).abs() < 1e-3,
            "orbit advanced azimuth"
        );
        apply(&mut app, 1, AgentCommand::Zoom { factor: 0.5 });
        assert!(
            app.camera_mut().distance < dist0,
            "zoom-in reduced distance"
        );

        // Non-finite inputs are no-ops (guarded), not panics.
        let az_now = app.camera_mut().azimuth_deg;
        apply(
            &mut app,
            1,
            AgentCommand::Orbit {
                dx_deg: f32::NAN,
                dy_deg: 0.0,
            },
        );
        assert_eq!(
            app.camera_mut().azimuth_deg,
            az_now,
            "NaN orbit delta ignored"
        );
    }

    #[test]
    fn frame_all_with_nothing_loaded_posts_warn() {
        // With no mesh/STL loaded, `frame_all` leaves the camera and posts a
        // warn note (never a panic).
        let mut app = ValenxApp::default();
        let dir = isolate_cmd_dir(&mut app, "frameall_empty");
        app.assistant
            .set_feed_path_for_test(dir.join("assistant_feed.jsonl"));
        assert!(app.mesh.is_none() && app.stl.is_none());
        apply(&mut app, 1, AgentCommand::FrameAll);
        let feed_path = crate::assistant_workbench::unit_feed_path(&app, 1);
        let body = std::fs::read_to_string(&feed_path).expect("unit-1 feed written");
        let warned = body
            .lines()
            .filter_map(|l| serde_json::from_str::<serde_json::Value>(l).ok())
            .any(|v| {
                v.get("kind").and_then(|k| k.as_str()) == Some("warn")
                    && v.get("detail")
                        .and_then(|d| d.as_str())
                        .is_some_and(|d| d.contains("no mesh or STL loaded"))
            });
        assert!(
            warned,
            "frame_all with nothing loaded warns; feed = {body:?}"
        );
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn add_sketch_point_appends_a_vertex_through_real_poll_path() {
        // End-to-end through the poll path: two `add_sketch_point` commands grow
        // the in-house CAD sketch anchor count from 0 → 2 (verified via
        // `sketch_points().len()` before/after).
        let mut app = ValenxApp::default();
        app.wb_agent_counter = 1;
        let dir = isolate_cmd_dir(&mut app, "sketchpoint");
        let path = cmd_path(&app, 1);
        assert_eq!(app.cad.sketch_points().len(), 0, "sketch starts empty");
        std::fs::write(
            &path,
            "{\"cmd\":\"add_sketch_point\",\"x\":0.0,\"y\":0.0}\n{\"cmd\":\"add_sketch_point\",\"x\":1.0,\"y\":0.0}\n",
        )
        .unwrap();
        poll_and_apply_agent_commands(&mut app);
        assert_eq!(
            app.cad.sketch_points().len(),
            2,
            "two bridge points appended two anchors"
        );
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn extrude_sketch_nonpositive_height_warns_no_panic() {
        // ExtrudeSketch with height <= 0 must post a `warn` note and not panic.
        // Seed a valid 3-anchor profile first so the height check is what fails.
        let mut app = ValenxApp::default();
        let dir = isolate_cmd_dir(&mut app, "extrude_badheight");
        app.assistant
            .set_feed_path_for_test(dir.join("assistant_feed.jsonl"));
        app.cad.sketch_add_point(0.0, 0.0);
        app.cad.sketch_add_point(1.0, 0.0);
        app.cad.sketch_add_point(1.0, 1.0);
        assert_eq!(app.cad.sketch_points().len(), 3);

        apply(&mut app, 1, AgentCommand::ExtrudeSketch { height: 0.0 });
        apply(&mut app, 1, AgentCommand::ExtrudeSketch { height: -2.0 });

        let feed_path = crate::assistant_workbench::unit_feed_path(&app, 1);
        let body = std::fs::read_to_string(&feed_path).expect("unit-1 feed written");
        let warned = body
            .lines()
            .filter_map(|l| serde_json::from_str::<serde_json::Value>(l).ok())
            .filter(|v| {
                v.get("kind").and_then(|k| k.as_str()) == Some("warn")
                    && v.get("detail")
                        .and_then(|d| d.as_str())
                        .is_some_and(|d| d.contains("height must be > 0"))
            })
            .count();
        assert!(
            warned >= 2,
            "both non-positive heights warned; feed = {body:?}"
        );
    }

    #[test]
    fn extrude_sketch_valid_height_flags_rebuild() {
        // A valid extrude (>0 height, ≥3 anchors) requests a viewport rebuild —
        // the same effect the panel button has. We can't read the private
        // `rebuild_request`, so assert success indirectly via the ack note.
        let mut app = ValenxApp::default();
        let dir = isolate_cmd_dir(&mut app, "extrude_ok");
        app.assistant
            .set_feed_path_for_test(dir.join("assistant_feed.jsonl"));
        app.cad.sketch_add_point(0.0, 0.0);
        app.cad.sketch_add_point(2.0, 0.0);
        app.cad.sketch_add_point(2.0, 2.0);
        app.cad.sketch_add_point(0.0, 2.0);

        apply(&mut app, 1, AgentCommand::ExtrudeSketch { height: 1.5 });

        let feed_path = crate::assistant_workbench::unit_feed_path(&app, 1);
        let body = std::fs::read_to_string(&feed_path).expect("unit-1 feed written");
        let acked = body
            .lines()
            .filter_map(|l| serde_json::from_str::<serde_json::Value>(l).ok())
            .any(|v| {
                v.get("kind").and_then(|k| k.as_str()) == Some("result")
                    && v.get("detail")
                        .and_then(|d| d.as_str())
                        .is_some_and(|d| d.contains("extruded sketch"))
            });
        assert!(acked, "valid extrude posts an ack note; feed = {body:?}");
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn add_2d_line_increments_drawing_entity_count_through_real_poll_path() {
        // End-to-end: an `add_2d_line` command grows the 2-D drawing's entity
        // count by one (verified via `entity_count()` before/after).
        let mut app = ValenxApp::default();
        app.wb_agent_counter = 1;
        let dir = isolate_cmd_dir(&mut app, "add2dline");
        let path = cmd_path(&app, 1);
        let before = app.draft2d.entity_count();
        std::fs::write(
            &path,
            "{\"cmd\":\"add_2d_line\",\"x1\":0.0,\"y1\":0.0,\"x2\":10.0,\"y2\":10.0}\n",
        )
        .unwrap();
        poll_and_apply_agent_commands(&mut app);
        assert_eq!(
            app.draft2d.entity_count(),
            before + 1,
            "bridge line added exactly one entity"
        );
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn add_2d_circle_nonpositive_radius_warns_no_change() {
        // Add2dCircle with r <= 0 posts a `warn` and adds nothing.
        let mut app = ValenxApp::default();
        let dir = isolate_cmd_dir(&mut app, "add2dcircle_badr");
        app.assistant
            .set_feed_path_for_test(dir.join("assistant_feed.jsonl"));
        let before = app.draft2d.entity_count();
        apply(
            &mut app,
            1,
            AgentCommand::Add2dCircle {
                cx: 0.0,
                cy: 0.0,
                r: 0.0,
            },
        );
        assert_eq!(
            app.draft2d.entity_count(),
            before,
            "non-positive radius added nothing"
        );
        let feed_path = crate::assistant_workbench::unit_feed_path(&app, 1);
        let body = std::fs::read_to_string(&feed_path).expect("unit-1 feed written");
        let warned = body
            .lines()
            .filter_map(|l| serde_json::from_str::<serde_json::Value>(l).ok())
            .any(|v| {
                v.get("kind").and_then(|k| k.as_str()) == Some("warn")
                    && v.get("detail")
                        .and_then(|d| d.as_str())
                        .is_some_and(|d| d.contains("radius must be > 0"))
            });
        assert!(warned, "bad-radius circle warns; feed = {body:?}");
        // A valid radius DOES add one.
        apply(
            &mut app,
            1,
            AgentCommand::Add2dCircle {
                cx: 0.0,
                cy: 0.0,
                r: 3.0,
            },
        );
        assert_eq!(
            app.draft2d.entity_count(),
            before + 1,
            "valid circle added one entity"
        );
        let _ = std::fs::remove_dir_all(&dir);
    }

    // ---- Generic accessibility-name bridge (invoke_named / list_buttons / read_text) ----

    #[test]
    fn parses_invoke_named_list_buttons_read_text() {
        // The three new generic commands round-trip from their wire form.
        let inv: AgentCommand =
            serde_json::from_str(r#"{"cmd":"invoke_named","name":"▶ Compute"}"#).unwrap();
        assert_eq!(
            inv,
            AgentCommand::InvokeNamed {
                name: "▶ Compute".into()
            }
        );
        let lb: AgentCommand = serde_json::from_str(r#"{"cmd":"list_buttons"}"#).unwrap();
        assert_eq!(lb, AgentCommand::ListButtons);
        let rt: AgentCommand = serde_json::from_str(r#"{"cmd":"read_text"}"#).unwrap();
        assert_eq!(rt, AgentCommand::ReadText);
    }

    #[test]
    fn probe_resolves_a_named_button_to_an_accesskit_node() {
        // The load-bearing claim: the headless probe of the ACTIVE workbench
        // emits an accesskit tree whose nodes carry the button captions, so a
        // name → NodeId resolution is possible with no per-workbench wiring.
        // Open the Springs workbench (its "▶ Analyze" button is a stable caption).
        let mut app = ValenxApp::default();
        app.tab_bar.open(TabKind::Springs);
        project_tabs::sync_active(&mut app);

        let nodes = probe_active_workbench(&mut app).expect("active workbench yields a tree");
        let analyze = nodes
            .iter()
            .find(|(_, n)| is_clickable(n) && n.name() == Some("▶ Analyze"));
        assert!(
            analyze.is_some(),
            "the Springs '▶ Analyze' button is a clickable, named node in the probe tree; \
             captions seen = {:?}",
            nodes
                .iter()
                .filter(|(_, n)| is_clickable(n))
                .filter_map(|(_, n)| n.name())
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn invoke_named_queues_an_action_and_acks() {
        // `invoke_named` on a real workbench button resolves the name and queues
        // exactly one accesskit Default action for raw_input_hook to inject, and
        // posts a `result` ack — no panic, no app-state mutation beyond the queue.
        let mut app = ValenxApp::default();
        let dir = isolate_cmd_dir(&mut app, "invoke_named_ok");
        app.assistant
            .set_feed_path_for_test(dir.join("assistant_feed.jsonl"));
        app.tab_bar.open(TabKind::Springs);
        project_tabs::sync_active(&mut app);
        assert!(app.pending_accesskit_actions.is_empty());

        apply(
            &mut app,
            1,
            AgentCommand::InvokeNamed {
                name: "▶ Analyze".into(),
            },
        );

        assert_eq!(
            app.pending_accesskit_actions.len(),
            1,
            "exactly one Default action was queued for next-frame injection"
        );
        assert_eq!(
            app.pending_accesskit_actions[0].1,
            egui::accesskit::Action::Default
        );
        let feed_path = crate::assistant_workbench::unit_feed_path(&app, 1);
        let body = std::fs::read_to_string(&feed_path).expect("unit-1 feed written");
        assert!(
            body.contains("queued click"),
            "a successful invoke acks with a result note; feed = {body:?}"
        );
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn invoke_named_is_case_insensitive_fallback() {
        // An agent that lower-cased the caption still hits via the
        // case-insensitive fallback (the '▶' is preserved; only ASCII case is
        // folded, which is enough for the "analyze" word).
        let mut app = ValenxApp::default();
        app.tab_bar.open(TabKind::Springs);
        project_tabs::sync_active(&mut app);
        apply(
            &mut app,
            0,
            AgentCommand::InvokeNamed {
                name: "▶ analyze".into(),
            },
        );
        assert_eq!(
            app.pending_accesskit_actions.len(),
            1,
            "case-insensitive fallback resolved '▶ analyze' to '▶ Analyze'"
        );
    }

    #[test]
    fn invoke_named_unknown_button_warns_and_queues_nothing() {
        // A caption with no matching clickable node is a fail-loud `warn` that
        // also lists what IS available — and queues nothing (never a panic).
        let mut app = ValenxApp::default();
        let dir = isolate_cmd_dir(&mut app, "invoke_named_bad");
        app.assistant
            .set_feed_path_for_test(dir.join("assistant_feed.jsonl"));
        app.tab_bar.open(TabKind::Springs);
        project_tabs::sync_active(&mut app);

        apply(
            &mut app,
            1,
            AgentCommand::InvokeNamed {
                name: "No Such Button".into(),
            },
        );

        assert!(
            app.pending_accesskit_actions.is_empty(),
            "no action queued for an unknown caption"
        );
        let feed_path = crate::assistant_workbench::unit_feed_path(&app, 1);
        let body = std::fs::read_to_string(&feed_path).expect("unit-1 feed written");
        let warned = body
            .lines()
            .filter_map(|l| serde_json::from_str::<serde_json::Value>(l).ok())
            .any(|v| {
                v.get("kind").and_then(|k| k.as_str()) == Some("warn")
                    && v.get("detail")
                        .and_then(|d| d.as_str())
                        .is_some_and(|d| d.contains("no clickable button named"))
            });
        assert!(warned, "unknown button warns; feed = {body:?}");
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn invoke_named_no_active_workbench_warns() {
        // With no active workbench tab there is nothing to probe → a `warn`, not
        // a panic, and nothing queued.
        let mut app = ValenxApp::default();
        let dir = isolate_cmd_dir(&mut app, "invoke_named_noactive");
        app.assistant
            .set_feed_path_for_test(dir.join("assistant_feed.jsonl"));
        // No tab opened → active_kind() is None.
        apply(
            &mut app,
            1,
            AgentCommand::InvokeNamed {
                name: "▶ Analyze".into(),
            },
        );
        assert!(app.pending_accesskit_actions.is_empty());
        let feed_path = crate::assistant_workbench::unit_feed_path(&app, 1);
        let body = std::fs::read_to_string(&feed_path).expect("unit-1 feed written");
        assert!(
            body.contains("no active workbench"),
            "no-active-workbench warns; feed = {body:?}"
        );
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn probe_is_side_effect_free_on_visibility() {
        // The probe force-sets the active workbench's `show_*` flag to render its
        // panel, then MUST restore it — a hidden workbench stays hidden after a
        // probe (so list_buttons / read_text don't pop panels open as a side
        // effect). Springs is active but its `show_*` flag is left FALSE here
        // (we don't call sync_active), so the probe must leave it false.
        let mut app = ValenxApp::default();
        app.tab_bar.open(TabKind::Springs);
        // Deliberately do NOT sync_active → show_springs_workbench stays false.
        assert!(!app.show_springs_workbench);
        let _ = probe_active_workbench(&mut app);
        assert!(
            !app.show_springs_workbench,
            "probe restored the prior (hidden) visibility — no side effect"
        );
    }

    #[test]
    fn list_buttons_lists_the_active_workbench_captions() {
        // `list_buttons` posts a `result` note enumerating the active workbench's
        // clickable captions (so an agent can discover the invoke_named names).
        let mut app = ValenxApp::default();
        let dir = isolate_cmd_dir(&mut app, "list_buttons");
        app.assistant
            .set_feed_path_for_test(dir.join("assistant_feed.jsonl"));
        app.tab_bar.open(TabKind::Springs);
        project_tabs::sync_active(&mut app);

        apply(&mut app, 1, AgentCommand::ListButtons);

        let feed_path = crate::assistant_workbench::unit_feed_path(&app, 1);
        let body = std::fs::read_to_string(&feed_path).expect("unit-1 feed written");
        assert!(
            body.contains("buttons (") && body.contains("Analyze"),
            "list_buttons enumerates the Springs captions incl. Analyze; feed = {body:?}"
        );
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn raw_input_hook_injects_and_drains_queued_actions() {
        // The injection half: a queued action is pushed into the frame's
        // RawInput as an AccessKitActionRequest and the queue is drained (each
        // action fires on exactly one frame).
        use eframe::App as _;
        let mut app = ValenxApp::default();
        // Any NodeId works here — raw_input_hook only forwards it verbatim into
        // the frame's RawInput (resolution from a name happens earlier, in
        // invoke_named). `NodeId` is a public newtype over a u64.
        let target = egui::accesskit::NodeId(0xC0FFEE);
        app.pending_accesskit_actions
            .push((target, egui::accesskit::Action::Default));

        let ctx = egui::Context::default();
        let mut raw = egui::RawInput::default();
        app.raw_input_hook(&ctx, &mut raw);

        assert!(
            app.pending_accesskit_actions.is_empty(),
            "the queue is drained after injection"
        );
        let injected = raw.events.iter().any(|e| {
            matches!(
                e,
                egui::Event::AccessKitActionRequest(req)
                    if req.action == egui::accesskit::Action::Default && req.target == target
            )
        });
        assert!(
            injected,
            "the queued action was injected as an AccessKitActionRequest"
        );
    }
}
