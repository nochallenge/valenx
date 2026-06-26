# Driving valenx from an external AI

valenx is **AI-drivable-first**: an external agent can open tabs, switch
workbenches, set tool inputs, run solves, and read results back **without a
mouse, a screenshot, or window focus**. It does this over a robust **file
bridge** that valenx polls every frame and applies through the *same* vetted
code paths a user click would. The implementation lives in
[`crates/valenx-app/src/agent_commands.rs`](../crates/valenx-app/src/agent_commands.rs);
this document is the operator's guide.

The quickest way in is the bundled helper
[`scripts/valenx-drive.ps1`](../scripts/valenx-drive.ps1):

```powershell
# one command, then it prints what valenx posts back
./scripts/valenx-drive.ps1 open_workbench thermo
./scripts/valenx-drive.ps1 set_control "Temperature [K]" 350
./scripts/valenx-drive.ps1 run_command thermo.compute
./scripts/valenx-drive.ps1 read_readout thermo
```

---

## 1. The channels (files on disk)

Everything is newline-delimited JSON (`.jsonl`), one command/feed entry per
line. Three file *kinds* live side by side in one **base directory**:

| File | Direction | Purpose |
|------|-----------|---------|
| `valenx_chat_cmd.jsonl`      | agent -> valenx | **GLOBAL command channel** (drive everything) |
| `valenx_chat_feed.jsonl`     | valenx -> agent | **GLOBAL feed** (acks, readouts, warnings) |
| `valenx_chat_cmd_u{n}.jsonl` | agent -> valenx | per-unit command channel for Workbench+Agent unit `n` |
| `valenx_chat_feed_u{n}.jsonl`| valenx -> agent | per-unit feed for unit `n` |

**You almost always want the GLOBAL pair** (`valenx_chat_cmd.jsonl` /
`valenx_chat_feed.jsonl`). The per-unit channels exist for the in-app
"Workbench + Agent" tiles, where each tile is its own conversation.

### Where the base directory is

Resolved the same way valenx resolves it (see `agent_commands::global_cmd_path`
and `assistant_workbench::assistant_feed_path`):

- **Command base dir** = the directory of `$VALENX_ASSISTANT_INBOX` if that env
  var is set, otherwise the per-OS **state dir**.
- **Global feed file** = `$VALENX_ASSISTANT_FEED` if set, otherwise
  `<state-dir>/assistant_feed.jsonl`.
- **State dir**: `%APPDATA%\valenx` (Windows), `~/Library/Application
  Support/valenx` (macOS), `$XDG_STATE_HOME/valenx` or `~/.local/state/valenx`
  (Linux).

The recommended setup is to launch valenx with both env vars pointing at a
known temp location, so the bridge is unambiguous:

```
VALENX_ASSISTANT_INBOX = C:/Users/<you>/AppData/Local/Temp/valenx_chat_inbox.jsonl
VALENX_ASSISTANT_FEED  = C:/Users/<you>/AppData/Local/Temp/valenx_chat_feed.jsonl
```

Then the global command file is `…/valenx_chat_cmd.jsonl` in that same `Temp`
directory and the global feed is the `…_feed.jsonl` you named.

### Write rules (important)

- **Append, never rewrite.** valenx keeps a cursor and only applies *newly
  appended* lines; rewriting the file can replay or skip commands.
- **BOM-free UTF-8.** A UTF-16 or BOM-prefixed line will fail to parse and be
  skipped. `valenx-drive.ps1` uses `[System.IO.File]::AppendAllText` with a
  no-BOM `UTF8Encoding($false)` -- do the same if you write the file yourself.
- **One JSON object per line.** A half-written final line is tolerated (skipped)
  but a malformed line is silently dropped.
- Stale command files are **wiped once at launch**, so the first poll that sees
  a freshly created file runs every line in it from the top.

### How fast it applies

valenx polls the files about **once per second** and, while the bridge is
active, schedules a ~6 fps repaint so the poll fires promptly even when the
window is unfocused or in the background. "Active" means either env channel is
set **or** a global command file already exists on disk
(`agent_commands::bridge_active`). A background wake thread also pokes the event
loop a few times a second. Net effect: a command you append is normally applied
within ~1 s without touching the window. (Caveat: a **minimized** window may not
repaint on some platforms; restore it and the queued commands flush
immediately.)

---

## 2. The global channel drives everything

The global `valenx_chat_cmd.jsonl` honours the **full** command vocabulary
(handled by `apply_global`), applied against the **active tab / app**:

- tab ops: `new_tab`, `open_workbench`, `focus_tab`, `rename_tab`, `close_tab`
- inputs & actions: `set_control`, `run_command`
- discovery & read-back: `list_controls`, `list_commands`, `read_readout`
- camera: `set_view`, `orbit`, `zoom`, `frame_all`
- in-house CAD sketch: `add_sketch_point`, `add_sketch_arc`, `extrude_sketch`
- in-house 2-D drafting: `add_2d_line`, `add_2d_circle`
- narration: `note`
- playback: `animate` (carries its own target unit `n`)
- bootstrap: `new_unit` (opens a fresh Workbench+Agent unit; **global-only**)

Acks, readouts, and warnings for global-channel commands are posted to the
**global feed** (`valenx_chat_feed.jsonl`). So you can drive the whole app from
the single global pair with **no unit bootstrap**.

The inherently per-unit *product* renders -- `show_product`, `show_3d`,
`show_2d` -- target a specific unit's `workspace:<n>` pane, so they are **no-ops
on the global channel**. To use them, open a unit with `new_unit` first, then
drive that unit's `valenx_chat_cmd_u{n}.jsonl` channel.

---

## 3. Command vocabulary

Each line is a flat JSON object internally tagged on `"cmd"`. Fields marked
*(opt)* may be omitted. Source of truth: the `AgentCommand` enum in
`agent_commands.rs`.

| `cmd` | Fields | Effect |
|-------|--------|--------|
| `new_tab` | `name`, `workbench` *(opt)* | Open a new project tab; bind a workbench by id if given (else blank), make it active. |
| `open_workbench` | `id` | Switch the **active** tab's workbench to id. |
| `focus_tab` | `name` | Activate the first tab titled `name`. |
| `rename_tab` | `name` | Rename the active tab. |
| `close_tab` | `name` *(opt)* | Open the "Close tab?" confirm for the named tab (or the active tab). |
| `set_control` | `name`, `value`, `workbench` *(opt)* | Set a labelled control by its **user-visible caption** to a typed value, on the named workbench (else the active tab). |
| `run_command` | `id` | Run a command-palette action by its stable id (e.g. `view.front`, `thermo.compute`). |
| `list_controls` | `workbench` *(opt)* | Post the settable control captions of a workbench to the feed. |
| `list_commands` | -- | Post every runnable `run_command` id to the feed. |
| `read_readout` | `workbench` *(opt)* | Post a workbench's **computed result** text back to the feed. |
| `set_view` | `dir` | Snap the central camera: `front`/`back`/`left`/`right`/`top`/`bottom`/`iso`. |
| `orbit` | `dx_deg`, `dy_deg` | Orbit the central camera by a degree delta (elevation clamped to +-89.9). |
| `zoom` | `factor` | Dolly the central camera (positive = in, negative = out). |
| `frame_all` | -- | Frame the loaded model's bounding box. |
| `add_sketch_point` | `x`, `y` | Add a Line-tool vertex to the in-house CAD sketch. |
| `add_sketch_arc` | `start`, `via`, `end` (each `[x,y]`) | Add a 3-point arc to the CAD sketch. |
| `extrude_sketch` | `height` | Extrude the current CAD sketch profile by `height` (> 0). |
| `add_2d_line` | `x1`, `y1`, `x2`, `y2` | Add a line to the in-house 2-D drawing. |
| `add_2d_circle` | `cx`, `cy`, `r` | Add a circle to the 2-D drawing (`r` > 0). |
| `note` | `text`, `kind` *(opt)* | Post a narration line to the feed (`kind`: `build`/`result`/`ship`/`warn`, default `ship`). |
| `animate` | `n` *(opt)*, `play` *(opt)*, `speed` *(opt)* | Drive an animated product's clock (Play/Pause, speed 0..4). On the global channel pass `n` to pick the unit. |
| `new_unit` | `kind` *(opt)*, `title` *(opt)*, `note` *(opt)*, `group` *(opt)* | **Global-only.** Open a fresh Workbench+Agent unit, optionally building product `kind` and filing the tab into a coloured `group` band. |
| `show_product` | `title`, `lines` *(opt)* | Per-unit only. Render a text result card into `workspace:<n>`. |
| `show_3d` | `kind`, `n` *(opt)* | Per-unit only. Render a live 3-D mesh product into `workspace:<n>`. |
| `show_2d` | `kind`, `n` *(opt)* | Per-unit only. Render a 2-D engineering drawing into `workspace:<n>`. |

### `set_control` values

`value` is untagged JSON: `true`/`false` for toggles, a whole number for integer
controls (`256`), a real for floats (`0.55`), a string for enum-by-name or text
(`"linear"`). A value of the wrong type for the named control posts a `warn`
note and changes nothing -- it never panics. The caption in `name` is the exact
text the user sees next to the widget (the same string the accessibility tree
exposes via `labelled_by`).

### Bad input is always safe

Unknown `cmd` tags, unknown workbench ids, unknown captions, unknown command
ids, and out-of-range values are all **no-ops with a `warn` feed note** -- never
a crash. So an agent can probe and self-correct.

---

## 4. The read-back loop (self-verify without a screenshot)

The point of the bridge is a closed loop: **set -> run -> read**, all over the
feed. A typical sequence on the global channel:

1. `open_workbench thermo` -- bring the tool up on the active tab.
2. `list_controls thermo` -- discover the exact captions you can set.
   *(Reads back, e.g.:* `controls: Fluid, EOS model, Temperature [K], ...`*)*
3. `set_control "Temperature [K]" 350` -- write an input.
4. `list_commands` -- discover the runnable action ids.
   *(Reads back, e.g.:* `commands (NN): view.front, ..., thermo.compute, ...`*)*
5. `run_command thermo.compute` -- fire the solve (acks with a one-line status).
6. `read_readout thermo` -- read the **computed answer** back as a `result`
   note.

Because each step's output lands in the global feed, an agent never needs to
see the screen. `valenx-drive.ps1` prints the new feed lines after every command
so you watch the conversation as you drive.

### Reading the feed yourself

Each feed line is `{"title": "...", "detail": "...", "kind": "..."}`. The
`detail` is the human-readable payload (a readout, an ack like `ran <id>`, or a
`warn` reason). Poll the file for new lines after you append a command;
`read_readout` results carry `kind: "result"`, failures carry `kind: "warn"`.

---

## 5. `valenx-drive.ps1` reference

```
valenx-drive.ps1 <command> [args...]   # shorthand
valenx-drive.ps1 -Raw '<json>'         # arbitrary command object
valenx-drive.ps1 ... -Wait <seconds>   # pause before reading the feed (default 1.5)
valenx-drive.ps1 ... -Tail <n>         # how many feed lines to print (default 8)
```

Shorthands map positional args to fields and emit numbers/booleans as JSON
scalars, strings otherwise:

```powershell
valenx-drive.ps1 set_control "a1 (coeff on x1)" 2.5 uq   # name value [workbench]
valenx-drive.ps1 open_workbench fem
valenx-drive.ps1 run_command view.iso
valenx-drive.ps1 read_readout uq
valenx-drive.ps1 list_controls uq
valenx-drive.ps1 list_commands
valenx-drive.ps1 note "starting the sweep" build
valenx-drive.ps1 new_tab "Hull v2" cad
valenx-drive.ps1 focus_tab "Hull v2"
valenx-drive.ps1 set_view iso
valenx-drive.ps1 frame_all
valenx-drive.ps1 new_unit gear "Reducer"
```

Anything not in the shorthand table (e.g. `orbit`, `zoom`, `add_2d_circle`,
`extrude_sketch`) is sent with `-Raw`:

```powershell
valenx-drive.ps1 -Raw '{"cmd":"orbit","dx_deg":30,"dy_deg":-10}'
valenx-drive.ps1 -Raw '{"cmd":"add_2d_circle","cx":0,"cy":0,"r":25}'
```

The helper resolves the channel paths from the same env vars valenx uses, so if
you launched valenx with `VALENX_ASSISTANT_INBOX`/`VALENX_ASSISTANT_FEED` set,
run the script in a shell with the same values exported and it targets the live
instance automatically.
