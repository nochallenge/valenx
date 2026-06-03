# How to Test the Native Desktop App

You keep asking "how do we test / view / read it during development" —
this doc answers exactly that. No web browsers involved. No localhost.
No dev server. No preview tab.

> **Critical safety note.** A blanket `cargo test --workspace` is
> **forbidden** in this repo — some crates open native file dialogs or
> require a display, and a workspace-wide run hangs the CI agent
> indefinitely. See [docs/QA.md](./docs/QA.md) for the rationale.
> Use the scoped QA harness instead:
>
> ```powershell
> bash scripts/qa.sh          # Unix / WSL / Git-Bash
> pwsh scripts/qa.ps1         # Windows PowerShell
> ```
>
> Use `cargo run -p <crate>` rather than `cargo run` to launch the app
> so the workspace's binary is unambiguous.

---

## The short answer

```powershell
cargo run -p valenx-app
```

That single command compiles the native app and **opens a window on your
desktop** with Valenx running. You interact with it like any other
desktop app — click buttons, drag panels, resize the window, use
keyboard shortcuts.

To close it, click the X like any other window.

---

## The development loop

### For every code change during development

```powershell
# Fast path — debug build, compiles in ~2-10 sec after the first build
cargo run -p valenx-app

# Release build — slower compile, runs at native speed
cargo run -p valenx-app --release
```

The first build takes a few minutes because it compiles all the
dependencies. After that, `cargo` only recompiles what changed, and
launches are <10 seconds.

### For watching changes live

```powershell
# Install once
cargo install cargo-watch

# Then run this — every time you save a .rs file, it rebuilds and relaunches
cargo watch -x "run -p valenx-app"
```

This is the closest thing to "hot reload" you'll have. It's not as
instant as a web page auto-refresh, but it's fast enough that you can
iterate on UI changes in seconds.

### For just one crate

The workspace has many crates (`valenx-app`, `valenx-core`, `valenx-viz`, ...).
You can work on one in isolation:

```powershell
cargo run -p valenx-app           # run just the app
cargo test -p valenx-cfd-native   # test just the CFD crate
cargo check -p valenx-geo         # type-check the CAD crate (fastest)
```

---

## How to "view" / "read" it

You see the app the same way you see FreeCAD: as a window on your
desktop. There's no web browser, no localhost URL, no preview pane.

If you're sharing a screenshot:
- Windows: `Win+Shift+S` then paste, or `Print Screen`
- macOS: `Cmd+Shift+4` for region, `Cmd+Shift+3` for full screen
- Linux: `gnome-screenshot -a` or `spectacle` or use system tools

If you want to record a demo:
- OBS Studio is the open-source standard
- Windows: built-in `Win+G` Game Bar for quick recordings
- macOS: `Cmd+Shift+5` built-in screen recorder

---

## Automated tests (what runs in CI)

### Unit tests — via the scoped QA harness

```powershell
# Scoped QA harness — preferred path; runs the per-crate tests CI runs.
bash scripts/qa.sh                # Unix / WSL / Git-Bash
pwsh scripts/qa.ps1               # Windows PowerShell

# Per-crate (when you want to iterate on one crate without the full run):
cargo test -p valenx-core         # one crate
cargo test -p valenx-core --release   # at release optimization
```

NEVER run `cargo test --workspace` — some crates require a display or
open native file dialogs that would hang non-interactive runs. The
harness explicitly excludes those.

These are pure-Rust tests — no window opens, runs in a terminal, takes
seconds.

### Integration tests

In `tests/` at the workspace root:

```powershell
cargo test --test workflow_integration
```

These test that:
- Opening a project file works
- Running a simple OpenFOAM adapter case produces expected output
- Workflow DAG executes correctly
- File I/O round-trips (.valenx save/load)

They don't open a window — they drive the core library directly.

### UI tests (visual regression)

Tougher but doable. Two options:

**Option A — snapshot testing with `egui_kittest`:**

```rust
#[test]
fn case_tree_renders_correctly() {
    let harness = Harness::new(|ui| render_case_tree(ui, &sample_case()));
    harness.snapshot("case_tree");     // saves reference image on first run
}
```

Subsequent test runs compare against the saved image and fail if pixels
differ beyond a threshold.

**Option B — headless screenshot-based tests:**

Launch the app in a headless window manager (xvfb on Linux, or hidden
window on Win/Mac), take a screenshot, compare with reference.

### Benchmark tests

For the solvers we write natively:

```powershell
cargo bench -p valenx-cfd-native
```

Uses `criterion.rs` to measure performance over time. Catches regressions
where a change makes something 2x slower.

### Validation tests (physics correctness)

`crates/valenx-cfd-native/src/benchmark.rs` (and sibling
`crates/valenx-fem/`) host the canonical validation benchmarks:
- Ghia cavity (lid-driven, Re=100, 1000, 3200) — in `valenx-cfd-native`
- Driver-Seegmiller backward-facing step
- NAFEMS linear elasticity suite — in `valenx-fem`
- GRI-Mech 3.0 adiabatic flame temperature

These run our native solvers and compare against published reference
data. They take minutes to hours; normally run nightly in CI.

```powershell
cargo test -p valenx-cfd-native --release -- --test-threads=1 ghia
cargo test -p valenx-fem --release -- --test-threads=1 nafems
```

---

## Live debugging

### Print debugging

Rust has `println!` and `dbg!` — work the same as anywhere:

```rust
dbg!(&case.parameters);    // prints filename, line, and value to terminal
println!("u_max = {}", u_max);
```

Output goes to the terminal you ran `cargo run -p valenx-app` in.

### Real debugger

Rust has first-class debugger support via LLDB (macOS) or CodeLLDB
extension in VSCode/JetBrains.

Set breakpoints in the IDE, press F5, app launches under the debugger.
Inspect variables, step through code, all the usual.

### Rust analyzer

Every serious Rust dev uses `rust-analyzer` (an IDE plugin). Gives you:
- Inline error messages (no need to compile to see type errors)
- Hover for type info
- Jump to definition
- Rename refactoring
- Auto-imports

Install: `rustup component add rust-analyzer`, then add the extension to
your editor.

---

## Logging + tracing

The native app uses `tracing` for structured logs:

```rust
tracing::info!(case_id = %case.id, "starting solve");
tracing::error!(error = ?e, "solve failed");
```

By default logs print to the terminal. In production, they go to a log
file at `~/.valenx/logs/valenx-YYYY-MM-DD.log`.

Turn up verbosity:

```powershell
$env:RUST_LOG = "valenx=debug,valenx_cfd=trace"
cargo run -p valenx-app
```

---

## What I personally do when developing

```powershell
# Terminal 1 — watch + auto-rebuild
cargo watch -x "run -p valenx-app"

# Terminal 2 — run the QA harness on file change (NEVER cargo test --workspace)
cargo watch -s "bash scripts/qa.sh"

# Terminal 3 — tail logs (production only)
Get-Content "$env:USERPROFILE\.valenx\logs\valenx.log" -Wait
```

Three terminals open, save a file, see the window update + tests pass.
Tight feedback loop, no browser anywhere.

---

## "But how do YOU see the app during my development session?"

Honest answer: I can't directly see a native window running on your
machine the way I can navigate a browser preview. What I can do:

1. **You describe what you see** — this is the most reliable. Screenshot + description.
2. **Automated snapshot tests** — I add `egui_kittest` tests that save
   reference images; you paste the saved image and I can reason about
   what it looks like.
3. **Screenshot via `screencapture` / `gnome-screenshot` / Windows `nircmd`** — I can script a capture and you share the file.
4. **Log the UI state to JSON** — for debugging layout issues, I can
   have the app dump its widget tree to stdout; I read that.
5. **Remote debugging via `tracing` over a socket** — optional, if
   you want me to be able to inspect UI state live.

For day-to-day development I think option 1 is most productive; for
reproducible UI testing, option 2 is best.
