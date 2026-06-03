# AGENTS.md — guide for AI coding agents

**If you are an AI agent (Claude, Codex, Cursor, Copilot, Aider, Devin, or
similar) working on Valenx: read this first.** Human contributors should read
[CONTRIBUTING.md](./CONTRIBUTING.md) — this file is the agent-specific
addendum.

Valenx **actively encourages AI-assisted development.** Using AI to move faster
is welcome and expected. The bar is not "was a human or an AI faster" — the bar
is **the change is correct, tested, and lands green.** That responsibility is on
you (and the human who reviews/merges your work), regardless of how it was
written.

---

## The Prime Directive: loop until zero failures

**Do not submit work that does not build and pass tests.** Before you open a PR
or hand a change back, run this loop and keep going until it is fully green:

```sh
cargo fmt --all
cargo build --workspace            # must compile, 0 errors
cargo test  --workspace            # must pass, 0 failed
cargo clippy --workspace --all-targets -- -D warnings   # 0 warnings
```

> **Faster tests:** on this ~150-crate workspace, prefer
> `cargo nextest run --workspace` over the `cargo test` line above — it runs
> every test across all crates in one parallel pool (filling all cores) instead
> of one binary at a time, so it is dramatically faster here. Install once:
> `cargo install cargo-nextest --version 0.9.114 --locked` (0.9.114 is the last
> release supporting the repo's rustc 1.88; newer needs rustc 1.91). Project
> config: `.config/nextest.toml`. The dev profile is already optimized
> (`opt-level = 1`, deps at `3`) so the numerical tests don't run unoptimized.
>
> **Resource-aware — maximize or minimize to fit the machine.** The default
> scales to all cores (same as `cargo test`), fastest on a capable box. On a
> modest machine, or to stay responsive while tests run in the background, use
> `cargo nextest run --profile gentle` (caps parallelism), and pass `-j N` to
> `cargo build` / `cargo nextest run` to throttle compile/run concurrency (e.g.
> `-j 2` on a low-RAM laptop). Nothing here is unbounded — it scales with cores,
> and you can always dial it down.

The loop, explicitly:

1. Make your change.
2. Run **build → test → clippy** (commands above; scope to the crates you
   touched first if the workspace is slow, then widen).
3. If **anything** fails — a compile error, a failing test, a clippy warning —
   **fix it and go back to step 2.**
4. Repeat until there are **zero** failures and warnings.
5. Only then commit / open the PR.

"It mostly works" is not done. "The part I changed works but I didn't run the
suite" is not done. **Green build + green tests + green clippy = done.** If a
pre-existing failure is unrelated to your change, say so explicitly in the PR
rather than silently leaving it.

---

## Non-negotiables (scientific software)

Valenx computes numbers people make real decisions from. Therefore:

- **Wrong numbers are worse than a crash.** If you cannot make a path correct,
  make it **fail loud** with a clear message — never emit a plausible-but-wrong
  result. Several native paths intentionally refuse to run when they are not yet
  calibrated/validated; preserve that behavior.
- **Verify the math, don't trust it.** When you touch a numerical/scientific
  routine, add or update a test that pins the correct value (ideally against a
  published benchmark). A test that encodes a *wrong* expectation is worse than
  no test.
- **Never delete features, crates, or whole workbenches** to make something
  compile or to "simplify." If a deletion seems necessary, stop and ask in the
  PR/issue first.
- **Keep changes focused.** One logical change per PR. Don't bundle unrelated
  refactors. Smaller diffs review faster and break less.

---

## House style

- **Language:** Rust everywhere. C only for unavoidable FFI. No C++/JS/Python
  written in-tree (Python is user scripting only).
- **Native-first:** prefer extending the in-house Rust engines
  (`crates/valenx-*`) over adding a dependency on an external binary. External
  tools belong behind adapter crates (`crates/valenx-adapters/...`), ideally
  with a native fallback.
- **Tests live next to code** (`#[cfg(test)] mod tests`), plus headless
  egui-logic tests for UI panels (draw without panic + input validation).
- **Commits:** conventional style (`feat(scope): …`, `fix(scope): …`,
  `docs: …`). Explain *why* in the body for non-obvious changes.
- **Docs:** update rustdoc comments and any user-facing docs your change
  affects. Don't let docs drift from behavior.

## Getting oriented

1. [ARCHITECTURE.md](./ARCHITECTURE.md) — how the workspace fits together.
2. [CONTRIBUTING.md](./CONTRIBUTING.md) — full contribution process + RFCs.
3. [ROADMAP.md](./ROADMAP.md) — where the project is going; good source of work.
4. Open issues labeled `good first issue` / `help wanted`.

## What we especially want help with

The native engines cover a lot, but these are **not yet in-house** and are
high-value targets (see the README roadmap): a native **electromagnetics**
solver, **parametric CAD** history, native **unstructured meshing**,
**industrial/turbulent CFD**, **DFT** in `valenx-qchem`, and a first-class
**FEA workbench** UI over the existing native `valenx-fem` solvers.

By contributing you agree your work is dual-licensed under MIT OR Apache-2.0.
