# valenx — salvaged work (2026-06-22)

Rescued during a disk cleanup that deleted 135 stale agent worktrees under
`valenx/.claude/worktrees/`. Those worktrees were throwaway checkouts from
completed workflow runs, but they still held work that had **never been
committed or merged anywhere**. This folder is that work. Total: 1.3 MB.

## `unmerged-crates/` — 13 crates, 11,342 lines of Rust

These crates existed **only** as untracked files inside the deleted worktrees.
They are not in `master`, not on GitHub, not on any branch. Verified missing
from the workspace before salvaging.

| Crate | Lines | Crate | Lines |
|---|---|---|---|
| valenx-buoyancy | 1069 | valenx-radiation | 1021 |
| valenx-doseresponse | 1039 | valenx-cablesag | 994 |
| valenx-hardyweinberg | 1003 | valenx-epicyclic | 964 |
| valenx-rlcresonance | 946 | valenx-nozzle | 830 |
| valenx-curvedbeam | 743 | valenx-fin | 730 |
| valenx-gyroscope | 722 | valenx-thevenin | 674 |
| valenx-thickcylinder | 607 | | |

Each is a normal crate (`Cargo.toml` + `src/`). To revive one: copy it into
`crates/`, add `"crates/<name>"` to the workspace `members` in the root
`Cargo.toml`, then `cargo test -p <name>`. They were never reviewed or CI-tested,
so treat them as drafts.

## `patches/` — 134 `git diff` files

One per dirty worktree, capturing every **tracked-file** modification (mostly
workspace `Cargo.toml`/`Cargo.lock` member additions, plus ~35 modified
`crates/valenx-app/src/*.rs` files that wired new crates into the app).
Apply with `git apply <file>.patch` from a repo root if ever needed.

## `MANIFEST.txt`

Full `git status --porcelain` of all 135 worktrees at deletion time — the record
of exactly what was there.

## Not included (deliberately)

114 other untracked crates found in the worktrees were **already merged into
`master`**, so copies were redundant and were not kept. Build artifacts
(`target/`) were deleted as pure regenerable cache.
