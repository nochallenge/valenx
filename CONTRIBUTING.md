# Contributing to Valenx

Thanks for considering contributing. Valenx is a long-term open-source
project (see [ROADMAP.md](./ROADMAP.md) — 20-year horizon), and everything
gets better with more hands. This doc covers how to get involved.

## TL;DR

1. Read the [ROADMAP.md](./ROADMAP.md) and pick something that interests you
2. For bug fixes / small features: open a PR directly
3. For bigger changes: **file an RFC first** (see [rfcs/README.md](./rfcs/README.md))
4. Be respectful — see [CODE_OF_CONDUCT.md](./CODE_OF_CONDUCT.md)
5. All contributions are dual-licensed under MIT OR Apache-2.0; by contributing you agree

---

## AI-assisted contributions welcome

We **encourage** using AI coding agents (Claude, Copilot, Cursor, Aider, …) to
move faster — a large part of Valenx is built that way. We don't care whether a
human or an AI wrote a line; we care that it is **correct, tested, and lands
green**.

If you're driving an AI agent — or you *are* one — follow **[AGENTS.md](./AGENTS.md)**.
The rule that matters most: **loop build + tests + clippy until they are fully
green before you submit** — never open a PR with a failing build or test. And
because Valenx is scientific software, **wrong numbers are worse than a crash**:
if a path can't be made correct, make it fail loud rather than emit a
plausible-but-wrong result.

---

## Ways to contribute

### Code
- Rust (primary — almost everything)
- C (rarely, for FFI wrappers to OSS libraries)
- No C++, JavaScript, Python written in-tree (Python only as user scripting)

### Documentation
- Inline doc comments (`rustdoc`)
- mdBook user manual under `docs/`
- RFCs under `rfcs/`
- Tutorials + examples

### Validation + physics
- Write regression tests against published benchmarks
- Port existing validation cases from the Python legacy
- Review solver output for physical correctness
- Publish validation papers (credit stays with author)

### UI / UX
- Iconography, themes, interaction design
- Accessibility audits
- Internationalization (i18n) strings

### Infrastructure
- CI / release pipelines
- Installer packaging (MSI, .app, AppImage)
- Docker / container images for headless runs

### Integrations
- Write new adapter crates for OSS tools not yet integrated
  (see the [tool registry](./ROADMAP.md#06-integrated-tool-registry))

---

## Development setup

### Prerequisites

- **Rust stable** (MSRV: current minus 3 releases)
- **Git**
- **VS Code** with rust-analyzer extension (recommended)
- For builds involving OSS tools:
  - OpenFOAM (Linux/WSL — current Valenx-Tools/ setup works)
  - OpenCASCADE dev headers (`sudo apt install libocct-*-dev` or equivalent)
  - gmsh binaries on PATH

### First build

```powershell
git clone https://github.com/<your-org>/valenx.git
cd valenx
cargo build --workspace
# QA via the scoped harness — never `cargo test --workspace`
# (see docs/QA.md for why a workspace-wide cargo test is forbidden).
bash scripts/qa.sh           # Unix / WSL / Git-Bash
pwsh scripts/qa.ps1          # Windows PowerShell
cargo run -p valenx-app
```

The first build compiles all dependencies — takes 5-15 min. Incremental
builds after that are seconds.

### Running tests

```powershell
# Scoped QA harness — preferred path; matches what CI runs.
# NEVER `cargo test --workspace` (see docs/QA.md for the rationale).
bash scripts/qa.sh             # Unix / WSL / Git-Bash
pwsh scripts/qa.ps1            # Windows PowerShell

# Per-crate iteration:
cargo test -p valenx-core         # one crate
cargo bench -p valenx-cfd-native  # benchmarks (criterion)
cargo test -p valenx-cfd-native --release -- --ignored ghia  # validation suite
```

See [TESTING.md](./TESTING.md) for the full testing guide.

### Pre-commit checks

Before each commit, run these locally — the same commands CI runs:

```powershell
cargo fmt --all --check
cargo clippy --workspace --all-targets -- -D warnings
# QA via the scoped harness — never `cargo test --workspace`.
# See docs/QA.md for the rationale (some crates require a display
# or open native file dialogs that would hang non-interactive runs).
bash scripts/qa.sh           # or: pwsh scripts/qa.ps1
```

Catching issues locally saves round-trips through CI. Two ways to make
this automatic:

**Option A — git pre-commit hook, via the `pre-commit` framework**
([pre-commit.com](https://pre-commit.com/)):

```powershell
pip install pre-commit
pre-commit install
```

A `.pre-commit-config.yaml` will be added alongside the workspace
scaffold; until then, use Option B.

**Option B — a shell script in `.git/hooks/pre-commit`**

Two lines: run the three commands above, exit non-zero on failure.
Simplest, zero deps.

Either way, the authoritative check is CI — hooks are for your
convenience.

---

## Contribution workflow

### Small changes (bug fixes, docs, minor features)

1. Fork the repo
2. Create a feature branch: `git checkout -b fix/residual-parser`
3. Make your change
4. Run `cargo fmt && cargo clippy --workspace --all-targets -- -D warnings && bash scripts/qa.sh`
   (NEVER `cargo test --workspace` — see docs/QA.md)
5. Commit with a clear message (see style below)
6. Open a PR against `master`
7. Respond to review feedback
8. A maintainer merges when green

### Large changes (new features, architecture shifts, new tools)

1. **Write an RFC first** — see [rfcs/README.md](./rfcs/README.md)
2. File the RFC as a PR
3. Discuss in the PR comments
4. RFC gets merged when consensus reached (or BDFL approves)
5. Only then start implementation

Reasons to RFC:
- Adding a new crate to the workspace
- Changing the plugin API
- Adding a new integrated OSS tool
- Changing the `.valenx` file format
- Adding / removing dependencies
- Performance-critical algorithm changes
- Breaking API changes

### Commit message style

Follow [Conventional Commits](https://www.conventionalcommits.org/):

```
<type>(<scope>): <short description>

<optional body>

<optional footer>
```

Types: `feat`, `fix`, `docs`, `style`, `refactor`, `perf`, `test`, `chore`

Examples:
```
feat(openfoam): wrap adjointOptimisationFoam
fix(viz): STL loader crashes on binary files > 2GB
docs(rfc): add 0005-adjoint-api
refactor(core): extract workflow DAG into own crate
perf(cfd): vectorize k-omega source term (30% speedup)
```

### PR requirements

- CI green (formatting, clippy, tests, docs build)
- At least one maintainer review approval
- No merge conflicts with `master`
- New features include tests
- Breaking changes include migration notes in the PR description
- New dependencies justified in PR description

---

## Code style

- **Formatting**: `cargo fmt` with default settings
- **Lints**: `cargo clippy -- -D warnings` must pass
- **No unsafe without justification** — `// SAFETY: ...` comment required
- **Prefer explicit over clever** — this is a long-lived codebase
- **Document public API** — every `pub` item needs a doc comment
- **Errors**: `thiserror` for libraries, `anyhow` for app code
- **Logging**: `tracing` macros, not `println!`
- **Serialization**: `serde` with explicit derives; don't rely on struct-order
- **Async**: `tokio` for I/O; `rayon` for CPU-parallel; don't mix
- **Durable writes**: never write durable state with a raw
  `std::fs::write`, `File::create`, or a writing `OpenOptions` chain.
  Route every persisted file (settings, manifests, project files,
  results, subprocess input decks, in-memory-then-write exports) through
  `valenx_core::io_caps::atomic_write_str` / `atomic_write_bytes` (or
  `atomic_write_streaming` to stream a large file without buffering it
  whole). These do sidecar → fsync → atomic-rename, so a crash or a
  concurrent writer can never leave a torn/truncated file on disk.

Full style guide: `docs/src/style.md` (to be written as the codebase
matures).

### The no-raw-fs-write guard

The durable-write rule above is enforced automatically, not by review:

```bash
cargo test -p valenx-core --test no_raw_fs_write_guard
```

This test parses every `crates/**/src/**/*.rs` file with `syn` (a real
Rust AST — `#[cfg(test)]` modules / `#[test]` fns are excluded
structurally, and comments / string literals never produce false
hits) and **fails** if any *production* call raw-writes durable state
via `fs::write`, `File::create` / `File::create_new`, or an
`OpenOptions::new()....open(..)` builder chain. A parse failure also
fails the test loudly — a file the guard can't parse is a file it can't
police.

Legitimate exceptions exist and are listed in the test's `ALLOWLIST`
(matched by `(path_suffix, fn_name)`), each with a one-line reason.
Today's categories:

- **Streaming exports** — STL / OBJ / PLY mesh writers stream
  element-by-element; buffering a multi-GB mesh just to atomic-rename it
  is worse, and a torn export is a re-runnable artifact, not durable
  state.
- **Subprocess stdout/stderr redirect sinks** — the `File` handle is
  moved into `Stdio::from(..)` / `cmd.stdout(..)`; the OS writes the
  child's output, not us.
- **Append-only logs** — e.g. the audit hash-chain, where an
  atomic-rename would clobber the prior chain instead of extending it.
- **The canonical/local atomic-write implementations themselves.**

If you add a genuinely-can't-be-atomic write, add an `ALLOWLIST` entry
with a reason (a future reviewer is then forced to look at it) rather
than weakening the scan. If you're persisting durable state, migrate to
`atomic_write_*` instead — that is almost always the right answer.

---

## Testing expectations

- **Unit tests** for every module. Aim for 80%+ coverage on pure-logic crates.
- **Integration tests** in `tests/` at workspace root.
- **Validation tests** for physics code — compare against published benchmarks.
- **Snapshot tests** for UI — `egui_kittest` or equivalent.
- **Benchmarks** (`criterion`) for hot paths.

PRs that touch physics code must include a validation reference (paper,
reference dataset, known-good value) or the PR won't be accepted.

---

## Governance

- **BDFL** (project lead) has final say through Year 2
- **Technical Steering Committee** forms in Year 2 (elected by contributors)
- **Foundation model** by Year 10 (Apache Software Foundation or similar)
- **RFCs** are how major decisions get made
- **Consensus-seeking** with tie-breaks by BDFL/TSC as needed

See [ROADMAP.md](./ROADMAP.md#6-governance) for the full governance plan.

---

## Becoming a maintainer

Maintainer status is offered to contributors who:
- Have landed 5+ substantial PRs
- Demonstrated sustained engagement (6+ months)
- Reviewed others' PRs constructively
- Follow the code of conduct

Maintainers get:
- Merge rights on specific subsystems
- A vote in RFC decisions
- Listed in MAINTAINERS.md

Invitation is by current maintainer consensus. Email or file an issue to
express interest.

---

## Licensing + CLA

By contributing, you agree your code is dual-licensed under MIT OR
Apache-2.0 (same as the project), without any additional terms or
conditions. No CLA required; the Developer Certificate of Origin
applies (your commits imply DCO agreement).

If you are contributing on behalf of an employer, make sure you have
authorization.

**Do not copy code from proprietary sources.** If you previously saw the
ANSYS or Siemens NX source, don't work on comparable subsystems — stay
clean-room. Use public-domain literature (papers, textbooks) as your
reference.

---

## Getting help

- GitHub Discussions — general questions
- GitHub Issues — bug reports, feature requests
- `#valenx-dev` on Matrix/Discord (set up by Year 1) — chat
- `security@valenx.org` — security reports (see [SECURITY.md](./SECURITY.md))
- Maintainer email for sensitive topics

We're patient with new contributors. Don't worry about "bothering" us —
questions are welcome.
