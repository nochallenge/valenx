# Project Policies

Binding commitments for how Valenx evolves. These apply once 1.0 ships;
during pre-alpha (now through ~Year 1) everything is subject to change,
and breaking changes are announced in the CHANGELOG but not gated by
the deprecation cycles below.

---

## Versioning

Valenx follows **Semantic Versioning 2.0.0** for public APIs and file
formats:

- **MAJOR** — incompatible API or file-format changes
- **MINOR** — additive, backward-compatible features
- **PATCH** — backward-compatible bug fixes

Additionally:

- **Pre-1.0** uses `0.MINOR.PATCH` with `0.MINOR` bumps for breaks
- **Nightly** builds use `X.Y.Z-nightly.YYYYMMDD`
- **LTS** tags use `X.Y.Z-lts` and are marked in the release notes

### What "public API" covers

SemVer applies to:

- **Rust crate APIs** exported by `valenx-*` crates tagged `pub` in `lib.rs`
- **Plugin ABI** (WIT interfaces exposed to WASM plugins)
- **Python/Lua scripting API** (`valenx.scripting` module)
- **`.valenx` project file format**
- **Adapter contract** between `valenx-core` and `valenx-adapters-*`
- **CLI flags** on the `valenx` binary

SemVer does **not** apply to:

- Internal modules marked `pub(crate)`
- Debug-only features behind `#[cfg(debug_assertions)]`
- Experimental items marked `#[doc(hidden)]` or gated by a feature flag
  prefixed with `unstable-`
- Log output format (parse the structured fields, not the text)
- File layout inside `~/.valenx/` (use the API, not direct access)

---

## Deprecation policy

We deprecate rather than delete. Once 1.0 ships:

1. **Announce** the deprecation in a minor release with:
   - `#[deprecated]` attribute on the Rust item
   - Note in the CHANGELOG under `### Deprecated`
   - Migration path documented
2. **Maintain** the deprecated item for at least **two minor versions**
   (at our quarterly minor cadence, this is roughly six months
   minimum — often longer)
3. **Remove** in the next major release, with a migration note in the
   upgrade guide

For file formats, the deprecation window is longer: **one full major
version**. A `.valenx` file written by 2.x must still load in 3.x (it
may be upgraded on load).

---

## LTS (Long-Term Support)

Once we reach 1.0, we will designate an LTS release roughly every
**18 months**. LTS releases:

- Receive **security fixes for 24 months** from tag date
- Receive **critical bug fixes for 12 months**
- Get **no new features** — stability over novelty
- Track a pinned `tools.lock` (third-party solver versions frozen)

The LTS commitment is what institutional users need to justify adoption:
they can standardize on an LTS, validate their workflows against it, and
know it won't shift under them.

---

## Release cadence

| Track | Cadence | Audience |
|-------|---------|----------|
| Nightly | daily (from `master`) | developers, early adopters |
| Stable | **quarterly minor** / patch as needed | everyday users |
| LTS | every **18 months** | institutional users, teaching |

Release dates are published on the website and in the mailing list.

---

## Supported versions

Security and bug fixes are backported according to this matrix:

| Version | Security fixes | Bug fixes | Notes |
|---------|----------------|-----------|-------|
| Current stable (`N`) | Yes — all | Yes — all | Patched in next point release |
| Previous stable (`N-1`) | Yes — critical only | No | For ~6 months after `N` ships |
| Current LTS | Yes — all | Yes — critical only | 24 months from LTS tag |
| Nightly / dev | Rolling — fixed on `master` | Rolling | No backports |
| Anything older | No | No | Upgrade path always documented |

Until 1.0 ships, `master` is the only supported branch; apply fixes by
pulling the latest commit. Vulnerability-reporting mechanics (where to
send reports, response SLAs, coordinated disclosure) live in
[SECURITY.md](./SECURITY.md).

---

## Backward compatibility for project files

`.valenx` project files are SemVer-versioned independently (`format = "1.2"`
in the file header). Load behavior:

| File version | Current app | Behavior |
|--------------|-------------|----------|
| Same major | Any | Load, possibly with forward compat warnings |
| Older major | Newer app | Load with auto-migration; write back in new format if edited |
| Newer major | Older app | Refuse to load, point user to upgrade |

A file written today on version 1.2 must still open in version 1.x
forever, and in version 2.x after one-time migration.

---

## Plugin ABI stability

Plugins compiled for version `X.Y` are expected to work on any version
`X.Z` where `Z >= Y`. A major version bump may require recompilation,
but the **WIT interface** changes follow the deprecation cycle above.

If a plugin becomes incompatible, the app shows a clear error pointing
at the specific WIT interface that changed and where to find the
migration doc.

---

## Supported platforms

Tier 1 — supported targets the project commits to keeping green:

- Windows 10+ x86_64
- macOS 12+ x86_64 and arm64
- Ubuntu 22.04+ x86_64

Tier 2 — best-effort:

- Other Linux (Fedora, Arch, openSUSE)
- Linux arm64
- Windows 11 arm64

Tier 3 — community-supported:

- FreeBSD
- Other Unix-likes

**CI status (2026-05-24):** the project's GitHub Actions workflows
(`ci.yml`, `ci-nightly.yml`, `release.yml`) are currently
**`workflow_dispatch`-only** — they are not auto-triggered on
push / PR / tag. The tiers above describe the targets the
maintainers manually verify on each release, not a CI-blocking
guarantee. See `docs/CI.md` for the rationale (workspace size +
flaky cross-platform display dependencies during the
pre-alpha period). The intent is to re-enable push-triggered CI
once the workspace stabilises; until then, `bash scripts/qa.sh`
is the local source of truth.

We don't support 32-bit platforms. We don't support Windows 7/8 or macOS
< 12.

---

## Minimum Supported Rust Version (MSRV)

The MSRV is the oldest stable Rust toolchain we commit to building on.

- **MSRV floor is `stable - 3`** — at any given time, our `rust-version`
  sits three six-week cycles (~18 weeks) behind the current stable
  Rust release.
- The current MSRV is the **specific version** named in both
  [`rust-toolchain.toml`](./rust-toolchain.toml) (`[toolchain] channel`)
  and the workspace [`Cargo.toml`](./Cargo.toml) (`rust-version`).
  Those two values must match — a mismatch fails CI.
- An **MSRV bump is a minor-version bump** for our published crates, and
  is called out in the CHANGELOG under `### Changed`. It is not a major
  bump, because raising a floor doesn't break callers on newer
  toolchains — but downstreams pinned to an older Rust need to know.
- MSRV bumps are justified in the PR description: a feature stabilized,
  a dep now requires it, a soundness fix, etc. "Nice to have" doesn't
  qualify on its own.
- Once 1.0 ships, an LTS release's MSRV is frozen at LTS-tag time and
  does not bump for the life of that LTS.

We don't promise to work on every Rust version forever — that would
freeze the language out of useful features indefinitely. The
`stable - 3` window is long enough that distro packagers and
conservative shops have time to catch up.

---

## Dependency policy

### Adding a new dependency

New direct dependencies require review and justification in the PR:

- Why this crate over alternatives?
- Who maintains it? How active?
- Licence compatible with Apache 2.0? (see below)
- Does it add substantial binary size / compile time?
- Does it pull in C/C++ code? (if yes, extra scrutiny)

We prefer:

- Maintained crates with > 1M downloads where possible
- Crates from rust-lang, tokio-rs, dtolnay, servo, or other vetted sources
- Pure-Rust over C-binding when performance is comparable

### Licence compatibility for dependencies

| Licence | Can use? |
|---------|----------|
| Apache 2.0, MIT, BSD, ISC, Unlicense, 0BSD | Yes |
| MPL 2.0 | Yes (file-level copyleft) |
| LGPL 2.1 / 3.0 | Yes for dynamic linking only |
| GPL 2.0 / 3.0 | No, except for **subprocess-only** use |
| AGPL | No |
| Commercial / proprietary | No |

GPL solvers (OpenFOAM, gmsh, Code_Aster, etc.) are isolated via the
subprocess firewall — we never link against them, we exec them. This is
documented per tool in the adapter crates.

---

## Code ownership and copyright

Valenx is copyright its contributors, collectively. The `LICENSE-MIT`
and `LICENSE-APACHE` files at the root govern use — Valenx is
dual-licensed under MIT OR Apache-2.0, and downstream consumers may
choose either license. Individual file headers are optional; if
present, they must not contradict the dual MIT/Apache-2.0 license.

No CLA is required. Contributions are accepted under the Developer
Certificate of Origin — by committing to this repo you certify that
your contribution is yours to give (or you have the rights to give it).

---

## Changes to these policies

Policy changes go through the RFC process (see [rfcs/README.md](./rfcs/README.md)).

A policy change never applies retroactively — if we tighten the
deprecation window, anything already deprecated stays on the old
schedule.
