# `vendor/vtkio` — local patches vs upstream `vtkio = "0.6.3"`

This is a vendored copy of [`vtkio`](https://crates.io/crates/vtkio)
0.6.3 with workspace-local patches. The patched version number is
`0.6.3+valenx-patch.1`. The `+valenx-patch.N` build-metadata suffix
keeps cargo-deny / SPDX tooling happy while letting us bump the suffix
on each subsequent local change.

## Patches applied

### 1. `lz4_flex` bump 0.7 → 0.11

The upstream `vtkio = "0.6.3"` `Cargo.toml` pinned `lz4 = "0.7"` (where
`lz4` is actually `package = "lz4_flex"`, the pure-Rust LZ4 implementation
maintained by `PSeitz`). That pin pulled `lz4_flex 0.7.x` into our
workspace.

`lz4_flex 0.7.x` was flagged by [`RUSTSEC-2026-0041`][rustsec] — an
out-of-bounds read in the decompressor when fed maliciously-crafted
input. The fix landed in `lz4_flex 0.11`, which is otherwise wire-
compatible with 0.7 for the subset of the API `vtkio` uses
(`decompress_into`, `compress_prepend_size`).

This vendored copy pins `lz4_flex = "0.11"` and rebuilds — no source
changes inside `src/` were required.

### 2. `[lints]` table

Added a [lints] table to Cargo.toml to silence lint warnings in this vendored crate (no source changes).

### Verifying the zero-source-diff invariant

To prove this vendored copy has **no source changes** vs upstream
0.6.3 — i.e. the patch is purely a `Cargo.toml`-level dep bump — run:

```bash
# From the repo root. The grep filter is conservative: only Cargo.toml
# differs (the dep bump + the +valenx-patch.N version suffix).
git diff --no-index --stat \
    vendor/vtkio/src \
    <(curl -sL https://crates.io/api/v1/crates/vtkio/0.6.3/download | \
      tar xz -C /tmp && echo /tmp/vtkio-0.6.3/src)
# Expected output: empty (no files differ).
```

If `src/` ever diverges from upstream, a separate `### N.` entry above
must document each functional patch, and the `+valenx-patch.N` suffix
in `vendor/vtkio/Cargo.toml` must bump.

[rustsec]: https://rustsec.org/advisories/RUSTSEC-2026-0041.html

## How to update

When upstream `vtkio` cuts a release that picks up the `lz4_flex 0.11`
bump on its own (or whatever future RustSec advisory lands next), drop
this vendor copy entirely and switch back to the crates.io version in
the workspace `Cargo.toml`. Until then, every patch we make here gets
a new entry above and a bumped `+valenx-patch.N` suffix in
`vendor/vtkio/Cargo.toml`'s `version` field.
