# Cross-crate integration tests

Tests here exercise more than one workspace crate together. Per-crate
unit tests live inside each crate's `tests/` or `#[cfg(test)]` modules.

## Layout

- `fixtures/` — sample `.valenx` projects, reference meshes, and other
  inputs used by integration tests. Fixtures are small (< 100 KB each)
  and checked in.
- `workflow_integration.rs` *(to be written in Phase 1)* — end-to-end
  tests that drive the workflow DAG with a stub adapter, verifying
  canonical types flow as designed.
- `project_roundtrip.rs` *(to be written in Phase 1)* — open `.valenx`
  projects from `fixtures/`, write them back, diff.

Running:

```powershell
# Unit + integration tests (the scoped QA harness — never
# `cargo test --workspace`; see ../docs/QA.md for the rationale).
bash scripts/qa.sh           # Unix / WSL / Git-Bash
pwsh scripts/qa.ps1          # Windows PowerShell

# Single integration suite (still scoped — bypasses the file-dialog
# crates that hang a workspace-wide run).
cargo test --test workflow_integration
```

See [TESTING.md](../TESTING.md) for the full testing guide.
