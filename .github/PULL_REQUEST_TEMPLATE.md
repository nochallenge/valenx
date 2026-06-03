## Summary

What does this PR do? 1-3 sentences.

## Related

- Fixes # (issue)
- Implements RFC # (if applicable)
- See also: …

## Type of change

- [ ] Bug fix (non-breaking)
- [ ] New feature (non-breaking)
- [ ] Breaking change
- [ ] Documentation only
- [ ] Refactor (no behavior change)
- [ ] CI / tooling / build
- [ ] RFC

## Checklist

- [ ] `cargo fmt --all --check` passes
- [ ] `cargo clippy --workspace --all-targets -- -D warnings` passes
- [ ] `bash scripts/qa.sh` (or `pwsh scripts/qa.ps1`) passes — the scoped QA harness, NEVER `cargo test --workspace` (see `docs/QA.md` for why)
- [ ] Added / updated tests for changed behavior
- [ ] Added / updated docs (rustdoc, mdBook, or CHANGELOG as appropriate)
- [ ] For new deps: justified in PR description, license compatible
- [ ] For breaking changes: migration notes included below

## UI change checklist (only for PRs that change user-facing UI)

- [ ] Mockup link or screenshot attached
- [ ] Snapshot test covering the new / changed state
- [ ] Only tokens used — no hard-coded colors, spacings, motion, shadows
- [ ] Keyboard-reachable; focus order reasonable
- [ ] All strings use a localization key (no hard-coded English)
- [ ] References the relevant design principle when making a contested choice
  (see [DESIGN_PRINCIPLES.md](../DESIGN_PRINCIPLES.md))

## Test plan

How did you verify this works? Commands run, cases exercised,
screenshots, etc.

## Migration notes (if breaking)

What do existing users / contributors need to do?
