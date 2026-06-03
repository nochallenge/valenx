# Test fixtures

Minimal, checked-in inputs used by cross-crate integration tests.
Keep them small — bigger reference data belongs out-of-tree.

## Contents

### `minimal.valenx/`

The smallest possible `.valenx` directory that conforms to RFC 0001.
Used by:

- **project-file load / save round-trip** (in-progress; lands with
  `valenx-core::project`)
- **tools.lock parsing**
- **ASCII STL import smoke test**

Structure:

```
minimal.valenx/
├── project.toml          # manifest
├── tools.lock            # pinned tool versions
├── geometry/
│   └── box.stl           # 2-triangle ASCII STL
└── cases/
    └── cfd-steady/
        └── case.toml     # minimal CFD case definition
```

## Conventions

- No results/thumbnails/caches checked in — tests that need outputs
  produce them at runtime in a temp dir.
- Every fixture file under 10 KB. Larger fixtures go in CI caches
  or `git-lfs` once we need them (not today).
- Every fixture is linked from at least one test; orphaned fixtures
  get removed in a cleanup PR.
