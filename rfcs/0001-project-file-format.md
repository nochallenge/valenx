# RFC 0001: `.valenx` Project File Format

- **Status:** Accepted (initial design; revisions expected during early implementation)
- **Author(s):** BDFL
- **Created:** 2026-04-21
- **Discussion PR:** (this commit)
- **Tracking issue:** TBD

---

## Summary

Define the on-disk structure of a Valenx project: a directory with a
`.valenx` extension containing a TOML manifest, a human-readable data
folder, binary assets, and results cache. The manifest is
SemVer-versioned independently from the application.

---

## Motivation

Every serious simulation suite has a project file. CAD has STEP and its
own native saves; CFD has case directories; FEA has `.inp` or `.fmu`.
Valenx needs its own — but designed once, carefully, so we don't end up
with the legacy-parser hell that ANSYS and Siemens live with.

Requirements pulled from the vision and roadmap:

1. **Reproducible** — a `.valenx` project opened in 2030 by someone else
   must produce the same results, given the same pinned tools
2. **Version-aware** — the file has its own SemVer; old files open in
   new apps via migration
3. **Diff-friendly** — version-controllable; humans can read the
   manifest
4. **Plays well with binary assets** — CAD, meshes, results can be large
5. **Partially loadable** — opening a project shouldn't require reading
   every results file
6. **Cross-platform** — no platform-specific path formats, line endings,
   or case-sensitivity assumptions
7. **Self-documenting** — someone who's never used Valenx can look at
   the files and understand the project structure

The anti-requirements matter too:

- Not a single monolithic binary blob (that's FreeCAD's FCStd and it
  makes diffs impossible)
- Not raw JSON as the manifest (too verbose for human editing)
- Not tied to Valenx version at load time (old files open forever)

---

## Guide-level explanation

A Valenx project lives in a directory named `<something>.valenx/`. From
the user's POV, they double-click it and the app opens the project;
they see one "thing" in Finder/Explorer. Internally it's a directory:

```
my-project.valenx/
├── project.toml              # manifest — always the entry point
├── tools.lock                # pinned tool versions for reproducibility
├── geometry/                 # CAD, STEP imports, BRep caches
│   ├── airfoil.step
│   └── airfoil.brep
├── mesh/                     # generated meshes, keyed by config hash
│   └── a3f92c/
│       ├── mesh.msh
│       └── metadata.toml
├── cases/                    # physics cases
│   ├── cfd-steady/
│   │   ├── case.toml         # this case's config (relative to manifest)
│   │   ├── inputs/           # things the user authored (BCs, initial fields)
│   │   └── results/          # things the solver wrote (can be deleted)
│   └── thermal/
│       └── case.toml
├── scripts/                  # user Python/Lua scripts
│   └── postprocess.py
├── thumbnails/               # UI cache — screenshots of viewport
│   └── scene.png
└── notes.md                  # optional user notes
```

### The manifest

`project.toml`:

```toml
[project]
format = "1.0"                            # RFC 0001 version
name = "airfoil-study"
valenx_min = "0.8.0"                      # minimum app version that can open this
created = "2026-04-21T09:30:00Z"
modified = "2026-04-22T17:12:44Z"
author = "Jane Doe <jane@example.com>"
description = "NACA 0012, Re=6e6, sweep of AoA from 0 to 15"

[units]
length = "m"
mass = "kg"
time = "s"
temperature = "K"

[geometry]
entries = [
    { id = "airfoil", source = "geometry/airfoil.step", format = "step-ap242" }
]

[mesh.default]
source = "mesh/a3f92c/mesh.msh"
config_hash = "a3f92c..."                 # so we know when the mesh is stale

[cases]
order = ["cfd-steady", "thermal"]         # UI display order
```

### Tool pinning

`tools.lock` (TOML, but machine-maintained — users don't hand-edit it):

```toml
format = "1.0"
generated_by = "valenx 0.8.1 on 2026-04-21T09:30:00Z"

[[tool]]
name = "openfoam"
version = "v2406"
checksum = "sha256:abc123..."
channel = "stable"

[[tool]]
name = "gmsh"
version = "4.12.2"
checksum = "sha256:def456..."
channel = "stable"
```

Opening a project pins those tools for the lifetime of the session.
Re-solving with a different tool version produces a warning; user can
explicitly update the lock.

### Case files

Each case is a subdirectory with its own `case.toml`:

```toml
[case]
format = "1.0"
name = "cfd-steady"
physics = "cfd"
solver = "openfoam.simpleFoam"
mesh = "default"                         # references [mesh.default] in project.toml

[flow]
turbulence = "kOmegaSST"
schemes = "upwind-first-order"
fluid = { name = "air", rho = 1.225, nu = 1.5e-5 }

[boundaries.inlet]
type = "velocity-inlet"
velocity = [50, 0, 0]
turbulence_intensity = 0.05

[boundaries.outlet]
type = "pressure-outlet"
pressure = 0

[boundaries.walls]
type = "no-slip"

[solve]
iterations = 2000
residual_target = 1e-5
```

The schema for `[flow]`, `[boundaries]`, `[solve]` is specific per
physics type and versioned by the `format` in the case file. Adapters
know how to translate from this canonical form into the solver's
native input deck.

---

## Reference-level explanation

### File naming

- A project is any directory ending in `.valenx`
- Platform-neutral: case-sensitive, UTF-8 paths, forward-slash
  internally even if the host OS uses backslashes
- No reserved filenames inside the project (i.e., a case can be called
  `AUX` on Windows because we validate on save)
- Max path length inside the project: 200 characters (to stay within
  Windows MAX_PATH when the project is placed in a deep user folder)

### Manifest parsing rules

- `project.toml` is **the** entry point. No other file in the project
  can be loaded until the manifest is parsed.
- `format` in `[project]` is mandatory. Everything else has defaults.
- Unknown keys are **preserved on load and written back on save**. This
  is the forward-compat escape hatch: a newer app writes a key, an
  older app doesn't know what it means but doesn't delete it.
- Invalid TOML → fail to load with a precise error (line, column, what
  was expected)

### Version handling

On load:

| File format | App format | Action |
|-------------|-----------|--------|
| Equal | Equal | Normal load |
| File minor < App minor (same major) | Same major | Load, apply forward-compat defaults for missing fields |
| File minor > App minor (same major) | Same major | Load with warning; unknown fields preserved |
| File major < App major | Newer app | Auto-migrate on open, write-back on next save |
| File major > App major | Older app | Refuse to load; prompt user to upgrade Valenx |

Migration paths are code in `valenx-core::project::migrate::vN_to_vM`.
Migrations are always idempotent and always reversible (we emit a
"previous format" file alongside on first migration).

### Relative vs absolute paths

**All paths inside the project are relative to `project.toml`.**
Never absolute. Never `~/`. Never OS-specific. This is a hard rule —
the project must be movable by zipping the directory and unzipping it
elsewhere.

Paths **outside** the project (tool binaries from `tools.lock`, system
libraries) live in `tools.lock` and in Valenx's user settings, not in
the project.

### Mesh caching and invalidation

Each mesh gets a directory under `mesh/` named after its
**configuration hash**: a SHA-256 over the mesh-generation inputs
(source geometry file hash + mesh sizing config + mesher tool version).

When the user re-runs mesh generation:

- If the hash matches an existing dir, we reuse
- If not, generate a new dir; don't delete the old (user may want to
  compare)

Garbage collection is an explicit user action: *File → Clean
unreferenced caches.*

### Results layout

`cases/<name>/results/` is entirely "owned" by the solver. Adapters
write whatever structure their underlying tool produces; the adapter
also writes a `cases/<name>/results/manifest.toml` summarizing:

- Solver used + version
- When it ran + how long
- Exit code
- Residual / convergence history (sampled)
- Pointers to result files (`fields.vtk`, `forces.csv`, etc.)
- Hash of inputs that produced this result

Results can be safely deleted — the case re-runs from `inputs/`.

### Thumbnails and UI state

`thumbnails/` is a cache; anything in there may be regenerated. UI
state that's **user intent** (which panels are open, viewport camera,
selected objects) lives in `project.toml` under a `[ui]` section.

### Binary assets

Large binary files (mesh, STEP, results) are stored as-is; not
serialized into the manifest. We do not use git-lfs-style pointers
because a Valenx project is expected to be self-contained.

For version control, users should `.gitignore` the `results/` and
`thumbnails/` directories.

### Schema validation

Shipped with the app: a JSON Schema for the project and case TOML
structures (auto-generated from the Rust types via `schemars`).
Surfaced in the UI as real-time validation; used in CI for test
fixtures.

---

## Drawbacks

- **Directory-as-file is awkward on some platforms.** macOS handles it
  well (bundle behavior); Windows shows a folder unless we use a shell
  extension; Linux shows a folder. We can't hide this without shipping
  per-platform shell extensions.
- **More complex than a single-file blob.** We accept this for diff
  friendliness; power users can still `tar` a project to share.
- **Paths are a perennial source of bugs.** Relative-paths-only is a
  discipline; it requires all code paths to go through a single
  `ProjectPath` wrapper type.
- **TOML has no schema language of its own.** We use `serde` + JSON
  Schema for validation, which is fine but means two sources of truth.

---

## Rationale and alternatives

**Single-file ZIP (like FCStd, .xlsx):**
Rejected. Diffs impossible, partial load impossible, binary-blob
corruption risk.

**Database-backed (SQLite):**
Rejected. Solves the "many-files" problem but loses diff-friendliness;
schema migration is harder than file migration; opaque to scripting.

**JSON instead of TOML:**
Rejected. Humans edit the manifest sometimes (advanced users); TOML is
easier on the eyes. Cost: TOML lacks some features (no real references,
no null), but we don't need them at this layer.

**Stick with what OpenFOAM case directories look like:**
Rejected as the top-level. Adapters still use OpenFOAM-native layouts
inside `cases/<name>/inputs/` — the canonical layer above is ours, but
we don't dictate what the adapter writes below.

**Don't version the file format at all, tie to app version:**
Rejected. This is how FreeCAD got stuck — the file format drifts
silently and old files break mysteriously. Explicit SemVer forces us
to think about back-compat on every change.

---

## Prior art

- **STEP (AP242)** — canonical for CAD interchange; we import from it
  but don't use it as our own project format (too narrow, no simulation)
- **CGNS** — CFD community standard for results; we use it internally
  but it's results-only
- **SimScale** project format — proprietary, web-hosted, not applicable
- **ANSYS `.wbpj`** — binary, version-locked, a cautionary tale
- **Rhino `.3dm`** — the "file is a directory" pattern, works fine on macOS
- **Xcode `.xcodeproj`** — also a directory bundle; has taught users
  the pattern
- **Rust crates / Cargo projects** — directory with `Cargo.toml` as
  entry point; `Cargo.lock` for reproducibility. The layout of our
  project file is **directly inspired** by Cargo.

---

## Unresolved questions

- **Multi-user / collaborative projects.** Out of scope for v1.0. We'll
  need another RFC for this — probably merge-conflict resolution on
  the manifest + CRDTs for UI state.
- **Versioning the case schema alongside the project schema.** Current
  plan: case `format = "1.0"` is independent from project `format =
  "1.0"`. That has upsides (cases evolve per physics) and downsides
  (more versions to track). Revisit after 6 months of real use.
- **Encryption / sensitive data.** Some industrial users will want
  encrypted project files. Defer to a separate RFC.
- **Partial-network projects** (cases stored on NAS, some geometry in
  cloud). Defer.

---

## Future possibilities

- A `.valenx-template` format — a manifest with placeholders, useful
  for case templates ("NACA airfoil study")
- A `.valenx-bundle` format — a zipped, signed, read-only archive for
  publishing validation cases
- Schema-driven UI form generation — the case editor's panels could be
  auto-generated from the TOML schema
- Full-text search inside a project — indexable thanks to TOML +
  markdown notes being plain text
