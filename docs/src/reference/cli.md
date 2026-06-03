# CLI reference

Valenx ships a suite of headless command-line tools so projects,
runs, and reports can be driven from CI pipelines without booting
the GUI. Every tool follows the same conventions:

- `-h` / `--help` prints usage and exits 0.
- `-V` / `--version` prints `<name> v<CARGO_PKG_VERSION>` and exits 0.
- `--format json` (where applicable) emits a single parseable JSON
  document for downstream tooling.
- Exit codes are pinned at four levels: `0` success, `1` content
  error (parse / structural), `2` usage error, `3` I/O error. A few
  tools add a domain-specific code (e.g. `4` for a quality-gate
  failure in `valenx-mesh-info`) — those are documented per-tool.

All six binaries live under `crates/<crate>/src/bin/` and are
buildable via `cargo run --bin <name>` without extra setup. They
share a tiny argument-parsing convention (one `match` over the
arg list) — there's no `clap` dep, so each binary stays small.

The three inspectors that consume a JSON file
(`valenx-mesh-info` / `valenx-results` / `valenx-report`) accept
`-` as the path argument to read from stdin per Unix convention:

```bash
cat workdir/results.json | valenx-results -
kubectl get pod ... -o json | valenx-results - --format json
gmsh -3 case.geo -o - | <pipeline> | valenx-mesh-info - --check max-skew=0.9
```

## `valenx-init`

Scaffold a fresh `.valenx` project from a template.

```bash
valenx-init my-cfd-case --template cfd
valenx-init experiments/run-42 --template fea --name 'cantilever-beam'
valenx-init --list-templates       # discovery — prints the catalogue
```

Templates ship with a runnable `case.toml` and (where the adapter
expects a text input) a sample input file alongside. After a
successful scaffold the tool prints two follow-up commands the
user can run immediately (`valenx-validate <dir>` and `$EDITOR
<cases>`).

**Templates** (canonical name; aliases in parens; default case dir):

| Template | Description | Case dir |
|---|---|---|
| `empty` | minimal skeleton — no per-physics block | `case-1` |
| `cfd` (`openfoam`) | OpenFOAM simpleFoam, RANS over a box | `cavity` |
| `fea` (`calculix`, `structural`) | CalculiX linear-static cantilever | `cantilever` |
| `chemistry` (`cantera`, `chem`) | Cantera equilibrium-HP, methane / air | `ch4-equilibrium` |
| `su2` (`compressible`, `aero`) | SU2 NACA 0012 airfoil starter | `naca0012` |
| `openradioss` (`crash`, `impact`) | OpenRadioss explicit dynamics | `drop-test` |
| `code-aster` (`aster`, `thermomechanical`) | Code_Aster `as_run` on a `.export` | `static-beam` |
| `netgen` (`meshing`, `mesh`) | Netgen CSG meshing — unit cube | `csg-box` |
| `meep` (`fdtd`, `photonics`) | Meep FDTD ring-resonator | `ring-resonator` |
| `gromacs` (`md`, `molecular-dynamics`) | GROMACS `gmx mdrun` on a `.tpr` | `lysozyme` |
| `gmsh` (`delaunay`) | gmsh procedural box mesh | `box-mesh` |
| `lammps` (`lj`, `classical-md`) | LAMMPS LJ FCC fluid (NVE) | `lj-fluid` |
| `elmer-heat` (`elmer`, `heat`) | Elmer steady heat conduction | `heat-cube` |

## `valenx-validate`

Structural pre-flight on a `.valenx` project. Walks the project
directory, loads the manifest, optional `tools.lock`, and every
`case.toml` listed in `[cases].order`, then reports a punch list
of any structural issues it finds.

```bash
valenx-validate path/to/project.valenx
valenx-validate path/to/project.valenx --format json
```

Designed for CI pre-flight gates: exits 0 on a clean project, 1 on
a structural issue, so the recipe can wire it in front of any
adapter run. Doesn't depend on the full adapter zoo, so the binary
stays small.

## `valenx-mesh-info`

Inspect a canonical mesh JSON (the `mesh.canonical.json` produced
by `gmsh.collect()` / `netgen.collect()`). Prints quality
statistics + AR / skewness histograms in either text or JSON.

```bash
valenx-mesh-info mesh.json
valenx-mesh-info mesh.json --format json
valenx-mesh-info mesh.json --check max-skew=0.9 --check inverted=0
```

Quality-gate flag (`--check METRIC=THRESHOLD`) is repeatable. Exits
4 if any check fails — useful for CI that wants a hard gate on
mesh quality. Supported metrics:

| Metric | Direction |
|---|---|
| `max-skew=<float>` | upper bound on `max_skewness` |
| `max-aspect=<float>` | upper bound on `max_aspect_ratio` |
| `inverted=<int>` | upper bound on `inverted_count` |
| `min-orthogonality=<float>` | lower bound on `min_orthogonality` |

## `valenx-audit`

Inspect Valenx's append-only audit log (`audit.log.jsonl`). Two
subcommands: `verify` checks the SHA-256 chain integrity; `tail`
prints recent entries.

```bash
valenx-audit verify ~/.local/state/valenx/audit.log.jsonl
valenx-audit tail -n 50 ~/.local/state/valenx/audit.log.jsonl
valenx-audit tail --since 2026-04-28T00:00:00Z log.jsonl
valenx-audit tail --json log.jsonl | jq '.[].action'
```

`tail --since <ISO-8601>` filters entries by timestamp before the
ring-buffer truncation runs, so combining with `-n` keeps the last
N entries on or after the cutoff.

## `valenx-results`

Headless inspector for the `results.json` sidecar that ValenxApp
writes next to every finished run. Lists fields (with timestep
counts), scalars (with sample counts), and artifacts; the
provenance header includes adapter / tool versions and the run
UUID.

```bash
valenx-results /path/to/workdir/results.json
valenx-results /path/to/workdir/results.json --format json
```

The text output caps the per-name scalar list at 12 entries to keep
a 10k-row PyBaMM time series from blasting the terminal — drop to
`--format json` for the firehose.

## `valenx-report`

Headless HTML / Markdown / CSV exporter for a `results.json`.
Wraps the `valenx-export` library's `write_html_report`,
`write_markdown_report`, and `write_scalars_csv` helpers and
refuses to no-op silently — at least one of `--html` /
`--markdown` / `--csv` must be supplied.

```bash
valenx-report results.json --html report.html
valenx-report results.json --markdown summary.md
valenx-report results.json --csv scalars.csv
valenx-report results.json --html out.html --markdown out.md --csv out.csv
```

The Markdown shape is GitHub-flavoured (pipe-delimited tables, code
spans for identifiers) — ideal for posting as a pull-request
comment or dropping into release notes.

## CI-pipeline recipe

The full workflow loop without the GUI:

```bash
# 1. Generate a project from a template.
valenx-init pr-12345 --template cfd

# 2. Pre-flight: confirm the project loads cleanly.
valenx-validate pr-12345 --format json | jq '.ok'  # true

# 3. (Run the solver from the GUI or via direct adapter invocation.)

# 4. Post-run: inspect what the solver produced.
valenx-results pr-12345/cases/cavity/workdir/results.json

# 5. Attach a Markdown report to the PR.
valenx-report pr-12345/cases/cavity/workdir/results.json \
  --markdown comments/run-summary.md

# 6. Confirm the audit trail is intact.
valenx-audit verify ~/.local/state/valenx/audit.log.jsonl
```

Every step is exit-code-driven so a failure short-circuits the
recipe. Combine with `valenx-mesh-info --check max-skew=0.9` to
enforce mesh-quality gates on the meshing step.
