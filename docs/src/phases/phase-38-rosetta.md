# Phase 38 — Rosetta family

**Status:** 🟢 Live — Rosetta + PyRosetta open the
**first Rosetta protein-modeling family** in Valenx alongside the
Phase 17 / 17.5 / 18 / 25 / 27 / 27.5 / 27.6 / 28 / 29 / 30 / 31 /
32 / 34 / 35 / 36 biology + structure-prediction + protein-design +
RNA-structure + population-genetics + phylogenetics + read-simulator
+ systems-biology + docking + CRISPR-design + cryo-EM beachheads
and the Phase 24 cheminformatics expansion.

## Goal

Open the canonical Rosetta protein-modeling suite in Valenx with
the two most-used entry points into the RosettaCommons code base —
`rosetta_scripts` (the XML-driven protocol runner that's the
de-facto Rosetta entry point in production: every `relax` / `dock`
/ `abinitio` / FastDesign / enzyme-design pipeline lives as an XML
protocol fed to this binary) and PyRosetta (Python bindings
exposing the same C++ core through a Pythonic API for users who
prefer scripting Rosetta from `.py` rather than authoring XML
protocols). Rosetta = single-binary CLI shape (sister to Phase 18
BWA): the user supplies an XML protocol, an input PDB, the
`-database` path Rosetta data needs at runtime, an output basename,
and the number of decoys (`-nstruct`); the adapter composes the
positional `rosetta_scripts -database <path> -parser:protocol <xml>
-in:file:s <pdb> -out:prefix <basename> -nstruct <N> [extras...]`
invocation. PyRosetta = Python-script subprocess shape (sister to
Phase 17 Biopython): the user authors a `.py` driver, the adapter
stages script + optional input PDB + writes
`valenx_params.json`, and `run()` invokes `python <script>`. Both
adapters surface the RosettaCommons license accurately via
`tool_license = "Rosetta-License"` (a custom non-OSS license — not
a recognised SPDX identifier) and emit a probe warning whenever
the binary / bindings are detected, with the literal `"academic"`
string baked into the warning as a stable anchor for license-aware
filters and tests. Phase 38 sits numerically after Phase 36 cryo-EM
and ships chronologically right after Phase 29 population genetics.

## Capability inventory

### Live adapters (2)

- **Rosetta** — RosettaCommons' flagship modeling suite (custom
  Rosetta-License — academic / non-commercial use only without a
  separate commercial agreement). Rosetta drives protein design,
  structure prediction, docking, ligand binding, and a long tail of
  related modeling tasks through its `rosetta_scripts` binary,
  which reads an XML protocol describing the modeling pipeline
  (filters, movers, scorefunctions) and applies it to an input
  `.pdb`. Single-binary subprocess shape (sister to Phase 18 BWA):
  the CLI is `rosetta_scripts -database <path> -parser:protocol
  <protocol> -in:file:s <input_pdb> -out:prefix <output_basename>
  -nstruct <N> [extras...]`. The `database` knob is required —
  every `rosetta_scripts` invocation needs `-database <path>`
  pointing at the Rosetta data directory (energy tables, fragment
  libraries, etc.) which is bundled with the source distribution
  but isn't on PATH. Schema knobs: `protocol` (XML protocol script;
  required), `input_pdb` (input PDB; required), `output_basename`
  (filename stem the binary uses to label output decoys —
  `<basename>_0001.pdb` etc.; required, non-empty), `nstruct`
  (number of independent decoys to generate; `u32`, ≥ 1),
  `database` (path to the Rosetta `database/` directory; required),
  `extra_args`. `prepare()` resolves the protocol, input PDB, and
  database paths against the case directory when relative,
  validates the protocol + PDB exist on disk (returns `InvalidCase`
  with a helpful message when missing), and composes the positional
  invocation. `run()` streams Rosetta's `protocols.jd2` startup
  banner / `apply` per-mover lines / `Finished` / `successfully
  completed` end-of-run sentinels into progress hints. `collect()`
  walks the workdir for `<output_basename>*.pdb` (`Native`,
  "Rosetta designed structure") plus the canonical `score.sc`
  scorefile (`Tabular`, "Rosetta scores"). Probe via
  `find_on_path(&["rosetta_scripts", "rosetta_scripts.linuxgccrelease",
  "rosetta_scripts.macosclangrelease"])` — Rosetta source builds
  emit platform-suffixed names by default, conda / packaged
  distributions install a bare `rosetta_scripts` shim, and the
  probe covers all three. Version range `3.13.0..4.0.0` (the stable
  3.x line landed at 3.13 in 2021; upper bound 4.0 reserves room
  for an eventual major bump). **Academic-license-only** — probe
  always pushes a `"academic"`-keyworded warning into
  `ProbeReport.warnings` whenever Rosetta is detected, and
  `tool_license` surfaces as `"Rosetta-License"` rather than
  mislabeling the custom RosettaCommons terms as a recognised SPDX
  identifier. `bio.rosetta.protocol` ribbon capability.
- **PyRosetta** — Python bindings to the Rosetta C++ core
  (Rosetta-License — inherits the same academic / non-commercial
  use terms as the upstream Rosetta distribution). PyRosetta
  exposes the entire Rosetta modeling pipeline (movers, filters,
  scorefunctions, task-operations) through a Pythonic API, letting
  users drive Rosetta from regular `.py` scripts rather than
  authoring XML protocols. Python-script subprocess shape (sister
  to Phase 17 Biopython): the user supplies a `.py` driver
  referenced from `[bio.pyrosetta].script` in `case.toml`, plus an
  optional input PDB and a required output basename. Schema knobs:
  `script` (path to user-authored Python script; required),
  `python` (interpreter name; default `"python3"`), `input_pdb`
  (optional input PDB the script will operate on — None when the
  script generates structures de novo; surfaced in
  `valenx_params.json` so the script can read it without re-parsing
  case.toml), `output_basename` (filename stem; required,
  non-empty). `prepare()` stages the script (and PDB, when present)
  into the workdir under their original filenames so the script can
  resolve them via relative paths, then writes a flat
  `valenx_params.json` with `input_pdb` (staged filename or literal
  `null`) and `output_basename`. `run()` invokes `python <script>`
  via the shared subprocess runner; the script can emit a sentinel
  `[valenx] pyrosetta done` line on stdout to signal completion
  before exit (lifted to a 95% progress tick). `collect()` walks
  the workdir for `<output_basename>*.pdb` (`Native`, "PyRosetta
  designed structure") and `*.sc` files (`Tabular`, "PyRosetta
  scores"). Probe via Python on PATH with an `import pyrosetta`
  check — when the import fails the probe still returns `ok = true`
  with a warning so users with PyRosetta installed under a
  different interpreter (referenced via the case-level `python`
  override) aren't blocked. Version range `4.0.0..5.0.0` (the
  modern release line is the 4.x series with weekly nightly drops
  post-2017; upper bound 5.0 reserves room for an eventual major
  bump). **Academic-license-only** — probe always pushes a
  `"academic"`-keyworded warning into `ProbeReport.warnings`
  whenever Python is detected (regardless of whether `pyrosetta`
  itself is importable, since the user is either about to install
  it or has it installed and needs reminding), and `tool_license`
  surfaces as `"Rosetta-License"` rather than mislabeling the
  inherited RosettaCommons terms as MIT / BSD. `bio.pyrosetta.script`
  ribbon capability.

### Canonical types

**No new canonical types.** Both adapters consume user-supplied
inputs (Rosetta XML protocols + input PDBs + `database/` data
directories, PyRosetta Python scripts + optional input PDBs) and
emit user-readable artifacts (`<basename>*.pdb` design decoys,
`score.sc` / `*.sc` scorefiles) that the unchanged
`Results.artifacts` collection model surfaces directly. The
existing `valenx_bio::format::pdb` reader inspects collected PDB
artifacts for chain / residue / atom counts. A first-class Rosetta
canonical type — a generic protocol + scorefile pair spanning both
back-ends, parsed into a typed scorefile model with
per-decoy energy terms — defers to a future phase along with
score-distribution visualizers and per-mutation Δ-energy heatmap
viewers.

### Headless CLIs

**No new CLIs.** Rosetta's `score.sc` tab-separated scorefiles and
PyRosetta's `*.sc` files are inspectable in any editor or through
the user's downstream Python pipeline (`pandas`, `numpy`); design
decoy `.pdb` files are inspectable through the existing Phase 17
`valenx-pdb-info` CLI. A canonical Rosetta-aware CLI — scorefile
diffing, top-N decoy ranking, per-residue energy-decomposition
inspection — defers to a future phase along with
score-distribution visualizers and per-mutation Δ-energy heatmap
viewers.

## Domain milestone

Phase 38 is the **first Rosetta protein-modeling family** to land
in Valenx. The biology adapter family started with Phase 17
(foundation — sequence / structure / trajectory canonical types +
classical MD + cheminformatics scripts) and expanded through Phase
17.5 / 18 / 18.5 / 18.6 / 19 / 19.5 / 20 / 22 / 23 / 24 / 25 / 27
/ 27.5 / 27.6 / 28 / 29 / 30 / 31 / 32 / 34 / 35 / 36 to cover
sequence prediction, alignment, RNA-seq, variant calling,
single-cell, transcript quantification, workflow orchestration,
molecular viewers, cheminformatics, quantum chemistry, protein
design, EvolutionaryScale models, RNA structure, population
genetics, phylogenetics, sequencing read simulation, systems
biology, small-molecule docking, CRISPR design, and cryo-EM
reconstruction — but until Phase 38 the canonical Rosetta surface
(XML-protocol-driven modeling, Python-bindings access to the core)
was absent. Phase 38 closes that gap with the two most-used
entry points into the RosettaCommons code base — `rosetta_scripts`
(production XML-driven protocol runner) and PyRosetta
(Python bindings to the same C++ core).

## What landed early

The implementation rode subagent-driven-development across 4
discrete implementation commits (2 adapters, the registry rollup,
the init-template extension) plus this docs pass — each landing one
adapter, the registry rollup, the init-template extension, or the
documentation pass. Every commit kept workspace clippy + rustdoc
clean.

## Acceptance checklist

- [x] `valenx-adapter-rosetta` adapter ships with case-input parser
      + 4 lib tests + 5 case-input tests + 1 academic-warning test
      covering parses-minimal / rejects-empty-protocol / rejects-
      empty-input-pdb / rejects-empty-output-basename / rejects-
      bad-nstruct, plus the single-binary subprocess shape that
      composes `rosetta_scripts -database <path> -parser:protocol
      <xml> -in:file:s <pdb> -out:prefix <basename> -nstruct <N>
      [extras...]` and the mandatory `"academic"`-anchored license
      warning surfaced via `ProbeReport.warnings`
- [x] `valenx-adapter-pyrosetta` adapter ships with case-input
      parser + 4 lib tests + 3 case-input tests + 1 academic-
      warning test covering parses-minimal / parses-with-input-pdb
      / rejects-empty-script / rejects-empty-output-basename, plus
      the Python-script subprocess shape that stages script +
      optional input PDB + writes `valenx_params.json` with
      `input_pdb` (staged filename or literal `null`) and
      `output_basename`, and the mandatory `"academic"`-anchored
      license warning surfaced via `ProbeReport.warnings`
      regardless of whether `pyrosetta` is importable in the
      detected Python interpreter
- [x] Both adapters wired into `valenx-app::init_registry` —
      live adapter count moves from 84 to **89** (alongside the
      Phase 29 population-genetics trio), opening the first
      canonical Rosetta protein-modeling family to ship in Valenx
- [x] 2 Rosetta-family templates in `valenx-init` (`rosetta`,
      `pyrosetta`), all round-tripping through `valenx-validate`
      (cross-binary roundtrip now sweeps **85 templates** clean)
- [x] STATUS.md / README.md / ROADMAP.md / CHANGELOG.md /
      QUICKSTART.md updated
- [ ] **Phase 38.5** — specific Rosetta app adapters that wrap
      individual Rosetta apps (`relax`, `dock_protocol`,
      `AbinitioRelax`, `enzyme_design`, `loopmodel`) directly
      rather than going through `rosetta_scripts` XML protocols.
      Out of scope for this beachhead — `rosetta_scripts` covers
      all of them via XML protocols. Rosetta@home / RosettaCommons
      cluster execution (different shape, defer to a future phase).

## Success metrics

| Metric                                                | Target          |
|-------------------------------------------------------|-----------------|
| New Rosetta-family adapter (template + tests)         | 1 day per       |
| `rosetta_scripts` XML protocol → decoy ensemble loop  | < tool baseline |
| PyRosetta Python script → decoy ensemble loop         | < tool baseline |

## Leads into

Phase 38 opens the canonical Rosetta protein-modeling family that
the user's bio / chemistry spec called out alongside the Phase 17
/ 17.5 biology + structure-prediction stack and the Phase 27 /
27.5 / 27.6 protein-design beachhead. Combined with the existing
design-guides → predict-off-targets → search-off-targets →
simulate-reads → align → quantify → call-variants → predict-
structure → fold-RNA → infer-tree → simulate-pathway → reconstruct-
3D → simulate-popgen → analyze-trees → validate loop, the
**XML-protocol-design → Python-bindings-design → predict-structure
→ design-sequence → fold-RNA → infer-tree → simulate-pathway →
reconstruct-3D → simulate-popgen → analyze-trees → validate** loop
now spans two Rosetta entry points (`rosetta_scripts`, PyRosetta)
feeding into the existing Phase 27 / 27.5 / 27.6 protein-design
adapters (RFdiffusion, ProteinMPNN, Chroma, ESM-IF, RFantibody,
ESM3, ESM Cambrian), the Phase 17.5 prediction stack (ESMFold,
OpenFold, AlphaFold 2/3, ColabFold), the Phase 28 RNA-structure
tools (ViennaRNA, RNAstructure, NUPACK), the Phase 29 population-
genetics trio (SLiM, msprime, tskit), the Phase 30 phylogenetic-
tree builders (IQ-TREE, RAxML-NG, FastTree), the Phase 31 read
simulators (ART, wgsim, Badread), the Phase 32 systems-biology
surface (COPASI, BioNetGen, PhysiCell), the Phase 34 docking pair
(AutoDock Vina, AutoDock 4), the Phase 35 CRISPR-design tools
(CHOPCHOP, CRISPOR, Cas-OFFinder), and the Phase 36 cryo-EM
reconstruction tools (RELION, EMAN2, CTFFIND) — all in one Valenx
shell with no glue code beyond the existing case-toml / prepare /
run / collect path.

The natural follow-up is **Phase 38.5** — the deferred
Rosetta-family work called out above (specific Rosetta app
adapters wrapping individual apps like `relax`, `dock_protocol`,
`AbinitioRelax`, `enzyme_design`, `loopmodel` directly rather than
going through `rosetta_scripts` XML protocols), slotting in
alongside the existing Rosetta + PyRosetta adapters with the same
single-binary subprocess shape (Rosetta sister tools) or
Python-script subprocess shape (PyRosetta sister tools).
Rosetta@home / RosettaCommons cluster execution sits in a separate
phase — the data shape is different enough (distributed work-unit
dispatch rather than single-host modeling) to warrant a sister
phase rather than 38.5 expansion. See the out-of-scope section of
`docs/superpowers/plans/2026-04-30-rosetta.md` for the full
follow-up phase list.
