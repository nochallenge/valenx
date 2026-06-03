# Phase 39 — DNA structural geometry

**Status:** 🟢 Live — X3DNA + Curves+ + DSSR open the
**first DNA structural geometry domain** in Valenx alongside the
Phase 17 / 17.5 / 18 / 25 / 27 / 27.5 / 27.6 / 28 / 29 / 30 / 30.5
/ 31 / 32 / 34 / 35 / 36 / 38 biology + structure-prediction +
protein-design + RNA-structure + population-genetics +
phylogenetics + Bayesian-phylogenetics + read-simulator + systems-
biology + docking + CRISPR-design + cryo-EM + Rosetta-family
beachheads and the Phase 24 cheminformatics expansion.

## Goal

Open the DNA / RNA structural-geometry analysis domain in Valenx
with three established open-source tools that span the structural-
geometry tradeoff space — base-pair / base-step parameter
calculation (X3DNA, the de-facto reference for canonical helical
parameters: twist, roll, tilt, slide, shift, rise plus per-base
intra-pair parameters buckle, propeller, opening, shear, stretch,
stagger), helical-axis curvature analysis with groove-geometry
characterisation (Curves+, the modern Curves successor — the
canonical tool for "is this DNA bent, and if so, how" questions in
protein-DNA / drug-DNA structural studies), and structural-feature
annotation as a single machine-readable JSON summary (DSSR,
Dissecting the Spatial Structure of RNA — X3DNA's modern Python-
fronted tool that enumerates every detected feature: base pairs,
multiplets, helices, stems, hairpin / internal / junction loops,
kissing loops, A-minor motifs, ribose zippers, pseudoknots,
splayed-apart conformations). All three are single-binary CLIs
(sister to Phase 18 BWA): X3DNA's `analyze` is positional-only
(input PDB → derived output filenames), Curves+'s `Cur+` takes
its parameters via stdin (a Fortran-style `&inp ... &end`
namelist block followed by strand / axis residue cards — fed via
the Phase 36 CTFFIND-style stdin-feed pattern with
`Stdio::from(file)`), and DSSR's `x3dna-dssr` takes flag-form
arguments (`-i=<pdb> -o=<json> --json`). All three are
**academic-license-flagged** — the X3DNA family (X3DNA + DSSR)
ships under custom non-OSS academic terms, Curves+ ships under a
custom non-OSS academic license — Valenx surfaces these accurately
via `tool_license = "X3DNA-License"` / `"Curves-License"` /
`"DSSR-License"` and a mandatory `"academic"`-keyworded probe
warning whenever each binary is detected. Phase 39 sits
numerically after Phase 38 Rosetta and ships chronologically right
after Phase 30.5 Bayesian phylogenetics.

## Capability inventory

### Live adapters (3)

- **X3DNA** — Wilma Olson and Xiang-Jun Lu's reference toolkit for
  DNA / RNA structural-geometry analysis (custom X3DNA-License —
  academic / non-commercial use only). X3DNA reads a nucleic-acid
  PDB, identifies base pairs, and computes the canonical helical-
  step parameters (twist, roll, tilt, slide, shift, rise) plus
  per-base intra-pair parameters (buckle, propeller, opening,
  shear, stretch, stagger). It is the workhorse behind structural-
  bioinformatics pipelines that need quantitative DNA geometry —
  bending studies, drug-DNA / protein-DNA complex analysis, RNA
  tertiary-structure annotation. Single-binary subprocess shape
  (sister to Phase 18 BWA): the CLI is `analyze <input_pdb>
  [extras...]`. `analyze` is positional-only — it derives every
  output filename from the input basename, so the adapter just
  hands it the PDB and any user-supplied extras. Schema knobs:
  `input_pdb` (input PDB; required), `output_basename` (filename
  stem the user expects X3DNA to produce — surfaced here so
  `collect()` can label artefacts uniformly without scraping
  `analyze`'s filename heuristics; required, non-empty),
  `extra_args`. `prepare()` resolves the input PDB against the
  case directory when relative, validates it exists on disk
  (returns `InvalidCase` with a helpful message when missing), and
  composes the positional invocation. `collect()` walks the
  workdir for `<output_basename>*.par` (`Tabular`, "X3DNA base-
  step parameters") and `*.out` (`Log`, the per-run log `analyze`
  writes alongside). Probe via `find_on_path(&["analyze"])`
  (X3DNA's main analysis binary is literally named `analyze`).
  Version range `2.4.0..3.0.0` (X3DNA 2.4 (2020) is the modern
  stable release and the floor we test against; upper bound 3.0
  reserves room for an eventual major bump). **Academic-license-
  only** — probe always pushes an `"academic"`-keyworded warning
  into `ProbeReport.warnings` whenever the binary is detected, and
  `tool_license` surfaces as `"X3DNA-License"` rather than
  mislabeling the custom X3DNA terms as a recognised SPDX
  identifier. `bio.x3dna.analyze` ribbon capability.
- **Curves+** — Richard Lavery's reference toolkit for DNA
  helical-axis analysis (custom Curves-License — academic /
  non-commercial use only). Curves+ fits a curvilinear helical
  axis through a nucleic-acid structure and reports per-base axis-
  curvature, base-pair parameters relative to that axis, and a
  `.cda` file describing the axis itself for downstream
  visualisation. It is the canonical tool for "is this DNA bent,
  and if so, how" questions in protein-DNA / drug-DNA structural
  studies. Single-binary subprocess shape with stdin-piped
  parameters (sister to Phase 36 CTFFIND): Curves+ takes its
  parameters as a Fortran-style `&inp ... &end` namelist block on
  stdin followed by strand / axis residue cards. The adapter
  authors that block at `prepare()` time and pipes it into `Cur+`'s
  stdin at `run()` time via `Stdio::from(file)` — the shared
  `subprocess::run` helper closes stdin which makes Curves+ read
  EOF before parsing its first parameter and exit, so the custom
  `run()` opens the parameters file with `File::open()` and hands
  its FD to the child via `Stdio::from(file)` (the custom run path
  mirrors the MAFFT stdout-redirect pattern but for stdin, same
  shape Phase 36 CTFFIND uses). Schema knobs: `input_pdb` (input
  PDB; required), `output_basename` (filename stem Curves+ uses
  for outputs — `<basename>.lis`, `<basename>.cda`, etc.;
  required, non-empty), `first_residue` (first inclusive residue
  index in the strand to analyse; `u32`), `last_residue` (last
  inclusive residue index in the strand; `u32`, ≥
  `first_residue` — a reverse range is rejected up front with a
  helpful message), `extra_args`. `prepare()` resolves the input
  PDB against the case directory when relative, validates it
  exists on disk, writes `curves_params.txt` containing the
  namelist body + residue-range cards, stashes the filename under
  the sentinel env var `VALENX_CURVES_PARAMS_FILE`, and the
  custom `run()` recovers the filename, strips the sentinel from
  the env table so Curves+ doesn't see it, opens the params file,
  and pipes its contents into the child. `collect()` walks the
  workdir for `<output_basename>*.lis` (`Log`, "Curves+ helical
  analysis") and `<output_basename>*.cda` (`Tabular`, "Curves+
  axis curve data"). Probe via `find_on_path(&["Cur+"])` (the
  binary name uses a literal `+`). Version range `2.0.0..3.0.0`
  (Curves+ 2.x is the modern stable line; 2.0 is the floor; upper
  bound 3.0 reserves room for an eventual major bump).
  **Academic-license-only** — probe always pushes an
  `"academic"`-keyworded warning into `ProbeReport.warnings`
  whenever the binary is detected, and `tool_license` surfaces as
  `"Curves-License"` rather than mislabeling the custom Curves+
  terms as a recognised SPDX identifier. `bio.curves.analyze`
  ribbon capability.
- **DSSR** — Dissecting the Spatial Structure of RNA / DNA, the
  modern Python-fronted X3DNA-family tool (custom DSSR-License —
  academic / non-commercial use only). DSSR reads a nucleic-acid
  PDB and emits a single JSON file enumerating every detected
  feature: base pairs (Watson-Crick, Hoogsteen, sugar-edge,
  ...), multiplets, double helices, stems, hairpin / internal /
  junction loops, kissing loops, A-minor motifs, ribose zippers,
  pseudoknots, splayed-apart conformations, and more. It is the
  standard machine-readable feature-extraction step in modern
  RNA-structure pipelines. Single-binary subprocess shape (sister
  to X3DNA): the CLI is `x3dna-dssr -i=<input_pdb>
  -o=<output_json> --json [extras...]` (DSSR uses `key=value`
  flag form on its short-form options — no space between flag and
  value). Schema knobs: `input_pdb` (input PDB; required),
  `output_json` (output JSON path; required), `extra_args`.
  `prepare()` resolves the input PDB against the case directory
  when relative, scopes the output JSON path to the workdir when
  relative, validates the input exists on disk, and composes the
  flag-form invocation. `collect()` reports the configured
  `output_json` file as a single `Tabular` artifact ("DSSR
  analysis (JSON)") — DSSR's JSON is the canonical machine-
  readable summary; tagged `Tabular` rather than `Native` so
  downstream serdes can key off a consistent kind. Probe via
  `find_on_path(&["x3dna-dssr"])`. Version range `2.0.0..3.0.0`
  (DSSR 2.x is the modern stable line that ships with X3DNA 2.4+;
  upper bound 3.0 reserves room for an eventual major bump).
  **Academic-license-only** — probe always pushes an
  `"academic"`-keyworded warning into `ProbeReport.warnings`
  whenever the binary is detected, and `tool_license` surfaces as
  `"DSSR-License"` rather than mislabeling the inherited X3DNA-
  family terms as a recognised SPDX identifier.
  `bio.dssr.analyze` ribbon capability.

### Canonical types

**No new canonical types.** All three adapters consume user-
supplied inputs (X3DNA / Curves+ / DSSR all take nucleic-acid
PDBs, plus the Curves+ residue-range knobs) and emit user-
readable artifacts (X3DNA `.par` base-step parameter tables and
`.out` per-run logs, Curves+ `.lis` helical-analysis logs and
`.cda` axis-curve data files, DSSR JSON structural-feature
summaries) that the unchanged `Results.artifacts` collection
model surfaces directly. The existing `valenx_bio::format::pdb`
reader inspects collected PDB inputs for chain / residue / atom
counts. A first-class DNA-geometry canonical type — a typed
helical-parameter representation spanning all three back-ends,
with parsed per-step parameter tables (twist, roll, tilt, slide,
shift, rise) and a typed structural-feature summary (DSSR JSON
parsed into a typed feature graph) — defers to a future phase
along with helical-axis visualizers and per-feature interactive
overlays.

### Headless CLIs

**No new CLIs.** X3DNA's `.par` parameter tables and `.out` run
logs, Curves+'s `.lis` analysis logs and `.cda` axis data files,
and DSSR's JSON summary are all standard formats inspectable in
any editor or through the user's downstream Python pipeline
(`pandas`, `numpy`, `json`); the input nucleic-acid PDBs are
inspectable through the existing Phase 17 `valenx-pdb-info` CLI.
A canonical DNA-geometry CLI — per-step parameter extraction,
helical-axis comparison, structural-feature diffing — defers to
a future phase along with helical-axis visualizers and per-
feature interactive overlays.

## Domain milestone

Phase 39 is the **first DNA structural-geometry domain** to land
in Valenx. The biology adapter family started with Phase 17
(foundation — sequence / structure / trajectory canonical types +
classical MD + cheminformatics scripts) and expanded through Phase
17.5 / 18 / 18.5 / 18.6 / 19 / 19.5 / 20 / 22 / 23 / 24 / 25 / 27
/ 27.5 / 27.6 / 28 / 29 / 30 / 30.5 / 31 / 32 / 34 / 35 / 36 / 38
to cover sequence prediction, alignment, RNA-seq, variant calling,
single-cell, transcript quantification, workflow orchestration,
molecular viewers, cheminformatics, quantum chemistry, protein
design, EvolutionaryScale models, RNA structure, population
genetics, phylogenetics, Bayesian phylogenetics, sequencing read
simulation, systems biology, small-molecule docking, CRISPR
design, cryo-EM reconstruction, and Rosetta protein modeling — but
until Phase 39 the DNA structural-geometry surface (canonical
helical parameters, helical-axis analysis, structural-feature
annotation) was absent. Phase 39 closes that gap with three
established open-source tools spanning the structural-geometry
tradeoff space — X3DNA at the canonical helical-parameter end,
Curves+ for the helical-axis curvature side, and DSSR as the
machine-readable structural-feature annotator that feeds modern
RNA-structure pipelines.

## What landed early

The implementation rode subagent-driven-development across 5
discrete implementation commits (3 adapters, the registry rollup,
the init-template extension) plus this docs pass — each landing
one adapter, the registry rollup, the init-template extension, or
the documentation pass. Every commit kept workspace clippy +
rustdoc clean.

## Acceptance checklist

- [x] `valenx-adapter-x3dna` adapter ships with case-input parser
      + 5 lib tests (incl. probe-warning) + 3 case-input tests
      covering parses-minimal / rejects-empty-input-pdb / rejects-
      empty-output-basename, plus the single-binary subprocess
      shape that composes `analyze <input_pdb> [extras...]` with
      no flag prefixes (`analyze` is positional-only and derives
      every output filename from the input basename), and the
      mandatory `"academic"`-anchored license warning surfaced via
      `ProbeReport.warnings`
- [x] `valenx-adapter-curves` adapter ships with case-input parser
      + 5 lib tests (incl. probe-warning) + 4 case-input tests
      covering parses-minimal / rejects-empty-input-pdb / rejects-
      empty-output-basename / rejects-reverse-residue-range, plus
      the stdin-feed subprocess shape that writes
      `curves_params.txt` with the namelist body + residue-range
      cards in `prepare()`, stashes the filename under
      `VALENX_CURVES_PARAMS_FILE`, and the custom `run()` opens
      the file and pipes its contents into `Cur+`'s stdin via
      `Stdio::from(file)`, plus the mandatory `"academic"`-
      anchored license warning surfaced via `ProbeReport.warnings`
- [x] `valenx-adapter-dssr` adapter ships with case-input parser
      + 5 lib tests (incl. probe-warning) + 3 case-input tests
      covering parses-minimal / rejects-empty-input-pdb / rejects-
      empty-output-json, plus the single-binary subprocess shape
      that composes `x3dna-dssr -i=<input_pdb> -o=<output_json>
      --json [extras...]` with the DSSR-style `key=value` flag
      form (no space between flag and value), and the mandatory
      `"academic"`-anchored license warning surfaced via
      `ProbeReport.warnings`
- [x] All 3 adapters wired into `valenx-app::init_registry` —
      live adapter count moves from 89 to **94** (alongside the
      Phase 30.5 Bayesian-phylogenetics pair), opening the first
      DNA structural-geometry domain to ship in Valenx
- [x] 3 DNA-geometry templates in `valenx-init` (`x3dna` with
      alias `3dna`, `curves` with alias `curves+`, `dssr`), all
      round-tripping through `valenx-validate` (cross-binary
      roundtrip now sweeps **90 templates** clean)
- [x] STATUS.md / README.md / ROADMAP.md / CHANGELOG.md /
      QUICKSTART.md updated
- [ ] **Phase 39.5** — sister-adapter expansion of Phase 39:
      web3DNA (web-fronted X3DNA companion; defer to sister-
      adapter expansion phase), 3D-DART (DNA-axis transformation
      tool; defer), MC-Sym / MC-Annotate (Major / Cedergren
      group's RNA structural annotation; defer to 39.5), Madbend
      (DNA-bending statistical analysis; defer), DNAtools (Berman
      lab's DNA-conformation toolkit; defer). Out of scope for
      this beachhead.

## Success metrics

| Metric                                                | Target          |
|-------------------------------------------------------|-----------------|
| New DNA-geometry adapter (template + tests)           | 1 day per       |
| Helical-parameter + axis-analysis + feature-extract loop across 3 tools | < tool baseline |

## Leads into

Phase 39 opens the DNA structural-geometry domain that the user's
bio / chemistry spec called out alongside the Phase 17 / 17.5
biology + structure-prediction stack and the Phase 28 RNA-
structure beachhead. Combined with the existing align → quantify
→ predict → infer-tree-ML → infer-tree-Bayesian → simulate-popgen
→ analyze-trees → simulate-pathway → reconstruct-3D → fold-RNA →
design-protein → validate loop, the **predict-structure → fold-
RNA → analyze-DNA-geometry → infer-tree-ML → infer-tree-Bayesian
→ simulate-popgen → analyze-trees → simulate-pathway →
reconstruct-3D → design-protein → validate** loop now spans three
DNA structural-geometry tools (X3DNA, Curves+, DSSR) feeding into
the existing Phase 17 / 17.5 prediction stack (ESMFold, OpenFold,
AlphaFold 2/3, ColabFold), the Phase 28 RNA-structure tools
(ViennaRNA, RNAstructure, NUPACK), the Phase 29 population-
genetics trio (SLiM, msprime, tskit), the Phase 30 phylogenetic-
tree builders (IQ-TREE, RAxML-NG, FastTree), the Phase 30.5
Bayesian-phylogenetics pair (BEAST 2, MrBayes), the Phase 32
systems-biology surface (COPASI, BioNetGen, PhysiCell), the Phase
34 docking pair (AutoDock Vina, AutoDock 4), the Phase 35 CRISPR-
design tools (CHOPCHOP, CRISPOR, Cas-OFFinder), the Phase 36
cryo-EM reconstruction tools (RELION, EMAN2, CTFFIND), and the
Phase 38 Rosetta-family adapters (Rosetta, PyRosetta) — all in
one Valenx shell with no glue code beyond the existing case-toml
/ prepare / run / collect path.

The natural follow-up is **Phase 39.5** — the deferred DNA-
geometry work called out above (web3DNA as a web-fronted X3DNA
companion sister to X3DNA, 3D-DART for DNA-axis transformation,
MC-Sym / MC-Annotate for the Major / Cedergren group's RNA
structural annotation, Madbend for DNA-bending statistical
analysis, DNAtools for the Berman lab's DNA-conformation
toolkit), slotting in alongside the existing X3DNA + Curves+ +
DSSR adapters with the same single-binary subprocess shape (X3DNA
/ DSSR sister tools) or stdin-feed subprocess shape (Curves+
sister tools). See the out-of-scope section of
`docs/superpowers/plans/2026-04-30-dna-geometry.md` for the full
follow-up phase list.
