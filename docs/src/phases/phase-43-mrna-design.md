# Phase 43 — mRNA design

**Status:** 🟢 Live — DNA Chisel + LinearDesign + iCodon open the
**first mRNA / vaccine therapeutic design domain** in Valenx
alongside the Phase 5.5 / 5.6 / 5.7 / 17 / 17.5 / 17.7 / 18 / 18.5
/ 18.6 / 18.7 / 19 / 19.5 / 19.6 / 20 / 22 / 22.5 / 23 / 24 / 25 /
27 / 27.5 / 27.6 / 28 / 29 / 30 / 30.5 / 31 / 32 / 32.5 / 33 / 34 /
35 / 36 / 38 / 39 / 40 / 41 / 42 biology / biotech / chemistry
beachheads.

## Goal

Open the **mRNA / vaccine therapeutic design** corner of the bio
surface in Valenx with three established open-source tools that
span the codon-optimization + mRNA-design tradeoff space — the
Edinburgh Genome Foundry's general-purpose constraint-driven codon
optimizer that's the de-facto Python choice for synthetic-gene
design (DNA Chisel, the MIT-licensed library that handles codon
optimization, restriction-site avoidance, repeat scanning, GC-
content tuning, and the long tail of constraint-driven sequence
design programmatically), Baidu Research's joint codon +
secondary-structure mRNA design tool that landed as the modern
mRNA-vaccine design workhorse since the 2021 _Nature_ paper
(LinearDesign, the Apache-2.0 single-binary CLI that jointly
optimizes codon usage and mRNA secondary-structure stability under
a tunable Lagrangian tradeoff parameter), and the Vejnar lab's
codon-level mRNA stability predictor (iCodon, the GPL-3.0 R-based
tool that scores per-position codon contributions to mRNA
half-life given a target organism). DNA Chisel + iCodon follow
the established Phase 17 Biopython + Phase 19.6 Seurat subprocess
patterns: the user supplies a `.py` / `.R` script that imports the
upstream package and reads `valenx_params.json` for the parsed
knobs. LinearDesign follows the established Phase 18 BWA single-
binary CLI pattern: protein FASTA in, optimized mRNA out. Phase 43
sits numerically after Phase 42 web visualization and ships
chronologically right after Phase 42 — same chronological-vs-
numerical convention used for Phase 17.5 / 24 / 28 / 31 / 35 / 39
/ 5.5 / 5.6 / 5.7 / 32.5 / 40 / 41 / 22.5 / 42.

## Capability inventory

### Live adapters (3)

- **DNA Chisel** — the Edinburgh Genome Foundry's constraint-
  driven codon-optimization library (MIT). DNA Chisel handles the
  long tail of synthetic-gene design programmatically — codon
  optimization against a target organism's codon usage table,
  restriction-site avoidance (one of the canonical pain points
  in cloning workflows), repeat scanning, GC-content tuning,
  forbidden-pattern matching, and arbitrary user-defined
  constraints composed into a single optimization objective. The
  library is the de-facto Python choice for end-to-end synthetic-
  gene design pipelines feeding into Phase 33 j5 assembly and
  Phase 41 pydna cloning workflows. Wrapped via the Python-script
  subprocess pattern (sister to Phase 17 Biopython, Phase 19.5
  Scanpy, Phase 33 pySBOL, Phase 41 pydna, Phase 42 Mol* / NGL
  Viewer). The user supplies a `.py` script referenced from
  `[bio.dnachisel].script` in `case.toml` that imports `dnachisel`
  and reads `valenx_params.json` for the parsed knobs. Schema
  knobs: `script` (path to user-supplied Python script; required,
  `.py` enforced), `python` (interpreter name; default
  `"python3"`), `input_fasta` (`Option<PathBuf>` — optional input
  FASTA the script can use as the starting sequence — `.fa` /
  `.fasta`; `None` when the script generates the sequence from
  scratch or fetches it inline), `output_basename` (filename stem
  the user's script uses for outputs — surfaced here so collect()
  can label artifacts uniformly; required, non-empty).
  `prepare()` enforces the `.py` extension on the script,
  resolves `script` and the optional `input_fasta` against the
  case directory when relative, **routes both staged paths through
  `confined_join`** so the staged filenames cannot escape the
  workdir, then writes a flat hand-rolled `valenx_params.json`
  containing `output_basename` always plus `input_fasta` (staged
  filename) only when set — the key is omitted entirely when
  `None` rather than emitted as `null`, matching the hand-rolled
  JSON convention the rest of the bio adapters use (Phase 19.6
  Seurat / AnnData, Phase 27.5 ESM-IF, Phase 41 pydna, Phase 42
  Mol* / NGL). `collect()` walks the workdir for
  `<output_basename>*.fasta` (`Native`, "DNA Chisel optimized
  FASTA"), `<output_basename>*.gb` / `.genbank` (`Native`, "DNA
  Chisel GenBank"), `<output_basename>*.json` (`Tabular`, "DNA
  Chisel constraint report" — the canonical machine-readable
  per-constraint-pass / per-constraint-fail report DNA Chisel
  emits for downstream automation), `<output_basename>*.png`
  (`Native`, "DNA Chisel plot"), and `*.log` (`Log`). Probe via
  Python on PATH then `<python> -c "import dnachisel"` — when the
  `import dnachisel` check fails the probe still returns `ok =
  true` with a targeted `"probe found python on PATH but could
  not import dnachisel — install with pip install dnachisel"`
  warning so users with Python ready but no `dnachisel` package
  see the install hint without failing the probe (sister to the
  Phase 19.5 scanpy / scvi / Phase 19.6 AnnData / Phase 5.6
  HOOMD-blue / Phase 5.7 MDTraj / Phase 41 pydna / Phase 42 Mol*
  / NGL probe convention). Version range `3.0.0..4.0.0` (DNA
  Chisel 3.x is the modern stable line shipping the contemporary
  constraint API + multi-objective optimizer; upper bound 4.0
  reserves room for the next major bump). `bio.dnachisel.optimize`
  ribbon capability.
- **LinearDesign** — Baidu Research's joint codon + secondary-
  structure mRNA design tool (Apache-2.0). LinearDesign landed as
  the modern mRNA-vaccine design workhorse following the 2021
  _Nature_ paper that demonstrated dramatic stability /
  expression gains for mRNA vaccines designed under the joint
  codon-usage + minimum-free-energy objective. The tool consumes
  a target protein sequence and emits an optimized mRNA sequence
  that maximizes a Lagrangian tradeoff between codon-adaptation-
  index (CAI) and predicted mRNA secondary-structure stability
  (MFE), tunable via the `lambda_param` knob — `lambda_param =
  0.0` collapses to pure MFE-optimal design, large `lambda_param`
  collapses to pure CAI-optimal design, intermediate values
  (default 1.0) hit the joint sweet spot demonstrated in the
  paper. Single-binary CLI — sister to Phase 18 BWA, Phase 32.5
  Smoldyn, Phase 5 GROMACS — with `lineardesign --aa <protein>
  --lambda <lambda_param> --codon_usage <codon_usage>
  --output_basename <basename> [extras...]`. Schema knobs:
  `protein` (path to protein FASTA; required — read in place, no
  staging), `output_basename` (filename stem LinearDesign uses
  for outputs; required, non-empty), `lambda_param` (`f64`,
  finite and ≥ 0.0; default 1.0; **note** the Rust field is
  `lambda_param` because `lambda` is a Rust reserved keyword —
  the CLI emits `--lambda <value>` regardless), `codon_usage`
  (target organism codon-usage table name; default `"human"` —
  selectable from the LinearDesign-shipped set: `"human"` /
  `"mouse"` / `"yeast"` / `"ecoli"` / etc.), `extra_args`.
  `prepare()` validates `lambda_param` is finite and ≥ 0.0
  (returns `InvalidCase` when the value is negative or NaN —
  LinearDesign's optimizer would either crash or silently
  collapse to MFE-only design on invalid input), resolves
  `protein` against the case directory when relative, validates
  the file exists on disk, and composes the invocation. Note
  that the `protein` FASTA is read in place rather than staged
  into the workdir via `confined_join` — same shape as Phase 18
  BWA's reference genome, Phase 35 Cas-OFFinder's reference, the
  Phase 18.5 / 18.6 / 18.7 aligner reference handling pattern.
  `collect()` walks the workdir for `<output_basename>*.fasta`
  (`Native`, "LinearDesign optimized mRNA"),
  `<output_basename>*.txt` (`Tabular`, "LinearDesign report" —
  the canonical per-design summary report LinearDesign writes
  alongside the FASTA output), and `*.log` (`Log`). Probe via
  `find_on_path(&["lineardesign"])` — when the `lineardesign`
  binary isn't found but Python is on PATH the probe surfaces a
  targeted `"probe found python on PATH but could not find
  lineardesign — clone https://github.com/LinearDesignSoftware/
  LinearDesign and add the bin directory to PATH"` warning so
  users see the install hint immediately. Version range
  `1.0.0..2.0.0` (LinearDesign 1.x is the modern stable line
  shipping the contemporary CAI + MFE joint-optimization
  algorithm; upper bound 2.0 reserves room for an eventual major
  bump). `bio.lineardesign.design` ribbon capability.
- **iCodon** — the Vejnar lab's codon-level mRNA stability
  prediction tool (GPL-3.0). iCodon predicts per-position codon
  contributions to mRNA half-life given a target organism — the
  canonical readout for "given this mRNA sequence, which codons
  are dragging down stability and where would a codon swap help".
  iCodon is the canonical R-based mRNA stability predictor and
  ships as a `devtools::install_github('santiago1234/iCodon')`
  R package. Wrapped via the Rscript subprocess pattern (sister
  to Phase 19.6 Seurat) — the user supplies an `.R` script
  referenced from `[bio.icodon].script` in `case.toml` that
  loads `library(iCodon)` and reads `valenx_params.json` for
  the parsed knobs via `jsonlite::fromJSON`. Schema knobs:
  `script` (path to user-supplied R script; required, `.R`
  enforced), `rscript` (R interpreter name; default
  `"Rscript"`), `input_fasta` (`Option<PathBuf>` — optional
  input mRNA FASTA the script can use as the sequence to score;
  `None` when the script generates the sequence inline or
  reads it from a different source), `output_basename`
  (filename stem the user's script uses for outputs; required,
  non-empty). `prepare()` enforces the `.R` extension on the
  script, resolves `script` and the optional `input_fasta`
  against the case directory when relative, routes both staged
  paths through `confined_join` so the staged filenames cannot
  escape the workdir, then writes a flat hand-rolled
  `valenx_params.json` containing `output_basename` always plus
  `input_fasta` (staged filename) only when set — the key is
  omitted entirely when `None` rather than emitted as `null`
  (same hand-rolled JSON convention as DNA Chisel above and
  every other Phase 19.6+ adapter that takes an optional input
  file). `collect()` walks the workdir for
  `<output_basename>*.csv` / `*.tsv` (`Tabular`, "iCodon
  stability table"), `<output_basename>*.rds` (`Native`,
  "iCodon R object (RDS)" — canonical R-serialised iCodon model
  output consumed by every downstream R-side stability /
  visualization pipeline), `<output_basename>*.png` (`Native`,
  "iCodon plot"), and `*.log` (`Log`). Probe via
  `find_on_path(&["Rscript"])` — does not attempt to confirm
  iCodon itself is installed because that would require running
  R, an expensive multi-second startup at probe time (same
  shape as Phase 19.6 Seurat); the `ToolNotInstalled` install
  hint mentions the canonical
  `devtools::install_github('santiago1234/iCodon')` install
  path. Version range `1.0.0..2.0.0` (iCodon 1.x is the modern
  stable line; upper bound 2.0 reserves room for an eventual
  major bump). `bio.icodon.predict` ribbon capability.

### Canonical types

**No new canonical types.** All three adapters consume user-
supplied inputs (Python / R scripts that import `dnachisel` /
`iCodon`, plus optional starting `.fa` / `.fasta` FASTA files for
DNA Chisel + iCodon, plus the LinearDesign `--aa <protein>` FASTA
input) and emit user-readable artifacts (DNA Chisel's optimized
`.fasta` / `.gb` GenBank + per-constraint `.json` report + plot
`.png`, LinearDesign's optimized `.fasta` mRNA + per-design
`.txt` report, iCodon's `.csv` / `.tsv` stability tables + `.rds`
R-serialised model output + plot `.png`) that the unchanged
`Results.artifacts` collection model surfaces directly. The
existing `valenx_bio::format::fasta` reader already inspects
collected FASTA inputs / outputs for sequence count + alphabet.
A first-class mRNA-design canonical type — a typed codon-
optimization / stability-prediction representation spanning DNA
Chisel constraint reports, LinearDesign per-design summaries,
and iCodon per-position stability scores, with parsed
constraint-pass / CAI / MFE / per-codon-half-life graphs —
defers to a future phase along with cross-tool design-comparison
viewers and per-codon stability-trace inspection CLIs.

### Headless CLIs

**No new CLIs.** DNA Chisel's `.fasta` / `.gb` GenBank outputs +
constraint `.json` reports + plot `.png` images, LinearDesign's
`.fasta` mRNA + `.txt` per-design summaries, and iCodon's `.csv`
/ `.tsv` stability tables + `.rds` R objects + plot `.png`
images are all standard formats inspectable in any editor or
through the user's downstream Python / R pipeline. Input FASTAs
are inspectable through the existing Phase 17 `valenx-fasta`
CLI. A canonical mRNA-design CLI — codon-optimization-pass
diffing, CAI / MFE / stability-score comparison, per-codon
stability-trace inspection — defers to a future phase along with
the canonical type.

## Domain milestone

Phase 43 opens the **first mRNA / vaccine therapeutic design
domain** in Valenx — the codon-optimization + joint-design half
of the mRNA workflow that the existing RNA structure-prediction
stack (Phase 28 ViennaRNA / RNAstructure / NUPACK) and synthetic-
biology composition stack (Phase 33 pySBOL / j5 / Cello + Phase
41 pydna) leave incomplete. Phase 28 covers RNA secondary-
structure prediction over a given mRNA / non-coding RNA
sequence; Phase 33 + 41 cover synthetic-gene composition and
plasmid-design programmatically. But until Phase 43 the
**codon-optimization + joint-design** surface (where the input
is a target protein sequence and the output is an optimized
mRNA encoding that protein under codon-usage / stability /
constraint objectives) was absent. Phase 43 closes the gap with
three established open-source tools spanning the mRNA-design
tradeoff space — DNA Chisel at the constraint-driven Python
end (codon optimization + restriction-site avoidance + repeat
scanning + GC tuning + arbitrary user constraints), LinearDesign
at the joint codon + secondary-structure end (the canonical
modern mRNA-vaccine design workhorse since the 2021 _Nature_
paper), and iCodon at the stability-prediction end (per-codon
mRNA half-life scoring for downstream codon-swap optimization).
This opens the canonical mRNA-vaccine and synthetic-gene
therapeutic design pipeline in Valenx — protein → optimized
mRNA → predicted secondary structure → predicted stability →
plasmid composition → assembly — all in one Valenx shell.

## What landed early

The implementation rode subagent-driven-development across 4
discrete implementation commits (3 adapters plus the registry +
init-template rollup) plus this docs pass — each landing one
adapter, the registry rollup, the init-template extension, or
the documentation pass. Every commit kept workspace clippy +
rustdoc clean.

## Acceptance checklist

- [x] `valenx-adapter-dnachisel` adapter ships with case-input
      parser + 4 lib tests + 3 case-input tests covering parses-
      minimal / parses-with-input-fasta / rejects-non-py-script,
      plus the Python-script subprocess shape that enforces
      `.py`, routes script + optional input_fasta through
      `confined_join`, writes `valenx_params.json` with
      `output_basename` always plus `input_fasta` (staged
      filename) only when set — key omitted entirely when `None`
      rather than emitted as `null`, matching the hand-rolled
      JSON convention the rest of the bio adapters use, plus the
      Python on PATH + `import dnachisel` probe with `"probe
      found python on PATH but could not import dnachisel"`
      warning when the import fails
- [x] `valenx-adapter-lineardesign` adapter ships with case-
      input parser + 4 lib tests + 4 case-input tests covering
      parses-minimal / parses-with-overrides / rejects-negative-
      lambda / rejects-empty-protein, plus the single-binary
      subprocess shape that composes `lineardesign --aa
      <protein> --lambda <lambda_param> --codon_usage
      <codon_usage> --output_basename <basename> [extras...]`
      with `protein` resolved against the case directory and
      validated on disk (read in place, no staging), and the
      `find_on_path(["lineardesign"])` probe with the
      `"clone https://github.com/LinearDesignSoftware/
      LinearDesign and add the bin directory to PATH"` warning
      whenever Python is on PATH but `lineardesign` is missing
- [x] `valenx-adapter-icodon` adapter ships with case-input
      parser + 4 lib tests + 3 case-input tests covering parses-
      minimal / parses-with-input-fasta / rejects-non-R-script,
      plus the Rscript subprocess shape that enforces `.R`,
      routes script + optional input_fasta through
      `confined_join`, writes `valenx_params.json` with the same
      hand-rolled JSON shape as DNA Chisel + Seurat (key omitted
      when `None`), and the `find_on_path(["Rscript"])` probe
      with the `devtools::install_github('santiago1234/iCodon')`
      install hint surfaced via `ToolNotInstalled`
- [x] All 3 adapters wired into `valenx-app::init_registry` —
      live adapter count moves from 124 to **127**, opening the
      first mRNA / vaccine therapeutic design domain to ship in
      Valenx
- [x] 3 mRNA-design templates in `valenx-init` (`dnachisel`,
      `lineardesign`, `icodon`), all round-tripping through
      `valenx-validate` (cross-binary roundtrip now sweeps **123
      templates** clean)
- [x] STATUS.md / README.md / ROADMAP.md / CHANGELOG.md /
      QUICKSTART.md updated
- [ ] **Phase 43.5** — sister-adapter expansion of Phase 43:
      OptiPyMer (Roche's Python codon-optimizer with a slightly
      narrower constraint surface than DNA Chisel; defer until
      upstream activity stabilises), CodonW (the legacy CAI /
      ENC / GC3 calculator from John Peden's PhD work; defer —
      modern DNA Chisel + iCodon cover the same ground), DEGRON
      (mRNA stability prediction sister to iCodon but Python-
      based; defer until upstream activity resumes), MRNA-Stab
      (Stanford's deep-learning mRNA stability predictor; defer
      pending model-checkpoint licensing review), CodonOpt
      (commercial codon-optimization service; out of scope as a
      hosted service). Out of scope for this beachhead.

## Success metrics

| Metric                                                | Target          |
|-------------------------------------------------------|-----------------|
| New mRNA-design adapter (template + tests)            | 1 day per       |
| Optimize-codons → design-mRNA → predict-stability loop across 3 tools | < tool baseline |

## Leads into

Phase 43 opens the **mRNA / vaccine therapeutic design** domain
that the user's bio / chemistry spec called out alongside the
existing Phase 28 RNA structure-prediction stack (ViennaRNA /
RNAstructure / NUPACK — every mRNA design deserves a secondary-
structure check) and the Phase 33 synthetic-biology composition
stack (pySBOL / j5 / Cello + Phase 41 pydna — every codon-
optimized mRNA needs to be composed into a plasmid /
construct). Combined with the existing render-Mol* → render-NGL
→ run-Galaxy-workflow → run-WDL → run-CWL → run-Nextflow → run-
Snakemake → design-plasmid → view-alignment → process-image →
segment-cells → classify-pixels → simulate-pathway → expand-
rules → grow-tissue → diffuse-particles → trace-MCell-
trajectories → simulate-MD → analyze-trajectory → reweight-
free-energy → fit-ENM → run-cpptraj-script → predict-structure
→ fold-RNA → analyze-DNA-geometry → infer-tree-ML → infer-tree-
Bayesian → simulate-popgen → analyze-trees → reconstruct-3D →
design-protein → validate loop, the **optimize-codons → design-
mRNA → predict-stability → render-Mol* → render-NGL → run-
Galaxy-workflow → run-WDL → run-CWL → run-Nextflow → run-
Snakemake → design-plasmid → view-alignment → process-image →
segment-cells → classify-pixels → simulate-pathway → expand-
rules → grow-tissue → diffuse-particles → trace-MCell-
trajectories → simulate-MD → analyze-trajectory → reweight-
free-energy → fit-ENM → run-cpptraj-script → predict-structure
→ fold-RNA → analyze-DNA-geometry → infer-tree-ML → infer-
tree-Bayesian → simulate-popgen → analyze-trees → reconstruct-
3D → design-protein → validate** loop now spans three mRNA-
design tools (DNA Chisel, LinearDesign, iCodon) feeding into
the existing Phase 28 RNA structure-prediction stack
(ViennaRNA, RNAstructure, NUPACK), the Phase 33 synthetic-
biology composition stack (pySBOL, j5, Cello), and the Phase
41 pydna plasmid-design library — closing the **canonical
mRNA-vaccine and synthetic-gene therapeutic design pipeline**:
target protein → DNA Chisel codon-optimized DNA / LinearDesign
joint-optimized mRNA → ViennaRNA / RNAstructure / NUPACK
secondary-structure check → iCodon stability prediction →
pydna / j5 plasmid composition → Cello logic-circuit
integration — all in one Valenx shell with no glue code beyond
the existing case-toml / prepare / run / collect path.

This **opens a new domain on top of the bio-ecosystem-complete
milestone reached at Phase 22.5 + 42** — Phase 43 layers mRNA /
vaccine therapeutic design on top of the bio surface that
already spans alignment, sequence editors, cheminformatics,
cryo-EM, CRISPR, DNA geometry, docking, MD analysis, MD engines,
microscopy, phylogenetics, population genetics, protein design,
quantum chemistry, RNA structure, sequence read simulators,
single-cell genomics, spatial stochastic simulation, structure
prediction, structure search, synthetic biology, systems biology,
variant calling, viewers (desktop + web), web visualization, and
workflow managers — **108 bio adapters across 39 biology /
biotech / chemistry phases**, all in one Valenx shell with no
glue code beyond the existing case-toml / prepare / run /
collect path.

The natural follow-up is **Phase 43.5** — the deferred mRNA-
design work called out above (OptiPyMer as a Python codon-
optimizer sister to DNA Chisel, CodonW for the legacy CAI / ENC
/ GC3 calculation surface if a use case emerges, DEGRON if
upstream activity resumes, MRNA-Stab pending model-checkpoint
licensing review), slotting in alongside the existing DNA
Chisel + LinearDesign + iCodon adapters with the same Python-
script subprocess shape (OptiPyMer / DEGRON / MRNA-Stab sister
tools), single-binary CLI shape (CodonW sister tool), or
Rscript subprocess shape (further iCodon-family R sisters). See
the out-of-scope section of `docs/superpowers/plans/2026-05-04-
mrna-design.md` for the full follow-up phase list.
