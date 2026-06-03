# Phase 29 — Population genetics

**Status:** 🟢 Live — SLiM + msprime + tskit open the
**first population-genetics / evolutionary-simulation domain** in
Valenx alongside the Phase 17 / 17.5 / 18 / 25 / 27 / 27.5 / 27.6
/ 28 / 30 / 31 / 32 / 34 / 35 / 36 / 38 biology + structure-
prediction + protein-design + RNA-structure + phylogenetics +
read-simulator + systems-biology + docking + CRISPR-design +
cryo-EM + Rosetta-family beachheads and the Phase 24
cheminformatics expansion.

## Goal

Open the population-genetics / evolutionary-simulation domain in
Valenx with three established open-source tools that span the
population-genetics tradeoff space — forward-time individual-based
simulation under arbitrary selection / demography / mating-system
specifications (SLiM, the de-facto forward simulator), coalescent
backward-time simulation of sample ancestries under configurable
demographies (msprime, the de-facto coalescent simulator and
companion to SLiM), and tree-sequence analysis / statistics on the
succinct tree-sequence outputs both simulators emit (tskit, the
canonical analysis library). SLiM = single-binary CLI shape
(sister to Phase 18 BWA): the user supplies a `.slim` Eidos script
referenced from `[bio.slim].script`, optionally a `seed` (passed
via `slim -s <N>`) and Eidos-constant overrides via
`extra_args`, and the adapter composes `slim [-s <seed>] [extras...]
<script>`. msprime + tskit = Python-script subprocess shape (sister
to Phase 17 Biopython): the user authors a `.py` driver, the
adapter stages script + writes `valenx_params.json`, and `run()`
invokes `python <script>`. msprime's `valenx_params.json` carries
the demographic knobs (`population_size`, `num_samples`,
`recombination_rate`, `mutation_rate`, `output_basename`) the
script reads back via `json.load(open("valenx_params.json"))`;
tskit's carries the input `.trees` filename and output basename
for collected statistics tables. All three feed into each other —
SLiM emits `.trees` via `treeSeqOutput()`, msprime emits `.trees`
via `dump()`, tskit consumes both for downstream
population-genetics statistics (π, Tajima's D, Fst, site-frequency
spectra). Phase 29 sits numerically between Phase 28 RNA structure
and Phase 30 phylogenetics, and ships chronologically right after
Phase 35 CRISPR design.

## Capability inventory

### Live adapters (3)

- **SLiM** — Philipp Messer's forward-time population-genetics
  simulator (GPL-3.0). SLiM evolves a finite-population model
  generation by generation under a user-defined Eidos script:
  mutation rates, selection coefficients, recombination maps,
  demographic events, migrations, mating systems. The state is
  sampled at any generation the script asks for, and tree-sequence
  recording (the `treeSeqOutput()` family) feeds straight into
  tskit / msprime downstream. Single-binary subprocess shape
  (sister to Phase 18 BWA): the CLI is `slim [-s <seed>] [extras...]
  <script>`. Schema knobs: `script` (`.slim` Eidos model file;
  required), `seed` (optional `u64`; passed via `slim -s <N>`
  when present, otherwise SLiM picks its own seed and prints it on
  the run banner), `output_basename` (filename stem the user's
  script uses for outputs — surfaced here so `collect()` can label
  artefacts uniformly even though SLiM scripts choose their own
  output paths; required, non-empty), `extra_args` (additional CLI
  arguments appended after the script path; `-d KEY=VALUE` is the
  canonical way to inject Eidos constants from outside the script).
  `prepare()` resolves the script against the case directory when
  relative, validates it exists on disk, and composes the
  invocation with `seed` injected before any extras and the script
  positional last. `run()` streams SLiM's `// Initial random seed`
  banner / periodic `// generation N` lines / `// Run finished`
  end-of-run sentinel into progress hints. `collect()` walks the
  workdir for any `<output_basename>*.trees` (`Native`, "SLiM
  tree sequence") and `<output_basename>*.log` (`Log`) the script
  emitted — the adapter doesn't try to predict the exact filenames
  the script will write, since SLiM scripts choose their own
  output paths via `writeFile()` / `treeSeqOutput()` calls. Probe
  via `find_on_path(&["slim"])` (conda-forge / source / Homebrew
  all install under the canonical lowercase `slim` name); the
  generic version detector tries both the conventional `--version`
  and the SLiM-native `-version` form. Version range
  `4.0.0..5.0.0` (the modern release line is the 4.x series from
  2022+; 4.0 introduced the streamlined Eidos type system and the
  `treeSeqOutput()` helpers we rely on for tskit interop; upper
  bound 5.0 reserves room for an eventual major bump).
  `bio.slim.simulate` ribbon capability.
- **msprime** — Jerome Kelleher's coalescent backwards-in-time
  population-genetics simulator (GPL-3.0). msprime simulates the
  ancestry of a sample under a configurable demography and
  recombination map, then layers mutations onto the resulting tree
  sequence. It is the speed-of-light coalescent simulator
  (millions of samples per minute on a workstation) and the
  canonical companion to SLiM (forward-time) and tskit (tree-
  sequence analysis). Python-script subprocess shape (sister to
  Phase 17 Biopython): the user authors a `simulate.py` referenced
  from `[bio.msprime].script` in `case.toml`. Schema knobs:
  `script` (path to user-authored Python script; required),
  `python` (interpreter name; default `"python3"`),
  `population_size` (`u32`, ≥ 1), `num_samples` (`u32`, ≥ 1),
  `recombination_rate` (`f64`, ≥ 0.0 and finite — per-site per-
  generation rate), `mutation_rate` (`f64`, ≥ 0.0 and finite —
  per-site per-generation rate), `output_basename` (filename stem;
  required, non-empty). `prepare()` stages the script into the
  workdir under its original filename and writes a flat
  `valenx_params.json` containing `population_size`, `num_samples`,
  `recombination_rate` (emitted via `{:e}` so Python's `json.load`
  parses it back as a float), `mutation_rate` (same), and
  `output_basename`. `run()` invokes `python <script>` via the
  shared subprocess runner. `collect()` walks the workdir for
  `<output_basename>.trees` (`Native`, "msprime tree sequence"),
  `<output_basename>.vcf` (`Tabular`, "msprime VCF"), and
  `<output_basename>.csv` (`Tabular`, "msprime per-sample
  summary") — user scripts emit any combination of these via
  msprime / tskit's tabular APIs. Probe via Python on PATH with an
  `import msprime` check — when the import fails the probe still
  returns `ok = true` with a warning so users with msprime
  installed under a different interpreter (referenced via the case-
  level `python` override) aren't blocked. Version range
  `1.3.0..2.0.0` (the modern `sim_ancestry()` / `sim_mutations()`
  split landed in 1.3 in 2024, paired with the tskit 0.5+ tree-
  sequence format we surface in collect(); upper bound 2.0
  reserves room for an eventual major bump). `bio.msprime.simulate`
  ribbon capability.
- **tskit** — the canonical tree-sequence analysis library
  (MIT), built around the succinct tree-sequence data structure
  pioneered by msprime. tskit computes population-genetics
  statistics (π, Tajima's D, Fst, site-frequency spectra,
  IBD shares), exposes per-tree iteration across the genome,
  converts between tree-sequence and VCF / Newick / table formats,
  and renders phylogenetic plots. It's the workhorse downstream of
  every Phase 29 simulator — msprime emits `.trees`, SLiM emits
  `.trees`, tskit consumes them. Python-script subprocess shape
  (sister to msprime). Schema knobs: `script` (path to
  user-authored Python script; required), `python` (interpreter
  name; default `"python3"`), `input_trees` (`.trees` file from
  SLiM or msprime; required), `output_basename` (filename stem;
  required, non-empty). `prepare()` stages script + tree-sequence
  file into the workdir under their original filenames so the
  script can resolve them via relative paths, then writes a flat
  `valenx_params.json` containing `input_trees` (staged filename)
  and `output_basename`. `run()` invokes `python <script>` via the
  shared subprocess runner. `collect()` walks the workdir for
  `<output_basename>*.csv` / `<output_basename>*.tsv` (`Tabular`,
  "tskit statistics") and `*.png` (`Native`, "tskit plot") — user
  scripts emit any combination of statistics tables and rendered
  plots. Probe via Python on PATH with an `import tskit` check —
  same `ok = true` + warning fallback as msprime. Version range
  `0.5.0..1.0.0` (tskit 0.5+ ships the modern `Statistics` API
  surface and the v3 tree-sequence file format msprime 1.3+ writes;
  upper bound 1.0 reserves room for the long-promised 1.0
  release). `bio.tskit.analyze` ribbon capability.

### Canonical types

**No new canonical types.** All three adapters consume user-
supplied inputs (SLiM `.slim` Eidos scripts; msprime + tskit
Python scripts; tskit input `.trees` tree-sequence files emitted
by SLiM or msprime) and emit user-readable artifacts (`.trees`
tree sequences, `.vcf` genotype calls, `.csv` / `.tsv` statistics
tables, `.png` rendered plots, `.log` run logs) that the unchanged
`Results.artifacts` collection model surfaces directly. A first-
class population-genetics canonical type — a typed tree-sequence
representation spanning all three back-ends, with parsed per-tree
edge / node / mutation tables and a typed statistics-table
representation — defers to a future phase along with tree-sequence
visualizers and per-population allele-frequency-spectrum viewers.

### Headless CLIs

**No new CLIs.** SLiM's `.trees` outputs and `.log` files,
msprime's `.trees` / `.vcf` / `.csv` outputs, and tskit's `.csv`
/ `.tsv` statistics tables and `.png` plots are all standard
formats inspectable in any editor or through the user's downstream
Python pipeline (`pandas`, `numpy`, `tskit`, `msprime`); VCF
outputs are inspectable through the existing Phase 19
`valenx-vcf-info` CLI. A canonical population-genetics CLI —
tree-sequence diffing, per-tree statistics extraction, demographic-
inference pipeline driver — defers to a future phase along with
tree-sequence visualizers.

## Domain milestone

Phase 29 is the **first population-genetics / evolutionary-
simulation domain** to land in Valenx. The biology adapter family
started with Phase 17 (foundation — sequence / structure /
trajectory canonical types + classical MD + cheminformatics
scripts) and expanded through Phase 17.5 / 18 / 18.5 / 18.6 / 19 /
19.5 / 20 / 22 / 23 / 24 / 25 / 27 / 27.5 / 27.6 / 28 / 30 / 31 /
32 / 34 / 35 / 36 / 38 to cover sequence prediction, alignment,
RNA-seq, variant calling, single-cell, transcript quantification,
workflow orchestration, molecular viewers, cheminformatics,
quantum chemistry, protein design, EvolutionaryScale models, RNA
structure, phylogenetics, sequencing read simulation, systems
biology, small-molecule docking, CRISPR design, cryo-EM
reconstruction, and Rosetta protein modeling — but until Phase 29
the population-genetics surface (forward-time individual-based
evolutionary simulation, coalescent backward-time ancestry
simulation, tree-sequence analysis) was absent. Phase 29 closes
that gap with three established open-source tools spanning the
population-genetics tradeoff space — SLiM at the forward-time
individual-based end, msprime as the canonical coalescent
companion, and tskit as the analysis workhorse downstream of both
simulators.

## What landed early

The implementation rode subagent-driven-development across 5
discrete implementation commits (3 adapters, the registry rollup,
the init-template extension) plus this docs pass — each landing one
adapter, the registry rollup, the init-template extension, or the
documentation pass. Every commit kept workspace clippy + rustdoc
clean.

## Acceptance checklist

- [x] `valenx-adapter-slim` adapter ships with case-input parser
      + 4 lib tests + 4 case-input tests covering parses-minimal /
      parses-with-seed / rejects-empty-script / rejects-empty-
      output-basename / rejects-negative-seed, plus the
      single-binary subprocess shape that composes `slim
      [-s <seed>] [extras...] <script>` with the script positional
      last so SLiM treats it as the model file rather than the
      value of an earlier flag
- [x] `valenx-adapter-msprime` adapter ships with case-input
      parser + 4 lib tests + 4 case-input tests covering parses-
      minimal / rejects-empty-script / rejects-zero-population-
      size / rejects-negative-recombination-rate / rejects-non-
      finite-mutation-rate, plus the Python-script subprocess
      shape that stages script + writes `valenx_params.json` with
      the parsed demographic knobs (`population_size`,
      `num_samples`, `recombination_rate` + `mutation_rate`
      emitted as `{:e}` so Python `json.load` parses them as
      floats, `output_basename`)
- [x] `valenx-adapter-tskit` adapter ships with case-input parser
      + 4 lib tests + 3 case-input tests covering parses-minimal
      / rejects-empty-script / rejects-empty-input-trees /
      rejects-empty-output-basename, plus the Python-script
      subprocess shape that stages script + input `.trees` file +
      writes `valenx_params.json` with the staged tree-sequence
      filename so the script can resolve it via the bare filename
- [x] All 3 adapters wired into `valenx-app::init_registry` —
      live adapter count moves from 84 to **89** (alongside the
      Phase 38 Rosetta-family pair), opening the first population-
      genetics / evolutionary-simulation domain to ship in Valenx
- [x] 3 population-genetics templates in `valenx-init` (`slim`,
      `msprime`, `tskit`), all round-tripping through
      `valenx-validate` (cross-binary roundtrip now sweeps **85
      templates** clean)
- [x] STATUS.md / README.md / ROADMAP.md / CHANGELOG.md /
      QUICKSTART.md updated
- [ ] **Phase 29.5** — sister-adapter expansion of Phase 29:
      fwdpy11 (Kevin Thornton's Python-driven forward-time
      population-genetics simulator; sister to SLiM with a
      different scripting surface; defer to sister-adapter
      expansion phase), simuPOP (Bo Peng's Python forward
      simulator with a long history; defer), pyslim (the SLiM
      tree-sequence Python interop layer that bridges SLiM
      `.trees` outputs into msprime / tskit; defer to 29.5),
      stdpopsim (the standardised population-genetics simulation
      library that wraps msprime / SLiM under a catalog of
      curated demographic models; defer), demes (the human-
      readable demographic-model specification format; defer to
      29.5). Out of scope for this beachhead.

## Success metrics

| Metric                                                | Target          |
|-------------------------------------------------------|-----------------|
| New population-genetics adapter (template + tests)    | 1 day per       |
| Forward-sim → coalescent-sim → analyze loop across 3 tools | < tool baseline |

## Leads into

Phase 29 opens the population-genetics / evolutionary-simulation
domain that the user's bio / chemistry spec called out alongside
the Phase 17 / 17.5 biology + structure-prediction stack and the
Phase 30 phylogenetics beachhead. Combined with the existing
design-guides → predict-off-targets → search-off-targets →
simulate-reads → align → quantify → call-variants → predict-
structure → fold-RNA → infer-tree → simulate-pathway → reconstruct-
3D → validate loop, the **forward-sim → coalescent-sim → analyze-
trees → simulate-reads → align → quantify → call-variants →
predict-structure → fold-RNA → infer-tree → simulate-pathway →
reconstruct-3D → design-protein → validate** loop now spans three
population-genetics tools (SLiM, msprime, tskit) feeding into the
existing Phase 31 read simulators (ART, wgsim, Badread), the eleven
Phase 18 / 18.5 / 18.6 alignment tools (BWA, Bowtie2, HISAT2,
STAR, minimap2, MAFFT, MUSCLE, HMMER, samtools, MMseqs2, DIAMOND),
the two Phase 20 transcript quantifiers (Salmon, Kallisto), the
three Phase 19 variant callers (bcftools, GATK, DeepVariant), the
Phase 17 / 17.5 prediction stack (ESMFold, OpenFold, AlphaFold
2/3, ColabFold), the Phase 28 RNA-structure tools (ViennaRNA,
RNAstructure, NUPACK), the Phase 30 phylogenetic-tree builders
(IQ-TREE, RAxML-NG, FastTree), the Phase 32 systems-biology
surface (COPASI, BioNetGen, PhysiCell), the Phase 35 CRISPR-design
tools (CHOPCHOP, CRISPOR, Cas-OFFinder), the Phase 36 cryo-EM
reconstruction tools (RELION, EMAN2, CTFFIND), and the Phase 38
Rosetta-family adapters (Rosetta, PyRosetta) — all in one Valenx
shell with no glue code beyond the existing case-toml / prepare /
run / collect path.

The natural follow-up is **Phase 29.5** — the deferred population-
genetics work called out above (fwdpy11 as a Python-driven
forward-time simulator sister to SLiM, simuPOP for the long-
history Python forward-simulator surface, pyslim as the SLiM
tree-sequence Python interop layer bridging SLiM `.trees` outputs
into msprime / tskit, stdpopsim as the standardised population-
genetics simulation library wrapping msprime / SLiM under a
catalog of curated demographic models, demes as the human-readable
demographic-model specification format), slotting in alongside the
existing Phase 29 adapters with the same single-binary subprocess
shape (SLiM sister tools) or Python-script subprocess shape
(msprime / tskit sister tools). See the out-of-scope section of
`docs/superpowers/plans/2026-04-30-population-genetics.md` for
the full follow-up phase list.
