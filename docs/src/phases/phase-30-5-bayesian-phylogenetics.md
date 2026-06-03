# Phase 30.5 — Bayesian phylogenetics

**Status:** 🟢 Live — BEAST 2 + MrBayes round out the molecular
phylogenetics surface alongside the Phase 30 IQ-TREE + RAxML-NG +
FastTree maximum-likelihood beachhead.

## Goal

Sister-domain expansion of Phase 30. Add the two de-facto Bayesian
phylogenetic inference engines — **BEAST 2** (Bayesian Evolutionary
Analysis by Sampling Trees v2 — the cross-platform XML-driven MCMC
framework with a sprawling package ecosystem covering tip-dated
trees, relaxed molecular clocks, coalescent demographic models,
birth-death speciation models, and the BDSKY / MASCOT / BEASTling /
StarBEAST3 universe of extensions) and **MrBayes** (the long-
standing Bayesian MCMC tree inference engine that remains the
de-facto choice alongside BEAST 2 for posterior tree sampling
across nucleotide / amino-acid / morphological datasets, with its
own NEXUS-embedded model-and-mcmc command language and built-in
Metropolis-coupled MCMC ("MC^3") chain swapping). Both adapters
follow the established Phase 18 BWA single-binary CLI pattern: a
user-authored model description (BEAST 2 XML or MrBayes NEXUS
file) in, posterior tree + parameter samples out — the same
shape the Phase 30 ML tools use, with the inputs swapped from
multiple-sequence alignments to MCMC model files. Phase 30.5 sits
numerically after Phase 30 and ships chronologically right after
Phase 38 Rosetta — same chronological-vs-numerical convention used
for Phase 17.5 / 24 / 28 / 31 / 35.

## Capability inventory

### Live adapters (2)

- **BEAST 2** — the cross-platform Bayesian Evolutionary Analysis
  by Sampling Trees v2 engine (LGPL-2.1). BEAST 2 is the canonical
  Bayesian MCMC framework for time-calibrated phylogenetics: tip-
  dated trees, relaxed molecular clocks, coalescent demographic
  models, birth-death speciation models, and the ever-growing
  universe of BEAST 2 packages (BDSKY, MASCOT, BEASTling,
  StarBEAST3, ...). It complements the maximum-likelihood Phase 30
  tools (IQ-TREE, RAxML-NG, FastTree) with a full posterior over
  tree topologies and parameters. Single-binary subprocess shape
  (sister to Phase 18 BWA): the CLI is `beast [-seed <N>] -threads
  <N> [-overwrite] <xml> [extras...]`. The user authors / generates
  the model XML (typically through BEAUti) and references it from
  `[bio.beast2].xml` — the adapter doesn't generate XML. Schema
  knobs: `xml` (BEAUti-generated XML model file; required), `seed`
  (optional `u64`; passed via `beast -seed <N>` when present,
  otherwise BEAST picks its own seed and prints it on the run
  banner), `threads` (number of threads BEAST uses for tree-
  likelihood evaluation; `u32`, ≥ 1, default 1), `overwrite`
  (default `false`; toggles `-overwrite` so an existing output set
  from a previous run is replaced rather than triggering a fail),
  `extra_args`. `prepare()` resolves the XML against the case
  directory when relative, validates it exists on disk, composes
  the invocation with `seed` injected before `-threads` and the XML
  positional last so BEAST treats it as the model file rather than
  the value of an earlier flag. `run()` streams BEAST's
  `Random number seed` / `BEAST v2` startup banner / periodic
  `Sample` / `posterior` chain-status lines / `End likelihood` /
  `Total calculation time` end-of-run sentinels into progress
  hints. `collect()` walks the workdir for `*.log` (`Log`, "BEAST 2
  trace log" — the parameter trace Tracer reads) and `*.trees`
  (`Native`, "BEAST 2 sampled trees" — the sampled tree posterior
  TreeAnnotator / DensiTree consumes); the adapter doesn't try to
  predict the exact filenames since BEAST writes whatever the XML's
  `<log fileName="...">` sites configure. Probe via
  `find_on_path(&["beast"])`; the generic version detector tries
  the conventional `--version` and BEAST's own `-version` form.
  Version range `2.7.0..3.0.0` (the modern stable line is the 2.7.x
  series from 2022+ that introduced modern threading + the package
  manager; upper bound 3.0 reserves room for an eventual major
  bump). `bio.beast2.mcmc` ribbon capability.
- **MrBayes** — the long-standing Bayesian MCMC phylogenetic
  inference engine (GPL-3.0). MrBayes is the historic workhorse for
  Bayesian phylogenetics: alongside BEAST 2 it remains the de-facto
  choice for posterior tree sampling across nucleotide / amino-
  acid / morphological datasets, with its own NEXUS-embedded
  model-and-mcmc command language and built-in Metropolis-coupled
  MCMC ("MC^3") chain swapping. Single-binary subprocess shape
  (sister to BEAST 2): the CLI is `mb [-i] <nexus> [extras...]`.
  The binary is literally named `mb` (the project's own
  convention). The user authors a NEXUS file with a DATA block plus
  a MRBAYES block embedding the model / MCMC parameters and `mcmc`
  command and references it from `[bio.mrbayes].nexus`. Schema
  knobs: `nexus` (NEXUS data file with embedded MRBAYES block;
  required), `batch` (default `false`; toggles `-i` so MrBayes runs
  the embedded commands non-interactively and exits cleanly rather
  than waiting on stdin at the prompt — the right default for non-
  interactive automation), `extra_args`. `prepare()` resolves the
  NEXUS path against the case directory when relative, validates it
  exists on disk, and composes the invocation with the NEXUS
  positional last so MrBayes treats it as the model file rather
  than the value of an earlier flag. `run()` streams MrBayes's
  `MrBayes v` / `Initializing` startup banner / periodic
  `Generation NNNN` / `Avg standard deviation of split frequencies`
  chain-status lines / `Analysis completed` / `Continue with
  analysis` end-of-run sentinels into progress hints. `collect()`
  walks the workdir for `*.t` (`Native`, "MrBayes tree samples"),
  `*.p` (`Tabular`, "MrBayes parameter samples"), and `*.con.tre`
  (`Native`, "MrBayes consensus tree"). Probe via
  `find_on_path(&["mb"])`. Version range `3.2.0..4.0.0` (the long-
  running stable 3.2.x line that every distro ships covers every
  release through 3.2.7; upper bound 4.0 reserves room for an
  eventual major bump). `bio.mrbayes.mcmc` ribbon capability.

### Canonical types

**No new canonical types.** Both adapters consume user-supplied
inputs (BEAST 2 XML model files, MrBayes NEXUS data files with
embedded MRBAYES blocks) and emit user-readable artifacts (BEAST 2
`.log` trace files + `.trees` sampled-tree posteriors, MrBayes
`.t` tree samples + `.p` parameter samples + `.con.tre` consensus
trees) that the unchanged `Results.artifacts` collection model
surfaces directly. A first-class Bayesian-phylogenetics canonical
type — a typed posterior representation spanning both back-ends,
with parsed per-generation parameter traces and per-sample tree
topologies plus convergence-diagnostic helpers (effective sample
size, Gelman-Rubin) — defers to a future phase along with trace
visualizers, tree-density plots, and consensus-tree viewers.

### Headless CLIs

**No new CLIs.** BEAST 2's `.log` parameter traces and `.trees`
sampled-tree posteriors plus MrBayes's `.t` / `.p` / `.con.tre`
families are all standard formats inspectable in the canonical
downstream tools (Tracer for `.log`, DensiTree / TreeAnnotator /
FigTree for `.trees` / `.t` / `.con.tre`, R `coda` / Python
`arviz` for parameter traces). A canonical Bayesian-phylogenetics
CLI — convergence-diagnostic computation, posterior-tree
summarisation, log-file diffing — defers to a future phase along
with trace visualizers and consensus-tree viewers.

## Domain milestone

Phase 30.5 closes the Bayesian half of the molecular phylogenetics
surface that Phase 30 opened from the maximum-likelihood side.
Phase 30 (IQ-TREE + RAxML-NG + FastTree) covered the ML tradeoff
space — modern ML with ModelFinder + UFBoot bootstrap, the
next-generation RAxML rewrite, and approximate-ML for very large
trees. Phase 30.5 (BEAST 2 + MrBayes) covers the Bayesian MCMC
side: BEAST 2 for time-calibrated trees with relaxed clocks,
demographic priors, and the sprawling package ecosystem; MrBayes
for the historic NEXUS-embedded MCMC workhorse with built-in MC^3
chain swapping. Together the five tools span the full
phylogenetics tradeoff space — point-estimate ML at the fast-and-
approximate end (FastTree), modern ML with bootstrap (IQ-TREE,
RAxML-NG), and full posterior MCMC (BEAST 2, MrBayes) at the
expensive-but-rigorous end.

## What landed early

The implementation landed across 4
discrete implementation commits (2 adapters, the registry rollup,
the init-template extension) plus this docs pass — each landing
one adapter, the registry rollup, the init-template extension, or
the documentation pass. Every commit kept workspace clippy +
rustdoc clean.

## Acceptance checklist

- [x] `valenx-adapter-beast2` adapter ships with case-input parser
      + 4 lib tests + 4 case-input tests covering parses-minimal /
      parses-with-seed / rejects-empty-xml / rejects-zero-threads /
      rejects-bad-extra-args, plus the single-binary subprocess
      shape that composes `beast [-seed <N>] -threads <N>
      [-overwrite] <xml> [extras...]` with the XML positional last
      so BEAST treats it as the model file rather than the value
      of an earlier flag
- [x] `valenx-adapter-mrbayes` adapter ships with case-input
      parser + 4 lib tests + 3 case-input tests covering parses-
      minimal / parses-with-batch / rejects-empty-nexus / rejects-
      bad-extra-args, plus the single-binary subprocess shape that
      composes `mb [-i if batch] <nexus> [extras...]` with the
      NEXUS positional last
- [x] Both adapters wired into `valenx-app::init_registry` —
      live adapter count moves from 89 to **94** (alongside the
      Phase 39 DNA structural geometry trio), rounding out the
      molecular phylogenetics surface that Phase 30 opened from
      the ML side
- [x] 2 Bayesian-phylogenetics templates in `valenx-init`
      (`beast2` with alias `beast`, `mrbayes` with alias `mb`),
      all round-tripping through `valenx-validate` (cross-binary
      roundtrip now sweeps **90 templates** clean)
- [x] STATUS.md / README.md / ROADMAP.md / CHANGELOG.md /
      QUICKSTART.md updated
- [ ] **Phase 30.6** — sister-adapter expansion of Phase 30.5:
      RevBayes (Sebastian Höhna's flexible probabilistic graphical-
      model phylogenetics framework with a Rev language; sister to
      BEAST 2 / MrBayes with a different scripting surface; defer
      to sister-adapter expansion phase), BEAST 1.x (the long-
      lived predecessor / sibling line to BEAST 2; defer), PhyloBayes
      (Nicolas Lartillot's site-heterogeneous CAT model Bayesian
      phylogenetics; defer), MIGRATE-N (Peter Beerli's Bayesian
      population-genetics inference of migration rates and effective
      population sizes; sits adjacent to Phase 29 population
      genetics; defer). Out of scope for this beachhead.

## Success metrics

| Metric                                                | Target          |
|-------------------------------------------------------|-----------------|
| New Bayesian-phylogenetics adapter (template + tests) | 1 day per       |
| Bayesian MCMC tree-inference loop across 2 tools      | < tool baseline |

## Leads into

Phase 30.5 rounds out the molecular phylogenetics surface that
Phase 30 opened from the maximum-likelihood side. Combined with
the existing align → quantify → predict → infer-tree → validate
loop, the **align → quantify → predict → infer-tree-ML → infer-
tree-Bayesian → validate** loop now spans eleven alignment / search
tools (BWA, Bowtie2, HISAT2, STAR, minimap2, MAFFT, MUSCLE, HMMER,
samtools, MMseqs2, DIAMOND), two transcript quantifiers (Salmon,
Kallisto), five prediction tools (ESMFold, OpenFold, AlphaFold 2,
AlphaFold 3, ColabFold), three maximum-likelihood phylogenetic-
tree builders (IQ-TREE, RAxML-NG, FastTree), and two Bayesian MCMC
tree-inference engines (BEAST 2, MrBayes) — all in one Valenx
shell with no glue code beyond the existing case-toml / prepare /
run / collect path.

The natural follow-up is **Phase 30.6** — the deferred Bayesian-
phylogenetics work called out above (RevBayes for the flexible
probabilistic graphical-model framework with a Rev language, BEAST
1.x for the long-lived predecessor / sibling line to BEAST 2,
PhyloBayes for the site-heterogeneous CAT model, MIGRATE-N for
the Bayesian population-genetics inference of migration rates that
sits adjacent to the Phase 29 population-genetics surface),
slotting in alongside the existing BEAST 2 + MrBayes adapters with
the same single-binary subprocess shape (RevBayes / BEAST 1.x /
PhyloBayes / MIGRATE-N sister tools).
