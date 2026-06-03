# Phase 30 — Phylogenetics

**Status:** 🟢 Live — IQ-TREE + RAxML-NG + FastTree open the molecular
phylogenetics domain in Valenx alongside the Phase 17 / 17.5 prediction
stack and the Phase 18 / 18.5 / 18.6 / 20 alignment + quantification
beachhead.

## Goal

Open the molecular phylogenetics domain in Valenx with the three most-
used maximum-likelihood tree inference tools: **IQ-TREE** (the de-facto
modern ML tree builder — ModelFinder + UFBoot ultrafast bootstrap),
**RAxML-NG** (the next-generation RAxML rewrite — successor to
classical `raxmlHPC`), and **FastTree** (approximate-ML, optimized for
very large trees — sub-quadratic in alignment size). All three follow
the established Phase 18 BWA single-binary CLI pattern: input alignment
in, tree out. No two-stage index step (the alignment is the input; the
tree is the output). Phase 30 sits numerically after Phase 27.5 and
ships chronologically right after the Phase 27.5 protein-design
expansion beachhead.

## Capability inventory

### Live adapters (3)

- **IQ-TREE** — Bui Quang Minh & Robert Lanfear's de-facto modern
  maximum-likelihood phylogenetic tree builder (GPL-2.0). Single-binary
  subprocess shape; alignment in, tree out. Schema knobs: `alignment`
  (FASTA / PHYLIP / NEXUS / CLUSTAL; required), `model` (default
  `"MFP"` — `"TEST"` / `"MFP"` trigger ModelFinder's automatic model
  selection; otherwise pass e.g. `"GTR+G"` / `"WAG+I+G"` verbatim;
  required non-empty), `bootstrap` (UFBoot ultrafast bootstrap
  replicates; `0` disables, default `1000`), `threads` (default
  `"AUTO"` — IQ-TREE's auto-detect, otherwise an integer count;
  validated against `^(AUTO|\d+)$`), `prefix` (output file prefix;
  required non-empty), `extra_args`. `prepare()` builds `iqtree2 -s
  <alignment> -m <model> -B <bootstrap> -T <threads> --prefix <prefix>
  [extras...]` (omitting `-B` when `bootstrap == 0`). `collect()`
  walks for `<prefix>.treefile` (`Native`, `"IQ-TREE ML tree"`),
  `<prefix>.iqtree` (`Log`), `<prefix>.log` (`Log`). Probe via
  `find_on_path(&["iqtree2", "iqtree"])` (newer 2.x ships as
  `iqtree2`; older 1.x as `iqtree`). `bio.iqtree.tree` ribbon
  capability.
- **RAxML-NG** — Alexey Kozlov's next-generation RAxML rewrite (AGPL-
  3.0). Successor to classical `raxmlHPC`. Single-binary subprocess
  shape with mode dispatch. Schema knobs: `alignment` (required),
  `model` (substitution model — `"GTR+G"` / `"WAG+I+G"` / etc.;
  required non-empty), `mode` (`"search"` single-tree ML, `"all"`
  search + bootstrap, or `"bootstrap"` bootstrap-only on existing
  tree), `bootstrap` (replicates — required ≥ 1 when `mode ∈ {all,
  bootstrap}`, ignored otherwise), `threads` (≥ 1, default 1),
  `prefix` (required non-empty), `extra_args`. `prepare()` builds
  `raxml-ng --<mode> --msa <alignment> --model <model> --threads <N>
  --prefix <prefix> [--bs-trees <bootstrap> if mode in {all,
  bootstrap}] [extras...]`. `collect()` walks for
  `<prefix>.raxml.bestTree` (`Native`, `"RAxML-NG ML tree"`),
  `<prefix>.raxml.support` (`Native`), `<prefix>.raxml.log` (`Log`).
  Probe via `find_on_path(&["raxml-ng"])`. `bio.raxml-ng.tree`
  ribbon capability.
- **FastTree** — Morgan Price's approximate-ML phylogenetic inference
  tool (GPL-2.0). Optimized for very large trees: sub-quadratic in
  alignment size. Single-binary subprocess shape; writes Newick to
  stdout (the MAFFT-style stdout-redirect pattern captures stdout to
  the `output` path). Schema knobs: `alignment` (required), `output`
  (Newick tree path; required), `seq_type` (`"nt"` nucleotide or
  `"aa"` amino acid), `use_gtr` (default `true` — uses GTR for
  nucleotides, ignored for amino acid; FastTree's default is JC
  without this flag), `gamma` (gamma rate-variation model toggle;
  default `false`), `extra_args`. `prepare()` builds — nucleotide:
  `FastTree [-nt] [-gtr if use_gtr] [-gamma if gamma] <alignment>` →
  stdout; amino-acid: `FastTree [-gamma if gamma] <alignment>` →
  stdout. `collect()` reports `output` as a `Native` artifact
  `"FastTree Newick tree"`. Probe via `find_on_path(&["FastTree",
  "fasttree"])` (binary name varies by distro). `bio.fasttree.tree`
  ribbon capability.

### Canonical types

**No new canonical types.** All three adapters consume the existing
Phase 18 / 20 multiple-sequence alignment inputs (FASTA / PHYLIP /
NEXUS / CLUSTAL) and emit Newick-format tree files that the unchanged
`Results.artifacts` collection model surfaces directly through the
`Native` artifact kind. A first-class `Tree` canonical type with a
Newick reader as a Valenx CLI defers to a future phase along with
visualization integrations.

### Headless CLIs

**No new CLIs.** Newick trees are short text files that any tree
viewer (FigTree, iTOL, Dendroscope) can ingest directly; the
adapter-emitted `<prefix>.iqtree` / `<prefix>.raxml.log` log files
are human-readable and don't need a dedicated inspector.

## What landed early

The implementation landed across 6 discrete
commits, each landing one adapter, the registry rollup, the init-
template extension, or the documentation pass. Every commit kept
workspace clippy + rustdoc clean.

## Acceptance checklist

- [x] `valenx-adapter-iqtree` adapter ships with case-input parser
      + 4 lib tests + 5 case-input tests
- [x] `valenx-adapter-raxml-ng` adapter ships with case-input parser
      + 4 lib tests + 5 case-input tests
- [x] `valenx-adapter-fasttree` adapter ships with case-input parser
      + 4 lib tests + 4 case-input tests
- [x] All 3 adapters wired into `valenx-app::init_registry` —
      live adapter count moves from 61 to 64
- [x] 3 phylogenetics templates in `valenx-init` (`iqtree` with
      alias `iqtree2`, `raxml-ng` with alias `raxml`, `fasttree`),
      all round-tripping through `valenx-validate` (cross-binary
      roundtrip now sweeps 60 templates clean)
- [x] STATUS.md / README.md / ROADMAP.md / CHANGELOG.md /
      QUICKSTART.md updated
- [ ] **Phase 30.5** — BEAST 2 / MrBayes / RevBayes (Bayesian
      phylogenetics — different shape, MCMC convergence-monitoring
      story), PhyML (niche; defer to user demand), ModelTest /
      jModelTest (model selection — workflow-orchestration concern,
      slot into the workflow-manager surface), tree visualization
      (FigTree, Dendroscope, TreeViewer — viewer concern, slot into
      Phase 23.5 if user demand surfaces). Out of scope for this
      beachhead.

## Success metrics

| Metric                                            | Target          |
|---------------------------------------------------|-----------------|
| New phylogenetics adapter (template + tests)      | 1 day per       |
| ML tree inference across 3 tools                  | < tool baseline |

## Leads into

Phase 30 opens the molecular phylogenetics domain that the user's bio
spec called out alongside the Phase 17 / 17.5 prediction stack and
the Phase 18 / 18.5 / 18.6 / 20 alignment + quantification beachhead.
Combined with the existing align → quantify → predict → validate
loop, the **align → quantify → predict → infer-tree → validate** loop
now spans eleven alignment / search tools (BWA, Bowtie2, HISAT2,
STAR, minimap2, MAFFT, MUSCLE, HMMER, samtools, MMseqs2, DIAMOND),
two transcript quantifiers (Salmon, Kallisto), five prediction tools
(ESMFold, OpenFold, AlphaFold 2, AlphaFold 3, ColabFold), and three
phylogenetic-tree builders (IQ-TREE, RAxML-NG, FastTree) — all in
one Valenx shell with no glue code beyond the existing case-toml /
prepare / run / collect path.

The natural follow-up is **Phase 30.5** — Bayesian phylogenetics
(BEAST 2 / MrBayes / RevBayes) with the MCMC convergence-monitoring
shape, plus model-selection orchestration (ModelTest / jModelTest)
that sits adjacent to Phase 22's workflow-manager surface.
