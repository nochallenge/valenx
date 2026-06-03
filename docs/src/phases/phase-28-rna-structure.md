# Phase 28 — RNA structure

**Status:** 🟢 Live — ViennaRNA + RNAstructure + NUPACK open the RNA
secondary-structure-prediction domain in Valenx alongside the Phase 17
biology stack and the Phase 17.5 / 27 / 27.5 / 30 protein-design,
prediction, and phylogenetics beachheads.

## Goal

Open the RNA structure prediction domain in Valenx with three
established tools: **ViennaRNA** (the most-cited RNA secondary-structure
suite — `RNAfold` minimum-free-energy folding), **RNAstructure**
(Mathews lab's classic toolkit, BSD-3-Clause licensed; `Fold` is the
flagship), and **NUPACK** (Caltech's nucleic-acid package — academic-
license-only, surfaces an awareness warning à la VMD / AlphaFold 3).
ViennaRNA follows the MAFFT-style stdout-redirect pattern (RNAfold
writes to stdout); RNAstructure follows the BWA single-binary CLI
pattern with explicit `-o`-style output; NUPACK follows the OpenMM /
Scanpy Python-script-subprocess pattern (NUPACK 4 is Python-driven; the
3.x CLI is deprecated). Phase 28 sits numerically between Phase 27.5
and Phase 30 and ships chronologically right after the Phase 30
phylogenetics beachhead — the same chronological-vs-numerical
convention as Phase 17.5 sits numerically between Phase 17 and Phase 18.

## Capability inventory

### Live adapters (3)

- **ViennaRNA** — the most-cited RNA secondary-structure suite (Ivo
  Hofacker et al., custom non-commercial / academic-use license). Single-
  binary subprocess shape with stdout-redirect (RNAfold writes the
  dot-bracket structure to stdout; the MAFFT-style stdout-capture
  pattern routes it to `output`). Schema knobs: `input` (FASTA file
  containing the sequence(s) to fold; required), `output` (dot-bracket
  output filename relative to workdir; required), `temperature`
  (Celsius; default 37.0; finite), `partition_function` (default false
  — toggles `-p` for partition function + base-pair probabilities),
  `allow_gu` (default true — `--noGU` disables GU pairs), `extra_args`.
  `prepare()` builds `RNAfold -i <input> -T <temperature> [-p]
  [--noGU if !allow_gu] [extras...]` → stdout captured to `output`.
  `collect()` reports `output` as a `Native` artifact "ViennaRNA
  secondary structure". Probe via `find_on_path(&["RNAfold"])` (capital
  R-N-A; that's ViennaRNA's binary name). **License callout:**
  ViennaRNA is licensed for non-commercial / academic use only — the
  probe surfaces an `"academic"`-keyworded license-awareness warning
  in `ProbeReport.warnings` reminding users to confirm their use case
  before redistributing folds or derived data. The init aliases
  `vienna` and `rnafold` resolve to the same template.
  `bio.viennarna.fold` ribbon capability.
- **RNAstructure** — Mathews lab's classic RNA folding toolkit (BSD-3-
  Clause). Single-binary subprocess shape with `-o`-style output
  (binary literally named `Fold`, mirroring the capital-F naming of
  `FastTree` and `STAR`). Schema knobs: `input` (FASTA or `.seq`
  RNAstructure-native format; required), `output` (`.ct` connection-
  table file; required), `max_structures` (number of structures to
  predict; default 20; ≥ 1), `max_percent` (energy difference cutoff
  as % of MFE; default 10; in `0..=100`), `temperature` (Kelvin —
  RNAstructure's convention; default 310.15; finite, > 0.0),
  `extra_args`. `prepare()` builds `Fold <input> <output> -m
  <max_structures> -p <max_percent> -t <temperature> [extras...]`.
  `collect()` reports `output` as a `Native` artifact "RNAstructure
  connectivity table". Probe via `find_on_path(&["Fold"])`.
  `bio.rnastructure.fold` ribbon capability.
- **NUPACK** — Caltech's nucleic-acid package (Niles Pierce lab,
  custom academic-only license). Python-script subprocess shape:
  NUPACK 4 is Python-driven (the 3.x CLI is deprecated), so the user
  supplies a Python script that imports `nupack` and reads
  `valenx_params.json` for the config knobs. Schema knobs: `script`
  (required Python file), `python` (default `python3`), `input`
  (optional FASTA / `.npc` NUPACK config), `output_basename` (script
  reads from `valenx_params.json` and writes outputs prefixed with
  this; required non-empty), `temperature` (Celsius; default 37.0;
  finite), `sodium` (salt concentration in molar — NUPACK's `sodium`
  parameter; default 1.0; > 0.0 and finite). `prepare()` stages the
  script + optional input, writes `valenx_params.json` with the staged
  filename / `output_basename` / `temperature` / `sodium`, builds
  `native_command = [python, script]`. `collect()` walks
  `<output_basename>*` (`Native`, "NUPACK output") and `.npc` /
  `.json` files (`Tabular` / `Log`). Probe via
  `find_on_path(&["python3", "python"])` then
  `python -c "import nupack; print(nupack.__version__)"` — surfaces
  an install hint when Python is on PATH but `nupack` isn't
  importable. **License callout:** NUPACK is licensed for non-
  commercial / academic use only — the probe surfaces an
  `"academic"`-keyworded license-awareness warning in
  `ProbeReport.warnings` reminding users that Caltech's NUPACK license
  restricts redistribution + commercial use; confirm your use case
  complies before publishing analyses. `bio.nupack.analyze` ribbon
  capability.

### Canonical types

**No new canonical types.** All three adapters consume the existing
Phase 17 FASTA inputs and emit user-readable artifacts (dot-bracket
structures, connection tables, NUPACK output files) that the
unchanged `Results.artifacts` collection model surfaces directly. A
first-class RNA-secondary-structure canonical type with a dot-bracket
or `.ct` reader as a Valenx CLI defers to a future phase along with
visualization integrations.

### Headless CLIs

**No new CLIs.** Dot-bracket structures and `.ct` connection tables
are short text files that any RNA structure viewer (RNAcanvas, VARNA,
forna) can ingest directly; NUPACK script outputs are user-supplied
shapes that the existing `Results.artifacts` collection covers.

## License callouts

Two of the three adapters wrap tools that ship under non-commercial /
academic-use-only licenses. Both surface a license-awareness warning
through their `probe()` call so the user sees it in the registry
status before they ship folds or analyses downstream:

- **ViennaRNA** — custom non-commercial license. Free for academic
  and non-commercial use; commercial use requires a separate license
  from the University of Vienna. The probe pushes an
  `"academic"`-keyworded warning into `ProbeReport.warnings` whenever
  `RNAfold` is on PATH.
- **NUPACK** — custom Caltech academic-only license. Free for non-
  commercial / academic use; commercial use requires a separate
  license from Caltech's Office of Technology Transfer. The probe
  pushes an `"academic"`-keyworded warning into `ProbeReport.warnings`
  whenever Python (and the importable `nupack` module) are present.

This mirrors the existing VMD (Phase 23) / AlphaFold 3 (Phase 17.5) /
ChimeraX (Phase 27.5 expansion) probe-warning pattern. RNAstructure
ships under BSD-3-Clause and needs no analogous callout.

## What landed early

The implementation landed across 6 discrete
commits, each landing one adapter, the registry rollup, the init-
template extension, or the documentation pass. Every commit kept
workspace clippy + rustdoc clean.

## Acceptance checklist

- [x] `valenx-adapter-viennarna` adapter ships with case-input parser
      + 4 lib tests + 4 case-input tests + 1 probe-warning test
      asserting the `"academic"` keyword surfaces
- [x] `valenx-adapter-rnastructure` adapter ships with case-input
      parser + 4 lib tests + 4 case-input tests
- [x] `valenx-adapter-nupack` adapter ships with case-input parser
      + 4 lib tests + 4 case-input tests + 1 probe-warning test
      asserting the `"academic"` keyword surfaces
- [x] All 3 adapters wired into `valenx-app::init_registry` —
      live adapter count moves from 64 to 67
- [x] 3 RNA-structure templates in `valenx-init` (`viennarna` with
      aliases `vienna` / `rnafold`, `rnastructure`, `nupack`), all
      round-tripping through `valenx-validate` (cross-binary
      roundtrip now sweeps 63 templates clean)
- [x] ViennaRNA's + NUPACK's academic-only licenses surface in the
      probe warnings so users see them before they ship folds or
      analyses downstream
- [x] STATUS.md / README.md / ROADMAP.md / CHANGELOG.md /
      QUICKSTART.md updated
- [ ] **Phase 28.5** — ContraFold / IPknot / ProbKnot (sub-tools of
      RNA suites — niche enough to defer until user demand surfaces),
      LocARNA (alignment-based RNA structure prediction — different
      shape, separate phase), SimRNA (3D RNA structure — different
      shape, would slot alongside the Phase 17.5 protein-prediction
      stack), mfold / UNAFold (predecessor to RNAstructure /
      superseded — defer). Out of scope for this beachhead.

## Success metrics

| Metric                                            | Target          |
|---------------------------------------------------|-----------------|
| New RNA-structure adapter (template + tests)      | 1 day per       |
| RNA secondary-structure prediction across 3 tools | < tool baseline |

## Leads into

Phase 28 opens the RNA structure prediction domain that the user's
bio spec called out alongside the Phase 17 / 17.5 protein-prediction
stack and the Phase 18 / 18.5 / 18.6 / 20 alignment + quantification
beachhead. Combined with the existing fold → analyze → predict →
infer-tree → validate loop, the **align → quantify → predict →
fold-RNA → infer-tree → validate** loop now spans eleven alignment /
search tools (BWA, Bowtie2, HISAT2, STAR, minimap2, MAFFT, MUSCLE,
HMMER, samtools, MMseqs2, DIAMOND), two transcript quantifiers
(Salmon, Kallisto), five prediction tools (ESMFold, OpenFold,
AlphaFold 2, AlphaFold 3, ColabFold), three RNA-structure tools
(ViennaRNA, RNAstructure, NUPACK), and three phylogenetic-tree
builders (IQ-TREE, RAxML-NG, FastTree) — all in one Valenx shell
with no glue code beyond the existing case-toml / prepare / run /
collect path.

The natural follow-up is **Phase 28.5** — the deferred RNA-structure
work called out above (ContraFold, IPknot, ProbKnot, LocARNA, SimRNA,
mfold / UNAFold), slotting in alongside the existing RNA-structure
adapters with the same single-binary or Python-script-subprocess
shape.
