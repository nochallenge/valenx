# Phase 35.5 — Base + prime editing design

**Status:** 🟢 Live — BE-Designer + BE-Hive + PrimeDesign + pegFinder
round out the **CRISPR editing-design surface** that Phase 35
CHOPCHOP / CRISPOR / Cas-OFFinder opened, alongside the Phase 5.5
/ 5.6 / 5.7 / 17 / 17.5 / 17.7 / 18 / 18.5 / 18.6 / 18.7 / 19 /
19.5 / 19.6 / 20 / 22 / 22.5 / 23 / 24 / 25 / 27 / 27.5 / 27.6 /
28 / 29 / 30 / 30.5 / 31 / 32 / 32.5 / 33 / 34 / 36 / 38 / 39 / 40
/ 41 / 42 / 43 / 44.5 biology / biotech / chemistry beachheads.

## Goal

Sister-adapter expansion of the existing Phase 35 CRISPR design
trio (CHOPCHOP / CRISPOR / Cas-OFFinder). Round out the CRISPR
guide-RNA design surface with four canonical **base + prime
editing** design tools that the existing Cas9-cut-focused adapters
don't cover — the Komor / Liu lab base-editor design / outcome-
prediction / pegRNA-finder family that defines the modern non-
cleavage CRISPR editing surface. **BE-Designer** (Komor lab base-
editor guide design — the de-facto first stop for "I want to make
a C→T or A→G base change at this position, what guides should I
order"), **BE-Hive** (Liu lab base-editing outcome predictor —
predicts the per-position editing efficiency given a guide / chassis
/ window combination, the canonical readout for "given my designed
guide, what fraction of edits will land on the right base"),
**PrimeDesign** (Liu lab prime editing design — the canonical
pegRNA designer for the Anzalone / Liu prime-editing system, given
a target edit it designs the pegRNA + nicking guide pair), and
**pegFinder** (Komor lab alternative pegRNA finder — sister tool
to PrimeDesign with a different scoring model emphasizing pegRNA
secondary-structure stability + RT template length tradeoffs). All
four ride the established Python-script subprocess pattern (sister
to Phase 17 Biopython, Phase 19.5 Scanpy, Phase 33 pySBOL, Phase
35 CHOPCHOP / CRISPOR, Phase 41 pydna, Phase 43 DNA Chisel / iCodon)
— the user supplies a `.py` script that imports the upstream
package and reads `valenx_params.json` for the parsed knobs. Phase
35.5 sits numerically adjacent to Phase 35 and ships chronologically
right after Phase 44.5 RNA folding expansion — same chronological-
vs-numerical convention used for Phase 17.5 / 24 / 28 / 31 / 35 /
39 / 5.5 / 5.6 / 5.7 / 32.5 / 40 / 41 / 22.5 / 42 / 44.5.

## Capability inventory

### Live adapters (4)

- **BE-Designer** — Komor lab's base-editor guide design tool
  (MIT). BE-Designer is the de-facto first stop in modern base-
  editing workflows: given a target genome region with a desired
  C→T or A→G base change, it enumerates candidate gRNAs, scores
  each against the editing window of the requested base-editor
  variant (BE3, BE4max, ABE7.10, ABEmax, etc.), and emits a guide
  table with per-guide editing-window predictions, off-target
  predictions, and PAM-compatibility filtering. Python-script
  subprocess shape (sister to Phase 35 CHOPCHOP / CRISPOR, Phase
  17 Biopython, Phase 41 pydna, Phase 43 DNA Chisel / iCodon).
  The user supplies a `.py` script referenced from
  `[bio.be-designer].script` in `case.toml` that imports
  `bedesigner` and reads `valenx_params.json` for the parsed
  knobs. Schema knobs: `script` (path to user-supplied Python
  script; required, `.py` enforced), `python` (interpreter name;
  default `"python3"`), `input_fasta` (`Option<PathBuf>` —
  optional target sequence FASTA; `None` when the script
  generates the target inline or fetches it from a database),
  `output_basename` (filename stem; required, non-empty).
  `prepare()` enforces the `.py` extension, routes script +
  optional input_fasta through `confined_join` to stage them
  safely in the workdir, then writes a flat hand-rolled
  `valenx_params.json` containing `output_basename` always plus
  `input_fasta` (staged filename) only when set — the key is
  omitted entirely when `None` rather than emitted as `null`,
  matching the hand-rolled JSON convention the rest of the bio
  adapters use (Phase 19.6 Seurat / AnnData, Phase 27.5 ESM-IF,
  Phase 41 pydna, Phase 42 Mol* / NGL, Phase 43 DNA Chisel /
  iCodon, Phase 44.5 EternaFold). `collect()` walks the workdir
  for `<output_basename>*.csv` (`Tabular`, "BE-Designer guide
  table" — the canonical per-guide ranking with PAM, editing
  window position, and off-target preview), `<output_basename>*
  .fasta` (`Native`, "BE-Designer designed sequences" — the
  designed gRNA sequences ready for synthesis), and `*.log`
  (`Log`). Probe via Python on PATH then `<python> -c "import
  bedesigner"` — when the import fails the probe still returns
  `ok = true` with a targeted `"probe found python on PATH but
  could not import bedesigner — install via pip install
  bedesigner"` warning so users with Python ready but no
  `bedesigner` package see the install hint without failing the
  probe (sister to the Phase 19.5 scanpy / scvi / Phase 19.6
  AnnData / Phase 5.6 HOOMD-blue / Phase 5.7 MDTraj / Phase 41
  pydna / Phase 42 Mol* / NGL / Phase 43 DNA Chisel / Phase 44.5
  EternaFold probe convention). Version range `1.0.0..2.0.0`
  (BE-Designer 1.x is the modern stable line; upper bound 2.0
  reserves room for an eventual major bump). `bio.be_designer
  .design` ribbon capability.
- **BE-Hive** — Liu lab's base-editing outcome predictor (MIT).
  BE-Hive answers the canonical sister question to "what guides
  should I order" — given a designed guide / target / base-editor
  chassis, **what fraction of edits will land on the right base**,
  what fraction will have unintended bystander edits, and where
  in the editing window will the edits concentrate. The Liu lab's
  large-scale base-editing outcome dataset trains a CNN-style
  model that predicts per-position editing efficiency for every
  major base-editor chassis (the BE3 / BE4 / ABE / CBE families).
  Python-script subprocess shape (sister to BE-Designer, Phase
  17 Biopython, Phase 35 CHOPCHOP / CRISPOR). The user supplies
  a `.py` script referenced from `[bio.be-hive].script` in
  `case.toml` that imports `be_predict` (the canonical Python
  module name for the Liu lab's BE-Hive predictor — the canonical
  pip-installable package on PyPI is `be_predict`) and reads
  `valenx_params.json` for the parsed knobs. Schema knobs: `script`
  (`.py` enforced) / `python` (default `"python3"`) / `input_fasta`
  (`Option<PathBuf>`) / `output_basename`. `prepare()` follows the
  same Python-script + `valenx_params.json` shape as BE-Designer
  (key omitted when `None`). `collect()` walks the workdir for
  `<output_basename>*.csv` (`Tabular`, "BE-Hive efficiency
  predictions" — per-position editing efficiency + bystander
  predictions), `<output_basename>*.png` (`Native`, "BE-Hive
  plot" — per-position bar chart of predicted editing
  probability), and `*.log` (`Log`). Probe via Python on PATH
  then `<python> -c "import be_predict"` — same `ok = true` +
  warning fallback as BE-Designer for users whose `be_predict`
  install lives in a non-standard environment. Version range
  `1.0.0..2.0.0`. `bio.be_hive.predict` ribbon capability.
- **PrimeDesign** — Liu lab's prime-editing design tool (MIT).
  PrimeDesign is the canonical pegRNA designer for the Anzalone /
  Liu prime-editing system: given a desired edit at a target
  locus (point mutation, small insertion, small deletion, or
  combinations thereof up to a few dozen base pairs), it designs
  the pegRNA — the chimeric guide RNA with a 3′ extension
  encoding the reverse-transcriptase template + primer-binding
  site — plus the optional secondary nicking guide that boosts
  editing efficiency under PE3 / PE3b. Python-script subprocess
  shape (sister to BE-Designer / BE-Hive). The user supplies a
  `.py` script referenced from `[bio.primedesign].script` in
  `case.toml` that imports `primedesign` and reads
  `valenx_params.json` for the parsed knobs. Schema knobs:
  `script` (`.py` enforced) / `python` / `input_fasta` / 
  `output_basename` — same shape as BE-Designer / BE-Hive.
  `prepare()` follows the same Python-script + `valenx_params
  .json` shape (key omitted when `None`). `collect()` walks the
  workdir for `<output_basename>*.csv` (`Tabular`, "PrimeDesign
  pegRNA table" — per-pegRNA design with RT template length, PBS
  length, predicted edit, and nicking-guide pairing),
  `<output_basename>*.txt` (`Tabular`, "PrimeDesign report" —
  human-readable summary of the recommended pegRNA + nicking
  guide), and `*.log` (`Log`). Probe via Python on PATH then
  `<python> -c "import primedesign"` — same `ok = true` +
  warning fallback as BE-Designer / BE-Hive. Version range
  `1.0.0..2.0.0`. `bio.primedesign.design` ribbon capability.
- **pegFinder** — Komor lab's alternative pegRNA finder (MIT).
  pegFinder is sister to PrimeDesign with a different scoring
  model that emphasizes pegRNA secondary-structure stability +
  RT-template-length tradeoffs — the Komor lab's analysis showed
  pegRNA secondary structure substantially impacts prime-editing
  efficiency, and pegFinder explicitly scores candidate pegRNAs
  on predicted RT-template + PBS thermodynamic stability. Having
  both PrimeDesign and pegFinder lets users cross-check pegRNA
  recommendations between the Liu and Komor labs' design
  philosophies — the canonical use case in modern prime-editing
  workflows. Python-script subprocess shape (sister to PrimeDesign).
  The user supplies a `.py` script referenced from
  `[bio.pegfinder].script` in `case.toml` that imports `pegfinder`
  and reads `valenx_params.json` for the parsed knobs. Schema
  knobs: `script` (`.py` enforced) / `python` / `input_fasta` /
  `output_basename` — same shape as PrimeDesign. `prepare()`
  follows the same Python-script + `valenx_params.json` shape
  (key omitted when `None`). `collect()` walks the workdir for
  `<output_basename>*.csv` (`Tabular`, "pegFinder pegRNA
  candidates" — per-pegRNA design with secondary-structure
  scores), `<output_basename>*.txt` (`Tabular`, "pegFinder
  summary"), and `*.log` (`Log`). Probe via Python on PATH then
  `<python> -c "import pegfinder"` — same `ok = true` + warning
  fallback. Version range `1.0.0..2.0.0`. `bio.pegfinder.design`
  ribbon capability.

### Canonical types

**No new canonical types.** All four adapters consume user-
supplied inputs (`.py` scripts that import `bedesigner` /
`be_predict` / `primedesign` / `pegfinder`, plus optional starting
FASTA target-sequence files for each) and emit user-readable
artifacts (BE-Designer guide tables + designed gRNA sequences,
BE-Hive editing-efficiency tables + per-position bar charts,
PrimeDesign pegRNA tables + recommendation reports, pegFinder
pegRNA candidate tables + summary reports) that the unchanged
`Results.artifacts` collection model surfaces directly. The
existing `valenx_bio::format::fasta` reader already inspects
target FASTAs for sequence count + identifiers + alphabets. A
first-class CRISPR-editing canonical type — a typed guide /
pegRNA / editing-window / efficiency-prediction representation
spanning all seven CRISPR design + editing tools (Phase 35
CHOPCHOP / CRISPOR / Cas-OFFinder + Phase 35.5 BE-Designer /
BE-Hive / PrimeDesign / pegFinder) — defers to a future phase
along with cross-tool guide-comparison viewers and editing-
window heatmap viewers.

### Headless CLIs

**No new CLIs.** BE-Designer / BE-Hive / PrimeDesign / pegFinder
emit standard CSV / TSV / PNG / TXT formats inspectable in any
editor or through the user's downstream Python pipeline (`pandas`,
`numpy`). Target FASTAs are inspectable through the existing
Phase 17 `valenx-fasta` CLI. A canonical CRISPR-editing CLI —
guide-vs-pegRNA cross-comparison, editing-window heatmap,
bystander-edit visualisation — defers to a future phase along
with the canonical type.

## Domain expansion

Phase 35.5 is a **sister-adapter expansion of the Phase 35 CRISPR
design trio** (CHOPCHOP / CRISPOR / Cas-OFFinder) — the same
CRISPR-design surface broadened with four more established tools
that cover the **non-cleavage editing** corner Phase 35 doesn't
reach. Phase 35 covers Cas9 / Cas12a / Cas13 cut-style guide
design + off-target prediction; Phase 35.5 covers the modern
non-cleavage editing tools — base editors (BE-Designer + BE-Hive)
and prime editors (PrimeDesign + pegFinder). With Phase 35.5 the
CRISPR-design surface in Valenx covers the full design tradeoff
space — Cas9 cut design (CHOPCHOP), comprehensive off-target
prediction (CRISPOR), pure off-target searching (Cas-OFFinder),
base-editor guide design (BE-Designer), base-editing outcome
prediction (BE-Hive), prime-editing design via the Liu lab's
PrimeDesign, and prime-editing design via the Komor lab's
pegFinder. The Liu and Komor labs are the two dominant labs in
the modern non-cleavage CRISPR space — Phase 35.5 ships
representatives from both.

## What landed early

The implementation landed across 5
discrete implementation commits (4 adapters plus the registry +
init-template rollup) plus this docs pass — each landing one
adapter, the registry rollup, the init-template extension, or
the documentation pass. Every commit kept workspace clippy +
rustdoc clean.

## Acceptance checklist

- [x] `valenx-adapter-be-designer` adapter ships with case-input
      parser + lib tests + case-input tests covering parses-
      minimal / parses-with-input-fasta / rejects-non-py-script,
      plus the Python-script subprocess shape that enforces `.py`,
      routes script + optional input_fasta through `confined_join`,
      writes `valenx_params.json` with `output_basename` always
      plus `input_fasta` (staged filename) only when set — key
      omitted entirely when `None` rather than emitted as `null`,
      matching the hand-rolled JSON convention the rest of the
      bio adapters use, plus the Python on PATH + `import
      bedesigner` probe with `"could not import bedesigner"`
      warning when the import fails
- [x] `valenx-adapter-be-hive` adapter ships with case-input
      parser + lib tests + case-input tests, plus the Python-
      script subprocess shape (sister to BE-Designer with
      `import be_predict` probe and `<output_basename>*.png`
      collection for the per-position bar chart)
- [x] `valenx-adapter-primedesign` adapter ships with case-input
      parser + lib tests + case-input tests, plus the Python-
      script subprocess shape (sister to BE-Designer with `import
      primedesign` probe and `<output_basename>*.txt` recommended-
      pegRNA report collection)
- [x] `valenx-adapter-pegfinder` adapter ships with case-input
      parser + lib tests + case-input tests, plus the Python-
      script subprocess shape (sister to PrimeDesign with
      `import pegfinder` probe and `<output_basename>*.txt`
      summary report collection)
- [x] All 4 adapters wired into `valenx-app::init_registry` —
      live adapter count moves from 131 to **135** as part of the
      13-adapter / 4-phase rollup that takes the headline live-
      adapter total to **141**, rounding out the CRISPR editing-
      design surface that Phase 35 CHOPCHOP / CRISPOR / Cas-
      OFFinder opened
- [x] 4 base + prime editing templates in `valenx-init`
      (`be-designer`, `be-hive`, `primedesign`, `pegfinder`),
      all round-tripping through `valenx-validate` (cross-binary
      roundtrip now sweeps **136 templates** clean alongside the
      Phase 44.5 / 35.6 / 45 sister rollups)
- [x] STATUS.md / README.md / ROADMAP.md / CHANGELOG.md /
      QUICKSTART.md updated
- [ ] **Phase 35.7** — sister-adapter expansion of Phase 35 / 35.5:
      DeepCRISPR (deep-learning Cas9 efficiency predictor; defer),
      Azimuth 2 (Microsoft Research's classic on-target efficiency
      scorer; defer pending upstream activity), TIDE / TIDER
      (sequencing-trace edit deconvolution; defer to a future
      genomics-analysis phase), CRISPResso2 (Pinello lab's amplicon-
      sequencing edit-outcome quantifier; defer — already partially
      reachable through Phase 19 GATK / DeepVariant on amplicon
      data). Out of scope for this beachhead.

## Success metrics

| Metric                                                | Target          |
|-------------------------------------------------------|-----------------|
| New base/prime-editing adapter (template + tests)     | 1 day per       |
| Design-base-edit → predict-base-outcome → design-prime-edit → score-pegRNA loop across 4 tools | < tool baseline |

## Leads into

Phase 35.5 rounds out the CRISPR-design surface that the user's
bio / chemistry spec called out alongside the existing Phase 35
trio (CHOPCHOP / CRISPOR / Cas-OFFinder). Combined with the
existing fold-mfold → fold-EternaFold → fold-LinearFold →
optimize-codons → design-mRNA → predict-stability → render-Mol*
→ render-NGL → run-Galaxy-workflow → run-WDL → run-CWL → run-
Nextflow → run-Snakemake → design-plasmid → view-alignment →
process-image → segment-cells → classify-pixels → simulate-
pathway → expand-rules → grow-tissue → diffuse-particles →
trace-MCell-trajectories → simulate-MD → analyze-trajectory →
reweight-free-energy → fit-ENM → run-cpptraj-script → predict-
structure → fold-RNA → analyze-DNA-geometry → infer-tree-ML →
infer-tree-Bayesian → simulate-popgen → analyze-trees →
reconstruct-3D → design-protein → validate loop, the **design-
base-edit → predict-base-outcome → design-prime-edit → score-
pegRNA → fold-mfold → fold-EternaFold → fold-LinearFold →
optimize-codons → design-mRNA → predict-stability → render-Mol*
→ render-NGL → run-Galaxy-workflow → run-WDL → run-CWL → run-
Nextflow → run-Snakemake → design-plasmid → view-alignment →
process-image → segment-cells → classify-pixels → simulate-
pathway → expand-rules → grow-tissue → diffuse-particles →
trace-MCell-trajectories → simulate-MD → analyze-trajectory →
reweight-free-energy → fit-ENM → run-cpptraj-script → predict-
structure → fold-RNA → analyze-DNA-geometry → infer-tree-ML →
infer-tree-Bayesian → simulate-popgen → analyze-trees →
reconstruct-3D → design-protein → validate** loop now spans
seven CRISPR-design + editing tools (the Phase 35 CHOPCHOP /
CRISPOR / Cas-OFFinder cut-style trio plus Phase 35.5
BE-Designer / BE-Hive / PrimeDesign / pegFinder editing
quartet) feeding into the existing Phase 33 / 41 synthetic-
biology composition stacks (pySBOL / j5 / Cello + pydna), the
Phase 43 mRNA design stack (DNA Chisel / LinearDesign / iCodon),
and the Phase 44.5 RNA folding expansion (mfold / EternaFold /
LinearFold) — all in one Valenx shell with no glue code beyond
the existing case-toml / prepare / run / collect path.

The natural follow-up is **Phase 35.7** — the deferred CRISPR-
editing work called out above (DeepCRISPR / Azimuth 2 for the
ML-based efficiency-predictor surface, TIDE / TIDER / CRISPResso2
for the sequencing-trace edit-outcome quantification surface),
slotting in alongside the existing BE-Designer / BE-Hive /
PrimeDesign / pegFinder adapters with the same Python-script
subprocess shape.
