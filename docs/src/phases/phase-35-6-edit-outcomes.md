# Phase 35.6 — Edit-outcome prediction

**Status:** 🟢 Live — inDelphi + FORECasT + AlphaMissense + CRISPRitz
close the **CRISPR design → predict-outcome → off-target loop**
that Phase 35 CHOPCHOP / CRISPOR / Cas-OFFinder opened and Phase
35.5 BE-Designer / BE-Hive / PrimeDesign / pegFinder broadened
into base + prime editing, alongside the Phase 5.5 / 5.6 / 5.7 /
17 / 17.5 / 17.7 / 18 / 18.5 / 18.6 / 18.7 / 19 / 19.5 / 19.6 /
20 / 22 / 22.5 / 23 / 24 / 25 / 27 / 27.5 / 27.6 / 28 / 29 / 30 /
30.5 / 31 / 32 / 32.5 / 33 / 34 / 36 / 38 / 39 / 40 / 41 / 42 /
43 / 44.5 biology / biotech / chemistry beachheads.

## Goal

Sister-adapter expansion of the existing Phase 35 CRISPR design
(CHOPCHOP / CRISPOR / Cas-OFFinder) and Phase 35.5 base + prime
editing design (BE-Designer / BE-Hive / PrimeDesign / pegFinder).
Round out the CRISPR-design loop with four canonical **edit-outcome
predictors** that close the gap between "I designed an edit" and
"I know what the edit will actually do" — **inDelphi** (Liu lab's
Cas9-cut indel pattern predictor — given a guide / target /
chassis combination, predicts the per-indel-pattern frequency in
the resulting edited cell population, the canonical readout for
Cas9 cut outcomes), **FORECasT** (Sanger Institute's alternative
indel predictor — different model, different training corpus, used
to cross-check inDelphi predictions), **AlphaMissense** (DeepMind's
missense-effect predictor — given a designed protein-coding edit
producing a missense change, predicts the pathogenicity score
under the AlphaFold-grounded model; the canonical readout for
"will this missense edit cause disease"; CC-BY-NC-SA-4.0 / academic
non-commercial weights), and **CRISPRitz** (Pinello lab's off-
target genome-wide search — sister to Phase 35 Cas-OFFinder with
a different scoring model + variant-aware off-target searching
that accounts for population SNVs / SNPs in the off-target landscape).
All four ride the established Python-script subprocess pattern
(sister to Phase 17 Biopython, Phase 19.5 Scanpy, Phase 33 pySBOL,
Phase 35 CHOPCHOP / CRISPOR, Phase 35.5 BE-Designer / BE-Hive /
PrimeDesign / pegFinder, Phase 41 pydna, Phase 43 DNA Chisel /
iCodon). AlphaMissense's CC-BY-NC-SA-4.0 academic-only weights
follow the established AlphaFold 3 academic-license-warning
pattern. Phase 35.6 sits numerically adjacent to Phase 35.5 and
ships chronologically right after Phase 35.5 — same chronological-
vs-numerical convention used for Phase 17.5 / 24 / 28 / 31 / 35 /
39 / 5.5 / 5.6 / 5.7 / 32.5 / 40 / 41 / 22.5 / 42 / 44.5 / 35.5.

## Capability inventory

### Live adapters (4)

- **inDelphi** — Liu lab's Cas9-cut indel pattern predictor (MIT).
  inDelphi is the de-facto first stop in modern Cas9-cut-outcome
  workflows: given a designed gRNA + target site + cell-line
  chassis (mESC / U2OS / HCT116 / HEK293 etc.), it predicts the
  per-indel-pattern frequency distribution that will result from
  the Cas9 double-strand break and the cell's own end-joining
  repair machinery (NHEJ / MMEJ). The canonical readout is "given
  my designed guide, what indel patterns will I see in the
  edited cell population, and what fraction of edits will land
  in the desired frame-shift / clean-deletion / target-precise
  category". Python-script subprocess shape (sister to Phase 35
  CHOPCHOP / CRISPOR, Phase 35.5 BE-Designer / BE-Hive / Phase
  17 Biopython, Phase 41 pydna, Phase 43 DNA Chisel / iCodon).
  The user supplies a `.py` script referenced from
  `[bio.indelphi].script` in `case.toml` that imports `inDelphi`
  and reads `valenx_params.json` for the parsed knobs. Schema
  knobs: `script` (path to user-supplied Python script; required,
  `.py` enforced), `python` (interpreter name; default
  `"python3"`), `input_fasta` (`Option<PathBuf>` — optional
  target sequence FASTA; `None` when the script generates the
  target inline or fetches it from a database), `output_basename`
  (filename stem; required, non-empty). `prepare()` enforces the
  `.py` extension, routes script + optional input_fasta through
  `confined_join` to stage them safely in the workdir, then
  writes a flat hand-rolled `valenx_params.json` containing
  `output_basename` always plus `input_fasta` (staged filename)
  only when set — the key is omitted entirely when `None` rather
  than emitted as `null`, matching the hand-rolled JSON
  convention the rest of the bio adapters use (Phase 19.6 Seurat /
  AnnData, Phase 27.5 ESM-IF, Phase 41 pydna, Phase 42 Mol* /
  NGL, Phase 43 DNA Chisel / iCodon, Phase 44.5 EternaFold,
  Phase 35.5 BE-Designer / BE-Hive / PrimeDesign / pegFinder).
  `collect()` walks the workdir for `<output_basename>*.csv`
  (`Tabular`, "inDelphi indel predictions" — per-indel-pattern
  predicted frequency table), `<output_basename>*.png` (`Native`,
  "inDelphi plot" — predicted indel-pattern bar chart), and
  `*.log` (`Log`). Probe via Python on PATH then `<python> -c
  "import inDelphi"` — when the import fails the probe still
  returns `ok = true` with a targeted `"probe found python on
  PATH but could not import inDelphi — install via pip install
  inDelphi"` warning so users with Python ready but no
  `inDelphi` package see the install hint without failing the
  probe (sister to the Phase 19.5 scanpy / scvi / Phase 19.6
  AnnData / Phase 5.6 HOOMD-blue / Phase 5.7 MDTraj / Phase 41
  pydna / Phase 42 Mol* / NGL / Phase 43 DNA Chisel / Phase
  44.5 EternaFold / Phase 35.5 BE-Designer probe convention).
  Version range `1.0.0..2.0.0` (inDelphi 1.x is the modern
  stable line; upper bound 2.0 reserves room for an eventual
  major bump). `bio.indelphi.predict` ribbon capability.
- **FORECasT** — Sanger Institute's alternative Cas9-cut indel
  predictor (Apache-2.0). FORECasT is the alternative predictor
  to inDelphi — different model architecture, trained on a
  different corpus (Felicity Allen's Sanger SelfTarget library),
  validated against a different assay design — that lets users
  cross-check Cas9-cut outcome predictions across two independent
  groups' models. The canonical use case is "inDelphi predicted
  X% frame-shift; what does FORECasT predict for the same
  guide / target / chassis combination, and how do the two
  predictors agree on the dominant indel pattern". The Python
  module is published under the SelfTarget name (`import
  selftarget` rather than `import forecast`) — FORECasT is the
  predictor's published name in the original 2018 _Nature
  Biotechnology_ paper but the GitHub repository the canonical
  pip-install lives in is named `SelfTarget` after Allen's
  data-collection assay. Python-script subprocess shape (sister
  to inDelphi). The user supplies a `.py` script referenced from
  `[bio.forecast].script` in `case.toml` that imports `selftarget`
  and reads `valenx_params.json` for the parsed knobs. Schema
  knobs: `script` (`.py` enforced) / `python` / `input_fasta` /
  `output_basename` — same shape as inDelphi. `prepare()` follows
  the same Python-script + `valenx_params.json` shape (key
  omitted when `None`). `collect()` walks the workdir for
  `<output_basename>*.csv` (`Tabular`, "FORECasT indel
  predictions"), `<output_basename>*.txt` (`Tabular`, "FORECasT
  summary" — text summary of the dominant indel patterns + per-
  pattern frequency), and `*.log` (`Log`). Probe via Python on
  PATH then `<python> -c "import selftarget"` — same `ok = true`
  + warning fallback as inDelphi for users whose `selftarget`
  install lives in a non-standard environment. Version range
  `1.0.0..2.0.0`. `bio.forecast.predict` ribbon capability.
- **AlphaMissense** — DeepMind's missense-effect predictor (CC-BY-
  NC-SA-4.0 / academic non-commercial). AlphaMissense extends
  the AlphaFold structural-prediction lineage to score per-
  position missense-mutation pathogenicity: given a protein
  sequence + a missense change, it predicts the pathogenicity
  score on a 0-to-1 scale that calibrates against ClinVar
  pathogenic / benign labels. The canonical readout in modern
  CRISPR-editing workflows is "given my designed CRISPR edit
  causing a missense change at this protein position, what's
  the predicted pathogenicity — should I expect this edit to
  produce a benign protein variant, an ambiguous variant, or a
  loss-of-function / pathogenic variant". Python-script
  subprocess shape (sister to inDelphi / FORECasT). The user
  supplies a `.py` script referenced from `[bio.alphamissense]
  .script` in `case.toml` that imports `alphamissense` and reads
  `valenx_params.json` for the parsed knobs. Schema knobs:
  `script` (`.py` enforced) / `python` / `input_fasta` /
  `output_basename` — same shape as inDelphi. `prepare()` follows
  the same Python-script + `valenx_params.json` shape (key
  omitted when `None`). `collect()` walks the workdir for
  `<output_basename>*.csv` (`Tabular`, "AlphaMissense
  pathogenicity scores"), `<output_basename>*.tsv` (`Tabular`,
  "AlphaMissense pathogenicity scores" — alternate extension
  AlphaMissense's exporters write under), `<output_basename>*
  .png` (`Native`, "AlphaMissense plot" — per-position
  pathogenicity bar chart), and `*.log` (`Log`). Probe via
  Python on PATH then `<python> -c "import alphamissense"` —
  when the import fails the probe still returns `ok = true`
  with a warning. **License callout:** AlphaMissense's model
  weights are released under the Creative Commons Attribution-
  NonCommercial-ShareAlike 4.0 International license, which
  restricts commercial use of the predictor's outputs. The
  probe surfaces an `"academic"` / `"non-commercial"`-keyworded
  license-awareness warning in `ProbeReport.warnings` whenever
  Python is detected — sister to the AlphaFold 3 academic
  warning pattern (Phase 17.5 AlphaFold 3 / Phase 28 ViennaRNA /
  NUPACK / Phase 23 VMD / Phase 5.6 NAMD / Phase 44.5 mfold).
  The license warning ALWAYS surfaces when probe finds Python
  on PATH, **regardless of whether `import alphamissense`
  succeeds** — sister to AlphaFold 3's mandatory probe warning.
  Version range `1.0.0..2.0.0` (AlphaMissense 1.x is the modern
  stable line; upper bound 2.0 reserves room for an eventual
  major bump). `bio.alphamissense.predict` ribbon capability.
- **CRISPRitz** — Pinello lab's variant-aware off-target genome-
  wide search (MIT). CRISPRitz is sister to Phase 35 Cas-OFFinder
  with a different scoring model and the distinguishing
  property of **variant-aware off-target searching** — given a
  reference genome plus a population VCF (1000 Genomes, gnomAD,
  etc.), CRISPRitz searches for off-target sites that exist
  only in specific haplotypes / specific population sub-groups
  rather than just the reference assembly. This is the canonical
  readout for "is my CRISPR guide safe across the human
  population, or are there off-targets that would only manifest
  in certain genetic backgrounds". Python-script subprocess
  shape (sister to inDelphi / FORECasT / AlphaMissense). The
  user supplies a `.py` script referenced from
  `[bio.crispritz].script` in `case.toml` that imports
  `crispritz` and reads `valenx_params.json` for the parsed
  knobs. Schema knobs: `script` (`.py` enforced) / `python` /
  `input_fasta` / `output_basename` — same shape as inDelphi.
  `prepare()` follows the same Python-script + `valenx_params
  .json` shape (key omitted when `None`). `collect()` walks the
  workdir for `<output_basename>*.txt` (`Tabular`, "CRISPRitz
  off-target table"), `<output_basename>*.bed` (`Tabular`,
  "CRISPRitz off-target BED" — UCSC BED-format off-target hits
  ready for genome-browser visualisation), and `*.log` (`Log`).
  Probe via Python on PATH then `<python> -c "import crispritz"`
  — same `ok = true` + warning fallback as inDelphi / FORECasT
  for users whose `crispritz` install lives in a non-standard
  environment. Version range `2.6.0..3.0.0` (CRISPRitz 2.6 is
  the modern stable line shipping the contemporary variant-
  aware mode; upper bound 3.0 reserves room for an eventual
  major bump). `bio.crispritz.search` ribbon capability.

### Canonical types

**No new canonical types.** All four adapters consume user-
supplied inputs (`.py` scripts that import `inDelphi` /
`selftarget` / `alphamissense` / `crispritz`, plus optional
target / protein FASTA inputs for each) and emit user-readable
artifacts (inDelphi indel-pattern frequency tables + per-pattern
plots, FORECasT indel-pattern tables + summary text, AlphaMissense
per-position pathogenicity scores in CSV / TSV + bar-chart plots,
CRISPRitz tabular off-target tables + UCSC BED-format off-target
hits) that the unchanged `Results.artifacts` collection model
surfaces directly. The existing `valenx_bio::format::fasta` reader
already inspects target / protein FASTAs for sequence count +
identifiers + alphabets. A first-class CRISPR edit-outcome
canonical type — a typed indel-pattern / pathogenicity / off-
target representation spanning all eleven CRISPR design + editing
+ outcome tools (Phase 35 CHOPCHOP / CRISPOR / Cas-OFFinder + Phase
35.5 BE-Designer / BE-Hive / PrimeDesign / pegFinder + Phase
35.6 inDelphi / FORECasT / AlphaMissense / CRISPRitz) — defers to
a future phase along with cross-tool outcome-comparison viewers
and population-aware off-target heatmap viewers.

### Headless CLIs

**No new CLIs.** inDelphi / FORECasT / AlphaMissense / CRISPRitz
emit standard CSV / TSV / TXT / BED / PNG formats inspectable in
any editor or through the user's downstream Python pipeline
(`pandas`, `numpy`, `pybedtools`). Target FASTAs are inspectable
through the existing Phase 17 `valenx-fasta` CLI; off-target BED
files are consumable through any standard genomics tool. A
canonical CRISPR-edit-outcome CLI — predictor-vs-predictor indel-
pattern diffing, missense-pathogenicity comparison against
ClinVar, variant-aware off-target prevalence inspection — defers
to a future phase along with the canonical type.

## Domain expansion

Phase 35.6 is a **sister-adapter expansion of the Phase 35 CRISPR
design trio + Phase 35.5 base + prime editing quartet** — the
same CRISPR-design surface broadened with four more established
tools that close the **prediction loop** between design and
outcome. Phase 35 covers Cas9 / Cas12a / Cas13 cut-style guide
design + off-target prediction; Phase 35.5 covers base + prime
editing design via BE-Designer / BE-Hive / PrimeDesign / pegFinder;
Phase 35.6 closes the loop with the modern outcome-prediction
tools — Cas9 indel pattern predictors (inDelphi + FORECasT for
cross-checking between two independent labs), missense pathogenicity
prediction (AlphaMissense for "will this missense edit cause
disease"), and variant-aware off-target searching (CRISPRitz for
population-genome off-target safety). With Phase 35.6 the
CRISPR-design surface in Valenx covers the full design → predict-
outcome → off-target loop with eleven adapters spanning every
canonical question in modern CRISPR workflows. AlphaMissense
inherits the CC-BY-NC-SA-4.0 academic-non-commercial-weights
license caveat from DeepMind, surfaced via the established
academic-license probe-warning pattern (sister to AlphaFold 3 /
ViennaRNA / NUPACK / VMD / NAMD / mfold).

## What landed early

The implementation landed across 5
discrete implementation commits (4 adapters plus the registry +
init-template rollup) plus this docs pass — each landing one
adapter, the registry rollup, the init-template extension, or
the documentation pass. Every commit kept workspace clippy +
rustdoc clean.

## Acceptance checklist

- [x] `valenx-adapter-indelphi` adapter ships with case-input
      parser + lib tests + case-input tests covering parses-
      minimal / parses-with-input-fasta / rejects-non-py-script,
      plus the Python-script subprocess shape that enforces `.py`,
      routes script + optional input_fasta through `confined_join`,
      writes `valenx_params.json` with `output_basename` always
      plus `input_fasta` (staged filename) only when set — key
      omitted entirely when `None` rather than emitted as `null`,
      matching the hand-rolled JSON convention the rest of the
      bio adapters use, plus the Python on PATH + `import
      inDelphi` probe with `"could not import inDelphi"` warning
      when the import fails
- [x] `valenx-adapter-forecast` adapter ships with case-input
      parser + lib tests + case-input tests, plus the Python-
      script subprocess shape (sister to inDelphi with `import
      selftarget` probe — note the probe targets the `selftarget`
      module name even though the adapter is `forecast`, since
      the upstream PyPI distribution lives under the SelfTarget
      data-collection-assay name — and `<output_basename>*.txt`
      summary collection)
- [x] `valenx-adapter-alphamissense` adapter ships with case-
      input parser + lib tests + case-input tests, plus the
      Python-script subprocess shape (sister to inDelphi with
      `import alphamissense` probe and `<output_basename>*.tsv`
      / `*.csv` pathogenicity score collection + `<output_basename>*
      .png` per-position bar chart collection), and the **mandatory
      academic / non-commercial license-awareness warning** that
      surfaces in `ProbeReport.warnings` whenever Python is
      detected — sister to the Phase 17.5 AlphaFold 3 mandatory
      probe-warning pattern, regardless of whether `import
      alphamissense` succeeds
- [x] `valenx-adapter-crispritz` adapter ships with case-input
      parser + lib tests + case-input tests, plus the Python-
      script subprocess shape (sister to inDelphi with `import
      crispritz` probe and `<output_basename>*.bed` UCSC BED-
      format off-target hits collection alongside the
      `<output_basename>*.txt` off-target table)
- [x] All 4 adapters wired into `valenx-app::init_registry` —
      live adapter count moves from 135 to **139** as part of the
      13-adapter / 4-phase rollup that takes the headline live-
      adapter total to **141**, closing the design → predict-
      outcome → off-target loop that Phase 35 CHOPCHOP / CRISPOR
      / Cas-OFFinder + Phase 35.5 BE-Designer / BE-Hive /
      PrimeDesign / pegFinder opened
- [x] 4 edit-outcome templates in `valenx-init` (`indelphi`,
      `forecast`, `alphamissense`, `crispritz`), all round-
      tripping through `valenx-validate` (cross-binary roundtrip
      now sweeps **136 templates** clean alongside the Phase 44.5
      / 35.5 / 45 sister rollups)
- [x] STATUS.md / README.md / ROADMAP.md / CHANGELOG.md /
      QUICKSTART.md updated
- [ ] **Phase 35.7** — sister-adapter expansion of Phase 35 / 35.5
      / 35.6: SpCas9-HF / eSpCas9 high-fidelity off-target
      predictors (defer pending upstream activity), DeepBE
      (deep-learning base-editing outcome predictor sister to
      BE-Hive; defer pending licensing review of the model
      checkpoint), CRISPRme (population-aware variant-aware off-
      target sister to CRISPRitz; defer — CRISPRitz already
      covers the variant-aware shape), Plant-edit (plant-genome-
      specific off-target searching; defer pending domain
      coverage decision). Out of scope for this beachhead.

## Success metrics

| Metric                                                | Target          |
|-------------------------------------------------------|-----------------|
| New edit-outcome adapter (template + tests)           | 1 day per       |
| Predict-cas9-indels → predict-base-bystander → predict-pathogenicity → search-off-targets loop across 4 tools | < tool baseline |

## Leads into

Phase 35.6 closes the CRISPR design → predict-outcome → off-target
loop that the user's bio / chemistry spec called out alongside
the existing Phase 35 trio (CHOPCHOP / CRISPOR / Cas-OFFinder)
and the Phase 35.5 base + prime editing quartet (BE-Designer /
BE-Hive / PrimeDesign / pegFinder). Combined with the existing
design-base-edit → predict-base-outcome → design-prime-edit →
score-pegRNA → fold-mfold → fold-EternaFold → fold-LinearFold →
optimize-codons → design-mRNA → predict-stability → render-Mol*
→ render-NGL → run-Galaxy-workflow → run-WDL → run-CWL → run-
Nextflow → run-Snakemake → design-plasmid → view-alignment →
process-image → segment-cells → classify-pixels → simulate-
pathway → expand-rules → grow-tissue → diffuse-particles →
trace-MCell-trajectories → simulate-MD → analyze-trajectory →
reweight-free-energy → fit-ENM → run-cpptraj-script → predict-
structure → fold-RNA → analyze-DNA-geometry → infer-tree-ML →
infer-tree-Bayesian → simulate-popgen → analyze-trees →
reconstruct-3D → design-protein → validate loop, the **predict-
cas9-indels → cross-check-indels → predict-missense-
pathogenicity → search-population-off-targets → design-base-
edit → predict-base-outcome → design-prime-edit → score-pegRNA
→ fold-mfold → fold-EternaFold → fold-LinearFold → optimize-
codons → design-mRNA → predict-stability → render-Mol* →
render-NGL → run-Galaxy-workflow → run-WDL → run-CWL → run-
Nextflow → run-Snakemake → design-plasmid → view-alignment →
process-image → segment-cells → classify-pixels → simulate-
pathway → expand-rules → grow-tissue → diffuse-particles →
trace-MCell-trajectories → simulate-MD → analyze-trajectory →
reweight-free-energy → fit-ENM → run-cpptraj-script → predict-
structure → fold-RNA → analyze-DNA-geometry → infer-tree-ML →
infer-tree-Bayesian → simulate-popgen → analyze-trees →
reconstruct-3D → design-protein → validate** loop now spans
eleven CRISPR-design + editing + outcome tools (the Phase 35
CHOPCHOP / CRISPOR / Cas-OFFinder cut-style trio plus Phase 35.5
BE-Designer / BE-Hive / PrimeDesign / pegFinder editing quartet
plus Phase 35.6 inDelphi / FORECasT / AlphaMissense / CRISPRitz
outcome-prediction quartet) feeding into the existing Phase 33 /
41 synthetic-biology composition stacks (pySBOL / j5 / Cello +
pydna), the Phase 43 mRNA design stack (DNA Chisel /
LinearDesign / iCodon), and the Phase 44.5 RNA folding expansion
(mfold / EternaFold / LinearFold) — closing the **canonical
CRISPR-editing safety + outcome-prediction pipeline**: target
guide → CHOPCHOP / CRISPOR ranking → Cas-OFFinder / CRISPRitz
off-target → inDelphi / FORECasT indel prediction → BE-Hive
base-editing outcome → AlphaMissense missense pathogenicity → all
in one Valenx shell with no glue code beyond the existing case-
toml / prepare / run / collect path.

The natural follow-up is **Phase 35.7** — the deferred CRISPR-
outcome work called out above (SpCas9-HF / eSpCas9 high-fidelity
off-target predictors, DeepBE deep-learning base-editing outcome
sister to BE-Hive, CRISPRme variant-aware off-target sister to
CRISPRitz, Plant-edit plant-genome-specific off-target searching),
slotting in alongside the existing inDelphi / FORECasT /
AlphaMissense / CRISPRitz adapters with the same Python-script
subprocess shape.
