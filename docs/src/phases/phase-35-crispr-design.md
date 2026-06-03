# Phase 35 — CRISPR design

**Status:** 🟢 Live — CHOPCHOP + CRISPOR + Cas-OFFinder open the
**first CRISPR guide-RNA design domain** in Valenx alongside the
Phase 17 / 17.5 / 18 / 25 / 27 / 27.5 / 28 / 30 / 31 / 32 / 34 / 36
biology + structure-prediction + protein-design + RNA-structure +
phylogenetics + read-simulator + systems-biology + docking + cryo-EM
beachheads and the Phase 24 cheminformatics expansion.

## Goal

Open the CRISPR guide-RNA design + off-target searching domain in
Valenx with three established open-source tools that span the
CRISPR-design tradeoff space — popular ranked guide design with
off-target scoring (CHOPCHOP), comprehensive guide design plus
rigorous off-target prediction across many enzymes (CRISPOR), and
pure off-target searching used as a primitive by most other
CRISPR-design web services and pipelines (Cas-OFFinder). CHOPCHOP +
CRISPOR follow the established Phase 17 Biopython Python-script
subprocess shape: the user supplies a Python script that imports
the upstream package and reads `valenx_params.json` for the
parsed knobs, the adapter stages script + target FASTA, and
`run()` invokes `python <script>`. Cas-OFFinder follows the
established Phase 18 BWA single-binary CLI pattern: input file +
backend selector + output path in, ranked off-target hit table out.
Phase 35 sits numerically before Phase 36 but ships chronologically
right after Phase 31 read simulators — same chronological-vs-
numerical convention used for Phase 17.5 / 24 / 28 / 31.

## Capability inventory

### Live adapters (3)

- **CHOPCHOP** — University of Bergen's web-and-script CRISPR
  guide-RNA design tool (MIT). The de-facto first stop for "I have a
  gene, what should I cut" in academic CRISPR workflows: scores
  candidate gRNAs against a target sequence under a configurable
  nuclease (Cas9, Cas12a, Cas13) or TALEN design pass, ranks by
  efficiency / specificity / off-target risk, and emits both a
  guide-ranking TSV and a guide-location BED. Python-script
  subprocess shape (sister to Phase 17 Biopython): the user supplies
  a Python script referenced from `[bio.chopchop].script` in
  `case.toml` that imports `chopchop` (or invokes `chopchop.py`) and
  reads `valenx_params.json` for the parsed knobs. Schema knobs:
  `script` (path to user-supplied Python script; required), `python`
  (interpreter name; default `"python3"`), `target` (target sequence
  FASTA; required), `genome` (CHOPCHOP-installed genome name —
  `"hg38"` / `"mm10"` / etc.; required), `cas_variant` (one of
  `"Cas9"` / `"Cas12a"` / `"Cas13"` / `"TALEN"`; required), `pam`
  (PAM sequence — `"NGG"` for Cas9, `"TTTV"` for Cas12a, etc.;
  required), `output_basename` (filename stem; required, non-empty).
  `prepare()` stages the script + target FASTA into the workdir,
  writes a flat `valenx_params.json` containing `target` (staged
  filename), `genome`, `cas_variant`, `pam`, `output_basename`, and
  composes `python <script_filename>` as the native command.
  `collect()` walks the workdir for `<output_basename>*.tsv`
  (`Tabular`, "CHOPCHOP guide rankings") and `<output_basename>*.bed`
  (`Tabular`, "CHOPCHOP guide locations"). Probe via Python on PATH
  with an `import chopchop` check — when the import fails the probe
  still returns `ok = true` with a warning so users with CHOPCHOP
  installed under a non-standard module name (`crispr_chopchop`,
  `chopchop_v3`) or invoked via the user-supplied script aren't
  blocked. Version range `3.0.0..4.0.0` (the modern web / script
  split landed in 3.0; upper bound 4.0 reserves room for an eventual
  major bump). `bio.chopchop.design` ribbon capability.
- **CRISPOR** — Maximilian Haeussler's CRISPR guide-RNA design +
  off-target prediction tool (GPL-3.0). CRISPOR's distinguishing
  feature is the rigorous off-target pass: scores candidate guides
  against a reference genome assembly with the CFD scoring model
  and reports an MIT-style specificity score per guide. It powers
  the public crispor.org service and is also distributed as a
  standalone Python script for batch / pipeline use, supporting
  many more enzymes / PAMs than CHOPCHOP. Python-script subprocess
  shape (sister to CHOPCHOP): the user supplies a Python script
  referenced from `[bio.crispor].script` in `case.toml`. Schema
  knobs: `script` (path to user-supplied Python script; required),
  `python` (interpreter name; default `"python3"`), `target` (target
  sequence FASTA; required), `genome` (CRISPOR-supported genome
  name; required), `pam` (PAM motif — `"NGG"` / `"NG"` / `"TTTV"` /
  etc.; required), `batch_id` (optional string; CRISPOR caches
  partial results by batch so passing the same `batch_id` resumes a
  previously-interrupted run), `output_basename` (filename stem;
  required, non-empty). `prepare()` stages the script + target
  FASTA, writes a flat `valenx_params.json` containing `target`
  (staged filename), `genome`, `pam`, `batch_id` (JSON string or
  literal `null`), `output_basename`, and composes `python
  <script_filename>` as the native command. `collect()` walks the
  workdir for `<output_basename>*.tsv` (`Tabular`, "CRISPOR guide
  rankings") and `<output_basename>*.txt` (`Log`). Probe via Python
  on PATH with an `import crispor` check (same `ok = true` + warning
  fallback as CHOPCHOP for users with non-standard installs). Version
  range `5.0.0..6.0.0` (the modern Python 3 / batch-mode rewrite
  landed in 5.0; upper bound 6.0 reserves room for an eventual
  major bump). `bio.crispor.design` ribbon capability.
- **Cas-OFFinder** — Bae / Park / Kim group's CRISPR off-target
  searching tool from Hanyang / Seoul National University
  (BSD-3-Clause). Cas-OFFinder is a fast, OpenCL-accelerated scanner:
  given a list of guide sequences + PAM patterns + mismatch budget
  in a plain-text input file, it walks a reference genome and reports
  every position whose sequence matches one of the guides within the
  configured Hamming distance. It's the workhorse off-target scanner
  sitting under most CRISPR design web services (CRISPOR,
  CRISPRdirect, …) and pipelines. Single-binary subprocess shape
  (sister to Phase 18 BWA): the CLI is fixed-shape `cas-offinder
  <input> {C|G|A} <output> [extras...]`. `<input>` is a 3+-line
  text file with the reference genome path, the PAM pattern, and one
  guide-sequence row per query. The middle positional argument
  selects the OpenCL device class — `C` (CPU), `G` (GPU), or `A`
  (auto-pick fastest at runtime). Schema knobs: `input`
  (Cas-OFFinder input file; required), `output` (output text file;
  required), `backend` (one of `"C"` / `"G"` / `"A"`; required),
  `extra_args`. `prepare()` resolves both paths against the case
  directory (when relative) and composes the invocation positionally
  — no `-i` / `-o` flags, the order is fixed. `collect()` reports
  the configured `output` file as a single `Tabular` artifact
  (`"Cas-OFFinder off-target hits"`). Probe via
  `find_on_path(&["cas-offinder"])`. Version range `2.4.0..3.0.0`
  (the modern OpenCL device-selection CLI stabilised at 2.4; upper
  bound 3.0 reserves room for an eventual major bump). The init
  alias `cas-off` resolves to the same template as the canonical
  `cas-offinder`. `bio.cas-offinder.search` ribbon capability.

### Canonical types

**No new canonical types.** All three adapters consume user-supplied
inputs (CHOPCHOP + CRISPOR target FASTAs and design Python scripts,
Cas-OFFinder fixed-shape input files specifying the genome path +
PAM + per-guide query lines) and emit user-readable artifacts
(CHOPCHOP guide-ranking TSV + guide-location BED, CRISPOR
guide-ranking TSV + log TXT, Cas-OFFinder ranked off-target hit
TSV) that the unchanged `Results.artifacts` collection model
surfaces directly. The existing `valenx_bio::format::fasta` reader
already inspects target FASTAs for sequence count + identifiers +
alphabets. A first-class CRISPR-design canonical type — a generic
guide / off-target / scoring type spanning all three back-ends —
defers to a future phase along with guide-ranking visualizers and
off-target heatmap viewers.

### Headless CLIs

**No new CLIs.** CHOPCHOP's guide-ranking TSVs + BED location
files, CRISPOR's TSV / log files, and Cas-OFFinder's tabular
off-target hit text file are all standard tabular formats
inspectable in any editor or through the user's downstream Python
pipeline (`pandas`, `numpy`). Genome FASTAs and target sequences
are inspectable through the existing Phase 17 `valenx-fasta` CLI.
A canonical CRISPR-design CLI defers to a future phase along with
guide-comparison and on-target-vs-off-target diffing integrations.

## Domain milestone

Phase 35 is the **first CRISPR guide-RNA design domain** to land in
Valenx. The biology adapter family started with Phase 17 (foundation
— sequence / structure / trajectory canonical types + classical MD
+ cheminformatics scripts) and expanded through Phase 17.5 / 18 /
18.5 / 18.6 / 19 / 19.5 / 20 / 22 / 23 / 24 / 25 / 27 / 27.5 / 27.6
/ 28 / 30 / 31 / 32 / 34 / 36 to cover sequence prediction,
alignment, RNA-seq, variant calling, single-cell, transcript
quantification, workflow orchestration, molecular viewers,
cheminformatics, quantum chemistry, protein design, EvolutionaryScale
models, RNA structure, phylogenetics, sequencing read simulation,
systems biology, small-molecule docking, and cryo-EM reconstruction
— but until Phase 35 the CRISPR-design surface (guide-RNA scoring,
off-target prediction, off-target searching) was absent. Phase 35
closes that gap with three established open-source tools spanning
the CRISPR-design tradeoff space — CHOPCHOP at the popular
ranked-design end, CRISPOR for comprehensive design + rigorous
off-target prediction across many enzymes, and Cas-OFFinder as the
canonical off-target-searching primitive.

## What landed early

The implementation landed across 5
discrete implementation commits (3 adapters, the registry rollup,
the init-template extension) plus this docs pass — each landing one
adapter, the registry rollup, the init-template extension, or the
documentation pass. Every commit kept workspace clippy + rustdoc
clean.

## Acceptance checklist

- [x] `valenx-adapter-chopchop` adapter ships with case-input parser
      + 4 lib tests + 4 case-input tests covering parses-minimal /
      rejects-empty-script / rejects-bad-cas-variant / rejects-
      empty-pam / rejects-empty-output-basename, plus the
      Python-script subprocess shape that stages script + target
      FASTA + writes `valenx_params.json` with the parsed knobs
- [x] `valenx-adapter-crispor` adapter ships with case-input parser
      + 4 lib tests + 4 case-input tests covering parses-minimal /
      parses-with-batch-id / rejects-empty-target / rejects-empty-
      genome / rejects-empty-output-basename, plus the optional
      `batch_id` shape that emits the JSON literal `null` when
      omitted (so user scripts can always do `params["batch_id"]`
      without an `in` check)
- [x] `valenx-adapter-cas-offinder` adapter ships with case-input
      parser + 4 lib tests + 4 case-input tests covering parses-
      minimal / rejects-bad-backend / rejects-empty-input / rejects-
      empty-output, plus the fixed-shape positional CLI that emits
      `cas-offinder <input> <backend> <output> [extras...]` with
      no flag prefixes
- [x] All 3 adapters wired into `valenx-app::init_registry` —
      live adapter count moves from 81 to **84**, opening the first
      CRISPR guide-RNA design domain to ship in Valenx
- [x] 3 CRISPR-design templates in `valenx-init` (`chopchop`,
      `crispor`, `cas-offinder` with alias `cas-off`), all round-
      tripping through `valenx-validate` (cross-binary roundtrip
      now sweeps **80 templates** clean)
- [x] STATUS.md / README.md / ROADMAP.md / CHANGELOG.md /
      QUICKSTART.md updated
- [ ] **Phase 35.5** — CRISPRitz (in-silico off-target search with
      variant-aware scoring; sister to Cas-OFFinder; defer to
      sister-adapter expansion phase), FlashFry (high-throughput
      guide design + scoring; defer), E-CRISP (Boutros lab guide
      design with conservation scoring; defer), CRISPRdirect
      (Naito lab web-service guide selector; defer), Guidescan
      (off-target enumeration via specificity scoring; defer to
      35.5). Out of scope for this beachhead. CRISPResso2 (CRISPR
      editing analysis from sequencing data — scope is
      **post-editing analysis** rather than guide design; different
      shape; defer to a future phase).

## Success metrics

| Metric                                                | Target          |
|-------------------------------------------------------|-----------------|
| New CRISPR-design adapter (template + tests)          | 1 day per       |
| Design + off-target-predict + off-target-search loop across 3 tools | < tool baseline |

## Leads into

Phase 35 opens the CRISPR guide-RNA design domain that the user's
bio / chemistry spec called out alongside the Phase 17 / 17.5
biology + structure-prediction stack and the Phase 27 / 27.5 /
27.6 protein-design beachhead. Combined with the existing
simulate-reads → align → quantify → call-variants → predict-
structure → fold-RNA → infer-tree → simulate-pathway → reconstruct-
3D → validate loop, the **design-guides → predict-off-targets →
search-off-targets → simulate-reads → align → call-variants →
predict-structure → fold-RNA → infer-tree → simulate-pathway →
reconstruct-3D → validate** loop now spans three CRISPR-design
tools (CHOPCHOP, CRISPOR, Cas-OFFinder) feeding into the existing
Phase 31 read simulators (ART, wgsim, Badread), the eleven Phase
18 / 18.5 / 18.6 alignment tools (BWA, Bowtie2, HISAT2, STAR,
minimap2, MAFFT, MUSCLE, HMMER, samtools, MMseqs2, DIAMOND), the
two Phase 20 transcript quantifiers (Salmon, Kallisto), the three
Phase 19 variant callers (bcftools, GATK, DeepVariant), the Phase
17 / 17.5 prediction stack (ESMFold, OpenFold, AlphaFold 2/3,
ColabFold), the Phase 28 RNA-structure tools (ViennaRNA,
RNAstructure, NUPACK), the Phase 30 phylogenetic-tree builders
(IQ-TREE, RAxML-NG, FastTree), the Phase 32 systems-biology
surface (COPASI, BioNetGen, PhysiCell), and the Phase 36 cryo-EM
reconstruction tools (RELION, EMAN2, CTFFIND) — all in one Valenx
shell with no glue code beyond the existing case-toml / prepare /
run / collect path.

The natural follow-up is **Phase 35.5** — the deferred CRISPR-design
work called out above (CRISPRitz as an in-silico off-target search
with variant-aware scoring sister to Cas-OFFinder, FlashFry for
high-throughput guide design + scoring, E-CRISP for the Boutros
lab guide design with conservation scoring, CRISPRdirect for the
Naito lab web-service guide selector, Guidescan for off-target
enumeration via specificity scoring), slotting in alongside the
existing CRISPR-design adapters with the same Python-script
subprocess shape (CHOPCHOP / CRISPOR sister tools) or the
single-binary subprocess shape (Cas-OFFinder sister tools).
CRISPResso2 (post-editing analysis from sequencing data) sits in
a separate phase — the data shape is different enough (mapped
read alignment + indel calling rather than guide-RNA design) to
warrant a sister phase rather than 35.5 expansion.
