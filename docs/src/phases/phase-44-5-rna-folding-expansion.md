# Phase 44.5 — RNA folding expansion

**Status:** 🟢 Live — mfold + EternaFold + LinearFold round out the
**RNA secondary-structure folding surface** that Phase 28 ViennaRNA /
RNAstructure / NUPACK opened, alongside the Phase 5.5 / 5.6 / 5.7 /
17 / 17.5 / 17.7 / 18 / 18.5 / 18.6 / 18.7 / 19 / 19.5 / 19.6 / 20 /
22 / 22.5 / 23 / 24 / 25 / 27 / 27.5 / 27.6 / 29 / 30 / 30.5 / 31 /
32 / 32.5 / 33 / 34 / 35 / 36 / 38 / 39 / 40 / 41 / 42 / 43 biology
/ biotech / chemistry beachheads.

## Goal

Sister-adapter expansion of the existing Phase 28 RNA structure
trio (ViennaRNA / RNAstructure / NUPACK). Round out the RNA
secondary-structure folding surface with three more canonical RNA
folders that span the modern tradeoff space — Michael Zuker's
classic dynamic-programming RNA folder that defined the field
(mfold/UNAFold, the academic-license workhorse the entire RNA-
structure literature was built on), the Eterna game's ML-aware
folder reachable via the `arnie` Python wrapper (EternaFold, MIT —
trained on a half-decade of crowd-sourced Eterna gameplay puzzles
plus thermodynamic + machine-learning corpora), and Baidu / Oregon
State's beam-search linear-time folder (LinearFold, the Apache-2.0
sister to LinearDesign that scales to long sequences — viral
genomes, CRISPR-prep mRNAs, full transcripts — at constant per-
nucleotide cost). mfold + LinearFold follow the established Phase
18 BWA single-binary CLI pattern: sequence in, structure +
connect-table / dot-bracket out. EternaFold follows the established
Phase 17 Biopython Python-script subprocess shape: the user
supplies a Python script that imports `arnie` (the canonical EternaFold
front-end) and reads `valenx_params.json` for the parsed knobs.
Phase 44.5 sits numerically adjacent to Phase 44 (the registry-
palette wire-in landed in commit `a29cdeb`) and ships chronologically
right after Phase 43 mRNA design — same chronological-vs-numerical
convention used for Phase 17.5 / 24 / 28 / 31 / 35 / 39 / 5.5 /
5.6 / 5.7 / 32.5 / 40 / 41 / 22.5 / 42.

## Capability inventory

### Live adapters (3)

- **mfold/UNAFold** — Michael Zuker's classic RNA / DNA secondary-
  structure folder (academic-use license). mfold is the original
  dynamic-programming Zuker / Stiegler RNA folder that defined the
  field: minimum-free-energy folding plus a configurable suboptimal-
  structure ensemble within an energy budget of the MFE optimum,
  using the canonical Turner / Mathews thermodynamic parameters
  loaded by temperature. The modern `mfold` / `UNAFold.pl` driver
  consumes a single-sequence FASTA / `.seq` input and writes a
  classic connect-table `.ct` describing the predicted base
  pairings, plus a PostScript / PDF structure plot and a per-run
  `.out` log. Single-binary subprocess shape (sister to Phase 18
  BWA): the CLI uses mfold's `KEY=VALUE`-style invocation `mfold
  SEQ=<sequence> NA=RNA T=<temperature> [extras...]` rather than
  POSIX-style `--seq` / `-T` flags. Schema knobs: `sequence`
  (`.fa` / `.fasta` / `.seq` single-sequence input; required),
  `output_basename` (filename stem mfold uses for `.ct` / `.ps` /
  `.out` outputs; required, non-empty), `temperature` (`f64`,
  finite; folding temperature in Celsius matching mfold's `T=`
  convention; default 37.0 — physiological / standard reference
  temperature), `extra_args` (additional `KEY=VALUE` pairs
  appended after the canonical three so users can pin
  `MAX_LP=...`, `IONS=...`, suboptimal-energy budget overrides,
  etc.). `prepare()` resolves `sequence` against the case
  directory when relative, validates the file exists on disk
  (returns `InvalidCase` with a helpful message when missing),
  and composes the `KEY=VALUE` invocation. `collect()` walks the
  workdir for `*.ct` (`Tabular`, "mfold connect-table" — the
  canonical mfold-format base-pair file consumed by every
  downstream RNA secondary-structure analyzer), `*.ps` / `*.pdf`
  (`Native`, "mfold structure plot" — PostScript / PDF rendering
  of the predicted secondary structure), and `*.out` (`Log`,
  "mfold log" — the per-run mfold log). **License callout:**
  mfold is licensed for academic / non-commercial use only —
  commercial use requires a separate license from RNA Structure
  Software, Inc. The probe surfaces an `"academic"` /
  `"non-commercial"`-keyworded license-awareness warning in
  `ProbeReport.warnings` reminding users to confirm their use
  case before redistributing folds or derived data, sister to the
  Phase 28 ViennaRNA / NUPACK and Phase 5.6 NAMD / Phase 23 VMD
  pattern. Probe via `find_on_path(&["mfold", "UNAFold.pl"])`
  (the modern UNAFold distribution renames the launcher to
  `UNAFold.pl` while keeping the original `mfold` symlink for
  backwards compat). Version range `3.8.0..4.0.0` (mfold 3.8 is
  the modern stable line; upper bound 4.0 reserves room for an
  eventual major bump). `bio.mfold.fold` ribbon capability.
- **EternaFold** — the Eterna project's ML-aware RNA folder
  (MIT). EternaFold rides on a half-decade of crowd-sourced
  Eterna gameplay-puzzle data plus the modern thermodynamic +
  machine-learning corpora to train a maximum-expected-accuracy
  (MEA) predictor that's competitive with ViennaRNA + RNAstructure
  on benchmark sets while doing better on the long tail of
  experimentally-characterised RNA puzzle structures the
  classical thermodynamic models miss. The canonical interface is
  the `arnie` Python wrapper (Das lab) — a unified RNA-folder
  front-end that proxies to ViennaRNA / RNAstructure / NUPACK /
  EternaFold / LinearFold under a single Python API — so
  EternaFold ships in Valenx via the standard Python-script
  subprocess pattern (sister to Phase 17 Biopython, Phase 19.5
  Scanpy, Phase 28 NUPACK, Phase 41 pydna, Phase 43 DNA Chisel).
  The user supplies a `.py` script referenced from
  `[bio.eternafold].script` in `case.toml` that imports `arnie`
  and reads `valenx_params.json` for the parsed knobs. Schema
  knobs: `script` (path to user-supplied Python script; required,
  `.py` enforced), `python` (interpreter name; default
  `"python3"`), `input_fasta` (`Option<PathBuf>` — optional input
  FASTA the script can fold; `None` when the script generates
  the sequence inline), `output_basename` (filename stem;
  required, non-empty). `prepare()` enforces the `.py` extension,
  routes script + optional input_fasta through `confined_join`
  to stage them safely in the workdir under their original
  filenames, then writes a flat hand-rolled `valenx_params.json`
  containing `output_basename` always plus `input_fasta` (staged
  filename) only when set — the key is omitted entirely when
  `None` rather than emitted as `null`, matching the hand-rolled
  JSON convention the rest of the bio adapters use (Phase 19.6
  Seurat / AnnData, Phase 27.5 ESM-IF, Phase 41 pydna, Phase 42
  Mol* / NGL, Phase 43 DNA Chisel / iCodon). `collect()` walks
  the workdir for `<output_basename>*.ct` (`Tabular`, "EternaFold
  connect-table"), `<output_basename>*.dot` (`Native`,
  "EternaFold dot-bracket" — the canonical Vienna-style dot-
  bracket secondary-structure notation), `<output_basename>*.csv`
  (`Tabular`, "EternaFold MEA / probabilities" — per-base
  pairing probabilities for downstream MEA / probabilistic
  decoding), and `*.log` (`Log`). Probe via Python on PATH then
  `<python> -c "import arnie"` — when the `import arnie` check
  fails the probe still returns `ok = true` with a targeted
  `"probe found python on PATH but could not import arnie —
  install via pip install arnie"` warning so users with Python
  ready but no `arnie` package see the install hint without
  failing the probe (sister to the Phase 19.5 scanpy / scvi /
  Phase 19.6 AnnData / Phase 5.6 HOOMD-blue / Phase 5.7 MDTraj
  / Phase 41 pydna / Phase 42 Mol* / NGL / Phase 43 DNA Chisel
  probe convention). Version range `1.3.0..2.0.0` (EternaFold
  1.3 is the modern stable line shipping the contemporary MEA
  predictor; upper bound 2.0 reserves room for an eventual
  major bump). `bio.eternafold.fold` ribbon capability.
- **LinearFold** — Baidu Research / Oregon State's beam-search
  linear-time RNA folder (Apache-2.0). LinearFold is the folding-
  only sister to Phase 43 LinearDesign — same beam-search core
  from the same group, applied to the inverse problem of "given
  a sequence, find the secondary structure" rather than "given
  a target protein, find the optimized mRNA". The decisive
  property is the linear (`O(N)` or `O(N · beam_size)`) per-
  nucleotide complexity, contrasted with the cubic `O(N^3)` cost
  of the classical Zuker / RNAstructure / ViennaRNA dynamic-
  programming folders — LinearFold scales to viral-genome-
  length sequences (~30 kb SARS-CoV-2, ~9 kb HIV, full-length
  pre-mRNAs) without the classical folders' polynomial blowup.
  LinearFold ships two model back-ends — `C` (CONTRAfold-style
  ML-trained scoring) and `V` (Vienna-style thermodynamic
  scoring) — selectable through the `model` knob. Single-binary
  subprocess shape (sister to LinearDesign / Phase 18 BWA / Phase
  32.5 Smoldyn / Phase 5 GROMACS) with a non-standard stdin
  contract: LinearFold reads the sequence from stdin and writes
  the predicted dot-bracket / energy lines to stdout. Schema
  knobs: `sequence` (`.fa` / `.fasta` / `.seq` single-sequence
  input; required — read in place, no staging), `output_basename`
  (filename stem; required, non-empty), `model` (default `"C"` —
  CONTRAfold-style; `"V"` selects ViennaRNA-style scoring),
  `beam_size` (`u32`, ≥ 1; beam-search width controlling the
  speed / accuracy tradeoff; default 100 — LinearFold's
  recommended default), `extra_args`. `prepare()` resolves
  `sequence` against the case directory when relative, validates
  the file exists on disk, and composes `linearfold -V
  <beam_size> [extras...]` (or `-C` when `model == "C"`) with
  stdin redirected from `sequence` and stdout redirected to
  `<output_basename>.txt`. `collect()` walks the workdir for
  `<output_basename>*.txt` (`Tabular`, "LinearFold structure
  output" — LinearFold's per-line `<sequence> <dot-bracket>
  (<energy>)` capture) and `*.log` (`Log`). Probe via
  `find_on_path(&["linearfold"])` — when `linearfold` isn't on
  PATH but Python is the probe surfaces a targeted `"LinearFold
  not found on PATH; clone https://github.com/LinearFold/
  LinearFold and add the bin directory to PATH"` warning so
  users see the install hint immediately. Version range
  `1.0.0..2.0.0` (LinearFold 1.x is the modern stable line
  shipping the contemporary `-V` / `-C` model selection; upper
  bound 2.0 reserves room for an eventual major bump).
  `bio.linearfold.fold` ribbon capability.

### Canonical types

**No new canonical types.** All three adapters consume user-
supplied inputs (mfold + LinearFold single-sequence FASTA / `.seq`
files, EternaFold Python scripts that import `arnie` plus optional
input FASTA) and emit user-readable artifacts (mfold's classic
connect-table `.ct` + structure-plot `.ps` / `.pdf` + per-run
`.out` log, EternaFold's connect-table `.ct` + dot-bracket `.dot`
+ MEA / probability `.csv` table, LinearFold's per-line
`<sequence> <dot-bracket> (<energy>)` `.txt` output) that the
unchanged `Results.artifacts` collection model surfaces directly.
The existing `valenx_bio::format::fasta` reader already inspects
sequence inputs for sequence count + identifiers + alphabet. A
first-class RNA-secondary-structure canonical type — a typed
connect-table / dot-bracket / per-base-probability representation
spanning all six RNA folders (Phase 28 ViennaRNA / RNAstructure /
NUPACK plus Phase 44.5 mfold / EternaFold / LinearFold) with
parsed base-pair / probability graphs — defers to a future phase
along with cross-folder structure-comparison viewers and
suboptimal-ensemble inspection CLIs.

### Headless CLIs

**No new CLIs.** mfold's `.ct` connect-tables + `.ps` / `.pdf`
plots, EternaFold's `.ct` / `.dot` / `.csv` outputs, and
LinearFold's per-line `.txt` output are all standard tabular /
plain-text formats inspectable in any editor or through the
user's downstream Python pipeline (`pandas`, `numpy`, `Biopython`,
`arnie`). Input FASTAs are inspectable through the existing
Phase 17 `valenx-fasta` CLI. A canonical RNA-folding CLI —
folder-comparison diffing, MEA / partition-function probability
inspection, suboptimal-ensemble traversal — defers to a future
phase along with the canonical type.

## Domain expansion

Phase 44.5 is a **sister-adapter expansion of the Phase 28 RNA
structure trio** (ViennaRNA / RNAstructure / NUPACK) — the same
RNA secondary-structure prediction surface broadened with three
more established folders that cover the corners Phase 28 doesn't
reach. ViennaRNA is the most-cited classical thermodynamic folder;
RNAstructure is the Mathews-lab BSD-3-Clause classic; NUPACK is
Caltech's Python-driven academic-license suite. With Phase 44.5
the RNA-folding surface in Valenx covers all three canonical
shapes — single-binary CLI thermodynamic / DP folders (mfold +
LinearFold), Python-package ML-aware folders via `arnie`
(EternaFold), and the original Phase 28 trio — six folders total
spanning the academic Zuker classic, the BSD modern Mathews-lab
classic, the Caltech Python-driven academic suite, the modern
thermodynamic web-cited workhorse (ViennaRNA), the ML-aware
crowd-sourced Eterna predictor, and the linear-time beam-search
folder for full-genome sequences. The Phase 28 ViennaRNA's
academic-use license caveat is inherited by mfold (and remains
contrasted with RNAstructure's BSD / EternaFold's MIT /
LinearFold's Apache-2.0).

## What landed early

The implementation rode subagent-driven-development across 4
discrete implementation commits (3 adapters plus the registry +
init-template rollup) plus this docs pass — each landing one
adapter, the registry rollup, the init-template extension, or
the documentation pass. Every commit kept workspace clippy +
rustdoc clean.

## Acceptance checklist

- [x] `valenx-adapter-mfold` adapter ships with case-input parser
      + 4 lib tests + 3 case-input tests covering parses-minimal /
      parses-with-overrides / rejects-empty-sequence, plus the
      single-binary subprocess shape that composes `mfold
      SEQ=<sequence> NA=RNA T=<temperature> [extras...]` with
      `sequence` resolved against the case directory and validated
      on disk, the `find_on_path(["mfold", "UNAFold.pl"])` probe,
      and the `"academic"` / `"non-commercial"`-keyworded
      license-awareness warning surfaced in `ProbeReport.warnings`
- [x] `valenx-adapter-eternafold` adapter ships with case-input
      parser + 4 lib tests + 3 case-input tests covering parses-
      minimal / parses-with-input-fasta / rejects-non-py-script,
      plus the Python-script subprocess shape that enforces `.py`,
      routes script + optional input_fasta through `confined_join`,
      writes `valenx_params.json` with `output_basename` always
      plus `input_fasta` (staged filename) only when set — key
      omitted entirely when `None` rather than emitted as `null`,
      matching the hand-rolled JSON convention the rest of the
      bio adapters use, plus the Python on PATH + `import arnie`
      probe with `"could not import arnie"` warning when the
      import fails
- [x] `valenx-adapter-linearfold` adapter ships with case-input
      parser + 4 lib tests + 3 case-input tests covering parses-
      minimal / parses-with-overrides / rejects-bad-beam-size,
      plus the single-binary subprocess shape with non-standard
      stdin contract — `sequence` resolved against the case
      directory and validated on disk (read in place, no staging),
      stdin redirected from the sequence file, stdout captured to
      `<output_basename>.txt`, `linearfold -V <beam_size>` (or
      `-C`) composed with `[extras...]`, and the
      `find_on_path(["linearfold"])` probe with the `"clone
      https://github.com/LinearFold/LinearFold and add the bin
      directory to PATH"` warning when Python is on PATH but
      `linearfold` is missing
- [x] All 3 adapters wired into `valenx-app::init_registry` —
      live adapter count moves from 128 to **131** as part of the
      13-adapter / 4-phase rollup that takes the headline live-
      adapter total to **141**, rounding out the RNA secondary-
      structure folding surface that Phase 28 ViennaRNA /
      RNAstructure / NUPACK opened
- [x] 3 RNA-folding templates in `valenx-init` (`mfold`,
      `eternafold`, `linearfold`), all round-tripping through
      `valenx-validate` (cross-binary roundtrip now sweeps **136
      templates** clean alongside the Phase 35.5 / 35.6 / 45
      sister rollups)
- [x] STATUS.md / README.md / ROADMAP.md / CHANGELOG.md /
      QUICKSTART.md updated
- [ ] **Phase 44.6** — sister-adapter expansion of Phase 44.5:
      RNAcontextfold (the contextfold ML predictor sister to
      EternaFold; defer pending upstream activity), CONTRAfold
      (the original DP-based ML scoring model EternaFold + LinearFold
      build on; defer — modern EternaFold + LinearFold-C cover the
      same ground), CentroidFold (centroid-style decoder defer),
      RNAfold-MEA (the ViennaRNA `RNAfold --MEA` mode that's already
      reachable through Phase 28's `extra_args`; defer — the user
      can already pass `--MEA`), MXfold2 (deep-learning RNA folder;
      defer pending licensing review of the model checkpoint),
      RNAstructure's `partition` / `MaxExpect` modes (already
      reachable through the existing Phase 28 RNAstructure adapter's
      `extra_args`; out of scope as a separate adapter). Out of
      scope for this beachhead.

## Success metrics

| Metric                                                | Target          |
|-------------------------------------------------------|-----------------|
| New RNA-folding adapter (template + tests)            | 1 day per       |
| Fold-mfold → fold-EternaFold → fold-LinearFold loop across 3 tools | < tool baseline |

## Leads into

Phase 44.5 rounds out the RNA secondary-structure folding surface
that the user's bio / chemistry spec called out alongside the
existing Phase 28 trio (ViennaRNA / RNAstructure / NUPACK).
Combined with the existing optimize-codons → design-mRNA →
predict-stability → render-Mol* → render-NGL → run-Galaxy-
workflow → run-WDL → run-CWL → run-Nextflow → run-Snakemake →
design-plasmid → view-alignment → process-image → segment-cells
→ classify-pixels → simulate-pathway → expand-rules → grow-
tissue → diffuse-particles → trace-MCell-trajectories → simulate-
MD → analyze-trajectory → reweight-free-energy → fit-ENM → run-
cpptraj-script → predict-structure → fold-RNA → analyze-DNA-
geometry → infer-tree-ML → infer-tree-Bayesian → simulate-popgen
→ analyze-trees → reconstruct-3D → design-protein → validate
loop, the **fold-mfold → fold-EternaFold → fold-LinearFold →
optimize-codons → design-mRNA → predict-stability → render-Mol*
→ render-NGL → run-Galaxy-workflow → run-WDL → run-CWL → run-
Nextflow → run-Snakemake → design-plasmid → view-alignment →
process-image → segment-cells → classify-pixels → simulate-
pathway → expand-rules → grow-tissue → diffuse-particles →
trace-MCell-trajectories → simulate-MD → analyze-trajectory →
reweight-free-energy → fit-ENM → run-cpptraj-script → predict-
structure → fold-RNA → analyze-DNA-geometry → infer-tree-ML →
infer-tree-Bayesian → simulate-popgen → analyze-trees →
reconstruct-3D → design-protein → validate** loop now spans six
RNA folders (the Phase 28 ViennaRNA / RNAstructure / NUPACK trio
plus Phase 44.5 mfold / EternaFold / LinearFold) feeding into the
existing Phase 43 mRNA design stack (DNA Chisel, LinearDesign,
iCodon — every codon-optimized mRNA deserves a multi-folder
secondary-structure check), the Phase 33 synthetic-biology
composition stack (pySBOL / j5 / Cello + Phase 41 pydna), and the
Phase 35 / 35.5 / 35.6 CRISPR / editing design stacks — all in
one Valenx shell with no glue code beyond the existing case-toml /
prepare / run / collect path.

The natural follow-up is **Phase 44.6** — the deferred RNA-
folding work called out above (RNAcontextfold pending upstream
activity, CONTRAfold if the modern ML predictors don't already
cover the use case, CentroidFold for the centroid-style decoder
surface, MXfold2 pending licensing review of the deep-learning
model checkpoint), slotting in alongside the existing mfold /
EternaFold / LinearFold adapters with the same single-binary CLI
shape (mfold / LinearFold sister tools), Python-script subprocess
shape (EternaFold sister tools), or new shape if upstream tools
require something novel. See the out-of-scope section of
`docs/superpowers/plans/2026-05-04-rna-folding-expansion.md` for
the full follow-up phase list.
