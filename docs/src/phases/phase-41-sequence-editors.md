# Phase 41 — Sequence editors

**Status:** 🟢 Live — pydna + Jalview open the **first plasmid-
design / alignment-viewer domain** in Valenx alongside the
Phase 5.5 / 5.6 / 5.7 / 17 / 17.5 / 17.7 / 18 / 25 / 27 / 27.5
/ 27.6 / 28 / 29 / 30 / 30.5 / 31 / 32 / 32.5 / 33 / 34 / 35 /
36 / 38 / 39 / 40 MD-engines + MD-analysis-expansion + biology
+ structure-prediction + protein-design + RNA-structure +
population-genetics + phylogenetics + Bayesian-phylogenetics +
read-simulator + systems-biology + spatial-stochastic +
synthetic-biology + docking + CRISPR-design + cryo-EM +
Rosetta-family + DNA-geometry + microscopy beachheads and the
Phase 24 cheminformatics expansion.

## Goal

Open the **plasmid / clone-design + alignment-viewer** domain
in Valenx with two established open-source tools that span the
sequence-editor tradeoff space — a Python plasmid-design
library that handles PCR primer design, restriction-enzyme
digests, and Gibson / Golden-Gate assembly programmatically
(pydna, Bjorn Johansson's BSD-3-Clause library that's the de-
facto Python choice for cloning automation), and the canonical
Java alignment viewer with a headless mode for batch image /
format conversion (Jalview, the Barton group's GPL-3.0 viewer
that's been the reference alignment viewer in molecular biology
labs since the 2000s and supports headless operation for
unattended pipeline integration). pydna follows the established
Phase 17 Biopython Python-script subprocess shape: the user
supplies a Python script that imports the upstream package and
reads `valenx_params.json` for the parsed knobs. Jalview is
JAR-distributed (no `jalview` launcher binary on PATH); the
user supplies the absolute path to the JAR via case input, and
we probe `java` itself rather than the JAR — same JAR-
distribution shape Phase 33 j5 / Cello use. Phase 41 sits
numerically after Phase 40 microscopy and ships chronologically
right after Phase 40 — same chronological-vs-numerical
convention used for Phase 17.5 / 24 / 28 / 31 / 35 / 39 / 5.5 /
5.6 / 5.7 / 32.5 / 40.

## Capability inventory

### Live adapters (2)

- **pydna** — Bjorn Johansson's Python plasmid / clone-design
  library (BSD-3-Clause). pydna handles the long tail of
  cloning operations programmatically: PCR primer design,
  restriction-enzyme digests, Gibson assembly, Golden-Gate
  assembly, sequence-overlap detection, ligation simulation,
  primer Tm calculation, cloning-strategy validation. The
  canonical use case is "I need to assemble these N parts into
  this target construct — what primers should I order, what
  enzymes should I cut with, and does the end-product match my
  target sequence?" — pydna replaces hours of manual work in
  ApE / SnapGene / Vector NTI with a few dozen lines of Python.
  Python-script subprocess shape (sister to Phase 17 Biopython,
  Phase 19.5 Scanpy, Phase 33 pySBOL): the user supplies a
  Python script referenced from `[bio.pydna].script` in
  `case.toml` that imports `pydna` and reads `valenx_params.json`
  for the parsed knobs. Schema knobs: `script` (path to user-
  supplied Python script; required, `.py` enforced), `python`
  (interpreter name; default `"python3"`), `input_genbank`
  (`Option<PathBuf>` — optional starting GenBank file the
  script can use as the parent / template construct; `None`
  when the script generates the design from scratch),
  `output_basename` (filename stem the user's script uses for
  outputs — surfaced here so collect() can label artifacts
  uniformly; required, non-empty). `prepare()` enforces the
  `.py` extension on the script, resolves `script` and the
  optional `input_genbank` against the case directory when
  relative, stages both into the workdir under their original
  filenames so the script can resolve them via relative paths,
  then writes a flat `valenx_params.json` containing
  `output_basename` always plus `input_genbank` (staged
  filename) only when set — the key is omitted entirely when
  `None` rather than emitted as `null`, matching the hand-
  rolled JSON convention the rest of the bio adapters use
  (Phase 19.6 Seurat / AnnData, Phase 27.5 ESM-IF). `collect()`
  walks the workdir for `<output_basename>*.gb` (`Native`,
  "pydna GenBank file" — the canonical GenBank `.gb` annotated-
  sequence format every downstream tool reads),
  `<output_basename>*.genbank` (`Native`, "pydna GenBank file"
  — same kind, alternate extension),
  `<output_basename>*.fasta` (`Native`, "pydna FASTA" — for
  scripts that emit FASTA-format primer / fragment lists),
  `<output_basename>*.csv` (`Tabular`, "pydna table" — for
  scripts that emit per-fragment / per-primer summary tables),
  and `*.log` (`Log`). Probe via Python on PATH with an `import
  pydna` check — when the import fails the probe still returns
  `ok = true` with a warning so users with pydna installed
  under a non-standard module name aren't blocked (sister to
  the Phase 19.5 scanpy / scvi / Phase 19.6 AnnData / Phase 5.6
  HOOMD-blue / Phase 5.7 MDTraj probe convention). Version
  range `5.0.0..7.0.0` (pydna 5.x is the modern stable line
  pairing with Biopython 1.8x; pydna 6.x ships in 2024 with the
  contemporary Gibson / Golden-Gate assembly improvements;
  upper bound 7.0 reserves room for an eventual major bump).
  `bio.pydna.design` ribbon capability.
- **Jalview** — the Barton group's Java alignment viewer
  (GPL-3.0). Jalview has been the reference alignment viewer in
  molecular biology labs since the 2000s — multiple-sequence
  alignment viewing + editing, conservation / consensus /
  occupancy plots, structural overlays via Jmol / Chimera links,
  per-column annotations, the canonical interactive front-end
  for the MSA outputs every Phase 18 / 18.5 / 18.6 / 18.7
  aligner (BWA / minimap2 / MAFFT / MUSCLE / Clustal Omega /
  T-Coffee) emits. Crucially Jalview ships a **headless mode**
  (`-nodisplay`) for batch image / format conversion: the user
  feeds an alignment in, picks an output format (PNG image,
  HTML report, SVG vector graphic, FASTA / Clustal alignment
  re-export), and Jalview writes the requested artifact without
  opening its GUI. JAR-distributed single-binary subprocess
  shape (sister to Phase 33 j5 / Cello): no `jalview` launcher
  binary on PATH; the user supplies the absolute path to the
  JAR via `[bio.jalview].jar` in `case.toml`. The CLI is
  `java -jar <jar> -nodisplay -open <input> -<output_format>
  <basename>.<ext> [extras...]`. Schema knobs: `jar` (absolute
  path to the Jalview jar; required), `input` (alignment input
  file — `.fa` / `.aln` / `.clustal` / `.stockholm` and friends
  Jalview reads natively; required), `output_basename`
  (filename stem the adapter pins as the Jalview output target;
  required, non-empty), `output_format` (string; default
  `"png"` — selectable from `"png"` / `"html"` / `"svg"` /
  `"fasta"` / `"clustal"` for the canonical headless output
  formats Jalview ships), `extra_args` (additional CLI
  arguments appended after the output target). `prepare()`
  resolves `jar` and `input` against the case directory when
  relative, validates each file exists on disk (returns
  `InvalidCase` with a helpful message when missing), derives
  the output extension from `output_format` (png → `.png`, html
  → `.html`, svg → `.svg`, fasta → `.fasta`, clustal → `.aln`,
  default → use the format string itself as the extension),
  composes the `java -jar` invocation. `collect()` walks the
  workdir for `<output_basename>*.png` (`Native`, "Jalview
  alignment image" — the canonical headless PNG render Jalview
  writes for unattended pipelines that want a static image of
  the alignment), `<output_basename>*.svg` (`Native`, "Jalview
  SVG" — vector alternative for downstream LaTeX / Inkscape
  workflows), `<output_basename>*.html` (`Native`, "Jalview
  HTML" — interactive HTML report), `<output_basename>*.fasta`
  (`Native`, "Jalview FASTA" — re-exported FASTA alignment
  for downstream phylogenetics / variant-calling pipelines),
  `<output_basename>*.aln` (`Tabular`, "Jalview alignment" —
  re-exported Clustal-format alignment), and `*.log` (`Log`).
  Probe via `find_on_path(&["java"])` — Jalview's version comes
  from the jar itself, not from `java`, so we surface no
  version here; the user pins the Jalview release implicitly by
  the jar they point at (same shape as Phase 33 j5 / Cello).
  Version range `2.11.0..3.0.0` (Jalview 2.11 (2022) is the
  modern stable line shipping the contemporary headless
  improvements; upper bound 3.0 reserves room for an eventual
  major bump). `bio.jalview.view` ribbon capability.

### Canonical types

**No new canonical types.** Both adapters consume user-supplied
inputs (pydna Python composition scripts + optional starting
GenBank files, Jalview alignment files in the formats Jalview
reads natively + the Jalview jar) and emit user-readable
artifacts (pydna `.gb` / `.genbank` annotated-sequence files +
`.fasta` primer / fragment lists + `.csv` per-fragment summary
tables, Jalview PNG / SVG / HTML alignment renders + re-
exported FASTA / Clustal alignments) that the unchanged
`Results.artifacts` collection model surfaces directly. The
existing `valenx_bio::format::fasta` reader already inspects
collected FASTA outputs for sequence counts; the existing
`valenx_bio::format::pdb` reader is structurally adjacent for
the cross-pipeline cases where pydna feeds into structure
prediction. A first-class plasmid-design canonical type — a
typed annotated-sequence representation spanning pydna GenBank
output + Jalview alignment edits, with parsed feature /
restriction-site / primer-binding-site graphs — defers to a
future phase along with plasmid-map visualizers and per-
restriction-site interactive overlays.

### Headless CLIs

**No new CLIs.** pydna's `.gb` / `.genbank` annotated-sequence
files, `.fasta` primer / fragment lists, and `.csv` per-fragment
summary tables are all standard formats inspectable through the
user's downstream Python pipeline (`Bio.SeqIO`,
`pydna.parsers.parse_primers`, `pandas`); Jalview's PNG / SVG /
HTML / FASTA / Clustal outputs are inspectable in any image
viewer or alignment viewer. A canonical sequence-editor CLI —
plasmid-construct inspection, alignment-edit diffing,
restriction-site comparison — defers to a future phase along
with the canonical type.

## Domain milestone

Phase 41 is the **first plasmid-design / alignment-viewer
domain** to land in Valenx. The biology adapter family started
with Phase 17 (foundation — sequence / structure / trajectory
canonical types + classical MD + cheminformatics scripts) and
expanded through Phase 5.5 / 5.6 / 5.7 / 17.5 / 17.7 / 18 /
18.5 / 18.6 / 18.7 / 19 / 19.5 / 19.6 / 20 / 22 / 23 / 24 / 25
/ 27 / 27.5 / 27.6 / 28 / 29 / 30 / 30.5 / 31 / 32 / 32.5 / 33
/ 34 / 35 / 36 / 38 / 39 / 40 to cover MD-trajectory analysis
expansion, bio MD engines, MDTraj, sequence prediction,
structure-tools expansion, alignment, RNA-seq, alignment-
toolkit expansion, variant calling, single-cell, single-cell
expansion, transcript quantification, workflow orchestration,
molecular viewers, cheminformatics, quantum chemistry, protein
design, protein-design expansion, EvolutionaryScale models, RNA
structure, population genetics, phylogenetics, Bayesian
phylogenetics, sequencing read simulation, systems biology,
spatial-stochastic simulators, synthetic biology, small-
molecule docking, CRISPR design, cryo-EM reconstruction,
Rosetta protein modeling, DNA structural geometry, and
microscopy — but until Phase 41 the plasmid-design / alignment-
viewer surface (Python plasmid / clone-design library, headless
alignment viewer for batch image / format conversion) was
absent. Phase 41 closes that gap with two established open-
source tools spanning the sequence-editor tradeoff space —
pydna at the Python plasmid-design end, and Jalview as the
canonical headless-mode-capable alignment viewer that closes
the loop on the entire Phase 18 / 18.5 / 18.6 / 18.7 alignment
beachhead by giving the MSA outputs a publication-quality
visual front end.

## What landed early

The implementation landed across 3
discrete implementation commits (2 adapters plus the registry +
init-template rollup) plus this docs pass — each landing one
adapter, the registry rollup, the init-template extension, or
the documentation pass. Every commit kept workspace clippy +
rustdoc clean.

## Acceptance checklist

- [x] `valenx-adapter-pydna` adapter ships with case-input
      parser + 4 lib tests + 3 case-input tests covering parses-
      minimal / parses-with-input-genbank / rejects-non-py-
      script, plus the Python-script subprocess shape that
      enforces `.py`, stages script + optional input_genbank,
      writes `valenx_params.json` with `output_basename` always
      plus `input_genbank` (staged filename) only when set —
      key omitted entirely when `None` rather than emitted as
      `null`, matching the hand-rolled JSON convention the rest
      of the bio adapters use
- [x] `valenx-adapter-jalview` adapter ships with case-input
      parser + 4 lib tests + 3 case-input tests covering parses-
      minimal / parses-with-format-overrides / rejects-empty-
      jar, plus the JAR-distributed single-binary subprocess
      shape that probes `java` on PATH and composes `java -jar
      <jar> -nodisplay -open <input> -<output_format>
      <basename>.<ext> [extras...]` with the output extension
      derived from `output_format` (png → .png, html → .html,
      svg → .svg, fasta → .fasta, clustal → .aln, default →
      format string itself as extension)
- [x] Both adapters wired into `valenx-app::init_registry` —
      live adapter count moves from 117 to **119** (alongside
      the Phase 32.5 spatial-stochastic pair and Phase 40
      microscopy trio that bring the total to 119), opening the
      first plasmid-design / alignment-viewer domain to ship in
      Valenx
- [x] 2 sequence-editor templates in `valenx-init` (`pydna`,
      `jalview`), all round-tripping through `valenx-validate`
      (cross-binary roundtrip now sweeps **115 templates** clean
      alongside the Phase 32.5 spatial-stochastic pair and
      Phase 40 microscopy trio)
- [x] STATUS.md / README.md / ROADMAP.md / CHANGELOG.md /
      QUICKSTART.md updated
- [ ] **Phase 41.5** — sister-adapter expansion of Phase 41:
      ApE (open-source plasmid editor; macOS / Windows GUI
      app with no headless mode — defer until upstream ships
      one), SnapGene Viewer (free read-only GUI; commercial
      editor; no headless mode — defer), Benchling (web-based
      collaborative cloning; out of scope as a web service),
      MEGA (canonical phylogenetics / sequence-editing GUI;
      `megacc` headless CLI exists — defer to a sister
      adapter), CIPRES (web-portal phylogenetics; out of scope
      as a web service), AliView (lightweight Java alignment
      viewer sister to Jalview; defer to sister-adapter
      expansion), SeaView (alternative alignment viewer with
      tree-editing; defer). Out of scope for this beachhead.

## Success metrics

| Metric                                                | Target          |
|-------------------------------------------------------|-----------------|
| New sequence-editor adapter (template + tests)        | 1 day per       |
| Design-plasmid → view-alignment loop across 2 tools   | < tool baseline |

## Leads into

Phase 41 opens the plasmid-design / alignment-viewer domain
that the user's bio / chemistry spec called out alongside the
Phase 17 / 17.5 / 17.7 biology + structure-prediction stack,
the Phase 18 / 18.5 / 18.6 / 18.7 alignment + search beachhead,
and the Phase 33 synthetic-biology trio (pySBOL / j5 / Cello —
the canonical SBOL composition + DNA assembly + Verilog →
circuit compilation surface that pydna sits naturally adjacent
to on the cloning-automation side). Combined with the existing
process-image → segment-cells → classify-pixels → simulate-
pathway → expand-rules → grow-tissue → diffuse-particles →
trace-MCell-trajectories → simulate-MD → analyze-trajectory →
reweight-free-energy → fit-ENM → run-cpptraj-script → predict-
structure → fold-RNA → analyze-DNA-geometry → infer-tree-ML →
infer-tree-Bayesian → simulate-popgen → analyze-trees →
reconstruct-3D → design-protein → validate loop, the **design-
plasmid → view-alignment → process-image → segment-cells →
classify-pixels → simulate-pathway → expand-rules → grow-tissue
→ diffuse-particles → trace-MCell-trajectories → simulate-MD →
analyze-trajectory → reweight-free-energy → fit-ENM → run-
cpptraj-script → predict-structure → fold-RNA → analyze-DNA-
geometry → infer-tree-ML → infer-tree-Bayesian → simulate-
popgen → analyze-trees → reconstruct-3D → design-protein →
validate** loop now spans two sequence-editor tools (pydna,
Jalview) feeding into the existing Phase 40 microscopy trio
(Fiji, CellProfiler, Ilastik), the Phase 5 / 5.6 GROMACS /
LAMMPS / NAMD / sander / HOOMD-blue MD engines, the Phase 5.5
/ 5.7 / 17 PLUMED / ProDy / cpptraj / MDTraj / MDAnalysis post-
MD analysis stack, the Phase 17 / 17.5 / 17.7 prediction stack
(ESMFold, OpenFold, AlphaFold 2/3, ColabFold, RoseTTAFold,
OmegaFold, FoldSeek), the Phase 18 / 18.5 / 18.6 / 18.7
alignment + search surface (BWA, minimap2, MAFFT, MUSCLE,
HMMER, samtools, Bowtie2, MMseqs2, DIAMOND, HISAT2, STAR,
BLAST+, Clustal Omega, T-Coffee — Jalview is the natural
viewer for every MSA the latter five emit), the Phase 28
RNA-structure tools (ViennaRNA, RNAstructure, NUPACK), the
Phase 29 population-genetics trio (SLiM, msprime, tskit), the
Phase 30 phylogenetic-tree builders (IQ-TREE, RAxML-NG,
FastTree), the Phase 30.5 Bayesian-phylogenetics pair (BEAST 2,
MrBayes), the Phase 32 systems-biology surface (COPASI,
BioNetGen, PhysiCell), the Phase 32.5 spatial-stochastic pair
(Smoldyn, MCell), the Phase 33 synthetic-biology trio (pySBOL,
j5, Cello — pydna sits naturally adjacent on the cloning-
automation side), the Phase 34 docking pair (AutoDock Vina,
AutoDock 4), the Phase 35 CRISPR-design tools (CHOPCHOP,
CRISPOR, Cas-OFFinder), the Phase 36 cryo-EM reconstruction
tools (RELION, EMAN2, CTFFIND), the Phase 38 Rosetta-family
adapters (Rosetta, PyRosetta), and the Phase 39 DNA-structural-
geometry tools (X3DNA, Curves+, DSSR) — all in one Valenx
shell with no glue code beyond the existing case-toml / prepare
/ run / collect path.

The natural follow-up is **Phase 41.5** — the deferred
sequence-editor work called out above (MEGA / `megacc` for the
canonical phylogenetics + sequence-editing CLI sister to
Jalview, AliView / SeaView as additional alignment viewers
sister to Jalview on the GUI side, ApE / SnapGene Viewer if
upstream ever ships headless modes, Benchling / CIPRES if the
web-service shape becomes acceptable for the registry pattern),
slotting in alongside the existing pydna + Jalview adapters
with the same Python-script subprocess shape (pydna sister
tools), JAR-distributed single-binary subprocess shape (Jalview
/ MEGA / AliView / SeaView sister tools), or app-launcher
subprocess shape (ApE / SnapGene sister tools if headless modes
land).
