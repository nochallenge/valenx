# Phase 17.7 — Structure tools expansion

**Status:** 🟢 Live — RoseTTAFold + OmegaFold + FoldSeek round
out the **protein structure prediction + structure search
surface** that Phase 17.5 ESMFold / OpenFold / AlphaFold 2 /
AlphaFold 3 + Phase 17 ColabFold opened.

## Goal

Sister-adapter expansion of the existing Phase 17.5 structure-
prediction beachhead and the Phase 17 ColabFold adapter. Round
out the protein structure prediction + structure search surface
with three more foundational open-source tools — Baker lab's
original 3-track structure-prediction network (RoseTTAFold, the
canonical pre-AlphaFold-3 sibling that established the 3-track
SE(3)-equivariant attention pattern), HelixonAI's single-sequence
structure predictor (OmegaFold, MSA-free like ESMFold but with a
larger pre-trained transformer backbone), and Steinegger lab's
3D-structure search tool (FoldSeek, the protein-3D analogue of
the Phase 18.5 MMseqs2 sequence search — both from the Steinegger
lab, both built on the same fast many-vs-many search core, but
FoldSeek encodes the per-residue 3D geometry as the canonical
"3Di alphabet" for sequence-search-style 3D matching at sequence-
search speed). RoseTTAFold + OmegaFold follow established Python-
script subprocess shapes (RoseTTAFold sister to Phase 17.5
ESMFold / OpenFold; OmegaFold ships its own CLI binary with
Python fallback). FoldSeek follows the established Phase 18.5
MMseqs2 single-binary CLI pattern. Phase 17.7 sits numerically
after Phase 17.5 (Phase 17.6 was reserved for a deferred
confidence-aware structure-ranking + mmCIF-reader work that
hasn't shipped yet) and ships chronologically right after Phase
5.7 MDTraj — same chronological-vs-numerical convention used for
Phase 17.5 / 24 / 28 / 31 / 35 / 39 / 5.5 / 5.6 / 5.7.

## Capability inventory

### Live adapters (3)

- **RoseTTAFold** — Baker lab's original 3-track structure-
  prediction network (MIT). RoseTTAFold is the canonical pre-
  AlphaFold-3 sibling that established the 3-track SE(3)-
  equivariant attention pattern — three concurrent attention
  tracks over the 1D sequence, the 2D pair-distance map, and the
  3D Cartesian backbone, with cross-track message passing
  refining all three jointly. RoseTTAFold drove the original
  Baker-lab + Anand-lab structure-prediction surge (2021), and
  remains in production use as one of the few open-source
  alternatives to AlphaFold for structure prediction at academic
  scale. Python-script subprocess shape (sister to Phase 17.5
  ESMFold / Phase 17 Biopython): the user supplies a Python
  script referenced from `[bio.rosettafold].script` in
  `case.toml` that imports the upstream RoseTTAFold inference
  module and reads `valenx_params.json` for the parsed knobs.
  Schema knobs: `script` (path to user-supplied Python script;
  required, `.py` enforced), `python` (interpreter name; default
  `"python3"`), `fasta` (input FASTA query sequence; required),
  `output_basename` (filename stem RoseTTAFold uses for the
  predicted PDB + confidence arrays; required, non-empty).
  `prepare()` enforces the `.py` extension, resolves `script`
  and `fasta` against the case directory when relative, stages
  both into the workdir under their original filenames so the
  script can resolve them via relative paths, then writes a flat
  `valenx_params.json` containing `output_basename` and the bare
  `fasta` filename, and builds `<python> <staged_script>`.
  `collect()` walks the workdir for `<output_basename>*.pdb`
  (`Native`, "RoseTTAFold predicted structure" — the predicted
  3D backbone with pLDDT-style per-residue confidence written
  into the B-factor column, lifted by the existing Phase 17 PDB
  reader without any structure-prediction-specific code path),
  `<output_basename>*.npz` (`Native`, "RoseTTAFold confidence
  arrays" — the per-residue + per-pair-distance confidence
  arrays the Baker lab ships alongside the PDB), and `*.log`
  (`Log`). Probe via `find_on_path(&["python3", "python"])` —
  deliberately doesn't try `import rosettafold` (RoseTTAFold is
  not a pip package, it's a clone-from-GitHub install with heavy
  ML deps). Pushes a probe warning into `ProbeReport.warnings`
  whenever Python is detected: "RoseTTAFold model weights +
  dependencies not bundled — clone
  https://github.com/RosettaCommons/RoseTTAFold and follow the
  install README". Version range `1.0.0..3.0.0` (RoseTTAFold 1.x
  is the original 2021 release; RoseTTAFold 2 (RoseTTAFold All-
  Atom) is the late-2023 follow-up that adds nucleic-acid + small-
  molecule support; upper bound 3.0 reserves room for an eventual
  major). `bio.rosettafold.predict` ribbon capability.
- **OmegaFold** — HelixonAI's single-sequence protein-structure
  predictor (Apache-2.0). OmegaFold is MSA-free like ESMFold but
  uses a larger pre-trained transformer backbone trained on a
  much wider sequence corpus — it doesn't need an MSA, so it
  works on single sequences (sister to ESMFold from Phase 17.5
  in that respect) but routinely matches AlphaFold-2-with-MSA
  quality on orphan / synthetic / fast-evolving sequences where
  MSA-based methods struggle. OmegaFold ships its **own CLI
  binary** (`omegafold <fasta> <output_dir>`) and falls back to
  `<python> -m omegafold ...` when the CLI launcher isn't on
  PATH but Python is. Schema knobs: `fasta` (input FASTA query
  sequence; required), `output_basename` (workdir-relative
  output directory name OmegaFold writes the predicted PDBs
  under; required, non-empty), `python` (interpreter name;
  default `"python3"`; used only as fallback when the OmegaFold
  CLI isn't on PATH), `model_dir` (`Option<PathBuf>` — optional
  pre-downloaded model checkpoint directory; OmegaFold defaults
  to ~/.cache/omegafold_ckpt when omitted). `prepare()` builds
  `omegafold <fasta> <output_basename> [--model <model_dir>]`
  with the FASTA passed by absolute path (NOT staged into the
  workdir — OmegaFold reads it once then writes everything else
  into the output directory). `collect()` walks one level deep
  into the `<output_basename>/` subdirectory for `*.pdb`
  (`Native`, "OmegaFold predicted structure" — the predicted 3D
  backbone with per-residue confidence in the B-factor column,
  lifted by the existing Phase 17 PDB reader) and `*.json`
  (`Log`, "OmegaFold metadata" — the per-prediction metadata
  sidecar OmegaFold writes alongside each PDB), plus the
  workdir-top-level `*.log` (`Log`). Probe via
  `find_on_path(&["omegafold", "python3", "python"])` — surfaces
  a warning if `omegafold` itself isn't on PATH but Python is
  ("OmegaFold CLI not found on PATH; install via pip install
  git+https://github.com/HeliXonProtein/OmegaFold.git"). No
  academic-license caveat (Apache-2.0). Version range
  `1.0.0..2.0.0` (OmegaFold 1.x is the modern release line;
  upper bound 2.0 reserves room for a future major).
  `bio.omegafold.predict` ribbon capability.
- **FoldSeek** — Steinegger lab's protein-structure search via
  the 3Di alphabet (GPL-3.0). FoldSeek is the **3D analogue of
  Phase 18.5 MMseqs2** — both from the Steinegger lab, both
  built on the same fast many-vs-many search core, but FoldSeek
  encodes the per-residue 3D geometry as a custom "3Di
  alphabet" (a 20-letter alphabet over local backbone geometry
  patterns, designed so structural matches have high 3Di
  alphabet identity) and runs 3Di-vs-3Di comparisons at
  sequence-search speed. The result is a structure search tool
  that finds structural homologs at PDB-scale in seconds rather
  than the hours / days HMM-based or geometry-based search tools
  take. FoldSeek single-binary subprocess shape (sister to Phase
  18.5 MMseqs2 / Phase 18 BWA): the CLI is `foldseek easy-search
  <query> <database> <basename>.m8 tmp_<basename> --threads <N>
  [extras...]` where `tmp_<basename>` is a per-run temp
  directory FoldSeek requires. Schema knobs: `query` (`.pdb` /
  `.cif` query structure; required), `database` (FoldSeek
  database path prefix — the user supplies the path stem and
  FoldSeek resolves the `<prefix>_*` sidecar files itself;
  required), `output_basename` (filename stem FoldSeek uses for
  the `.m8` BLAST-style hit table; required, non-empty),
  `threads` (`u32`, default 1), `extra_args`. `prepare()`
  resolves `query` against the case directory when relative,
  validates the `database` parent directory exists on disk (the
  database files themselves use the prefix convention so we
  cannot validate them by name — same shape as Phase 18.7
  BLAST+'s `database` validation), and composes the invocation
  with the per-run temp directory pinned to `tmp_<basename>`.
  `collect()` walks the workdir for `<output_basename>.m8`
  (`Tabular`, "FoldSeek search results" — the canonical BLAST-
  style M8 hit table format every downstream FoldSeek pipeline
  reads) and `*.log` (`Log`). The temp directory is not surfaced
  in artifacts — it's intermediate. Probe via
  `find_on_path(&["foldseek"])`. Version range `8.0.0..10.0.0`
  (FoldSeek 8.x is the modern stable line that ships the easy-
  search command surface and the 3Di alphabet; upper bound 10.0
  reserves room for a future major). `bio.foldseek.search`
  ribbon capability.

### Canonical types

**No new canonical types.** All three adapters consume user-
supplied inputs (RoseTTAFold + OmegaFold FASTA queries, FoldSeek
PDB / CIF structures + database prefixes) and emit user-readable
artifacts (RoseTTAFold + OmegaFold `.pdb` predicted structures
with per-residue confidence in the B-factor column, RoseTTAFold
`.npz` confidence arrays, OmegaFold `.json` metadata sidecars,
FoldSeek `.m8` BLAST-style hit tables) that the unchanged
`Results.artifacts` collection model surfaces directly. The
predicted PDBs flow naturally through the existing
`valenx_bio::format::pdb` reader — same shape as Phase 17.5
ESMFold / OpenFold / AlphaFold 2 / AlphaFold 3 — and the per-
residue confidence rides the existing `Atom.b_factor` field. A
first-class structure-search canonical type — a typed structural-
homolog hit-list spanning FoldSeek + future structural-search
tools (DALI, TM-align, US-align) — defers to a future phase along
with cross-tool hit-table converters.

### Headless CLIs

**No new CLIs.** RoseTTAFold + OmegaFold predicted structures are
PDB files inspectable through the existing Phase 17
`valenx-pdb-info` CLI (chain / residue / atom counts + element
tally + B-factor lifted as confidence); FoldSeek's `.m8` BLAST-
style hit tables are tabular text inspectable in any editor or
through the user's downstream pipeline (BLAST `.m8` / `outfmt 6`
parsers). RoseTTAFold's `.npz` confidence arrays are inspectable
through `numpy.load`. A canonical structure-search CLI defers to
a future phase along with the canonical type.

## Domain expansion

Phase 17.7 is a **sister expansion of the Phase 17.5 structure-
prediction beachhead** plus a structure-search adapter on the side
— the same protein-structure-prediction surface broadened with
two more established prediction tools (RoseTTAFold, OmegaFold)
and the canonical structure-search tool (FoldSeek) that the Phase
17.5 ESMFold / OpenFold / AlphaFold 2 / AlphaFold 3 set deferred.
ESMFold (Phase 17.5) is the canonical MSA-free Meta language-model
predictor; OpenFold (Phase 17.5) is the canonical PyTorch AF2
reimplementation; AlphaFold 2 (Phase 17.5) is the DeepMind
reference implementation; AlphaFold 3 (Phase 17.5) is the all-atom
complex predictor; ColabFold (Phase 17) is the de-facto MSA-driven
prediction front-end; RoseTTAFold (Phase 17.7) is the original
Baker-lab 3-track network; OmegaFold (Phase 17.7) is the MSA-free
HelixonAI larger-backbone predictor; FoldSeek (Phase 17.7) is the
Steinegger-lab structure-search analogue of MMseqs2. With Phase
17.7 the protein structure prediction + search surface in Valenx
covers all eight canonical shapes — MSA-driven (ColabFold),
MSA-free language-model (ESMFold), PyTorch AF2 reimplementation
(OpenFold), DeepMind reference (AlphaFold 2 / AlphaFold 3), Baker-
lab 3-track (RoseTTAFold), MSA-free transformer (OmegaFold), and
3Di-alphabet structure search (FoldSeek).

## What landed early

The implementation landed across 5
discrete implementation commits (3 adapters, the registry rollup,
the init-template extension) plus this docs pass — each landing
one adapter, the registry rollup, the init-template extension, or
the documentation pass. Every commit kept workspace clippy +
rustdoc clean.

## Acceptance checklist

- [x] `valenx-adapter-rosettafold` adapter ships with case-input
      parser + 4 lib tests + 3 case-input tests covering parses-
      minimal / parses-with-overrides / rejects-non-py-script,
      plus the Python-script subprocess shape that enforces
      `.py`, stages script + fasta, writes `valenx_params.json`
      with `output_basename` + `fasta` (staged filename), and
      composes `<python> <staged_script>`; probe pushes the
      "model weights + dependencies not bundled — clone
      https://github.com/RosettaCommons/RoseTTAFold and follow
      the install README" warning whenever Python is detected
- [x] `valenx-adapter-omegafold` adapter ships with case-input
      parser + 4 lib tests + 3 case-input tests covering parses-
      minimal / parses-with-overrides / rejects-empty-output-
      basename, plus the single-binary CLI subprocess shape that
      composes `omegafold <fasta> <output_basename> [--model
      <model_dir>]` with the FASTA passed by absolute path (not
      staged); probe via `find_on_path(&["omegafold", "python3",
      "python"])` surfaces a "OmegaFold CLI not found on PATH;
      install via pip install
      git+https://github.com/HeliXonProtein/OmegaFold.git"
      warning when `omegafold` is absent but Python is present
- [x] `valenx-adapter-foldseek` adapter ships with case-input
      parser + 4 lib tests + 3 case-input tests covering parses-
      minimal / parses-with-overrides / rejects-empty-database,
      plus the single-binary subprocess shape that composes
      `foldseek easy-search <query> <database> <basename>.m8
      tmp_<basename> --threads <N> [extras...]` with the database
      parent directory validated on disk
- [x] All 3 adapters wired into `valenx-app::init_registry` —
      live adapter count moves from 109 to **112** (alongside the
      Phase 5.6 MD-engine trio and Phase 5.7 MDTraj single-
      adapter), rounding out the protein structure prediction +
      structure search surface that Phase 17.5 ESMFold / OpenFold
      / AlphaFold 2 / AlphaFold 3 + Phase 17 ColabFold opened
- [x] 3 structure-tool templates in `valenx-init` (`rosettafold`,
      `omegafold`, `foldseek`), all round-tripping through
      `valenx-validate` (cross-binary roundtrip now sweeps **108
      templates** clean alongside the Phase 5.6 MD-engine trio
      and Phase 5.7 MDTraj single template)
- [x] STATUS.md / README.md / ROADMAP.md / CHANGELOG.md /
      QUICKSTART.md updated
- [ ] **Future structure-tool work** — DALI (the original
      structural-alignment tool sister to FoldSeek; defer),
      TM-align (per-pair structural alignment with TM-score
      output; defer), US-align (the universal structural
      alignment successor to TM-align; defer), Phyre2 (web-
      service-only, doesn't fit the local-binary adapter
      pattern), I-TASSER (academic-only with a complex install
      pipeline; defer), boltz-1 (the open boltz reimplementation
      of an AlphaFold-3-style all-atom predictor; defer to a
      future structure-prediction expansion), Chai-1 (the chai
      lab's open AlphaFold-3 successor; defer). Out of scope for
      this expansion.

## Success metrics

| Metric                                                | Target          |
|-------------------------------------------------------|-----------------|
| New structure-tool adapter (template + tests)         | 1 day per       |
| Protein structure prediction + structure search surface across 8 tools (Phase 17 ColabFold + Phase 17.5 ESMFold / OpenFold / AlphaFold 2 / AlphaFold 3 + Phase 17.7 RoseTTAFold / OmegaFold / FoldSeek) | < tool baseline |

## Leads into

Phase 17.7 rounds out the protein structure prediction + structure
search surface alongside the Phase 17.5 ESMFold / OpenFold /
AlphaFold 2 / AlphaFold 3 and Phase 17 ColabFold beachheads.
Combined with the existing predict-structure → fold-RNA → analyze-
DNA-geometry → infer-tree-ML → infer-tree-Bayesian → simulate-
popgen → analyze-trees → simulate-pathway → reconstruct-3D →
design-protein → validate loop, the **predict-structure-
ESMFold → predict-structure-OpenFold → predict-structure-AF2 →
predict-structure-AF3 → predict-structure-RoseTTAFold → predict-
structure-OmegaFold → search-structure-FoldSeek → fold-RNA →
analyze-DNA-geometry → infer-tree-ML → infer-tree-Bayesian →
simulate-popgen → analyze-trees → simulate-pathway → reconstruct-
3D → design-protein → validate** loop now spans seven structure
predictors (Phase 17 ColabFold + Phase 17.5 ESMFold / OpenFold /
AlphaFold 2 / AlphaFold 3 + Phase 17.7 RoseTTAFold / OmegaFold)
and one structure-search tool (Phase 17.7 FoldSeek) feeding into
the Phase 27 / 27.5 / 27.6 protein-design stack, the Phase 28 RNA-
structure tools (ViennaRNA, RNAstructure, NUPACK), the Phase 29
population-genetics trio (SLiM, msprime, tskit), the Phase 30
phylogenetic-tree builders (IQ-TREE, RAxML-NG, FastTree), the
Phase 30.5 Bayesian-phylogenetics pair (BEAST 2, MrBayes), the
Phase 32 systems-biology surface (COPASI, BioNetGen, PhysiCell),
the Phase 33 synthetic-biology trio (pySBOL, j5, Cello), the
Phase 34 docking pair (AutoDock Vina, AutoDock 4), the Phase 35
CRISPR-design tools (CHOPCHOP, CRISPOR, Cas-OFFinder), the Phase
36 cryo-EM reconstruction tools (RELION, EMAN2, CTFFIND), the
Phase 38 Rosetta-family adapters (Rosetta, PyRosetta), and the
Phase 39 DNA-structural-geometry tools (X3DNA, Curves+, DSSR) —
all in one Valenx shell with no glue code beyond the existing
case-toml / prepare / run / collect path.

The natural follow-up is **Phase 17.8** — the deferred structure-
tool work called out above (DALI / TM-align / US-align as
structural-alignment sister tools to FoldSeek, boltz-1 and Chai-1
as open AlphaFold-3-style all-atom predictors), slotting in
alongside the existing RoseTTAFold / OmegaFold / FoldSeek adapters
with the same Python-script subprocess shape (RoseTTAFold sister
tools) or single-binary subprocess shape (FoldSeek / OmegaFold
sister tools). Phase 17.6 (confidence-aware structure ranking +
mmCIF reader) remains a separate deferred docs-cleanup phase.
