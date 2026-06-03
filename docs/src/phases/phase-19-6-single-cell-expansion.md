# Phase 19.6 — Single-cell expansion

**Status:** 🟢 Live — Seurat + AnnData round out the **single-
cell genomics surface** that Phase 19.5 Scanpy + scVI opened with
the dominant R-runtime single-cell library (Seurat) and the
canonical Python data-container library (AnnData) every other
scverse tool reads and writes.

## Goal

Sister-adapter expansion of the existing Phase 19.5 single-cell
Scanpy + scVI beachhead. Round out the single-cell genomics
surface with the two most-requested tools that Phase 19.5
explicitly deferred — the dominant R-based single-cell analysis
toolkit (**Seurat**, the Satija lab's reference single-cell
library that drives roughly half of single-cell papers
worldwide) and the canonical Python data-container library
(**AnnData**, the scverse foundation library that scanpy / scvi /
scirpy / squidpy / muon all read and write). Phase 19.6
explicitly extends Phase 19.5 with R + AnnData: Seurat
introduces the **Rscript subprocess pattern** to Valenx (the R
analogue of the Python-script pattern that Phase 17 Biopython /
RDKit / OpenMM, Phase 19.5 Scanpy / scVI, and Phase 33 pySBOL
established), and AnnData reuses the existing Python-script
pattern for the canonical container that ties the Phase 19.5
scanpy / scvi adapters to every downstream scverse tool. Phase
19.6 sits numerically after Phase 19.5 and ships chronologically
right after it — same convention used for Phase 17.5 / 18.5 /
18.6 / 18.7 / 24 / 28 / 31 / 35.

## Capability inventory

### Live adapters (2)

- **Seurat** — Satija lab's R-based single-cell analysis toolkit
  (MIT). Seurat is the dominant single-cell analysis library on
  the R side of the bioinformatics ecosystem — clustering,
  dimensionality reduction, integration, marker discovery,
  spatial transcriptomics, plus the broader Satija lab tooling
  (Azimuth, signac, BPCells). Phase 19.6 introduces the
  **Rscript subprocess pattern** to Valenx — the R analogue of
  the Python-script pattern that Phase 17 Biopython, Phase 19.5
  Scanpy, and Phase 33 pySBOL established. The user supplies an
  `.R` script referenced from `[bio.seurat].script` in
  `case.toml` that loads `library(Seurat)` and reads
  `valenx_params.json` for the parsed knobs via
  `jsonlite::fromJSON`. Schema knobs: `script` (path to user-
  supplied R script; required, must end in `.R`), `rscript`
  (binary name; default `"Rscript"`), `input_data`
  (`Option<PathBuf>` — optional input matrix that the script
  loads; supports `.h5` / `.mtx` / `.rds` so users can drop in
  10x HDF5, sparse Matrix Market, or pre-saved Seurat object
  formats), `output_basename` (filename stem the script uses
  for outputs — surfaced here so collect() can label artifacts
  uniformly; required, non-empty). `prepare()` enforces the
  `.R` extension, stages script + optional input_data into the
  workdir under their original filenames so the script can
  resolve them via relative paths, then writes a flat
  `valenx_params.json` containing `output_basename` and
  `input_data` (staged filename when set; the key is **omitted
  entirely** when `input_data` is `None` rather than emitted as
  `null`, matching the hand-rolled JSON convention the rest of
  the bio adapters use). Builds `native_command = [rscript,
  script]`. `collect()` walks the workdir for
  `<output_basename>*.rds` (`Native`, "Seurat object (RDS)" —
  the canonical R-serialised Seurat object format consumed by
  every downstream Seurat / signac / Azimuth pipeline),
  `<output_basename>*.csv` (`Tabular`, "Seurat output table"),
  `<output_basename>*.png` (`Native`, "Seurat plot"), and
  `*.log` (`Log`). Probe via `find_on_path(&["Rscript"])` —
  the probe surfaces an `ok = true` warning when Rscript is
  missing rather than failing, same shape as the Phase 17
  Biopython probe; the probe deliberately does **not** attempt
  to confirm Seurat itself is installed because that would
  require running R (an expensive multi-second startup at probe
  time that conflicts with the rest of the registry's snappy
  PATH-lookup probes). Version range `4.0.0..6.0.0` (Seurat
  4.x is the modern stable line that every recent paper cites;
  Seurat 5 (2024) is the current major; upper bound 6.0
  reserves room for an eventual major bump).
  `bio.seurat.analyze` ribbon capability.
- **AnnData** — scverse's Python single-cell HDF5 data
  container library (BSD-3-Clause). AnnData is the canonical
  container that ties the entire scverse Python ecosystem
  together — scanpy reads and writes `.h5ad`, scvi reads and
  writes `.h5ad`, scirpy reads and writes `.h5ad`, squidpy
  reads and writes `.h5ad`, muon multi-omics frames `.h5mu`
  files around `.h5ad` per modality. The `.h5ad` HDF5 format
  is the de-facto single-cell interchange format; `AnnData` as
  a standalone adapter lets users preprocess, convert, or
  inspect `.h5ad` files independent of the analysis pipelines
  that consume them. Python-script subprocess shape (sister to
  Phase 19.5 Scanpy / scVI): the user supplies a Python script
  referenced from `[bio.anndata].script` in `case.toml` that
  imports `anndata` and reads `valenx_params.json` for the
  parsed knobs. Schema knobs: `script` (path to user-supplied
  Python script; required, must end in `.py`), `python`
  (interpreter name; default `"python3"`), `input_h5ad`
  (`Option<PathBuf>` — optional input single-cell file the
  script loads; supports `.h5ad` (the canonical AnnData
  format) and `.h5` (10x HDF5)), `output_basename` (filename
  stem the script uses for outputs; required, non-empty).
  `prepare()` enforces the `.py` extension, stages script +
  optional input_h5ad into the workdir under their original
  filenames, then writes a flat `valenx_params.json` with the
  same hand-rolled shape as Seurat — `output_basename` plus
  `input_h5ad` (staged filename when set, key omitted when
  `None`). Builds `native_command = [python, script]`.
  `collect()` walks the workdir for `<output_basename>*.h5ad`
  (`Native`, "AnnData h5ad file" — the canonical AnnData
  HDF5-backed format every scverse tool reads and writes),
  `<output_basename>*.csv` (`Tabular`, "AnnData output
  table"), `<output_basename>*.png` (`Native`, "AnnData
  plot"), and `*.log` (`Log`). Probe via `find_on_path(&
  ["python3", "python"])` then `<python> -c "import anndata"`
  — on import failure surface as a `ProbeReport.warnings`
  entry (not error) so non-standard installs aren't blocked,
  same shape as the Scanpy / scVI / Biopython probes. Version
  range `0.9.0..1.0.0` (AnnData 0.9 is the modern stable line
  that pairs with scanpy 1.9 / scvi-tools 1.x; upper bound 1.0
  reserves room for the eventual 1.0 stabilisation).
  `bio.anndata.process` ribbon capability.

### Canonical types

**No new canonical types.** Both adapters consume user-supplied
inputs (Seurat takes an `.R` script plus optional `.h5` / `.mtx`
/ `.rds` input data; AnnData takes a `.py` script plus optional
`.h5ad` / `.h5` input file) and emit user-readable artifacts
(Seurat `.rds` Seurat objects + `.csv` tables + `.png` plots,
AnnData `.h5ad` containers + `.csv` tables + `.png` plots) that
the unchanged `Results.artifacts` collection model surfaces
directly. The first-class AnnData canonical type — a typed
`.h5ad` reader exposing the obs / var / X / layers / obsm /
varm / uns hierarchy as parsed Rust structs — defers to a
future phase along with the `hdf5` crate dependency it requires
(a non-trivial C-library dep that the Phase 19.5 docs already
flagged as out of scope), Seurat-object inspection beyond the
existing artifact-collection model, and per-cell / per-gene
inspection viewers.

### Headless CLIs

**No new CLIs.** Seurat's `.rds` outputs are inspectable
through the user's downstream R pipeline (`readRDS` plus the
broader Seurat / signac / Azimuth toolkit); AnnData's `.h5ad`
outputs are inspectable through the user's downstream Python
pipeline (`anndata.read_h5ad`, `scanpy.read`, `scvi.data`).
The Phase 19.5 deferral note about a canonical AnnData reader
as a Valenx CLI continues to apply — it slots in alongside the
`.h5ad` canonical-type work in a future phase when the `hdf5`
crate dependency lands. CSV / PNG outputs surface directly
through the existing `Tabular` / `Native` artifact kinds.

## Domain expansion

Phase 19.6 is a **sister-adapter expansion of the Phase 19.5
single-cell Scanpy + scVI beachhead** — the same single-cell
genomics surface broadened with two more established tools that
cover the two corners Phase 19.5 explicitly deferred. Seurat is
the dominant R-based single-cell library that drives roughly
half of single-cell papers worldwide (the other half ride on
the Phase 19.5 scverse Python ecosystem); the Rscript
subprocess pattern that Seurat introduces opens the door to
every other R-based bioinformatics tool (Bioconductor, edgeR,
DESeq2, limma, monocle3, scran, BiocManager-installed
packages) for future phases. AnnData is the canonical Python
container library that ties the entire scverse ecosystem
together — scanpy / scvi / scirpy / squidpy / muon all read
and write `.h5ad` — so the standalone AnnData adapter lets
users preprocess, convert, or inspect `.h5ad` files without
needing to drive a full scanpy / scvi pipeline. With Phase
19.6 the single-cell surface in Valenx covers all four
canonical shapes — Python analysis library (Scanpy + scVI
from Phase 19.5), R analysis library (Seurat), and the
canonical Python data-container library (AnnData) — across
the two language ecosystems the field is split between.

## What landed early

The implementation rode subagent-driven-development across 5
discrete implementation commits (2 adapters, the registry
rollup, the init-template extension, plus a follow-up commit
for Seurat) plus this docs pass — each landing one adapter,
the registry rollup, the init-template extension, or the
documentation pass. Every commit kept workspace clippy +
rustdoc clean. Phase 19.6 introduces the **Rscript subprocess
pattern** to the workspace for the first time alongside the
Seurat adapter — the runtime infrastructure (R-script staging,
`valenx_params.json` shape compatible with `jsonlite::
fromJSON`, `.R` extension validation, `Rscript` binary probe)
ships in the Seurat adapter and is reusable from any future
R-based bioinformatics adapter without further plumbing.

## Acceptance checklist

- [x] `valenx-adapter-seurat` adapter ships with case-input
      parser + 4 lib tests + 3 case-input tests, plus the
      Rscript subprocess shape that enforces the `.R`
      extension, stages script + optional input_data, and
      writes `valenx_params.json` with `output_basename` and
      `input_data` (staged filename when set, key omitted
      when `None`)
- [x] `valenx-adapter-anndata` adapter ships with case-input
      parser + 4 lib tests + 3 case-input tests, plus the
      Python-script subprocess shape that enforces the `.py`
      extension, stages script + optional input_h5ad, and
      writes `valenx_params.json` with the same hand-rolled
      shape as Seurat
- [x] Both adapters wired into `valenx-app::init_registry` —
      live adapter count moves from 103 to **105** (alongside
      the Phase 18.7 alignment-toolkit-expansion trio that
      brings the total to **105**), rounding out the single-
      cell genomics surface that Phase 19.5 Scanpy + scVI
      opened
- [x] 2 single-cell-expansion templates in `valenx-init`
      (`seurat`, `anndata`), both round-tripping through
      `valenx-validate` (cross-binary roundtrip now sweeps
      **101 templates** clean alongside the Phase 18.7
      alignment-toolkit-expansion trio)
- [x] STATUS.md / README.md / ROADMAP.md / CHANGELOG.md /
      QUICKSTART.md updated
- [ ] **Phase 19.7** — sister-adapter expansion of Phase
      19.5 / 19.6: scirpy (single-cell immune-receptor
      analysis sister to Scanpy in the scverse ecosystem;
      defer), squidpy (single-cell spatial-omics analysis
      sister to Scanpy; defer), muon (single-cell multi-omics
      analysis library that wraps multiple AnnData containers
      via `.h5mu` files; defer), dandelion (single-cell BCR /
      TCR repertoire analysis sister to scirpy; defer), the
      first-class `.h5ad` canonical type backed by the `hdf5`
      crate (non-trivial C-library dep deferred from Phase
      19.5; defer to its own infrastructure phase), Bioconductor
      / edgeR / DESeq2 R-based bulk RNA-seq analysis tools
      that ride the Phase 19.6 Rscript subprocess pattern
      (defer). Out of scope for this beachhead.

## Success metrics

| Metric                                                | Target          |
|-------------------------------------------------------|-----------------|
| New single-cell adapter (template + tests)            | 1 day per       |
| Seurat `CreateSeuratObject → NormalizeData → FindVariableFeatures → RunPCA → FindClusters` loop | < tool baseline |
| AnnData `read_h5ad → preprocess → write_h5ad` round-trip | < tool baseline |

## Leads into

Phase 19.6 rounds out the single-cell genomics surface that
the user's bio spec called out alongside the Phase 19.5
Scanpy + scVI beachhead. Combined with the existing predict-
structure → fold-RNA → analyze-DNA-geometry → infer-tree-ML
→ infer-tree-Bayesian → simulate-popgen → analyze-trees →
simulate-pathway → reconstruct-3D → design-protein → validate
loop, the **load-data → cluster → integrate → annotate →
visualise** single-cell loop now spans four single-cell tools
across both language ecosystems (the Phase 19.5 Scanpy / scVI
Python pair plus the Phase 19.6 Seurat R + AnnData Python
container pair) feeding into the existing Phase 17 / 17.5 /
18 / 18.5 / 18.6 / 18.7 / 19 / 20 / 22 / 23 / 24 / 25 / 27 /
27.5 / 27.6 / 28 / 29 / 30 / 30.5 / 31 / 32 / 33 / 34 / 35 /
36 / 38 / 39 biology / biotech / chemistry stack — all in one
Valenx shell with no glue code beyond the existing case-toml
/ prepare / run / collect path. The Rscript subprocess
pattern that Seurat introduces also opens the door to every
other R-based bioinformatics tool for future phases.

The natural follow-up is **Phase 19.7** — the deferred
single-cell work called out above (scirpy / squidpy / muon /
dandelion as the remaining scverse Python tools sister to the
Phase 19.5 / 19.6 single-cell adapters, the first-class
`.h5ad` canonical type backed by the `hdf5` crate, and the
broader Bioconductor / edgeR / DESeq2 R-based bulk RNA-seq
analysis tools that ride the Phase 19.6 Rscript subprocess
pattern), slotting in alongside the existing Scanpy / scVI /
Seurat / AnnData adapters with the same Python-script
subprocess shape (scverse sister tools) or Rscript
subprocess shape (Bioconductor sister tools). See the
out-of-scope section of `docs/superpowers/plans/2026-05-02-
single-cell-expansion.md` for the full follow-up phase list.
