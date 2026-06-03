# Phase 19.5 — Single-cell genomics

**Status:** 🟢 Live — Scanpy + scVI open the single-cell genomics
domain alongside the Phase 19 variant-calling stack.

## Goal

Open the single-cell genomics domain in Valenx with the two
most-used Python tools: **Scanpy** (the de-facto single-cell
analysis library — clustering, dimensionality reduction, marker
discovery) and **scVI** (probabilistic deep-learning models for
single-cell data via the `scvi-tools` package). Both adapters
follow the established Phase 17 Biopython / RDKit pattern —
Python-script subprocess where the user's script imports
`scanpy` or `scvi` and reads `valenx_params.json` (auto-written
by the adapter) for config knobs. Outputs (typically `.h5ad`
AnnData files + plots + clustering tables) get walked on
`collect()`. Phase 19.5 sits numerically before Phase 22 but
ships chronologically after Phase 22 — same convention as Phase
17.5 / 24.

## Capability inventory

### Live adapters (2)

- **Scanpy** — de-facto Python single-cell analysis library
  (BSD-3-Clause). Python-script subprocess shape: the user's
  script imports `scanpy as sc` and reads `valenx_params.json`
  for the config knobs. Schema knobs: `script` (required
  Python file), `python` (default `python3`), `input_h5ad`
  (input AnnData file or 10x mtx directory; staged into the
  workdir if relative), `output_h5ad` (output AnnData filename
  the script should write), `n_top_genes` (highly variable gene
  count, default 2000, ≥ 1), `n_pcs` (PCA components, default
  50, ≥ 1), `n_neighbors` (k-NN graph neighbors, default 15,
  ≥ 1), `resolution` (Leiden resolution, default 1.0, > 0 and
  finite). `prepare()` stages the script + input `.h5ad`,
  writes `valenx_params.json` with the staged filename and the
  knob values, builds `native_command = [python, script]`.
  `collect()` walks for `.h5ad` (kind `Native`, label
  `"Scanpy AnnData output"`), `.png` / `.pdf` (kind `Native`,
  label `"Scanpy plot"`), `.csv` / `.tsv` (kind `Tabular`,
  label `"Scanpy table"`). Probe via `find_on_path(&["python3",
  "python"])` then `python -c "import scanpy"` — surfaces an
  install hint when Python is on PATH but `scanpy` isn't
  importable. `bio.scanpy.analyse` ribbon capability.
- **scVI** — probabilistic deep-learning models for single-cell
  data via the `scvi-tools` package (BSD-3-Clause). Python-
  script subprocess shape mirroring Scanpy: the user's script
  imports `scvi` and reads `valenx_params.json` for config
  knobs. Schema knobs: `script` (required Python file),
  `python` (default `python3`), `input_h5ad` (input AnnData;
  staged into the workdir if relative), `output_h5ad` (output
  AnnData filename the script should write), `model` (one of
  `scvi` / `scanvi` / `totalvi` / `linear-scvi`; default
  `scvi`), `n_latent` (latent-space dimensionality, default 10,
  ≥ 1), `n_layers` (encoder/decoder layer count, default 2,
  ≥ 1), `max_epochs` (training epoch budget, default 400,
  ≥ 1), `batch_key` (optional categorical column name in
  `adata.obs` used for batch correction). Same stage / probe /
  collect shape as Scanpy. The init alias `scvi-tools` resolves
  to the same template. `bio.scvi.train` ribbon capability.

### Canonical types

**No new canonical types.** Both adapters consume `.h5ad`
(HDF5-backed AnnData) and emit `.h5ad` plus user-readable
plots / tables. The `.h5ad` reader-as-canonical-type story
needs the `hdf5` crate (a non-trivial C-library dep) and is
deferred to Phase 19.6 along with the Seurat R-runtime work.

### Headless CLIs

**No new CLIs.** Single-cell outputs (`.h5ad`, `.png` / `.pdf`
plots, `.csv` clustering tables) flow through the unchanged
`Results.artifacts` collection model directly; no Valenx-side
inspector is wanted at this beachhead.

## What landed early

The implementation rode subagent-driven-development across 5
discrete commits, each landing one adapter, the registry rollup,
the init-template extension, or the documentation pass. Every
commit kept workspace clippy + rustdoc clean.

## Acceptance checklist

- [x] `valenx-adapter-scanpy` adapter ships with case-input
      parser + 4 lib tests + 4 case-input tests
- [x] `valenx-adapter-scvi` adapter ships with case-input
      parser + 4 lib tests + 4 case-input tests
- [x] Both adapters wired into `valenx-app::init_registry` —
      live adapter count moves from 49 to 51
- [x] 2 single-cell templates in `valenx-init` (`scanpy`,
      `scvi` with alias `scvi-tools`), both round-tripping
      through `valenx-validate` (cross-binary roundtrip now
      sweeps 47 templates clean)
- [x] STATUS.md / README.md / ROADMAP.md / CHANGELOG.md /
      QUICKSTART.md updated
- [ ] **Phase 19.6** — sister-tool expansion: Seurat (the
      dominant R single-cell library — needs an R-runtime
      adapter pattern (`Rscript`-based subprocess) that's new
      to Valenx; the runtime infrastructure ships alongside
      the adapter), AnnData reader as a canonical type
      (`.h5ad` is HDF5-backed; the reader needs the `hdf5`
      crate, a non-trivial C-library dep), `scanpy-spatial`
      (niche enough to defer), CellxGene visualization (viewer
      concern; slot into Phase 23.5). Out of scope for this
      beachhead.

## Success metrics

| Metric                                            | Target          |
|---------------------------------------------------|-----------------|
| New single-cell adapter (template + tests)        | 1 day per       |
| Scanpy `pp.normalize_total → pca → leiden` loop   | < tool baseline |

## Leads into

Phase 19.5 opens the single-cell genomics domain in Valenx
alongside the existing Phase 17 / 17.5 / 18 / 19 / 22 / 23 / 24 /
27 / 34 biology stack. Single-cell workflows now drive through
the same case-toml / prepare / run / collect shell as any other
biology adapter — the user's Python script imports `scanpy` or
`scvi`, reads the auto-written `valenx_params.json`, and writes
`.h5ad` + plots / tables that flow through `collect()`.

The natural follow-up is **Phase 19.6** — the deferred single-
cell work called out above: Seurat (needs the new R-runtime
adapter pattern; ships the `Rscript`-based subprocess
infrastructure for the first time, opening the door to other
R-based bioinformatics tools), AnnData reader as a canonical
type (`.h5ad` HDF5-backed format; needs the `hdf5` crate and
its C-library dep). See the future-phases table at the end of
`docs/superpowers/plans/2026-04-30-single-cell-genomics.md`
for the full follow-up phase list.
