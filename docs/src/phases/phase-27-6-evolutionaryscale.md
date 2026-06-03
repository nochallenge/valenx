# Phase 27.6 — EvolutionaryScale models

**Status:** 🟢 Live — ESM3 + ESM Cambrian (ESMC) close out the
EvolutionaryScale open-source ESM lineup alongside the Phase 17.5
ESMFold and Phase 27.5 ESM-IF beachheads.

## Goal

Complete EvolutionaryScale's open-source ESM lineup in Valenx. Phase
17.5 + 27.5 already shipped **ESMFold** (single-sequence structure
prediction) and **ESM-IF** (GVP-based inverse-folding sequence
design). Phase 27.6 adds the remaining two: **ESM3**
(EvolutionaryScale's flagship generative multi-modal protein model
— joint reasoning over sequence + structure + function tracks, with
modes `design` / `inverse-fold` / `scaffold` / `predict`) and **ESM
Cambrian / ESMC** (the smaller-faster protein representation model
for embedding-driven downstream ML, two open checkpoints
`esmc-300m` / `esmc-600m`). Both follow the established Phase 17.5
ESMFold / Phase 27.5 ESM-IF pattern — Python-script subprocess
where the user's script imports the relevant package (the same
EvolutionaryScale `esm` package both ESMFold and ESM-IF already
pull in) and reads `valenx_params.json` (auto-written by the
adapter) for config knobs. No new infrastructure. Phase 27.6 sits
numerically after Phase 27.5 and ships chronologically right after
Phase 25 quantum chemistry — same convention as Phase 17.5 sits
between Phase 17 and Phase 18 numerically.

## Capability inventory

### Live adapters (2)

- **ESM3** — EvolutionaryScale's flagship generative multi-modal
  protein model (Cambrian-Open-License — open weights for the
  smaller checkpoints, non-commercial for the largest Forge-only
  variants). Python-script subprocess shape: the user's script
  imports `esm` and reads `valenx_params.json` for the config
  knobs. Where ESMFold and ESM-IF each tackle a single direction
  (sequence → structure, structure → sequence), ESM3 reasons
  jointly over sequence, structure, and function tracks and can be
  conditioned on any subset to fill in the rest. Schema knobs:
  `script` (required Python file), `python` (default `python3`),
  `model_variant` ∈ `{open, open-multimer, small}` (the open-weight
  variants — larger Forge-only variants are not supported by this
  adapter), `mode` ∈ `{design, inverse-fold, scaffold, predict}`,
  `num_samples` (default 4, ≥ 1), `input_pdb` (optional PDB —
  required for `inverse-fold` and `scaffold`), `input_fasta`
  (optional FASTA — required for `predict`), `temperature` (default
  1.0, > 0 and finite), `output_basename` (designs land at
  `<output_basename>*.{pdb,fa}`). `prepare()` stages the script +
  optional PDB / FASTA, writes `valenx_params.json` with the staged
  filenames and the knob values, builds
  `native_command = [python, script]`. `collect()` walks for
  `<output_basename>*.pdb` (kind `Native`, label
  `"ESM3 generated structure"`) and `<output_basename>*.fa` (kind
  `Tabular`, label `"ESM3 generated sequence"`). Probe via
  `find_on_path(&["python3", "python"])` then
  `python -c "import esm; print(esm.__version__)"` — surfaces an
  install hint when Python is on PATH but `esm` isn't importable
  (ESMFold / ESM-IF / ESMC convention). `bio.esm3.generate` ribbon
  capability.
- **ESM Cambrian / ESMC** — EvolutionaryScale's open-weight
  protein representation model (Cambrian-Open-License). Same
  Python-script subprocess shape as ESM3; the user's script imports
  `esm` and reads `valenx_params.json` for the config knobs. Where
  ESM3 is generative and ESMFold targets structure prediction, ESMC
  is the workhorse: it produces high-quality per-residue (or
  pooled) embeddings that downstream classifiers and regressors
  consume directly. Schema knobs: `script` (required Python file),
  `python` (default `python3`), `input_fasta` (sequences to embed,
  required), `model_variant` ∈ `{esmc-300m, esmc-600m}` (the two
  open release sizes — 300M fits on a consumer GPU, 600M for
  larger / better representations), `pooling` ∈
  `{per-residue, mean}`, `output_basename` (embeddings land at
  `<output_basename>.{npy,npz,parquet}`). `prepare()` stages the
  script + input FASTA, writes `valenx_params.json` with the staged
  filename and the knob values. `collect()` walks for
  `<output_basename>.npy` / `<output_basename>.npz` /
  `<output_basename>.parquet` (kind `Tabular`, label
  `"ESMC embeddings"`). Probe identical to ESM3 (`import esm`).
  The init alias `esm-cambrian` resolves to the same template as
  the canonical `esmc` name. `bio.esmc.embed` ribbon capability.

### Canonical types

**No new canonical types.** Both adapters consume the existing
Phase 17 PDB / FASTA inputs and emit user-readable artifacts (PDB
backbones, FASTA sequences, NumPy `.npy` / `.npz` and Parquet
embedding tables) that the unchanged `Results.artifacts` collection
model surfaces directly. The PDB ↔ FASTA pairing flows naturally
through the existing `valenx_bio::format` readers — same shape as
Phase 17.5 ESMFold and Phase 27.5 ESM-IF. The embedding tables
ESMC writes are tabular numerical arrays best inspected through
the user's downstream Python pipeline rather than a Valenx CLI.

### Headless CLIs

**No new CLIs.** ESM3 PDB outputs are already inspectable through
the Phase 17 `valenx-pdb-info` CLI; ESM3 FASTA outputs and ESMC's
embedding sidecars (when expressed as Parquet) are inspectable
through the user's downstream Python tooling — the existing
`valenx-fasta` CLI covers the FASTA side. A canonical embeddings
CLI defers to a future phase along with HDF5 / Arrow reader work.

## EvolutionaryScale lineup complete

Phase 27.6 closes out the **open-source EvolutionaryScale ESM
lineup** at 4 of 4 tools:

- **ESMFold** (Phase 17.5) — single-sequence structure prediction.
- **ESM-IF** (Phase 27.5) — GVP-based inverse-folding sequence
  design (alternative to ProteinMPNN).
- **ESM3** (Phase 27.6) — generative multi-modal joint reasoning
  over sequence + structure + function with `design` /
  `inverse-fold` / `scaffold` / `predict` modes.
- **ESM Cambrian / ESMC** (Phase 27.6) — protein representation
  embeddings for downstream ML (`esmc-300m` / `esmc-600m`).

All four ride the same EvolutionaryScale `esm` Python package
under the hood — installing one installs them all — so the user
sees a single Python-side install hint across every probe and a
unified "ESM is importable" gate via the shared
`detect_esm_version` helper.

## What landed early

The implementation rode subagent-driven-development across 4
discrete implementation commits (2 adapters, the registry rollup,
the init-template extension) plus this docs pass — each landing
one adapter, the registry rollup, the init-template extension, or
the documentation pass. Every commit kept workspace clippy +
rustdoc clean.

## Acceptance checklist

- [x] `valenx-adapter-esm3` adapter ships with case-input parser
      + 4 lib tests + 5 case-input tests covering parses-design-
      minimal / parses-inverse-fold-with-pdb / rejects-inverse-
      fold-without-pdb / rejects-predict-without-fasta / rejects-
      unknown-mode plus the `model_variant` ∈
      `{open, open-multimer, small}` and `mode` ∈
      `{design, inverse-fold, scaffold, predict}` whitelists
- [x] `valenx-adapter-esmc` adapter ships with case-input parser
      + 4 lib tests + 4 case-input tests covering parses-minimal /
      parses-with-600m-and-mean-pooling / rejects-unknown-variant /
      rejects-unknown-pooling plus the `model_variant` ∈
      `{esmc-300m, esmc-600m}` and `pooling` ∈
      `{per-residue, mean}` whitelists
- [x] Both adapters wired into `valenx-app::init_registry` —
      live adapter count moves from 70 to **72**, completing the
      open-source EvolutionaryScale ESM lineup at 4 of 4 tools
- [x] 2 EvolutionaryScale templates in `valenx-init` (`esm3`,
      `esmc` with alias `esm-cambrian`), all round-tripping
      through `valenx-validate` (cross-binary roundtrip now sweeps
      **68 templates** clean)
- [x] STATUS.md / README.md / ROADMAP.md / CHANGELOG.md /
      QUICKSTART.md updated
- [ ] **Future phases** — ESM3 commercial Forge API (would need
      an HTTP-API client, not subprocess; out of scope for the
      open-weights adapter), other EvolutionaryScale Biohub items
      (CELL×GENE, CryoET Data Portal — different shape; tracked as
      future work). Out of scope for this lineup-completion phase.

## Success metrics

| Metric                                          | Target          |
|-------------------------------------------------|-----------------|
| New EvolutionaryScale adapter (template + tests) | 1 day per       |
| Generative + embedding loop across 2 tools       | < tool baseline |

## Leads into

Phase 27.6 closes out the open-source EvolutionaryScale ESM lineup
at 4 of 4 tools (ESMFold + ESM-IF + ESM3 + ESMC). Combined with the
Phases 17 + 17.5 prediction stack and the Phase 27 + 27.5 design
stack, the **design → predict → embed → infer → validate** loop now
spans seven design tools (RFdiffusion, Chroma, RFantibody for
backbone / antibody design; ProteinMPNN, ESM-IF for sequence design;
ESM3 in either generative or scaffold-fill mode), five prediction
tools (ColabFold, ESMFold, OpenFold, AlphaFold 2, AlphaFold 3, plus
ESM3 in `predict` mode), and one embedding workhorse (ESMC) — all
in one Valenx shell with no glue code beyond the existing case-toml
/ prepare / run / collect path.

The natural follow-up is broadening past the EvolutionaryScale
lineup itself — the deferred work called out above (ESM3 commercial
Forge API requiring HTTP-API rather than subprocess shape, other
EvolutionaryScale Biohub items like CELL×GENE / CryoET Data Portal
on a different shape). See the out-of-scope section of
`docs/superpowers/plans/2026-04-30-evolutionaryscale.md` for the
full follow-up phase list.
