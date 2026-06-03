# Phase 27 — Protein design

**Status:** 🟢 Live — de novo protein-design beachhead landed.

## Goal

Pair the structure-prediction adapters Valenx already ships
(Phases 17 + 17.5: ColabFold, ESMFold, OpenFold, AlphaFold 2/3)
with their de novo design counterparts. Phase 27 ships
RFdiffusion (GPU-driven protein backbone generation) and
ProteinMPNN (sequence design from backbone). Together with
the prediction stack, this gives Valenx the complete
**design → predict → validate** loop in one shell. Both
adapters follow the established Phase 17 ColabFold shape:
Python-script subprocess, FASTA / PDB in, PDB / FASTA out,
user-managed model weights and GPU runtime.

## Capability inventory

### Live adapters (2)

- **RFdiffusion** — GPU-driven protein backbone generation.
  Drives off a user-supplied Python entry script that imports
  `rfdiffusion` and reads `valenx_params.json` (written by the
  adapter into the workdir) for config knobs. Supports four
  modes via the `mode` field — `motif` (motif scaffolding),
  `binder` (binder design against a target context),
  `unconditional` (free generation), and `partial-diffusion`
  (refinement of an input structure). `num_designs` controls
  the sample count (default 8) and `diffusion_steps` the
  schedule depth (default 50, RFdiffusion's recommended value).
  Sampled designs land at `<output_basename>_0.pdb`,
  `<output_basename>_1.pdb`, … and are surfaced as typed
  `Native` artifacts via `valenx_bio::format::pdb::read`
  (RFdiffusion writes pLDDT into the B-factor column too —
  same as the prediction tools). BSD-3-Clause licensed.
  `bio.rfdiffusion.design` ribbon capability.
- **ProteinMPNN** — sequence design from a backbone PDB. Same
  Python-script-subprocess pattern as RFdiffusion; takes a
  backbone PDB and emits FASTA sequences (one per design)
  with per-residue probabilities. Three model variants via
  `model_variant` — `vanilla` (the published model),
  `soluble` (soluble-protein bias), and `ca-only` (Cα-only
  backbone input). `temperature` (default 0.1) controls
  sampling diversity; `num_seq_per_target` (default 8)
  controls the per-design sample count. Output FASTA lands
  at `<output_basename>.fa` and is parsed via
  `valenx_bio::format::fasta::read_str` for a richer
  `"ProteinMPNN · N sequences"` artifact label.
  MIT licensed. `bio.proteinmpnn.design` ribbon capability.

### Canonical types

**No new canonical types.** Both adapters consume the existing
Phase 17 PDB inputs and emit user-readable artifacts (PDB
backbones, FASTA sequences) that the unchanged
`Results.artifacts` collection model surfaces directly.
The PDB ↔ FASTA pairing flows naturally through the existing
`valenx_bio::format` readers.

### Headless CLIs

**No new CLIs.** RFdiffusion's PDB outputs and ProteinMPNN's
FASTA outputs are already inspectable through the Phase 17
`valenx-pdb-info` and `valenx-fasta` CLIs respectively — the
existing tooling covers the design loop without further work.

## What landed early

The implementation rode subagent-driven-development across 5
discrete commits, each landing one adapter, the registry
rollup, the init-template extension, or the documentation
pass. Every commit kept workspace clippy + rustdoc clean.

## Acceptance checklist

- [x] `valenx-adapter-rfdiffusion` adapter ships with case-input
      parser + 4 lib tests + 4 case-input tests
- [x] `valenx-adapter-proteinmpnn` adapter ships with case-input
      parser + 4 lib tests + 4 case-input tests
- [x] Both adapters wired into `valenx-app::init_registry` —
      live adapter count moves from 40 to 42
- [x] 2 design templates in `valenx-init` (`rfdiffusion` with
      `rfd` alias, `proteinmpnn` with `mpnn` alias), both
      round-tripping through `valenx-validate` (cross-binary
      roundtrip now sweeps 38 templates clean)
- [x] STATUS.md / README.md / ROADMAP.md / CHANGELOG.md /
      QUICKSTART.md updated
- [ ] **Phase 27.5** — sister-adapter expansion: Chroma,
      framediff, Genie (alternative diffusion-based design
      models), RFantibody (antibody-specialised RFdiffusion),
      ESM-IF (inverse folding from ESM). Out of scope for
      this beachhead; next plan covers it.

## Success metrics

| Metric                                        | Target          |
|-----------------------------------------------|-----------------|
| New design adapter (template + tests)         | 1 day per       |
| Backbone generation / sequence design loop    | < tool baseline |

## Leads into

Phase 27 paired with Phases 17 + 17.5 closes the
**design → predict → validate** loop: RFdiffusion generates
backbones → ProteinMPNN designs sequences → ColabFold /
ESMFold / OpenFold / AlphaFold 2 / AlphaFold 3 fold the
sequences back to validate the design. The complete cycle
runs in one Valenx shell with no glue code beyond the
existing case-toml / prepare / run / collect path.

The natural follow-up is **Phase 27.5** — sister-adapter
expansion: Chroma, framediff, Genie, RFantibody, and ESM-IF
slot in alongside RFdiffusion + ProteinMPNN with the same
Python-script-subprocess shape. See the future-phases table at
the end of `docs/superpowers/plans/2026-04-30-biology-foundation.md`
for the full follow-up phase list (Phases 19.5 → 43 cover the
remaining ~190 tools from the user's spec).
