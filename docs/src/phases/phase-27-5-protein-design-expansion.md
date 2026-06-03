# Phase 27.5 — Protein design expansion

**Status:** 🟢 Live — Chroma + ESM-IF + RFantibody round out the
de novo design surface alongside the Phase 27 RFdiffusion +
ProteinMPNN beachhead.

## Goal

Sister-adapter expansion of Phase 27. Add three more open-source
protein design tools to round out the de novo design surface:
**Chroma** (Generate Biomedicines' diffusion model — backbone +
sequence joint), **ESM-IF** (Meta's inverse-folding sequence
designer — alternative to ProteinMPNN), and **RFantibody**
(RosettaCommons antibody-specific RFdiffusion fork). All three
follow the established Phase 27 RFdiffusion / ProteinMPNN pattern
— Python-script subprocess where the user's script imports the
relevant package and reads `valenx_params.json` (auto-written by
the adapter) for config knobs. No new infrastructure. Phase 27.5
sits numerically after Phase 27 and ships chronologically right
after Phase 27 — same convention as Phase 17.5 sits between
Phase 17 and Phase 18 numerically.

## Capability inventory

### Live adapters (3)

- **Chroma** — Generate Biomedicines' joint backbone-and-sequence
  diffusion model (Apache-2.0). Python-script subprocess shape:
  the user's script imports `chroma` and reads
  `valenx_params.json` for the config knobs. Schema knobs:
  `script` (required Python file), `python` (default `python3`),
  `num_samples` (sample count, default 4, ≥ 1), `length`
  (residue count to design, ≥ 1), `temperature` (sampling
  temperature for the diffusion, default 1.0, > 0 and finite),
  `output_basename` (designs land at
  `<output_basename>_N.{pdb,fa}`). `prepare()` stages the script,
  writes `valenx_params.json` with the knob values, builds
  `native_command = [python, script]`. `collect()` walks for
  `<output_basename>*.pdb` (kind `Native`, label
  `"Chroma design"`) and `<output_basename>*.fa` (kind
  `Tabular`, label `"Chroma sequence"`). Probe via
  `find_on_path(&["python3", "python"])` then
  `python -c "import chroma"` — surfaces an install hint when
  Python is on PATH but `chroma` isn't importable.
  `bio.chroma.design` ribbon capability.
- **ESM-IF** — Meta's GVP-based inverse-folding sequence designer
  via the `fair-esm` package (MIT). Same Python-script subprocess
  shape; the user's script imports `esm` (the same package as
  ESMFold) and reads `valenx_params.json` for config knobs.
  Schema knobs: `script` (required Python file), `python`
  (default `python3`), `input_pdb` (input PDB; staged into the
  workdir if relative), `model` (default
  `esm_if1_gvp4_t16_142M_UR50` — non-empty; not whitelisted
  because ESM-IF model identifiers evolve fast and the upstream
  package validates), `temperature` (default 1.0, > 0 and
  finite), `num_samples` (default 8, ≥ 1), `output_basename`.
  `prepare()` stages the script + input PDB, writes
  `valenx_params.json` with the staged filename and the knob
  values. `collect()` walks for `<output_basename>.fa` — parsed
  via `valenx_bio::format::fasta::read` for a richer
  `"ESM-IF · N sequences"` label with the actual count, falling
  back to `"ESM-IF designed sequences"` on parse failure
  (ProteinMPNN pattern). The init aliases `esmif` and
  `inverse-folding` resolve to the same template.
  `bio.esm-if.design` ribbon capability.
- **RFantibody** — RosettaCommons antibody-specific fork of
  RFdiffusion (BSD-3-Clause). Same Python-script subprocess
  shape; the user's script imports `rfantibody` and reads
  `valenx_params.json` for config knobs. Adds antibody-aware
  modes and a CDR-loop-focused sampling protocol. Schema knobs:
  `script` (required Python file), `python` (default `python3`),
  `framework_pdb` (antibody framework PDB; staged into the
  workdir if relative), `target_pdb` (target antigen PDB; staged
  the same way), `design_loops` (CDR loops to design, non-empty
  subset of `["H1", "H2", "H3", "L1", "L2", "L3"]`),
  `num_designs` (default 8, ≥ 1), `diffusion_steps` (default 50,
  ≥ 1), `output_basename`. `prepare()` stages the script + both
  PDBs, writes `valenx_params.json` with the staged filenames and
  the knob values. `collect()` walks for `<output_basename>*.pdb`
  (kind `Native`, label `"RFantibody design"`). The init alias
  `rfab` resolves to the same template.
  `bio.rfantibody.design` ribbon capability.

### Canonical types

**No new canonical types.** All three adapters consume the
existing Phase 17 PDB inputs and emit user-readable artifacts
(PDB backbones, FASTA sequences) that the unchanged
`Results.artifacts` collection model surfaces directly. The
PDB ↔ FASTA pairing flows naturally through the existing
`valenx_bio::format` readers — same shape as Phase 27.

### Headless CLIs

**No new CLIs.** Chroma + RFantibody PDB outputs and ESM-IF FASTA
outputs are already inspectable through the Phase 17
`valenx-pdb-info` and `valenx-fasta` CLIs respectively — the
existing tooling covers the expanded design surface without
further work.

## What landed early

The implementation landed across 6
discrete commits, each landing one adapter, the registry rollup,
the init-template extension, or the documentation pass. Every
commit kept workspace clippy + rustdoc clean.

## Acceptance checklist

- [x] `valenx-adapter-chroma` adapter ships with case-input
      parser + 4 lib tests + 4 case-input tests
- [x] `valenx-adapter-esm-if` adapter ships with case-input
      parser + 4 lib tests + 4 case-input tests
- [x] `valenx-adapter-rfantibody` adapter ships with case-input
      parser + 4 lib tests + 4 case-input tests
- [x] All 3 adapters wired into `valenx-app::init_registry` —
      live adapter count moves from 51 to 54
- [x] 3 design templates in `valenx-init` (`chroma`, `esm-if`
      with aliases `esmif` / `inverse-folding`, `rfantibody`
      with alias `rfab`), all round-tripping through
      `valenx-validate` (cross-binary roundtrip now sweeps 50
      templates clean)
- [x] STATUS.md / README.md / ROADMAP.md / CHANGELOG.md /
      QUICKSTART.md updated
- [ ] **Phase 27.6** — further design-tool expansion: framediff
      and Genie (alternative diffusion-based design models —
      niche enough to defer until user demand surfaces),
      AlphaFold-Multimer-Design (different shape, would need
      direct AlphaFold integration), Hallucination /
      TrDesign-style design (different shape, separate phase).
      Out of scope for this expansion.

## Success metrics

| Metric                                            | Target          |
|---------------------------------------------------|-----------------|
| New design adapter (template + tests)             | 1 day per       |
| Backbone / sequence design loop across 3 tools    | < tool baseline |

## Leads into

Phase 27.5 broadens the de novo design surface alongside the
Phase 27 RFdiffusion + ProteinMPNN beachhead. Combined with the
Phases 17 + 17.5 prediction stack, the **design → predict →
validate** loop now spans five design tools (RFdiffusion, Chroma,
RFantibody for backbone / antibody design; ProteinMPNN, ESM-IF
for sequence design) feeding into five prediction tools
(ColabFold, ESMFold, OpenFold, AlphaFold 2, AlphaFold 3) — all in
one Valenx shell with no glue code beyond the existing case-toml
/ prepare / run / collect path.

The natural follow-up is **Phase 27.6** — the deferred design
work called out above (framediff, Genie, AlphaFold-Multimer-
Design, Hallucination / TrDesign), slotting in alongside the
existing design adapters with the same Python-script-subprocess
shape.
