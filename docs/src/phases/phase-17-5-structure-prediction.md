# Phase 17.5 — Structure prediction expansion

**Status:** 🟢 Live — structure-prediction beachhead expanded.

## Goal

Round out the Phase 17 protein-structure-prediction beachhead from
the single ColabFold adapter into the full set of open-source
sibling tools that take a FASTA query (or AF3-style JSON job spec)
and produce ranked PDB models with pLDDT scores already living in
the B-factor column. Phase 17.5 ships 4 new first-class adapters
covering the de-facto open-source structure-prediction stack
(ESMFold, OpenFold, AlphaFold 2, AlphaFold 3). No new canonical
types — pLDDT rides the existing `Atom.b_factor` field that the
Phase 17 PDB reader already lifts cleanly.

## Capability inventory

### Canonical types (`valenx-bio`)
No new canonical types. Each adapter writes a standard PDB ATOM
record with pLDDT in the B-factor column; the Phase 17
`valenx_bio::format::pdb` reader already round-trips that shape
through `Structure` / `Atom`.

### Format readers
No new format readers. AF3's JSON job spec is consumed by the
underlying tool — the adapter stages the file as-is rather than
parsing it into a canonical type (the tool owns the schema and
DeepMind iterates on it independently).

### Live adapters (4)
- **ESMFold** — Meta's protein language model for single-sequence
  structure prediction. No MSA, no separate database step. Probes
  via `python -c "import esm"`. Adapter shape mirrors ColabFold.
- **OpenFold** — PyTorch reimplementation of AlphaFold 2 with the
  full preset family (`model_1` through `model_5_multimer_v3`)
  validated at the case-input layer. Optional `use_templates`
  knob; default `num_recycles = 3` matches OpenFold's default.
- **AlphaFold 2** — DeepMind's reference AF2 implementation via
  `run_alphafold.py`. Validates `model_preset ∈ {monomer,
  monomer_ptm, multimer}` and that `max_template_date` matches
  `\d{4}-\d{2}-\d{2}`. The MSA / template database stays
  user-provided (path passed through).
- **AlphaFold 3** — DeepMind's all-atom complex predictor (protein
  + nucleic acid + ligand). Consumes AF3's JSON job-spec format
  rather than a raw FASTA. The probe surfaces a non-commercial
  warning into `ProbeReport.warnings` because AF3's model weights
  are released under CC-BY-NC-4.0 — the adapter's own invocation
  mode stays `LicenseMode::Subprocess` but the per-tool license
  surfaces in the registry.

### Headless CLIs (0)
No new CLIs. Predicted structures are PDB files that
`valenx-pdb-info` (Phase 17) already inspects; pLDDT in the
B-factor column is summarised by that CLI's element / residue
tally without any structure-prediction-specific code path.

## What landed early

The implementation landed across 6
discrete commits (one per adapter, one app-level register pass,
one init-template pass). Every commit kept workspace clippy +
rustdoc clean. The Phase 17.5 work piggybacks on Phase 18 —
no new canonical types, no new CLIs, just adapters + templates
+ docs.

## Acceptance checklist

- [x] 4 structure-prediction adapter crates under
      `crates/valenx-adapters/bio/`
- [x] All 4 adapters wired into `valenx-app::init_registry`
- [x] 4 init templates (`esmfold` / `openfold` / `alphafold2` /
      `alphafold3`) with aliases (`af2` / `af3`), all
      round-tripping through `valenx-validate`
- [x] Cross-binary roundtrip test sweeps 30 templates clean
- [x] AF3 probe pushes a non-commercial-weights warning into
      `ProbeReport.warnings`
- [x] STATUS.md / README.md / ROADMAP.md / CHANGELOG.md /
      QUICKSTART.md updated
- [ ] **Phase 17.6** — confidence-aware structure ranking and
      mmCIF reader (deferred from Phase 17). Out of scope for
      this expansion; the next plan covers it.

## Success metrics

| Metric                                        | Target          |
|-----------------------------------------------|-----------------|
| New structure-prediction adapter (template + tests) | 1 day per |
| AF3 non-commercial warning surfaces in registry UI  | 100% of probes  |
| pLDDT round-trips PDB → `Atom.b_factor` → PDB        | exact          |

## Leads into

Phase 17.6 — confidence-aware structure ranking (rank PDBs by
mean pLDDT) and mmCIF reader (the Phase 17 deferral). Phase 18.5
continues the alignment toolkit with BAM (binary BGZF) reading
plus Bowtie2 / HISAT2 / STAR / MMseqs2 / DIAMOND adapters that
share the BWA / minimap2 shape.
