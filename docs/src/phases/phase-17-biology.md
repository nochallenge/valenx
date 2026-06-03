# Phase 17 — Biology + biotech foundation

**Status:** 🟢 Live — beachhead landed.

## Goal

Bring the biology / biotech tool ecosystem under the same shell
as Valenx's existing physics-domain coverage. Phase 17 ships the
canonical types (Sequence / Structure / Trajectory) plus 7 first-
class adapters covering the most-used workflows.

## Capability inventory

### Canonical types (`valenx-bio`)
- `Sequence` — DNA / RNA / Protein with IUPAC alphabet validation,
  case-insensitive normalisation, serde round-trip.
- `Structure` — atom / residue / chain hierarchy for proteins,
  nucleic acids, small molecules. Round-trips through PDB ATOM
  records.
- `Trajectory` — per-frame atomic coordinates from MD output.
  Validated for cross-frame atom-count consistency.

### Format readers
- FASTA reader + writer (60-char body wrap on output)
- PDB reader (ATOM / HETATM records, v3 column layout)
- DCD reader (NAMD / CHARMM / VMD interchange — 3-atom / 2-frame
  synthesised tests prove the wire format works on OpenMM output)
- mmCIF reader stub (full impl deferred to Phase 17.5)

### Live adapters (7)
- **Biopython** — Python-script subprocess. User script imports
  `Bio` and writes outputs the adapter classifies on collect.
- **RDKit** — same pattern, imports `rdkit`. Optional inline
  SMILES list in case.toml.
- **OpenMM** — Python-native MD. Adapter parses the output PDB
  via `valenx_bio::format::pdb::read` and lists the DCD
  trajectory as a typed artifact.
- **ChimeraX** — `.cxc` command scripts. `--nogui` mode renders
  PNG / .cxs session files for headless visualisation.
- **oxDNA** — coarse-grained DNA / RNA MD on `input.dat`. Stages
  optional explicit topology.
- **MDAnalysis** — Python-script wrapper. collect() parses any
  produced DCD via the canonical reader, surfacing
  "trajectory · N frames · M atoms" labels.
- **ColabFold** — protein structure prediction from FASTA. Walks
  the result/ subdir for predicted PDB models.

### Headless CLIs (3)
- `valenx-fasta` — inspect / validate / extract on FASTA files.
  Subcommands + text/JSON output + stdin via `-`.
- `valenx-pdb-info` — structural summary (chain count / residue
  range / element tally) from a PDB file.
- `valenx-blast` — thin BLAST+ wrapper. Auto-detects alphabet from
  the query (DNA → blastn / Protein → blastp). Tests cover help /
  version / missing-binary / empty-query paths; real BLAST runs
  stay out of CI.

## What landed early

The implementation rode subagent-driven-development across 18
discrete commits, each landing one canonical-type, format reader,
adapter, or CLI. Every commit kept workspace clippy + rustdoc
clean.

## Acceptance checklist

- [x] `valenx-bio` crate at clean baseline (95+ tests)
- [x] 7 bio adapters wired into `valenx-app::init_registry`
- [x] 7 bio templates in `valenx-init`, all round-tripping
      through `valenx-validate`
- [x] STATUS.md / README.md / ROADMAP.md / CHANGELOG.md /
      QUICKSTART.md updated
- [ ] **Phase 17.5** — full mmCIF reader, structured oxDNA
      energy CSV → ScalarRecord parsing, pLDDT extraction in
      ColabFold collect(). Out of scope for this beachhead;
      next plan covers it.

## Success metrics

| Metric                                       | Target          |
|----------------------------------------------|-----------------|
| `valenx-fasta inspect <bundled>` time        | < 100 ms        |
| `valenx-bio` test suite wall time            | < 5 s           |
| New bio adapter (template + tests)           | 1 day per       |

## Leads into

Phase 18 — Sequence editors / cloning / plasmid design (~12
adapters: ApE, SerialCloner, UGENE, pLannotate, pydna, …). See
the future-phases table at the end of `docs/superpowers/plans/2026-04-30-biology-foundation.md`
for the full follow-up phase list (Phases 17.5 → 43 cover the
remaining ~190 tools from the user's spec).
