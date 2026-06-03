# Phase 34 — Molecular docking

**Status:** 🟢 Live — small-molecule docking beachhead landed.

## Goal

Add the de-facto open-source small-molecule docking pair to
Valenx's biology / chemistry stack: AutoDock Vina (the modern
single-binary docker) and AutoDock 4 (the older two-stage
`autogrid4 → autodock4` workflow that's still widely used in
pharma teaching + tutorials). Both adapters follow the
established Phase 18 BWA shape — single-action subprocess,
file in / file out, no GPU required. AutoDock 4 has the only
twist: a two-stage `prepare()` (autogrid4 writes the grid maps,
then autodock4 reads the grid + ligand) that mirrors BWA's
`bwa index` → `bwa mem` pattern.

## Capability inventory

### Live adapters (2)

- **AutoDock Vina** — modern single-binary small-molecule
  docker. Takes a receptor PDBQT + ligand PDBQT (prepared
  upstream via `prepare_receptor4.py` / Open Babel — adapter
  does NOT manage prep) and writes ranked binding poses to a
  user-named output PDBQT. The search space is defined by a
  centre `[x, y, z]` and a size `[x, y, z]` in Å, both required
  inputs. `exhaustiveness` (default 8, range 1..=32) tunes
  search depth; `num_modes` (default 9) controls the number of
  poses surfaced; `energy_range` (default 3.0 kcal/mol) bounds
  the energy window between best and worst returned poses;
  `cpu` (default 0 = auto-detect) selects thread count. Output
  PDBQT is collected as a `Native` artifact with label
  `"AutoDock Vina docked poses"`. Apache-2.0 licensed.
  `bio.vina.dock` ribbon capability.
- **AutoDock 4** — the older two-stage docker still common in
  teaching + tutorials. Stage 1: `autogrid4 -p <receptor>.gpf
  -l <grid_log>` writes the grid maps. Stage 2: `autodock4
  -p <ligand>.dpf -l <dock_log>` reads the maps + the docking
  parameter file and runs the docking. Adapter mirrors BWA's
  two-stage shape: stage 1 runs synchronously inside
  `prepare()`, stage 2 lands as the `PreparedJob.native_command`
  for the shared subprocess runner. `skip_grid` (default
  `false`) lets users reuse pre-generated grid maps;
  `grid_log` (default `"autogrid4.glg"`) and `dock_log`
  (default `"autodock4.dlg"`) name the per-stage log files
  inside the workdir. Probe surfaces a warning if `autogrid4`
  is missing from PATH while `autodock4` is present (since the
  full workflow needs both binaries unless `skip_grid` is on).
  GPL-2.0-or-later licensed. `bio.autodock4.dock` ribbon
  capability.

### Canonical types

**No new canonical types.** Both adapters consume PDBQT
inputs (receptor + ligand) and write PDBQT / `.dlg` / `.glg`
outputs that the unchanged `Results.artifacts` collection
model surfaces directly. PDBQT is a PDB-extension format the
existing `valenx_bio::format::pdb` reader can already inspect
for atom counts; ranked poses are user-readable as plain text.

### Headless CLIs

**No new CLIs.** Vina's docked-pose PDBQT outputs and
AutoDock 4's `.dlg` / `.pdbqt` outputs are already inspectable
through the Phase 17 `valenx-pdb-info` CLI (PDBQT is a PDB
superset) — the existing tooling covers the docking loop
without further work.

## What landed early

The implementation rode subagent-driven-development across 5
discrete commits, each landing one adapter, the registry
rollup, the init-template extension, or the documentation
pass. Every commit kept workspace clippy + rustdoc clean.

## Acceptance checklist

- [x] `valenx-adapter-vina` adapter ships with case-input
      parser + 4 lib tests + 5 case-input tests
- [x] `valenx-adapter-autodock4` adapter ships with case-input
      parser + 4 lib tests + 4 case-input tests
- [x] Both adapters wired into `valenx-app::init_registry` —
      live adapter count moves from 42 to 44
- [x] 2 docking templates in `valenx-init` (`vina` with
      `autodock-vina` alias, `autodock4` with `ad4` alias),
      both round-tripping through `valenx-validate`
      (cross-binary roundtrip now sweeps 40 templates clean)
- [x] STATUS.md / README.md / ROADMAP.md / CHANGELOG.md /
      QUICKSTART.md updated
- [ ] **Phase 34.5** — sister-adapter expansion: smina, GNINA
      (deep-learning rescoring fork of smina), DiffDock
      (diffusion-based blind docking), HADDOCK (information-
      driven protein-protein + protein-ligand), rDock (open-
      source docking with cavity detection). The AutoDock-GPU
      fork slots in here too. Out of scope for this beachhead;
      next plan covers it.

## Success metrics

| Metric                                        | Target          |
|-----------------------------------------------|-----------------|
| New docking adapter (template + tests)        | 1 day per       |
| Receptor + ligand → ranked poses loop         | < tool baseline |

## Leads into

Phase 34 paired with the Phase 17 RDKit + ChimeraX adapters
gives Valenx the upstream → docking → visualisation chain in
one shell: RDKit prepares the ligand → AutoDock Vina /
AutoDock 4 docks against the receptor → ChimeraX or PyMOL
(Phase 23) renders the ranked poses. Receptor preparation
itself (`prepare_receptor4.py` chains, Open Babel
preprocessing) lands in **Phase 24** — out of scope here.

The natural follow-up is **Phase 34.5** — sister-adapter
expansion: smina, GNINA, DiffDock, HADDOCK, and rDock slot in
alongside Vina + AutoDock 4 with the same single-binary or
two-stage subprocess shape. The AutoDock-GPU fork (different
binary from CPU Vina, GPU-accelerated) slots in here too.
See the future-phases table at the end of
`docs/superpowers/plans/2026-04-30-biology-foundation.md`
for the full follow-up phase list.
