# Phase 5.7 — MDTraj

**Status:** 🟢 Live — MDTraj rounds out the **MD trajectory
analysis surface** as a single-adapter sister to Phase 17
MDAnalysis, alongside the Phase 5.5 PLUMED + ProDy + cpptraj
analysis trio.

## Goal

Sister-adapter expansion of the existing Phase 17 MDAnalysis
adapter and the Phase 5.5 PLUMED / ProDy / cpptraj analysis trio.
Round out the post-MD analysis surface with the second-most-used
Python MD trajectory analyzer — **MDTraj** (Pande / VanderSpoel /
Beauchamp lab, LGPL-2.1), which has wider format support than
MDAnalysis (`.xtc` / `.dcd` / `.h5` / `.nc` / `.trr` / `.binpos` /
`.lh5` / `.amber` / `.gromacs` / etc.) and deeper integration with
the OpenMM ecosystem (the Pande / Beauchamp lab is co-located with
the OpenMM developers — MDTraj's HDF5 trajectory format is
OpenMM's native streaming output). Phase 5.7 follows the
established Phase 17 Biopython Python-script subprocess shape —
the user supplies a Python script that imports the upstream
package and reads `valenx_params.json` for the parsed knobs. **This
is a single-adapter phase** — there are precedents for one-adapter
phases when an established tool fills a clearly-defined corner of
an existing surface (no new infrastructure, no new canonical
types, no new CLIs — just one more Python-script adapter with the
trajectory + topology pair pinned at prepare time so the script can
resolve them via relative paths). Phase 5.7 sits numerically
adjacent to Phase 5.5 + Phase 5.6 and ships chronologically right
after Phase 5.6 — same chronological-vs-numerical convention used
for Phase 17.5 / 24 / 28 / 31 / 35 / 39 / 5.5 / 5.6.

## Capability inventory

### Live adapters (1)

- **MDTraj** — Pande / VanderSpoel / Beauchamp lab's Python MD
  trajectory analysis library (LGPL-2.1). MDTraj is the second-
  most-used Python MD trajectory analyzer alongside MDAnalysis —
  wider format support (`.xtc` / `.dcd` / `.h5` / `.nc` / `.trr`
  / `.binpos` / `.lh5` / `.amber` / `.gromacs` / etc.), deeper
  integration with the OpenMM ecosystem (the Pande / Beauchamp
  lab is co-located with the OpenMM developers — MDTraj's HDF5
  trajectory format is OpenMM's native streaming output), and a
  pandas-friendly per-frame property API. Python-script
  subprocess shape (sister to Phase 17 Biopython, Phase 5.5
  ProDy, Phase 17 OpenMM): the user supplies a Python script
  referenced from `[bio.mdtraj].script` in `case.toml` that
  imports `mdtraj` and reads `valenx_params.json` for the parsed
  knobs. Schema knobs: `script` (path to user-supplied Python
  script; required, `.py` enforced), `python` (interpreter name;
  default `"python3"`), `trajectory` (`.xtc` / `.dcd` / `.h5` /
  `.nc` / `.trr` / `.binpos` / `.lh5` MDTraj-supported
  trajectory; required), `topology` (`.pdb` / `.prmtop` / `.gro`
  / `.psf` topology MDTraj uses for atom + residue + chain
  metadata; required), `output_basename` (filename stem the
  user's script uses for outputs — surfaced here so collect()
  can label artefacts uniformly even though the Python script
  chooses its own output paths; required, non-empty).
  `prepare()` enforces the `.py` extension on the script,
  resolves all three input paths against the case directory when
  relative, stages script + trajectory + topology into the
  workdir under their original filenames so the script can
  resolve them via relative paths, then writes a flat
  `valenx_params.json` containing `output_basename`, the bare
  `trajectory` filename, and the bare `topology` filename, and
  builds `<python> <staged_script>`. `collect()` walks the
  workdir for `<output_basename>*.csv` (`Tabular`, "MDTraj
  analysis table" — the canonical pandas-friendly per-frame
  output every downstream MDTraj pipeline reads),
  `<output_basename>*.npz` (`Native`, "MDTraj numpy archive" —
  the per-frame arrays MDTraj writes when downstream consumers
  prefer NumPy over pandas), `<output_basename>*.h5` (`Native`,
  "MDTraj HDF5 output" — the OpenMM-native HDF5 streaming format
  MDTraj writes for re-emission of processed trajectories),
  `<output_basename>*.png` (`Native`, "MDTraj plot" —
  visualisation outputs from matplotlib / seaborn calls inside
  the user's analysis script), and `*.log` (`Log`). Probe via
  `find_on_path(&["python3", "python"])` then `<python> -c
  "import mdtraj"` — on import failure surface as a
  `ProbeReport.warnings` entry (not error) so non-standard
  installs aren't blocked (sister to the Phase 19.5 scanpy /
  scvi / Phase 19.6 AnnData / Phase 5.6 HOOMD-blue probe
  convention). Version range `1.9.0..2.0.0` (MDTraj 1.9 (2022)
  is the modern stable line that pairs with OpenMM 8.x; upper
  bound 2.0 reserves room for the eventual 2.0 stabilisation).
  `bio.mdtraj.analyze` ribbon capability.

### Canonical types

**No new canonical types.** MDTraj consumes user-supplied inputs
(MDTraj Python analysis scripts + trajectories + topologies) and
emits user-readable artifacts (MDTraj `.csv` per-frame analysis
tables + `.npz` numpy archives + `.h5` HDF5 trajectories + `.png`
plots) that the unchanged `Results.artifacts` collection model
surfaces directly. The existing `valenx_bio::format::pdb` reader
already inspects collected PDB topologies; the existing
`valenx_bio::format::dcd` reader already inspects collected DCD
trajectories. A first-class MD-analysis canonical type — a typed
collective-variable / normal-mode / per-frame-statistics
representation spanning all five back-ends (MDAnalysis, MDTraj,
PLUMED, ProDy, cpptraj) — defers to a future phase along with
COLVAR plotters, normal-mode visualizers, and per-statistic time-
series viewers (deferred from Phase 5.5 to the same future phase).

### Headless CLIs

**No new CLIs.** MDTraj's `.csv` per-frame tables are inspectable
in any editor or through the user's downstream Python pipeline
(`pandas.read_csv`); `.npz` numpy archives are inspectable through
`numpy.load`; `.h5` HDF5 trajectories are inspectable through the
user's downstream Python pipeline (`mdtraj.load`); `.png` plots
are inspectable through any image viewer. Trajectory topologies
(PDB) are inspectable through the existing Phase 17
`valenx-pdb-info` CLI. A canonical MD-analysis CLI defers to a
future phase along with the canonical type (deferred from
Phase 5.5 to the same future phase).

## Domain expansion

Phase 5.7 is a **single-adapter expansion of the Phase 17
MDAnalysis adapter and the Phase 5.5 PLUMED / ProDy / cpptraj
analysis trio** — the same MD-trajectory analysis surface
broadened with one more established Python tool that covers the
corners MDAnalysis doesn't reach (wider format support, deeper
OpenMM integration). MDAnalysis is the de-facto Python library
for trajectory I/O + per-frame property calculation (Phase 17);
PLUMED is the de-facto plug-in for biased / enhanced-sampling
work and free-energy reweighting (Phase 5.5); ProDy is the
de-facto Python library for elastic-network / normal-mode protein
dynamics (Phase 5.5); cpptraj is the canonical AmberTools
trajectory analysis CLI (Phase 5.5); MDTraj is the second-most-
used Python MD trajectory analyzer with wider format support and
OpenMM-native HDF5 streaming (Phase 5.7). With Phase 5.7 the
post-MD analysis surface in Valenx covers all five canonical
shapes — Python library API (MDAnalysis + ProDy + MDTraj),
enhanced-sampling plug-in CLI (PLUMED), and AmberTools domain-
language CLI (cpptraj). One-adapter phases are a precedent in
Valenx — when an established tool fills a clearly-defined corner
of an existing surface without requiring new infrastructure, the
phase ships as a single adapter (precedent: a future
docs-cleanup phase could combine multiple such single-adapter
expansions into a roll-up note, but each lands cleanly on its own
at adapter-shipped time).

## What landed early

The implementation landed across 3
discrete implementation commits (1 adapter, the registry rollup,
the init-template extension) plus this docs pass — each landing
the adapter, the registry rollup, the init-template extension, or
the documentation pass. Every commit kept workspace clippy +
rustdoc clean.

## Acceptance checklist

- [x] `valenx-adapter-mdtraj` adapter ships with case-input
      parser + 4 lib tests + 3 case-input tests covering parses-
      minimal / parses-with-overrides / rejects-non-py-script,
      plus the Python-script subprocess shape that enforces
      `.py`, resolves trajectory + topology + script against the
      case directory, stages all three under their original
      filenames, writes `valenx_params.json` with
      `output_basename` + bare `trajectory` filename + bare
      `topology` filename, and composes
      `<python> <staged_script>`
- [x] Adapter wired into `valenx-app::init_registry` — live
      adapter count moves from 108 to **109** (alongside the
      Phase 5.6 MD-engine trio and Phase 17.7 structure-tools
      trio that bring the total to **112**), rounding out the
      post-MD analysis surface that Phase 17 MDAnalysis + Phase
      5.5 PLUMED / ProDy / cpptraj opened
- [x] 1 MD-analysis template in `valenx-init` (`mdtraj`), round-
      tripping through `valenx-validate` (cross-binary roundtrip
      now sweeps **108 templates** clean alongside the Phase 5.6
      MD-engine trio and Phase 17.7 structure-tools trio)
- [x] STATUS.md / README.md / ROADMAP.md / CHANGELOG.md /
      QUICKSTART.md updated
- [ ] **Future MD-analysis work** — nMOLDYN (neutron-scattering
      observables from MD trajectories; defer), CHARMM-GUI (web-
      fronted CHARMM input generator; defer), Lomap2 (alchemical
      free-energy GPU engines sister to PLUMED's reweighting;
      defer to a docking / free-energy phase). Out of scope for
      this expansion.

## Success metrics

| Metric                                                | Target          |
|-------------------------------------------------------|-----------------|
| New MD-analysis adapter (template + tests)            | 1 day per       |
| Post-MD analysis surface across 5 tools (Phase 17 MDAnalysis + Phase 5.5 PLUMED / ProDy / cpptraj + Phase 5.7 MDTraj) | < tool baseline |

## Leads into

Phase 5.7 rounds out the post-MD analysis surface alongside the
Phase 17 MDAnalysis and Phase 5.5 PLUMED / ProDy / cpptraj
beachheads. Combined with the existing simulate-MD → analyze-
trajectory → reweight-free-energy → fit-ENM → run-cpptraj-script
loop, the **simulate-MD → analyze-trajectory-MDAnalysis →
analyze-trajectory-MDTraj → reweight-free-energy → fit-ENM →
run-cpptraj-script → predict-structure → fold-RNA → analyze-DNA-
geometry → infer-tree-ML → infer-tree-Bayesian → simulate-popgen
→ analyze-trees → simulate-pathway → reconstruct-3D → design-
protein → validate** loop now spans five MD-analysis tools (the
Phase 17 MDAnalysis adapter plus PLUMED, ProDy, cpptraj from
Phase 5.5 plus MDTraj from Phase 5.7) feeding into the Phase 5
GROMACS / LAMMPS and Phase 5.6 NAMD / sander / HOOMD-blue MD
engines and the entire Phase 17 → 39 biology / biotech /
chemistry expansion — all in one Valenx shell with no glue code
beyond the existing case-toml / prepare / run / collect path.

The natural follow-up is **a docs-cleanup phase or a future MD-
analysis expansion** — the deferred work called out above (nMOLDYN
for neutron-scattering observables, CHARMM-GUI for the web-fronted
CHARMM input generator, Lomap2 for alchemical free-energy GPU
engines), slotting in alongside the existing PLUMED / ProDy /
cpptraj / MDTraj adapters with the same Python-script subprocess
shape (MDTraj sister tools) or single-binary subprocess shape
(PLUMED / cpptraj sister tools).
