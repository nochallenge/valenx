# Phase 5.6 — Bio MD engines

**Status:** 🟢 Live — NAMD + AmberTools sander + HOOMD-blue round
out the **all-atom + GPU-native MD engine surface** that Phase 5
GROMACS / LAMMPS opened, alongside the Phase 17 OpenMM Python-
native engine and the Phase 5.5 / 17 PLUMED / ProDy / cpptraj /
MDAnalysis post-MD analysis stack.

## Goal

Sister-domain expansion of the existing Phase 5 GROMACS / LAMMPS
MD engine beachhead. Round out the all-atom MD-engine surface with
three more established open-source engines that span the corners
GROMACS / LAMMPS / OpenMM don't reach — the canonical UIUC
academic NAMD all-atom engine (NAMD, the de-facto choice in
biomolecular MD pedagogy and a workhorse on every academic HPC
cluster), the AmberTools OSS portion of AMBER's MD engine
(`sander`, the OSS sibling of cpptraj that Phase 5.5 already wraps
— sister tool from the same release), and the Glotzer-lab GPU-
native particle simulator (HOOMD-blue, the canonical Python-
scripted GPU-first engine for soft-matter / coarse-grained
particle systems). NAMD + sander follow the established Phase 18
BWA single-binary CLI pattern: configuration file in, trajectory
+ logs out. HOOMD-blue follows the established Phase 17 OpenMM
Python-script subprocess shape: the user supplies a Python script
that imports the upstream package and reads `valenx_params.json`
for the `output_basename` knob the adapter writes for every run.
NAMD is **academic-license-flagged** (sister to ChimeraX, VMD,
ViennaRNA, NUPACK, CTFFIND, Rosetta, X3DNA, Curves+, DSSR — the
established academic-only-tool convention surfaces a `tool_license
= "NAMD-License"` and pushes an `"academic"`-keyworded warning
into `ProbeReport.warnings` whenever the binary is detected). Phase
5.6 sits numerically adjacent to Phase 5.5 and ships
chronologically right after Phase 17.7 structure tools — same
chronological-vs-numerical convention used for Phase 17.5 / 24 /
28 / 31 / 35 / 39 / 5.5.

## Capability inventory

### Live adapters (3)

- **NAMD** — UIUC's flagship all-atom MD engine (custom NAMD-
  License — academic / non-commercial use only). NAMD is the
  de-facto choice in biomolecular MD pedagogy and a workhorse on
  every academic HPC cluster — NAMD 2.x ships an SMP-threaded
  CHARMM-style integrator, NAMD 3.x adds GPU-resident kernels.
  Single-binary subprocess shape (sister to Phase 5 LAMMPS /
  GROMACS): the CLI is `<binary> +p<processors> <config>
  [extras...]` where `<binary>` is `namd2` or `namd3` (the probe
  accepts either) and `+p<N>` (with no space — NAMD's own flag
  syntax) is NAMD's threading flag (multi-threaded SMP build uses
  it for thread count, MPI build uses MPI-rank count from the
  launcher). Schema knobs: `config` (path to NAMD `.namd` /
  `.conf` configuration file; required), `processors` (`u32`,
  default 1; emitted as the single OsString `+p<N>` so the flag
  and value travel together exactly as NAMD parses them),
  `extra_args` (additional CLI arguments appended after the
  config). `prepare()` resolves `config` against the case
  directory when relative and composes the `<binary> +p<N>
  <config> [extras...]` invocation. `collect()` walks the workdir
  for `*.dcd` (`Native`, "NAMD trajectory (DCD)" — the canonical
  NAMD / CHARMM / VMD interchange trajectory format), `*.coor`
  (`Native`, "NAMD coordinates" — restart binary), `*.vel`
  (`Native`, "NAMD velocities" — restart binary), `*.xsc`
  (`Tabular`, "NAMD extended system" — periodic-cell parameters),
  and `*.log` (`Log`). Probe via `find_on_path(&["namd2",
  "namd3"])` — accepting either NAMD 2.14 or NAMD 3.x; pushes an
  `"academic"`-keyworded warning containing both `"academic"` and
  `"non-commercial"` substrings into `ProbeReport.warnings`
  whenever the binary is detected, and `tool_license` surfaces as
  `"NAMD-License"` rather than mislabeling the custom UIUC NAMD
  terms as a recognised SPDX identifier. Version range
  `2.14.0..4.0.0` (NAMD 2.14 (2020) is the long-stable line still
  in production use; NAMD 3.x (2022+) is the GPU-resident
  rewrite; upper bound 4.0 reserves room for a future major).
  `bio.namd.simulate` ribbon capability.
- **AmberTools sander** — AMBER's OSS MD engine portion of
  AmberTools (GPL-3.0 — sander itself is OSS; the proprietary
  `pmemd.cuda` GPU engine is NOT wrapped here and would require
  the per-site AMBER license). sander reads an Amber `.prmtop` /
  `.parm7` topology, an `.inpcrd` / `.rst7` coordinate file, and
  a `.in` / `.mdin` simulation control file, runs the integrator
  for the requested number of steps, and emits an `.out` mdout
  log + `.rst` restart + `.nc` NetCDF trajectory. Sister to
  Phase 5.5 cpptraj (also AmberTools — installing AmberTools
  installs both). Single-binary subprocess shape (sister to
  Phase 18 BWA): the CLI is `sander -O -i <config> -p <topology>
  -c <coordinates> -o <basename>.out -r <basename>.rst -x
  <basename>.nc [extras...]`. The `-O` flag overwrites existing
  outputs (the standard re-run convention). Schema knobs:
  `topology` (`.prmtop` / `.parm7` Amber topology; required),
  `coordinates` (`.inpcrd` / `.rst7` Amber coordinate file;
  required), `config` (`.in` / `.mdin` Amber simulation control
  file; required), `output_basename` (filename stem the adapter
  pins for the three sander output flags so collect() walks
  deterministically; required, non-empty), `extra_args`.
  `prepare()` resolves all three input paths against the case
  directory when relative, validates each file exists on disk
  (returns `InvalidCase` with a helpful message when missing),
  and composes the invocation with the three output paths pinned
  to `<output_basename>.{out,rst,nc}`. `collect()` walks the
  workdir for `<output_basename>*.out` (`Log`, "sander mdout" —
  the AmberTools mdout log every downstream Amber tool reads),
  `<output_basename>*.nc` (`Native`, "sander NetCDF trajectory"
  — the AmberTools NetCDF trajectory format cpptraj already
  consumes), `<output_basename>*.rst` (`Native`, "sander restart
  coordinates" — the per-step restart file), and
  `<output_basename>*.mdinfo` (`Log`, "sander mdinfo" — the
  per-step performance / progress sidecar). Probe via
  `find_on_path(&["sander"])` — no academic-license caveat
  (sander itself is GPL-3.0 OSS; the proprietary `pmemd.cuda`
  variant is not part of this adapter's surface). Version range
  `22.0.0..26.0.0` (AmberTools 22 (2022) is the floor we test
  against; upper bound 26.0 reserves room for the next several
  AmberTools releases). `bio.sander.simulate` ribbon capability.
- **HOOMD-blue** — Glotzer lab's GPU-native particle simulator
  (BSD-3-Clause). HOOMD-blue v3+ is fully Python-scripted (no
  native CLI) — the user supplies a `.py` script that does
  `import hoomd` and runs the simulation; HOOMD's GPU-resident
  kernels handle the per-step force evaluation transparently
  underneath. HOOMD-blue is the canonical engine for soft-matter
  / coarse-grained particle systems, polymers, colloids, and
  rigid-body assemblies — sister to LAMMPS in the particle-MD
  surface but GPU-first by design. Python-script subprocess
  shape (sister to Phase 17 OpenMM): the user supplies a Python
  script referenced from `[bio.hoomd].script` in `case.toml`
  that imports `hoomd` and reads `valenx_params.json` for the
  `output_basename` knob the adapter writes for every run.
  Schema knobs: `script` (path to user-supplied Python script;
  required, `.py` enforced), `python` (interpreter name; default
  `"python3"`), `output_basename` (filename stem the user's
  script uses for outputs — surfaced here so collect() can label
  artefacts uniformly even though the Python script chooses its
  own output paths; required, non-empty). `prepare()` enforces
  the `.py` extension, stages the script into the workdir under
  its original filename so the script can resolve it via
  relative paths, then writes a flat `valenx_params.json`
  containing `output_basename`, and builds `<python>
  <staged_script>`. `collect()` walks the workdir for
  `<output_basename>*.gsd` (`Native`, "HOOMD trajectory (GSD)" —
  the canonical HOOMD trajectory format consumed by every
  downstream HOOMD / freud / ovito pipeline),
  `<output_basename>*.h5` (`Native`, "HOOMD HDF5 output" —
  HOOMD's tabular HDF5 sidecar), and `*.log` (`Log`). Probe via
  `find_on_path(&["python3", "python"])` then `<python> -c
  "import hoomd"` — on import failure surface as a
  `ProbeReport.warnings` entry (not error) so non-standard
  installs aren't blocked (sister to the Phase 19.5 scanpy /
  scvi / Phase 19.6 AnnData probe convention). Version range
  `3.0.0..6.0.0` (HOOMD-blue 3.x (2022) is the modern Python-
  first rewrite; HOOMD-blue 4.x / 5.x are the rolling current
  releases; upper bound 6.0 reserves room for a future major).
  `bio.hoomd.simulate` ribbon capability.

### Canonical types

**No new canonical types.** All three adapters consume user-
supplied inputs (NAMD configuration files, sander topology +
coordinates + simulation control files, HOOMD-blue Python
scripts) and emit user-readable artifacts (NAMD `.dcd`
trajectories + `.coor` / `.vel` restarts + `.xsc` extended-system
tables + `.log` logs, sander `.out` mdout logs + `.nc` NetCDF
trajectories + `.rst` restarts + `.mdinfo` sidecars, HOOMD-blue
`.gsd` trajectories + `.h5` HDF5 sidecars + `.log` logs) that the
unchanged `Results.artifacts` collection model surfaces directly.
The existing `valenx_bio::format::dcd` reader already inspects
collected NAMD DCD trajectories. A first-class MD-engine canonical
type — a typed force-field-and-integrator-state representation
spanning all six engines (Phase 5 LAMMPS / GROMACS, Phase 17
OpenMM, Phase 5.6 NAMD / sander / HOOMD-blue) — defers to a future
phase along with cross-engine restart-file converters and topology
swappers.

### Headless CLIs

**No new CLIs.** NAMD's `.dcd` trajectories are already
inspectable through the existing Phase 17 DCD reader; sander's
`.nc` NetCDF trajectories are inspectable through the user's
downstream Phase 5.5 cpptraj pipeline (cpptraj reads `.nc`
natively); HOOMD-blue's `.gsd` trajectories are inspectable
through the user's downstream HOOMD / freud / ovito pipeline.
NAMD `.xsc` extended-system + sander `.mdinfo` performance
sidecars are tabular text inspectable in any editor. A canonical
MD-engine CLI — cross-engine trajectory inspection, restart-file
diffing, topology cross-walks — defers to a future phase along
with the canonical type.

## Domain expansion

Phase 5.6 is a **sister-domain expansion of the Phase 5 GROMACS /
LAMMPS MD engine beachhead** — the same all-atom MD engine
surface broadened with three more established engines that cover
the corners GROMACS / LAMMPS don't reach. GROMACS is the de-facto
academic + industrial all-atom engine for biomolecular MD; LAMMPS
is the canonical materials-science particle-MD engine; OpenMM
(Phase 17) is the Python-native MD engine for embedded and
scriptable workflows; NAMD is the UIUC academic all-atom engine
for biomolecular MD pedagogy + HPC clusters; AmberTools sander is
the OSS portion of AMBER's MD engine + the canonical companion to
the Phase 5.5 cpptraj analyzer; HOOMD-blue is the Glotzer-lab
GPU-native particle simulator for soft-matter / coarse-grained
work. With Phase 5.6 the all-atom + particle MD-engine surface in
Valenx covers all six canonical shapes — single-binary CLI
(GROMACS, LAMMPS, NAMD, sander), Python-scripted (OpenMM,
HOOMD-blue), and the academic-license caveat (NAMD).

## What landed early

The implementation landed across 5
discrete implementation commits (3 adapters, the registry rollup,
the init-template extension) plus this docs pass — each landing
one adapter, the registry rollup, the init-template extension, or
the documentation pass. Every commit kept workspace clippy +
rustdoc clean.

## Acceptance checklist

- [x] `valenx-adapter-namd` adapter ships with case-input parser
      + 4 lib tests + 3 case-input tests covering parses-minimal /
      parses-with-overrides / rejects-empty-config plus the
      single-binary subprocess shape that composes
      `<binary> +p<processors> <config> [extras...]` with the
      `+p<N>` flag emitted as a single OsString (no space) so
      NAMD's flag-and-value parser reads them together
- [x] `valenx-adapter-amber-sander` adapter ships with case-input
      parser + 4 lib tests + 3 case-input tests covering parses-
      minimal / rejects-empty-output-basename / rejects-bad-paths,
      plus the single-binary subprocess shape that composes
      `sander -O -i <config> -p <topology> -c <coordinates> -o
      <basename>.out -r <basename>.rst -x <basename>.nc
      [extras...]` with the three input paths resolved against the
      case directory and validated on disk
- [x] `valenx-adapter-hoomd` adapter ships with case-input parser
      + 4 lib tests + 3 case-input tests covering parses-minimal /
      parses-with-overrides / rejects-non-py-script, plus the
      Python-script subprocess shape that enforces `.py`, stages
      the script, writes `valenx_params.json` with
      `output_basename`, and composes `<python> <staged_script>`
- [x] All 3 adapters wired into `valenx-app::init_registry` —
      live adapter count moves from 105 to **108** (alongside the
      Phase 5.7 MDTraj single-adapter and Phase 17.7 structure-
      tools trio that bring the total to **112**), rounding out
      the all-atom + GPU-native MD-engine surface that Phase 5
      GROMACS / LAMMPS opened
- [x] 3 MD-engine templates in `valenx-init` (`namd`, `sander`,
      `hoomd`), all round-tripping through `valenx-validate`
      (cross-binary roundtrip now sweeps **108 templates** clean
      alongside the Phase 5.7 single template and the Phase 17.7
      structure-tools trio)
- [x] NAMD probe pushes an `"academic"`-keyworded warning into
      `ProbeReport.warnings` containing both `"academic"` and
      `"non-commercial"` substrings, and `tool_license` surfaces
      as `"NAMD-License"` rather than mislabeling the custom UIUC
      NAMD terms as a recognised SPDX identifier
- [x] STATUS.md / README.md / ROADMAP.md / CHANGELOG.md /
      QUICKSTART.md updated
- [ ] **Future MD-engine work** — AMBER `pmemd.cuda` (proprietary
      AMBER per-site license; would require an HTTP-based license-
      check shape rather than subprocess; defer to a future
      proprietary-MD-engine phase), CHARMM (similar academic
      license + bespoke build process; defer), GENESIS (RIKEN's
      MD engine; defer), Tinker / Tinker-HP (polarisable force
      fields; defer), Desmond (Schrödinger's commercial MD engine
      with academic version; defer). Out of scope for this
      beachhead.

## Success metrics

| Metric                                                | Target          |
|-------------------------------------------------------|-----------------|
| New MD-engine adapter (template + tests)              | 1 day per       |
| All-atom + GPU-native MD-engine surface across 6 tools (Phase 5 GROMACS / LAMMPS + Phase 17 OpenMM + Phase 5.6 NAMD / sander / HOOMD-blue) | < tool baseline |

## Leads into

Phase 5.6 rounds out the all-atom + GPU-native MD-engine surface
that the user's bio / chemistry spec called out alongside the
Phase 5 GROMACS / LAMMPS beachhead. Combined with the existing
simulate-MD → analyze-trajectory → reweight-free-energy → fit-ENM
→ run-cpptraj-script → predict-structure → fold-RNA → analyze-DNA-
geometry → infer-tree-ML → infer-tree-Bayesian → simulate-popgen →
analyze-trees → simulate-pathway → reconstruct-3D → design-protein
→ validate loop, the **simulate-MD-NAMD → simulate-MD-sander →
simulate-MD-HOOMD → analyze-trajectory-MDAnalysis → analyze-
trajectory-MDTraj → reweight-free-energy → fit-ENM → run-cpptraj-
script → predict-structure → fold-RNA → analyze-DNA-geometry →
infer-tree-ML → infer-tree-Bayesian → simulate-popgen → analyze-
trees → simulate-pathway → reconstruct-3D → design-protein →
validate** loop now spans six MD engines (the Phase 5 GROMACS /
LAMMPS pair plus the Phase 17 OpenMM Python-native engine plus
the Phase 5.6 NAMD / sander / HOOMD-blue trio) feeding into the
Phase 5.5 / 17 PLUMED / ProDy / cpptraj / MDAnalysis (and the
Phase 5.7 MDTraj) post-MD analysis stack and the entire Phase 17
→ 39 biology / biotech / chemistry expansion — all in one Valenx
shell with no glue code beyond the existing case-toml / prepare /
run / collect path.

The natural follow-up is **Phase 5.8** — the deferred MD-engine
work called out above (AMBER `pmemd.cuda` requiring proprietary
license-check shape, CHARMM with a similar academic license +
bespoke build process, GENESIS for the RIKEN MD engine, Tinker /
Tinker-HP for polarisable force fields, Desmond for the
Schrödinger commercial MD engine), slotting in alongside the
existing NAMD / sander / HOOMD-blue adapters with the same
single-binary subprocess shape (sander / NAMD sister tools) or
Python-script subprocess shape (HOOMD-blue sister tools).
