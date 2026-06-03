# Phase 5.5 — MD analysis expansion

**Status:** 🟢 Live — PLUMED + ProDy + cpptraj round out the
**MD-trajectory analysis surface** that Phase 17 MDAnalysis opened
alongside the Phase 17 / 17.5 / 18 / 25 / 27 / 27.5 / 27.6 / 28 /
29 / 30 / 30.5 / 31 / 32 / 33 / 34 / 35 / 36 / 38 / 39 biology +
structure-prediction + protein-design + RNA-structure + population-
genetics + phylogenetics + Bayesian-phylogenetics + read-simulator +
systems-biology + synthetic-biology + docking + CRISPR-design +
cryo-EM + Rosetta-family + DNA-geometry beachheads and the Phase 24
cheminformatics expansion.

## Goal

Sister-adapter expansion of the existing Phase 17 MDAnalysis adapter.
Round out the MD-trajectory analysis surface with three more
established open-source tools that span the post-MD analysis
tradeoff space — enhanced-sampling collective-variable evaluation +
free-energy reweighting (PLUMED, the de-facto plug-in that wraps
every major MD engine for biased-simulation / reweighting work),
protein-dynamics elastic-network / normal-mode analysis (ProDy, the
canonical Python toolkit for ENM / GNM / ANM and ensemble PCA), and
canonical AmberTools trajectory analysis via cpptraj's domain
language (cpptraj, the reference workhorse for `rms` / `radgyr` /
`hbond` / `clustering` over Amber-format trajectories). PLUMED +
cpptraj follow the established Phase 18 BWA single-binary CLI
pattern: trajectory + script in, analysis tables out. ProDy follows
the established Phase 17 Biopython Python-script subprocess shape:
the user supplies a Python script that imports the upstream package
and reads `valenx_params.json` for the parsed knobs. Phase 5.5 sits
numerically adjacent to the original Phase 5 MD beachhead and
ships chronologically right after Phase 39 DNA structural geometry —
same chronological-vs-numerical convention used for Phase 17.5 / 24
/ 28 / 31 / 35 / 39.

## Capability inventory

### Live adapters (3)

- **PLUMED** — the de-facto enhanced-sampling and free-energy
  plug-in that wraps every major MD engine (GROMACS, LAMMPS, AMBER,
  NAMD, OpenMM) (LGPL-3.0). PLUMED defines collective variables
  (RMSD, dihedrals, distances, contact maps), biases (metadynamics,
  well-tempered metad, umbrella sampling, ABF), and a reweighting
  framework that turns biased trajectories back into unbiased
  free-energy surfaces. The `plumed driver` sub-command runs PLUMED
  standalone over a pre-computed trajectory: read frames from
  `--mf_xtc <traj>`, evaluate the collective variables defined in
  `--plumed <plumed.dat>`, write COLVAR / bias / HILLS files into
  the workdir. Single-binary subprocess shape (sister to Phase 18
  BWA): the CLI is `plumed driver --plumed <plumed_dat> --mf_xtc
  <trajectory> --kt <kt> [extras...]`. Schema knobs: `plumed_dat`
  (PLUMED input file describing the collective variables and bias
  to compute; required), `trajectory` (XTC trajectory; required —
  users running DCD / TRR can swap to `--mf_dcd` / `--mf_trr` via
  `extra_args`), `output_basename` (filename stem PLUMED uses for
  COLVAR / bias outputs; required, non-empty), `kt` (`f64`, > 0.0
  and finite; PLUMED's `k_B T` in its energy units — kJ/mol by
  default; default 2.494 = room temperature 300 K; a zero or NaN
  `kt` would crash PLUMED's reweighting on the first frame),
  `extra_args`. `prepare()` resolves both paths against the case
  directory when relative, validates each file exists on disk
  (returns `InvalidCase` with a helpful message when missing), and
  composes the invocation. `collect()` walks the workdir for
  `<output_basename>*.dat` (`Tabular`, "PLUMED COLVAR output") and
  `<output_basename>*.bias` (`Tabular`, "PLUMED bias"). Probe via
  `find_on_path(&["plumed"])`. Version range `2.9.0..3.0.0` (PLUMED
  2.9 (2023) is the modern stable line — the `driver` sub-command,
  the metadynamics / OPES bias family, and the Python interface are
  all mature; upper bound 3.0 reserves room for the long-promised
  next major). `bio.plumed.analyze` ribbon capability.
- **ProDy** — Bahar lab's canonical Python library for protein
  dynamics (MIT). ProDy ships elastic-network models (ENM / GNM /
  ANM), normal-mode analysis, ensemble PCA, the NMD trajectory
  format consumed by VMD's NMWiz plug-in, and integrations with the
  BLAST / DALI / PDB databases. Python-script subprocess shape
  (sister to Phase 17 Biopython): the user supplies a Python script
  referenced from `[bio.prody].script` in `case.toml` that imports
  `prody` and reads `valenx_params.json` for the parsed knobs.
  Schema knobs: `script` (path to user-supplied Python script;
  required), `python` (interpreter name; default `"python3"`),
  `input_pdb` (input PDB; required), `output_basename` (filename
  stem ProDy uses for ENM / mode / NMD outputs; required, non-
  empty), `num_modes` (`u32`, ≥ 1; number of normal modes to
  compute; default 20), `cutoff` (`f64`, > 0.0 and finite; ENM
  contact cutoff in Å; default 15.0). `prepare()` stages the script
  + input PDB into the workdir under their original filenames so
  the script can resolve them via relative paths, then writes a
  flat `valenx_params.json` containing `input_pdb` (staged
  filename), `output_basename`, `num_modes`, and `cutoff`.
  `collect()` walks the workdir for `<output_basename>*.npz`
  (`Native`, "ProDy ENM modes"), `<output_basename>*.nmd` (`Native`,
  "ProDy NMD trajectory" — the NMD format consumed by VMD's NMWiz
  plug-in for normal-mode visualisation), and
  `<output_basename>*.csv` (`Tabular`, "ProDy table"). Probe via
  Python on PATH with an `import prody` check — when the import
  fails the probe still returns `ok = true` with a warning so users
  with ProDy installed under a non-standard module name aren't
  blocked. Version range `2.4.0..3.0.0` (ProDy 2.x is the modern
  stable line; 2.4 is the floor we test against; upper bound 3.0
  reserves room for an eventual major bump). `bio.prody.analyze`
  ribbon capability.
- **cpptraj** — AmberTools' canonical trajectory analysis tool
  (GPL-3.0). cpptraj reads Amber `.prmtop` / `.parm7` topologies
  plus `.nc` / `.dcd` / `.mdcrd` trajectories, runs an analysis
  script authored in cpptraj's domain language (`trajin`, `rms`,
  `radgyr`, `hbond`, `volume`, `clustering`, ...), and writes
  results into the workdir as `.dat` (per-frame tables), `.agr`
  (XmGrace plot data), or `.gnu` (gnuplot scripts). Single-binary
  subprocess shape (sister to Phase 18 BWA): the CLI is `cpptraj -p
  <topology> -i <script> [extras...]`. Schema knobs: `script`
  (`.ptraj` / `.cpptraj` analysis script; required), `topology`
  (Amber `.prmtop` / `.parm7`; required), `extra_args`. `prepare()`
  resolves both paths against the case directory when relative,
  validates each file exists on disk, and composes the invocation.
  `collect()` walks the workdir for `*.dat` (`Tabular`, "cpptraj
  analysis output"), `*.agr` (`Tabular`, "cpptraj XmGrace plot"),
  and `*.gnu` (`Log`, "cpptraj gnuplot script"). Probe via
  `find_on_path(&["cpptraj"])`. Version range `6.0.0..7.0.0`
  (cpptraj 6.x is the modern stable line shipped with AmberTools 23+
  (2023); upper bound 7.0 reserves room for the next major bump).
  `bio.cpptraj.analyze` ribbon capability.

### Canonical types

**No new canonical types.** All three adapters consume user-
supplied inputs (PLUMED collective-variable scripts + MD
trajectories, ProDy Python analysis scripts + input PDBs, cpptraj
domain-language scripts + Amber topologies + trajectories) and emit
user-readable artifacts (PLUMED COLVAR `.dat` / `.bias` tables,
ProDy `.npz` ENM-mode arrays + `.nmd` NMD-format trajectories +
`.csv` analysis tables, cpptraj `.dat` per-frame tables + `.agr`
XmGrace plots + `.gnu` gnuplot scripts) that the unchanged
`Results.artifacts` collection model surfaces directly. The
existing `valenx_bio::format::pdb` reader already inspects collected
PDB inputs for chain / residue / atom counts. A first-class
MD-analysis canonical type — a typed collective-variable / normal-
mode / per-frame-statistics representation spanning all three
back-ends and the existing Phase 17 MDAnalysis adapter — defers to
a future phase along with COLVAR plotters, normal-mode visualizers,
and per-statistic time-series viewers.

### Headless CLIs

**No new CLIs.** PLUMED's COLVAR / bias `.dat` files, ProDy's
`.npz` mode arrays + `.nmd` NMD-format trajectories + `.csv`
analysis tables, and cpptraj's `.dat` / `.agr` / `.gnu` outputs are
all standard tabular / NumPy formats inspectable in any editor or
through the user's downstream Python pipeline (`pandas`, `numpy`,
`prody.parseNMD`). Input PDBs are inspectable through the existing
Phase 17 `valenx-pdb-info` CLI. A canonical MD-analysis CLI —
COLVAR comparison, normal-mode diffing, per-statistic trace
inspection — defers to a future phase along with the canonical type.

## Domain expansion

Phase 5.5 is a **sister-adapter expansion of the Phase 17
MDAnalysis adapter** — the same MD-trajectory analysis surface
broadened with three more established tools that cover the corners
MDAnalysis doesn't reach. MDAnalysis is the de-facto Python library
for trajectory I/O + per-frame property calculation; PLUMED is the
de-facto plug-in for biased / enhanced-sampling work and free-
energy reweighting; ProDy is the de-facto Python library for
elastic-network / normal-mode protein dynamics; cpptraj is the
canonical AmberTools trajectory analysis CLI for the long tail of
`rms` / `radgyr` / `hbond` / `clustering` analyses. With Phase 5.5
the post-MD analysis surface in Valenx covers all four canonical
shapes — Python library API (MDAnalysis + ProDy), enhanced-sampling
plug-in CLI (PLUMED), and AmberTools domain-language CLI (cpptraj).

## What landed early

The implementation landed across 5
discrete implementation commits (3 adapters, the registry rollup,
the init-template extension) plus this docs pass — each landing one
adapter, the registry rollup, the init-template extension, or the
documentation pass. Every commit kept workspace clippy + rustdoc
clean.

## Acceptance checklist

- [x] `valenx-adapter-plumed` adapter ships with case-input parser
      + 4 lib tests + 4 case-input tests covering parses-minimal /
      parses-with-overrides / rejects-zero-kt / rejects-empty-
      trajectory, plus the single-binary subprocess shape that
      composes `plumed driver --plumed <plumed_dat> --mf_xtc
      <trajectory> --kt <kt> [extras...]` with both files resolved
      against the case directory and validated on disk
- [x] `valenx-adapter-prody` adapter ships with case-input parser
      + 4 lib tests + 4 case-input tests covering parses-minimal /
      parses-with-overrides / rejects-bad-num-modes / rejects-bad-
      cutoff, plus the Python-script subprocess shape that stages
      script + input PDB and writes `valenx_params.json` with
      `input_pdb` (staged filename), `output_basename`, `num_modes`,
      and `cutoff`
- [x] `valenx-adapter-cpptraj` adapter ships with case-input parser
      + 4 lib tests + 3 case-input tests covering parses-minimal /
      rejects-empty-script / rejects-empty-topology, plus the
      single-binary subprocess shape that composes `cpptraj -p
      <topology> -i <script> [extras...]` with both files resolved
      against the case directory and validated on disk
- [x] All 3 adapters wired into `valenx-app::init_registry` —
      live adapter count moves from 94 to **97** (alongside the
      Phase 33 synthetic-biology trio that brings the total to
      **100**), rounding out the post-MD analysis surface that
      Phase 17 MDAnalysis opened
- [x] 3 MD-analysis templates in `valenx-init` (`plumed`, `prody`,
      `cpptraj`), all round-tripping through `valenx-validate`
      (cross-binary roundtrip now sweeps **96 templates** clean
      alongside the Phase 33 synthetic-biology trio)
- [x] STATUS.md / README.md / ROADMAP.md / CHANGELOG.md /
      QUICKSTART.md updated
- [ ] **Phase 5.6** — sister-adapter expansion of Phase 5.5:
      AMBER `cpptraj`'s GPU-accelerated sibling `pmemd.cuda`
      (already part of AmberTools but a different shape — a real
      MD engine rather than a trajectory analyzer; defer to a
      future MD-engine phase), MDTraj (sister Python trajectory
      analyzer to MDAnalysis; defer), nMOLDYN (neutron-scattering
      observables from MD trajectories; defer), OPLS-AA-style
      free-energy GPU engines (Lomap2, alchemical perturbations;
      defer to a docking / free-energy phase), CHARMM-GUI
      (web-fronted CHARMM input generator; defer), HOOMD-blue
      (GPU-native particle simulator; defer to a future MD-engine
      phase). Out of scope for this beachhead.

## Success metrics

| Metric                                                | Target          |
|-------------------------------------------------------|-----------------|
| New MD-analysis adapter (template + tests)            | 1 day per       |
| Enhanced-sampling reweighting + ENM normal-mode + cpptraj per-frame analysis loop across 3 tools | < tool baseline |

## Leads into

Phase 5.5 rounds out the post-MD analysis surface that the user's
bio / chemistry spec called out alongside the Phase 17 MDAnalysis
beachhead. Combined with the existing predict-structure → fold-RNA
→ analyze-DNA-geometry → infer-tree-ML → infer-tree-Bayesian →
simulate-popgen → analyze-trees → simulate-pathway → reconstruct-
3D → design-protein → validate loop, the **simulate-MD → analyze-
trajectory → reweight-free-energy → fit-ENM → run-cpptraj-script
→ predict-structure → fold-RNA → analyze-DNA-geometry → infer-
tree-ML → infer-tree-Bayesian → simulate-popgen → analyze-trees →
simulate-pathway → reconstruct-3D → design-protein → validate**
loop now spans four MD-analysis tools (the Phase 17 MDAnalysis
adapter plus PLUMED, ProDy, cpptraj) feeding into the existing Phase
5 GROMACS / LAMMPS MD engines, the Phase 17 / 17.5 prediction stack
(ESMFold, OpenFold, AlphaFold 2/3, ColabFold), the Phase 28 RNA-
structure tools (ViennaRNA, RNAstructure, NUPACK), the Phase 29
population-genetics trio (SLiM, msprime, tskit), the Phase 30
phylogenetic-tree builders (IQ-TREE, RAxML-NG, FastTree), the Phase
30.5 Bayesian-phylogenetics pair (BEAST 2, MrBayes), the Phase 32
systems-biology surface (COPASI, BioNetGen, PhysiCell), the Phase
33 synthetic-biology trio (pySBOL, j5, Cello), the Phase 34 docking
pair (AutoDock Vina, AutoDock 4), the Phase 35 CRISPR-design tools
(CHOPCHOP, CRISPOR, Cas-OFFinder), the Phase 36 cryo-EM
reconstruction tools (RELION, EMAN2, CTFFIND), the Phase 38 Rosetta-
family adapters (Rosetta, PyRosetta), and the Phase 39 DNA-
structural-geometry tools (X3DNA, Curves+, DSSR) — all in one
Valenx shell with no glue code beyond the existing case-toml /
prepare / run / collect path.

The natural follow-up is **Phase 5.6** — the deferred MD-analysis
work called out above (MDTraj as a sister Python trajectory
analyzer, nMOLDYN for neutron-scattering observables, alchemical
free-energy engines like Lomap2 sister to PLUMED's reweighting,
CHARMM-GUI for the web-fronted input generator), slotting in
alongside the existing PLUMED / ProDy / cpptraj adapters with the
same Python-script subprocess shape (ProDy sister tools) or single-
binary subprocess shape (PLUMED / cpptraj sister tools).
