# Phase 32.5 — Spatial stochastic

**Status:** 🟢 Live — Smoldyn + MCell round out the
**systems-biology / multiscale modeling surface** that Phase 32
COPASI / BioNetGen / PhysiCell opened, alongside the Phase 5.5 /
5.6 / 5.7 / 17 / 17.5 / 17.7 / 18 / 25 / 27 / 27.5 / 27.6 / 28 /
29 / 30 / 30.5 / 31 / 33 / 34 / 35 / 36 / 38 / 39 MD-engines +
MD-analysis-expansion + biology + structure-prediction +
protein-design + RNA-structure + population-genetics +
phylogenetics + Bayesian-phylogenetics + read-simulator +
synthetic-biology + docking + CRISPR-design + cryo-EM +
Rosetta-family + DNA-geometry beachheads and the Phase 24
cheminformatics expansion.

## Goal

Sister-adapter expansion of the existing Phase 32 systems-biology
trio (COPASI well-mixed CRN / BioNetGen rule-based / PhysiCell
agent-based). Round out the systems-biology / multiscale modeling
surface with the two canonical **spatial stochastic / cell-scale
reaction-diffusion** simulators that Phase 32 explicitly deferred
— **Smoldyn** (Andrews lab, the de-facto particle-based spatial
stochastic simulator that resolves individual molecules diffusing
and reacting in continuous 3D space) and **MCell** (Salk
Institute / Stiles, Bartol — the canonical Monte Carlo cell-
scale spatial stochastic simulator built around the MDL model
description language). Both adapters follow the established Phase
18 BWA single-binary CLI pattern: model file in, reaction-data +
trajectory artifacts out. Phase 32.5 sits numerically adjacent to
Phase 32 and ships chronologically right after Phase 17.7
structure tools — same chronological-vs-numerical convention used
for Phase 17.5 / 24 / 28 / 31 / 35 / 39 / 5.5 / 5.6 / 5.7.

## Capability inventory

### Live adapters (2)

- **Smoldyn** — Steve Andrews's spatial stochastic reaction-
  diffusion simulator (LGPL-2.1). Smoldyn resolves individual
  molecules as particles diffusing and reacting in continuous 3D
  space (no lattice discretisation), the canonical choice when
  the question is "where does each molecule actually end up over
  time" rather than "what is the well-mixed concentration vs. t"
  that Phase 32 COPASI's ODE / SSA covers. The configuration is
  a plain-text Smoldyn config file describing the simulation
  geometry (boundaries, surfaces, compartments), molecule species
  + diffusion coefficients, and per-pair / per-surface reactions.
  Single-binary subprocess shape (sister to Phase 18 BWA): the
  CLI is `smoldyn <config> [extras...]`. Schema knobs: `config`
  (Smoldyn `.txt` configuration file; required), `extra_args`
  (additional CLI arguments appended after the config path).
  `prepare()` resolves `config` against the case directory when
  relative, validates it exists on disk (returns `InvalidCase`
  with a helpful message when missing), and composes the
  invocation. `collect()` walks the workdir for `*.txt`
  (`Tabular`, "Smoldyn output table" — Smoldyn's per-step
  particle / reaction tables), `*.dat` (`Tabular`, "Smoldyn
  data" — reaction-event / molecule-position dumps the config
  may direct here), and `*.log` (`Log`, "Smoldyn log"). Probe
  via `find_on_path(&["smoldyn"])`. Version range
  `2.70.0..3.0.0` (Smoldyn 2.70 (2023) is the modern stable line
  shipping the contemporary surface-reaction model + lattice-Monte-
  Carlo mode; upper bound 3.0 reserves room for the long-promised
  next major). `bio.smoldyn.simulate` ribbon capability.
- **MCell** — Salk Institute / Stiles, Bartol's cell-scale
  Monte Carlo spatial stochastic simulator (GPL-2.0). MCell
  walks the user's `.mdl` (Model Description Language) model —
  geometry built from triangle meshes, molecule species with
  diffusion coefficients, surface / volume reactions, release
  patterns, observation counts — and runs Brownian-dynamics
  particle trajectories with Monte Carlo reaction sampling. The
  canonical use case is sub-cellular signaling (synaptic
  transmission, calcium dynamics, receptor binding) where the
  geometry is intricate enough that Smoldyn's continuous-space
  mode would be overkill but a well-mixed COPASI / BioNetGen
  treatment misses the spatial structure. Single-binary
  subprocess shape (sister to Smoldyn): the CLI is
  `mcell [-seed <N>] <mdl> [extras...]` with `-seed` as a
  separate token from its integer argument. Schema knobs: `mdl`
  (`.mdl` MCell model description file; required), `seed`
  (`Option<u32>` — when `Some(n)` the adapter emits `-seed` and
  `<n>` as TWO separate args; when `None` MCell picks its own
  seed and prints it on the run banner — same shape as the Phase
  29 SLiM `-s` and Phase 30.5 BEAST 2 `-seed` knobs),
  `extra_args` (additional CLI arguments appended after the MDL
  path). `prepare()` resolves `mdl` against the case directory
  when relative, validates it exists on disk, and composes the
  invocation, threading `-seed` + `<n>` as two separate OsStrings
  into the argv only when `seed` is `Some(_)`. `collect()` walks
  the workdir for `*.dat` (`Tabular`, "MCell reaction data" —
  per-observation count tables MCell writes from the model's
  REACTION_DATA_OUTPUT block), `*.dx` (`Native`, "MCell
  visualization data" — DReAMM / OpenDX visualization frames),
  and `*.log` (`Log`, "MCell log"). Probe via
  `find_on_path(&["mcell"])`. Version range `4.0.0..5.0.0`
  (MCell 4.0 (2022) is the modern Python-friendly C++ rewrite;
  the older 3.x line is deprecated; upper bound 5.0 reserves
  room for a future major). `bio.mcell.simulate` ribbon
  capability.

### Canonical types

**No new canonical types.** Both adapters consume user-supplied
inputs (Smoldyn `.txt` config files, MCell `.mdl` model
description files plus an optional integer random seed) and emit
user-readable artifacts (Smoldyn `.txt` per-step particle /
reaction tables, `.dat` reaction-event / molecule-position dumps,
`.log` run logs; MCell `.dat` per-observation count tables, `.dx`
DReAMM / OpenDX visualization frames, `.log` run logs) that the
unchanged `Results.artifacts` collection model surfaces directly.
A first-class spatial-stochastic canonical type — a typed
particle-trajectory + reaction-event-stream representation
spanning both back-ends — defers to a future phase along with
particle-trajectory visualizers and per-species per-time
concentration-field plotters.

### Headless CLIs

**No new CLIs.** Smoldyn's `.txt` / `.dat` outputs are tabular
text files inspectable in any editor or through the user's
downstream Python pipeline (`pandas`, `numpy`); MCell's `.dat`
per-observation count tables are similarly tabular text, and the
`.dx` OpenDX visualization frames are consumed by DReAMM, VMD,
or ParaView. A canonical spatial-stochastic CLI — particle-
trajectory inspection, per-species count diffing, reaction-
event-stream comparison — defers to a future phase along with
the canonical type.

## Domain expansion

Phase 32.5 is a **sister-adapter expansion of the Phase 32
systems-biology trio** — the same systems-biology / multiscale
modeling surface broadened with two more established tools that
cover the spatial-stochastic / cell-scale corner Phase 32 left
out. Phase 32 COPASI is the de-facto deterministic biochemical-
pathway / ODE simulator (well-mixed CRN); Phase 32 BioNetGen is
the de-facto rule-based combinatorial-signaling-network
expander; Phase 32 PhysiCell is the de-facto agent-based off-
lattice multicellular tissue simulator; Phase 32.5 Smoldyn is
the de-facto particle-based continuous-space spatial-stochastic
simulator (sub-cellular molecular trajectories); Phase 32.5
MCell is the de-facto Monte Carlo cell-scale spatial-stochastic
simulator with the MDL model-description language (synaptic /
sub-cellular signaling on intricate triangle-mesh geometry).
With Phase 32.5 the systems-biology / multiscale modeling
surface in Valenx covers all five canonical shapes — well-mixed
ODE / SSA (COPASI), rule-based combinatorial expansion
(BioNetGen), agent-based tissue (PhysiCell), particle-based
continuous-space spatial stochastic (Smoldyn), and Monte Carlo
cell-scale spatial stochastic (MCell).

## What landed early

The implementation landed across 3
discrete implementation commits (2 adapters plus the registry +
init-template rollup) plus this docs pass — each landing one
adapter, the registry rollup, the init-template extension, or
the documentation pass. Every commit kept workspace clippy +
rustdoc clean.

## Acceptance checklist

- [x] `valenx-adapter-smoldyn` adapter ships with case-input
      parser + 4 lib tests + 3 case-input tests covering parses-
      minimal / parses-with-extras / rejects-empty-config, plus
      the single-binary subprocess shape that composes
      `smoldyn <config> [extras...]` with `config` resolved
      against the case directory and validated on disk
- [x] `valenx-adapter-mcell` adapter ships with case-input
      parser + 4 lib tests + 3 case-input tests covering parses-
      minimal / parses-with-seed / rejects-empty-mdl, plus the
      single-binary subprocess shape that composes
      `mcell [-seed <N>] <mdl> [extras...]` with `-seed` and
      `<N>` emitted as two separate OsStrings only when `seed`
      is `Some(_)`
- [x] Both adapters wired into `valenx-app::init_registry` —
      live adapter count moves from 112 to 114 (alongside the
      Phase 40 microscopy trio and Phase 41 sequence-editors
      pair that bring the total to **119**), rounding out the
      systems-biology / multiscale modeling surface that Phase
      32 COPASI / BioNetGen / PhysiCell opened
- [x] 2 spatial-stochastic templates in `valenx-init` (`smoldyn`,
      `mcell`), all round-tripping through `valenx-validate`
      (cross-binary roundtrip now sweeps **115 templates** clean
      alongside the Phase 40 microscopy trio and Phase 41
      sequence-editors pair)
- [x] STATUS.md / README.md / ROADMAP.md / CHANGELOG.md /
      QUICKSTART.md updated
- [ ] **Phase 32.6** — sister-adapter expansion of Phase 32 +
      32.5: Tellurium / libRoadRunner (Python library, fits the
      Biopython subprocess pattern; defer to sister-adapter
      expansion phase), VCell (Java GUI app; `vcell-cli` exists
      but workflow is heavy; defer), E-Cell / Morpheus /
      CompuCell3D (niche cell / tissue simulators; defer),
      StochPy / libSBML / PySB (Python libraries; defer to a
      future Python-systems-biology expansion), URDME / STEPS
      (sister spatial-stochastic simulators with mesh-based
      reaction-diffusion solvers; defer alongside Smoldyn /
      MCell). Out of scope for this beachhead.

## Success metrics

| Metric                                                | Target          |
|-------------------------------------------------------|-----------------|
| New spatial-stochastic adapter (template + tests)     | 1 day per       |
| Particle-based + Monte-Carlo spatial stochastic loop across 2 tools | < tool baseline |

## Leads into

Phase 32.5 rounds out the systems-biology / multiscale modeling
surface that the user's bio / chemistry spec called out alongside
the Phase 32 COPASI / BioNetGen / PhysiCell beachhead. Combined
with the existing simulate-pathway → expand-rules → grow-tissue
→ simulate-MD → analyze-trajectory → reweight-free-energy →
fit-ENM → run-cpptraj-script → predict-structure → fold-RNA →
analyze-DNA-geometry → infer-tree-ML → infer-tree-Bayesian →
simulate-popgen → analyze-trees → reconstruct-3D → design-protein
→ validate loop, the **simulate-pathway → expand-rules →
grow-tissue → diffuse-particles → trace-MCell-trajectories →
simulate-MD → analyze-trajectory → reweight-free-energy →
fit-ENM → run-cpptraj-script → predict-structure → fold-RNA →
analyze-DNA-geometry → infer-tree-ML → infer-tree-Bayesian →
simulate-popgen → analyze-trees → reconstruct-3D → design-protein
→ validate** loop now spans five systems-biology / multiscale
modeling tools (the Phase 32 COPASI / BioNetGen / PhysiCell trio
plus the Phase 32.5 Smoldyn / MCell pair) feeding into the
existing Phase 5 / 5.6 GROMACS / LAMMPS / NAMD / sander / HOOMD-
blue MD engines, the Phase 5.5 / 5.7 / 17 PLUMED / ProDy /
cpptraj / MDTraj / MDAnalysis post-MD analysis stack, the Phase
17 / 17.5 / 17.7 prediction stack (ESMFold, OpenFold, AlphaFold
2/3, ColabFold, RoseTTAFold, OmegaFold, FoldSeek), the Phase 28
RNA-structure tools (ViennaRNA, RNAstructure, NUPACK), the Phase
29 population-genetics trio (SLiM, msprime, tskit), the Phase 30
phylogenetic-tree builders (IQ-TREE, RAxML-NG, FastTree), the
Phase 30.5 Bayesian-phylogenetics pair (BEAST 2, MrBayes), the
Phase 33 synthetic-biology trio (pySBOL, j5, Cello), the Phase
34 docking pair (AutoDock Vina, AutoDock 4), the Phase 35
CRISPR-design tools (CHOPCHOP, CRISPOR, Cas-OFFinder), the Phase
36 cryo-EM reconstruction tools (RELION, EMAN2, CTFFIND), the
Phase 38 Rosetta-family adapters (Rosetta, PyRosetta), and the
Phase 39 DNA-structural-geometry tools (X3DNA, Curves+, DSSR) —
all in one Valenx shell with no glue code beyond the existing
case-toml / prepare / run / collect path.

The natural follow-up is **Phase 32.6** — the deferred systems-
biology / spatial-stochastic work called out above (Tellurium /
libRoadRunner as a Biopython-style Python-script subprocess,
VCell with its `vcell-cli` batch mode, E-Cell / Morpheus /
CompuCell3D as additional spatial / cell simulators alongside
PhysiCell + Smoldyn + MCell, StochPy / libSBML / PySB as Python
libraries fitting the Phase 24 cheminformatics-style adapter
pattern, URDME / STEPS as additional spatial-stochastic
simulators with mesh-based reaction-diffusion solvers), slotting
in alongside the existing Smoldyn / MCell adapters with the same
single-binary subprocess shape (Smoldyn / MCell / VCell sister
tools) or the Python-script subprocess shape (Tellurium / StochPy
sister tools).
