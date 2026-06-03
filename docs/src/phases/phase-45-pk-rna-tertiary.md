# Phase 45 — Pharmacokinetics + RNA tertiary structure

**Status:** 🟢 Live — PK-Sim + SimRNA open the **first
pharmacokinetics / pharmacodynamics (PK/PD) modeling domain in
Valenx** plus the **first RNA tertiary (3D) structure prediction
domain in Valenx**, alongside the Phase 5.5 / 5.6 / 5.7 / 17 /
17.5 / 17.7 / 18 / 18.5 / 18.6 / 18.7 / 19 / 19.5 / 19.6 / 20 /
22 / 22.5 / 23 / 24 / 25 / 27 / 27.5 / 27.6 / 28 / 29 / 30 / 30.5
/ 31 / 32 / 32.5 / 33 / 34 / 35 / 35.5 / 35.6 / 36 / 38 / 39 / 40
/ 41 / 42 / 43 / 44.5 biology / biotech / chemistry beachheads.

## Goal

Open two NEW domains the existing bio surface doesn't reach with
two canonical open-source tools — **PK/PD pharmacokinetics**
(PK-Sim, the Open Systems Pharmacology suite's GPL-2.0 physiologically-
based PK simulator that's the de-facto open-source PBPK modeling
tool) and **RNA tertiary (3D) structure** (SimRNA, the Bujnicki
group's GPL-3.0 Monte Carlo RNA 3D folder that's the canonical
open-source coarse-grained RNA tertiary-structure predictor). The
two domains are different but Phase 45 ships them together because
each is a single canonical adapter and each opens a never-before-
seen domain in Valenx. PK-Sim is the **first PK/PD modeling category
in Valenx** — distinct from the Phase 32 / 32.5 systems-biology
ODE modeling which covers cellular-scale biochemical pathways
rather than whole-body absorption / distribution / metabolism /
excretion (ADME) drug pharmacokinetics. SimRNA is the **first RNA
tertiary structure prediction category in Valenx** — distinct from
the Phase 28 ViennaRNA / RNAstructure / NUPACK and Phase 44.5
mfold / EternaFold / LinearFold which all cover RNA secondary
(2D base-pairing) structure rather than the full 3D Cartesian
backbone tertiary structure SimRNA predicts. Both adapters follow
the established Phase 18 BWA single-binary CLI pattern: model file
in, simulation outputs out. Phase 45 sits numerically after Phase
44.5 RNA folding expansion and ships chronologically right after
Phase 35.6 edit-outcome prediction — same chronological-vs-
numerical convention used for Phase 17.5 / 24 / 28 / 31 / 35 / 39
/ 5.5 / 5.6 / 5.7 / 32.5 / 40 / 41 / 22.5 / 42 / 44.5 / 35.5 /
35.6.

## Capability inventory

### Live adapters (2)

- **PK-Sim** — Open Systems Pharmacology suite's physiologically-
  based PK (PBPK) simulator (GPL-2.0). PK-Sim is the de-facto
  open-source PBPK modeling tool, descended from the Bayer
  internal pharmacokinetic simulator opened to the community via
  the Open Systems Pharmacology Initiative. It models whole-body
  drug absorption / distribution / metabolism / excretion (ADME)
  using a physiologically-grounded compartmental representation —
  every major organ (liver, kidney, lung, gut, adipose, muscle,
  brain) is a compartment with its own blood flow, volume,
  partition coefficient, and metabolic capacity, and the user
  supplies a `.pksim5` project file (XML-based, authored in the
  PK-Sim GUI or programmatically through the OSP Python API)
  describing the drug + dosing protocol + simulated population.
  The headless `pksim` / `PKSim.CLI` CLI runs the simulation in
  batch mode and writes per-compartment concentration-time profiles
  + simulation metadata. Single-binary subprocess shape (sister
  to Phase 18 BWA / Phase 32.5 Smoldyn / Phase 5 GROMACS / Phase
  43 LinearDesign): the CLI is `pksim --project <project> --output
  <output_basename> [extras...]`. Schema knobs: `project`
  (`.pksim5` project file containing the PBPK model + dosing
  protocol + population specification; required — read in place
  from the case directory, no staging via `confined_join` since
  the PK-Sim binary reads the project file directly rather than
  the adapter staging it), `output_basename` (filename stem
  PK-Sim uses for the output CSVs + JSON metadata; required, non-
  empty), `extra_args` (additional CLI arguments appended after
  the canonical `--project` / `--output` pair so users can pin
  population sizes, simulation horizons, output-variable
  selections through PK-Sim's own flag surface). `prepare()`
  resolves `project` against the case directory when relative,
  validates the file exists on disk (returns `InvalidCase` with a
  helpful message when missing), and composes the invocation.
  `collect()` walks the workdir for `<output_basename>*.csv`
  (`Tabular`, "PK-Sim simulation results" — the canonical per-
  compartment concentration-time table PK-Sim writes for
  downstream Tlag / Cmax / AUC analysis), `<output_basename>*
  .json` (`Tabular`, "PK-Sim metadata" — simulation metadata
  including dosing protocol, simulated population statistics,
  per-compartment volume / blood-flow assumptions, model version),
  and `*.log` (`Log`). Probe via `find_on_path(&["pksim",
  "PKSim.CLI"])` (the modern OSP distribution ships both the
  generic `pksim` launcher and the .NET-style `PKSim.CLI` invoker
  for Windows installs). Version range `11.0.0..13.0.0` (PK-Sim
  11.x is the modern stable line shipping the contemporary OSP
  Python API + headless CLI; upper bound 13.0 reserves room for
  the next OSP major release). `bio.pksim.simulate` ribbon
  capability.
- **SimRNA** — Bujnicki group's coarse-grained Monte Carlo RNA
  tertiary-structure predictor (GPL-3.0). SimRNA predicts the
  full 3D Cartesian backbone of an RNA from its sequence — NOT
  just the 2D base-pairing pattern that the Phase 28 ViennaRNA /
  RNAstructure / NUPACK and Phase 44.5 mfold / EternaFold /
  LinearFold folders predict, but the actual 3D structure with
  per-residue Cartesian coordinates suitable for ChimeraX / PyMOL
  visualisation, structural alignment against PDB, MD simulation
  via Phase 5 GROMACS, etc. The model represents each nucleotide
  as five coarse-grained beads (one per phosphate / sugar / base
  ring) and runs replica-exchange Monte Carlo over the resulting
  reduced-coordinate energy landscape, sampling the conformational
  ensemble around the predicted minimum. Single-binary subprocess
  shape (sister to PK-Sim / Phase 18 BWA / Phase 32.5 Smoldyn):
  the CLI is `SimRNA -c <config> -s <sequence> -o <output_basename>
  -R <n_replicas> [extras...]`. Schema knobs: `config` (SimRNA
  configuration file specifying force-field parameters,
  temperature schedule, exchange acceptance criteria, simulation
  duration; required — read in place from the case directory,
  no staging since the SimRNA binary reads the config directly),
  `sequence` (`.seq` SimRNA-format sequence file; required —
  read in place, no staging), `output_basename` (filename stem
  SimRNA uses for the predicted PDB + trajectory + energy
  outputs; required, non-empty), `n_replicas` (`u32`, ≥ 1; number
  of replica-exchange replicas; default 1 — replica exchange
  helps escape local minima in the rugged RNA tertiary landscape;
  larger values produce more thorough sampling at proportionally
  larger compute cost), `extra_args`. `prepare()` resolves both
  `config` and `sequence` against the case directory when relative,
  validates each file exists on disk, and composes the invocation.
  `collect()` walks the workdir for `<output_basename>*.pdb`
  (`Native`, "SimRNA tertiary structure" — the predicted 3D
  Cartesian backbone in PDB format, lifted by the existing Phase
  17 PDB reader without any RNA-3D-specific code path),
  `<output_basename>*.trafl` (`Native`, "SimRNA trajectory" —
  the per-MC-step replica trajectory in SimRNA's compressed
  `.trafl` flat-array format, consumable by the Bujnicki group's
  `SimRNA_trafl2pdbs.py` post-processor for per-frame PDB
  extraction), `<output_basename>*.txt` (`Tabular`, "SimRNA
  energy log" — per-step energy + acceptance statistics for
  trajectory-quality assessment), and `*.log` (`Log`). Probe via
  `find_on_path(&["SimRNA", "simrna"])` (the modern Bujnicki
  distribution ships the binary as `SimRNA` with a lowercase
  `simrna` symlink for backwards compat). Version range
  `3.20.0..4.0.0` (SimRNA 3.20 is the modern stable line shipping
  the contemporary replica-exchange + force-field improvements;
  upper bound 4.0 reserves room for the next major bump). `bio
  .simrna.fold` ribbon capability.

### Canonical types

**No new canonical types.** Both adapters consume user-supplied
inputs (PK-Sim `.pksim5` PBPK project files + SimRNA configuration
+ sequence files) and emit user-readable artifacts (PK-Sim
concentration-time CSVs + simulation metadata JSON, SimRNA
predicted PDB tertiary structures + replica-exchange `.trafl`
trajectories + energy `.txt` logs) that the unchanged
`Results.artifacts` collection model surfaces directly. The
existing `valenx_bio::format::pdb` reader already inspects SimRNA's
predicted PDB outputs for chain / residue / atom counts — RNA
tertiary structures inherit the established protein-PDB
inspection path without any RNA-3D-specific code. A first-class
PK/PD canonical type — a typed compartmental concentration-time
representation with parsed Tlag / Cmax / AUC / clearance — defers
to a future PK/PD beachhead-expansion phase, as does a first-
class RNA-tertiary-structure canonical type with parsed base-
pair geometry / glycosidic-bond angles.

### Headless CLIs

**No new CLIs.** PK-Sim's `.csv` concentration-time tables +
`.json` simulation metadata, and SimRNA's `.pdb` tertiary
structures + `.trafl` trajectories + `.txt` energy logs are all
standard formats inspectable in any editor or through the user's
downstream Python pipeline (`pandas`, `numpy`, `MDAnalysis`,
`Biopython`). PDB outputs are inspectable through the existing
Phase 17 `valenx-pdb-info` CLI. A canonical PK/PD CLI — Tlag /
Cmax / AUC computation, dose-response curve plotting — and a
canonical RNA-tertiary CLI — base-pair geometry + glycosidic-
bond-angle inspection — defer to future PK/PD and RNA-tertiary
beachhead-expansion phases along with their canonical types.

## Domain expansion

Phase 45 opens **two new categories** in Valenx — PK/PD
pharmacokinetics and RNA tertiary structure — neither of which
existed in the bio surface before. The two domains are unrelated
biologically but ship together because each is a single canonical
adapter (PK-Sim is the dominant open-source PBPK tool; SimRNA is
the dominant open-source RNA-3D folder) and each is a structural
copy of the same Phase 18 BWA single-binary CLI shape. PK-Sim is
the **first pharmacokinetics / pharmacodynamics (PK/PD) modeling
category in Valenx** — the existing Phase 32 / 32.5 systems-
biology modeling tools (COPASI / BioNetGen / PhysiCell / Smoldyn /
MCell) cover cellular-scale biochemical pathway / signaling /
spatial-stochastic simulation but stop short of the whole-body
absorption / distribution / metabolism / excretion (ADME)
pharmacokinetics PK-Sim's PBPK compartmental representation
delivers. SimRNA is the **first RNA tertiary (3D) structure
prediction category in Valenx** — the existing Phase 28 ViennaRNA
/ RNAstructure / NUPACK trio and the Phase 44.5 mfold /
EternaFold / LinearFold expansion all cover RNA secondary (2D
base-pairing) structure but stop short of the full 3D Cartesian
backbone SimRNA predicts. Phase 45 ships them as a coupled
beachhead because pairing PBPK whole-body simulation with RNA-3D
structure prediction enables the modern mRNA-vaccine
pharmacokinetics workflow — predict the mRNA's 3D structure (Phase
45 SimRNA), check its 2D folding (Phase 28 / 44.5), simulate its
absorption / lipid-nanoparticle distribution / hepatic metabolism
in a target population (Phase 45 PK-Sim), and integrate with the
existing Phase 43 mRNA-design + Phase 44.5 RNA-folding loop. Each
domain remains its own beachhead with its own future expansion
phase (Phase 45.5 for PK/PD, Phase 28.5 for RNA tertiary) but the
single-adapter packing keeps the rollout concise.

## What landed early

The implementation landed across 3
discrete implementation commits (2 adapters plus the registry +
init-template rollup) plus this docs pass — each landing one
adapter, the registry rollup, the init-template extension, or
the documentation pass. Every commit kept workspace clippy +
rustdoc clean.

## Acceptance checklist

- [x] `valenx-adapter-pksim` adapter ships with case-input parser
      + 4 lib tests + 3 case-input tests covering parses-minimal /
      parses-with-extra-args / rejects-empty-project, plus the
      single-binary subprocess shape that composes `pksim --project
      <project> --output <output_basename> [extras...]` with
      `project` resolved against the case directory and validated
      on disk (read in place, no staging — same shape as Phase 18
      BWA's reference genome and Phase 43 LinearDesign's `protein`
      FASTA), and the `find_on_path(["pksim", "PKSim.CLI"])` probe
- [x] `valenx-adapter-simrna` adapter ships with case-input parser
      + 4 lib tests + 3 case-input tests covering parses-minimal /
      parses-with-overrides / rejects-bad-n-replicas, plus the
      single-binary subprocess shape that composes `SimRNA -c
      <config> -s <sequence> -o <output_basename> -R <n_replicas>
      [extras...]` with both `config` and `sequence` resolved
      against the case directory and validated on disk (read in
      place, no staging), and the `find_on_path(["SimRNA",
      "simrna"])` probe
- [x] All 2 adapters wired into `valenx-app::init_registry` —
      live adapter count moves from 139 to **141** as part of the
      13-adapter / 4-phase rollup, opening the first PK/PD
      modeling and RNA tertiary structure prediction categories
      in Valenx
- [x] 2 PK/PD + RNA-tertiary templates in `valenx-init` (`pksim`,
      `simrna`), all round-tripping through `valenx-validate`
      (cross-binary roundtrip now sweeps **136 templates** clean
      alongside the Phase 44.5 / 35.5 / 35.6 sister rollups)
- [x] STATUS.md / README.md / ROADMAP.md / CHANGELOG.md /
      QUICKSTART.md updated
- [ ] **Phase 45.5** — PK/PD beachhead expansion: MoBi (Open
      Systems Pharmacology's QSP / mechanistic-model sister to
      PK-Sim; defer pending integration with PK-Sim's project
      file format), PoPy (population PK / Bayesian fitting;
      defer), Pumas (Julia-based PK/PD; defer pending domain-
      coverage decision on the Julia ecosystem). Out of scope
      for this beachhead.
- [ ] **Phase 28.5** — RNA tertiary beachhead expansion: ROSIE
      ("Rosetta Online Server that Includes Everyone" — the
      Rosetta-family RNA 3D module; defer — already partially
      reachable via Phase 38 PyRosetta), Vfold3D (alternative
      3D RNA folder sister to SimRNA; defer pending upstream
      activity), RNAComposer (template-based 3D-folding alternative
      to the de-novo SimRNA approach; defer — proprietary in
      parts). Out of scope for this beachhead.

## Success metrics

| Metric                                                | Target          |
|-------------------------------------------------------|-----------------|
| New PK/PD adapter (template + tests)                  | 1 day per       |
| New RNA-tertiary adapter (template + tests)           | 1 day per       |
| Simulate-PBPK → predict-RNA-3D loop across 2 tools    | < tool baseline |

## Leads into

Phase 45 opens **two new domains** that the user's bio / chemistry
spec called out as missing categories — PK/PD pharmacokinetics
and RNA tertiary structure prediction. Combined with the existing
predict-cas9-indels → cross-check-indels → predict-missense-
pathogenicity → search-population-off-targets → design-base-edit
→ predict-base-outcome → design-prime-edit → score-pegRNA →
fold-mfold → fold-EternaFold → fold-LinearFold → optimize-codons
→ design-mRNA → predict-stability → render-Mol* → render-NGL →
run-Galaxy-workflow → run-WDL → run-CWL → run-Nextflow → run-
Snakemake → design-plasmid → view-alignment → process-image →
segment-cells → classify-pixels → simulate-pathway → expand-
rules → grow-tissue → diffuse-particles → trace-MCell-
trajectories → simulate-MD → analyze-trajectory → reweight-
free-energy → fit-ENM → run-cpptraj-script → predict-structure
→ fold-RNA → analyze-DNA-geometry → infer-tree-ML → infer-
tree-Bayesian → simulate-popgen → analyze-trees → reconstruct-
3D → design-protein → validate loop, the **simulate-PBPK →
predict-RNA-3D → predict-cas9-indels → cross-check-indels →
predict-missense-pathogenicity → search-population-off-targets
→ design-base-edit → predict-base-outcome → design-prime-edit
→ score-pegRNA → fold-mfold → fold-EternaFold → fold-LinearFold
→ optimize-codons → design-mRNA → predict-stability → render-
Mol* → render-NGL → run-Galaxy-workflow → run-WDL → run-CWL →
run-Nextflow → run-Snakemake → design-plasmid → view-alignment
→ process-image → segment-cells → classify-pixels → simulate-
pathway → expand-rules → grow-tissue → diffuse-particles →
trace-MCell-trajectories → simulate-MD → analyze-trajectory →
reweight-free-energy → fit-ENM → run-cpptraj-script → predict-
structure → fold-RNA → analyze-DNA-geometry → infer-tree-ML →
infer-tree-Bayesian → simulate-popgen → analyze-trees →
reconstruct-3D → design-protein → validate** loop now spans
PK/PD pharmacokinetics (PK-Sim) feeding into the existing Phase
43 mRNA-design + Phase 44.5 RNA-folding + Phase 35 / 35.5 / 35.6
CRISPR design + editing + outcome stacks, and RNA tertiary 3D
structure prediction (SimRNA) feeding into the existing Phase 23
PyMOL / VMD / IGV + Phase 42 Mol* / NGL Viewer 3D-structure
visualisation stack — opening the **canonical mRNA-vaccine
pharmacokinetics + tertiary-structure-aware design pipeline**:
target protein → optimized mRNA → predicted 2D fold → predicted
3D structure → simulated whole-body PBPK distribution → all in
one Valenx shell with no glue code beyond the existing case-toml /
prepare / run / collect path.

The natural follow-up is **Phase 45.5** for the PK/PD domain
(MoBi as the QSP / mechanistic-model sister to PK-Sim, PoPy for
population PK / Bayesian fitting, Pumas for the Julia-ecosystem
sister) and **Phase 28.5** for the RNA tertiary domain (ROSIE
for the Rosetta-family RNA 3D module sister to PyRosetta,
Vfold3D as an alternative de-novo folder, RNAComposer for the
template-based folding sister), slotting in alongside the
existing PK-Sim + SimRNA adapters with the same single-binary
CLI shape (PK-Sim sister tools) or new shape if upstream tools
require something novel.
