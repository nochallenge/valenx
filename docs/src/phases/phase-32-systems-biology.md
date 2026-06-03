# Phase 32 — Systems biology

**Status:** 🟢 Live — COPASI + BioNetGen + PhysiCell open the **first
systems-biology / multiscale modeling domain** in Valenx alongside the
Phase 17 / 17.5 / 18 / 25 / 27 / 27.5 / 28 / 30 / 34 biology +
structure-prediction + protein-design + RNA-structure + phylogenetics +
docking beachheads and the Phase 24 cheminformatics expansion.

## Goal

Open the systems biology / multiscale modeling domain in Valenx with
three established open-source tools that span the systems-biology
tradeoff space — biochemical pathway / ODE simulation at the
deterministic end (COPASI), rule-based modeling for combinatorially-
complex signaling networks in the middle (BioNetGen), and agent-based
multicellular tissue simulation at the spatial / multiscale end
(PhysiCell). All three follow the established Phase 18 BWA single-
binary CLI pattern: model file in, results out. COPASI reads `.cps`
native or SBML `.xml` and runs the tasks defined inside; BioNetGen
reads `.bngl` rule-based models and emits `<basename>.net` /
`.gdat` / `.cdat` outputs through its Perl driver `BNG2.pl`;
PhysiCell models compile per-project to a project-specific C++
binary, so the adapter takes both the user's compiled `binary` path
and the run-time XML configuration. Phase 32 sits numerically after
Phase 30 and ships chronologically right after Phase 25 quantum
chemistry — same chronological-vs-numerical convention used for
Phase 17.5 / 24 / 28.

## Capability inventory

### Live adapters (3)

- **COPASI** — the COmplex PAthway SImulator (Artistic-2.0). The
  de-facto desktop suite for biochemical pathway and ODE-based
  systems-biology models, descended from the Gepasi lineage. Single-
  binary subprocess shape: COPASI's headless CLI is `CopasiSE`
  (capital `C-S-E`, "Self-Executing"), a task runner that reads a
  COPASI native `.cps` archive (or an SBML `.xml`) and executes the
  simulation / scan / fitting tasks defined inside. Schema knobs:
  `model` (`.cps` / `.sbml` / `.xml`; required), `report` (optional
  `--save <report>` target so the run output lands at a known path
  collect() can find without walking), `run_all` (default `false`;
  when `true` adds `--scheduled`, executing every task in the file
  rather than just the primary one), `extra_args`. `prepare()`
  composes `CopasiSE [--save <report>] <model> [--scheduled]
  [extras...]`. `collect()` reports the explicit `report` path when
  supplied (`Tabular`, "COPASI report") and otherwise walks the
  workdir top-level for `.csv` / `.txt` files (COPASI's tabular
  outputs). Probe via `find_on_path(&["CopasiSE"])`. Version range
  `4.40.0..5.0.0` (4.x is the long-running stable line; 4.40 is a
  recent floor that ships SBML L3v2 + the task scheduler).
  `bio.copasi.simulate` ribbon capability.
- **BioNetGen** — the rule-based modeling language and tool suite
  for combinatorially-complex signaling networks (MIT). The user
  writes BNGL (BioNetGen Language) files describing molecular
  species, sites, and reaction *rules*, and BioNetGen expands the
  rules into the underlying reaction network and (optionally)
  integrates it deterministically (ODE) or stochastically (SSA).
  Single-binary subprocess shape: `BNG2.pl` is the canonical Perl
  driver. Schema knobs: `model` (`.bngl`; required), `output_basename`
  (required; becomes the `-o` prefix every output file inherits so
  collect() walks deterministically), `generate_only` (default
  `false`; when `true` adds `--no-execute`, skipping simulate / scan
  / fitting actions and emitting just the expanded reaction network),
  `extra_args`. `prepare()` builds `BNG2.pl [--no-execute if
  generate_only] -o <output_basename> <model> [extras...]`.
  `collect()` walks the workdir top-level for `<output_basename>*.net`
  (`Native`, "BioNetGen reaction network"), `<output_basename>*.gdat`
  (`Tabular`, "BioNetGen species trajectories"), and
  `<output_basename>*.cdat` (`Tabular`, "BioNetGen concentrations") —
  `parameter_scan` per-trial variants share the basename prefix (e.g.
  `<basename>_001.gdat`) so the prefix-restricted walk picks them up
  too. Probe via `find_on_path(&["BNG2.pl"])`. Version range
  `2.8.0..3.0.0`. The `valenx-init` template ships with the alias
  `bng` alongside the canonical `bionetgen`.
  `bio.bionetgen.simulate` ribbon capability.
- **PhysiCell** — Paul Macklin's agent-based, off-lattice
  multicellular simulator (BSD-3-Clause). PhysiCell models tens to
  hundreds of thousands of individual cells (each an agent with
  state, mechanics, secretion, and phenotype) coupled to a reaction-
  diffusion microenvironment for substrates like oxygen and drugs.
  The canonical use case is tumour growth and immunology. Unlike a
  typical CLI tool, PhysiCell models compile to a project-specific
  C++ executable: the user clones the framework, edits the project's
  `custom_modules/` source, runs `make`, and ends up with e.g.
  `./project` next to the project directory. The adapter therefore
  takes both a `binary` path and the run-time XML configuration.
  Schema knobs: `binary` (required; the per-project compiled
  executable), `config` (required; the `.xml` settings file
  PhysiCell binaries accept as a positional argument), `extra_args`.
  `prepare()` validates `binary` and `config` exist on disk (returns
  `InvalidCase` with a helpful "PhysiCell models compile per-project
  — clone the framework, edit the project's `custom_modules/` source,
  run `make`, and point this field at the resulting executable."
  message if missing), then builds `<binary> <config> [extras...]`.
  `collect()` walks `output/`. PhysiCell drops a stack of per-
  snapshot files there: `output<N>.xml` (manifest), `output<N>_*.mat`
  (cell + microenvironment state in MATLAB v4 binary), and optional
  `*.csv` scalar summaries — typed `Native` ("PhysiCell tissue
  snapshot") for `.xml` / `.mat` and `Tabular` ("PhysiCell scalar
  table") for `.csv`. Probe via `find_on_path(&["physicell"])` —
  most installs won't have a generic `physicell` binary on PATH (the
  per-project build pattern means there isn't a canonical one), so
  the probe returns `ok = true` either way and attaches a warning
  that PhysiCell models compile per-project; the real validation
  happens in `prepare()` against the user's `binary` field. Version
  range `1.13.0..2.0.0`. `bio.physicell.simulate` ribbon capability.

### Canonical types

**No new canonical types.** All three adapters consume user-supplied
inputs (COPASI `.cps` / SBML `.xml` archives, BioNetGen `.bngl` rule-
based models, PhysiCell per-project compiled binaries + XML config)
and emit user-readable artifacts (CSV / TXT tabular reports,
reaction-network `.net` files, species-trajectory `.gdat` /
concentration `.cdat` tabular files, `.xml` / `.mat` per-snapshot
tissue state, per-cell scalar `.csv` summaries) that the unchanged
`Results.artifacts` collection model surfaces directly. A first-class
systems-biology canonical type — a generic SBML / BNGL / per-cell
state type spanning all three back-ends — defers to a future phase
along with SBML readers and tissue-snapshot visualizers.

### Headless CLIs

**No new CLIs.** COPASI's `.csv` / `.txt` reports, BioNetGen's
`.gdat` / `.cdat` species trajectories, and PhysiCell's per-cell
`.csv` summaries are tabular text files inspectable in any editor or
through the user's downstream Python pipeline (`pandas`, `numpy`).
PhysiCell's `.mat` per-snapshot tissue state is a MATLAB v4 binary
that the user reads through `scipy.io.loadmat` or the official
PhysiCell Python loader. A canonical systems-biology CLI defers to a
future phase along with SBML / BNGL / `.mat` reader work and
visualization integrations.

## Domain milestone

Phase 32 is the **first systems biology / multiscale modeling
domain** to land in Valenx. The biology adapter family started with
Phase 17 (foundation — sequence / structure / trajectory canonical
types + classical MD + cheminformatics scripts) and expanded through
Phase 17.5 / 18 / 18.5 / 18.6 / 19 / 19.5 / 20 / 22 / 23 / 24 / 25 /
27 / 27.5 / 27.6 / 28 / 30 / 34 to cover sequence prediction,
alignment, RNA-seq, variant calling, single-cell, transcript
quantification, workflow orchestration, molecular viewers,
cheminformatics, quantum chemistry, protein design,
EvolutionaryScale models, RNA structure, phylogenetics, and small-
molecule docking — but until Phase 32 the systems-biology /
multiscale-modeling surface (biochemical pathway ODE simulation,
rule-based signaling networks, agent-based multicellular tissue
simulation) was absent. Phase 32 closes that gap with three
established open-source tools that span the systems-biology tradeoff
space — COPASI at the deterministic-pathway end, BioNetGen for
rule-based combinatorial signaling, and PhysiCell at the spatial /
agent-based / multiscale end.

## What landed early

The implementation landed across 5
discrete implementation commits (3 adapters, the registry rollup,
the init-template extension) plus this docs pass — each landing one
adapter, the registry rollup, the init-template extension, or the
documentation pass. Every commit kept workspace clippy + rustdoc
clean.

## Acceptance checklist

- [x] `valenx-adapter-copasi` adapter ships with case-input parser
      + 4 lib tests + 3 case-input tests covering parses-minimal /
      parses-with-report / rejects-empty-model
- [x] `valenx-adapter-bionetgen` adapter ships with case-input parser
      + 4 lib tests + 3 case-input tests covering the
      `output_basename` + `generate_only` knob shape
- [x] `valenx-adapter-physicell` adapter ships with case-input parser
      + 4 lib tests + 3 case-input tests, plus the per-project-binary
      probe shape that returns `ok = true` with a warning when
      `physicell` isn't on PATH (the typical case, since PhysiCell
      models compile per-project) and validates `binary` + `config`
      at prepare time
- [x] All 3 adapters wired into `valenx-app::init_registry` —
      live adapter count moves from 72 to **75**, opening the first
      systems-biology / multiscale modeling domain to ship in Valenx
- [x] 3 systems-biology templates in `valenx-init` (`copasi`,
      `bionetgen` with alias `bng`, `physicell`), all round-tripping
      through `valenx-validate` (cross-binary roundtrip now sweeps
      **71 templates** clean)
- [x] STATUS.md / README.md / ROADMAP.md / CHANGELOG.md /
      QUICKSTART.md updated
- [ ] **Phase 32.5** — Tellurium / libRoadRunner (Python library,
      fits the Biopython subprocess pattern; defer to sister-adapter
      expansion phase), VCell (Java GUI app; `vcell-cli` exists but
      workflow is heavy; defer), E-Cell / Morpheus / CompuCell3D
      (niche; defer), Smoldyn / MCell (particle-based simulators;
      different shape; defer), StochPy / libSBML / PySB (Python
      libraries; future Phase 32.5). Out of scope for this beachhead.

## Success metrics

| Metric                                          | Target          |
|-------------------------------------------------|-----------------|
| New systems-biology adapter (template + tests)  | 1 day per       |
| Pathway + rule-based + agent-based loop across 3 tools | < tool baseline |

## Leads into

Phase 32 opens the systems-biology / multiscale modeling domain that
the user's bio / chemistry spec called out alongside the Phase 17 /
17.5 / 27 / 27.5 / 27.6 biology + protein-design stack and the Phase
24 / 25 cheminformatics + quantum-chemistry expansion. Combined with
the existing fold → analyze → predict → infer-tree → validate loop,
the **build-network → simulate-pathway → expand-rules → grow-tissue →
predict-structure → fold-RNA → infer-tree → validate** loop now spans
three systems-biology tools (COPASI, BioNetGen, PhysiCell) feeding
into the Phase 24 / 25 cheminformatics + quantum-chemistry surface
(DeepChem, Open Babel, Avogadro 2, Psi4, NWChem, xTB), the Phase 17 /
17.5 prediction stack (ESMFold, OpenFold, AlphaFold 2/3, ColabFold),
the Phase 28 RNA-structure tools (ViennaRNA, RNAstructure, NUPACK),
and the Phase 30 phylogenetic-tree builders (IQ-TREE, RAxML-NG,
FastTree) — all in one Valenx shell with no glue code beyond the
existing case-toml / prepare / run / collect path.

The natural follow-up is **Phase 32.5** — the deferred systems-
biology work called out above (Tellurium / libRoadRunner as a
Biopython-style Python-script subprocess, VCell with its `vcell-cli`
batch mode, E-Cell / Morpheus / CompuCell3D as additional spatial
simulators alongside PhysiCell, Smoldyn / MCell as particle-based
simulators, StochPy / libSBML / PySB as Python libraries fitting the
Phase 24 cheminformatics-style adapter pattern), slotting in
alongside the existing systems-biology adapters with the same single-
binary subprocess shape (or the OpenMM / Scanpy Python-script-
subprocess shape for the Python-library tools).
