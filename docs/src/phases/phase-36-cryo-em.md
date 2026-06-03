# Phase 36 — Cryo-EM

**Status:** 🟢 Live — RELION + EMAN2 + CTFFIND open the **first
cryo-electron microscopy reconstruction domain** in Valenx alongside
the Phase 17 / 17.5 / 18 / 25 / 27 / 27.5 / 28 / 30 / 32 / 34 biology
+ structure-prediction + protein-design + RNA-structure +
phylogenetics + systems-biology + docking beachheads and the Phase 24
cheminformatics expansion.

## Goal

Open the cryo-electron microscopy reconstruction domain in Valenx
with three established open-source tools that span the cryo-EM
pipeline — Bayesian 3D reconstruction at the core (RELION), broad-
spectrum image processing across the full cryo-EM workflow (EMAN2),
and per-micrograph contrast transfer function (CTF) estimation as
the canonical preprocessing step (CTFFIND). All three follow the
established Phase 18 BWA single-binary CLI pattern: the user supplies
particle stacks / reference maps / micrographs in, the adapter
streams progress out and walks the workdir for typed artifacts on
collect. RELION + EMAN2 each have many subcommands; this phase wraps
the most-used entry point per tool (RELION's `relion_refine`,
EMAN2's `e2refine_easy.py`, CTFFIND's `ctffind`); sister adapters
extending coverage defer to Phase 36.5. Phase 36 sits numerically
after Phase 34 and ships chronologically right after Phase 32 systems
biology — same chronological-vs-numerical convention used for Phase
17.5 / 24 / 28.

## Capability inventory

### Live adapters (3)

- **RELION** — Sjors Scheres' REgularised LIkelihood OptimisatioN
  suite (GPL-2.0). The de-facto Bayesian 3D reconstruction workhorse
  in cryo-EM facilities worldwide — particle classification, 3D
  refinement, CTF correction, post-processing. Single-binary
  subprocess shape with optional MPI wrapping: `relion_refine` for
  the single-process path, `mpirun -n <N> relion_refine_mpi` for
  multi-rank runs (RELION ships these as separate `_mpi`-suffixed
  binaries so the launcher knows which transport to use). Schema
  knobs: `particles` (`*_data.star` particle STAR file; required),
  `reference` (initial reference map `.mrc`; required),
  `output_basename` (becomes the `--o` prefix every output inherits
  so collect() walks deterministically; required), `angpix` (pixel
  size in Angstroms; required, > 0.0 and finite), `mpi_procs`
  (default 1, ≥ 1; > 1 switches to the MPI binary), `threads`
  (OpenMP threads per MPI rank, default 1, ≥ 1), `extra_args`.
  `prepare()` dispatches on `mpi_procs`: single-rank composes
  `relion_refine --i <particles> --ref <reference> --o
  <output_basename> --angpix <angpix> --j <threads> [extras...]`;
  multi-rank prepends `mpirun -n <mpi_procs> relion_refine_mpi ...`
  and surfaces a helpful install-hint `InvalidCase` ("install
  OpenMPI (`apt install openmpi-bin`, `brew install open-mpi`) or
  MPICH (`apt install mpich`) to enable multi-rank RELION runs") if
  `mpirun` isn't on PATH. `collect()` walks the workdir recursively
  for `<output_basename>*_class*.mrc` (`Native`, "RELION class
  average"), `<output_basename>*_data.star` (`Tabular`, "RELION
  particle assignments"), and `<output_basename>*_model.star`
  (`Log`, "RELION model summary"). Probe via
  `find_on_path(&["relion_refine"])`. Version range `4.0.0..6.0.0`
  (4.0 is the current stable line, predecessor 3.1; upper bound
  6.0 reserves room for the next major). `bio.relion.refine`
  ribbon capability.
- **EMAN2** — Steve Ludtke's broad-spectrum cryo-EM image-processing
  package (BSD-3-Clause). The "Swiss army knife" of single-particle
  cryo-EM: particle picking, 2D classification, initial-model
  building, 3D refinement (CTF corrected, simultaneous tilt-pair
  handling), and a sprawling Python toolkit (`e2*.py`) for
  everything in between. Single-binary subprocess shape: the adapter
  wraps `e2refine_easy.py`, EMAN2's high-level orchestrator that
  drives the rest of the toolkit. Schema knobs: `particles`
  (particle stack `.bdb` / `.hdf` / `.mrcs`; required), `model`
  (initial 3D model `.hdf` / `.mrc`; required), `output_basename`
  (becomes the `--path` argument; EMAN2 turns this into a
  `<basename>_NN/` results directory under the workdir; required),
  `target_resolution` (target resolution in Angstroms; required, >
  0.0 and finite), `symmetry` (point group — `"c1"` / `"d2"` /
  `"icos"` / etc.; required, default `"c1"`), `threads` (default 1,
  ≥ 1), `extra_args`. `prepare()` builds `e2refine_easy.py --input
  <particles> --model <model> --path <output_basename> --targetres
  <target_resolution> --sym <symmetry> --threads <threads>
  [extras...]`. `collect()` walks recursively for
  `<output_basename>_*/threed_*.hdf` (`Native`, "EMAN2
  reconstruction") and `<output_basename>_*/log.txt` (`Log`, "EMAN2
  log"). Probe via `find_on_path(&["e2refine_easy.py"])`. Version
  range `2.99.0..3.0.0` (the 2.99 line is the current pre-3.0
  stable release; upper bound 3.0 reserves room for the long-
  rumoured 3.x line). The `valenx-init` template ships with the
  alias `eman` alongside the canonical `eman2`. `bio.eman2.refine`
  ribbon capability.
- **CTFFIND** — Niko Grigorieff's contrast transfer function
  estimation tool (Janelia Research Campus non-commercial /
  academic-only license). The gold standard for fitting per-
  micrograph CTF parameters (defocus, astigmatism, phase shift) in
  single-particle cryo-EM workflows; RELION, cryoSPARC, EMAN2, and
  most automated pipelines all wrap CTFFIND under the hood as a
  preprocessing step. Single-binary subprocess shape with stdin-
  piped parameters: CTFFIND's CLI is interactive and prompts the
  user line-by-line for each microscope parameter on startup, so
  the adapter writes a parameters text file in the workdir during
  `prepare()` and uses a custom `run()` that pipes the file into
  the child's stdin via `Stdio::from(file)` (the shared
  `subprocess::run` helper closes stdin which makes CTFFIND read
  EOF before its first prompt and exit; the custom run path
  mirrors the MAFFT stdout-redirect pattern but for stdin). Schema
  knobs: `micrograph` (input micrograph `.mrc`; required),
  `output_diagnostic` (output diagnostic image `.mrc`; required),
  `output_txt` (output text file with CTF parameters; required),
  `pixel_size` (Angstroms; required, > 0.0 and finite), `voltage`
  (acceleration voltage in kV; default 300.0, > 0.0), `cs`
  (spherical aberration in mm; default 2.7, > 0.0),
  `amplitude_contrast` (fraction; required, in `0.0..=1.0` — 0.07
  typical for cryo, 0.1 for negative stain), `extra_args`.
  `prepare()` writes `ctffind_params.txt` containing one parameter
  per CTFFIND-v4.1 prompt in order (input image, output diagnostic,
  pixel size, voltage, Cs, amplitude contrast, plus standard
  defaults for box size / min res / max res / defocus search /
  expert sub-prompts) and stashes the filename under a sentinel env
  var (`VALENX_CTFFIND_PARAMS_FILE`). The custom `run()` recovers
  the filename, strips the sentinel from the env table so CTFFIND
  doesn't see it, opens the params file with `File::open()`, and
  hands the FD to the child — CTFFIND sees a pipe pre-loaded with
  one parameter per prompt and responds as if a human had typed
  each line. `collect()` reports `output_diagnostic` (`Native`,
  "CTFFIND diagnostic image") and `output_txt` (`Tabular`,
  "CTFFIND parameters"). Probe via `find_on_path(&["ctffind"])`.
  Version range `4.1.0..5.0.0` (CTFFIND4 is the long-running
  stable line; upper bound 5.0 reserves room for the announced
  CTFFIND5 line). `bio.ctffind.estimate` ribbon capability.
  **Academic-license-only** — probe pushes the literal string
  `"academic"` (the asserted anchor) into `ProbeReport.warnings`
  with the full reminder: "CTFFIND is licensed for non-commercial /
  academic use only. Confirm your use case complies with the
  Janelia license before redistributing CTF estimates or derived
  data." Tool license surfaces as `Janelia-License` rather than
  mislabeling as MIT / BSD.

### Canonical types

**No new canonical types.** All three adapters consume user-supplied
inputs (RELION particle STAR files + reference MRC volumes, EMAN2
particle stacks `.bdb` / `.hdf` / `.mrcs` + initial 3D models,
CTFFIND micrograph `.mrc` files) and emit user-readable artifacts
(RELION class-average MRC volumes, particle-assignment STAR files,
model-summary STAR files, EMAN2 `threed_*.hdf` reconstructions plus
per-run log files, CTFFIND diagnostic-image MRC plus per-micrograph
parameter text files) that the unchanged `Results.artifacts`
collection model surfaces directly. A first-class cryo-EM canonical
type — a generic `.mrc` volume / particle-stack / micrograph type
spanning all three back-ends — defers to a future phase along with
MRC readers and reconstruction visualizers.

### Headless CLIs

**No new CLIs.** RELION's `*_data.star` / `*_model.star` files,
EMAN2's `log.txt`, and CTFFIND's parameter text file are tabular
text inspectable in any editor or through the user's downstream
Python pipeline (`pandas`, `starfile`, `mrcfile`). The `.mrc`
volumes (RELION class averages, EMAN2 `threed_*.hdf`, CTFFIND
diagnostic image) and `.hdf` reconstructions are binary volumetric
data the user reads through `mrcfile` / `EMAN2.EMData()` / `h5py`
or visualises in Chimera / ChimeraX (the Phase 17 ChimeraX adapter
already covers headless rendering). A canonical cryo-EM CLI defers
to a future phase along with `.mrc` / `.hdf` reader work and
reconstruction-visualization integrations.

## Domain milestone

Phase 36 is the **first cryo-electron microscopy reconstruction
domain** to land in Valenx. The biology adapter family started with
Phase 17 (foundation — sequence / structure / trajectory canonical
types + classical MD + cheminformatics scripts) and expanded through
Phase 17.5 / 18 / 18.5 / 18.6 / 19 / 19.5 / 20 / 22 / 23 / 24 / 25 /
27 / 27.5 / 27.6 / 28 / 30 / 32 / 34 to cover sequence prediction,
alignment, RNA-seq, variant calling, single-cell, transcript
quantification, workflow orchestration, molecular viewers,
cheminformatics, quantum chemistry, protein design, EvolutionaryScale
models, RNA structure, phylogenetics, systems biology, and small-
molecule docking — but until Phase 36 the cryo-electron microscopy
reconstruction surface (Bayesian 3D refinement, broad-spectrum image
processing, CTF estimation) was absent. Phase 36 closes that gap with
three established open-source tools spanning the cryo-EM pipeline —
RELION at the Bayesian-reconstruction core, EMAN2 across the broader
single-particle workflow, and CTFFIND as the canonical CTF-estimation
preprocessing step.

## What landed early

The implementation rode subagent-driven-development across 5
discrete implementation commits (3 adapters, the registry rollup,
the init-template extension) plus this docs pass — each landing one
adapter, the registry rollup, the init-template extension, or the
documentation pass. Every commit kept workspace clippy + rustdoc
clean.

## Acceptance checklist

- [x] `valenx-adapter-relion` adapter ships with case-input parser
      + 4 lib tests + 5 case-input tests covering parses-minimal /
      parses-with-mpi-and-threads / rejects-zero-angpix / rejects-
      zero-mpi / rejects-zero-threads, plus the MPI-dispatch shape
      that prepends `mpirun -n <N>` and switches to
      `relion_refine_mpi` when `mpi_procs > 1`
- [x] `valenx-adapter-eman2` adapter ships with case-input parser
      + 4 lib tests + 4 case-input tests covering the `--path` /
      `--targetres` / `--sym` / `--threads` knob shape; init alias
      `eman` resolves to the same template as the canonical `eman2`
- [x] `valenx-adapter-ctffind` adapter ships with case-input parser
      + 4 lib tests + 5 case-input tests + 1 extra
      (`probe_warning_mentions_academic`) for 10 total, plus the
      stdin-piped-parameters shape that writes `ctffind_params.txt`
      in `prepare()` and routes the file into the child via
      `Stdio::from(file)` in a custom `run()` (the shared
      `subprocess::run` helper would close stdin and break CTFFIND's
      interactive-prompt CLI)
- [x] All 3 adapters wired into `valenx-app::init_registry` —
      live adapter count moves from 75 to **78**, opening the first
      cryo-electron microscopy reconstruction domain to ship in
      Valenx
- [x] 3 cryo-EM templates in `valenx-init` (`relion`, `eman2` with
      alias `eman`, `ctffind` with inline academic-license note),
      all round-tripping through `valenx-validate` (cross-binary
      roundtrip now sweeps **74 templates** clean)
- [x] STATUS.md / README.md / ROADMAP.md / CHANGELOG.md /
      QUICKSTART.md updated
- [ ] **Phase 36.5** — cisTEM (full single-particle cryo-EM
      pipeline UI; defer to sister-adapter expansion phase),
      SPHIRE (TransPhire / SPHIRE pipeline framework; defer),
      IMOD (cryo-ET reconstruction; different shape; defer),
      Bsoft (broad cryo-EM + electron crystallography toolkit;
      defer), Scipion (full cryo-EM pipeline orchestrator akin to
      Nextflow / Snakemake but cryo-EM-specific; defer), Frealign
      (predecessor to cisTEM; defer), motion correction (MotionCor2
      / RELION's own `relion_motioncorr`; defer to 36.5), particle
      picking (Topaz / crYOLO; defer), tomography (TomoBEAR /
      nextPYP; different shape — slot under cryo-ET separately).
      Out of scope for this beachhead.

## Success metrics

| Metric                                                | Target          |
|-------------------------------------------------------|-----------------|
| New cryo-EM adapter (template + tests)                | 1 day per       |
| Reconstruct + image-process + CTF-estimate loop across 3 tools | < tool baseline |

## Leads into

Phase 36 opens the cryo-electron microscopy reconstruction domain
that the user's bio / chemistry spec called out alongside the Phase
17 / 17.5 / 27 / 27.5 / 27.6 biology + protein-design stack and the
Phase 24 / 25 cheminformatics + quantum-chemistry expansion.
Combined with the existing fold → analyze → predict → infer-tree →
validate loop, the **estimate-CTF → process-images → reconstruct-3D
→ predict-structure → fold-RNA → infer-tree → simulate-pathway →
validate** loop now spans three cryo-EM tools (RELION, EMAN2,
CTFFIND) feeding into the Phase 32 systems-biology surface (COPASI,
BioNetGen, PhysiCell), the Phase 24 / 25 cheminformatics + quantum-
chemistry surface (DeepChem, Open Babel, Avogadro 2, Psi4, NWChem,
xTB), the Phase 17 / 17.5 prediction stack (ESMFold, OpenFold,
AlphaFold 2/3, ColabFold), the Phase 28 RNA-structure tools
(ViennaRNA, RNAstructure, NUPACK), and the Phase 30 phylogenetic-
tree builders (IQ-TREE, RAxML-NG, FastTree) — all in one Valenx
shell with no glue code beyond the existing case-toml / prepare /
run / collect path.

The natural follow-up is **Phase 36.5** — the deferred cryo-EM work
called out above (cisTEM as a full single-particle pipeline UI,
SPHIRE / TransPhire as pipeline frameworks, IMOD for cryo-ET
reconstruction, Bsoft for broad cryo-EM + electron-crystallography
coverage, Scipion as a full cryo-EM pipeline orchestrator akin to
Nextflow / Snakemake but cryo-EM-specific, Frealign as the cisTEM
predecessor, MotionCor2 / `relion_motioncorr` for motion correction,
Topaz / crYOLO for particle picking), slotting in alongside the
existing cryo-EM adapters with the same single-binary subprocess
shape (or the OpenMM / Scanpy Python-script-subprocess shape for
the Python-library tools). Tomography (TomoBEAR / nextPYP) sits in
a separate cryo-electron-tomography (cryo-ET) phase — the data
shape is different enough (tilt series rather than single
particles) to warrant a sister phase rather than 36.5 expansion.
See the out-of-scope section of
`docs/superpowers/plans/2026-04-30-cryo-em.md` for the full follow-
up phase list.
