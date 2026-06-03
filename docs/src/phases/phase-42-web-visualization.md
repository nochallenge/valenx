# Phase 42 — Web visualization

**Status:** 🟢 Live — Mol* + NGL Viewer open the **first modern
web 3D molecular visualization domain** in Valenx alongside the
Phase 5.5 / 5.6 / 5.7 / 17 / 17.5 / 17.7 / 18 / 18.5 / 18.6 / 18.7
/ 19 / 19.5 / 19.6 / 20 / 22 / 22.5 / 23 / 24 / 25 / 27 / 27.5 /
27.6 / 28 / 29 / 30 / 30.5 / 31 / 32 / 32.5 / 33 / 34 / 35 / 36 /
38 / 39 / 40 / 41 biology / biotech / chemistry beachheads.

## Goal

Open the **modern web 3D molecular visualization** corner of the
bio surface in Valenx with two established open-source WebGL
viewers that span the web-visualization tradeoff space — the
canonical PDBe / RCSB modern viewer that powers the structural-
biology web (Mol*, the EMBL-EBI / RCSB-led MIT-licensed WebGL
toolkit that has become the de-facto modern molecular viewer
embedded in the PDB / PDBe / AlphaFold DB / ESM Atlas web
properties since the late 2010s), and the Rose lab's WebGL
framework (NGL Viewer, Alexander Rose's MIT-licensed
high-performance WebGL framework for molecular visualization that
predated Mol* and still powers a large fraction of the
Jupyter-friendly notebook visualization ecosystem via its
`nglview` Python binding). Both are JavaScript browser libraries
in their primary distribution form; we wrap them via their
**Python bindings** (`molstar` and `nglview`, the official PyPI-
distributed Python interfaces both projects maintain) so they
slot into the existing Python-script subprocess pattern Phase 17
Biopython / Phase 19.5 Scanpy / Phase 33 pySBOL / Phase 41 pydna
established. The user supplies a `.py` script that imports the
binding to compose a state file (Mol* `.molj` JSON state),
rendered output (PNG image, HTML viewer page), or notebook-
embedded view; the adapter stages the script, writes a
`valenx_params.json` with the parsed knobs, and runs `<python>
<staged_script>`. Phase 42 sits numerically after Phase 41
sequence editors and ships chronologically right after Phase 22.5
workflow expansion — same chronological-vs-numerical convention
used for Phase 17.5 / 24 / 28 / 31 / 35 / 39 / 5.5 / 5.6 / 5.7 /
32.5 / 40 / 41 / 22.5.

## Capability inventory

### Live adapters (2)

- **Mol*** — the EMBL-EBI / RCSB-led modern WebGL molecular
  viewer (MIT). Mol* (pronounced "Mol-Star") has become the de-
  facto modern molecular viewer embedded in the PDB / PDBe /
  AlphaFold DB / ESM Atlas web properties since the late 2010s,
  superseding the older NGL / LiteMol / PV / 3Dmol viewers as
  the canonical browser-side renderer for protein / RNA /
  small-molecule structures, electron-density maps, cryo-EM
  volumes, and AlphaFold per-residue confidence overlays. Wrapped
  via the `molstar` Python binding (the official PyPI-distributed
  Python interface) so it slots into the existing Python-script
  subprocess pattern (sister to Phase 17 Biopython, Phase 19.5
  Scanpy, Phase 33 pySBOL, Phase 41 pydna). The user supplies a
  `.py` script referenced from `[bio.molstar].script` in
  `case.toml` that imports `molstar` and reads
  `valenx_params.json` for the parsed knobs. Schema knobs:
  `script` (path to user-supplied Python script; required, `.py`
  enforced), `python` (interpreter name; default `"python3"`),
  `input_structure` (`Option<PathBuf>` — optional input
  structure file the script can use as the molecule to render —
  `.pdb` / `.cif` / `.mmcif`; `None` when the script fetches
  from the PDB / generates the structure inline), `output_basename`
  (filename stem the user's script uses for outputs — surfaced
  here so collect() can label artifacts uniformly; required,
  non-empty). `prepare()` enforces the `.py` extension on the
  script, resolves `script` and the optional `input_structure`
  against the case directory when relative, stages both into the
  workdir under their original filenames so the script can
  resolve them via relative paths, then writes a flat hand-rolled
  `valenx_params.json` containing `output_basename` always plus
  `input_structure` (staged filename) only when set — the key is
  omitted entirely when `None` rather than emitted as `null`,
  matching the hand-rolled JSON convention the rest of the bio
  adapters use (Phase 19.6 Seurat / AnnData, Phase 27.5 ESM-IF,
  Phase 41 pydna). `collect()` walks the workdir for
  `<output_basename>*.html` (`Native`, "Mol* viewer HTML" — the
  canonical interactive HTML viewer page Mol* writes for
  embeddable / archival viewing), `<output_basename>*.molj`
  (`Native`, "Mol* state file" — the JSON state format that
  captures the entire viewer state, including representations /
  colourings / camera angles, for reproducible replay in any
  Mol*-conformant viewer), `<output_basename>*.png` (`Native`,
  "Mol* rendered image" — for scripts that emit a static
  publication-quality image render), and `*.log` (`Log`).
  Probe via Python on PATH then `<python> -c "import molstar"`
  — when the `import molstar` check fails the probe still
  returns `ok = true` with a targeted `"probe found python on
  PATH but could not import molstar — install with pip install
  molstar"` warning so users with Python ready but no `molstar`
  package see the install hint without failing the probe (sister
  to the Phase 19.5 scanpy / scvi / Phase 19.6 AnnData / Phase
  5.6 HOOMD-blue / Phase 5.7 MDTraj / Phase 41 pydna probe
  convention). Version range `3.0.0..5.0.0` (Mol* 3.x is the
  modern post-redesign line shipping the contemporary
  representation API + state file format; upper bound 5.0
  reserves room for the next two majors of this actively-
  maintained EMBL-EBI / RCSB project). `bio.molstar.view`
  ribbon capability.
- **NGL Viewer** — the Rose lab's high-performance WebGL
  framework for molecular visualization (MIT). NGL Viewer
  predated Mol* and still powers a large fraction of the
  Jupyter-friendly notebook visualization ecosystem via its
  `nglview` Python binding — the canonical choice for embedded-
  in-notebook 3D molecular views, especially for trajectory-
  scrubbing workflows where the user wants a MDAnalysis /
  MDTraj / Bio.PDB Universe / Trajectory directly visualized in
  a notebook cell. Wrapped via the `nglview` Python binding (the
  Jupyter-friendly Python interface the project maintains) so it
  slots into the existing Python-script subprocess pattern
  (sister to Mol*). The user supplies a `.py` script referenced
  from `[bio.ngl].script` in `case.toml` that imports `nglview`
  and reads `valenx_params.json` for the parsed knobs. Schema
  knobs: `script` (path to user-supplied Python script; required,
  `.py` enforced), `python` (interpreter name; default
  `"python3"`), `input_structure` (`Option<PathBuf>` — optional
  input structure file the script can use as the molecule to
  render — `.pdb` / `.cif` / `.mmcif` / and the wider format set
  NGL accepts via its readers; `None` when the script fetches
  from the PDB / generates the structure inline),
  `output_basename` (filename stem; required, non-empty).
  `prepare()` enforces the `.py` extension on the script,
  resolves `script` and the optional `input_structure` against
  the case directory when relative, stages both into the workdir
  under their original filenames so the script can resolve them
  via relative paths, then writes a flat hand-rolled
  `valenx_params.json` containing `output_basename` always plus
  `input_structure` (staged filename) only when set — the key is
  omitted entirely when `None` rather than emitted as `null`
  (same hand-rolled JSON convention as Mol*, Phase 19.6 Seurat /
  AnnData, Phase 27.5 ESM-IF, Phase 41 pydna). `collect()` walks
  the workdir for `<output_basename>*.html` (`Native`, "NGL
  viewer HTML" — the canonical interactive HTML viewer page NGL
  writes for embeddable / archival viewing),
  `<output_basename>*.png` (`Native`, "NGL rendered image" —
  for scripts that emit a static publication-quality image
  render), `<output_basename>*.json` (`Tabular`, "NGL state
  JSON" — JSON state captures of the viewer configuration), and
  `*.log` (`Log`). Probe via Python on PATH then `<python> -c
  "import nglview"` — when the `import nglview` check fails the
  probe still returns `ok = true` with a targeted `"probe found
  python on PATH but could not import nglview — install with
  pip install nglview"` warning so users with Python ready but
  no `nglview` package see the install hint without failing the
  probe (sister to the Mol* probe convention). Version range
  `3.0.0..5.0.0` (NGL 3.x is the modern actively-maintained
  line shipping the contemporary `nglview` Jupyter integration;
  upper bound 5.0 reserves room for the next two majors).
  `bio.ngl.view` ribbon capability.

### Canonical types

**No new canonical types.** Both adapters consume user-supplied
inputs (Python scripts that import `molstar` / `nglview`, plus
optional starting `.pdb` / `.cif` / `.mmcif` structure files) and
emit user-readable artifacts (HTML viewer pages, PNG image
renders, Mol* `.molj` state files, NGL state JSON) that the
unchanged `Results.artifacts` collection model surfaces directly.
The existing `valenx_bio::format::pdb` reader already inspects
collected PDB inputs for chain / residue / atom counts. A first-
class web-visualization canonical type — a typed viewer state
representation spanning Mol* `.molj` + NGL state JSON, with
parsed representation / colouring / camera-angle graphs — defers
to a future phase along with viewer-state diff tools and per-
representation visual-inspection CLIs.

### Headless CLIs

**No new CLIs.** Mol*'s `.html` viewer pages, `.molj` state
files, and `.png` image renders are all standard formats
inspectable in any browser, JSON viewer, or image viewer; NGL
Viewer's `.html` viewer pages, `.png` image renders, and
`.json` state captures are the same. Input PDB / CIF structures
are inspectable through the existing Phase 17 `valenx-pdb-info`
CLI. A canonical web-visualization CLI — Mol* / NGL viewer-state
diffing, per-representation comparison, headless screenshot
verification — defers to a future phase along with the canonical
type.

## Domain milestone

Phase 42 is the **first modern web 3D molecular visualization
domain** to land in Valenx. The biology adapter family started
with Phase 17 (foundation — sequence / structure / trajectory
canonical types + classical MD + cheminformatics scripts +
ChimeraX as the first script-driven molecular renderer) and
expanded through Phase 5.5 / 5.6 / 5.7 / 17.5 / 17.7 / 18 / 18.5
/ 18.6 / 18.7 / 19 / 19.5 / 19.6 / 20 / 22 / 22.5 / 23 / 24 / 25
/ 27 / 27.5 / 27.6 / 28 / 29 / 30 / 30.5 / 31 / 32 / 32.5 / 33 /
34 / 35 / 36 / 38 / 39 / 40 / 41 to cover MD-trajectory analysis
expansion, bio MD engines, MDTraj, sequence prediction,
structure-tools expansion, alignment, RNA-seq, alignment-toolkit
expansion, variant calling, single-cell, single-cell expansion,
transcript quantification, workflow orchestration, workflow
expansion, molecular viewers (Phase 23 — desktop-app sister to
Phase 17 ChimeraX with PyMOL / VMD / IGV `igvtools`),
cheminformatics, quantum chemistry, protein design, protein-
design expansion, EvolutionaryScale models, RNA structure,
population genetics, phylogenetics, Bayesian phylogenetics,
sequencing read simulation, systems biology, spatial-stochastic
simulators, synthetic biology, small-molecule docking, CRISPR
design, cryo-EM reconstruction, Rosetta protein modeling, DNA
structural geometry, microscopy, and sequence editing — but
until Phase 42 the **modern web 3D molecular visualization**
surface (browser-based WebGL viewers wrapped via their Python
bindings for headless artifact emission) was absent. Phase 23
covered the desktop / Tcl / `.pml` script side of molecular
visualization (PyMOL / VMD / IGV), while Phase 17 ChimeraX
covered the `.cxc` script-driven canonical desktop viewer; Phase
42 closes the web-viewer gap with two established open-source
WebGL toolkits spanning the web-visualization tradeoff space —
Mol* at the modern PDBe / RCSB embedded-viewer end, and NGL
Viewer as the Jupyter-friendly notebook-visualization framework
that closes the loop on the entire Phase 17 / 17.5 / 17.7
prediction stack and the Phase 5.5 / 5.7 trajectory analysis
stack by giving their structure / trajectory outputs a modern
WebGL visualization front end.

## What landed early

The implementation rode subagent-driven-development across 4
discrete implementation commits (2 adapters plus the registry +
init-template rollup) plus this docs pass — each landing one
adapter, the registry rollup, the init-template extension, or
the documentation pass. Every commit kept workspace clippy +
rustdoc clean.

## Acceptance checklist

- [x] `valenx-adapter-molstar` adapter ships with case-input
      parser + 4 lib tests + 3 case-input tests covering parses-
      minimal / parses-with-input-structure / rejects-non-py-
      script, plus the Python-script subprocess shape that
      enforces `.py`, stages script + optional input_structure,
      writes `valenx_params.json` with `output_basename` always
      plus `input_structure` (staged filename) only when set —
      key omitted entirely when `None` rather than emitted as
      `null`, matching the hand-rolled JSON convention the rest
      of the bio adapters use, plus the Python on PATH + `import
      molstar` probe with `"probe found python on PATH but could
      not import molstar"` warning when the import fails
- [x] `valenx-adapter-ngl` adapter ships with case-input parser
      + 4 lib tests + 3 case-input tests covering parses-minimal
      / parses-with-input-structure / rejects-non-py-script,
      plus the Python-script subprocess shape that mirrors Mol*
      with `import nglview` probe, the same hand-rolled
      `valenx_params.json` shape (key omitted when `None`), and
      the NGL-specific collect() filter set
- [x] Both adapters wired into `valenx-app::init_registry` —
      live adapter count moves from 122 to **124** alongside the
      Phase 22.5 workflow-expansion trio that brings the total
      to **124**, opening the first modern web 3D molecular
      visualization domain to ship in Valenx
- [x] 2 web-visualization templates in `valenx-init` (`molstar`,
      `ngl`), all round-tripping through `valenx-validate`
      (cross-binary roundtrip now sweeps **120 templates** clean
      alongside the Phase 22.5 workflow-expansion trio)
- [x] STATUS.md / README.md / ROADMAP.md / CHANGELOG.md /
      QUICKSTART.md updated
- [ ] **Phase 42.5** — sister-adapter expansion of Phase 42:
      3Dmol.js (sister WebGL viewer with a more compact API
      surface; the `py3Dmol` Python binding gives the same
      Python-script entry point — defer to sister-adapter
      expansion), LiteMol (Mol*'s direct predecessor; now in
      maintenance-only mode — defer until upstream activity
      resumes), PV (formerly the canonical lightweight WebGL
      viewer; now superseded by Mol* / NGL — defer), MolView
      (web-app rather than embeddable library; out of scope as
      a hosted service), Web3DMol (Chinese-academic alternative;
      defer), MoleculeKit (Acellera's PyTorch-friendly viewer;
      defer to a future phase). Out of scope for this beachhead.

## Success metrics

| Metric                                                | Target          |
|-------------------------------------------------------|-----------------|
| New web-visualization adapter (template + tests)      | 1 day per       |
| Render-Mol* → render-NGL loop across 2 tools          | < tool baseline |

## Leads into

Phase 42 opens the modern web 3D molecular visualization domain
that the user's bio / chemistry spec called out alongside the
Phase 17 / 23 desktop viewer beachheads, the Phase 17 / 17.5 /
17.7 prediction stack (every predicted structure deserves a
viewer), and the Phase 5.5 / 5.7 trajectory analysis stack (every
analyzed trajectory deserves a scrubbable visualization).
Combined with the existing run-Galaxy-workflow → run-WDL → run-
CWL → run-Nextflow → run-Snakemake → design-plasmid → view-
alignment → process-image → segment-cells → classify-pixels →
simulate-pathway → expand-rules → grow-tissue → diffuse-particles
→ trace-MCell-trajectories → simulate-MD → analyze-trajectory →
reweight-free-energy → fit-ENM → run-cpptraj-script → predict-
structure → fold-RNA → analyze-DNA-geometry → infer-tree-ML →
infer-tree-Bayesian → simulate-popgen → analyze-trees →
reconstruct-3D → design-protein → validate loop, the **render-
Mol* → render-NGL → run-Galaxy-workflow → run-WDL → run-CWL →
run-Nextflow → run-Snakemake → design-plasmid → view-alignment
→ process-image → segment-cells → classify-pixels → simulate-
pathway → expand-rules → grow-tissue → diffuse-particles →
trace-MCell-trajectories → simulate-MD → analyze-trajectory →
reweight-free-energy → fit-ENM → run-cpptraj-script → predict-
structure → fold-RNA → analyze-DNA-geometry → infer-tree-ML →
infer-tree-Bayesian → simulate-popgen → analyze-trees →
reconstruct-3D → design-protein → validate** loop now spans two
modern web 3D molecular visualization tools (Mol*, NGL Viewer)
feeding into every structural output Valenx's biology stack
emits — predicted structures from the Phase 17 / 17.5 / 17.7
prediction stack (ESMFold, OpenFold, AlphaFold 2/3, ColabFold,
RoseTTAFold, OmegaFold, FoldSeek), trajectories from the Phase
5 / 5.6 GROMACS / LAMMPS / NAMD / sander / HOOMD-blue MD engines
analyzed by the Phase 5.5 / 5.7 / 17 PLUMED / ProDy / cpptraj /
MDTraj / MDAnalysis post-MD analysis stack, designed proteins
from the Phase 27 / 27.5 / 27.6 design stack (RFdiffusion,
ProteinMPNN, Chroma, ESM-IF, RFantibody, ESM3, ESM Cambrian),
docked poses from the Phase 34 docking pair (AutoDock Vina,
AutoDock 4), reconstructed cryo-EM volumes from the Phase 36
reconstruction tools (RELION, EMAN2, CTFFIND), Rosetta-modeled
structures from the Phase 38 family (Rosetta, PyRosetta), and
the DNA structural-geometry outputs from the Phase 39 tools
(X3DNA, Curves+, DSSR) — all in one Valenx shell with no glue
code beyond the existing case-toml / prepare / run / collect
path.

**This completes the bio ecosystem from the original /review
list** — every major category the user originally called out
is now covered by a live adapter in Valenx. The bio surface
spans alignment, sequence editors, cheminformatics, cryo-EM,
CRISPR, DNA geometry, docking, MD analysis, MD engines,
microscopy, phylogenetics, population genetics, protein design,
quantum chemistry, RNA structure, sequence read simulators,
single-cell genomics, spatial stochastic simulation, structure
prediction, structure search, synthetic biology, systems
biology, variant calling, viewers (desktop + web), and workflow
managers — **105 bio adapters across 38 biology / biotech /
chemistry phases**, all in one Valenx shell with no glue code
beyond the existing case-toml / prepare / run / collect path.

The natural follow-up is **Phase 42.5** — the deferred web-
visualization work called out above (3Dmol.js / `py3Dmol` as a
sister WebGL viewer with a more compact API surface, LiteMol if
upstream activity resumes, PV if the lightweight-viewer niche
re-opens, MoleculeKit for the PyTorch-friendly visualization
shape), slotting in alongside the existing Mol* + NGL adapters
with the same Python-script subprocess shape (3Dmol.js / py3Dmol
sister tools) or web-service shape (MolView / Web3DMol if the
hosted-service pattern becomes acceptable for the registry
pattern). See the out-of-scope section of `docs/superpowers/
plans/2026-05-04-web-visualization.md` for the full follow-up
phase list.
