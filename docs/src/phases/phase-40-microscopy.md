# Phase 40 — Microscopy

**Status:** 🟢 Live — Fiji + CellProfiler + Ilastik open the
**first microscopy / bioimage analysis domain** in Valenx
alongside the Phase 5.5 / 5.6 / 5.7 / 17 / 17.5 / 17.7 / 18 /
25 / 27 / 27.5 / 27.6 / 28 / 29 / 30 / 30.5 / 31 / 32 / 32.5 /
33 / 34 / 35 / 36 / 38 / 39 MD-engines + MD-analysis-expansion +
biology + structure-prediction + protein-design + RNA-structure
+ population-genetics + phylogenetics + Bayesian-phylogenetics +
read-simulator + systems-biology + spatial-stochastic +
synthetic-biology + docking + CRISPR-design + cryo-EM +
Rosetta-family + DNA-geometry beachheads and the Phase 24
cheminformatics expansion.

## Goal

Open the **microscopy / bioimage analysis** domain in Valenx
with three established open-source tools that span the bioimage
analysis tradeoff space — script-driven general-purpose image
processing in headless mode (Fiji, the canonical ImageJ
distribution that's the de-facto first stop in bioimage
analysis), pipeline-driven cell segmentation + measurement
(CellProfiler, the Broad Institute pipeline-driven workhorse
that powers most high-content-screening assays), and
interactive-ML pixel / object classification (Ilastik, the
Hamprecht lab tool that leans on user-trained random-forest
classifiers for hard segmentation tasks where rule-based
pipelines struggle). All three run in headless mode for batch
processing — Fiji + Ilastik via app-launcher binaries that ship
in the upstream distribution, CellProfiler via its Python CLI.
No new canonical types — image inputs and outputs are all
standard formats (TIFF, PNG, HDF5, CSV) that the existing
`Results.artifacts` collection model surfaces directly. Phase 40
sits numerically between Phase 39 DNA structural geometry and
Phase 41 sequence editors and ships chronologically right after
Phase 32.5 spatial stochastic — same chronological-vs-numerical
convention used for Phase 17.5 / 24 / 28 / 31 / 35 / 39 / 5.5 /
5.6 / 5.7 / 32.5.

## Capability inventory

### Live adapters (3)

- **Fiji** — the [Fiji Is Just ImageJ](https://fiji.sc/)
  distribution of NIH ImageJ (Schindelin et al, GPL-3.0). Fiji
  bundles ImageJ2 + a curated set of plugins for biological
  image processing — channel splitting, thresholding, particle
  analysis, deconvolution, registration, segmentation, the
  TrackMate single-particle tracker, the BoneJ trabecular-bone
  toolkit, and the entire ImageJ macro / Jython / scripting
  surface. Fiji ships per-platform launcher binaries
  (`ImageJ-linux64`, `ImageJ.exe`, `Contents/MacOS/ImageJ-
  macosx`); the `--headless --console -macro` invocation runs a
  user-authored `.ijm` macro file without GUI and prints to
  stdout. App-launcher subprocess shape (sister to Phase 36
  RELION / EMAN2): the user supplies the absolute path to the
  per-platform Fiji launcher via `[bio.fiji].fiji_app` in
  `case.toml`. The CLI is `<fiji_app> --headless --console
  -macro <macro_file> [extras...]`. Schema knobs: `fiji_app`
  (absolute path to the per-platform Fiji launcher; required),
  `macro_file` (`.ijm` Fiji macro file describing the image-
  processing pipeline; required), `input_image`
  (`Option<PathBuf>` — optional input image the macro will
  operate on, typically picked up via `getArgument()`; the macro
  is responsible for opening it), `output_basename` (filename
  stem the user's macro uses for outputs — surfaced here so
  collect() can label artifacts uniformly; required, non-empty),
  `extra_args` (additional CLI arguments appended after the
  macro path — typically the `getArgument()` payload Fiji passes
  to the macro). `prepare()` resolves `fiji_app`, `macro_file`,
  and the optional `input_image` against the case directory when
  relative, validates each existing file exists on disk (returns
  `InvalidCase` with a helpful message when missing), and
  composes the invocation. `collect()` walks the workdir for
  `<output_basename>*.tif` (`Native`, "Fiji image (TIFF)" — the
  canonical multi-page TIFF format Fiji writes for processed
  image stacks), `<output_basename>*.tiff` (`Native`, "Fiji
  image (TIFF)" — same kind, alternate extension), `<output_
  basename>*.png` (`Native`, "Fiji image (PNG)" — single-channel
  / RGB PNG output from `Save As → PNG` calls in the macro),
  `<output_basename>*.csv` (`Tabular`, "Fiji measurements" —
  per-particle / per-ROI measurement tables from the
  `Analyze → Measure` family), and `*.log` (`Log`). Probe via
  `find_on_path(&["ImageJ-linux64", "ImageJ-macosx",
  "ImageJ.exe", "fiji"])` — accepting any of the per-platform
  launcher names; surfaces a Java-fallback warning whenever
  Fiji isn't on PATH but `java` is, since the user can still
  launch Fiji via `java -jar` if they prefer to bypass the
  launcher script. Version range `2.0.0..3.0.0` (Fiji's
  ImageJ2-based 2.x line is the modern stable shipping the
  curated bundle; upper bound 3.0 reserves room for an eventual
  major bump). `bio.fiji.process` ribbon capability.
- **CellProfiler** — Broad Institute's pipeline-driven cell
  segmentation + measurement suite (BSD-3-Clause). CellProfiler
  is the canonical tool for high-content screening: the user
  authors a `.cppipe` pipeline in the GUI (a chain of modules
  like LoadImages → IdentifyPrimaryObjects → MeasureObjectShape
  → ExportToSpreadsheet) and the CLI runs that pipeline over a
  directory of input images, emitting per-object measurement
  CSVs + segmented label-image overlays. Python-CLI subprocess
  shape (sister to Phase 19.5 Scanpy on the import-with-warning
  side, but with a bundled `cellprofiler` launcher as the
  primary entry point and `<python> -m cellprofiler ...` as the
  fallback when the launcher isn't on PATH but Python is). The
  CLI is `cellprofiler -c -r -p <pipeline> -i <input_dir> -o
  <basename> [extras...]` (`-c` = run without GUI, `-r` = run
  pipeline immediately, `-p` = pipeline file, `-i` = input
  directory, `-o` = output directory). Schema knobs: `pipeline`
  (`.cppipe` / `.cpproj` pipeline file; required), `input_dir`
  (directory containing input images; required — the adapter
  validates it is a directory at prepare time), `output_basename`
  (filename stem the adapter pins as the `-o` output directory;
  required, non-empty), `python` (interpreter name; default
  `"python3"` — used for the `<python> -m cellprofiler ...`
  fallback when the launcher isn't on PATH), `extra_args`
  (additional CLI arguments appended after the `-o` flag).
  `prepare()` resolves `pipeline` and `input_dir` against the
  case directory when relative, validates `pipeline` exists on
  disk and `input_dir` is a directory (returns `InvalidCase`
  with helpful messages when missing or not-a-directory), looks
  up the `cellprofiler` binary on PATH first then falls back to
  `<python> -m cellprofiler ...`, and composes the invocation.
  `collect()` walks **one level deep** into the
  `<output_basename>/` subdirectory for `*.csv` (`Tabular`,
  "CellProfiler measurements" — the per-object / per-image
  measurement tables ExportToSpreadsheet writes), `*.tif` /
  `*.tiff` (`Native`, "CellProfiler segmented image" — label-
  image overlays SaveImages writes), and `*.png` (`Native`,
  "CellProfiler plot" — diagnostic plots from the
  DisplayDataOnImage / DisplayHistogram family); plus the
  workdir top-level for `*.log` (`Log`, "CellProfiler log"). The
  one-level-deep walk into `<output_basename>/` mirrors Phase
  36 EMAN2's `<basename>_NN/` pattern — output-directory-rooted
  walks let collect() pick up everything the pipeline writes
  without needing to know the per-module filename conventions.
  Probe via `find_on_path(&["cellprofiler", "python3",
  "python"])`; surfaces a warning when `cellprofiler` itself
  isn't on PATH but Python is ("CellProfiler not found on PATH;
  install via pip install cellprofiler or download from
  https://cellprofiler.org/releases"). Version range
  `4.0.0..5.0.0` (CellProfiler 4.x is the modern Python 3 line;
  the older 3.x Python 2 line is deprecated; upper bound 5.0
  reserves room for an eventual major bump).
  `bio.cellprofiler.segment` ribbon capability.
- **Ilastik** — Hamprecht lab's interactive-ML pixel / object
  classification suite (GPL-3.0). Ilastik leans on user-trained
  random-forest classifiers — the user paints a few foreground /
  background strokes per image in the GUI to teach the
  classifier, saves the resulting `.ilp` project file, and then
  runs the headless CLI to apply that trained classifier to a
  batch of new images. The canonical use case is hard
  segmentation tasks where rule-based pipelines (CellProfiler) or
  threshold-driven macros (Fiji) struggle — light-sheet imagery,
  tissue cross-sections, anything with low contrast or
  irregular textures. App-launcher subprocess shape (sister to
  Fiji): the user supplies the absolute path to the per-platform
  Ilastik launcher (`ilastik`, `run_ilastik.sh`, or
  `ilastik.exe`) via `[bio.ilastik].ilastik_app` in `case.toml`.
  The CLI is `<ilastik_app> --headless --project=<project>
  --output_filename_format=<basename>_{nickname}.h5
  <input_images...> [extras...]`. The `--project=<project>` and
  `--output_filename_format=<basename>_{nickname}.h5` flags are
  emitted as single OsString args each (so `=` and the value
  travel together, matching Ilastik's own argv parser). The
  literal `{nickname}` substring in the format string is
  Ilastik's per-image nickname placeholder — it must reach
  Ilastik unmodified for per-input-image output disambiguation.
  Schema knobs: `ilastik_app` (absolute path to the per-platform
  Ilastik launcher; required), `project` (`.ilp` Ilastik project
  file containing the trained classifier; required),
  `input_images` (`Vec<PathBuf>` — one or more input images;
  required, must contain ≥ 1 entry — the adapter rejects an
  empty vector at prepare time), `output_basename` (filename
  stem; required, non-empty), `workflow` (string; default
  `"Pixel Classification"` — selectable workflow name from
  Ilastik's set: `"Pixel Classification"`, `"Object
  Classification"`, etc.), `extra_args`. `prepare()` resolves
  `ilastik_app`, `project`, and each `input_images` entry
  against the case directory when relative, validates
  `input_images` is non-empty, and composes the invocation.
  `collect()` walks the workdir for `<output_basename>*.h5`
  (`Native`, "Ilastik probability map (HDF5)" — Ilastik's
  canonical per-pixel class-probability output in HDF5),
  `<output_basename>*.tif` (`Native`, "Ilastik segmentation" —
  alternate TIFF segmentation output for downstream pipelines
  that don't read HDF5), and `*.log` (`Log`). Probe via
  `find_on_path(&["ilastik", "run_ilastik.sh",
  "ilastik.exe"])`; surfaces a warning when nothing matches
  ("Ilastik not found on PATH; download from
  https://www.ilastik.org/download.html and add bin to PATH")
  but still returns `ok = true` since the user can supply the
  launcher path via `case.toml` even when nothing on PATH
  resolves at probe time (sister to Phase 32 PhysiCell's
  per-project-binary probe convention). Version range
  `1.4.0..2.0.0` (Ilastik 1.4 is the modern stable line shipping
  the contemporary headless mode + workflow set; upper bound
  2.0 reserves room for an eventual major bump).
  `bio.ilastik.classify` ribbon capability.

### Canonical types

**No new canonical types.** All three adapters consume user-
supplied inputs (Fiji `.ijm` macros + optional input images,
CellProfiler `.cppipe` pipelines + input image directories,
Ilastik `.ilp` project files + one-or-more input images) and
emit user-readable artifacts (Fiji TIFF / PNG processed image
stacks + CSV measurement tables, CellProfiler CSV per-object
measurement tables + TIFF label-image segmentation overlays +
PNG diagnostic plots, Ilastik HDF5 per-pixel class-probability
maps + TIFF segmentation overlays) that the unchanged
`Results.artifacts` collection model surfaces directly. A
first-class bioimage canonical type — a typed image-stack +
per-channel + per-frame representation spanning all three back-
ends, with parsed TIFF / HDF5 metadata and a typed segmentation-
mask representation — defers to a future phase along with
image-overlay viewers and per-segmentation-class histogram
plotters.

### Headless CLIs

**No new CLIs.** Fiji's TIFF / PNG / CSV outputs, CellProfiler's
CSV / TIFF / PNG outputs, and Ilastik's HDF5 / TIFF outputs are
all standard image / tabular formats inspectable in any image
viewer (Fiji itself, ImageJ, napari, Imaris, Bio-Formats) or
through the user's downstream Python pipeline (`scikit-image`,
`tifffile`, `h5py`, `pandas`, `numpy`). A canonical bioimage
analysis CLI — image-stack inspection, segmentation-mask
diffing, per-channel histogram comparison — defers to a future
phase along with the canonical type.

## Domain milestone

Phase 40 is the **first microscopy / bioimage analysis domain**
to land in Valenx. The biology adapter family started with Phase
17 (foundation — sequence / structure / trajectory canonical
types + classical MD + cheminformatics scripts) and expanded
through Phase 5.5 / 5.6 / 5.7 / 17.5 / 17.7 / 18 / 18.5 / 18.6
/ 18.7 / 19 / 19.5 / 19.6 / 20 / 22 / 23 / 24 / 25 / 27 / 27.5
/ 27.6 / 28 / 29 / 30 / 30.5 / 31 / 32 / 32.5 / 33 / 34 / 35 /
36 / 38 / 39 to cover MD-trajectory analysis expansion, bio MD
engines, MDTraj, sequence prediction, structure-tools expansion,
alignment, RNA-seq, alignment-toolkit expansion, variant
calling, single-cell, single-cell expansion, transcript
quantification, workflow orchestration, molecular viewers,
cheminformatics, quantum chemistry, protein design, protein-
design expansion, EvolutionaryScale models, RNA structure,
population genetics, phylogenetics, Bayesian phylogenetics,
sequencing read simulation, systems biology, spatial-stochastic
simulators, synthetic biology, small-molecule docking, CRISPR
design, cryo-EM reconstruction, Rosetta protein modeling, and
DNA structural geometry — but until Phase 40 the microscopy /
bioimage analysis surface (general-purpose script-driven image
processing, pipeline-driven cell segmentation + measurement,
interactive-ML pixel / object classification) was absent. Phase
40 closes that gap with three established open-source tools
spanning the bioimage analysis tradeoff space — Fiji at the
script-driven general-purpose end, CellProfiler for the
canonical pipeline-driven cell-segmentation surface, and
Ilastik for the interactive-ML hard-segmentation corner where
rule-based pipelines struggle.

## What landed early

The implementation rode subagent-driven-development across 4
discrete implementation commits (3 adapters plus the registry +
init-template rollup) plus this docs pass — each landing one
adapter, the registry rollup, the init-template extension, or
the documentation pass. Every commit kept workspace clippy +
rustdoc clean.

## Acceptance checklist

- [x] `valenx-adapter-fiji` adapter ships with case-input parser
      + 4 lib tests + 3 case-input tests covering parses-minimal
      / parses-with-input-image / rejects-empty-output-basename,
      plus the app-launcher subprocess shape that composes
      `<fiji_app> --headless --console -macro <macro_file>
      [extras...]` with `fiji_app` + `macro_file` + optional
      `input_image` resolved against the case directory and
      validated on disk
- [x] `valenx-adapter-cellprofiler` adapter ships with case-
      input parser + 4 lib tests + 3 case-input tests covering
      parses-minimal / rejects-non-directory-input / rejects-
      empty-output-basename, plus the Python-CLI subprocess
      shape that composes `cellprofiler -c -r -p <pipeline> -i
      <input_dir> -o <basename> [extras...]` with `<python> -m
      cellprofiler ...` fallback when the launcher isn't on
      PATH but Python is, and the one-level-deep walk into
      `<output_basename>/` for collected artifacts
- [x] `valenx-adapter-ilastik` adapter ships with case-input
      parser + 4 lib tests + 3 case-input tests covering parses-
      minimal / parses-with-multiple-images / rejects-empty-
      input-images, plus the app-launcher subprocess shape that
      composes `<ilastik_app> --headless --project=<project>
      --output_filename_format=<basename>_{nickname}.h5
      <input_images...> [extras...]` with the `--project=` and
      `--output_filename_format=` flags emitted as single
      OsString args each (so the `=` and value travel together)
      and the literal `{nickname}` substring preserved unmodified
      in the format string for Ilastik's per-image
      disambiguation
- [x] All 3 adapters wired into `valenx-app::init_registry` —
      live adapter count moves from 114 to 117 (alongside the
      Phase 32.5 spatial-stochastic pair and Phase 41
      sequence-editors pair that bring the total to **119**),
      opening the first microscopy / bioimage analysis domain
      to ship in Valenx
- [x] 3 microscopy templates in `valenx-init` (`fiji`,
      `cellprofiler`, `ilastik`), all round-tripping through
      `valenx-validate` (cross-binary roundtrip now sweeps
      **115 templates** clean alongside the Phase 32.5 spatial-
      stochastic pair and Phase 41 sequence-editors pair)
- [x] STATUS.md / README.md / ROADMAP.md / CHANGELOG.md /
      QUICKSTART.md updated
- [ ] **Phase 40.5** — sister-adapter expansion of Phase 40:
      napari (Python-script subprocess sister to Fiji on the
      Python side, fits the Phase 17 Biopython subprocess
      pattern; defer), QuPath (Java GUI app for whole-slide
      digital-pathology image analysis with a `qupath script`
      headless mode; defer), DeepCell / StarDist / Cellpose
      (deep-learning-based cell-segmentation tools sister to
      Ilastik on the ML side; defer to a future ML-bioimage
      expansion phase), Bioformats CLI (canonical bioimage
      format converter sister to Open Babel for chemistry; defer
      alongside the canonical-type work), Imaris (commercial
      with academic version + Python API; defer). Out of scope
      for this beachhead.

## Success metrics

| Metric                                                | Target          |
|-------------------------------------------------------|-----------------|
| New microscopy adapter (template + tests)             | 1 day per       |
| Process-image → segment-cells → classify-pixels loop across 3 tools | < tool baseline |

## Leads into

Phase 40 opens the microscopy / bioimage analysis domain that
the user's bio / chemistry spec called out alongside the Phase
17 / 17.5 / 17.7 biology + structure-prediction stack, the
Phase 32 / 32.5 systems-biology + spatial-stochastic surface,
and the Phase 36 cryo-EM reconstruction tools (RELION, EMAN2,
CTFFIND — also image-processing-adjacent but at the per-particle
electron-microscopy 3D-reconstruction end of the spectrum
rather than the cell-scale light-microscopy 2D / 3D segmentation
end Phase 40 covers). Combined with the existing simulate-
pathway → expand-rules → grow-tissue → diffuse-particles →
trace-MCell-trajectories → simulate-MD → analyze-trajectory →
reweight-free-energy → fit-ENM → run-cpptraj-script → predict-
structure → fold-RNA → analyze-DNA-geometry → infer-tree-ML →
infer-tree-Bayesian → simulate-popgen → analyze-trees →
reconstruct-3D → design-protein → validate loop, the **process-
image → segment-cells → classify-pixels → simulate-pathway →
expand-rules → grow-tissue → diffuse-particles → trace-MCell-
trajectories → simulate-MD → analyze-trajectory → reweight-
free-energy → fit-ENM → run-cpptraj-script → predict-structure
→ fold-RNA → analyze-DNA-geometry → infer-tree-ML → infer-tree-
Bayesian → simulate-popgen → analyze-trees → reconstruct-3D →
design-protein → validate** loop now spans three microscopy /
bioimage analysis tools (Fiji, CellProfiler, Ilastik) feeding
into the existing Phase 5 / 5.6 GROMACS / LAMMPS / NAMD /
sander / HOOMD-blue MD engines, the Phase 5.5 / 5.7 / 17 PLUMED
/ ProDy / cpptraj / MDTraj / MDAnalysis post-MD analysis stack,
the Phase 17 / 17.5 / 17.7 prediction stack (ESMFold, OpenFold,
AlphaFold 2/3, ColabFold, RoseTTAFold, OmegaFold, FoldSeek), the
Phase 28 RNA-structure tools (ViennaRNA, RNAstructure, NUPACK),
the Phase 29 population-genetics trio (SLiM, msprime, tskit),
the Phase 30 phylogenetic-tree builders (IQ-TREE, RAxML-NG,
FastTree), the Phase 30.5 Bayesian-phylogenetics pair (BEAST 2,
MrBayes), the Phase 32 systems-biology surface (COPASI,
BioNetGen, PhysiCell), the Phase 32.5 spatial-stochastic pair
(Smoldyn, MCell), the Phase 33 synthetic-biology trio (pySBOL,
j5, Cello), the Phase 34 docking pair (AutoDock Vina, AutoDock
4), the Phase 35 CRISPR-design tools (CHOPCHOP, CRISPOR, Cas-
OFFinder), the Phase 36 cryo-EM reconstruction tools (RELION,
EMAN2, CTFFIND), the Phase 38 Rosetta-family adapters (Rosetta,
PyRosetta), and the Phase 39 DNA-structural-geometry tools
(X3DNA, Curves+, DSSR) — all in one Valenx shell with no glue
code beyond the existing case-toml / prepare / run / collect
path.

The natural follow-up is **Phase 40.5** — the deferred
microscopy work called out above (napari as a Python-script
subprocess sister to Fiji, QuPath for whole-slide digital
pathology, DeepCell / StarDist / Cellpose as deep-learning-
based cell-segmentation tools sister to Ilastik on the ML side,
Bioformats CLI as the canonical bioimage format converter
alongside the canonical-type work, Imaris with its academic
version + Python API), slotting in alongside the existing Fiji
/ CellProfiler / Ilastik adapters with the same app-launcher
subprocess shape (Fiji / Ilastik / QuPath / Imaris sister tools),
Python-CLI subprocess shape (CellProfiler sister tools), or
Python-script subprocess shape (napari / DeepCell / StarDist /
Cellpose sister tools). See the out-of-scope section of
`docs/superpowers/plans/2026-05-03-microscopy.md` for the full
follow-up phase list.
