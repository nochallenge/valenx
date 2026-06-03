# Phase 22.5 — Workflow expansion

**Status:** 🟢 Live — planemo + Cromwell + cwltool round out the
**bio workflow-orchestration surface** that Phase 22 Nextflow +
Snakemake opened alongside the Phase 5.5 / 5.6 / 5.7 / 17 / 17.5 /
17.7 / 18 / 18.5 / 18.6 / 18.7 / 19 / 19.5 / 19.6 / 20 / 22 / 23 /
24 / 25 / 27 / 27.5 / 27.6 / 28 / 29 / 30 / 30.5 / 31 / 32 / 32.5 /
33 / 34 / 35 / 36 / 38 / 39 / 40 / 41 biology / biotech / chemistry
beachheads.

## Goal

Sister-adapter expansion of the existing Phase 22 Nextflow +
Snakemake workflow-manager pair. Round out the bio workflow-
orchestration surface with three more canonical workflow tools that
cover the corners Nextflow / Snakemake don't reach — Galaxy
ecosystem CLI for tool development + workflow execution outside a
full Galaxy server (planemo, the official Galaxy command-line
companion that lints tool wrappers, runs Galaxy workflow tests, and
executes `.ga` / `.gxwf.yml` workflows without requiring a Galaxy
server stand-up), Broad Institute WDL workflow engine (Cromwell,
the canonical Workflow Description Language runner powering most
production GATK + Terra pipelines), and the Common Workflow Language
reference runner (cwltool, the cross-tool standard implementation
for describing analytical workflows in YAML / JSON). planemo +
cwltool follow the established Phase 22 Nextflow / Snakemake single-
binary CLI pattern: workflow file in, run artifacts out. Cromwell
is JAR-distributed (no `cromwell` launcher binary on PATH); the
user supplies the absolute path to `cromwell-<version>.jar` via
case input, and we probe `java` itself rather than the JAR — same
JAR-distribution shape Phase 33 j5 / Cello and Phase 41 Jalview
use. Phase 22.5 sits numerically adjacent to Phase 22 and ships
chronologically right after Phase 41 sequence editors — same
chronological-vs-numerical convention used for Phase 17.5 / 24 /
28 / 31 / 35 / 39 / 5.5 / 5.6 / 5.7 / 32.5 / 40.

## Capability inventory

### Live adapters (3)

- **planemo** — the Galaxy project's official command-line
  companion for tool development + workflow execution outside a
  full Galaxy server (AFL-3.0). The same `planemo` binary lints
  tool wrappers, runs Galaxy workflow tests, and executes `.ga`
  workflow files / `.gxwf.yml` Galaxy-flavoured workflows; the
  `action` knob picks which sub-command to invoke. Single-binary
  subprocess shape (sister to Phase 22 Nextflow / Snakemake): the
  CLI is `planemo <action> <workflow> [inputs] [extras...]`.
  Schema knobs: `workflow` (`.ga` Galaxy workflow file or
  `.gxwf.yml` Galaxy-flavoured workflow; required), `inputs`
  (`Option<PathBuf>` — optional inputs JSON; `None` for workflows
  that take no inputs or have only defaulted inputs),
  `output_basename` (filename stem `collect()` uses to filter HTML
  reports; required, non-empty), `action` (string; default
  `"run"` — `prepare()` rejects values other than `run` / `test`
  / `lint` at parse time so the adapter doesn't forward
  unsupported sub-commands), `extra_args` (additional CLI
  arguments). `prepare()` resolves both `workflow` and the
  optional `inputs` against the case directory when relative,
  validates each file exists on disk (returns `InvalidCase` with
  a helpful message when missing), and composes the invocation.
  `collect()` walks the workdir for `<output_basename>*.html`
  (`Native`, "Planemo report" — Galaxy-style HTML run reports
  with per-step status), `*.json` (`Tabular`, "Planemo run JSON" —
  machine-readable workflow status + provenance), and `*.log`
  (`Log`, "Planemo log"). Probe via `find_on_path(&["planemo"])`.
  Version range `0.75.0..1.0.0` (planemo 0.75 (early 2023) is the
  floor where the modern `workflow_run` semantics and Galaxy
  23.0+ compatibility stabilised; upper bound 1.0 reserves room
  for the long-promised next major). `bio.planemo.run` ribbon
  capability.
- **Cromwell** — the Broad Institute's canonical Workflow
  Description Language (WDL) workflow engine (BSD-3-Clause).
  Cromwell powers most production GATK + Terra pipelines + a
  large fraction of academic genomics workflows, parsing the WDL
  language (a workflow DSL with explicit `task` + `workflow`
  blocks, typed inputs / outputs, and per-task Docker /
  Singularity backends) and dispatching tasks to local-shell /
  SLURM / SGE / Google Cloud / AWS Batch backends. **JAR-
  distributed** — no `cromwell` launcher binary on PATH; the user
  supplies the absolute path to `cromwell-<version>.jar` via
  `[bio.cromwell].jar` in `case.toml`. Single-binary subprocess
  shape (sister to Phase 33 j5 / Cello, Phase 41 Jalview) with
  `java -jar <jar> <action> <workflow> [-i <inputs>]
  [extras...]`. Schema knobs: `jar` (absolute path to
  `cromwell-<version>.jar`; required), `workflow` (`.wdl`
  workflow file; required), `inputs` (`Option<PathBuf>` —
  optional inputs JSON; emitted as **two separate args** `-i`
  + `<inputs>` only when `Some` — the flag is suppressed entirely
  when `None` rather than emitted with an empty value),
  `output_basename` (filename stem `collect()` uses to filter
  metadata JSON files; required, non-empty), `action` (string;
  default `"run"` — `prepare()` rejects values other than `run`
  / `submit` / `validate` at parse time), `extra_args`.
  `prepare()` resolves `jar`, `workflow`, and the optional
  `inputs` against the case directory when relative, validates
  each file exists on disk, and composes the `java -jar`
  invocation. `collect()` walks **the top level only** of the
  workdir for `<output_basename>*.json` (`Tabular`, "Cromwell
  metadata" — the per-run metadata JSON Cromwell writes
  alongside the workflow root) and `*.log` (`Log`, "Cromwell
  log"); per-task subdirectories under `cromwell-executions/`
  are out of scope for this adapter. Probe via
  `find_on_path(&["java"])` — Cromwell's version comes from the
  jar itself, not from `java`, so we surface no version here;
  the user pins the Cromwell release implicitly by the jar they
  point at (same shape as Phase 33 j5 / Cello and Phase 41
  Jalview). Probe attaches a warning that the real validation
  happens at prepare time. Version range `80.0.0..100.0.0`
  (Cromwell 80 (2023) is the modern stable line shipping the
  contemporary WDL 1.0 + WDL 1.1 grammar; upper bound 100
  reserves room for several majors of this actively maintained
  Broad project). `bio.cromwell.run` ribbon capability.
- **cwltool** — the reference implementation of the Common
  Workflow Language (Apache-2.0). CWL is the cross-tool standard
  for describing analytical workflows in YAML / JSON, with first-
  class support across Galaxy, Arvados, Toil, Cromwell, and
  cwltool itself; cwltool is the CWL specification's official
  Python reference implementation, used as the conformance test
  driver for every other CWL runner. Single-binary subprocess
  shape (sister to Phase 22 Snakemake) with `cwltool --outdir
  <output_dir> [extras...] <workflow> [inputs]`. cwltool itself
  stages tools (in-process or via Docker / Singularity / podman
  per the workflow's `DockerRequirement`) and writes the
  workflow's declared `outputs` into `<output_dir>/`. Schema
  knobs: `workflow` (`.cwl` tool / workflow document; required),
  `inputs` (`Option<PathBuf>` — optional CWL input-object
  document in JSON or YAML; `None` when the workflow takes no
  inputs or has only defaulted inputs; cwltool accepts either
  YAML or JSON for the input object), `output_dir` (`--outdir`
  target subdirectory under the workdir; required, non-empty),
  `extra_args`. `prepare()` resolves `workflow` and the optional
  `inputs` against the case directory when relative, validates
  each file exists on disk, and composes the invocation; the
  `--outdir <output_dir>` flag is workdir-relative — the
  subprocess runner's cwd is the workdir, so cwltool resolves
  the basename correctly and writes into `<workdir>/<output_dir>/`.
  Probe prefers the `cwltool` console-script entry-point and
  falls back to a Python-on-PATH detection — when only Python is
  reachable, the probe still returns `ok = true` with a targeted
  warning ("cwltool not found on PATH; install via `pip install
  cwltool`") so users with a Python environment ready but no
  `cwltool` package see the install hint without failing the
  probe (sister to the Phase 40 CellProfiler probe convention).
  `collect()` walks **one level deep** into `<output_dir>/`
  surfacing every file as a `Native` artifact labeled `"cwltool
  output"` (CWL's per-workflow output declarations are typed by
  the workflow itself, so the adapter intentionally doesn't try
  to second-guess the file taxonomy), plus the workdir-top-level
  `*.log` (`Log`, "cwltool log"). Probe via
  `find_on_path(&["cwltool"])`. Version range `3.1.0..4.0.0`
  (cwltool 3.1 (2023) is the modern stable line shipping the
  contemporary CWL 1.2 conformance + Docker / Singularity
  backends; upper bound 4.0 reserves room for the next major
  bump). `bio.cwltool.run` ribbon capability.

### Canonical types

**No new canonical types.** Workflow managers are meta-orchestrators
— they don't produce a single canonical artifact of their own. The
pipelines they invoke produce whatever the underlying tools do
(BAM via BWA, VCF via bcftools, FASTA via ColabFold, structures
via AlphaFold / RoseTTAFold, …), and the unchanged
`Results.artifacts` collection model surfaces them through their
respective adapters' canonical types. Phase 22.5's adapters report
the workflow-manager-level metadata (Galaxy reports + per-run
status JSON for planemo, Cromwell metadata JSON, CWL declared
outputs for cwltool) only; the per-tool typing happens in the
adapters underneath. A first-class workflow-manager canonical type
— a typed per-step status / per-task provenance representation
spanning all five back-ends (Phase 22 Nextflow + Snakemake, Phase
22.5 planemo + Cromwell + cwltool) — defers to a future phase
along with workflow-DAG visualizers and per-step status
inspectors.

### Headless CLIs

**No new CLIs.** planemo's Galaxy HTML reports + JSON status
files, Cromwell's metadata JSON + per-task logs, and cwltool's
declared outputs are all standard formats inspectable in any
browser, JSON viewer, or through the user's downstream genomics
pipeline. The per-pipeline outputs that all three orchestrators
produce are already covered by Valenx's existing CLIs (BAM via
`valenx-sam-info`, VCF via `valenx-vcf-info`, FASTA via
`valenx-fasta`, PDB via `valenx-pdb-info`). A canonical workflow-
manager CLI — per-step status diffing, per-task wall-clock
inspection, cross-orchestrator workflow comparison — defers to a
future phase along with the canonical type.

## Domain expansion

Phase 22.5 is a **sister-adapter expansion of the Phase 22
Nextflow + Snakemake workflow-manager pair** — the same bio
workflow-orchestration surface broadened with three more
established tools that cover the corners Nextflow / Snakemake
don't reach. Nextflow is the de-facto DSL-driven pipeline
language behind nf-core; Snakemake is the de-facto Python-
flavoured rule-based orchestrator; planemo is the official Galaxy
command-line companion that brings the Galaxy tool-development +
workflow-execution surface to a headless CLI; Cromwell is the
canonical Workflow Description Language runner powering most
production GATK + Terra pipelines; cwltool is the reference
implementation of the cross-tool Common Workflow Language
standard. With Phase 22.5 the bio workflow-orchestration surface
in Valenx covers all five canonical languages — Nextflow DSL
(Phase 22 Nextflow), Python-rule-based (Phase 22 Snakemake),
Galaxy `.ga` / `.gxwf.yml` (Phase 22.5 planemo), Workflow
Description Language `.wdl` (Phase 22.5 Cromwell), and Common
Workflow Language `.cwl` (Phase 22.5 cwltool).

## What landed early

The implementation rode subagent-driven-development across 5
discrete implementation commits (3 adapters, the registry rollup,
the init-template extension) plus this docs pass — each landing
one adapter, the registry rollup, the init-template extension, or
the documentation pass. Every commit kept workspace clippy +
rustdoc clean.

## Acceptance checklist

- [x] `valenx-adapter-planemo` adapter ships with case-input
      parser + 4 lib tests + 3 case-input tests covering parses-
      minimal / parses-with-inputs / rejects-bad-action, plus the
      single-binary subprocess shape that composes `planemo
      <action> <workflow> [inputs] [extras...]` with both files
      resolved against the case directory and validated on disk
      and `action ∈ {run, test, lint}` enforced at parse time
- [x] `valenx-adapter-cromwell` adapter ships with case-input
      parser + 4 lib tests + 3 case-input tests covering parses-
      minimal / parses-with-inputs / rejects-bad-action, plus the
      JAR-distributed single-binary subprocess shape that probes
      `java` on PATH and composes `java -jar <jar> <action>
      <workflow> [-i <inputs>] [extras...]` with all three paths
      resolved against the case directory and `action ∈ {run,
      submit, validate}` enforced at parse time, and `-i` /
      `<inputs>` emitted as TWO separate args only when `Some`
- [x] `valenx-adapter-cwltool` adapter ships with case-input
      parser + 4 lib tests + 3 case-input tests covering parses-
      minimal / parses-with-inputs / rejects-empty-output-dir,
      plus the single-binary subprocess shape that composes
      `cwltool --outdir <output_dir> [extras...] <workflow>
      [inputs]` with both files resolved against the case
      directory and validated on disk; probe prefers `cwltool`
      console-script and falls back to a Python-on-PATH detection
      with a "cwltool not found on PATH; install via `pip install
      cwltool`" warning
- [x] All 3 adapters wired into `valenx-app::init_registry` —
      live adapter count moves from 119 to **122** alongside the
      Phase 42 web-visualization pair that brings the total to
      **124**, rounding out the bio workflow-orchestration
      surface that Phase 22 Nextflow + Snakemake opened
- [x] 3 workflow-expansion templates in `valenx-init` (`planemo`,
      `cromwell`, `cwltool`), all round-tripping through
      `valenx-validate` (cross-binary roundtrip now sweeps **120
      templates** clean alongside the Phase 42 web-visualization
      pair)
- [x] STATUS.md / README.md / ROADMAP.md / CHANGELOG.md /
      QUICKSTART.md updated
- [ ] **Phase 22.6** — sister-adapter expansion of Phase 22.5:
      Toil (sister CWL / WDL runner with strong distributed
      backends; defer to a future workflow-manager phase),
      WorkflowHub (registry-style CWL workflow discovery; defer
      as a registry-not-runner shape), Tibanna (AWS-Lambda CWL
      runner; defer as a cloud-only shape), Arvados (CWL +
      Crunch; defer as a multi-component platform), Galaxy
      itself (web-server pipeline runner; out of scope as a web
      service rather than a single-binary CLI). Out of scope for
      this beachhead.

## Success metrics

| Metric                                                | Target          |
|-------------------------------------------------------|-----------------|
| New workflow-expansion adapter (template + tests)     | 1 day per       |
| run-Galaxy-workflow → run-WDL → run-CWL loop across 3 tools | < tool baseline |

## Leads into

Phase 22.5 rounds out the bio workflow-orchestration surface that
the user's bio / chemistry spec called out alongside the Phase 22
Nextflow + Snakemake beachhead. Combined with the existing
design-plasmid → view-alignment → process-image → segment-cells →
classify-pixels → simulate-pathway → expand-rules → grow-tissue →
diffuse-particles → trace-MCell-trajectories → simulate-MD →
analyze-trajectory → reweight-free-energy → fit-ENM → run-cpptraj-
script → predict-structure → fold-RNA → analyze-DNA-geometry →
infer-tree-ML → infer-tree-Bayesian → simulate-popgen → analyze-
trees → reconstruct-3D → design-protein → validate loop, the
**run-Galaxy-workflow → run-WDL → run-CWL → run-Nextflow → run-
Snakemake → design-plasmid → view-alignment → process-image →
segment-cells → classify-pixels → simulate-pathway → expand-rules
→ grow-tissue → diffuse-particles → trace-MCell-trajectories →
simulate-MD → analyze-trajectory → reweight-free-energy → fit-ENM
→ run-cpptraj-script → predict-structure → fold-RNA → analyze-
DNA-geometry → infer-tree-ML → infer-tree-Bayesian → simulate-
popgen → analyze-trees → reconstruct-3D → design-protein →
validate** loop now spans five workflow orchestrators (Phase 22
Nextflow + Snakemake plus Phase 22.5 planemo + Cromwell + cwltool)
that drive the existing per-tool adapters underneath — Phase 5 /
5.6 GROMACS / LAMMPS / NAMD / sander / HOOMD-blue MD engines, the
Phase 5.5 / 5.7 / 17 PLUMED / ProDy / cpptraj / MDTraj /
MDAnalysis post-MD analysis stack, the Phase 17 / 17.5 / 17.7
prediction stack (ESMFold, OpenFold, AlphaFold 2/3, ColabFold,
RoseTTAFold, OmegaFold, FoldSeek), the Phase 18 / 18.5 / 18.6 /
18.7 alignment + search surface (BWA, minimap2, MAFFT, MUSCLE,
HMMER, samtools, Bowtie2, MMseqs2, DIAMOND, HISAT2, STAR, BLAST+,
Clustal Omega, T-Coffee), the Phase 19 / 19.5 / 19.6 variant +
single-cell stack (bcftools, GATK, DeepVariant, Scanpy, scVI,
Seurat, AnnData), the Phase 20 transcript-quantification pair
(Salmon, Kallisto), the Phase 23 / 24 / 25 viewer + cheminformatics
+ quantum-chemistry tools (PyMOL, VMD, IGV, DeepChem, Open Babel,
Avogadro 2, Psi4, NWChem, xTB), the Phase 27 / 27.5 / 27.6 protein-
design stack (RFdiffusion, ProteinMPNN, Chroma, ESM-IF, RFantibody,
ESM3, ESM Cambrian), the Phase 28 RNA-structure tools (ViennaRNA,
RNAstructure, NUPACK), the Phase 29 / 30 / 30.5 population-
genetics + phylogenetics surface (SLiM, msprime, tskit, IQ-TREE,
RAxML-NG, FastTree, BEAST 2, MrBayes), the Phase 31 read-simulator
trio (ART, wgsim, Badread), the Phase 32 / 32.5 systems-biology +
spatial-stochastic surface (COPASI, BioNetGen, PhysiCell, Smoldyn,
MCell), the Phase 33 synthetic-biology trio (pySBOL, j5, Cello),
the Phase 34 docking pair (AutoDock Vina, AutoDock 4), the Phase
35 CRISPR-design tools (CHOPCHOP, CRISPOR, Cas-OFFinder), the
Phase 36 cryo-EM reconstruction tools (RELION, EMAN2, CTFFIND),
the Phase 38 Rosetta-family adapters (Rosetta, PyRosetta), the
Phase 39 DNA-structural-geometry tools (X3DNA, Curves+, DSSR), the
Phase 40 microscopy trio (Fiji, CellProfiler, Ilastik), and the
Phase 41 sequence-editor pair (pydna, Jalview) — all in one Valenx
shell with no glue code beyond the existing case-toml / prepare /
run / collect path.

The natural follow-up is **Phase 22.6** — the deferred workflow-
expansion work called out above (Toil as a sister CWL / WDL
runner with strong distributed backends, WorkflowHub for
registry-style CWL discovery, Tibanna / Arvados as cloud-shape
runners, Galaxy itself if the web-service shape becomes
acceptable for the registry pattern), slotting in alongside the
existing planemo / Cromwell / cwltool adapters with the same
single-binary CLI subprocess shape (Toil sister tools), JAR-
distributed shape (Arvados sister tools), or web-service shape
(Galaxy / WorkflowHub if the registry pattern shifts). See the
out-of-scope section of `docs/superpowers/plans/2026-05-04-
workflow-expansion.md` for the full follow-up phase list.
