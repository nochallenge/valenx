# Phase 22 — Workflow managers

**Status:** 🟢 Live — Nextflow + Snakemake round out the bio
workflow-orchestration surface.

## Goal

Add the two de-facto bioinformatics workflow orchestrators to
Valenx. Phase 22 ships two sister adapters: **Nextflow** (the
DSL-driven pipeline language behind nf-core; `nextflow run
<pipeline>` shape) and **Snakemake** (the Python-flavoured
rule-based orchestrator; `snakemake -s <Snakefile> --cores N`
shape). Unlike the per-tool adapters Phase 17 / 18 / 19 / 23 / 24
ship (BWA, samtools, bcftools, …), these are **meta-tools** —
they invoke pipelines that themselves call other bio adapters'
underlying binaries. The Valenx adapter just orchestrates the
orchestrator, keeping the rest of the registry useful underneath.
Both follow the established Phase 18 BWA single-binary CLI shape:
probe / prepare / run / collect, output-in-workdir. The semantic
difference: where BWA's `out.sam` is the artifact, a workflow
manager's `work/` directory contains intermediate outputs from
arbitrarily many sub-pipelines — the adapter doesn't try to
introspect those, it just reports the workdir and lets the user
walk results. Phase 22 sits numerically before Phase 23 but ships
chronologically after Phase 24 — same convention as Phase 17.5.

## Capability inventory

### Live adapters (2)

- **Nextflow** — DSL-driven pipeline orchestrator behind nf-core
  and most modern bioinformatics workflows. Single-binary CLI
  shape: `nextflow run <pipeline> [-c <config>] [-profile
  <profile>] [-resume] [--<key> <value> for each param]
  [extras...]`. The `pipeline` field accepts a local `.nf`
  filename, a relative/absolute path, or a registry identifier
  like `nf-core/rnaseq`. `params` maps as `--<key> <value>` on
  the command line (values stringified; numeric / bool conversions
  happen in the Nextflow DSL). `profile` selects a config profile
  (`-profile <name>`); `resume` toggles `-resume` for incremental
  re-runs; optional `config` path passes through `-c <file>`.
  The native_command lives at the workdir's parent so Nextflow
  writes its `work/` and `.nextflow/` directories there.
  `collect()` reports the workdir as a `Native` artifact with
  label `"Nextflow run workdir"` and walks for `report.html` /
  `timeline.html` / `dag.svg` (Nextflow's standard observability
  outputs) surfacing them as `Log` artifacts. Apache-2.0 licensed.
  `bio.nextflow.run` ribbon capability.
- **Snakemake** — Python-flavoured rule-based pipeline orchestrator
  used heavily across academic genomics. Single-binary CLI shape:
  `snakemake -s <snakefile> --cores N [--use-conda] [-n]
  [--configfile <path>] [<targets>...] [extras...]`. `snakefile`
  points at the canonical `Snakefile` (default name; relative to
  the case dir or absolute); `targets` lists specific rules to
  build (empty = all default targets); `cores` (default 1, must be
  ≥ 1) sets `--cores N` for parallel-rule execution; `use_conda`
  toggles `--use-conda` for managed environments; `dry_run` toggles
  `-n` for plan-only inspection; optional `config_file` passes
  through `--configfile`. `collect()` reports the workdir as a
  `Native` artifact with label `"Snakemake run workdir"` and walks
  `.snakemake/log/*.log` if present, surfacing the most-recent log
  file as a `Log` artifact. MIT licensed. `bio.snakemake.run`
  ribbon capability.

### Canonical types

**No new canonical types.** Workflow managers are
meta-orchestrators — they don't produce a single canonical artifact
of their own. The pipelines they invoke produce whatever the
underlying tools do (BAM via BWA, VCF via bcftools, FASTA via
ColabFold, …), and the unchanged `Results.artifacts` collection
model surfaces them through their respective adapters' canonical
types. Phase 22's adapters report the workdir + observability logs
only; the per-tool typing happens in the adapters underneath.

### Headless CLIs

**No new CLIs.** Workflow-manager runs produce per-pipeline outputs
that Valenx's existing CLIs already cover (BAM via `valenx-sam-info`,
VCF via `valenx-vcf-info`, FASTA via `valenx-fasta`, PDB via
`valenx-pdb-info`, etc.). The workflow report HTML / DAG SVG that
Nextflow emits and the `.snakemake/log/` files Snakemake emits are
plain user-readable artifacts that don't need a Valenx-side
inspector.

## What landed early

The implementation landed across 5
discrete commits, each landing one adapter, the registry rollup,
the init-template extension, or the documentation pass. Every
commit kept workspace clippy + rustdoc clean.

## Acceptance checklist

- [x] `valenx-adapter-nextflow` adapter ships with case-input
      parser + 4 lib tests + 4 case-input tests
- [x] `valenx-adapter-snakemake` adapter ships with case-input
      parser + 4 lib tests + 4 case-input tests
- [x] Both adapters wired into `valenx-app::init_registry` —
      live adapter count moves from 47 to 49
- [x] 2 workflow-manager templates in `valenx-init` (`nextflow`
      with `nf` alias, `snakemake` with `smk` alias), both
      round-tripping through `valenx-validate` (cross-binary
      roundtrip now sweeps 45 templates clean)
- [x] STATUS.md / README.md / ROADMAP.md / CHANGELOG.md /
      QUICKSTART.md updated
- [ ] **Phase 22.5** — sister-orchestrator expansion: Galaxy
      (separate XML-tool-definition shape; doesn't fit the BWA
      pattern), Cromwell + WDL (Broad's pipeline language —
      similar shape to Nextflow but with a Java JVM dependency),
      CWL (Common Workflow Language; runs through `cwltool` /
      `toil` runners), nf-core pipeline-catalog browsing UI
      (Nextflow runs them via `nextflow run nf-core/<pipeline>`;
      no catalog UI in this beachhead), pipeline-graph
      visualization (DAG-rendering follow-up — Phase 23.5
      embedded-viewer territory). Out of scope for this beachhead.

## Success metrics

| Metric                                            | Target          |
|---------------------------------------------------|-----------------|
| New workflow-manager adapter (template + tests)   | 1 day per       |
| `nextflow run nf-core/<pipeline>` smoke loop      | < tool baseline |

## Leads into

Phase 22 closes the orchestration link in Valenx's bio toolchain:
where Phases 17 / 18 / 19 / 23 / 24 / 27 / 34 ship the per-tool
adapters that user-side `.nf` / `Snakefile` pipelines call into,
Phase 22's two adapters let users invoke an entire pipeline (and
its sub-tools) through the same case-toml / prepare / run / collect
shell as any single-tool adapter. Together with the underlying
BWA / samtools / bcftools / GATK / DeepVariant / RDKit / … chain,
Valenx now drives both the per-step and the whole-pipeline
patterns from one shell.

The natural follow-up is **Phase 22.5** — the deferred workflow
managers called out above (Galaxy, Cromwell + WDL, CWL via
`cwltool` / `toil`, nf-core catalog browsing, pipeline-graph
visualization).
