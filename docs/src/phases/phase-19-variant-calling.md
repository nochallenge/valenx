# Phase 19 — Variant calling toolkit

**Status:** 🟢 Live — variant calling beachhead landed.

## Goal

Take the Phase 18 alignment toolkit (BWA / minimap2 / samtools — already
shipped) and close the next link of the genomics workflow: variant
calling from aligned reads. Phase 19 ships a new `Vcf` / `VcfRecord`
canonical type plus a minimal VCF text reader, a headless
`valenx-vcf-info` CLI, and 3 first-class adapters wrapping the
conventional + ML-driven variant-calling stack (bcftools, GATK
HaplotypeCaller, DeepVariant).

## Capability inventory

### Canonical types (`valenx-bio`)
- `Vcf` — VCF v4.x file as parsed: `##` header lines verbatim,
  `samples` extracted from the `#CHROM` column header (empty when no
  per-sample columns), and a list of `VcfRecord`s. Serde round-trip;
  default-constructible to the empty file.
- `VcfRecord` — single variant row. `chrom` / `pos` / optional
  `id` / `ref_allele` / comma-split `alt` list / optional
  Phred-scaled `qual` / `;`-split `filter` list / raw `info` string /
  optional `format` plus per-sample column strings in the order of
  `Vcf.samples`. `is_pass()` recognises both `["PASS"]` and the
  `"."` (= unfiltered) convention; `has_alt()` rejects the
  `"."` ALT.

### Format readers
- Minimal VCF text reader (`valenx_bio::format::vcf::read_str`) —
  handles 8 mandatory columns + optional FORMAT + per-sample
  columns. Plain-text only; BCF (binary) and bgzf-compressed VCF
  are out of scope (convert with `bcftools view` first).

### Live adapters (3)
- **bcftools** — VCF/BCF multitool covering the four subcommands
  the variant-calling workflow most often touches: `view` (output
  filtering / format conversion), `call` (sensible-default
  multiallelic-caller variant calling on a BAM), `filter` (post-
  call filtering), and `concat` (joining per-region call outputs).
  Per-action dispatch validates that `input` / `inputs` /
  `reference` are populated only where the action requires them.
- **GATK** — Broad Institute reference variant caller. This phase
  focuses on the most-used tool (`HaplotypeCaller`) — joint
  genotyping / GVCF workflows are deferred to a follow-up. The
  adapter wraps `gatk --java-options "-Xmx<heap>" HaplotypeCaller`
  with reference + sorted-indexed BAM staging and an optional
  intervals (BED) restriction. Java heap validated against the
  conventional `8g` / `16g` style suffix.
- **DeepVariant** — Google's deep-learning variant caller. Wraps
  `run_deepvariant` with a typed `model_type` ∈ `{WGS, WES, PACBIO,
  ONT_R104, HYBRID_PACBIO_ILLUMINA}` and `num_shards` knob. Probe
  hint mentions both the direct binary and the Docker / Singularity
  wrapper paths — the adapter does not manage container runtimes,
  the user brings their own.

### Headless CLIs (1)
- `valenx-vcf-info` — VCF summary (header-line count, sample count,
  total records, PASS / FAIL split, no-ALT count) from a VCF file or
  stdin via `-`. Text + `--format json` modes mirror `valenx-sam-info`.

## What landed early

The implementation rode subagent-driven-development across 7 discrete
commits, each landing one canonical-type, format reader, adapter, CLI,
or registry / template wiring. Every commit kept workspace clippy +
rustdoc clean.

## Acceptance checklist

- [x] `valenx-bio` extends with `Vcf` + `VcfRecord` types
- [x] Minimal VCF text reader in `valenx_bio::format::vcf`
- [x] `valenx-vcf-info` CLI ships with tests
- [x] 3 variant-calling adapters wired into `valenx-app::init_registry`
- [x] 3 variant-calling templates in `valenx-init`, all round-tripping
      through `valenx-validate` (cross-binary roundtrip now sweeps 33
      templates clean)
- [x] STATUS.md / README.md / ROADMAP.md / CHANGELOG.md /
      QUICKSTART.md updated
- [ ] **Phase 19.5** — BCF (binary VCF) reader, plus the variant
      callers that share bcftools / GATK's adapter shape
      (Strelka2, FreeBayes, DELLY, Manta, Pindel, vcftools).
      Out of scope for this beachhead; next plan covers it.

## Success metrics

| Metric                                        | Target          |
|-----------------------------------------------|-----------------|
| `valenx-vcf-info <bundled>` time              | < 100 ms        |
| New variant-calling adapter (template + tests)| 1 day per       |

## Leads into

Phase 19.5 — BCF (binary VCF) reading, plus Strelka2 / FreeBayes /
DELLY / Manta / Pindel / vcftools adapters that share the
bcftools / GATK shape; joint genotyping / GVCF workflows; variant
annotation (SnpEff / VEP) follows in Phase 43. See the future-phases
table at the end of
`docs/superpowers/plans/2026-04-30-sequence-alignment-toolkit.md`
for the full follow-up phase list (Phases 19.5 → 43 cover the
remaining ~190 tools from the user's spec).
