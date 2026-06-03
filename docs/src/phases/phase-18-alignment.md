# Phase 18 — Sequence alignment toolkit

**Status:** 🟢 Live — alignment beachhead landed.

## Goal

Round out Valenx's biology beachhead with the second-most-used
bio toolset after Biopython / OpenMM: sequence alignment + read
mapping. Phase 18 ships two new canonical types (`FastqRecord` /
`Alignment`) plus 6 first-class adapters wrapping the de-facto
open-source CLIs (BWA, minimap2, MAFFT, MUSCLE, HMMER, samtools).

## Capability inventory

### Canonical types (`valenx-bio`)
- `FastqRecord` — sequence + per-base quality scores. Validates
  that sequence and quality strings have matching lengths;
  serde round-trip; happy-path constructors for ASCII Phred+33
  quality.
- `Alignment` — multiple-sequence alignment as a list of named
  gapped sequences with shared length. Validates that every
  row has the same alignment length; preserves row insertion
  order.

### Format readers
- FASTQ reader + writer (4-line format, robust against trailing
  whitespace and `\r\n` line endings)
- Minimal SAM-text reader — parses header (`@HD` / `@SQ` / `@PG`)
  + records (QNAME / FLAG / RNAME / POS / CIGAR + …) sufficient
  for summary inspection. BAM (binary BGZF) deferred to Phase 18.5.

### Live adapters (6)
- **BWA** — short-read alignment via `bwa mem` / `bwa aln`. Stages
  the reference index reference; collect() classifies the produced
  SAM file as a typed artifact.
- **minimap2** — long-read + spliced + asm-vs-asm alignment.
  Selectable preset (`map-ont` / `map-pb` / `splice` / `asm5` /
  …) per case.
- **MAFFT** — multiple-sequence alignment. Wraps `mafft --auto`
  by default; algorithm overrides exposed via `[seqalign.mafft]`.
- **MUSCLE** — alternate MSA back-end. Wraps `muscle -align` /
  `-super5` for v5+; -in / -out flags for v3 fallback.
- **HMMER** — profile-HMM search. Covers `hmmbuild` / `hmmsearch`
  / `phmmer` / `jackhmmer` via a single dispatching adapter.
- **samtools** — SAM / BAM utilities. Initial coverage:
  `flagstat` / `view` / `sort` / `index` / `stats` against a
  user-provided file.

### Headless CLIs (2)
- `valenx-fastq` — inspect / validate FASTQ files. Subcommands +
  text/JSON output + stdin via `-`.
- `valenx-sam-info` — alignment summary (record count, mapped /
  unmapped tally, reference list, average MAPQ) from a SAM file.

## What landed early

The implementation landed across 12
discrete commits, each landing one canonical-type, format reader,
adapter, or CLI. Every commit kept workspace clippy + rustdoc
clean.

## Acceptance checklist

- [x] `valenx-bio` extends with `FastqRecord` + `Alignment` types
- [x] Minimal SAM-text reader in `valenx_bio::format::sam`
- [x] FASTQ reader + writer in `valenx_bio::format::fastq`
- [x] 6 alignment adapters wired into `valenx-app::init_registry`
- [x] 6 alignment templates in `valenx-init`, all round-tripping
      through `valenx-validate`
- [x] `valenx-fastq` + `valenx-sam-info` CLIs ship with tests
- [x] STATUS.md / README.md / ROADMAP.md / CHANGELOG.md /
      QUICKSTART.md updated
- [ ] **Phase 18.5** — BAM (BGZF) reader, BAM-aware
      `valenx-sam-info`, plus aligners that share BWA / minimap2's
      adapter shape (Bowtie2, HISAT2, STAR, MMseqs2, DIAMOND).
      Out of scope for this beachhead; next plan covers it.

## Success metrics

| Metric                                        | Target          |
|-----------------------------------------------|-----------------|
| `valenx-fastq inspect <bundled>` time         | < 100 ms        |
| `valenx-sam-info <bundled>` time              | < 100 ms        |
| New alignment adapter (template + tests)      | 1 day per       |

## Leads into

Phase 18.5 — BAM (binary BGZF) reading, plus Bowtie2 / HISAT2 /
STAR / MMseqs2 / DIAMOND adapters that share the BWA / minimap2
shape. Phase 19 picks up variant calling (GATK, bcftools,
DeepVariant, FreeBayes, Strelka2).
