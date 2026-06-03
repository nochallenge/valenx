# Phase 18.6 — RNA-seq alignment

**Status:** 🟢 Live — HISAT2 + STAR add the two de-facto
splice-aware RNA-seq aligners alongside the Phase 18 BWA /
minimap2 / MAFFT / MUSCLE / HMMER / samtools and Phase 18.5
Bowtie2 / MMseqs2 / DIAMOND alignment beachhead.

## Goal

Sister-adapter expansion of Phase 18 / 18.5. Add the two de-facto
RNA-seq aligners to Valenx: **HISAT2** (Daehwan Kim's graph-based
splice-aware aligner — successor to TopHat) and **STAR**
(Alex Dobin's most-used spliced aligner, the reference RNA-seq
mapper backing GTEx / TCGA / ENCODE pipelines and the only Phase
18.x aligner that doubles as a chromatin-conformation tool).
Both are spliced extensions of Phase 18's BWA / Phase 18.5's
Bowtie2 — they handle reads that span exon-exon junctions, where
the linear short-read aligners would soft-clip or misalign. Both
adapters mirror the established Phase 18 BWA two-stage shape —
single-binary CLI subprocess, file in / file out, `index → align`
pipeline. STAR has a heavier index step (genomic + splice-junction
database) but the same overall shape. No new infrastructure.
Phase 18.6 sits numerically after Phase 18.5 and ships
chronologically right after it — same convention as Phase 17.5
sits between Phase 17 and Phase 18 numerically.

## Capability inventory

### Live adapters (2)

- **HISAT2** — Daehwan Kim's graph-based splice-aware RNA-seq
  aligner (GPL-3.0). Single-binary subprocess shape with a
  two-stage `hisat2-build → hisat2` pipeline that mirrors BWA's
  `bwa index → bwa mem` pattern and Bowtie2's `bowtie2-build →
  bowtie2` pattern. Schema knobs: `reference` (FASTA; required),
  `reads` (1 or 2 entries — single-end / paired-end FASTQ),
  `threads` (≥ 1, default 1), `skip_index` (default `false`;
  set `true` to reuse a pre-built HFM index), `strandness`
  (default `"unstranded"`; whitelist
  `["unstranded", "F", "R", "FR", "RF"]` — F/R variants match
  Illumina TruSeq stranded library prep conventions),
  `extra_args`. `prepare()` synchronously runs
  `hisat2-build <reference> <reference_basename>` unless
  `skip_index` is set, then composes `hisat2 -x <ref_basename>
  -p <threads> -S out.sam [--rna-strandness <strandness>]
  [-U <single-read> | -1 <r1> -2 <r2>] [extras...]`. The
  `--rna-strandness` flag is omitted when
  `strandness = "unstranded"` because HISAT2 treats unstranded
  data as the default. `collect()` walks for `out.sam`
  (`Tabular`, `"HISAT2 aligned reads"`). Probe via
  `find_on_path(&["hisat2"])`. `bio.hisat2.align` ribbon
  capability. The init alias `hisat` resolves to the same
  template.
- **STAR** — Alex Dobin's spliced RNA-seq aligner (MIT). Note
  the capitalized binary name — `find_on_path(&["STAR"])`, not
  `star`. Single-binary subprocess shape with a two-stage
  `--runMode genomeGenerate → --runMode alignReads` pipeline.
  STAR's index step is heavier than BWA / Bowtie2 / HISAT2 — it
  builds a suffix-array-indexed genome under `genome_dir/` and
  optionally a splice-junction database from a GTF — but the
  adapter shape is the same. Schema knobs: `genome_dir` (the
  pre-built STAR index directory, or where the adapter writes
  one if `skip_index = false`; required), `reference` (FASTA;
  required only when generating the index), `reads` (1 or 2
  entries), `threads` (≥ 1, default 1), `skip_index` (default
  `false`), `output_type` (default `"BAM_SortedByCoordinate"`;
  whitelist `["BAM_Unsorted", "BAM_SortedByCoordinate", "SAM"]`),
  `sjdb_gtf` (optional GTF for splice-junction-database-aware
  indexing), `extra_args`. The `output_type` underscore-
  delimited canonical names map to STAR's two-arg
  `--outSAMtype` form: `"BAM_Unsorted"` →
  `--outSAMtype BAM Unsorted`, `"BAM_SortedByCoordinate"` →
  `--outSAMtype BAM SortedByCoordinate`, `"SAM"` →
  `--outSAMtype SAM`. `prepare()` synchronously runs
  `STAR --runMode genomeGenerate --genomeDir <genome_dir>
  --genomeFastaFiles <reference> --runThreadN N
  [--sjdbGTFfile <sjdb_gtf>] [--sjdbOverhang 100]` unless
  `skip_index` is set, then composes `STAR --runMode alignReads
  --genomeDir <genome_dir> --readFilesIn <reads...>
  --runThreadN N --outSAMtype <output_type spec>
  --outFileNamePrefix star_ [extras...]`. `collect()` walks
  for `star_Aligned.out.{bam,sam}` (`Tabular` for SAM,
  `Native` for BAM, `"STAR aligned reads"`) and
  `star_Log.final.out` (`Log`, `"STAR alignment summary"`).
  Validation: when `skip_index == false`, `reference` is
  required. `bio.star.align` ribbon capability.

### Canonical types

**No new canonical types.** Both adapters consume the existing
Phase 17 FASTA and Phase 18 FASTQ inputs and emit SAM / BAM
outputs that the unchanged `Results.artifacts` collection model
surfaces directly.

### Headless CLIs

**No new CLIs.** HISAT2's SAM output is already inspectable
through the Phase 18 `valenx-sam-info` CLI. STAR's BAM output
needs the existing samtools adapter (`samtools view`) to convert
to SAM before `valenx-sam-info` can read it; STAR's
`star_Log.final.out` is plain text and surfaces directly through
the `Log` artifact kind.

## What landed early

The implementation landed across 5
discrete commits, each landing one adapter, the registry rollup,
the init-template extension, or the documentation pass. Every
commit kept workspace clippy + rustdoc clean.

## Acceptance checklist

- [x] `valenx-adapter-hisat2` adapter ships with case-input
      parser + 4 lib tests + 5 case-input tests
- [x] `valenx-adapter-star` adapter ships with case-input
      parser + 4 lib tests + 5 case-input tests
- [x] Both adapters wired into `valenx-app::init_registry` —
      live adapter count moves from 57 to 59
- [x] 2 RNA-seq alignment templates in `valenx-init` (`hisat2`
      with alias `hisat`, `star`), all round-tripping through
      `valenx-validate` (cross-binary roundtrip now sweeps 55
      templates clean)
- [x] STATUS.md / README.md / ROADMAP.md / CHANGELOG.md /
      QUICKSTART.md updated
- [ ] **Phase 20** — transcript-quantification (Salmon, Kallisto)
      — different shape (k-mer-based pseudoalignment, no genome
      index), deserves its own phase. TopHat (deprecated; HISAT2
      is the successor) skipped. Cufflinks (assembly, not
      alignment) defers to Phase 20.5 if user demand surfaces.
      Out of scope for this expansion.

## Success metrics

| Metric                                            | Target          |
|---------------------------------------------------|-----------------|
| New aligner adapter (template + tests)            | 1 day per       |
| RNA-seq splice-aware align across 2 tools         | < tool baseline |

## Leads into

Phase 18.6 closes the RNA-seq alignment gap that Phase 18 and
Phase 18.5 explicitly deferred. Combined with the Phase 17
ColabFold + Phase 17.5 ESMFold / OpenFold / AlphaFold 2 /
AlphaFold 3 prediction stack and the Phase 18 / 18.5 alignment
beachhead, the **search → align → predict → validate** loop now
spans eleven alignment / search tools (BWA, Bowtie2, HISAT2,
STAR, minimap2, MAFFT, MUSCLE, HMMER, samtools, MMseqs2, DIAMOND)
feeding into five prediction tools — all in one Valenx shell with
no glue code beyond the existing case-toml / prepare / run /
collect path.

The natural follow-up is **Phase 20** — transcript-quantification
adapters (Salmon, Kallisto) — sitting downstream of the Phase 18.6
RNA-seq aligners with a different shape (k-mer-based
pseudoalignment, no genome index).
