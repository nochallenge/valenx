# Phase 20 — Transcript quantification

**Status:** 🟢 Live — Salmon + Kallisto add the two de-facto
transcript-level quantification tools alongside the Phase 18 BWA /
minimap2 / MAFFT / MUSCLE / HMMER / samtools, Phase 18.5 Bowtie2 /
MMseqs2 / DIAMOND, and Phase 18.6 HISAT2 / STAR alignment beachhead.

## Goal

Sister-domain expansion of Phase 18 / 18.5 / 18.6. Add the two
de-facto transcript-level quantification tools to Valenx: **Salmon**
(Rob Patro's quasi-mapping plus two-phase EM transcript-level
quantification — the GTEx / TCGA / nf-core / GENCODE reference
quantifier) and **Kallisto** (Lior Pachter's pseudoalignment-based
quantifier — the original "skip the alignment" approach). Both
pseudo-align reads to a transcriptome and report TPM / count per
transcript without producing intermediate SAM / BAM files; they are
shape-distinct from the Phase 18.6 RNA-seq aligners (HISAT2, STAR)
because they emit per-transcript abundance tables rather than aligned
reads. Both adapters mirror the established Phase 18 BWA two-stage
shape — single-binary CLI subprocess, file in / file out, `index →
quant` pipeline. No new infrastructure. Phase 20 sits numerically
after Phase 19.5 and ships chronologically right after the Phase 18.6
RNA-seq alignment beachhead — same convention as Phase 18.6 sits
between Phase 18.5 and Phase 19 numerically.

## Capability inventory

### Live adapters (2)

- **Salmon** — Rob Patro's quasi-mapping plus two-phase EM
  transcript-level quantification tool (GPL-3.0). Single-binary
  subprocess shape with a two-stage `salmon index → salmon quant`
  pipeline that mirrors BWA's `bwa index → bwa mem` pattern, Bowtie2's
  `bowtie2-build → bowtie2` pattern, and HISAT2's `hisat2-build →
  hisat2` pattern. Schema knobs: `transcriptome` (FASTA; required —
  used to build the index when `skip_index = false`), `index_dir`
  (the index directory salmon writes to / reads from; required),
  `reads` (1 or 2 entries — single-end / paired-end FASTQ),
  `output_dir` (`salmon quant -o`; required), `threads` (≥ 1, default
  1), `skip_index` (default `false`; set `true` to reuse a pre-built
  salmon index), `libtype` (default `"A"`; library-type DSL — `"A"`
  auto-detects orientation, `"U"` unstranded, `"ISF"` / `"ISR"`
  paired-end stranded forward / reverse, `"IU"` paired-end
  unstranded — left non-whitelisted because Salmon's libtype DSL has
  many valid combos), `extra_args`. `prepare()` synchronously runs
  `salmon index -t <transcriptome> -i <index_dir> -p <threads>`
  unless `skip_index` is set, then composes the quant command:
  single-end: `salmon quant -i <index_dir> -l <libtype> -p <threads>
  -o <output_dir> -r <reads[0]> [extras...]`; paired-end: `salmon
  quant -i <index_dir> -l <libtype> -p <threads> -o <output_dir> -1
  <reads[0]> -2 <reads[1]> [extras...]`. `collect()` walks
  `<output_dir>` for `quant.sf` (`Tabular`, `"Salmon transcript
  quantification"`) and `cmd_info.json` (`Log`, `"Salmon command
  info"`). Probe via `find_on_path(&["salmon"])`. `bio.salmon.quant`
  ribbon capability.
- **Kallisto** — Lior Pachter's pseudoalignment-based transcript
  quantifier (BSD-2-Clause). Single-binary subprocess shape with a
  two-stage `kallisto index → kallisto quant` pipeline. Kallisto's
  index is a single `.idx` file (not a directory) — the only shape
  difference from Salmon. Schema knobs: `transcriptome` (FASTA;
  required), `index` (single `.idx` file path — kallisto convention),
  `reads` (1 or 2 entries), `output_dir` (required), `threads` (≥ 1,
  default 1), `skip_index` (default `false`), `fragment_length`
  (optional `f64` — required for single-end reads only; `kallisto
  quant -l`), `fragment_sd` (optional `f64` — required for single-end
  reads only; `kallisto quant -s`), `extra_args`. Validation: when
  `reads.len() == 1`, both `fragment_length` and `fragment_sd` must
  be present, finite, and `> 0.0` (kallisto auto-detects fragment
  statistics from paired-end reads but cannot for single-end).
  `prepare()` synchronously runs `kallisto index -i <index>
  <transcriptome>` unless `skip_index` is set, then composes the
  quant command: paired-end: `kallisto quant -i <index> -o
  <output_dir> -t <threads> <reads[0]> <reads[1]> [extras...]`;
  single-end: `kallisto quant -i <index> -o <output_dir> -t <threads>
  --single -l <fragment_length> -s <fragment_sd> <reads[0]>
  [extras...]`. `collect()` walks `<output_dir>` for `abundance.tsv`
  (`Tabular`, `"Kallisto transcript abundance"`), `abundance.h5`
  (`Native`, `"Kallisto HDF5 abundance"`), and `run_info.json`
  (`Log`, `"Kallisto run info"`). Probe via
  `find_on_path(&["kallisto"])`. `bio.kallisto.quant` ribbon
  capability.

### Canonical types

**No new canonical types.** Both adapters consume the existing
Phase 17 FASTA and Phase 18 FASTQ inputs and emit per-transcript
abundance tables (Salmon `quant.sf`, Kallisto `abundance.tsv`)
that the unchanged `Results.artifacts` collection model surfaces
directly through the `Tabular` artifact kind.

### Headless CLIs

**No new CLIs.** Salmon's `quant.sf` and Kallisto's `abundance.tsv`
are tab-separated text and surface directly through the `Tabular`
artifact kind; no dedicated inspector CLI is needed because the
schema is human-readable. Kallisto's `abundance.h5` HDF5 sidecar
needs an external tool (`h5dump` / Python `h5py`) to inspect — the
canonical H5 reader as a Valenx CLI defers to Phase 19.6 along with
the Seurat / AnnData R-runtime work.

## What landed early

The implementation rode subagent-driven-development across 5
discrete commits, each landing one adapter, the registry rollup,
the init-template extension, or the documentation pass. Every
commit kept workspace clippy + rustdoc clean.

## Acceptance checklist

- [x] `valenx-adapter-salmon` adapter ships with case-input
      parser + 4 lib tests + 5 case-input tests
- [x] `valenx-adapter-kallisto` adapter ships with case-input
      parser + 4 lib tests + 5 case-input tests
- [x] Both adapters wired into `valenx-app::init_registry` —
      live adapter count moves from 59 to 61
- [x] 2 transcript quantification templates in `valenx-init`
      (`salmon`, `kallisto`), all round-tripping through
      `valenx-validate` (cross-binary roundtrip now sweeps 57
      templates clean)
- [x] STATUS.md / README.md / ROADMAP.md / CHANGELOG.md /
      QUICKSTART.md updated
- [ ] **Phase 19.6** — AnnData / `.h5ad` HDF5 reader as a canonical
      type (needs the `hdf5` crate, a non-trivial C-library dep)
      and Seurat R-runtime adapter pattern. Out of scope for this
      expansion. StringTie / Cufflinks (transcript assembly —
      different workflow shape) defers to Phase 20.5 if user demand
      surfaces. Tximport / DESeq2 / edgeR (downstream differential
      expression — R-runtime territory) slot into Phase 19.6.

## Success metrics

| Metric                                            | Target          |
|---------------------------------------------------|-----------------|
| New quantifier adapter (template + tests)         | 1 day per       |
| Transcript-level quant across 2 tools             | < tool baseline |

## Leads into

Phase 20 closes the transcript-level quantification gap that Phase
18.6 explicitly deferred ("Salmon / Kallisto — transcript
quantification, not alignment — different shape, k-mer-based
pseudoalignment, no genome index; defer to Phase 20"). Combined with
the Phase 17 ColabFold + Phase 17.5 ESMFold / OpenFold / AlphaFold 2
/ AlphaFold 3 prediction stack and the Phase 18 / 18.5 / 18.6
alignment beachhead, the **search → align → quantify → predict →
validate** loop now spans eleven alignment / search tools (BWA,
Bowtie2, HISAT2, STAR, minimap2, MAFFT, MUSCLE, HMMER, samtools,
MMseqs2, DIAMOND) plus two transcript quantifiers (Salmon, Kallisto)
feeding into five prediction tools — all in one Valenx shell with no
glue code beyond the existing case-toml / prepare / run / collect
path.

The natural follow-up is **Phase 19.6** — AnnData / `.h5ad` canonical
type plus the Seurat R-runtime adapter pattern (Tximport / DESeq2 /
edgeR slot into the same R-runtime infrastructure). See the
future-phases section of
`docs/superpowers/plans/2026-04-30-transcript-quantification.md` for
the full follow-up phase list.
