# Phase 31 — Sequencing read simulators

**Status:** 🟢 Live — ART + wgsim + Badread open the **first
sequencing read-simulation domain** in Valenx alongside the Phase
17 / 17.5 / 18 / 18.5 / 18.6 / 19 / 19.5 / 20 / 25 / 27 / 27.5 /
27.6 / 28 / 30 / 32 / 34 / 36 biology + structure-prediction +
alignment + variant-calling + transcript-quantification +
protein-design + RNA-structure + phylogenetics + systems-biology
+ docking + cryo-EM beachheads.

## Goal

Open the sequencing read-simulation domain in Valenx with three
established open-source tools that span all three major
sequencing-technology classes — **ART** (the de-facto Illumina
short-read simulator with per-platform empirical error profiles),
**wgsim** (the classic short-read simulator that ships bundled with
samtools, deliberately simple under a uniform error model), and
**Badread** (the Nanopore long-read simulator with realistic per-
platform error profiles including chimeras, adapters, and identity
drift). All three follow the established Phase 18 BWA single-binary
CLI pattern: reference FASTA in, simulated FASTQ(s) out. ART writes
its output via the `-o` prefix flag, wgsim takes positional output
arguments after the reference, and Badread writes to stdout (handled
via the MAFFT-style stdout-redirect-to-file pattern). Phase 31 sits
numerically before Phase 32 but ships chronologically right after
Phase 36 cryo-EM — same chronological-vs-numerical convention used
for Phase 17.5 / 24 / 28.

## Capability inventory

### Live adapters (3)

- **ART** — Weichun Huang's NIEHS Illumina-platform read simulator
  (GPL-3.0). The de-facto choice for synthesising FASTQs that match
  the empirical error profile of a given Illumina sequencing system
  (HiSeq 2500, HiSeq X, MiSeq v3, NextSeq 500, MiniSeq) so
  downstream pipelines can be validated against a known-truth
  reference at controlled coverage and read length. Single-binary
  subprocess shape: the adapter wraps `art_illumina`, the workhorse
  of the ART family (companion `art_454` / `art_SOLiD` binaries
  cover platforms this adapter does not surface). Schema knobs:
  `reference` (FASTA; required), `output_prefix` (filename stem;
  ART writes `<prefix>.fq` for single-end or `<prefix>1.fq` +
  `<prefix>2.fq` for paired-end; required, non-empty),
  `sequencing_system` (one of `"HS25"` / `"HSXt"` / `"MSv3"` /
  `"NS50"` / `"MinS"`; required), `read_length` (≥ 1; required),
  `fold_coverage` (> 0.0; required), `paired_end` (default `false`),
  `fragment_mean` (mean insert size for paired-end; default 200.0,
  > 0.0 when `paired_end`), `fragment_sd` (insert-size stddev for
  paired-end; default 10.0, > 0.0 when `paired_end`), `extra_args`.
  `prepare()` builds `art_illumina -ss <sequencing_system> -i
  <reference> -l <read_length> -f <fold_coverage> -o
  <output_prefix> [-p -m <fragment_mean> -s <fragment_sd> if
  paired_end] [extras...]`. `collect()` walks the workdir top-level
  for `<output_prefix>*.fq` (`Tabular`, "ART simulated reads") and
  `<output_prefix>*.aln` (`Log`, "ART alignment record" — the per-
  read alignment record ART writes alongside the FASTQ, useful for
  validating aligner accuracy against the simulated truth). Probe
  via `find_on_path(&["art_illumina"])`. Version range
  `2.5.0..3.0.0` (the long-running ChocolateCherryCake `2.5.x`
  series since 2016; Bioconda + Homebrew ship 2.5.8). The
  `valenx-init` template ships with the alias `art-illumina`
  alongside the canonical `art`. `bio.art.simulate` ribbon
  capability.
- **wgsim** — Heng Li's classic Whole-Genome SIMulator that ships
  alongside samtools (MIT). Always paired-end, always position-
  uniform, deliberately simple under a uniform sequencing-error
  model with configurable insert size, read length, and per-base
  error rate. Unlike ART (which models per-platform empirical error
  profiles), wgsim is the canonical "small + classic" simulator for
  fast smoke-testing of mappers and variant callers when realistic
  error spectra are not required. Single-binary subprocess shape:
  `wgsim` takes the reference and both output FASTQs as positional
  arguments (no stdout-redirect needed). Schema knobs: `reference`
  (FASTA; required), `output1` (FASTQ for read 1; required, non-
  empty), `output2` (FASTQ for read 2; required, non-empty —
  wgsim is paired-end only), `num_pairs` (≥ 1; required), `length1`
  (read 1 length, default 70, ≥ 1), `length2` (read 2 length,
  default 70, ≥ 1), `fragment_size` (outer fragment length, default
  500, > 0), `error_rate` (per-base error rate in `0.0..=1.0`,
  default 0.02 — typical Illumina baseline), `extra_args`.
  `prepare()` builds `wgsim -N <num_pairs> -1 <length1> -2
  <length2> -d <fragment_size> -e <error_rate> <reference>
  <output1> <output2> [extras...]`. `collect()` reports `output1`
  and `output2` as `Tabular` artifacts ("wgsim simulated reads").
  Probe via `find_on_path(&["wgsim"])`. Version range `1.0.0..2.0.0`
  (wgsim is versioned alongside the parent samtools 1.x line; the
  binary historically prints the matching samtools tag on startup).
  `bio.wgsim.simulate` ribbon capability.
- **Badread** — Ryan Wick's long-read simulator with realistic
  Nanopore (and PacBio CLR) error profiles (GPL-3.0). Badread's
  per-platform error models are calibrated against actual sequencer
  output: random / chimeric / adapter / glitch read types, junk-read
  injection, identity drift, and length distributions that match
  what users see from a live flowcell. The de-facto choice for
  stress-testing long-read pipelines under realistic conditions.
  Single-binary subprocess shape with stdout-redirect: Badread
  writes its simulated FASTQ to stdout (no `-o` flag), so `run()`
  borrows MAFFT's stdout-redirect-to-file pattern — spawn the
  child directly, attach stdout to a `File` via `Stdio::from(file)`,
  stream stderr through the line handler. Schema knobs: `reference`
  (FASTA; required), `output` (FASTQ output path; required, non-
  empty), `quantity` (Badread's `--quantity` literal — one or more
  decimal digits followed by an optional `K` / `M` / `G` / `T` SI
  suffix, e.g. `"100M"` for 100 megabases or `"5G"` for 5 gigabases;
  validated against the per-platform `is_valid_quantity` helper),
  `error_model` (one of `"nanopore2018"` / `"nanopore2020"` /
  `"nanopore2023"` / `"pacbio2016"`; required — selects the per-
  platform error profile baked into the Badread distribution),
  `identity_mean` (read identity mean as a percentage in
  `0.0..=100.0`; default 87.5), `length_mean` (read length mean in
  bases; default 15000.0, > 0.0), `length_sd` (read length stddev
  in bases; default 13000.0, > 0.0), `extra_args`. `prepare()`
  builds `badread simulate --reference <reference> --quantity
  <quantity> --error_model <error_model> --identity <identity_mean>
  --length <length_mean>,<length_sd> [extras...]` → stdout, captured
  to `output` via the MAFFT-style stdout-redirect pattern.
  `collect()` reports `output` as a single `Tabular` artifact
  ("Badread simulated reads"). Probe via `find_on_path(&["badread"])`.
  Version range `0.4.0..1.0.0` (the long-running 0.4.x stable
  series; a 1.0 cut hasn't happened yet but the upper bound reserves
  room for it). `bio.badread.simulate` ribbon capability.

### Canonical types

**No new canonical types.** All three adapters consume user-supplied
reference FASTAs (the existing `valenx_bio::format::fasta` reader
already inspects sequence count + identifiers + alphabets) and emit
FASTQ files that the existing Phase 18 `valenx-fastq` CLI inspects
for record count, base quality distributions, and read-length
statistics. The unchanged `Results.artifacts` collection model
surfaces every emitted FASTQ + ART alignment record directly. A
first-class read-simulation provenance type — recording which
simulator produced which FASTQ under which error model — defers to
a future phase along with simulator-aware pipeline stitching.

### Headless CLIs

**No new CLIs.** ART's `<prefix>*.fq` simulated reads, wgsim's
paired `output1` / `output2` FASTQs, and Badread's stdout-captured
FASTQ are all standard four-line FASTQ records that the Phase 18
`valenx-fastq` CLI already inspects for record count, mean quality,
length distributions, and per-base validity. ART's `<prefix>*.aln`
per-read alignment record is a tabular text file inspectable in
any editor or through the user's downstream Python pipeline
(`pandas`, `numpy`). A canonical read-simulator CLI defers to a
future phase along with simulator-comparison and ground-truth-
diffing integrations.

## Domain milestone

Phase 31 is the **first sequencing read-simulation domain** to land
in Valenx. The biology adapter family started with Phase 17
(foundation — sequence / structure / trajectory canonical types +
classical MD + cheminformatics scripts) and expanded through Phase
17.5 / 18 / 18.5 / 18.6 / 19 / 19.5 / 20 / 22 / 23 / 24 / 25 / 27 /
27.5 / 27.6 / 28 / 30 / 32 / 34 / 36 to cover sequence prediction,
alignment, RNA-seq, variant calling, single-cell, transcript
quantification, workflow orchestration, molecular viewers,
cheminformatics, quantum chemistry, protein design,
EvolutionaryScale models, RNA structure, phylogenetics, systems
biology, small-molecule docking, and cryo-EM reconstruction — but
until Phase 31 the read-simulation surface (synthetic FASTQ
generation across the three major sequencing-technology classes)
was absent. Phase 31 closes that gap with three established open-
source tools spanning the full simulator tradeoff space — ART for
empirical-error-profile Illumina short reads, wgsim for the simple-
uniform-error short-read baseline, and Badread for realistic-error-
profile Nanopore long reads.

## What landed early

The implementation rode subagent-driven-development across 5
discrete implementation commits (3 adapters, the registry rollup,
the init-template extension) plus this docs pass — each landing one
adapter, the registry rollup, the init-template extension, or the
documentation pass. Every commit kept workspace clippy + rustdoc
clean.

## Acceptance checklist

- [x] `valenx-adapter-art` adapter ships with case-input parser
      + 4 lib tests + 5 case-input tests covering parses-minimal /
      parses-paired-end / rejects-bad-sequencing-system / rejects-
      zero-read-length / rejects-zero-fold-coverage, plus the
      paired-end dispatch shape that adds `-p -m <mean> -s <sd>`
      to the `art_illumina` invocation when `paired_end = true`
- [x] `valenx-adapter-wgsim` adapter ships with case-input parser
      + 4 lib tests + 4 case-input tests covering the
      `-N <num_pairs> -1 <length1> -2 <length2> -d <fragment_size>
      -e <error_rate>` knob shape (positional output arguments
      after the reference, no stdout-redirect)
- [x] `valenx-adapter-badread` adapter ships with case-input parser
      + 4 lib tests + 5 case-input tests covering parses-minimal /
      parses-with-extras / rejects-bad-quantity / rejects-bad-error-
      model / rejects-zero-length, plus the stdout-redirect-to-file
      shape that mirrors MAFFT's custom `run()` (Badread writes its
      simulated FASTQ to stdout rather than via `-o`)
- [x] All 3 adapters wired into `valenx-app::init_registry` —
      live adapter count moves from 78 to **81**, opening the first
      sequencing read-simulation domain to ship in Valenx
- [x] 3 read-simulator templates in `valenx-init` (`art` with alias
      `art-illumina`, `wgsim`, `badread`), all round-tripping
      through `valenx-validate` (cross-binary roundtrip now sweeps
      **77 templates** clean)
- [x] STATUS.md / README.md / ROADMAP.md / CHANGELOG.md /
      QUICKSTART.md updated
- [ ] **Phase 31.5** — DWGSIM (Nils Homer's wgsim fork with
      structural-variant injection; sister to wgsim; defer to
      sister-adapter expansion phase), pIRS (BGI's profile-based
      Illumina simulator; sister to ART; defer), InSilicoSeq
      (HMM-based Illumina + ONT simulator; defer), Mason (UCSC's
      single-binary read simulator covering Illumina + 454; defer),
      CuReSim (PCR / amplicon-aware simulator; defer), pbsim2 /
      pbsim3 (PacBio HiFi simulator, sister to Badread; defer to
      31.5), NanoSim (Nanopore simulator with model training;
      different shape — requires a per-flowcell training step; defer
      to a future phase). Out of scope for this beachhead.

## Success metrics

| Metric                                                | Target          |
|-------------------------------------------------------|-----------------|
| New read-simulator adapter (template + tests)         | 1 day per       |
| Simulate-Illumina + simulate-classic + simulate-Nanopore loop across 3 tools | < tool baseline |

## Leads into

Phase 31 opens the sequencing read-simulation domain that the user's
bio / chemistry spec called out alongside the Phase 18 / 18.5 /
18.6 alignment beachhead and the Phase 19 / 20 variant-calling +
transcript-quantification stack. Combined with the existing
simulate → align → quantify → call-variants → validate loop, the
**simulate-reads → align → quantify → call-variants → predict-
structure → fold-RNA → infer-tree → simulate-pathway → reconstruct-
3D → validate** loop now spans three read simulators (ART, wgsim,
Badread) feeding into the eleven Phase 18 / 18.5 / 18.6 alignment
tools (BWA, Bowtie2, HISAT2, STAR, minimap2, MAFFT, MUSCLE, HMMER,
samtools, MMseqs2, DIAMOND), the two Phase 20 transcript quantifiers
(Salmon, Kallisto), the three Phase 19 variant callers (bcftools,
GATK, DeepVariant), and through them into the existing Phase 17 /
17.5 prediction stack, the Phase 28 RNA-structure tools, the Phase
30 phylogenetic-tree builders, the Phase 32 systems-biology surface,
and the Phase 36 cryo-EM reconstruction tools — all in one Valenx
shell with no glue code beyond the existing case-toml / prepare /
run / collect path.

The natural follow-up is **Phase 31.5** — the deferred read-simulator
work called out above (DWGSIM as a wgsim fork with structural-variant
injection, pIRS as a profile-based Illumina alternative to ART,
InSilicoSeq as an HMM-based Illumina + ONT simulator, Mason for the
single-binary Illumina + 454 surface, CuReSim for PCR / amplicon-
aware simulation, pbsim2 / pbsim3 for the PacBio HiFi sister to
Badread), slotting in alongside the existing read simulators with
the same single-binary subprocess shape. NanoSim (Nanopore simulator
with model training) sits in a separate phase — the per-flowcell
training step is a different shape than the prepare / run / collect
path the Phase 31 simulators ride. See the out-of-scope section of
`docs/superpowers/plans/2026-04-30-read-simulators.md` for the full
follow-up phase list.
