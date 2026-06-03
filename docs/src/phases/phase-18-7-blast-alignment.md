# Phase 18.7 — Alignment toolkit expansion

**Status:** 🟢 Live — BLAST+ + Clustal Omega + T-Coffee round out
the **sequence-alignment surface** that Phase 18 BWA / minimap2 /
MAFFT / MUSCLE / HMMER / samtools opened, alongside the Phase 18.5
Bowtie2 / MMseqs2 / DIAMOND aligners-expansion trio and the Phase
18.6 HISAT2 / STAR RNA-seq aligners.

## Goal

Sister-adapter expansion of the existing Phase 18 / 18.5 / 18.6
alignment beachhead. Round out the foundational sequence-alignment
surface with three more established open-source tools that the
existing BWA / MAFFT / MUSCLE adapters explicitly left out —
canonical sequence-database search via NCBI BLAST+ (`blastn` /
`blastp` / `blastx` / `tblastn` / `tblastx`, the seminal nucleotide
/ protein search tool every wet-lab and bioinformatics pipeline
reaches for first), modern progressive multiple-sequence alignment
via Clustal Omega (Sievers / Higgins' HMM-driven successor to
ClustalW, the de-facto choice for routine MSA work alongside MAFFT
and MUSCLE), and consensus / library-based multiple-sequence
alignment via T-Coffee (Notredame / Higgins' library-based aligner
that combines pairwise alignments from many sources into a single
consistency-weighted MSA, the canonical choice for difficult
distantly-related sequences). All three follow the established
Phase 18 BWA single-binary CLI pattern: file in, alignment table /
search hits out. No new infrastructure. Phase 18.7 sits numerically
after Phase 18.6 and ships chronologically right after it — same
chronological-vs-numerical convention used for Phase 17.5 / 24 /
28 / 31 / 35 / 39.

## Capability inventory

### Live adapters (3)

- **BLAST+** — NCBI's seminal sequence-database search tool
  (Public Domain — US government work). BLAST+ ships five user-
  facing search programs covering every nucleotide / protein
  search direction: `blastn` (nucleotide query against nucleotide
  database), `blastp` (protein query against protein database),
  `blastx` (translated nucleotide query against protein database),
  `tblastn` (protein query against translated nucleotide
  database), and `tblastx` (translated nucleotide query against
  translated nucleotide database). Single-binary subprocess shape
  (sister to Phase 18 BWA): the per-program CLI is `<program>
  -query <query> -db <database> -out blast_results.txt -evalue
  <evalue> -outfmt <outfmt> -num_threads <threads> [extras...]`.
  Schema knobs: `program` (one of the five — required),
  `query` (FASTA query file; required), `database` (BLAST
  database path prefix — the user supplies the path stem and
  BLAST resolves the `<prefix>.nhr` / `<prefix>.phr` / etc.
  sidecars itself; required), `evalue` (`f64`, default 10.0 —
  BLAST's own default; the canonical "report a hit if its
  expected number under random chance is ≤ this value" cutoff),
  `outfmt` (`u8`, default 0 — BLAST's pairwise-alignment text
  output; 6 is the tab-separated format every downstream tool
  parses, 7 is the same with comment lines, 11 is the BLAST
  archive format), `threads` (`usize`, default 1), `extra_args`.
  `prepare()` resolves `query` against the case directory when
  relative, validates the database parent directory exists on
  disk (the database files themselves use the prefix
  convention so we cannot validate them by name), looks up the
  per-program binary via `find_on_path(&[&input.program])`, and
  composes the invocation with the output filename pinned to
  `blast_results.txt`. `collect()` walks the workdir for
  `blast_results.txt` (`Tabular`, "BLAST search results") and
  `*.log` (`Log`). Probe via `find_on_path(&["blastn", "blastp"])`
  — at least one BLAST+ binary on PATH counts as installed; the
  per-program lookup happens at `prepare()` based on the
  `program` field. Version range `2.10.0..3.0.0` (BLAST+ 2.10
  (2019) is the modern stable line that ships every NCBI release;
  upper bound 3.0 reserves room for an eventual major bump).
  `bio.blast.search` ribbon capability.
- **Clustal Omega** — Sievers / Higgins' modern HMM-driven
  multiple-sequence aligner (GPL-2.0). Clustal Omega is the
  modern successor to ClustalW: scales to thousands of
  sequences using HMM-based progressive alignment, the de-facto
  choice for routine MSA work alongside MAFFT and MUSCLE.
  Single-binary subprocess shape (sister to Phase 18 MAFFT):
  the CLI is `clustalo -i <input> -o <basename>.<ext>
  --outfmt=<outfmt> --threads=<N> [extras...]`. Schema knobs:
  `input` (FASTA multi-sequence input; required),
  `output_basename` (filename stem the user expects Clustal
  Omega to produce — e.g. `"alignment"` resolves to
  `"alignment.aln"`; required, non-empty), `outfmt` (default
  `"clustal"`; whitelist follows Clustal Omega's `--outfmt` set
  — `clustal` / `fasta` / `phylip` / `vienna` / `nexus`),
  `threads` (`usize`, default 1), `extra_args`. `prepare()`
  resolves `input` against the case directory when relative,
  validates it exists on disk, derives `<ext>` from `outfmt`
  (`clustal` → `.aln`, `fasta` → `.fasta`, `phylip` → `.phy`,
  `vienna` → `.vie`, `nexus` → `.nex`, default `.aln`), and
  composes the invocation. `collect()` walks the workdir for
  `<output_basename>*` (`Tabular`, "Clustal Omega alignment")
  and `*.log` (`Log`). Probe via `find_on_path(&["clustalo"])`.
  Version range `1.2.0..2.0.0` (Clustal Omega 1.2 is the modern
  stable release line; upper bound 2.0 reserves room for an
  eventual major bump). `bio.clustalo.align` ribbon capability.
- **T-Coffee** — Notredame / Higgins' library-based multiple-
  sequence aligner (GPL-2.0). T-Coffee combines pairwise
  alignments from many sources into a single consistency-
  weighted MSA, the canonical choice for difficult distantly-
  related sequences where progressive aligners (MAFFT, Clustal
  Omega, MUSCLE) lose accuracy. The library approach also
  supports specialised modes — `expresso` (structure-informed
  alignment via PDB lookups), `psicoffee` (PSI-BLAST profile-
  driven alignment), `mcoffee` (meta-aligner combining many
  back-ends) — selectable through the `mode` knob. Single-
  binary subprocess shape (sister to Clustal Omega): the CLI is
  `t_coffee <input> -output=<outfmt> -outfile=<basename>.aln
  [-mode=<mode>] [extras...]`. Schema knobs: `input` (FASTA
  multi-sequence input; required), `output_basename` (filename
  stem T-Coffee uses for outputs — output is always written as
  `<basename>.aln`; required, non-empty), `outfmt` (default
  `"clustalw"` — T-Coffee's own naming follows ClustalW
  conventions: `clustalw` / `fasta_aln` / `phylip` / `msf`),
  `mode` (`Option<String>` — omit to use T-Coffee's default
  progressive mode; set to `expresso` / `psicoffee` / `mcoffee`
  / etc. to opt into a specialised back-end), `extra_args`.
  Note T-Coffee's `=`-style flag form (`-output=<value>`,
  `-outfile=<value>`, `-mode=<value>`) instead of the more
  common space-separated form. `prepare()` resolves `input`
  against the case directory when relative, validates it exists
  on disk, and composes the invocation with the output always
  pinned to `.aln`. `collect()` walks the workdir for
  `<output_basename>*` (`Tabular`, "T-Coffee alignment"),
  `*.dnd` (`Native`, "T-Coffee guide tree" — T-Coffee writes
  the guide tree it constructs as a sibling Newick `.dnd` file
  for downstream phylogenetics consumption), and `*.log`
  (`Log`). Probe via `find_on_path(&["t_coffee"])` — note the
  underscore: T-Coffee installs as `t_coffee` (underscore, not
  hyphen), the project's own convention. Version range
  `13.0.0..14.0.0` (T-Coffee 13.x is the modern stable line;
  upper bound 14.0 reserves room for an eventual major bump).
  `bio.tcoffee.align` ribbon capability.

### Canonical types

**No new canonical types.** All three adapters consume the
existing Phase 17 FASTA inputs (BLAST+ takes a FASTA query plus
a BLAST-formatted database the user pre-builds with `makeblastdb`;
Clustal Omega and T-Coffee take a FASTA multi-sequence input)
and emit user-readable artifacts (BLAST+ tabular search results,
Clustal Omega alignment in the user-selected format, T-Coffee
alignment plus guide-tree Newick file) that the unchanged
`Results.artifacts` collection model surfaces directly. The
existing `valenx_bio::format::fasta` reader inspects collected
FASTA inputs for sequence counts; the existing
`valenx_bio::format::sam` reader covers the SAM outputs
neighbouring aligners produce. A first-class search-results
canonical type — typed BLAST hit tables with parsed query /
subject / e-value / bit-score / alignment-coordinate fields —
defers to a future phase along with MSA visualisation overlays
and per-residue conservation viewers.

### Headless CLIs

**No new CLIs.** BLAST+ tabular output (`outfmt = 6` / `7`) is a
plain `\t`-separated table inspectable in any text editor or
through the user's downstream Python pipeline (`pandas`,
`Bio.Blast.NCBIXML`); Clustal Omega and T-Coffee MSA outputs are
the standard ClustalW / FASTA / PHYLIP / NEXUS formats every
phylogenetics tool downstream of Phase 30 IQ-TREE / RAxML-NG /
FastTree consumes. Input FASTA files remain inspectable through
the existing Phase 17 `valenx-fasta` CLI; the existing Phase 17
`valenx-blast` CLI continues to wrap the auto-routing
`blastp` / `blastn` shorthand for users who don't need the full
adapter-driven pipeline. A canonical alignment-results CLI —
BLAST hit-table comparison, MSA conservation diffing, per-
residue alignment-quality inspection — defers to a future phase
along with the canonical type.

## Domain expansion

Phase 18.7 is a **sister-adapter expansion of the Phase 18 /
18.5 / 18.6 alignment beachhead** — the same sequence-alignment
surface broadened with three more established tools that cover
the corners the existing BWA / minimap2 / MAFFT / MUSCLE / HMMER
/ samtools / Bowtie2 / MMseqs2 / DIAMOND / HISAT2 / STAR set
doesn't reach. BLAST+ is the seminal database-search tool every
wet-lab and bioinformatics pipeline reaches for first (the
canonical pairwise / batch sequence-database search workflow
that MMseqs2 and DIAMOND accelerate but don't replace). Clustal
Omega is the modern HMM-driven progressive aligner that scales
to thousands of sequences alongside MAFFT and MUSCLE. T-Coffee
is the consistency-weighted library-based aligner that combines
pairwise alignments from many sources into a single MSA — the
canonical choice for difficult distantly-related sequences where
the progressive aligners lose accuracy. With Phase 18.7 the
alignment surface in Valenx covers all four canonical shapes —
short-read alignment (BWA, Bowtie2), long-read / spliced
alignment (minimap2, HISAT2, STAR), profile / database search
(BLAST+, HMMER, MMseqs2, DIAMOND), and multiple-sequence
alignment (MAFFT, MUSCLE, Clustal Omega, T-Coffee).

## What landed early

The implementation rode subagent-driven-development across 5
discrete implementation commits (3 adapters, the registry rollup,
the init-template extension) plus this docs pass — each landing
one adapter, the registry rollup, the init-template extension, or
the documentation pass. Every commit kept workspace clippy +
rustdoc clean.

## Acceptance checklist

- [x] `valenx-adapter-blast` adapter ships with case-input parser
      + 4 lib tests + 3 case-input tests, plus the single-binary
      subprocess shape that composes `<program> -query <query>
      -db <database> -out blast_results.txt -evalue <evalue>
      -outfmt <outfmt> -num_threads <threads> [extras...]` with
      `query` resolved against the case directory and the
      database parent directory validated on disk
- [x] `valenx-adapter-clustalo` adapter ships with case-input
      parser + 4 lib tests + 3 case-input tests, plus the
      single-binary subprocess shape that composes `clustalo -i
      <input> -o <basename>.<ext> --outfmt=<outfmt>
      --threads=<N> [extras...]` with `<ext>` mapped from
      `outfmt` (`clustal` → `.aln`, `fasta` → `.fasta`,
      `phylip` → `.phy`, `vienna` → `.vie`, `nexus` → `.nex`)
- [x] `valenx-adapter-tcoffee` adapter ships with case-input
      parser + 4 lib tests + 3 case-input tests, plus the
      single-binary subprocess shape that composes `t_coffee
      <input> -output=<outfmt> -outfile=<basename>.aln
      [-mode=<mode>] [extras...]` with the `=`-style flag form
      and the output always pinned to `.aln`
- [x] All 3 adapters wired into `valenx-app::init_registry` —
      live adapter count moves from 100 to **103** (alongside
      the Phase 19.6 single-cell-expansion pair that brings the
      total to **105**), rounding out the foundational
      sequence-alignment surface that Phase 18 / 18.5 / 18.6
      opened
- [x] 3 alignment-toolkit-expansion templates in `valenx-init`
      (`blast`, `clustalo`, `tcoffee`), all round-tripping
      through `valenx-validate` (cross-binary roundtrip now
      sweeps **101 templates** clean alongside the Phase 19.6
      single-cell-expansion pair)
- [x] STATUS.md / README.md / ROADMAP.md / CHANGELOG.md /
      QUICKSTART.md updated
- [ ] **Phase 18.8** — sister-adapter expansion of Phase 18.7:
      `makeblastdb` as a first-class adapter (currently the
      user pre-builds the BLAST database outside Valenx and
      points `database` at the prefix; defer to its own phase
      since the database-construction shape differs from the
      search shape), Kalign (lightweight progressive aligner
      sister to Clustal Omega; defer), ProbCons (alternative
      consistency-based aligner sister to T-Coffee; defer),
      MUSCLE 5 (the modern MUSCLE rewrite — the existing Phase
      18 MUSCLE adapter wraps the v3 / v5 fallback path; defer
      a v5-only adapter to a future phase if user demand
      surfaces), PRANK (phylogeny-aware aligner; defer), the
      remaining BLAST sister tools (`psiblast`, `rpsblast`,
      `deltablast`, `rpstblastn`, `legacy_blast.pl`; defer to
      18.8). Out of scope for this beachhead.

## Success metrics

| Metric                                                | Target          |
|-------------------------------------------------------|-----------------|
| New alignment-toolkit adapter (template + tests)      | 1 day per       |
| BLAST search + Clustal Omega MSA + T-Coffee MSA loop across 3 tools | < tool baseline |

## Leads into

Phase 18.7 rounds out the foundational sequence-alignment surface
that the user's bio / chemistry spec called out alongside the
Phase 18 / 18.5 / 18.6 alignment beachhead. Combined with the
existing predict-structure → fold-RNA → analyze-DNA-geometry →
infer-tree-ML → infer-tree-Bayesian → simulate-popgen → analyze-
trees → simulate-pathway → reconstruct-3D → design-protein →
validate loop, the **search → align → predict → validate** loop
now spans fourteen alignment / search tools (BWA, Bowtie2,
HISAT2, STAR, minimap2, MAFFT, MUSCLE, HMMER, samtools, MMseqs2,
DIAMOND, BLAST+, Clustal Omega, T-Coffee) feeding into the five
Phase 17 + 17.5 prediction tools (ColabFold, ESMFold, OpenFold,
AlphaFold 2, AlphaFold 3) — all in one Valenx shell with no glue
code beyond the existing case-toml / prepare / run / collect
path.

The natural follow-up is **Phase 18.8** — the deferred
alignment-toolkit work called out above (`makeblastdb` as the
first-class database-construction adapter, Kalign / ProbCons as
sister progressive / consistency-based MSA tools, the remaining
BLAST sister tools `psiblast` / `rpsblast` / `deltablast` /
`rpstblastn`), slotting in alongside the existing BLAST+ /
Clustal Omega / T-Coffee adapters with the same single-binary
subprocess shape. See the out-of-scope section of
`docs/superpowers/plans/2026-05-02-blast-alignment.md` for the
full follow-up phase list.
