# Phase 18.5 — Aligners expansion

**Status:** 🟢 Live — Bowtie2 + MMseqs2 + DIAMOND broaden the
alignment + protein-search surface alongside the Phase 18
BWA / minimap2 / MAFFT / MUSCLE / HMMER / samtools beachhead.

## Goal

Sister-adapter expansion of Phase 18. Add three more open-source
aligners covering distinct user-facing use cases: **Bowtie2**
(Langmead & Salzberg's gapped short-read aligner — alternative to
BWA for RNA-seq / ChIP-seq / bisulfite pipelines), **MMseqs2**
(Söding lab's "many vs. many" protein search + clustering toolkit —
fast alternative to BLAST and the prefilter behind ColabFold's
MSA generation), and **DIAMOND** (Buchfink, Reuter & Drost's
ultra-fast BLAST-protocol-compatible protein aligner — two to three
orders of magnitude faster than BLASTP / BLASTX for whole-metagenome
and UniRef-scale searches). All three follow the established
Phase 18 BWA pattern — single-binary CLI subprocess, file in / file
out. Bowtie2 mirrors BWA's two-stage `index → align` shape; MMseqs2
and DIAMOND dispatch per-action through the bcftools-style
`build_command(...) -> Result<Vec<OsString>, AdapterError>`
helper. No new infrastructure. Phase 18.5 sits numerically after
Phase 18 and ships chronologically right after Phase 27.5 — same
convention as Phase 17.5 sits between Phase 17 and Phase 18
numerically.

## Capability inventory

### Live adapters (3)

- **Bowtie2** — Langmead & Salzberg's gapped FM-index-based short-
  read aligner (GPL-3.0). Single-binary subprocess shape with a
  two-stage `bowtie2-build → bowtie2` pipeline that mirrors BWA's
  `bwa index → bwa mem` pattern. Schema knobs: `reference` (FASTA;
  required), `reads` (1 or 2 entries — single-end / paired-end
  FASTQ), `threads` (≥ 1, default 1), `skip_index` (default
  `false`; set `true` to reuse a pre-built FM-index), `preset`
  (default `"sensitive"`; whitelist
  `["very-fast", "fast", "sensitive", "very-sensitive"]` —
  Bowtie2's end-to-end preset family), `extra_args`. `prepare()`
  synchronously runs `bowtie2-build <reference> <reference>`
  unless `skip_index` is set, then composes
  `bowtie2 -x <ref_basename> --<preset> -p <threads> -S out.sam
  [-U <single-read> | -1 <r1> -2 <r2>] [extras...]`. `collect()`
  walks for `out.sam` (`Tabular`, `"Bowtie2 aligned reads"`) and
  any `.log` files (`Log`, `"Bowtie2 log"`). Probe via
  `find_on_path(&["bowtie2"])`. `bio.bowtie2.align` ribbon
  capability. The init alias `bt2` resolves to the same template.
- **MMseqs2** — Söding lab's high-throughput protein search +
  clustering toolkit (MIT). Single-binary subprocess; the user
  picks one of three high-level easy-* workflows via `action` in
  `[bio.mmseqs2]`. Schema knobs: `action` (whitelist
  `["easy-search", "easy-cluster", "easy-linsearch"]`), `query`
  (required), `target` (required for `easy-search` / `easy-
  linsearch`; ignored for `easy-cluster`), `output` (required),
  `sensitivity` (default 7.5 = max sensitivity; range
  `1.0..=7.5`, finite-checked), `threads` (≥ 1, default 1),
  `extra_args`. Per-action dispatch lives in
  `build_command(...) -> Result<Vec<OsString>, AdapterError>`
  (post-fix bcftools shape — `InvalidCase` on schema drift, never
  panics): `easy-search` →
  `mmseqs easy-search <query> <target> <output> tmp -s
  <sensitivity> --threads N [extras...]`; `easy-linsearch` →
  `mmseqs easy-linsearch <query> <target> <output> tmp --threads
  N [extras...]`; `easy-cluster` → `mmseqs easy-cluster <query>
  <output_prefix> tmp -s <sensitivity> --threads N [extras...]`.
  `collect()` reports the `output` path as `Tabular` with a
  per-action label (`"MMseqs2 easy-search hits"` /
  `"MMseqs2 easy-linsearch hits"` / `"MMseqs2 easy-cluster
  output"`). Probe via `find_on_path(&["mmseqs"])` — the on-disk
  binary is just `mmseqs` (no `2` suffix). MMseqs2 versions are
  git-hash-tagged (e.g. `14-7e284`); `version_range` spans
  `14.0.0..17.0.0` to cover the current major lines.
  `bio.mmseqs2.search` ribbon capability. The init alias `mmseqs`
  resolves to the same template.
- **DIAMOND** — Buchfink, Reuter & Drost's ultra-fast protein
  aligner (GPL-3.0). Single-binary subprocess; the user picks the
  mode via `action` in `[bio.diamond]`. Schema knobs: `action`
  (whitelist `["blastp", "blastx", "makedb"]`), `query`,
  `database`, `output` (all required), `sensitivity` (whitelist
  `["default", "fast", "sensitive", "more-sensitive",
  "very-sensitive", "ultra-sensitive"]`), `threads` (≥ 1),
  `extra_args`. Per-action dispatch in `build_command(...) ->
  Result<...>`: `blastp` / `blastx` → `diamond <action> -q
  <query> -d <database> -o <output> --<sensitivity> -p <threads>
  [extras...]` (the `--default` flag is omitted when
  `sensitivity = "default"` because DIAMOND's out-of-the-box
  default has no flag); `makedb` → `diamond makedb --in <query>
  -d <database> -p <threads> [extras...]`. In `makedb` mode the
  schema field roles flip — `query` is the input FASTA and
  `database` is the output DB basename (DIAMOND appends
  `.dmnd`); the adapter mirrors the upstream CLI directly so the
  schema names stay stable across actions. `collect()` reports
  the `output` for `blastp` / `blastx` (`Tabular`, `"DIAMOND
  <action> hits"`); for `makedb` reports `<database>.dmnd`
  (`Native`, `"DIAMOND .dmnd database"`). Probe via
  `find_on_path(&["diamond"])`. `bio.diamond.search` ribbon
  capability. The init alias `dmnd` resolves to the same
  template.

### Canonical types

**No new canonical types.** All three adapters consume the
existing Phase 17 FASTA and Phase 18 FASTQ inputs and emit
SAM (Bowtie2) or tabular hit-table (MMseqs2 / DIAMOND BLAST
format-8) outputs that the unchanged `Results.artifacts`
collection model surfaces directly.

### Headless CLIs

**No new CLIs.** Bowtie2's SAM output is already inspectable
through the Phase 18 `valenx-sam-info` CLI. MMseqs2 + DIAMOND
emit plain-text tabular hit tables; the existing `Tabular`
artifact kind surfaces them in `Results.artifacts` without
further tooling.

## What landed early

The implementation rode subagent-driven-development across 6
discrete commits, each landing one adapter, the registry rollup,
the init-template extension, or the documentation pass. Every
commit kept workspace clippy + rustdoc clean.

## Acceptance checklist

- [x] `valenx-adapter-bowtie2` adapter ships with case-input
      parser + 4 lib tests + 5 case-input tests
- [x] `valenx-adapter-mmseqs2` adapter ships with case-input
      parser + 4 lib tests + 5 case-input tests
- [x] `valenx-adapter-diamond` adapter ships with case-input
      parser + 4 lib tests + 5 case-input tests
- [x] All 3 adapters wired into `valenx-app::init_registry` —
      live adapter count moves from 54 to 57
- [x] 3 alignment templates in `valenx-init` (`bowtie2` with
      alias `bt2`, `mmseqs2` with alias `mmseqs`, `diamond` with
      alias `dmnd`), all round-tripping through `valenx-validate`
      (cross-binary roundtrip now sweeps 53 templates clean)
- [x] STATUS.md / README.md / ROADMAP.md / CHANGELOG.md /
      QUICKSTART.md updated
- [ ] **Phase 18.6** — RNA-seq-specific aligners: HISAT2 and STAR
      (different shape — splice-aware genome index + transcript
      annotation feed, deserves its own phase). LAST (niche
      pairwise aligner) defers further. BAM (binary) reader needs
      BGZF + same scope reason as Phase 18. Out of scope for this
      expansion.

## Success metrics

| Metric                                            | Target          |
|---------------------------------------------------|-----------------|
| New aligner adapter (template + tests)            | 1 day per       |
| Short-read align / protein-search loop across 3   | < tool baseline |

## Leads into

Phase 18.5 broadens the alignment + protein-search surface
alongside the Phase 18 BWA / minimap2 / MAFFT / MUSCLE / HMMER /
samtools beachhead. Combined with the Phase 17 ColabFold + Phase
17.5 ESMFold / OpenFold / AlphaFold 2 / AlphaFold 3 prediction
stack, the **search → align → predict → validate** loop now spans
nine alignment / search tools (BWA, Bowtie2, minimap2, MAFFT,
MUSCLE, HMMER, samtools, MMseqs2, DIAMOND) feeding into five
prediction tools — all in one Valenx shell with no glue code
beyond the existing case-toml / prepare / run / collect path.

The natural follow-up is **Phase 18.6** — the deferred RNA-seq
aligner work called out above (HISAT2, STAR), slotting in
alongside the existing alignment adapters with the same single-
binary-subprocess shape but with splice-aware index + transcript
annotation feed. See the future-phases table at the end of
`docs/superpowers/plans/2026-04-30-aligners-expansion.md` for the
full follow-up phase list.
