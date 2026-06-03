# valenx-e2e-tests

End-to-end integration tests that spawn real upstream tools through the
Valenx adapter trait. Every test wires `prepare()` → `run()` →
`collect()` against a freshly installed binary (BWA, samtools, MAFFT,
ViennaRNA, bcftools, HMMER, minimap2, MUSCLE, Bowtie2, MMseqs2,
DIAMOND, HISAT2, BLAST+, Clustal Omega, T-Coffee, IQ-TREE, RAxML-NG,
FastTree, MrBayes, RNAstructure, ART, wgsim, Badread, Snakemake,
Nextflow, cwltool, Smoldyn, Open Babel, xTB, Cas-OFFinder, ...) and
verifies the adapter's subprocess machinery — argument layout, working
directory, stdout / stderr handling, artifact discovery — survives
contact with the real tool. Unit tests cover schema parsing and trait
wiring; this crate closes the "does it actually run?" gap.

Tests that find their upstream binary missing on `PATH` print a
`Skipping E2E test — ...` line and exit cleanly. That keeps the crate
useful on developer machines without 30+ bio tools installed locally.

The crate currently covers **30** adapters end-to-end.

## Running locally

Install just the tools you need for the adapter you're testing, then
run that test target:

```bash
# All 30 current tools via conda-forge / bioconda:
conda install -c bioconda -c conda-forge \
    bwa samtools mafft viennarna bcftools hmmer \
    minimap2 muscle bowtie2 mmseqs2 diamond hisat2 \
    blast clustalo t-coffee iqtree raxml-ng fasttree \
    mrbayes rnastructure art wgsim badread snakemake \
    nextflow cwltool smoldyn openbabel xtb cas-offinder

# Or a single one:
conda install -c bioconda bwa

# Run one test target:
cargo test -p valenx-e2e-tests --test bwa

# Run the whole crate (skips anything not installed):
cargo test -p valenx-e2e-tests -- --nocapture
```

## Running in CI

`.github/workflows/ci-nightly.yml` installs every adapter's tool via
conda and runs `cargo test -p valenx-e2e-tests` on Ubuntu every night
at 04:00 UTC. Manual runs are available via the **Actions → CI nightly
(E2E) → Run workflow** button on GitHub.

## Adding a new adapter test

1. Add the adapter as a `[dev-dependencies]` entry in
   `crates/valenx-e2e-tests/Cargo.toml`:

   ```toml
   valenx-adapter-yourtool = { path = "../valenx-adapters/<domain>/valenx-adapter-yourtool" }
   ```

2. Add a new `tests/<adapter>.rs` following the template
   (`tests/bwa.rs` is the canonical example).

3. Build a minimal real `case.toml` plus any input files inside the
   test body. Keep inputs tiny (inline byte strings; no fixture files
   in version control). The point is to confirm the binary runs, not
   to verify scientifically correct output — a 30-bp FASTA, a single
   VCF record, or a one-hairpin RNA are all plenty.

4. Add the upstream tool to the conda install line in
   `.github/workflows/ci-nightly.yml`.

5. PR it.

## Why tests skip when the binary is missing

Each test starts with:

```rust
fn skip_if_missing(adapter: &dyn Adapter) -> bool {
    match adapter.probe() {
        Ok(report) if report.ok => false,
        _ => {
            eprintln!(
                "Skipping E2E test — `{}` upstream binary not installed on PATH.",
                adapter.info().id
            );
            true
        }
    }
}
```

Without that check, the crate would fail loudly on every developer
machine that doesn't have whichever niche bioinformatics tool the test
targets. With it, the crate's a no-op on machines that aren't set up
for the specific adapter — and a full integration test on machines
(and CI runners) that are.

The trade-off: a `0 failed; N ignored` result on a developer machine
**does not** prove the adapter works. It just proves nothing exploded.
The CI nightly workflow is the load-bearing run; this is the
"developer machines stay quiet" half of the deal.

## Adapters not currently E2E-tested

A handful of conda-installable adapters were deliberately left out of
this crate. Each entry below pairs the adapter with the reason it
isn't yet wired up:

- **`star`** (Phase 18.7) — STAR refuses to run on references smaller
  than ~1 MB and wants 1–30 GB of RAM for its smallest meaningful
  index. A genuine smoke test would need fixture files > 1 KB and
  several seconds of memory allocation; out of scope for an inline-
  byte-string-only crate.
- **`salmon`** (Phase 18.6) — Same issue as STAR: salmon's quant
  workflow requires a transcriptome large enough that the smallest
  realistic test exceeds the 1 KB fixture cap.
- **`kallisto`** (Phase 18.6) — Same issue as salmon — kallisto's
  pseudoalignment needs a transcriptome much larger than a single-
  line inline FASTA can supply.
- **`foldseek`** (Phase 18.5) — FoldSeek's `easy-search` needs either
  a pre-built structure database or a PDB query large enough to
  contain real backbone geometry. A minimum-viable PDB exceeds the
  1 KB fixture cap.
- **`beast2`** (Phase 30) — BEAST 2 takes a BEAUti-generated XML
  describing the phylogenetic model. Authoring such an XML inline
  reliably is non-trivial and the BEAUti-equivalent boilerplate would
  exceed the 1 KB cap.
- **`gatk`** / **`deepvariant`** (Phase 19) — GATK requires a 1+ GB
  jar download plus JVM tuning that's outside a clean conda install;
  DeepVariant is Docker-only.
- **`cromwell`** (Phase 22.5) — Cromwell is a Java JAR rather than a
  console-script. The adapter requires the user to supply an absolute
  jar path, which means there's no single conda binary the test can
  exercise without adding fixture files for the JAR.
- **`nupack`** (Phase 28) — Academic license restricts redistribution
  via conda or any other package channel.
- **`cpptraj`** (Phase 5.5 / 32.5) — cpptraj ships as part of
  AmberTools, which is a 1+ GB conda install. The conda CI runner
  budget doesn't justify pulling it in just for a single adapter
  smoke test.
- **`psi4`** / **`nwchem`** (Phase 25) — Both are heavyweight quantum-
  chemistry packages (Psi4 alone is ~500 MB after conda install).
  Skipped to keep the CI runner image lean; the xTB adapter covers
  the "small molecule QM" smoke-test slot.
- **`vina`** (Phase 34) — AutoDock Vina needs PDBQT files for both
  receptor and ligand. Even the smallest realistic receptor exceeds
  the 1 KB fixture cap.
- **`chopchop`** / **`crispor`** (Phase 35) — Both are pip-only Python
  packages. The CI workflow uses conda exclusively for tool install,
  so neither has a single-command bioconda path.
- **`planemo`** (Phase 22 / 22.5) — Same story as chopchop: pip-only,
  no bioconda package.
