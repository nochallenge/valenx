# Valenx QA runbook

This is the operating manual for validating the Valenx workspace locally. It
covers **what the QA process is**, the **scoped-safe test rules and why they
exist**, **what each crate's tests cover**, and **how to run the one-command
harness**.

> TL;DR — run `./scripts/qa.sh` (or `scripts\qa.ps1` on Windows). It runs the
> entire safe validation suite and exits non-zero on any failure. It never runs
> `cargo test --workspace`, never launches the app, never opens a file dialog.

---

## 1. Why QA is scoped — the file-dialog crash history

**A blanket `cargo test --workspace` is forbidden in this repo.**

The `valenx-app` crate is the desktop application. Its library has UI-coupled
unit tests, and some of them exercise the open/save code paths, which call
[`rfd::FileDialog`](https://docs.rs/rfd). `rfd::FileDialog::save_file` /
`pick_file` pop a **native OS file-picker window** and block the calling thread
until the user clicks something. In a headless / unattended run there is no
user — the test hangs forever. Earlier in the project this wedged a developer
machine hard enough to need a reboot.

So the iron rule, also recorded in `docs/INCREMENTAL_GOAL.md`:

- **NEVER** `cargo test --workspace`.
- **NEVER** `cargo test` with no `-p` filter.
- **NEVER** `cargo test -p valenx-app` *unfiltered* — its lib unit tests spawn
  the dialog.
- **NEVER** `cargo run`, `cargo bench`, or launch any `target/**/valenx*`
  binary.
- **NEVER** run anything that calls `rfd::FileDialog`.

What **is** safe, and what the harness does:

| Allowed command | Why it is safe |
|---|---|
| `cargo test -p <pure-crate>` (one at a time) | The 20 pure computational crates are pure-Rust algorithm libraries — they do not depend on `rfd`, open no window, spawn no subprocess. |
| `cargo test -p valenx-app headless_ui_tests` | A **name filter**. Every windowless egui-logic test lives in a module named `headless_ui_tests`; the filter selects only those and excludes every file-dialog test. The selected tests drive panel logic in an off-screen `egui::Context` and never call `rfd`. |
| `cargo test -p valenx-app --test pipeline_e2e` | Compiles and runs **only** the one integration-test file `crates/valenx-app/tests/pipeline_e2e.rs`. That file contains the cross-crate end-to-end suite and the adapter smoke test — none of it touches `rfd`. |
| `cargo check --workspace` | Build-only. Never executes a test or `main()`. |
| `cargo clippy --workspace --all-targets -- -D warnings` | Lint-only. Never executes anything. |
| `cargo doc --workspace --no-deps` | Doc-build only. |

`cargo check` / `clippy` / `doc` are safe workspace-wide because they compile
but never *run* code. Only `cargo test` / `run` / `bench` execute, and only
those are scoped.

---

## 2. The one-command harness

Two equivalent entry points live in `scripts/`:

- **`scripts/qa.sh`** — Bash (Linux, macOS, Git-Bash / WSL on Windows).
- **`scripts/qa.ps1`** — PowerShell (Windows-native).

### Run everything

```bash
./scripts/qa.sh
```
```powershell
.\scripts\qa.ps1
```

It runs, in order:

1. `cargo test -p <crate>` for each of the 20 pure computational crates;
2. `cargo test -p valenx-app headless_ui_tests` (the workbench UI-logic tests);
3. `cargo test -p valenx-app --test pipeline_e2e` (the cross-crate e2e suite);
4. `cargo check --workspace`;
5. `cargo clippy --workspace --all-targets -- -D warnings`;
6. `cargo doc --workspace --no-deps`;
7. `cargo deny check` (license + vulnerability + advisory gate — matches the
   round-8 SECURITY.md claim that the local QA script enforces what CI no
   longer auto-runs).

Each step prints `PASS` / `FAIL`; the script exits `0` only if every step
passed, and prints the list of failures otherwise.

### Partial runs

```bash
./scripts/qa.sh --tests    # only steps 1-3 (the scoped test runs)
./scripts/qa.sh --gates    # only steps 4-7 (the workspace gates incl. cargo deny check)
./scripts/qa.sh --help
```
```powershell
.\scripts\qa.ps1 -Tests
.\scripts\qa.ps1 -Gates
```

### Expected baseline

A healthy `master` exits `ALL QA STEPS PASSED`. The **only** accepted noise is
~5 pre-existing rustdoc warnings in the untouched `valenx-solvespace-3d` crate
— a known, documented baseline that predates this harness. Everything else
green means green.

### Runtime

A full run is dominated by compilation on a cold `target/`. With a warm cache
the test phase is a few minutes; the `valenx-aero` crate is the slowest
(a real 3-D Navier-Stokes solver — its test suite runs on deliberately coarse
grids to keep the wall-clock feasible).

---

## 3. What each crate's tests cover

### The 20 pure computational crates

These are pure-Rust algorithm libraries (Round 6 computational-biology blocks
plus the native physics solvers). Each has a `cargo test -p <crate>` suite of
unit tests plus, where a published reference exists, reference-value tests.

| Crate | Domain | Test coverage focus |
|---|---|---|
| `valenx-bioseq` | Biological sequences | FASTA/FASTQ/GenBank/EMBL I/O, reverse-complement, transcription, translation across NCBI code tables, ORF finding, restriction digest, primer design, ProtParam protein properties. |
| `valenx-align` | Sequence alignment | Needleman-Wunsch / Gotoh affine / Smith-Waterman, banded + Hirschberg, BLOSUM/PAM matrices, k-mer + FM-index search, **progressive MSA (affine-gap profile alignment)**, profile + pair HMMs. |
| `valenx-phylo` | Phylogenetics | Distance methods (p / JC69 / K80 / TN93) with exact closed-form reference tests, UPGMA / NJ / BIONJ, parsimony, Felsenstein likelihood, Robinson-Foulds / quartet distances, consensus. |
| `valenx-popgen` | Population genetics | Allele/genotype frequencies, Hardy-Weinberg, F-statistics, LD, the site-frequency spectrum, neutrality + selection scans, coalescent simulation. |
| `valenx-rnastruct` | RNA secondary structure | Nussinov + Zuker MFE folding against the **full Turner-2004 parameter set** (reference-value tests vs the analytic Turner sum and ViennaRNA), McCaskill partition function, suboptimal + multiloop folding, ct/bpseq I/O. |
| `valenx-md` | Molecular dynamics | Bonded forces (bond / angle / dihedral / improper) verified against finite differences, Lennard-Jones + PME electrostatics, thermostats / barostats, PBC, integrators, RDF. |
| `valenx-cheminf` | Cheminformatics | SMILES / SMARTS / MOL-SDF / InChI parsing, ring perception + aromaticity, descriptors (logP, TPSA, Lipinski, Veber), Morgan / MACCS fingerprints, Gasteiger PEOE charges, periodic-table reference data. |
| `valenx-biostruct` | Macromolecular structure | PDB / mmCIF I/O, geometry (distances, torsions, Ramachandran, radius of gyration, SASA), DSSP secondary structure, Kabsch / quaternion superposition, TM-score, base-pair parameters. |
| `valenx-qchem` | Quantum chemistry | One- and two-electron integrals (McMurchie-Davidson), RHF / UHF SCF with DIIS, MP2, extended Hückel, Mulliken / Löwdin populations — reference-value energy tests on small molecules. |
| `valenx-genomics` | Genomics | Read QC, assembly, variant normalisation / calling / stats, CRISPR guide design + PAM scanning + off-target scoring, amplicon analysis. |
| `valenx-sysbio` | Systems biology | Reaction-network model, mass-action / Michaelis-Menten / Hill kinetics, RK4 / RK45 / BDF ODE integration, Gillespie SSA + tau-leaping, the from-scratch simplex / FBA layer, conserved-moiety analysis. |
| `valenx-dock-screen` | Docking & virtual screening | Vina + AutoDock4 scoring functions, affinity-grid precompute, Lamarckian GA / Monte-Carlo search, pose clustering, ensemble + consensus, MM-GBSA rescoring, ligand prep (protonation + torsion tree). |
| `valenx-genediting` | Gene editing | CRISPR nuclease database + guide design, NHEJ / HDR / base-editing / prime-editing design, mRNA construct design, codon optimisation, off-target safety aggregation. |
| `valenx-structpredict` | Structure prediction | Classical homology + ab-initio modelling, CCD loop closure, fragment assembly, statistical potentials, rotamer repacking, Kabsch RMSD + GDT, classical cryo-EM reconstruction (CTF, back-projection, FSC). |
| `valenx-rnadesign` | Synthetic RNA design | Inverse folding, riboswitch design, the multi-objective simulated-annealing optimiser, ensemble-defect metric, in-silico validation, DNA-template + IVT-plan generation. |
| `valenx-aero` | 3-D external aerodynamics | The immersed-boundary RANS solver, cut-cell vs staircased walls, k-ε / k-ω SST turbulence, drag / lift / moment coefficients, wake survey, AoA sweep — tested on coarse grids for feasible runtime. |
| `valenx-cfd-native` | 2-D CFD | The SIMPLE staggered-grid Navier-Stokes solver, lid-driven-cavity recirculation, channel-flow parabolic profile, k-ε turbulence, transient time-marching. |
| `valenx-fem` | Finite-element analysis | Linear-static / modal / thermal / nonlinear Tet4 solvers, von Mises plasticity, penalty contact, Newmark dynamics, eigenvalue buckling — verified against analytic beam / conduction references. |
| `valenx-pathtrace` | Path tracing | The Monte-Carlo path tracer, BVH traversal, BSDF sampling, MIS, dielectric refraction, the a-trous denoiser, volumetric rendering — white-furnace + analytic-lighting tests. |
| `valenx-render-bridge` | Render bridge / PBR | The Cook-Torrance PBR library, the split-sum IBL, irradiance-volume GI, the WGSL shader (naga-validated), HDR environment decoding, render-job persistence. |

### Workbench UI-logic tests — `cargo test -p valenx-app headless_ui_tests`

151 headless tests covering the desktop panels: the 14 Genetics-workbench
panels, the 8-section Wind Tunnel workbench, the CAD workbenches, and the GPU
render path. Each panel is drawn in a windowless `egui::Context` across
representative states (no panic), has its Run action driven against the real
backend crate (sane result), and is fed bad input (graceful error). No window
is opened, no file dialog is shown.

### Cross-crate e2e tests — `cargo test -p valenx-app --test pipeline_e2e`

The `crates/valenx-app/tests/pipeline_e2e.rs` integration suite. It runs
realistic multi-crate workflows start-to-finish and asserts the final result is
physically / biologically sane — catching integration bugs (a type-shape
mismatch at a crate seam, a unit error, a convention clash) that a single
crate's unit tests cannot see:

1. **Comparative genomics** — FASTA parse → pairwise + MSA alignment → distance
   matrix → neighbor-joining tree; asserts the tree recovers the known clades.
2. **Drug-design prep** — SMILES → cheminf descriptors → dock-screen ligand
   preparation (protonation + torsion tree).
3. **Protein structure** — PDB parse → geometry (Rg) + DSSP secondary
   structure → Kabsch superposition.
4. **Gene expression** — DNA → ORF finding → translation → ProtParam protein
   properties.
5. **Systems biology** — a hand-built reaction network → ODE time-course AND
   Gillespie SSA; cross-checks the two engines agree in the large-N limit.
6. **Virtual wind tunnel** — a cube mesh → the aero immersed-boundary RANS
   solver → a drag coefficient in the bluff-body band.
7. **Quantum chemistry** — an H2 geometry → RHF/STO-3G → total energy +
   properties against the textbook reference.

Plus the pre-existing `gmsh → OpenFOAM` adapter smoke test, which tolerates the
tools being absent.

---

## 4. Continuous integration

**The local `scripts/qa.sh` harness is the source of truth.** Run it before
pushing; that is the gate.

The repo *does* ship `.github/workflows/ci.yml`, `ci-nightly.yml`, and
`release.yml`, but they are gated to **`workflow_dispatch`** (manual runs)
plus the explicit `release` tag — the maintainer manages Actions minutes
tightly and the workflows do not fire on every push. The active CI also
invokes `scripts/qa.sh --tests` rather than a blanket `cargo test
--workspace`, mirroring the local rules.

The template below documents the scoped-rule recipe in stand-alone form
for forks / mirrors that want a push-triggered version.

```yaml
# .github/workflows/qa.yml  —  TEMPLATE ONLY, not active in this repo.
# Mirrors scripts/qa.sh. Note the scoped per-crate test matrix: there is
# deliberately no `cargo test --workspace` step — see docs/QA.md section 1.
name: qa
on:
  workflow_dispatch:        # manual only — never on every push
jobs:
  scoped-tests:
    runs-on: ubuntu-latest
    strategy:
      fail-fast: false
      matrix:
        crate:
          - valenx-bioseq
          - valenx-align
          - valenx-phylo
          - valenx-popgen
          - valenx-rnastruct
          - valenx-md
          - valenx-cheminf
          - valenx-biostruct
          - valenx-qchem
          - valenx-genomics
          - valenx-sysbio
          - valenx-dock-screen
          - valenx-genediting
          - valenx-structpredict
          - valenx-rnadesign
          - valenx-aero
          - valenx-cfd-native
          - valenx-fem
          - valenx-pathtrace
          - valenx-render-bridge
    steps:
      - uses: actions/checkout@v4
      - uses: dtolnay/rust-toolchain@stable
      - run: cargo test -p ${{ matrix.crate }}
  app-scoped-tests:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: dtolnay/rust-toolchain@stable
      # NAME-FILTERED + single-file only — never unfiltered valenx-app tests.
      - run: cargo test -p valenx-app headless_ui_tests
      - run: cargo test -p valenx-app --test pipeline_e2e
  gates:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: dtolnay/rust-toolchain@stable
        with: { components: clippy }
      - run: cargo check --workspace
      - run: cargo clippy --workspace --all-targets -- -D warnings
      - run: cargo doc --workspace --no-deps
```

---

## 5. Test-coverage measurement

Line / region coverage of the pure crates is measured with
[`cargo-llvm-cov`](https://github.com/taisukef/cargo-llvm-cov) (a dev tool):

```bash
cargo install cargo-llvm-cov            # one-time
cargo llvm-cov --no-cfg-coverage -p valenx-rnastruct --summary-only
```

Run it **scoped per crate**, exactly like `cargo test` — the same lockdown
applies. As of the coverage pass that introduced this runbook the 19+ pure
crates sit at roughly **88-96 % line coverage**; see `docs/VALIDATION.md` for
the per-crate numbers and the gap-filling history.
