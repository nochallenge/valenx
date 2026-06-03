# Valenx

> A native, open-source desktop simulation suite — written in Rust, no
> browser, no subscription, no vendor lock-in. One app spanning:
> **engineering** — aerospace (launch-vehicle ascent, orbital mechanics,
> re-entry), CFD, FEA, electromagnetics, multibody dynamics, thermal,
> battery, and a parametric CAD/CAM suite; **chemistry & materials** —
> molecular dynamics, quantum chemistry, reaction dynamics, and
> cheminformatics; and **computational biology** — genomics, sequence
> alignment, phylogenetics, population genetics, protein structure &
> design, RNA/mRNA design, CRISPR gene-editing, cryo-EM, and systems
> biology. Native Rust solvers plus 141 open-source tool integrations.

**Status:** `0.1.0-alpha.1` — pre-release. The workflow loop is
usable end-to-end (load project, click a case, **Prepare**, **Run**,
inspect results) but real-world testing is just beginning. Expect
sharp edges; please file issues.

## What it does

Valenx unifies a stack of [open-source simulation tools](#supported-solvers)
behind one Rust-native desktop shell so you can:

- Run **CFD / FEA / EM / battery / multibody** simulations from one UI
  (OpenFOAM, SU2, CalculiX, Code_Aster, Elmer, OpenRadioss, openEMS,
  Meep, PyBaMM, MuJoCo, preCICE, …)
- Drive **molecular dynamics** (GROMACS, LAMMPS, OpenMM, NAMD,
  AmberTools, HOOMD-blue)
- Predict **protein structure** (AlphaFold 2/3, ESMFold, OpenFold,
  RoseTTAFold, OmegaFold, ColabFold)
- Design **proteins** (RFdiffusion, ProteinMPNN, ESM-IF, Chroma,
  RFantibody, ESM3)
- Design **CRISPR guides** + analyze edit outcomes (CHOPCHOP, CRISPOR,
  Cas-OFFinder, inDelphi, FORECasT, BE-Designer, PrimeDesign)
- Fold **RNA** + design **mRNA vaccines** (ViennaRNA, NUPACK, mfold,
  LinearFold, DNA Chisel, LinearDesign)
- Reconstruct **cryo-EM** volumes (RELION, EMAN2, CTFFIND)
- Pipe everything through **reproducible workflows** (Snakemake,
  Nextflow, Cromwell, cwltool)

**The idea:** native Rust engines for the core science, one app, your laptop —
plus optional adapters to 141 external tools when you want them. No cloud, no
API keys, your data never leaves your machine.

## Native engines — included, nothing to install

A large part of Valenx runs **without any external tool** — native in-house
Rust solvers ship inside the app and work out of the box:

- **Computational biology** — a 14-panel **Genetics Workbench**: sequence
  analysis, pairwise + multiple alignment, phylogenetics, population genetics,
  RNA secondary structure, RNA/mRNA design, molecular dynamics,
  cheminformatics, macromolecular structure (PDB/mmCIF, DSSP, superposition),
  quantum chemistry (Hartree–Fock / MP2), genomics, systems biology, docking,
  and CRISPR / gene-edit design — all native (`valenx-bioseq`, `valenx-align`,
  `valenx-phylo`, `valenx-rnastruct`, `valenx-md`, `valenx-qchem`, …).
- **Engineering** — native 3D external-aerodynamics (a **wind-tunnel
  workbench** with k-ω SST turbulence), 2D laminar CFD, finite-element analysis
  (**8 native solvers**: static, modal, thermal, nonlinear, plasticity, beam),
  a CAD geometry kernel (primitives, booleans, fillets, NURBS surfaces), and a
  launch / orbital workbench.

The external tools below are **optional** — reach for them when you want a
reference implementation, a GPU/ML model, or a domain not yet native. A few
domains are still external-only and on the roadmap: a native
**electromagnetics** solver, **parametric CAD** history, native **unstructured
meshing**, **industrial/turbulent CFD**, and **DFT**. **Contributors welcome —
AI-assisted included** (see [CONTRIBUTING.md](./CONTRIBUTING.md) +
[AGENTS.md](./AGENTS.md)).

## Install

| Platform | How |
| --- | --- |
| Windows | Download `Valenx-<ver>-x86_64.msi` from [Releases][releases] and run it |
| macOS   | Download `Valenx-<ver>.dmg` from [Releases][releases] and drag to Applications |
| Linux (Debian/Ubuntu) | `sudo apt install ./valenx_<ver>_amd64.deb` |
| Linux (Fedora/RHEL/openSUSE) | `sudo dnf install ./valenx-<ver>.x86_64.rpm` |

Full installer guide with per-OS pinning recipes:
[docs/INSTALLER.md](./docs/INSTALLER.md).

[releases]: https://github.com/nochallenge/valenx/releases

## Build from source

```sh
# Linux / macOS
git clone https://github.com/nochallenge/valenx
cd valenx
cargo build --release -p valenx-app
./target/release/valenx
```

```powershell
# Windows (PowerShell)
git clone https://github.com/nochallenge/valenx
cd valenx
cargo build --release -p valenx-app
.\target\release\valenx.exe
```

Requires Rust 1.88+. See [CONTRIBUTING.md](./CONTRIBUTING.md) for
full dev setup.

## Documentation

- **[QUICKSTART.md](./QUICKSTART.md)** — five-minute walkthrough
- **[docs/INSTALLER.md](./docs/INSTALLER.md)** — install + pin to dock/taskbar
- **[STATUS.md](./STATUS.md)** — what's done, scaffolded, documentation-only
- **[ARCHITECTURE.md](./ARCHITECTURE.md)** — how the pieces fit together
- **[ROADMAP.md](./ROADMAP.md)** — 20-year plan
- **[CHANGELOG.md](./CHANGELOG.md)** — release history
- **[docs/CI.md](./docs/CI.md)** — CI policy (manual-trigger only; here's why)
- **[CONTRIBUTING.md](./CONTRIBUTING.md)** — how to contribute + dev setup
- **[SECURITY.md](./SECURITY.md)** — vulnerability disclosure
- **[CODE_OF_CONDUCT.md](./CODE_OF_CONDUCT.md)** — Contributor Covenant v2.1
- **[POLICIES.md](./POLICIES.md)** — SemVer, deprecation, MSRV, dependency rules
- **[rfcs/](./rfcs/)** — RFC process + accepted proposals

## What's actually working

| | |
|---|---|
| Live adapters | **141 fully live** with real prepare/run/collect across CFD / FEA / EM / chemistry / MD / battery / multibody / coupling + a 123-adapter biology stack (structure prediction, alignment, variant calling, single-cell, workflow managers, viewers, cheminformatics, quantum chemistry, protein design, RNA structure, phylogenetics, CRISPR, cryo-EM, microscopy, mRNA design, and more) |
| Tests | 10,000+ passing, 0 failed, 0 clippy warnings, 0 rustdoc warnings (workspace-wide) |
| Workflow loop | load project → click case → Run / Prepare / Open workdir / Run-from-prepared, all live |
| Sweeps + export | grid / Latin-hypercube / gradient-descent optimisers, materialise → run (sync or threaded) → assemble dataset (csv / npy / npz / manifest.json) |
| Mesh quality | Equiangle skewness · cell-face orthogonality · aspect-ratio + skewness histograms · cell-face + edge adjacency. Surfaced in the GUI Quality panel and via the headless `valenx-mesh-info` CLI (text / JSON / `--check METRIC=THRESHOLD` for CI gates) |
| HPC | `valenx-executor-slurm` does end-to-end remote-cluster submission: rsync staging + ssh-dispatched `sbatch` / `squeue` / `sacct` / `scancel`, sacct fallback for terminal state, GPU (`--gres`) + multi-rank `srun`, post-completion `fetch_results` |
| Visual results | OpenFOAM + Elmer + CalculiX + SU2 populate `Results.fields` with parsed `.vtu` / `.vtk` legacy / `.frd` / appended-binary data; auto-routed via vtk_dispatch. Vector + tensor fields render as magnitude / Frobenius norm. Viewport paints filled triangles + wireframe by field value with colour-bar legend + time-scrubbing |
| Compliance | RBAC (3 roles + per-project override) + append-only audit log with SHA-256 chain. `valenx-audit verify` and `valenx-audit tail` CLIs for offline integrity checks + headless inspection |
| CLI tooling | `valenx-init` scaffolds a project. `valenx-validate` runs a structural pre-flight on a project bundle. `valenx-results` inspects the `results.json` sidecar a finished run leaves on disk. `valenx-report` writes a self-contained HTML report and/or a flat scalar history CSV from a finished run. All four are exit-code-driven; the inspectors offer JSON output for downstream CI tooling |

## Supported solvers

Beyond its [native engines](#native-engines--included-nothing-to-install),
Valenx integrates **141 external open-source tools** as optional adapters — for
reference implementations, GPU/ML models, and domains not yet native. Each tool
below is wrapped by a Valenx adapter crate
(`crates/valenx-adapters/<domain>/valenx-adapter-<tool>`) that handles probing,
case preparation, subprocess execution, and result collection. You install only
the ones you actually use.

**CFD / FEA / EM / multiphysics**
OpenFOAM, SU2, gmsh, Netgen, FreeCAD, CalculiX, Elmer, Code_Aster,
OpenRadioss, Cantera, openEMS, Meep, PyBaMM, MuJoCo, preCICE.

**Molecular dynamics**
GROMACS, LAMMPS, OpenMM, NAMD, AmberTools sander, HOOMD-blue.
Trajectory analysis: MDAnalysis, MDTraj, PLUMED, ProDy, cpptraj.

**Protein structure prediction + design**
AlphaFold 2, AlphaFold 3, ESMFold, OpenFold, RoseTTAFold, OmegaFold,
ColabFold, FoldSeek. Design: RFdiffusion, ProteinMPNN, ESM-IF,
Chroma, RFantibody, ESM3, ESM Cambrian.

**Sequence alignment + search + variant calling**
BWA, minimap2, Bowtie2, HISAT2, STAR, MMseqs2, DIAMOND, BLAST+,
MAFFT, MUSCLE, Clustal Omega, T-Coffee, HMMER, samtools, bcftools,
GATK HaplotypeCaller, DeepVariant.

**Transcript quantification + single-cell**
Salmon, Kallisto, Scanpy, scVI, Seurat, AnnData.

**Cheminformatics + quantum chemistry**
RDKit, DeepChem, Open Babel, Avogadro 2, Psi4, NWChem, xTB.

**Molecular docking + viewers**
AutoDock Vina, AutoDock 4. PyMOL, VMD, ChimeraX, IGV, Mol\*, NGL.

**CRISPR design + edit-outcome prediction**
CHOPCHOP, CRISPOR, Cas-OFFinder. Base/prime: BE-Designer, BE-Hive,
PrimeDesign, pegFinder. Outcomes: inDelphi, FORECasT, AlphaMissense,
CRISPRitz.

**RNA structure + mRNA design**
ViennaRNA, RNAstructure, NUPACK, mfold, EternaFold, LinearFold.
mRNA: DNA Chisel, LinearDesign, iCodon. Tertiary: SimRNA.

**Cryo-EM + microscopy**
RELION, EMAN2, CTFFIND. Bioimage: Fiji, CellProfiler, Ilastik.

**Phylogenetics + population genetics**
IQ-TREE, RAxML-NG, FastTree, BEAST 2, MrBayes. PopGen: SLiM,
msprime, tskit.

**Systems / synthetic biology**
COPASI, BioNetGen, PhysiCell, Smoldyn, MCell. SynBio: pySBOL, j5,
Cello, pydna, Jalview.

**Sequencing read simulators + Rosetta family**
ART, wgsim, Badread. Rosetta, PyRosetta.

**DNA structural geometry + pharmacokinetics**
X3DNA, Curves+, DSSR. PK-Sim.

**Workflow managers**
Snakemake, Nextflow, planemo, Cromwell, cwltool.

Several upstream tools ship under academic / non-commercial terms
(AlphaFold 3, NAMD, VMD, CTFFIND, Rosetta, ViennaRNA, NUPACK, mfold,
AlphaMissense, X3DNA family, Curves+). Valenx surfaces these via
mandatory `"academic"`-keyworded probe warnings — you choose whether
to install the binary based on your use case.

Full per-tool status and capability matrix: [STATUS.md](./STATUS.md).

## Design principles

1. **Native desktop app.** No browser. No localhost. Downloads like
   FreeCAD, runs like FreeCAD.
2. **Open source, dual MIT / Apache-2.0.** Free forever, commercial
   friendly, ecosystem-standard.
3. **Rust first.** Safe, fast, modern tooling.
4. **Bundle existing solvers.** Reuse decades of validated physics.
   Don't rewrite what already works.
5. **Replace where it matters.** The UX, the integration, the
   workflow — that's what we build new.

## License

Dual-licensed under either of:

- **MIT License** ([LICENSE-MIT](./LICENSE-MIT))
- **Apache License, Version 2.0** ([LICENSE-APACHE](./LICENSE-APACHE))

at your option.

This is the [Rust ecosystem
standard](https://rust-lang.github.io/api-guidelines/necessities.html#crate-and-its-dependencies-have-a-permissive-license-c-permissive)
— the same dual layout used by rustc, Tokio, Serde, Bevy, and most
of the Rust crate ecosystem.

### Contribution

Unless you explicitly state otherwise, any contribution intentionally
submitted for inclusion in Valenx by you, as defined in the
Apache-2.0 license, shall be dual-licensed as above, without any
additional terms or conditions.

See [CONTRIBUTING.md](./CONTRIBUTING.md) for guidelines.

## Contribute

We're pre-alpha. If you want to help — shell scaffolding, adapters,
physics solvers, documentation — start here:

- Read [CONTRIBUTING.md](./CONTRIBUTING.md) for dev setup and workflow
- Big ideas → file an RFC (see [rfcs/README.md](./rfcs/README.md))
- Small fixes → open a PR directly
- Security issues → [SECURITY.md](./SECURITY.md) (please don't file
  public issues)
- Be respectful → [CODE_OF_CONDUCT.md](./CODE_OF_CONDUCT.md)
