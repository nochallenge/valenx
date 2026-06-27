# Valenx — AI-controlled design and simulation

> A native, open-source desktop simulation suite — written in Rust, no
> browser, no subscription, no vendor lock-in — and **AI-drivable end to end**
> (every workbench exposes named accessible widgets + an agent-command bridge,
> so an AI agent can operate the whole app, not just a human). One app spanning:
> **engineering** — aerospace (launch-vehicle ascent, orbital mechanics,
> re-entry), CFD (plus a GPU-compute path), FEA, electromagnetics, multibody
> dynamics, thermal, **thermodynamics** (equations of state), **acoustics**,
> **optics**, **topology optimization**, **large-deformation (MPM)**, battery,
> and a parametric CAD/CAM suite (with a B-Rep solid kernel); **chemistry &
> materials** —
> molecular dynamics, quantum chemistry, reaction dynamics, and
> cheminformatics; and **computational biology** — genomics, sequence
> alignment, phylogenetics, population genetics, protein structure &
> design, RNA/mRNA design, CRISPR gene-editing, cryo-EM, and systems
> biology; and **neuroengineering** — neural-interface / BCI stimulation
> (extracellular fields, Hodgkin–Huxley cable models, bioheat); and
> **gravitational physics** — an in-house general-relativity / black-hole
> engine (Kerr–Newman spacetimes, geodesics, Hawking thermodynamics, and a
> shadow ray-tracer); **quantum computing** — a state-vector circuit
> simulator; and **defense / mission planning** — a map-based tactical planner
> (OSM basemap, MIL-STD-2525 symbology, A\* routing, line-of-sight), analysis /
> planning only. **A complete in-house suite** — its own native Rust solvers
> across every domain above, validated against published ground truth — with
> *optional* adapters to external tools when you want them.

**Status:** `0.1.0-alpha.1` — pre-release. The workflow loop is
usable end-to-end (load project, click a case, **Prepare**, **Run**,
inspect results) but real-world testing is just beginning. Expect
sharp edges; please file issues.

## 📺 Demo

[![Watch the Valenx demo](https://img.youtube.com/vi/n1Je-qk2es0/maxresdefault.jpg)](https://www.youtube.com/watch?v=n1Je-qk2es0)

A walkthrough of Valenx — generating and simulating ~190 products from one app,
the 3-D editor, and the in-app AI agent. It's an **early developer build, not a
commercial release** — **contributors wanted.**

## What it does

Valenx is a **complete, native simulation suite** — its own from-scratch Rust
solvers across engineering, physics, chemistry, and biology, in one desktop app,
on your laptop. The native engines are the product, and they're
[validated against published ground truth](#validation):

- **Aerospace** — launch-vehicle ascent to orbit, 3-D orbital mechanics (J2),
  a 6-DOF attitude core, impulsive + Lambert-rendezvous maneuvers
- **CFD / FEA / EM / multibody / thermal** — a 3-D wind-tunnel aero workbench,
  k-ω SST CFD, 8 FEA solvers, topology optimization, MPM large-deformation
- **Chemistry & materials** — molecular dynamics, quantum chemistry
  (Hartree–Fock / MP2 / DFT), reaction dynamics (AIMD), cheminformatics
- **Computational biology** — a 14-panel Genetics Workbench: alignment,
  phylogenetics, population genetics, RNA/mRNA design, protein structure,
  docking, CRISPR gene-editing
- **More native engines** — general relativity / black-holes, neural-interface
  / BCI, geometric optics, a quantum-circuit simulator, and a parametric
  CAD / CAM / CAE suite with a B-Rep solid kernel

Everything is **AI-drivable end to end** — every workbench exposes named
accessible widgets + an agent-command bridge, so an AI agent can set parameters,
run, and read results across the whole app, not just a human.

**Optional interop:** when you specifically want a reference implementation, a
GPU/ML model (e.g. AlphaFold), or a domain not yet native, Valenx can *also*
drive external tools through optional adapters — a convenience layer, never a
requirement. No cloud, no API keys, your data never leaves your machine.

## Native engines — included, nothing to install

A large part of Valenx runs **without any external tool** — native in-house
Rust solvers ship inside the app and work out of the box:

- **Computational biology** — a 14-panel **Genetics Workbench**: sequence
  analysis, pairwise + multiple alignment, phylogenetics, population genetics,
  RNA secondary structure, RNA/mRNA design, molecular dynamics,
  cheminformatics, macromolecular structure (PDB/mmCIF, DSSP, superposition),
  quantum chemistry (Hartree–Fock / MP2 / Kohn–Sham DFT), genomics, systems biology, docking,
  and CRISPR / gene-edit design — all native (`valenx-bioseq`, `valenx-align`,
  `valenx-phylo`, `valenx-rnastruct`, `valenx-md`, `valenx-qchem`, …).
- **Aerospace / GN&C** — a launch-vehicle + orbital engine: multi-stage
  gravity-turn ascent to orbit, 3-D orbital mechanics with **J2** oblateness, a
  **6-DOF rigid-body** attitude core, and impulsive maneuvers (Hohmann,
  bi-elliptic, plane-change, **Lambert rendezvous**) — validated against
  textbook/analytic results (`valenx-astro`; see [Validation](#validation)).
- **Fluids & structures** — a 3-D external-aerodynamics **wind-tunnel
  workbench** and 2-D CFD with Menter **k-ω SST** turbulence, plus
  finite-element analysis (**8 native solvers**: static, modal, thermal,
  nonlinear, plasticity, beam, buckling, transient dynamics).
- **CAD / CAM / CAE** — a geometry kernel (primitives, booleans, fillets, NURBS
  surfaces), a **parametric feature-tree history**, CAM toolpaths + G-code
  post-processing, technical drafting (DXF / SVG / PDF), assemblies with
  interference detection, architectural / structural modelling (IFC4, Eurocode),
  surface modelling (NURBS fitting + blends), and JT / STEP-AP242 / IGES interop.
- **Engineering design workbenches** — a growing set of native parametric design
  calculators, each a right-side panel (View menu) with a form → compute → readout
  and, where the geometry suits it, a live wireframe preview: helical **springs**
  (spring index, axial stiffness, Wahl factor, developed wire length + a 3-D helix
  preview), involute **gears** (pitch / base / addendum / dedendum diameters, mesh
  ratio + a tooth-profile preview), **piping** (NPS bore, flow + metal area, wetted +
  external surface), structural **frames** (section area + perimeter), **sheet-metal**
  bends (bend allowance / deduction), **fasteners** (ISO 4017 bolt dimensions),
  **collision** (AABB volume / surface / overlap), **geodesics** (haversine + rhumb
  distance, bearings, cross-track), and **field statistics** — all closed-form and
  computed in-process (`valenx-springs`, `valenx-gears`, `valenx-frames`, …). These
  are research-grade design readouts, not a commercial parametric-modelling
  replacement (the roadmap targets FreeCAD-class parity, not Fusion 360 / SolidWorks).
- **Reaction dynamics & graphics** — a **reaction-dynamics / AIMD** simulator
  (velocity-Verlet on quantum-chemistry forces) and a physically-based **path
  tracer** (light-tree importance sampling, bidirectional path tracing,
  subsurface scattering).
- **Gravitational physics / general relativity** — an in-house GR & black-hole
  engine (`valenx-relativity`): Schwarzschild / Kerr / Reissner–Nordström /
  **Kerr–Newman** spacetimes, curvature tensors via automatic differentiation,
  null + timelike **geodesics**, observables (photon sphere, ISCO, shadow
  radius), **Hawking thermodynamics**, and a **black-hole shadow ray-tracer** —
  checked against closed-form results (Schwarzschild photon sphere = 3M, shadow
  radius = √27·M).
- **Neuroengineering / BCI** — a neural-interface stimulation **and recording**
  suite (`valenx-neuro`): an implanted electrode's **extracellular FEM field**
  (reusing the FEA solver), **Hodgkin–Huxley** cable axons, the **Rattay
  activating function** coupling the two (an electrode recruits nearby
  neurons), **bioheat** tissue heating, and **electrode–tissue impedance**,
  plus an **unconditionally-stable implicit cable solver**, **myelinated
  saltatory fibers** (conduction velocity matched to the empirical 6·D rule),
  **strength–duration** curves, an **anisotropic-tensor FEM field**, and
  **multi-contact current steering**, and **extracellular recording** (the
  forward-EAP read-out model) — each validated against a textbook or
  closed-form result (see [Validation](#validation)).
- **Thermo-fluids & acoustics** — fluid **thermodynamics** (Peng–Robinson / SRK
  cubic equations of state, compressibility, saturation pressure — validated vs
  NIST), free-field **acoustic radiation** (monopole / dipole directivity), and a
  GPU-compute pressure-Poisson sweep for the CFD core (`valenx-thermo`,
  `valenx-acoustics`, `valenx-cfd-native`).
- **Mechanics depth** — **topology optimization** (SIMP, generative
  minimum-weight structures), the **Material Point Method** (2-D MLS-MPM for
  large-deformation / failure where mesh FEM inverts), and **porous-media /
  groundwater** flow (Darcy + unsaturated Richards) on the FEM core
  (`valenx-topopt`, `valenx-mpm`, `valenx-fem`).
- **Optics, quantum & signals** — geometric **optics** (thin-lens + paraxial
  ABCD ray-transfer), a **quantum-circuit** state-vector simulator
  (H/X/Y/Z/S/T/CNOT/CZ + Born-rule measurement, Bell + GHZ validated), and a
  **digital-waveform** (VCD) viewer (`valenx-optics`, `valenx-quantum`,
  `valenx-waveform`).
- **Defense / mission planning** — a map-based tactical planner: real OSM basemap
  tiles, **MIL-STD-2525 / APP-6** affiliation symbology, **A\*** routing over a
  terrain cost grid, and terrain-masked **line-of-sight**. Analysis / planning
  posture only — no weapons-lethality / targeting.
- **AI-drivable everywhere** — every workbench is a reactive panel with **named,
  accessible widgets** (an AccessKit tree) plus an **agent-command bridge**, so an
  AI agent — or any MCP client, via the in-tree `valenx-mcp` server — can set
  parameters, run, and read results across the whole app. This is the
  "AI-controlled" half of the name.

The external tools below are **optional** — reach for them when you want a
reference implementation, a GPU/ML model, or a domain not yet native. A few
domains are still external-only or only partly native and on the roadmap: a
native **electromagnetics** solver, fuller **unstructured meshing** (a native
voxel→hex mesher ships now; general tet meshing is next), and **large-scale
3-D / industrial CFD** (a GPU-compute pressure-Poisson path ships now; full
industrial scale is next). **Contributors welcome —
AI-assisted included** (see [CONTRIBUTING.md](./CONTRIBUTING.md) +
[AGENTS.md](./AGENTS.md)).

## Validation

Native solvers are checked against published references or analytic results —
the figures below are reproduced by the test suite, not asserted. Full detail —
every per-crate validation suite and the ~200 bugs surfaced and fixed by
running them — lives in [docs/VALIDATION.md](./docs/VALIDATION.md).

Or verify it yourself in seconds — **`valenx --self-test`** drives all **56
products'** real compute and prints a machine-parseable `PASS / FAIL / value`
report (**53 PASS · 0 FAIL · 3 SKIP**; 50 of them deep known-value checks). It
needs no GUI and exits non-zero on any failure, so it doubles as a CI gate — and
its line-oriented output is consumable by any agent or script, not just a human
(see [Self-validation](./docs/VALIDATION.md#self-validation--valenx---self-test)).

| Solver | Benchmark / reference | Result |
| --- | --- | --- |
| Orbital — `valenx-astro` | LEO→GEO Hohmann Δv vs textbook ~3.9 km/s | **3,892 m/s** total (2,425 + 1,467), ToF 5.27 h — [test](./crates/valenx-astro/src/maneuver.rs) |
| Orbital — `valenx-astro` | J2 secular nodal regression vs closed-form rate | within **5%** of the analytic dΩ/dt; `a`, `i` secularly unchanged — [test](./crates/valenx-astro/src/orbit3d.rs) |
| CFD — `valenx-cfd-native` | Lid-driven cavity vs Ghia, Ghia & Shin 1982 | centerline MAE **0.035 / 0.016 / 0.024** at Re 100 / 400 / 1000 |
| CFD — `valenx-cfd-native` | Poiseuille channel vs analytic | u_max **1.4949 vs 1.5000** (0.34% error) |
| CFD — `valenx-cfd-native` | Backward-facing step vs Armaly / Gartling | reattachment x_r/h ≈ **4.5**, inside the published envelope |
| FEA — `valenx-fem` | Constant-strain patch test | satisfied to **~1e-9**; Hex8 90% / Tet10 112% of the Euler–Bernoulli tip |
| MD — `valenx-md` | Argon NVE energy conservation | std-to-mean **< 1%**, no secular drift; lattice sum < 0.5% vs analytic |
| Quantum chem — `valenx-qchem` | Hydrogen atom, exact Kohn–Sham | reproduces **−0.212742 Ha**; LDA → uniform-gas and PBE → LDA limits hold |
| Docking — `valenx-dock-screen` | 1HVR / 3PTB / 1STP redocking | RMSD **0.305 / 0.263 / 0.139 Å** (mean 0.236), 100% < 2 Å |
| Neuro — `valenx-neuro` | Hodgkin–Huxley axon | textbook action potential (~100 mV overshoot, threshold, refractory) + propagating cable |
| Neuro — `valenx-neuro` | field + activating function | point source φ = I/(4πσr); cathodic depolarizes, anodic hyperpolarizes (Rattay) |
| Neuro — `valenx-neuro` | electrode recruitment | threshold **rises with electrode–fiber distance**; recruitment monotone in current |
| Neuro — `valenx-neuro` | bioheat + impedance | ΔT = Q/(4πk r); electrode access resistance R = 1/(4σa) |
| Neuro — `valenx-neuro` | implicit cable + myelinated fiber | unconditionally-stable implicit solver (stable at 100 µm, fixes the v1 RK4 limit); myelinated **CV ≈ 6·D** within ~6% (57 / 113 m/s at 10 / 20 µm), ∝ D not √D |
| Neuro — `valenx-neuro` | strength–duration | rheobase + **chronaxie 1.65 ms** (≈ ½ membrane τ); Lapicque constant-charge law holds to **< 1%** at short widths |
| Neuro — `valenx-neuro` | anisotropic FEM + steering | tensor-σ point source matches the **closed-form** anisotropic Green's function within ~10%; multi-contact **current steering** shifts the focus |
| Neuro — `valenx-neuro` | extracellular recording (forward-EAP) | **biphasic** spike with a dominant-negative (sink) phase; amplitude falls off with electrode distance; membrane currents conserve charge (Σ ≈ 0) |

## Download & Install

> **Pre-alpha (`0.1.0-alpha.1`).** Pre-built portable binaries for Windows and
> Linux are published on the [Releases page][releases]. They're *tested-build*
> (they compile + link in CI on every release) but real-world runtime testing is
> just beginning — expect sharp edges. You can always
> [build from source](#build-from-source) instead.
>
> **GPU drivers required.** Valenx renders with `wgpu`, so you need working GPU
> drivers: **Vulkan or DX12 on Windows, Vulkan or OpenGL on Linux**.

Grab the archive for your OS from the **[latest release][releases-latest]**:

| Platform | Download | Status |
| --- | --- | --- |
| **Windows** (x86_64) | `valenx-<ver>-x86_64-pc-windows-msvc.zip` | tested-build |
| **Linux** (x86_64) | `valenx-<ver>-x86_64-unknown-linux-gnu.tar.gz` | tested-build |
| **macOS** (Apple Silicon) | `valenx-<ver>-aarch64-apple-darwin.tar.gz` | **experimental** |

**Windows**

1. Download the `.zip` and unzip it.
2. Run `valenx.exe`.

> SmartScreen may warn — the binary is unsigned in this pre-alpha. Choose
> *More info → Run anyway*.

**Linux**

```sh
tar -xzf valenx-<ver>-x86_64-unknown-linux-gnu.tar.gz
cd valenx-<ver>-x86_64-unknown-linux-gnu
chmod +x valenx
./valenx
```

**macOS — experimental**

The macOS build **compiles in CI but the GUI is not yet verified on macOS**, so
treat it as experimental (and the archive may be absent if that best-effort build
didn't complete). To try it:

```sh
tar -xzf valenx-<ver>-aarch64-apple-darwin.tar.gz
cd valenx-<ver>-aarch64-apple-darwin
xattr -dr com.apple.quarantine valenx   # clear Gatekeeper quarantine (unsigned)
chmod +x valenx
./valenx
```

Each archive also bundles `LICENSE-MIT`, `LICENSE-APACHE`, and a short
`README.txt`. Signed native installers (`.msi` / `.dmg` / `.deb` / `.rpm`) are
the planned `1.0` distribution — see [docs/INSTALLER.md](./docs/INSTALLER.md).

[releases]: https://github.com/nochallenge/valenx/releases
[releases-latest]: https://github.com/nochallenge/valenx/releases/latest

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
- **[docs/CI.md](./docs/CI.md)** — CI policy (auto-runs on every push + PR; the `CI OK` gate)
- **[CONTRIBUTING.md](./CONTRIBUTING.md)** — how to contribute + dev setup
- **[SECURITY.md](./SECURITY.md)** — vulnerability disclosure
- **[CODE_OF_CONDUCT.md](./CODE_OF_CONDUCT.md)** — Contributor Covenant v2.1
- **[POLICIES.md](./POLICIES.md)** — SemVer, deprecation, MSRV, dependency rules
- **[rfcs/](./rfcs/)** — RFC process + accepted proposals

## What's actually working

| | |
|---|---|
| Native solvers | **150+ in-house workbenches** across aerospace, CFD/FEA, MD, quantum chemistry, biology, neuro, GR, optics, and CAD/CAM — written from scratch in Rust and **validated against published ground truth** (see [Validation](#validation)) |
| Optional adapters | **141** external open-source tools *also* integrate (CFD / FEA / MD / chemistry / protein / RNA / CRISPR / cryo-EM / workflows) for reference implementations + GPU-ML models — a convenience layer, not the core |
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
4. **Native-first.** The core science is our own from-scratch Rust solvers,
   validated against published ground truth — that's the product, not a wrapper.
5. **Optional interop, not dependency.** External tools integrate through
   optional adapters (reference implementations, GPU/ML models, not-yet-native
   domains) — a convenience layer you reach for only when you want it.
6. **AI-drivable first.** Every workbench exposes named accessible
   widgets + an agent-command bridge, so AI agents operate the app as
   first-class users — not just humans.

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
