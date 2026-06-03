# Valenx — 20 Year Roadmap

A native, open-source desktop simulation suite that unifies CAD, meshing,
CFD, FEA, EM, chemistry, molecular dynamics, battery modeling, robotics,
and multi-physics coupling into one download. Written in Rust. No browser.
No web server. No electron. Ships like FreeCAD does — you download an
installer, you run the app, it opens a window.

**Planning horizon:** 20 years. The capability inventory in Section 0.5
describes the target scope; matching it takes two decades even with a
talented team. The plan below paces that honestly.

**Supporting documents.** This roadmap sets direction; the "how" lives in:

- [ARCHITECTURE.md](./ARCHITECTURE.md) — how the system fits together
- [DESIGN.md](./DESIGN.md) — design plan (principles, system, screens, timeline)
- [LANGUAGES.md](./LANGUAGES.md) — language + library choices
- [TESTING.md](./TESTING.md) — how development + testing works
- [CONTRIBUTING.md](./CONTRIBUTING.md) — contributor workflow
- [MAINTAINERS.md](./MAINTAINERS.md) — who to reach
- [POLICIES.md](./POLICIES.md) — SemVer, deprecation, LTS, MSRV, dep rules
- [SECURITY.md](./SECURITY.md) — vulnerability disclosure
- [CODE_OF_CONDUCT.md](./CODE_OF_CONDUCT.md) — community standard
- [rfcs/](./rfcs/) — design RFCs. Initial set:
  [0001 project file format](./rfcs/0001-project-file-format.md),
  [0002 adapter contract](./rfcs/0002-adapter-contract.md),
  [0003 plugin API](./rfcs/0003-plugin-api.md)
- [CHANGELOG.md](./CHANGELOG.md) — what shipped when

---

## Decades at a glance

| Decade | What it looks like |
|---|---|
| **Year 0–5** | Native Rust app wrapping existing solvers. Fusion-360-quality UX. Users get everything the underlying OSS tools do, unified. First native Rust physics solvers (CFD, FEA basics) emerge. |
| **Year 5–10** | Native solvers reach validation parity with OpenFOAM/CalculiX/openEMS for core physics. Full CAD kernel + meshing native. Verticals mature (aerospace, automotive, biomed, civil). Industry adoption begins. |
| **Year 10–15** | Advanced multi-physics (FSI, rotating machinery, adjoint optimization) native. ML/AI integration (surrogate models, physics-informed NNs). Cloud + HPC mature. Plugin marketplace thrives. |
| **Year 15–20** | Feature parity with ANSYS Workbench across 5+ verticals. Industry certifications (ISO, FDA, DO-178C for safety-critical). Self-sustaining via consortium + support. Infrastructure for the next generation. |

---

## 0. Vision

### The pitch
One download that gives you **every scientific-computing capability the
integrated open-source tools offer** — CAD, meshing, CFD, FEA, EM,
chemistry, molecular dynamics, battery modeling, multi-physics — behind a
single beautifully-integrated native desktop application. Free, forever.

### Why this is worth 20 years of work
Every existing open-source scientific tool is world-class at its physics but
painful at its edges. FreeCAD ships a capable CAD kernel behind a 2003-era
UI. OpenFOAM is 20 years of validated CFD with no usable GUI. Code_Aster
is used by nuclear regulators but has a Salome interface from another
dimension. Buying ANSYS to escape this costs $50K/seat/yr. We fix the
problem by building the missing native integration layer — the thing all of
these tools share with none of them.

### What "done" looks like
A researcher, engineer, or hobbyist downloads one installer. Whatever they
want to simulate — a turbine blade, a battery cell, a burning jet, a drug
molecule, a wind-loaded building, a reacting flow in a catalytic reactor,
a cracked pressure vessel — **is a workflow inside Valenx**. They didn't
touch a shell, didn't edit a text file, didn't install Python packages. The
underlying solvers are the same peer-reviewed codes the industry trusts;
only the experience is new.

### UX reference: Fusion 360, not FreeCAD
This is the north star for how Valenx should feel:

- **Ribbon toolbar at the top**, tab-grouped by workflow (Geometry,
  Simulate, Mesh, Results, Render, Drawing, Manufacture)
- **Browser tree on the left** — bodies, components, sketches, features,
  assemblies, simulation studies
- **Timeline at the bottom** — scrubbable feature history, undo across time
- **Central viewport** with a proper ViewCube, orbit/pan/zoom gestures,
  smooth animated transitions, context-sensitive right-click menus
- **Command palette** (`Cmd/Ctrl+K`) for instant access to any operation
- **Properties panel on the right** that adapts to what's selected
- **Polished icons, animations, typography** — it should feel like a
  premium tool, not a lab-grown one
- **Multi-document tabs** — a project can have multiple CAD + simulation
  documents open at once
- **Light + dark themes**, both designed first-class
- **Measure overlays** that render in the viewport, not in a dialog

We are **not copying FreeCAD's docked-toolbox interface**. Fusion 360's
ribbon + timeline + ViewCube is the pattern worth stealing.

### Non-goals
- Not a browser app
- Not competing with ANSYS on marketing; competing on merit
- Not rewriting validated physics from scratch where OSS versions exist
- Not proprietary — Apache 2.0, forever
- Not a UI copy of FreeCAD (explicitly)

---

## 0.5 Capability inventory — what "fully inclusive" actually means

Every underlying OSS tool has decades of features. The integrated app must
expose **all of them** over time — not just the subset we implemented in
earlier iterations. This section is the scope target. Year-1 coverage is
"adapter-level" (subprocess wrappers expose the capability); long-term
coverage is "native-level" (rewritten into Rust with Valenx-native UX).

### Geometry / CAD (FreeCAD + OpenCASCADE)
- **Parametric solid modeling** — sketches with constraints, extrude, revolve, sweep, loft, helix
- **Boolean operations** — union, cut, intersect, compound, split
- **Feature operations** — fillet, chamfer, draft, shell, hole, pocket, rib, pattern (linear, polar, mirror), thicken
- **Surface modeling** — NURBS, lofted/swept surfaces, blends, filling, offset
- **Assemblies** — both constraint-based (assembly4) and joint-based (A2+), rigid + flexible
- **Sheet metal** — unfold, bend, flatten, punch, louver, seam
- **2D drafting** — dimensions, technical drawings, orthographic views, section views, detail views, bill of materials
- **Architecture / BIM** — walls, floors, roofs, windows, doors, structural elements, IFC import/export
- **Ship design** — hull forms, naval architecture workbench
- **CAM / toolpaths** — 2.5D / 3D milling, drilling, G-code output, post-processors for 20+ CNC dialects
- **Mesh editing** — STL cleanup, remeshing, decimation, boolean on tri meshes
- **Point cloud** — import, filtering, reconstruction
- **Import/export** — STEP, IGES, BREP, STL, OBJ, DXF, DWG, IFC, 3DS, X3D, DAE, VRML, PLY, Rhino 3DM, glTF, Collada
- **Parametric spreadsheets** — drive geometry from tables
- **Python scripting** — full API exposure
- **Addon ecosystem** — workbenches like Curves, Defeaturing, A2+, Fasteners, Gear, Symbols, Spreadsheet

### Meshing (gmsh + snappyHexMesh + cfMesh)
- **Structured hex** — blockMesh-style for simple geometries
- **Unstructured tet** — Delaunay, frontal, MeshAdapt algorithms
- **Unstructured hex** — snappyHexMesh-style, octree refinement
- **Quad/hex recombination** — for tet-dominant regions
- **Surface meshing** — conforming to STL surfaces
- **Boundary layers** — prism extrusion for wall-resolved CFD
- **Mixed meshes** — tet + prism + pyramid transitions
- **Adaptive refinement** — error-indicator-driven, solution-adaptive
- **Periodic meshes** — for turbomachinery, RVE analysis
- **Polyhedral meshes** — via polyDualMesh
- **Mesh quality metrics** — skewness, orthogonality, aspect ratio, determinant, volume ratio
- **Mesh morphing** — for shape optimization
- **Import/export** — .msh, .vtu, .foam, .unv, .cgns, .med, .nas, .inp

### CFD (OpenFOAM + SU2)
- **Incompressible**
  - Steady: `simpleFoam`, `SRFSimpleFoam` (single rotating frame), `porousSimpleFoam`
  - Transient: `pisoFoam`, `icoFoam`, `pimpleFoam`, `pimpleDyMFoam` (dynamic mesh)
- **Compressible**
  - Density-based: `rhoCentralFoam` (strong shocks), `rhoPimpleFoam`
  - Pressure-based: `rhoSimpleFoam`, `sonicFoam`, `sonicFoamAuto`
- **Turbulence**
  - **RANS**: k-ε (standard, RNG, realizable), k-ω (standard, SST, EARSM), v²-f, SpalartAllmaras, LRR, LaunderSharma, Chien, ~20 models total
  - **LES**: Smagorinsky, WALE, dynamic Smagorinsky, dynamic k-equation, Vreman, ~10 models
  - **Hybrid**: DES, DDES, IDDES, SAS
- **Multiphase**
  - VOF: `interFoam`, `compressibleInterFoam`, `interIsoFoam`
  - Euler-Euler: `twoPhaseEulerFoam`, `multiphaseEulerFoam`, `reactingTwoPhaseEulerFoam`
  - Lagrangian particles: `DPMFoam`, `MPPICFoam`, `coalChemistryFoam`
- **Combustion**
  - `reactingFoam`, `XiFoam`, `fireFoam`, `PDRFoam`, `chemFoam`
- **Heat transfer**
  - Boussinesq: `buoyantBoussinesqSimpleFoam`, `buoyantBoussinesqPimpleFoam`
  - Full compressible: `buoyantSimpleFoam`, `buoyantPimpleFoam`
  - Conjugate (CHT): `chtMultiRegionFoam` with solid regions
  - Radiation: P1, fvDOM, view factors
- **Rotating machinery**
  - MRF (multiple reference frame), AMI (arbitrary mesh interface), sliding mesh, overset mesh
- **Aeroacoustics** — FW-H surface noise, Curle's analogy, Lighthill tensor
- **Dynamic mesh / 6-DOF** — floating bodies, wing flutter, sloshing
- **Adjoint** — continuous + discrete for shape optimization
- **Electromagnetic coupling** — `mhdFoam` magnetohydrodynamics
- **Atmospheric boundary layer** — atmBL inlet, urban-wind templates
- **Boundary conditions** — 100+ types: cyclicAMI, overset, mapped, fan, totalPressure, waveTransmissive, ...
- **SU2 specifically** — continuous adjoint optimization at scale, multi-disciplinary FSI

### FEA (Code_Aster + CalculiX + Elmer)
- **Linear static** — stress, strain, deformation under load
- **Nonlinear static**
  - Geometric nonlinearity (large deformation, large rotation)
  - Material nonlinearity (plasticity — J2, Drucker-Prager, Mohr-Coulomb, Hill anisotropic; creep — Norton, Blackburn; viscoelasticity; hyperelasticity — Mooney-Rivlin, Neo-Hookean, Ogden, Gent)
  - Contact — penalty, Lagrangian, mortar, friction (Coulomb, regularized)
- **Dynamics**
  - Modal analysis (eigenfrequencies, mode shapes)
  - Harmonic / frequency response
  - Transient (implicit Newmark, explicit central difference)
  - Seismic / base excitation
  - Rotordynamics (unbalance, Campbell diagrams)
- **Thermal** — steady, transient, coupled thermal-mechanical
- **Fluid-structure interaction (FSI)** — via preCICE coupling
- **Fatigue** — Miner's rule, rainflow counting, Dang Van criterion, multi-axial
- **Fracture mechanics** — J-integral, G-theta, XFEM, cohesive zones
- **Acoustics** — Helmholtz, vibro-acoustic coupling
- **Piezoelectric / electro-mechanical coupling**
- **Composite layups** — ply-by-ply, failure criteria (Tsai-Wu, Puck, Hashin)
- **Buckling** — linear eigenvalue, nonlinear arc-length
- **Soil mechanics** — Mohr-Coulomb, Cam-Clay, consolidation
- **Pressure-vessel** — shakedown, limit analysis
- **Element library** — solid (tet, hex, wedge, pyramid), shell, plate, beam, truss, membrane, gap, spring
- **NAFEMS benchmark coverage** for validation

### Electromagnetics (openEMS + Elmer + SU2-Maxwell)
- **FDTD time-domain** — Maxwell solver (openEMS)
- **Antenna design** — dipoles, patch, horn, array; farfield radiation patterns; VSWR, S-parameters, efficiency
- **Microwave circuits** — filters, couplers, waveguides
- **Radar cross-section** — monostatic, bistatic, stealth analysis
- **SAR calculation** — body exposure, compliance with ICNIRP/IEEE C95.1
- **Periodic structures** — metamaterials, frequency-selective surfaces
- **Moving boundaries** — Doppler shift simulation
- **Electrostatics / magnetostatics** — Elmer
- **Time-harmonic** — eddy currents, coil design
- **Magneto-hydrodynamics** — OpenFOAM `mhdFoam`

### Chemistry (Cantera)
- **Equilibrium** — Gibbs free-energy minimization at (T,P), (H,P), (U,V), (T,V)
- **Kinetics** — homogeneous + heterogeneous (surface) reactions
- **Reactors**
  - 0-D: constant-P, constant-V, plug-flow, well-stirred (PSR, PFR, batch)
  - 1-D flames: premixed free, premixed burner-stabilized, counterflow diffusion, partially premixed
  - Shock tube / detonation
- **Transport** — mixture-averaged, multi-component, Soret/Dufour
- **Thermo** — NASA polynomials, Shomate, constant-Cp
- **Mechanism analysis** — sensitivity coefficients, reaction path analysis, eigenvalue flame decomposition
- **Mechanism library** — GRI-Mech 3.0, USC Mech II, Li/Dryer, AramcoMech, hundreds of published schemes
- **Coupling** — with OpenFOAM (`reactingFoam`) for CFD + kinetics

### Molecular dynamics (LAMMPS + GROMACS)
- **LAMMPS — materials + general MD**
  - Force fields: LJ, Coulomb, EAM, MEAM, Tersoff, REBO, AIREBO, ReaxFF, COMPASS, Buckingham, Stillinger-Weber, Morse, custom tabulated
  - Ensembles: NVE, NVT (Nose-Hoover, Berendsen, Langevin), NPT, NPH, isostress, isoenthalpy
  - Methods: MD, minimization, Monte Carlo, DPD, SPH, peridynamics, rigid-body, granular (DEM)
  - Accelerated: replica exchange, metadynamics, NEB (nudged elastic band), parallel tempering, hyperdynamics
  - Reactive: ReaxFF, COMB, bond/angle/dihedral creation + breaking
  - Coupling: QM/MM, ML potentials (DeepMD, HDNNP), SPPARKS kinetic MC
  - Boundary conditions: periodic, shrink-wrapped, mixed
  - Free energy: thermodynamic integration, umbrella sampling, metadynamics
  - Analysis: RDF, MSD, VACF, stress, energy decomposition
- **GROMACS — biomolecular MD**
  - Force fields: AMBER, CHARMM, GROMOS, OPLS, Martini (coarse-grained)
  - Solvation: TIP3P, TIP4P, SPC, SPC/E
  - Particle Mesh Ewald long-range electrostatics
  - Constraints: LINCS, SHAKE, SETTLE
  - Free energy: BAR, thermodynamic integration
  - Replica exchange
  - Protein folding, membrane dynamics, drug-target binding
  - GPU acceleration (CUDA, SYCL, OpenCL)

### Battery (PyBaMM)
- **Electrochemical models**
  - Doyle-Fuller-Newman (DFN) — full porous electrode
  - Single Particle Model (SPM)
  - Single Particle + Electrolyte (SPMe)
  - Multi-particle + electrolyte (SPMecc)
- **Degradation**
  - SEI layer growth
  - Lithium plating
  - Loss of active material (mechanical, stress-driven)
  - Electrolyte decomposition
  - Particle cracking (mechanical)
- **Pack modeling** — series/parallel cell assemblies, thermal coupling
- **Parameter estimation** — fit models to experimental data
- **Cycling protocols** — CC-CV, pulse, drive cycles, fast-charge, aging tests
- **Parameter sets** — LG M50, Chen2020, Ai2020, Ecker2015, Marquis2019, Mohtat2020, NCA_Kim2011, O'Kane2022, ORegan2022, Prada2013, Ramadass2004, Xu2019 + custom
- **Open-circuit voltage curves** — hysteresis, path-dependent

### Robotics + multibody (MuJoCo)
- **Rigid body dynamics** with continuous contact
- **Soft body** — cable, cloth, tendons
- **Articulated systems** — robot arms, humanoids, quadrupeds, grippers
- **Actuators** — position, velocity, torque, muscle models
- **Sensors** — joint, IMU, touch, rangefinder, camera
- **Inverse dynamics, inverse kinematics**
- **Differentiable simulation** — for RL + trajectory optimization
- **Domain randomization** hooks

### Multi-physics coupling (preCICE)
- **Fluid-structure interaction** — OpenFOAM ↔ CalculiX / Code_Aster / Elmer
- **Conjugate heat transfer** (beyond single-solver CHT)
- **Particle-fluid coupling**
- **Partitioned + monolithic coupling schemes**
- **Quasi-Newton acceleration** (IQN-ILS, IQN-IMVJ)

### Post-processing (VTK + ParaView + native)
- **Scalar / vector / tensor field rendering**
- **Iso-surfaces** with marching cubes/tetrahedra
- **Streamlines, pathlines, streaklines**
- **Vector glyphs** — arrows, cones, ellipsoids for tensors
- **Clipping planes** — scalar-defined + arbitrary planes
- **Slice planes** with on-slice field overlays
- **Contour plots** — 2D field cross-sections
- **Probe lines / points** — field extraction at arbitrary locations
- **Volume rendering** — transfer functions, ray-marched
- **Surface integration** — forces, heat flux, mass flow
- **Line integration** — circulation, work, path integrals
- **Statistics** — mean, variance, FFT, PSD, autocorrelation
- **Animation** — time-series playback, camera keyframes, export to MP4/GIF
- **Export** — PNG, EPS, SVG, glTF, USD, tabular CSV, VTU

### Optimization
- **Adjoint-based shape optimization** (OpenFOAM, SU2)
- **Design-of-experiments** — full factorial, Latin hypercube, OAT
- **Response surfaces** — polynomial, radial basis, Gaussian process
- **Multi-objective** — NSGA-II, SPEA2
- **Gradient-based** — SLSQP, L-BFGS
- **Topology optimization** — density-based, level-set
- **Uncertainty quantification** — polynomial chaos, Monte Carlo, Sobol

### HPC + distribution
- **Shared-memory parallel** — rayon for Rust-native solvers
- **MPI** — rsmpi wrapper for distributed solvers
- **HPC batch** — Slurm / PBS / LSF / SGE script generation + submission
- **Cloud dispatch** — AWS Batch, Azure Batch, Google Cloud, RunPod, Modal
- **Local multi-process** — job queue with `num_cpus`-aware scheduling
- **GPU offload** — where supported by the underlying solver (GROMACS, some OpenFOAM)

### Provenance + project management
- **Case history** — every run recorded in SQLite with inputs, outputs, hashes
- **Comparison** — diff two runs side-by-side (params + results)
- **Reproducibility** — every run emits a `.valenx` archive that regenerates it exactly
- **Tagging + notes** — organize runs by project, tag, keyword
- **Search** — full-text over case notes + parameters
- **AiiDA-style provenance graph** — every data node tracks origin
- **Report generation** — HTML + PDF with charts, tables, screenshots
- **Collaboration** — shared cases (local file + git integration)

### Scripting + automation
- **Python** — embedded interpreter (PyO3); full API exposure
- **Lua** — alternative (mlua); smaller footprint
- **Macros** — record user actions, replay as script
- **Plugin SDK** — WIT-based; community can add case types, post-processors, importers

### Sustainability / community
- Apache 2.0 everywhere
- Plugin marketplace
- Validation suite visible to all users
- mdBook documentation site
- Tutorial videos
- Built-in example gallery (50+ cases across domains)

---

## 0.6 Integrated tool registry

Every open-source tool Valenx integrates, grouped by domain. Licensing has
been verified against our policy:

- **Bundle**: license is permissive (BSD/MIT/Apache/zlib/public domain) —
  we ship binaries, link statically or dynamically, no contamination
- **Dynamic link**: LGPL tools — we link against their `.dll`/`.so`/`.dylib`
  at runtime; Valenx itself stays Apache 2.0
- **Subprocess**: GPL tools — we invoke the binary as a separate process;
  GPL does not cross process boundaries; Valenx stays Apache 2.0

### CAD / Geometry
| Tool | License | Mode | Phase | Purpose |
|---|---|---|---|---|
| **OpenCASCADE (OCCT)** | LGPL-2.1 | Dynamic link | 1 | Primary CAD kernel — BRep, STEP/IGES/BREP I/O |
| **FreeCAD (as reference)** | LGPL-2.1 | — | 1 | UX reference; not linked directly (we use OCCT directly) |
| **SolveSpace** | GPL-3 | Subprocess | 2 | 2D/3D constraint solver for sketcher |
| **OpenSCAD** | GPL-2 | Subprocess | later | Code-based CAD for scripting users |
| **Blender** | GPL-3 | Subprocess | later | 3D rendering + visualization (CLI mode) |
| **MeshLab** | GPL-3 | Subprocess | later | Triangle mesh cleanup + remeshing |
| **CGAL** | GPL-3 (dual) | Subprocess | later | Computational geometry algorithms |
| **fornjot / truck** | Apache/MIT | Bundle | later | Rust-native CAD kernels (consider contributing upstream) |

### Meshing
| Tool | License | Mode | Phase | Purpose |
|---|---|---|---|---|
| **gmsh** | GPL-2 | Subprocess / FFI | 1 | Primary mesher — tet, hex, quad |
| **NetGen** | LGPL-2.1 | Dynamic link | 2 | Alternative tet/prism mesher |
| **MMG** | LGPL (lib) / GPL (exe) | Mixed | 3 | Mesh adaptation + remeshing |
| **cfMesh** | GPL-3 | Subprocess | 1 | Polyhedral meshing (comes with OpenFOAM) |
| **snappyHexMesh** | GPL-3 | Subprocess | 1 | Refinement zones (comes with OpenFOAM) |
| **TetGen** | AGPL-3 | **Skip** | — | AGPL conflicts with our policy |
| **Triangle** | Academic-only | **Skip** | — | License restricts commercial use |

### CFD
| Tool | License | Mode | Phase | Purpose |
|---|---|---|---|---|
| **OpenFOAM** | GPL-3 | Subprocess | 1 | Primary CFD — incompressible + compressible + multiphase + combustion |
| **SU2** | LGPL-2.1 | Dynamic link | 1 | CFD + continuous adjoint optimization |
| **Code_Saturne** | GPL-3 | Subprocess | 2 | EDF's industrial CFD; alternative backend |
| **Nek5000 / NekRS** | BSD-3 | Bundle | later | High-order spectral-element CFD |
| **Basilisk** | GPL-3 | Subprocess | later | Adaptive-mesh-refinement CFD |
| **MFIX** | public domain | Bundle | later | Multiphase + DEM (fluidized beds) |
| **Gerris** | GPL-2 | Subprocess | later | Multigrid adaptive CFD |
| **Palabos / OpenLB** | GPL-3 | Subprocess | later | Lattice Boltzmann (alternative paradigm) |
| **SpectralDNS** | BSD-3 | Bundle | later | High-order DNS for research |
| **PeleC / PeleLM** | BSD-3 | Bundle | later | Combustion CFD (AMReX-based) |

### FEA / Mechanics
| Tool | License | Mode | Phase | Purpose |
|---|---|---|---|---|
| **Code_Aster** | GPL-3 | Subprocess | 1 | Industrial-grade nonlinear FEA (EDF France) |
| **CalculiX** | GPL-2 | Subprocess | 1 | Linear + nonlinear FEA (well-known) |
| **Elmer FEM** | LGPL-2.1+ | Dynamic link | 1 | Multiphysics FEM (CSC Finland) |
| **MFEM** | BSD-3 | Bundle | 2 | High-order FEM (LLNL, GPU-ready) |
| **FEniCS / DOLFINx** | LGPL-3 | Dynamic link | 2 | Symbolic FEM (automatic weak-form) |
| **deal.II** | LGPL-2.1 | Dynamic link | 2 | Mature FEM library |
| **scikit-fem** | BSD-3 | Bundle | 1 | Pure-Python FEM for quick cases |
| **Moose Framework** | LGPL-2.1 | Dynamic link | later | Multi-physics PDE (Idaho NL, reactor) |
| **NGSolve** | LGPL-2.1 | Dynamic link | later | High-order FEM |
| **Peridigm** | BSD-3 | Bundle | later | Peridynamics (Sandia) |

### Electromagnetics
| Tool | License | Mode | Phase | Purpose |
|---|---|---|---|---|
| **openEMS** | GPL-3 | Subprocess | 1 | FDTD time-domain EM |
| **Meep** | GPL-2 | Subprocess | 2 | MIT's FDTD (alternative) |
| **Palace** | BSD-3 | Bundle | 3 | AWS Annapurna modern Maxwell solver |
| **FDTD++** | GPL-3 | Subprocess | later | Alternative FDTD |

### Chemistry / Reactions
| Tool | License | Mode | Phase | Purpose |
|---|---|---|---|---|
| **Cantera** | BSD-3 | Bundle | 1 | Kinetics, equilibrium, flames, reactors |
| **RDKit** | BSD-3 | Bundle | 2 | Chemistry informatics, molecular properties |
| **Open Babel** | GPL-2 | Subprocess | 2 | Molecular format converter |
| **Quantum ESPRESSO** | GPL-2 | Subprocess | later | Ab-initio DFT, condensed matter |
| **CP2K** | GPL-2 | Subprocess | later | DFT + classical MD |
| **NWChem** | Apache-2.0 | Bundle | later | Quantum chemistry (PNNL) |
| **PSI4** | LGPL-3 | Dynamic link | later | Modern quantum chemistry |

### Molecular dynamics
| Tool | License | Mode | Phase | Purpose |
|---|---|---|---|---|
| **LAMMPS** | GPL-2 | Subprocess | 1 | Primary MD — materials, general |
| **GROMACS** | LGPL-2.1 | Dynamic link | 1 | Biomolecular MD |
| **OpenMM** | MIT | Bundle | 2 | GPU-ready MD (Stanford) |
| **ASE** (Atomic Simulation Environment) | LGPL-2.1 | Dynamic link | 2 | Glue for atomistic tools |
| **ESPResSo** | GPL-3 | Subprocess | later | Soft-matter MD (colloids, polymers) |
| **HOOMD-blue** | BSD-3 | Bundle | later | GPU MD for polymers |
| **MDAnalysis** | GPL-2 | Subprocess | later | Trajectory analysis |
| **MDTraj** | LGPL-2.1 | Dynamic link | later | Trajectory analysis |

### Battery / Electrochemistry
| Tool | License | Mode | Phase | Purpose |
|---|---|---|---|---|
| **PyBaMM** | BSD-3 | Bundle | 1 | DFN / SPM battery models |

### Robotics / Multibody
| Tool | License | Mode | Phase | Purpose |
|---|---|---|---|---|
| **MuJoCo** | Apache-2.0 | Bundle | 1 | Rigid + soft body dynamics |
| **Bullet Physics** | zlib | Bundle | later | Rigid-body / game physics |
| **Drake** | BSD-3 | Bundle | later | MIT/TRI robotics + optimization |
| **Gazebo / Ignition** | Apache-2.0 | Bundle | later | Robot simulator |
| **ROS 2** | Apache-2.0 | Bundle | later | Robot Operating System integration |

### Nuclear / Radiation (new vertical, Year 2+)
| Tool | License | Mode | Phase | Purpose |
|---|---|---|---|---|
| **OpenMC** | MIT | Bundle | 2 | Monte Carlo neutronics |
| **Shift** | BSD-3 | Bundle | later | ORNL reactor physics |

### Building / Energy (new vertical, Year 2+)
| Tool | License | Mode | Phase | Purpose |
|---|---|---|---|---|
| **EnergyPlus** | BSD-3 | Bundle | 2 | Building energy simulation (NREL) |
| **OpenStudio** | LGPL-3 | Dynamic link | 2 | Workflow layer over EnergyPlus |
| **Radiance** | BSD-like | Bundle | 2 | Daylighting / lighting simulation (LBNL) |

### Atmospheric / Climate
| Tool | License | Mode | Phase | Purpose |
|---|---|---|---|---|
| **WRF** | public domain | Bundle | later | Weather research & forecasting |
| **MPAS** | BSD-3 | Bundle | later | Global ocean + atmosphere |
| **Delft3D** | LGPL-3 | Dynamic link | later | Coastal hydrodynamics |
| **GEOS-Chem** | MIT | Bundle | later | Atmospheric chemistry |

### Seismology / Geosciences
| Tool | License | Mode | Phase | Purpose |
|---|---|---|---|---|
| **SPECFEM3D** | GPL-3 | Subprocess | later | Seismic wave propagation |
| **PyLith** | MIT | Bundle | later | Crustal deformation |

### Power systems
| Tool | License | Mode | Phase | Purpose |
|---|---|---|---|---|
| **OpenDSS** | BSD-3 | Bundle | later | Distribution system analysis |
| **GridLAB-D** | BSD-3 | Bundle | later | Power distribution |
| **PyPSA** | MIT | Bundle | later | Python power system analysis |

### Multi-physics coupling
| Tool | License | Mode | Phase | Purpose |
|---|---|---|---|---|
| **preCICE** | LGPL-3 | Dynamic link | 1 | FSI / CHT / particle-fluid coupling |

### Optimization
| Tool | License | Mode | Phase | Purpose |
|---|---|---|---|---|
| **Ipopt** | EPL-2.0 | Dynamic link | 2 | Nonlinear optimization |
| **NLopt** | LGPL-2.1 | Dynamic link | 2 | Optimization library |
| **CasADi** | LGPL-3 | Dynamic link | 3 | Nonlinear optimization (optimal control) |
| **Pagmo / pygmo** | LGPL-3 | Dynamic link | 3 | Parallel + multi-objective optimization |
| **Pyomo** | BSD-3 | Bundle | 3 | Optimization modeling language (Sandia) |

### Provenance / workflow
| Tool | License | Mode | Phase | Purpose |
|---|---|---|---|---|
| **AiiDA** | MIT | Bundle | 1 | Provenance graph |

### Visualization / Post
| Tool | License | Mode | Phase | Purpose |
|---|---|---|---|---|
| **VTK** | BSD-3 | Bundle | 1 | Data-structure + rendering backbone |
| **ParaView** | BSD-3 | Bundle | 1 | Results viewer (embedded or companion) |
| **VisIt** | BSD-3 | Bundle | later | LLNL's alternative to ParaView |
| **yt** | BSD-3 | Bundle | later | Volumetric astrophysics |

### Numerical infrastructure (used internally by our solvers + some tools above)
| Library | License | Mode | Phase | Purpose |
|---|---|---|---|---|
| **PETSc** | BSD-2 | Bundle | 2 | Sparse linear + nonlinear solvers |
| **SUNDIALS** | BSD-3 | Bundle | 2 | ODE / DAE integrators |
| **hypre** | Apache-2.0 | Bundle | 2 | Scalable multigrid solvers |
| **Kokkos** | Apache-2.0 | Bundle | 2 | GPU abstraction (performance portability) |
| **ParMETIS** | BSD-style | Bundle | 2 | Parallel graph partitioning |
| **SCOTCH** | CeCILL-C (≈LGPL) | Dynamic link | 2 | Graph partitioning |
| **HDF5** | BSD-like | Bundle | 1 | Hierarchical scientific data format |
| **NetCDF** | MIT-like | Bundle | 1 | Array-oriented data format |
| **CGNS** | zlib-like | Bundle | 1 | CFD data format standard |
| **Eigen** | MPL-2.0 | Bundle | 1 | Dense linear algebra C++ |
| **FFTW** | GPL-2 | Subprocess / replace with FFTS (BSD) | 2 | Fast Fourier transforms — use BSD alternative |
| **BLAS / LAPACK** | BSD-like | Bundle | 1 | Fundamental linear algebra |
| **OpenBLAS** | BSD-3 | Bundle | 1 | Optimized BLAS |

### ML / AI (Phase 13)
| Tool | License | Mode | Phase | Purpose |
|---|---|---|---|---|
| **PyTorch** | BSD-3 | Bundle | 13 | Deep learning, surrogate models |
| **JAX** | Apache-2.0 | Bundle | 13 | Differentiable simulation |
| **ONNX** | MIT | Bundle | 13 | Model interchange |
| **scikit-learn** | BSD-3 | Bundle | 13 | Classical ML |
| **DeepXDE** | LGPL-2.1 | Dynamic link | 13 | Physics-informed neural networks |

### Commercial / licensed — NOT included
| Tool | Reason for skip |
|---|---|
| ANSYS / Siemens NX / COMSOL / Abaqus | Proprietary commercial — we compete with these |
| NAMD | Non-commercial free only; license conflict |
| AMBER | Commercial license required |
| MCNP | Export-controlled (US government) |
| SERPENT | Not fully open |
| FLASH | Academic-only license |
| IsaacGym / IsaacLab | Closed source (NVIDIA) |
| TetGen | AGPL — viral license |
| Triangle | Academic-only |
| DL_POLY | Academic only |

---

**Total integrated tools:**
- Phase 1 (Year 1): **~25 tools** — all the previously-wrapped ones plus essential numerical infrastructure
- Phase 2-3 (Years 2-3): **+15 tools** — MFEM, FEniCS, deal.II, OpenMC, EnergyPlus, NetGen, MMG, OpenMM, ASE, RDKit, Palace, Meep, Code_Saturne, OpenStudio, Radiance
- Phase 4-10 (Years 4-10): **+20 tools** — alternative backends, atmospheric/climate/seismic/power verticals, advanced optimization, Blender/MeshLab/OpenSCAD integration
- Phase 11+ (Year 10+): **+10 tools** — long-tail specialty tools as verticals mature
- Phase 13 (Year 10+): **+5 ML/AI tools**

**~75 tools** total in the mature suite. Every one verified license-compatible with Valenx's Apache 2.0 distribution model.

---

## 1. Architecture

### The layer cake

```
┌─ valenx-app  ─ native window · egui · wgpu 3D ─────────────────┐
│  Tabbed workbench: Geometry · Mesh · Physics · Solve · Post    │
│  Docking panels · keyboard shortcuts · command palette         │
├─ valenx-core  ─ orchestration brain ───────────────────────────┤
│  Project / Case / Workflow data model                          │
│  DAG engine with incremental re-run (nix/bazel-like)           │
│  Provenance graph · SQLite run history · plugin loader         │
├─ valenx-geo  ─ CAD kernel ─────────────────────────────────────┤
│  OpenCASCADE via FFI (LGPL-safe) — primary forever             │
│  Optional: contribute to fornjot/truck as pure-Rust alternative│
├─ valenx-mesh  ─ meshing ───────────────────────────────────────┤
│  gmsh FFI + NetGen dynamic link + cfMesh/snappy subprocess     │
│  Native tet + hex + prism layer added over time (Phase 3)      │
├─ valenx-viz  ─ results viewer · wgpu ──────────────────────────┤
│  Iso-surfaces · streamlines · clipping planes · vector glyphs  │
│  Reads VTK · EnSight · native Valenx format                    │
├─ valenx-solvers  ─ native physics (gradually built) ───────────┤
│  cfd · fea · em · chem · md · battery                          │
├─ valenx-adapters  ─ subprocess + FFI wrappers for 75+ tools ───┤
│  Year 1:   OpenFOAM · FreeCAD/OCCT · Code_Aster · CalculiX ·   │
│            Elmer · Cantera · SU2 · openEMS · LAMMPS · GROMACS ·│
│            PyBaMM · gmsh · preCICE · MuJoCo · scikit-fem · ... │
│  Year 2-3: MFEM · FEniCS · deal.II · OpenMC · EnergyPlus ·     │
│            NetGen · MMG · OpenMM · ASE · RDKit · Meep · Palace │
│  Year 4-10: 20+ more — climate, seismic, power, nuclear,       │
│            advanced optimization, alternative backends         │
│  Full registry in Section 0.6                                  │
└─ Bundled solver binaries ─ shipped per version manifest ───────┘
```

### Why this shape works

**Adapters first, native second.** Year 1 ships a product by wrapping
existing solvers. Users get something real. Native Rust solvers replace
adapters one at a time over years 3-8.

**Subprocess firewall for GPL solvers.** OpenFOAM, Code_Aster, CalculiX,
openEMS, LAMMPS, gmsh are GPL. We invoke them as external processes; no
static linking. That keeps Valenx Apache 2.0 regardless of what we bundle.

**LGPL tools can be linked.** FreeCAD/OCCT, Elmer, SU2, GROMACS, preCICE
ship as dynamic libraries we bind via FFI. Faster, cleaner integration.

**BSD/Apache tools are a free lunch.** Cantera, PyBaMM, scikit-fem, MuJoCo,
OpenMC, EnergyPlus, MFEM, and many others — bundle them however we want,
no strings.

### Cargo workspace layout

```
valenx/
├── Cargo.toml                workspace root
├── crates/
│   ├── valenx-app            native UI entrypoint
│   ├── valenx-core           project model · DAG · workflow
│   ├── valenx-geo            CAD kernel / OCCT bindings
│   ├── valenx-mesh           meshing
│   ├── valenx-viz            3D results viewer (wgpu)
│   ├── valenx-cfd            native CFD
│   ├── valenx-fea            native FEA
│   ├── valenx-em             native EM
│   ├── valenx-chem           native chemistry
│   ├── valenx-md             native molecular dynamics
│   ├── valenx-battery        native battery models
│   ├── valenx-scripting      PyO3 + Lua bindings
│   ├── valenx-plugin-api     WIT for third-party plugins
│   ├── valenx-bench          validation + benchmarks
│   └── valenx-adapters/
│       ├── openfoam
│       ├── freecad
│       ├── calculix
│       ├── elmer
│       ├── cantera
│       ├── pybamm
│       ├── gmsh
│       ├── su2
│       ├── openems
│       ├── lammps
│       ├── gromacs
│       └── precice
├── docs/                     mdBook
├── assets/                   icons · themes · default cases
├── installers/               msi · pkg · appimage configs
├── third-party/              vendored LGPL/BSD tool source (build from)
└── tests/                    integration + regression
```

---

## 2. The 20-year phase plan

Phases overlap — later phases start while earlier ones continue. Dates are
calendar-year from project start.

### Phase 0 — Foundation · months 0–6
Lock in decisions the next 20 years won't regret.

- License (Apache 2.0), governance, RFC process, domain, GitHub org
- CI on Win/macOS/Linux; nightly release pipeline
- Cargo workspace skeleton; pre-commit; MSRV policy (Rust stable)
- Pick UI framework (recommend egui; alternative Slint)
- Pick 3D backend (wgpu via egui or directly)
- Docs site (mdBook)
- Contributor guide, code of conduct
- Define the `.valenx` project file format (TOML-based)
- RFC-0001: Architecture overview
- RFC-0002: Plugin interface (WIT)
- RFC-0003: Adapter contract

**Deliverable:** `valenx 0.0.1` — empty window that says "Valenx" and closes cleanly.

### Phase 1 — Shell + adapters · months 3–18
The native app that matches everything the old web Valenx did.

- `valenx-app` main window, Fusion-360-style ribbon toolbar, docking panels,
  case tree, properties panel, timeline (feature history)
- ViewCube, orbit/pan/zoom gestures, smooth animated camera transitions
- Command palette (`Cmd/Ctrl+K`)
- Light + dark themes, both first-class
- `valenx-viz` baseline — STL, glTF, point cloud, axis helpers
- `valenx-core` project file I/O, workflow DAG, provenance
- Port every adapter from the legacy Python code to Rust
- Case templates: airfoil · cavity · channel · BFS · cylinder · heat
  transfer · multiphase · compressible · MRF · CHT · 3D wing · reactor · lung · more
- Validation suite: 8+ benchmarks with graded A/B/C scorecards
- Run history DB (SQLite) + comparison view
- `.valenx` archive export/import
- HPC dispatch (Slurm/PBS/LSF scripts)
- PDF/HTML report generator
- Installer: MSI / .app / AppImage (signed)
- Multi-document tabs

**Deliverable:** `valenx 0.1` — native app, downloadable, does everything
the legacy Valenx did, looks like Fusion 360 feels. Bundled size ~1 GB.

### Phase 2 — Geometry kernel · months 9–36
Replace FreeCAD as the CAD frontend.

- Architecture RFC: half-edge BRep data structures, serde-serializable
- OCCT FFI bindings crate (rewrite `opencascade-rs` if needed)
- Sketcher: 2D primitives, constraint solver (SolveSpace-style)
- Feature tree: extrude, revolve, sweep, loft, helix, fillet, chamfer, draft, shell, hole, pocket, rib, pattern (linear/polar/mirror), thicken
- CSG ops (union, subtract, intersect, compound, split) via OCCT initially, native later
- Surface modeling: NURBS, lofted/swept, blends, filling, offset
- Assembly: constraint-based + joint-based
- Parametric re-execution on parameter change
- STEP / IGES / BREP / STL / OBJ / DXF / DWG / IFC / glTF / USD import/export
- 100+ example parts that work end-to-end
- Evaluate fornjot and truck — prefer adopting/contributing over rewriting

**Deliverable:** `valenx-geo 0.1` — parametric parts and assemblies in a
native Rust UI. Feature-complete enough for 80% of engineering CAD.

### Phase 3 — Meshing · months 12–36
Replace gmsh + snappyHexMesh as the mesher.

- Structured hex meshing (blockMesh-equivalent)
- Surface-conforming tet meshing (bundle gmsh FFI short-term, native mid-term)
- Unstructured hex, octree refinement
- Prism boundary layer extrusion
- Quad/hex recombination
- Polyhedral meshes (polyDualMesh)
- Interactive refinement zone placement in the 3D viewport
- Quality metrics: skewness, orthogonality, aspect ratio heatmaps
- Adaptive mesh refinement API (error-indicator-driven)
- Periodic meshes, mesh morphing
- Import/export: .msh, .vtu, .foam, .unv, .cgns, .med, .nas, .inp

**Deliverable:** `valenx-mesh 0.1` — drag-to-mesh workflow that beats both
snappyHexMesh and gmsh in UX.

### Phase 4 — Native CFD basics · months 18–48
Incompressible RANS first; prove the architecture.

- Finite-volume framework in `valenx-cfd`
- Incompressible Navier-Stokes: SIMPLE, PISO, PIMPLE
- Turbulence: k-ε, k-ω SST, SpalartAllmaras (RANS)
- Parallelization: rayon (shared memory) + rsmpi (MPI)
- Published-benchmark regression suite (Ghia cavity, Driver-Seegmiller, NACA Mason, periodic hill)

**Deliverable:** `valenx-cfd 0.1` — incompressible RANS graded A/B on 6+ benchmarks, within 5% of OpenFOAM accuracy.

### Phase 5 — Native CFD advanced · months 30–72
Compressible + multiphase + combustion + rotating machinery.

- Density-based compressible (rhoCentralFoam-equivalent, strong shocks)
- Pressure-based compressible (low-subsonic to transonic)
- LES: Smagorinsky, WALE, dynamic-k, Vreman
- Hybrid RANS-LES: DES, DDES, IDDES, SAS
- VOF multiphase + Euler-Euler + Lagrangian particles
- Combustion: reactingFoam-equivalent, coupled with Cantera mechanisms
- Conjugate heat transfer (native CHT, no preCICE needed)
- MRF + AMI + sliding mesh + overset (rotating machinery, turbomachinery)
- Dynamic mesh + 6-DOF
- Aeroacoustics (FW-H, Curle's, Lighthill tensor)
- Atmospheric BL, urban wind templates
- 100+ boundary condition types

**Deliverable:** `valenx-cfd 0.5` — parity with ~80% of OpenFOAM solver families.

### Phase 6 — Native FEA basics · months 24–60
Linear + nonlinear static, basic dynamics.

- Linear elasticity (tet + hex + wedge + pyramid elements, Lagrange shape functions)
- Geometric nonlinearity (large deformation, large rotation)
- Material nonlinearity: J2 plasticity, Drucker-Prager, Mohr-Coulomb, Hill anisotropic
- Creep: Norton, Blackburn
- Viscoelasticity: generalized Maxwell, Prony series
- Hyperelasticity: Mooney-Rivlin, Neo-Hookean, Ogden, Gent
- Contact: penalty, Lagrangian, mortar, Coulomb friction
- Transient dynamics: implicit Newmark, explicit central difference
- Thermal steady + transient
- Thermomechanical coupling
- Modal analysis, harmonic response, buckling
- Element library: solid, shell, plate, beam, truss, membrane, gap, spring
- NAFEMS benchmark validation

**Deliverable:** `valenx-fea 0.1` — covers ~90% of CalculiX; ~60% of Code_Aster.

### Phase 7 — Native FEA advanced · months 48–96
Where Code_Aster really earns its reputation.

- Fatigue: Miner's rule, rainflow counting, Dang Van, multi-axial
- Fracture mechanics: J-integral, G-theta, XFEM, cohesive zones
- Rotordynamics: unbalance, Campbell diagrams, critical speeds
- Acoustics: Helmholtz, vibro-acoustic coupling
- Piezoelectric / electro-mechanical coupling
- Composite layups: ply-by-ply, Tsai-Wu, Puck, Hashin failure
- Seismic / base excitation, response spectrum
- Soil mechanics: Mohr-Coulomb, Cam-Clay, consolidation
- Pressure-vessel: shakedown, limit analysis

**Deliverable:** `valenx-fea 0.5` — ~80% of Code_Aster capability.

### Phase 8 — Other physics natives · months 30–84
Parallel tracks, each led ideally by a domain expert.

- `valenx-em` FDTD (Yee grid, PML boundaries) — replaces openEMS
  - Antenna design, radar cross-section, SAR, metamaterials, moving boundaries
- `valenx-chem` kinetic mechanisms + equilibrium — replaces Cantera
  - 0-D reactors, 1-D flames, shock tubes, sensitivity analysis, coupling to CFD
- `valenx-battery` DFN + SPM + SPMe — replaces PyBaMM
  - Degradation models (SEI, plating, LAM, cracking), pack modeling
- `valenx-md` basic MD — complements LAMMPS
  - Verlet + LJ + Coulomb + EAM; ensembles NVE/NVT/NPT; minimization
- `valenx-robot` — rigid body dynamics with contact (MuJoCo-inspired)

**Deliverable:** Five new solvers at 0.1-0.3 capability each — enough to
handle common cases in each domain natively.

### Phase 9 — Advanced optimization + UQ · months 72–132
Stuff commercial tools charge extra for.

- Native adjoint optimization in CFD and FEA
- Topology optimization (density-based, level-set)
- Multi-objective: NSGA-II, SPEA2, MOEA/D
- Gradient-based: SLSQP, L-BFGS, trust-region
- Response surfaces: polynomial, RBF, Gaussian process, neural net
- Design of experiments: full factorial, Latin hypercube, Taguchi, OAT
- Uncertainty quantification: polynomial chaos, Monte Carlo, Sobol, MLMC
- Robust design, reliability-based design optimization (RBDO)
- Adjoint-based sensitivities for shape + topology
- Surrogate modeling integration

**Deliverable:** `valenx-opt 0.1` — optimization and UQ exposed as first-class
workflows in the UI.

### Phase 10 — Advanced multi-physics · months 84–156
Where we stop needing preCICE.

- Native FSI (partitioned quasi-Newton, monolithic for small problems)
- Native multi-region CHT (beyond single-solver CHT)
- Native particle-fluid (fluid + DEM, fluid + Lagrangian)
- Native reactive flow (CFD + chemistry integrated, not subprocess-coupled)
- Electromagneto-structural (coil deformation under Lorentz force)
- Thermo-hydraulic-mechanical (for nuclear, geothermal)
- Electrochemistry + fluid (for flow batteries, electrolyzers)
- Multi-phase + phase-change (boiling, condensation, freezing)

**Deliverable:** `valenx-coupling 0.1` — most multi-physics cases run natively
without external coupling libraries.

### Phase 11 — Verticals maturity · years 7–14
Each vertical gets curated depth.

- **Aerospace**: airfoils, wings, propellers, rockets, jet engines, stability derivatives, flutter
- **Automotive**: external aero, underhood thermal, HVAC, crash (explicit dynamics), NVH
- **Biomedical**: blood flow, lung airflow, drug delivery, orthopedic FEA, medical device CFD
- **Civil / architecture**: urban wind, tall-building loads, HVAC, BIM integration, structural steel
- **Energy**: wind turbines, solar panels, battery packs, fuel cells, nuclear thermal-hydraulics
- **Defense**: radar cross-section, blast modeling, ballistic impact (non-export-controlled)
- **Manufacturing**: welding thermal, additive manufacturing thermal-mechanical, forming, casting
- **Electronics**: CHT for PCB/chip cooling, EMC/EMI, solder-joint fatigue
- **Marine**: hull resistance, cavitation, propeller design, VOF seakeeping

Each vertical gets:
- Packaged templates + example gallery
- Reference benchmarks validated against published data
- Industry-partner validation (where possible)
- Specialized BC/material libraries
- Publications documenting the validation

**Deliverable:** `valenx-verticals` — each domain has a polished workbench
with 10+ example cases.

### Phase 12 — Cloud + HPC + collaboration · years 8–15
From single-workstation to infrastructure.

- Native HPC dispatch: Slurm, PBS, LSF, SGE with result streaming
- Cloud compute: AWS Batch, Azure Batch, GCP Batch, RunPod, Modal
  - Transparent pricing display, cost estimates per case
  - Spot-instance handling, checkpoint-restart
- GPU offload: CUDA, ROCm, Metal — where solver supports it
- Collaboration: shared cases via git-LFS, comments, review workflow
- Live collaborative editing (optional) — CRDT-based
- Institutional deployments: site license, LDAP/SSO, audit logging
- Artifact storage: S3-compatible backends for large result archives

**Deliverable:** Valenx runs anywhere from laptop to 10k-core cluster with
the same UI.

### Phase 13 — ML/AI integration · years 10–17
The biggest differentiator from classical tools.

- Surrogate models: train neural nets on CFD/FEA outputs, use for real-time
  preview + optimization warm-start
- Physics-informed neural networks (PINNs) for PDE approximation
- ML-accelerated meshing (learn from good meshes)
- Foundation models for simulation (e.g., ALPHAFOLD-style for materials, turbulence)
- Anomaly detection in results (auto-flag non-physical outputs)
- Reinforcement learning for shape + topology optimization
- Generative CAD (text → 3D model via diffusion or next-generation)
- Natural-language case setup (LLM-driven — builds on existing pattern matcher)
- Auto-validation: compare new runs against historical database, flag drift
- Differentiable simulation for end-to-end gradient-based design

**Deliverable:** `valenx-ml` — optional but well-integrated AI-assisted workflows.

### Phase 14 — Ecosystem maturity · years 12–18
The platform becomes infrastructure.

- Plugin marketplace (1000s of third-party plugins, curated + reviewed)
- Certification path: ISO 9001, ISO 27001
- Safety-critical certifications where relevant (FDA for medical CFD,
  DO-178C / DO-254 for avionics, IEC 61508 for industrial control)
- Teaching infrastructure: university course curricula, reference textbooks,
  built-in pedagogical modes
- Commercial support ecosystem: multiple vendors offering Valenx consulting
- Annual conference (ValenxCon — physical + virtual)
- Language translations: docs + UI in Spanish, Chinese, French, German,
  Japanese, Hindi, Portuguese
- Paper count: 500+ peer-reviewed publications citing Valenx in methods

**Deliverable:** Valenx is known as "the FreeCAD of everything" — the
obvious OSS default for scientific computing.

### Phase 15 — Infrastructure for the next generation · years 15–20+
The project outlives its original maintainers.

- Apache Software Foundation or Eclipse Foundation home (if not already)
- Endowed development positions at partner universities
- Long-term-stable API (backward compatibility for 10+ years)
- Immutable reproducibility — a Valenx case from 2026 runs identically in 2046
- Archival of historical solver versions (everyone can re-run old papers)
- Deprecation policy for adapters as native solvers fully supplant
- Student + junior-developer onboarding program
- Valenx-as-a-library: embed it in commercial products via well-defined API

**Deliverable:** Valenx continues without any single individual; the tooling
outlives the founders.

---

## 🎯 Bio ecosystem · COMPLETE ✅ + Phases 43 / 44.5 / 35.5 / 35.6 / 45

The bio adapter family below (Phase 5.5 / 5.6 / 5.7 / 17 → 45)
**completes the bio ecosystem from the user's original /review
list** as of the Phase 22.5 + Phase 42 pair — every major category
called out is now covered — **plus Phase 43 layers mRNA / vaccine
therapeutic design on top, Phase 44.5 sister-expands the Phase 28
RNA secondary-structure trio with mfold + EternaFold + LinearFold,
Phase 35.5 sister-expands the Phase 35 CRISPR design trio with the
Liu and Komor labs' base / prime editing tools, Phase 35.6 closes
the design → predict-outcome → off-target loop with inDelphi /
FORECasT / AlphaMissense / CRISPRitz, and Phase 45 opens two new
domains — the first PK/PD pharmacokinetics modeling category in
Valenx (PK-Sim) plus the first RNA tertiary 3D structure
prediction category in Valenx (SimRNA)**. The bio surface now
spans **121 bio adapters across 43 biology / biotech / chemistry
phases**, covering alignment / base + prime editing /
cheminformatics / CRISPR / cryo-EM / DNA geometry / docking /
edit-outcome prediction / MD analysis / MD engines / microscopy /
mRNA design / pharmacokinetics / phylogenetics / population
genetics / protein design / quantum chemistry / RNA structure (2D
+ 3D) / sequence editors / sequence read simulators / single-cell
genomics / spatial stochastic simulation / structure prediction /
structure search / synthetic biology / systems biology / variant
calling / viewers (desktop + web) / web visualization / workflow
managers — all in one Valenx shell with no glue code beyond the
existing case-toml / prepare / run / collect path.

---

### Phase 5.5 — MD analysis expansion · live
Sister-adapter expansion of the existing Phase 17 MDAnalysis adapter.
Round out the post-MD analysis surface that Phase 17 MDAnalysis
opened with three more established open-source tools that span the
post-MD analysis tradeoff space — enhanced-sampling collective-
variable evaluation + free-energy reweighting (PLUMED, the de-facto
plug-in that wraps every major MD engine for biased-simulation /
reweighting work; LGPL-3.0; defines collective variables (RMSD,
dihedrals, distances, contact maps), biases (metadynamics, well-
tempered metad, umbrella sampling, ABF), and a reweighting
framework that turns biased trajectories back into unbiased free-
energy surfaces; the `plumed driver` sub-command runs PLUMED
standalone over a pre-computed trajectory: read frames, evaluate
the collective variables defined in `plumed.dat`, write COLVAR /
bias / HILLS files; single-binary subprocess shape sister to Phase
18 BWA), protein-dynamics elastic-network / normal-mode analysis
(ProDy, the canonical Python toolkit for ENM / GNM / ANM and
ensemble PCA; MIT; ships elastic-network models, normal-mode
analysis, ensemble PCA, the NMD trajectory format consumed by
VMD's NMWiz plug-in, and integrations with the BLAST / DALI / PDB
databases; Python-script subprocess shape sister to Phase 17
Biopython), and canonical AmberTools trajectory analysis via
cpptraj's domain language (cpptraj, the reference workhorse for
`rms` / `radgyr` / `hbond` / `clustering` over Amber-format
trajectories; GPL-3.0; reads Amber `.prmtop` / `.parm7` topologies
plus `.nc` / `.dcd` / `.mdcrd` trajectories, runs an analysis
script authored in cpptraj's domain language, and writes results
into the workdir as `.dat` per-frame tables, `.agr` XmGrace plot
data, or `.gnu` gnuplot scripts; single-binary subprocess shape
sister to PLUMED). PLUMED + cpptraj follow the established Phase
18 BWA single-binary CLI pattern: trajectory + script in, analysis
tables out. ProDy follows the established Phase 17 Biopython
Python-script subprocess shape: the user supplies a Python script
that imports `prody` and reads `valenx_params.json` for the parsed
knobs. Phase 5.5 sits numerically adjacent to the original Phase 5
MD beachhead and ships chronologically right after Phase 39 DNA
structural geometry — same chronological-vs-numerical convention
used for Phase 17.5 / 24 / 28 / 31 / 35 / 39.

- 3 new adapter crates land under `crates/valenx-adapters/bio/`:
  PLUMED (single-binary subprocess sister to Phase 18 BWA;
  `case.toml` knobs `plumed_dat` (PLUMED input file describing the
  collective variables and bias to compute; required), `trajectory`
  (XTC trajectory; required — users running DCD / TRR can swap to
  `--mf_dcd` / `--mf_trr` via `extra_args`), `output_basename`
  (filename stem PLUMED uses for COLVAR / bias outputs; required,
  non-empty), `kt` (`f64`, > 0.0 and finite; PLUMED's `k_B T` in
  its energy units — kJ/mol by default; default 2.494 = room
  temperature 300 K; a zero or NaN `kt` would crash PLUMED's
  reweighting on the first frame), `extra_args`; `prepare()`
  resolves both paths against the case directory when relative,
  validates each file exists on disk (returns `InvalidCase` with a
  helpful message when missing), and composes `plumed driver
  --plumed <plumed_dat> --mf_xtc <trajectory> --kt <kt>
  [extras...]`; `run()` streams PLUMED's `PLUMED: PLUMED is
  starting` startup banner / periodic per-frame status / `PLUMED:
  Finishing` end-of-run sentinels into progress hints; collects
  `<output_basename>*.dat` (`Tabular`, "PLUMED COLVAR output") and
  `<output_basename>*.bias` (`Tabular`, "PLUMED bias"); probe via
  `find_on_path(&["plumed"])`; LGPL-3.0 licensed; version range
  `2.9.0..3.0.0` (PLUMED 2.9 (2023) is the modern stable line —
  the `driver` sub-command, the metadynamics / OPES bias family,
  and the Python interface are all mature); `bio.plumed.analyze`
  ribbon capability), ProDy (Python-script subprocess sister to
  Phase 17 Biopython; `case.toml` knobs `script` (path to user-
  supplied Python script; required), `python` (interpreter name;
  default `"python3"`), `input_pdb` (input PDB; required),
  `output_basename` (filename stem ProDy uses for ENM / mode / NMD
  outputs; required, non-empty), `num_modes` (`u32`, ≥ 1; number
  of normal modes to compute; default 20), `cutoff` (`f64`, > 0.0
  and finite; ENM contact cutoff in Å; default 15.0); `prepare()`
  stages script + input PDB into the workdir under their original
  filenames so the script can resolve them via relative paths,
  then writes a flat `valenx_params.json` containing `input_pdb`
  (staged filename), `output_basename`, `num_modes`, and `cutoff`;
  collects `<output_basename>*.npz` (`Native`, "ProDy ENM modes"),
  `<output_basename>*.nmd` (`Native`, "ProDy NMD trajectory" — the
  NMD format consumed by VMD's NMWiz plug-in for normal-mode
  visualisation), and `<output_basename>*.csv` (`Tabular`, "ProDy
  table"); probe via Python on PATH with an `import prody` check
  (returns `ok = true` with a warning when import fails so non-
  standard installs aren't blocked); MIT licensed; version range
  `2.4.0..3.0.0` (ProDy 2.x is the modern stable line; 2.4 is the
  floor we test against); `bio.prody.analyze` ribbon capability),
  and cpptraj (single-binary subprocess sister to PLUMED;
  `case.toml` knobs `script` (`.ptraj` / `.cpptraj` analysis
  script; required), `topology` (Amber `.prmtop` / `.parm7`;
  required), `extra_args`; `prepare()` resolves both paths against
  the case directory when relative, validates each file exists on
  disk, and composes `cpptraj -p <topology> -i <script>
  [extras...]`; collects `*.dat` (`Tabular`, "cpptraj analysis
  output"), `*.agr` (`Tabular`, "cpptraj XmGrace plot"), and
  `*.gnu` (`Log`, "cpptraj gnuplot script"); probe via
  `find_on_path(&["cpptraj"])`; GPL-3.0 licensed; version range
  `6.0.0..7.0.0` (cpptraj 6.x is the modern stable line shipped
  with AmberTools 23+ (2023)); `bio.cpptraj.analyze` ribbon
  capability). Each wired into `valenx-app::init_registry`.
- No new canonical types, no new format readers, no new CLIs. All
  three adapters consume user-supplied inputs (PLUMED collective-
  variable scripts + MD trajectories, ProDy Python analysis
  scripts + input PDBs, cpptraj domain-language scripts + Amber
  topologies + trajectories) and emit user-readable artifacts
  (PLUMED COLVAR `.dat` / `.bias` tables, ProDy `.npz` ENM-mode
  arrays + `.nmd` NMD-format trajectories + `.csv` analysis
  tables, cpptraj `.dat` per-frame tables + `.agr` XmGrace plots
  + `.gnu` gnuplot scripts) that the unchanged `Results.artifacts`
  collection model surfaces directly. The existing
  `valenx_bio::format::pdb` reader inspects collected PDB inputs
  for chain / residue / atom counts. A first-class MD-analysis
  canonical type — a typed collective-variable / normal-mode /
  per-frame-statistics representation spanning all three back-ends
  and the existing Phase 17 MDAnalysis adapter — defers to a
  future phase along with COLVAR plotters, normal-mode
  visualizers, and per-statistic time-series viewers.
- 3 new `valenx-init` templates ship: `plumed` (`plumed-analyze`),
  `prody` (`prody-analyze`), and `cpptraj` (`cpptraj-analyze`).
  Cross-binary roundtrip test sweeps all 96 templates clean.

The full per-phase shape lives in
[docs/src/phases/phase-5-5-md-analysis.md](./docs/src/phases/phase-5-5-md-analysis.md);
the implementation plan lives at
[docs/superpowers/plans/2026-04-30-md-analysis.md](./docs/superpowers/plans/2026-04-30-md-analysis.md).

**Deliverable:** Adapter inventory at **100 of 101** fully live
after this phase (3 new MD-analysis-expansion adapters added
alongside the Phase 33 synthetic-biology trio on top of the Phase
17 + Phase 17.5 + Phase 18 + Phase 18.5 + Phase 18.6 + Phase 19 +
Phase 19.5 + Phase 20 + Phase 22 + Phase 23 + Phase 24 + Phase 25
+ Phase 27 + Phase 27.5 + Phase 27.6 + Phase 28 + Phase 29 + Phase
30 + Phase 30.5 + Phase 31 + Phase 32 + Phase 34 + Phase 35 +
Phase 36 + Phase 38 + Phase 39 totals); Phase 5.5 rounds out the
post-MD analysis surface that Phase 17 MDAnalysis opened,
broadening the simulate-MD → analyze-trajectory loop to four
post-MD analysis tools (the Phase 17 MDAnalysis adapter plus
PLUMED, ProDy, cpptraj) feeding into the existing Phase 5 GROMACS
/ LAMMPS MD engines and the entire Phase 17 → 39 biology /
biotech / chemistry expansion. **Crosses the 100-adapter
milestone** alongside the Phase 33 synthetic-biology trio.

### Phase 5.6 — Bio MD engines · live
Sister-domain expansion of the existing Phase 5 GROMACS / LAMMPS
MD engine beachhead. Round out the all-atom + GPU-native MD-engine
surface with three more established open-source engines that span
the corners GROMACS / LAMMPS / OpenMM don't reach — the canonical
UIUC academic NAMD all-atom engine (NAMD, the de-facto choice in
biomolecular MD pedagogy and a workhorse on every academic HPC
cluster; custom NAMD-License — academic / non-commercial use only,
flagged via mandatory `"academic"`-keyworded probe warning), the
AmberTools OSS portion of AMBER's MD engine (`sander`, the OSS
sibling of cpptraj that Phase 5.5 already wraps; GPL-3.0), and
the Glotzer-lab GPU-native particle simulator (HOOMD-blue, the
canonical Python-scripted GPU-first engine for soft-matter /
coarse-grained particle systems; BSD-3-Clause). NAMD + sander
follow the established Phase 18 BWA single-binary CLI pattern;
HOOMD-blue follows the established Phase 17 OpenMM Python-script
subprocess shape. Phase 5.6 sits numerically adjacent to Phase 5.5
and ships chronologically right after Phase 17.7 — same
chronological-vs-numerical convention used for Phase 17.5 / 24 /
28 / 31 / 35 / 39 / 5.5.

- 3 new adapter crates land under `crates/valenx-adapters/bio/`:
  NAMD (single-binary subprocess sister to Phase 5 LAMMPS / GROMACS;
  `case.toml` knobs `config` (NAMD `.namd` / `.conf` configuration;
  required), `processors` (`u32`, default 1; emitted as a single
  OsString `+p<N>`), `extra_args`; `prepare()` resolves `config`
  against the case directory when relative and composes
  `<binary> +p<N> <config> [extras...]` where `<binary>` is `namd2`
  or `namd3`; collects `*.dcd` (`Native`, "NAMD trajectory (DCD)"),
  `*.coor` (`Native`, "NAMD coordinates"), `*.vel` (`Native`,
  "NAMD velocities"), `*.xsc` (`Tabular`, "NAMD extended system"),
  `*.log` (`Log`); probe via `find_on_path(&["namd2", "namd3"])`
  pushes an `"academic"`-keyworded warning containing both
  `"academic"` and `"non-commercial"` substrings; custom NAMD-
  License surfaced as `tool_license = "NAMD-License"`; version
  range `2.14.0..4.0.0`; `bio.namd.simulate` ribbon capability),
  AmberTools sander (single-binary subprocess sister to Phase 18
  BWA; `case.toml` knobs `topology` (`.prmtop` / `.parm7`;
  required), `coordinates` (`.inpcrd` / `.rst7`; required),
  `config` (`.in` / `.mdin`; required), `output_basename`,
  `extra_args`; `prepare()` resolves all three input paths against
  the case directory when relative, validates each file exists on
  disk, composes `sander -O -i <config> -p <topology> -c
  <coordinates> -o <basename>.out -r <basename>.rst -x
  <basename>.nc [extras...]`; collects `<basename>*.out` (`Log`,
  "sander mdout"), `<basename>*.nc` (`Native`, "sander NetCDF
  trajectory"), `<basename>*.rst` (`Native`, "sander restart
  coordinates"), `<basename>*.mdinfo` (`Log`, "sander mdinfo");
  GPL-3.0; sister to Phase 5.5 cpptraj — installing AmberTools
  installs both; version range `22.0.0..26.0.0`;
  `bio.sander.simulate` ribbon capability), and HOOMD-blue
  (Python-script subprocess sister to Phase 17 OpenMM; `case.toml`
  knobs `script` (`.py` enforced), `python` (default `"python3"`),
  `output_basename`; `prepare()` enforces the `.py` extension,
  stages the script into the workdir under its original filename,
  writes a flat `valenx_params.json` containing `output_basename`,
  builds `<python> <staged_script>`; collects `<basename>*.gsd`
  (`Native`, "HOOMD trajectory (GSD)"), `<basename>*.h5` (`Native`,
  "HOOMD HDF5 output"), `*.log`; probe via Python on PATH with an
  `import hoomd` check (returns `ok = true` with a warning when
  import fails so non-standard installs aren't blocked);
  BSD-3-Clause; version range `3.0.0..6.0.0`;
  `bio.hoomd.simulate` ribbon capability). Each wired into
  `valenx-app::init_registry`.
- No new canonical types, no new format readers, no new CLIs. All
  three adapters consume user-supplied inputs (NAMD configuration
  files, sander topology + coordinates + simulation control files,
  HOOMD-blue Python scripts) and emit user-readable artifacts
  (NAMD `.dcd` / `.coor` / `.vel` / `.xsc`, sander `.out` / `.nc`
  / `.rst` / `.mdinfo`, HOOMD-blue `.gsd` / `.h5` / `.log`) that
  the unchanged `Results.artifacts` collection model surfaces
  directly.
- 3 new `valenx-init` templates ship: `namd` (`namd-simulate`),
  `sander` (`sander-simulate`), `hoomd` (`hoomd-simulate`).

The full per-phase shape lives in
[docs/src/phases/phase-5-6-md-engines.md](./docs/src/phases/phase-5-6-md-engines.md);
the implementation plan lives at
[docs/superpowers/plans/2026-05-02-md-engines.md](./docs/superpowers/plans/2026-05-02-md-engines.md).

**Deliverable:** Adapter inventory at **108 of 109** fully live
after this phase (3 new bio MD-engine adapters added on top of
the existing 105-adapter total). Phase 5.6 rounds out the all-
atom + GPU-native MD-engine surface alongside the Phase 5 GROMACS
/ LAMMPS pair and the Phase 17 OpenMM Python-native engine,
broadening the simulate-MD → analyze-trajectory loop to six MD
engines (Phase 5 GROMACS / LAMMPS + Phase 17 OpenMM + Phase 5.6
NAMD / sander / HOOMD-blue).

### Phase 5.7 — MDTraj · live
Single-adapter sister to the Phase 17 MDAnalysis adapter and the
Phase 5.5 PLUMED / ProDy / cpptraj analysis trio. Round out the
post-MD analysis surface with the second-most-used Python MD
trajectory analyzer — **MDTraj** (Pande / VanderSpoel / Beauchamp
lab, LGPL-2.1), which has wider format support than MDAnalysis
(`.xtc` / `.dcd` / `.h5` / `.nc` / `.trr` / `.binpos` / `.lh5` /
`.amber` / `.gromacs`) and deeper integration with the OpenMM
ecosystem (the Pande / Beauchamp lab is co-located with the
OpenMM developers — MDTraj's HDF5 trajectory format is OpenMM's
native streaming output). Single-adapter phases are a precedent
in Valenx — when an established tool fills a clearly-defined
corner of an existing surface without requiring new infrastructure,
the phase ships as a single adapter. Phase 5.7 sits numerically
adjacent to Phase 5.5 + Phase 5.6 and ships chronologically right
after Phase 5.6 — same chronological-vs-numerical convention used
for Phase 17.5 / 24 / 28 / 31 / 35 / 39 / 5.5 / 5.6.

- 1 new adapter crate lands under `crates/valenx-adapters/bio/`:
  MDTraj (Python-script subprocess sister to Phase 17 Biopython,
  Phase 5.5 ProDy, Phase 17 OpenMM; `case.toml` knobs `script`
  (`.py` enforced), `python` (default `"python3"`), `trajectory`
  (`.xtc` / `.dcd` / `.h5` / `.nc` / `.trr` / `.binpos` / `.lh5`
  MDTraj-supported trajectory; required), `topology` (`.pdb` /
  `.prmtop` / `.gro` / `.psf` topology MDTraj uses for atom +
  residue + chain metadata; required), `output_basename`;
  `prepare()` enforces the `.py` extension on the script, resolves
  all three input paths against the case directory when relative,
  stages script + trajectory + topology into the workdir under
  their original filenames so the script can resolve them via
  relative paths, then writes a flat `valenx_params.json`
  containing `output_basename`, the bare `trajectory` filename,
  and the bare `topology` filename, builds `<python>
  <staged_script>`; collects `<output_basename>*.csv` (`Tabular`,
  "MDTraj analysis table"), `<output_basename>*.npz` (`Native`,
  "MDTraj numpy archive"), `<output_basename>*.h5` (`Native`,
  "MDTraj HDF5 output"), `<output_basename>*.png` (`Native`,
  "MDTraj plot"), `*.log`; probe via Python on PATH with an
  `import mdtraj` check (returns `ok = true` with a warning when
  import fails so non-standard installs aren't blocked); LGPL-2.1;
  version range `1.9.0..2.0.0`; `bio.mdtraj.analyze` ribbon
  capability). Wired into `valenx-app::init_registry`.
- No new canonical types, no new format readers, no new CLIs.
  MDTraj consumes user-supplied inputs (MDTraj Python analysis
  scripts + trajectories + topologies) and emits user-readable
  artifacts (MDTraj `.csv` / `.npz` / `.h5` / `.png` outputs)
  that the unchanged `Results.artifacts` collection model surfaces
  directly.
- 1 new `valenx-init` template ships: `mdtraj` (`mdtraj-analyze`).

The full per-phase shape lives in
[docs/src/phases/phase-5-7-mdtraj.md](./docs/src/phases/phase-5-7-mdtraj.md);
the implementation plan lives at
[docs/superpowers/plans/2026-05-02-mdtraj.md](./docs/superpowers/plans/2026-05-02-mdtraj.md).

**Deliverable:** Adapter inventory at **109 of 110** fully live
after this phase (1 new MDTraj adapter added on top of the
existing 108-adapter total alongside the Phase 5.6 bio MD-engine
trio). Phase 5.7 rounds out the post-MD analysis surface alongside
the Phase 17 MDAnalysis and Phase 5.5 PLUMED / ProDy / cpptraj
beachheads — broadening the analyze-trajectory loop to five MD-
analysis tools (Phase 17 MDAnalysis + Phase 5.5 PLUMED / ProDy /
cpptraj + Phase 5.7 MDTraj).

### Phase 17 — Biology + biotech foundation · live
The first non-physics domain ships in Valenx — a beachhead that brings
the biology / biotech tool ecosystem under the same shell as the existing
physics-domain coverage.

- New `valenx-bio` crate hosts canonical types: `Sequence` (DNA / RNA /
  Protein with IUPAC alphabet validation), `Structure` (atom / residue /
  chain hierarchy round-tripping through PDB ATOM records), and
  `Trajectory` (per-frame atomic coordinates from MD output).
- Format readers: FASTA (with 60-char body wrap on output), PDB (ATOM /
  HETATM, v3 column layout), DCD (NAMD / CHARMM / VMD interchange).
  mmCIF reader stub deferred to Phase 17.5.
- 7 first-class adapters covering the most-used workflows: Biopython,
  RDKit, OpenMM (Python-native MD), ChimeraX (3D viz), oxDNA (CG DNA),
  MDAnalysis (trajectory analysis), ColabFold (protein structure
  prediction). Each wired into `valenx-app::init_registry`.
- 3 new headless CLIs round out the workflow loop: `valenx-fasta`,
  `valenx-pdb-info`, `valenx-blast` (alphabet auto-detection routing
  to blastp / blastn).
- 7 new `valenx-init` templates (one per bio adapter); the cross-binary
  roundtrip test sweeps all 20 templates clean.

The full per-phase shape lives in
[docs/src/phases/phase-17-biology.md](./docs/src/phases/phase-17-biology.md);
the implementation plan plus the future-phases scope-out is at
[docs/superpowers/plans/2026-04-30-biology-foundation.md](./docs/superpowers/plans/2026-04-30-biology-foundation.md).

**Deliverable:** Adapter inventory at 24 of 25 fully live (only `occt`
remains stub-only pending an `occt-sys` C++ FFI shim); biology workflows
land alongside the existing physics domains in the same case-toml /
prepare / run / collect shell.

### Phase 17.5 — Structure prediction expansion · live
The Phase 17 ColabFold adapter expands into the full set of open-source
sibling tools that take a FASTA query (or AF3-style JSON job spec) and
produce ranked PDB models with pLDDT scores in the B-factor column.

- 4 new first-class adapters: ESMFold (Meta protein language model —
  single-sequence prediction, no MSA), OpenFold (PyTorch AF2
  reimplementation with full preset family validated at the case-input
  layer), AlphaFold 2 (DeepMind reference `run_alphafold.py`;
  `monomer` / `monomer_ptm` / `multimer` presets), AlphaFold 3
  (DeepMind all-atom complex predictor — JSON job spec, not FASTA).
  Each wired into `valenx-app::init_registry`.
- AlphaFold 3's probe pushes a non-commercial-weights warning into
  `ProbeReport.warnings` because AF3's model weights are released
  under CC-BY-NC-4.0; the per-tool license surfaces in the registry
  UI even though the adapter's invocation mode stays
  `LicenseMode::Subprocess`.
- No new canonical types, no new format readers, no new CLIs. pLDDT
  rides the existing `Atom.b_factor` field that the Phase 17 PDB
  reader already lifts cleanly. AF3's JSON job spec is consumed as-is
  by the underlying tool.
- 4 new `valenx-init` templates (`esmfold` / `openfold` / `alphafold2`
  with `af2` alias / `alphafold3` with `af3` alias); the cross-binary
  roundtrip test sweeps all 30 templates clean.

The full per-phase shape lives in
[docs/src/phases/phase-17-5-structure-prediction.md](./docs/src/phases/phase-17-5-structure-prediction.md);
the implementation plan plus the future-phases scope-out is at
[docs/superpowers/plans/2026-04-30-structure-prediction-expansion.md](./docs/superpowers/plans/2026-04-30-structure-prediction-expansion.md).

**Deliverable:** Adapter inventory at 34 of 35 fully live after this
phase (4 new structure-prediction adapters added on top of the
Phase 17 + Phase 18 totals — Phase 17.5 ships chronologically after
Phase 18 even though it sits numerically between Phase 17 and
Phase 18); structure-prediction workflows land alongside the existing
biology adapters in the same case-toml / prepare / run / collect shell.

### Phase 17.7 — Structure tools expansion · live
Sister-adapter expansion of the existing Phase 17.5 structure-
prediction beachhead and the Phase 17 ColabFold adapter. Round out
the protein structure prediction + structure search surface with
three more foundational open-source tools — Baker lab's original
3-track structure-prediction network (RoseTTAFold, the canonical
pre-AlphaFold-3 sibling that established the 3-track SE(3)-
equivariant attention pattern; MIT), HelixonAI's single-sequence
structure predictor (OmegaFold, MSA-free like ESMFold but with a
larger pre-trained transformer backbone; Apache-2.0), and
Steinegger lab's 3D-structure search tool (FoldSeek, the protein-
3D analogue of the Phase 18.5 MMseqs2 sequence search; GPL-3.0).
RoseTTAFold + OmegaFold follow established Python-script subprocess
shapes (RoseTTAFold sister to Phase 17.5 ESMFold / OpenFold;
OmegaFold ships its own CLI binary with Python fallback). FoldSeek
follows the established Phase 18.5 MMseqs2 single-binary CLI
pattern. Phase 17.7 sits numerically after Phase 17.5 (Phase 17.6
was reserved for a deferred confidence-aware structure-ranking +
mmCIF-reader work that hasn't shipped yet) and ships chronologically
right after Phase 5.7 MDTraj — same chronological-vs-numerical
convention used for Phase 17.5 / 24 / 28 / 31 / 35 / 39 / 5.5 /
5.6 / 5.7.

- 3 new adapter crates land under `crates/valenx-adapters/bio/`:
  RoseTTAFold (Python-script subprocess sister to Phase 17.5
  ESMFold / Phase 17 Biopython; `case.toml` knobs `script` (`.py`
  enforced), `python` (default `"python3"`), `fasta` (input FASTA
  query sequence; required), `output_basename`; `prepare()`
  enforces the `.py` extension, resolves `script` and `fasta`
  against the case directory when relative, stages both into the
  workdir under their original filenames, writes a flat
  `valenx_params.json` containing `output_basename` and the bare
  `fasta` filename, builds `<python> <staged_script>`; collects
  `<output_basename>*.pdb` (`Native`, "RoseTTAFold predicted
  structure" — pLDDT-style per-residue confidence in the B-factor
  column, lifted by the existing Phase 17 PDB reader without any
  structure-prediction-specific code path), `<output_basename>*.npz`
  (`Native`, "RoseTTAFold confidence arrays"), `*.log`; probe via
  `find_on_path(&["python3", "python"])` — deliberately doesn't try
  `import rosettafold` (RoseTTAFold is not a pip package, it's a
  clone-from-GitHub install with heavy ML deps); pushes a probe
  warning whenever Python is detected: "RoseTTAFold model weights
  + dependencies not bundled — clone
  https://github.com/RosettaCommons/RoseTTAFold and follow the
  install README"; MIT; version range `1.0.0..3.0.0` (RoseTTAFold
  1.x is the original 2021 release; RoseTTAFold 2 / RoseTTAFold
  All-Atom is the late-2023 follow-up that adds nucleic-acid +
  small-molecule support); `bio.rosettafold.predict` ribbon
  capability), OmegaFold (single-binary CLI subprocess with
  Python fallback; `case.toml` knobs `fasta` (input FASTA query
  sequence; required), `output_basename` (workdir-relative output
  directory name OmegaFold writes the predicted PDBs under),
  `python` (interpreter name; default `"python3"`; used only as
  fallback when the OmegaFold CLI isn't on PATH), `model_dir`
  (`Option<PathBuf>` — optional pre-downloaded model checkpoint
  directory; OmegaFold defaults to `~/.cache/omegafold_ckpt` when
  omitted); `prepare()` builds `omegafold <fasta> <output_basename>
  [--model <model_dir>]` with the FASTA passed by absolute path
  (NOT staged into the workdir — OmegaFold reads it once then
  writes everything else into the output directory); collects
  walks one level deep into the `<output_basename>/` subdirectory
  for `*.pdb` (`Native`, "OmegaFold predicted structure" — per-
  residue confidence in the B-factor column) and `*.json` (`Log`,
  "OmegaFold metadata"), plus the workdir-top-level `*.log`; probe
  via `find_on_path(&["omegafold", "python3", "python"])` —
  surfaces a warning if `omegafold` itself isn't on PATH but
  Python is ("OmegaFold CLI not found on PATH; install via pip
  install git+https://github.com/HeliXonProtein/OmegaFold.git");
  Apache-2.0; version range `1.0.0..2.0.0`;
  `bio.omegafold.predict` ribbon capability), and FoldSeek
  (single-binary subprocess sister to Phase 18.5 MMseqs2 / Phase
  18 BWA; `case.toml` knobs `query` (`.pdb` / `.cif` query
  structure; required), `database` (FoldSeek database path prefix
  — the user supplies the path stem and FoldSeek resolves the
  `<prefix>_*` sidecar files itself; required), `output_basename`,
  `threads` (`u32`, default 1), `extra_args`; `prepare()` resolves
  `query` against the case directory when relative, validates the
  `database` parent directory exists on disk (the database files
  themselves use the prefix convention so we cannot validate them
  by name — same shape as Phase 18.7 BLAST+'s `database`
  validation), composes `foldseek easy-search <query> <database>
  <basename>.m8 tmp_<basename> --threads <N> [extras...]` (the
  `tmp_<basename>` is a per-run temp directory FoldSeek requires);
  collects `<output_basename>.m8` (`Tabular`, "FoldSeek search
  results" — the canonical BLAST-style M8 hit table format every
  downstream FoldSeek pipeline reads) and `*.log`; the temp
  directory is not surfaced in artifacts — it's intermediate;
  probe via `find_on_path(&["foldseek"])`; GPL-3.0; version range
  `8.0.0..10.0.0`; `bio.foldseek.search` ribbon capability). Each
  wired into `valenx-app::init_registry`.
- No new canonical types, no new format readers, no new CLIs. All
  three adapters consume user-supplied inputs (RoseTTAFold +
  OmegaFold FASTA queries, FoldSeek PDB / CIF structures +
  database prefixes) and emit user-readable artifacts (RoseTTAFold
  + OmegaFold `.pdb` predicted structures with per-residue
  confidence in the B-factor column flowing through the existing
  `valenx_bio::format::pdb` reader, RoseTTAFold `.npz` confidence
  arrays, OmegaFold `.json` metadata sidecars, FoldSeek `.m8`
  BLAST-style hit tables) that the unchanged `Results.artifacts`
  collection model surfaces directly.
- 3 new `valenx-init` templates ship: `rosettafold`
  (`rosettafold-predict`), `omegafold` (`omegafold-predict`),
  `foldseek` (`foldseek-search`).

The full per-phase shape lives in
[docs/src/phases/phase-17-7-structure-tools.md](./docs/src/phases/phase-17-7-structure-tools.md);
the implementation plan lives at
[docs/superpowers/plans/2026-05-02-structure-tools.md](./docs/superpowers/plans/2026-05-02-structure-tools.md).

**Deliverable:** Adapter inventory at **112 of 113** fully live
after this phase (3 new structure-tools-expansion adapters added
on top of the existing 109-adapter total alongside the Phase 5.6
bio MD-engine trio and Phase 5.7 MDTraj single-adapter); Phase
17.7 rounds out the protein structure prediction + structure
search surface alongside the Phase 17.5 ESMFold / OpenFold /
AlphaFold 2 / AlphaFold 3 and Phase 17 ColabFold beachheads,
broadening the predict-structure → search-structure loop to seven
structure predictors (Phase 17 ColabFold + Phase 17.5 ESMFold /
OpenFold / AlphaFold 2 / AlphaFold 3 + Phase 17.7 RoseTTAFold /
OmegaFold) and one structure-search tool (Phase 17.7 FoldSeek).

### Phase 18 — Sequence alignment toolkit · live
The second non-physics phase ships in Valenx — extending the biology
beachhead with the most-used alignment + read-mapping toolset.

- `valenx-bio` extends with two new canonical types: `FastqRecord`
  (sequence + per-base quality, length-validated) and `Alignment`
  (multiple-sequence alignment as a list of named gapped sequences
  with shared length).
- New format readers: FASTQ (4-line) and a minimal SAM-text reader
  (header + records sufficient for summary inspection). BAM (binary
  BGZF) deferred to Phase 18.5.
- 6 first-class adapters covering the user-visible workflows: BWA
  (short-read alignment), minimap2 (long-read + spliced + asm-vs-asm),
  MAFFT and MUSCLE (multiple-sequence alignment), HMMER (profile-HMM
  search), samtools (SAM/BAM utilities). Each wired into
  `valenx-app::init_registry`.
- 2 new headless CLIs round out the workflow loop: `valenx-fastq`
  (FASTQ inspect / validate) and `valenx-sam-info` (SAM alignment
  summary).
- 6 new `valenx-init` templates (one per alignment adapter); the
  cross-binary roundtrip test sweeps all 26 templates clean.

The full per-phase shape lives in
[docs/src/phases/phase-18-alignment.md](./docs/src/phases/phase-18-alignment.md);
the implementation plan plus the future-phases scope-out is at
[docs/superpowers/plans/2026-04-30-sequence-alignment-toolkit.md](./docs/superpowers/plans/2026-04-30-sequence-alignment-toolkit.md).

**Deliverable:** Adapter inventory at 30 of 31 fully live (only `occt`
remains stub-only); sequence-alignment workflows land alongside the
Phase 17 biology adapters in the same case-toml / prepare / run /
collect shell.

### Phase 18.5 — Aligners expansion · live
Sister-adapter expansion of Phase 18. Add three more aligners covering
distinct user-facing use cases: Bowtie2 (Langmead & Salzberg's gapped
FM-index short-read aligner — alternative to BWA for RNA-seq /
ChIP-seq / bisulfite pipelines), MMseqs2 (Söding lab's many-vs-many
protein search + clustering toolkit — fast alternative to BLAST and
the prefilter behind ColabFold's MSA generation), and DIAMOND
(Buchfink, Reuter & Drost's ultra-fast BLAST-protocol-compatible
protein aligner — two to three orders of magnitude faster than
BLASTP / BLASTX for whole-metagenome and UniRef-scale searches). All
three follow the established Phase 18 BWA pattern — single-binary CLI
subprocess, file in / file out. Bowtie2 mirrors BWA's two-stage
`index → align` shape; MMseqs2 and DIAMOND dispatch per-action via
the bcftools-style `build_command(...) -> Result<Vec<OsString>,
AdapterError>` helper. Phase 18.5 sits numerically after Phase 18 but
ships chronologically after Phase 27.5 — the same chronological-vs-
numerical convention used for Phase 17.5.

- 3 new adapter crates land under `crates/valenx-adapters/bio/`:
  Bowtie2 (single-binary subprocess; two-stage `bowtie2-build →
  bowtie2` pipeline; `case.toml` knobs `reference` (FASTA, required),
  `reads` (1 entry single-end, 2 entries paired-end), `threads`
  default 1, `skip_index` default false, `preset` default
  `"sensitive"` whitelist `["very-fast", "fast", "sensitive",
  "very-sensitive"]`, `extra_args`; collects `out.sam` typed
  `Tabular` "Bowtie2 aligned reads" and any `.log` files;
  `bio.bowtie2.align` ribbon capability; GPL-3.0 licensed),
  MMseqs2 (single-binary subprocess via `mmseqs <action>`; per-action
  dispatch on `action ∈ ["easy-search", "easy-cluster",
  "easy-linsearch"]`; `case.toml` knobs `query`, `target` (required
  for search modes, ignored for cluster), `output`, `sensitivity`
  default 7.5 range `1.0..=7.5` finite-checked, `threads` default 1,
  `extra_args`; binary is `mmseqs` no `2` suffix; collects `output`
  typed `Tabular` with per-action label like "MMseqs2 easy-search
  hits"; `bio.mmseqs2.search` ribbon capability; MIT licensed),
  DIAMOND (single-binary subprocess via `diamond <action>`; per-
  action dispatch on `action ∈ ["blastp", "blastx", "makedb"]`;
  `case.toml` knobs `query`, `database`, `output`, `sensitivity`
  whitelist `["default", "fast", "sensitive", "more-sensitive",
  "very-sensitive", "ultra-sensitive"]`, `threads`, `extra_args`;
  `--default` flag is omitted when `sensitivity = "default"` because
  DIAMOND's out-of-the-box default has no flag; in `makedb` mode the
  field roles flip — `query` is the input FASTA and `database` is
  the output `.dmnd` basename; collects `output` typed `Tabular`
  "DIAMOND <action> hits" for blastp/blastx, `<database>.dmnd` typed
  `Native` "DIAMOND .dmnd database" for makedb;
  `bio.diamond.search` ribbon capability; GPL-3.0 licensed). Each
  wired into `valenx-app::init_registry`.
- No new canonical types, no new format readers, no new CLIs. All
  three adapters consume the existing Phase 17 FASTA + Phase 18
  FASTQ inputs and emit SAM (Bowtie2) or tabular hit-table
  (MMseqs2 / DIAMOND BLAST format-8) outputs that the unchanged
  `Results.artifacts` collection model surfaces directly. Bowtie2
  SAM outputs are inspectable through the existing
  `valenx-sam-info` CLI; MMseqs2 + DIAMOND tabular hits surface
  via the existing `Tabular` artifact kind.
- 3 new `valenx-init` templates ship: `bowtie2` with alias `bt2`
  (`bowtie2-align`), `mmseqs2` with alias `mmseqs`
  (`mmseqs2-search`), and `diamond` with alias `dmnd`
  (`diamond-search`). Cross-binary roundtrip test sweeps all 53
  templates clean.

The full per-phase shape lives in
[docs/src/phases/phase-18-5-aligners.md](./docs/src/phases/phase-18-5-aligners.md);
the implementation plan lives at
[docs/superpowers/plans/2026-04-30-aligners-expansion.md](./docs/superpowers/plans/2026-04-30-aligners-expansion.md).

**Deliverable:** Adapter inventory at 57 of 58 fully live after this
phase (3 new aligners-expansion adapters added on top of the Phase 17
+ Phase 17.5 + Phase 18 + Phase 19 + Phase 19.5 + Phase 22 + Phase 23
+ Phase 24 + Phase 27 + Phase 27.5 + Phase 34 totals); aligners-
expansion workflows land alongside the existing biology adapters in
the same case-toml / prepare / run / collect shell, broadening the
search → align → predict → validate loop to nine alignment / search
tools (BWA, Bowtie2, minimap2, MAFFT, MUSCLE, HMMER, samtools,
MMseqs2, DIAMOND) feeding into the five Phase 17 + 17.5 prediction
tools.

### Phase 18.6 — RNA-seq alignment · live
Sister-adapter expansion of Phase 18 / 18.5 closing the splice-aware
RNA-seq alignment gap that Phase 18 and Phase 18.5 explicitly deferred.
Add the two de-facto RNA-seq aligners to Valenx: HISAT2 (Daehwan Kim's
graph-based splice-aware aligner — successor to TopHat) and STAR
(Alex Dobin's most-used spliced aligner, the reference RNA-seq mapper
backing GTEx / TCGA / ENCODE pipelines and the only Phase 18.x aligner
that doubles as a chromatin-conformation tool). Both are spliced
extensions of Phase 18's BWA / Phase 18.5's Bowtie2 — they handle
reads that span exon-exon junctions, where the linear short-read
aligners would soft-clip or misalign. Both adapters mirror the
established Phase 18 BWA two-stage shape — single-binary CLI
subprocess, file in / file out, `index → align` pipeline. STAR has a
heavier index step (genomic + splice-junction database) but the same
overall shape. No new infrastructure. Phase 18.6 sits numerically
after Phase 18.5 and ships chronologically right after it.

- 2 new adapter crates land under `crates/valenx-adapters/bio/`:
  HISAT2 (single-binary subprocess; two-stage `hisat2-build →
  hisat2` pipeline; `case.toml` knobs `reference` (FASTA, required),
  `reads` (1 entry single-end, 2 entries paired-end), `threads`
  default 1, `skip_index` default false, `strandness` default
  `"unstranded"` whitelist `["unstranded", "F", "R", "FR", "RF"]`
  (F/R variants match Illumina TruSeq stranded library prep
  conventions), `extra_args`; the `--rna-strandness` flag is
  omitted when `strandness = "unstranded"` because HISAT2 treats
  unstranded data as the default; collects `out.sam` typed
  `Tabular` "HISAT2 aligned reads"; `bio.hisat2.align` ribbon
  capability; GPL-3.0 licensed), STAR (single-binary subprocess
  with capitalized binary name — `find_on_path(&["STAR"])`, not
  `star`; two-stage `--runMode genomeGenerate → --runMode
  alignReads` pipeline; STAR's index step is heavier than BWA /
  Bowtie2 / HISAT2 — it builds a suffix-array-indexed genome under
  `genome_dir/` and optionally a splice-junction database from a
  GTF — but the adapter shape is the same; `case.toml` knobs
  `genome_dir` (the pre-built STAR index directory, or where the
  adapter writes one if `skip_index = false`; required),
  `reference` (FASTA; required only when generating the index),
  `reads` (1 or 2 entries), `threads` default 1, `skip_index`
  default false, `output_type` default `"BAM_SortedByCoordinate"`
  whitelist `["BAM_Unsorted", "BAM_SortedByCoordinate", "SAM"]`
  (underscore canonical names map to STAR's two-arg `--outSAMtype`
  form, e.g. `"BAM_SortedByCoordinate"` →
  `--outSAMtype BAM SortedByCoordinate`), `sjdb_gtf` optional GTF
  for splice-junction-database-aware indexing, `extra_args`;
  collects `star_Aligned.out.{bam,sam}` typed `Tabular` for SAM /
  `Native` for BAM "STAR aligned reads" and `star_Log.final.out`
  typed `Log` "STAR alignment summary"; `bio.star.align` ribbon
  capability; MIT licensed). Each wired into
  `valenx-app::init_registry`.
- No new canonical types, no new format readers, no new CLIs. Both
  adapters consume the existing Phase 17 FASTA + Phase 18 FASTQ
  inputs and emit SAM (HISAT2 + STAR `SAM` mode) or BAM (STAR `BAM_*`
  modes) outputs that the unchanged `Results.artifacts` collection
  model surfaces directly. HISAT2 SAM outputs are inspectable
  through the existing `valenx-sam-info` CLI; STAR BAM outputs need
  the existing samtools adapter (`samtools view`) to convert to SAM
  before `valenx-sam-info` can read them; STAR's
  `star_Log.final.out` is plain text and surfaces directly through
  the `Log` artifact kind.
- 2 new `valenx-init` templates ship: `hisat2` with alias `hisat`
  (`hisat2-align`) and `star` (`star-align`). Cross-binary roundtrip
  test sweeps all 55 templates clean.

The full per-phase shape lives in
[docs/src/phases/phase-18-6-rna-seq.md](./docs/src/phases/phase-18-6-rna-seq.md);
the implementation plan lives at
[docs/superpowers/plans/2026-04-30-rna-seq-alignment.md](./docs/superpowers/plans/2026-04-30-rna-seq-alignment.md).

**Deliverable:** Adapter inventory at 59 of 60 fully live after this
phase (2 new RNA-seq alignment adapters added on top of the Phase 17
+ Phase 17.5 + Phase 18 + Phase 18.5 + Phase 19 + Phase 19.5 + Phase 22
+ Phase 23 + Phase 24 + Phase 27 + Phase 27.5 + Phase 34 totals);
RNA-seq alignment workflows land alongside the existing biology
adapters in the same case-toml / prepare / run / collect shell,
broadening the search → align → predict → validate loop to eleven
alignment / search tools (BWA, Bowtie2, HISAT2, STAR, minimap2, MAFFT,
MUSCLE, HMMER, samtools, MMseqs2, DIAMOND) feeding into the five
Phase 17 + 17.5 prediction tools.

### Phase 18.7 — Alignment toolkit expansion · live ✅
Sister-adapter expansion of Phase 18 / 18.5 / 18.6 rounding out the
foundational sequence-alignment surface with three more established
open-source tools the existing BWA / minimap2 / MAFFT / MUSCLE adapters
explicitly left out: BLAST+ (NCBI's seminal sequence-database search
tool — five user-facing search programs `blastn` / `blastp` / `blastx`
/ `tblastn` / `tblastx` covering every nucleotide / protein search
direction; Public Domain — US government work), Clustal Omega
(Sievers / Higgins' modern HMM-driven progressive multiple-sequence
aligner — modern successor to ClustalW that scales to thousands of
sequences; GPL-2.0), and T-Coffee (Notredame / Higgins' library-based
consistency-weighted multiple-sequence aligner — combines pairwise
alignments from many sources into a single consistency-weighted MSA,
the canonical choice for difficult distantly-related sequences;
GPL-2.0). All three follow the established Phase 18 BWA single-binary
CLI pattern: file in, alignment table / search hits out. No new
infrastructure. Phase 18.7 sits numerically after Phase 18.6 and ships
chronologically right after it — same convention used for Phase 17.5 /
24 / 28 / 31 / 35 / 39.

- 3 new adapter crates land under `crates/valenx-adapters/bio/`:
  BLAST+ (single-binary subprocess sister to Phase 18 BWA; per-program
  CLI `<program> -query <query> -db <database> -out blast_results.txt
  -evalue <evalue> -outfmt <outfmt> -num_threads <threads> [extras...]`;
  `case.toml` knobs `program` (one of the five — required), `query`
  (FASTA query file; required), `database` (BLAST database path prefix
  — the user supplies the path stem and BLAST resolves the
  `<prefix>.nhr` / `<prefix>.phr` etc. sidecars itself; required),
  `evalue` (`f64`, default 10.0 — BLAST's own default), `outfmt`
  (`u8`, default 0 — pairwise text; 6 = tabular), `threads` (`usize`,
  default 1), `extra_args`; `prepare()` resolves `query` against the
  case directory when relative, validates the database parent directory
  exists on disk (the database files themselves use the prefix
  convention so we can't validate them by name), looks up the
  per-program binary via `find_on_path(&[&input.program])`, and pins
  the output filename to `blast_results.txt`; collects
  `blast_results.txt` (`Tabular`, "BLAST search results") + `*.log`
  (`Log`); probe via `find_on_path(&["blastn", "blastp"])` — at least
  one BLAST+ binary on PATH counts as installed; `bio.blast.search`
  ribbon capability; Public Domain licensed), Clustal Omega
  (single-binary subprocess sister to Phase 18 MAFFT; CLI `clustalo
  -i <input> -o <basename>.<ext> --outfmt=<outfmt> --threads=<N>
  [extras...]`; `case.toml` knobs `input` (FASTA multi-sequence input;
  required), `output_basename` (filename stem; required), `outfmt`
  (default `"clustal"`; whitelist `clustal` / `fasta` / `phylip` /
  `vienna` / `nexus`), `threads` (`usize`, default 1), `extra_args`;
  `prepare()` derives `<ext>` from `outfmt` (`clustal` → `.aln`,
  `fasta` → `.fasta`, `phylip` → `.phy`, `vienna` → `.vie`, `nexus`
  → `.nex`, default `.aln`); collects `<output_basename>*` (`Tabular`,
  "Clustal Omega alignment") + `*.log` (`Log`); probe via
  `find_on_path(&["clustalo"])`; `bio.clustalo.align` ribbon
  capability; GPL-2.0 licensed), T-Coffee (single-binary subprocess
  sister to Clustal Omega; CLI `t_coffee <input> -output=<outfmt>
  -outfile=<basename>.aln [-mode=<mode>] [extras...]` (note T-Coffee's
  `=`-style flag form); `case.toml` knobs `input`, `output_basename`,
  `outfmt` (default `"clustalw"`; T-Coffee's own naming follows
  ClustalW conventions: `clustalw` / `fasta_aln` / `phylip` / `msf`),
  `mode` (`Option<String>` — omit for default progressive mode; set
  to `expresso` / `psicoffee` / `mcoffee` / etc. to opt into a
  specialised back-end), `extra_args`; output always pinned to `.aln`;
  collects `<output_basename>*` (`Tabular`, "T-Coffee alignment") +
  `*.dnd` (`Native`, "T-Coffee guide tree" — Newick guide tree
  consumed by downstream phylogenetics) + `*.log` (`Log`); probe via
  `find_on_path(&["t_coffee"])` — note the underscore: T-Coffee
  installs as `t_coffee`, not `t-coffee`; `bio.tcoffee.align` ribbon
  capability; GPL-2.0 licensed). Each wired into
  `valenx-app::init_registry`.
- No new canonical types, no new format readers, no new CLIs. All
  three adapters consume the existing Phase 17 FASTA inputs and emit
  standard tabular search-hit / MSA outputs that the unchanged
  `Results.artifacts` collection model surfaces directly. The existing
  Phase 17 `valenx-blast` CLI continues to wrap the auto-routing
  `blastp` / `blastn` shorthand for users who don't need the full
  adapter-driven pipeline.
- 3 new `valenx-init` templates ship: `blast` (`blast-search`),
  `clustalo` (`clustalo-align`), and `tcoffee` (`tcoffee-align`).
  Cross-binary roundtrip test sweeps all 101 templates clean alongside
  the Phase 19.6 single-cell-expansion pair.

The full per-phase shape lives in
[docs/src/phases/phase-18-7-blast-alignment.md](./docs/src/phases/phase-18-7-blast-alignment.md);
the implementation plan lives at
[docs/superpowers/plans/2026-05-02-blast-alignment.md](./docs/superpowers/plans/2026-05-02-blast-alignment.md).

**Deliverable:** Adapter inventory at 103 of 104 fully live after this
phase (3 new alignment-toolkit-expansion adapters added on top of the
Phase 17 + Phase 17.5 + Phase 18 + Phase 18.5 + Phase 18.6 + Phase 19
+ Phase 19.5 + Phase 20 + Phase 22 + Phase 23 + Phase 24 + Phase 25 +
Phase 27 + Phase 27.5 + Phase 27.6 + Phase 28 + Phase 29 + Phase 30 +
Phase 30.5 + Phase 31 + Phase 32 + Phase 33 + Phase 34 + Phase 35 +
Phase 36 + Phase 38 + Phase 39 + Phase 5.5 totals); alignment-toolkit
workflows land alongside the existing biology adapters in the same
case-toml / prepare / run / collect shell, broadening the search →
align → predict → validate loop to fourteen alignment / search tools
(BWA, Bowtie2, HISAT2, STAR, minimap2, MAFFT, MUSCLE, HMMER, samtools,
MMseqs2, DIAMOND, BLAST+, Clustal Omega, T-Coffee) feeding into the
five Phase 17 + 17.5 prediction tools.

### Phase 19 — Variant calling toolkit · live
The next link of the genomics workflow after Phase 18 read mapping:
variant calling from aligned reads. Adds the conventional + ML-driven
variant-calling stack on top of the Phase 18 alignment beachhead.

- `valenx-bio` extends with two new canonical types: `Vcf` (file-level:
  `##` header lines, sample IDs from the `#CHROM` column header, list
  of records) and `VcfRecord` (single variant row — chrom / pos /
  optional id / ref / comma-split alt / optional Phred qual / `;`-split
  filter / raw info / optional format + per-sample columns). `is_pass()`
  recognises both `["PASS"]` and the `"."` (= unfiltered) convention.
- New format reader: minimal VCF text reader (`valenx_bio::format::vcf`).
  Plain-text only; BCF (binary) and bgzf-compressed VCF deferred to
  Phase 19.5 — convert with `bcftools view` first.
- 3 first-class adapters covering the user-visible workflows: bcftools
  (VCF/BCF multitool — `view` / `call` / `filter` / `concat` per-action
  dispatch), GATK HaplotypeCaller (Broad Institute reference variant
  caller — Java heap validated against the conventional `8g` / `16g`
  suffix; optional intervals (BED) restriction), DeepVariant (Google
  ML-based caller — typed `model_type` ∈ `{WGS, WES, PACBIO,
  ONT_R104, HYBRID_PACBIO_ILLUMINA}` and `num_shards` knob; probe hint
  mentions both the direct binary and the Docker / Singularity wrapper
  paths, the adapter does not manage container runtimes). Each wired
  into `valenx-app::init_registry`.
- 1 new headless CLI rounds out the workflow loop: `valenx-vcf-info`
  (header-line count, sample count, total records, PASS / FAIL split,
  no-ALT count from a VCF file or stdin via `-`; text + JSON output).
- 3 new `valenx-init` templates ship: `bcftools` (`bcftools-call`),
  `gatk` (alias `hc`, `gatk-haplotype`), and `deepvariant` (alias
  `dv`, `deepvariant-call`). Cross-binary roundtrip test sweeps all
  33 templates clean.

What's not in this phase: BCF (binary VCF) reading, joint genotyping /
GVCF workflows, and the remaining variant callers that share
bcftools / GATK's adapter shape (Strelka2, FreeBayes, DELLY, Manta,
Pindel, vcftools). These land in Phase 19.5; variant annotation
(SnpEff / VEP) follows in Phase 43.

The full per-phase shape lives in
[docs/src/phases/phase-19-variant-calling.md](./docs/src/phases/phase-19-variant-calling.md);
the implementation plan plus the future-phases scope-out is at
[docs/superpowers/plans/2026-04-30-variant-calling-toolkit.md](./docs/superpowers/plans/2026-04-30-variant-calling-toolkit.md).

**Deliverable:** Adapter inventory at 37 of 38 fully live after this
phase (3 new variant-calling adapters added on top of the Phase 17 +
Phase 17.5 + Phase 18 totals); variant-calling workflows land
alongside the existing biology adapters in the same case-toml /
prepare / run / collect shell.

### Phase 19.5 — Single-cell genomics · live
Open the single-cell genomics domain in Valenx with the two most-used
Python tools: Scanpy (the de-facto single-cell analysis library —
clustering, dimensionality reduction, marker discovery) and scVI
(probabilistic deep-learning models for single-cell data via the
`scvi-tools` package). Both adapters follow the established Phase 17
Biopython / RDKit pattern — Python-script subprocess where the user's
script imports `scanpy` or `scvi` and reads `valenx_params.json`
(auto-written by the adapter) for config knobs. Phase 19.5 sits
numerically before Phase 22 but ships chronologically after Phase 22 —
the same chronological-vs-numerical convention used for Phase 17.5.

- 2 new adapter crates land under `crates/valenx-adapters/bio/`:
  Scanpy (de-facto Python single-cell analysis library, BSD-3-Clause;
  `valenx_params.json` knobs `input_h5ad` / `output_h5ad` /
  `n_top_genes` (default 2000) / `n_pcs` (default 50) / `n_neighbors`
  (default 15) / `resolution` (default 1.0); collects `.h5ad`
  ("Scanpy AnnData output"), `.png` / `.pdf` ("Scanpy plot"),
  `.csv` / `.tsv` (`Tabular`, "Scanpy table"); `bio.scanpy.analyse`
  ribbon capability) and scVI (probabilistic deep-learning models
  for single-cell data via `scvi-tools`, BSD-3-Clause;
  `valenx_params.json` knobs `input_h5ad` / `output_h5ad` / typed
  `model` ∈ `{scvi, scanvi, totalvi, linear-scvi}` (default `scvi`)
  / `n_latent` (default 10) / `n_layers` (default 2) / `max_epochs`
  (default 400) / optional `batch_key`; collects the same set as
  Scanpy under "scVI" labels; `bio.scvi.train` ribbon capability).
  Each wired into `valenx-app::init_registry`.
- No new canonical types, no new format readers, no new CLIs.
  AnnData reader as a canonical type (`.h5ad` is HDF5-backed and
  needs the `hdf5` crate, a non-trivial C-library dep) is deferred
  to Phase 19.6 along with the Seurat R-runtime work.
- 2 new `valenx-init` templates ship: `scanpy` (`scanpy-analyse`)
  and `scvi` with alias `scvi-tools` (`scvi-train`). Cross-binary
  roundtrip test sweeps all 47 templates clean.

The full per-phase shape lives in
[docs/src/phases/phase-19-5-single-cell.md](./docs/src/phases/phase-19-5-single-cell.md);
the implementation plan lives at
[docs/superpowers/plans/2026-04-30-single-cell-genomics.md](./docs/superpowers/plans/2026-04-30-single-cell-genomics.md).

**Deliverable:** Adapter inventory at 51 of 52 fully live after this
phase (2 new single-cell genomics adapters added on top of the
Phase 17 + Phase 17.5 + Phase 18 + Phase 19 + Phase 22 + Phase 23 +
Phase 24 + Phase 27 + Phase 34 totals — Phase 19.5 ships
chronologically after Phase 22 even though it sits numerically before
Phase 22); single-cell workflows land alongside the existing biology
adapters in the same case-toml / prepare / run / collect shell.

### Phase 19.6 — Single-cell expansion · live ✅
Sister-adapter expansion of the Phase 19.5 single-cell Scanpy + scVI
beachhead. Round out the single-cell genomics surface with the two
most-requested tools that Phase 19.5 explicitly deferred — the
dominant R-based single-cell analysis toolkit (Seurat, the Satija
lab's reference single-cell library that drives roughly half of
single-cell papers worldwide) and the canonical Python data-container
library (AnnData, the scverse foundation library that scanpy / scvi /
scirpy / squidpy / muon all read and write). Phase 19.6 explicitly
extends Phase 19.5 with R + AnnData: Seurat introduces the **Rscript
subprocess pattern** to Valenx (the R analogue of the Python-script
pattern that Phase 17 Biopython / RDKit / OpenMM, Phase 19.5 Scanpy /
scVI, and Phase 33 pySBOL established), and AnnData reuses the
existing Python-script pattern for the canonical container that ties
the Phase 19.5 scanpy / scvi adapters to every downstream scverse
tool. Phase 19.6 sits numerically after Phase 19.5 and ships
chronologically right after it — same convention used for Phase 17.5
/ 18.5 / 18.6 / 18.7 / 24 / 28 / 31 / 35.

- 2 new adapter crates land under `crates/valenx-adapters/bio/`:
  Seurat (introduces the **Rscript subprocess pattern** to Valenx;
  user supplies an `.R` script referenced from `[bio.seurat].script`
  in `case.toml` that loads `library(Seurat)` and reads
  `valenx_params.json` for the parsed knobs via `jsonlite::fromJSON`;
  `case.toml` knobs `script` (path to user-supplied R script;
  required, must end in `.R`), `rscript` (binary name; default
  `"Rscript"`), `input_data` (`Option<PathBuf>` — optional input
  matrix that the script loads; supports `.h5` / `.mtx` / `.rds` so
  users can drop in 10x HDF5, sparse Matrix Market, or pre-saved
  Seurat object formats), `output_basename`; `prepare()` enforces
  the `.R` extension, stages script + optional input_data into the
  workdir under their original filenames, then writes a flat
  `valenx_params.json` containing `output_basename` and `input_data`
  (staged filename when set; the key is **omitted entirely** when
  `input_data` is `None` rather than emitted as `null`, matching the
  hand-rolled JSON convention the rest of the bio adapters use);
  builds `native_command = [rscript, script]`; collects
  `<output_basename>*.rds` (`Native`, "Seurat object (RDS)" — the
  canonical R-serialised Seurat object format consumed by every
  downstream Seurat / signac / Azimuth pipeline),
  `<output_basename>*.csv` (`Tabular`, "Seurat output table"),
  `<output_basename>*.png` (`Native`, "Seurat plot"), and `*.log`
  (`Log`); probe via `find_on_path(&["Rscript"])` — surfaces a
  warning when Rscript is missing rather than failing (same shape as
  Phase 17 Biopython probe); the probe deliberately does **not**
  attempt to confirm Seurat itself is installed because that would
  require running R (an expensive multi-second startup at probe time
  that conflicts with the rest of the registry's snappy PATH-lookup
  probes); `bio.seurat.analyze` ribbon capability; MIT licensed) and
  AnnData (Python-script subprocess shape sister to Phase 19.5
  Scanpy / scVI; user supplies a Python script referenced from
  `[bio.anndata].script` in `case.toml` that imports `anndata` and
  reads `valenx_params.json` for the parsed knobs; `case.toml` knobs
  `script` (path to user-supplied Python script; required, must end
  in `.py`), `python` (interpreter name; default `"python3"`),
  `input_h5ad` (`Option<PathBuf>` — optional input single-cell file
  the script loads; supports `.h5ad` (the canonical AnnData format)
  and `.h5` (10x HDF5)), `output_basename`; `prepare()` enforces the
  `.py` extension, stages script + optional input_h5ad, writes
  `valenx_params.json` with the same hand-rolled shape as Seurat
  (`output_basename` plus `input_h5ad` (staged filename when set,
  key omitted when `None`)); builds `native_command = [python,
  script]`; collects `<output_basename>*.h5ad` (`Native`, "AnnData
  h5ad file" — the canonical AnnData HDF5-backed format every
  scverse tool reads and writes), `<output_basename>*.csv`
  (`Tabular`, "AnnData output table"), `<output_basename>*.png`
  (`Native`, "AnnData plot"), and `*.log` (`Log`); probe via
  `find_on_path(&["python3", "python"])` then `<python> -c "import
  anndata"` — on import failure surface as a `ProbeReport.warnings`
  entry (not error) so non-standard installs aren't blocked, same
  shape as the Scanpy / scVI / Biopython probes;
  `bio.anndata.process` ribbon capability; BSD-3-Clause licensed).
  Each wired into `valenx-app::init_registry`.
- No new canonical types, no new format readers, no new CLIs. The
  first-class `.h5ad` canonical type backed by the `hdf5` crate (a
  non-trivial C-library dep that Phase 19.5 already flagged as out
  of scope) defers to a future phase along with Seurat-object
  inspection beyond the existing artifact-collection model. Phase
  19.6 introduces the **Rscript subprocess pattern** to the
  workspace for the first time alongside the Seurat adapter — the
  runtime infrastructure (R-script staging, `valenx_params.json`
  shape compatible with `jsonlite::fromJSON`, `.R` extension
  validation, `Rscript` binary probe) ships in the Seurat adapter
  and is reusable from any future R-based bioinformatics adapter
  without further plumbing.
- 2 new `valenx-init` templates ship: `seurat` (`seurat-analyze`)
  and `anndata` (`anndata-process`). Cross-binary roundtrip test
  sweeps all 101 templates clean alongside the Phase 18.7
  alignment-toolkit-expansion trio.

The full per-phase shape lives in
[docs/src/phases/phase-19-6-single-cell-expansion.md](./docs/src/phases/phase-19-6-single-cell-expansion.md);
the implementation plan lives at
[docs/superpowers/plans/2026-05-02-single-cell-expansion.md](./docs/superpowers/plans/2026-05-02-single-cell-expansion.md).

**Deliverable:** Adapter inventory at 105 of 106 fully live after this
phase (2 new single-cell-expansion adapters added on top of the Phase
17 + Phase 17.5 + Phase 18 + Phase 18.5 + Phase 18.6 + Phase 18.7 +
Phase 19 + Phase 19.5 + Phase 20 + Phase 22 + Phase 23 + Phase 24 +
Phase 25 + Phase 27 + Phase 27.5 + Phase 27.6 + Phase 28 + Phase 29 +
Phase 30 + Phase 30.5 + Phase 31 + Phase 32 + Phase 33 + Phase 34 +
Phase 35 + Phase 36 + Phase 38 + Phase 39 + Phase 5.5 totals);
single-cell-expansion workflows land alongside the existing biology
adapters in the same case-toml / prepare / run / collect shell,
broadening the load-data → cluster → integrate → annotate → visualise
single-cell loop to four single-cell tools across both language
ecosystems (the Phase 19.5 Scanpy / scVI Python pair plus the Phase
19.6 Seurat R + AnnData Python container pair). The Rscript subprocess
pattern that Seurat introduces also opens the door to every other
R-based bioinformatics tool for future phases.

### Phase 20 — Transcript quantification · live
Sister-domain expansion of Phase 18 / 18.5 / 18.6 closing the
transcript-level quantification gap that Phase 18.6 explicitly
deferred. Add the two de-facto transcript-level quantification tools
to Valenx: Salmon (Rob Patro's quasi-mapping plus two-phase EM
quantifier — the GTEx / TCGA / nf-core / GENCODE reference
quantifier) and Kallisto (Lior Pachter's pseudoalignment-based
quantifier — the original "skip the alignment" approach). Both
pseudo-align reads to a transcriptome and report TPM / count per
transcript without producing intermediate SAM / BAM files; they are
shape-distinct from the Phase 18.6 RNA-seq aligners (HISAT2, STAR)
because they emit per-transcript abundance tables rather than aligned
reads. Both adapters mirror the established Phase 18 BWA two-stage
shape — single-binary CLI subprocess, file in / file out, `index →
quant` pipeline. No new infrastructure. Phase 20 sits numerically
after Phase 19.5 and ships chronologically right after the Phase 18.6
RNA-seq alignment beachhead.

- 2 new adapter crates land under `crates/valenx-adapters/bio/`:
  Salmon (single-binary subprocess; two-stage `salmon index → salmon
  quant` pipeline; `case.toml` knobs `transcriptome` (FASTA, required
  — used to build the index when `skip_index = false`), `index_dir`
  (required), `reads` (1 entry single-end, 2 entries paired-end),
  `output_dir` (required), `threads` default 1, `skip_index` default
  false, `libtype` default `"A"` (Salmon's libtype DSL — `"A"` auto-
  detects orientation, `"U"` unstranded, `"ISF"` / `"ISR"` paired-end
  stranded forward / reverse, `"IU"` paired-end unstranded — left
  non-whitelisted because the libtype DSL has many valid combos),
  `extra_args`; collects `quant.sf` typed `Tabular` "Salmon transcript
  quantification" and `cmd_info.json` typed `Log` "Salmon command
  info"; `bio.salmon.quant` ribbon capability; GPL-3.0 licensed),
  Kallisto (single-binary subprocess; two-stage `kallisto index →
  kallisto quant` pipeline; index is a single `.idx` file (kallisto
  convention — not a directory like Salmon / STAR); `case.toml` knobs
  `transcriptome` (FASTA, required), `index` (single `.idx` file,
  required), `reads` (1 or 2 entries), `output_dir` (required),
  `threads` default 1, `skip_index` default false, `fragment_length`
  + `fragment_sd` optional `f64` pair (both required when
  `reads.len() == 1` — kallisto auto-detects fragment statistics from
  paired-end reads but cannot for single-end; both must be finite and
  `> 0.0`), `extra_args`; collects `abundance.tsv` typed `Tabular`
  "Kallisto transcript abundance", `abundance.h5` typed `Native`
  "Kallisto HDF5 abundance", and `run_info.json` typed `Log`
  "Kallisto run info"; `bio.kallisto.quant` ribbon capability;
  BSD-2-Clause licensed). Each wired into
  `valenx-app::init_registry`.
- No new canonical types, no new format readers, no new CLIs. Both
  adapters consume the existing Phase 17 FASTA + Phase 18 FASTQ
  inputs and emit per-transcript abundance tables that the unchanged
  `Results.artifacts` collection model surfaces directly through the
  `Tabular` artifact kind. Kallisto's `abundance.h5` HDF5 sidecar
  needs an external tool (`h5dump` / Python `h5py`) to inspect — the
  canonical H5 reader as a Valenx CLI defers to Phase 19.6 along with
  the Seurat / AnnData R-runtime work.
- 2 new `valenx-init` templates ship: `salmon` (`salmon-quant`) and
  `kallisto` (`kallisto-quant`). Cross-binary roundtrip test sweeps
  all 57 templates clean.

The full per-phase shape lives in
[docs/src/phases/phase-20-transcript-quantification.md](./docs/src/phases/phase-20-transcript-quantification.md);
the implementation plan lives at
[docs/superpowers/plans/2026-04-30-transcript-quantification.md](./docs/superpowers/plans/2026-04-30-transcript-quantification.md).

**Deliverable:** Adapter inventory at 61 of 62 fully live after this
phase (2 new transcript-quantification adapters added on top of the
Phase 17 + Phase 17.5 + Phase 18 + Phase 18.5 + Phase 18.6 + Phase 19
+ Phase 19.5 + Phase 22 + Phase 23 + Phase 24 + Phase 27 + Phase 27.5
+ Phase 34 totals); transcript-quantification workflows land alongside
the existing biology adapters in the same case-toml / prepare / run /
collect shell, broadening the search → align → quantify → predict →
validate loop to eleven alignment / search tools (BWA, Bowtie2,
HISAT2, STAR, minimap2, MAFFT, MUSCLE, HMMER, samtools, MMseqs2,
DIAMOND) plus two transcript quantifiers (Salmon, Kallisto) feeding
into the five Phase 17 + 17.5 prediction tools.

### Phase 22 — Workflow managers · live
Add the two de-facto bioinformatics workflow orchestrators to Valenx:
Nextflow (the DSL-driven pipeline language behind nf-core) and
Snakemake (the Python-flavoured rule-based orchestrator). Unlike the
per-tool adapters Phase 17 / 18 / 19 / 23 / 24 ship (BWA, samtools,
bcftools, …), these are meta-tools — they invoke pipelines that
themselves call other bio adapters' underlying binaries. The Valenx
adapter just orchestrates the orchestrator, keeping the rest of the
registry useful underneath. Both follow the established Phase 18 BWA
single-binary CLI shape: probe / prepare / run / collect, output-in-
workdir. Phase 22 sits numerically before Phase 23 but ships
chronologically after Phase 24 — the same chronological-vs-numerical
convention used for Phase 17.5.

- 2 new adapter crates land under `crates/valenx-adapters/bio/`:
  Nextflow (DSL-driven pipeline orchestrator behind nf-core; single-
  binary CLI shape `nextflow run <pipeline> [-c <config>] [-profile
  <profile>] [-resume] [--<key> <value>...] [extras...]`; `pipeline`
  accepts a local `.nf`, an absolute path, or a registry identifier
  like `nf-core/rnaseq`; `params` map maps to `--<key> <value>`;
  collects the workdir as a `Native` artifact `"Nextflow run workdir"`
  and walks `report.html` / `timeline.html` / `dag.svg` as `Log`
  artifacts; Apache-2.0 licensed; `bio.nextflow.run` ribbon
  capability) and Snakemake (Python-flavoured rule-based pipeline
  orchestrator; single-binary CLI shape `snakemake -s <snakefile>
  --cores N [--use-conda] [-n] [--configfile <path>] [<targets>...]`;
  `targets` lists specific rules to build (empty = all default
  targets); `cores` (≥ 1) sets parallelism; `use_conda` toggles
  managed envs; `dry_run` toggles plan-only `-n`; collects the
  workdir as a `Native` artifact `"Snakemake run workdir"` and walks
  `.snakemake/log/*.log`; MIT licensed; `bio.snakemake.run` ribbon
  capability). Each wired into `valenx-app::init_registry`.
- No new canonical types, no new format readers, no new CLIs.
  Workflow managers are meta-orchestrators — they don't produce a
  single canonical artifact of their own; the pipelines they invoke
  produce whatever the underlying tools do (BAM via BWA, VCF via
  bcftools, FASTA via ColabFold, …), and the unchanged
  `Results.artifacts` collection model surfaces them through their
  respective adapters' canonical types.
- 2 new `valenx-init` templates ship: `nextflow` with alias `nf`
  (`nextflow-pipeline`), and `snakemake` with alias `smk`
  (`snakemake-pipeline`). Cross-binary roundtrip test sweeps all 45
  templates clean.

The full per-phase shape lives in
[docs/src/phases/phase-22-workflow-managers.md](./docs/src/phases/phase-22-workflow-managers.md);
the implementation plan lives at
[docs/superpowers/plans/2026-04-30-workflow-managers.md](./docs/superpowers/plans/2026-04-30-workflow-managers.md).

**Deliverable:** Adapter inventory at 49 of 50 fully live after this
phase (2 new workflow-manager adapters added on top of the Phase 17 +
Phase 17.5 + Phase 18 + Phase 19 + Phase 23 + Phase 24 + Phase 27 +
Phase 34 totals — Phase 22 ships chronologically after Phase 24 even
though it sits numerically before Phase 23); workflow-manager
pipelines drive the existing biology adapters' underlying binaries
through the same case-toml / prepare / run / collect shell.

### Phase 22.5 — Workflow expansion · live ✅
Sister-adapter expansion of the Phase 22 Nextflow + Snakemake
workflow-manager pair. Round out the bio workflow-orchestration
surface that Phase 22 opened with three more canonical workflow
tools that cover the corners Nextflow / Snakemake don't reach —
Galaxy ecosystem CLI for tool development + workflow execution
outside a full Galaxy server (planemo, the official Galaxy command-
line companion that lints tool wrappers, runs Galaxy workflow
tests, and executes `.ga` / `.gxwf.yml` workflows; AFL-3.0; single-
binary subprocess shape sister to Phase 22 Nextflow / Snakemake),
Broad Institute Workflow Description Language engine (Cromwell, the
canonical WDL runner powering most production GATK + Terra
pipelines; BSD-3-Clause; **JAR-distributed** sister to Phase 33 j5
/ Cello / Phase 41 Jalview), and Common Workflow Language reference
runner (cwltool, the cross-tool standard implementation for
describing analytical workflows in YAML / JSON; Apache-2.0; single-
binary subprocess shape sister to Phase 22 Snakemake). planemo +
cwltool follow the established Phase 22 Nextflow / Snakemake
single-binary CLI pattern: workflow file in, run artifacts out.
Cromwell is JAR-distributed; the user supplies the absolute path to
`cromwell-<version>.jar` via case input, and we probe `java`
itself rather than the JAR. Phase 22.5 sits numerically adjacent to
Phase 22 and ships chronologically right after Phase 41 sequence
editors — same chronological-vs-numerical convention used for Phase
17.5 / 24 / 28 / 31 / 35 / 39 / 5.5 / 5.6 / 5.7 / 32.5 / 40 / 41.

- 3 new adapter crates land under `crates/valenx-adapters/bio/`:
  planemo (single-binary subprocess sister to Phase 22 Nextflow /
  Snakemake; `case.toml` knobs `workflow` (`.ga` / `.gxwf.yml`;
  required) / `inputs` (`Option<PathBuf>` — optional inputs JSON;
  `None` for workflows that take no inputs) / `output_basename`
  (filename stem `collect()` uses to filter HTML reports;
  required, non-empty) / `action` (string; default `"run"` —
  `prepare()` rejects values other than `run` / `test` / `lint`
  at parse time so the adapter doesn't forward unsupported sub-
  commands) / `extra_args`; `prepare()` resolves both `workflow`
  and the optional `inputs` against the case directory when
  relative, validates each file exists on disk, builds `planemo
  <action> <workflow> [inputs] [extras...]`; `collect()` walks
  the workdir for `<output_basename>*.html` (`Native`, "Planemo
  report"), `*.json` (`Tabular`, "Planemo run JSON"), `*.log`
  (`Log`); probe via `find_on_path(&["planemo"])`; AFL-3.0
  licensed; version range `0.75.0..1.0.0`; `bio.planemo.run`
  ribbon capability), Cromwell (JAR-distributed single-binary
  subprocess sister to Phase 33 j5 / Cello / Phase 41 Jalview;
  `case.toml` knobs `jar` (absolute path to
  `cromwell-<version>.jar`; required) / `workflow` (`.wdl`
  workflow file; required) / `inputs` (`Option<PathBuf>` —
  optional inputs JSON; emitted as **two separate args** `-i` +
  `<inputs>` only when `Some` — the flag is suppressed entirely
  when `None`) / `output_basename` / `action` (string; default
  `"run"` — `prepare()` rejects values other than `run` /
  `submit` / `validate` at parse time) / `extra_args`;
  `prepare()` resolves all three paths against the case
  directory when relative, validates each file exists on disk,
  builds `java -jar <jar> <action> <workflow> [-i <inputs>]
  [extras...]`; `collect()` walks **the top level only** of the
  workdir for `<output_basename>*.json` (`Tabular`, "Cromwell
  metadata"), `*.log` (`Log`); probe via
  `find_on_path(&["java"])` — Cromwell's version comes from the
  jar itself, not from `java`, so we surface no version here
  (same shape as Phase 33 j5 / Cello / Phase 41 Jalview);
  BSD-3-Clause licensed; version range `80.0.0..100.0.0`;
  `bio.cromwell.run` ribbon capability), and cwltool (single-
  binary subprocess sister to Phase 22 Snakemake; `case.toml`
  knobs `workflow` (`.cwl` tool / workflow document; required) /
  `inputs` (`Option<PathBuf>` — optional CWL input-object
  document in JSON or YAML) / `output_dir` (`--outdir` target
  subdirectory under the workdir; required, non-empty) /
  `extra_args`; `prepare()` resolves both files against the case
  directory when relative, validates each file exists on disk,
  builds `cwltool --outdir <output_dir> [extras...] <workflow>
  [inputs]`; `collect()` walks **one level deep** into
  `<output_dir>/` for any file (`Native`, "cwltool output"),
  top-level `*.log` (`Log`); probe prefers the `cwltool`
  console-script entry-point and falls back to a Python-on-PATH
  detection with `"cwltool not found on PATH; install via pip
  install cwltool"` warning when only Python is present;
  Apache-2.0 licensed; version range `3.1.0..4.0.0`;
  `bio.cwltool.run` ribbon capability). Each wired into
  `valenx-app::init_registry`.
- No new canonical types, no new format readers, no new CLIs.
  Workflow managers are meta-orchestrators — the per-tool
  outputs are surfaced through their respective adapters'
  canonical types.
- 3 new `valenx-init` templates ship: `planemo`, `cromwell`,
  `cwltool`. Cross-binary roundtrip test sweeps all 120 templates
  clean.

The full per-phase shape lives in
[docs/src/phases/phase-22-5-workflow-expansion.md](./docs/src/phases/phase-22-5-workflow-expansion.md);
the implementation plan lives at
[docs/superpowers/plans/2026-05-04-workflow-expansion.md](./docs/superpowers/plans/2026-05-04-workflow-expansion.md).

**Deliverable:** Adapter inventory at 122 of 123 fully live after
this phase trio (3 new workflow-expansion adapters added on top of
the prior set, alongside the Phase 42 web-visualization pair that
brings the total to **124 of 125 fully live**); Phase 22.5 sister-
expands the Phase 22 workflow-manager pair with three more
workflow languages (Galaxy / WDL / CWL).

### Phase 23 — Molecular viewers · live
Round out the visualization surface for everything Valenx's biology
stack produces. Phase 17 shipped ChimeraX as the first script-driven
molecular renderer; Phase 23 ships its three most-used siblings,
following the same shape (script-driven subprocess, headless mode,
output-in-workdir).

- 3 new adapter crates land under `crates/valenx-adapters/bio/`:
  PyMOL (open-source build — the Schrödinger fork is proprietary;
  drives off `.pml` Python-style command files; defaults to
  `pymol -c -q <script>` for headless quiet rendering), VMD
  (Tcl-scripted MD trajectory viewer — `vmd -dispdev text -e
  <script>`; optional `structure` field stages a `.pdb` / `.gro` /
  `.psf` as a positional arg), and IGV (`igvtools` wrapper for
  headless BAM / VCF / WIG indexing — per-action dispatch on
  `action ∈ {index, count, sort, tile}`; the companion GUI
  viewer is out of scope, this is the headless-tooling adapter
  only). Each wired into `valenx-app::init_registry`.
- VMD's adapter pushes a license-awareness warning into
  `ProbeReport.warnings` because VMD ships under a custom
  non-OSS-but-free-for-academic-use license. The warning contains
  the keyword `"academic"` so users see it before they ship
  renders or derived data downstream.
- No new canonical types, no new format readers, no new CLIs.
  All three viewers consume existing Phase 17 / 18 / 19 inputs
  (`Structure` / PDB, BAM, trajectories) and emit user-readable
  artifacts (PNG / PSE / TGA renders, BAI / IDX index sidecars)
  that the unchanged `Results.artifacts` collection model surfaces
  directly.
- 3 new `valenx-init` templates ship: `pymol` (`pymol-render`),
  `vmd` (`vmd-render`), and `igv` with alias `igvtools` (`igv-
  index`). Cross-binary roundtrip test sweeps all 36 templates
  clean.

The full per-phase shape lives in
[docs/src/phases/phase-23-viewers.md](./docs/src/phases/phase-23-viewers.md);
the implementation plan lives at
[docs/superpowers/plans/2026-04-30-molecular-viewers.md](./docs/superpowers/plans/2026-04-30-molecular-viewers.md).

**Deliverable:** Adapter inventory at 40 of 41 fully live after this
phase (3 new molecular-viewer adapters added on top of the Phase 17 +
Phase 17.5 + Phase 18 + Phase 19 totals); molecular-viewer workflows
land alongside the existing biology adapters in the same case-toml /
prepare / run / collect shell.

### Phase 24 — Cheminformatics expansion · live
Round out the cheminformatics surface that Phase 17's RDKit adapter
started. Phase 24 ships three sister adapters (DeepChem + Open Babel +
Avogadro 2) that together with RDKit (Phase 17) and the Phase 34
docking pair give Valenx the complete small-molecule + cheminformatics
stack. All three follow established patterns: DeepChem mirrors RDKit's
Python-script subprocess shape, Open Babel uses BWA's single-binary
CLI shape (`obabel <in> -O <out>`), and Avogadro 2 mirrors ChimeraX's
script-driven-headless pattern. Phase 24 sits numerically between
Phase 23 and Phase 27 but ships chronologically after Phase 34 — the
same chronological-vs-numerical convention used for Phase 17.5.

- 3 new adapter crates land under `crates/valenx-adapters/bio/`:
  DeepChem (PyTorch-backed deep-learning cheminformatics; user-
  provided Python script + `valenx_params.json` knobs — optional
  inline `smiles` list, optional `dataset_csv`, optional `checkpoint`;
  walks workdir for `.csv` (`Tabular`, "DeepChem analysis output"),
  `.png` ("DeepChem plot"), `.pkl` / `.pt` ("DeepChem model
  checkpoint"); MIT licensed; `bio.deepchem.script` ribbon
  capability), Open Babel (de-facto open-source chemistry-format
  converter, ~120 formats; `obabel <input> -O <output>` single-
  binary CLI with explicit `input_format` / `output_format`
  overrides, `gen_3d` toggles `--gen3D` for 2D → 3D coord generation,
  `add_hydrogens` toggles `-h`; collects converted output as a
  `Native` artifact "Open Babel converted file"; GPL-2.0 licensed;
  `bio.openbabel.convert` ribbon capability), and Avogadro 2
  (Python-scriptable chemistry editor + small-molecule rendering
  pipeline; `avogadro2 --script <script.py>` with optional
  `structure` (`.cml` / `.mol` / `.xyz` / `.pdb`) staged as a
  positional arg, `headless` (default true) toggles `--no-gui`;
  collects `.png` ("Avogadro 2 render") and `.cml` / `.mol` /
  `.xyz` ("Avogadro 2 exported structure"); GPL-2.0-or-later
  licensed; `bio.avogadro.render` ribbon capability). Each wired
  into `valenx-app::init_registry`.
- No new canonical types, no new format readers, no new CLIs.
  All three adapters consume existing Phase 17 / 18 / 19 inputs
  (PDB / SDF / SMILES / CSV) and emit user-readable artifacts
  (CSV tables, PNG renders, exported chemistry structures, PyTorch /
  Pickle model checkpoints) that the unchanged `Results.artifacts`
  collection model surfaces directly.
- 3 new `valenx-init` templates ship: `deepchem` with alias `dc`
  (`deepchem-screen`), `openbabel` with alias `obabel`
  (`openbabel-convert`), and `avogadro` with alias `avogadro2`
  (`avogadro-render`). Cross-binary roundtrip test sweeps all 43
  templates clean.

The full per-phase shape lives in
[docs/src/phases/phase-24-cheminformatics.md](./docs/src/phases/phase-24-cheminformatics.md);
the implementation plan lives at
[docs/superpowers/plans/2026-04-30-cheminformatics-expansion.md](./docs/superpowers/plans/2026-04-30-cheminformatics-expansion.md).

**Deliverable:** Adapter inventory at 47 of 48 fully live after this
phase (3 new cheminformatics adapters added on top of the Phase 17 +
Phase 17.5 + Phase 18 + Phase 19 + Phase 23 + Phase 27 + Phase 34
totals — Phase 24 ships chronologically after Phase 34 even though it
sits numerically between Phase 23 and Phase 27); cheminformatics
workflows land alongside the existing biology adapters in the same
case-toml / prepare / run / collect shell.

### Phase 25 — Quantum chemistry · live
Open the **first quantum-chemistry domain** in Valenx with three
established open-source tools that span the quantum-chemistry tradeoff
space — semiempirical at the fast-and-approximate end, general-purpose
HF/DFT/post-HF in the middle, massively-parallel ab initio at the high
end: Psi4 (HF/DFT/post-HF general-purpose ab initio quantum chemistry,
LGPL-3.0; Psithon-scriptable input), NWChem (Pacific Northwest National
Lab's massively-parallel ab initio + plane-wave DFT package, ECL-2.0;
its own `.nw` input format with optional `mpirun` launcher when
`mpi_procs > 1`), and xTB (Stefan Grimme's extended tight-binding
semiempirical method, LGPL-3.0; reads `.xyz` coordinates directly with
all options on the CLI). Psi4 follows the Phase 18 BWA single-binary
CLI shape with explicit `-i`/`-o` arguments; NWChem follows the same
BWA shape with optional MPI wrapping plus the MAFFT-style stdout-
redirect pattern (NWChem writes its run report to stdout, captured to
`output`); xTB follows the BWA shape with stdout captured to `xtb.log`
via the same MAFFT-style stdout-redirect pattern. Phase 25 sits
numerically between Phase 24 and Phase 27 but ships chronologically
right after Phase 28 — the same chronological-vs-numerical convention
used for Phase 17.5 / 24 / 28.

- 3 new adapter crates land under `crates/valenx-adapters/bio/`:
  Psi4 (open-source HF/DFT/post-HF quantum chemistry; `case.toml`
  knobs `input` (`.in` / `.dat` Psithon script; required), `output`
  (output filename relative to workdir; required), `threads` (default
  1; ≥ 1), `memory` (default `"1 gb"`; matches
  `^\d+\s*(mb|gb|MB|GB)$` via `is_valid_memory` helper), `extra_args`;
  `prepare()` builds `psi4 -i <input> -o <output> -n <threads>
  [-m <memory>] [extras...]` (`-m` only emitted when `memory` is non-
  default to preserve Psi4's own `"500 mb"` internal default);
  collects `output` as `Log` artifact "Psi4 output" plus `.fchk`
  (`Native`, "Psi4 formatted checkpoint") and `.molden` (`Native`,
  "Psi4 Molden orbital data") files; probe via
  `find_on_path(&["psi4"])`; LGPL-3.0 licensed; `bio.psi4.compute`
  ribbon capability), NWChem (Pacific Northwest National Lab's
  massively-parallel ab initio + plane-wave DFT package; `case.toml`
  knobs `input` (`.nw` NWChem-format script; required), `output`
  (required), `mpi_procs` (default 1; ≥ 1), `extra_args`; `prepare()`
  builds — serial: `nwchem [extras...] <input>`; parallel:
  `mpirun -n <mpi_procs> [extras...] nwchem <input>`; when
  `mpi_procs > 1`, prepare resolves `mpirun` via `find_on_path` and
  fails with a helpful install-hint `InvalidCase` if it's missing
  (`apt install openmpi-bin` / `apt install mpich`) rather than
  letting the child fail later with a less obvious "command not
  found"; output path is stashed in
  `PreparedJob.environment[VALENX_NWCHEM_OUTPUT]` so `run()` can
  redirect stdout to it via the MAFFT-style stdout-redirect pattern
  without re-parsing the case TOML; collects `output` as `Log`
  artifact "NWChem output"; probe via `find_on_path(&["nwchem"])`;
  ECL-2.0 licensed; `bio.nwchem.compute` ribbon capability), and xTB
  (Stefan Grimme's extended tight-binding semiempirical quantum
  chemistry; `case.toml` knobs `input` (`.xyz` geometry; required),
  `mode` ∈ `{single-point, opt, ohess, hess, md}` (default
  `"single-point"` — xTB's default run type, no flag emitted; every
  other mode maps to `--<mode>`), `charge` (electron-balance `i32`;
  default 0), `uhf` (xTB's multiplicity convention — number of
  unpaired electrons, `u32`; default 0), `gfn` ∈ `{0, 1, 2}` (GFN
  method; default 2 — GFN2-xTB), `solvent` (optional ALPB solvent
  name e.g. `"water"` / `"thf"`; `None` = gas phase), `extra_args`;
  `prepare()` builds `xtb <input> --gfn <gfn> --chrg <charge> --uhf
  <uhf> [--<mode> if mode != "single-point"] [--alpb <solvent> if
  Some] [extras...]`; collects `xtb.log` as `Log` artifact "xTB
  stdout log" plus `xtbopt.xyz` (`Native`, "xTB optimized geometry"),
  `xtbopt.log` (`Log`), `gradient` / `hessian` files (`Native`);
  probe via `find_on_path(&["xtb"])`; LGPL-3.0 licensed;
  `bio.xtb.compute` ribbon capability). Each wired into
  `valenx-app::init_registry`.
- No new canonical types, no new format readers, no new CLIs. All
  three adapters consume user-supplied input files (Psithon `.in` /
  NWChem `.nw` / xyz `.xyz` coordinates) and emit user-readable
  artifacts (text output reports, formatted checkpoints, Molden
  orbital data, optimized xyz geometries, gradient / hessian files)
  that the unchanged `Results.artifacts` collection model surfaces
  directly. A first-class quantum-chemistry canonical type — a
  generic energy / geometry / orbital data type spanning all three
  back-ends — defers to a future phase along with `.fchk` / Molden /
  `.cube` reader CLIs and visualization integrations.
- 3 new `valenx-init` templates ship: `psi4` (`psi4-compute`),
  `nwchem` (`nwchem-compute`), and `xtb` (`xtb-compute`). Canonical
  names only — no aliases beyond the canonical names themselves.
  Cross-binary roundtrip test sweeps all 66 templates clean.

The full per-phase shape lives in
[docs/src/phases/phase-25-quantum-chemistry.md](./docs/src/phases/phase-25-quantum-chemistry.md);
the implementation plan lives at
[docs/superpowers/plans/2026-04-30-quantum-chemistry.md](./docs/superpowers/plans/2026-04-30-quantum-chemistry.md).

**Deliverable:** Adapter inventory at 70 of 71 fully live after this
phase (3 new quantum-chemistry adapters added on top of the Phase 17 +
Phase 17.5 + Phase 18 + Phase 18.5 + Phase 18.6 + Phase 19 + Phase
19.5 + Phase 20 + Phase 22 + Phase 23 + Phase 24 + Phase 27 + Phase
27.5 + Phase 28 + Phase 30 + Phase 34 totals — Phase 25 ships
chronologically after Phase 28 even though it sits numerically between
Phase 24 and Phase 27); quantum-chemistry workflows land alongside the
existing biology + chemistry adapters in the same case-toml / prepare
/ run / collect shell, opening the **first quantum-chemistry domain**
in Valenx and broadening the build-geometry → optimize → compute-
energy → predict-structure → fold-RNA → infer-tree → validate loop
with three quantum-chemistry tools (Psi4, NWChem, xTB) feeding into
the Phase 24 cheminformatics surface (DeepChem, Open Babel, Avogadro
2), the Phase 17 / 17.5 prediction stack (ESMFold, OpenFold, AlphaFold
2/3, ColabFold), the Phase 28 RNA-structure tools (ViennaRNA,
RNAstructure, NUPACK), and the Phase 30 phylogenetic-tree builders
(IQ-TREE, RAxML-NG, FastTree).

### Phase 27 — Protein design · live
Pair the structure-prediction adapters Valenx already ships
(Phases 17 + 17.5: ColabFold, ESMFold, OpenFold, AlphaFold 2/3) with
their de novo design counterparts. Phase 27 ships RFdiffusion (GPU-
driven protein backbone generation) and ProteinMPNN (sequence design
from backbone). Together with the prediction stack, this gives Valenx
the complete **design → predict → validate** loop in one shell.

- 2 new adapter crates land under `crates/valenx-adapters/bio/`:
  RFdiffusion (Python-script subprocess; `valenx_params.json` for
  config knobs; `mode ∈ {motif, binder, unconditional, partial-
  diffusion}`; `num_designs` and `diffusion_steps` defaults of 8
  and 50 respectively; collects `<output_basename>_*.pdb` typed via
  `valenx_bio::format::pdb::read` — RFdiffusion writes pLDDT into
  the B-factor column too; BSD-3-Clause licensed; `bio.rfdiffusion.
  design` ribbon capability) and ProteinMPNN (same Python-script-
  subprocess shape; `model_variant ∈ {vanilla, soluble, ca-only}`;
  `temperature` default 0.1 and `num_seq_per_target` default 8;
  collects `<output_basename>.fa` parsed via `valenx_bio::format::
  fasta::read_str` for a richer `"ProteinMPNN · N sequences"`
  artifact label; MIT licensed; `bio.proteinmpnn.design` ribbon
  capability). Each wired into `valenx-app::init_registry`.
- No new canonical types, no new format readers, no new CLIs.
  Both adapters consume the existing Phase 17 PDB inputs and emit
  user-readable artifacts (PDB backbones, FASTA sequences) that the
  unchanged `Results.artifacts` collection model surfaces directly.
  RFdiffusion's PDB outputs are inspectable through the existing
  `valenx-pdb-info` CLI; ProteinMPNN's FASTA outputs through the
  existing `valenx-fasta` CLI.
- 2 new `valenx-init` templates ship: `rfdiffusion` with alias `rfd`
  (`rfdiffusion-design`), and `proteinmpnn` with alias `mpnn`
  (`proteinmpnn-design`). Cross-binary roundtrip test sweeps all
  38 templates clean.

The full per-phase shape lives in
[docs/src/phases/phase-27-protein-design.md](./docs/src/phases/phase-27-protein-design.md);
the implementation plan lives at
[docs/superpowers/plans/2026-04-30-protein-design.md](./docs/superpowers/plans/2026-04-30-protein-design.md).

**Deliverable:** Adapter inventory at 42 of 43 fully live after this
phase (2 new protein-design adapters added on top of the Phase 17 +
Phase 17.5 + Phase 18 + Phase 19 + Phase 23 totals); protein-design
workflows land alongside the existing biology adapters in the same
case-toml / prepare / run / collect shell, closing the design →
predict → validate loop with no glue code.

### Phase 27.5 — Protein design expansion · live
Sister-adapter expansion of Phase 27. Add three more open-source
protein design tools to round out the de novo design surface:
Chroma (Generate Biomedicines' joint backbone + sequence diffusion
model), ESM-IF (Meta's GVP-based inverse-folding sequence designer
— alternative to ProteinMPNN), and RFantibody (RosettaCommons
antibody-specific RFdiffusion fork). All three follow the
established Phase 27 RFdiffusion / ProteinMPNN pattern — Python-
script subprocess with `valenx_params.json` for config knobs. No
new infrastructure. Phase 27.5 sits numerically after Phase 27
and ships chronologically right after Phase 27 — same convention
as Phase 17.5 sits between Phase 17 and Phase 18 numerically.

- 3 new adapter crates land under `crates/valenx-adapters/bio/`:
  Chroma (Python-script subprocess; `valenx_params.json` knobs
  `num_samples` default 4, `length`, `temperature` default 1.0,
  `output_basename`; collects `<output_basename>*.pdb` typed
  `Native` "Chroma design" and `<output_basename>*.fa` typed
  `Tabular` "Chroma sequence"; Apache-2.0 licensed;
  `bio.chroma.design` ribbon capability), ESM-IF (same Python-
  script-subprocess shape; `valenx_params.json` knobs `input_pdb`,
  `model` default `esm_if1_gvp4_t16_142M_UR50`, `temperature`
  default 1.0, `num_samples` default 8, `output_basename`;
  collects `<output_basename>.fa` parsed via
  `valenx_bio::format::fasta::read` for a richer `"ESM-IF · N
  sequences"` artifact label; MIT licensed via the `fair-esm`
  package; `bio.esm-if.design` ribbon capability), and RFantibody
  (same Python-script-subprocess shape; `valenx_params.json` knobs
  `framework_pdb`, `target_pdb`, `design_loops` ⊆ `{H1, H2, H3,
  L1, L2, L3}`, `num_designs` default 8, `diffusion_steps` default
  50, `output_basename`; collects `<output_basename>*.pdb` typed
  `Native` "RFantibody design"; BSD-3-Clause licensed;
  `bio.rfantibody.design` ribbon capability). Each wired into
  `valenx-app::init_registry`.
- No new canonical types, no new format readers, no new CLIs.
  All three adapters consume the existing Phase 17 PDB inputs
  and emit user-readable artifacts (PDB backbones, FASTA
  sequences) that the unchanged `Results.artifacts` collection
  model surfaces directly. Chroma + RFantibody PDB outputs are
  inspectable through the existing `valenx-pdb-info` CLI; ESM-IF
  FASTA outputs through the existing `valenx-fasta` CLI.
- 3 new `valenx-init` templates ship: `chroma` (`chroma-design`),
  `esm-if` with aliases `esmif` / `inverse-folding`
  (`esm-if-design`), and `rfantibody` with alias `rfab`
  (`rfantibody-design`). Cross-binary roundtrip test sweeps all
  50 templates clean.

The full per-phase shape lives in
[docs/src/phases/phase-27-5-protein-design-expansion.md](./docs/src/phases/phase-27-5-protein-design-expansion.md);
the implementation plan lives at
[docs/superpowers/plans/2026-04-30-protein-design-expansion.md](./docs/superpowers/plans/2026-04-30-protein-design-expansion.md).

**Deliverable:** Adapter inventory at 54 of 55 fully live after
this phase (3 new protein-design-expansion adapters added on top
of the Phase 17 + Phase 17.5 + Phase 18 + Phase 19 + Phase 19.5 +
Phase 22 + Phase 23 + Phase 24 + Phase 27 + Phase 34 totals);
protein-design-expansion workflows land alongside the existing
biology adapters in the same case-toml / prepare / run / collect
shell, broadening the design → predict → validate loop to five
design tools (RFdiffusion, Chroma, RFantibody for backbone /
antibody design; ProteinMPNN, ESM-IF for sequence design) feeding
into the five Phase 17 + 17.5 prediction tools.

### Phase 27.6 — EvolutionaryScale models · live
Complete EvolutionaryScale's open-source ESM lineup in Valenx.
Phase 17.5 + 27.5 already shipped ESMFold (single-sequence
structure prediction) and ESM-IF (GVP-based inverse-folding
sequence design). Phase 27.6 adds the remaining two: ESM3
(EvolutionaryScale's flagship generative multi-modal protein model
— joint reasoning over sequence + structure + function tracks,
modes `design` / `inverse-fold` / `scaffold` / `predict`) and ESM
Cambrian / ESMC (the smaller-faster protein representation model
for embedding-driven downstream ML, two open checkpoints
`esmc-300m` / `esmc-600m`). Both follow the established Phase
17.5 ESMFold / Phase 27.5 ESM-IF pattern — Python-script
subprocess where the user's script imports the EvolutionaryScale
`esm` package (the same one ESMFold + ESM-IF already pull in) and
reads `valenx_params.json` (auto-written by the adapter) for
config knobs. No new infrastructure. Phase 27.6 sits numerically
after Phase 27.5 and ships chronologically right after Phase 25
quantum chemistry — same convention as Phase 17.5 sits between
Phase 17 and Phase 18 numerically.

- 2 new adapter crates land under `crates/valenx-adapters/bio/`:
  ESM3 (Python-script subprocess; `valenx_params.json` knobs
  `model_variant` ∈ `{open, open-multimer, small}`, `mode` ∈
  `{design, inverse-fold, scaffold, predict}`, `num_samples`
  default 4, optional `input_pdb` (required for `inverse-fold` /
  `scaffold`), optional `input_fasta` (required for `predict`),
  `temperature` default 1.0, `output_basename`; collects
  `<output_basename>*.pdb` typed `Native` "ESM3 generated
  structure" and `<output_basename>*.fa` typed `Tabular` "ESM3
  generated sequence"; Cambrian-Open-License — open weights for
  the smaller checkpoints, non-commercial for the largest Forge-
  only variants; `bio.esm3.generate` ribbon capability), ESMC
  (same Python-script-subprocess shape; `valenx_params.json` knobs
  `input_fasta` (required), `model_variant` ∈
  `{esmc-300m, esmc-600m}`, `pooling` ∈ `{per-residue, mean}`,
  `output_basename`; collects `<output_basename>.{npy,npz,parquet}`
  typed `Tabular` "ESMC embeddings"; Cambrian-Open-License;
  `bio.esmc.embed` ribbon capability). Each wired into
  `valenx-app::init_registry`. Both ride the same EvolutionaryScale
  `esm` Python package as ESMFold / ESM-IF — installing one
  installs them all, and the probe surfaces a single unified "esm
  is importable" gate via the shared `detect_esm_version` helper.
- No new canonical types, no new format readers, no new CLIs.
  Both adapters consume the existing Phase 17 PDB / FASTA inputs
  and emit user-readable artifacts (PDB backbones, FASTA
  sequences, NumPy `.npy` / `.npz` and Parquet embedding tables)
  that the unchanged `Results.artifacts` collection model surfaces
  directly. ESM3 PDB outputs are inspectable through the existing
  `valenx-pdb-info` CLI; ESM3 FASTA outputs through the existing
  `valenx-fasta` CLI; ESMC's embedding sidecars feed straight into
  the user's downstream Python pipeline.
- 2 new `valenx-init` templates ship: `esm3` (`esm3-generate`)
  and `esmc` with alias `esm-cambrian` (`esmc-embed`). Cross-
  binary roundtrip test sweeps all 68 templates clean.

The full per-phase shape lives in
[docs/src/phases/phase-27-6-evolutionaryscale.md](./docs/src/phases/phase-27-6-evolutionaryscale.md);
the implementation plan lives at
[docs/superpowers/plans/2026-04-30-evolutionaryscale.md](./docs/superpowers/plans/2026-04-30-evolutionaryscale.md).

**Deliverable:** Adapter inventory at 72 of 73 fully live after
this phase (2 new EvolutionaryScale-model adapters added on top
of the Phase 17 + Phase 17.5 + Phase 18 + Phase 18.5 + Phase 18.6
+ Phase 19 + Phase 19.5 + Phase 20 + Phase 22 + Phase 23 +
Phase 24 + Phase 25 + Phase 27 + Phase 27.5 + Phase 28 + Phase 30
+ Phase 34 totals); Phase 27.6 closes out the open-source
EvolutionaryScale ESM lineup at 4 of 4 tools (ESMFold + ESM-IF +
ESM3 + ESMC), broadening the design → predict → embed → infer →
validate loop to seven design tools (RFdiffusion, Chroma,
RFantibody for backbone / antibody design; ProteinMPNN, ESM-IF
for sequence design; ESM3 in either generative or scaffold-fill
mode), five prediction tools (ColabFold, ESMFold, OpenFold,
AlphaFold 2, AlphaFold 3, plus ESM3 in `predict` mode), and one
embedding workhorse (ESMC).

### Phase 28 — RNA structure · live
Open the RNA secondary-structure prediction domain in Valenx with three
established tools: ViennaRNA (the most-cited RNA secondary-structure
suite — `RNAfold` minimum-free-energy folding; custom non-commercial /
academic-use license, flagged via probe warning), RNAstructure
(Mathews lab's classic RNA folding toolkit — `Fold` is the flagship;
BSD-3-Clause), and NUPACK (Caltech's nucleic-acid package — academic-
license-only, flagged via probe warning the same way as VMD /
AlphaFold 3 / ChimeraX). ViennaRNA follows the MAFFT-style stdout-
redirect pattern (RNAfold writes to stdout); RNAstructure follows the
BWA single-binary CLI pattern with explicit `-o`-style output; NUPACK
follows the OpenMM / Scanpy Python-script-subprocess pattern (NUPACK 4
is Python-driven — the 3.x CLI is deprecated). Phase 28 sits
numerically between Phase 27.5 and Phase 30 and ships chronologically
right after the Phase 30 phylogenetics beachhead — the same
chronological-vs-numerical convention as Phase 17.5 sits numerically
between Phase 17 and Phase 18.

- 3 new adapter crates land under `crates/valenx-adapters/bio/`:
  ViennaRNA (single-binary subprocess with stdout-redirect — RNAfold
  writes the dot-bracket structure to stdout, captured to `output`
  via the MAFFT-style stdout-capture pattern; `case.toml` knobs
  `input` (FASTA), `output` (dot-bracket file), `temperature` Celsius
  default 37.0, `partition_function` default false (toggles `-p` for
  partition function + base-pair probabilities), `allow_gu` default
  true (`--noGU` disables GU pairs), `extra_args`; `prepare()` builds
  `RNAfold -i <input> -T <temperature> [-p] [--noGU if !allow_gu]
  [extras...]` → stdout captured to `output`; collects `output` as
  `Native` artifact "ViennaRNA secondary structure"; probe via
  `find_on_path(&["RNAfold"])` (capital R-N-A); custom non-commercial
  / academic license — probe pushes an `"academic"`-keyworded
  warning into `ProbeReport.warnings`; `bio.viennarna.fold` ribbon
  capability), RNAstructure (single-binary subprocess; binary
  literally named `Fold`; `case.toml` knobs `input` (FASTA or `.seq`
  RNAstructure-native format), `output` (`.ct` connection-table file),
  `max_structures` default 20 (≥ 1), `max_percent` % of MFE default
  10 (in `0..=100`), `temperature` Kelvin default 310.15 (> 0.0,
  finite), `extra_args`; `prepare()` builds `Fold <input> <output>
  -m <max_structures> -p <max_percent> -t <temperature> [extras...]`;
  collects `output` as `Native` artifact "RNAstructure connectivity
  table"; probe via `find_on_path(&["Fold"])`; BSD-3-Clause licensed;
  `bio.rnastructure.fold` ribbon capability), and NUPACK (Python-
  script subprocess shape — NUPACK 4 is Python-driven, the 3.x CLI is
  deprecated; user-supplied Python script imports `nupack` and reads
  `valenx_params.json` for config knobs; `case.toml` knobs `script`
  (required Python file), `python` default `python3`, `input`
  (optional FASTA / `.npc`), `output_basename`, `temperature` Celsius
  default 37.0, `sodium` molar default 1.0; `prepare()` stages the
  script + optional input, writes `valenx_params.json` with the
  staged filename / `output_basename` / `temperature` / `sodium`,
  builds `native_command = [python, script]`; collects
  `<output_basename>*` (`Native`, "NUPACK output") and `.npc` /
  `.json` files (`Tabular` / `Log`); probe via
  `find_on_path(&["python3", "python"])` then
  `python -c "import nupack"` — surfaces an install hint when Python
  is on PATH but `nupack` isn't importable; custom Caltech academic-
  only license — probe pushes an `"academic"`-keyworded warning into
  `ProbeReport.warnings`; `bio.nupack.analyze` ribbon capability).
  Each wired into `valenx-app::init_registry`.
- No new canonical types, no new format readers, no new CLIs. All
  three adapters consume the existing Phase 17 FASTA inputs and emit
  user-readable artifacts (dot-bracket structures, `.ct` connection
  tables, NUPACK output files) that the unchanged
  `Results.artifacts` collection model surfaces directly. A first-
  class RNA-secondary-structure canonical type with a dot-bracket or
  `.ct` reader as a Valenx CLI defers to a future phase along with
  visualization integrations.
- 3 new `valenx-init` templates ship: `viennarna` with aliases
  `vienna` / `rnafold` (`viennarna-fold`), `rnastructure`
  (`rnastructure-fold`), and `nupack` (`nupack-analyze`). Cross-
  binary roundtrip test sweeps all 63 templates clean.

The full per-phase shape lives in
[docs/src/phases/phase-28-rna-structure.md](./docs/src/phases/phase-28-rna-structure.md);
the implementation plan lives at
[docs/superpowers/plans/2026-04-30-rna-structure.md](./docs/superpowers/plans/2026-04-30-rna-structure.md).

**Deliverable:** Adapter inventory at 67 of 68 fully live after this
phase (3 new RNA-structure adapters added on top of the Phase 17 +
Phase 17.5 + Phase 18 + Phase 18.5 + Phase 18.6 + Phase 19 + Phase
19.5 + Phase 20 + Phase 22 + Phase 23 + Phase 24 + Phase 27 + Phase
27.5 + Phase 30 + Phase 34 totals); RNA-structure workflows land
alongside the existing biology adapters in the same case-toml /
prepare / run / collect shell, broadening the align → quantify →
predict → fold-RNA → infer-tree → validate loop to eleven alignment /
search tools (BWA, Bowtie2, HISAT2, STAR, minimap2, MAFFT, MUSCLE,
HMMER, samtools, MMseqs2, DIAMOND), two transcript quantifiers
(Salmon, Kallisto), five prediction tools (ESMFold, OpenFold,
AlphaFold 2, AlphaFold 3, ColabFold), three RNA-structure tools
(ViennaRNA, RNAstructure, NUPACK), and three phylogenetic-tree
builders (IQ-TREE, RAxML-NG, FastTree).

### Phase 29 — Population genetics · live
Open the **first population-genetics / evolutionary-simulation
domain** in Valenx with three established open-source tools that
span the population-genetics tradeoff space: SLiM (Philipp Messer's
forward-time individual-based population-genetics simulator —
evolves a finite-population model generation by generation under a
user-defined Eidos script (mutation rates, selection coefficients,
recombination maps, demographic events, migrations, mating
systems); GPL-3.0; single-binary subprocess shape sister to Phase
18 BWA with `slim [-s <seed>] [extras...] <script>`), msprime
(Jerome Kelleher's coalescent backwards-in-time population-genetics
simulator — speed-of-light coalescent simulator, millions of
samples per minute on a workstation; the canonical companion to
SLiM and tskit; GPL-3.0; Python-script subprocess shape sister to
Phase 17 Biopython), and tskit (the canonical tree-sequence
analysis library, MIT — built around the succinct tree-sequence
data structure pioneered by msprime; computes population-genetics
statistics (π, Tajima's D, Fst, site-frequency spectra, IBD
shares); the workhorse downstream of every Phase 29 simulator —
msprime emits `.trees`, SLiM emits `.trees`, tskit consumes them).
SLiM follows the established Phase 18 BWA single-binary CLI
pattern (script positional last so SLiM treats it as the model file
rather than the value of an earlier flag). msprime + tskit follow
the established Phase 17 Biopython Python-script subprocess pattern
(user authors a `.py` driver, the adapter stages script + writes
`valenx_params.json`, run() invokes `python <script>`). Phase 29
sits numerically between Phase 28 and Phase 30 and ships
chronologically right after Phase 35 CRISPR design — same
chronological-vs-numerical convention as Phase 17.5 / 24 / 28 / 31
/ 35.

- 3 new adapter crates land under `crates/valenx-adapters/bio/`:
  SLiM (single-binary subprocess sister to Phase 18 BWA;
  `case.toml` knobs `script` (`.slim` Eidos model file; required),
  `seed` (optional `u64`; passed via `slim -s <N>` when present,
  otherwise SLiM picks its own seed and prints it on the run
  banner), `output_basename` (filename stem the user's script uses
  for outputs — surfaced here so collect() can label artefacts
  uniformly even though SLiM scripts choose their own output paths;
  required, non-empty), `extra_args` (additional CLI arguments
  appended after the script path; `-d KEY=VALUE` is the canonical
  way to inject Eidos constants from outside the script);
  `prepare()` resolves the script against the case directory when
  relative, validates it exists on disk, and composes the
  invocation with `seed` injected before any extras and the script
  positional last; collects `<output_basename>*.trees` (`Native`,
  "SLiM tree sequence") and `<output_basename>*.log` (`Log`); probe
  via `find_on_path(&["slim"])` (conda-forge / source / Homebrew
  all install under the canonical lowercase `slim` name); GPL-3.0
  licensed; version range `4.0.0..5.0.0` (the modern release line
  is the 4.x series from 2022+; 4.0 introduced the streamlined
  Eidos type system and the `treeSeqOutput()` helpers we rely on
  for tskit interop; upper bound 5.0 reserves room for an eventual
  major bump); `bio.slim.simulate` ribbon capability), msprime
  (Python-script subprocess sister to Phase 17 Biopython;
  `case.toml` knobs `script` (path to user-authored Python script;
  required), `python` (interpreter name; default `"python3"`),
  `population_size` (`u32`, ≥ 1), `num_samples` (`u32`, ≥ 1),
  `recombination_rate` (`f64`, ≥ 0.0 and finite — per-site per-
  generation rate), `mutation_rate` (`f64`, ≥ 0.0 and finite),
  `output_basename` (filename stem; required, non-empty);
  `prepare()` stages the script into the workdir under its
  original filename and writes a flat `valenx_params.json`
  containing `population_size`, `num_samples`, `recombination_rate`
  (emitted via `{:e}` so Python's `json.load` parses it back as a
  float), `mutation_rate` (same), and `output_basename`; collects
  `<output_basename>.trees` (`Native`, "msprime tree sequence"),
  `<output_basename>.vcf` (`Tabular`, "msprime VCF"), and
  `<output_basename>.csv` (`Tabular`, "msprime per-sample
  summary"); probe via Python on PATH with an `import msprime`
  check (returns `ok = true` with a warning when import fails so
  non-standard installs aren't blocked); GPL-3.0 licensed; version
  range `1.3.0..2.0.0` (the modern `sim_ancestry()` /
  `sim_mutations()` split landed in 1.3 in 2024, paired with the
  tskit 0.5+ tree-sequence format we surface in collect(); upper
  bound 2.0 reserves room for an eventual major bump);
  `bio.msprime.simulate` ribbon capability), and tskit (Python-
  script subprocess sister to msprime; `case.toml` knobs `script`
  (path to user-authored Python script; required), `python`
  (interpreter name; default `"python3"`), `input_trees` (`.trees`
  file from SLiM or msprime; required), `output_basename`
  (filename stem; required, non-empty); `prepare()` stages script
  + tree-sequence file into the workdir under their original
  filenames so the script can resolve them via relative paths,
  then writes a flat `valenx_params.json` containing
  `input_trees` (staged filename) and `output_basename`; collects
  `<output_basename>*.csv` / `<output_basename>*.tsv` (`Tabular`,
  "tskit statistics") and `*.png` (`Native`, "tskit plot"); probe
  via Python on PATH with an `import tskit` check (same
  `ok = true` + warning fallback as msprime); MIT licensed;
  version range `0.5.0..1.0.0` (tskit 0.5+ ships the modern
  `Statistics` API surface and the v3 tree-sequence file format
  msprime 1.3+ writes; upper bound 1.0 reserves room for the
  long-promised 1.0 release); `bio.tskit.analyze` ribbon
  capability). Each wired into `valenx-app::init_registry`.
- No new canonical types, no new format readers, no new CLIs. All
  three adapters consume user-supplied inputs (SLiM `.slim` Eidos
  scripts; msprime + tskit Python scripts; tskit input `.trees`
  tree-sequence files emitted by SLiM or msprime) and emit user-
  readable artifacts (`.trees` tree sequences, `.vcf` genotype
  calls, `.csv` / `.tsv` statistics tables, `.png` rendered plots,
  `.log` run logs) that the unchanged `Results.artifacts`
  collection model surfaces directly. A first-class population-
  genetics canonical type — a typed tree-sequence representation
  spanning all three back-ends, with parsed per-tree edge / node /
  mutation tables and a typed statistics-table representation —
  defers to a future phase along with tree-sequence visualizers
  and per-population allele-frequency-spectrum viewers.
- 3 new `valenx-init` templates ship: `slim` (`slim-simulate`),
  `msprime` (`msprime-simulate`), and `tskit` (`tskit-analyze`).
  Cross-binary roundtrip test sweeps all 85 templates clean.

The full per-phase shape lives in
[docs/src/phases/phase-29-population-genetics.md](./docs/src/phases/phase-29-population-genetics.md);
the implementation plan lives at
[docs/superpowers/plans/2026-04-30-population-genetics.md](./docs/superpowers/plans/2026-04-30-population-genetics.md).

**Deliverable:** Adapter inventory at 89 of 90 fully live after this
phase (3 new population-genetics adapters added on top of the Phase
17 + Phase 17.5 + Phase 18 + Phase 18.5 + Phase 18.6 + Phase 19 +
Phase 19.5 + Phase 20 + Phase 22 + Phase 23 + Phase 24 + Phase 25 +
Phase 27 + Phase 27.5 + Phase 27.6 + Phase 28 + Phase 30 + Phase 31
+ Phase 32 + Phase 34 + Phase 35 + Phase 36 + Phase 38 totals);
Phase 29 opens the first population-genetics / evolutionary-
simulation domain to ship in Valenx, broadening the design-guides →
predict-off-targets → search-off-targets → simulate-reads → align →
quantify → call-variants → predict-structure → fold-RNA → infer-tree
→ simulate-pathway → reconstruct-3D → validate loop to three
population-genetics tools (SLiM, msprime, tskit) feeding into the
existing Phase 31 read simulators, the Phase 18 / 18.5 / 18.6
alignment surface, the Phase 19 / 20 variant-calling + transcript-
quantification stack, the Phase 17 / 17.5 prediction stack, the
Phase 28 RNA-structure tools, the Phase 30 phylogenetic-tree
builders, the Phase 32 systems-biology surface, the Phase 35
CRISPR-design tools, the Phase 36 cryo-EM reconstruction tools, and
the Phase 38 Rosetta-family adapters.

### Phase 30 — Phylogenetics · live
Open the molecular phylogenetics domain in Valenx with the three most-
used maximum-likelihood tree inference tools: IQ-TREE (Bui Quang Minh
& Robert Lanfear's de-facto modern ML tree builder — ModelFinder +
UFBoot ultrafast bootstrap), RAxML-NG (Alexey Kozlov's next-generation
RAxML rewrite — successor to classical `raxmlHPC`), and FastTree
(Morgan Price's approximate-ML inference, optimized for very large
trees — sub-quadratic in alignment size). All three follow the
established Phase 18 BWA single-binary CLI pattern: input alignment
in, tree out. No two-stage index step (the alignment is the input;
the tree is the output). Phase 30 sits numerically after Phase 27.5
and ships chronologically right after the Phase 27.5 protein-design
expansion beachhead.

- 3 new adapter crates land under `crates/valenx-adapters/bio/`:
  IQ-TREE (single-binary subprocess; `case.toml` knobs `alignment`
  (FASTA / PHYLIP / NEXUS / CLUSTAL; required), `model` default
  `"MFP"` (`"TEST"` / `"MFP"` trigger ModelFinder's automatic model
  selection; otherwise pass e.g. `"GTR+G"` / `"WAG+I+G"` verbatim),
  `bootstrap` UFBoot ultrafast bootstrap replicates default 1000
  (`0` disables), `threads` default `"AUTO"` validated against
  `^(AUTO|\d+)$`, `prefix` required, `extra_args`; `prepare()`
  builds `iqtree2 -s <alignment> -m <model> -B <bootstrap> -T
  <threads> --prefix <prefix> [extras...]` (omitting `-B` when
  `bootstrap == 0`); collects `<prefix>.treefile` (`Native`,
  "IQ-TREE ML tree") + `.iqtree` / `.log` (`Log`); probe via
  `find_on_path(&["iqtree2", "iqtree"])`; GPL-2.0 licensed;
  `bio.iqtree.tree` ribbon capability), RAxML-NG (single-binary
  subprocess with mode dispatch; `case.toml` knobs `alignment`,
  `model`, `mode` ∈ `{search, all, bootstrap}`, `bootstrap` (≥ 1
  required when `mode ∈ {all, bootstrap}`), `threads` ≥ 1 default
  1, `prefix`, `extra_args`; `prepare()` builds `raxml-ng --<mode>
  --msa <alignment> --model <model> --threads <N> --prefix
  <prefix> [--bs-trees <bootstrap> if mode in {all, bootstrap}]
  [extras...]`; collects `<prefix>.raxml.bestTree` (`Native`,
  "RAxML-NG ML tree") + `<prefix>.raxml.support` (`Native`) +
  `<prefix>.raxml.log` (`Log`); probe via `find_on_path(&["raxml-
  ng"])`; AGPL-3.0 licensed; `bio.raxml-ng.tree` ribbon
  capability), and FastTree (single-binary subprocess; writes
  Newick to stdout — captured to `output` via the MAFFT-style
  stdout-redirect pattern; `case.toml` knobs `alignment`, `output`
  (Newick path), `seq_type` ∈ `{nt, aa}`, `use_gtr` default `true`
  (uses GTR for nucleotides — FastTree's default is JC without
  this flag; ignored for amino acid), `gamma` default `false`
  (gamma rate-variation toggle), `extra_args`; `prepare()` builds
  — nucleotide: `FastTree [-nt] [-gtr if use_gtr] [-gamma if
  gamma] <alignment>` → stdout; amino-acid: `FastTree [-gamma if
  gamma] <alignment>` → stdout; collects `output` as `Native`
  artifact "FastTree Newick tree"; probe via
  `find_on_path(&["FastTree", "fasttree"])` (binary name varies by
  distro); GPL-2.0 licensed; `bio.fasttree.tree` ribbon
  capability). Each wired into `valenx-app::init_registry`.
- No new canonical types, no new format readers, no new CLIs. All
  three adapters consume the existing Phase 18 / 20 multiple-
  sequence alignment inputs (FASTA / PHYLIP / NEXUS / CLUSTAL) and
  emit Newick-format tree files that the unchanged
  `Results.artifacts` collection model surfaces directly through
  the `Native` artifact kind. A first-class `Tree` canonical type
  with a Newick reader as a Valenx CLI defers to a future phase
  along with visualization integrations.
- 3 new `valenx-init` templates ship: `iqtree` with alias
  `iqtree2` (`iqtree-build`), `raxml-ng` with alias `raxml`
  (`raxml-ng-build`), and `fasttree` (`fasttree-build`). Cross-
  binary roundtrip test sweeps all 60 templates clean.

The full per-phase shape lives in
[docs/src/phases/phase-30-phylogenetics.md](./docs/src/phases/phase-30-phylogenetics.md);
the implementation plan lives at
[docs/superpowers/plans/2026-04-30-phylogenetics.md](./docs/superpowers/plans/2026-04-30-phylogenetics.md).

**Deliverable:** Adapter inventory at 64 of 65 fully live after this
phase (3 new phylogenetics adapters added on top of the Phase 17 +
Phase 17.5 + Phase 18 + Phase 18.5 + Phase 18.6 + Phase 19 + Phase
19.5 + Phase 20 + Phase 22 + Phase 23 + Phase 24 + Phase 27 + Phase
27.5 + Phase 34 totals); phylogenetics workflows land alongside the
existing biology adapters in the same case-toml / prepare / run /
collect shell, broadening the align → quantify → predict → infer-
tree → validate loop to eleven alignment / search tools (BWA,
Bowtie2, HISAT2, STAR, minimap2, MAFFT, MUSCLE, HMMER, samtools,
MMseqs2, DIAMOND), two transcript quantifiers (Salmon, Kallisto),
five prediction tools (ESMFold, OpenFold, AlphaFold 2, AlphaFold 3,
ColabFold), and three phylogenetic-tree builders (IQ-TREE, RAxML-NG,
FastTree).

### Phase 30.5 — Bayesian phylogenetics · live
Sister-domain expansion of Phase 30. Round out the molecular
phylogenetics surface with the two de-facto Bayesian phylogenetic
inference engines: BEAST 2 (the cross-platform Bayesian Evolutionary
Analysis by Sampling Trees v2 engine — canonical Bayesian MCMC
framework for time-calibrated phylogenetics: tip-dated trees,
relaxed molecular clocks, coalescent demographic models, birth-death
speciation models, and the ever-growing universe of BEAST 2
packages (BDSKY, MASCOT, BEASTling, StarBEAST3, ...); LGPL-2.1;
single-binary subprocess shape sister to Phase 18 BWA), and MrBayes
(the long-standing Bayesian MCMC phylogenetic inference engine —
historic workhorse for Bayesian phylogenetics; alongside BEAST 2
the de-facto choice for posterior tree sampling across nucleotide
/ amino-acid / morphological datasets, with its own NEXUS-embedded
model-and-mcmc command language and built-in Metropolis-coupled
MCMC ("MC^3") chain swapping; GPL-3.0; single-binary subprocess
shape sister to BEAST 2). Both adapters share the established Phase
18 BWA single-binary CLI pattern: a user-authored model description
(BEAST 2 XML or MrBayes NEXUS file) in, posterior tree + parameter
samples out — the same shape the Phase 30 ML tools use, with the
inputs swapped from multiple-sequence alignments to MCMC model
files. Phase 30.5 sits numerically after Phase 30 and ships
chronologically right after Phase 38 Rosetta — same chronological-
vs-numerical convention used for Phase 17.5 / 24 / 28 / 31 / 35.

- 2 new adapter crates land under `crates/valenx-adapters/bio/`:
  BEAST 2 (single-binary subprocess sister to Phase 18 BWA;
  `case.toml` knobs `xml` (BEAUti-generated XML model file;
  required), `seed` (optional `u64`; passed via `beast -seed <N>`
  when present, otherwise BEAST picks its own seed and prints it
  on the run banner), `threads` (`u32`, ≥ 1, default 1; maps to
  `-threads N` for tree-likelihood evaluation parallelism),
  `overwrite` (default `false`; toggles `-overwrite` so an
  existing output set from a previous run is replaced rather than
  triggering a fail), `extra_args`; `prepare()` resolves the XML
  against the case directory when relative, validates it exists
  on disk (returns `InvalidCase` with a helpful message when
  missing), and composes `beast [-seed <N>] -threads <N>
  [-overwrite] <xml> [extras...]` with the XML positional last so
  BEAST treats it as the model file rather than the value of an
  earlier flag; `run()` streams BEAST's `Random number seed` /
  `BEAST v2` startup banner / periodic `Sample` / `posterior`
  chain-status lines / `End likelihood` / `Total calculation
  time` end-of-run sentinels into progress hints; collects
  `*.log` (`Log`, "BEAST 2 trace log") and `*.trees` (`Native`,
  "BEAST 2 sampled trees"); the adapter doesn't try to predict
  the exact filenames since BEAST writes whatever the XML's
  `<log fileName="...">` sites configure; probe via
  `find_on_path(&["beast"])` (the generic version detector tries
  the conventional `--version` and BEAST's own `-version` form);
  LGPL-2.1 licensed; version range `2.7.0..3.0.0` (the modern
  stable line is the 2.7.x series from 2022+ that introduced
  modern threading + the package manager); `bio.beast2.mcmc`
  ribbon capability), and MrBayes (single-binary subprocess
  sister to BEAST 2; `case.toml` knobs `nexus` (NEXUS data file
  with embedded MRBAYES block driving the run; required), `batch`
  (default `false`; toggles `-i` so MrBayes runs the embedded
  commands non-interactively and exits cleanly rather than
  waiting on stdin at the prompt — the right default for non-
  interactive automation), `extra_args`; `prepare()` resolves the
  NEXUS path against the case directory when relative, validates
  it exists on disk, and composes `mb [-i if batch] <nexus>
  [extras...]` with the NEXUS positional last so MrBayes treats
  it as the model file rather than the value of an earlier flag;
  `run()` streams MrBayes's `MrBayes v` / `Initializing` startup
  banner / periodic `Generation NNNN` / `Avg standard deviation
  of split frequencies` chain-status lines / `Analysis completed`
  / `Continue with analysis` end-of-run sentinels into progress
  hints; collects `*.t` (`Native`, "MrBayes tree samples"), `*.p`
  (`Tabular`, "MrBayes parameter samples"), and `*.con.tre`
  (`Native`, "MrBayes consensus tree"); probe via
  `find_on_path(&["mb"])` (the binary is literally named `mb` —
  the project's own convention); GPL-3.0 licensed; version range
  `3.2.0..4.0.0` (the long-running stable 3.2.x line that every
  distro ships covers every release through 3.2.7);
  `bio.mrbayes.mcmc` ribbon capability). Each wired into
  `valenx-app::init_registry`.
- No new canonical types, no new format readers, no new CLIs.
  Both adapters consume user-supplied inputs (BEAST 2 XML model
  files, MrBayes NEXUS data files with embedded MRBAYES blocks)
  and emit user-readable artifacts (BEAST 2 `.log` trace files +
  `.trees` sampled-tree posteriors, MrBayes `.t` tree samples +
  `.p` parameter samples + `.con.tre` consensus trees) that the
  unchanged `Results.artifacts` collection model surfaces
  directly. A first-class Bayesian-phylogenetics canonical type —
  a typed posterior representation spanning both back-ends, with
  parsed per-generation parameter traces and per-sample tree
  topologies plus convergence-diagnostic helpers (effective
  sample size, Gelman-Rubin) — defers to a future phase along
  with trace visualizers, tree-density plots, and consensus-tree
  viewers.
- 2 new `valenx-init` templates ship: `beast2` with alias `beast`
  (`beast2-mcmc`) and `mrbayes` with alias `mb` (`mrbayes-mcmc`).
  Cross-binary roundtrip test sweeps all 90 templates clean.

The full per-phase shape lives in
[docs/src/phases/phase-30-5-bayesian-phylogenetics.md](./docs/src/phases/phase-30-5-bayesian-phylogenetics.md);
the implementation plan lives at
[docs/superpowers/plans/2026-04-30-bayesian-phylogenetics.md](./docs/superpowers/plans/2026-04-30-bayesian-phylogenetics.md).

**Deliverable:** Adapter inventory at 94 of 95 fully live after this
phase (2 new Bayesian-phylogenetics adapters added alongside the
Phase 39 DNA-structural-geometry trio on top of the Phase 17 +
Phase 17.5 + Phase 18 + Phase 18.5 + Phase 18.6 + Phase 19 + Phase
19.5 + Phase 20 + Phase 22 + Phase 23 + Phase 24 + Phase 25 + Phase
27 + Phase 27.5 + Phase 27.6 + Phase 28 + Phase 29 + Phase 30 +
Phase 31 + Phase 32 + Phase 34 + Phase 35 + Phase 36 + Phase 38
totals); Phase 30.5 rounds out the molecular phylogenetics surface
that Phase 30 opened from the maximum-likelihood side, broadening
the align → quantify → predict → infer-tree-ML → infer-tree-
Bayesian → validate loop to two Bayesian MCMC tree-inference
engines (BEAST 2, MrBayes) feeding into the existing Phase 30
maximum-likelihood phylogenetic-tree builders (IQ-TREE, RAxML-NG,
FastTree), the Phase 18 / 18.5 / 18.6 alignment surface, the Phase
20 transcript-quantification stack, the Phase 17 / 17.5 prediction
stack, the Phase 28 RNA-structure tools, and the Phase 29
population-genetics trio.

### Phase 31 — Sequencing read simulators · live
Open the **first sequencing read-simulation domain** in Valenx with
three established open-source tools that span all three major
sequencing-technology classes: ART (Weichun Huang's NIEHS Illumina-
platform read simulator — de-facto choice for synthesising FASTQs
that match per-platform empirical error profiles for HiSeq 2500 /
HiSeq X / MiSeq v3 / NextSeq 500 / MiniSeq; GPL-3.0; single-binary
`art_illumina`), wgsim (Heng Li's classic Whole-Genome SIMulator
that ships alongside samtools — always paired-end, always position-
uniform, deliberately simple under a uniform sequencing-error model;
the canonical "small + classic" simulator for fast smoke-testing of
mappers and variant callers when realistic error spectra are not
required; MIT; single-binary `wgsim` with positional output
arguments after the reference), and Badread (Ryan Wick's long-read
simulator with realistic Nanopore + PacBio CLR error profiles —
random / chimeric / adapter / glitch read types, junk-read injection,
identity drift, length distributions calibrated against actual
sequencer output; GPL-3.0; single-binary `badread simulate` that
writes its FASTQ to stdout via the MAFFT-style stdout-redirect-to-
file pattern). All three follow the established Phase 18 BWA single-
binary CLI pattern: reference FASTA in, simulated FASTQ(s) out.
Phase 31 sits numerically before Phase 32 but ships chronologically
right after Phase 36 cryo-EM — same chronological-vs-numerical
convention used for Phase 17.5 / 24 / 28.

- 3 new adapter crates land under `crates/valenx-adapters/bio/`:
  ART (single-binary subprocess wrapping `art_illumina`; `case.toml`
  knobs `reference` (FASTA; required), `output_prefix` (filename
  stem; ART writes `<prefix>.fq` single-end or `<prefix>1.fq` +
  `<prefix>2.fq` paired-end; required, non-empty), `sequencing_system`
  ∈ `{HS25, HSXt, MSv3, NS50, MinS}` (HiSeq 2500 / HiSeq X TruSeq /
  MiSeq v3 / NextSeq 500 / MiniSeq; required), `read_length`
  (≥ 1; required), `fold_coverage` (> 0.0; required), `paired_end`
  (default `false`), `fragment_mean` (mean insert size for paired-
  end; default 200.0, > 0.0 when `paired_end`), `fragment_sd`
  (insert-size stddev for paired-end; default 10.0, > 0.0 when
  `paired_end`), `extra_args`; `prepare()` builds `art_illumina -ss
  <sequencing_system> -i <reference> -l <read_length> -f
  <fold_coverage> -o <output_prefix> [-p -m <fragment_mean> -s
  <fragment_sd> if paired_end] [extras...]`; collects
  `<output_prefix>*.fq` (`Tabular`, "ART simulated reads") and
  `<output_prefix>*.aln` (`Log`, "ART alignment record" — the per-
  read alignment record ART writes alongside, useful for validating
  aligner accuracy against the simulated truth); probe via
  `find_on_path(&["art_illumina"])`; GPL-3.0 licensed; version range
  `2.5.0..3.0.0` (the long-running ChocolateCherryCake `2.5.x`
  series since 2016; Bioconda + Homebrew ship 2.5.8); the
  `valenx-init` template ships with the alias `art-illumina`
  alongside the canonical `art`; `bio.art.simulate` ribbon
  capability), wgsim (single-binary subprocess; `case.toml` knobs
  `reference` (FASTA; required), `output1` (FASTQ for read 1;
  required, non-empty), `output2` (FASTQ for read 2; required, non-
  empty — wgsim is paired-end only), `num_pairs` (≥ 1; required),
  `length1` (read 1 length, default 70, ≥ 1), `length2` (read 2
  length, default 70, ≥ 1), `fragment_size` (outer fragment length,
  default 500, > 0), `error_rate` (per-base error rate in
  `0.0..=1.0`, default 0.02 — typical Illumina baseline),
  `extra_args`; `prepare()` builds `wgsim -N <num_pairs> -1
  <length1> -2 <length2> -d <fragment_size> -e <error_rate>
  <reference> <output1> <output2> [extras...]`; collects `output1`
  and `output2` as `Tabular` artifacts ("wgsim simulated reads");
  probe via `find_on_path(&["wgsim"])`; MIT licensed; version range
  `1.0.0..2.0.0` (wgsim is versioned alongside the parent samtools
  1.x line); `bio.wgsim.simulate` ribbon capability), and Badread
  (single-binary subprocess with stdout-redirect — Badread writes
  its simulated FASTQ to stdout, captured to `output` via the
  MAFFT-style stdout-redirect-to-file pattern (spawn the child
  directly, attach stdout to a `File` via `Stdio::from(file)`,
  stream stderr through the line handler); `case.toml` knobs
  `reference` (FASTA; required), `output` (FASTQ output path;
  required, non-empty), `quantity` (Badread `--quantity` literal —
  one or more decimal digits with optional `K` / `M` / `G` / `T`
  SI suffix, e.g. `"100M"` for 100 megabases or `"5G"` for 5
  gigabases; validated via the `is_valid_quantity` helper),
  `error_model` ∈ `{nanopore2018, nanopore2020, nanopore2023,
  pacbio2016}` (per-platform error profile baked into the Badread
  distribution), `identity_mean` (read identity mean as a percentage
  in `0.0..=100.0`, default 87.5), `length_mean` (read length mean
  in bases, default 15000.0, > 0.0), `length_sd` (read length stddev
  in bases, default 13000.0, > 0.0), `extra_args`; `prepare()`
  builds `badread simulate --reference <reference> --quantity
  <quantity> --error_model <error_model> --identity <identity_mean>
  --length <length_mean>,<length_sd> [extras...]` → stdout, captured
  to `output` via the MAFFT-style stdout-redirect pattern; collects
  `output` as a single `Tabular` artifact ("Badread simulated
  reads"); probe via `find_on_path(&["badread"])`; GPL-3.0 licensed;
  version range `0.4.0..1.0.0` (the long-running 0.4.x stable
  series; a 1.0 cut hasn't happened yet but the upper bound reserves
  room for it); `bio.badread.simulate` ribbon capability). Each
  wired into `valenx-app::init_registry`.
- No new canonical types, no new format readers, no new CLIs. All
  three adapters consume user-supplied reference FASTAs (the existing
  `valenx_bio::format::fasta` reader already inspects sequence count
  + identifiers + alphabets) and emit FASTQ files that the existing
  Phase 18 `valenx-fastq` CLI inspects for record count, base
  quality distributions, and read-length statistics. The unchanged
  `Results.artifacts` collection model surfaces every emitted FASTQ
  + ART alignment record directly. A first-class read-simulation
  provenance type — recording which simulator produced which FASTQ
  under which error model — defers to a future phase along with
  simulator-aware pipeline stitching.
- 3 new `valenx-init` templates ship: `art` with alias `art-illumina`
  (`art-simulate`), `wgsim` (`wgsim-simulate`), and `badread`
  (`badread-simulate`). Cross-binary roundtrip test sweeps all 77
  templates clean.

The full per-phase shape lives in
[docs/src/phases/phase-31-read-simulators.md](./docs/src/phases/phase-31-read-simulators.md);
the implementation plan lives at
[docs/superpowers/plans/2026-04-30-read-simulators.md](./docs/superpowers/plans/2026-04-30-read-simulators.md).

**Deliverable:** Adapter inventory at 81 of 82 fully live after this
phase (3 new read-simulator adapters added on top of the Phase 17 +
Phase 17.5 + Phase 18 + Phase 18.5 + Phase 18.6 + Phase 19 + Phase
19.5 + Phase 20 + Phase 22 + Phase 23 + Phase 24 + Phase 25 + Phase
27 + Phase 27.5 + Phase 27.6 + Phase 28 + Phase 30 + Phase 32 +
Phase 34 + Phase 36 totals); Phase 31 opens the first sequencing
read-simulation domain to ship in Valenx, broadening the simulate-
reads → align → quantify → call-variants → predict-structure →
fold-RNA → infer-tree → simulate-pathway → reconstruct-3D → validate
loop to three read simulators (ART, wgsim, Badread) feeding into the
existing Phase 18 / 18.5 / 18.6 alignment surface, the Phase 19 / 20
variant-calling + transcript-quantification stack, the Phase 17 /
17.5 prediction stack, the Phase 28 RNA-structure tools, the Phase
30 phylogenetic-tree builders, the Phase 32 systems-biology surface,
and the Phase 36 cryo-EM reconstruction tools.

### Phase 32 — Systems biology · live
Open the **first systems-biology / multiscale modeling domain** in
Valenx with three established open-source tools that span the systems-
biology tradeoff space: COPASI (the COmplex PAthway SImulator — de-
facto biochemical pathway / ODE-based systems-biology suite descended
from the Gepasi lineage; Artistic-2.0; reads `.cps` native or SBML
`.xml`), BioNetGen (rule-based modeling language + tool suite for
combinatorially-complex signaling networks; MIT; Perl driver `BNG2.pl`
reads `.bngl` rule-based models and emits `<basename>.net` /
`<basename>.gdat` / `<basename>.cdat` outputs), and PhysiCell (Paul
Macklin's agent-based, off-lattice multicellular simulator — tens to
hundreds of thousands of individual cells coupled to a reaction-
diffusion microenvironment for substrates like oxygen and drugs;
canonical use case is tumour growth + immunology; BSD-3-Clause; models
compile per-project to a project-specific C++ binary, so the adapter
takes both the user's compiled `binary` path and the run-time XML
configuration). All three follow the established Phase 18 BWA single-
binary CLI pattern: model file in, results out. Phase 32 sits
numerically after Phase 30 and ships chronologically right after
Phase 25 quantum chemistry — same chronological-vs-numerical
convention used for Phase 17.5 / 24 / 28.

- 3 new adapter crates land under `crates/valenx-adapters/bio/`:
  COPASI (single-binary subprocess; `case.toml` knobs `model` (`.cps`
  or SBML `.xml`; required), `report` (optional `--save <report>`
  target so collect() finds the run output deterministically without
  walking), `run_all` (default false; toggles `--scheduled` to execute
  every defined task in the file rather than just the primary one),
  `extra_args`; `prepare()` builds `CopasiSE [--save <report>]
  <model> [--scheduled if run_all] [extras...]`; collects the
  explicit `report` path as `Tabular` ("COPASI report") when supplied,
  else walks the workdir top-level for `.csv` / `.txt` files; probe
  via `find_on_path(&["CopasiSE"])` (capital `C-S-E` "Self-Executing"
  CLI is the canonical headless name); Artistic-2.0 licensed; version
  range `4.40.0..5.0.0` (4.x is the long-running stable line; 4.40 is
  a recent floor that ships SBML L3v2 + the task scheduler);
  `bio.copasi.simulate` ribbon capability), BioNetGen (single-binary
  Perl-driver subprocess; `case.toml` knobs `model` (`.bngl`;
  required), `output_basename` (required — becomes the `-o` prefix
  every output file inherits so collect() walks deterministically),
  `generate_only` (default false; adds `--no-execute` to skip simulate
  / scan / fitting actions and emit just the expanded reaction
  network), `extra_args`; `prepare()` builds `BNG2.pl [--no-execute
  if generate_only] -o <output_basename> <model> [extras...]`;
  collects `<output_basename>*.net` (`Native`, "BioNetGen reaction
  network"), `<output_basename>*.gdat` (`Tabular`, "BioNetGen species
  trajectories"), `<output_basename>*.cdat` (`Tabular`, "BioNetGen
  concentrations") — `parameter_scan` per-trial variants share the
  basename prefix (e.g. `<basename>_001.gdat`) so the prefix-
  restricted walk picks them up too; probe via
  `find_on_path(&["BNG2.pl"])`; MIT licensed; version range
  `2.8.0..3.0.0`; `bio.bionetgen.simulate` ribbon capability), and
  PhysiCell (per-project-binary subprocess; `case.toml` knobs
  `binary` (the per-project compiled executable; required), `config`
  (the `.xml` settings file PhysiCell binaries accept as a positional
  argument; required), `extra_args`; `prepare()` validates `binary`
  and `config` exist on disk and returns a helpful "PhysiCell models
  compile per-project — clone the framework, edit the project's
  `custom_modules/` source, run `make`, and point this field at the
  resulting executable." `InvalidCase` if the binary is missing,
  rather than letting the child fail later with a less obvious
  "command not found"; builds `<binary> <config> [extras...]`;
  collects `output/*.xml` and `output/*.mat` (`Native`, "PhysiCell
  tissue snapshot") plus `output/*.csv` (`Tabular`, "PhysiCell scalar
  table"); probe via `find_on_path(&["physicell"])` returns
  `ok = true` either way (most installs won't have a generic
  `physicell` binary on PATH — the per-project build pattern means
  there isn't a canonical one) and attaches a warning that the real
  validation happens at prepare time; BSD-3-Clause licensed; version
  range `1.13.0..2.0.0`; `bio.physicell.simulate` ribbon capability).
  Each wired into `valenx-app::init_registry`.
- No new canonical types, no new format readers, no new CLIs. All
  three adapters consume user-supplied inputs (COPASI `.cps` / SBML
  `.xml` archives, BioNetGen `.bngl` rule-based models, PhysiCell
  per-project compiled binaries + XML config) and emit user-readable
  artifacts (CSV / TXT tabular reports, reaction-network `.net`
  files, species-trajectory `.gdat` / concentration `.cdat` tabular
  files, `.xml` / `.mat` per-snapshot tissue state, per-cell scalar
  `.csv` summaries) that the unchanged `Results.artifacts`
  collection model surfaces directly. A first-class systems-biology
  canonical type — a generic SBML / BNGL / per-cell state type
  spanning all three back-ends — defers to a future phase along with
  SBML readers and tissue-snapshot visualizers.
- 3 new `valenx-init` templates ship: `copasi` (`copasi-simulate`),
  `bionetgen` with alias `bng` (`bionetgen-simulate`), and
  `physicell` (`physicell-simulate`). Cross-binary roundtrip test
  sweeps all 71 templates clean.

The full per-phase shape lives in
[docs/src/phases/phase-32-systems-biology.md](./docs/src/phases/phase-32-systems-biology.md);
the implementation plan lives at
[docs/superpowers/plans/2026-04-30-systems-biology.md](./docs/superpowers/plans/2026-04-30-systems-biology.md).

**Deliverable:** Adapter inventory at 75 of 76 fully live after this
phase (3 new systems-biology adapters added on top of the Phase 17 +
Phase 17.5 + Phase 18 + Phase 18.5 + Phase 18.6 + Phase 19 + Phase
19.5 + Phase 20 + Phase 22 + Phase 23 + Phase 24 + Phase 25 +
Phase 27 + Phase 27.5 + Phase 27.6 + Phase 28 + Phase 30 + Phase 34
totals); Phase 32 opens the first systems-biology / multiscale
modeling domain to ship in Valenx, broadening the build-network →
simulate-pathway → expand-rules → grow-tissue → predict-structure →
fold-RNA → infer-tree → validate loop to three systems-biology tools
(COPASI, BioNetGen, PhysiCell) feeding into the existing Phase 17 /
17.5 prediction stack, the Phase 24 / 25 cheminformatics + quantum-
chemistry surface, the Phase 28 RNA-structure tools, and the Phase 30
phylogenetic-tree builders.

### Phase 32.5 — Spatial stochastic · live ✅
Sister-adapter expansion of the Phase 32 systems-biology trio
(COPASI / BioNetGen / PhysiCell). Round out the systems-biology /
multiscale modeling surface with the two canonical **spatial
stochastic / cell-scale reaction-diffusion** simulators that
Phase 32 explicitly deferred: Smoldyn (Steve Andrews's spatial
stochastic reaction-diffusion simulator, LGPL-2.1 — resolves
individual molecules as particles diffusing and reacting in
continuous 3D space, the canonical choice when the question is
"where does each molecule actually end up over time" rather
than "what is the well-mixed concentration vs. t" Phase 32
COPASI's ODE / SSA covers; single-binary subprocess shape
sister to Phase 18 BWA), MCell (Salk Institute / Stiles, Bartol's
cell-scale Monte Carlo spatial stochastic simulator, GPL-2.0 —
walks the user's `.mdl` (Model Description Language) model and
runs Brownian-dynamics particle trajectories with Monte Carlo
reaction sampling on intricate triangle-mesh geometry; canonical
use case is sub-cellular signaling (synaptic transmission,
calcium dynamics, receptor binding); single-binary subprocess
shape sister to Smoldyn). Both adapters follow the established
Phase 18 BWA single-binary CLI pattern: model file in, reaction-
data + trajectory artifacts out. Phase 32.5 sits numerically
adjacent to Phase 32 and ships chronologically right after Phase
17.7 structure tools — same chronological-vs-numerical
convention used for Phase 17.5 / 24 / 28 / 31 / 35 / 39 / 5.5 /
5.6 / 5.7.

- 2 new adapter crates land under `crates/valenx-adapters/bio/`:
  Smoldyn (single-binary subprocess sister to Phase 18 BWA;
  `case.toml` knobs `config` (Smoldyn `.txt` configuration
  describing simulation geometry — boundaries, surfaces,
  compartments — plus molecule species + diffusion coefficients
  and per-pair / per-surface reactions; required), `extra_args`;
  `prepare()` resolves `config` against the case directory when
  relative, validates it exists on disk, builds `smoldyn
  <config> [extras...]`; `collect()` walks for `*.txt` (`Tabular`,
  "Smoldyn output table"), `*.dat` (`Tabular`, "Smoldyn data"),
  `*.log` (`Log`, "Smoldyn log"); probe via
  `find_on_path(&["smoldyn"])`; LGPL-2.1 licensed; version range
  `2.70.0..3.0.0`; `bio.smoldyn.simulate` ribbon capability),
  and MCell (single-binary subprocess sister to Smoldyn;
  `case.toml` knobs `mdl` (`.mdl` MCell model description file;
  required), `seed` (`Option<u32>` — when `Some(n)` the adapter
  emits `-seed` and `<n>` as TWO separate args; when `None`
  MCell picks its own seed and prints it on the run banner —
  same shape as the Phase 29 SLiM `-s` and Phase 30.5 BEAST 2
  `-seed` knobs), `extra_args`; `prepare()` resolves `mdl`
  against the case directory when relative, validates it exists
  on disk, threads `-seed` + `<n>` as two separate OsStrings
  into the argv only when `seed` is `Some(_)`, builds `mcell
  [-seed <N>] <mdl> [extras...]`; `collect()` walks for `*.dat`
  (`Tabular`, "MCell reaction data"), `*.dx` (`Native`, "MCell
  visualization data"), `*.log`; probe via
  `find_on_path(&["mcell"])`; GPL-2.0 licensed; version range
  `4.0.0..5.0.0`; `bio.mcell.simulate` ribbon capability). Each
  wired into `valenx-app::init_registry`.
- No new canonical types, no new format readers, no new CLIs.
- 2 new `valenx-init` templates ship: `smoldyn`, `mcell`. Cross-
  binary roundtrip test sweeps all 115 templates clean.

The full per-phase shape lives in
[docs/src/phases/phase-32-5-spatial-stochastic.md](./docs/src/phases/phase-32-5-spatial-stochastic.md);
the implementation plan lives at
[docs/superpowers/plans/2026-05-03-spatial-stochastic.md](./docs/superpowers/plans/2026-05-03-spatial-stochastic.md).

**Deliverable:** Adapter inventory at 119 of 120 fully live after
this phase trio (2 new spatial-stochastic adapters added on top
of the prior set, alongside the Phase 40 microscopy trio and
Phase 41 sequence-editors pair); Phase 32.5 sister-expands Phase
32 systems biology with cell-scale spatial stochastic simulators,
broadening the simulate-pathway → expand-rules → grow-tissue
loop to also cover diffuse-particles → trace-MCell-trajectories
on intricate sub-cellular geometry.

### Phase 33 — Synthetic biology · live
Open the **first synthetic biology / genetic-circuit design domain**
in Valenx with three established open-source tools that span the
synthetic-biology tradeoff space: pySBOL (the Python implementation
(pySBOL3) of the Synthetic Biology Open Language standard, Apache-
2.0; SBOL captures components, sequences, interactions, constraints,
and the full provenance of a synthetic design as RDF/XML or JSON-LD
that round-trips with every SBOL-conformant tool — j5, Cello,
SynBioHub, iBioSim; Python-script subprocess shape sister to Phase
17 Biopython), j5 (JBEI's canonical DNA-assembly automation tool,
BSD-3-Clause; consumes a target circuit design (CSV row per
cassette) plus a parts library (CSV row per part / oligo), then
plans the optimal Gibson / Golden-Gate / SLIC / SLIM assembly
strategy and writes the per-step protocol + GenBank construct
files; **JAR-distributed** — single-binary subprocess shape sister
to Phase 18 BWA but the user supplies the absolute path to the JAR
via case input and we probe `java` itself rather than the JAR), and
Cello (CIDAR's canonical genetic-circuit DNA compiler — Cello v2,
BSD-3-Clause; consumes a Verilog netlist describing the desired
logic function plus a triplet of JSON constraint files (a user
constraint file pinning the chassis / library, an input sensor file
pinning the input promoters, an output device file pinning the
reporter) and emits a fully assembled DNA construct that implements
the logic in a living cell, running a simulated-annealing
optimization over the gate-assignment problem and outputting a
Graphviz `.dot` netlist + circuit diagram PNG + human-readable
report; **JAR-distributed** — single-binary subprocess shape sister
to j5). pySBOL follows the established Phase 17 Biopython Python-
script subprocess shape: the user supplies a Python script that
imports `sbol3` and reads `valenx_params.json` for the parsed
knobs. j5 + Cello are both JAR-distributed: no `j5` / `cello`
launcher binary on PATH; the user supplies the absolute path to the
JAR via case input, and we probe `java` itself rather than the JAR
— different sites pin different j5 / Cello releases under different
paths, so the j5 / Cello version implicit in the jar is the
authoritative pin. Phase 33 sits numerically between Phase 32
systems biology and Phase 34 molecular docking and ships
chronologically right after Phase 39 DNA structural geometry —
same chronological-vs-numerical convention used for Phase 17.5 /
24 / 28 / 31 / 35 / 39.

- 3 new adapter crates land under `crates/valenx-adapters/bio/`:
  pySBOL (Python-script subprocess sister to Phase 17 Biopython;
  `case.toml` knobs `script` (path to user-supplied Python script;
  required), `python` (interpreter name; default `"python3"`),
  `input_sbol` (optional starting SBOL XML document; `None` when
  the script generates the design from scratch), `output_basename`
  (filename stem the user's script uses for outputs — surfaced
  here so collect() can label artifacts uniformly; required, non-
  empty); `prepare()` stages script + optional input SBOL into the
  workdir under their original filenames so the script can resolve
  them via relative paths, then writes a flat `valenx_params.json`
  containing `input_sbol` (staged filename or literal `null`) and
  `output_basename`; collects `<output_basename>*.xml` (`Tabular`,
  "pySBOL document") and `<output_basename>*.json` (`Log`, "pySBOL
  composition log"); probe via Python on PATH with an `import
  sbol3` check (returns `ok = true` with a warning when import
  fails so non-standard installs aren't blocked); Apache-2.0
  licensed; version range `3.0.0..4.0.0` (pySBOL3 is the modern
  Python rewrite — the older 2.x line is deprecated; 3.0 is the
  floor); the init alias `sbol` resolves to the same template as
  the canonical `pysbol` name; `bio.pysbol.compose` ribbon
  capability), j5 (JAR-distributed single-binary subprocess sister
  to Phase 18 BWA; `case.toml` knobs `jar` (absolute path to
  `j5.jar`; required), `design_csv` (j5 design CSV with parts,
  oligos, target; required), `parts_csv` (parts list CSV;
  required), `output_basename` (filename stem the user expects j5
  to produce; required, non-empty), `extra_args`; `prepare()`
  resolves all three input paths against the case directory when
  relative, validates each file exists on disk (returns
  `InvalidCase` with a helpful message when missing), and composes
  `java -jar <jar> -d <design_csv> -p <parts_csv> -o
  <output_basename> [extras...]`; collects `<output_basename>*.csv`
  (`Tabular`, "j5 assembly plan") and `<output_basename>*.gb`
  (`Native`, "j5 GenBank output"); probe via
  `find_on_path(&["java"])` — j5's version comes from the jar
  itself, not from `java`, so we surface no version here (the user
  pins the j5 release implicitly by the jar they point at); BSD-
  3-Clause licensed; version range `1.0.0..2.0.0` (j5 has been on
  a 1.x line for over a decade); `bio.j5.assemble` ribbon
  capability), and Cello (JAR-distributed single-binary subprocess
  sister to j5; `case.toml` knobs `jar` (absolute path to the
  Cello jar; required), `verilog` (`.v` Verilog circuit
  description; required), `user_constraints` (`.UCF` user
  constraints file pinning the chassis / library; required),
  `input_sensors` (`.input.json` pinning the input promoters;
  required), `output_devices` (`.output.json` pinning the
  reporter; required), `output_basename` (filename stem Cello uses
  for the output directory; required, non-empty), `extra_args`;
  `prepare()` resolves all five input paths against the case
  directory when relative, validates each file exists on disk, and
  composes `java -jar <jar> -inputNetlist <verilog>
  -targetDataFile <user_constraints> -inputSensorFile
  <input_sensors> -outputDeviceFile <output_devices> -outputDir
  <output_basename> [extras...]`; collects `<output_basename>*.txt`
  (`Log`, "Cello report"), `<output_basename>*.png` (`Native`,
  "Cello circuit diagram"), and `<output_basename>*.dot` (`Native`,
  "Cello Graphviz netlist"); probe via `find_on_path(&["java"])`
  (same JAR-versioning shape as j5); BSD-3-Clause licensed;
  version range `2.0.0..3.0.0` (Cello v2 is the modern Java
  rewrite (2020+); the v1 line was Python and is deprecated);
  `bio.cello.compile` ribbon capability). Each wired into
  `valenx-app::init_registry`.
- No new canonical types, no new format readers, no new CLIs. All
  three adapters consume user-supplied inputs (pySBOL Python
  composition scripts + optional starting SBOL XML, j5 design +
  parts CSVs + the j5 jar, Cello Verilog netlist + UCF + input-
  sensor + output-device JSONs + the Cello jar) and emit user-
  readable artifacts (pySBOL XML / JSON SBOL documents, j5 CSV
  assembly plans + GenBank `.gb` constructs, Cello Graphviz `.dot`
  netlists + circuit-diagram PNGs + human-readable text reports)
  that the unchanged `Results.artifacts` collection model surfaces
  directly. A first-class synthetic-biology canonical type — a
  typed SBOL-document representation spanning pySBOL output + j5
  GenBank + Cello netlists, with parsed component / sequence /
  interaction graphs — defers to a future phase along with
  circuit-diagram visualizers and per-construct interactive
  overlays.
- 3 new `valenx-init` templates ship: `pysbol` with alias `sbol`
  (`pysbol-compose`), `j5` (`j5-assemble`), and `cello`
  (`cello-compile`). Cross-binary roundtrip test sweeps all 96
  templates clean.

The full per-phase shape lives in
[docs/src/phases/phase-33-synthetic-biology.md](./docs/src/phases/phase-33-synthetic-biology.md);
the implementation plan lives at
[docs/superpowers/plans/2026-04-30-synthetic-biology.md](./docs/superpowers/plans/2026-04-30-synthetic-biology.md).

**Deliverable:** Adapter inventory at **100 of 101** fully live
after this phase (3 new synthetic-biology adapters added alongside
the Phase 5.5 MD-analysis-expansion trio on top of the Phase 17 +
Phase 17.5 + Phase 18 + Phase 18.5 + Phase 18.6 + Phase 19 + Phase
19.5 + Phase 20 + Phase 22 + Phase 23 + Phase 24 + Phase 25 +
Phase 27 + Phase 27.5 + Phase 27.6 + Phase 28 + Phase 29 + Phase
30 + Phase 30.5 + Phase 31 + Phase 32 + Phase 34 + Phase 35 +
Phase 36 + Phase 38 + Phase 39 totals); Phase 33 opens the first
synthetic biology / genetic-circuit design domain to ship in
Valenx, broadening the design-protein → simulate-pathway → grow-
tissue → predict-structure loop to three synthetic-biology tools
(pySBOL, j5, Cello) feeding into the existing Phase 17 / 17.5
prediction stack, the Phase 32 systems-biology surface, the Phase
27 / 27.5 / 27.6 protein-design beachhead, and the entire Phase 28
→ 39 biology / biotech expansion. **Crosses the 100-adapter
milestone** alongside the Phase 5.5 MD-analysis-expansion trio.

### Phase 34 — Molecular docking · live
Add the de-facto open-source small-molecule docking pair to Valenx's
biology / chemistry stack: AutoDock Vina (the modern single-binary
docker) and AutoDock 4 (the older two-stage `autogrid4 → autodock4`
workflow that's still widely used in pharma teaching + tutorials).
Both adapters follow the established Phase 18 BWA shape — single-
action subprocess, file in / file out, no GPU required. AutoDock 4's
two-stage `prepare()` mirrors BWA's `bwa index` → `bwa mem` pattern.

- 2 new adapter crates land under `crates/valenx-adapters/bio/`:
  AutoDock Vina (single-binary subprocess; receptor PDBQT + ligand
  PDBQT in, ranked-pose PDBQT out; required `center [x, y, z]` and
  `size [x, y, z]` Å search-space knobs; `exhaustiveness` default 8
  range 1..=32; `num_modes` default 9; `energy_range` default 3.0
  kcal/mol; `cpu` default 0 = auto-detect; collects the output
  PDBQT as a `Native` artifact `"AutoDock Vina docked poses"`;
  Apache-2.0 licensed; `bio.vina.dock` ribbon capability) and
  AutoDock 4 (two-stage subprocess: `autogrid4 -p <gpf> -l
  <grid_log>` runs synchronously inside `prepare()`, `autodock4
  -p <dpf> -l <dock_log>` lands as the `PreparedJob.native_command`;
  `skip_grid` toggle reuses pre-generated grid maps; `grid_log`
  defaults `"autogrid4.glg"`, `dock_log` defaults `"autodock4.dlg"`;
  probe warns if `autogrid4` missing from PATH; collects `.dlg` /
  `.pdbqt` outputs alongside the dock log; GPL-2.0-or-later
  licensed; `bio.autodock4.dock` ribbon capability). Each wired
  into `valenx-app::init_registry`.
- No new canonical types, no new format readers, no new CLIs.
  PDBQT is a PDB-extension format the existing
  `valenx_bio::format::pdb` reader can already inspect for atom
  counts; ranked poses are user-readable as plain text. Vina's
  docked-pose PDBQT outputs and AutoDock 4's `.dlg` / `.pdbqt`
  outputs are inspectable through the existing `valenx-pdb-info`
  CLI.
- 2 new `valenx-init` templates ship: `vina` with alias
  `autodock-vina` (`vina-dock`), and `autodock4` with alias `ad4`
  (`autodock4-dock`). Cross-binary roundtrip test sweeps all 40
  templates clean.

The full per-phase shape lives in
[docs/src/phases/phase-34-docking.md](./docs/src/phases/phase-34-docking.md);
the implementation plan lives at
[docs/superpowers/plans/2026-04-30-docking.md](./docs/superpowers/plans/2026-04-30-docking.md).

**Deliverable:** Adapter inventory at 44 of 45 fully live after this
phase (2 new molecular-docking adapters added on top of the Phase 17
+ Phase 17.5 + Phase 18 + Phase 19 + Phase 23 + Phase 27 totals);
small-molecule docking workflows land alongside the existing biology
adapters in the same case-toml / prepare / run / collect shell, with
no glue code beyond the standard subprocess shape.

### Phase 35 — CRISPR design · live
Open the **first CRISPR guide-RNA design domain** in Valenx with
three established open-source tools that span the CRISPR-design
tradeoff space: CHOPCHOP (University of Bergen's web-and-script
CRISPR guide-RNA design tool — de-facto first stop in academic
CRISPR workflows; scores candidate gRNAs against a target sequence
under a configurable nuclease (Cas9, Cas12a, Cas13) or TALEN design
pass; MIT; Python-script subprocess shape sister to Phase 17
Biopython), CRISPOR (Maximilian Haeussler's CRISPR guide-RNA design
+ off-target prediction tool behind the public crispor.org service —
distinguishing feature is the rigorous off-target pass via the CFD
scoring model; supports many more enzymes / PAMs than CHOPCHOP;
GPL-3.0; Python-script subprocess shape sister to CHOPCHOP), and
Cas-OFFinder (Bae / Park / Kim group's CRISPR off-target searching
tool from Hanyang / Seoul National University — fast, OpenCL-
accelerated scanner that's the workhorse off-target scanner sitting
under most CRISPR design web services and pipelines; BSD-3-Clause;
single-binary subprocess shape sister to Phase 18 BWA with fixed-
shape positional CLI `cas-offinder <input> {C|G|A} <output>
[extras...]` — no `-i` / `-o` flags, the order is fixed). All three
follow established Phase 17 / 18 patterns: target FASTA / Python
script / Cas-OFFinder input file in, ranked guide / hit table out.
Phase 35 sits numerically before Phase 36 but ships chronologically
right after Phase 31 read simulators — same chronological-vs-
numerical convention used for Phase 17.5 / 24 / 28 / 31.

- 3 new adapter crates land under `crates/valenx-adapters/bio/`:
  CHOPCHOP (Python-script subprocess sister to Phase 17 Biopython;
  `case.toml` knobs `script` (path to user-supplied Python script
  that imports `chopchop` and reads `valenx_params.json`; required),
  `python` (interpreter name; default `"python3"`), `target` (target
  sequence FASTA; required), `genome` (CHOPCHOP-installed genome
  name — `"hg38"` / `"mm10"` / etc.; required), `cas_variant` ∈
  `{Cas9, Cas12a, Cas13, TALEN}` (required), `pam` (PAM sequence —
  `"NGG"` for Cas9, `"TTTV"` for Cas12a, etc.; required),
  `output_basename` (filename stem; required, non-empty);
  `prepare()` stages script + target FASTA into the workdir, writes
  a flat `valenx_params.json` containing `target` (staged filename) /
  `genome` / `cas_variant` / `pam` / `output_basename`, and composes
  `python <script_filename>` as the native command; collects
  `<output_basename>*.tsv` (`Tabular`, "CHOPCHOP guide rankings")
  and `<output_basename>*.bed` (`Tabular`, "CHOPCHOP guide
  locations"); probe via Python on PATH with an `import chopchop`
  check (returns `ok = true` with a warning when import fails so
  non-standard installs aren't blocked); MIT licensed; version range
  `3.0.0..4.0.0` (the modern web / script split landed in 3.0;
  upper bound 4.0 reserves room for an eventual major bump);
  `bio.chopchop.design` ribbon capability), CRISPOR (Python-script
  subprocess sister to CHOPCHOP; `case.toml` knobs `script` /
  `python` / `target` / `genome` / `pam` / `batch_id` (optional —
  CRISPOR caches partial results by batch so passing the same
  `batch_id` resumes a previously-interrupted run) / `output_basename`;
  `prepare()` stages script + target FASTA, writes a flat
  `valenx_params.json` containing `target` (staged filename) /
  `genome` / `pam` / `batch_id` (JSON string or literal `null` so
  user scripts can always do `params["batch_id"]` without an `in`
  check) / `output_basename`, and composes `python <script_filename>`
  as the native command; collects `<output_basename>*.tsv`
  (`Tabular`, "CRISPOR guide rankings") and `<output_basename>*.txt`
  (`Log`); probe via Python on PATH with an `import crispor` check
  (same `ok = true` + warning fallback as CHOPCHOP); GPL-3.0
  licensed; version range `5.0.0..6.0.0` (the modern Python 3 /
  batch-mode rewrite landed in 5.0; upper bound 6.0 reserves room
  for an eventual major bump); `bio.crispor.design` ribbon
  capability), and Cas-OFFinder (single-binary subprocess sister to
  Phase 18 BWA; `case.toml` knobs `input` (Cas-OFFinder input file —
  3+-line text file with reference genome path, PAM pattern, and
  one guide-sequence row per query; required), `output` (output text
  file; required), `backend` ∈ `{C, G, A}` (OpenCL device class —
  CPU / GPU / auto-pick fastest at runtime; required), `extra_args`;
  `prepare()` resolves both paths against the case directory when
  relative and composes `cas-offinder <input> <backend> <output>
  [extras...]` — Cas-OFFinder's CLI is purely positional, no `-i` /
  `-o` flags, the order is fixed; collects the configured `output`
  file as a single `Tabular` artifact ("Cas-OFFinder off-target
  hits"); probe via `find_on_path(&["cas-offinder"])`; the init
  alias `cas-off` resolves to the same template as the canonical
  `cas-offinder`; BSD-3-Clause licensed; version range
  `2.4.0..3.0.0` (the modern OpenCL device-selection CLI stabilised
  at 2.4; upper bound 3.0 reserves room for an eventual major
  bump); `bio.cas-offinder.search` ribbon capability). Each wired
  into `valenx-app::init_registry`.
- No new canonical types, no new format readers, no new CLIs. All
  three adapters consume user-supplied inputs (CHOPCHOP + CRISPOR
  target FASTAs and design Python scripts, Cas-OFFinder fixed-shape
  input files specifying the genome path + PAM + per-guide query
  lines) and emit user-readable artifacts (CHOPCHOP guide-ranking
  TSV + guide-location BED, CRISPOR guide-ranking TSV + log TXT,
  Cas-OFFinder ranked off-target hit TSV) that the unchanged
  `Results.artifacts` collection model surfaces directly. A first-
  class CRISPR-design canonical type — a generic guide / off-target
  / scoring type spanning all three back-ends — defers to a future
  phase along with guide-ranking visualizers and off-target heatmap
  viewers.
- 3 new `valenx-init` templates ship: `chopchop` (`chopchop-design`),
  `crispor` (`crispor-design`), and `cas-offinder` with alias
  `cas-off` (`cas-offinder-search`). Cross-binary roundtrip test
  sweeps all 80 templates clean.

The full per-phase shape lives in
[docs/src/phases/phase-35-crispr-design.md](./docs/src/phases/phase-35-crispr-design.md);
the implementation plan lives at
[docs/superpowers/plans/2026-04-30-crispr-design.md](./docs/superpowers/plans/2026-04-30-crispr-design.md).

**Deliverable:** Adapter inventory at 84 of 85 fully live after this
phase (3 new CRISPR-design adapters added on top of the Phase 17 +
Phase 17.5 + Phase 18 + Phase 18.5 + Phase 18.6 + Phase 19 + Phase
19.5 + Phase 20 + Phase 22 + Phase 23 + Phase 24 + Phase 25 + Phase
27 + Phase 27.5 + Phase 27.6 + Phase 28 + Phase 30 + Phase 31 +
Phase 32 + Phase 34 + Phase 36 totals); Phase 35 opens the first
CRISPR guide-RNA design domain to ship in Valenx, broadening the
design-guides → predict-off-targets → search-off-targets →
simulate-reads → align → quantify → call-variants → predict-
structure → fold-RNA → infer-tree → simulate-pathway → reconstruct-
3D → validate loop to three CRISPR-design tools (CHOPCHOP, CRISPOR,
Cas-OFFinder) feeding into the existing Phase 31 read simulators,
the Phase 18 / 18.5 / 18.6 alignment surface, the Phase 19 / 20
variant-calling + transcript-quantification stack, the Phase 17 /
17.5 prediction stack, the Phase 28 RNA-structure tools, the Phase
30 phylogenetic-tree builders, the Phase 32 systems-biology surface,
and the Phase 36 cryo-EM reconstruction tools.

### Phase 36 — Cryo-EM · live
Open the **first cryo-electron microscopy reconstruction domain** in
Valenx with three established open-source tools that span the cryo-EM
pipeline: RELION (Sjors Scheres' REgularised LIkelihood OptimisatioN
suite — de-facto Bayesian 3D reconstruction workhorse in cryo-EM
facilities worldwide; GPL-2.0; single-binary `relion_refine` for the
single-process path or `mpirun -n <N> relion_refine_mpi` for multi-
rank, since RELION ships separate `_mpi`-suffixed binaries),
EMAN2 (Steve Ludtke's broad-spectrum cryo-EM image-processing package
— "Swiss army knife" of single-particle cryo-EM, BSD-3-Clause; high-
level driver `e2refine_easy.py` orchestrates particle picking, 2D
classification, initial-model building, and 3D refinement across the
sprawling `e2*.py` toolkit), and CTFFIND (Niko Grigorieff's contrast
transfer function estimation tool — gold standard for fitting per-
micrograph CTF parameters that RELION, cryoSPARC, EMAN2, and most
automated pipelines all wrap as a preprocessing step; Janelia non-
commercial / academic-only license, surfaced as `Janelia-License`
and flagged via mandatory `"academic"`-keyworded probe warning;
single-binary `ctffind` with stdin-piped parameters since the CLI is
interactive). All three follow the established Phase 18 BWA single-
binary CLI pattern: particles / micrographs / reference maps in,
typed artifacts out. Phase 36 sits numerically after Phase 34 and
ships chronologically right after Phase 32 systems biology — same
chronological-vs-numerical convention used for Phase 17.5 / 24 / 28.

- 3 new adapter crates land under `crates/valenx-adapters/bio/`:
  RELION (single-binary subprocess with optional MPI wrapping;
  `case.toml` knobs `particles` (`*_data.star` particle STAR file;
  required), `reference` (initial reference map `.mrc`; required),
  `output_basename` (becomes the `--o` prefix every output inherits
  so collect() walks deterministically; required), `angpix` (pixel
  size in Angstroms; required, > 0.0 and finite), `mpi_procs`
  (default 1, ≥ 1; > 1 switches to the MPI binary), `threads`
  (OpenMP threads per MPI rank, default 1, ≥ 1), `extra_args`;
  `prepare()` dispatches on `mpi_procs`: single-rank composes
  `relion_refine --i <particles> --ref <reference> --o
  <output_basename> --angpix <angpix> --j <threads> [extras...]`;
  multi-rank prepends `mpirun -n <mpi_procs> relion_refine_mpi ...`
  and surfaces a helpful install-hint `InvalidCase` ("install
  OpenMPI (`apt install openmpi-bin`, `brew install open-mpi`) or
  MPICH (`apt install mpich`) to enable multi-rank RELION runs") if
  `mpirun` isn't on PATH; collects `<output_basename>*_class*.mrc`
  (`Native`, "RELION class average"), `<output_basename>*_data.star`
  (`Tabular`, "RELION particle assignments"),
  `<output_basename>*_model.star` (`Log`, "RELION model summary");
  probe via `find_on_path(&["relion_refine"])`; GPL-2.0 licensed;
  version range `4.0.0..6.0.0` (4.0 is the current stable line,
  predecessor 3.1; upper bound 6.0 reserves room for the next
  major); `bio.relion.refine` ribbon capability), EMAN2 (single-
  binary subprocess wrapping `e2refine_easy.py`, EMAN2's high-level
  orchestrator that drives the rest of the toolkit; `case.toml`
  knobs `particles` (`.bdb` / `.hdf` / `.mrcs`; required), `model`
  (initial 3D model `.hdf` / `.mrc`; required), `output_basename`
  (becomes the `--path` argument; EMAN2 turns this into a
  `<basename>_NN/` results directory under the workdir; required),
  `target_resolution` (Å; required, > 0.0 and finite), `symmetry`
  (point group — `"c1"` / `"d2"` / `"icos"` / etc.; required,
  default `"c1"`), `threads` (default 1, ≥ 1), `extra_args`;
  `prepare()` builds `e2refine_easy.py --input <particles> --model
  <model> --path <output_basename> --targetres <target_resolution>
  --sym <symmetry> --threads <threads> [extras...]`; collects
  `<output_basename>_*/threed_*.hdf` (`Native`, "EMAN2
  reconstruction") and `<output_basename>_*/log.txt` (`Log`, "EMAN2
  log"); probe via `find_on_path(&["e2refine_easy.py"])`; init alias
  `eman` resolves to the same template as canonical `eman2`;
  BSD-3-Clause licensed; version range `2.99.0..3.0.0` (2.99 is the
  current pre-3.0 stable release; upper bound 3.0 reserves room for
  the long-rumoured 3.x line); `bio.eman2.refine` ribbon
  capability), and CTFFIND (single-binary subprocess with stdin-
  piped parameters; `case.toml` knobs `micrograph` (input `.mrc`;
  required), `output_diagnostic` (output diagnostic `.mrc`;
  required), `output_txt` (output text file with CTF parameters;
  required), `pixel_size` (Å; required, > 0.0 and finite), `voltage`
  (kV; default 300.0, > 0.0), `cs` (spherical aberration mm;
  default 2.7, > 0.0), `amplitude_contrast` (fraction in
  `0.0..=1.0`; required — 0.07 typical for cryo, 0.1 for negative
  stain), `extra_args`; `prepare()` writes `ctffind_params.txt`
  containing one parameter per CTFFIND-v4.1 prompt in order (input
  image, output diagnostic, pixel size, voltage, Cs, amplitude
  contrast, plus standard defaults for box size / min res / max
  res / defocus search / expert sub-prompts) and stashes the
  filename under a sentinel env var (`VALENX_CTFFIND_PARAMS_FILE`);
  the custom `run()` recovers the filename, strips the sentinel
  from the env table so CTFFIND doesn't see it, opens the params
  file with `File::open()`, and hands the FD to the child via
  `Stdio::from(file)` — CTFFIND sees a pipe pre-loaded with one
  parameter per prompt and responds as if a human had typed each
  line (the shared `subprocess::run` helper closes stdin which
  makes CTFFIND read EOF before its first prompt and exit; the
  custom run path mirrors the MAFFT stdout-redirect pattern but
  for stdin); collects `output_diagnostic` (`Native`, "CTFFIND
  diagnostic image") and `output_txt` (`Tabular`, "CTFFIND
  parameters"); probe via `find_on_path(&["ctffind"])` pushes the
  literal string `"academic"` into `ProbeReport.warnings` with the
  full reminder ("CTFFIND is licensed for non-commercial /
  academic use only. Confirm your use case complies with the
  Janelia license before redistributing CTF estimates or derived
  data."); tool license surfaces as `Janelia-License` rather than
  mislabeling as MIT / BSD; version range `4.1.0..5.0.0` (CTFFIND4
  is the long-running stable line; upper bound 5.0 reserves room
  for the announced CTFFIND5 line); `bio.ctffind.estimate` ribbon
  capability). Each wired into `valenx-app::init_registry`.
- No new canonical types, no new format readers, no new CLIs. All
  three adapters consume user-supplied inputs (RELION particle STAR
  files + reference MRC volumes, EMAN2 particle stacks + initial 3D
  models, CTFFIND micrograph `.mrc` files) and emit user-readable
  artifacts (RELION class-average MRC volumes, particle-assignment
  STAR files, model-summary STAR files, EMAN2 `threed_*.hdf`
  reconstructions plus per-run log files, CTFFIND diagnostic-image
  MRC plus per-micrograph parameter text files) that the unchanged
  `Results.artifacts` collection model surfaces directly. A first-
  class cryo-EM canonical type — a generic `.mrc` volume / particle-
  stack / micrograph type spanning all three back-ends — defers to
  a future phase along with MRC readers and reconstruction
  visualizers.
- 3 new `valenx-init` templates ship: `relion` (`relion-refine`),
  `eman2` with alias `eman` (`eman2-refine`), and `ctffind`
  (`ctffind-estimate`; carries an inline academic-license note in
  the scaffolded `case.toml`). Cross-binary roundtrip test sweeps
  all 74 templates clean.

The full per-phase shape lives in
[docs/src/phases/phase-36-cryo-em.md](./docs/src/phases/phase-36-cryo-em.md);
the implementation plan lives at
[docs/superpowers/plans/2026-04-30-cryo-em.md](./docs/superpowers/plans/2026-04-30-cryo-em.md).

**Deliverable:** Adapter inventory at 78 of 79 fully live after this
phase (3 new cryo-EM adapters added on top of the Phase 17 + Phase
17.5 + Phase 18 + Phase 18.5 + Phase 18.6 + Phase 19 + Phase 19.5 +
Phase 20 + Phase 22 + Phase 23 + Phase 24 + Phase 25 + Phase 27 +
Phase 27.5 + Phase 27.6 + Phase 28 + Phase 30 + Phase 32 + Phase 34
totals); Phase 36 opens the first cryo-electron microscopy
reconstruction domain to ship in Valenx, broadening the estimate-CTF
→ process-images → reconstruct-3D → predict-structure → fold-RNA →
infer-tree → simulate-pathway → validate loop to three cryo-EM tools
(RELION, EMAN2, CTFFIND) feeding into the existing Phase 32 systems-
biology surface, the Phase 17 / 17.5 prediction stack, the Phase 24
/ 25 cheminformatics + quantum-chemistry surface, the Phase 28 RNA-
structure tools, and the Phase 30 phylogenetic-tree builders.

### Phase 38 — Rosetta family · live
Open the **first Rosetta protein-modeling family** in Valenx with
the two most-used entry points into the RosettaCommons code base:
Rosetta (RosettaCommons' flagship modeling suite — drives protein
design, structure prediction, docking, ligand binding, and a long
tail of related modeling tasks through `rosetta_scripts`, the
XML-driven protocol runner that's the de-facto Rosetta entry point
in production: every `relax` / `dock` / `abinitio` / FastDesign /
enzyme-design pipeline lives as an XML protocol fed to this binary;
custom Rosetta-License — academic / non-commercial use only,
surfaced as `Rosetta-License` and flagged via mandatory
`"academic"`-keyworded probe warning; single-binary subprocess
shape sister to Phase 18 BWA), and PyRosetta (Python bindings to
the Rosetta C++ core — exposes the entire Rosetta modeling
pipeline (movers, filters, scorefunctions, task-operations) through
a Pythonic API, letting users drive Rosetta from regular `.py`
scripts rather than authoring XML protocols; inherits the same
academic / non-commercial use terms; Python-script subprocess shape
sister to Phase 17 Biopython). Both adapters surface the
RosettaCommons license accurately via `tool_license =
"Rosetta-License"` (a custom non-OSS license — not a recognised
SPDX identifier) and emit a probe warning whenever the binary /
bindings are detected, with the literal `"academic"` string baked
into the warning as a stable anchor for license-aware filters and
tests. Phase 38 sits numerically after Phase 36 cryo-EM and ships
chronologically right after Phase 29 population genetics.

- 2 new adapter crates land under `crates/valenx-adapters/bio/`:
  Rosetta (single-binary subprocess sister to Phase 18 BWA;
  `case.toml` knobs `protocol` (XML protocol script driving
  `rosetta_scripts`; required), `input_pdb` (input PDB Rosetta will
  operate on; required), `output_basename` (filename stem the
  binary uses to label output decoys — `<basename>_0001.pdb` etc.;
  required, non-empty), `nstruct` (number of independent decoys to
  generate; `u32`, ≥ 1), `database` (path to the Rosetta
  `database/` directory — required because every `rosetta_scripts`
  invocation needs `-database <path>` pointing at the energy
  tables / fragment libraries / etc. bundled with the source
  distribution), `extra_args`; `prepare()` resolves the protocol,
  input PDB, and database paths against the case directory when
  relative, validates the protocol + PDB exist on disk (returns
  `InvalidCase` with a helpful message when missing), and composes
  `rosetta_scripts -database <path> -parser:protocol <protocol>
  -in:file:s <input_pdb> -out:prefix <output_basename> -nstruct <N>
  [extras...]`; collects `<output_basename>*.pdb` (`Native`,
  "Rosetta designed structure") plus the canonical `score.sc`
  scorefile (`Tabular`, "Rosetta scores"); probe via
  `find_on_path(&["rosetta_scripts",
  "rosetta_scripts.linuxgccrelease",
  "rosetta_scripts.macosclangrelease"])` — Rosetta source builds
  emit platform-suffixed names by default, conda / packaged
  distributions install a bare `rosetta_scripts` shim, and the
  probe covers all three; **academic-license-only** — probe always
  pushes an `"academic"`-keyworded warning into
  `ProbeReport.warnings` whenever the binary is detected, and
  `tool_license` surfaces as `"Rosetta-License"` rather than
  mislabeling the custom RosettaCommons terms as a recognised SPDX
  identifier; version range `3.13.0..4.0.0` (the stable 3.x line
  landed at 3.13 in 2021; upper bound 4.0 reserves room for an
  eventual major bump); `bio.rosetta.protocol` ribbon capability),
  and PyRosetta (Python-script subprocess sister to Phase 17
  Biopython; `case.toml` knobs `script` (path to user-authored
  Python script; required), `python` (interpreter name; default
  `"python3"`), `input_pdb` (optional input PDB the script will
  operate on — None when the script generates structures de novo;
  surfaced in `valenx_params.json` so the script can read it
  without re-parsing case.toml), `output_basename` (filename stem;
  required, non-empty); `prepare()` stages the script (and PDB,
  when present) into the workdir under their original filenames so
  the script can resolve them via relative paths, then writes a
  flat `valenx_params.json` with `input_pdb` (staged filename or
  literal `null` so user scripts can always do
  `params["input_pdb"]` without an `in` check) and
  `output_basename`; collects `<output_basename>*.pdb` (`Native`,
  "PyRosetta designed structure") and `*.sc` files (`Tabular`,
  "PyRosetta scores"); probe via Python on PATH with an
  `import pyrosetta` check (returns `ok = true` with a warning
  when import fails so non-standard installs aren't blocked);
  **academic-license-only** — probe always pushes an
  `"academic"`-keyworded warning into `ProbeReport.warnings`
  whenever Python is detected (regardless of whether `pyrosetta`
  itself is importable, since the user is either about to install
  it or has it installed and needs reminding); version range
  `4.0.0..5.0.0` (the modern release line is the 4.x series with
  weekly nightly drops post-2017; upper bound 5.0 reserves room
  for an eventual major bump); `bio.pyrosetta.script` ribbon
  capability). Each wired into `valenx-app::init_registry`.
- No new canonical types, no new format readers, no new CLIs.
  Both adapters consume user-supplied inputs (Rosetta XML
  protocols + input PDBs + `database/` data directories,
  PyRosetta Python scripts + optional input PDBs) and emit user-
  readable artifacts (`<basename>*.pdb` design decoys, `score.sc`
  / `*.sc` scorefiles) that the unchanged `Results.artifacts`
  collection model surfaces directly. The existing
  `valenx_bio::format::pdb` reader inspects collected PDB
  artifacts for chain / residue / atom counts. A first-class
  Rosetta canonical type — a generic protocol + scorefile pair
  spanning both back-ends, parsed into a typed scorefile model
  with per-decoy energy terms — defers to a future phase along
  with score-distribution visualizers and per-mutation Δ-energy
  heatmap viewers.
- 2 new `valenx-init` templates ship: `rosetta`
  (`rosetta-protocol`) and `pyrosetta` (`pyrosetta-script`).
  Cross-binary roundtrip test sweeps all 85 templates clean.

The full per-phase shape lives in
[docs/src/phases/phase-38-rosetta.md](./docs/src/phases/phase-38-rosetta.md);
the implementation plan lives at
[docs/superpowers/plans/2026-04-30-rosetta.md](./docs/superpowers/plans/2026-04-30-rosetta.md).

**Deliverable:** Adapter inventory at 89 of 90 fully live after this
phase (2 new Rosetta-family adapters added on top of the Phase 17 +
Phase 17.5 + Phase 18 + Phase 18.5 + Phase 18.6 + Phase 19 + Phase
19.5 + Phase 20 + Phase 22 + Phase 23 + Phase 24 + Phase 25 + Phase
27 + Phase 27.5 + Phase 27.6 + Phase 28 + Phase 29 + Phase 30 +
Phase 31 + Phase 32 + Phase 34 + Phase 35 + Phase 36 totals); Phase
38 opens the first Rosetta protein-modeling family to ship in
Valenx, broadening the design-guides → predict-off-targets →
search-off-targets → simulate-reads → align → quantify → call-
variants → predict-structure → fold-RNA → infer-tree → simulate-
pathway → reconstruct-3D → simulate-popgen → analyze-trees →
validate loop to two Rosetta entry points (Rosetta, PyRosetta)
feeding into the existing Phase 27 / 27.5 / 27.6 protein-design
adapters (RFdiffusion, ProteinMPNN, Chroma, ESM-IF, RFantibody,
ESM3, ESM Cambrian), the Phase 17 / 17.5 prediction stack, the
Phase 28 RNA-structure tools, the Phase 29 population-genetics
trio, the Phase 30 phylogenetic-tree builders, the Phase 32
systems-biology surface, the Phase 34 docking pair, the Phase 35
CRISPR-design tools, and the Phase 36 cryo-EM reconstruction tools.

### Phase 39 — DNA structural geometry · live
Open the **first DNA structural-geometry domain** in Valenx with
three established open-source tools that span the structural-
geometry tradeoff space: X3DNA (Wilma Olson and Xiang-Jun Lu's
reference toolkit for DNA / RNA structural-geometry analysis —
de-facto reference for canonical helical-step parameters (twist,
roll, tilt, slide, shift, rise) plus per-base intra-pair
parameters (buckle, propeller, opening, shear, stretch, stagger);
custom X3DNA-License — academic / non-commercial use only,
surfaced as `X3DNA-License` and flagged via mandatory
`"academic"`-keyworded probe warning; single-binary `analyze`
positional CLI), Curves+ (Richard Lavery's reference toolkit for
DNA helical-axis analysis — fits a curvilinear helical axis
through a nucleic-acid structure and reports per-base axis-
curvature, base-pair parameters relative to that axis, and a
`.cda` file describing the axis itself; the canonical tool for
"is this DNA bent, and if so, how" questions in protein-DNA /
drug-DNA structural studies; custom Curves-License — academic /
non-commercial use only; single-binary `Cur+` CLI with stdin-
piped namelist parameters since the binary takes its parameters
as a Fortran-style `&inp ... &end` block on stdin), and DSSR
(Dissecting the Spatial Structure of RNA / DNA — the modern
Python-fronted X3DNA-family tool; reads a nucleic-acid PDB and
emits a single JSON file enumerating every detected structural
feature: base pairs, multiplets, helices, stems, hairpin /
internal / junction loops, kissing loops, A-minor motifs, ribose
zippers, pseudoknots; the standard machine-readable feature-
extraction step in modern RNA-structure pipelines; custom
DSSR-License — academic / non-commercial use only; single-binary
`x3dna-dssr` CLI). All three follow the established Phase 18 BWA
single-binary CLI pattern: nucleic-acid PDB in, typed structural-
geometry artifacts out (X3DNA / DSSR positional or flag-form CLIs;
Curves+ stdin-fed via the Phase 36 CTFFIND-style stdin-feed
pattern with `Stdio::from(file)`). All three are **academic-
license-flagged** — Valenx surfaces the X3DNA-family + Curves+
licenses accurately via `tool_license = "X3DNA-License"` /
`"Curves-License"` / `"DSSR-License"` and a mandatory
`"academic"`-keyworded probe warning whenever each binary is
detected. Phase 39 sits numerically after Phase 38 cryo-EM
Rosetta and ships chronologically right after Phase 30.5 Bayesian
phylogenetics — same chronological-vs-numerical convention used
for Phase 17.5 / 24 / 28 / 31 / 35.

- 3 new adapter crates land under `crates/valenx-adapters/bio/`:
  X3DNA (single-binary subprocess sister to Phase 18 BWA;
  `case.toml` knobs `input_pdb` (input PDB; required),
  `output_basename` (filename stem the user expects X3DNA to
  produce — surfaced here so collect() can label artefacts
  uniformly without scraping `analyze`'s filename heuristics;
  required, non-empty), `extra_args`; `prepare()` resolves the
  input PDB against the case directory when relative, validates
  it exists on disk (returns `InvalidCase` with a helpful message
  when missing), and composes `analyze <input_pdb> [extras...]`
  (`analyze` is positional-only — derives every output filename
  from the input basename, so the adapter just hands it the PDB
  and any user-supplied extras); collects `<output_basename>*.par`
  (`Tabular`, "X3DNA base-step parameters") and `*.out` (`Log`,
  the per-run log `analyze` writes alongside); probe via
  `find_on_path(&["analyze"])` (X3DNA's main analysis binary is
  literally named `analyze`); **academic-license-only** — probe
  always pushes an `"academic"`-keyworded warning into
  `ProbeReport.warnings` whenever the binary is detected, and
  `tool_license` surfaces as `"X3DNA-License"` rather than
  mislabeling the custom X3DNA terms as a recognised SPDX
  identifier; version range `2.4.0..3.0.0` (X3DNA 2.4 (2020) is
  the modern stable release and the floor we test against);
  `bio.x3dna.analyze` ribbon capability), Curves+ (single-binary
  subprocess with stdin-piped parameters sister to Phase 36
  CTFFIND; `case.toml` knobs `input_pdb` (input PDB; required),
  `output_basename` (filename stem Curves+ uses for outputs —
  `<basename>.lis`, `<basename>.cda`, etc.; required, non-empty),
  `first_residue` (`u32` — first inclusive residue index in the
  strand to analyse; required), `last_residue` (`u32`, ≥
  `first_residue` — a reverse range is rejected up front with a
  helpful message; required), `extra_args`; `prepare()` resolves
  the input PDB against the case directory when relative,
  validates it exists on disk, writes `curves_params.txt`
  containing the namelist body + residue-range cards, stashes the
  filename under the sentinel env var `VALENX_CURVES_PARAMS_FILE`,
  and the custom `run()` recovers the filename, strips the
  sentinel from the env table so Curves+ doesn't see it, opens
  the params file with `File::open()`, and hands the FD to the
  child via `Stdio::from(file)` — the shared `subprocess::run`
  helper closes stdin which makes Curves+ read EOF before parsing
  its first parameter and exit; the custom run path mirrors the
  MAFFT stdout-redirect pattern but for stdin (same shape Phase
  36 CTFFIND uses); collects `<output_basename>*.lis` (`Log`,
  "Curves+ helical analysis") and `<output_basename>*.cda`
  (`Tabular`, "Curves+ axis curve data"); probe via
  `find_on_path(&["Cur+"])` (the binary name uses a literal `+`);
  **academic-license-only** — probe always pushes an
  `"academic"`-keyworded warning into `ProbeReport.warnings` and
  `tool_license` surfaces as `"Curves-License"`; version range
  `2.0.0..3.0.0` (Curves+ 2.x is the modern stable line; 2.0 is
  the floor); `bio.curves.analyze` ribbon capability), and DSSR
  (single-binary subprocess sister to X3DNA; `case.toml` knobs
  `input_pdb` (input PDB; required), `output_json` (output JSON
  path; required), `extra_args`; `prepare()` resolves the input
  PDB against the case directory when relative, scopes the output
  JSON path to the workdir when relative, validates the input
  exists on disk, and composes `x3dna-dssr -i=<input_pdb>
  -o=<output_json> --json [extras...]` (DSSR uses `key=value`
  flag form on its short-form options — no space between flag and
  value); collects the configured `output_json` file as a single
  `Tabular` artifact ("DSSR analysis (JSON)") — DSSR's JSON is
  the canonical machine-readable summary; tagged `Tabular` rather
  than `Native` so downstream serdes can key off a consistent
  kind; probe via `find_on_path(&["x3dna-dssr"])`; **academic-
  license-only** — probe always pushes an `"academic"`-keyworded
  warning into `ProbeReport.warnings` and `tool_license` surfaces
  as `"DSSR-License"` rather than mislabeling the inherited
  X3DNA-family terms as a recognised SPDX identifier; version
  range `2.0.0..3.0.0` (DSSR 2.x is the modern stable line that
  ships with X3DNA 2.4+); `bio.dssr.analyze` ribbon capability).
  Each wired into `valenx-app::init_registry`.
- No new canonical types, no new format readers, no new CLIs.
  All three adapters consume user-supplied inputs (X3DNA / Curves+
  / DSSR all take nucleic-acid PDBs, plus the Curves+ residue-
  range knobs) and emit user-readable artifacts (X3DNA `.par`
  base-step parameter tables and `.out` per-run logs, Curves+
  `.lis` helical-analysis logs and `.cda` axis-curve data files,
  DSSR JSON structural-feature summaries) that the unchanged
  `Results.artifacts` collection model surfaces directly. The
  existing `valenx_bio::format::pdb` reader inspects collected PDB
  inputs for chain / residue / atom counts. A first-class
  DNA-geometry canonical type — a typed helical-parameter
  representation spanning all three back-ends, with parsed per-
  step parameter tables and a typed structural-feature summary —
  defers to a future phase along with helical-axis visualizers
  and per-feature interactive overlays.
- 3 new `valenx-init` templates ship: `x3dna` with alias `3dna`
  (`x3dna-analyze`), `curves` with alias `curves+` (`curves-
  analyze`; carries an inline academic-license note in the
  scaffolded `case.toml`), and `dssr` (`dssr-analyze`; carries an
  inline academic-license note in the scaffolded `case.toml`).
  Cross-binary roundtrip test sweeps all 90 templates clean.

The full per-phase shape lives in
[docs/src/phases/phase-39-dna-geometry.md](./docs/src/phases/phase-39-dna-geometry.md);
the implementation plan lives at
[docs/superpowers/plans/2026-04-30-dna-geometry.md](./docs/superpowers/plans/2026-04-30-dna-geometry.md).

**Deliverable:** Adapter inventory at 94 of 95 fully live after this
phase (3 new DNA-structural-geometry adapters added alongside the
Phase 30.5 Bayesian-phylogenetics pair on top of the Phase 17 +
Phase 17.5 + Phase 18 + Phase 18.5 + Phase 18.6 + Phase 19 + Phase
19.5 + Phase 20 + Phase 22 + Phase 23 + Phase 24 + Phase 25 + Phase
27 + Phase 27.5 + Phase 27.6 + Phase 28 + Phase 29 + Phase 30 +
Phase 31 + Phase 32 + Phase 34 + Phase 35 + Phase 36 + Phase 38
totals); Phase 39 opens the first DNA structural-geometry domain to
ship in Valenx, broadening the predict-structure → fold-RNA →
analyze-DNA-geometry → infer-tree-ML → infer-tree-Bayesian →
simulate-popgen → analyze-trees → simulate-pathway → reconstruct-
3D → design-protein → validate loop to three DNA structural-
geometry tools (X3DNA, Curves+, DSSR) feeding into the existing
Phase 17 / 17.5 prediction stack, the Phase 28 RNA-structure
tools, the Phase 29 population-genetics trio, the Phase 30
phylogenetic-tree builders, the Phase 30.5 Bayesian-
phylogenetics pair, the Phase 32 systems-biology surface, the
Phase 34 docking pair, the Phase 35 CRISPR-design tools, the
Phase 36 cryo-EM reconstruction tools, and the Phase 38 Rosetta-
family adapters.

### Phase 40 — Microscopy · live ✅
Open the **first microscopy / bioimage analysis domain** in
Valenx with three established open-source tools that span the
bioimage analysis tradeoff space — script-driven general-purpose
image processing in headless mode (Fiji, the canonical ImageJ
distribution that's the de-facto first stop in bioimage
analysis), pipeline-driven cell segmentation + measurement
(CellProfiler, the Broad Institute pipeline-driven workhorse
that powers most high-content-screening assays), and
interactive-ML pixel / object classification (Ilastik, the
Hamprecht lab tool that leans on user-trained random-forest
classifiers for hard segmentation tasks where rule-based
pipelines struggle). All three run in headless mode for batch
processing — Fiji + Ilastik via app-launcher binaries that ship
in the upstream distribution, CellProfiler via its Python CLI.
No new canonical types — image inputs and outputs are all
standard formats (TIFF, PNG, HDF5, CSV) that the existing
`Results.artifacts` collection model surfaces directly. Phase 40
sits numerically between Phase 39 DNA structural geometry and
Phase 41 sequence editors and ships chronologically right after
Phase 32.5 spatial stochastic — same chronological-vs-numerical
convention used for Phase 17.5 / 24 / 28 / 31 / 35 / 39 / 5.5 /
5.6 / 5.7 / 32.5.

- 3 new adapter crates land under `crates/valenx-adapters/bio/`:
  Fiji (app-launcher subprocess sister to Phase 36 RELION /
  EMAN2; `case.toml` knobs `fiji_app` (absolute path to the
  per-platform Fiji launcher — `ImageJ-linux64`, `ImageJ.exe`,
  `Contents/MacOS/ImageJ-macosx`; required), `macro_file` (`.ijm`
  Fiji macro file describing the image-processing pipeline;
  required), `input_image` (`Option<PathBuf>` — optional input
  image the macro will operate on, typically picked up via
  `getArgument()`), `output_basename` (filename stem the user's
  macro uses for outputs — surfaced here so collect() can label
  artifacts uniformly; required, non-empty), `extra_args`;
  `prepare()` resolves `fiji_app`, `macro_file`, and the optional
  `input_image` against the case directory when relative,
  validates each existing file exists on disk, builds `<fiji_app>
  --headless --console -macro <macro_file> [extras...]`;
  `collect()` walks for `<output_basename>*.tif` / `.tiff`
  (`Native`, "Fiji image (TIFF)"), `<output_basename>*.png`
  (`Native`, "Fiji image (PNG)"), `<output_basename>*.csv`
  (`Tabular`, "Fiji measurements"), `*.log`; probe via
  `find_on_path(&["ImageJ-linux64", "ImageJ-macosx",
  "ImageJ.exe", "fiji"])` — surfaces a Java-fallback warning
  whenever Fiji isn't on PATH but `java` is; GPL-3.0 licensed;
  version range `2.0.0..3.0.0`; `bio.fiji.process` ribbon
  capability), CellProfiler (Python-CLI subprocess; `case.toml`
  knobs `pipeline` (`.cppipe` / `.cpproj` pipeline file;
  required), `input_dir` (directory containing input images;
  required — the adapter validates it is a directory at prepare
  time), `output_basename` (filename stem the adapter pins as the
  `-o` output directory; required, non-empty), `python`
  (interpreter name; default `"python3"` — used for the `<python>
  -m cellprofiler ...` fallback when the launcher isn't on PATH),
  `extra_args`; `prepare()` looks up the `cellprofiler` binary on
  PATH first then falls back to `<python> -m cellprofiler ...`,
  builds `cellprofiler -c -r -p <pipeline> -i <input_dir> -o
  <basename> [extras...]`; `collect()` walks **one level deep**
  into `<output_basename>/` for `*.csv` (`Tabular`, "CellProfiler
  measurements"), `*.tif` / `.tiff` (`Native`, "CellProfiler
  segmented image"), `*.png` (`Native`, "CellProfiler plot");
  top-level `*.log` (`Log`, "CellProfiler log"); probe via
  `find_on_path(&["cellprofiler", "python3", "python"])` with
  warning when `cellprofiler` itself isn't on PATH but Python is;
  BSD-3-Clause licensed; version range `4.0.0..5.0.0`;
  `bio.cellprofiler.segment` ribbon capability), and Ilastik
  (app-launcher subprocess sister to Fiji; `case.toml` knobs
  `ilastik_app` (absolute path to per-platform Ilastik launcher
  — `ilastik`, `run_ilastik.sh`, or `ilastik.exe`; required),
  `project` (`.ilp` Ilastik project file containing the trained
  classifier; required), `input_images` (`Vec<PathBuf>` — must
  contain ≥ 1 entry; the adapter rejects an empty vector at
  prepare time), `output_basename` (filename stem; required,
  non-empty), `workflow` (string; default `"Pixel Classification"`
  — selectable from Ilastik's set: `"Pixel Classification"`,
  `"Object Classification"`, etc.), `extra_args`; `prepare()`
  resolves `ilastik_app`, `project`, and each `input_images`
  entry against the case directory when relative, validates
  `input_images` is non-empty, builds `<ilastik_app> --headless
  --project=<project>
  --output_filename_format=<basename>_{nickname}.h5
  <input_images...> [extras...]` with `--project=` and
  `--output_filename_format=` flags emitted as single OsString
  args each (so `=` and the value travel together) and the
  literal `{nickname}` substring preserved unmodified for
  Ilastik's per-image disambiguation; `collect()` walks for
  `<output_basename>*.h5` (`Native`, "Ilastik probability map
  (HDF5)"), `<output_basename>*.tif` (`Native`, "Ilastik
  segmentation"), `*.log`; probe via `find_on_path(&["ilastik",
  "run_ilastik.sh", "ilastik.exe"])` with warning when nothing
  matches but still returns `ok = true` since the user can
  supply the launcher via `case.toml` (sister to Phase 32
  PhysiCell's per-project-binary probe convention); GPL-3.0
  licensed; version range `1.4.0..2.0.0`;
  `bio.ilastik.classify` ribbon capability). Each wired into
  `valenx-app::init_registry`.
- No new canonical types, no new format readers, no new CLIs.
- 3 new `valenx-init` templates ship: `fiji`, `cellprofiler`,
  `ilastik`. Cross-binary roundtrip test sweeps all 115 templates
  clean.

The full per-phase shape lives in
[docs/src/phases/phase-40-microscopy.md](./docs/src/phases/phase-40-microscopy.md);
the implementation plan lives at
[docs/superpowers/plans/2026-05-03-microscopy.md](./docs/superpowers/plans/2026-05-03-microscopy.md).

**Deliverable:** Adapter inventory at 119 of 120 fully live after
this phase trio (3 new microscopy adapters added on top of the
prior set, alongside the Phase 32.5 spatial-stochastic pair and
Phase 41 sequence-editors pair); Phase 40 opens the first
microscopy / bioimage analysis domain in Valenx.

### Phase 41 — Sequence editors · live ✅
Open the **first plasmid-design / alignment-viewer domain** in
Valenx with two established open-source tools that span the
sequence-editor tradeoff space — a Python plasmid-design library
that handles PCR primer design, restriction-enzyme digests, and
Gibson / Golden-Gate assembly programmatically (pydna, Bjorn
Johansson's BSD-3-Clause library that's the de-facto Python
choice for cloning automation), and the canonical Java alignment
viewer with a headless mode for batch image / format conversion
(Jalview, the Barton group's GPL-3.0 viewer that's been the
reference alignment viewer in molecular biology labs since the
2000s and supports headless operation for unattended pipeline
integration). pydna follows the established Phase 17 Biopython
Python-script subprocess shape: the user supplies a Python script
that imports the upstream package and reads `valenx_params.json`
for the parsed knobs. Jalview is JAR-distributed (no `jalview`
launcher binary on PATH); the user supplies the absolute path
to the JAR via case input, and we probe `java` itself rather
than the JAR — same JAR-distribution shape Phase 33 j5 / Cello
use. Phase 41 sits numerically after Phase 40 microscopy and
ships chronologically right after Phase 40 — same chronological-
vs-numerical convention used for Phase 17.5 / 24 / 28 / 31 / 35
/ 39 / 5.5 / 5.6 / 5.7 / 32.5 / 40.

- 2 new adapter crates land under `crates/valenx-adapters/bio/`:
  pydna (Python-script subprocess sister to Phase 17 Biopython
  / Phase 19.5 Scanpy / Phase 33 pySBOL; `case.toml` knobs
  `script` (path to user-supplied Python script; required, `.py`
  enforced), `python` (interpreter name; default `"python3"`),
  `input_genbank` (`Option<PathBuf>` — optional starting GenBank
  file the script can use as the parent / template construct;
  `None` when the script generates the design from scratch),
  `output_basename` (filename stem the user's script uses for
  outputs — surfaced here so collect() can label artifacts
  uniformly; required, non-empty); `prepare()` enforces the `.py`
  extension on the script, resolves `script` and the optional
  `input_genbank` against the case directory when relative,
  stages both into the workdir under their original filenames,
  writes a flat `valenx_params.json` containing `output_basename`
  always plus `input_genbank` (staged filename) only when set —
  key omitted entirely when `None` rather than emitted as `null`,
  matching the hand-rolled JSON convention the rest of the bio
  adapters use (Phase 19.6 Seurat / AnnData, Phase 27.5 ESM-IF),
  builds `<python> <staged_script>`; `collect()` walks for
  `<output_basename>*.gb` / `.genbank` (`Native`, "pydna GenBank
  file"), `<output_basename>*.fasta` (`Native`, "pydna FASTA"),
  `<output_basename>*.csv` (`Tabular`, "pydna table"), `*.log`;
  probe via Python on PATH with `import pydna` check (returns
  `ok = true` with a warning when import fails so non-standard
  installs aren't blocked); BSD-3-Clause licensed; version range
  `5.0.0..7.0.0`; `bio.pydna.design` ribbon capability), and
  Jalview (JAR-distributed single-binary subprocess sister to
  Phase 33 j5 / Cello; `case.toml` knobs `jar` (absolute path to
  the Jalview jar; required), `input` (alignment input — `.fa` /
  `.aln` / `.clustal` / `.stockholm` and friends Jalview reads
  natively; required), `output_basename` (filename stem the
  adapter pins as the Jalview output target; required, non-empty),
  `output_format` (string; default `"png"` — selectable from
  `"png"` / `"html"` / `"svg"` / `"fasta"` / `"clustal"` for the
  canonical headless output formats Jalview ships), `extra_args`;
  `prepare()` resolves `jar` and `input` against the case
  directory when relative, validates each file exists on disk,
  derives the output extension from `output_format` (png →
  `.png`, html → `.html`, svg → `.svg`, fasta → `.fasta`, clustal
  → `.aln`, default → use the format string itself as the
  extension), builds `java -jar <jar> -nodisplay -open <input>
  -<output_format> <basename>.<ext> [extras...]`; `collect()`
  walks for `<output_basename>*.png` (`Native`, "Jalview
  alignment image"), `<output_basename>*.svg` (`Native`, "Jalview
  SVG"), `<output_basename>*.html` (`Native`, "Jalview HTML"),
  `<output_basename>*.fasta` (`Native`, "Jalview FASTA"),
  `<output_basename>*.aln` (`Tabular`, "Jalview alignment"),
  `*.log`; probe via `find_on_path(&["java"])` — Jalview's
  version comes from the jar itself, not from `java`, so we
  surface no version here (same shape as Phase 33 j5 / Cello);
  GPL-3.0 licensed; version range `2.11.0..3.0.0`;
  `bio.jalview.view` ribbon capability). Each wired into
  `valenx-app::init_registry`.
- No new canonical types, no new format readers, no new CLIs.
- 2 new `valenx-init` templates ship: `pydna`, `jalview`. Cross-
  binary roundtrip test sweeps all 115 templates clean.

The full per-phase shape lives in
[docs/src/phases/phase-41-sequence-editors.md](./docs/src/phases/phase-41-sequence-editors.md);
the implementation plan lives at
[docs/superpowers/plans/2026-05-03-sequence-editors.md](./docs/superpowers/plans/2026-05-03-sequence-editors.md).

**Deliverable:** Adapter inventory at 119 of 120 fully live after
this phase trio (2 new sequence-editor adapters added on top of
the prior set, alongside the Phase 32.5 spatial-stochastic pair
and Phase 40 microscopy trio); Phase 41 opens the first plasmid-
design / alignment-viewer domain in Valenx and crosses the
**100-bio-adapter milestone** with the bio-adapter count now
totaling **100 across 36 biology / biotech / chemistry phases**.

### Phase 42 — Web visualization · live ✅
Open the **first modern web 3D molecular visualization domain**
in Valenx with two established open-source WebGL viewers that span
the web-visualization tradeoff space — the canonical PDBe / RCSB
modern viewer that powers the structural-biology web (Mol*, the
EMBL-EBI / RCSB-led MIT-licensed WebGL toolkit that has become the
de-facto modern molecular viewer embedded in the PDB / PDBe /
AlphaFold DB / ESM Atlas web properties since the late 2010s), and
the Rose lab's WebGL framework that predated Mol* and still powers
a large fraction of the Jupyter-friendly notebook visualization
ecosystem (NGL Viewer, Alexander Rose's MIT-licensed high-
performance WebGL framework). Both are JavaScript browser libraries
in their primary distribution form; we wrap them via their
**Python bindings** (`molstar` / `nglview`, the official PyPI-
distributed Python interfaces both projects maintain) so they slot
into the existing Python-script subprocess pattern (sister to
Phase 17 Biopython, Phase 19.5 Scanpy, Phase 33 pySBOL, Phase 41
pydna). The user supplies a `.py` script that imports the binding
to compose a state file, rendered output, or notebook-embedded
view. Phase 42 sits numerically after Phase 41 sequence editors
and ships chronologically right after Phase 22.5 workflow
expansion — same chronological-vs-numerical convention used for
Phase 17.5 / 24 / 28 / 31 / 35 / 39 / 5.5 / 5.6 / 5.7 / 32.5 / 40
/ 41 / 22.5.

- 2 new adapter crates land under `crates/valenx-adapters/bio/`:
  Mol* (Python-script subprocess sister to Phase 17 Biopython /
  Phase 19.5 Scanpy / Phase 33 pySBOL / Phase 41 pydna; `case.toml`
  knobs `script` (path to user-supplied Python script; required,
  `.py` enforced) / `python` (interpreter name; default
  `"python3"`) / `input_structure` (`Option<PathBuf>` — optional
  `.pdb` / `.cif` / `.mmcif` structure file; `None` when the
  script fetches from the PDB / generates the structure inline) /
  `output_basename` (filename stem; required, non-empty);
  `prepare()` enforces `.py`, stages script + optional
  input_structure, writes `valenx_params.json` with
  `output_basename` always plus `input_structure` (staged
  filename) only when set — key omitted entirely when `None`
  rather than emitted as `null`, matching the hand-rolled JSON
  convention the rest of the bio adapters use, builds `<python>
  <staged_script>`; `collect()` walks for `<output_basename>*.html`
  (`Native`, "Mol* viewer HTML"), `<output_basename>*.molj`
  (`Native`, "Mol* state file"), `<output_basename>*.png`
  (`Native`, "Mol* rendered image"), `*.log`; probe via Python on
  PATH then `<python> -c "import molstar"` — on import failure
  surface as a `ProbeReport.warnings` entry, not error (sister to
  the Phase 19.5 scanpy / scvi / Phase 19.6 AnnData / Phase 5.6
  HOOMD-blue / Phase 5.7 MDTraj / Phase 41 pydna probe
  convention); MIT licensed; version range `3.0.0..5.0.0`;
  `bio.molstar.view` ribbon capability), and NGL Viewer (Python-
  script subprocess sister to Mol*; `case.toml` knobs `script`
  / `python` / `input_structure` / `output_basename` mirror
  Mol* exactly; `prepare()` mirrors Mol* shape exactly; `collect()`
  walks for `<output_basename>*.html` (`Native`, "NGL viewer
  HTML"), `<output_basename>*.png` (`Native`, "NGL rendered
  image"), `<output_basename>*.json` (`Tabular`, "NGL state
  JSON"), `*.log`; probe via Python on PATH then `<python> -c
  "import nglview"` — sister probe to Mol* with `nglview`-
  specific warning; MIT licensed; version range `3.0.0..5.0.0`;
  `bio.ngl.view` ribbon capability). Each wired into
  `valenx-app::init_registry`.
- No new canonical types, no new format readers, no new CLIs.
- 2 new `valenx-init` templates ship: `molstar`, `ngl`. Cross-
  binary roundtrip test sweeps all 120 templates clean.

The full per-phase shape lives in
[docs/src/phases/phase-42-web-visualization.md](./docs/src/phases/phase-42-web-visualization.md);
the implementation plan lives at
[docs/superpowers/plans/2026-05-04-web-visualization.md](./docs/superpowers/plans/2026-05-04-web-visualization.md).

**Deliverable:** Adapter inventory at **124 of 125 fully live**
after this phase pair (2 new web-visualization adapters added on
top of the prior set, alongside the Phase 22.5 workflow-expansion
trio); Phase 42 opens the first modern web 3D molecular
visualization domain in Valenx and **completes the bio ecosystem
from the user's original /review list** — every major category
called out is now covered, with the bio-adapter count totaling
**105 bio adapters across 38 biology / biotech / chemistry
phases**. **🎯 Bio ecosystem · COMPLETE ✅**.

### Phase 43 — mRNA design · live ✅
Open the **first mRNA / vaccine therapeutic design domain** in
Valenx with three established open-source tools that span the
codon-optimization + joint-design tradeoff space — the Edinburgh
Genome Foundry's general-purpose constraint-driven codon optimizer
that's the de-facto Python choice for synthetic-gene design (DNA
Chisel, the MIT-licensed library that handles codon optimization,
restriction-site avoidance, repeat scanning, GC-content tuning,
forbidden-pattern matching, and arbitrary user-defined constraints
composed into a single optimization objective), Baidu Research's
joint codon + secondary-structure mRNA design tool that landed as
the modern mRNA-vaccine design workhorse since the 2021 _Nature_
paper (LinearDesign, the Apache-2.0 single-binary CLI that jointly
optimizes codon usage and mRNA secondary-structure stability under
a tunable Lagrangian tradeoff parameter), and the Vejnar lab's
codon-level mRNA stability predictor (iCodon, the GPL-3.0 R-based
tool that scores per-position codon contributions to mRNA half-
life given a target organism). DNA Chisel + iCodon follow the
established Phase 17 Biopython + Phase 19.6 Seurat subprocess
patterns: the user supplies a `.py` / `.R` script that imports the
upstream package and reads `valenx_params.json` for the parsed
knobs. LinearDesign follows the established Phase 18 BWA single-
binary CLI pattern: protein FASTA in, optimized mRNA out. Phase 43
sits numerically after Phase 42 web visualization and ships
chronologically right after Phase 42 — same chronological-vs-
numerical convention used for Phase 17.5 / 24 / 28 / 31 / 35 / 39
/ 5.5 / 5.6 / 5.7 / 32.5 / 40 / 41 / 22.5 / 42.

- 3 new adapter crates land under `crates/valenx-adapters/bio/`:
  DNA Chisel (Python-script subprocess sister to Phase 17
  Biopython / Phase 19.5 Scanpy / Phase 33 pySBOL / Phase 41
  pydna / Phase 42 Mol* / NGL Viewer; `case.toml` knobs `script`
  (path to user-supplied Python script; required, `.py` enforced)
  / `python` (interpreter name; default `"python3"`) /
  `input_fasta` (`Option<PathBuf>` — optional starting `.fa` /
  `.fasta` FASTA; `None` when the script generates the sequence
  from scratch or fetches it inline) / `output_basename`
  (filename stem; required, non-empty); `prepare()` enforces
  `.py`, routes script + optional input_fasta through
  `confined_join` to stage them safely in the workdir, writes
  `valenx_params.json` with `output_basename` always plus
  `input_fasta` (staged filename) only when set — key omitted
  entirely when `None` rather than emitted as `null`, matching
  the hand-rolled JSON convention the rest of the bio adapters
  use, builds `<python> <staged_script>`; `collect()` walks for
  `<output_basename>*.fasta` (`Native`, "DNA Chisel optimized
  FASTA"), `<output_basename>*.gb` / `.genbank` (`Native`, "DNA
  Chisel GenBank"), `<output_basename>*.json` (`Tabular`, "DNA
  Chisel constraint report" — the canonical machine-readable
  per-constraint-pass / per-constraint-fail report DNA Chisel
  emits for downstream automation), `<output_basename>*.png`
  (`Native`, "DNA Chisel plot"), `*.log`; probe via Python on
  PATH then `<python> -c "import dnachisel"` — on import failure
  surface as a `ProbeReport.warnings` entry, not error (sister
  to the Phase 19.5 scanpy / scvi / Phase 19.6 AnnData / Phase
  5.6 HOOMD-blue / Phase 5.7 MDTraj / Phase 41 pydna / Phase 42
  Mol* / NGL probe convention); MIT licensed; version range
  `3.0.0..4.0.0`; `bio.dnachisel.optimize` ribbon capability),
  LinearDesign (single-binary CLI subprocess sister to Phase 18
  BWA / Phase 32.5 Smoldyn / Phase 5 GROMACS; `case.toml` knobs
  `protein` (path to protein FASTA; required — read in place,
  no staging) / `output_basename` / `lambda_param` (`f64`,
  finite and ≥ 0.0; default 1.0; the Rust field is
  `lambda_param` because `lambda` is a Rust reserved keyword —
  the CLI emits `--lambda <value>` regardless; tunable
  Lagrangian tradeoff between codon-adaptation-index and
  predicted mRNA secondary-structure stability — `0.0` = pure
  MFE-optimal, large = pure CAI-optimal, intermediate values
  like default `1.0` hit the joint sweet spot demonstrated in
  the paper) / `codon_usage` (target organism codon-usage table
  name; default `"human"` — selectable from the LinearDesign-
  shipped set: `"human"` / `"mouse"` / `"yeast"` / `"ecoli"` /
  etc.) / `extra_args`; `prepare()` validates `lambda_param`
  is finite and ≥ 0.0 (returns `InvalidCase` when negative or
  NaN), resolves `protein` against the case directory when
  relative, validates the file exists on disk (read in place,
  no staging — same shape as Phase 18 BWA's reference genome),
  builds `lineardesign --aa <protein> --lambda <lambda_param>
  --codon_usage <codon_usage> --output_basename <basename>
  [extras...]`; `collect()` walks for `<output_basename>*.fasta`
  (`Native`, "LinearDesign optimized mRNA"),
  `<output_basename>*.txt` (`Tabular`, "LinearDesign report" —
  the canonical per-design summary report LinearDesign writes
  alongside the FASTA output), `*.log`; probe via
  `find_on_path(&["lineardesign"])` — when the `lineardesign`
  binary isn't found but Python is on PATH the probe surfaces a
  targeted `"clone https://github.com/LinearDesignSoftware/
  LinearDesign and add the bin directory to PATH"` warning;
  Apache-2.0 licensed; version range `1.0.0..2.0.0`;
  `bio.lineardesign.design` ribbon capability), and iCodon
  (Rscript subprocess pattern sister to Phase 19.6 Seurat;
  `case.toml` knobs `script` (path to user-supplied R script;
  required, `.R` enforced) / `rscript` (R interpreter name;
  default `"Rscript"`) / `input_fasta` (`Option<PathBuf>` —
  optional input mRNA FASTA the script can score; `None` when
  the script generates the sequence inline or reads it from a
  different source) / `output_basename` (filename stem;
  required, non-empty); `prepare()` enforces the `.R` extension,
  routes script + optional input_fasta through `confined_join`
  to stage them safely in the workdir, writes
  `valenx_params.json` with the same hand-rolled JSON shape as
  DNA Chisel (key omitted when `None`); `collect()` walks for
  `<output_basename>*.csv` / `*.tsv` (`Tabular`, "iCodon
  stability table"), `<output_basename>*.rds` (`Native`,
  "iCodon R object (RDS)" — canonical R-serialised iCodon model
  output consumed by every downstream R-side stability /
  visualization pipeline), `<output_basename>*.png` (`Native`,
  "iCodon plot"), `*.log`; probe via
  `find_on_path(&["Rscript"])` — does not attempt to confirm
  iCodon itself is installed because that would require running
  R, an expensive multi-second startup at probe time (same
  shape as Phase 19.6 Seurat); the `ToolNotInstalled` install
  hint mentions the canonical
  `devtools::install_github('santiago1234/iCodon')` install
  path; GPL-3.0 licensed; version range `1.0.0..2.0.0`;
  `bio.icodon.predict` ribbon capability). Each wired into
  `valenx-app::init_registry`.
- No new canonical types, no new format readers, no new CLIs.
- 3 new `valenx-init` templates ship: `dnachisel`,
  `lineardesign`, `icodon`. Cross-binary roundtrip test sweeps
  all 123 templates clean.

The full per-phase shape lives in
[docs/src/phases/phase-43-mrna-design.md](./docs/src/phases/phase-43-mrna-design.md);
the implementation plan lives at
[docs/superpowers/plans/2026-05-04-mrna-design.md](./docs/superpowers/plans/2026-05-04-mrna-design.md).

**Deliverable:** Adapter inventory at **127 fully live** after
this phase trio (3 new mRNA-design adapters added on top of the
bio-ecosystem-complete milestone reached at Phase 22.5 + 42);
Phase 43 opens the first mRNA / vaccine therapeutic design domain
in Valenx and layers it on top of the bio-ecosystem-complete
milestone, with the bio-adapter count totaling **108 bio adapters
across 39 biology / biotech / chemistry phases**.

### Phase 44.5 — RNA folding expansion · live ✅
Sister-adapter expansion of the existing Phase 28 RNA structure
trio (ViennaRNA / RNAstructure / NUPACK). Round out the RNA
secondary-structure folding surface with three more canonical
folders that span the modern tradeoff space — **mfold/UNAFold**
(Michael Zuker's classic dynamic-programming Zuker / Stiegler RNA
folder; academic-license; minimum-free-energy folding plus
suboptimal-structure ensembles via the canonical Turner / Mathews
thermodynamic parameters; single-binary subprocess shape sister to
Phase 18 BWA with mfold's `KEY=VALUE`-style invocation `mfold SEQ=
<sequence> NA=RNA T=<temperature>`; surfaces `"academic"` /
`"non-commercial"`-keyworded license-awareness warning sister to
ViennaRNA / NUPACK / VMD / NAMD), **EternaFold** (the Eterna
project's MIT-licensed ML-aware folder via the Das lab's `arnie`
Python wrapper; trained on a half-decade of crowd-sourced Eterna
gameplay puzzles plus thermodynamic + ML corpora; Python-script
subprocess shape sister to Phase 17 Biopython / Phase 28 NUPACK),
**LinearFold** (Baidu / Oregon State's Apache-2.0 beam-search
linear-time folder — folding-only sister to Phase 43 LinearDesign
with the same beam-search core from the same group applied to the
inverse problem; linear `O(N · beam_size)` complexity scales to
viral-genome-length sequences without the cubic blowup of classical
DP folders; non-standard stdin contract — sequence on stdin,
structure on stdout). Each wired into `valenx-app::init_registry`.
3 new `valenx-init` templates (`mfold`, `eternafold`, `linearfold`).

### Phase 35.5 — Base + prime editing design · live ✅
Sister-adapter expansion of the existing Phase 35 CRISPR design
trio (CHOPCHOP / CRISPOR / Cas-OFFinder). Round out the CRISPR
guide-RNA design surface with four canonical non-cleavage editing
tools the Phase 35 Cas9-cut-focused adapters don't cover — **BE-
Designer** (Komor lab base-editor guide design, MIT — de-facto
first stop for "I want to make a C→T or A→G base change"),
**BE-Hive** (Liu lab base-editing outcome predictor, MIT — the
canonical Python module is `be_predict`), **PrimeDesign** (Liu
lab pegRNA designer for the Anzalone / Liu prime-editing system,
MIT), **pegFinder** (Komor lab alternative pegRNA finder, MIT —
sister to PrimeDesign with a different scoring model emphasising
pegRNA secondary-structure stability + RT-template-length
tradeoffs). All four ride the established Python-script subprocess
pattern (sister to Phase 17 Biopython, Phase 35 CHOPCHOP /
CRISPOR, Phase 41 pydna, Phase 43 DNA Chisel / iCodon). Each
wired into `valenx-app::init_registry`. 4 new `valenx-init`
templates (`be-designer`, `be-hive`, `primedesign`, `pegfinder`).

### Phase 35.6 — Edit-outcome prediction · live ✅
Sister-adapter expansion of the existing Phase 35 + 35.5 CRISPR
design + editing surface. Close the design → predict-outcome →
off-target loop with four canonical outcome predictors —
**inDelphi** (Liu lab's MIT-licensed Cas9-cut indel pattern
predictor), **FORECasT** (Sanger Institute's Apache-2.0
alternative indel predictor — the Python module is `selftarget`
named after Allen's data-collection assay rather than the
predictor's published name), **AlphaMissense** (DeepMind's
missense-effect predictor extending the AlphaFold lineage to
score per-position pathogenicity calibrated against ClinVar
pathogenic / benign labels; **CC-BY-NC-SA-4.0 / academic non-
commercial weights** — probe pushes a mandatory `"academic"` /
`"non-commercial"`-keyworded warning whenever Python is on PATH,
regardless of whether `import alphamissense` succeeds, sister to
AlphaFold 3's mandatory probe-warning pattern), **CRISPRitz**
(Pinello lab's MIT-licensed variant-aware off-target genome-wide
search — sister to Phase 35 Cas-OFFinder with the distinguishing
property of walking population VCFs (1000 Genomes, gnomAD) for
off-target sites that exist only in specific haplotypes). All
four ride the established Python-script subprocess pattern. Each
wired into `valenx-app::init_registry`. 4 new `valenx-init`
templates (`indelphi`, `forecast`, `alphamissense`, `crispritz`).

### Phase 45 — Pharmacokinetics + RNA tertiary structure · live ✅
Open **two new domains** in Valenx with two single-canonical-
adapter beachheads — **PK-Sim** (Open Systems Pharmacology
suite's GPL-2.0 physiologically-based PK / PBPK simulator, the
de-facto open-source PBPK modeling tool descended from the Bayer
internal pharmacokinetic simulator opened to the community via
the OSP Initiative; models whole-body drug ADME using a
physiologically-grounded compartmental representation; consumes a
`.pksim5` XML project file authored in PK-Sim's GUI or
programmatically through the OSP Python API; single-binary
subprocess shape sister to Phase 18 BWA with `pksim --project
<project> --output <output_basename> [extras...]`; collects per-
compartment concentration-time CSVs + simulation metadata JSON;
opens the **first PK/PD pharmacokinetics modeling category in
Valenx** distinct from Phase 32 / 32.5 cellular-scale ODE /
spatial-stochastic modeling), and **SimRNA** (Bujnicki group's
GPL-3.0 coarse-grained Monte Carlo RNA tertiary-structure
predictor — predicts the **full 3D Cartesian backbone** of an RNA
from its sequence using five-bead per-nucleotide coarse-graining
and replica-exchange Monte Carlo over the reduced-coordinate
energy landscape; single-binary subprocess shape sister to PK-Sim
/ Phase 18 BWA / Phase 32.5 Smoldyn with `SimRNA -c <config>
-s <sequence> -o <output_basename> -R <n_replicas> [extras...]`;
collects predicted PDB tertiary structures + replica-exchange
`.trafl` trajectories + per-step energy `.txt` logs; opens the
**first RNA tertiary 3D structure prediction category in Valenx**
distinct from Phase 28 + 44.5 RNA secondary 2D folders). Each
wired into `valenx-app::init_registry`. 2 new `valenx-init`
templates (`pksim`, `simrna`).

**Deliverable:** Adapter inventory at **141 fully live** after
the four-phase rollup (3 RNA-folding-expansion adapters + 4
base + prime editing design adapters + 4 edit-outcome-prediction
adapters + 2 pharmacokinetics + RNA-tertiary adapters added on
top of the Phase 43 mRNA design beachhead and the bio-ecosystem-
complete milestone reached at Phase 22.5 + 42); the bio-adapter
count now totals **121 bio adapters across 43 biology / biotech
/ chemistry phases**, covering every category from the user's
original /review list plus the four new sister-domain expansions
and two new domain beachheads.

### Future phases — 18.5 → 45
Phases 17.5 → 43 cover the remaining ~190 biology / biotech tools from
the user's spec (sequence editors, alignment + search, variant calling,
workflow managers, additional MD engines, DNA origami design, protein
design, RNA structure, population genetics, phylogenetics, whole-cell
systems biology, synthetic biology, molecular docking, CRISPR design,
quantum chemistry, cryo-EM reconstruction, ML / AI foundation models,
Rosetta family, and more). The full table — phase number, title, adapter
list, canonical-type extensions, estimated size — lives at the bottom of
[docs/superpowers/plans/2026-04-30-biology-foundation.md](./docs/superpowers/plans/2026-04-30-biology-foundation.md).
Cumulative scope after all phases: ~210 adapters, ~1300 tasks ≈ 650 hours
of focused implementation work, on top of the existing physics-domain
phases.

---

## 3. Staffing and funding

### Honest team sizing (to full 20-year feature depth)

| Team | Year to Phase 1 (v0.1) | Year to Phase 5 (native CFD) | Year to Phase 11 (verticals mature) | Year to Phase 15 (infrastructure) |
|---|---|---|---|---|
| 1 dev (solo) | 2 | 10+ | unlikely | unlikely |
| 2-3 contributors | 1 | 5-7 | 12-15 | 20 |
| 5 full-time | 0.75 | 3-4 | 8-10 | 15-18 |
| 10 full-time + community | 0.5 | 2-3 | 6-8 | 12-15 |
| 25+ core + broader OSS community | 0.25 | 1-2 | 5 | 10-12 |

**Takeaway:** a two-decade roadmap is realistic for a team of 5-10 by
year 15-20. Solo or 2-3 person teams can reach v0.1 and native-CFD but
won't cover the full vertical depth without more contributors arriving.

### Funding stages

| Stage | Years | Source | Amount |
|---|---|---|---|
| 1 | 0-2 | Self-fund + GitHub Sponsors + small grants | ~$50-200K/yr |
| 2 | 2-5 | **NSF POSE** ($500K-1M), **EU Horizon**, **DOE Office of Science**, **DARPA** | $500K-2M/yr |
| 3 | 4-10 | **Industry consortium** — automotive, aero, energy members at $50-100K/yr | $300K-1.5M/yr |
| 4 | 6-15 | **Commercial support** (Red Hat / Databricks model — paid support for companies using free software) | $500K-5M/yr |
| 5 | 10-20 | **Foundation endowment** (Apache, Eclipse, or purpose-built Valenx Foundation) | $2-10M/yr |
| 6 | 15-20+ | **Long-term infrastructure funding** — NIH R24, institutional core, federal science infrastructure | $1-5M/yr sustaining |

### Total cost to 20-year maturity

Cumulative: roughly **$50M-$150M** over 20 years depending on team size
and vertical depth. Competing ANSYS alone spends $800M/yr on R&D — we
are 1-2% of their budget spread across two decades, competing through
focus + community contribution + the fact that we don't pay license-fee
overhead.

### Sustainability model

By Year 10-15, the project should be self-sustaining via:
- Industrial consortium dues (voting rights + early-access features)
- Commercial support contracts for large deployments
- Endowed university positions (core maintainers on tenure track)
- Foundation governance and ops (similar to how Python Software Foundation funds CPython)
- Marketplace revenue share on paid plugins (split with plugin authors)

---

## 4. Risks

| Risk | Likelihood | Mitigation |
|---|---|---|
| Solver physics bugs destroy credibility | High | Validation-first: every solver ships with published-benchmark regression tests before "ready" label |
| Scope creep — rewrite everything at once | Near certain | Ruthless prioritization; adapters keep app usable while natives mature |
| Lone-maintainer burnout | High (years 1-5) | Document architecture early; recruit 2-3 co-maintainers in Year 1; RFC-based governance |
| Low contributor influx to Rust systems code | Medium | Rust community is strong; good RFC docs + mentoring; target domain experts willing to learn Rust |
| GPL license contagion accidentally | Medium | Every upstream dep reviewed; adapters are subprocess-only crates with clear licensing boundaries |
| Industry won't switch from ANSYS without validation | High (years 1-4) | Position as complementary early; publish validation papers; target universities/startups first |
| Grant funding doesn't materialize | Medium | Keep burn rate low (1-2 FTE) until Phase 1 ships; grants look much better with a real product |
| Rust numerical ecosystem gaps | Low-Medium | `nalgebra`, `ndarray`, `faer-rs`, `rsmpi` are mature; BLAS/LAPACK via bindings for hot paths |
| **Founder departure before Year 10** | High over 20yr | Governance transition to TSC by Year 3; no-bus-factor-1 rule; endowed positions at universities by Year 8 |
| **Rust loses popularity over 20 years** | Low-Medium | Rust is in Linux kernel, adopted broadly; migration path to any future systems language is mechanical |
| **UI framework (egui) abandoned** | Low | Rust UI ecosystem has multiple mature alternatives (Slint, Iced); internal abstraction layer limits blast radius |
| **Upstream OSS tool stagnates (e.g., FreeCAD development slows)** | Medium | We fork + maintain + contribute back; adapters are replaceable |
| **Upstream OSS tool license changes** | Low-Medium | We vendor specific versions; LGPL/Apache permissive licenses don't typically tighten |
| **ANSYS / Siemens aggressive legal response** | Low | Prior art is on our side; OSS alternative is protected speech; we don't copy proprietary algorithms |
| **Commercial support competition undermines sustainability** | Medium | Federated support model — multiple vendors offer support; foundation prevents capture |
| **National regulatory changes for export-controlled physics** (e.g., nuclear, ITAR) | Medium | Keep export-controlled physics OUT of core; provide interfaces for users to bring their own |
| **Long-tail feature coverage stagnates** | Medium | Plugin marketplace incentivizes community; vertical sponsors fund their own capability gaps |
| **Bus factor of 1 for key subsystems** | High initially | Mentorship program; every subsystem needs 2+ maintainers by Year 5 |
| **Security incidents** (since we ship binaries) | Medium over 20yr | Signed releases; SLSA provenance; reproducible builds; security@valenx address from Day 1 |
| **Competing OSS project emerges and fragments community** | Low-Medium | Welcome collaboration; merge if scope aligns; differentiate on focus + polish |

---

## 5. Success metrics

### Year 1
- `valenx 0.1` shipped (Phase 1 complete)
- 500 GitHub stars; 100 monthly active users
- All legacy adapters ported to Rust
- Native app feels like Fusion 360 in polish

### Year 3
- Phase 2 (geometry) + Phase 3 (meshing) native
- Phase 4 (native CFD basics) underway
- 5000 stars; 1000 MAU; 3 companies using in production
- First NSF POSE grant secured

### Year 5
- Native CFD graded A/B on 10+ published benchmarks (Phase 4 done)
- Native FEA linear + nonlinear static shipped (Phase 6 partial)
- 20000 stars; industry consortium formed with 5+ members
- First peer-reviewed validation paper published
- First vertical (aerospace) fully matured

### Year 8
- Phase 5 (advanced CFD: compressible, multiphase, rotating) shipped
- Phase 7 (advanced FEA) shipped
- Phase 8 (EM, chemistry, battery, MD natives) ~50% done
- 50000 stars; consortium at 10+ members
- 5 companies using in production daily
- Teaching universities adopt Valenx for simulation courses

### Year 10
- All core physics natives at 0.5+ capability
- Phase 9 (optimization + UQ) shipped
- Phase 11 (verticals) — 3 verticals (aero, automotive, electronics cooling) at production-grade
- Self-sustaining via grants + consortium + commercial support
- 100K+ installs; 20+ peer-reviewed Valenx-based publications

### Year 15
- Feature parity with ANSYS Workbench in 5+ verticals
- Phase 13 (ML/AI integration) mature
- Phase 12 (cloud + HPC) production-grade; deployed at multiple national labs
- 250K+ installs; 500K+ downloads/yr
- Plugin marketplace with 500+ plugins
- Industry certifications in place where relevant (FDA medical, ISO structural)
- Foundation governance in place; project outlives any single maintainer

### Year 20
- Known as "the FreeCAD of everything" — the obvious OSS default across
  scientific computing
- Feature depth rivals ANSYS / Siemens NX / COMSOL for most use cases
- Endowed development positions at 3+ partner universities
- 1M+ installs; 2000+ papers citing Valenx in methods
- Marketplace with 1000s of plugins; thriving contributor ecosystem
- Reference implementation for reproducible simulation (runs 20-year-old
  `.valenx` cases identically)
- Phase 15 complete — the project is infrastructure

---

## 6. Governance

### Model evolution over 20 years

| Years | Governance |
|---|---|
| 0-2 | **BDFL** (you). Fast decisions. RFC process for major changes. |
| 2-5 | **BDFL + Technical Steering Committee** — elected by contributors. BDFL retains tie-break vote. |
| 5-10 | **TSC-only** — BDFL steps down from tie-breaks. Foundation affiliate (Apache or similar) in progress. |
| 10-20 | **Full foundation governance** — Valenx Foundation (or Apache/Eclipse home). Board elected by members. Development funded by consortium + grants + endowment. |

### Core operating rules

- **License:** Apache 2.0 (permissive, industry-friendly, enables
  commercial adoption without GPL stigma). Bundled tools retain their
  original licenses with clear notices.
- **RFCs** follow Rust's format — markdown docs in
  [`rfcs/`](./rfcs/), discussion in a PR, merge = accepted.
  Required for any user-facing or API change.
  See [rfcs/README.md](./rfcs/README.md) for the full process.
- **Versioning:** SemVer; 0.x until native CFD ships; 1.0 = Phase 5
  complete (~Year 6). 2.0 = Phase 11 (verticals mature, ~Year 12).
  Full SemVer contract in [POLICIES.md](./POLICIES.md).
- **Release cadence:** nightly from `main`; quarterly minor, patch as
  needed on Stable; LTS every ~18 months, supported 24 months.
- **Code of Conduct:** Contributor Covenant v2.1 — see
  [CODE_OF_CONDUCT.md](./CODE_OF_CONDUCT.md).
- **Contributor workflow:** see [CONTRIBUTING.md](./CONTRIBUTING.md).
- **Backward compatibility:** `.valenx` project files forward-compatible
  forever (one-time migrate on major bumps); APIs deprecated over
  two minor versions before removal — details in
  [POLICIES.md](./POLICIES.md).
- **Reproducibility guarantee:** a case from version X must run
  identically in version Y for at least 10 years, given the same
  pinned `tools.lock` (Section 9 and [RFC 0001](./rfcs/0001-project-file-format.md)).
- **Bus-factor rule:** no single maintainer owns a subsystem alone
  past Year 5 — at least 2 maintainers per crate.
- **Security:** signed releases from Day 1; SLSA provenance;
  reproducible builds; published security policy with response SLA
  in [SECURITY.md](./SECURITY.md).

---

## 7. Decisions already made

- **Primary language:** Rust
- **Secondary:** C (for FFI to OCCT, VTK, BLAS), Python (embedded via PyO3 for scripting)
- **UI framework:** egui (recommendation; alternative Slint if declarative feel matters more)
- **3D:** wgpu
- **CAD kernel:** OpenCASCADE via Rust FFI — **wrap forever**, not rewrite. Pure-Rust alternatives (fornjot, truck) are upstream contributions, not a goal.
- **Physics solvers:** wrap existing validated OSS tools (OpenFOAM, Code_Aster, etc.) — **wrap forever, not rewrite**. Native Rust solvers are optional, added only where a clear reason exists.
- **License:** Apache 2.0 for Valenx itself. Bundled tools retain their licenses (BSD/LGPL/GPL handled per integration mode in Section 0.6).
- **Distribution:** GitHub Releases; per-OS native installers; no SaaS, no subscription.
- **UX reference:** Fusion 360. Not FreeCAD.
- **Version policy:** every Valenx release pins exact tool versions via a lock manifest; reproducibility guaranteed 10+ years (Section 9).
- **Install experience:** small installer + first-run tool-picker wizard; modules downloaded on demand (Section 10).

## 8. Decisions pending

1. Project name: keep Valenx vs. rebrand?
2. Should macOS ship from Day 1 or Windows-first?
3. Should we support Python scripting from Day 1 (PyO3) or wait for Phase 2?
4. Does the v0.1 include any native Rust physics (e.g., valenx-cfd linear solver) or pure adapters only? (Recommend pure adapters for v0.1.)
5. Do we run a hosted documentation site + community forum, or GitHub Discussions only?
6. egui vs. Slint — final pick (recommend egui, but Slint's declarative approach may age better over 20 years)
7. Single mono-installer vs. per-platform variants (MSI/dmg/deb/AppImage all at once, or stagger)
8. When we hit Phase 5, name the Valenx 1.0 release after what milestone? (Native CFD validation? Vertical launch?)

---

## 9. Tool versioning strategy

**Every Valenx release pins every integrated tool to an exact version** in
a lock manifest (`tools.lock`) shipped inside the release. Tools are
upgraded in new Valenx releases after validation — never silently under
users' feet.

### Release channels
| Channel | Tool updates | Audience |
|---|---|---|
| **Stable** | Quarterly, tested | Most users |
| **LTS** | Frozen 2 years, security backports only | Enterprises, teaching |
| **Nightly** | Tool bumps as soon as they pass CI | Early adopters |

### Per-tool update policy
- **OpenFOAM / Code_Aster / Code_Saturne (major breaking releases)**: pin per Valenx major version; dual-support two majors during transition
- **OCCT / Cantera / gmsh (minor semver)**: pin minor, take patches freely
- **CalculiX / openEMS (infrequent)**: pin exact version
- **Numerical libs (PETSc/SUNDIALS/hypre/Kokkos)**: pin minor, take patches

### Reproducibility guarantee
A `.valenx` case file records which tool versions were used. Opening in
Valenx 10 years later uses those exact versions (cached on release
server). A case run in 2026 runs identically in 2036+.

### User overrides
Set `VALENX_{TOOL}_PATH=/custom/path` in settings to use your own install
instead of the bundled version. Useful for HPC, custom-patched tools,
pre-release testing.

### Maintainer workflow
1. Upstream tool releases new version
2. CI branch created with manifest bump
3. Full validation benchmark suite runs
4. If no grade regression → merge; ship in next Valenx stable
5. If regression → investigate + patch upstream + hold bump

---

## 10. Install experience

**Design goals:** Fusion-360 ease-of-install. No Python packages, no
`apt install`, no PATH editing, no `.bashrc` tweaks. Download → double-click
→ running in under a minute. First-run wizard handles physics-module
selection and per-user download.

### Platform installers
| Platform | Primary | Alternatives |
|---|---|---|
| **Windows** | `.exe` NSIS (signed) | `.msi`, winget, Chocolatey |
| **macOS** | `.dmg` (notarized — drag to Applications) | Homebrew cask, `.pkg` |
| **Linux** | `.AppImage` (no install, just chmod+x + run) | `.deb`, `.rpm`, Flatpak, Snap, AUR |

### Installer size
- **Installer itself:** ~250 MB (shell + core + CAD + viewer + embedded Python + essential numerical libs)
- **Full first-run download:** +1.5–2.5 GB for typical physics selection
- **Offline bundle** (separate download, ~3.5 GB): everything pre-included for air-gapped environments

### First-run setup wizard
On first launch:
1. Welcome screen
2. Tool-picker checklist by domain (CFD / FEA / Chemistry / Battery / EM / MD / etc.)
3. Default selection: "Recommended" (CFD + FEA + Chemistry + Battery ≈ 2 GB)
4. Alternative choices: "Core only" (skip all), "Everything" (all ~3 GB), "Advanced: use my own tools"
5. Progress bar with per-tool download + integrity check (SHA256)
6. When done: Valenx main window opens

### Post-install: Tool Manager
`Settings → Tools` shows every supported tool, its install status, disk
usage, and version. Users can install, uninstall, update, or point at a
custom install path.

### Updates
- Valenx checks for updates weekly on launch (configurable)
- When a new Valenx release drops, installer patches the shell + bumps
  all pinned tool versions in the manifest (downloads only what changed)
- User's projects and preferences preserved across updates

### Signing
- Signed via a CA code-signing certificate (DigiCert, Sectigo — budget
  $300-500/yr) to avoid Windows SmartScreen warnings
- macOS: Apple Developer certificate for notarization ($99/yr)
- Linux: GPG-signed `.deb` / `.rpm` / checksums for AppImage

---

## 11. What we ship this month

**Already landed** (repository foundation, Phase 0):

- Cleaned up the old web-app codebase; preserved reference material in
  [legacy-reference/](./legacy-reference/)
- Planning + architecture docs: [ROADMAP.md](./ROADMAP.md),
  [ARCHITECTURE.md](./ARCHITECTURE.md), [DESIGN.md](./DESIGN.md),
  [DESIGN_PRINCIPLES.md](./DESIGN_PRINCIPLES.md),
  [LANGUAGES.md](./LANGUAGES.md), [TESTING.md](./TESTING.md)
- Governance scaffold: [CONTRIBUTING.md](./CONTRIBUTING.md),
  [MAINTAINERS.md](./MAINTAINERS.md),
  [CODE_OF_CONDUCT.md](./CODE_OF_CONDUCT.md),
  [SECURITY.md](./SECURITY.md), [POLICIES.md](./POLICIES.md)
- RFCs: [0001 project file format](./rfcs/0001-project-file-format.md),
  [0002 adapter contract](./rfcs/0002-adapter-contract.md),
  [0003 plugin API](./rfcs/0003-plugin-api.md),
  [0004 results + fields](./rfcs/0004-results-and-fields.md),
  [0005 design principles](./rfcs/0005-design-principles.md),
  [0006 token system](./rfcs/0006-token-system.md)
- CI skeleton ([`.github/workflows/ci.yml`](./.github/workflows/ci.yml))
  + cargo-deny policy ([`deny.toml`](./deny.toml)) + issue / PR templates
  (now with a UI-PR checklist)
- Toolchain + editor configs: `rust-toolchain.toml`, `rustfmt.toml`,
  `.editorconfig`
- **Rust workspace scaffold** — root `Cargo.toml` plus stub crates:
  `valenx-app`, `valenx-core`, `valenx-geo`, `valenx-mesh`,
  `valenx-fields`, `valenx-viz`, `valenx-icons`, `valenx-fonts`,
  `valenx-design-tokens` (with real token JSON + build-time Rust
  codegen), and a first adapter crate stub
  `valenx-adapter-openfoam`
- **mdBook skeleton** at `docs/` with `book.toml`, `src/SUMMARY.md`,
  intro chapter
- **Design sources** at `docs/design/` — README, icon inventory,
  mockups and patterns placeholders
- [CHANGELOG.md](./CHANGELOG.md)

**Next** (Phase 0 → Phase 1 handoff):

1. ✅ Clean-slate Rust workspace scaffolded at `crates/`
2. ✅ Canonical types landed (`valenx-fields`, `valenx-geo`,
   `valenx-mesh`, `valenx-core`) with unit tests
3. ✅ `valenx-core::Adapter` trait + registry matches
   [RFC 0002](./rfcs/0002-adapter-contract.md)
4. ✅ `.valenx` project loader + writer
   (`valenx-core::project`) matching RFC 0001 — load, validate,
   round-trip, save atomically
5. ✅ Workflow DAG (`valenx-core::workflow`) — typed ports,
   cycle detection, topological ordering
6. ◑ `valenx-adapter-openfoam` probe implemented; prepare / run /
   collect return honest not-implemented errors — next step is real
   `prepare()` against an OpenFOAM case directory template
7. ☐ `valenx-app` opens a native egui window with the Workspace
   shell layout from [DESIGN.md § 6](./DESIGN.md)
8. ☐ `valenx-viz` can load and display an STL with the ViewCube
9. ☐ End-to-end: click "New Airfoil Case" → pick NACA → run
   simpleFoam → plot residuals in the native window — no web
   browser anywhere

Everything beyond that follows the phases. 20 years is a long time,
but we build it one phase at a time — and every phase ships something
usable. Nothing waits 20 years for its first release.

---

## Appendix A — 20-year Gantt summary

```
Year      0────2────4────6────8────10───12───14───16───18───20
Phase 0   ██
Phase 1   ████
Phase 2    ██████████
Phase 3     ██████████
Phase 4       ████████████
Phase 5         ████████████████
Phase 6        ████████████████
Phase 7             ████████████████
Phase 8          ████████████████████████
Phase 9                   ████████████████████████
Phase 10                      ████████████████████████
Phase 11                         ████████████████████████████
Phase 12                             ████████████████████████████
Phase 13                                  ████████████████████████████
Phase 14                                       ████████████████████████
Phase 15                                             ████████████████████+
```

Phases overlap heavily — every phase beyond 0/1 runs in parallel with
several others. By Year 10 we're running 6+ phases simultaneously on
different parts of the codebase, each with their own maintainer group.

## Appendix B — What each year looks like from the user's perspective

- **Year 1:** "Oh nice, a native Valenx. Feels like Fusion 360. Wraps
  OpenFOAM cleanly. I can run CFD without opening a terminal."
- **Year 3:** "I can now do CAD natively without switching to FreeCAD.
  The meshing UI is the best I've used."
- **Year 5:** "Valenx has its own CFD solver now. It's validated against
  Ghia and Driver-Seegmiller. We're using it in production."
- **Year 8:** "All of our FEA cases run in Valenx now. Native EM solver
  is up. Chemistry and battery are starting to replace Cantera/PyBaMM."
- **Year 10:** "We don't even install OpenFOAM separately anymore. Valenx
  handles the full incompressible + compressible + multiphase + combustion
  stack natively. Optimization workflows are built in."
- **Year 15:** "Our university runs all sim courses on Valenx. AI-assisted
  case setup is like having a grad student who read every paper. Cloud
  dispatch to AWS is one click. Cheaper than ANSYS by 99%."
- **Year 20:** "Valenx is what you reach for when you need to simulate
  something. Like NumPy is for arrays. It's just there."
