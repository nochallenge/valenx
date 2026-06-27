# valenx Molecular Simulation + Visualization — build queue

**Date:** 2026-06-24
**Source:** user request — survey https://github.com/topics/molecular-dynamics?l=rust, add the worthwhile
Rust MD/viewer projects. Two named priorities: `David-OConnor/molchanica` (viewer) and
`lumol-org/lumol` (MD engine).

## License + gap analysis (honest)

valenx **already has `valenx-md`** (classical MD, force fields, thermostats) and a wgpu 3-D viewport +
`genetics::molecule_view`. So the rule is **fill gaps / add new capability**, not redundantly rewrite.
License-gated: permissive → port the code with attribution in `THIRD-PARTY-NOTICES`; GPL → clean-room
from published methods only.

| Project | Stars | License | Verdict |
|---|---|---|---|
| `lumol-org/lumol` | 209 | **BSD-3 ✓** | port the gaps: Ewald/Wolf electrostatics, Metropolis MC, NPT |
| `David-OConnor/molchanica` | 144 | **MIT ✓** | port viewer representations: cartoon/ribbon, SAS surface, density isosurface |
| `yesint/molar` | 54 | (verify) | analysis lib — optional, later |
| `seatonullberg/velvet` | 11 | (verify) | MD engine — overlaps valenx-md, **skip** |
| `ma3ke/molly`, `dnlbauer/xdrfile` | ~8 | (verify) | GROMACS XTC/TRR trajectory I/O — worth a pure-Rust reader |
| `caltechmsc/dreid-*`, `ForblazeProject/uff-relax` | ~8 | (verify) | DREIDING/UFF force-field typers — optional |
| `caltechmsc/cheq` | 8 | (verify) | QEq partial charges — optional, pairs with electrostatics |

## The tracks

### V1 — Molecular viewer representations (from molchanica, MIT)  ← user priority, IN PROGRESS
Enhance valenx's molecule viewer (in `valenx-app`, reusing the existing wgpu mesh rendering; no new
crate, reimplement marching-cubes — no dep): **ball-and-stick/sticks**, **cartoon/ribbon** (spline
through the Cα backbone), **SAS molecular surface** (marching-cubes on a union-of-spheres / Gaussian
density), and **density isosurface** (marching-cubes on a scalar grid, adjustable iso-level). Reactive
representation picker + AI-drivable (named widget / agent command), per the standing gate.

### V2 — valenx-md long-range electrostatics + Monte Carlo (from lumol, BSD-3)  ← IN PROGRESS
Extend `valenx-md` (no new dep) with the gaps lumol fills: **Ewald summation** (real + reciprocal +
self energy/forces), the **Wolf** damped-shifted-force method, and a **Metropolis Monte Carlo (NVT)**
driver; an **NPT barostat** only if genuinely absent. Pinned to the **NaCl Madelung constant
(≈1.747565)** for Ewald — the standard validation.

### V3 — MD trajectory I/O (from molly/xdrfile, license-permitting)
A pure-Rust **XTC/TRR trajectory reader** (the GROMACS format is documented; `molly`/`chemfiles` as
reference) so valenx can **load + play back MD trajectories** in the new viewer. License-verify the
references; reimplement the format (don't bind the LGPL C `xdrfile`).

### Optional / later (verify licenses first)
UFF/DREIDING force-field typers (assign FF params to an arbitrary molecule — feeds valenx-md), QEq
partial charges (`cheq`), and GROMACS-analysis utilities (order parameters, etc.).

## Loop integration
V1 (`valenx-app`) and V2 (`valenx-md`) are disjoint from the connective-layer + defense crates, so they
run concurrently (no new dep → no `Cargo.lock` contention with the one new-crate slot). V3 and the
optional tracks follow as slots free. Each: build subagent (TDD + scoped-gated green) → commit
email-safe (real files only, leak 0) → review to zero findings → next. GitHub HELD (local).
