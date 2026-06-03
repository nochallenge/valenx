# Phase 25 — Quantum chemistry

**Status:** 🟢 Live — Psi4 + NWChem + xTB open the **first quantum-
chemistry domain** in Valenx alongside the Phase 17 / 17.5 / 27 /
27.5 / 28 / 30 biology + structure-prediction + protein-design + RNA-
structure + phylogenetics beachheads and the Phase 24 cheminformatics
expansion.

## Goal

Open the quantum chemistry domain in Valenx with three established
open-source tools: **Psi4** (HF/DFT/post-HF general-purpose ab initio
quantum chemistry — Psithon-scriptable input, LGPL-3.0), **NWChem**
(massively-parallel ab initio quantum chemistry + plane-wave DFT — its
own `.nw` input format, optional `mpirun` launcher when
`mpi_procs > 1`, ECL-2.0), and **xTB** (Grimme group's extended tight-
binding semiempirical method — fast geometry optimization +
property screening on systems too big for full ab initio, LGPL-3.0).
Psi4 reads `.in` Python-scriptable input files; NWChem reads its own
`.nw` input format with optional `mpirun` wrapping for parallel runs;
xTB reads `.xyz` coordinates directly with all options on the CLI and
captures stdout to `xtb.log` via the MAFFT-style stdout-redirect
pattern. Phase 25 sits numerically between Phase 24 and Phase 27 but
ships chronologically right after Phase 28 — the same chronological-
vs-numerical convention used for Phase 17.5 / 24 / 28.

## Capability inventory

### Live adapters (3)

- **Psi4** — open-source HF/DFT/post-HF quantum chemistry (Justin
  Turney et al., LGPL-3.0). Single-binary subprocess shape: Psi4 reads
  a Psithon (Python-scriptable) input file and writes its run report
  to a user-named output file. Schema knobs: `input` (`.in` / `.dat`
  Psithon script; required), `output` (output filename relative to
  workdir; required), `threads` (default 1; ≥ 1), `memory` (default
  `"1 gb"`; matches `^\d+\s*(mb|gb|MB|GB)$` via the
  `is_valid_memory` helper), `extra_args`. `prepare()` builds
  `psi4 -i <input> -o <output> -n <threads> [-m <memory>] [extras...]`.
  The `-m` flag is only emitted when the user asked for something
  other than the documented `"1 gb"` default — passing `-m` every
  time would override Psi4's internal `"500 mb"` default with our
  fixed value even when the user didn't ask for one. `collect()`
  reports `output` as a `Log` artifact "Psi4 output" and walks the
  workdir for `.fchk` (`Native`, "Psi4 formatted checkpoint") and
  `.molden` (`Native`, "Psi4 Molden orbital data") files. Probe via
  `find_on_path(&["psi4"])`. `bio.psi4.compute` ribbon capability.
- **NWChem** — Pacific Northwest National Laboratory's massively-
  parallel ab initio + plane-wave DFT package (ECL-2.0). Single-
  binary subprocess shape with optional MPI wrapping: NWChem reads
  its own `.nw` input format and writes its run report to stdout —
  captured to a user-named output file via the MAFFT-style stdout-
  redirect pattern. Schema knobs: `input` (`.nw` NWChem-format
  script; required), `output` (output filename relative to workdir;
  required), `mpi_procs` (default 1; ≥ 1), `extra_args`. `prepare()`
  builds — serial: `nwchem [extras...] <input>`; parallel:
  `mpirun -n <mpi_procs> [extras...] nwchem <input>`. When
  `mpi_procs > 1`, prepare resolves `mpirun` via `find_on_path` and
  fails with a helpful install-hint `InvalidCase` if it's missing
  (`apt install openmpi-bin` / `apt install mpich`) rather than
  letting the child fail later with a less obvious "command not
  found". Output path is stashed in
  `PreparedJob.environment[VALENX_NWCHEM_OUTPUT]` so `run()` can
  redirect stdout to it without re-parsing the case TOML. `collect()`
  reports `output` as a `Log` artifact "NWChem output". Probe via
  `find_on_path(&["nwchem"])`. `bio.nwchem.compute` ribbon
  capability.
- **xTB** — Stefan Grimme's extended tight-binding semiempirical
  quantum chemistry package (LGPL-3.0). Single-binary subprocess
  shape with stdout-redirect: xTB reads `.xyz` coordinates directly
  and writes its run report to stdout — captured to `xtb.log` via
  the MAFFT-style stdout-redirect pattern. Schema knobs: `input`
  (`.xyz` geometry; required), `mode` ∈
  `{single-point, opt, ohess, hess, md}` (default `"single-point"`),
  `charge` (electron-balance `i32`; default 0), `uhf` (xTB's
  multiplicity convention — number of unpaired electrons, `u32`;
  default 0), `gfn` ∈ `{0, 1, 2}` (GFN method; default 2 — GFN2-xTB
  is the modern default), `solvent` (optional ALPB solvent name e.g.
  `"water"` / `"thf"`; `None` = gas phase), `extra_args`. `prepare()`
  builds `xtb <input> --gfn <gfn> --chrg <charge> --uhf <uhf>
  [--<mode> if mode != "single-point"] [--alpb <solvent> if Some]
  [extras...]`. `single-point` is xTB's default run type so it gets
  no flag; every other mode maps to `--<mode>`. Charge, multiplicity,
  and the GFN parameter set are always emitted so the invocation is
  unambiguous regardless of whether xTB's own defaults match ours.
  `collect()` reports `xtb.log` as a `Log` artifact "xTB stdout log"
  and walks the workdir for `xtbopt.xyz` (`Native`, "xTB optimized
  geometry"), `xtbopt.log` (`Log`), `gradient` / `hessian` files
  (`Native`). Probe via `find_on_path(&["xtb"])`. `bio.xtb.compute`
  ribbon capability.

### Canonical types

**No new canonical types.** All three adapters consume user-
supplied input files (Psithon `.in` / NWChem `.nw` / xyz `.xyz`
coordinates) and emit user-readable artifacts (text output reports,
formatted checkpoints, Molden orbital data, optimized xyz geometries,
gradient / hessian files) that the unchanged `Results.artifacts`
collection model surfaces directly. A first-class quantum-chemistry
canonical type — a generic energy / geometry / orbital data type
spanning all three back-ends — defers to a future phase along with
`.fchk` / Molden / `.cube` reader CLIs and visualization integrations.

### Headless CLIs

**No new CLIs.** Psi4 / NWChem / xTB output files are short text
reports that the user can inspect in any editor; richer formats
(`.fchk` formatted checkpoints, `.molden` orbital data, `.cube`
volumetric grids) are best inspected with the dedicated viewers
(Avogadro 2 from Phase 24, ChimeraX from Phase 17, PyMOL from
Phase 23) that already ship in the Valenx adapter zoo.

## Domain milestone

Phase 25 is the **first quantum-chemistry domain** to land in Valenx.
The biology adapter family started with Phase 17 (foundation —
sequence / structure / trajectory canonical types + classical MD +
cheminformatics scripts), expanded through Phase 17.5 / 18 / 18.5 /
18.6 / 19 / 19.5 / 20 / 22 / 23 / 24 / 27 / 27.5 / 28 / 30 / 34 to
cover sequence prediction, alignment, RNA-seq, variant calling,
single-cell, transcript quantification, workflow orchestration,
molecular viewers, cheminformatics, protein design, RNA structure,
phylogenetics, and small-molecule docking — but until Phase 25 the
quantum-mechanics surface (HF / DFT / post-HF / semiempirical
methods) was absent. Phase 25 closes that gap with three established
open-source tools that span the quantum-chemistry tradeoff space —
xTB at the fast-and-approximate end, Psi4 in the general-purpose
HF/DFT/post-HF middle, and NWChem at the massively-parallel ab initio
end.

## What landed early

The implementation landed across 5 discrete
implementation commits (3 adapters, the registry rollup, the init-
template extension) plus this docs pass — each landing one adapter,
the registry rollup, the init-template extension, or the
documentation pass. Every commit kept workspace clippy + rustdoc
clean.

## Acceptance checklist

- [x] `valenx-adapter-psi4` adapter ships with case-input parser
      + 4 lib tests + 4 case-input tests covering parse-minimal /
      parse-with-overrides / reject-zero-threads / reject-invalid-
      memory plus the `is_valid_memory` helper unit tests
- [x] `valenx-adapter-nwchem` adapter ships with case-input parser
      + 4 lib tests + 3 case-input tests, plus an `InvalidCase` path
      asserting `mpi_procs > 1` without `mpirun` on PATH surfaces a
      helpful install-hint error rather than a downstream "command
      not found"
- [x] `valenx-adapter-xtb` adapter ships with case-input parser
      + 4 lib tests + 5 case-input tests covering the
      `mode ∈ {single-point, opt, ohess, hess, md}` and
      `gfn ∈ {0, 1, 2}` whitelists
- [x] All 3 adapters wired into `valenx-app::init_registry` —
      live adapter count moves from 67 to **70**, opening the first
      quantum-chemistry domain to ship in Valenx
- [x] 3 quantum-chemistry templates in `valenx-init` (`psi4`,
      `nwchem`, `xtb` — canonical names only, no aliases beyond
      themselves), all round-tripping through `valenx-validate`
      (cross-binary roundtrip now sweeps **66 templates** clean)
- [x] STATUS.md / README.md / ROADMAP.md / CHANGELOG.md /
      QUICKSTART.md updated
- [ ] **Phase 25.5** — CP2K / Quantum ESPRESSO / GAMESS-US (different
      shape — plane-wave / massively parallel; defer to a sister-
      adapter expansion phase), DFTB+ / ABINIT / Octopus (niche;
      defer), ORCA (proprietary binary, free-tier-only license; would
      need a separate license-mode flag like AlphaFold 3), PySCF
      (Python library, fits Phase 24 cheminformatics pattern; can be
      added there). Out of scope for this beachhead.

## Success metrics

| Metric                                        | Target          |
|-----------------------------------------------|-----------------|
| New quantum-chemistry adapter (template + tests) | 1 day per       |
| Quantum-chemistry compute across 3 tools       | < tool baseline |

## Leads into

Phase 25 opens the quantum-chemistry domain that the user's bio /
chemistry spec called out alongside the Phase 17 / 17.5 / 27 / 27.5
biology + protein-design stack and the Phase 24 cheminformatics
expansion. Combined with the existing fold → analyze → predict →
infer-tree → validate loop, the **build-geometry → optimize →
compute-energy → predict-structure → fold-RNA → infer-tree →
validate** loop now spans three quantum-chemistry tools (Psi4,
NWChem, xTB) feeding into the Phase 24 cheminformatics surface
(DeepChem, Open Babel, Avogadro 2), the Phase 17 / 17.5 prediction
stack (ESMFold, OpenFold, AlphaFold 2/3, ColabFold), the Phase 28
RNA-structure tools (ViennaRNA, RNAstructure, NUPACK), and the
Phase 30 phylogenetic-tree builders (IQ-TREE, RAxML-NG, FastTree) —
all in one Valenx shell with no glue code beyond the existing case-
toml / prepare / run / collect path.

The natural follow-up is **Phase 25.5** — the deferred quantum-
chemistry work called out above (CP2K, Quantum ESPRESSO, GAMESS-US,
DFTB+, ABINIT, Octopus, ORCA with its proprietary license-mode
handling, PySCF as a Phase 24 cheminformatics-style Python adapter),
slotting in alongside the existing quantum-chemistry adapters with
the same single-binary subprocess shape (or the OpenMM / Scanpy
Python-script-subprocess shape for PySCF).
