# Phase 24 — Cheminformatics expansion

**Status:** 🟢 Live — DeepChem + Open Babel + Avogadro 2 round
out the cheminformatics surface that Phase 17's RDKit started.

## Goal

Round out the cheminformatics surface that Phase 17's RDKit
adapter started. Phase 24 ships three sister adapters:
**DeepChem** (PyTorch-backed deep-learning cheminformatics
library — sister to RDKit's classical chemistry), **Open Babel**
(the de-facto open-source chemistry-format converter — handles
~120 chemistry file formats), and **Avogadro 2** (Python-
scriptable chemistry editor with a small-molecule rendering
pipeline). Together with RDKit (already shipped) and the
Phase 34 docking adapters, Valenx now has the complete
small-molecule + cheminformatics stack. All three follow
established Phase 17 / Phase 18 patterns: DeepChem mirrors
RDKit's Python-script subprocess shape, Open Babel uses BWA's
single-binary CLI shape (`obabel <in> -O <out>`), and
Avogadro 2 mirrors ChimeraX's script-driven-headless pattern.

## Capability inventory

### Live adapters (3)

- **DeepChem** — PyTorch-backed deep-learning cheminformatics
  library; sister to RDKit's classical chemistry. Drives off a
  user-provided Python script that imports `deepchem` and
  reads `valenx_params.json` (written by the adapter into the
  workdir) for config knobs: optional inline `smiles` list
  (passed through the params file for the script to consume),
  optional `dataset_csv` (staged into the workdir), optional
  `checkpoint` model path. Output classification walks the
  workdir for `.csv` (kind `Tabular`, label
  `"DeepChem analysis output"`), `.png` (kind `Native`, label
  `"DeepChem plot"`), and `.pkl` / `.pt` (kind `Native`, label
  `"DeepChem model checkpoint"`). MIT licensed.
  `bio.deepchem.script` ribbon capability.
- **Open Babel** — the de-facto open-source chemistry-format
  converter; `obabel` translates between ~120 file formats
  (SMILES, MOL, MOL2, PDB, SDF, XYZ, …). Single-binary CLI
  shape: `obabel <input> -O <output> [-i <input_format>]
  [-o <output_format>] [--gen3D] [-h] [extras…]`. `gen_3d`
  (default `false`) toggles `--gen3D` for 2D → 3D coordinate
  generation; `add_hydrogens` (default `false`) toggles the
  `-h` hydrogen-adding flag; explicit `input_format` /
  `output_format` overrides let users pin a format that the
  extension would mis-detect. Output collected as a `Native`
  artifact with label `"Open Babel converted file"`. GPL-2.0
  licensed. `bio.openbabel.convert` ribbon capability.
- **Avogadro 2** — Python-scriptable chemistry editor with a
  small-molecule rendering pipeline. Drives off a user-supplied
  Python script via `avogadro2 --script <script.py>`; an
  optional `structure` field (`.cml` / `.mol` / `.xyz` /
  `.pdb`) gets staged + passed as a positional arg so the
  script doesn't need to know the path. `headless` (default
  `true`) toggles `--no-gui` for batch / CI use. Output
  classification walks the workdir for `.png` (label
  `"Avogadro 2 render"`), `.cml` / `.mol` / `.xyz` (label
  `"Avogadro 2 exported structure"`). GPL-2.0-or-later
  licensed. `bio.avogadro.render` ribbon capability.

### Canonical types

**No new canonical types.** All three adapters consume existing
Phase 17 / 18 / 19 inputs (PDB / SDF / SMILES / CSV) and emit
user-readable artifacts (CSV tables, PNG renders, exported
chemistry structures, PyTorch / Pickle model checkpoints) that
the unchanged `Results.artifacts` collection model surfaces
directly. Open Babel's converted outputs land in whatever
format the user requested; the existing PDB reader inspects
the PDB targets, the existing FASTA reader handles SMILES /
sequence outputs.

### Headless CLIs

**No new CLIs.** DeepChem's CSV outputs are inspectable through
standard Unix tooling; Open Babel's converted files land in
formats the existing `valenx-pdb-info` (PDB / PDBQT) and
`valenx-fasta` (SMILES / sequence) CLIs already cover;
Avogadro 2's PNG / CML / MOL / XYZ outputs are user-readable
directly.

## What landed early

The implementation landed across 6
discrete commits, each landing one adapter, the registry
rollup, the init-template extension, or the documentation
pass. Every commit kept workspace clippy + rustdoc clean.

## Acceptance checklist

- [x] `valenx-adapter-deepchem` adapter ships with case-input
      parser + 4 lib tests + 3 case-input tests
- [x] `valenx-adapter-openbabel` adapter ships with case-input
      parser + 4 lib tests + 4 case-input tests
- [x] `valenx-adapter-avogadro` adapter ships with case-input
      parser + 4 lib tests + 3 case-input tests
- [x] All 3 adapters wired into `valenx-app::init_registry` —
      live adapter count moves from 44 to 47
- [x] 3 cheminformatics templates in `valenx-init` (`deepchem`
      with `dc` alias, `openbabel` with `obabel` alias,
      `avogadro` with `avogadro2` alias), all round-tripping
      through `valenx-validate` (cross-binary roundtrip now
      sweeps 43 templates clean)
- [x] STATUS.md / README.md / ROADMAP.md / CHANGELOG.md /
      QUICKSTART.md updated
- [ ] **Phase 24.5** — sister-adapter expansion: OEChem (commercial-
      use restrictions; defer to enterprise extensions), datamol
      (Python lib already covered transitively by the RDKit
      adapter), Open Babel Python bindings (the `pybel` module —
      route through the Biopython adapter for users who want it),
      ML model training pipelines for DeepChem (user-script
      territory). Out of scope for this beachhead.

## Success metrics

| Metric                                            | Target          |
|---------------------------------------------------|-----------------|
| New cheminformatics adapter (template + tests)    | 1 day per       |
| SMILES → 3D structure → render loop               | < tool baseline |

## Leads into

Phase 24 paired with the Phase 17 RDKit adapter and the
Phase 34 Vina + AutoDock 4 docking adapters gives Valenx the
complete small-molecule cheminformatics chain in one shell:
RDKit / DeepChem prepares the ligand → Open Babel converts
between formats and generates 3D coordinates → Avogadro 2
edits / renders the structure → AutoDock Vina or AutoDock 4
docks against a receptor → ChimeraX / PyMOL / VMD (Phases 17
+ 23) renders the ranked poses. Sequence-driven prep
(Biopython / RDKit Python scripts) stays user-side; the
adapters cover the conversion, ML, and rendering links.

The natural follow-up is **Phase 24.5** — the deferred
cheminformatics tools called out above (OEChem behind the
enterprise paywall, the pybel Python bindings via Biopython,
ML training pipelines as user scripts).
