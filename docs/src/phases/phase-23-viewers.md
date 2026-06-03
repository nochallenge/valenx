# Phase 23 — Molecular viewers

**Status:** 🟢 Live — molecular-viewer beachhead landed.

## Goal

Round out the visualization surface for everything Valenx's biology
stack produces. Phase 17 shipped the ChimeraX adapter as the first
script-driven molecular renderer; Phase 23 ships its three most-used
siblings — PyMOL (Python-scriptable, ubiquitous in structural
biology), VMD (Tcl-scriptable, the de-facto MD trajectory viewer),
and IGV (genome browser — `igvtools` headless mode for indexing +
tile generation). All three follow the established Phase 17 ChimeraX
shape: script-driven subprocess, headless mode, output-in-workdir.

## Capability inventory

### Live adapters (3)

- **PyMOL** — open-source PyMOL build (the
  Schrödinger fork is proprietary; we wrap the BSD-licensed open-
  source line). Drives off `.pml` (Python-style) command files.
  Defaults to `-c` (no GUI) + `-q` (quiet) — i.e. `pymol -c -q
  <script.pml>` — and collects the `.png` / `.pse` / `.cif` / `.pdb`
  outputs the script generates. `bio.pymol.render` ribbon
  capability.
- **VMD** — Tcl-scripted MD trajectory viewer
  (`vmd -dispdev text -e <script>`). Optional `structure` field in
  the case-input schema lets a `.pdb` / `.gro` / `.psf` get pre-
  loaded as a positional arg without the script having to know
  the path. Collects `.png` / `.tga` / `.bmp` renders, `.pdb` /
  `.gro` exported structures, and `.dat` analysis output. **Note:
  VMD ships under a custom non-OSS-but-free-for-academic-use
  license**; the adapter's probe pushes a license-awareness
  warning into `ProbeReport.warnings` reminding users to confirm
  their use case before redistributing renders or derived data.
  `bio.vmd.render` ribbon capability.
- **IGV** — `igvtools` wrapper for headless BAM / VCF / WIG
  indexing + tile generation. Per-action dispatch on
  `action ∈ {index, count, sort, tile}` — `index` writes the
  `.bai` / `.idx` sidecar next to the input (no `output` field),
  the other three actions consume an explicit `output` path.
  `count` exposes the conventional 25-bp default `window_size`.
  The companion GUI viewer (interactive genome browser) is out
  of scope — this is the headless-tooling adapter only.
  `bio.igv.index` ribbon capability.

### Canonical types

**No new canonical types.** All three viewers consume the existing
Phase 17 / Phase 18 / Phase 19 inputs (`Structure` / PDB, BAM,
trajectories) and emit user-readable artifacts (PNG / PSE / TGA
renders, BAI / IDX index sidecars) that the unchanged
`Results.artifacts` collection model surfaces directly.

### Headless CLIs

**No new CLIs.** The viewers' outputs are images and binary indices
that are best inspected by their native tooling, not summarised
through a Valenx text/JSON inspector.

## What landed early

The implementation rode subagent-driven-development across 6 discrete
commits, each landing one adapter, the registry rollup, the init-
template extension, or the documentation pass. Every commit kept
workspace clippy + rustdoc clean.

## Acceptance checklist

- [x] `valenx-adapter-pymol` adapter ships with case-input parser
      + 4 lib tests + 3 case-input tests
- [x] `valenx-adapter-vmd` adapter ships with case-input parser +
      4 lib tests + 3 case-input tests + 1 probe-warning test
      asserting the `"academic"` keyword surfaces
- [x] `valenx-adapter-igv` adapter ships with case-input parser
      (per-action dispatch) + 4 lib tests + 5 case-input tests
- [x] All 3 adapters wired into `valenx-app::init_registry` —
      live adapter count moves from 37 to 40
- [x] 3 viewer templates in `valenx-init` (`pymol` / `vmd` / `igv`
      with `igvtools` alias), all round-tripping through
      `valenx-validate` (cross-binary roundtrip now sweeps 36
      templates clean)
- [x] VMD's academic-only license surfaces in the probe warnings
      so users see it before they ship renders downstream
- [x] STATUS.md / README.md / ROADMAP.md / CHANGELOG.md /
      QUICKSTART.md updated
- [ ] **Phase 23.5** — in-app embedded viewers (Mol* / NGL Viewer /
      3Dmol.js webview integration), Avogadro 2 (cheminformatics —
      different shape; slots into Phase 24), LiteMol / iCn3D
      (niche enough to wait on user demand). Out of scope for this
      beachhead; next plan covers it.

## Success metrics

| Metric                                        | Target          |
|-----------------------------------------------|-----------------|
| New viewer adapter (template + tests)         | 1 day per       |
| Headless render / index round-trip            | < tool baseline |

## Leads into

Phase 23.5 — in-app embedded viewers (Mol* / NGL Viewer /
3Dmol.js webview integration) so the user can see the renders
without leaving the app; Avogadro 2 cheminformatics adapter
slots into Phase 24. See the future-phases table at the end of
`docs/superpowers/plans/2026-04-30-biology-foundation.md`
for the full follow-up phase list (Phases 19.5 → 43 cover the
remaining ~190 tools from the user's spec).
