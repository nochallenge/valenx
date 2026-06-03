# Phase 2 — CAD + meshing integration

**Status:** 🟢 In progress — gmsh adapter live; FreeCAD/OCCT/Netgen still scaffolded.

## Goal

Make Valenx useful for people who have geometry in STEP/IGES and want
a mesh — without leaving the app.

## Capability inventory

- STEP (AP242 + AP214) and IGES import via FreeCAD's headless
  kernel.
- OpenCASCADE (OCCT) dynamic-linked for in-process BRep operations
  (booleans, fillets, chamfers).
- Parametric feature tree in the browser — edit a parameter, the
  model rebuilds.
- gmsh-driven unstructured tet/hex meshing with prism layers for
  boundary-layer resolution.
- Netgen as an alternate mesher for curved-geometry biased work.
- Mesh quality metrics: aspect ratio, skewness, orthogonality,
  non-orthogonality angle.
- Boundary group naming preserved across CAD → mesh → solver.

## Integrated tools graduating to Implemented

| Tool       | Adapter crate                     | Status this phase      |
|------------|-----------------------------------|------------------------|
| FreeCAD    | `valenx-adapter-freecad`          | 🟢 probe + prepare + run + collect (import + export; parameter rebuild still pending) |
| OCCT       | `valenx-adapter-occt`             | 🔲 probe + dynamic-linked translate |
| gmsh       | `valenx-adapter-gmsh`             | 🟢 probe + prepare + run + collect (live) |
| Netgen     | `valenx-adapter-netgen`           | 🔲 probe + prepare + run + collect |

## Acceptance checklist

- [x] Import a STEP file via FreeCAD — the `valenx-adapter-freecad`
      adapter generates a deterministic `FreeCADCmd` Python script
      that opens the STEP, walks `ActiveDocument.Objects`, exports
      STL/BREP/STEP, and emits `summary.json` with parts, feature
      tree, bounding box, volume, and area. Wiring the feature
      tree into the Valenx browser is the next UX step.
- [ ] Edit a linear-pattern parameter, model rebuilds in < 500 ms.
- [x] Generate an unstructured tet mesh via gmsh (procedural Box /
      Sphere domains, or Merged STL) — `.geo` generator +
      `.msh` v4.1 parser both shipped and tested.
- [x] Prism-layer generation: the gmsh `.geo` writer now emits a
      `BoundaryLayer` field (`hwall_n` + `ratio` + `thickness` from
      the typed `[mesh.boundary_layer]` section). Users specify
      first-cell thickness + growth rate + layer count in the case
      file; the generator computes the total stack thickness from
      the geometric series and stamps the field onto the domain's
      Physical Surface. y+ calibration to a specific flow still
      lives in the user's hands (needs friction velocity).
- [x] Mesh quality metrics in `valenx-mesh::quality` (signed
      volume / area / length, aspect ratio, inverted-element
      count). Displaying them in the inspector UI is the next UX
      step.
- [x] Run the OpenFOAM simpleFoam pipeline against this mesh:
      `valenx-adapter-openfoam::prepare()` now detects a
      `mesh.msh` in either the case directory or the workdir and
      invokes `gmshToFoam` to materialise `constant/polyMesh/`
      before emitting the solver dicts. Structured
      `ToolNotInstalled` / `Run` errors surface when gmshToFoam
      is missing or rejects the mesh. The end-to-end
      "converges on a real geometry" demo still needs OpenFOAM
      installed to actually execute, but every adapter-side link
      in that chain is now in place.

## Success metrics

| Metric                                              | Target   |
|-----------------------------------------------------|----------|
| STEP import (5 MB, ~50 parts)                       | < 3 s    |
| gmsh tet mesh for a 500k-cell domain                | < 30 s   |
| Feature-tree rebuild latency on parameter change     | < 500 ms |
| Mesh quality thresholds surfaced in UI               | yes      |

## Leads into

[Phase 3 — Finite-element analysis](./phase-03-fea.md).
