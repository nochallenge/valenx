# Introduction

Welcome to **Valenx** — a native desktop simulation suite. One download
gives you CAD, meshing, CFD, FEA, electromagnetics, chemistry, molecular
dynamics, battery modelling, and multi-physics coupling in a single
unified window. Apache 2.0, free forever, no subscription.

Valenx does not reinvent the physics. Under the hood it drives decades
of proven open-source solvers — OpenFOAM, Code_Aster, CalculiX, Elmer,
SU2, Cantera, openEMS, PyBaMM, and many more. What Valenx builds is the
workflow: a consistent native UI, a common project format, canonical
result types, and an adapter layer that lets you stop learning one
solver per physics.

---

## This manual

- **Getting started** — install the app, complete the first-run setup,
  run your first case in ten minutes
- **Concepts** — projects, physics, adapters, the tool registry, units
- **Physics** — per-vertical guides (CFD, FEA, …)
- **Workflow** — the tasks you do in any project (geometry → mesh →
  boundaries → solve → post-process → export)
- **Reference** — keyboard shortcuts, settings, file formats, CLI
- **Scripting** — Python and Lua inside the app
- **Plugins** — installing and authoring third-party extensions
- **Tutorials** — end-to-end, copy-along examples
- **Validation gallery** — canonical reference cases with published
  benchmarks

## Status

Valenx is **pre-alpha**. The architecture and plan are set; the Rust
workspace scaffold is in place; the native app is being built Phase by
Phase (see the [roadmap][roadmap] for the 20-year plan). Sections of
this manual marked *TBD* describe intent rather than current behaviour.

## Getting help

- **Discussions** — open-ended questions
- **Issues** — bug reports, feature requests
- **Security** — see [SECURITY.md][security]

See the [Contributing][contributing] chapter for how to get involved.

[roadmap]: https://github.com/nochallenge/valenx/blob/main/ROADMAP.md
[security]: https://github.com/nochallenge/valenx/blob/main/SECURITY.md
[contributing]: ./contributing.md
