# Valenx Architecture

This is the "how it fits together" doc. Read this before you read code
or file an RFC — it's the shared mental model everything else refers to.

Nothing here is a specification; specs live in the [RFCs](./rfcs/).
This is the overview that makes the specs make sense. Direct-reference
map from this doc to the specs:

- Project file format → [RFC 0001](./rfcs/0001-project-file-format.md)
- Adapter contract → [RFC 0002](./rfcs/0002-adapter-contract.md)
- Plugin (WASM) API → [RFC 0003](./rfcs/0003-plugin-api.md)
- Results and fields data model → [RFC 0004](./rfcs/0004-results-and-fields.md)
- Design principles → [RFC 0005](./rfcs/0005-design-principles.md)
- Design token system and pipeline → [RFC 0006](./rfcs/0006-token-system.md)

---

## 1. The one-paragraph summary

Valenx is a native desktop app written in Rust. Its job is to give
simulation engineers and life-science researchers a unified frontend
over ~141 open-source scientific-computing tools — CAD, meshing, CFD,
FEA, EM, chemistry, MD, battery, multi-physics, plus biology
(alignment / CRISPR design / structure prediction / variant calling
/ single-cell / phylogenetics / RNA folding / pharmacokinetics, with
biology now the largest surface at 123 of 141 adapters) — without
rewriting the science those tools already do well. The app owns the
UX, the project model, the workflow engine, the 3D viewer, and thin
Rust adapters that translate between a canonical case representation
and each tool's native inputs. The tools themselves run beside (or
inside) the app according to their license: permissively-licensed
code is bundled, LGPL code is
dynamically linked, GPL code runs as a child process. Reproducibility
is guaranteed by version-pinning every external tool per Valenx
release.

---

## 2. The layer cake

```
┌──────────────────────────────────────────────────────────────────┐
│                          valenx-app                              │
│  Ribbon · Browser tree · Timeline · Viewport · Command palette   │
│                        (egui + wgpu)                             │
├──────────────────────────────────────────────────────────────────┤
│                         valenx-core                              │
│  Project model · Workflow DAG engine · Adapter registry ·        │
│  Unit system · Case validation · Settings · Logging · IPC        │
├──────────────────────────────────────────────────────────────────┤
│         valenx-geo    valenx-mesh    valenx-fields               │
│         (canonical geometry · mesh · results types)              │
├──────────────────────────────────────────────────────────────────┤
│                       valenx-adapters-*                          │
│  One crate per integrated OSS tool. Implements RFC 0002.         │
│  Translate canonical Case → native inputs, invoke tool,          │
│  collect native outputs → canonical Results.                     │
├──────────────────────────────────────────────────────────────────┤
│                       External tools                             │
│  OpenFOAM · FreeCAD/OCCT · gmsh · Code_Aster · CalculiX · Elmer  │
│  Cantera · SU2 · openEMS · LAMMPS · PyBaMM · MuJoCo · preCICE …  │
│                                                                  │
│  Integration modes:                                              │
│     Bundled   — compiled into Valenx (permissive licenses)       │
│     DynamicLinked — .so / .dll at runtime (LGPL)                 │
│     Subprocess   — child process, license-isolated (GPL)         │
└──────────────────────────────────────────────────────────────────┘
```

**Reading from the top down:** The app never talks directly to an
external tool. Every cross-tool interaction goes through the adapter
layer and the canonical types in `valenx-geo` / `valenx-mesh` /
`valenx-fields`. That's what keeps 141 tools from becoming 141 special
cases in the UI.

**Reading from the bottom up:** Each external tool has exactly one
adapter crate. Each adapter exposes the same trait (`Adapter`,
defined in RFC 0002). The registry in `valenx-core` walks them all on
startup, probes, and tells the UI what's available.

---

## 3. The Cargo workspace

```
valenx/
├── Cargo.toml                          ← workspace manifest
├── crates/
│   ├── valenx-app/                     ← the binary (egui + wgpu frontend)
│   ├── valenx-core/                    ← registry · workflow · project model
│   ├── valenx-geo/                     ← canonical geometry types
│   ├── valenx-mesh/                    ← canonical mesh types
│   ├── valenx-fields/                  ← canonical field / results types
│   ├── valenx-viz/                     ← wgpu 3D viewer
│   ├── valenx-scripting/               ← embedded Python (pyo3) + Lua (mlua)
│   ├── valenx-plugins/                 ← WASM plugin host (wasmtime)
│   ├── valenx-bench/                   ← validation suite
│   └── valenx-adapters/
│       ├── cad/
│       ├── mesh/
│       ├── cfd/
│       ├── fea/
│       ├── em/
│       ├── chem/
│       ├── md/
│       ├── battery/
│       ├── robotics/
│       ├── coupling/         (preCICE)
│       ├── opt/
│       └── viz/              (ParaView, VTK)
├── tests/                              ← cross-crate integration tests
└── installers/                         ← per-platform packaging
```

Every adapter under `crates/valenx-adapters/` has the same internal
shape (see RFC 0002 § Reference-level explanation). The uniform shape
is what makes reviewing 141 adapters tractable: diffs are always in
`translate/` and `collect/`, never in the scaffolding around them.

---

## 4. Canonical types — the shared vocabulary

Four crates carry data between adapters:

| Crate | Types | Produced by | Consumed by |
|-------|-------|-------------|-------------|
| `valenx-geo` | `Geometry`, `BRep`, `Shape` | CAD adapters (FreeCAD, OCCT, STEP) | Meshers |
| `valenx-mesh` | `Mesh`, `Region`, `BoundaryGroup` | Meshers (gmsh, cfMesh, Netgen) | Solvers |
| `valenx-core::case` | `Case`, `Boundary`, `Material` | UI (user-authored) | Solvers |
| `valenx-fields` | `Field<Scalar/Vector/Tensor>`, `Results` | Solvers | Viewer, post-processors |

Adapters **never see each other's native types**. An OpenFOAM
adapter's internal `FoamDict` representation never crosses the
crate boundary; only `Mesh + Case` goes in, `Results` comes out. Swap
the solver, the rest of the workflow stays intact.

---

## 5. The adapter registry — runtime organization

```
┌───────────────────────── App startup ──────────────────────────┐
│                                                                │
│  1. inventory::iter() collects every adapter crate's           │
│     Adapter trait object                                       │
│                                                                │
│  2. rayon::par_iter() runs probe() on each, in parallel        │
│       ~200-500 ms total; target < 500 ms on warm install       │
│                                                                │
│  3. Classify: Ready / Missing / Outdated / Broken / Disabled   │
│                                                                │
│  4. Publish the registry to valenx-core                        │
│                                                                │
│  5. Emit ProbeComplete event; UI enables menus accordingly     │
│                                                                │
└────────────────────────────────────────────────────────────────┘
```

The registry is the **single source of truth** for "what physics does
this install of Valenx support?" It's indexed three ways:

- **By tool id** — `"openfoam"`, `"gmsh"`, ...
- **By physics** — `Physics::Cfd`, `Physics::Fea`, ...
- **By capability** — `Capability::Meshing2D`, `Capability::AdjointOptimization`, ...

Every UI piece that asks "should this be enabled?" queries the
registry. The user experiences this as: "I installed OpenFOAM but
not SU2, so the SU2 path disappears from my New Project dialog," not
as a generic error.

---

## 6. The workflow DAG — orchestration

A user action like "run this case" is really a DAG of adapter calls
over canonical types:

```
   Geometry       Mesh          Case         Results        VizScene
      │             │             │             │              │
      ▼             ▼             ▼             ▼              ▼
 ┌─────────┐  ┌─────────┐  ┌─────────┐  ┌─────────┐  ┌─────────────┐
 │ freecad │→ │  gmsh   │→ │openfoam │→ │ collect │→ │   viz /     │
 │ (CAD)   │  │ (mesh)  │  │ (solve) │  │ (VTK)   │  │  paraview   │
 └─────────┘  └─────────┘  └─────────┘  └─────────┘  └─────────────┘
```

Nodes are adapter invocations. Edges are typed data from the canonical
crates. The DAG is serialized into the `.valenx` project file (per RFC
0001), which is how a project re-runs deterministically years later
with the same pinned tools.

Swap `gmsh` → `cfMesh`: both consume `Geometry`, both produce `Mesh` —
the CFD node doesn't know. Swap `OpenFOAM` → `SU2`: both consume
`Mesh + Case`, both produce `Results`. The post-processor doesn't know.

---

## 7. Multi-physics coupling

Coupled problems (FSI, CHT, reacting flow, electro-thermal) aren't a
separate subsystem. `valenx-adapter-precice` is a **meta-adapter**:

```
                   ┌──────────────────┐
                   │    preCICE       │
                   │    (coupler)     │
                   └────┬────────┬────┘
                        │ traction│displacement
                        │         │
                   ┌────▼──┐  ┌──▼─────┐
                   │openfoam│  │calculix│
                   │  (CFD)│  │ (FEA)  │
                   └───────┘  └────────┘
```

preCICE handles timestep synchronization, data mapping between
non-matching meshes, implicit coupling iteration, and convergence.
From `valenx-core`'s perspective the coupled subgraph looks like one
node; underneath, preCICE is orchestrating the participating solvers
concurrently. Adding a new coupled physics combination is mostly a
preCICE configuration exercise, not a new orchestration layer.

---

## 8. License isolation — the firewall

The only reason 141 heterogeneously-licensed tools can live under an
Apache 2.0 project is that we treat license mode as a first-class
structural concern, not an afterthought.

Every adapter declares a `LicenseMode`:

| Mode | Meaning | Enforcement |
|------|---------|-------------|
| `Bundled` | Tool compiled into the Valenx binary (statically or via dynamic lib we ship) | Licenses limited to Apache / BSD / MIT / ISC / MPL — see `deny.toml` |
| `DynamicLinked` | Linked at runtime to tool's `.so` / `.dll` / `.dylib` | LGPL compatible; link honored dynamically; user can swap in their own build |
| `Subprocess` | Tool runs as a child process, communicating via args + files | GPL tools live here; we never link against them; they cannot pollute our binary |

**Three enforcement layers:**

1. **CI (`cargo-deny`)** — workspace-level `deny.toml` plus crate-local
   rules. A Subprocess adapter with a link dep on the tool's `-sys`
   crate fails CI.
2. **Runtime (`valenx-core::process`)** — the subprocess-spawning
   helper cross-checks the adapter's declared `LicenseMode`; an
   adapter declared `Bundled` cannot accidentally shell out.
3. **OS (process isolation)** — GPL tools run in their own process
   space. Crashes there don't take down Valenx; more importantly,
   they don't share an address space, so no copyleft obligation
   attaches.

Result: the Valenx binary is Apache 2.0 top to bottom, even though
the user may have GPL binaries on disk that Valenx invokes.

---

## 9. Reproducibility — tools.lock

Every Valenx release ships a `tools.lock` at the repo root pinning
the exact version and checksum of every integrated external tool for
that release. When the user opens a project, the project's *own*
`tools.lock` (inside the `.valenx` bundle, per RFC 0001) is compared
against the installed versions. Running the same project on the same
Valenx release produces bit-identical results — that's the
reproducibility guarantee, and it's what makes Valenx viable for
certification work.

Three channels, different freeze cadences:

- **Stable** — `tools.lock` frozen for ~3 months between minor releases
- **LTS** — frozen for the 24-month life of the LTS tag
- **Nightly** — tracks upstream tools as they release

A user can override individual tools via `~/.valenx/tool-overrides.toml`,
at which point the reproducibility guarantee detaches for that tool and
a visible "unlocked" indicator appears in the UI.

---

## 10. The UI layer

The UI speaks the canonical `Case`. Adapters translate up and down.
The user never types a solver's native input deck unless they
explicitly ask to.

- **Ribbon** reshapes per physics. Adapters contribute ribbon entries
  through `Adapter::capabilities()` so enabling a new tool expands the
  ribbon without editing `valenx-app`.
- **Browser tree** — uniform shape regardless of which solver's
  behind the case: Project → Geometry → Mesh → Cases → Results.
- **Timeline** — ordered user actions, replayable; ~Fusion-style.
- **Log viewer** — streams the underlying tool's stdout/stderr verbatim,
  tabbed by tool name, so power users see exactly what OpenFOAM or
  CalculiX said.
- **"Reveal native inputs"** — expert escape hatch exposing the
  adapter's workdir with the translated input deck, for debugging or
  hand-tuning.

The UI is one window, rendered by `egui` through `wgpu`. No browser,
no localhost, no dev server, no Electron. See [TESTING.md](./TESTING.md)
for how to launch it.

---

## 11. Extensibility

Two extension surfaces, for different trust levels:

**In-tree adapters** (RFC 0002). First-party Rust crates under
`crates/valenx-adapters/`. Compiled with the app; go through CI and
code review. This is the path for maintainers and trusted
contributors integrating a new OSS tool.

**Plugins** (RFC 0003). Third-party `.wasm` files loaded at runtime
via `wasmtime`. Sandboxed; capabilities granted by the user at
install. This is the path for community extensions and proprietary
user code that the main project doesn't ship.

Both surfaces target the same canonical types (`Geometry`, `Mesh`,
`Case`, `Results`). A plugin adapter looks like an in-tree adapter
from the registry's point of view; it's just sandboxed.

---

## 12. Key numbers (targets, not current state)

| Metric | Target |
|--------|--------|
| Cold start to main window | < 1.5 s |
| Adapter probe (all adapters) | < 500 ms |
| STEP import, 100 MB file | < 10 s |
| Mesh visualization, 5M cells | interactive 30+ fps |
| Binary size (base installer) | ~250 MB |
| Compiled Rust (debug) | ~3 GB workspace target |
| Tool downloads (optional) | up to several GB depending on selection |

These are design constraints informing the architecture, not
promises about current code.

---

## 13. How this document stays accurate

- Structural changes (new layer, new canonical type, new integration
  mode) require an RFC that also updates this file
- Minor clarifications are normal PRs
- If the ROADMAP and this document disagree, the ROADMAP is
  forward-looking intent; this document is the current system design

See [rfcs/](./rfcs/) for the specifications this overview is built on,
and [ROADMAP.md](./ROADMAP.md) for the 20-year direction this
architecture is meant to carry.
