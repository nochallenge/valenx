# Languages & Tools

Every language in the project, what it's used for, why.

---

## Primary: **Rust** (~95% of code we write)

### Why Rust
- **Performance:** C/C++ speed, no garbage collection pauses
- **Safety:** memory safety without a runtime; prevents an entire
  category of bugs that plagues C++ numerics code
- **Cross-platform:** one codebase ships on Win, macOS, Linux
- **Modern tooling:** `cargo` is the best package manager in any systems
  language; rustfmt, clippy, rust-analyzer all top-tier
- **Concurrency:** fearless parallelism via `rayon`, `tokio`, `crossbeam`;
  huge for solvers
- **FFI is first-class:** binding to OpenCASCADE, VTK, BLAS, MPI is clean
- **Growing scientific ecosystem:** `nalgebra`, `ndarray`, `faer-rs`,
  `rsmpi`, `sprs`, `peroxide` — not as mature as Python's, but catching up fast
- **No garbage collector** means real-time-capable if we ever want
  live-interaction physics (soft body deformation, haptics, etc.)

### Where Rust lives in the stack

- The entire native app (`valenx-app`)
- Orchestration, project model, workflow engine (`valenx-core`)
- 3D results viewer (`valenx-viz` via `wgpu`)
- Adapters for every external solver (subprocess drivers)
- All native physics solvers we write (`valenx-cfd-native`, `valenx-fem`, etc.)
- Plugin API (via WIT / `wit-bindgen` / `wasmtime`)

### Rust crates we depend on

| Purpose | Crate |
|---|---|
| **UI framework** | `egui` (immediate-mode; default choice) |
| **3D graphics** | `wgpu` (Vulkan/Metal/DX12/OpenGL abstraction) |
| **Linear algebra** | `nalgebra` (dense), `sprs` (sparse), `faer-rs` (fastest) |
| **Arrays** | `ndarray` (NumPy-like) |
| **Parallelism** | `rayon` (shared-memory), `rsmpi` (distributed MPI) |
| **Serialization** | `serde`, `toml`, `bincode`, `ciborium` (CBOR) |
| **Persistence** | `sqlx` (SQLite backend for run history) |
| **Async** | `tokio` (for job queue, subprocess management) |
| **Logging** | `tracing` + `tracing-subscriber` |
| **Errors** | `thiserror` (library errors), `anyhow` (app errors) |
| **Testing** | `cargo test`, `criterion` (benchmarks), `proptest`, `insta` |
| **Python bridge** | `pyo3` (embed Python interpreter for scripting) |
| **Lua bridge** | `mlua` (alternative scripting) |
| **Subprocess** | `tokio::process` (for adapter calls to OpenFOAM etc.) |
| **Plugin host** | `wasmtime` + `wit-bindgen` (sandboxed WASM plugins) |

### MSRV (Minimum Supported Rust Version)

Stable channel, current minus 3 releases. We update ~quarterly.

---

## Secondary: **C** (~2% via FFI)

### Where
- `opencascade-rs` — binding to OCCT (OpenCASCADE is C++ but exposes a C
  API wrapper)
- `vtk-sys` — binding to VTK for results file parsing
- `blas-sys` / `lapack-sys` — binding to optimized numerical libraries
- `gmsh` C API — if we wrap it as a library rather than subprocess
- `cantera-sys` — if we want tighter integration than subprocess

### Why C, not C++
Rust can FFI to C cleanly. C++ FFI is painful because of name mangling,
ABI drift, and the complexity of C++ types. For the OSS tools we wrap,
they all expose C APIs — that's what we bind to. We never write new C
code; we only use it to talk to existing C libraries.

### We do NOT use C++ directly
No new C++ is written in this project. If we ever need to extend OCCT
or VTK, we contribute upstream in their language; we don't embed C++
into our workspace.

---

## Build glue: **TOML** (config)

- `Cargo.toml` — every crate's manifest
- `.valenx` project files — user-facing project saves
- Rust's ecosystem convention; trivial to parse with `serde`

---

## Scripting: **Python** (optional, for users)

### Why Python at all
Scientific computing is Python-heavy. Researchers already know it.
Offering Python scripting inside Valenx lets them:
- Write batch jobs
- Post-process results programmatically
- Define custom case templates
- Call NumPy/SciPy/matplotlib from inside the app

### How we embed it
`pyo3` embeds a Python 3.11+ interpreter inside the Rust binary. No
system Python required. A `valenx.scripting` module exposes the app's
internal data model to Python scripts.

### Alternative: Lua
Smaller footprint, faster startup. Option for users who find Python
overkill. Provided via `mlua`. Both languages will be supported in
parallel — users choose.

---

## Web (minimal usage)

### Docs site — Markdown → HTML
The **user-facing documentation site** uses `mdBook` (Rust-native static
site generator). Markdown source, HTML output. Hosted on GitHub Pages.
No runtime web server in the app itself.

### The OpenAPI docs page (`docs/index.html`)
A single static HTML file that was useful when we had a web backend.
**Deleted in this phase** — the native app replaces it with in-app
context-sensitive help.

---

## Build tooling

- `cargo` — primary build tool; manages the workspace
- `just` or `make` — optional task runner (Rust alternative: `cargo xtask`)
- **GitHub Actions** for CI/CD (YAML)
- **Pre-commit** hooks (Rust: `cargo fmt`, `cargo clippy --workspace --all-targets -- -D warnings`, `bash scripts/qa.sh` — NEVER `cargo test --workspace`; see docs/QA.md)

---

## Languages explicitly NOT used

| Language | Why not |
|---|---|
| **JavaScript / TypeScript** | No web frontend — the UI is pure Rust / egui |
| **C++** | See FFI section — we bind to C APIs only |
| **Fortran** | We bind to BLAS/LAPACK (their Fortran is wrapped by C) but don't write new Fortran |
| **Go** | No reason to mix Go into a Rust codebase |
| **Swift / Kotlin** | Native iOS/Android is Phase 9+ if at all |

---

## Language breakdown estimate (when project is mature)

| Language | % of LOC | Where |
|---|---|---|
| Rust | 95% | All new code |
| C (FFI) | 2% | Thin wrappers around OCCT, VTK, BLAS |
| TOML | 2% | Config files |
| Python (examples) | 1% | Tutorial scripts users will write |
| Shell (build scripts) | <1% | CI helpers |
| Markdown | huge, but not LOC | Docs |

---

## Quick answers to anticipated questions

**"Why not TypeScript like the old frontend was?"**
The old frontend was a browser app. We're building a desktop app.
TypeScript is great for web UIs; for native, Rust + a Rust UI toolkit
is faster, smaller, and doesn't depend on Chromium being installed.

**"Why not Go or Zig or Nim?"**
- Go: garbage-collected; pauses hurt solver performance
- Zig: too new; small ecosystem; no mature numerical libraries
- Nim: small community; less stable release cadence
- Rust has the best combination of maturity + performance + safety in
  a modern language in 2026.

**"Why not just use C++?"**
C++ is fine for performance but worse for:
- Build tooling (Cargo is genuinely decade-ahead of CMake)
- Memory safety (we skip entire classes of bugs)
- Dependency management (Rust's strict versioning vs C++'s fragile ABIs)
- Parallelism (Rust's ownership model prevents data races at compile time)

For a 10-year project with high contributor turnover, Rust's safety is
worth a lot. C++ code written in 2015 with 5 authors is often a
maintenance nightmare by 2025; Rust code tends to stay reviewable.

**"What if Rust loses popularity?"**
Unlikely on a 10-year horizon — it's in the Linux kernel, widely adopted
at Google/Microsoft/Meta. If it did wane, the Rust → C++ migration path
is mechanical; there are bindings for the heavy lifting. But we're not
betting against Rust.
