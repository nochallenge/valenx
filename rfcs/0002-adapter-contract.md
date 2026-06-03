# RFC 0002: Adapter Contract

- **Status:** Accepted (initial design)
- **Author(s):** BDFL
- **Created:** 2026-04-21
- **Discussion PR:** (this commit)
- **Tracking issue:** TBD

---

## Summary

Define the Rust trait and lifecycle contract that every integrated
open-source tool (OpenFOAM, gmsh, Code_Aster, Cantera, …) must
implement to be driven by Valenx. Adapters are **in-tree Rust crates**
that translate Valenx's canonical case format into the tool's native
input, invoke the tool (subprocess or dynamic-link depending on
license), collect results, and surface errors.

---

## Motivation

Valenx's entire value proposition is: **one native frontend, many
validated solvers underneath**. That only works if:

1. Adding a new tool is a well-defined engineering task, not spelunking
2. Every adapter behaves the same from the core's point of view
   (progress reporting, cancellation, error classification, result
   metadata)
3. Solver-specific quirks stay inside the adapter, never leaking into
   `valenx-core`
4. License isolation is mechanical and auditable: a GPL subprocess
   adapter cannot accidentally link against the GPL tool

Without a contract, every adapter is a snowflake — exactly the mess we
want to avoid for a 20-year codebase.

---

## Guide-level explanation

### What an adapter is

An adapter is a Rust crate named `valenx-adapter-<tool>` that
implements the `Adapter` trait from `valenx-core::adapter`:

```rust
pub trait Adapter: Send + Sync {
    fn info(&self) -> AdapterInfo;
    fn probe(&self) -> Result<ProbeReport, AdapterError>;

    fn prepare(
        &self,
        case: &Case,
        workdir: &Path,
    ) -> Result<PreparedJob, AdapterError>;

    fn run(
        &self,
        job: &PreparedJob,
        ctx: &mut RunContext,
    ) -> Result<RunReport, AdapterError>;

    fn collect(
        &self,
        job: &PreparedJob,
    ) -> Result<Results, AdapterError>;
}
```

Each adapter crate also exports a plain constructor:

```rust
pub fn adapter() -> Box<dyn Adapter> { Box::new(OpenFoamAdapter::new()) }
```

which `valenx-core` calls at startup to register the adapter.

### Anatomy of a run

Five phases, every adapter, no exceptions:

1. **Info** — static metadata: name, supported physics, license class,
   minimum tool version
2. **Probe** — can we actually use this tool on this machine right now?
   (Check binary exists, version matches `tools.lock`, required
   dependencies present.)
3. **Prepare** — given a `Case`, write the solver's native input deck
   into a temporary working directory. Deterministic: same case → same
   files, byte for byte.
4. **Run** — invoke the tool. Report progress. Respect cancellation.
   Stream logs back to the UI.
5. **Collect** — parse outputs into the canonical `Results` structure;
   point at the raw files so the UI can load them for visualization.

The core orchestrates these phases in order. The adapter never talks to
the UI directly; all communication goes through `RunContext`.

### License modes

Every adapter declares its **license mode**, baked into its
`AdapterInfo`:

```rust
pub enum LicenseMode {
    /// Tool shipped as part of our binary; linked statically or dynamically
    /// with a permissive license (MIT, BSD, Apache).
    Bundled,

    /// Tool linked at runtime; license permits dynamic linking (LGPL).
    DynamicLinked,

    /// Tool runs as a child process; we never link against it. For GPL
    /// tools, this is the only permitted mode.
    Subprocess,
}
```

The `Subprocess` mode is enforced at the code level: a subprocess-mode
adapter crate **must not** have the GPL tool as a direct dependency in
`Cargo.toml`. CI checks this.

### Progress and cancellation

Adapters report progress through `RunContext`:

```rust
pub struct RunContext<'a> {
    pub cancel: &'a CancellationToken,
    pub progress: Box<dyn ProgressSink + 'a>,
    pub log: Box<dyn LogSink + 'a>,
}

impl<'a> RunContext<'a> {
    pub fn check_cancel(&self) -> Result<(), Cancelled> { … }
    pub fn report(&self, pct: f32, message: &str) { … }
    pub fn log_line(&self, level: LogLevel, line: &str) { … }
}
```

Cancellation must be checked at least every few seconds. For subprocess
adapters, `Drop` on the `PreparedJob` kills the child process.

### Error classification

Adapter errors are not just strings — they're a structured enum:

```rust
pub enum AdapterError {
    ToolNotInstalled { name: &'static str, hint: String },
    ToolVersionMismatch { expected: SemVer, found: SemVer },
    InvalidCase { case_path: PathBuf, reason: String },
    Translate(TranslateError),   // case → native input
    Run { exit_code: i32, stderr: String, phase: RunPhase },
    ParseOutput { file: PathBuf, reason: String },
    Cancelled,
    IO(std::io::Error),
    Other(anyhow::Error),
}
```

This matters because the UI renders them differently:

- `ToolNotInstalled` → modal with "install now" button
- `InvalidCase` → inline validation error on the offending parameter
- `Run` → log viewer showing the tool's own stderr
- etc.

---

## Reference-level explanation

### Directory structure

```
crates/
└── valenx-adapters/
    ├── valenx-adapter-openfoam/        # Subprocess, GPL isolation
    │   ├── Cargo.toml
    │   └── src/lib.rs
    ├── valenx-adapter-gmsh/            # Subprocess; also has optional C-API cfg-gate
    ├── valenx-adapter-freecad/         # Dynamic-linked (LGPL via OCCT)
    ├── valenx-adapter-cantera/         # Bundled (BSD)
    ├── valenx-adapter-code-aster/      # Subprocess
    ├── valenx-adapter-calculix/        # Subprocess
    ├── valenx-adapter-elmer/           # Dynamic-linked (LGPL)
    ├── valenx-adapter-su2/             # Dynamic-linked (LGPL)
    ├── valenx-adapter-preciyce/        # Dynamic-linked (LGPL), for coupling
    └── ...
```

Every adapter crate is independent, has its own `Cargo.toml`, and can
be built in isolation. Adapters don't depend on each other; they all
depend on `valenx-core::adapter`.

### The `Adapter` trait in full

```rust
use std::path::{Path, PathBuf};
use std::sync::Arc;

pub trait Adapter: Send + Sync {
    /// Static metadata. Called once at registration.
    fn info(&self) -> AdapterInfo;

    /// Is the tool usable on this system?
    fn probe(&self) -> Result<ProbeReport, AdapterError>;

    /// Translate a canonical Case into the tool's native input deck
    /// inside `workdir`. Must be deterministic.
    fn prepare(
        &self,
        case: &Case,
        workdir: &Path,
    ) -> Result<PreparedJob, AdapterError>;

    /// Run the tool. Long-running. Must respect ctx.check_cancel().
    fn run(
        &self,
        job: &PreparedJob,
        ctx: &mut RunContext,
    ) -> Result<RunReport, AdapterError>;

    /// Parse outputs into the canonical Results.
    fn collect(
        &self,
        job: &PreparedJob,
    ) -> Result<Results, AdapterError>;

    /// Optional: validate a case before queueing it, fast.
    /// Default implementation returns Ok(()).
    fn validate(&self, case: &Case) -> Result<(), AdapterError> {
        let _ = case;
        Ok(())
    }

    /// Optional: provide a capability map so the UI can show / hide features.
    /// Default returns Capabilities::default().
    fn capabilities(&self) -> Capabilities {
        Capabilities::default()
    }
}
```

### `AdapterInfo`

```rust
pub struct AdapterInfo {
    pub id: &'static str,                    // "openfoam", "gmsh", ...
    pub display_name: &'static str,          // "OpenFOAM"
    pub version_range: SemVerRange,          // "v2306..v2506"
    pub physics: &'static [Physics],         // &[Physics::CFD]
    pub license_mode: LicenseMode,
    pub tool_license: &'static str,          // SPDX identifier
    pub docs_url: &'static str,
    pub homepage_url: &'static str,
}
```

### `ProbeReport`

```rust
pub struct ProbeReport {
    pub ok: bool,
    pub found_version: Option<SemVer>,
    pub binary_path: Option<PathBuf>,
    pub warnings: Vec<String>,
    pub required_env: Vec<(&'static str, String)>,   // e.g. FOAM_ETC=/opt/openfoam/etc
}
```

### `PreparedJob`

```rust
pub struct PreparedJob {
    pub workdir: PathBuf,
    pub native_command: Vec<OsString>,       // argv for subprocess adapters
    pub environment: Vec<(OsString, OsString)>,
    pub estimated_runtime: Option<Duration>,
    pub kill_on_drop: bool,
}
```

### `RunReport` and `Results`

```rust
pub struct RunReport {
    pub exit_code: i32,
    pub wall_time: Duration,
    pub converged: Option<bool>,             // Some(true) / Some(false) / None (unknown)
    pub residual_history: Vec<ResidualSample>,
    pub warnings: Vec<String>,
}

pub struct Results {
    pub fields: Vec<FieldRef>,               // pointers to VTK / CGNS files
    pub scalars: HashMap<String, f64>,       // drag coefficient, etc.
    pub artifacts: Vec<Artifact>,            // forces.csv, probes.dat, ...
    pub manifest_path: PathBuf,              // cases/<name>/results/manifest.toml
}
```

### Determinism requirements

`prepare()` must be **deterministic** — given the same `Case`, the
files written to `workdir` must be byte-identical. This is how we make
`config_hash` (from RFC 0001) work. In practice this means:

- Sort map entries before writing TOML
- Use a fixed line ending
- Don't embed timestamps in inputs (timestamps go in `results/manifest.toml`)
- Don't use `HashMap` iteration order — use `BTreeMap` or explicitly sorted keys

We ship a test harness: `valenx-adapter-test` which verifies
determinism by running `prepare()` twice and diffing.

### Error taxonomy recap

```rust
pub enum AdapterError {
    ToolNotInstalled { name: &'static str, hint: String },
    ToolVersionMismatch { expected: SemVerRange, found: SemVer },
    InvalidCase { case_path: PathBuf, reason: String },
    Translate(TranslateError),
    Run { exit_code: i32, stderr: String, phase: RunPhase },
    ParseOutput { file: PathBuf, reason: String },
    Cancelled,
    IO(std::io::Error),
    Other(anyhow::Error),
}

pub enum RunPhase {
    Startup, MeshRead, Solve, Output, Shutdown,
}

#[derive(Debug, thiserror::Error)]
pub enum TranslateError {
    #[error("unsupported feature {feature:?}")]
    Unsupported { feature: String },
    #[error("required field {field:?} missing")]
    MissingField { field: String },
    #[error("value {value} out of range [{min}, {max}] for {field}")]
    OutOfRange { field: String, value: f64, min: f64, max: f64 },
}
```

### Subprocess-mode rules (hard constraints)

A subprocess-mode adapter:

- **Must not** list the GPL tool in `Cargo.toml` as a dependency
- **Must not** `dlopen` or `LoadLibrary` the tool's shared objects
- **Must** communicate only via: standard input/output, arguments, and
  files in its `workdir`
- **Must** kill its child process if dropped
- **Must** treat everything the tool prints as untrusted text (no
  unchecked parsing that could OOM)

Enforced by:

- A CI check (`cargo deny` policy file per adapter) that denies any
  link against forbidden-license crates
- Code review — subprocess adapters go through extra scrutiny
- `#![forbid(unsafe_code)]` at the adapter crate level unless
  specifically justified

### Tool version pinning

Each adapter crate declares its supported version range:

```rust
fn info(&self) -> AdapterInfo {
    AdapterInfo {
        version_range: SemVerRange::from("v2306..v2506"),
        ...
    }
}
```

At startup, Valenx runs `probe()` on every registered adapter. The
`tools.lock` entry's version must fall inside the adapter's range; if
not, the adapter is marked unavailable and the UI explains why.

### Cross-platform subprocess invocation

Subprocess adapters must handle:

- **Windows** — child processes need a job object so killing the parent
  kills children (uses `valenx-core::process::Job`)
- **macOS** — code-signing sometimes blocks spawning bundled binaries;
  probe detects and reports this
- **Linux** — typically straightforward; containerized runs (Docker)
  need PID 1 awareness for signal handling

Helpers live in `valenx-core::process`, not in each adapter.

### Progress conventions

- Progress `pct` is 0.0 to 100.0; inclusive bounds
- Report at least every 2 seconds or every 5% change, whichever first
- Messages are user-facing text, < 80 chars
- Stage boundaries are reported as `progress(x, "Meshing...")`,
  `progress(y, "Solving iteration 1200/2000")`, etc.

---

## Drawbacks

- **Contract is wide.** Five phases plus capabilities/validate means
  adapters take real work to write correctly. Counterpoint: if we had
  fewer phases, we'd force them all to cram logic into one function.
- **Enforcing license-mode boundaries is partly honor-system.** `cargo
  deny` catches direct deps but not transitive weirdness; code review
  handles the rest.
- **Determinism is hard.** Floating-point output in intermediate files,
  timestamps in headers, locale-dependent number formatting — many
  gotchas. We accept the test cost.
- **Error taxonomy is opinionated.** Some real-world failures will
  squish awkwardly into `Run { phase, exit_code }`. Over time we'll
  extend with more specific variants.

---

## Rationale and alternatives

**One trait per physics (CfdAdapter, FeaAdapter, …):**
Rejected. Many tools span physics (OpenFOAM does CFD + CHT; FreeCAD
does CAD + meshing via plugins). A single trait with a `Physics` flag
in `AdapterInfo` is simpler.

**Adapters as plugins (WASM):**
Rejected for built-in adapters. Performance and FFI complexity aren't
worth it for in-tree code. WASM plugins are a separate concept,
covered in RFC 0003 — for user-contributed, sandboxed extensions.

**gRPC / IPC instead of subprocess stdout parsing:**
Considered. Most GPL solvers don't expose an IPC interface, so we'd be
patching them upstream. For tools that *do* expose IPC (Cantera's C
API, MuJoCo's C API), the `DynamicLinked` mode already covers it.

**No error taxonomy, just `anyhow`:**
Rejected. The UI needs to react differently to different errors;
string matching `anyhow::Error` messages is brittle.

---

## Prior art

- **Salome-Meca** — `salome_run` is the closest analog: a launcher
  that drives Code_Aster, Code_Saturne, others. Lessons learned: each
  tool gets a separate adapter module; Salome's problem is that
  adapters leaked solver concepts up into the UI. We avoid that with
  the canonical `Case`/`Results` types.
- **CAELinux** — bundles many solvers; no unified API, just a Linux
  live CD. Not a reference for the contract, just for the bundling.
- **Rust Analyzer's proc-macro server** — subprocess isolation for
  *compiled* code. Demonstrated that subprocess can be fast enough for
  interactive loops.
- **Tauri's IPC layer** — subprocess with structured JSON messaging.
  Works well for their use case; overkill for ours.

---

## Unresolved questions

- **Where do adapter-specific config options live?** Global (machine-wide
  preferences), per-user, or per-project? Current plan: per-project,
  in `project.toml` under `[adapter.<id>]`. Revisit after a year of use.
- **How to handle adapters for tools with exotic license terms** — e.g.,
  a tool that's open-source but requires an academic-use-only license
  for certain modules. Defer to a case-by-case RFC.
- **Remote adapters** — running the solver on a cluster or cloud. The
  contract above assumes local execution; remote mode is a future RFC
  that subclasses or wraps this one.
- **Adapter sandboxing** — even in-tree adapters could benefit from
  some OS-level sandbox (seccomp, AppContainer). Defer.

---

## Future possibilities

- **Adapter marketplace** — once the WASM plugin API (RFC 0003) is
  stable, users could ship third-party adapters. The built-in crates
  become the reference implementation.
- **Adapter benchmarking harness** — run a fixed validation case
  through every CFD adapter, compare results and runtime
- **Auto-generated adapters** — for tools with OpenAPI or similar
  schemas, a code generator could produce the Rust stub
- **Cross-physics orchestration** (FSI, CHT) via preCICE is already
  scoped as its own adapter; this RFC's contract supports it.
