# RFC 0004: Results and Fields Data Model

- **Status:** Accepted (initial design; revisions expected as physics verticals land)
- **Author(s):** BDFL
- **Created:** 2026-04-23
- **Discussion PR:** (this commit)
- **Tracking issue:** TBD

---

## Summary

Define the canonical in-memory and on-disk representation of
simulation **Results** — the fields (scalar / vector / tensor),
derived scalars, provenance, and units — that every adapter emits
and every downstream consumer (viewer, plots, tables, report
generator, export) reads. This is the shared vocabulary that lets
one UI speak across 75 solvers.

---

## Motivation

Without a canonical results type, every adapter emits its own
shape, every UI pane reads N different shapes, and every export
path is custom per solver. That's how commercial suites end up
with brittle, duplicated rendering and export code. Valenx
commits to the opposite: one schema, every adapter converts into
it, every consumer reads it.

Design requirements:

1. **Typed fields.** Pressure, velocity, stress, temperature — we
   know what kind of data each is, its rank, its units.
2. **Dimension-agnostic.** The same types work for 1-D (battery
   cell), 2-D (cross-sections, slabs), and 3-D (full CFD / FEA).
3. **Unit-aware.** Every numeric value has units; conversions are
   explicit; unit mismatches are errors, not silent bugs.
4. **Streamable.** Fields can be hundreds of MB; the model must
   support lazy loading / memory-mapped access / on-disk chunks.
5. **Serializable to open formats.** VTK / CGNS / HDF5 are the
   non-negotiable export paths.
6. **Derived scalars and metadata cheap to query.** "What was the
   drag coefficient?" is answerable without loading any field.
7. **Reproducible.** Every `Results` object carries enough
   provenance to explain where every value came from — which
   adapter, which solver version, which input hash.

Anti-requirements:

- Not tied to any one solver's convention (OpenFOAM's `fields`
  layout, CalculiX `frd`, etc.)
- Not a full CAE data framework — we borrow from CGNS conceptually
  but don't adopt it wholesale; too much scope for v1
- Not a database — just typed containers with clear serialization

---

## Guide-level explanation

A `Results` object is what every adapter's `collect()` returns (per
RFC 0002). It holds:

```rust
pub struct Results {
    pub meta: ResultMeta,
    pub fields: FieldCatalog,
    pub scalars: ScalarCatalog,
    pub artifacts: Vec<Artifact>,
    pub provenance: Provenance,
}
```

- **`meta`** — who produced this, when, what case hash, mesh hash
- **`fields`** — the big-data fields (pressure, velocity, stress)
- **`scalars`** — small derived quantities (drag coefficient,
  max stress, wall-time, iteration count)
- **`artifacts`** — pointers to raw files on disk (VTK, CSVs, raw
  solver outputs) for power-user inspection
- **`provenance`** — chain of adapter + solver versions that
  produced this

### Fields

A **field** is a named array of values defined on a mesh region,
with a rank (scalar / vector / tensor) and units.

```rust
pub struct Field {
    pub name: String,             // "pressure", "velocity_x", "stress_von_mises"
    pub kind: FieldKind,          // Scalar | Vector{dim} | Tensor{rank}
    pub location: Location,       // OnNode | OnCell | OnFace | OnEdge
    pub region: RegionRef,        // mesh region ID (body, boundary, volume subset)
    pub units: Units,             // SI dimension tuple
    pub time: TimeKey,            // Steady | Transient{t} | Iteration{k}
    pub data: FieldData,          // eager Vec or lazy handle
    pub range: Option<(f64, f64)>, // cached min/max for viz
}
```

Five locations (node / cell / face / edge / region-constant)
cover 99% of practical CAE data.

### Scalars

Small derived values, indexed by name:

```rust
pub struct ScalarRecord {
    pub name: String,             // "cd", "cl", "max_stress", "frequency_hz"
    pub value: f64,
    pub units: Units,
    pub time: TimeKey,
    pub description: Option<String>,
    pub source: ScalarSource,     // Computed | Extracted | UserDefined
}
```

### Units

Represented as a signed 7-tuple of base SI dimensions (length,
mass, time, current, temperature, amount-of-substance, luminous
intensity) plus a scale factor for display:

```rust
pub struct Units {
    pub dims: [i8; 7],            // powers of base SI
    pub scale: f64,               // for display (km vs m)
    pub display: &'static str,    // "Pa", "m/s", "K", "mol/L"
}
```

Arithmetic on `Units` is checked — `Pa + m/s` is a compile-ish
error (runtime, but caught before any numerics happen). A common
constants table (`units::PASCAL`, `units::METER_PER_SECOND`) is
shipped.

### Time

```rust
pub enum TimeKey {
    Steady,
    Iteration(u64),                      // SIMPLE/coupled iteration count
    Time { value: f64, units: Units },   // for transient runs
}
```

Transient results store a time series — many `Field` records each
with a different `TimeKey::Time`, keyed by a shared `name`.

### Provenance

```rust
pub struct Provenance {
    pub adapter: String,                // "openfoam"
    pub adapter_version: SemVer,
    pub tool: String,                   // "OpenFOAM-v2406"
    pub tool_version: String,
    pub case_hash: Sha256,
    pub mesh_hash: Sha256,
    pub input_hash: Sha256,
    pub tools_lock_hash: Sha256,
    pub run_id: Uuid,
    pub wall_time: Duration,
    pub completed_at: DateTime<Utc>,
    pub ancestors: Vec<ProvenanceRef>,  // for derived results
}
```

The `ancestors` list is how we track derived results: a
post-processed field carries a reference to the raw field it came
from. Scientific reproducibility demands this.

---

## Reference-level explanation

### Crate layout

Lives in `crates/valenx-fields/`:

```
valenx-fields/
├── Cargo.toml
├── src/
│   ├── lib.rs
│   ├── field.rs              // Field, FieldKind, Location, FieldData
│   ├── scalar.rs             // ScalarRecord, ScalarCatalog
│   ├── catalog.rs            // FieldCatalog (indexed by name + time)
│   ├── units.rs              // Units, arithmetic, canonical constants
│   ├── time.rs               // TimeKey, time-series queries
│   ├── region.rs             // RegionRef, region resolution against Mesh
│   ├── provenance.rs         // Provenance, ancestors, hashing
│   ├── meta.rs               // ResultMeta
│   ├── artifact.rs           // Artifact (path + format tag + checksum)
│   ├── results.rs            // Results struct, queries
│   ├── data.rs               // FieldData: eager + lazy variants
│   └── io/
│       ├── vtk.rs            // VTK XML legacy + unstructured writers
│       ├── cgns.rs           // CGNS for CFD
│       ├── hdf5.rs           // generic HDF5
│       └── json.rs           // metadata-only JSON for CI gates
└── tests/                    // round-trip + unit-arithmetic tests
```

### `FieldData` — eager vs. lazy

```rust
pub enum FieldData {
    Dense(Array),                 // ndarray::ArrayD<f64> in memory
    Mmap(MmapHandle),             // memory-mapped file-backed
    Chunked(ChunkedHandle),       // loaded on demand per chunk
    External(PathBuf),             // not loaded; user decides when
}
```

Consumers use trait methods (`.len()`, `.slice(...)`,
`.as_slice()`) that work uniformly across the variants. A
memory-constrained visualizer can hold `Mmap` fields and OS page
cache handles paging.

### `FieldCatalog` indexing

```rust
impl FieldCatalog {
    pub fn by_name(&self, name: &str) -> Iter<'_, Field>;
    pub fn at_time(&self, name: &str, t: TimeKey) -> Option<&Field>;
    pub fn names(&self) -> impl Iterator<Item = &str>;
    pub fn time_series(&self, name: &str) -> Vec<TimeKey>;
    pub fn insert(&mut self, field: Field);
}
```

Typical usage:

```rust
let results: Results = adapter.collect(&job)?;
let cd = results.scalars.get("cd")?;                 // drag coeff
let p  = results.fields.at_time("pressure", TimeKey::Steady)?;
for t in results.fields.time_series("velocity") {
    let v = results.fields.at_time("velocity", t)?;
    // plot, render, export
}
```

### Unit arithmetic

```rust
impl std::ops::Mul for Units { ... }  // dims add
impl std::ops::Div for Units { ... }  // dims subtract
impl Units {
    pub fn is_compatible(&self, other: &Units) -> bool;
    pub fn convert(&self, to: &Units) -> Option<f64>; // scale factor
}
```

A library of constants:

```rust
pub mod units {
    pub const METER: Units;
    pub const KILOGRAM: Units;
    pub const SECOND: Units;
    pub const KELVIN: Units;
    pub const PASCAL: Units;               // kg·m⁻¹·s⁻²
    pub const METER_PER_SECOND: Units;
    pub const NEWTON: Units;
    // ... ~60 common ones
}
```

### Serialization

**VTK.** Default for CFD results. We use `vtkio` for legacy `.vtk`
and XML `.vtu` / `.vtm`. Every `Field` maps to a VTK
`PointData` / `CellData` array; `FieldCatalog` maps to a VTK
multiblock.

**CGNS.** For CFD interchange with other suites. HDF5-backed,
industry-standard, richer than VTK for boundary metadata.

**HDF5.** Generic container for anything that doesn't fit VTK or
CGNS — FEA modal results, battery state-of-charge over cycles,
chemistry time histories.

**JSON.** Metadata-only mode (everything but field arrays). For
CI assertions, diffing two runs, headless validation.

Every format goes through a trait:

```rust
pub trait ResultsWriter {
    fn write(&self, results: &Results, path: &Path) -> Result<(), IoError>;
    fn format_name(&self) -> &'static str;
    fn supported_fields(&self) -> FieldCapabilityMask;
}
```

And the reverse for round-trip tests:

```rust
pub trait ResultsReader {
    fn read(&self, path: &Path) -> Result<Results, IoError>;
}
```

### On-disk layout inside a `.valenx` project

Per RFC 0001, `cases/<name>/results/`:

```
results/
├── manifest.toml         # Results-format v1 manifest
├── scalars.csv           # human-readable scalar dump
├── provenance.json
├── fields/
│   ├── pressure-steady.vtu
│   ├── velocity-steady.vtu
│   └── temperature-t=1.00.vtu
└── artifacts/
    ├── solver.log
    └── residual-history.csv
```

`manifest.toml` is the authoritative index — it lists every field
and scalar with its metadata, avoiding a full directory scan.

### Derived results

A post-processing operation (e.g., computing Q-criterion from
velocity gradient) produces a new `Field` whose `Provenance.ancestors`
points at the source field. The viewer and report generator can
follow that chain all the way back to the solver run.

### Comparisons and diffs

```rust
impl Results {
    pub fn compare(&self, other: &Results) -> ComparisonReport;
}
```

Used by the "A/B results" UI (Year 3+) and by CI validation cases
(compare today's run against the pinned reference). The report
includes per-field L² and L∞ norms, per-scalar absolute and
relative differences, and highlights fields present in one but not
the other.

---

## Drawbacks

- **More scope than a single-RFC minimum.** We could ship with
  just `HashMap<String, Vec<f64>>` and iterate. But that leaks
  into every consumer and becomes painful to change later. Paying
  the design cost now is cheaper than paying it in year 5.
- **Unit system is non-trivial.** SI-tuple dimension algebra is a
  known pattern but adds about 1000 LOC. Worth it; unit bugs are
  the worst bugs in engineering code.
- **Lazy `FieldData` variants add complexity** to consumers that
  need random access. Mitigated by the trait-based accessor API.
- **CGNS / HDF5 dependencies** pull in C code. HDF5 in particular
  is a big dep. Made optional via Cargo features; only built when
  the user wants those formats.

---

## Rationale and alternatives

**Adopt CGNS wholesale.**
Rejected. CGNS is CFD-centric; doesn't fit FEA modal, chemistry
time series, battery results. Our schema is inspired by CGNS's
node/cell-data distinction and boundary metadata but stays
physics-agnostic.

**`HashMap<String, Vec<f64>>` everywhere.**
Rejected. No types, no units, no provenance. Will become a
maintenance nightmare within a year.

**Use Apache Arrow.**
Considered. Strong type system, good columnar perf, first-class
Python interop. But Arrow isn't designed for meshes; fields on a
mesh are not columnar. We stay with domain-appropriate containers.

**Store everything in HDF5 by default.**
Rejected. HDF5 is a large dep and opaque on disk. VTK + CSV +
TOML per-case is diff-friendly, inspectable, and works offline.

---

## Prior art

- **CGNS (CFD General Notation System)** — the industry standard
  for CFD. Our node/cell-data separation and boundary-metadata
  approach is directly inspired.
- **VTK** — the de-facto visualization format. We always export to
  it.
- **OpenFOAM's `fields` directory** — per-time-step, per-field
  files. We mirror the principle (time series as separate files).
- **Ansys CDB, Abaqus ODB, Nastran OP2** — proprietary binary
  results formats. We interoperate via export-only, not import.
- **xarray (Python)** — labeled n-d arrays with metadata. Close to
  what we want but Python-specific; Rust's `ndarray` + our wrapper
  types achieve similar ergonomics.
- **Apache Arrow** — columnar store, excellent for tabular; not
  right for meshed fields.

---

## Unresolved questions

- **Complex-valued fields.** Frequency-domain EM, eigenmode
  analysis. Probably a separate `ComplexField` variant; deferred
  until EM adapter lands.
- **Adaptive mesh refinement.** Results on meshes that changed
  mid-run. The `region` reference handles this via mesh-hash
  binding but the UX for it is undefined.
- **Per-material-point data** (plasticity state, damage
  variables). Extension of `Location`; not yet specified.
- **Streaming results during a running solve.** Do we mutate the
  `Results` as the solver reports progress, or snapshot on
  completion? Current plan: on completion; progress comes through
  `RunContext` in RFC 0002.

---

## Future possibilities

- **Lineage graph UI** — visualize the provenance chain across
  multiple derived results; useful for debugging and reports
- **Results diffing in the browser tree** — two runs side by side
  with highlighted scalar differences
- **Probe-time plots** — click a point in the viewport, see its
  field value over time as a sparkline in a hover card
- **Streaming / incremental load** for >10 GB results on HPC
  systems
- **Integration with Common Data Model (CDM)** for climate /
  atmospheric data when those adapters land
