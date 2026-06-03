# RFC 0012 — ML training-data export

| | |
|---|---|
| **Status** | Draft |
| **Phase target** | 13 (ML / surrogate models) |
| **Type** | Architecture |
| **Author** | Maintainers |
| **Discussion** | open in PR — not yet merged |
| **Related** | RFC 0011 (parameter sweeps), RFC 0004 (results + fields) |

## Summary

Define a stable, framework-agnostic export format for turning a sweep
of completed Valenx runs (RFC 0011) into a tensor dataset for
training neural-network surrogates, regression models, or any other
ML pipeline. Output is a single directory of numpy `.npz` files
(plus a JSON manifest) that PyTorch / TensorFlow / scikit-learn can
load with three lines of code.

This RFC is intentionally framework-neutral. We don't pick a Python
framework; we don't ship a training loop; we don't define a model
zoo. Those are the user's call. Valenx's job ends at "here's a
clean dataset — go train."

## Motivation

The point of running 200 parameter-sweep cases isn't 200 runs — it's
a model that can predict the 201st without running it. Today turning
Valenx run output into ML-trainable tensors requires:

- Iterating each run's `Results.fields` / `Results.scalars` by hand.
- Dealing with per-adapter quirks (OpenFOAM stores velocity in a
  flat `[ux, uy, uz, ux, uy, uz, …]` buffer; CalculiX stores
  displacement the same way but tagged differently).
- Re-aligning fields onto a common grid for training (different
  cases mesh differently if they involve geometry sweeps).
- Splitting train/val/test, normalising, batching.

We don't solve all of those, but we get the user as close to "a
torch.utils.data.Dataset" as we can while staying honest about
units, provenance, and field semantics.

## Design

### CLI / app entry point

```bash
valenx export --sweep <sweep-dir> --format npz --output <out-dir>
```

Or via Run menu → "Export training dataset…". Either way, the user
points at a sweep result directory and gets back a tensor dump.

### Output layout

```
<out-dir>/
├── manifest.json
├── inputs/
│   ├── sample_0001.npz
│   ├── sample_0002.npz
│   └── …
└── outputs/
    ├── sample_0001.npz
    ├── sample_0002.npz
    └── …
```

`manifest.json`:

```json
{
  "valenx_export_version": "0.1",
  "sweep_source": "../sweep-airfoil-aoa-2026-04-25/",
  "sample_count": 200,
  "inputs": {
    "schema": [
      { "name": "aoa_deg",       "shape": [1], "units": "deg",   "dtype": "f32" },
      { "name": "inlet_velocity",  "shape": [3], "units": "m/s", "dtype": "f32" }
    ]
  },
  "outputs": {
    "schema": [
      { "name": "drag_coefficient", "shape": [1], "units": "1",       "dtype": "f32" },
      { "name": "pressure",         "shape": [N], "units": "Pa",      "dtype": "f32",
        "location": "OnNode", "mesh_id": "common-airfoil" }
    ],
    "common_mesh": {
      "id": "common-airfoil",
      "node_count": 12345,
      "source": "interpolated to baseline mesh"
    }
  }
}
```

Each `sample_NNNN.npz` is a numpy zip-of-arrays whose keys match
the manifest's `name` entries. Loading from PyTorch:

```python
import numpy as np
data = np.load("sample_0001.npz")
aoa = data["aoa_deg"]      # shape (1,)
p   = data["pressure"]     # shape (N,)
```

### Field interpolation

When samples have different meshes (geometry sweep), the exporter
interpolates each output Field onto a common mesh declared by the
user. Default: the mesh from the first run. Configurable via
`--common-mesh <case>`. The interpolation algorithm is
nearest-neighbour for scalar / vector fields on point data;
proper RBF / k-d-tree-weighted interpolation is a follow-up RFC.

When samples share a mesh (parameter-only sweep, no geometry
change), the interpolation step is a no-op and we emit the field
verbatim.

### Scalar inputs

`Results.scalars` from each run becomes the input vector. Sweep-
generated parameters (the `[[sweep.parameter]]` declarations from
RFC 0011) are added as inputs automatically; user can opt extra
scalars in/out via the export config.

### Provenance

The manifest carries the sweep's git commit hash, the Valenx
version, every adapter version, and a SHA-256 of every source
case.toml. No way to "lose" the connection between training
data and the simulations that produced it.

### Train/val/test split

A simple `--split 0.7,0.15,0.15` flag produces three subdirectories.
Sample assignment is deterministic (hash of sample_id) so a re-export
gives the same split. Stratification by parameter ranges is a
follow-up.

## Drawbacks

- npz isn't the only ML serialisation format users want. PyTorch
  prefers `.pt`; TensorFlow prefers `.tfrecord`. We ship `.npz`
  as the lowest-common-denominator (everyone can read it); per-
  framework converters can land later.
- Interpolation onto a common mesh loses information when meshes
  are very different. The naive nearest-neighbour algo will give
  visible artefacts in surrogate predictions for fine-feature
  geometry sweeps. A real RBF interpolant is needed for production
  use.
- Dataset versioning is on the user — if the export schema changes,
  models trained against the old schema break. We provide
  `valenx_export_version` in the manifest so loaders can detect
  this; we don't auto-migrate.

## Alternatives considered

- **Direct PyTorch DataLoader integration** — would tie us to
  PyTorch's API, lock in CUDA assumptions, force a Python dep on
  the Valenx binary. Hard no.
- **Parquet / Arrow output** — modern, columnar, language-agnostic,
  but most ML frameworks don't have first-class loaders. Could
  ship as a second format (`--format arrow`) once `--format npz`
  is established.
- **HDF5 dump of the Results catalog verbatim** — works but burdens
  the user with re-aligning per-sample meshes themselves; the
  point of the exporter is to do that for them.

## Migration path

Phase 13.0 (this RFC): land `valenx-export` crate with the npz
writer + manifest schema + common-mesh interpolation (nearest-
neighbour only).

Phase 13.1: real RBF interpolant + per-framework converters
(`--format pt`, `--format tfrecord`).

Phase 13.2: a "training quickstart" Python module under
`tools/train-quickstart/` that demonstrates loading the exported
data into PyTorch + a small MLP for regression. Documentation only,
not a Valenx-shipped artifact.

Phase 13.3: surrogate-assisted optimization — feed the trained
surrogate back into RFC 0011's `GradientDescentOptimizer` so
expensive simulations only run when the surrogate's uncertainty
is high.

## Open questions

1. **Time-series export.** Today's design assumes one snapshot per
   sample. Transient runs have hundreds; do we export each
   timestep as its own sample (huge dataset), or pick a single
   "interesting" timestep (which one?), or emit a per-sample
   time-series tensor (changes the schema)?
2. **Categorical inputs.** Parameter sweeps over discrete choices
   (`turbulence = ["kEpsilon", "kOmegaSST"]`) need one-hot encoding.
   Do we emit the raw string + a separate encoder, or pre-encode?
3. **Missing values.** A failed run can't contribute output but its
   inputs are still meaningful for "predicting failure." Out of
   scope for v0; the exporter skips failed runs entirely today.
4. **Streaming export for huge sweeps.** A 10,000-run sweep produces
   GB-scale output. Today's design loads each run into memory
   before writing the npz; that's fine to ~1k samples but not
   beyond. Streaming variant is a follow-up.

## References

- RFC 0011 — parameter sweeps (the upstream that makes this RFC's
  output worth having).
- RFC 0004 — results + fields (the canonical types we're exporting).
- numpy `.npz` format spec.
- PyTorch `torch.utils.data.Dataset` reference.
- HDF5 SciML — competing format we considered.
