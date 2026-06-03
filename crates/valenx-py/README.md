# valenx-py

Phase 11 of the FreeCAD-parity roadmap: a Python scripting bridge that
exposes the `valenx-cad`, `valenx-sketch`, and `valenx-feature-tree`
crates as a `valenx` Python module via PyO3.

## Scope

Honestly small but real, not a full FreeCAD API clone. v1 covers the
authoring path a typical user reaches for from a script:

- `valenx.cad` — primitives (`box`, `cylinder`, `sphere`, `cone`,
  `torus`, `prism`), booleans (`union`, `difference`, `intersection`),
  tessellation (`solid_to_mesh`)
- `valenx.sketch` — 2D parametric sketcher with constraints and
  extrude-back-to-solid
- `valenx.feature_tree` — ordered Pad / Pocket / Revolve / Mirror /
  Pattern / Fillet / Chamfer / ImportedSolid tree with `replay()`
- `valenx.mesh` — read-only mesh handle returned by tessellation

## Building the wheel

`valenx-py` ships as both an `rlib` (so the workspace's `cargo check`
keeps working) and a `cdylib` (the actual loadable Python extension).
The recommended toolchain is [`maturin`](https://www.maturin.rs/):

```bash
pip install maturin
cd crates/valenx-py
maturin develop --release
```

This produces a `valenx_py` extension importable as `import valenx`
(the `[lib]` name maps to the Python module name).

## Usage example

```python
import valenx

# Build a punched cube via the cad submodule.
cube = valenx.cad.box(10, 10, 10)
hole = valenx.cad.sphere(4).translated(5, 5, 5)
punched = valenx.cad.difference(cube, hole)
print(repr(punched))                       # Solid(faces=..., edges=..., vertices=...)

# Tessellate for downstream rendering / STL.
mesh = valenx.cad.solid_to_mesh(punched, 0.5)
print(f"{mesh.triangle_count} triangles")

# Build a parametric sketch and extrude.
sk = valenx.sketch.Sketch()
a = sk.add_point(0, 0)
b = sk.add_point(1, 0)
c = sk.add_point(1, 1)
d = sk.add_point(0, 1)
sk.add_line(a, b)
sk.add_line(b, c)
sk.add_line(c, d)
sk.add_line(d, a)
sk.add_constraint("horizontal", {"line": 5})  # line a-b
sk.add_constraint("vertical",   {"line": 6})  # line b-c
prism = sk.extrude(2.0)

# Feature tree: pad a profile, pocket a hole, replay.
tree = valenx.feature_tree.FeatureTree()
sketch_idx = tree.add_sketch(sk)
pad = tree.add_feature(
    "pad",
    {"sketch": sketch_idx, "depth": 2.0, "direction_positive": True},
    "Base Pad",
)
solid = tree.replay()
```

## v1 limitations

- **No fillet / chamfer wrappers in `valenx.cad`.** The Phase 3
  mesh-domain fillet doesn't round-trip through booleans cleanly;
  expose later under `valenx.mesh_ops` when the BRep fillet ships.
- **No persistence wrappers.** RON serialisation lives on the Rust
  side. Phase 11.5 will add `valenx.cad.save_solid` /
  `valenx.feature_tree.save_tree` etc.
- **No solver invocation from Python.** Constraints are added but
  `valenx_sketch::solver::solve` isn't yet wrapped. Extrude works on
  the as-typed coordinates; deferred to Phase 11.5.

## Why submodules

Mirrors FreeCAD's `Part`, `Sketcher`, `PartDesign` split — users
coming from FreeCAD already know to reach for `valenx.cad.box(...)`,
not a flat `valenx.box(...)`.

## Smoke test

`tests/smoke.rs` is gated behind the `embed-python` feature so it
only runs when explicitly invoked with a linkable Python on the build
machine:

```bash
cargo test -p valenx-py --features embed-python -- --ignored
```

Without `embed-python`, `cargo test -p valenx-py` runs zero tests
(the smoke test is `#[ignore]`-tagged).
