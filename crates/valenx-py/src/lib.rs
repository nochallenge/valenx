//! # valenx-py
//!
//! Python scripting bridge for Valenx — Phase 11 of the FreeCAD-parity
//! roadmap. Exposes the [`valenx_cad`], [`valenx_sketch`], and
//! [`valenx_feature_tree`] crates as a single `valenx` Python module via
//! PyO3.
//!
//! Scope
//! =====
//!
//! Honestly small but real, not a full FreeCAD API clone:
//!
//! - **cad** — primitives (`box`, `cylinder`, `sphere`, `cone`,
//!   `torus`, `prism`), booleans (`union`, `difference`,
//!   `intersection`), tessellation (`solid_to_mesh`). Mirrors
//!   `valenx_cad`'s public API one-for-one.
//! - **sketch** — [`crate::sketch::PySketch`] wraps
//!   [`valenx_sketch::Sketch`] with `add_point`, `add_line`,
//!   `add_circle`, `add_constraint` (constraint kind as a string),
//!   and `extrude`.
//! - **feature_tree** — [`crate::feature_tree::PyFeatureTree`] wraps
//!   [`valenx_feature_tree::FeatureTree`] with `add_sketch`,
//!   `add_feature` (kind + params as a dict), and `replay`.
//!
//! ## Usage from Python
//!
//! ```python
//! import valenx
//!
//! box    = valenx.cad.box(10, 20, 30)
//! sphere = valenx.cad.sphere(5).translated(5, 10, 15)
//! result = valenx.cad.difference(box, sphere)
//! mesh   = valenx.cad.solid_to_mesh(result, 0.5)
//! print(f"mesh has {mesh.triangle_count} triangles")
//! ```
//!
//! ## Why submodules
//!
//! Mirrors FreeCAD's `Part`, `Sketcher`, `PartDesign` split — Python
//! users coming from FreeCAD already know to reach for
//! `valenx.cad.box(...)`, not a flat `valenx.box(...)`. Each submodule
//! is a [`pyo3::types::PyModule`] registered under the top-level
//! `valenx` module.
//!
//! ## v1 limitations
//!
//! - **No fillet / chamfer wrappers in `valenx.cad`.** Phase 3 lives
//!   in `valenx_fillet` and writes mesh-backed solids that don't round
//!   trip cleanly through booleans; expose later under a
//!   `valenx.mesh_ops` submodule when the BRep fillet ships.
//! - **No persistence wrappers.** `SketchFile` / `FeatureTreeFile` /
//!   `ValenxProject` RON envelopes aren't bridged — Python users can
//!   construct the trees and immediately tessellate without writing
//!   project files for now.
//! - **No constraint solver invocation from Python.** Users build
//!   constraints and add them to the sketch, but
//!   [`valenx_sketch::solver::solve`] is not yet wrapped. Phase 11.5
//!   adds it.

#![forbid(unsafe_code)]
#![warn(missing_docs)]
// PyO3 0.22's `#[pymethods]` proc macro emits `.into()` calls on the
// return value of every fallible method, which clippy flags as
// `useless_conversion` when the function already returns `PyResult<T>`
// (the macro-emitted `Into::<PyErr>::into` is a no-op for an
// already-`PyErr` value). The conversion lives inside the macro
// expansion and is not addressable from user code, so the lint is
// silenced crate-wide. See pyo3#3743.
#![allow(clippy::useless_conversion)]

pub mod cad;
pub mod error;
pub mod feature_tree;
pub mod mesh;
pub mod sketch;

use pyo3::prelude::*;

/// PyO3 entry point. `maturin develop` builds the cdylib and Python
/// loads it as the `valenx_py` module; users `import valenx_py as valenx`
/// or set `from valenx_py import cad, sketch, feature_tree`.
///
/// Each submodule is registered via its own `register` function so
/// the per-module wiring stays close to the classes being added.
#[pymodule]
fn valenx_py(m: &Bound<'_, PyModule>) -> PyResult<()> {
    use pyo3::types::PyModule;

    let py = m.py();

    m.add("__version__", env!("CARGO_PKG_VERSION"))?;

    // ----- Submodule: valenx_py.cad -----
    let cad_mod = PyModule::new(py, "cad")?;
    cad::register(&cad_mod)?;
    m.add_submodule(&cad_mod)?;

    // ----- Submodule: valenx_py.sketch -----
    let sketch_mod = PyModule::new(py, "sketch")?;
    sketch::register(&sketch_mod)?;
    m.add_submodule(&sketch_mod)?;

    // ----- Submodule: valenx_py.feature_tree -----
    let ft_mod = PyModule::new(py, "feature_tree")?;
    feature_tree::register(&ft_mod)?;
    m.add_submodule(&ft_mod)?;

    // ----- Submodule: valenx_py.mesh (PyMesh return type) -----
    let mesh_mod = PyModule::new(py, "mesh")?;
    mesh::register(&mesh_mod)?;
    m.add_submodule(&mesh_mod)?;

    Ok(())
}
