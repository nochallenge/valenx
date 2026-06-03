//! Python bindings for [`valenx_cad`] primitives, booleans, and
//! tessellation.
//!
//! Each function is a thin `#[pyfunction]` wrapper that calls the
//! upstream Rust API and maps its typed error into a Python
//! exception via [`crate::error::cad_err`].
//!
//! Names match `valenx_cad` exactly except for `box_solid` → `box`
//! (Python users expect the unqualified name; `box` isn't a Python
//! reserved word so it's safe). The `_solid` suffix on the Rust side
//! exists only to avoid colliding with the language keyword in Rust.

use pyo3::prelude::*;
use pyo3::types::PyModule;

use crate::error::cad_err;
use crate::mesh::PyMesh;

/// Opaque wrapper around a [`valenx_cad::Solid`] held in a Python
/// object. The Solid stays on the Rust heap; Python sees only the
/// wrapper handle, which it can pass back into other cad.* calls.
///
/// Solids are immutable from Python — every op returns a fresh one.
/// Internally we hold an `Arc` so the same `Solid` can be shared
/// between multiple [`PySolid`] handles without cloning the BRep.
#[pyclass(name = "Solid", module = "valenx.cad", frozen)]
#[derive(Clone)]
pub struct PySolid {
    inner: std::sync::Arc<valenx_cad::Solid>,
}

impl PySolid {
    /// Wrap a Rust [`valenx_cad::Solid`] for Python consumption.
    pub fn wrap(solid: valenx_cad::Solid) -> Self {
        Self {
            inner: std::sync::Arc::new(solid),
        }
    }
    /// Borrow the underlying solid — used by inter-module bridges
    /// (sketch.extrude, feature_tree.replay).
    pub fn inner(&self) -> &valenx_cad::Solid {
        &self.inner
    }
}

#[pymethods]
impl PySolid {
    /// Number of BRep faces. Returns 0 for mesh-backed solids
    /// (see [`valenx_cad::Solid`] docs).
    fn faces(&self) -> usize {
        self.inner.faces()
    }
    /// Number of distinct BRep edges. Returns 0 for mesh-backed solids.
    fn edges(&self) -> usize {
        self.inner.edges()
    }
    /// Number of distinct BRep vertices. Returns 0 for mesh-backed solids.
    fn vertices(&self) -> usize {
        self.inner.vertices()
    }
    /// Translate the solid by `(dx, dy, dz)`. Returns a fresh handle.
    /// Round-6: returns `PyValueError` when any of `dx`, `dy`, `dz` is
    /// non-finite (NaN or ±inf) — the previous infallible signature
    /// silently propagated bad inputs into the BRep matrix.
    fn translated(&self, dx: f64, dy: f64, dz: f64) -> PyResult<Self> {
        self.inner
            .translated(dx, dy, dz)
            .map(Self::wrap)
            .map_err(cad_err)
    }
    fn __repr__(&self) -> String {
        format!(
            "Solid(faces={}, edges={}, vertices={})",
            self.inner.faces(),
            self.inner.edges(),
            self.inner.vertices()
        )
    }
}

// ----- Primitives -----

/// Axis-aligned box from origin to `(dx, dy, dz)`.
#[pyfunction]
#[pyo3(name = "box")]
fn py_box(dx: f64, dy: f64, dz: f64) -> PyResult<PySolid> {
    valenx_cad::box_solid(dx, dy, dz)
        .map(PySolid::wrap)
        .map_err(cad_err)
}

/// Right circular cylinder, base disk in XY plane, axis along +Z.
#[pyfunction]
fn cylinder(radius: f64, height: f64) -> PyResult<PySolid> {
    valenx_cad::cylinder(radius, height)
        .map(PySolid::wrap)
        .map_err(cad_err)
}

/// Sphere centred on the origin.
#[pyfunction]
fn sphere(radius: f64) -> PyResult<PySolid> {
    valenx_cad::sphere(radius)
        .map(PySolid::wrap)
        .map_err(cad_err)
}

/// Truncated cone (frustum). `top_radius=0` gives a pointed cone.
#[pyfunction]
fn cone(base_radius: f64, top_radius: f64, height: f64) -> PyResult<PySolid> {
    valenx_cad::cone(base_radius, top_radius, height)
        .map(PySolid::wrap)
        .map_err(cad_err)
}

/// Torus with major axis along Z. `minor_radius < major_radius`.
#[pyfunction]
fn torus(major_radius: f64, minor_radius: f64) -> PyResult<PySolid> {
    valenx_cad::torus(major_radius, minor_radius)
        .map(PySolid::wrap)
        .map_err(cad_err)
}

/// Extrude a closed polygon profile (list of `(x, y)` tuples) by
/// `height` along +Z. Profile must have at least 3 points; the
/// closing edge is implicit.
#[pyfunction]
fn prism(profile_xy: Vec<(f64, f64)>, height: f64) -> PyResult<PySolid> {
    valenx_cad::prism(&profile_xy, height)
        .map(PySolid::wrap)
        .map_err(cad_err)
}

// ----- Booleans -----

/// Union (A ∪ B).
#[pyfunction]
fn union(a: &PySolid, b: &PySolid) -> PyResult<PySolid> {
    valenx_cad::union(a.inner(), b.inner())
        .map(PySolid::wrap)
        .map_err(cad_err)
}

/// Intersection (A ∩ B).
#[pyfunction]
fn intersection(a: &PySolid, b: &PySolid) -> PyResult<PySolid> {
    valenx_cad::intersection(a.inner(), b.inner())
        .map(PySolid::wrap)
        .map_err(cad_err)
}

/// Difference (A − B).
#[pyfunction]
fn difference(a: &PySolid, b: &PySolid) -> PyResult<PySolid> {
    valenx_cad::difference(a.inner(), b.inner())
        .map(PySolid::wrap)
        .map_err(cad_err)
}

// ----- Tessellation -----

/// Tessellate a Solid into a triangle mesh with the given chord-error
/// tolerance. Smaller tolerance = denser mesh.
#[pyfunction]
fn solid_to_mesh(solid: &PySolid, tolerance: f64) -> PyResult<PyMesh> {
    valenx_cad::solid_to_mesh(solid.inner(), tolerance)
        .map(PyMesh::wrap)
        .map_err(cad_err)
}

/// Register every cad.* function on the submodule.
pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PySolid>()?;
    m.add_function(wrap_pyfunction!(py_box, m)?)?;
    m.add_function(wrap_pyfunction!(cylinder, m)?)?;
    m.add_function(wrap_pyfunction!(sphere, m)?)?;
    m.add_function(wrap_pyfunction!(cone, m)?)?;
    m.add_function(wrap_pyfunction!(torus, m)?)?;
    m.add_function(wrap_pyfunction!(prism, m)?)?;
    m.add_function(wrap_pyfunction!(union, m)?)?;
    m.add_function(wrap_pyfunction!(intersection, m)?)?;
    m.add_function(wrap_pyfunction!(difference, m)?)?;
    m.add_function(wrap_pyfunction!(solid_to_mesh, m)?)?;
    Ok(())
}
