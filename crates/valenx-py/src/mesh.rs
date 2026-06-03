//! Python bindings for [`valenx_mesh::Mesh`].
//!
//! Read-only wrapper — Python callers receive [`PyMesh`] from
//! `cad.solid_to_mesh` and use it for inspection (triangle count,
//! vertex count, bounds) or for round-tripping into a downstream
//! Rust API. There's no constructor; meshes come from CAD ops.
//!
//! Why no `to_stl` here
//! --------------------
//!
//! STL / OBJ export lives in `valenx-export` and `valenx-mesh::format`.
//! We deliberately keep the Python `valenx.mesh` surface tiny in v1 —
//! Python users who need to write meshes can build the toolpath / RON
//! envelope on the Rust side and let the desktop UI handle export,
//! or pipe the mesh out through numpy via a future buffer-protocol
//! wrapper.

use pyo3::prelude::*;
use pyo3::types::PyModule;

/// Opaque wrapper around a [`valenx_mesh::Mesh`].
#[pyclass(name = "Mesh", module = "valenx.mesh", frozen)]
#[derive(Clone)]
pub struct PyMesh {
    inner: std::sync::Arc<valenx_mesh::Mesh>,
}

impl PyMesh {
    /// Wrap a Rust [`valenx_mesh::Mesh`] for Python consumption.
    pub fn wrap(mesh: valenx_mesh::Mesh) -> Self {
        Self {
            inner: std::sync::Arc::new(mesh),
        }
    }
    /// Borrow the inner mesh — for future inter-module bridges.
    #[allow(dead_code)]
    pub fn inner(&self) -> &valenx_mesh::Mesh {
        &self.inner
    }
}

#[pymethods]
impl PyMesh {
    /// Number of distinct vertex positions.
    #[getter]
    fn vertex_count(&self) -> usize {
        self.inner.nodes.len()
    }
    /// Total triangle count across every element block. Quad / n-gon
    /// blocks are NOT counted — the BRep tessellator only emits Tri3.
    #[getter]
    fn triangle_count(&self) -> usize {
        self.inner.total_elements()
    }
    /// Mesh identifier — usually "cad" for the BRep tessellator's
    /// output. Mirrors `valenx_mesh::Mesh::id`.
    #[getter]
    fn id(&self) -> String {
        self.inner.id.clone()
    }
    fn __repr__(&self) -> String {
        format!(
            "Mesh(id='{}', vertices={}, triangles={})",
            self.inner.id,
            self.inner.nodes.len(),
            self.inner.total_elements()
        )
    }
}

/// Register every mesh.* item on the submodule.
pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyMesh>()?;
    Ok(())
}
