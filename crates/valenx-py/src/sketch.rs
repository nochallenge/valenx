//! Python bindings for [`valenx_sketch`].
//!
//! [`PySketch`] wraps [`valenx_sketch::Sketch`] with the four most
//! useful authoring methods: `add_point`, `add_line`, `add_circle`,
//! and `add_constraint`. `extrude(depth)` is the bridge back into
//! [`valenx_cad::Solid`] for downstream CAD ops.
//!
//! Constraint encoding
//! -------------------
//!
//! [`valenx_sketch::constraint::Constraint`] is a Rust enum that
//! Python can't construct directly without a verbose `Constraint::new_*`
//! enumeration on the Rust side. Instead, `PySketch.add_constraint`
//! takes a `kind` string plus a dict of per-kind parameters. Each
//! supported kind has its expected keys:
//!
//! - `"coincident"`        — `{a, b}` (entity IDs)
//! - `"horizontal"`        — `{line}` (entity ID)
//! - `"vertical"`          — `{line}` (entity ID)
//! - `"parallel"`          — `{a, b}` (line IDs)
//! - `"perpendicular"`     — `{a, b}` (line IDs)
//! - `"tangent"`           — `{a, b}` (line-or-circle ID + circle ID)
//! - `"equal_length"`      — `{a, b}` (line IDs)
//! - `"distance"`          — `{a, b, target}` (point IDs + float target)
//! - `"angle"`             — `{a, b, target}` (line IDs + radians target)
//! - `"radius"`            — `{entity, target}` (circle/arc ID + float radius)
//!
//! Unknown kinds and missing keys both raise [`pyo3::exceptions::PyValueError`].

use pyo3::prelude::*;
use pyo3::types::{PyDict, PyModule};

use valenx_sketch::constraint::Constraint;
use valenx_sketch::geom::EntityId;

use crate::cad::PySolid;
use crate::error::{cad_err, sketch_err, value_err};

/// Wrapped [`valenx_sketch::Sketch`]. Held by value, not Arc — every
/// `add_*` mutates in place, so concurrent Python references to the
/// same sketch would race. PyO3's default `#[pyclass]` puts the inner
/// behind a `RefCell`, which is exactly what we want here.
#[pyclass(name = "Sketch", module = "valenx.sketch")]
#[derive(Default)]
pub struct PySketch {
    pub(crate) inner: valenx_sketch::Sketch,
}

impl PySketch {
    /// Borrow the inner sketch for use by feature_tree.add_sketch.
    pub(crate) fn clone_inner(&self) -> valenx_sketch::Sketch {
        self.inner.clone()
    }
}

#[pymethods]
impl PySketch {
    /// Empty sketch.
    #[new]
    fn new() -> Self {
        Self::default()
    }

    /// Add a point at `(x, y)`. Returns the entity ID (1-based).
    fn add_point(&mut self, x: f64, y: f64) -> usize {
        self.inner.add_point(x, y).0
    }

    /// Add a line between two existing point IDs. Returns the new
    /// entity ID.
    fn add_line(&mut self, start: usize, end: usize) -> PyResult<usize> {
        self.inner
            .add_line(EntityId(start), EntityId(end))
            .map(|id| id.0)
            .map_err(sketch_err)
    }

    /// Add a circle around an existing centre point with the given
    /// initial radius.
    fn add_circle(&mut self, center: usize, radius: f64) -> PyResult<usize> {
        self.inner
            .add_circle(EntityId(center), radius)
            .map(|id| id.0)
            .map_err(sketch_err)
    }

    /// Add an arc around an existing centre point. Angles in radians.
    fn add_arc(
        &mut self,
        center: usize,
        radius: f64,
        start_angle: f64,
        end_angle: f64,
    ) -> PyResult<usize> {
        self.inner
            .add_arc(EntityId(center), radius, start_angle, end_angle)
            .map(|id| id.0)
            .map_err(sketch_err)
    }

    /// Add a geometric constraint. See module docs for the supported
    /// `kind` strings and their parameter dicts.
    fn add_constraint(&mut self, kind: &str, params: &Bound<'_, PyDict>) -> PyResult<()> {
        let c = parse_constraint(kind, params)?;
        self.inner.add_constraint(c);
        Ok(())
    }

    /// Total number of constraints added so far.
    #[getter]
    fn constraint_count(&self) -> usize {
        self.inner.constraints.len()
    }

    /// Total number of entities added so far.
    #[getter]
    fn entity_count(&self) -> usize {
        self.inner.entities.len()
    }

    /// Extrude this sketch's closed profile by `depth` along +Z.
    /// Returns a Solid handle usable in cad.* ops.
    fn extrude(&self, depth: f64) -> PyResult<PySolid> {
        // Map the sketch error first (it covers profile-extraction
        // failures), then the cad error from the actual truck call.
        // `extrude` returns SketchError, which already wraps cad
        // errors in `SketchError::CadFailed`, so a single mapping
        // suffices.
        self.inner
            .extrude(depth)
            .map(PySolid::wrap)
            .map_err(sketch_err)
    }
}

/// Parse a `(kind, params)` pair into a typed
/// [`valenx_sketch::constraint::Constraint`].
fn parse_constraint(kind: &str, p: &Bound<'_, PyDict>) -> PyResult<Constraint> {
    // Helper closures — fetch a required key, surface a clear error
    // when missing.
    let id = |key: &str| -> PyResult<EntityId> {
        match p.get_item(key)? {
            Some(v) => v
                .extract::<usize>()
                .map(EntityId)
                .map_err(|e| value_err(format!("constraint key '{key}': {e}"))),
            None => Err(value_err(format!(
                "constraint '{kind}' missing required key '{key}'"
            ))),
        }
    };
    let scalar = |key: &str| -> PyResult<f64> {
        match p.get_item(key)? {
            Some(v) => v
                .extract::<f64>()
                .map_err(|e| value_err(format!("constraint key '{key}': {e}"))),
            None => Err(value_err(format!(
                "constraint '{kind}' missing required key '{key}'"
            ))),
        }
    };

    Ok(match kind {
        "coincident" => Constraint::Coincident {
            a: id("a")?,
            b: id("b")?,
        },
        "horizontal" => Constraint::Horizontal(id("line")?),
        "vertical" => Constraint::Vertical(id("line")?),
        "parallel" => Constraint::Parallel {
            a: id("a")?,
            b: id("b")?,
        },
        "perpendicular" => Constraint::Perpendicular {
            a: id("a")?,
            b: id("b")?,
        },
        "tangent" => Constraint::Tangent {
            line_or_circle_a: id("a")?,
            circle_b: id("b")?,
        },
        "equal_length" => Constraint::EqualLength {
            a: id("a")?,
            b: id("b")?,
        },
        "distance" => Constraint::Distance {
            a: id("a")?,
            b: id("b")?,
            target: scalar("target")?,
        },
        "angle" => Constraint::Angle {
            a: id("a")?,
            b: id("b")?,
            target: scalar("target")?,
        },
        "radius" => Constraint::Radius {
            circle_or_arc: id("entity")?,
            target: scalar("target")?,
        },
        other => {
            return Err(value_err(format!(
                "unknown constraint kind '{other}' — supported: coincident, \
                 horizontal, vertical, parallel, perpendicular, tangent, \
                 equal_length, distance, angle, radius"
            )))
        }
    })
}

// Silence dead-code on the cad_err import — kept here for symmetry
// with the other modules and so a future constraint-parser failure
// path can drop a Cad error directly without re-importing.
#[allow(dead_code)]
fn _suppress_unused_cad_err() -> fn(valenx_cad::CadError) -> PyErr {
    cad_err
}

/// Register every sketch.* item on the submodule.
pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PySketch>()?;
    Ok(())
}
