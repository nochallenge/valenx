//! Python bindings for [`valenx_feature_tree`].
//!
//! [`PyFeatureTree`] wraps [`valenx_feature_tree::FeatureTree`] with
//! `add_sketch` (consumes a [`crate::sketch::PySketch`] handle),
//! `add_feature` (kind string + params dict), and `replay` (returns the
//! final solid).
//!
//! Feature encoding
//! ----------------
//!
//! Each [`valenx_feature_tree::Feature`] variant has its own params
//! struct. To keep the Python surface flat (no per-variant pyclass
//! explosion), `add_feature` takes a `kind` string and a parameter
//! dict whose keys depend on the kind:
//!
//! - `"pad"`     — `{sketch: int, depth: float, direction_positive: bool}`
//! - `"pocket"`  — `{sketch: int, depth: float, direction_positive: bool}`
//! - `"revolve"` — `{sketch: int, axis_origin: (x,y,z), axis_direction: (x,y,z), angle: float}`
//! - `"mirror"`  — `{target: int, plane_origin: (x,y,z), plane_normal: (x,y,z), keep_original: bool}`
//! - `"linear_pattern"`   — `{target: int, direction: (x,y,z), count: int, spacing: float}`
//! - `"circular_pattern"` — `{target: int, axis_origin: (x,y,z), axis_direction: (x,y,z), count: int, total_angle: float}`
//! - `"fillet"`           — `{target: int, radius: float, threshold_deg: float}`
//! - `"chamfer"`          — `{target: int, distance: float, threshold_deg: float}`
//! - `"imported_solid"`   — `{source_path: str}`
//!
//! The `name` argument is the display label shown in the tree UI.
//!
//! Unknown kinds and missing keys raise [`pyo3::exceptions::PyValueError`].

use nalgebra::Vector3;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyModule};

use valenx_feature_tree::feature::{
    ChamferParams, CircularPatternParams, FeatureId, FilletParams, ImportedSolidParams,
    LinearPatternParams, MirrorParams, PadParams, PocketParams, RevolveParams, SketchRef,
};
use valenx_feature_tree::{Feature, FeatureTree};

use crate::cad::PySolid;
use crate::error::{feature_err, value_err};
use crate::sketch::PySketch;

/// Wrapped [`valenx_feature_tree::FeatureTree`]. Held in a RefCell
/// implicitly via `#[pyclass]` since `add_*` methods take `&mut self`.
#[pyclass(name = "FeatureTree", module = "valenx.feature_tree")]
#[derive(Default)]
pub struct PyFeatureTree {
    inner: FeatureTree,
}

#[pymethods]
impl PyFeatureTree {
    /// Empty tree.
    #[new]
    fn new() -> Self {
        Self::default()
    }

    /// Append a sketch into the tree's sketch table. Returns the
    /// 0-based `SketchRef` index that features can reference.
    fn add_sketch(&mut self, sketch: &PySketch) -> usize {
        // Take a snapshot of the sketch so the tree owns it. Python
        // can keep using the same PySketch afterwards without
        // affecting the tree's copy.
        let s = sketch.clone_inner();
        self.inner.add_sketch(s).0
    }

    /// Append a feature. `kind` is one of the strings in the module
    /// docs; `params` carries the per-kind keys; `name` is the
    /// display label. Returns the 0-based `FeatureId`.
    fn add_feature(
        &mut self,
        kind: &str,
        params: &Bound<'_, PyDict>,
        name: &str,
    ) -> PyResult<usize> {
        let feat = parse_feature(kind, params)?;
        Ok(self.inner.add_feature(feat, name).0)
    }

    /// Total number of features.
    #[getter]
    fn feature_count(&self) -> usize {
        self.inner.features.len()
    }

    /// Total number of sketches.
    #[getter]
    fn sketch_count(&self) -> usize {
        self.inner.sketches.len()
    }

    /// Walk the tree and return the final Solid. Returns None if the
    /// tree is empty or every feature is suppressed.
    fn replay(&self) -> PyResult<Option<PySolid>> {
        valenx_feature_tree::replay(&self.inner)
            .map(|opt| opt.map(PySolid::wrap))
            .map_err(feature_err)
    }
}

/// Parse a `(kind, params)` pair into a [`Feature`].
fn parse_feature(kind: &str, p: &Bound<'_, PyDict>) -> PyResult<Feature> {
    let int = |key: &str| -> PyResult<usize> {
        match p.get_item(key)? {
            Some(v) => v
                .extract::<usize>()
                .map_err(|e| value_err(format!("feature '{kind}' key '{key}': {e}"))),
            None => Err(value_err(format!("feature '{kind}' missing key '{key}'"))),
        }
    };
    let scalar = |key: &str| -> PyResult<f64> {
        match p.get_item(key)? {
            Some(v) => v
                .extract::<f64>()
                .map_err(|e| value_err(format!("feature '{kind}' key '{key}': {e}"))),
            None => Err(value_err(format!("feature '{kind}' missing key '{key}'"))),
        }
    };
    let boolean = |key: &str| -> PyResult<bool> {
        match p.get_item(key)? {
            Some(v) => v
                .extract::<bool>()
                .map_err(|e| value_err(format!("feature '{kind}' key '{key}': {e}"))),
            None => Err(value_err(format!("feature '{kind}' missing key '{key}'"))),
        }
    };
    let vec3 = |key: &str| -> PyResult<Vector3<f64>> {
        match p.get_item(key)? {
            Some(v) => {
                let t: (f64, f64, f64) = v
                    .extract()
                    .map_err(|e| value_err(format!("feature '{kind}' key '{key}': {e}")))?;
                Ok(Vector3::new(t.0, t.1, t.2))
            }
            None => Err(value_err(format!("feature '{kind}' missing key '{key}'"))),
        }
    };
    let count_u32 = |key: &str| -> PyResult<u32> {
        match p.get_item(key)? {
            Some(v) => v
                .extract::<u32>()
                .map_err(|e| value_err(format!("feature '{kind}' key '{key}': {e}"))),
            None => Err(value_err(format!("feature '{kind}' missing key '{key}'"))),
        }
    };
    let string = |key: &str| -> PyResult<String> {
        match p.get_item(key)? {
            Some(v) => v
                .extract::<String>()
                .map_err(|e| value_err(format!("feature '{kind}' key '{key}': {e}"))),
            None => Err(value_err(format!("feature '{kind}' missing key '{key}'"))),
        }
    };

    Ok(match kind {
        "pad" => Feature::Pad(PadParams {
            sketch: SketchRef(int("sketch")?),
            depth: scalar("depth")?.into(),
            direction_positive: boolean("direction_positive")?,
        }),
        "pocket" => Feature::Pocket(PocketParams {
            sketch: SketchRef(int("sketch")?),
            depth: scalar("depth")?.into(),
            direction_positive: boolean("direction_positive")?,
        }),
        "revolve" => Feature::Revolve(RevolveParams {
            sketch: SketchRef(int("sketch")?),
            axis_origin: vec3("axis_origin")?,
            axis_direction: vec3("axis_direction")?,
            angle: scalar("angle")?.into(),
        }),
        "mirror" => Feature::Mirror(MirrorParams {
            target: FeatureId(int("target")?),
            plane_origin: vec3("plane_origin")?,
            plane_normal: vec3("plane_normal")?,
            keep_original: boolean("keep_original")?,
        }),
        "linear_pattern" => Feature::LinearPattern(LinearPatternParams {
            target: FeatureId(int("target")?),
            direction: vec3("direction")?,
            count: count_u32("count")?,
            spacing: scalar("spacing")?,
        }),
        "circular_pattern" => Feature::CircularPattern(CircularPatternParams {
            target: FeatureId(int("target")?),
            axis_origin: vec3("axis_origin")?,
            axis_direction: vec3("axis_direction")?,
            count: count_u32("count")?,
            total_angle: scalar("total_angle")?,
        }),
        "fillet" => Feature::Fillet(FilletParams {
            target: FeatureId(int("target")?),
            radius: scalar("radius")?,
            threshold_deg: scalar("threshold_deg")?,
            // Phase 14: Python callers omit `edge_indices` for the
            // backward-compatible "auto by angle threshold" behavior.
            // Explicit indices are reserved for a future Python API
            // expansion when the BRep path graduates from soft-error
            // fall-through to live BRep substitution.
            edge_indices: None,
        }),
        "chamfer" => Feature::Chamfer(ChamferParams {
            target: FeatureId(int("target")?),
            distance: scalar("distance")?,
            threshold_deg: scalar("threshold_deg")?,
            edge_indices: None,
        }),
        "imported_solid" => Feature::ImportedSolid(ImportedSolidParams {
            source_path: string("source_path")?,
        }),
        other => {
            return Err(value_err(format!(
                "unknown feature kind '{other}' — supported: pad, pocket, \
                 revolve, mirror, linear_pattern, circular_pattern, fillet, \
                 chamfer, imported_solid"
            )))
        }
    })
}

/// Register every feature_tree.* item on the submodule.
pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyFeatureTree>()?;
    Ok(())
}
