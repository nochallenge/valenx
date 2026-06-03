//! Drawing tree — HeeksCAD's object tree with Layer grouping.

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

/// A 2D sketch entity — line or arc within a sketch.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum SketchEntity {
    /// Line.
    Line {
        /// Start point.
        a: [f64; 2],
        /// End point.
        b: [f64; 2],
    },
    /// Arc.
    Arc {
        /// Centre.
        centre: [f64; 2],
        /// Radius.
        radius: f64,
        /// Start angle in radians.
        start: f64,
        /// Sweep in radians.
        sweep: f64,
    },
}

/// HeeksCAD primitive object.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum HeeksObject {
    /// 2D sketch.
    Sketch {
        /// Object name.
        name: String,
        /// Sketch plane normal direction (e.g. (0,0,1)).
        plane_normal: Vector3<f64>,
        /// Entities making up the sketch.
        entities: Vec<SketchEntity>,
    },
    /// Pad (extrude a sketch upward).
    Pad {
        /// Object name.
        name: String,
        /// Referenced sketch by name.
        sketch_ref: String,
        /// Extrusion height.
        height: f64,
    },
    /// Pocket (subtract a sketch's interior to a depth).
    Pocket {
        /// Object name.
        name: String,
        /// Referenced sketch.
        sketch_ref: String,
        /// Depth.
        depth: f64,
    },
    /// Drill — a series of holes at named positions.
    Drill {
        /// Object name.
        name: String,
        /// Hole positions in XY.
        positions: Vec<[f64; 2]>,
        /// Depth.
        depth: f64,
        /// Diameter.
        diameter: f64,
    },
}

impl HeeksObject {
    /// Object name.
    pub fn name(&self) -> &str {
        match self {
            HeeksObject::Sketch { name, .. }
            | HeeksObject::Pad { name, .. }
            | HeeksObject::Pocket { name, .. }
            | HeeksObject::Drill { name, .. } => name,
        }
    }
}

/// Layer container.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Layer {
    /// Layer name.
    pub name: String,
    /// Visibility flag.
    pub visible: bool,
    /// Layer-owned objects.
    pub objects: Vec<HeeksObject>,
}

impl Layer {
    /// New visible empty layer.
    pub fn new(name: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            visible: true,
            objects: Vec::new(),
        }
    }
}

/// Top-level drawing.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct Drawing {
    /// Layers.
    pub layers: Vec<Layer>,
}

impl Drawing {
    /// New drawing with a default layer named `"0"`.
    pub fn new() -> Self {
        Self {
            layers: vec![Layer::new("0")],
        }
    }

    /// Mutable access to a layer by name.
    pub fn layer_mut(&mut self, name: &str) -> Option<&mut Layer> {
        self.layers.iter_mut().find(|l| l.name == name)
    }

    /// Read-only access to a layer by name.
    pub fn layer(&self, name: &str) -> Option<&Layer> {
        self.layers.iter().find(|l| l.name == name)
    }

    /// Total object count across all layers.
    pub fn total_objects(&self) -> usize {
        self.layers.iter().map(|l| l.objects.len()).sum()
    }
}
