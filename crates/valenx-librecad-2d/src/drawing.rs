//! Drawing data model — Drawing2D + Layer + Block + Entity2D enum.
//!
//! Maps 1:1 onto DXF concepts: a drawing has many layers, blocks
//! (reusable groups), and entities; an entity belongs to a layer; an
//! INSERT entity references a block.

use serde::{Deserialize, Serialize};

/// Top-level drawing.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct Drawing2D {
    /// All layers (always non-empty after [`Drawing2D::new`] — the
    /// `"0"` default layer is created automatically).
    pub layers: Vec<Layer>,
    /// All blocks.
    pub blocks: Vec<Block>,
    /// All entities at the drawing root (not inside any block).
    pub entities: Vec<Entity2D>,
}

impl Drawing2D {
    /// Fresh drawing with a default layer named `"0"`.
    pub fn new() -> Self {
        Self {
            layers: vec![Layer::default_layer()],
            blocks: Vec::new(),
            entities: Vec::new(),
        }
    }

    /// Add an entity (panics if `layer_name` is unknown — callers
    /// should add layers first).
    pub fn add(&mut self, entity: Entity2D) {
        self.entities.push(entity);
    }
}

/// A DXF layer — name + colour index + linetype + visibility.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Layer {
    /// Layer name.
    pub name: String,
    /// AutoCAD colour index (1..=255, 7=white default).
    pub color: u8,
    /// Linetype string (`"CONTINUOUS"`, `"DASHED"`, …).
    pub linetype: String,
    /// Visibility flag — `false` hides the layer.
    pub visible: bool,
}

impl Layer {
    /// The default `"0"` layer.
    pub fn default_layer() -> Self {
        Self {
            name: "0".into(),
            color: 7,
            linetype: "CONTINUOUS".into(),
            visible: true,
        }
    }
}

/// A reusable group of entities, instantiated via [`Entity2D::Insert`].
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Block {
    /// Block name.
    pub name: String,
    /// Insertion-point reference (block-local origin).
    pub origin: [f64; 2],
    /// Entities making up the block.
    pub entities: Vec<Entity2D>,
}

/// 2D entity zoo.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub enum Entity2D {
    /// Straight line segment.
    Line {
        /// Owning layer name.
        layer: String,
        /// Start point.
        a: [f64; 2],
        /// End point.
        b: [f64; 2],
    },
    /// Full circle.
    Circle {
        /// Owning layer name.
        layer: String,
        /// Centre.
        centre: [f64; 2],
        /// Radius.
        radius: f64,
    },
    /// Circular arc — angle units are degrees, CCW from +X.
    Arc {
        /// Owning layer name.
        layer: String,
        /// Centre.
        centre: [f64; 2],
        /// Radius.
        radius: f64,
        /// Start angle (degrees).
        start_angle_deg: f64,
        /// End angle (degrees).
        end_angle_deg: f64,
    },
    /// Light-weight polyline (LWPOLYLINE in DXF).
    Polyline {
        /// Owning layer name.
        layer: String,
        /// Vertices.
        vertices: Vec<[f64; 2]>,
        /// Closed flag.
        closed: bool,
    },
    /// B-spline (control-point representation; v1 stores only control
    /// points and a degree).
    Spline {
        /// Owning layer name.
        layer: String,
        /// Control points.
        control_points: Vec<[f64; 2]>,
        /// Degree (defaults to 3 for cubic).
        degree: u8,
    },
    /// Solid hatch — bounded by a polyline.
    Hatch {
        /// Owning layer name.
        layer: String,
        /// Outer boundary loop.
        boundary: Vec<[f64; 2]>,
        /// Hatch pattern name (`"SOLID"` default).
        pattern: String,
    },
    /// Single-line text.
    Text {
        /// Owning layer name.
        layer: String,
        /// Insertion point.
        position: [f64; 2],
        /// Text height (drawing units).
        height: f64,
        /// Text content.
        text: String,
    },
    /// Multi-line text (MTEXT in DXF) — wraps to `width`.
    MText {
        /// Owning layer name.
        layer: String,
        /// Insertion point.
        position: [f64; 2],
        /// Text height (drawing units).
        height: f64,
        /// Wrap width.
        width: f64,
        /// Text content (may contain newlines).
        text: String,
    },
    /// Linear dimension.
    Dimension {
        /// Owning layer name.
        layer: String,
        /// Definition point 1.
        a: [f64; 2],
        /// Definition point 2.
        b: [f64; 2],
        /// Text location.
        text_pos: [f64; 2],
        /// Measured-value label.
        text: String,
    },
    /// Block reference (INSERT entity).
    Insert {
        /// Owning layer name.
        layer: String,
        /// Referenced block name.
        block: String,
        /// Insertion point.
        position: [f64; 2],
        /// Uniform scale.
        scale: f64,
        /// Rotation in degrees CCW.
        rotation_deg: f64,
    },
}

impl Entity2D {
    /// Stable kebab-cased label, used by [`crate::dxf`] dispatcher
    /// and diagnostics.
    pub fn kind(&self) -> &'static str {
        match self {
            Entity2D::Line { .. } => "LINE",
            Entity2D::Circle { .. } => "CIRCLE",
            Entity2D::Arc { .. } => "ARC",
            Entity2D::Polyline { .. } => "LWPOLYLINE",
            Entity2D::Spline { .. } => "SPLINE",
            Entity2D::Hatch { .. } => "HATCH",
            Entity2D::Text { .. } => "TEXT",
            Entity2D::MText { .. } => "MTEXT",
            Entity2D::Dimension { .. } => "DIMENSION",
            Entity2D::Insert { .. } => "INSERT",
        }
    }

    /// Owning layer.
    pub fn layer(&self) -> &str {
        match self {
            Entity2D::Line { layer, .. }
            | Entity2D::Circle { layer, .. }
            | Entity2D::Arc { layer, .. }
            | Entity2D::Polyline { layer, .. }
            | Entity2D::Spline { layer, .. }
            | Entity2D::Hatch { layer, .. }
            | Entity2D::Text { layer, .. }
            | Entity2D::MText { layer, .. }
            | Entity2D::Dimension { layer, .. }
            | Entity2D::Insert { layer, .. } => layer,
        }
    }
}
