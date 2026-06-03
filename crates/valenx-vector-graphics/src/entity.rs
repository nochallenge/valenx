//! SVG-style 2D entity types.

use nalgebra::Vector2;
use serde::{Deserialize, Serialize};

/// One path command — mirrors the SVG `d` attribute mini-language.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum PathSegment {
    /// `M x y` — move to.
    MoveTo(Vector2<f64>),
    /// `L x y` — line to.
    LineTo(Vector2<f64>),
    /// `C cx1 cy1 cx2 cy2 x y` — cubic Bezier.
    CurveTo {
        /// First control point.
        c1: Vector2<f64>,
        /// Second control point.
        c2: Vector2<f64>,
        /// Endpoint.
        end: Vector2<f64>,
    },
    /// `Q cx cy x y` — quadratic Bezier.
    QuadTo {
        /// Control point.
        c: Vector2<f64>,
        /// Endpoint.
        end: Vector2<f64>,
    },
    /// `A rx ry rot large_arc_flag sweep_flag x y` — elliptical arc.
    Arc {
        /// X radius.
        rx: f64,
        /// Y radius.
        ry: f64,
        /// X-axis rotation in degrees.
        x_axis_rotation_deg: f64,
        /// Large-arc flag.
        large_arc: bool,
        /// Sweep flag.
        sweep: bool,
        /// Endpoint.
        end: Vector2<f64>,
    },
    /// `Z` — close current sub-path.
    Close,
}

/// A vector entity.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum VectorEntity {
    /// Straight line segment.
    Line {
        /// Start.
        a: Vector2<f64>,
        /// End.
        b: Vector2<f64>,
    },
    /// SVG-style path.
    Path(Vec<PathSegment>),
    /// Axis-aligned rectangle.
    Rect {
        /// Origin (lowest corner).
        origin: Vector2<f64>,
        /// Size (width, height).
        size: Vector2<f64>,
    },
    /// Ellipse (centre + radii).
    Ellipse {
        /// Centre.
        centre: Vector2<f64>,
        /// X radius.
        rx: f64,
        /// Y radius.
        ry: f64,
    },
    /// Closed polygon (vertex loop).
    Polygon(Vec<Vector2<f64>>),
    /// Text label.
    Text {
        /// Anchor point.
        anchor: Vector2<f64>,
        /// Font size in user units.
        font_size: f64,
        /// String content.
        text: String,
    },
}

impl VectorEntity {
    /// Tag string.
    pub fn kind(&self) -> &'static str {
        match self {
            VectorEntity::Line { .. } => "line",
            VectorEntity::Path(_) => "path",
            VectorEntity::Rect { .. } => "rect",
            VectorEntity::Ellipse { .. } => "ellipse",
            VectorEntity::Polygon(_) => "polygon",
            VectorEntity::Text { .. } => "text",
        }
    }
}
