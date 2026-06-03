//! Draft entity variants.
//!
//! Each variant is a 2D primitive expressed in the local frame of the
//! parent [`super::WorkingPlane`]. There is no parametric solver —
//! coordinates are stored directly. Direct-modelling, suitable for
//! construction lines, annotations, and reference geometry.

use serde::{Deserialize, Serialize};

/// One drawable element of a [`super::DraftDocument`].
///
/// All numeric fields are 2D coordinates in the local frame of the
/// document's working plane. The UI / renderer projects them to
/// world space via [`super::WorkingPlane::local_to_world`].
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum DraftEntity {
    /// A straight line segment between two endpoints.
    Line {
        /// Start endpoint (x, y).
        start: [f64; 2],
        /// End endpoint (x, y).
        end: [f64; 2],
    },

    /// An open or closed sequence of straight segments through `points`.
    /// `closed = true` connects the last point back to the first.
    Polyline {
        /// Ordered vertices of the polyline.
        points: Vec<[f64; 2]>,
        /// When true, the renderer closes the loop by drawing a final
        /// segment from `points.last()` back to `points.first()`.
        closed: bool,
    },

    /// A circular arc with explicit start / end angles in radians.
    /// Angles increase counter-clockwise; sweep direction is
    /// `sign(end_angle - start_angle)`.
    Arc {
        /// Centre of the arc (x, y).
        center: [f64; 2],
        /// Radius (must be positive at render time; zero is treated
        /// as a no-op).
        radius: f64,
        /// Start angle in radians.
        start_angle: f64,
        /// End angle in radians.
        end_angle: f64,
    },

    /// A full circle centred at `center` with the given `radius`.
    Circle {
        /// Centre of the circle.
        center: [f64; 2],
        /// Radius.
        radius: f64,
    },

    /// An axis-aligned rectangle defined by its `min` (bottom-left) and
    /// `max` (top-right) corners in the local frame.
    Rectangle {
        /// Bottom-left corner (x_min, y_min).
        min: [f64; 2],
        /// Top-right corner (x_max, y_max).
        max: [f64; 2],
    },

    /// A regular polygon inscribed in a circle of `radius` centred at
    /// `center`, with `sides` vertices. The first vertex sits at
    /// angle 0 (i.e. +x of the local frame) and subsequent vertices
    /// are placed counter-clockwise.
    Polygon {
        /// Centre of the inscribing circle.
        center: [f64; 2],
        /// Inscribing radius.
        radius: f64,
        /// Vertex count (3 or more is sensible).
        sides: u32,
    },

    /// A linear dimension between two points, drawn parallel to the
    /// `from`→`to` axis at perpendicular distance `offset`.
    LinearDimension {
        /// First measured endpoint.
        from: [f64; 2],
        /// Second measured endpoint.
        to: [f64; 2],
        /// Perpendicular offset of the dimension line. Positive
        /// places the dimension on the +Y side of the from→to axis.
        offset: f64,
    },

    /// A text label placed at `position` with the given content and
    /// nominal `size` (in local-frame units).
    Text {
        /// Anchor position (text grows to the +x direction).
        position: [f64; 2],
        /// Text content.
        content: String,
        /// Nominal character size in local-frame units.
        size: f64,
    },
}

impl DraftEntity {
    /// Short human-readable kind tag (used by the UI entity list).
    pub fn kind(&self) -> &'static str {
        match self {
            DraftEntity::Line { .. } => "Line",
            DraftEntity::Polyline { .. } => "Polyline",
            DraftEntity::Arc { .. } => "Arc",
            DraftEntity::Circle { .. } => "Circle",
            DraftEntity::Rectangle { .. } => "Rectangle",
            DraftEntity::Polygon { .. } => "Polygon",
            DraftEntity::LinearDimension { .. } => "LinearDim",
            DraftEntity::Text { .. } => "Text",
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn constructs_each_variant_and_reports_kind() {
        let cases: Vec<(DraftEntity, &'static str)> = vec![
            (
                DraftEntity::Line {
                    start: [0.0, 0.0],
                    end: [1.0, 1.0],
                },
                "Line",
            ),
            (
                DraftEntity::Polyline {
                    points: vec![[0.0, 0.0], [1.0, 0.0]],
                    closed: false,
                },
                "Polyline",
            ),
            (
                DraftEntity::Arc {
                    center: [0.0, 0.0],
                    radius: 1.0,
                    start_angle: 0.0,
                    end_angle: std::f64::consts::PI,
                },
                "Arc",
            ),
            (
                DraftEntity::Circle {
                    center: [0.0, 0.0],
                    radius: 2.0,
                },
                "Circle",
            ),
            (
                DraftEntity::Rectangle {
                    min: [0.0, 0.0],
                    max: [1.0, 1.0],
                },
                "Rectangle",
            ),
            (
                DraftEntity::Polygon {
                    center: [0.0, 0.0],
                    radius: 1.0,
                    sides: 6,
                },
                "Polygon",
            ),
            (
                DraftEntity::LinearDimension {
                    from: [0.0, 0.0],
                    to: [10.0, 0.0],
                    offset: 1.0,
                },
                "LinearDim",
            ),
            (
                DraftEntity::Text {
                    position: [0.0, 0.0],
                    content: "hi".into(),
                    size: 1.0,
                },
                "Text",
            ),
        ];
        for (e, expected_kind) in cases {
            assert_eq!(e.kind(), expected_kind, "wrong kind for {e:?}");
        }
    }

    #[test]
    fn line_round_trips_via_ron() {
        let l = DraftEntity::Line {
            start: [1.0, 2.0],
            end: [3.0, 4.0],
        };
        let ron = ron::to_string(&l).unwrap();
        let back: DraftEntity = ron::from_str(&ron).unwrap();
        assert_eq!(l, back);
    }
}
