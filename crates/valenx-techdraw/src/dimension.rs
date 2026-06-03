//! Drawing dimensions — linear, angular, radial, diameter.
//!
//! A [`Dimension`] is purely data: endpoints, offset, computed value.
//! The export crates (SVG / PDF / DXF) own the rendering — they
//! decide how to draw extension lines, arrowheads, and text. Keeping
//! the dimension as data means the same drawing renders identically
//! across all three formats.

use serde::{Deserialize, Serialize};

/// Return type of [`Dimension::render_segments`]:
/// `(line segments, (label_x, label_y, label_text))`.
pub type RenderedDimension = (Vec<[(f64, f64); 2]>, (f64, f64, String));

/// One annotated dimension overlaid on a drawing view.
///
/// All coordinates are in the **drawing's sheet frame** (millimeters
/// from the bottom-left of the paper). The export layer doesn't have
/// to know which view a dimension "belongs to" — the values are
/// already in sheet space.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub enum Dimension {
    /// Linear distance between two points. `offset` is the
    /// perpendicular distance from the witness line to the dimension
    /// line (positive = above/right depending on the segment's
    /// orientation). `value` is the literal label, decoupled from
    /// the geometric distance so it survives unit conversion / drift
    /// (callers compute it once).
    Linear {
        /// First endpoint in mm.
        from: [f64; 2],
        /// Second endpoint in mm.
        to: [f64; 2],
        /// Perpendicular offset of the dimension line.
        offset: f64,
        /// Label value (mm).
        value: f64,
    },
    /// Angle between two rays sharing a vertex. `value` is in
    /// degrees so the renderer doesn't have to do unit conversion.
    Angular {
        /// Shared vertex of the two rays.
        vertex: [f64; 2],
        /// Point along the first ray.
        a: [f64; 2],
        /// Point along the second ray.
        b: [f64; 2],
        /// Perpendicular offset of the arc from `vertex` (mm).
        offset: f64,
        /// Label value (degrees).
        value: f64,
    },
    /// Radial dimension on a circle / arc.
    Radial {
        /// Center of the circle/arc.
        center: [f64; 2],
        /// Radius in mm.
        radius: f64,
        /// Position of the label text.
        label_pos: [f64; 2],
        /// Label value (mm) — usually = `radius`.
        value: f64,
    },
    /// Diameter dimension on a circle.
    Diameter {
        /// Center of the circle.
        center: [f64; 2],
        /// Radius in mm.
        radius: f64,
        /// Position of the label text.
        label_pos: [f64; 2],
        /// Label value (mm) — usually = `2 * radius`.
        value: f64,
    },
}

impl Dimension {
    /// Short label for UI display.
    pub fn kind(&self) -> &'static str {
        match self {
            Dimension::Linear { .. } => "Linear",
            Dimension::Angular { .. } => "Angular",
            Dimension::Radial { .. } => "Radial",
            Dimension::Diameter { .. } => "Diameter",
        }
    }

    /// Numeric value (mm or degrees as appropriate).
    pub fn value(&self) -> f64 {
        match self {
            Dimension::Linear { value, .. }
            | Dimension::Angular { value, .. }
            | Dimension::Radial { value, .. }
            | Dimension::Diameter { value, .. } => *value,
        }
    }

    /// Render placeholder — Task 17 in the plan keeps the dimension
    /// as pure data; the export crates (SVG / PDF / DXF) call back
    /// into [`Dimension::render_segments`] for the geometric primitives.
    ///
    /// This method exists so that future renderers (egui overlay in
    /// the viewport) have an API hook. For now it returns the same
    /// segments that the SVG / PDF exporters would draw.
    pub fn render(
        &self,
        segments: &mut Vec<[(f64, f64); 2]>,
        labels: &mut Vec<(f64, f64, String)>,
    ) {
        let (lines, label) = self.render_segments();
        segments.extend(lines);
        labels.push(label);
    }

    /// Compute the line segments + label position/text that depict
    /// this dimension. Used by every exporter; pulled out so the
    /// arrow / extension-line geometry stays in one place.
    pub fn render_segments(&self) -> RenderedDimension {
        let mut out: Vec<[(f64, f64); 2]> = Vec::new();
        match *self {
            Dimension::Linear {
                from,
                to,
                offset,
                value,
            } => {
                let dx = to[0] - from[0];
                let dy = to[1] - from[1];
                let len = (dx * dx + dy * dy).sqrt().max(1e-9);
                let nx = -dy / len;
                let ny = dx / len;
                let p0 = (from[0] + nx * offset, from[1] + ny * offset);
                let p1 = (to[0] + nx * offset, to[1] + ny * offset);
                // Dimension line.
                out.push([p0, p1]);
                // Witness lines (extension lines).
                out.push([(from[0], from[1]), p0]);
                out.push([(to[0], to[1]), p1]);
                // Arrowheads — small V at each end.
                let arrow_len = 2.0;
                let arrow_w = 0.6;
                let head = |p: (f64, f64), dir: f64, lines: &mut Vec<[(f64, f64); 2]>| {
                    let ux = (p1.0 - p0.0) / len * dir;
                    let uy = (p1.1 - p0.1) / len * dir;
                    let tip = (p.0 + ux * 0.0, p.1 + uy * 0.0);
                    let base = (p.0 - ux * arrow_len, p.1 - uy * arrow_len);
                    let l = (base.0 - nx * arrow_w, base.1 - ny * arrow_w);
                    let r = (base.0 + nx * arrow_w, base.1 + ny * arrow_w);
                    lines.push([tip, l]);
                    lines.push([tip, r]);
                };
                head(p0, -1.0, &mut out);
                head(p1, 1.0, &mut out);
                // Label centered on dimension line, offset perpendicular.
                let mid = ((p0.0 + p1.0) * 0.5, (p0.1 + p1.1) * 0.5);
                let label_pos = (mid.0 + nx * 2.0, mid.1 + ny * 2.0);
                (out, (label_pos.0, label_pos.1, format!("{value:.2}")))
            }
            Dimension::Angular {
                vertex,
                a,
                b,
                offset,
                value,
            } => {
                // Approximate the arc with a short polyline. Two
                // rays from vertex through a and b, plus a single
                // chord at `offset` distance from vertex.
                let da = vec_norm([a[0] - vertex[0], a[1] - vertex[1]]);
                let db = vec_norm([b[0] - vertex[0], b[1] - vertex[1]]);
                out.push([
                    (vertex[0], vertex[1]),
                    (
                        vertex[0] + da[0] * offset * 1.2,
                        vertex[1] + da[1] * offset * 1.2,
                    ),
                ]);
                out.push([
                    (vertex[0], vertex[1]),
                    (
                        vertex[0] + db[0] * offset * 1.2,
                        vertex[1] + db[1] * offset * 1.2,
                    ),
                ]);
                // Arc chord (poor man's arc — single straight line).
                let arc_a = (vertex[0] + da[0] * offset, vertex[1] + da[1] * offset);
                let arc_b = (vertex[0] + db[0] * offset, vertex[1] + db[1] * offset);
                out.push([arc_a, arc_b]);
                let mid = ((arc_a.0 + arc_b.0) * 0.5, (arc_a.1 + arc_b.1) * 0.5);
                (out, (mid.0, mid.1, format!("{value:.1}°")))
            }
            Dimension::Radial {
                center,
                radius,
                label_pos,
                value,
            } => {
                // Single radial line from center to a point on the
                // circle pointing toward `label_pos`.
                let dx = label_pos[0] - center[0];
                let dy = label_pos[1] - center[1];
                let len = (dx * dx + dy * dy).sqrt().max(1e-9);
                let edge = (center[0] + dx / len * radius, center[1] + dy / len * radius);
                out.push([(center[0], center[1]), edge]);
                out.push([edge, (label_pos[0], label_pos[1])]);
                (out, (label_pos[0], label_pos[1], format!("R{value:.2}")))
            }
            Dimension::Diameter {
                center,
                radius,
                label_pos,
                value,
            } => {
                // Two opposite radii through center toward label_pos.
                let dx = label_pos[0] - center[0];
                let dy = label_pos[1] - center[1];
                let len = (dx * dx + dy * dy).sqrt().max(1e-9);
                let ex = dx / len;
                let ey = dy / len;
                let a = (center[0] + ex * radius, center[1] + ey * radius);
                let b = (center[0] - ex * radius, center[1] - ey * radius);
                out.push([a, b]);
                out.push([a, (label_pos[0], label_pos[1])]);
                (out, (label_pos[0], label_pos[1], format!("⌀{value:.2}")))
            }
        }
    }
}

fn vec_norm(v: [f64; 2]) -> [f64; 2] {
    let l = (v[0] * v[0] + v[1] * v[1]).sqrt().max(1e-9);
    [v[0] / l, v[1] / l]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn linear_renders_dimension_witness_and_arrowheads() {
        let d = Dimension::Linear {
            from: [0.0, 0.0],
            to: [10.0, 0.0],
            offset: 5.0,
            value: 10.0,
        };
        let (segs, (_, _, label)) = d.render_segments();
        // Dimension line + 2 witness + 4 arrowhead lines = 7.
        assert_eq!(segs.len(), 7);
        assert!(label.contains("10.00"));
    }

    #[test]
    fn radial_label_includes_r_prefix() {
        let d = Dimension::Radial {
            center: [0.0, 0.0],
            radius: 3.5,
            label_pos: [10.0, 0.0],
            value: 3.5,
        };
        let (_, (_, _, label)) = d.render_segments();
        assert!(label.starts_with('R'));
        assert!(label.contains("3.50"));
    }

    #[test]
    fn diameter_label_includes_diameter_glyph() {
        let d = Dimension::Diameter {
            center: [0.0, 0.0],
            radius: 5.0,
            label_pos: [10.0, 0.0],
            value: 10.0,
        };
        let (_, (_, _, label)) = d.render_segments();
        assert!(label.starts_with('⌀'));
    }

    #[test]
    fn angular_label_includes_degree_glyph() {
        let d = Dimension::Angular {
            vertex: [0.0, 0.0],
            a: [10.0, 0.0],
            b: [0.0, 10.0],
            offset: 5.0,
            value: 90.0,
        };
        let (_, (_, _, label)) = d.render_segments();
        assert!(label.ends_with('°'));
    }

    #[test]
    fn kind_label_is_distinct_per_variant() {
        let cases = [
            Dimension::Linear {
                from: [0.0; 2],
                to: [0.0; 2],
                offset: 0.0,
                value: 0.0,
            }
            .kind(),
            Dimension::Angular {
                vertex: [0.0; 2],
                a: [0.0; 2],
                b: [0.0; 2],
                offset: 0.0,
                value: 0.0,
            }
            .kind(),
            Dimension::Radial {
                center: [0.0; 2],
                radius: 0.0,
                label_pos: [0.0; 2],
                value: 0.0,
            }
            .kind(),
            Dimension::Diameter {
                center: [0.0; 2],
                radius: 0.0,
                label_pos: [0.0; 2],
                value: 0.0,
            }
            .kind(),
        ];
        let set: std::collections::HashSet<_> = cases.iter().collect();
        assert_eq!(set.len(), cases.len());
    }
}
