//! InternalAlignment — a point lies on (and is "owned by") the major
//! or minor axis endpoint of an ellipse.
//!
//! Phase 12B Task 17. This is FreeCAD's mechanism for exposing an
//! ellipse's geometric handles (`Ellipse.MajorEnd`,
//! `Ellipse.MinorEnd`) as ordinary point entities the user can dimension.
//!
//! Two residuals: `(p.x - tx, p.y - ty)` where `(tx, ty)` is one of the
//! ellipse's axis endpoints (major +/- or minor +/-) computed from its
//! `center`, `major_axis` and `minor_radius`.

use crate::geom::EntityId;
use crate::sketch::Sketch;

/// Which ellipse axis endpoint the point is aligned with.
#[derive(Copy, Clone, Debug, serde::Serialize, serde::Deserialize, PartialEq, Eq)]
pub enum AlignmentKind {
    /// Major-axis +radius endpoint (`center + major_axis`).
    MajorEnd,
    /// Major-axis -radius endpoint (`center - major_axis`).
    MajorStart,
    /// Minor-axis +radius endpoint (`center + perp(major_unit) * minor_radius`).
    MinorEnd,
    /// Minor-axis -radius endpoint.
    MinorStart,
}

/// Residuals: (p.x - target.x, p.y - target.y).
pub fn residuals(
    sketch: &Sketch,
    point: EntityId,
    ellipse: EntityId,
    kind: AlignmentKind,
    out: &mut [f64],
) {
    let (Ok(p), Ok(e)) = (sketch.point_at(point), sketch.ellipse_at(ellipse)) else {
        out[0] = 0.0;
        out[1] = 0.0;
        return;
    };
    let (px, py) = p.read(&sketch.vars);
    let (cx, cy) = e.center_xy(&sketch.vars);
    let (mx, my) = e.major_axis(&sketch.vars);
    let a = (mx * mx + my * my).sqrt();
    let b = e.minor_radius(&sketch.vars);
    let (tx, ty) = match kind {
        AlignmentKind::MajorEnd => (cx + mx, cy + my),
        AlignmentKind::MajorStart => (cx - mx, cy - my),
        AlignmentKind::MinorEnd => {
            if a < 1e-15 {
                (cx, cy)
            } else {
                let ux = mx / a;
                let uy = my / a;
                let vx = -uy;
                let vy = ux;
                (cx + vx * b, cy + vy * b)
            }
        }
        AlignmentKind::MinorStart => {
            if a < 1e-15 {
                (cx, cy)
            } else {
                let ux = mx / a;
                let uy = my / a;
                let vx = -uy;
                let vy = ux;
                (cx - vx * b, cy - vy * b)
            }
        }
    };
    out[0] = px - tx;
    out[1] = py - ty;
}

/// Jacobian: finite-difference. 7 candidate variables per row.
pub fn jacobian(
    sketch: &Sketch,
    point: EntityId,
    ellipse: EntityId,
    kind: AlignmentKind,
    triplets: &mut Vec<(usize, usize, f64)>,
) {
    let (Ok(p), Ok(e)) = (sketch.point_at(point), sketch.ellipse_at(ellipse)) else {
        return;
    };
    let vars = vec![
        p.x_var,
        p.y_var,
        e.center.x_var,
        e.center.y_var,
        e.major_x_var,
        e.major_y_var,
        e.minor_radius_var,
    ];
    let mut base = vec![0.0; 2];
    residuals(sketch, point, ellipse, kind, &mut base);
    let h = 1e-7;
    let mut perturbed = sketch.clone();
    for var in vars {
        let saved = perturbed.vars[var];
        perturbed.vars[var] = saved + h;
        let mut r = vec![0.0; 2];
        residuals(&perturbed, point, ellipse, kind, &mut r);
        for row in 0..2 {
            let d = (r[row] - base[row]) / h;
            if d.abs() > 1e-15 {
                triplets.push((row, var, d));
            }
        }
        perturbed.vars[var] = saved;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn residual_zero_for_major_end_on_axis_aligned_ellipse() {
        let mut s = Sketch::new();
        let c = s.add_point(0.0, 0.0);
        let e = s.add_ellipse(c, (2.0, 0.0), 1.0).unwrap();
        let p = s.add_point(2.0, 0.0);
        let mut out = vec![0.0; 2];
        residuals(&s, p, e, AlignmentKind::MajorEnd, &mut out);
        assert!(out[0].abs() < 1e-12 && out[1].abs() < 1e-12);
    }
}
