//! 2-D Ellipse + Elliptical Arc primitives for the sketcher.
//!
//! Phase 12A.
//!
//! An [`Ellipse2`] is stored by:
//! - a `center` [`Point2`],
//! - a `major_axis_var` pair (one variable each for x and y components
//!   of the semi-major axis vector — its magnitude is the semi-major
//!   length and its direction sets the ellipse's rotation),
//! - a `minor_radius_var` scalar (the semi-minor magnitude).
//!
//! Parametric evaluation:
//! `P(t) = center + cos(t) * major_axis + sin(t) * minor_axis`
//! where `minor_axis = perp(major_axis_unit) * minor_radius`.

use serde::{Deserialize, Serialize};

use crate::geom::Point2;

/// Full ellipse on the sketch plane.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash, Serialize, Deserialize)]
pub struct Ellipse2 {
    /// Centre point.
    pub center: Point2,
    /// Index for the x component of the semi-major axis vector.
    pub major_x_var: usize,
    /// Index for the y component of the semi-major axis vector.
    pub major_y_var: usize,
    /// Index for the semi-minor magnitude.
    pub minor_radius_var: usize,
}

impl Ellipse2 {
    /// Read the (cx, cy) of the centre.
    pub fn center_xy(&self, vars: &[f64]) -> (f64, f64) {
        self.center.read(vars)
    }

    /// Semi-major axis vector (mx, my) in world units. Out-of-range
    /// handles read as `0.0` (R33 H1 defense-in-depth — see
    /// [`crate::geom::Point2::read`]).
    pub fn major_axis(&self, vars: &[f64]) -> (f64, f64) {
        (
            vars.get(self.major_x_var).copied().unwrap_or(0.0),
            vars.get(self.major_y_var).copied().unwrap_or(0.0),
        )
    }

    /// Semi-major length (magnitude of major axis).
    pub fn major_radius(&self, vars: &[f64]) -> f64 {
        let (mx, my) = self.major_axis(vars);
        (mx * mx + my * my).sqrt()
    }

    /// Semi-minor magnitude. Out-of-range handle reads as `0.0`
    /// (R33 H1 defense-in-depth — see [`crate::geom::Point2::read`]).
    pub fn minor_radius(&self, vars: &[f64]) -> f64 {
        vars.get(self.minor_radius_var).copied().unwrap_or(0.0)
    }

    /// Evaluate at parameter `t ∈ [0, 2π]`. Returns (x, y) in world.
    pub fn evaluate(&self, vars: &[f64], t: f64) -> (f64, f64) {
        let (cx, cy) = self.center_xy(vars);
        let (mx, my) = self.major_axis(vars);
        let m_mag = (mx * mx + my * my).sqrt();
        if m_mag < 1e-30 {
            return (cx, cy);
        }
        // Minor axis direction = +90° rotation of major axis unit vec.
        let ux = mx / m_mag;
        let uy = my / m_mag;
        let vx = -uy;
        let vy = ux;
        let minor = self.minor_radius(vars);
        let ct = t.cos();
        let st = t.sin();
        (
            cx + ct * mx + st * vx * minor,
            cy + ct * my + st * vy * minor,
        )
    }

    /// (start, end) for a full ellipse — both at t=0 (same point).
    pub fn endpoints(&self, vars: &[f64]) -> ((f64, f64), (f64, f64)) {
        let p = self.evaluate(vars, 0.0);
        (p, p)
    }
}

/// Arc segment of an ellipse — same as [`Ellipse2`] plus start/end
/// angle variables (radians, CCW from major axis).
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash, Serialize, Deserialize)]
pub struct EllipticalArc2 {
    /// Underlying ellipse.
    pub ellipse: Ellipse2,
    /// Start-angle variable index.
    pub start_angle_var: usize,
    /// End-angle variable index.
    pub end_angle_var: usize,
}

impl EllipticalArc2 {
    /// (start_angle, end_angle) in radians. Out-of-range handles read
    /// as `0.0` (R33 H1 defense-in-depth — see
    /// [`crate::geom::Point2::read`]).
    pub fn angles(&self, vars: &[f64]) -> (f64, f64) {
        (
            vars.get(self.start_angle_var).copied().unwrap_or(0.0),
            vars.get(self.end_angle_var).copied().unwrap_or(0.0),
        )
    }

    /// Sweep = end - start. May be negative.
    pub fn sweep(&self, vars: &[f64]) -> f64 {
        let (s, e) = self.angles(vars);
        e - s
    }

    /// Evaluate at parameter `t` (radians). Same as `ellipse.evaluate`.
    pub fn evaluate(&self, vars: &[f64], t: f64) -> (f64, f64) {
        self.ellipse.evaluate(vars, t)
    }

    /// (start, end) points at the start/end angles.
    pub fn endpoints(&self, vars: &[f64]) -> ((f64, f64), (f64, f64)) {
        let (sa, ea) = self.angles(vars);
        (
            self.ellipse.evaluate(vars, sa),
            self.ellipse.evaluate(vars, ea),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sketch::Sketch;
    use std::f64::consts::{PI, TAU};

    fn unit_circle_ellipse(s: &mut Sketch) -> Ellipse2 {
        let center = s.add_point(0.0, 0.0);
        let cpt = s.point_at(center).unwrap();
        let mx = s.alloc_var(1.0);
        let my = s.alloc_var(0.0);
        let mr = s.alloc_var(1.0);
        Ellipse2 {
            center: cpt,
            major_x_var: mx,
            major_y_var: my,
            minor_radius_var: mr,
        }
    }

    // R33 H1 defense-in-depth: out-of-range major/minor/angle handles
    // on the ellipse read path degrade to 0.0 instead of panicking.
    #[test]
    fn ellipse_major_axis_with_out_of_range_handle_does_not_panic() {
        let e = Ellipse2 {
            center: Point2 { x_var: 0, y_var: 0 },
            major_x_var: 999,
            major_y_var: 999,
            minor_radius_var: 999,
        };
        let vars = vec![0.0];
        assert_eq!(e.major_axis(&vars), (0.0, 0.0));
        assert_eq!(e.minor_radius(&vars), 0.0);
        // evaluate falls through the degenerate-major branch → center.
        assert_eq!(e.evaluate(&vars, 1.23), (0.0, 0.0));
    }

    #[test]
    fn elliptical_arc_angles_with_out_of_range_handle_does_not_panic() {
        let arc = EllipticalArc2 {
            ellipse: Ellipse2 {
                center: Point2 { x_var: 0, y_var: 0 },
                major_x_var: 0,
                major_y_var: 0,
                minor_radius_var: 0,
            },
            start_angle_var: 999,
            end_angle_var: 999,
        };
        let vars = vec![0.0];
        assert_eq!(arc.angles(&vars), (0.0, 0.0));
        assert_eq!(arc.sweep(&vars), 0.0);
    }

    #[test]
    fn ellipse_at_zero_angle_is_at_center_plus_major() {
        let mut s = Sketch::new();
        let e = unit_circle_ellipse(&mut s);
        let (x, y) = e.evaluate(&s.vars, 0.0);
        assert!((x - 1.0).abs() < 1e-12);
        assert!(y.abs() < 1e-12);
    }

    #[test]
    fn ellipse_at_pi_is_opposite_major() {
        let mut s = Sketch::new();
        let e = unit_circle_ellipse(&mut s);
        let (x, y) = e.evaluate(&s.vars, PI);
        assert!((x + 1.0).abs() < 1e-12);
        assert!(y.abs() < 1e-12);
    }

    #[test]
    fn ellipse_at_half_pi_is_at_minor_radius() {
        let mut s = Sketch::new();
        let e = unit_circle_ellipse(&mut s);
        let (x, y) = e.evaluate(&s.vars, PI / 2.0);
        assert!(x.abs() < 1e-12);
        assert!((y - 1.0).abs() < 1e-12);
    }

    #[test]
    fn elliptical_arc_sweep_is_end_minus_start() {
        let mut s = Sketch::new();
        let e = unit_circle_ellipse(&mut s);
        let sa = s.alloc_var(0.0);
        let ea = s.alloc_var(PI);
        let arc = EllipticalArc2 {
            ellipse: e,
            start_angle_var: sa,
            end_angle_var: ea,
        };
        assert!((arc.sweep(&s.vars) - PI).abs() < 1e-12);
    }

    #[test]
    fn elliptical_arc_endpoints_match_angles() {
        let mut s = Sketch::new();
        let e = unit_circle_ellipse(&mut s);
        let sa = s.alloc_var(0.0);
        let ea = s.alloc_var(PI);
        let arc = EllipticalArc2 {
            ellipse: e,
            start_angle_var: sa,
            end_angle_var: ea,
        };
        let ((sx, sy), (ex, ey)) = arc.endpoints(&s.vars);
        assert!((sx - 1.0).abs() < 1e-12 && sy.abs() < 1e-12);
        assert!((ex + 1.0).abs() < 1e-12 && ey.abs() < 1e-12);
        let _ = TAU; // silence unused
    }
}
