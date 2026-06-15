//! Planar points, forces, and the scalar moment `M = r x F`.
//!
//! Everything in 2D statics is built from two primitives, taken
//! directly from `nalgebra`:
//!
//! - a **point** `nalgebra::Point2<f64>` — a position in the `xy`
//!   plane;
//! - a **force** `nalgebra::Vector2<f64>` — a free vector `(Fx, Fy)`.
//!
//! An [`AppliedForce`] bundles a force vector with its point of
//! application so the moment computation needs no separate position
//! argument.
//!
//! The one non-obvious operation is the **moment** of a force about a
//! point. In a plane every position vector and every force lies in the
//! `xy` plane, so the cross product `r x F` points purely along `z`. We
//! therefore represent the moment as the signed scalar
//!
//! ```text
//! M = r_x * F_y - r_y * F_x
//! ```
//!
//! the `z`-component of `r x F`, where `r` is the vector from the
//! reference (pivot) point to the force's point of application. The
//! sign follows the right-hand rule: positive is counter-clockwise
//! (out of the page).

use nalgebra::{Point2, Vector2};

/// A force applied at a point in the plane.
///
/// `application` is the point of application; `vector` is the force
/// `(Fx, Fy)` in consistent units (e.g. newtons). The struct bundles
/// the two so a [`moment_about`](AppliedForce::moment_about) needs no
/// separate position argument.
#[derive(Debug, Clone, Copy, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct AppliedForce {
    /// Point of application in the `xy` plane.
    pub application: Point2<f64>,
    /// Force vector `(Fx, Fy)`.
    pub vector: Vector2<f64>,
}

impl AppliedForce {
    /// Build an [`AppliedForce`] from raw coordinates and components.
    ///
    /// `(x, y)` is the point of application; `(fx, fy)` is the force.
    #[must_use]
    pub fn new(x: f64, y: f64, fx: f64, fy: f64) -> Self {
        Self {
            application: Point2::new(x, y),
            vector: Vector2::new(fx, fy),
        }
    }

    /// `true` if every coordinate and component is finite
    /// (neither `NaN` nor infinite).
    #[must_use]
    pub fn is_finite(&self) -> bool {
        self.application.x.is_finite()
            && self.application.y.is_finite()
            && self.vector.x.is_finite()
            && self.vector.y.is_finite()
    }

    /// The scalar (out-of-plane, `z`) moment of this force about the
    /// point `pivot`.
    ///
    /// Computes `M = r_x * F_y - r_y * F_x` with `r = application -
    /// pivot`. Positive is counter-clockwise.
    #[must_use]
    pub fn moment_about(&self, pivot: Point2<f64>) -> f64 {
        let r = self.application - pivot;
        moment_z(r, self.vector)
    }
}

/// The `z`-component of the cross product `r x F` for planar vectors.
///
/// `moment_z(r, f) = r_x * f_y - r_y * f_x`. This is the scalar moment
/// of a force `f` whose line of action passes through the tip of the
/// position vector `r` measured from the pivot. Positive is
/// counter-clockwise (right-hand rule, `+z` out of the page).
#[must_use]
pub fn moment_z(r: Vector2<f64>, f: Vector2<f64>) -> f64 {
    r.x * f.y - r.y * f.x
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-12;

    #[test]
    fn moment_arm_unit_force() {
        // A unit upward force 3 m to the right of the pivot makes a
        // +3 N*m counter-clockwise moment (M = d * F = 3 * 1).
        let f = AppliedForce::new(3.0, 0.0, 0.0, 1.0);
        let m = f.moment_about(Point2::origin());
        assert!((m - 3.0).abs() < EPS, "got {m}");
    }

    #[test]
    fn moment_sign_is_right_hand_rule() {
        // Upward force to the +x side -> CCW (+). Downward force to the
        // +x side -> CW (-). Magnitudes equal, signs opposite.
        let up = AppliedForce::new(2.0, 0.0, 0.0, 5.0);
        let down = AppliedForce::new(2.0, 0.0, 0.0, -5.0);
        let m_up = up.moment_about(Point2::origin());
        let m_down = down.moment_about(Point2::origin());
        assert!((m_up - 10.0).abs() < EPS, "got {m_up}");
        assert!((m_down + 10.0).abs() < EPS, "got {m_down}");
    }

    #[test]
    fn force_through_pivot_has_zero_moment() {
        // A force whose line of action passes through the pivot exerts
        // no moment: collinear r and F -> cross product zero.
        let f = AppliedForce::new(4.0, 0.0, 7.0, 0.0); // along +x, applied on +x axis
        let m = f.moment_about(Point2::origin());
        assert!(m.abs() < EPS, "got {m}");
    }

    #[test]
    fn moment_about_arbitrary_pivot() {
        // Force (0, 10) at (5, 2); pivot at (1, 2).
        // r = (4, 0); M = 4*10 - 0*0 = 40.
        let f = AppliedForce::new(5.0, 2.0, 0.0, 10.0);
        let m = f.moment_about(Point2::new(1.0, 2.0));
        assert!((m - 40.0).abs() < EPS, "got {m}");
    }

    #[test]
    fn moment_z_matches_determinant() {
        let r = Vector2::new(3.0, 4.0);
        let f = Vector2::new(-2.0, 6.0);
        // 3*6 - 4*(-2) = 18 + 8 = 26.
        assert!((moment_z(r, f) - 26.0).abs() < EPS);
    }

    #[test]
    fn finiteness_detects_nan() {
        let good = AppliedForce::new(1.0, 2.0, 3.0, 4.0);
        assert!(good.is_finite());
        let bad = AppliedForce::new(f64::NAN, 0.0, 0.0, 0.0);
        assert!(!bad.is_finite());
    }
}
