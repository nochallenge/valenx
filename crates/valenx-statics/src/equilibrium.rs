//! Coplanar force systems and the three planar equilibrium equations.
//!
//! A [`ForceSystem`] is a bag of [`AppliedForce`]s acting in the `xy`
//! plane. From it you can read the **resultant**:
//!
//! - the net force `(sum Fx, sum Fy)` ([`ForceSystem::resultant_force`]),
//! - the net moment `sum M` about any chosen pivot
//!   ([`ForceSystem::resultant_moment`]).
//!
//! A rigid body is in **static equilibrium** when all three vanish
//! simultaneously:
//!
//! ```text
//! sum Fx = 0,    sum Fy = 0,    sum M = 0.
//! ```
//!
//! A key theorem of planar statics: if `sum Fx = sum Fy = 0` *and*
//! `sum M = 0` about **one** point, then `sum M = 0` about **every**
//! point. [`ForceSystem::is_in_equilibrium`] uses that fact — it tests
//! the moment about the origin only, which is sufficient once the force
//! sums are zero.

use crate::force::AppliedForce;
use nalgebra::{Point2, Vector2};

/// A collection of coplanar applied forces acting on one rigid body.
#[derive(Debug, Clone, Default, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct ForceSystem {
    /// The applied forces, in no particular order.
    pub forces: Vec<AppliedForce>,
}

impl ForceSystem {
    /// An empty system (no forces). Its resultant is the zero force and
    /// zero moment, so it is trivially in equilibrium.
    #[must_use]
    pub fn new() -> Self {
        Self { forces: Vec::new() }
    }

    /// Build a system from an existing list of forces.
    #[must_use]
    pub fn from_forces(forces: Vec<AppliedForce>) -> Self {
        Self { forces }
    }

    /// Append a force and return `&mut self` for chaining.
    pub fn push(&mut self, f: AppliedForce) -> &mut Self {
        self.forces.push(f);
        self
    }

    /// The resultant force `(sum Fx, sum Fy)` of the whole system.
    #[must_use]
    pub fn resultant_force(&self) -> Vector2<f64> {
        self.forces
            .iter()
            .fold(Vector2::zeros(), |acc, f| acc + f.vector)
    }

    /// `sum Fx` — the net horizontal force.
    #[must_use]
    pub fn sum_fx(&self) -> f64 {
        self.resultant_force().x
    }

    /// `sum Fy` — the net vertical force.
    #[must_use]
    pub fn sum_fy(&self) -> f64 {
        self.resultant_force().y
    }

    /// The resultant moment `sum M` of the system about `pivot`
    /// (counter-clockwise positive).
    #[must_use]
    pub fn resultant_moment(&self, pivot: Point2<f64>) -> f64 {
        self.forces.iter().map(|f| f.moment_about(pivot)).sum()
    }

    /// `sum M` about the origin — a convenience for the common case.
    #[must_use]
    pub fn sum_moment_origin(&self) -> f64 {
        self.resultant_moment(Point2::origin())
    }

    /// `true` if the system is in static equilibrium to within `tol`:
    /// `|sum Fx| <= tol`, `|sum Fy| <= tol`, and `|sum M| <= tol` about
    /// the origin.
    ///
    /// Testing the moment about the origin alone is sufficient: once the
    /// force resultant is zero, the moment is the same about every
    /// point, so a zero moment about one point implies a zero moment
    /// everywhere.
    #[must_use]
    pub fn is_in_equilibrium(&self, tol: f64) -> bool {
        let r = self.resultant_force();
        r.x.abs() <= tol && r.y.abs() <= tol && self.sum_moment_origin().abs() <= tol
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    #[test]
    fn empty_system_is_in_equilibrium() {
        let sys = ForceSystem::new();
        assert!(sys.is_in_equilibrium(EPS));
        assert!(sys.sum_fx().abs() < EPS);
        assert!(sys.sum_fy().abs() < EPS);
        assert!(sys.sum_moment_origin().abs() < EPS);
    }

    #[test]
    fn opposing_collinear_forces_balance() {
        // Two equal-and-opposite forces sharing a line of action:
        // sum F = 0 and sum M = 0 -> equilibrium.
        let mut sys = ForceSystem::new();
        sys.push(AppliedForce::new(0.0, 0.0, 10.0, 0.0))
            .push(AppliedForce::new(0.0, 0.0, -10.0, 0.0));
        assert!((sys.sum_fx()).abs() < EPS, "sum Fx = {}", sys.sum_fx());
        assert!(sys.is_in_equilibrium(EPS));
    }

    #[test]
    fn couple_has_zero_force_but_nonzero_moment() {
        // A force couple: two equal-and-opposite forces on parallel but
        // distinct lines. Net force is zero, but the moment is not, so
        // the body is NOT in equilibrium.
        let mut sys = ForceSystem::new();
        sys.push(AppliedForce::new(0.0, 1.0, 5.0, 0.0)) // +x at y = +1
            .push(AppliedForce::new(0.0, -1.0, -5.0, 0.0)); // -x at y = -1
        assert!(sys.sum_fx().abs() < EPS);
        assert!(sys.sum_fy().abs() < EPS);
        // Each contributes -5 N*m about origin: M = -(r_y * F_x).
        // Top: r=(0,1),F=(5,0) -> 0*0 - 1*5 = -5. Bottom likewise -5.
        let m = sys.sum_moment_origin();
        assert!((m + 10.0).abs() < EPS, "got {m}");
        assert!(!sys.is_in_equilibrium(EPS));
    }

    #[test]
    fn moment_independent_of_pivot_when_force_balanced() {
        // A pure couple has the same moment about every point. Verify
        // the moment about the origin equals the moment about (7, -3).
        let mut sys = ForceSystem::new();
        sys.push(AppliedForce::new(0.0, 1.0, 5.0, 0.0))
            .push(AppliedForce::new(0.0, -1.0, -5.0, 0.0));
        let m0 = sys.resultant_moment(Point2::origin());
        let m1 = sys.resultant_moment(Point2::new(7.0, -3.0));
        assert!((m0 - m1).abs() < EPS, "m0 {m0}, m1 {m1}");
    }

    #[test]
    fn three_force_equilibrium() {
        // Classic worked example: three concurrent forces summing to
        // zero. (3,4) + (-3,0) + (0,-4) = (0,0); all act at the origin
        // so every moment is zero.
        let mut sys = ForceSystem::new();
        sys.push(AppliedForce::new(0.0, 0.0, 3.0, 4.0))
            .push(AppliedForce::new(0.0, 0.0, -3.0, 0.0))
            .push(AppliedForce::new(0.0, 0.0, 0.0, -4.0));
        assert!(sys.is_in_equilibrium(EPS));
    }
}
