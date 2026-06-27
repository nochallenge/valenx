//! A simple axis-aligned box boundary (reflect / clamp).
//!
//! The cheapest workable container for a particle fluid: an axis-aligned
//! bounding box. After integration any particle that has left the box is
//! **clamped** back to the nearest wall, and the velocity component normal to
//! that wall is **reflected** and damped by a restitution coefficient
//! `0 ≤ e ≤ 1` (`e = 0` fully absorbs the normal velocity, `e = 1` is a
//! perfectly elastic bounce). This is a penalty-free positional boundary — it is
//! enough to keep a column of fluid in a tank for the hydrostatic and
//! conservation tests, but it is **not** a physical no-slip wall and does not
//! contribute boundary particles to the SPH density sum (a known limitation,
//! documented at the crate root).

use nalgebra::Vector3;

use crate::error::FluidError;
use crate::particle::Particle;

/// An axis-aligned box `[min, max]` with a wall restitution coefficient.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BoxBoundary {
    min: Vector3<f64>,
    max: Vector3<f64>,
    restitution: f64,
}

impl BoxBoundary {
    /// Build a box from its `min` and `max` corners and a restitution
    /// coefficient.
    ///
    /// # Errors
    /// - [`FluidError::NonFinite`] if any corner component or the restitution is
    ///   non-finite.
    /// - [`FluidError::InvalidDomain`] if `min` is not strictly less than `max`
    ///   on every axis.
    /// - [`FluidError::InvalidConfig`] if `restitution` is outside `[0, 1]`.
    pub fn new(min: Vector3<f64>, max: Vector3<f64>, restitution: f64) -> Result<Self, FluidError> {
        if !min.iter().chain(max.iter()).all(|c| c.is_finite()) {
            return Err(FluidError::NonFinite("box corner".into()));
        }
        if !restitution.is_finite() {
            return Err(FluidError::NonFinite("restitution".into()));
        }
        for axis in 0..3 {
            if min[axis] >= max[axis] {
                return Err(FluidError::InvalidDomain(format!(
                    "box min must be < max on every axis (axis {axis}: {} !< {})",
                    min[axis], max[axis]
                )));
            }
        }
        if !(0.0..=1.0).contains(&restitution) {
            return Err(FluidError::InvalidConfig(format!(
                "restitution must be in [0, 1], got {restitution}"
            )));
        }
        Ok(Self {
            min,
            max,
            restitution,
        })
    }

    /// The minimum corner.
    #[must_use]
    pub fn min(&self) -> Vector3<f64> {
        self.min
    }

    /// The maximum corner.
    #[must_use]
    pub fn max(&self) -> Vector3<f64> {
        self.max
    }

    /// The wall restitution coefficient.
    #[must_use]
    pub fn restitution(&self) -> f64 {
        self.restitution
    }

    /// Whether a position lies inside (or on the boundary of) the box.
    #[must_use]
    pub fn contains(&self, p: Vector3<f64>) -> bool {
        (0..3).all(|a| p[a] >= self.min[a] && p[a] <= self.max[a])
    }

    /// Clamp a particle back inside the box, reflecting and damping the normal
    /// velocity component on any wall it has crossed.
    ///
    /// Returns `true` if the particle was outside and got corrected.
    pub fn enforce(&self, particle: &mut Particle) -> bool {
        let mut corrected = false;
        for axis in 0..3 {
            if particle.position[axis] < self.min[axis] {
                particle.position[axis] = self.min[axis];
                // Reflect the inward-normal velocity only if it is still moving
                // outward (into the wall).
                if particle.velocity[axis] < 0.0 {
                    particle.velocity[axis] *= -self.restitution;
                }
                corrected = true;
            } else if particle.position[axis] > self.max[axis] {
                particle.position[axis] = self.max[axis];
                if particle.velocity[axis] > 0.0 {
                    particle.velocity[axis] *= -self.restitution;
                }
                corrected = true;
            }
        }
        corrected
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn unit_box() -> BoxBoundary {
        BoxBoundary::new(
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 1.0, 1.0),
            0.5,
        )
        .unwrap()
    }

    #[test]
    fn rejects_degenerate_box() {
        // min == max on an axis.
        assert!(BoxBoundary::new(
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 1.0),
            0.5
        )
        .is_err());
        // min > max.
        assert!(BoxBoundary::new(
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 1.0),
            0.5
        )
        .is_err());
        // restitution out of range.
        assert!(BoxBoundary::new(
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 1.0, 1.0),
            1.5
        )
        .is_err());
        // non-finite corner.
        assert!(BoxBoundary::new(
            Vector3::new(0.0, 0.0, f64::NAN),
            Vector3::new(1.0, 1.0, 1.0),
            0.5
        )
        .is_err());
    }

    #[test]
    fn contains_is_inclusive_of_walls() {
        let b = unit_box();
        assert!(b.contains(Vector3::new(0.5, 0.5, 0.5)));
        assert!(b.contains(Vector3::new(0.0, 1.0, 0.5)));
        assert!(!b.contains(Vector3::new(-0.01, 0.5, 0.5)));
        assert!(!b.contains(Vector3::new(0.5, 1.01, 0.5)));
    }

    #[test]
    fn enforce_clamps_and_reflects() {
        let b = unit_box();
        // Below the floor moving down.
        let mut p = Particle::new(Vector3::new(0.5, -0.2, 0.5), Vector3::new(0.0, -2.0, 0.0));
        assert!(b.enforce(&mut p));
        assert_eq!(p.position.y, 0.0); // clamped to the floor
        assert_eq!(p.velocity.y, 1.0); // reflected & damped: -0.5 * -2.0
    }

    #[test]
    fn enforce_leaves_interior_particles_untouched() {
        let b = unit_box();
        let mut p = Particle::new(Vector3::new(0.5, 0.5, 0.5), Vector3::new(1.0, 1.0, 1.0));
        let before = p;
        assert!(!b.enforce(&mut p));
        assert_eq!(p, before);
    }

    #[test]
    fn enforce_does_not_re_reflect_when_already_moving_inward() {
        // A particle pressed against the floor but already moving up should be
        // clamped (if outside) without flipping its (already inward) velocity.
        let b = unit_box();
        let mut p = Particle::new(Vector3::new(0.5, -0.1, 0.5), Vector3::new(0.0, 3.0, 0.0));
        b.enforce(&mut p);
        assert_eq!(p.position.y, 0.0);
        assert_eq!(p.velocity.y, 3.0); // unchanged — it was already moving inward
    }
}
