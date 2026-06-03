//! A single (position, orientation) pair used to instance a solid.

use nalgebra::{UnitQuaternion, Vector3};
use serde::{Deserialize, Serialize};

/// One instance pose. `position` is in world mm; `orientation` is the
/// rotation from the source solid's frame to the world.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Placement {
    /// World-space position.
    pub position: Vector3<f64>,
    /// World-space orientation.
    pub orientation: UnitQuaternion<f64>,
}

impl Placement {
    /// Identity placement at the origin.
    pub fn identity() -> Self {
        Self {
            position: Vector3::zeros(),
            orientation: UnitQuaternion::identity(),
        }
    }

    /// Placement at `position` with identity orientation.
    pub fn at(position: Vector3<f64>) -> Self {
        Self {
            position,
            orientation: UnitQuaternion::identity(),
        }
    }
}

impl Default for Placement {
    fn default() -> Self {
        Self::identity()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn identity_is_origin() {
        let p = Placement::identity();
        assert_eq!(p.position, Vector3::zeros());
    }

    #[test]
    fn at_uses_position() {
        let p = Placement::at(Vector3::new(1.0, 2.0, 3.0));
        assert_eq!(p.position, Vector3::new(1.0, 2.0, 3.0));
    }
}
