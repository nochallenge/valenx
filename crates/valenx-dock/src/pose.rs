//! Ligand pose: translation + orientation + per-bond torsion angles.

use nalgebra::{UnitQuaternion, Vector3};

/// Full pose: position of the root atom plus orientation plus torsions.
#[derive(Clone, Debug, PartialEq)]
pub struct Pose {
    /// Translation applied to the root group's reference frame (Å).
    pub translation: Vector3<f64>,
    /// Global orientation of the root group.
    pub orientation: UnitQuaternion<f64>,
    /// One angle (radians) per rotatable bond, indexed by the
    /// non-root-group order in [`crate::ligand::Ligand::groups`].
    pub torsions: Vec<f64>,
}

impl Pose {
    /// Construct the identity pose: no translation, no rotation,
    /// all torsions zero.
    pub fn identity(n_torsions: usize) -> Self {
        Self {
            translation: Vector3::zeros(),
            orientation: UnitQuaternion::identity(),
            torsions: vec![0.0; n_torsions],
        }
    }

    /// Degrees of freedom = 3 (trans) + 3 (rot, parameterized as axis-angle) + N (tors).
    pub fn n_dofs(&self) -> usize {
        6 + self.torsions.len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn identity_has_no_movement() {
        let p = Pose::identity(3);
        assert_eq!(p.translation, Vector3::zeros());
        assert_eq!(p.orientation, UnitQuaternion::identity());
        assert_eq!(p.torsions, vec![0.0, 0.0, 0.0]);
        assert_eq!(p.n_dofs(), 9);
    }
}
