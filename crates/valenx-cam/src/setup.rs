//! Setup struct (Phase 17F) — stock + fixture + work-origin
//! transform, used for multi-side machining workflows.
//!
//! v1 stores a list of [`Setup`] per part. Each setup has its own
//! stock orientation (so flipping in the vise for the 2nd side
//! machining is one `Setup`), and a 4×4 transform that maps
//! setup-local coordinates back to part-world.

use nalgebra::{Matrix4, Vector3};
use serde::{Deserialize, Serialize};

use crate::fixture::Fixture;
use crate::stock::Stock;

/// One machining setup — stock placed in the fixture at a fixed
/// orientation, ready for some subset of operations.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Setup {
    /// Free-form label (`"OP10 first-side"`, `"OP20 second-side"`).
    pub label: String,
    /// Stock block as seen by the operations in this setup.
    pub stock: Stock,
    /// Fixture geometry for this setup.
    pub fixture: Fixture,
    /// 4×4 transform from setup-local to part-world coordinates.
    /// Identity = setup matches the part's reference frame.
    pub work_origin_transform: [[f64; 4]; 4],
}

impl Default for Setup {
    fn default() -> Self {
        Self {
            label: "Setup 1".into(),
            stock: Stock::default(),
            fixture: Fixture::new(),
            work_origin_transform: Matrix4::identity().into(),
        }
    }
}

impl Setup {
    /// New setup with the given label and an identity work-origin
    /// transform.
    pub fn new(label: impl Into<String>, stock: Stock) -> Self {
        Self {
            label: label.into(),
            stock,
            fixture: Fixture::new(),
            work_origin_transform: Matrix4::identity().into(),
        }
    }

    /// Reconstruct the work-origin transform as an nalgebra matrix.
    pub fn transform(&self) -> Matrix4<f64> {
        Matrix4::from(self.work_origin_transform)
    }

    /// Transform a setup-local point to part-world coordinates.
    pub fn local_to_world(&self, p: Vector3<f64>) -> Vector3<f64> {
        let m = self.transform();
        let h = m * nalgebra::Vector4::new(p.x, p.y, p.z, 1.0);
        Vector3::new(h.x, h.y, h.z)
    }
}

/// Collection of setups for a part. Operations reference a
/// `setup_index` into this list when running multi-side workflows.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct SetupSet {
    /// Ordered list of setups.
    pub setups: Vec<Setup>,
}

impl SetupSet {
    /// Empty setup set.
    pub fn new() -> Self {
        Self::default()
    }

    /// Append a setup. Returns its index.
    pub fn push(&mut self, setup: Setup) -> usize {
        let idx = self.setups.len();
        self.setups.push(setup);
        idx
    }

    /// Number of setups.
    pub fn len(&self) -> usize {
        self.setups.len()
    }

    /// `true` when the set is empty.
    pub fn is_empty(&self) -> bool {
        self.setups.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_setup_identity_transform() {
        let s = Setup::default();
        let m = s.transform();
        let p = Vector3::new(1.0, 2.0, 3.0);
        let q = s.local_to_world(p);
        assert!((q - p).norm() < 1e-9);
        // Confirm identity literal too.
        for i in 0..4 {
            for j in 0..4 {
                assert!((m[(i, j)] - if i == j { 1.0 } else { 0.0 }).abs() < 1e-9);
            }
        }
    }

    #[test]
    fn setup_set_push_returns_index() {
        let mut ss = SetupSet::new();
        let i0 = ss.push(Setup::new("op10", Stock::default()));
        let i1 = ss.push(Setup::new("op20", Stock::default()));
        assert_eq!(i0, 0);
        assert_eq!(i1, 1);
        assert_eq!(ss.len(), 2);
    }
}
