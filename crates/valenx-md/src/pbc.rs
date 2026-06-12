//! Periodic boundary conditions — the simulation box.
//!
//! **Roadmap feature 11.** A [`SimBox`] holds the three lattice
//! vectors **a**, **b**, **c** of the periodic cell as the columns of
//! a 3×3 matrix `h`. Two geometries are supported:
//!
//! - **Orthorhombic** — `h` is diagonal, the common rectangular box.
//! - **Triclinic** — a general lower-triangular `h` in the GROMACS
//!   convention (`a` along x; `b` in the xy-plane; `c` anywhere), as
//!   produced by truncated-octahedron and rhombic-dodecahedron
//!   solvent boxes.
//!
//! The central operation is the **minimum-image convention**: given a
//! displacement vector, return the shortest equivalent vector under
//! the lattice. For an orthorhombic box this is the textbook
//! component-wise wrap. For a triclinic box this crate uses the
//! standard fractional-coordinate round followed by a short search
//! over the 27 neighbouring images — exact whenever the cutoff is
//! below half the smallest box width, which is the usual MD
//! requirement.
//!
//! A box may also be **non-periodic** ([`SimBox::none`]); then the
//! minimum image is the identity and the volume is infinite.

use nalgebra::{Matrix3, Vector3};
use serde::{Deserialize, Serialize};

use crate::error::{MdError, Result};

/// A simulation cell: three lattice vectors plus a periodicity flag.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SimBox {
    /// Lattice matrix `h` — column `i` is lattice vector `i`. For a
    /// non-periodic box this is the identity and `periodic` is false.
    h: Matrix3<f64>,
    /// Cached inverse of `h` (fractional-coordinate transform).
    h_inv: Matrix3<f64>,
    /// Whether the box wraps. A non-periodic box behaves as free
    /// space.
    periodic: bool,
}

impl SimBox {
    /// A non-periodic (open) box. Minimum-image is the identity;
    /// [`volume`](Self::volume) is infinite.
    pub fn none() -> Self {
        SimBox {
            h: Matrix3::identity(),
            h_inv: Matrix3::identity(),
            periodic: false,
        }
    }

    /// An orthorhombic (rectangular) periodic box with the given edge
    /// lengths (nm).
    ///
    /// # Errors
    /// [`MdError::Invalid`] if any edge is not strictly positive and
    /// finite.
    pub fn orthorhombic(lx: f64, ly: f64, lz: f64) -> Result<Self> {
        for (name, v) in [("lx", lx), ("ly", ly), ("lz", lz)] {
            if !(v.is_finite() && v > 0.0) {
                return Err(MdError::invalid(
                    "box",
                    format!("edge {name} must be finite and positive, got {v}"),
                ));
            }
        }
        let h = Matrix3::from_diagonal(&Vector3::new(lx, ly, lz));
        Ok(SimBox {
            h,
            h_inv: Matrix3::from_diagonal(&Vector3::new(1.0 / lx, 1.0 / ly, 1.0 / lz)),
            periodic: true,
        })
    }

    /// A cubic periodic box of edge `l` (nm).
    ///
    /// # Errors
    /// [`MdError::Invalid`] if `l` is not strictly positive.
    pub fn cubic(l: f64) -> Result<Self> {
        Self::orthorhombic(l, l, l)
    }

    /// A general triclinic box from three lattice vectors.
    ///
    /// The vectors become the columns of `h`. They must be linearly
    /// independent (non-zero determinant); a positive determinant is
    /// enforced so the cell is right-handed.
    ///
    /// # Errors
    /// [`MdError::Invalid`] if the three vectors are coplanar /
    /// degenerate.
    pub fn triclinic(a: Vector3<f64>, b: Vector3<f64>, c: Vector3<f64>) -> Result<Self> {
        let h = Matrix3::from_columns(&[a, b, c]);
        let det = h.determinant();
        if !det.is_finite() || det.abs() < 1e-12 {
            return Err(MdError::invalid(
                "box",
                format!("lattice vectors are degenerate (det = {det})"),
            ));
        }
        if det < 0.0 {
            return Err(MdError::invalid(
                "box",
                "lattice vectors form a left-handed cell; swap two vectors",
            ));
        }
        let h_inv = h
            .try_inverse()
            .ok_or_else(|| MdError::invalid("box", "lattice matrix is not invertible"))?;
        Ok(SimBox {
            h,
            h_inv,
            periodic: true,
        })
    }

    /// Whether the box is periodic.
    pub fn is_periodic(&self) -> bool {
        self.periodic
    }

    /// The lattice matrix `h` (columns are the lattice vectors).
    pub fn matrix(&self) -> &Matrix3<f64> {
        &self.h
    }

    /// Cell volume (nm³). [`f64::INFINITY`] for a non-periodic box.
    pub fn volume(&self) -> f64 {
        if self.periodic {
            self.h.determinant().abs()
        } else {
            f64::INFINITY
        }
    }

    /// The three edge lengths (norms of the lattice vectors), nm.
    pub fn edge_lengths(&self) -> [f64; 3] {
        [
            self.h.column(0).norm(),
            self.h.column(1).norm(),
            self.h.column(2).norm(),
        ]
    }

    /// `true` if `h` is diagonal to a tight tolerance — i.e. the box
    /// is orthorhombic and the fast minimum-image path applies.
    pub fn is_orthorhombic(&self) -> bool {
        let h = &self.h;
        let off = h.m12.abs() + h.m13.abs() + h.m21.abs() + h.m23.abs() + h.m31.abs() + h.m32.abs();
        off < 1e-9
    }

    /// Returns the minimum-image displacement equivalent to `d`.
    ///
    /// For a non-periodic box this is `d` unchanged. For an
    /// orthorhombic box it is the component-wise nearest-image wrap.
    /// For a triclinic box it is a fractional-coordinate round
    /// followed by a 3×3×3 image search — exact when the cutoff is
    /// below half the minimum box width.
    pub fn min_image(&self, d: Vector3<f64>) -> Vector3<f64> {
        if !self.periodic {
            return d;
        }
        if self.is_orthorhombic() {
            let lx = self.h.m11;
            let ly = self.h.m22;
            let lz = self.h.m33;
            return Vector3::new(
                d.x - lx * (d.x / lx).round(),
                d.y - ly * (d.y / ly).round(),
                d.z - lz * (d.z / lz).round(),
            );
        }
        // Triclinic: round in fractional space, then refine over the
        // 27 nearest images to catch the cases where the rounded
        // fractional vector is not quite the geometric minimum.
        let frac = self.h_inv * d;
        let base = frac - frac.map(|x| x.round());
        let mut best = self.h * base;
        let mut best_sq = best.norm_squared();
        for i in -1..=1 {
            for j in -1..=1 {
                for k in -1..=1 {
                    if i == 0 && j == 0 && k == 0 {
                        continue;
                    }
                    let shifted = base + Vector3::new(i as f64, j as f64, k as f64);
                    let cart = self.h * shifted;
                    let sq = cart.norm_squared();
                    if sq < best_sq {
                        best_sq = sq;
                        best = cart;
                    }
                }
            }
        }
        best
    }

    /// The minimum-image distance between two points.
    pub fn distance(&self, a: Vector3<f64>, b: Vector3<f64>) -> f64 {
        self.min_image(b - a).norm()
    }

    /// Wraps a position into the primary cell `[0, 1)³` in fractional
    /// coordinates and returns it back in Cartesian space. A
    /// non-periodic box returns `p` unchanged.
    pub fn wrap(&self, p: Vector3<f64>) -> Vector3<f64> {
        if !self.periodic {
            return p;
        }
        let frac = self.h_inv * p;
        let wrapped = frac.map(|x| x - x.floor());
        self.h * wrapped
    }

    /// Largest cutoff radius for which the minimum-image convention is
    /// unambiguous: half the smallest perpendicular box width.
    ///
    /// For a triclinic cell the perpendicular widths are
    /// `volume / |area of the opposite face|`.
    pub fn max_cutoff(&self) -> f64 {
        if !self.periodic {
            return f64::INFINITY;
        }
        let a = self.h.column(0).into_owned();
        let b = self.h.column(1).into_owned();
        let c = self.h.column(2).into_owned();
        let vol = self.volume();
        let wa = vol / b.cross(&c).norm();
        let wb = vol / c.cross(&a).norm();
        let wc = vol / a.cross(&b).norm();
        0.5 * wa.min(wb).min(wc)
    }

    /// Uniformly scales the box by `factor` (used by the barostats).
    ///
    /// # Errors
    /// [`MdError::Invalid`] if `factor` is not finite and positive.
    pub fn scaled(&self, factor: f64) -> Result<Self> {
        if !(factor.is_finite() && factor > 0.0) {
            return Err(MdError::invalid(
                "box-scale",
                format!("factor must be finite and positive, got {factor}"),
            ));
        }
        if !self.periodic {
            return Ok(self.clone());
        }
        let h = self.h * factor;
        Ok(SimBox {
            h,
            h_inv: self.h_inv / factor,
            periodic: true,
        })
    }

    /// Anisotropically scales the box by a per-axis factor (used by
    /// the Parrinello-Rahman barostat).
    ///
    /// # Errors
    /// [`MdError::Invalid`] if any factor is not finite and positive.
    pub fn scaled_aniso(&self, factors: Vector3<f64>) -> Result<Self> {
        for (axis, f) in ["x", "y", "z"].iter().zip(factors.iter()) {
            if !(f.is_finite() && *f > 0.0) {
                return Err(MdError::invalid(
                    "box-scale",
                    format!("factor for {axis} must be finite and positive, got {f}"),
                ));
            }
        }
        if !self.periodic {
            return Ok(self.clone());
        }
        let scale = Matrix3::from_diagonal(&factors);
        let h = scale * self.h;
        let h_inv = h
            .try_inverse()
            .ok_or_else(|| MdError::invalid("box-scale", "scaled box became singular"))?;
        Ok(SimBox {
            h,
            h_inv,
            periodic: true,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn orthorhombic_volume_and_edges() {
        let b = SimBox::orthorhombic(2.0, 3.0, 4.0).unwrap();
        assert!((b.volume() - 24.0).abs() < 1e-9);
        assert_eq!(b.edge_lengths(), [2.0, 3.0, 4.0]);
        assert!(b.is_orthorhombic());
    }

    #[test]
    fn rejects_degenerate_boxes() {
        assert!(SimBox::orthorhombic(0.0, 1.0, 1.0).is_err());
        assert!(SimBox::orthorhombic(1.0, -1.0, 1.0).is_err());
        assert!(SimBox::cubic(f64::NAN).is_err());
    }

    #[test]
    fn min_image_orthorhombic_wraps() {
        let b = SimBox::cubic(10.0).unwrap();
        // A 9 nm displacement in a 10 nm box should fold to -1 nm.
        let d = b.min_image(Vector3::new(9.0, 0.0, 0.0));
        assert!((d.x - (-1.0)).abs() < 1e-9, "got {}", d.x);
        // Already minimal: unchanged.
        let d = b.min_image(Vector3::new(2.0, -3.0, 1.0));
        assert!((d - Vector3::new(2.0, -3.0, 1.0)).norm() < 1e-9);
    }

    #[test]
    fn min_image_distance_is_symmetric() {
        let b = SimBox::cubic(5.0).unwrap();
        let p = Vector3::new(0.5, 0.5, 0.5);
        let q = Vector3::new(4.7, 0.5, 0.5);
        // Across the boundary, the short way is 0.8 nm not 4.2 nm.
        assert!((b.distance(p, q) - 0.8).abs() < 1e-9);
        assert!((b.distance(p, q) - b.distance(q, p)).abs() < 1e-12);
    }

    #[test]
    fn triclinic_volume_matches_determinant() {
        // A simple sheared cell.
        let a = Vector3::new(4.0, 0.0, 0.0);
        let b = Vector3::new(1.0, 4.0, 0.0);
        let c = Vector3::new(1.0, 1.0, 4.0);
        let bx = SimBox::triclinic(a, b, c).unwrap();
        assert!((bx.volume() - 64.0).abs() < 1e-9);
        assert!(!bx.is_orthorhombic());
    }

    #[test]
    fn triclinic_min_image_is_shortest() {
        let a = Vector3::new(4.0, 0.0, 0.0);
        let b = Vector3::new(1.0, 4.0, 0.0);
        let c = Vector3::new(0.5, 0.5, 4.0);
        let bx = SimBox::triclinic(a, b, c).unwrap();
        // Pick a displacement near a lattice vector; the min image
        // must be at least as short as any 27-image candidate.
        let d = Vector3::new(3.6, 3.7, 0.2);
        let mi = bx.min_image(d);
        let h = bx.matrix();
        for i in -1..=1 {
            for j in -1..=1 {
                for k in -1..=1 {
                    let cand = d - h * Vector3::new(i as f64, j as f64, k as f64);
                    assert!(mi.norm() <= cand.norm() + 1e-9);
                }
            }
        }
    }

    #[test]
    fn wrap_lands_in_primary_cell() {
        let b = SimBox::cubic(10.0).unwrap();
        let w = b.wrap(Vector3::new(23.0, -4.0, 10.0));
        for c in [w.x, w.y, w.z] {
            assert!((0.0..10.0).contains(&c) || c.abs() < 1e-9, "c = {c}");
        }
    }

    #[test]
    fn non_periodic_is_identity() {
        let b = SimBox::none();
        assert!(!b.is_periodic());
        assert_eq!(b.volume(), f64::INFINITY);
        let d = Vector3::new(100.0, 200.0, 300.0);
        assert_eq!(b.min_image(d), d);
        assert_eq!(b.wrap(d), d);
    }

    #[test]
    fn scaling_changes_volume_cubically() {
        let b = SimBox::cubic(2.0).unwrap();
        let s = b.scaled(1.5).unwrap();
        assert!((s.volume() - 8.0 * 1.5_f64.powi(3)).abs() < 1e-9);
    }

    #[test]
    fn max_cutoff_is_half_min_width() {
        let b = SimBox::orthorhombic(4.0, 6.0, 8.0).unwrap();
        assert!((b.max_cutoff() - 2.0).abs() < 1e-9);
    }
}
