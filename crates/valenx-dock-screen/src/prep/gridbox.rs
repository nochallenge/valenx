//! Feature 1 — receptor grid-box definition.
//!
//! A [`GridBox`] is the rectangular search volume a docking run
//! explores: a centre, three edge lengths and a grid spacing. It is
//! the docking analogue of AutoGrid's `npts` / `gridcenter` /
//! `spacing` triple. The grid maps in [`crate::score`] are sampled on
//! the lattice this box defines.

use nalgebra::Vector3;

use crate::error::{DockScreenError, Result};

/// A rectangular docking search box: centre, edge lengths and the
/// lattice spacing of the affinity grids sampled inside it.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GridBox {
    /// Centre of the box in receptor coordinates (Å).
    pub center: Vector3<f64>,
    /// Edge lengths along x / y / z (Å). All strictly positive.
    pub size: Vector3<f64>,
    /// Lattice spacing of the affinity grids (Å). AutoGrid's default
    /// is `0.375`; Vina uses the same value.
    pub spacing: f64,
}

impl GridBox {
    /// AutoGrid / Vina default grid spacing in Å.
    pub const DEFAULT_SPACING: f64 = 0.375;

    /// Build a grid box from a centre and edge lengths, using the
    /// default spacing. Returns [`DockScreenError::Invalid`] if any
    /// edge is non-positive or non-finite.
    pub fn new(center: [f64; 3], size: [f64; 3]) -> Result<Self> {
        Self::with_spacing(center, size, Self::DEFAULT_SPACING)
    }

    /// Build a grid box with an explicit spacing.
    pub fn with_spacing(center: [f64; 3], size: [f64; 3], spacing: f64) -> Result<Self> {
        for (axis, v) in ["x", "y", "z"].iter().zip(size.iter()) {
            if !v.is_finite() || *v <= 0.0 {
                return Err(DockScreenError::invalid(
                    "size",
                    format!("box edge `{axis}` must be positive and finite, got {v}"),
                ));
            }
        }
        for (axis, v) in ["x", "y", "z"].iter().zip(center.iter()) {
            if !v.is_finite() {
                return Err(DockScreenError::invalid(
                    "center",
                    format!("box centre `{axis}` must be finite, got {v}"),
                ));
            }
        }
        if !spacing.is_finite() || spacing <= 0.0 {
            return Err(DockScreenError::invalid(
                "spacing",
                format!("grid spacing must be positive and finite, got {spacing}"),
            ));
        }
        Ok(GridBox {
            center: Vector3::new(center[0], center[1], center[2]),
            size: Vector3::new(size[0], size[1], size[2]),
            spacing,
        })
    }

    /// A cubic box of the given edge length centred on `center`.
    pub fn cubic(center: [f64; 3], edge: f64) -> Result<Self> {
        Self::new(center, [edge, edge, edge])
    }

    /// The corner of the box with the smallest x / y / z — the grid
    /// origin (`center - size/2`).
    pub fn origin(&self) -> Vector3<f64> {
        self.center - self.size / 2.0
    }

    /// The corner of the box with the largest x / y / z.
    pub fn max_corner(&self) -> Vector3<f64> {
        self.center + self.size / 2.0
    }

    /// Lattice dimensions `(nx, ny, nz)` — the number of grid points
    /// along each axis. Matches [`valenx_dock::DockConfig::grid_dims`]:
    /// `ceil(size / spacing) + 1`, with a floor of 2 points per axis.
    pub fn dims(&self) -> (usize, usize, usize) {
        let n = |edge: f64| ((edge / self.spacing).ceil() as usize + 1).max(2);
        (n(self.size.x), n(self.size.y), n(self.size.z))
    }

    /// Total number of grid points (`nx * ny * nz`).
    pub fn point_count(&self) -> usize {
        let (nx, ny, nz) = self.dims();
        nx * ny * nz
    }

    /// Box volume in Å³.
    pub fn volume(&self) -> f64 {
        self.size.x * self.size.y * self.size.z
    }

    /// `true` if `p` lies inside the box (boundary inclusive).
    pub fn contains(&self, p: Vector3<f64>) -> bool {
        let lo = self.origin();
        let hi = self.max_corner();
        p.x >= lo.x && p.x <= hi.x && p.y >= lo.y && p.y <= hi.y && p.z >= lo.z && p.z <= hi.z
    }

    /// Build a box that snugly encloses a set of atom positions plus a
    /// `padding` margin on every side. The classic "blind docking"
    /// setup — derive the search volume from the receptor / a known
    /// ligand rather than guessing coordinates.
    ///
    /// Returns [`DockScreenError::Invalid`] if `points` is empty.
    pub fn enclosing(points: &[Vector3<f64>], padding: f64) -> Result<Self> {
        if points.is_empty() {
            return Err(DockScreenError::invalid(
                "points",
                "cannot build an enclosing box from zero points",
            ));
        }
        let pad = padding.max(0.0);
        let mut lo = points[0];
        let mut hi = points[0];
        for p in points {
            lo = lo.inf(p);
            hi = hi.sup(p);
        }
        let center = (lo + hi) / 2.0;
        // Each edge is at least 2*spacing wide so the box never
        // degenerates for a single-point input.
        let size = (hi - lo).map(|e| (e + 2.0 * pad).max(2.0 * Self::DEFAULT_SPACING));
        Self::new([center.x, center.y, center.z], [size.x, size.y, size.z])
    }

    /// Bridge to the `valenx-dock` config — fills the centre, size and
    /// spacing of a [`valenx_dock::DockConfig`], leaving the search
    /// knobs at their defaults.
    pub fn to_dock_config(&self) -> valenx_dock::DockConfig {
        valenx_dock::DockConfig {
            center: self.center,
            size: self.size,
            grid_spacing: self.spacing,
            ..valenx_dock::DockConfig::default()
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn new_rejects_bad_edges() {
        assert!(GridBox::new([0.0, 0.0, 0.0], [0.0, 20.0, 20.0]).is_err());
        assert!(GridBox::new([0.0, 0.0, 0.0], [20.0, -1.0, 20.0]).is_err());
        assert!(GridBox::new([0.0, 0.0, 0.0], [20.0, f64::NAN, 20.0]).is_err());
        assert!(GridBox::new([f64::INFINITY, 0.0, 0.0], [20.0, 20.0, 20.0]).is_err());
    }

    #[test]
    fn new_rejects_bad_spacing() {
        assert!(GridBox::with_spacing([0.0; 3], [20.0; 3], 0.0).is_err());
        assert!(GridBox::with_spacing([0.0; 3], [20.0; 3], -0.1).is_err());
    }

    #[test]
    fn origin_and_max_corner_straddle_center() {
        let gb = GridBox::cubic([10.0, 20.0, 30.0], 20.0).unwrap();
        assert_eq!(gb.origin(), Vector3::new(0.0, 10.0, 20.0));
        assert_eq!(gb.max_corner(), Vector3::new(20.0, 30.0, 40.0));
    }

    #[test]
    fn dims_match_dock_convention() {
        // 7.5 Å / 0.375 = 20, + 1 = 21 — same as DockConfig::grid_dims.
        let gb = GridBox::with_spacing([0.0; 3], [7.5, 7.5, 7.5], 0.375).unwrap();
        assert_eq!(gb.dims(), (21, 21, 21));
        assert_eq!(gb.point_count(), 21 * 21 * 21);
    }

    #[test]
    fn contains_is_boundary_inclusive() {
        let gb = GridBox::cubic([0.0, 0.0, 0.0], 10.0).unwrap();
        assert!(gb.contains(Vector3::new(0.0, 0.0, 0.0)));
        assert!(gb.contains(Vector3::new(5.0, 5.0, 5.0))); // exact corner
        assert!(!gb.contains(Vector3::new(5.1, 0.0, 0.0)));
    }

    #[test]
    fn enclosing_box_covers_points_with_padding() {
        let pts = vec![Vector3::new(0.0, 0.0, 0.0), Vector3::new(4.0, 6.0, 8.0)];
        let gb = GridBox::enclosing(&pts, 2.0).unwrap();
        // Centre is the midpoint.
        assert!((gb.center - Vector3::new(2.0, 3.0, 4.0)).norm() < 1e-9);
        // Every point + a 2 Å pad must sit inside.
        for p in &pts {
            assert!(gb.contains(*p), "point {p:?} fell outside enclosing box");
        }
        // The extent should be the span (4,6,8) plus 2*pad on each.
        assert!((gb.size.x - 8.0).abs() < 1e-9);
        assert!((gb.size.y - 10.0).abs() < 1e-9);
    }

    #[test]
    fn enclosing_rejects_empty() {
        assert!(GridBox::enclosing(&[], 2.0).is_err());
    }

    #[test]
    fn enclosing_single_point_does_not_degenerate() {
        let gb = GridBox::enclosing(&[Vector3::new(1.0, 2.0, 3.0)], 0.0).unwrap();
        // Even with zero padding the box keeps a minimum positive edge.
        assert!(gb.size.x > 0.0 && gb.size.y > 0.0 && gb.size.z > 0.0);
    }

    #[test]
    fn to_dock_config_carries_geometry() {
        let gb = GridBox::with_spacing([1.0, 2.0, 3.0], [10.0, 11.0, 12.0], 0.5).unwrap();
        let cfg = gb.to_dock_config();
        assert_eq!(cfg.center, Vector3::new(1.0, 2.0, 3.0));
        assert_eq!(cfg.size, Vector3::new(10.0, 11.0, 12.0));
        assert_eq!(cfg.grid_spacing, 0.5);
    }

    #[test]
    fn volume_is_product_of_edges() {
        let gb = GridBox::new([0.0; 3], [2.0, 3.0, 4.0]).unwrap();
        assert!((gb.volume() - 24.0).abs() < 1e-9);
    }
}
