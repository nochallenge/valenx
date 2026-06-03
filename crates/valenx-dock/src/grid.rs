//! 3-D scalar grid + trilinear interpolation.
//!
//! Vina precomputes one grid per ligand atom type: the value at each
//! grid point is the energy that a single atom of that type would
//! feel at that position due to the entire receptor. Pose evaluation
//! is then a per-atom trilinear lookup.

use nalgebra::Vector3;

use crate::atom_type::Ad4AtomType;
use crate::receptor::Receptor;
use crate::score::{pair_score, surface_distance};
use crate::PAIR_CUTOFF;

/// 3-D scalar grid in receptor coordinates.
#[derive(Clone, Debug)]
pub struct Grid3D {
    /// World-space coordinate of the grid origin (grid index (0,0,0)).
    pub origin: Vector3<f64>,
    /// Spacing in Å along each axis (Vina default: 0.375 Å).
    pub spacing: f64,
    /// Grid dimensions (nx, ny, nz).
    pub dims: (usize, usize, usize),
    /// Row-major data: data[ix + iy*nx + iz*nx*ny].
    pub data: Vec<f64>,
}

impl Grid3D {
    /// Build an empty grid (all zeros) with the given geometry.
    pub fn zeros(origin: Vector3<f64>, spacing: f64, dims: (usize, usize, usize)) -> Self {
        let total = dims.0 * dims.1 * dims.2;
        Self {
            origin,
            spacing,
            dims,
            data: vec![0.0; total],
        }
    }

    /// Linear index for (ix, iy, iz).
    pub fn index(&self, ix: usize, iy: usize, iz: usize) -> usize {
        ix + iy * self.dims.0 + iz * self.dims.0 * self.dims.1
    }

    /// World-space position of grid point (ix, iy, iz).
    pub fn position(&self, ix: usize, iy: usize, iz: usize) -> Vector3<f64> {
        self.origin + Vector3::new(ix as f64, iy as f64, iz as f64) * self.spacing
    }

    /// Trilinear interpolation. Returns 0.0 outside the grid (Vina
    /// treats out-of-box as zero potential; the search rarely strays).
    ///
    /// Defensive: also returns 0.0 if `spacing` is non-positive or
    /// non-finite — a degenerate grid behaves as if everything is
    /// out-of-box. Real callers should validate via
    /// [`crate::DockConfig::validate`] before constructing a grid;
    /// this guard exists so a logic error in one of those callers
    /// cannot produce NaN-poisoned scores that propagate silently.
    pub fn sample(&self, p: Vector3<f64>) -> f64 {
        if !self.spacing.is_finite() || self.spacing <= 0.0 {
            return 0.0;
        }
        let local = (p - self.origin) / self.spacing;
        let fx = local.x.floor();
        let fy = local.y.floor();
        let fz = local.z.floor();
        let ix = fx as isize;
        let iy = fy as isize;
        let iz = fz as isize;
        if ix < 0
            || iy < 0
            || iz < 0
            || ix + 1 >= self.dims.0 as isize
            || iy + 1 >= self.dims.1 as isize
            || iz + 1 >= self.dims.2 as isize
        {
            return 0.0;
        }
        let dx = local.x - fx;
        let dy = local.y - fy;
        let dz = local.z - fz;
        let ix = ix as usize;
        let iy = iy as usize;
        let iz = iz as usize;
        let c000 = self.data[self.index(ix, iy, iz)];
        let c100 = self.data[self.index(ix + 1, iy, iz)];
        let c010 = self.data[self.index(ix, iy + 1, iz)];
        let c001 = self.data[self.index(ix, iy, iz + 1)];
        let c110 = self.data[self.index(ix + 1, iy + 1, iz)];
        let c101 = self.data[self.index(ix + 1, iy, iz + 1)];
        let c011 = self.data[self.index(ix, iy + 1, iz + 1)];
        let c111 = self.data[self.index(ix + 1, iy + 1, iz + 1)];

        let c00 = c000 * (1.0 - dx) + c100 * dx;
        let c01 = c001 * (1.0 - dx) + c101 * dx;
        let c10 = c010 * (1.0 - dx) + c110 * dx;
        let c11 = c011 * (1.0 - dx) + c111 * dx;
        let c0 = c00 * (1.0 - dy) + c10 * dy;
        let c1 = c01 * (1.0 - dy) + c11 * dy;
        c0 * (1.0 - dz) + c1 * dz
    }
}

use std::collections::HashMap;

use rayon::prelude::*;

use crate::ligand::Ligand;

/// One grid per ligand atom type — the per-pose evaluation reads
/// from these.
#[derive(Clone, Debug)]
pub struct GridBundle {
    /// Keyed by ligand atom type.
    pub grids: HashMap<Ad4AtomType, Grid3D>,
}

impl GridBundle {
    /// Build a grid for every distinct atom type in the ligand.
    /// Grids are independent (each is a pure function of the receptor
    /// and probe atom type), so we precompute them in parallel via
    /// rayon — a typical 6-atom-type bundle drops from ~6× one
    /// grid-time to roughly one grid-time on a multicore host.
    pub fn build(
        receptor: &Receptor,
        ligand: &Ligand,
        origin: Vector3<f64>,
        spacing: f64,
        dims: (usize, usize, usize),
    ) -> Self {
        let mut types: std::collections::BTreeSet<Ad4AtomType> = Default::default();
        for a in &ligand.atoms {
            types.insert(a.ad4_type);
        }
        let grids: HashMap<Ad4AtomType, Grid3D> = types
            .into_par_iter()
            .map(|t| {
                (
                    t,
                    precompute_receptor_grid(receptor, t, origin, spacing, dims),
                )
            })
            .collect();
        Self { grids }
    }
}

/// Build a grid of the inter-molecular energy a single probe atom of
/// type `probe` would feel at each grid point due to all receptor
/// atoms within the Vina cutoff.
///
/// Defensive: if `spacing` is non-positive or non-finite the function
/// short-circuits to an all-zeros grid (the same behaviour
/// [`Grid3D::sample`] degrades to). Callers should validate via
/// [`crate::DockConfig::validate`] first; this guard is a
/// belt-and-suspenders for the case where the validate call is
/// skipped by a faulty caller.
pub fn precompute_receptor_grid(
    receptor: &Receptor,
    probe: Ad4AtomType,
    origin: Vector3<f64>,
    spacing: f64,
    dims: (usize, usize, usize),
) -> Grid3D {
    let mut grid = Grid3D::zeros(origin, spacing, dims);
    if !spacing.is_finite() || spacing <= 0.0 {
        return grid;
    }
    let probe_vdw = probe.props().vdw_radius;
    let cutoff_sq = PAIR_CUTOFF * PAIR_CUTOFF;
    for iz in 0..dims.2 {
        for iy in 0..dims.1 {
            for ix in 0..dims.0 {
                let p = grid.position(ix, iy, iz);
                let mut sum = 0.0;
                for ra in &receptor.atoms {
                    let r2 = (p - ra.position).norm_squared();
                    if r2 > cutoff_sq {
                        continue;
                    }
                    let d =
                        surface_distance(p, ra.position, probe_vdw, ra.ad4_type.props().vdw_radius);
                    sum += pair_score(probe, ra.ad4_type, d);
                }
                let idx = grid.index(ix, iy, iz);
                grid.data[idx] = sum;
            }
        }
    }
    grid
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empty_grid_sample_is_zero_at_origin() {
        let g = Grid3D::zeros(Vector3::zeros(), 1.0, (3, 3, 3));
        assert_eq!(g.sample(Vector3::new(1.0, 1.0, 1.0)), 0.0);
    }

    #[test]
    fn linear_field_recovered_exactly_by_trilinear() {
        // Fill grid so that data(x,y,z) = x (a linear function).
        // Trilinear interpolation should reproduce it exactly.
        let mut g = Grid3D::zeros(Vector3::zeros(), 1.0, (4, 4, 4));
        for iz in 0..4 {
            for iy in 0..4 {
                for ix in 0..4 {
                    let idx = g.index(ix, iy, iz);
                    g.data[idx] = ix as f64;
                }
            }
        }
        // Sample at (1.5, 1.5, 1.5) — expect 1.5.
        let v = g.sample(Vector3::new(1.5, 1.5, 1.5));
        assert!((v - 1.5).abs() < 1e-12);
        // Sample at (2.25, 0.5, 1.0) — expect 2.25.
        let v = g.sample(Vector3::new(2.25, 0.5, 1.0));
        assert!((v - 2.25).abs() < 1e-12);
    }

    #[test]
    fn out_of_box_returns_zero() {
        let g = Grid3D::zeros(Vector3::zeros(), 1.0, (3, 3, 3));
        assert_eq!(g.sample(Vector3::new(-1.0, 0.0, 0.0)), 0.0);
        assert_eq!(g.sample(Vector3::new(10.0, 0.0, 0.0)), 0.0);
    }

    use crate::atom_type::Ad4AtomType;
    use crate::receptor::{Receptor, ReceptorAtom};

    #[test]
    fn precompute_grid_has_negative_minimum_for_attractive_pair() {
        // Single C atom at origin. Grid covers (-3..3) at 1 Å spacing.
        // Probe type = C (hydrophobic). At grid points near the
        // receptor atom the score should be NEGATIVE (attractive).
        let receptor = Receptor {
            atoms: vec![ReceptorAtom {
                position: Vector3::zeros(),
                ad4_type: Ad4AtomType::C,
                partial_charge: 0.0,
            }],
        };
        let g = precompute_receptor_grid(
            &receptor,
            Ad4AtomType::C,
            Vector3::new(-3.0, -3.0, -3.0),
            1.0,
            (7, 7, 7),
        );
        let minv = g.data.iter().cloned().fold(f64::INFINITY, f64::min);
        assert!(minv < 0.0, "expected attractive min, got {minv}");
    }

    use crate::ligand::Ligand;

    #[test]
    fn grid_bundle_covers_every_ligand_type() {
        // Ligand has C and OA — bundle should have two grids.
        let pdbqt = "\
ROOT
ATOM      1  C1  LIG A   1       0.000   0.000   0.000  1.00  0.00     0.000 C 
ATOM      2  O1  LIG A   1       1.500   0.000   0.000  1.00  0.00     0.000 OA
ENDROOT
TORSDOF 0
";
        let lig = Ligand::from_pdbqt(pdbqt).unwrap();
        let receptor = Receptor {
            atoms: vec![ReceptorAtom {
                position: Vector3::new(10.0, 0.0, 0.0),
                ad4_type: Ad4AtomType::C,
                partial_charge: 0.0,
            }],
        };
        let bundle = GridBundle::build(
            &receptor,
            &lig,
            Vector3::new(-5.0, -5.0, -5.0),
            1.0,
            (11, 11, 11),
        );
        assert_eq!(bundle.grids.len(), 2);
        assert!(bundle.grids.contains_key(&Ad4AtomType::C));
        assert!(bundle.grids.contains_key(&Ad4AtomType::OA));
    }
}
