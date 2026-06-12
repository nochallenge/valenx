//! Features 7 & 8 — affinity-map precomputation and trilinear scoring.
//!
//! Re-evaluating every receptor/ligand atom pair for every pose is the
//! dominant cost of a docking search. The classical fix, used by both
//! AutoGrid and Vina, is to *precompute* the receptor's contribution
//! once: for each ligand atom type, build a 3-D grid whose value at
//! each lattice point is the energy a single probe atom of that type
//! would feel there. Pose scoring then becomes one trilinear lookup
//! per ligand atom — independent of receptor size.
//!
//! - **Feature 7** — [`AffinityMapSet::precompute`] builds one
//!   [`AffinityMap`] per distinct ligand atom type over a
//!   [`GridBox`]. Maps are independent (each is a pure function of the
//!   receptor and the probe type), so they are computed in parallel.
//! - **Feature 8** — [`AffinityMap::interpolate`] does the trilinear
//!   lookup, and [`score_ligand_on_maps`] sums it across a posed
//!   ligand for a fast grid-based score.
//!
//! Both the Vina-class and the AutoDock4-class scoring functions can
//! drive the precompute (see [`MapKind`]). The Vina maps reuse
//! [`valenx_dock`]'s `Grid3D` precompute directly; the AutoDock4 maps
//! sample [`crate::score::ad4`] — the electrostatic component is kept
//! as its own *charge-independent* grid (energy per unit ligand
//! charge) so a single grid serves ligand atoms of any charge.

use std::collections::BTreeSet;

use nalgebra::Vector3;
use rayon::prelude::*;

use valenx_dock::atom_type::Ad4AtomType;
use valenx_dock::receptor::Receptor;

use crate::error::{DockScreenError, Result};
use crate::prep::gridbox::GridBox;

/// Which scoring function an affinity map family was sampled from.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum MapKind {
    /// Maps sampled from the Vina-class scoring function.
    Vina,
    /// Maps sampled from the AutoDock4-class force field. The
    /// electrostatic contribution is stored separately (see
    /// [`AffinityMapSet::electrostatic`]).
    AutoDock4,
}

/// A single 3-D scalar affinity map: the energy a probe atom of one
/// type feels at each lattice point of a [`GridBox`].
#[derive(Clone, Debug, PartialEq)]
pub struct AffinityMap {
    /// World-space coordinate of lattice point `(0,0,0)`.
    pub origin: Vector3<f64>,
    /// Lattice spacing (Å).
    pub spacing: f64,
    /// Lattice dimensions `(nx, ny, nz)`.
    pub dims: (usize, usize, usize),
    /// Row-major scalar data: `data[ix + iy*nx + iz*nx*ny]`.
    pub data: Vec<f64>,
}

impl AffinityMap {
    /// An all-zeros map with the geometry of `grid`.
    pub fn zeros(grid: &GridBox) -> Self {
        let dims = grid.dims();
        AffinityMap {
            origin: grid.origin(),
            spacing: grid.spacing,
            dims,
            data: vec![0.0; dims.0 * dims.1 * dims.2],
        }
    }

    /// Linear index of lattice point `(ix, iy, iz)`.
    pub fn index(&self, ix: usize, iy: usize, iz: usize) -> usize {
        ix + iy * self.dims.0 + iz * self.dims.0 * self.dims.1
    }

    /// World-space position of lattice point `(ix, iy, iz)`.
    pub fn position(&self, ix: usize, iy: usize, iz: usize) -> Vector3<f64> {
        self.origin + Vector3::new(ix as f64, iy as f64, iz as f64) * self.spacing
    }

    /// Trilinear interpolation of the map at world position `p`.
    /// Returns `0.0` for points outside the lattice — the same
    /// out-of-box convention AutoGrid and Vina use.
    pub fn interpolate(&self, p: Vector3<f64>) -> f64 {
        if !self.spacing.is_finite() || self.spacing <= 0.0 {
            return 0.0;
        }
        let local = (p - self.origin) / self.spacing;
        let fx = local.x.floor();
        let fy = local.y.floor();
        let fz = local.z.floor();
        let (ix, iy, iz) = (fx as isize, fy as isize, fz as isize);
        if ix < 0
            || iy < 0
            || iz < 0
            || ix + 1 >= self.dims.0 as isize
            || iy + 1 >= self.dims.1 as isize
            || iz + 1 >= self.dims.2 as isize
        {
            return 0.0;
        }
        let (dx, dy, dz) = (local.x - fx, local.y - fy, local.z - fz);
        let (ix, iy, iz) = (ix as usize, iy as usize, iz as usize);
        let c = |x: usize, y: usize, z: usize| self.data[self.index(x, y, z)];
        let c000 = c(ix, iy, iz);
        let c100 = c(ix + 1, iy, iz);
        let c010 = c(ix, iy + 1, iz);
        let c001 = c(ix, iy, iz + 1);
        let c110 = c(ix + 1, iy + 1, iz);
        let c101 = c(ix + 1, iy, iz + 1);
        let c011 = c(ix, iy + 1, iz + 1);
        let c111 = c(ix + 1, iy + 1, iz + 1);
        let c00 = c000 * (1.0 - dx) + c100 * dx;
        let c01 = c001 * (1.0 - dx) + c101 * dx;
        let c10 = c010 * (1.0 - dx) + c110 * dx;
        let c11 = c011 * (1.0 - dx) + c111 * dx;
        let c0 = c00 * (1.0 - dy) + c10 * dy;
        let c1 = c01 * (1.0 - dy) + c11 * dy;
        c0 * (1.0 - dz) + c1 * dz
    }

    /// The most negative (most favourable) value anywhere on the map.
    pub fn min_value(&self) -> f64 {
        self.data.iter().copied().fold(f64::INFINITY, f64::min)
    }
}

/// A family of affinity maps — one per ligand atom type — plus, for
/// the AutoDock4 family, a charge-independent electrostatic map.
#[derive(Clone, Debug)]
pub struct AffinityMapSet {
    /// Which scoring function the maps were sampled from.
    pub kind: MapKind,
    /// `(atom type, map)` pairs, sorted by atom type.
    pub maps: Vec<(Ad4AtomType, AffinityMap)>,
    /// AutoDock4 only: the electrostatic potential per unit ligand
    /// charge. `None` for [`MapKind::Vina`] (the Vina function folds
    /// electrostatics into its per-type maps).
    pub electrostatic: Option<AffinityMap>,
}

impl AffinityMapSet {
    /// The map for a given atom type, if present.
    pub fn map_for(&self, t: Ad4AtomType) -> Option<&AffinityMap> {
        self.maps.iter().find(|(mt, _)| *mt == t).map(|(_, m)| m)
    }

    /// Number of per-type maps in the set.
    pub fn len(&self) -> usize {
        self.maps.len()
    }

    /// `true` if the set has no per-type maps.
    pub fn is_empty(&self) -> bool {
        self.maps.is_empty()
    }

    /// Feature 7 — precompute an affinity-map family covering every
    /// distinct atom type in `ligand_types`, over `grid`, using
    /// scoring function `kind`.
    ///
    /// Returns [`DockScreenError::Invalid`] if `ligand_types` is empty
    /// or the receptor has no atoms.
    pub fn precompute(
        receptor: &Receptor,
        ligand_types: &[Ad4AtomType],
        grid: &GridBox,
        kind: MapKind,
    ) -> Result<Self> {
        if receptor.atoms.is_empty() {
            return Err(DockScreenError::invalid_receptor("receptor has no atoms"));
        }
        let types: BTreeSet<Ad4AtomType> = ligand_types.iter().copied().collect();
        if types.is_empty() {
            return Err(DockScreenError::invalid(
                "ligand_types",
                "cannot precompute affinity maps for zero atom types",
            ));
        }
        let maps: Vec<(Ad4AtomType, AffinityMap)> = types
            .par_iter()
            .map(|&t| (t, precompute_one(receptor, t, grid, kind)))
            .collect();
        let electrostatic = match kind {
            MapKind::Vina => None,
            MapKind::AutoDock4 => Some(precompute_electrostatic(receptor, grid)),
        };
        Ok(AffinityMapSet {
            kind,
            maps,
            electrostatic,
        })
    }
}

/// Precompute a single per-type affinity map.
fn precompute_one(
    receptor: &Receptor,
    probe: Ad4AtomType,
    grid: &GridBox,
    kind: MapKind,
) -> AffinityMap {
    match kind {
        MapKind::Vina => precompute_vina(receptor, probe, grid),
        MapKind::AutoDock4 => precompute_ad4(receptor, probe, grid),
    }
}

/// Vina per-type map — delegates to `valenx-dock`'s grid precompute so
/// the values are bit-identical to a native `valenx-dock` run.
fn precompute_vina(receptor: &Receptor, probe: Ad4AtomType, grid: &GridBox) -> AffinityMap {
    let g = valenx_dock::grid::precompute_receptor_grid(
        receptor,
        probe,
        grid.origin(),
        grid.spacing,
        grid.dims(),
    );
    AffinityMap {
        origin: g.origin,
        spacing: g.spacing,
        dims: g.dims,
        data: g.data,
    }
}

/// AutoDock4 per-type map — samples the AD4 vdW + H-bond + desolvation
/// terms (electrostatics is its own map). The probe is treated as
/// uncharged so the per-type map is charge-independent; the
/// electrostatic contribution is added at lookup time from the
/// separate map scaled by the ligand atom's charge.
fn precompute_ad4(receptor: &Receptor, probe: Ad4AtomType, grid: &GridBox) -> AffinityMap {
    let mut map = AffinityMap::zeros(grid);
    let (nx, ny, nz) = map.dims;
    for iz in 0..nz {
        for iy in 0..ny {
            for ix in 0..nx {
                let p = map.position(ix, iy, iz);
                // Score a single uncharged probe atom of `probe` type
                // at this lattice point against the whole receptor.
                let terms = crate::score::ad4::score_complex(receptor, &[(p, probe, 0.0)], 0);
                // Drop electrostatics (zero anyway — uncharged probe)
                // and the torsional term (zero, n_torsions=0).
                let idx = map.index(ix, iy, iz);
                map.data[idx] = terms.vdw + terms.hbond + terms.desolvation;
            }
        }
    }
    map
}

/// AutoDock4 electrostatic map — the screened-Coulomb potential a unit
/// positive charge would feel at each lattice point. Scaled by the
/// ligand atom's actual charge at lookup time.
fn precompute_electrostatic(receptor: &Receptor, grid: &GridBox) -> AffinityMap {
    let mut map = AffinityMap::zeros(grid);
    let (nx, ny, nz) = map.dims;
    for iz in 0..nz {
        for iy in 0..ny {
            for ix in 0..nx {
                let p = map.position(ix, iy, iz);
                // Probe with unit positive charge; use a tiny VDW
                // atom type (HD) so only the electrostatic term
                // matters here.
                let terms =
                    crate::score::ad4::score_complex(receptor, &[(p, Ad4AtomType::HD, 1.0)], 0);
                let idx = map.index(ix, iy, iz);
                map.data[idx] = terms.electrostatic;
            }
        }
    }
    map
}

/// Feature 8 — fast grid-based score of a posed ligand.
///
/// `ligand_atoms` is `(world position, AD4 type, partial charge)` per
/// atom. Each atom's per-type map is sampled by trilinear
/// interpolation; for an AutoDock4 map set the electrostatic map is
/// also sampled and scaled by the atom's charge. Atoms whose type has
/// no map (should not happen for a set built from this ligand) score
/// `0`.
pub fn score_ligand_on_maps(
    maps: &AffinityMapSet,
    ligand_atoms: &[(Vector3<f64>, Ad4AtomType, f64)],
) -> f64 {
    let mut sum = 0.0;
    for &(p, t, q) in ligand_atoms {
        if let Some(m) = maps.map_for(t) {
            sum += m.interpolate(p);
        }
        if let Some(elec) = &maps.electrostatic {
            sum += q * elec.interpolate(p);
        }
    }
    sum
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_dock::receptor::ReceptorAtom;

    fn carbon_receptor() -> Receptor {
        Receptor {
            atoms: vec![ReceptorAtom {
                position: Vector3::zeros(),
                ad4_type: Ad4AtomType::C,
                partial_charge: -0.2,
            }],
        }
    }

    #[test]
    fn precompute_rejects_empty_inputs() {
        let grid = GridBox::cubic([0.0; 3], 6.0).unwrap();
        let empty = Receptor::default();
        assert!(
            AffinityMapSet::precompute(&empty, &[Ad4AtomType::C], &grid, MapKind::Vina).is_err()
        );
        let r = carbon_receptor();
        assert!(AffinityMapSet::precompute(&r, &[], &grid, MapKind::Vina).is_err());
    }

    #[test]
    fn vina_map_set_covers_every_type() {
        let grid = GridBox::cubic([0.0; 3], 6.0).unwrap();
        let r = carbon_receptor();
        let set = AffinityMapSet::precompute(
            &r,
            &[Ad4AtomType::C, Ad4AtomType::OA, Ad4AtomType::C],
            &grid,
            MapKind::Vina,
        )
        .unwrap();
        // Duplicate C collapses → two distinct maps.
        assert_eq!(set.len(), 2);
        assert!(set.map_for(Ad4AtomType::C).is_some());
        assert!(set.map_for(Ad4AtomType::OA).is_some());
        // Vina folds electrostatics in — no separate elec map.
        assert!(set.electrostatic.is_none());
    }

    #[test]
    fn ad4_map_set_has_separate_electrostatic_map() {
        let grid = GridBox::cubic([0.0; 3], 6.0).unwrap();
        let r = carbon_receptor();
        let set =
            AffinityMapSet::precompute(&r, &[Ad4AtomType::C], &grid, MapKind::AutoDock4).unwrap();
        assert_eq!(set.kind, MapKind::AutoDock4);
        assert!(set.electrostatic.is_some());
    }

    #[test]
    fn vina_map_has_a_favourable_minimum() {
        // A carbon probe over a carbon receptor — somewhere on the
        // grid the energy is attractive (negative).
        let grid = GridBox::with_spacing([0.0; 3], [6.0; 3], 1.0).unwrap();
        let r = carbon_receptor();
        let set = AffinityMapSet::precompute(&r, &[Ad4AtomType::C], &grid, MapKind::Vina).unwrap();
        let m = set.map_for(Ad4AtomType::C).unwrap();
        assert!(m.min_value() < 0.0, "Vina map min = {}", m.min_value());
    }

    #[test]
    fn trilinear_recovers_a_linear_field_exactly() {
        // Build a map whose value equals the x lattice index, then
        // interpolate — trilinear must reproduce a linear field.
        // The box is *centred* at the origin, so its grid origin (the
        // min corner) is center - size/2 = (-1.5, -1.5, -1.5). To
        // sample lattice-local coordinates (1.5, 0.5, 1.0) — half a
        // cell past x-index 1 — the world point is origin + that
        // offset = (0.0, -1.0, -0.5). The x-index-valued field there
        // interpolates to exactly 1.5.
        let grid = GridBox::with_spacing([0.0; 3], [3.0; 3], 1.0).unwrap();
        let mut m = AffinityMap::zeros(&grid);
        let (nx, ny, nz) = m.dims;
        for iz in 0..nz {
            for iy in 0..ny {
                for ix in 0..nx {
                    let idx = m.index(ix, iy, iz);
                    m.data[idx] = ix as f64;
                }
            }
        }
        let v = m.interpolate(Vector3::new(0.0, -1.0, -0.5));
        assert!((v - 1.5).abs() < 1e-12, "got {v}");
    }

    #[test]
    fn interpolate_out_of_box_is_zero() {
        let grid = GridBox::with_spacing([0.0; 3], [3.0; 3], 1.0).unwrap();
        let m = AffinityMap::zeros(&grid);
        assert_eq!(m.interpolate(Vector3::new(-5.0, 0.0, 0.0)), 0.0);
        assert_eq!(m.interpolate(Vector3::new(100.0, 0.0, 0.0)), 0.0);
    }

    #[test]
    fn grid_score_matches_direct_pairwise_within_tolerance() {
        // The whole point of grids: a grid score should track the
        // direct pairwise score. Build a fine Vina grid and compare.
        let grid = GridBox::with_spacing([0.0; 3], [10.0; 3], 0.25).unwrap();
        let r = carbon_receptor();
        let set = AffinityMapSet::precompute(&r, &[Ad4AtomType::C], &grid, MapKind::Vina).unwrap();
        let pose = vec![(Vector3::new(3.8, 0.0, 0.0), Ad4AtomType::C, 0.0)];
        let grid_score = score_ligand_on_maps(&set, &pose);
        // Direct Vina score of the same single-atom pose.
        let direct = crate::score::vina::score_complex(
            &r,
            &[(Vector3::new(3.8, 0.0, 0.0), Ad4AtomType::C)],
            0,
        )
        .intermolecular();
        assert!(
            (grid_score - direct).abs() < 0.05,
            "grid {grid_score} vs direct {direct}"
        );
    }

    #[test]
    fn ad4_grid_electrostatics_scale_with_ligand_charge() {
        let grid = GridBox::with_spacing([0.0; 3], [8.0; 3], 0.5).unwrap();
        let r = carbon_receptor(); // receptor carbon carries -0.2
        let set =
            AffinityMapSet::precompute(&r, &[Ad4AtomType::N], &grid, MapKind::AutoDock4).unwrap();
        // A +0.5 ligand atom and a -0.5 ligand atom at the same place
        // get opposite-sign electrostatic contributions.
        let p = Vector3::new(3.0, 0.0, 0.0);
        let pos = score_ligand_on_maps(&set, &[(p, Ad4AtomType::N, 0.5)]);
        let neg = score_ligand_on_maps(&set, &[(p, Ad4AtomType::N, -0.5)]);
        // The vdW part is identical; only electrostatics flips. So the
        // two scores must differ.
        assert!((pos - neg).abs() > 1e-9, "charge sign had no effect");
    }
}
