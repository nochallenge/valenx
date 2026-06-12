//! 3D conformer generation by distance geometry.
//!
//! [`embed_3d`] generates a 3D conformer the way RDKit's
//! `EmbedMolecule` does in outline:
//!
//! 1. **bounds matrix** — for every atom pair compute a lower and an
//!    upper distance bound. 1,2 pairs (bonded) get the covalent bond
//!    length; 1,3 pairs (a shared neighbour) get the distance implied
//!    by an idealised valence angle; everything else gets `[sum of vdW
//!    radii, large]` — then the bounds are *triangle-smoothed* so they
//!    are geometrically consistent;
//! 2. **metric-matrix embedding** — sample a trial distance for each
//!    pair uniformly within its bounds, build the Gram matrix, take
//!    its top three eigenvectors (`nalgebra` symmetric eigensolver) →
//!    an initial set of 3D coordinates;
//! 3. the raw embedding is then refined by [`super::super::forcefield`]
//!    so bond lengths and angles relax to chemically sensible values.
//!
//! A deterministic seedable PCG generator drives the trial-distance
//! sampling so a given `(molecule, seed)` always yields the same
//! conformer.
//!
//! **v1 simplifications:** the bounds use idealised tetrahedral /
//! trigonal / linear angles by hybridisation guess, not a full
//! experimental-geometry table; chirality and double-bond geometry are
//! not *enforced* during embedding (the force-field cleanup keeps the
//! geometry reasonable but a specified `@`/`E`/`Z` is not guaranteed in
//! the output) — that distance-geometry chirality constraint is a
//! documented limitation.

use crate::error::{CheminfError, Result};
use crate::molecule::{BondOrder, Molecule};
use crate::perceive::hydrogen::add_explicit_hydrogens;
use nalgebra::{DMatrix, SymmetricEigen};

/// A small deterministic PCG-XSH-RR random generator so conformer
/// embedding is reproducible for a given seed.
#[derive(Clone, Debug)]
pub struct Pcg32 {
    state: u64,
    inc: u64,
}

impl Pcg32 {
    /// Seed the generator.
    pub fn new(seed: u64) -> Self {
        let mut g = Pcg32 {
            state: 0,
            inc: (seed << 1) | 1,
        };
        g.next_u32();
        g.state = g.state.wrapping_add(seed);
        g.next_u32();
        g
    }

    /// Next pseudo-random `u32`.
    pub fn next_u32(&mut self) -> u32 {
        let old = self.state;
        self.state = old.wrapping_mul(6364136223846793005).wrapping_add(self.inc);
        let xorshifted = (((old >> 18) ^ old) >> 27) as u32;
        let rot = (old >> 59) as u32;
        xorshifted.rotate_right(rot)
    }

    /// Uniform `f64` in `[0, 1)`.
    pub fn next_f64(&mut self) -> f64 {
        (self.next_u32() as f64) / (u32::MAX as f64 + 1.0)
    }

    /// Uniform `f64` in `[lo, hi)`.
    pub fn uniform(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.next_f64()
    }
}

/// Generate a 3D conformer for `mol`, returning a **new** molecule with
/// explicit hydrogens and [`Molecule::coords`] populated
/// (`coords_3d = true`).
///
/// `seed` makes the embedding reproducible.
pub fn embed_3d(mol: &Molecule, seed: u64) -> Result<Molecule> {
    let mut work = add_explicit_hydrogens(mol);
    let n = work.atoms.len();
    if n == 0 {
        return Err(CheminfError::invalid(
            "molecule",
            "cannot embed an empty molecule",
        ));
    }
    if n == 1 {
        work.coords = vec![[0.0, 0.0, 0.0]];
        work.coords_3d = true;
        return Ok(work);
    }

    let (lower, upper) = bounds_matrix(&work);
    let smoothed = triangle_smooth(lower, upper);
    let coords = metric_embed(&smoothed, n, seed);
    work.coords = coords;
    work.coords_3d = true;
    // refine with the force field so bonds / angles relax
    crate::forcefield::clean_up_geometry(&mut work, 200);
    Ok(work)
}

/// Distance-geometry embedding followed by **production MMFF94**
/// cleanup (the new default). Equivalent to [`embed_3d`] but the
/// post-embedding relaxation uses the tabulated MMFF94 force field
/// from [`crate::forcefield_mmff94`] rather than the legacy reduced
/// force field. Recommended path for any caller that wants the
/// best-available conformer geometry from this crate.
pub fn embed_3d_mmff94(mol: &Molecule, seed: u64) -> Result<Molecule> {
    let mut work = add_explicit_hydrogens(mol);
    let n = work.atoms.len();
    if n == 0 {
        return Err(CheminfError::invalid(
            "molecule",
            "cannot embed an empty molecule",
        ));
    }
    if n == 1 {
        work.coords = vec![[0.0, 0.0, 0.0]];
        work.coords_3d = true;
        return Ok(work);
    }
    let (lower, upper) = bounds_matrix(&work);
    let smoothed = triangle_smooth(lower, upper);
    let coords = metric_embed(&smoothed, n, seed);
    work.coords = coords;
    work.coords_3d = true;
    crate::forcefield_mmff94::clean_up_geometry(&mut work, 300);
    Ok(work)
}

/// Build the `(lower, upper)` inter-atomic distance bound matrices.
fn bounds_matrix(mol: &Molecule) -> (DMatrix<f64>, DMatrix<f64>) {
    let n = mol.atoms.len();
    let mut lower = DMatrix::from_element(n, n, 0.0);
    let mut upper = DMatrix::from_element(n, n, 100.0);
    for i in 0..n {
        lower[(i, i)] = 0.0;
        upper[(i, i)] = 0.0;
    }
    // 1,2 — bonded pairs
    for b in &mol.bonds {
        let len = ideal_bond_length(mol, b.a, b.b, b.order);
        set_bound(&mut lower, &mut upper, b.a, b.b, len * 0.95, len * 1.05);
    }
    // 1,3 — pairs sharing a neighbour, from the central valence angle
    for c in 0..n {
        let nbrs = mol.neighbors(c);
        for ii in 0..nbrs.len() {
            for jj in ii + 1..nbrs.len() {
                let (a, b) = (nbrs[ii], nbrs[jj]);
                let l_ac = ideal_bond_length(mol, a, c, bond_order(mol, a, c));
                let l_bc = ideal_bond_length(mol, b, c, bond_order(mol, b, c));
                let angle = ideal_angle(mol, c);
                // law of cosines for the 1,3 distance
                let d = (l_ac * l_ac + l_bc * l_bc - 2.0 * l_ac * l_bc * angle.cos()).sqrt();
                set_bound(&mut lower, &mut upper, a, b, d * 0.93, d * 1.07);
            }
        }
    }
    // everything else: lower bound = sum of vdW radii, big upper bound
    for i in 0..n {
        for j in i + 1..n {
            if upper[(i, j)] >= 99.0 {
                let vdw =
                    vdw_radius(mol.atoms[i].atomic_number) + vdw_radius(mol.atoms[j].atomic_number);
                lower[(i, j)] = vdw * 0.8;
                lower[(j, i)] = vdw * 0.8;
                // upper stays large
            }
        }
    }
    (lower, upper)
}

fn set_bound(
    lower: &mut DMatrix<f64>,
    upper: &mut DMatrix<f64>,
    i: usize,
    j: usize,
    lo: f64,
    hi: f64,
) {
    lower[(i, j)] = lo;
    lower[(j, i)] = lo;
    upper[(i, j)] = hi;
    upper[(j, i)] = hi;
}

/// Ideal bond length = sum of covalent radii, shortened for multiple
/// bonds.
fn ideal_bond_length(mol: &Molecule, a: usize, b: usize, order: BondOrder) -> f64 {
    let ra = crate::element::covalent_radius(mol.atoms[a].atomic_number);
    let rb = crate::element::covalent_radius(mol.atoms[b].atomic_number);
    let base = ra + rb;
    match order {
        BondOrder::Double => base * 0.87,
        BondOrder::Triple => base * 0.78,
        BondOrder::Aromatic => base * 0.91,
        _ => base,
    }
}

fn bond_order(mol: &Molecule, a: usize, b: usize) -> BondOrder {
    mol.bond_between(a, b)
        .map(|bi| mol.bonds[bi].order)
        .unwrap_or(BondOrder::Single)
}

/// Ideal valence angle at `center` from a hybridisation guess.
fn ideal_angle(mol: &Molecule, center: usize) -> f64 {
    let a = &mol.atoms[center];
    let multiples = mol
        .bonds_on(center)
        .iter()
        .filter(|&&b| {
            matches!(
                mol.bonds[b].order,
                BondOrder::Double | BondOrder::Triple | BondOrder::Aromatic
            )
        })
        .count();
    let triple = mol
        .bonds_on(center)
        .iter()
        .any(|&b| mol.bonds[b].order == BondOrder::Triple);
    if triple {
        std::f64::consts::PI // 180° sp
    } else if multiples > 0 || a.aromatic {
        2.094 // 120° sp2
    } else {
        1.911 // 109.47° sp3
    }
}

/// van der Waals radius for the non-bonded lower bound.
fn vdw_radius(z: u8) -> f64 {
    match z {
        1 => 1.20,
        6 => 1.70,
        7 => 1.55,
        8 => 1.52,
        9 => 1.47,
        15 => 1.80,
        16 => 1.80,
        17 => 1.75,
        35 => 1.85,
        53 => 1.98,
        _ => 1.70,
    }
}

/// Triangle-inequality smoothing: tighten upper bounds and raise lower
/// bounds until they are geometrically consistent (a Floyd-Warshall
/// pass on the upper bounds, then a lower-bound pass).
fn triangle_smooth(
    mut lower: DMatrix<f64>,
    mut upper: DMatrix<f64>,
) -> (DMatrix<f64>, DMatrix<f64>) {
    let n = upper.nrows();
    // upper bounds: u[i][j] ≤ u[i][k] + u[k][j]
    for k in 0..n {
        for i in 0..n {
            for j in 0..n {
                let via = upper[(i, k)] + upper[(k, j)];
                if via < upper[(i, j)] {
                    upper[(i, j)] = via;
                }
            }
        }
    }
    // lower bounds: l[i][j] ≥ l[i][k] − u[k][j]
    for k in 0..n {
        for i in 0..n {
            for j in 0..n {
                let cand = lower[(i, k)] - upper[(k, j)];
                if cand > lower[(i, j)] {
                    lower[(i, j)] = cand;
                }
                // keep lower ≤ upper
                if lower[(i, j)] > upper[(i, j)] {
                    lower[(i, j)] = upper[(i, j)];
                }
            }
        }
    }
    (lower, upper)
}

/// Metric-matrix embedding: sample distances, build the Gram matrix,
/// take its top-3 eigenvectors as 3D coordinates.
fn metric_embed(bounds: &(DMatrix<f64>, DMatrix<f64>), n: usize, seed: u64) -> Vec<[f64; 3]> {
    let (lower, upper) = bounds;
    let mut rng = Pcg32::new(seed.wrapping_add(0x9E3779B9));
    // trial distance matrix d[i][j] sampled in [lower, upper]
    let mut d2 = DMatrix::from_element(n, n, 0.0);
    for i in 0..n {
        for j in i + 1..n {
            let lo = lower[(i, j)];
            let hi = upper[(i, j)].min(lo + 8.0).max(lo);
            let dist = rng.uniform(lo, hi.max(lo + 1e-3));
            d2[(i, j)] = dist * dist;
            d2[(j, i)] = dist * dist;
        }
    }
    // metric (Gram) matrix from squared distances:
    // g[i][j] = (d0i^2 + d0j^2 - dij^2)/2 with the centroid as origin
    // use the standard double-centering of the squared-distance matrix
    let mut row_mean = vec![0.0; n];
    let mut total = 0.0;
    for i in 0..n {
        let mut s = 0.0;
        for j in 0..n {
            s += d2[(i, j)];
        }
        row_mean[i] = s / n as f64;
        total += s;
    }
    let grand = total / (n as f64 * n as f64);
    let mut gram = DMatrix::from_element(n, n, 0.0);
    for i in 0..n {
        for j in 0..n {
            gram[(i, j)] = -0.5 * (d2[(i, j)] - row_mean[i] - row_mean[j] + grand);
        }
    }
    // symmetric eigendecomposition; top 3 eigenpairs → coordinates
    let eig = SymmetricEigen::new(gram);
    let mut pairs: Vec<(f64, usize)> = eig
        .eigenvalues
        .iter()
        .copied()
        .enumerate()
        .map(|(i, v)| (v, i))
        .collect();
    pairs.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap_or(std::cmp::Ordering::Equal));

    let mut coords = vec![[0.0f64; 3]; n];
    for (axis, &(eigval, col)) in pairs.iter().take(3).enumerate() {
        let scale = eigval.max(0.0).sqrt();
        for (atom, c) in coords.iter_mut().enumerate() {
            c[axis] = eig.eigenvectors[(atom, col)] * scale;
        }
    }
    coords
}

/// Root-mean-square deviation of `mol`'s actual bond lengths from their
/// ideal values — a conformer-quality probe.
pub fn bond_length_rmsd(mol: &Molecule) -> f64 {
    if mol.coords.is_empty() || mol.bonds.is_empty() {
        return 0.0;
    }
    let mut sq = 0.0;
    for b in &mol.bonds {
        let ideal = ideal_bond_length(mol, b.a, b.b, b.order);
        let p = mol.coords[b.a];
        let q = mol.coords[b.b];
        let d = ((p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2) + (p[2] - q[2]).powi(2)).sqrt();
        sq += (d - ideal).powi(2);
    }
    (sq / mol.bonds.len() as f64).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mol_from_smiles;

    #[test]
    fn pcg_is_deterministic() {
        let mut a = Pcg32::new(42);
        let mut b = Pcg32::new(42);
        for _ in 0..50 {
            assert_eq!(a.next_u32(), b.next_u32());
        }
    }

    #[test]
    fn embed_produces_3d_coords() {
        let m = mol_from_smiles("CCO").unwrap();
        let conf = embed_3d(&m, 1).unwrap();
        assert!(conf.coords_3d);
        assert_eq!(conf.coords.len(), conf.atom_count());
        // explicit hydrogens were added (ethanol C2H6O → 9 atoms)
        assert_eq!(conf.atom_count(), 9);
        // genuine 3D: at least one atom off the xy-plane
        let has_z = conf.coords.iter().any(|c| c[2].abs() > 1e-3);
        assert!(has_z, "embedding is flat");
    }

    #[test]
    fn embedding_is_reproducible() {
        let m = mol_from_smiles("CCN").unwrap();
        let a = embed_3d(&m, 7).unwrap();
        let b = embed_3d(&m, 7).unwrap();
        assert_eq!(a.coords, b.coords);
    }

    #[test]
    fn bond_lengths_are_reasonable() {
        let m = mol_from_smiles("CCCC").unwrap();
        let conf = embed_3d(&m, 3).unwrap();
        let rmsd = bond_length_rmsd(&conf);
        // after force-field cleanup bonds should be within ~0.3 A of ideal
        assert!(rmsd < 0.5, "bond-length RMSD too large: {rmsd}");
    }

    #[test]
    fn single_atom_embeds() {
        let m = mol_from_smiles("[Ne]").unwrap();
        let conf = embed_3d(&m, 1).unwrap();
        assert_eq!(conf.coords.len(), 1);
    }

    #[test]
    fn benzene_embeds() {
        let m = mol_from_smiles("c1ccccc1").unwrap();
        let conf = embed_3d(&m, 5).unwrap();
        // benzene C6H6 → 12 atoms
        assert_eq!(conf.atom_count(), 12);
        assert!(conf.coords_3d);
    }
}
