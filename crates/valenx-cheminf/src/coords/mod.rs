//! Molecular coordinates — 2D depiction and 3D conformer generation.
//!
//! - [`compute_2d_coords`] lays a molecule out for drawing (ring
//!   polygons + chain growth + repulsion relaxation);
//! - [`embed_3d`] generates a single conformer by distance geometry
//!   (bounds matrix → metric embedding → force-field cleanup);
//! - [`etkdg`] is the production-depth conformer generator —
//!   distance-geometry biased by the Riniker-Landrum experimental
//!   torsion library, with multi-conformer generation, MMFF94 cleanup
//!   and RMSD pruning.
//!
//! Both store their result in [`Molecule::coords`]; the
//! [`Molecule::coords_3d`] flag distinguishes a depiction from a
//! conformer.

pub mod embed3d;
pub mod etkdg;
pub mod layout2d;

pub use embed3d::{bond_length_rmsd, embed_3d, embed_3d_mmff94, Pcg32};
pub use etkdg::{
    etkdg_embed, generate_conformers, heavy_atom_rmsd, torsion_pref, EtkdgOptions, ScoredConformer,
    TorsionClass, TorsionPref, TORSION_LIBRARY,
};
pub use layout2d::{compute_2d_coords, mean_bond_length_2d};

use crate::molecule::Molecule;

/// Geometric centroid of a molecule's coordinates, or the origin if it
/// has none.
pub fn centroid(mol: &Molecule) -> [f64; 3] {
    if mol.coords.is_empty() {
        return [0.0, 0.0, 0.0];
    }
    let n = mol.coords.len() as f64;
    let mut c = [0.0, 0.0, 0.0];
    for p in &mol.coords {
        c[0] += p[0];
        c[1] += p[1];
        c[2] += p[2];
    }
    [c[0] / n, c[1] / n, c[2] / n]
}

/// Translate `mol` so its centroid sits at the origin (mutates the
/// coordinates in place).
pub fn center_at_origin(mol: &mut Molecule) {
    let c = centroid(mol);
    for p in &mut mol.coords {
        p[0] -= c[0];
        p[1] -= c[1];
        p[2] -= c[2];
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mol_from_smiles;

    #[test]
    fn centroid_and_centering() {
        let mut m = mol_from_smiles("CCO").unwrap();
        compute_2d_coords(&mut m);
        center_at_origin(&mut m);
        let c = centroid(&m);
        assert!(c[0].abs() < 1e-9 && c[1].abs() < 1e-9);
    }
}
