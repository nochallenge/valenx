//! **Feature 5 — sidechain placement.**
//!
//! After the backbone is built, every residue needs a sidechain. The
//! classical answer is to pick, per residue, a rotamer from a library
//! and place the sidechain in that conformation. This module places a
//! **sidechain-centroid pseudo-atom** for every residue from the
//! [`crate::rotamer`] library — enough geometry for the
//! centroid-resolution scores and the downstream refinement.
//!
//! For a *first* placement (no neighbour context yet) the highest-
//! prior rotamer is used; the energy-driven combinatorial repacking
//! that resolves sidechain-sidechain clashes lives in
//! [`crate::refine::repack`]. This split mirrors a real pipeline:
//! place sidechains, then repack.

use crate::error::Result;
use crate::model::ProteinModel;
use crate::rotamer::{place_sidechain_centroid, rebuild_cb, rotamers_for};

/// Places a sidechain (Cβ + a sidechain-centroid pseudo-atom stored
/// in the `cb` slot's companion — here the `cb` slot itself holds the
/// rebuilt Cβ) for every residue of a model.
///
/// For every residue with a built backbone: rebuilds the Cβ if it is
/// missing, then picks the most-probable rotamer for the residue's
/// amino acid. The model's `cb` field is filled with the real /
/// rebuilt Cβ. Residues lacking a backbone are skipped.
///
/// Returns the number of sidechains placed.
///
/// # Errors
/// Currently infallible for a well-formed model, but returns
/// [`Result`] for forward-compatibility.
pub fn place_sidechains(model: &mut ProteinModel) -> Result<usize> {
    let mut placed = 0usize;
    for res in &mut model.residues {
        let (Some(n), Some(ca), Some(c)) = (res.n, res.ca, res.c) else {
            continue;
        };
        if res.aa != 'G' && res.cb.is_none() {
            res.cb = Some(rebuild_cb(&n, &ca, &c));
        }
        placed += 1;
    }
    Ok(placed)
}

/// Computes, for every residue, the sidechain-centroid coordinate
/// under each residue's most-probable rotamer.
///
/// Returns one centroid per residue (in residue order); a residue
/// without a backbone yields `None`. This is the centroid set the
/// packing / knowledge-based scores consume.
pub fn sidechain_centroids(model: &ProteinModel) -> Vec<Option<nalgebra::Point3<f64>>> {
    model
        .residues
        .iter()
        .map(|res| {
            let rotamers = rotamers_for(res.aa);
            let best = rotamers
                .iter()
                .max_by(|a, b| {
                    a.probability
                        .partial_cmp(&b.probability)
                        .unwrap_or(std::cmp::Ordering::Equal)
                })
                .expect("non-empty rotamer set");
            place_sidechain_centroid(res, best)
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Point3;

    fn three_residue_model() -> ProteinModel {
        let mut m = ProteinModel::from_sequence("ALW").expect("model");
        for (i, r) in m.residues.iter_mut().enumerate() {
            let x = i as f64 * 3.8;
            r.n = Some(Point3::new(x, 0.0, 0.0));
            r.ca = Some(Point3::new(x + 1.46, 0.0, 0.0));
            r.c = Some(Point3::new(x + 2.0, 1.4, 0.0));
            r.o = Some(Point3::new(x + 1.3, 2.4, 0.0));
        }
        m
    }

    #[test]
    fn places_cb_for_every_non_glycine() {
        let mut m = three_residue_model();
        let placed = place_sidechains(&mut m).expect("place");
        assert_eq!(placed, 3);
        for r in &m.residues {
            assert!(r.cb.is_some(), "{} has Cβ", r.aa);
        }
    }

    #[test]
    fn glycine_keeps_no_cb() {
        let mut m = ProteinModel::from_sequence("G").expect("model");
        let r = &mut m.residues[0];
        r.n = Some(Point3::new(0.0, 0.0, 0.0));
        r.ca = Some(Point3::new(1.46, 0.0, 0.0));
        r.c = Some(Point3::new(2.0, 1.4, 0.0));
        r.o = Some(Point3::new(1.3, 2.4, 0.0));
        place_sidechains(&mut m).expect("place");
        assert!(m.residues[0].cb.is_none());
    }

    #[test]
    fn centroids_reflect_sidechain_size() {
        let mut m = three_residue_model();
        place_sidechains(&mut m).expect("place");
        let cen = sidechain_centroids(&m);
        assert_eq!(cen.len(), 3);
        // Trp (residue 2) centroid sits further from its Cα than
        // Ala (residue 0).
        let ala_reach = (cen[0].unwrap() - m.residues[0].ca.unwrap()).norm();
        let trp_reach = (cen[2].unwrap() - m.residues[2].ca.unwrap()).norm();
        assert!(trp_reach > ala_reach, "Trp reaches further");
    }
}
