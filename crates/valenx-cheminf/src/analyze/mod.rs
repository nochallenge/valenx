//! High-level molecular analysis.
//!
//! - [`pharmacophore`](mod@pharmacophore) — a pharmacophore feature
//!   model (donor / acceptor / aromatic / hydrophobic / charged
//!   feature points with 3D positions);
//! - [`qed`](mod@qed) — the Quantitative Estimate of Drug-likeness;
//! - [`report`] — the batch [`MoleculeReport`] bundle and SDF-batch
//!   processing utilities.

pub mod pharmacophore;
pub mod qed;
pub mod report;

pub use pharmacophore::{
    feature_count, feature_distances, pharmacophore, FeatureKind, PharmacophoreFeature,
};
pub use qed::{qed, qed_components, structural_alert_count, QedComponents};
pub use report::{analyze_batch, analyze_sdf, filter_sdf, summarize, BatchSummary, MoleculeReport};

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mol_from_smiles;

    #[test]
    fn analyze_surface_smoke() {
        let m = mol_from_smiles("CC(=O)Nc1ccccc1").unwrap(); // acetanilide
        let report = MoleculeReport::analyze(&m);
        assert!(report.qed > 0.0);
        let feats = pharmacophore(&m);
        assert!(!feats.is_empty());
        assert!((0.0..=1.0).contains(&qed(&m)));
    }
}
