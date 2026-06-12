//! Molecular fingerprints and similarity.
//!
//! Four fingerprint kinds, all producing a dense [`FingerprintBits`]:
//!
//! - [`ecfp`] / [`fcfp`] — Morgan circular fingerprints (ECFP /
//!   functional-class FCFP);
//! - [`maccs_keys`] — a MACCS-class structural-key set;
//! - [`path_fingerprint`] — a Daylight-style topological path
//!   fingerprint;
//!
//! plus the [`tanimoto`] / [`dice`] / [`cosine`] / [`tversky`]
//! similarity coefficients.

pub mod bitvec;
pub mod maccs;
pub mod morgan;
pub mod path;
pub mod similarity;

pub use bitvec::FingerprintBits;
pub use maccs::{maccs_key_count, maccs_keys};
pub use morgan::{ecfp, ecfp_features, fcfp};
pub use path::path_fingerprint;
pub use similarity::{cosine, dice, nearest_neighbor, tanimoto, tanimoto_distance, tversky};

use crate::molecule::Molecule;

/// The fingerprint families this crate can compute — a tag for
/// generic / config-driven callers.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum FingerprintKind {
    /// Morgan / ECFP circular fingerprint (carries its radius).
    Ecfp(usize),
    /// Functional-class FCFP circular fingerprint.
    Fcfp(usize),
    /// MACCS-class structural keys (fixed length).
    Maccs,
    /// Daylight-style topological path fingerprint (carries max path
    /// length).
    Path(usize),
}

/// Compute the fingerprint named by `kind`. `n_bits` is ignored for
/// [`FingerprintKind::Maccs`] (its length is fixed).
pub fn compute(mol: &Molecule, kind: FingerprintKind, n_bits: usize) -> FingerprintBits {
    match kind {
        FingerprintKind::Ecfp(r) => ecfp(mol, r, n_bits),
        FingerprintKind::Fcfp(r) => fcfp(mol, r, n_bits),
        FingerprintKind::Maccs => maccs_keys(mol),
        FingerprintKind::Path(l) => path_fingerprint(mol, l, n_bits),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mol_from_smiles;

    #[test]
    fn compute_dispatches() {
        let m = mol_from_smiles("c1ccccc1").unwrap();
        for kind in [
            FingerprintKind::Ecfp(2),
            FingerprintKind::Fcfp(2),
            FingerprintKind::Maccs,
            FingerprintKind::Path(7),
        ] {
            let fp = compute(&m, kind, 1024);
            assert!(fp.count_ones() > 0, "{kind:?} produced no bits");
        }
    }

    #[test]
    fn similarity_pipeline() {
        let a = mol_from_smiles("c1ccccc1").unwrap();
        let b = mol_from_smiles("c1ccccc1C").unwrap(); // toluene
        let fa = ecfp(&a, 2, 2048);
        let fb = ecfp(&b, 2, 2048);
        let sim = tanimoto(&fa, &fb);
        // similar but not identical
        assert!(sim > 0.0 && sim < 1.0, "benzene/toluene tanimoto = {sim}");
    }
}
