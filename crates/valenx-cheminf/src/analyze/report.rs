//! Batch molecule analysis — the [`MoleculeReport`] bundle and
//! SDF-batch utilities.
//!
//! [`MoleculeReport`] runs every headline analysis on one molecule —
//! formula, weights, descriptors, Lipinski / Veber verdicts, QED, the
//! canonical SMILES and a fingerprint — and bundles them into a single
//! struct. [`analyze_sdf`] / [`filter_sdf`] apply the report across an
//! SD file: the typical "score this library" workflow.

use crate::descriptors::{self, LipinskiResult, VeberResult};
use crate::error::Result;
use crate::fingerprint::{ecfp, FingerprintBits};
use crate::molecule::Molecule;
use crate::perceive::formula::{average_molecular_weight, molecular_formula, monoisotopic_mass};
use crate::smiles::write_canonical_smiles;

/// A full single-molecule analysis bundle.
#[derive(Clone, Debug, PartialEq)]
pub struct MoleculeReport {
    /// Molecule title (from SDF) — empty if unnamed.
    pub name: String,
    /// Canonical SMILES.
    pub canonical_smiles: String,
    /// Hill-system molecular formula.
    pub formula: String,
    /// Average molecular weight (g/mol).
    pub average_mw: f64,
    /// Monoisotopic (exact) mass (u).
    pub monoisotopic_mass: f64,
    /// Heavy-atom count.
    pub heavy_atoms: usize,
    /// Crippen logP.
    pub logp: f64,
    /// Topological polar surface area (Å²).
    pub tpsa: f64,
    /// Hydrogen-bond donor count.
    pub hbd: usize,
    /// Hydrogen-bond acceptor count.
    pub hba: usize,
    /// Rotatable-bond count.
    pub rotatable_bonds: usize,
    /// SSSR ring count.
    pub ring_count: usize,
    /// Aromatic-ring count.
    pub aromatic_rings: usize,
    /// Fraction of `sp³` carbons.
    pub fraction_csp3: f64,
    /// Net formal charge.
    pub formal_charge: i32,
    /// Lipinski rule-of-five verdict.
    pub lipinski: LipinskiResult,
    /// Veber rule verdict.
    pub veber: VeberResult,
    /// QED drug-likeness score (0..1).
    pub qed: f64,
    /// An ECFP4 (radius 2, 2048-bit) fingerprint.
    pub fingerprint: FingerprintBits,
}

impl MoleculeReport {
    /// Compute the full report for `mol`.
    ///
    /// `mol` should already have perception run ([`crate::mol_from_smiles`]
    /// does this); the report re-perceives defensively.
    pub fn analyze(mol: &Molecule) -> Self {
        let mut m = mol.clone();
        crate::perceive::perceive_all(&mut m);
        MoleculeReport {
            name: m.name.clone(),
            canonical_smiles: write_canonical_smiles(&m),
            formula: molecular_formula(&m),
            average_mw: average_molecular_weight(&m),
            monoisotopic_mass: monoisotopic_mass(&m),
            heavy_atoms: m.heavy_atom_count(),
            logp: descriptors::crippen_logp(&m),
            tpsa: descriptors::tpsa(&m),
            hbd: descriptors::hbd(&m),
            hba: descriptors::hba(&m),
            rotatable_bonds: descriptors::rotatable_bonds(&m),
            ring_count: descriptors::ring_count(&m),
            aromatic_rings: descriptors::aromatic_ring_count(&m),
            fraction_csp3: descriptors::fraction_csp3(&m),
            formal_charge: descriptors::formal_charge(&m),
            lipinski: descriptors::lipinski(&m),
            veber: descriptors::veber(&m),
            qed: super::qed::qed(&m),
            fingerprint: ecfp(&m, 2, 2048),
        }
    }

    /// Analyse a molecule directly from a SMILES string.
    pub fn from_smiles(smiles: &str) -> Result<Self> {
        let mol = crate::mol_from_smiles(smiles)?;
        Ok(MoleculeReport::analyze(&mol))
    }

    /// `true` if the molecule is "drug-like" — passes Lipinski (≤ 1
    /// violation) and Veber.
    pub fn is_drug_like(&self) -> bool {
        self.lipinski.passes() && self.veber.passes
    }

    /// A compact one-line text summary, handy for logs / CSV rows.
    pub fn summary_line(&self) -> String {
        format!(
            "{} | {} | MW {:.1} | logP {:.2} | HBD {} HBA {} | TPSA {:.1} | QED {:.3} | {}",
            if self.name.is_empty() {
                "(unnamed)"
            } else {
                &self.name
            },
            self.formula,
            self.average_mw,
            self.logp,
            self.hbd,
            self.hba,
            self.tpsa,
            self.qed,
            if self.is_drug_like() {
                "drug-like"
            } else {
                "non-drug-like"
            },
        )
    }
}

/// Analyse every molecule in an SD file, returning one
/// [`MoleculeReport`] per record.
pub fn analyze_sdf(sdf_text: &str) -> Result<Vec<MoleculeReport>> {
    let molecules = crate::molfile::read_sdf(sdf_text)?;
    Ok(molecules.iter().map(MoleculeReport::analyze).collect())
}

/// Analyse every molecule in a slice — the in-memory equivalent of
/// [`analyze_sdf`].
pub fn analyze_batch(molecules: &[Molecule]) -> Vec<MoleculeReport> {
    molecules.iter().map(MoleculeReport::analyze).collect()
}

/// Read an SD file and keep only the records whose report satisfies
/// `predicate` — the "filter a library" workflow. Returns the surviving
/// molecules (with their original SDF properties intact).
pub fn filter_sdf(
    sdf_text: &str,
    predicate: impl Fn(&MoleculeReport) -> bool,
) -> Result<Vec<Molecule>> {
    let molecules = crate::molfile::read_sdf(sdf_text)?;
    Ok(molecules
        .into_iter()
        .filter(|m| predicate(&MoleculeReport::analyze(m)))
        .collect())
}

/// Summary statistics over a batch of reports — the headline numbers a
/// library triage wants.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct BatchSummary {
    /// Number of molecules analysed.
    pub count: usize,
    /// Mean average-molecular-weight.
    pub mean_mw: f64,
    /// Mean Crippen logP.
    pub mean_logp: f64,
    /// Mean QED.
    pub mean_qed: f64,
    /// Number passing the drug-likeness filter.
    pub drug_like: usize,
}

/// Roll a slice of reports into a [`BatchSummary`].
pub fn summarize(reports: &[MoleculeReport]) -> BatchSummary {
    if reports.is_empty() {
        return BatchSummary::default();
    }
    let n = reports.len() as f64;
    BatchSummary {
        count: reports.len(),
        mean_mw: reports.iter().map(|r| r.average_mw).sum::<f64>() / n,
        mean_logp: reports.iter().map(|r| r.logp).sum::<f64>() / n,
        mean_qed: reports.iter().map(|r| r.qed).sum::<f64>() / n,
        drug_like: reports.iter().filter(|r| r.is_drug_like()).count(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::molfile::write_sdf;

    #[test]
    fn report_has_all_fields() {
        let r = MoleculeReport::from_smiles("CCO").unwrap();
        assert_eq!(r.formula, "C2H6O");
        assert!(r.average_mw > 45.0 && r.average_mw < 47.0);
        assert_eq!(r.hbd, 1);
        assert_eq!(r.hba, 1);
        assert!(r.fingerprint.count_ones() > 0);
        assert!((0.0..=1.0).contains(&r.qed));
    }

    #[test]
    fn drug_like_classification() {
        let small = MoleculeReport::from_smiles("c1ccccc1C(=O)NCCO").unwrap();
        assert!(small.is_drug_like());
        // a huge greasy chain fails Lipinski (logP / MW)
        let big = MoleculeReport::from_smiles("CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC").unwrap();
        assert!(!big.is_drug_like());
    }

    #[test]
    fn summary_line_is_readable() {
        let r = MoleculeReport::from_smiles("CCO").unwrap();
        let line = r.summary_line();
        assert!(line.contains("C2H6O"));
        assert!(line.contains("MW"));
    }

    #[test]
    fn sdf_batch_round_trip() {
        // build a 3-molecule SDF, analyse it back
        let mut mols = vec![
            crate::mol_from_smiles("CCO").unwrap(),
            crate::mol_from_smiles("c1ccccc1").unwrap(),
            crate::mol_from_smiles("CC(=O)O").unwrap(),
        ];
        for (i, m) in mols.iter_mut().enumerate() {
            m.name = format!("mol{i}");
        }
        let sdf = write_sdf(&mols);
        let reports = analyze_sdf(&sdf).unwrap();
        assert_eq!(reports.len(), 3);
        assert_eq!(reports[0].name, "mol0");
        assert_eq!(reports[1].formula, "C6H6");
    }

    #[test]
    fn filter_sdf_keeps_matching() {
        let mols = vec![
            crate::mol_from_smiles("CCO").unwrap(),
            crate::mol_from_smiles("CCCCCCCCCCCCCCCCCCCC").unwrap(),
        ];
        let sdf = write_sdf(&mols);
        // keep only molecules with ≤ 5 heavy atoms
        let kept = filter_sdf(&sdf, |r| r.heavy_atoms <= 5).unwrap();
        assert_eq!(kept.len(), 1);
        assert_eq!(kept[0].heavy_atom_count(), 3);
    }

    #[test]
    fn batch_summary() {
        let reports = analyze_batch(&[
            crate::mol_from_smiles("CCO").unwrap(),
            crate::mol_from_smiles("c1ccccc1").unwrap(),
        ]);
        let s = summarize(&reports);
        assert_eq!(s.count, 2);
        assert!(s.mean_mw > 0.0);
        assert!((0.0..=1.0).contains(&s.mean_qed));
    }
}
