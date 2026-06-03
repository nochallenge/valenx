//! Feature 23 — ADMET-lite library filtering for screening.
//!
//! Before docking a library it pays to *triage* it: a compound that
//! will never be a drug — too big, too greasy, riddled with reactive
//! groups — is not worth a docking slot. This is the
//! "lead-likeness / drug-likeness" pre-filter every virtual screen
//! applies.
//!
//! All the underlying chemistry lives in [`valenx_cheminf`]:
//! [`valenx_cheminf::descriptors::lipinski`] (rule of five),
//! [`valenx_cheminf::descriptors::veber`] (oral-bioavailability
//! rules) and [`valenx_cheminf::analyze::qed::structural_alert_count`]
//! (a PAINS-lite reactive-group / promiscuity flag). This module
//! reuses them and adds the *filter*: an [`AdmetFilter`] of tunable
//! thresholds, an [`AdmetVerdict`] per molecule, and a batch
//! [`admet_filter`] that partitions a library into pass / fail.

use valenx_cheminf::analyze::qed::{qed, structural_alert_count};
use valenx_cheminf::descriptors::{
    crippen_logp, hba, hbd, lipinski, rotatable_bonds, tpsa, veber,
};
use valenx_cheminf::molecule::Molecule;
use valenx_cheminf::perceive::formula::average_molecular_weight;

use crate::error::{DockScreenError, Result};

/// Tunable thresholds for the ADMET-lite filter. The defaults are the
/// classical drug-likeness limits (Lipinski + Veber + a structural
/// alert cap).
#[derive(Clone, Copy, Debug)]
pub struct AdmetFilter {
    /// Maximum allowed average molecular weight (g/mol).
    pub max_molecular_weight: f64,
    /// Maximum allowed Crippen logP.
    pub max_logp: f64,
    /// Maximum allowed hydrogen-bond donor count.
    pub max_hbd: usize,
    /// Maximum allowed hydrogen-bond acceptor count.
    pub max_hba: usize,
    /// Maximum allowed rotatable-bond count (Veber).
    pub max_rotatable_bonds: usize,
    /// Maximum allowed topological polar surface area (Å², Veber).
    pub max_tpsa: f64,
    /// Maximum allowed number of structural alerts (PAINS-lite).
    pub max_structural_alerts: usize,
    /// How many Lipinski rule-of-five violations to tolerate (the
    /// classical "drug-like" cutoff is 1).
    pub allowed_lipinski_violations: u8,
}

impl Default for AdmetFilter {
    fn default() -> Self {
        AdmetFilter {
            max_molecular_weight: 500.0,
            max_logp: 5.0,
            max_hbd: 5,
            max_hba: 10,
            max_rotatable_bonds: 10,
            max_tpsa: 140.0,
            max_structural_alerts: 0,
            allowed_lipinski_violations: 1,
        }
    }
}

impl AdmetFilter {
    /// A stricter "lead-like" filter (smaller, less lipophilic
    /// molecules — the Teague lead-likeness window).
    pub fn lead_like() -> Self {
        AdmetFilter {
            max_molecular_weight: 350.0,
            max_logp: 3.5,
            max_hbd: 3,
            max_hba: 6,
            max_rotatable_bonds: 7,
            max_tpsa: 120.0,
            max_structural_alerts: 0,
            allowed_lipinski_violations: 0,
        }
    }

    /// A permissive "fragment-like" filter (the rule of three).
    pub fn fragment_like() -> Self {
        AdmetFilter {
            max_molecular_weight: 300.0,
            max_logp: 3.0,
            max_hbd: 3,
            max_hba: 3,
            max_rotatable_bonds: 3,
            max_tpsa: 90.0,
            max_structural_alerts: 1,
            allowed_lipinski_violations: 1,
        }
    }
}

/// The ADMET-lite verdict for one molecule.
#[derive(Clone, Debug, PartialEq)]
pub struct AdmetVerdict {
    /// The molecule's index in the input library.
    pub molecule_index: usize,
    /// `true` if the molecule passed every filter criterion.
    pub passes: bool,
    /// Average molecular weight (g/mol).
    pub molecular_weight: f64,
    /// Crippen logP.
    pub logp: f64,
    /// Hydrogen-bond donor count.
    pub hbd: usize,
    /// Hydrogen-bond acceptor count.
    pub hba: usize,
    /// Rotatable-bond count.
    pub rotatable_bonds: usize,
    /// Topological polar surface area (Å²).
    pub tpsa: f64,
    /// Structural-alert (PAINS-lite) count.
    pub structural_alerts: usize,
    /// Number of Lipinski rule-of-five violations.
    pub lipinski_violations: u8,
    /// QED drug-likeness score `[0,1]` (informational — not a pass /
    /// fail criterion here).
    pub qed: f64,
    /// Human-readable reasons the molecule failed (empty if it passed).
    pub failures: Vec<String>,
}

/// The result of filtering a library.
#[derive(Clone, Debug)]
pub struct AdmetFilterResult {
    /// One verdict per input molecule, in input order.
    pub verdicts: Vec<AdmetVerdict>,
}

impl AdmetFilterResult {
    /// Indices of the molecules that passed.
    pub fn passing_indices(&self) -> Vec<usize> {
        self.verdicts
            .iter()
            .filter(|v| v.passes)
            .map(|v| v.molecule_index)
            .collect()
    }

    /// Number of molecules that passed.
    pub fn n_passed(&self) -> usize {
        self.verdicts.iter().filter(|v| v.passes).count()
    }

    /// Number of molecules that failed.
    pub fn n_failed(&self) -> usize {
        self.verdicts.iter().filter(|v| !v.passes).count()
    }
}

/// Evaluate one molecule against an ADMET-lite filter.
pub fn evaluate(index: usize, mol: &Molecule, filter: &AdmetFilter) -> AdmetVerdict {
    let mw = average_molecular_weight(mol);
    let lp = crippen_logp(mol);
    let d = hbd(mol);
    let a = hba(mol);
    let rb = rotatable_bonds(mol);
    let psa = tpsa(mol);
    let alerts = structural_alert_count(mol);
    let lip = lipinski(mol);
    let q = qed(mol);
    let _ = veber(mol); // exercised for parity; thresholds applied directly below

    let mut failures: Vec<String> = Vec::new();
    if mw > filter.max_molecular_weight {
        failures.push(format!(
            "molecular weight {mw:.1} > {:.1}",
            filter.max_molecular_weight
        ));
    }
    if lp > filter.max_logp {
        failures.push(format!("logP {lp:.2} > {:.2}", filter.max_logp));
    }
    if d > filter.max_hbd {
        failures.push(format!("HBD {d} > {}", filter.max_hbd));
    }
    if a > filter.max_hba {
        failures.push(format!("HBA {a} > {}", filter.max_hba));
    }
    if rb > filter.max_rotatable_bonds {
        failures.push(format!(
            "rotatable bonds {rb} > {}",
            filter.max_rotatable_bonds
        ));
    }
    if psa > filter.max_tpsa {
        failures.push(format!("TPSA {psa:.1} > {:.1}", filter.max_tpsa));
    }
    if alerts > filter.max_structural_alerts {
        failures.push(format!(
            "structural alerts {alerts} > {}",
            filter.max_structural_alerts
        ));
    }
    if lip.violations > filter.allowed_lipinski_violations {
        failures.push(format!(
            "Lipinski violations {} > {}",
            lip.violations, filter.allowed_lipinski_violations
        ));
    }

    AdmetVerdict {
        molecule_index: index,
        passes: failures.is_empty(),
        molecular_weight: mw,
        logp: lp,
        hbd: d,
        hba: a,
        rotatable_bonds: rb,
        tpsa: psa,
        structural_alerts: alerts,
        lipinski_violations: lip.violations,
        qed: q,
        failures,
    }
}

/// Feature 23 — filter a molecule library with an ADMET-lite filter.
///
/// Every molecule is evaluated; the result carries a pass / fail
/// verdict and the property values for each. Use
/// [`AdmetFilterResult::passing_indices`] to get the screening-ready
/// subset.
///
/// Returns [`DockScreenError::Invalid`] for an empty library.
pub fn admet_filter(library: &[Molecule], filter: &AdmetFilter) -> Result<AdmetFilterResult> {
    if library.is_empty() {
        return Err(DockScreenError::invalid(
            "library",
            "cannot filter an empty molecule library",
        ));
    }
    let verdicts: Vec<AdmetVerdict> = library
        .iter()
        .enumerate()
        .map(|(i, m)| evaluate(i, m, filter))
        .collect();
    Ok(AdmetFilterResult { verdicts })
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_cheminf::mol_from_smiles;

    #[test]
    fn rejects_empty_library() {
        assert!(admet_filter(&[], &AdmetFilter::default()).is_err());
    }

    #[test]
    fn a_small_drug_like_molecule_passes() {
        // Benzene — tiny, well within every drug-likeness limit.
        let benzene = mol_from_smiles("c1ccccc1").unwrap();
        let result = admet_filter(&[benzene], &AdmetFilter::default()).unwrap();
        assert!(result.verdicts[0].passes, "benzene should pass: {:?}", result.verdicts[0].failures);
        assert_eq!(result.n_passed(), 1);
    }

    #[test]
    fn a_huge_greasy_molecule_fails_the_filter() {
        // A long alkane — high molecular weight and very high logP.
        let big = mol_from_smiles("CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC").unwrap();
        let result = admet_filter(&[big], &AdmetFilter::default()).unwrap();
        assert!(!result.verdicts[0].passes, "a C40 alkane should fail");
        assert!(!result.verdicts[0].failures.is_empty());
        assert_eq!(result.n_failed(), 1);
    }

    #[test]
    fn lead_like_filter_is_stricter_than_default() {
        // A molecule that passes the default filter but is too large
        // for the lead-like window.
        let medium = mol_from_smiles("CCCCCCCCCCCCCCCCCCCCCCCC").unwrap(); // C24
        let by_default = admet_filter(&[medium.clone()], &AdmetFilter::default())
            .unwrap()
            .verdicts[0]
            .clone();
        let by_lead = admet_filter(&[medium], &AdmetFilter::lead_like())
            .unwrap()
            .verdicts[0]
            .clone();
        // The lead-like filter cannot pass more than the default does:
        // if it passed, the default must have passed too.
        assert!(!by_lead.passes || by_default.passes);
    }

    #[test]
    fn passing_indices_lists_only_passes() {
        let benzene = mol_from_smiles("c1ccccc1").unwrap();
        let big = mol_from_smiles("CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC").unwrap();
        let result = admet_filter(&[benzene, big], &AdmetFilter::default()).unwrap();
        let passing = result.passing_indices();
        // Benzene (index 0) passes; the C40 alkane (index 1) does not.
        assert!(passing.contains(&0));
        assert!(!passing.contains(&1));
    }

    #[test]
    fn verdict_carries_all_property_values() {
        let m = mol_from_smiles("CCO").unwrap();
        let v = evaluate(0, &m, &AdmetFilter::default());
        assert!(v.molecular_weight > 0.0);
        assert!((0.0..=1.0).contains(&v.qed));
        // Ethanol: 1 donor, 1 acceptor.
        assert_eq!(v.hbd, 1);
        assert_eq!(v.hba, 1);
    }

    #[test]
    fn fragment_filter_admits_only_small_fragments() {
        let benzene = mol_from_smiles("c1ccccc1").unwrap();
        let result = admet_filter(&[benzene], &AdmetFilter::fragment_like()).unwrap();
        // Benzene (78 g/mol) is well within the rule-of-three window.
        assert!(result.verdicts[0].passes);
    }
}
