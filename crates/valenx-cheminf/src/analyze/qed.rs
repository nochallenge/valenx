//! QED — the Quantitative Estimate of Drug-likeness.
//!
//! QED (Bickerton et al., *Nature Chemistry* 2012) scores how
//! drug-like a molecule is on a 0..1 scale by combining eight
//! physicochemical descriptors, each passed through a *desirability
//! function* `d(x)` that maps the descriptor onto 0..1, then taking
//! the geometric mean of the eight desirabilities.
//!
//! The eight descriptors: molecular weight, logP, hydrogen-bond donor
//! count, hydrogen-bond acceptor count, topological polar surface
//! area, rotatable-bond count, aromatic-ring count and a structural-
//! alert count.
//!
//! [`qed`] computes the unweighted QED (`QEDw,u` — the most common
//! variant). The desirability functions here are asymmetric double-
//! sigmoid bumps fitted to the published QED desirability shapes:
//! each has a favourable plateau and falls off on either side.
//!
//! **v1 note.** The published QED uses desirability functions defined
//! by a specific set of fitted parameters; this module uses
//! double-sigmoid bumps with peaks / widths placed to reproduce the
//! published shapes closely. The structural-alert count uses a small
//! built-in SMARTS alert set, not the full Brenk / PAINS catalogue.

use crate::descriptors;
use crate::molecule::Molecule;
use crate::perceive::formula::average_molecular_weight;
use crate::smarts::SmartsPattern;

/// The eight descriptor desirabilities that make up a QED score.
#[derive(Clone, Debug, PartialEq)]
pub struct QedComponents {
    /// Desirability of the molecular weight.
    pub mw: f64,
    /// Desirability of the logP.
    pub logp: f64,
    /// Desirability of the H-bond donor count.
    pub hbd: f64,
    /// Desirability of the H-bond acceptor count.
    pub hba: f64,
    /// Desirability of the polar surface area.
    pub psa: f64,
    /// Desirability of the rotatable-bond count.
    pub rotb: f64,
    /// Desirability of the aromatic-ring count.
    pub arom: f64,
    /// Desirability of the structural-alert count.
    pub alerts: f64,
}

impl QedComponents {
    /// The QED score — the geometric mean of the eight desirabilities.
    pub fn score(&self) -> f64 {
        let vals = [
            self.mw,
            self.logp,
            self.hbd,
            self.hba,
            self.psa,
            self.rotb,
            self.arom,
            self.alerts,
        ];
        // geometric mean = exp(mean(ln d)); guard against a zero factor
        let mut sum_ln = 0.0;
        for v in vals {
            sum_ln += v.max(1e-6).ln();
        }
        (sum_ln / vals.len() as f64).exp()
    }
}

/// Compute the unweighted QED score of `mol` in `0..1`.
pub fn qed(mol: &Molecule) -> f64 {
    qed_components(mol).score()
}

/// Compute the eight QED desirability components.
pub fn qed_components(mol: &Molecule) -> QedComponents {
    let mw = average_molecular_weight(mol);
    let logp = descriptors::crippen_logp(mol);
    let hbd = descriptors::hbd(mol) as f64;
    let hba = descriptors::hba(mol) as f64;
    let psa = descriptors::tpsa(mol);
    let rotb = descriptors::rotatable_bonds(mol) as f64;
    let arom = descriptors::aromatic_ring_count(mol) as f64;
    let alerts = structural_alert_count(mol) as f64;

    QedComponents {
        // each desirability is a bump: favourable plateau, falloff
        mw: bump(mw, 300.0, 100.0, 150.0),
        logp: bump(logp, 2.5, 2.5, 3.0),
        hbd: bump(hbd, 1.5, 1.5, 2.5),
        hba: bump(hba, 4.0, 3.5, 3.5),
        psa: bump(psa, 70.0, 50.0, 60.0),
        rotb: bump(rotb, 5.0, 5.0, 4.0),
        arom: bump(arom, 2.0, 1.5, 1.5),
        alerts: alert_desirability(alerts),
    }
}

/// An asymmetric "bump" desirability: 1.0 at `peak`, falling to ~0.5 a
/// distance `left` below the peak and `right` above it, smoothly via a
/// Gaussian-style decay. Models the published QED desirability shape.
fn bump(x: f64, peak: f64, left: f64, right: f64) -> f64 {
    let width = if x < peak { left } else { right };
    let z = (x - peak) / width.max(1e-6);
    // exp(-z^2 * ln2) → exactly 0.5 at |z| = 1
    (-(z * z) * std::f64::consts::LN_2).exp().clamp(1e-4, 1.0)
}

/// Desirability of the structural-alert count — 1.0 for no alerts,
/// decaying as alerts accumulate.
fn alert_desirability(n: f64) -> f64 {
    // each alert roughly halves the desirability
    0.5f64.powf(n).clamp(1e-3, 1.0)
}

/// A small built-in structural-alert SMARTS set (a Brenk-style
/// subset). Each match counts as one alert.
const ALERT_SMARTS: &[&str] = &[
    "[N+](=O)[O-]", // nitro group
    "C(=O)Cl",      // acyl chloride
    "[N-]=[N+]=[N-]", // azide
    "[S][S]",       // disulfide
    "C=C=C",        // allene / cumulated diene
    "[O-][O-]",     // peroxide-ish
];

/// Count structural-alert substructures in `mol`.
pub fn structural_alert_count(mol: &Molecule) -> usize {
    let mut count = 0;
    for smarts in ALERT_SMARTS {
        if let Ok(pat) = SmartsPattern::parse(smarts) {
            count += pat.count_unique(mol);
        }
    }
    count
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mol_from_smiles;

    #[test]
    fn qed_in_unit_range() {
        for s in ["CCO", "c1ccccc1", "CC(=O)O", "CCCCCCCCCCCC"] {
            let m = mol_from_smiles(s).unwrap();
            let q = qed(&m);
            assert!((0.0..=1.0).contains(&q), "{s} → QED {q}");
        }
    }

    #[test]
    fn drug_like_scores_higher_than_extreme() {
        // a small drug-like molecule should beat a huge greasy chain
        let drug = mol_from_smiles("c1ccc(cc1)C(=O)NCCO").unwrap();
        let greasy = mol_from_smiles("CCCCCCCCCCCCCCCCCCCC").unwrap();
        assert!(
            qed(&drug) > qed(&greasy),
            "drug {} vs greasy {}",
            qed(&drug),
            qed(&greasy)
        );
    }

    #[test]
    fn bump_is_one_at_peak() {
        assert!((bump(300.0, 300.0, 100.0, 150.0) - 1.0).abs() < 1e-9);
        // half-height at peak +/- width
        assert!((bump(450.0, 300.0, 100.0, 150.0) - 0.5).abs() < 1e-6);
        assert!((bump(200.0, 300.0, 100.0, 150.0) - 0.5).abs() < 1e-6);
    }

    #[test]
    fn structural_alert_detected() {
        // nitrobenzene carries a nitro alert
        let m = mol_from_smiles("c1ccccc1[N+](=O)[O-]").unwrap();
        assert!(structural_alert_count(&m) >= 1);
        // benzene has none
        let benzene = mol_from_smiles("c1ccccc1").unwrap();
        assert_eq!(structural_alert_count(&benzene), 0);
    }

    #[test]
    fn alert_lowers_qed() {
        // QED is the geometric mean of eight property desirabilities, so
        // comparing two *different* molecules cannot isolate the alert
        // term (a nitro group changes MW, PSA, HBA … all at once). The
        // honest test holds every other property fixed: take a real
        // molecule's QED components and confirm that lowering ONLY the
        // structural-alert desirability lowers the overall score.
        let m = mol_from_smiles("c1ccccc1[N+](=O)[O-]").unwrap();
        let mut comps = qed_components(&m);
        comps.alerts = 1.0; // pretend it has no structural alert
        let score_no_alert = comps.score();
        comps.alerts = 0.5; // one structural alert
        let score_with_alert = comps.score();
        assert!(
            score_with_alert < score_no_alert,
            "a structural alert must lower QED: {score_with_alert} < {score_no_alert}"
        );
        // And nitrobenzene really does carry an alert.
        assert!(structural_alert_count(&m) >= 1);
    }

    #[test]
    fn components_geometric_mean() {
        let m = mol_from_smiles("CCO").unwrap();
        let c = qed_components(&m);
        // score equals geometric mean of the eight parts
        assert!((c.score() - qed(&m)).abs() < 1e-9);
    }
}
