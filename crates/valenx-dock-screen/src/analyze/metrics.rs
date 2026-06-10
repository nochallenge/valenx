//! Post-docking efficiency metrics.

/// Ligand efficiency — the binding affinity normalised by ligand size:
/// `LE = −binding_score / heavy_atoms` (kcal·mol⁻¹ per heavy atom). Docking scores are
/// negative for favourable binding, so the negation makes a stronger binder give a larger
/// (more positive) LE; typical lead-like values are ~0.3–0.5. Returns `0.0` when there are
/// no heavy atoms or the score is non-finite.
pub fn ligand_efficiency(binding_score: f64, heavy_atoms: usize) -> f64 {
    if heavy_atoms == 0 || !binding_score.is_finite() {
        return 0.0;
    }
    -binding_score / heavy_atoms as f64
}

/// Lipophilic efficiency (LLE) — binding potency penalised by lipophilicity:
/// `LLE = pIC50 − logP` (dimensionless), rewarding high potency with low lipophilicity (a key
/// ADMET balance). `potency_pic50` is the negative log of the binding affinity (molar); `logp`
/// is the partition coefficient. Typical lead-like values are 5–7. Returns `0.0` when either
/// input is non-finite. Distinct from [`ligand_efficiency`] (which normalises by ligand size).
pub fn lipophilic_efficiency(potency_pic50: f64, logp: f64) -> f64 {
    if !potency_pic50.is_finite() || !logp.is_finite() {
        return 0.0;
    }
    potency_pic50 - logp
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ligand_efficiency_basic() {
        // −9.0 kcal/mol over 30 heavy atoms → 0.3; −4.5 over 9 → 0.5.
        assert!((ligand_efficiency(-9.0, 30) - 0.3).abs() < 1e-9);
        assert!((ligand_efficiency(-4.5, 9) - 0.5).abs() < 1e-9);
        // A favourable (negative) score yields a positive LE.
        assert!(ligand_efficiency(-6.0, 12) > 0.0);
        // Inverse scaling: at a fixed score, doubling the heavy-atom count halves LE.
        assert!(
            (ligand_efficiency(-10.0, 10) - 2.0 * ligand_efficiency(-10.0, 20)).abs() < 1e-9
        );
        // Guards: zero heavy atoms or a non-finite score → 0.0.
        assert_eq!(ligand_efficiency(-9.0, 0), 0.0);
        assert_eq!(ligand_efficiency(f64::NAN, 10), 0.0);
        assert_eq!(ligand_efficiency(f64::INFINITY, 10), 0.0);
    }

    #[test]
    fn lipophilic_efficiency_basic() {
        // LLE = pIC50 − logP. (8,3) → 5; (7,2) → 5.
        assert!((lipophilic_efficiency(8.0, 3.0) - 5.0).abs() < 1e-9);
        assert!((lipophilic_efficiency(7.0, 2.0) - 5.0).abs() < 1e-9);
        // +1 potency raises LLE by 1; +1 logP lowers it by 1.
        assert!(
            (lipophilic_efficiency(9.0, 3.0) - lipophilic_efficiency(8.0, 3.0) - 1.0).abs() < 1e-9
        );
        assert!(
            (lipophilic_efficiency(8.0, 4.0) - lipophilic_efficiency(8.0, 3.0) + 1.0).abs() < 1e-9
        );
        // Guards: non-finite inputs → 0.0.
        assert_eq!(lipophilic_efficiency(f64::NAN, 3.0), 0.0);
        assert_eq!(lipophilic_efficiency(8.0, f64::INFINITY), 0.0);
    }
}
