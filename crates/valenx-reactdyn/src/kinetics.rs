//! Chemical-kinetics scalars — the Arrhenius rate constant and the first-order half-life.

/// Ideal gas constant `R` in J/(mol·K) (CODATA 2018).
const GAS_CONSTANT_J_PER_MOL_K: f64 = 8.314_462_618;

/// Arrhenius rate constant `k = A·exp(−Eₐ/(R·T))` for pre-exponential factor `a`, activation
/// energy `ea_j_per_mol` (J/mol), and absolute temperature `temp_k` (K). Returns `0.0` for a
/// non-positive or non-finite temperature, or any non-finite input.
pub fn arrhenius_rate(a: f64, ea_j_per_mol: f64, temp_k: f64) -> f64 {
    if !a.is_finite() || !ea_j_per_mol.is_finite() || !temp_k.is_finite() || temp_k <= 0.0 {
        return 0.0;
    }
    a * (-ea_j_per_mol / (GAS_CONSTANT_J_PER_MOL_K * temp_k)).exp()
}

/// First-order half-life `t½ = ln2 / k`, where `k` is the [`arrhenius_rate`]. Returns `0.0` when
/// the rate is non-positive or non-finite (e.g. an invalid temperature).
pub fn arrhenius_half_life_1st_order(a: f64, ea_j_per_mol: f64, temp_k: f64) -> f64 {
    let k = arrhenius_rate(a, ea_j_per_mol, temp_k);
    if !k.is_finite() || k <= 0.0 {
        return 0.0;
    }
    std::f64::consts::LN_2 / k
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn arrhenius_rate_matches_formula_and_trends() {
        // Eₐ = 0 → k = A (temperature-independent).
        assert!((arrhenius_rate(5.0e12, 0.0, 298.15) - 5.0e12).abs() < 1.0e-3);
        // k increases with T and decreases with Eₐ.
        assert!(arrhenius_rate(1e13, 50_000.0, 373.15) > arrhenius_rate(1e13, 50_000.0, 298.15));
        assert!(arrhenius_rate(1e13, 200_000.0, 298.15) < arrhenius_rate(1e13, 10_000.0, 298.15));
        // Worked: A=1e13, Eₐ=50 kJ/mol, T=300 K → k = 1e13·exp(−50000/(8.314463·300)).
        let expected = 1e13 * (-50_000.0_f64 / (8.314_462_618 * 300.0)).exp();
        assert!((arrhenius_rate(1e13, 50_000.0, 300.0) - expected).abs() < 1e-6 * expected);
        // Guard: T ≤ 0 / non-finite → 0.
        assert_eq!(arrhenius_rate(1e13, 50_000.0, 0.0), 0.0);
        assert_eq!(arrhenius_rate(1e13, 50_000.0, -10.0), 0.0);
        assert_eq!(arrhenius_rate(f64::NAN, 50_000.0, 300.0), 0.0);
    }

    #[test]
    fn arrhenius_half_life_is_ln2_over_rate() {
        let k = arrhenius_rate(1e14, 50_000.0, 300.0);
        let t_half = arrhenius_half_life_1st_order(1e14, 50_000.0, 300.0);
        assert!((t_half - std::f64::consts::LN_2 / k).abs() < 1e-15 * t_half);
        // Higher T → shorter half-life.
        assert!(
            arrhenius_half_life_1st_order(1e13, 60_000.0, 373.15)
                < arrhenius_half_life_1st_order(1e13, 60_000.0, 298.15)
        );
        // Guard: invalid T → 0 (no panic, no division by zero).
        assert_eq!(arrhenius_half_life_1st_order(1e13, 50_000.0, 0.0), 0.0);
    }
}
