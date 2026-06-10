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

/// Thermodynamic equilibrium constant `K = exp(−ΔG°/(R·T))` from the standard Gibbs free-energy
/// change `delta_g_j_per_mol` (J/mol) and absolute temperature `temp_k` (K) — the rearrangement of
/// `ΔG° = −R·T·ln K`. Distinct from [`arrhenius_rate`] (a kinetic rate constant). Returns `0.0` for
/// a non-positive or non-finite temperature, or non-finite ΔG.
pub fn equilibrium_constant(delta_g_j_per_mol: f64, temp_k: f64) -> f64 {
    if !delta_g_j_per_mol.is_finite() || !temp_k.is_finite() || temp_k <= 0.0 {
        return 0.0;
    }
    (-delta_g_j_per_mol / (GAS_CONSTANT_J_PER_MOL_K * temp_k)).exp()
}

/// Standard Gibbs free-energy change `ΔG° = −R·T·ln K` (J/mol) recovered from a thermodynamic
/// equilibrium constant `k_eq` and absolute temperature `temp_k` (K) — the inverse of
/// [`equilibrium_constant`], for the common case where `K` is measured/known and `ΔG°` is wanted
/// (e.g. thermodynamic feasibility bounds in flux-balance analysis). Returns `0.0` for a
/// non-positive or non-finite `k_eq` (where `ln K` is undefined) or temperature.
pub fn gibbs_free_energy_from_equilibrium_constant(k_eq: f64, temp_k: f64) -> f64 {
    if !k_eq.is_finite() || k_eq <= 0.0 || !temp_k.is_finite() || temp_k <= 0.0 {
        return 0.0;
    }
    -GAS_CONSTANT_J_PER_MOL_K * temp_k * k_eq.ln()
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

    #[test]
    fn equilibrium_constant_from_gibbs_energy() {
        // ΔG = 0 → K = 1.
        assert!((equilibrium_constant(0.0, 298.15) - 1.0).abs() < 1e-12);
        // ΔG = −10000 J/mol @ 300 K → K ≈ 55.08 (exergonic, K > 1).
        let expected = (10_000.0_f64 / (8.314_462_618 * 300.0)).exp();
        assert!((equilibrium_constant(-10_000.0, 300.0) - expected).abs() < 1e-2 * expected);
        assert!(equilibrium_constant(-10_000.0, 300.0) > 1.0);
        assert!(equilibrium_constant(10_000.0, 300.0) < 1.0);
        // Inverse: K(ΔG)·K(−ΔG) = 1.
        assert!(
            (equilibrium_constant(-15_000.0, 298.15) * equilibrium_constant(15_000.0, 298.15)
                - 1.0)
                .abs()
                < 1e-12
        );
        // Guard: T ≤ 0 / non-finite → 0.
        assert_eq!(equilibrium_constant(-10_000.0, 0.0), 0.0);
        assert_eq!(equilibrium_constant(f64::NAN, 300.0), 0.0);
    }

    #[test]
    fn gibbs_from_keq_inverts_equilibrium_constant() {
        // Round-trip ΔG → K → ΔG recovers the original.
        let dg = -12_345.0;
        let t = 310.0;
        let k = equilibrium_constant(dg, t);
        let back = gibbs_free_energy_from_equilibrium_constant(k, t);
        assert!((back - dg).abs() < 1e-9 * dg.abs());
        // K = 1 ⇒ ΔG° = 0; K > 1 (exergonic) ⇒ ΔG° < 0; K < 1 ⇒ ΔG° > 0.
        assert!(gibbs_free_energy_from_equilibrium_constant(1.0, 298.15).abs() < 1e-9);
        assert!(gibbs_free_energy_from_equilibrium_constant(50.0, 298.15) < 0.0);
        assert!(gibbs_free_energy_from_equilibrium_constant(0.02, 298.15) > 0.0);
        // Worked: K=10, T=298.15 → ΔG = −R·T·ln10 ≈ −5708 J/mol.
        let expected = -8.314_462_618 * 298.15 * 10.0_f64.ln();
        assert!(
            (gibbs_free_energy_from_equilibrium_constant(10.0, 298.15) - expected).abs() < 1e-6
        );
        // Guard: non-positive K / T, or non-finite → 0 (no panic, no ln of ≤ 0).
        assert_eq!(gibbs_free_energy_from_equilibrium_constant(0.0, 298.15), 0.0);
        assert_eq!(gibbs_free_energy_from_equilibrium_constant(-1.0, 298.15), 0.0);
        assert_eq!(gibbs_free_energy_from_equilibrium_constant(10.0, 0.0), 0.0);
        assert_eq!(gibbs_free_energy_from_equilibrium_constant(f64::NAN, 298.15), 0.0);
    }
}
