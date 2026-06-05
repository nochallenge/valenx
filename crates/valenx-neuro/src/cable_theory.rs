//! Passive-cable constants (Rall cable theory).
//!
//! Before an axon spikes, it behaves as a leaky RC cable. Two constants set
//! how a subthreshold signal spreads and settles:
//!
//! - **Length (space) constant** `λ = sqrt(d·R_m / (4·R_i))` — the distance
//!   over which a steady voltage decays to `1/e`. Larger fibres, higher
//!   membrane resistance, and lower axoplasmic resistivity all let signals
//!   reach farther.
//! - **Membrane time constant** `τ = R_m·C_m` — how fast the membrane charges.
//!
//! and the dimensionless **electrotonic length** `L = ℓ/λ` rates a real
//! segment's length against `λ` (a compartment with `L ≪ 1` is electrically
//! compact; `L ≫ 1` is "long").
//!
//! Units are the cellular-physiology CGS-ish convention: `d` in cm, specific
//! membrane resistance `R_m` in Ω·cm², axial resistivity `R_i` in Ω·cm,
//! specific membrane capacitance `C_m` in F/cm² — giving `λ` in cm and `τ` in
//! seconds. These are the exact passive results, not a fit.

/// Length (space) constant `λ = sqrt(d·R_m / (4·R_i))` in cm, from fibre
/// diameter `d_cm` (cm), specific membrane resistance `r_m_ohm_cm2` (Ω·cm²),
/// and axial resistivity `r_i_ohm_cm` (Ω·cm).
pub fn length_constant_cm(d_cm: f64, r_m_ohm_cm2: f64, r_i_ohm_cm: f64) -> f64 {
    (d_cm * r_m_ohm_cm2 / (4.0 * r_i_ohm_cm)).sqrt()
}

/// Membrane time constant `τ = R_m·C_m` in seconds, from specific membrane
/// resistance `r_m_ohm_cm2` (Ω·cm²) and capacitance `c_m_f_cm2` (F/cm²).
pub fn time_constant_s(r_m_ohm_cm2: f64, c_m_f_cm2: f64) -> f64 {
    r_m_ohm_cm2 * c_m_f_cm2
}

/// Dimensionless electrotonic length `L = ℓ / λ` for a physical length
/// `length_cm` and a length constant `lambda_cm` (both cm). `lambda_cm` must
/// be positive.
pub fn electrotonic_length(length_cm: f64, lambda_cm: f64) -> f64 {
    length_cm / lambda_cm
}

/// Semi-infinite cable input resistance `R_∞ = √(r_m·r_a)` in ohms — the
/// resistance seen looking into one end of a long cylindrical fibre, from fibre
/// diameter `d_cm` (cm), specific membrane resistance `r_m_ohm_cm2` (Ω·cm²) and
/// axial resistivity `r_i_ohm_cm` (Ω·cm). Equivalently `R_∞ = r_a·λ`, the axial
/// resistance per unit length `r_a = 4·R_i/(π·d²)` times the length constant λ;
/// it falls as `d^(−3/2)`, so a thicker fibre presents a much lower input load.
pub fn semi_infinite_input_resistance(d_cm: f64, r_m_ohm_cm2: f64, r_i_ohm_cm: f64) -> f64 {
    // Per-unit-length membrane (Ω·cm) and axial (Ω/cm) resistances.
    let r_m = r_m_ohm_cm2 / (std::f64::consts::PI * d_cm);
    let r_a = 4.0 * r_i_ohm_cm / (std::f64::consts::PI * d_cm * d_cm);
    (r_m * r_a).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn length_constant_matches_the_formula_and_scales() {
        // λ = sqrt(d·R_m/(4·R_i)); d=4e-4, R_m=1e4, R_i=100 → sqrt(0.01) = 0.1 cm.
        assert!((length_constant_cm(4e-4, 1e4, 100.0) - 0.1).abs() < 1e-12);
        // λ ∝ √d: quadrupling the diameter doubles λ.
        let l1 = length_constant_cm(4e-4, 1e4, 100.0);
        let l4 = length_constant_cm(16e-4, 1e4, 100.0);
        assert!((l4 - 2.0 * l1).abs() < 1e-12);
        // λ ∝ 1/√(R_i): quadrupling axial resistivity halves λ.
        let lr = length_constant_cm(4e-4, 1e4, 400.0);
        assert!((lr - 0.5 * l1).abs() < 1e-12);
    }

    #[test]
    fn time_constant_is_rm_times_cm() {
        // τ = R_m·C_m; 1e4 Ω·cm² × 1e-6 F/cm² = 0.01 s (10 ms).
        assert!((time_constant_s(1e4, 1e-6) - 0.01).abs() < 1e-15);
        // Linear in each factor.
        assert!((time_constant_s(2e4, 1e-6) - 2.0 * time_constant_s(1e4, 1e-6)).abs() < 1e-15);
    }

    #[test]
    fn electrotonic_length_is_length_over_lambda() {
        assert!((electrotonic_length(0.2, 0.1) - 2.0).abs() < 1e-12);
        // A compartment much shorter than λ is electrically compact (L ≪ 1).
        assert!(electrotonic_length(0.01, 0.1) < 0.2);
    }

    #[test]
    fn semi_infinite_input_resistance_relates_to_length_constant() {
        // d=4e-4 cm (4 µm), R_m=1e4 Ω·cm², R_i=100 Ω·cm.
        let r_inf = semi_infinite_input_resistance(4e-4, 1e4, 100.0);
        assert!(r_inf > 0.0 && r_inf.is_finite());
        // R_∞ = r_a·λ — cross-check against the length constant.
        let lambda = length_constant_cm(4e-4, 1e4, 100.0);
        let r_a = 4.0 * 100.0 / (std::f64::consts::PI * 4e-4 * 4e-4);
        assert!((r_inf - r_a * lambda).abs() / r_inf < 1e-9, "R_∞ = r_a·λ: {r_inf}");
        // ~80 MΩ for this thin fibre.
        assert!((r_inf / 1e6 - 79.6).abs() < 1.0, "R_∞ {} MΩ", r_inf / 1e6);
        // R_∞ ∝ d^(−3/2): a 4× thicker fibre has 1/8 the input resistance.
        let thick = semi_infinite_input_resistance(16e-4, 1e4, 100.0);
        assert!((thick - r_inf / 8.0).abs() / r_inf < 1e-9, "d^(−3/2) scaling");
    }
}
