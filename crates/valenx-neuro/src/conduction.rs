//! Action-potential conduction velocity — closed-form estimates of how fast a
//! spike travels along an axon, as a function of fibre diameter.
//!
//! - **Myelinated** fibres obey Hursh's rule: velocity scales *linearly* with
//!   the outer (myelin-included) diameter, `CV ≈ k·d`, with `k ≈ 6` m/s per µm
//!   for mammalian fibres (Hursh 1939). Saltatory conduction — the spike
//!   jumping node to node — is what buys this steep diameter dependence.
//! - **Unmyelinated** fibres obey the cable-theory scaling `CV ∝ √d`: thicker
//!   axons have lower axial resistance, but only as the square root, so growing
//!   an unmyelinated fibre is a far less efficient way to go fast (Hodgkin 1954).
//!
//! These are scaling laws with empirical constants, not a compartmental solve —
//! the [`crate::myelinated`] module integrates the full saltatory cable when a
//! mechanistic velocity is needed. The proportionality constants are left as
//! parameters because they depend on species, temperature, and membrane
//! properties; [`HURSH_FACTOR_M_PER_S_PER_UM`] is the standard myelinated value.

/// Hursh's proportionality constant `k` for myelinated mammalian fibres:
/// ≈ 6 m/s of conduction velocity per micrometre of outer fibre diameter.
pub const HURSH_FACTOR_M_PER_S_PER_UM: f64 = 6.0;

/// Conduction velocity (m/s) of a **myelinated** fibre by Hursh's rule,
/// `CV = k·d`, for outer fibre diameter `diameter_um` (µm) and proportionality
/// `k_m_per_s_per_um` (m/s per µm; the classic value is
/// [`HURSH_FACTOR_M_PER_S_PER_UM`]).
pub fn myelinated_conduction_velocity(diameter_um: f64, k_m_per_s_per_um: f64) -> f64 {
    k_m_per_s_per_um * diameter_um
}

/// Conduction velocity (m/s) of an **unmyelinated** fibre by the cable-theory
/// scaling `CV = k·√d`, for axon diameter `diameter_um` (µm) and scale constant
/// `k_m_per_s_per_sqrt_um` (m/s per √µm, set by membrane and temperature).
pub fn unmyelinated_conduction_velocity(diameter_um: f64, k_m_per_s_per_sqrt_um: f64) -> f64 {
    k_m_per_s_per_sqrt_um * diameter_um.max(0.0).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn myelinated_follows_hursh_six_times_diameter() {
        // A 10 µm mammalian fibre conducts at ~60 m/s (Hursh's rule).
        let cv = myelinated_conduction_velocity(10.0, HURSH_FACTOR_M_PER_S_PER_UM);
        assert!((cv - 60.0).abs() < 1e-9, "10 µm fibre should be ~60 m/s, got {cv}");
        // Linear in diameter: doubling d doubles the velocity.
        assert!(
            (myelinated_conduction_velocity(20.0, 6.0)
                - 2.0 * myelinated_conduction_velocity(10.0, 6.0))
            .abs()
                < 1e-9
        );
    }

    #[test]
    fn unmyelinated_scales_with_sqrt_diameter() {
        // CV ∝ √d: quadrupling the diameter only doubles the velocity.
        let cv1 = unmyelinated_conduction_velocity(1.0, 2.0);
        let cv4 = unmyelinated_conduction_velocity(4.0, 2.0);
        assert!((cv4 - 2.0 * cv1).abs() < 1e-9, "√d scaling: {cv1} -> {cv4}");
        // At unit diameter CV = k.
        assert!((unmyelinated_conduction_velocity(1.0, 1.8) - 1.8).abs() < 1e-12);
        // Non-positive diameter is clamped to a zero velocity, not a NaN.
        assert_eq!(unmyelinated_conduction_velocity(-1.0, 2.0), 0.0);
    }

    #[test]
    fn myelination_is_faster_for_a_typical_fibre() {
        // A 10 µm myelinated fibre vastly outruns a 1 µm unmyelinated one —
        // the evolutionary payoff of myelin.
        let myelinated = myelinated_conduction_velocity(10.0, HURSH_FACTOR_M_PER_S_PER_UM);
        let unmyelinated = unmyelinated_conduction_velocity(1.0, 2.0);
        assert!(myelinated > 10.0 * unmyelinated, "{myelinated} vs {unmyelinated}");
    }
}
