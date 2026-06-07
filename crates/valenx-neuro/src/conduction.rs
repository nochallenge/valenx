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

/// The **conduction-velocity crossover diameter** `d* = (k_u/k_m)²` (µm) — the fibre
/// diameter at which a myelinated fibre (Hursh `CV = k_m·d`,
/// [`myelinated_conduction_velocity`]) and an unmyelinated one (cable-theory
/// `CV = k_u·√d`, [`unmyelinated_conduction_velocity`]) conduct at the *same* speed,
/// from the myelinated constant `k_myel_per_um` `k_m` (m/s per µm) and the
/// unmyelinated constant `k_unmyel_per_sqrt_um` `k_u` (m/s per √µm).
///
/// Setting `k_m·d = k_u·√d` and solving gives `d* = (k_u/k_m)²`. **Below** `d*` the
/// `√d` law wins and an *unmyelinated* axon is actually faster (and metabolically
/// cheaper), which is why the thinnest fibres — the bulk of a peripheral nerve — are
/// left bare; **above** `d*` myelin's linear scaling pulls ahead, the payoff that
/// makes thick fibres worth insulating. Returns `0` for a non-positive or non-finite
/// myelinated constant, or a negative / non-finite unmyelinated constant (the
/// crossover is then undefined).
pub fn myelination_crossover_diameter(k_myel_per_um: f64, k_unmyel_per_sqrt_um: f64) -> f64 {
    if !k_myel_per_um.is_finite()
        || k_myel_per_um <= 0.0
        || !k_unmyel_per_sqrt_um.is_finite()
        || k_unmyel_per_sqrt_um < 0.0
    {
        return 0.0;
    }
    (k_unmyel_per_sqrt_um / k_myel_per_um).powi(2)
}

/// The **axonal conduction delay** (latency) `t = distance / velocity` in seconds —
/// the time an action potential takes to propagate a distance `distance_m` (m) along
/// a fibre conducting at `velocity_m_per_s` (m/s, e.g. from
/// [`myelinated_conduction_velocity`]). It completes the conduction kinematic triple
/// with the velocity and the distance, and is the quantity a nerve-conduction study
/// actually measures: the latency over a known nerve segment yields the velocity.
/// Returns `0` for non-physical input (negative or non-finite distance, or a
/// non-positive / non-finite velocity, where the latency is undefined).
pub fn conduction_delay_s(distance_m: f64, velocity_m_per_s: f64) -> f64 {
    if !distance_m.is_finite()
        || !velocity_m_per_s.is_finite()
        || distance_m < 0.0
        || velocity_m_per_s <= 0.0
    {
        return 0.0;
    }
    distance_m / velocity_m_per_s
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn conduction_delay_is_distance_over_velocity() {
        // Worked point: 0.5 m at 50 m/s → 10 ms latency.
        assert!((conduction_delay_s(0.5, 50.0) - 0.01).abs() < 1e-12, "10 ms over 0.5 m @ 50 m/s");

        // Defining inverse: delay · velocity = distance.
        for &(d, v) in &[(1.0_f64, 72.0_f64), (0.3, 12.0), (2.0, 0.5)] {
            assert!((conduction_delay_s(d, v) * v - d).abs() < 1e-12, "delay·v = distance");
        }

        // Threads myelinated_conduction_velocity: a 12 µm fibre (Hursh v = 6·12 = 72 m/s)
        // has latency 1/v over a 1 m segment.
        let v = myelinated_conduction_velocity(12.0, HURSH_FACTOR_M_PER_S_PER_UM);
        assert!((conduction_delay_s(1.0, v) - 1.0 / v).abs() < 1e-12, "latency = 1/v(12 µm)");

        // Faster fibres have shorter latency.
        assert!(
            conduction_delay_s(1.0, 80.0) < conduction_delay_s(1.0, 40.0),
            "latency falls with conduction velocity"
        );

        // Non-physical input → 0.
        assert_eq!(conduction_delay_s(1.0, 0.0), 0.0); // velocity 0
        assert_eq!(conduction_delay_s(-1.0, 50.0), 0.0); // distance < 0
        assert_eq!(conduction_delay_s(1.0, f64::NAN), 0.0); // non-finite
    }

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

    #[test]
    fn myelination_crossover_diameter_equalises_the_two_conduction_velocities() {
        // k_m = 6 m/s/µm (Hursh), k_u = 2 m/s/√µm → d* = (k_u/k_m)² = (1/3)² = 1/9 µm.
        let d_star = myelination_crossover_diameter(6.0, 2.0);
        assert!((d_star - 1.0 / 9.0).abs() < 1e-12, "d* = (k_u/k_m)² = 1/9, got {d_star}");
        // STRONG dual cross-check: at d* the two CV laws agree exactly (threading BOTH
        // conduction-velocity functions); below d* the √d law is faster, above it the
        // linear (myelinated) law wins.
        for &(km, ku) in &[(6.0_f64, 2.0_f64), (5.0, 3.0), (10.0, 1.5), (4.0, 4.0)] {
            let d = myelination_crossover_diameter(km, ku);
            let cv_m = myelinated_conduction_velocity(d, km);
            let cv_u = unmyelinated_conduction_velocity(d, ku);
            assert!((cv_m - cv_u).abs() / cv_m < 1e-9, "CVs equal at d* for k_m={km}, k_u={ku}");
            assert!(
                unmyelinated_conduction_velocity(0.5 * d, ku) > myelinated_conduction_velocity(0.5 * d, km),
                "below d*: unmyelinated faster (k_m={km}, k_u={ku})"
            );
            assert!(
                myelinated_conduction_velocity(2.0 * d, km) > unmyelinated_conduction_velocity(2.0 * d, ku),
                "above d*: myelinated faster (k_m={km}, k_u={ku})"
            );
        }
        // Non-physical input → 0 (the crossover is undefined).
        assert_eq!(myelination_crossover_diameter(0.0, 2.0), 0.0); // k_m ≤ 0
        assert_eq!(myelination_crossover_diameter(-6.0, 2.0), 0.0); // k_m < 0
        assert_eq!(myelination_crossover_diameter(f64::NAN, 2.0), 0.0); // non-finite k_m
        assert_eq!(myelination_crossover_diameter(6.0, f64::INFINITY), 0.0); // non-finite k_u
        assert_eq!(myelination_crossover_diameter(6.0, -1.0), 0.0); // k_u < 0
    }
}
