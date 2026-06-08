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

/// The **myelinated fibre diameter from a measured conduction velocity** `d = CV/k` (µm) —
/// the inverse of Hursh's rule [`myelinated_conduction_velocity`] (`CV = k·d`), recovering
/// outer fibre diameter from a velocity `velocity_m_per_s` (m/s) and the proportionality
/// `k_m_per_s_per_um` (m/s per µm; classic [`HURSH_FACTOR_M_PER_S_PER_UM`]). This is the
/// standard nerve-conduction-study size estimate (a 60 m/s response ≈ a 10 µm fibre).
/// Returns `0` for non-physical input (negative or non-finite velocity, or non-positive /
/// non-finite `k`).
pub fn myelinated_fiber_diameter_um(velocity_m_per_s: f64, k_m_per_s_per_um: f64) -> f64 {
    if !velocity_m_per_s.is_finite()
        || velocity_m_per_s < 0.0
        || !k_m_per_s_per_um.is_finite()
        || k_m_per_s_per_um <= 0.0
    {
        return 0.0;
    }
    velocity_m_per_s / k_m_per_s_per_um
}

/// The **unmyelinated fibre diameter from a measured conduction velocity** `d = (CV/k)²`
/// (µm) — the inverse of the cable-theory `CV = k·√d` law
/// ([`unmyelinated_conduction_velocity`]), for a velocity `velocity_m_per_s` (m/s) and the
/// scale constant `k_m_per_s_per_sqrt_um` (m/s per √µm). Returns `0` for non-physical input
/// (negative or non-finite velocity, or non-positive / non-finite `k`).
pub fn unmyelinated_fiber_diameter_um(velocity_m_per_s: f64, k_m_per_s_per_sqrt_um: f64) -> f64 {
    if !velocity_m_per_s.is_finite()
        || velocity_m_per_s < 0.0
        || !k_m_per_s_per_sqrt_um.is_finite()
        || k_m_per_s_per_sqrt_um <= 0.0
    {
        return 0.0;
    }
    (velocity_m_per_s / k_m_per_s_per_sqrt_um).powi(2)
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

/// The **conduction propagation distance** `distance = velocity · time` in metres — how
/// far an action potential travels in time `time_s` (s) along a fibre conducting at
/// `velocity_m_per_s` (m/s, e.g. from [`myelinated_conduction_velocity`]). It is the
/// inverse of the [`conduction_delay_s`] latency and the third member of the conduction
/// kinematic triple (velocity, distance, delay). Returns `0` for non-physical input
/// (negative or non-finite velocity or time).
pub fn conduction_distance_m(velocity_m_per_s: f64, time_s: f64) -> f64 {
    if !velocity_m_per_s.is_finite() || !time_s.is_finite() || velocity_m_per_s < 0.0 || time_s < 0.0
    {
        return 0.0;
    }
    velocity_m_per_s * time_s
}

/// The **conduction velocity from a measured latency** `velocity = distance / time`
/// (m/s) — the speed inferred when an action potential covers a known distance
/// `distance_m` (m) in time `time_s` (s). This is exactly the clinical nerve-conduction
/// study measurement (a latency over a known nerve segment yields the velocity), and the
/// third permutation of the conduction kinematic triple — the inverse of both
/// [`conduction_delay_s`] (in the velocity) and [`conduction_distance_m`]. Distinct from
/// the diameter-based [`myelinated_conduction_velocity`]. Returns `0` for non-physical
/// input (negative or non-finite distance, or a non-positive / non-finite time, where the
/// velocity is undefined).
pub fn conduction_velocity_from_latency(distance_m: f64, time_s: f64) -> f64 {
    if !distance_m.is_finite() || !time_s.is_finite() || distance_m < 0.0 || time_s <= 0.0 {
        return 0.0;
    }
    distance_m / time_s
}

/// The **conduction velocity at the myelination crossover diameter** (m/s) — the
/// single speed at which a myelinated fibre (`v = k_m·d`) and an unmyelinated fibre
/// (`v = k_u·√d`) conduct equally. Below it, myelination gives no velocity benefit
/// (it is why the thinnest fibres are left unmyelinated); ≈ 1–2 m/s in mammals.
/// `k_myel_per_um` `k_m` and `k_unmyel_per_sqrt_um` `k_u` are the two velocity
/// constants (see [`myelinated_conduction_velocity`] /
/// [`unmyelinated_conduction_velocity`]). It is the myelinated velocity evaluated at
/// the [`myelination_crossover_diameter`], analytically `k_u²/k_m`. Returns `0` for
/// the same non-physical inputs as [`myelination_crossover_diameter`].
pub fn conduction_velocity_at_crossover(k_myel_per_um: f64, k_unmyel_per_sqrt_um: f64) -> f64 {
    let d = myelination_crossover_diameter(k_myel_per_um, k_unmyel_per_sqrt_um);
    myelinated_conduction_velocity(d, k_myel_per_um)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn conduction_velocity_from_latency_inverts_delay_and_distance() {
        // Double round-trip: recover v from the latency it produces (inverts
        // conduction_delay_s) and from the distance it covers (inverts conduction_distance_m).
        for &(d, v, t) in &[(0.05_f64, 50.0_f64, 0.001_f64), (1.0, 120.0, 0.05), (0.3, 2.5, 0.2)] {
            assert!(
                (conduction_velocity_from_latency(d, conduction_delay_s(d, v)) - v).abs() <= 1e-12 * v,
                "v = d/t inverts t = d/v"
            );
            assert!(
                (conduction_velocity_from_latency(conduction_distance_m(v, t), t) - v).abs()
                    <= 1e-12 * v,
                "v = d/t inverts d = v·t"
            );
        }

        // Worked: a 5 cm segment crossed in 1 ms → 50 m/s (normal myelinated NCV).
        assert!(
            (conduction_velocity_from_latency(0.05, 0.001) - 50.0).abs() < 1e-12,
            "5 cm / 1 ms = 50 m/s"
        );

        // Linear in distance, inverse in time.
        assert!(
            (conduction_velocity_from_latency(0.10, 0.001)
                - 2.0 * conduction_velocity_from_latency(0.05, 0.001))
            .abs()
                < 1e-12,
            "linear in distance"
        );
        assert!(
            (conduction_velocity_from_latency(0.05, 0.002)
                - 0.5 * conduction_velocity_from_latency(0.05, 0.001))
            .abs()
                < 1e-12,
            "inverse in time"
        );

        // Zero distance → zero velocity (valid); non-physical (time ≤ 0, negatives, NaN) → 0.
        assert_eq!(conduction_velocity_from_latency(0.0, 0.001), 0.0);
        assert_eq!(conduction_velocity_from_latency(0.05, 0.0), 0.0);
        assert_eq!(conduction_velocity_from_latency(-0.05, 0.001), 0.0);
        assert_eq!(conduction_velocity_from_latency(0.05, -0.001), 0.0);
        assert_eq!(conduction_velocity_from_latency(f64::NAN, 0.001), 0.0);
    }

    #[test]
    fn conduction_distance_inverts_the_delay() {
        // Round-trip: distance recovered from the latency it produces (the exact inverse
        // of conduction_delay_s).
        for &(d, v) in &[(0.05_f64, 50.0_f64), (1.0, 120.0), (0.3, 2.5)] {
            let recovered = conduction_distance_m(v, conduction_delay_s(d, v));
            assert!((recovered - d).abs() <= 1e-12 * d, "distance = v·t inverts t = d/v");
        }

        // Worked: a 50 m/s myelinated fibre travels 5 cm in 1 ms.
        assert!((conduction_distance_m(50.0, 0.001) - 0.05).abs() < 1e-12, "50 m/s × 1 ms = 5 cm");

        // Linear in velocity and time.
        assert!(
            (conduction_distance_m(100.0, 0.002) - 2.0 * conduction_distance_m(50.0, 0.002)).abs()
                < 1e-12,
            "linear in v"
        );
        assert!(
            (conduction_distance_m(50.0, 0.004) - 2.0 * conduction_distance_m(50.0, 0.002)).abs()
                < 1e-12,
            "linear in t"
        );

        // Zero velocity or time → zero distance (valid); negatives / NaN → 0 sentinel.
        assert_eq!(conduction_distance_m(0.0, 0.001), 0.0);
        assert_eq!(conduction_distance_m(50.0, 0.0), 0.0);
        assert_eq!(conduction_distance_m(-50.0, 0.001), 0.0);
        assert_eq!(conduction_distance_m(50.0, -0.001), 0.0);
        assert_eq!(conduction_distance_m(f64::NAN, 0.001), 0.0);
    }

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
    fn conduction_velocity_at_crossover_is_where_the_curves_meet() {
        // At the crossover diameter the myelinated and unmyelinated velocity curves
        // intersect, so the crossover velocity equals BOTH (the impl takes the
        // myelinated branch; this confirms the unmyelinated branch agrees) and the
        // closed form k_u²/k_m.
        for &(km, ku) in &[(6.0_f64, 4.0_f64), (5.0, 3.0), (6.0, 2.0)] {
            let v = conduction_velocity_at_crossover(km, ku);
            let d = myelination_crossover_diameter(km, ku);
            assert!(
                (v - unmyelinated_conduction_velocity(d, ku)).abs() / v < 1e-12,
                "crossover velocity = unmyelinated v at d*"
            );
            assert!((v - ku * ku / km).abs() / v < 1e-12, "= k_u²/k_m");
            assert!(
                (v - myelinated_conduction_velocity(d, km)).abs() / v < 1e-12,
                "= myelinated v at d*"
            );
            assert!(v > 0.0, "positive crossover velocity");
        }
        // Non-physical k_m → 0.
        assert_eq!(conduction_velocity_at_crossover(0.0, 4.0), 0.0);
        assert_eq!(conduction_velocity_at_crossover(-1.0, 4.0), 0.0);
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

    #[test]
    fn fiber_diameter_from_conduction_velocity_inverts_the_velocity_laws() {
        // (a) WORKED myelinated: v = 60 m/s at Hursh k = 6 → d = 10 µm.
        assert!(
            (myelinated_fiber_diameter_um(60.0, 6.0) - 10.0).abs() <= 1e-9 * 10.0,
            "myelinated d = CV/k = 10 µm"
        );

        // (b) ROUND-TRIP myelinated threading myelinated_conduction_velocity (both directions).
        for &(d, k) in &[(10.0_f64, 6.0_f64), (5.0, 4.5)] {
            let cv = myelinated_conduction_velocity(d, k);
            assert!(
                (myelinated_conduction_velocity(myelinated_fiber_diameter_um(cv, k), k) - cv).abs()
                    <= 1e-9 * cv,
                "CV(d(CV)) = CV"
            );
            assert!((myelinated_fiber_diameter_um(cv, k) - d).abs() <= 1e-9 * d, "d(CV(d)) = d");
        }

        // (c) WORKED + ROUND-TRIP unmyelinated: v = 2 m/s at k = 2 → d = (2/2)² = 1 µm.
        assert!(
            (unmyelinated_fiber_diameter_um(2.0, 2.0) - 1.0).abs() <= 1e-9,
            "unmyelinated d = (CV/k)² = 1 µm"
        );
        for &(v, k) in &[(2.0_f64, 2.0_f64), (1.5, 0.8)] {
            assert!(
                (unmyelinated_conduction_velocity(unmyelinated_fiber_diameter_um(v, k), k) - v)
                    .abs()
                    <= 1e-9 * v,
                "CV(d(CV)) = CV (unmyelinated)"
            );
        }

        // (d) HURSH constant thread: a 60 m/s response is a ~10 µm fibre.
        assert!(
            (myelinated_fiber_diameter_um(60.0, HURSH_FACTOR_M_PER_S_PER_UM) - 10.0).abs()
                <= 1e-9 * 10.0,
            "Hursh: 60 m/s → 10 µm"
        );

        // (e) GUARD: non-physical input → 0 (both fns).
        assert_eq!(myelinated_fiber_diameter_um(60.0, 0.0), 0.0);
        assert_eq!(myelinated_fiber_diameter_um(-1.0, 6.0), 0.0);
        assert_eq!(myelinated_fiber_diameter_um(f64::NAN, 6.0), 0.0);
        assert_eq!(unmyelinated_fiber_diameter_um(2.0, 0.0), 0.0);
        assert_eq!(unmyelinated_fiber_diameter_um(-1.0, 2.0), 0.0);
        assert_eq!(unmyelinated_fiber_diameter_um(f64::NAN, 2.0), 0.0);
    }
}
