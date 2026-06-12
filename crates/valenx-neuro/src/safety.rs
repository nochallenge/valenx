//! Charge-injection safety limits for stimulating electrodes.
//!
//! A stimulating electrode injures tissue if it delivers too much charge per
//! phase relative to its area. The **Shannon (1992)** criterion captures this
//! with a single empirical line in the log–log plane of charge density `D`
//! (µC/cm²/phase) versus charge per phase `Q` (µC/phase):
//!
//! ```text
//! log10(D) = k − log10(Q)        ⟺        k = log10(D) + log10(Q) = log10(Q² / A)
//! ```
//!
//! where `A` is the electrode's geometric area (cm²) and `D = Q / A`. An
//! operating point's **k-value** measures how aggressive it is: Shannon's
//! cortical-stimulation data place the injury threshold near `k ≈ 1.7–2.0`,
//! and `k ≲ 1.5` is conservatively considered safe. Points below the chosen
//! limit line are safe; points above it are expected to damage tissue.
//!
//! This is a deliberately first-order model — a single line fit to historical
//! data. McCreery and others refine it (pulse-rate, duty-cycle, and electrode
//! material dependence), and it is **research / education-grade only**, not a
//! clinical safety guarantee.

/// Charge density `D = Q / A` (µC/cm²/phase) for a charge of `q_uc`
/// (µC/phase) injected through a geometric electrode area `area_cm2` (cm²).
pub fn charge_density(q_uc: f64, area_cm2: f64) -> f64 {
    q_uc / area_cm2
}

/// **Charge per phase** `Q = I·t` (µC/phase) — the charge a single rectangular stimulus
/// phase of amplitude `current_ma` (mA) and duration `pulse_width_ms` (ms) injects. Since
/// `1 mA · 1 ms = 1 µC`, the product is already in µC/phase: the source quantity that feeds
/// the Shannon chain — the [`charge_density`] (`Q/A`), the [`shannon_k`] value, and the
/// [`is_safe`] check, bounded by [`max_safe_charge_per_phase`].
pub fn charge_per_phase_uc(current_ma: f64, pulse_width_ms: f64) -> f64 {
    current_ma * pulse_width_ms
}

/// Shannon **k-value** of an operating point:
/// `k = log10(D) + log10(Q) = log10(Q² / A)`.
///
/// Larger `k` means a more aggressive (closer to / past the damage line)
/// stimulus. Inputs must be positive (µC/phase and cm²).
pub fn shannon_k(q_uc: f64, area_cm2: f64) -> f64 {
    (q_uc * q_uc / area_cm2).log10()
}

/// Whether an operating point sits strictly **below** a Shannon limit line
/// `k_limit` — i.e. its [`shannon_k`] value is `< k_limit`.
pub fn is_safe(q_uc: f64, area_cm2: f64, k_limit: f64) -> bool {
    shannon_k(q_uc, area_cm2) < k_limit
}

/// Largest charge per phase (µC/phase) that stays at or below `k_limit` for an
/// electrode of area `area_cm2`: solving `k_limit = log10(Q² / A)` gives
/// `Q_max = sqrt(A · 10^k_limit)`.
pub fn max_safe_charge_per_phase(area_cm2: f64, k_limit: f64) -> f64 {
    (area_cm2 * 10f64.powf(k_limit)).sqrt()
}

/// Largest charge **density** (µC/cm²/phase) that stays at or below `k_limit` for an
/// electrode of area `area_cm2` — the [`max_safe_charge_per_phase`] spread over the
/// electrode area, `D_max = Q_max / A = sqrt(10^k_limit / A)`. This is the y-axis of the
/// Shannon–McCreery plot: the charge density is the quantity directly compared against the
/// tissue-damage threshold, so this is the safe ceiling a stimulator's [`charge_density`]
/// must stay under.
pub fn max_safe_charge_density(area_cm2: f64, k_limit: f64) -> f64 {
    max_safe_charge_per_phase(area_cm2, k_limit) / area_cm2
}

/// The **minimum safe electrode area** `A_min = Q² / 10^k_limit` (cm²) — the smallest
/// geometric electrode area that keeps a charge per phase of `q_uc` `Q` (µC/phase) at or
/// below the Shannon limit line `k_limit`, inverting [`shannon_k`] for the area
/// (`k = log10(Q²/A)` ⟹ `A = Q²/10^k`). It is the electrode-sizing companion to
/// [`max_safe_charge_per_phase`] (the charge inverse): a smaller electrode concentrates the
/// same charge into a higher density, so below `A_min` the operating point crosses the
/// damage line. Inputs are `Q` in µC/phase and the dimensionless `k_limit`.
pub fn min_safe_electrode_area(q_uc: f64, k_limit: f64) -> f64 {
    q_uc * q_uc / 10f64.powf(k_limit)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn k_value_matches_the_definition() {
        // Q = 1 µC on A = 1 cm² → D = 1 → k = log10(1) + log10(1) = 0.
        assert!(shannon_k(1.0, 1.0).abs() < 1e-12);
        // Q = 10 µC, A = 1 cm² → k = log10(100) = 2.
        assert!((shannon_k(10.0, 1.0) - 2.0).abs() < 1e-12);
        // Charge density is simply Q / A.
        assert!((charge_density(10.0, 2.0) - 5.0).abs() < 1e-12);
    }

    #[test]
    fn max_safe_charge_scales_with_sqrt_area() {
        // Q_max = sqrt(A · 10^k). A=1,k=2 → 10; A=4,k=2 → 20 (√4 = 2× the charge).
        let q1 = max_safe_charge_per_phase(1.0, 2.0);
        let q4 = max_safe_charge_per_phase(4.0, 2.0);
        assert!((q1 - 10.0).abs() < 1e-9, "q1 = {q1}");
        assert!((q4 - 20.0).abs() < 1e-9, "q4 = {q4}");
        // A point exactly at Q_max sits on the limit line (k == k_limit).
        assert!((shannon_k(q1, 1.0) - 2.0).abs() < 1e-9);
    }

    #[test]
    fn max_safe_charge_density_is_the_boundary_in_density_terms() {
        let (a, k) = (0.01, 1.85); // 0.01 cm² electrode, classic Shannon k = 1.85
        let q_max = max_safe_charge_per_phase(a, k);
        let d_max = max_safe_charge_density(a, k);

        // Threads max_safe_charge_per_phase + charge_density: D_max = Q_max / A.
        assert!(
            (d_max - q_max / a).abs() <= 1e-9 * d_max,
            "D_max = Q_max / A"
        );
        assert!(
            (d_max - charge_density(q_max, a)).abs() <= 1e-9 * d_max,
            "D_max = charge_density(Q_max, A)"
        );

        // Boundary round-trip via shannon_k: at the max charge, k == k_limit.
        assert!(
            (shannon_k(q_max, a) - k).abs() <= 1e-9 * k,
            "shannon_k(Q_max, A) = k"
        );

        // is_safe straddles the boundary: just inside is safe, just outside is not.
        assert!(
            is_safe(0.999 * q_max, a, k),
            "just inside the limit is safe"
        );
        assert!(
            !is_safe(1.001 * q_max, a, k),
            "just outside the limit is unsafe"
        );
    }

    #[test]
    fn min_safe_electrode_area_is_the_charge_boundary_in_area_terms() {
        // (a) WORKED: A_min = Q²/10^k. Q = 10 µC, k = 2 → 100/100 = 1 cm².
        assert!(
            (min_safe_electrode_area(10.0, 2.0) - 1.0).abs() <= 1e-9,
            "A_min = Q²/10^k = 1 cm²"
        );

        // (b) ROUND-TRIP threading max_safe_charge_per_phase (non-tautological): the area
        // that just admits Q lets exactly Q through as the max safe charge. (c) BOUNDARY
        // threading shannon_k: at (Q, A_min) the k-value is exactly k_limit.
        for &(q, k) in &[(5.0_f64, 1.85_f64), (20.0, 2.0), (0.8, 1.5)] {
            let a = min_safe_electrode_area(q, k);
            assert!(
                (max_safe_charge_per_phase(a, k) - q).abs() <= 1e-9 * q,
                "max_safe_charge_per_phase(A_min, k) = Q"
            );
            assert!(
                (shannon_k(q, a) - k).abs() <= 1e-9 * k,
                "shannon_k(Q, A_min) = k"
            );
        }

        // (d) MONOTONICITY: more charge needs more area; a stricter (lower) k needs more
        // area for the same charge.
        assert!(
            min_safe_electrode_area(20.0, 2.0) > min_safe_electrode_area(10.0, 2.0),
            "more charge → larger area"
        );
        assert!(
            min_safe_electrode_area(10.0, 1.5) > min_safe_electrode_area(10.0, 2.0),
            "stricter k → larger area"
        );
    }

    #[test]
    fn charge_per_phase_uc_feeds_the_shannon_chain() {
        // (a) WORKED: 2 mA × 0.5 ms → Q = 1.0 µC (since mA·ms = µC).
        assert!(
            (charge_per_phase_uc(2.0, 0.5) - 1.0).abs() <= 1e-9,
            "Q = I·t = 1 µC"
        );

        // (b) THREAD charge_density (non-tautological): D = Q/A.
        let (i, pw, a) = (2.0_f64, 0.5_f64, 0.01_f64);
        assert!(
            (charge_density(charge_per_phase_uc(i, pw), a) - i * pw / a).abs()
                <= 1e-9 * (i * pw / a),
            "D = (I·t)/A"
        );

        // (c) THREAD max_safe_charge_per_phase + is_safe (non-tautological): a 2 mA × 0.5 ms
        // pulse on a 0.01 cm² electrode delivers exactly the k = 2.0 limit charge (1 µC).
        assert!(
            (charge_per_phase_uc(2.0, 0.5) - max_safe_charge_per_phase(0.01, 2.0)).abs() <= 1e-9,
            "Q = Q_max at k = 2.0"
        );
        assert!(
            !is_safe(charge_per_phase_uc(2.0, 0.5), 0.01, 2.0),
            "at the boundary → unsafe"
        );
        assert!(
            is_safe(charge_per_phase_uc(2.0, 0.495), 0.01, 2.0),
            "just below → safe"
        );

        // (d) PROPORTIONALITY + ZERO: linear in I and t; no current → no charge.
        let base = charge_per_phase_uc(2.0, 0.5);
        assert!(
            (charge_per_phase_uc(4.0, 0.5) - 2.0 * base).abs() <= 1e-9 * 2.0 * base,
            "∝ I"
        );
        assert!(
            (charge_per_phase_uc(2.0, 1.0) - 2.0 * base).abs() <= 1e-9 * 2.0 * base,
            "∝ t"
        );
        assert_eq!(charge_per_phase_uc(0.0, 0.5), 0.0);
    }

    #[test]
    fn safe_below_and_unsafe_above_the_line() {
        // k_limit = 2.0: Q=5,A=1 → k ≈ 1.40 (safe); Q=20,A=1 → k ≈ 2.60 (unsafe).
        assert!(is_safe(5.0, 1.0, 2.0));
        assert!(!is_safe(20.0, 1.0, 2.0));
        // A larger electrode raises the safe charge: the same 20 µC that is
        // unsafe on a 1 cm² pad is safe on a 100 cm² one.
        assert!(is_safe(20.0, 100.0, 2.0));
    }
}
