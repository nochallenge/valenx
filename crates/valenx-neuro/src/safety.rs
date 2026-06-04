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
    fn safe_below_and_unsafe_above_the_line() {
        // k_limit = 2.0: Q=5,A=1 → k ≈ 1.40 (safe); Q=20,A=1 → k ≈ 2.60 (unsafe).
        assert!(is_safe(5.0, 1.0, 2.0));
        assert!(!is_safe(20.0, 1.0, 2.0));
        // A larger electrode raises the safe charge: the same 20 µC that is
        // unsafe on a 1 cm² pad is safe on a 100 cm² one.
        assert!(is_safe(20.0, 100.0, 2.0));
    }
}
