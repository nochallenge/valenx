//! The current–distance relationship for extracellular microstimulation.
//!
//! Stoney et al. (1968) and BeMent & Ranck established that the current needed
//! to fire a neuron with a point electrode grows **quadratically** with the
//! electrode–neuron distance:
//!
//! ```text
//! I_th(r) = I_0 + k · r²
//! ```
//!
//! - `I_0` (µA) — the at-electrode threshold (`r = 0`).
//! - `k` (µA/µm²) — the **current–distance constant** characterising
//!   excitability. Cortical microstimulation values are ~1000–4000 µA/mm²
//!   (i.e. `0.001–0.004` µA/µm²); a smaller `k` means a more excitable target
//!   reachable from farther away.
//! - `r` (µm) — distance.
//!
//! Inverting the law gives the **activation radius** — the radius within which
//! every fiber fires for a given stimulus current — which is the geometric
//! basis for stimulation selectivity. This is a phenomenological population
//! law (research/education-grade), complementary to the mechanistic
//! field → activating-function → cable path elsewhere in this crate.

/// Threshold current (µA) to activate a fiber at distance `r_um` (µm):
/// `I_th = I_0 + k · r²`.
pub fn threshold_current(i0_ua: f64, k_ua_per_um2: f64, r_um: f64) -> f64 {
    i0_ua + k_ua_per_um2 * r_um * r_um
}

/// Activation radius (µm) for a stimulus of `i_ua` (µA): the largest distance
/// at which a fiber still fires, from inverting `I = I_0 + k · r²`:
/// `r = sqrt((I − I_0) / k)`. Returns `0` when the current is at or below the
/// at-electrode threshold, or `k` is non-positive.
pub fn activation_radius(i0_ua: f64, k_ua_per_um2: f64, i_ua: f64) -> f64 {
    if k_ua_per_um2 <= 0.0 || i_ua <= i0_ua {
        return 0.0;
    }
    ((i_ua - i0_ua) / k_ua_per_um2).sqrt()
}

/// Fit the current–distance constant `k` (µA/µm²) from one measured
/// `(distance, threshold)` pair and a known at-electrode threshold `i0_ua`:
/// `k = (I_th − I_0) / r²`. Returns `0` for `r = 0` (no distance to fit).
pub fn fit_constant(i0_ua: f64, r_um: f64, i_th_ua: f64) -> f64 {
    if r_um == 0.0 {
        return 0.0;
    }
    (i_th_ua - i0_ua) / (r_um * r_um)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn threshold_is_quadratic_in_distance() {
        // I_th = I_0 + k·r²  with I_0 = 10 µA, k = 2 µA/µm².
        assert!((threshold_current(10.0, 2.0, 0.0) - 10.0).abs() < 1e-12);
        assert!((threshold_current(10.0, 2.0, 1.0) - 12.0).abs() < 1e-12);
        // Doubling the distance quadruples the excess over I_0 (2 → 8).
        let excess_1 = threshold_current(10.0, 2.0, 1.0) - 10.0;
        let excess_2 = threshold_current(10.0, 2.0, 2.0) - 10.0;
        assert!((excess_2 - 4.0 * excess_1).abs() < 1e-12);
    }

    #[test]
    fn activation_radius_inverts_the_threshold() {
        // The current that just fires a fiber at r = 3 µm activates exactly
        // out to r = 3 µm.
        let i = threshold_current(10.0, 2.0, 3.0);
        assert!((activation_radius(10.0, 2.0, i) - 3.0).abs() < 1e-9);
        // At or below the at-electrode threshold, nothing beyond r = 0 fires.
        assert_eq!(activation_radius(10.0, 2.0, 10.0), 0.0);
        assert_eq!(activation_radius(10.0, 2.0, 5.0), 0.0);
    }

    #[test]
    fn fit_recovers_the_constant() {
        let i_th = threshold_current(10.0, 2.0, 3.0); // 28 µA
        assert!((fit_constant(10.0, 3.0, i_th) - 2.0).abs() < 1e-12);
    }
}
