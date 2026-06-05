//! Synaptic conductance time courses.
//!
//! When a presynaptic spike arrives, transmitter binding opens postsynaptic
//! channels and the synaptic conductance rises sharply, then decays. The
//! **alpha function** (Rall, 1967) is the classic single-parameter waveform for
//! that transient — the simplest shape that captures the finite rise time the
//! cruder instantaneous-jump model misses, while needing only one time constant.

/// The **alpha-function** synaptic conductance `g(t) = g_max·(t/τ)·e^(1−t/τ)` at
/// time `t_s` (s) after a presynaptic event at `t = 0`, with peak conductance
/// `g_max` and time-to-peak `tau_s` (`τ`, s).
///
/// The `e^(1−t/τ)` normalisation makes the conductance reach **exactly** `g_max`
/// at `t = τ`: it is zero at the moment of the spike, rises to the peak at `τ`,
/// then decays smoothly back toward zero — the textbook single-parameter
/// postsynaptic conductance transient. (Its area is `g_max·τ·e`, the total
/// charge the synapse passes.) Returns `0` before the event (`t < 0`) and for
/// non-physical input (`τ ≤ 0`, or any non-finite argument).
pub fn alpha_synapse_conductance(g_max: f64, tau_s: f64, t_s: f64) -> f64 {
    if !g_max.is_finite()
        || !tau_s.is_finite()
        || tau_s <= 0.0
        || !t_s.is_finite()
        || t_s < 0.0
    {
        return 0.0;
    }
    let x = t_s / tau_s;
    g_max * x * (1.0 - x).exp()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn alpha_synapse_peaks_at_g_max_at_the_time_constant() {
        let (g_max, tau) = (5.0e-9, 0.002); // 5 nS peak, 2 ms time-to-peak
        // No conductance at the moment of the presynaptic spike.
        assert_eq!(alpha_synapse_conductance(g_max, tau, 0.0), 0.0);
        // The peak is exactly g_max at t = τ.
        assert!((alpha_synapse_conductance(g_max, tau, tau) - g_max).abs() < 1e-18, "peak at τ");
        // Below the peak on both sides — rising before, decaying after — and it
        // never exceeds g_max.
        let before = alpha_synapse_conductance(g_max, tau, tau / 2.0);
        let after = alpha_synapse_conductance(g_max, tau, 2.0 * tau);
        assert!(before < g_max && before > 0.0, "rising before the peak: {before}");
        assert!(after < g_max && after > 0.0, "decaying after the peak: {after}");
        // Matches the closed form at an arbitrary time.
        let t = 0.005;
        let x = t / tau;
        let expected = g_max * x * (1.0 - x).exp();
        assert!((alpha_synapse_conductance(g_max, tau, t) - expected).abs() < 1e-18);
        // Decays toward zero long after the event.
        assert!(alpha_synapse_conductance(g_max, tau, 50.0 * tau) < 0.001 * g_max);
        // Before the event and non-physical inputs → 0.
        assert_eq!(alpha_synapse_conductance(g_max, tau, -0.001), 0.0);
        assert_eq!(alpha_synapse_conductance(g_max, 0.0, tau), 0.0);
        assert_eq!(alpha_synapse_conductance(g_max, f64::NAN, tau), 0.0);
    }
}
