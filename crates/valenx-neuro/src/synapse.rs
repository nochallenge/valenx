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

/// The **double-exponential** ("bi-exponential", Exp2Syn) synaptic conductance
/// `g(t) = g_max·[e^(−t/τ_d) − e^(−t/τ_r)] / [e^(−t_p/τ_d) − e^(−t_p/τ_r)]` at
/// time `t_s` (s) after a presynaptic event, with peak conductance `g_max`, rise
/// time constant `tau_rise_s` (`τ_r`, s) and decay time constant `tau_decay_s`
/// (`τ_d`, s).
///
/// This is the general two-time-constant postsynaptic conductance — the standard
/// waveform for AMPA / GABA / NMDA synapses, with an independent fast rise and a
/// slower decay. It **generalises** [`alpha_synapse_conductance`], which is its
/// `τ_r = τ_d` limit (a single time constant): here the two are set separately,
/// so the rise and decay are shaped independently. The denominator normalises
/// the peak to **exactly** `g_max`, reached at the analytic peak time
/// `t_p = (τ_r·τ_d/(τ_d − τ_r))·ln(τ_d/τ_r)`. (By convention `τ_r < τ_d`, but the
/// expression is symmetric in the two constants.)
///
/// Returns `0` at and before the event (`g(0) = 0`, `t < 0`), for non-physical
/// input (any `τ ≤ 0`, or non-finite), and for the degenerate `τ_r = τ_d`
/// singularity — use [`alpha_synapse_conductance`] for equal time constants.
pub fn dual_exponential_synapse_conductance(
    g_max: f64,
    tau_rise_s: f64,
    tau_decay_s: f64,
    t_s: f64,
) -> f64 {
    if !g_max.is_finite()
        || !tau_rise_s.is_finite()
        || !tau_decay_s.is_finite()
        || tau_rise_s <= 0.0
        || tau_decay_s <= 0.0
        || !t_s.is_finite()
        || t_s < 0.0
    {
        return 0.0;
    }
    // The τ_r = τ_d limit is the alpha function — this form is 0/0 there.
    if (tau_rise_s - tau_decay_s).abs() <= f64::EPSILON * tau_rise_s.max(tau_decay_s) {
        return 0.0;
    }
    let t_peak =
        (tau_rise_s * tau_decay_s / (tau_decay_s - tau_rise_s)) * (tau_decay_s / tau_rise_s).ln();
    let norm = (-t_peak / tau_decay_s).exp() - (-t_peak / tau_rise_s).exp();
    let raw = (-t_s / tau_decay_s).exp() - (-t_s / tau_rise_s).exp();
    g_max * raw / norm
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

    #[test]
    fn dual_exponential_synapse_peaks_at_g_max_at_the_analytic_peak_time() {
        let (g_max, tau_r, tau_d) = (1.0e-9, 0.0005, 0.003); // 1 nS, τ_r = 0.5 ms, τ_d = 3 ms
        // No conductance at the moment of the presynaptic spike.
        assert_eq!(
            dual_exponential_synapse_conductance(g_max, tau_r, tau_d, 0.0),
            0.0
        );
        // The peak is exactly g_max at the analytic peak time t_p.
        let tp = (tau_r * tau_d / (tau_d - tau_r)) * (tau_d / tau_r).ln();
        assert!(
            (dual_exponential_synapse_conductance(g_max, tau_r, tau_d, tp) - g_max).abs() < 1e-18,
            "peak g_max at t_p"
        );
        // Rising before the peak, decaying after, never exceeding g_max.
        let before = dual_exponential_synapse_conductance(g_max, tau_r, tau_d, tp / 2.0);
        let after = dual_exponential_synapse_conductance(g_max, tau_r, tau_d, 2.0 * tp);
        assert!(before > 0.0 && before < g_max, "rising before the peak: {before}");
        assert!(after > 0.0 && after < g_max, "decaying after the peak: {after}");
        assert!(
            before < dual_exponential_synapse_conductance(g_max, tau_r, tau_d, 0.8 * tp),
            "monotonic rise toward t_p"
        );
        // Matches the closed form at an arbitrary time.
        let t = 0.002;
        let norm = (-tp / tau_d).exp() - (-tp / tau_r).exp();
        let expected = g_max * ((-t / tau_d).exp() - (-t / tau_r).exp()) / norm;
        assert!(
            (dual_exponential_synapse_conductance(g_max, tau_r, tau_d, t) - expected).abs() < 1e-18
        );
        // Decays toward zero long after the event (the slow τ_d dominates).
        assert!(dual_exponential_synapse_conductance(g_max, tau_r, tau_d, 50.0 * tau_d) < 0.001 * g_max);
        // Before the event, non-physical input, and the τ_r = τ_d singularity → 0.
        assert_eq!(dual_exponential_synapse_conductance(g_max, tau_r, tau_d, -0.001), 0.0);
        assert_eq!(dual_exponential_synapse_conductance(g_max, 0.0, tau_d, 0.001), 0.0);
        assert_eq!(dual_exponential_synapse_conductance(g_max, tau_r, -1.0, 0.001), 0.0);
        assert_eq!(dual_exponential_synapse_conductance(g_max, f64::NAN, tau_d, 0.001), 0.0);
        assert_eq!(dual_exponential_synapse_conductance(g_max, 0.002, 0.002, 0.001), 0.0); // τ_r = τ_d
    }
}
