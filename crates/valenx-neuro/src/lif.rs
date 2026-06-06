//! Leaky integrate-and-fire (LIF) neuron — the reduced spiking model.
//!
//! Where the [`crate::cable`] module integrates the full Hodgkin–Huxley
//! membrane, the LIF strips the spike down to its essence: a leaky `RC`
//! membrane that charges toward `R·I`, fires the instant it crosses a fixed
//! threshold, then resets and sits out a refractory gap. It is the standard
//! single-neuron building block of large network and firing-rate models, where
//! the cost of a full conductance-based spike per cell is prohibitive.

/// The **leaky integrate-and-fire firing rate** `f` (Hz) of a neuron driven by a
/// *constant* input current `current` `I` (A), in closed form:
///
/// ```text
///   f = 1 / (t_ref + τ_m·ln( R·I / (R·I − V_th) ))     for  R·I > V_th
///   f = 0                                              otherwise
/// ```
///
/// The leaky membrane `τ_m·dV/dt = −V + R·I` charges exponentially toward its
/// steady value `R·I`; the time to climb from reset (`0`) to the threshold
/// `v_threshold` `V_th` is `τ_m·ln(R·I/(R·I − V_th))`, and adding the refractory
/// dead time `t_refractory_s` `t_ref` gives the inter-spike interval whose
/// reciprocal is the rate. `resistance` is the input resistance `R` (Ω) and
/// `tau_m_s` the membrane time constant `τ_m` (s).
///
/// If the asymptotic drive `R·I` does not exceed the threshold the membrane
/// never reaches it and the cell is silent (`f = 0`) — the rheobase of the LIF.
/// As the current grows the logarithm shrinks toward `0` and the rate saturates
/// at the refractory ceiling `1/t_ref`. Returns `0` for non-physical input
/// (`R`, `τ_m`, or `V_th` non-positive, `t_ref` negative, or any non-finite
/// argument).
pub fn lif_firing_rate(
    current: f64,
    resistance: f64,
    tau_m_s: f64,
    v_threshold: f64,
    t_refractory_s: f64,
) -> f64 {
    if !t_refractory_s.is_finite() || t_refractory_s < 0.0 {
        return 0.0;
    }
    // The inter-spike interval is the refractory dead time plus the climb time
    // from reset to threshold; the rate is its reciprocal. A non-finite latency
    // means a subthreshold (or non-physical) drive — the cell is silent.
    let climb = lif_time_to_first_spike(current, resistance, tau_m_s, v_threshold);
    if !climb.is_finite() {
        return 0.0;
    }
    1.0 / (t_refractory_s + climb)
}

/// The **leaky integrate-and-fire subthreshold membrane potential** `V(t)` (V)
/// of a neuron driven from rest (`V(0) = 0`) by a *constant* input current
/// `current` `I` (A), the closed-form `RC` charging curve:
///
/// ```text
///   V(t) = R·I·(1 − e^(−t/τ_m))     for  t ≥ 0
/// ```
///
/// The leaky membrane `τ_m·dV/dt = −V + R·I` relaxes exponentially toward its
/// steady value `R·I` with time constant `tau_m_s` `τ_m` (s); `resistance` is the
/// input resistance `R` (Ω). This is the *time-course* companion to
/// [`lif_firing_rate`]: that gives the steady spike rate, this gives the membrane
/// trajectory leading up to a spike (the depolarisation that, once it reaches the
/// threshold, triggers one). `V(0) = 0` at the reset, `V(τ_m) = R·I·(1 − 1/e) ≈
/// 0.632·R·I`, and `V → R·I` as `t → ∞`. If the asymptote `R·I` is below the
/// firing threshold the cell never spikes and `V` simply saturates there (the LIF
/// rheobase); a hyperpolarising (negative) current drives `V` negative the same
/// way. Returns `0` before the stimulus (`t < 0`) and for non-physical input
/// (`R` or `τ_m` non-positive, or any non-finite argument).
pub fn lif_membrane_potential(current: f64, resistance: f64, tau_m_s: f64, t_s: f64) -> f64 {
    if !current.is_finite()
        || !resistance.is_finite()
        || resistance <= 0.0
        || !tau_m_s.is_finite()
        || tau_m_s <= 0.0
        || !t_s.is_finite()
        || t_s < 0.0
    {
        return 0.0;
    }
    resistance * current * (1.0 - (-t_s / tau_m_s).exp())
}

/// The **leaky integrate-and-fire time-to-first-spike** (response latency, s) —
/// the time `t₁ = τ_m·ln(R·I/(R·I − V_th))` for the leaky membrane, driven from
/// reset by a *constant* current `current` `I` (A), to climb to the threshold
/// `v_threshold` `V_th` and fire. `resistance` is the input resistance `R` (Ω)
/// and `tau_m_s` the membrane time constant `τ_m` (s).
///
/// This is the spike-timing companion to the rate and the trajectory: it is the
/// climb time [`lif_firing_rate`] adds to the refractory dead time to form the
/// inter-spike interval (`f = 1/(t_ref + t₁)`), and the exact instant at which
/// [`lif_membrane_potential`] reaches `V_th`. A stronger drive `R·I` shortens the
/// latency (the membrane charges past threshold sooner), and as `R·I` approaches
/// the threshold from above the latency diverges. Returns `f64::INFINITY` — the
/// cell never spikes — when the asymptotic drive does not exceed the threshold
/// (`R·I ≤ V_th`, the rheobase) or for non-physical input (`R`, `τ_m`, or `V_th`
/// non-positive, or any non-finite argument).
pub fn lif_time_to_first_spike(
    current: f64,
    resistance: f64,
    tau_m_s: f64,
    v_threshold: f64,
) -> f64 {
    if !current.is_finite()
        || !resistance.is_finite()
        || resistance <= 0.0
        || !tau_m_s.is_finite()
        || tau_m_s <= 0.0
        || !v_threshold.is_finite()
        || v_threshold <= 0.0
    {
        return f64::INFINITY;
    }
    let drive = resistance * current; // the asymptotic membrane voltage R·I
    if drive <= v_threshold {
        return f64::INFINITY; // subthreshold — the membrane never reaches threshold
    }
    tau_m_s * (drive / (drive - v_threshold)).ln()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lif_firing_rate_traces_the_f_i_curve() {
        let (r, tau, v_th, t_ref) = (1.0e8, 0.02, 0.015, 0.002); // 100 MΩ, 20 ms, 15 mV, 2 ms
        // Subthreshold: R·I ≤ V_th → the cell never fires (the LIF rheobase).
        assert_eq!(lif_firing_rate(1.4e-10, r, tau, v_th, t_ref), 0.0); // R·I = 14 mV < 15
        assert_eq!(lif_firing_rate(1.5e-10, r, tau, v_th, t_ref), 0.0); // R·I = 15 mV = V_th
        // Worked point: I = 0.2 nA → R·I = 20 mV → f = 1/(t_ref + τ·ln 4) ≈ 33.6 Hz.
        let f = lif_firing_rate(2.0e-10, r, tau, v_th, t_ref);
        let expected = 1.0 / (t_ref + tau * (0.02_f64 / 0.005).ln());
        assert!((f - expected).abs() < 1e-9, "closed form");
        assert!((f - 33.64).abs() < 0.1, "≈33.6 Hz, got {f}");
        // Monotone increasing with current.
        assert!(
            lif_firing_rate(5.0e-10, r, tau, v_th, t_ref) > f,
            "more current → higher rate"
        );
        // Saturates toward the refractory ceiling 1/t_ref as I → ∞ (the ln → 0).
        let f_huge = lif_firing_rate(1.0e-6, r, tau, v_th, t_ref);
        assert!(f_huge < 1.0 / t_ref, "stays below the 1/t_ref ceiling");
        assert!(f_huge > 0.99 / t_ref, "approaches the 1/t_ref ceiling, got {f_huge}");
        // Non-physical input → 0.
        assert_eq!(lif_firing_rate(2.0e-10, 0.0, tau, v_th, t_ref), 0.0);
        assert_eq!(lif_firing_rate(2.0e-10, r, 0.0, v_th, t_ref), 0.0);
        assert_eq!(lif_firing_rate(2.0e-10, r, tau, 0.0, t_ref), 0.0);
        assert_eq!(lif_firing_rate(2.0e-10, r, tau, v_th, -1.0), 0.0);
        assert_eq!(lif_firing_rate(f64::NAN, r, tau, v_th, t_ref), 0.0);
    }

    #[test]
    fn lif_membrane_potential_traces_the_rc_charging_curve() {
        use std::f64::consts::E;
        let (r, tau) = (1.0e8, 0.02); // 100 MΩ, 20 ms (the firing-rate fixture)
        let i = 2.0e-10; // 0.2 nA → R·I = 20 mV, the worked point
        let drive = r * i; // 0.02 V steady value
        // Starts from rest at the reset.
        assert_eq!(lif_membrane_potential(i, r, tau, 0.0), 0.0);
        // One time constant → R·I·(1 − 1/e) ≈ 0.632·R·I.
        let v_tau = lif_membrane_potential(i, r, tau, tau);
        assert!((v_tau - drive * (1.0 - 1.0 / E)).abs() / drive < 1e-12, "V(τ) = R·I·(1−1/e)");
        // Saturates toward the steady value R·I: at 10·τ_m the membrane is within
        // 0.1% of R·I but still strictly below it (e^(−10) ≈ 4.5e-5; a larger t
        // would round the exponential gap below f64 precision at this R·I).
        let v_inf = lif_membrane_potential(i, r, tau, 10.0 * tau);
        assert!(v_inf < drive && v_inf > 0.999 * drive, "V → R·I, got {v_inf} vs {drive}");
        // Monotone increasing with time for a depolarising current.
        assert!(lif_membrane_potential(i, r, tau, 2.0 * tau) > v_tau, "charges over time");
        // Cross-check tying the trajectory to lif_firing_rate (#173): at the climb
        // time t* = τ·ln(R·I/(R·I − V_th)) the membrane reaches EXACTLY V_th.
        let v_th = 0.015;
        let t_climb = lif_time_to_first_spike(i, r, tau, v_th);
        let v_at_climb = lif_membrane_potential(i, r, tau, t_climb);
        assert!((v_at_climb - v_th).abs() < 1e-9, "V(t*) = V_th, got {v_at_climb}");
        // A hyperpolarising (negative) current drives V negative.
        assert!(lif_membrane_potential(-i, r, tau, tau) < 0.0, "negative current → negative V");
        // Before the stimulus and non-physical input → 0.
        assert_eq!(lif_membrane_potential(i, r, tau, -0.001), 0.0);
        assert_eq!(lif_membrane_potential(i, 0.0, tau, tau), 0.0); // R ≤ 0
        assert_eq!(lif_membrane_potential(i, r, 0.0, tau), 0.0); // τ_m ≤ 0
        assert_eq!(lif_membrane_potential(f64::NAN, r, tau, tau), 0.0); // non-finite I
        assert_eq!(lif_membrane_potential(i, r, tau, f64::INFINITY), 0.0); // non-finite t
    }

    #[test]
    fn lif_time_to_first_spike_is_the_latency_to_threshold() {
        let (r, tau, v_th) = (1.0e8, 0.02, 0.015); // 100 MΩ, 20 ms, 15 mV
        let i = 2.0e-10; // 0.2 nA → R·I = 20 mV, the firing-rate worked point
        let t1 = lif_time_to_first_spike(i, r, tau, v_th);
        // Closed form τ·ln(R·I/(R·I − V_th)) = 0.02·ln(0.02/0.005) = 0.02·ln 4 ≈ 27.7 ms.
        assert!((t1 - tau * (0.02_f64 / 0.005).ln()).abs() < 1e-12, "closed form");
        assert!((t1 - 0.02773).abs() < 1e-4, "≈27.7 ms, got {t1}");
        // The membrane reaches threshold EXACTLY at this latency (ties to #191,
        // non-tautological): V(t1) = R·I·(1 − e^(−t1/τ)) = V_th.
        let v_at_t1 = lif_membrane_potential(i, r, tau, t1);
        assert!((v_at_t1 - v_th).abs() < 1e-12, "V(t1) = V_th, got {v_at_t1}");
        // It is the climb time the firing rate is built on: f = 1/(t_ref + t1).
        let t_ref = 0.002;
        let f = lif_firing_rate(i, r, tau, v_th, t_ref);
        assert!((f - 1.0 / (t_ref + t1)).abs() < 1e-9, "f = 1/(t_ref + t1)");
        // A stronger drive shortens the latency (charges past threshold sooner).
        assert!(
            lif_time_to_first_spike(5.0e-10, r, tau, v_th) < t1,
            "more current → shorter latency"
        );
        // Subthreshold (R·I ≤ V_th) → never spikes → infinite latency.
        assert!(lif_time_to_first_spike(1.4e-10, r, tau, v_th).is_infinite(), "subthreshold → ∞");
        // Non-physical input → ∞ (so the firing rate it feeds returns 0).
        assert!(lif_time_to_first_spike(i, 0.0, tau, v_th).is_infinite()); // R ≤ 0
        assert!(lif_time_to_first_spike(i, r, 0.0, v_th).is_infinite()); // τ_m ≤ 0
        assert!(lif_time_to_first_spike(i, r, tau, 0.0).is_infinite()); // V_th ≤ 0
        assert!(lif_time_to_first_spike(f64::NAN, r, tau, v_th).is_infinite()); // non-finite I
    }
}
