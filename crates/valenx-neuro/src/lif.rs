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

/// The **leaky integrate-and-fire rheobase current** `I_rh = V_th / R` (A) — the
/// minimum *constant* input current that makes the neuron fire at all. The leaky
/// membrane charges toward the steady value `R·I`; only if that asymptote exceeds
/// the threshold `v_threshold` `V_th` (V) does `V` ever reach it, so the firing
/// boundary is `R·I = V_th`, i.e. `I = V_th/R`. `resistance` is the input
/// resistance `R` (Ω).
///
/// This is the current-axis threshold the rest of the LIF family is defined
/// against: [`lif_firing_rate`] is `0` for `I ≤ I_rh` and positive above it,
/// [`lif_time_to_first_spike`] is `f64::INFINITY` for `I ≤ I_rh` and finite above
/// it, and [`lif_membrane_potential`] driven by exactly `I_rh` asymptotes to
/// `V_th` from below (the marginal drive that never quite fires). Returns
/// `f64::INFINITY` — no finite current elicits a spike — for non-physical input
/// (`R` or `V_th` non-positive, or either non-finite), mirroring the triad's
/// treatment of the subthreshold / non-physical regime.
pub fn lif_rheobase_current(resistance: f64, v_threshold: f64) -> f64 {
    if !resistance.is_finite()
        || resistance <= 0.0
        || !v_threshold.is_finite()
        || v_threshold <= 0.0
    {
        return f64::INFINITY;
    }
    v_threshold / resistance
}

/// The **leaky integrate-and-fire steady-state potential** `V∞ = R·I` (V) — the
/// passive subthreshold depolarisation the leaky membrane relaxes to under a
/// constant current `current` `I` (A) through input resistance `resistance` `R`
/// (Ω), measured from the `0` reset/rest level.
///
/// It is the `t → ∞` asymptote of the charging trajectory
/// [`lif_membrane_potential`] (`V(t) = R·I·(1 − e^(−t/τ)) → R·I`) and the quantity
/// the rheobase is defined against: firing needs `V∞ > V_th`, so the
/// [`lif_rheobase_current`] `I_rh = V_th/R` is exactly the current whose steady
/// state sits *at* the threshold. Linear in both the current and the resistance,
/// and sign-preserving (a hyperpolarising current gives a negative `V∞`). Returns
/// `0` for non-physical input (`R` non-positive, or either argument non-finite),
/// matching [`lif_membrane_potential`].
pub fn lif_steady_state_potential(current: f64, resistance: f64) -> f64 {
    if !current.is_finite() || !resistance.is_finite() || resistance <= 0.0 {
        return 0.0;
    }
    resistance * current
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

    #[test]
    fn lif_rheobase_current_is_the_firing_threshold() {
        let (r, tau, v_th, t_ref) = (1.0e8, 0.02, 0.015, 0.002); // the LIF fixture
        let i_rh = lif_rheobase_current(r, v_th);
        // Worked point: I_rh = V_th/R = 0.015 / 1e8 = 1.5e-10 A (0.15 nA) — exactly
        // the silent/fire boundary the f–I test brackets at 1.4e-10 / 1.5e-10.
        assert!((i_rh - 1.5e-10).abs() / i_rh < 1e-12, "I_rh = V_th/R = 0.15 nA, got {i_rh}");
        // Definitional identity: R·I_rh = V_th exactly (the asymptote sits AT threshold).
        assert!((r * i_rh - v_th).abs() / v_th < 1e-12, "R·I_rh = V_th");
        // Inversely proportional to R, linearly proportional to V_th.
        assert!((lif_rheobase_current(2.0 * r, v_th) - 0.5 * i_rh).abs() / i_rh < 1e-12, "∝ 1/R");
        assert!((lif_rheobase_current(r, 2.0 * v_th) - 2.0 * i_rh).abs() / i_rh < 1e-12, "∝ V_th");
        // STRONG non-tautological cross-check to the whole triad: just BELOW I_rh the
        // cell is silent (rate 0, infinite latency); just ABOVE it fires (rate > 0,
        // finite latency). The impl is V_th/R; the checks use the three triad fns.
        let (below, above) = (0.99 * i_rh, 1.01 * i_rh);
        assert_eq!(lif_firing_rate(below, r, tau, v_th, t_ref), 0.0, "below rheobase → silent");
        assert!(lif_time_to_first_spike(below, r, tau, v_th).is_infinite(), "below → ∞ latency");
        assert!(lif_firing_rate(above, r, tau, v_th, t_ref) > 0.0, "above rheobase → fires");
        assert!(lif_time_to_first_spike(above, r, tau, v_th).is_finite(), "above → finite latency");
        // Driven by EXACTLY the rheobase current the membrane asymptotes to V_th from
        // below (ties to #191): V(10·τ) is within 0.1% of V_th but strictly under it —
        // the marginal drive that never quite fires.
        let v_marginal = lif_membrane_potential(i_rh, r, tau, 10.0 * tau);
        assert!(
            v_marginal < v_th && v_marginal > 0.999 * v_th,
            "V(I_rh, 10τ) → V_th⁻, got {v_marginal}"
        );
        // Non-physical input → ∞ (no finite current elicits a spike).
        assert!(lif_rheobase_current(0.0, v_th).is_infinite()); // R ≤ 0
        assert!(lif_rheobase_current(-1.0e8, v_th).is_infinite()); // R < 0
        assert!(lif_rheobase_current(r, 0.0).is_infinite()); // V_th ≤ 0
        assert!(lif_rheobase_current(f64::NAN, v_th).is_infinite()); // non-finite R
        assert!(lif_rheobase_current(r, f64::INFINITY).is_infinite()); // non-finite V_th
    }

    #[test]
    fn lif_steady_state_potential_is_the_charging_asymptote_and_rheobase_basis() {
        // Worked point: I = 0.5 nA through R = 100 MΩ → V∞ = R·I = 0.05 V.
        let (i, r) = (0.5e-9, 100.0e6);
        let v_inf = lif_steady_state_potential(i, r);
        assert!((v_inf - 0.05).abs() < 1e-12, "V∞ = R·I = 0.05 V, got {v_inf}");
        // Linear in I and in R; sign-preserving (a hyperpolarising current → −V∞).
        assert!((lif_steady_state_potential(2.0 * i, r) - 2.0 * v_inf).abs() < 1e-12, "∝ I");
        assert!((lif_steady_state_potential(i, 3.0 * r) - 3.0 * v_inf).abs() < 1e-12, "∝ R");
        assert!((lif_steady_state_potential(-i, r) + v_inf).abs() < 1e-12, "sign-preserving");
        // STRONG cross-check (1): it is the t→∞ limit of the charging curve
        // lif_membrane_potential — at t = 20·τ the transient e^(−20) ≈ 2e-9 has died.
        let tau = 10.0e-3;
        let v_late = lif_membrane_potential(i, r, tau, 20.0 * tau);
        assert!(
            (v_inf - v_late).abs() / v_inf < 1e-7,
            "V∞ = lim_t lif_membrane_potential: {v_inf} vs {v_late}"
        );
        // STRONG cross-check (2): the rheobase identity V∞(I_rh) = V_th — the rheobase
        // current is exactly the one whose steady state sits at threshold.
        let v_th = 0.015; // 15 mV threshold
        let i_rh = lif_rheobase_current(r, v_th);
        assert!((lif_steady_state_potential(i_rh, r) - v_th).abs() < 1e-12, "V∞(I_rh) = V_th");
        // Non-physical input → 0 (matching lif_membrane_potential).
        assert_eq!(lif_steady_state_potential(i, 0.0), 0.0); // R ≤ 0
        assert_eq!(lif_steady_state_potential(i, -1.0e6), 0.0); // R < 0
        assert_eq!(lif_steady_state_potential(f64::NAN, r), 0.0); // non-finite I
        assert_eq!(lif_steady_state_potential(i, f64::INFINITY), 0.0); // non-finite R
    }
}
