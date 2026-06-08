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
/// postsynaptic conductance transient. (Its time-integral is `g_max·τ·e` — the
/// synaptic efficacy, exposed as [`alpha_synapse_conductance_integral`].)
/// Returns `0` before the event (`t < 0`) and for
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

/// The **time-integral of the alpha-function conductance**,
/// `∫₀^∞ g(t) dt = g_max·τ·e` (S·s) — the area under the
/// [`alpha_synapse_conductance`] transient, with peak conductance `g_max` and
/// time-to-peak `tau_s` (`τ`, s).
///
/// This single number is the synapse's **efficacy** (its "weight"): multiplied
/// by the synaptic driving force `(V − E_syn)` it gives the total charge the
/// event injects postsynaptically, so two synapses with the same integral
/// deliver the same charge however differently their `g_max` and `τ` trade off.
/// The closed form follows from `∫₀^∞ (t/τ)·e^(1−t/τ) dt = τ·e` (substitute
/// `u = t/τ`: `e·∫₀^∞ u·e^(−u) du = e·Γ(2) = e`), so it is linear in both `g_max`
/// and `τ`. Returns `0` for non-physical input (`τ ≤ 0`, or any non-finite
/// argument), matching [`alpha_synapse_conductance`].
pub fn alpha_synapse_conductance_integral(g_max: f64, tau_s: f64) -> f64 {
    if !g_max.is_finite() || !tau_s.is_finite() || tau_s <= 0.0 {
        return 0.0;
    }
    g_max * tau_s * std::f64::consts::E
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
    let t_peak = dual_exponential_peak_time(tau_rise_s, tau_decay_s);
    let norm = (-t_peak / tau_decay_s).exp() - (-t_peak / tau_rise_s).exp();
    let raw = (-t_s / tau_decay_s).exp() - (-t_s / tau_rise_s).exp();
    g_max * raw / norm
}

/// The **time-to-peak** of the double-exponential synaptic conductance,
/// `t_p = (τ_r·τ_d / (τ_d − τ_r))·ln(τ_d/τ_r)` (s) — the moment after the
/// presynaptic event at which [`dual_exponential_synapse_conductance`] reaches
/// its maximum `g_max`, given the rise time constant `tau_rise_s` (`τ_r`, s) and
/// decay time constant `tau_decay_s` (`τ_d`, s).
///
/// It is the synapse's **rise time** — the characteristic the single decay
/// constant cannot express, and what distinguishes a fast AMPA contact from a
/// slow NMDA one. It always lies strictly between `τ_r` and `τ_d`, and is
/// symmetric in the two constants. (For equal time constants the waveform is the
/// alpha function, whose peak is simply `τ` — see [`alpha_synapse_conductance`].)
/// Returns `0` for non-physical input (`τ ≤ 0` or non-finite) and for the
/// degenerate `τ_r = τ_d` singularity, matching
/// [`dual_exponential_synapse_conductance`].
pub fn dual_exponential_peak_time(tau_rise_s: f64, tau_decay_s: f64) -> f64 {
    if !tau_rise_s.is_finite()
        || !tau_decay_s.is_finite()
        || tau_rise_s <= 0.0
        || tau_decay_s <= 0.0
    {
        return 0.0;
    }
    if (tau_rise_s - tau_decay_s).abs() <= f64::EPSILON * tau_rise_s.max(tau_decay_s) {
        return 0.0; // the τ_r = τ_d limit is the alpha function (peak at τ)
    }
    (tau_rise_s * tau_decay_s / (tau_decay_s - tau_rise_s)) * (tau_decay_s / tau_rise_s).ln()
}

/// The **time-integral of the double-exponential conductance**,
/// `∫₀^∞ g(t) dt = g_max·(τ_d − τ_r)/norm` (S·s) — the synapse's **efficacy** (its
/// "weight"), the dual-exponential analogue of
/// [`alpha_synapse_conductance_integral`]. Multiplied by the synaptic driving
/// force `(V − E_syn)` it gives the total charge a single event injects
/// postsynaptically, so it is the size-of-the-PSP figure independent of how the
/// rise time constant `tau_rise_s` (`τ_r`, s) and decay time constant
/// `tau_decay_s` (`τ_d`, s) trade off; `g_max` is the (normalised) peak
/// conductance.
///
/// The bare difference of exponentials integrates to `τ_d − τ_r`
/// (`∫₀^∞ e^(−t/τ) dt = τ`), and the peak-normalisation factor `norm` (the same one
/// [`dual_exponential_synapse_conductance`] divides by, evaluated at the peak time
/// [`dual_exponential_peak_time`]) scales that to the `g_max`-peaked waveform — so
/// the efficacy is linear in `g_max`. Returns `0` for non-physical input (any
/// `τ ≤ 0` or non-finite) and for the degenerate `τ_r = τ_d` singularity (use
/// [`alpha_synapse_conductance_integral`], `g_max·τ·e`, for equal time constants),
/// matching [`dual_exponential_synapse_conductance`].
pub fn dual_exponential_conductance_integral(g_max: f64, tau_rise_s: f64, tau_decay_s: f64) -> f64 {
    if !g_max.is_finite()
        || !tau_rise_s.is_finite()
        || !tau_decay_s.is_finite()
        || tau_rise_s <= 0.0
        || tau_decay_s <= 0.0
    {
        return 0.0;
    }
    // The τ_r = τ_d limit is the alpha function — this form is 0/0 there.
    if (tau_rise_s - tau_decay_s).abs() <= f64::EPSILON * tau_rise_s.max(tau_decay_s) {
        return 0.0;
    }
    let t_peak = dual_exponential_peak_time(tau_rise_s, tau_decay_s);
    let norm = (-t_peak / tau_decay_s).exp() - (-t_peak / tau_rise_s).exp();
    g_max * (tau_decay_s - tau_rise_s) / norm
}

/// The **NMDA-receptor Mg²⁺ block** `B(V) = 1/(1 + ([Mg²⁺]/3.57)·e^(−0.062·V))` —
/// the fraction of NMDA channels *unblocked* (conducting) at membrane potential
/// `voltage_mv` `V` (mV) and extracellular magnesium concentration `mg_conc_mm`
/// `[Mg²⁺]` (mM), from the Jahr–Stevens (1990) fit.
///
/// At rest the channel pore is plugged by an extracellular Mg²⁺ ion, so the
/// receptor barely conducts; depolarising the membrane electrostatically expels
/// the ion and *relieves* the block. This voltage dependence is what makes the
/// NMDA receptor a molecular **coincidence detector** — it passes current only
/// when presynaptic transmitter (which opens the channel) and postsynaptic
/// depolarisation (which clears the Mg²⁺) arrive together, the biophysical basis
/// of Hebbian long-term potentiation. It is the orthogonal voltage *gate* that
/// multiplies the synaptic conductance *time course*
/// ([`alpha_synapse_conductance`], [`dual_exponential_synapse_conductance`]):
/// the full NMDA current is `g_max · g(t) · B(V) · (V − E_syn)`.
///
/// `B` rises monotonically from ≈0 (deeply hyperpolarised, blocked) toward `1`
/// (strongly depolarised, fully relieved); with no magnesium (`[Mg²⁺] = 0`) it
/// is `1` at every voltage (nothing to block). Returns `0` for non-physical
/// input (non-finite `V` or `[Mg²⁺]`, or negative `[Mg²⁺]`).
pub fn nmda_mg_block(voltage_mv: f64, mg_conc_mm: f64) -> f64 {
    if !voltage_mv.is_finite() || !mg_conc_mm.is_finite() || mg_conc_mm < 0.0 {
        return 0.0;
    }
    1.0 / (1.0 + (mg_conc_mm / 3.57) * (-0.062 * voltage_mv).exp())
}

/// The **NMDA Mg²⁺ half-block voltage** `V½ = ln([Mg²⁺]/3.57) / 0.062` (mV) — the
/// membrane potential at which the Jahr–Stevens block [`nmda_mg_block`] is exactly
/// half-relieved (`B = 0.5`), for extracellular magnesium `mg_conc_mm` `[Mg²⁺]` (mM).
/// It is the inflection of the receptor's sigmoidal voltage dependence and the single
/// number that summarises how depolarised the cell must be to unblock the channel:
/// ≈ −20.5 mV at the physiological 1 mM `[Mg²⁺]`, rising to `0` mV at 3.57 mM and more
/// positive as magnesium increases. Returns `NaN` for non-physical input (non-finite
/// or non-positive `[Mg²⁺]`, where there is no block to half-relieve).
pub fn nmda_mg_block_half_voltage(mg_conc_mm: f64) -> f64 {
    if !mg_conc_mm.is_finite() || mg_conc_mm <= 0.0 {
        return f64::NAN;
    }
    (mg_conc_mm / 3.57).ln() / 0.062
}

/// The **NMDA Mg²⁺ block voltage** `V = ln( (f/(1−f)) · [Mg²⁺]/3.57 ) / 0.062` (mV) —
/// the general inverse of the Jahr–Stevens block [`nmda_mg_block`]: the membrane
/// potential at which the *unblocked* fraction reaches an arbitrary target
/// `unblock_fraction` `B = f ∈ (0, 1)`, for extracellular magnesium `mg_conc_mm`
/// `[Mg²⁺]` (mM). It generalises [`nmda_mg_block_half_voltage`] (the `f = 0.5` case,
/// where `f/(1−f) = 1` collapses the formula to `ln([Mg²⁺]/3.57)/0.062`) to any
/// unblock level — it answers "how depolarised must the cell be for the NMDA
/// receptor to be `f`-relieved?".
///
/// `V` rises with both the demanded unblock fraction (more relief needs more
/// depolarisation) and the magnesium concentration (more block to clear). As
/// `f → 0` the channel is only fully blocked at `V → −∞`, and as `f → 1` full
/// relief needs `V → +∞`; both limits are non-physical. Returns `NaN` for a
/// fraction outside the open interval (`f ≤ 0` or `f ≥ 1`, or non-finite) or for
/// non-physical magnesium (non-finite or non-positive `[Mg²⁺]`).
pub fn nmda_mg_block_voltage(unblock_fraction: f64, mg_conc_mm: f64) -> f64 {
    if !unblock_fraction.is_finite()
        || unblock_fraction <= 0.0
        || unblock_fraction >= 1.0
        || !mg_conc_mm.is_finite()
        || mg_conc_mm <= 0.0
    {
        return f64::NAN;
    }
    ((unblock_fraction / (1.0 - unblock_fraction)) * (mg_conc_mm / 3.57)).ln() / 0.062
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn nmda_mg_block_half_voltage_relieves_half_the_block() {
        // Round-trip: at V½ the Jahr–Stevens block is exactly half-relieved (B = 0.5).
        for &mg in &[0.5_f64, 1.0, 1.5, 3.57] {
            let v_half = nmda_mg_block_half_voltage(mg);
            assert!((nmda_mg_block(v_half, mg) - 0.5).abs() < 1e-12, "B(V½) = 0.5 at [Mg]={mg}");
        }

        // Worked: V½ = 0 at [Mg] = 3.57 mM, ≈ −20.5 mV at the physiological 1 mM.
        assert!((nmda_mg_block_half_voltage(3.57) - 0.0).abs() < 1e-12, "V½(3.57) = 0");
        assert!((nmda_mg_block_half_voltage(1.0) - (-20.526)).abs() < 1e-2, "V½(1 mM) ≈ −20.5 mV");

        // More magnesium needs more depolarisation to half-relieve.
        assert!(
            nmda_mg_block_half_voltage(2.0) > nmda_mg_block_half_voltage(1.0),
            "V½ rises with [Mg²⁺]"
        );

        // Non-physical [Mg²⁺] → NaN (0 mV is a valid output, so NaN, not 0).
        assert!(nmda_mg_block_half_voltage(0.0).is_nan());
        assert!(nmda_mg_block_half_voltage(-1.0).is_nan());
        assert!(nmda_mg_block_half_voltage(f64::NAN).is_nan());
    }

    #[test]
    fn nmda_mg_block_voltage_inverts_the_block() {
        // (a) ROUND-TRIP (non-tautological): feeding V(f) back into the forward
        // block recovers the demanded unblock fraction, across (f, [Mg]) pairs.
        for &(f, mg) in &[(0.2_f64, 1.0_f64), (0.5, 1.0), (0.7, 2.0), (0.9, 0.5)] {
            let v = nmda_mg_block_voltage(f, mg);
            assert!(
                (nmda_mg_block(v, mg) - f).abs() < 1e-9,
                "B(V(f)) = f at f={f}, [Mg]={mg}"
            );
        }

        // (b) SPECIALIZATION: at f = 0.5 it collapses to the half-block voltage —
        // cross-checks the existing nmda_mg_block_half_voltage.
        for &mg in &[0.5_f64, 1.0, 3.57] {
            assert!(
                (nmda_mg_block_voltage(0.5, mg) - nmda_mg_block_half_voltage(mg)).abs() < 1e-12,
                "V(0.5, [Mg]) = V½([Mg]) at [Mg]={mg}"
            );
        }

        // (c) WORKED INDEPENDENT: at [Mg] = 3.57 the factor [Mg]/3.57 = 1, so
        // V = ln(f/(1−f))/0.062. f = 0.5 → ln(1)/0.062 = 0 mV; f = e/(1+e) →
        // f/(1−f) = e → V = ln(e)/0.062 = 1/0.062 ≈ 16.129 mV, both direct.
        assert!((nmda_mg_block_voltage(0.5, 3.57) - 0.0).abs() < 1e-12, "V(0.5, 3.57) = 0");
        let f_e = std::f64::consts::E / (1.0 + std::f64::consts::E);
        let v_e = nmda_mg_block_voltage(f_e, 3.57);
        assert!(
            (v_e - 1.0 / 0.062).abs() / (1.0 / 0.062) < 1e-9,
            "V(e/(1+e), 3.57) = 1/0.062 ≈ 16.129 mV, got {v_e}"
        );

        // (d) MONOTONICITY: a higher demanded unblock fraction needs more
        // depolarisation, and so does more magnesium at the same fraction.
        assert!(
            nmda_mg_block_voltage(0.7, 1.0) > nmda_mg_block_voltage(0.3, 1.0),
            "more unblock ⇒ more depolarised"
        );
        assert!(
            nmda_mg_block_voltage(0.5, 2.0) > nmda_mg_block_voltage(0.5, 1.0),
            "more [Mg²⁺] ⇒ more depolarised"
        );

        // (e) NaN guards: fraction must be strictly in (0, 1) and [Mg²⁺] > 0.
        assert!(nmda_mg_block_voltage(0.0, 1.0).is_nan());
        assert!(nmda_mg_block_voltage(1.0, 1.0).is_nan());
        assert!(nmda_mg_block_voltage(-0.1, 1.0).is_nan());
        assert!(nmda_mg_block_voltage(1.5, 1.0).is_nan());
        assert!(nmda_mg_block_voltage(f64::NAN, 1.0).is_nan());
        assert!(nmda_mg_block_voltage(0.5, 0.0).is_nan());
        assert!(nmda_mg_block_voltage(0.5, -1.0).is_nan());
        assert!(nmda_mg_block_voltage(0.5, f64::NAN).is_nan());
    }

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
    fn alpha_synapse_conductance_integral_is_g_max_tau_e() {
        use std::f64::consts::E;
        let (g_max, tau) = (5.0e-9, 2.0e-3); // 5 nS peak, 2 ms time-to-peak
        let area = alpha_synapse_conductance_integral(g_max, tau);
        // Closed form: the area under the alpha transient is exactly g_max·τ·e.
        assert!((area - g_max * tau * E).abs() / area < 1e-12, "closed form g_max·τ·e");
        // Numerical cross-check: midpoint-integrate the alpha conductance (#149)
        // out to 20τ — ties the closed form to the time course (non-tautological).
        let n = 100_000;
        let dt = 20.0 * tau / n as f64;
        let mut sum = 0.0;
        for k in 0..n {
            let t = (k as f64 + 0.5) * dt;
            sum += alpha_synapse_conductance(g_max, tau, t) * dt;
        }
        assert!((sum - area).abs() / area < 1e-3, "numerical ∫g dt {sum} ≈ {area}");
        // Linear in g_max and in τ.
        assert!(
            (alpha_synapse_conductance_integral(2.0 * g_max, tau) - 2.0 * area).abs() / area < 1e-12,
            "linear in g_max"
        );
        assert!(
            (alpha_synapse_conductance_integral(g_max, 3.0 * tau) - 3.0 * area).abs() / area < 1e-12,
            "linear in τ"
        );
        // Non-physical input → 0 (mirrors alpha_synapse_conductance).
        assert_eq!(alpha_synapse_conductance_integral(g_max, 0.0), 0.0);
        assert_eq!(alpha_synapse_conductance_integral(g_max, -1.0e-3), 0.0);
        assert_eq!(alpha_synapse_conductance_integral(f64::NAN, tau), 0.0);
        assert_eq!(alpha_synapse_conductance_integral(g_max, f64::INFINITY), 0.0);
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
        let tp = dual_exponential_peak_time(tau_r, tau_d);
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

    #[test]
    fn dual_exponential_peak_time_matches_the_closed_form() {
        let (tau_r, tau_d) = (0.0005, 0.003); // 0.5 ms rise, 3 ms decay
        let tp = dual_exponential_peak_time(tau_r, tau_d);
        // Closed form (τ_r·τ_d/(τ_d−τ_r))·ln(τ_d/τ_r).
        let expected = (tau_r * tau_d / (tau_d - tau_r)) * (tau_d / tau_r).ln();
        assert!((tp - expected).abs() < 1e-15, "closed form, got {tp}");
        // The rise time lies strictly between the two time constants.
        assert!(tp > tau_r && tp < tau_d, "τ_r < t_p < τ_d, got {tp}");
        // It is the time at which the conductance (#155) actually peaks at g_max.
        let g = dual_exponential_synapse_conductance(1.0, tau_r, tau_d, tp);
        assert!((g - 1.0).abs() < 1e-12, "g(t_p) = g_max, got {g}");
        // Symmetric in the two time constants (swapping τ_r and τ_d is invariant).
        assert!(
            (dual_exponential_peak_time(tau_d, tau_r) - tp).abs() < 1e-15,
            "symmetric in τ_r, τ_d"
        );
        // Degenerate / non-physical input → 0.
        assert_eq!(dual_exponential_peak_time(tau_r, tau_r), 0.0); // τ_r = τ_d (alpha limit)
        assert_eq!(dual_exponential_peak_time(0.0, tau_d), 0.0); // τ ≤ 0
        assert_eq!(dual_exponential_peak_time(tau_r, f64::NAN), 0.0); // non-finite
    }

    #[test]
    fn dual_exponential_conductance_integral_is_the_synaptic_efficacy() {
        let (g_max, tau_r, tau_d) = (1.0, 0.001, 0.005); // 1 ms rise, 5 ms decay
        let area = dual_exponential_conductance_integral(g_max, tau_r, tau_d);
        // Worked point: t_p ≈ 2.012 ms, norm ≈ 0.53499, ∫ = g_max·(τd−τr)/norm ≈ 7.4767e-3.
        assert!((area - 7.4767e-3).abs() / area < 1e-4, "efficacy ≈ 7.4767e-3 S·s, got {area}");
        // Linear in g_max (the peak conductance scales the whole waveform).
        assert!(
            (dual_exponential_conductance_integral(2.0 * g_max, tau_r, tau_d) - 2.0 * area).abs()
                / area
                < 1e-12,
            "linear in g_max"
        );
        // STRONG non-tautological cross-check: midpoint-integrate the actual
        // conductance #155 over [0, 60·τd] (the tail beyond is ~e^(−60), negligible).
        // The impl is a closed form; this is an independent Riemann sum of the fn.
        let dt = tau_d / 2000.0;
        let n = (60.0 * tau_d / dt) as usize;
        let mut numeric = 0.0;
        for k in 0..n {
            let t = (k as f64 + 0.5) * dt;
            numeric += dual_exponential_synapse_conductance(g_max, tau_r, tau_d, t) * dt;
        }
        assert!(
            (numeric - area).abs() / area < 1e-4,
            "∫g dt (numeric {numeric}) matches the closed form {area}"
        );
        // Degenerate τ_r = τ_d → 0 (use the alpha integral there).
        assert_eq!(dual_exponential_conductance_integral(g_max, tau_r, tau_r), 0.0);
        // Non-physical input → 0.
        assert_eq!(dual_exponential_conductance_integral(g_max, 0.0, tau_d), 0.0); // τ ≤ 0
        assert_eq!(dual_exponential_conductance_integral(g_max, tau_r, f64::NAN), 0.0); // non-finite τ
        assert_eq!(dual_exponential_conductance_integral(f64::NAN, tau_r, tau_d), 0.0); // non-finite g
    }

    #[test]
    fn nmda_mg_block_relieves_with_depolarisation() {
        let mg = 1.0; // 1 mM physiological extracellular magnesium
        // At rest (V = −80 mV) the channel is deeply blocked.
        let rest = nmda_mg_block(-80.0, mg);
        assert!(rest < 0.05, "blocked at rest: {rest}");
        assert!(
            (rest - 1.0 / (1.0 + (1.0 / 3.57) * (0.062_f64 * 80.0).exp())).abs() < 1e-12,
            "closed form at rest"
        );
        // Strong depolarisation (V = +40 mV) relieves the block.
        let depol = nmda_mg_block(40.0, mg);
        assert!(depol > 0.9, "relieved when depolarised: {depol}");
        // Monotonic in voltage: more depolarised → more unblocked.
        let mid = nmda_mg_block(0.0, mg);
        assert!(mid > rest && depol > mid, "monotone in V");
        // Always a valid open fraction in (0, 1].
        for v in [-100.0_f64, -60.0, -20.0, 0.0, 20.0, 60.0] {
            let b = nmda_mg_block(v, mg);
            assert!(b > 0.0 && b <= 1.0, "B in (0,1] at V={v}: {b}");
        }
        // No magnesium → no block: B = 1 at any voltage (from the formula itself).
        assert!((nmda_mg_block(-80.0, 0.0) - 1.0).abs() < 1e-12, "Mg=0 → open at rest");
        assert!((nmda_mg_block(40.0, 0.0) - 1.0).abs() < 1e-12, "Mg=0 → open depolarised");
        // Heavier magnesium deepens the block at a fixed voltage.
        assert!(
            nmda_mg_block(-30.0, 100.0) < nmda_mg_block(-30.0, 1.0),
            "more Mg → more block"
        );
        // Non-physical input → 0.
        assert_eq!(nmda_mg_block(f64::NAN, mg), 0.0);
        assert_eq!(nmda_mg_block(-80.0, -1.0), 0.0);
        assert_eq!(nmda_mg_block(-80.0, f64::INFINITY), 0.0);
    }
}
