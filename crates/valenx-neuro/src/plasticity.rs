//! Spike-timing-dependent plasticity (STDP).
//!
//! Where [`crate::synapse`] gives the *conductance* a synaptic event produces,
//! this module gives how the synaptic *weight itself* changes when pre- and
//! post-synaptic spikes are paired — the canonical activity-dependent learning
//! rule behind Hebbian potentiation and competitive synaptic refinement.
//!
//! The **pair-based STDP rule** (Bi & Poo 1998; Song, Miller & Abbott 2000)
//! makes the change depend on the spike-time difference `Δt = t_post − t_pre`
//! through two decaying exponentials:
//!
//! ```text
//! Δw(Δt) =  A₊·e^(−Δt/τ₊)   for Δt > 0   (pre before post → potentiation, LTP)
//! Δw(Δt) = −A₋·e^(+Δt/τ₋)   for Δt < 0   (post before pre → depression,  LTD)
//! ```
//!
//! The window is **causal and asymmetric**: a presynaptic spike that *precedes*
//! the postsynaptic one (it could have helped cause it) strengthens the synapse,
//! while the reverse order weakens it, with an influence that decays over the
//! tens-of-milliseconds correlation windows `τ₊`, `τ₋`.

/// Parameters of the **pair-based spike-timing-dependent plasticity (STDP)**
/// rule (Bi & Poo 1998; Song, Miller & Abbott 2000). See the
/// [module documentation](crate::plasticity) for the window equation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StdpRule {
    /// LTP amplitude `A₊` — the maximal synaptic *strengthening*, approached as
    /// `Δt → 0⁺` (pre just before post). Dimensionless weight units.
    pub a_plus: f64,
    /// LTD amplitude `A₋` (a positive magnitude) — the maximal *weakening*,
    /// approached as `Δt → 0⁻` (post just before pre).
    pub a_minus: f64,
    /// LTP time constant `τ₊` (s) — the pre-before-post correlation window over
    /// which potentiation decays (typically ≈ 20 ms).
    pub tau_plus_s: f64,
    /// LTD time constant `τ₋` (s) — the post-before-pre correlation window over
    /// which depression decays (typically ≈ 20 ms).
    pub tau_minus_s: f64,
}

impl StdpRule {
    /// The **STDP weight change** `Δw` for a pre/post spike pair separated by
    /// `dt_s` `Δt = t_post − t_pre` (s).
    ///
    /// Positive for `Δt > 0` (pre before post → LTP, decaying from `A₊` with
    /// time constant `τ₊`), negative for `Δt < 0` (post before pre → LTD,
    /// decaying from `−A₋` with `τ₋`), and `0` for exactly coincident spikes
    /// (`Δt = 0`, the discontinuous limit). Returns `0` for a non-finite `Δt` or
    /// non-physical parameters (`τ ≤ 0`, or any non-finite field).
    pub fn weight_change(&self, dt_s: f64) -> f64 {
        if !dt_s.is_finite()
            || !self.a_plus.is_finite()
            || !self.a_minus.is_finite()
            || !self.tau_plus_s.is_finite()
            || !self.tau_minus_s.is_finite()
            || self.tau_plus_s <= 0.0
            || self.tau_minus_s <= 0.0
        {
            return 0.0;
        }
        if dt_s > 0.0 {
            self.a_plus * (-dt_s / self.tau_plus_s).exp()
        } else if dt_s < 0.0 {
            -self.a_minus * (dt_s / self.tau_minus_s).exp()
        } else {
            0.0
        }
    }

    /// The **total potentiation area** `∫₀^∞ A₊·e^(−Δt/τ₊) dΔt = A₊·τ₊` — the
    /// integral of the LTP lobe of the [`weight_change`](Self::weight_change)
    /// window.
    ///
    /// With [`depression_integral`](Self::depression_integral) it fixes the
    /// rule's behaviour under *uncorrelated* pre/post activity (uniformly
    /// distributed `Δt`): the net drift is `A₊·τ₊ − A₋·τ₋`, so when
    /// `A₋·τ₋ > A₊·τ₊` uncorrelated firing nets *depression* — the
    /// Song–Miller–Abbott condition that keeps weights competitive and bounded.
    /// Returns `0` for non-physical parameters (`τ₊ ≤ 0` or non-finite).
    pub fn potentiation_integral(&self) -> f64 {
        if !self.a_plus.is_finite() || !self.tau_plus_s.is_finite() || self.tau_plus_s <= 0.0 {
            return 0.0;
        }
        self.a_plus * self.tau_plus_s
    }

    /// The **total depression area** `∫_{−∞}^0 A₋·e^(+Δt/τ₋) dΔt = A₋·τ₋` — the
    /// integral of the LTD lobe of the [`weight_change`](Self::weight_change)
    /// window. See [`potentiation_integral`](Self::potentiation_integral) for the
    /// net-balance / stability interpretation. Returns `0` for non-physical
    /// parameters (`τ₋ ≤ 0` or non-finite).
    pub fn depression_integral(&self) -> f64 {
        if !self.a_minus.is_finite() || !self.tau_minus_s.is_finite() || self.tau_minus_s <= 0.0 {
            return 0.0;
        }
        self.a_minus * self.tau_minus_s
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::E;

    fn rule() -> StdpRule {
        // Classic Song–Miller–Abbott / Bi–Poo parameters: 20 ms windows, with a
        // slightly stronger depression lobe (A₋ > A₊) for competitive stability.
        StdpRule {
            a_plus: 0.010,
            a_minus: 0.012,
            tau_plus_s: 0.020,
            tau_minus_s: 0.020,
        }
    }

    #[test]
    fn stdp_window_is_causal_and_asymmetric() {
        let r = rule();
        // Pre before post (Δt > 0) strengthens; post before pre (Δt < 0) weakens.
        assert!(r.weight_change(0.005) > 0.0, "pre-before-post potentiates");
        assert!(r.weight_change(-0.005) < 0.0, "post-before-pre depresses");
        // Coincident spikes: the discontinuous limit is conventionally 0.
        assert_eq!(r.weight_change(0.0), 0.0);
        // The window decays monotonically in |Δt| on each side.
        assert!(
            r.weight_change(0.002) > r.weight_change(0.010),
            "LTP decays with lag"
        );
        assert!(
            r.weight_change(-0.002) < r.weight_change(-0.010),
            "LTD decays with lag"
        );
    }

    #[test]
    fn stdp_peaks_approach_the_amplitudes_at_zero_lag() {
        let r = rule();
        // As Δt → 0⁺, Δw → A₊; as Δt → 0⁻, Δw → −A₋.
        assert!(
            (r.weight_change(1e-9) - r.a_plus).abs() < 1e-6,
            "LTP peak → A₊"
        );
        assert!(
            (r.weight_change(-1e-9) + r.a_minus).abs() < 1e-6,
            "LTD peak → −A₋"
        );
    }

    #[test]
    fn stdp_window_matches_the_exponential_form_at_one_time_constant() {
        // GROUND TRUTH: at Δt = τ the exponential window has decayed to exactly
        // 1/e of its zero-lag amplitude — e is universal, so this pins the
        // exponential shape against an independent constant, not the code.
        let r = rule();
        let ltp_at_tau = r.weight_change(r.tau_plus_s);
        assert!(
            (ltp_at_tau - r.a_plus / E).abs() < 1e-15,
            "Δw(τ₊) = A₊/e = {}, got {ltp_at_tau}",
            r.a_plus / E
        );
        let ltd_at_tau = r.weight_change(-r.tau_minus_s);
        assert!(
            (ltd_at_tau + r.a_minus / E).abs() < 1e-15,
            "Δw(−τ₋) = −A₋/e = {}, got {ltd_at_tau}",
            -r.a_minus / E
        );
        // Worked absolute values: A₊=0.010, τ₊=20 ms → Δw(20 ms) = 0.010/e =
        // 0.0036787944; Δw(10 ms) = 0.010·e^(−0.5) = 0.0060653066.
        assert!((r.weight_change(0.020) - 0.003_678_794_4).abs() < 1e-9);
        assert!((r.weight_change(0.010) - 0.006_065_306_6).abs() < 1e-9);
    }

    #[test]
    fn stdp_lobe_integrals_are_amplitude_times_tau() {
        // GROUND TRUTH: ∫₀^∞ A·e^(−Δt/τ) dΔt = A·τ (independent calculus), tied
        // to the window itself by a numerical Riemann sum — non-tautological.
        let r = rule();
        assert!(
            (r.potentiation_integral() - r.a_plus * r.tau_plus_s).abs() < 1e-18,
            "∫ LTP = A₊·τ₊"
        );
        assert!(
            (r.depression_integral() - r.a_minus * r.tau_minus_s).abs() < 1e-18,
            "∫ LTD = A₋·τ₋"
        );
        // Numerically integrate the LTP lobe out to 30·τ₊ and match A₊·τ₊.
        let n = 300_000;
        let dt = 30.0 * r.tau_plus_s / n as f64;
        let mut sum = 0.0;
        for k in 0..n {
            let t = (k as f64 + 0.5) * dt;
            sum += r.weight_change(t) * dt;
        }
        let exact = r.potentiation_integral();
        assert!(
            (sum - exact).abs() / exact < 1e-3,
            "numerical ∫ LTP {sum} ≈ A₊·τ₊ {exact}"
        );
    }

    #[test]
    fn stdp_net_balance_follows_the_lobe_integrals() {
        // Song–Miller–Abbott competitive stability: uncorrelated activity nets
        // depression when A₋·τ₋ > A₊·τ₊. This rule (A₋=0.012 > A₊=0.010, equal τ)
        // is net-depressing.
        let r = rule();
        assert!(
            r.depression_integral() > r.potentiation_integral(),
            "this rule net-depresses uncorrelated pairs"
        );
        // A balanced rule (A₊·τ₊ = A₋·τ₋) has exactly zero net drift.
        let balanced = StdpRule {
            a_minus: 0.010,
            ..r
        };
        assert!(
            (balanced.potentiation_integral() - balanced.depression_integral()).abs() < 1e-18,
            "balanced rule has zero net drift"
        );
    }

    #[test]
    fn stdp_guards_non_physical_input() {
        let r = rule();
        assert_eq!(r.weight_change(f64::NAN), 0.0);
        assert_eq!(r.weight_change(f64::INFINITY), 0.0);
        // Non-physical parameters → 0 (no panic, no NaN).
        let bad_tau_plus = StdpRule {
            tau_plus_s: 0.0,
            ..r
        };
        assert_eq!(bad_tau_plus.weight_change(0.005), 0.0);
        assert_eq!(bad_tau_plus.potentiation_integral(), 0.0);
        let bad_tau_minus = StdpRule {
            tau_minus_s: -1.0,
            ..r
        };
        assert_eq!(bad_tau_minus.weight_change(-0.005), 0.0);
        assert_eq!(bad_tau_minus.depression_integral(), 0.0);
        let nan_amp = StdpRule {
            a_plus: f64::NAN,
            ..r
        };
        assert_eq!(nan_amp.weight_change(0.005), 0.0);
    }
}
