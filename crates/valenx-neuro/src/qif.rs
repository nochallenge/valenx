//! Quadratic integrate-and-fire (QIF) neuron — the canonical Type-I model.
//!
//! The QIF is the normal form of the saddle-node-on-invariant-circle (SNIC)
//! bifurcation that marks the onset of repetitive firing in every Type-I
//! neuron. In normalized units the sub-threshold dynamics are
//!
//! ```text
//! τ dV/dt = V² + I,   reset:  V → V_reset  when  V ≥ V_peak
//! ```
//!
//! - **I < 0** — two fixed points at `V = ±√(−I)`: a stable rest at `−√(−I)`
//!   and an unstable threshold at `+√(−I)`. The cell is excitable but
//!   quiescent.
//! - **I = 0** — the two collide at `V = 0` (the saddle-node); the firing
//!   threshold.
//! - **I > 0** — no fixed point; `V` runs to `+∞` in finite time and the cell
//!   fires periodically, with the exact period
//!
//! ```text
//! T = τ ∫_{V_reset}^{V_peak} dV/(V²+I)
//!   = (τ/√I)·[ arctan(V_peak/√I) − arctan(V_reset/√I) ].
//! ```
//!
//! The firing rate `f = 1/T` therefore has a closed form, and as the thresholds
//! `V_peak, V_reset → ±∞` it approaches `f → √I/(τπ)` — the hallmark Type-I
//! `f ∝ √I` onset (a continuous f–I curve rising from zero at `I = 0`), unlike
//! the discontinuous Type-II onset of e.g. Hodgkin–Huxley.
//!
//! Reference: Ermentrout & Kopell (1986); Latham et al., *J. Neurophysiol.*
//! (2000); Izhikevich, *Dynamical Systems in Neuroscience* (2007), ch. 8 (the
//! QIF / theta-neuron normal form).
//!
//! # Honest scope
//!
//! Research/educational. A phenomenological canonical model in normalized
//! units — it captures the *qualitative* Type-I onset structure exactly, not a
//! specific cell's quantitative biophysics; use Hodgkin–Huxley
//! ([`crate::cable`]) for the latter.

/// Parameters of a normalized quadratic integrate-and-fire neuron.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct QifParams {
    /// Membrane time constant `τ` (`> 0`; sets the absolute time-scale).
    pub tau: f64,
    /// Reset potential `V_reset` applied after a spike.
    pub v_reset: f64,
    /// Peak / spike-detection threshold `V_peak` (the numerical stand-in for
    /// the `+∞` blow-up).
    pub v_peak: f64,
}

impl Default for QifParams {
    fn default() -> Self {
        Self {
            tau: 1.0,
            v_reset: -10.0,
            v_peak: 10.0,
        }
    }
}

/// The fixed points of the QIF at constant input `i`.
///
/// Returns `Some((rest, threshold)) = (−√(−i), +√(−i))` for `i < 0`, and `None`
/// for `i ≥ 0` — at `i = 0` they have merged at the origin (the saddle-node)
/// and the cell sits at the firing threshold; for `i > 0` no fixed point
/// exists.
pub fn qif_fixed_points(i: f64) -> Option<(f64, f64)> {
    if i < 0.0 {
        let r = (-i).sqrt();
        Some((-r, r))
    } else {
        None
    }
}

/// The exact analytic firing rate `f = 1/T` for constant input `i`, or `None`
/// when the cell does not fire repetitively (`i ≤ 0`) or the parameters are
/// degenerate (`τ ≤ 0`, `V_peak ≤ V_reset`).
pub fn qif_firing_rate(i: f64, params: &QifParams) -> Option<f64> {
    if i <= 0.0 || params.tau <= 0.0 || params.v_peak <= params.v_reset {
        return None;
    }
    let s = i.sqrt();
    let period = (params.tau / s) * ((params.v_peak / s).atan() - (params.v_reset / s).atan());
    Some(1.0 / period)
}

/// One QIF neuron's state.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct QifNeuron {
    /// Membrane potential `V`.
    pub v: f64,
    /// The neuron's parameters.
    pub params: QifParams,
}

impl QifNeuron {
    /// A neuron initialised at its reset potential `V_reset`.
    pub fn new(params: QifParams) -> Self {
        Self {
            v: params.v_reset,
            params,
        }
    }

    /// A neuron initialised at an explicit potential `v0`.
    pub fn at(v0: f64, params: QifParams) -> Self {
        Self { v: v0, params }
    }

    /// Advance the state by `dt` with constant input `current`, using a
    /// 4th-order Runge–Kutta step on `τ dV/dt = V² + I`. Returns `true` if the
    /// neuron crossed `V_peak` this step (and was reset to `V_reset`).
    pub fn step(&mut self, dt: f64, current: f64) -> bool {
        let tau = self.params.tau;
        let f = |v: f64| (v * v + current) / tau;
        let k1 = f(self.v);
        let k2 = f(self.v + 0.5 * dt * k1);
        let k3 = f(self.v + 0.5 * dt * k2);
        let k4 = f(self.v + dt * k3);
        self.v += (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
        if !self.v.is_finite() || self.v >= self.params.v_peak {
            self.v = self.params.v_reset;
            true
        } else {
            false
        }
    }
}

/// Simulate a QIF neuron from `v0` under constant `current` and return the spike
/// times. `dt` is the integration step, `duration` the total simulated time.
pub fn qif_spike_times(
    v0: f64,
    params: QifParams,
    current: f64,
    dt: f64,
    duration: f64,
) -> Vec<f64> {
    let mut neuron = QifNeuron::at(v0, params);
    let mut spikes = Vec::new();
    let steps = (duration / dt).round() as usize;
    for i in 0..steps {
        if neuron.step(dt, current) {
            spikes.push(i as f64 * dt);
        }
    }
    spikes
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn fixed_points_are_plus_minus_sqrt_minus_i() {
        let (rest, thr) = qif_fixed_points(-1.0).unwrap();
        assert!(close(rest, -1.0, 1e-12) && close(thr, 1.0, 1e-12));
        let (rest, thr) = qif_fixed_points(-4.0).unwrap();
        assert!(close(rest, -2.0, 1e-12) && close(thr, 2.0, 1e-12));
        // At and above the saddle-node there is no fixed point.
        assert!(qif_fixed_points(0.0).is_none());
        assert!(qif_fixed_points(1.0).is_none());
    }

    #[test]
    fn analytic_firing_rate_closed_form() {
        // Default params, I = 1: T = 2·arctan(10) = 2.942255, f = 0.339875.
        let f = qif_firing_rate(1.0, &QifParams::default()).unwrap();
        assert!(close(f, 0.339_875, 1e-5), "f = {f}");
        // No repetitive firing at or below the saddle-node.
        assert!(qif_firing_rate(0.0, &QifParams::default()).is_none());
        assert!(qif_firing_rate(-1.0, &QifParams::default()).is_none());
    }

    #[test]
    fn simulated_rate_matches_analytic() {
        // The RK4 integrator's steady-state rate must match the exact closed
        // form. Default params, I = 1, fine step, long enough for many spikes.
        let params = QifParams::default();
        let analytic = qif_firing_rate(1.0, &params).unwrap();
        let spikes = qif_spike_times(params.v_reset, params, 1.0, 1e-3, 90.0);
        assert!(spikes.len() >= 10, "only {} spikes", spikes.len());
        // Steady-state mean inter-spike interval (drops transient + edge bias).
        let n = spikes.len();
        let mean_isi = (spikes[n - 1] - spikes[0]) / (n - 1) as f64;
        let sim_rate = 1.0 / mean_isi;
        assert!(
            close(sim_rate, analytic, 0.01 * analytic),
            "sim {sim_rate} vs analytic {analytic}"
        );
    }

    #[test]
    fn type_one_sqrt_onset() {
        // With wide thresholds f → √I/(τπ): a continuous f∝√I onset from zero.
        let wide = QifParams {
            tau: 1.0,
            v_reset: -1.0e4,
            v_peak: 1.0e4,
        };
        let f_small = qif_firing_rate(0.01, &wide).unwrap();
        assert!(
            close(f_small, 0.1 / std::f64::consts::PI, 1e-4),
            "f(0.01) = {f_small}"
        );
        // Quadrupling I doubles √I, hence (in the wide limit) doubles the rate.
        let ratio = qif_firing_rate(0.16, &wide).unwrap() / qif_firing_rate(0.04, &wide).unwrap();
        assert!(close(ratio, 2.0, 1e-2), "f(4I)/f(I) = {ratio}");
    }

    #[test]
    fn excitable_below_threshold_is_quiescent_above_fires_once() {
        // I = -1 → rest -1, threshold +1. Started below threshold the cell
        // relaxes to rest and never fires; started above it fires exactly once
        // (the reset to -10 < rest leaves it on the resting branch).
        let params = QifParams::default();
        let quiet = qif_spike_times(0.5, params, -1.0, 1e-3, 50.0);
        assert!(
            quiet.is_empty(),
            "expected silence, got {} spikes",
            quiet.len()
        );
        let fired = qif_spike_times(2.0, params, -1.0, 1e-3, 50.0);
        assert_eq!(
            fired.len(),
            1,
            "expected a single spike, got {}",
            fired.len()
        );
    }

    #[test]
    fn degenerate_params_have_no_rate() {
        let bad_tau = QifParams {
            tau: 0.0,
            ..Default::default()
        };
        assert!(qif_firing_rate(1.0, &bad_tau).is_none());
        let bad_thr = QifParams {
            tau: 1.0,
            v_reset: 5.0,
            v_peak: 5.0,
        };
        assert!(qif_firing_rate(1.0, &bad_thr).is_none());
    }
}
