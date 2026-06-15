//! The Izhikevich (2003) simple spiking-neuron model.
//!
//! A two-variable reduction that is far cheaper than Hodgkin–Huxley yet
//! reproduces the major cortical firing patterns (regular spiking, bursting,
//! chattering, fast spiking, …) by changing four parameters. The dynamics are
//!
//! ```text
//! v' = 0.04 v² + 5 v + 140 − u + I
//! u' = a (b v − u)
//! if v ≥ 30 mV:  v ← c,  u ← u + d
//! ```
//!
//! with `v` the membrane potential (mV), `u` a recovery variable, and `I` the
//! input current. The peak at +30 mV is the spike; `c` is the reset potential
//! and `d` the after-spike recovery jump. Integration uses Izhikevich's own
//! numerically-stable scheme — two half-steps on `v` per step, then `u`.
//!
//! Reference: E. M. Izhikevich, "Simple Model of Spiking Neurons", *IEEE Trans.
//! Neural Networks* 14(6), 2003. The canonical parameter sets below are from
//! that paper's Fig. 2.
//!
//! # Honest scope
//!
//! Research/educational. This is a phenomenological reduced model: `v` is in
//! the paper's mV-like units and the quadratic coefficients are fixed
//! constants, not fitted to a specific cell. It captures spike *patterns* and
//! relative timing, not quantitative sub-threshold biophysics — use
//! Hodgkin–Huxley ([`crate::cable`]) for the latter.

/// The four Izhikevich parameters.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct IzhikevichParams {
    /// Recovery time-scale (larger = faster recovery).
    pub a: f64,
    /// Recovery sensitivity to sub-threshold `v`.
    pub b: f64,
    /// After-spike reset value of `v` (mV).
    pub c: f64,
    /// After-spike jump added to `u`.
    pub d: f64,
}

impl IzhikevichParams {
    /// Regular-spiking cortical excitatory neuron (RS): `a=0.02, b=0.2, c=−65, d=8`.
    pub fn regular_spiking() -> Self {
        Self {
            a: 0.02,
            b: 0.2,
            c: -65.0,
            d: 8.0,
        }
    }

    /// Intrinsically-bursting neuron (IB): `a=0.02, b=0.2, c=−55, d=4`.
    pub fn intrinsically_bursting() -> Self {
        Self {
            a: 0.02,
            b: 0.2,
            c: -55.0,
            d: 4.0,
        }
    }

    /// Chattering neuron (CH): `a=0.02, b=0.2, c=−50, d=2`.
    pub fn chattering() -> Self {
        Self {
            a: 0.02,
            b: 0.2,
            c: -50.0,
            d: 2.0,
        }
    }

    /// Fast-spiking inhibitory neuron (FS): `a=0.1, b=0.2, c=−65, d=2`.
    pub fn fast_spiking() -> Self {
        Self {
            a: 0.1,
            b: 0.2,
            c: -65.0,
            d: 2.0,
        }
    }

    /// Low-threshold-spiking neuron (LTS): `a=0.02, b=0.25, c=−65, d=2`.
    pub fn low_threshold_spiking() -> Self {
        Self {
            a: 0.02,
            b: 0.25,
            c: -65.0,
            d: 2.0,
        }
    }

    /// The stable resting potential (mV) with `I = 0`: the lower root of the
    /// sub-threshold balance `0.04 v² + (5 − b) v + 140 = 0` (using `u = b v`),
    /// or `None` when there is no real fixed point.
    pub fn resting_potential(&self) -> Option<f64> {
        let p = 5.0 - self.b;
        let discriminant = p * p - 4.0 * 0.04 * 140.0;
        if discriminant < 0.0 {
            return None;
        }
        // The more-negative (stable-node) root.
        Some((-p - discriminant.sqrt()) / (2.0 * 0.04))
    }
}

/// One Izhikevich neuron's state.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct IzhikevichNeuron {
    /// Membrane potential (mV).
    pub v: f64,
    /// Recovery variable.
    pub u: f64,
    /// The neuron's parameters.
    pub params: IzhikevichParams,
}

impl IzhikevichNeuron {
    /// A neuron initialised at its resting fixed point (`v = resting_potential`,
    /// falling back to `−70` mV if none exists; `u = b v`).
    pub fn at_rest(params: IzhikevichParams) -> Self {
        let v = params.resting_potential().unwrap_or(-70.0);
        Self {
            v,
            u: params.b * v,
            params,
        }
    }

    /// Advance the state by `dt_ms` with constant input current `current`,
    /// using Izhikevich's two-half-step scheme for `v`. Returns `true` if the
    /// neuron spiked this step (and was reset).
    pub fn step(&mut self, dt_ms: f64, current: f64) -> bool {
        let half = 0.5 * dt_ms;
        for _ in 0..2 {
            self.v += half * (0.04 * self.v * self.v + 5.0 * self.v + 140.0 - self.u + current);
        }
        self.u += dt_ms * self.params.a * (self.params.b * self.v - self.u);
        if self.v >= 30.0 {
            self.v = self.params.c;
            self.u += self.params.d;
            true
        } else {
            false
        }
    }
}

/// Simulate a single neuron under constant `current` and return the spike times
/// (ms). Starts from rest. `dt_ms` is the integration step (1.0 or 0.5 ms are
/// the usual choices); `duration_ms` is the total simulated time.
pub fn simulate_spike_times(
    params: IzhikevichParams,
    current: f64,
    dt_ms: f64,
    duration_ms: f64,
) -> Vec<f64> {
    let mut neuron = IzhikevichNeuron::at_rest(params);
    let mut spikes = Vec::new();
    let steps = (duration_ms / dt_ms).round() as usize;
    for i in 0..steps {
        if neuron.step(dt_ms, current) {
            spikes.push(i as f64 * dt_ms);
        }
    }
    spikes
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn regular_spiking_resting_potential_is_minus_70() {
        // Analytic lower root of 0.04 v² + 4.8 v + 140 = 0 is exactly −70.
        let rest = IzhikevichParams::regular_spiking()
            .resting_potential()
            .unwrap();
        assert!(
            (rest - (-70.0)).abs() < 1e-9,
            "RS rest {rest} should be −70"
        );
    }

    #[test]
    fn lts_resting_potential_matches_quadratic() {
        // b = 0.25 -> root [−4.75 − sqrt(4.75² − 22.4)] / 0.08 ≈ −64.40.
        let rest = IzhikevichParams::low_threshold_spiking()
            .resting_potential()
            .unwrap();
        assert!(
            (rest - (-64.40)).abs() < 0.05,
            "LTS rest {rest} should be ≈ −64.4"
        );
    }

    #[test]
    fn rest_is_a_true_fixed_point() {
        // At rest with no input, the neuron must not drift or spike.
        let mut n = IzhikevichNeuron::at_rest(IzhikevichParams::regular_spiking());
        let v0 = n.v;
        let mut spikes = 0;
        for _ in 0..2000 {
            if n.step(1.0, 0.0) {
                spikes += 1;
            }
        }
        assert_eq!(spikes, 0, "a resting neuron must not spike");
        assert!((n.v - v0).abs() < 1e-6, "v drifted from {v0} to {}", n.v);
    }

    #[test]
    fn supra_threshold_drive_produces_regular_spiking() {
        // RS neuron with I = 10 fires repeatedly over 1 s.
        let spikes = simulate_spike_times(IzhikevichParams::regular_spiking(), 10.0, 0.5, 1000.0);
        assert!(
            spikes.len() > 5,
            "RS at I=10 should spike repeatedly, got {}",
            spikes.len()
        );
    }

    #[test]
    fn firing_rate_is_monotonic_in_current() {
        let p = IzhikevichParams::regular_spiking();
        let n_lo = simulate_spike_times(p, 3.0, 0.5, 1000.0).len();
        let n_mid = simulate_spike_times(p, 8.0, 0.5, 1000.0).len();
        let n_hi = simulate_spike_times(p, 15.0, 0.5, 1000.0).len();
        assert!(
            n_hi >= n_mid && n_mid >= n_lo,
            "f-I should be monotonic: I=3 -> {n_lo}, I=8 -> {n_mid}, I=15 -> {n_hi}"
        );
    }

    #[test]
    fn spike_resets_v_to_c() {
        let mut n = IzhikevichNeuron::at_rest(IzhikevichParams::regular_spiking());
        // Drive hard until the first spike, then check the reset.
        for _ in 0..100_000 {
            if n.step(0.5, 20.0) {
                assert!(
                    (n.v - n.params.c).abs() < 1e-9,
                    "after a spike v should reset to c={}, got {}",
                    n.params.c,
                    n.v
                );
                return;
            }
        }
        panic!("neuron never spiked under strong drive");
    }
}
