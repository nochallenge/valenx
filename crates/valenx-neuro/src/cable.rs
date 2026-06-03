//! Hodgkin–Huxley membrane dynamics.
//!
//! A single-compartment model (this task) of an axon membrane, integrated
//! with RK4, using the classic Hodgkin–Huxley (1952) squid-giant-axon
//! parameters in the modern convention (resting potential −65 mV). The
//! multi-compartment cable is added on top in the next task.
//!
//! Native units are the standard HH system, which is internally consistent
//! so the membrane equations need no conversion: potential **mV**, time
//! **ms**, capacitance **µF/cm²**, conductance density **mS/cm²**, current
//! density **µA/cm²**.

// --- pinned HH 1952 constants (modern −65 mV convention) ----------------

/// Membrane capacitance (µF/cm²).
pub const C_M: f64 = 1.0;
/// Maximal sodium conductance (mS/cm²).
pub const G_NA: f64 = 120.0;
/// Maximal potassium conductance (mS/cm²).
pub const G_K: f64 = 36.0;
/// Leak conductance (mS/cm²).
pub const G_L: f64 = 0.3;
/// Sodium reversal potential (mV).
pub const E_NA: f64 = 50.0;
/// Potassium reversal potential (mV).
pub const E_K: f64 = -77.0;
/// Leak reversal potential (mV).
pub const E_L: f64 = -54.4;
/// Resting membrane potential (mV).
pub const V_REST: f64 = -65.0;

// --- gating rate functions (V in mV, rates in 1/ms) ---------------------
// The two `x / (1 − e^{−x/10})` forms are 0/0 at V = −40 (αm) and V = −55
// (αn); each is guarded with its analytic limit there.

fn alpha_m(v: f64) -> f64 {
    let x = v + 40.0;
    if x.abs() < 1e-6 {
        1.0 // limit of 0.1·x / (1 − e^{−x/10}) as x → 0
    } else {
        0.1 * x / (1.0 - (-x / 10.0).exp())
    }
}
fn beta_m(v: f64) -> f64 {
    4.0 * (-(v + 65.0) / 18.0).exp()
}
fn alpha_h(v: f64) -> f64 {
    0.07 * (-(v + 65.0) / 20.0).exp()
}
fn beta_h(v: f64) -> f64 {
    1.0 / (1.0 + (-(v + 35.0) / 10.0).exp())
}
fn alpha_n(v: f64) -> f64 {
    let x = v + 55.0;
    if x.abs() < 1e-6 {
        0.1 // limit of 0.01·x / (1 − e^{−x/10}) as x → 0
    } else {
        0.01 * x / (1.0 - (-x / 10.0).exp())
    }
}
fn beta_n(v: f64) -> f64 {
    0.125 * (-(v + 65.0) / 80.0).exp()
}

fn steady(alpha: f64, beta: f64) -> f64 {
    alpha / (alpha + beta)
}

/// A rectangular stimulating current pulse (intracellular, µA/cm²).
#[derive(Debug, Clone, Copy)]
pub struct StimPulse {
    /// Amplitude (µA/cm²).
    pub amp_ua_cm2: f64,
    /// Onset time (ms).
    pub start_ms: f64,
    /// Pulse width (ms).
    pub width_ms: f64,
}

impl StimPulse {
    /// A zero (off) pulse.
    pub fn off() -> Self {
        Self { amp_ua_cm2: 0.0, start_ms: 0.0, width_ms: 0.0 }
    }
    /// The stimulus current (µA/cm²) at time `t_ms`.
    fn current_at(&self, t_ms: f64) -> f64 {
        if t_ms >= self.start_ms && t_ms < self.start_ms + self.width_ms {
            self.amp_ua_cm2
        } else {
            0.0
        }
    }
}

/// A single Hodgkin–Huxley membrane compartment: potential `v` (mV) and the
/// three gating variables `m`, `h`, `n` (dimensionless, in [0, 1]).
#[derive(Debug, Clone, Copy)]
pub struct HhCompartment {
    /// Membrane potential (mV).
    pub v: f64,
    /// Sodium activation gate.
    pub m: f64,
    /// Sodium inactivation gate.
    pub h: f64,
    /// Potassium activation gate.
    pub n: f64,
}

impl HhCompartment {
    /// A compartment at rest: V = −65 mV with every gate at its steady-state
    /// value for that potential.
    pub fn at_rest() -> Self {
        let v = V_REST;
        Self {
            v,
            m: steady(alpha_m(v), beta_m(v)),
            h: steady(alpha_h(v), beta_h(v)),
            n: steady(alpha_n(v), beta_n(v)),
        }
    }

    /// Time derivatives `(dV, dm, dh, dn)` given an injected current
    /// `i_stim` (µA/cm²).
    fn derivs(v: f64, m: f64, h: f64, n: f64, i_stim: f64) -> (f64, f64, f64, f64) {
        let i_na = G_NA * m.powi(3) * h * (v - E_NA);
        let i_k = G_K * n.powi(4) * (v - E_K);
        let i_l = G_L * (v - E_L);
        let dv = (i_stim - (i_na + i_k + i_l)) / C_M;
        let dm = alpha_m(v) * (1.0 - m) - beta_m(v) * m;
        let dh = alpha_h(v) * (1.0 - h) - beta_h(v) * h;
        let dn = alpha_n(v) * (1.0 - n) - beta_n(v) * n;
        (dv, dm, dh, dn)
    }

    /// Advance one RK4 step of `dt` (ms) under injected current `i_stim`.
    fn step_rk4(&mut self, i_stim: f64, dt: f64) {
        let (v, m, h, n) = (self.v, self.m, self.h, self.n);
        let k1 = Self::derivs(v, m, h, n, i_stim);
        let k2 = Self::derivs(
            v + 0.5 * dt * k1.0,
            m + 0.5 * dt * k1.1,
            h + 0.5 * dt * k1.2,
            n + 0.5 * dt * k1.3,
            i_stim,
        );
        let k3 = Self::derivs(
            v + 0.5 * dt * k2.0,
            m + 0.5 * dt * k2.1,
            h + 0.5 * dt * k2.2,
            n + 0.5 * dt * k2.3,
            i_stim,
        );
        let k4 = Self::derivs(
            v + dt * k3.0,
            m + dt * k3.1,
            h + dt * k3.2,
            n + dt * k3.3,
            i_stim,
        );
        self.v += dt / 6.0 * (k1.0 + 2.0 * k2.0 + 2.0 * k3.0 + k4.0);
        self.m += dt / 6.0 * (k1.1 + 2.0 * k2.1 + 2.0 * k3.1 + k4.1);
        self.h += dt / 6.0 * (k1.2 + 2.0 * k2.2 + 2.0 * k3.2 + k4.2);
        self.n += dt / 6.0 * (k1.3 + 2.0 * k2.3 + 2.0 * k3.3 + k4.3);
    }

    /// Integrate for `duration_ms` at timestep `dt_ms` under a single `stim`
    /// pulse, returning the membrane-potential trace (mV) — one sample per
    /// step, including the initial value.
    pub fn run(&mut self, stim: StimPulse, duration_ms: f64, dt_ms: f64) -> Vec<f64> {
        self.run_two(stim, StimPulse::off(), duration_ms, dt_ms)
    }

    /// Like [`run`](Self::run) but with two stimulus pulses summed — used to
    /// probe the refractory period.
    pub fn run_two(
        &mut self,
        a: StimPulse,
        b: StimPulse,
        duration_ms: f64,
        dt_ms: f64,
    ) -> Vec<f64> {
        let n_steps = (duration_ms / dt_ms).round() as usize;
        let mut trace = Vec::with_capacity(n_steps + 1);
        trace.push(self.v);
        let mut t = 0.0;
        for _ in 0..n_steps {
            let i_stim = a.current_at(t) + b.current_at(t);
            self.step_rk4(i_stim, dt_ms);
            t += dt_ms;
            trace.push(self.v);
        }
        trace
    }
}

/// Count action potentials in a trace: the number of upward crossings of
/// `threshold` (mV).
pub fn count_spikes(trace: &[f64], threshold: f64) -> usize {
    let mut count = 0;
    for w in trace.windows(2) {
        if w[0] < threshold && w[1] >= threshold {
            count += 1;
        }
    }
    count
}

#[cfg(test)]
mod tests {
    use super::*;

    fn vmax(trace: &[f64]) -> f64 {
        trace.iter().cloned().fold(f64::MIN, f64::max)
    }

    #[test]
    fn suprathreshold_stimulus_fires_action_potential() {
        // A 50 µA/cm² × 0.5 ms pulse charges the membrane ~25 mV past rest —
        // well over the HH firing threshold — so a full AP fires.
        let mut c = HhCompartment::at_rest();
        let trace = c.run(
            StimPulse { amp_ua_cm2: 50.0, start_ms: 1.0, width_ms: 0.5 },
            20.0,
            0.005,
        );
        let peak = vmax(&trace);
        assert!(peak > 20.0, "AP should overshoot 0 mV; peak={peak}");
        assert!(peak < 60.0, "but not blow up; peak={peak}");
    }

    #[test]
    fn subthreshold_stimulus_does_not_fire() {
        // 5 µA/cm² × 0.5 ms ≈ 2.5 mV — far below threshold, no AP.
        let mut c = HhCompartment::at_rest();
        let trace = c.run(
            StimPulse { amp_ua_cm2: 5.0, start_ms: 1.0, width_ms: 0.5 },
            20.0,
            0.005,
        );
        let peak = vmax(&trace);
        assert!(peak < -40.0, "weak stim must not trigger an AP; peak={peak}");
    }

    #[test]
    fn resting_state_is_stable() {
        let mut c = HhCompartment::at_rest();
        let trace = c.run(StimPulse::off(), 10.0, 0.005);
        let drift = trace.last().unwrap() - V_REST;
        assert!(drift.abs() < 1.0, "rest must not drift; drift={drift}");
    }

    #[test]
    fn refractory_period_blocks_a_too_soon_second_spike() {
        // Two identical strong pulses. Close together (2 ms apart) the second
        // lands in the refractory period and is blocked → 1 spike. Far apart
        // (17 ms) the membrane has recovered and both fire → 2 spikes.
        let pulse = |start| StimPulse { amp_ua_cm2: 50.0, start_ms: start, width_ms: 0.5 };
        let mut c_close = HhCompartment::at_rest();
        let close = c_close.run_two(pulse(1.0), pulse(3.0), 30.0, 0.005);
        let mut c_far = HhCompartment::at_rest();
        let far = c_far.run_two(pulse(1.0), pulse(18.0), 30.0, 0.005);
        assert_eq!(count_spikes(&close, 0.0), 1, "second pulse blocked by refractory period");
        assert_eq!(count_spikes(&far, 0.0), 2, "well-separated pulses both fire");
    }
}
