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

pub(crate) fn alpha_m(v: f64) -> f64 {
    let x = v + 40.0;
    if x.abs() < 1e-6 {
        1.0 // limit of 0.1·x / (1 − e^{−x/10}) as x → 0
    } else {
        0.1 * x / (1.0 - (-x / 10.0).exp())
    }
}
pub(crate) fn beta_m(v: f64) -> f64 {
    4.0 * (-(v + 65.0) / 18.0).exp()
}
pub(crate) fn alpha_h(v: f64) -> f64 {
    0.07 * (-(v + 65.0) / 20.0).exp()
}
pub(crate) fn beta_h(v: f64) -> f64 {
    1.0 / (1.0 + (-(v + 35.0) / 10.0).exp())
}
pub(crate) fn alpha_n(v: f64) -> f64 {
    let x = v + 55.0;
    if x.abs() < 1e-6 {
        0.1 // limit of 0.01·x / (1 − e^{−x/10}) as x → 0
    } else {
        0.01 * x / (1.0 - (-x / 10.0).exp())
    }
}
pub(crate) fn beta_n(v: f64) -> f64 {
    0.125 * (-(v + 65.0) / 80.0).exp()
}

fn steady(alpha: f64, beta: f64) -> f64 {
    alpha / (alpha + beta)
}

/// The voltage-dependent steady states and time constants of the three
/// Hodgkin–Huxley gates at one membrane potential — the curves that define the
/// channel kinetics.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GatingKinetics {
    /// Na⁺ activation steady state `m∞(V)` (dimensionless, in `[0, 1]`).
    pub m_inf: f64,
    /// Na⁺ inactivation steady state `h∞(V)`.
    pub h_inf: f64,
    /// K⁺ activation steady state `n∞(V)`.
    pub n_inf: f64,
    /// Na⁺ activation time constant `τ_m(V)` (ms).
    pub tau_m_ms: f64,
    /// Na⁺ inactivation time constant `τ_h(V)` (ms).
    pub tau_h_ms: f64,
    /// K⁺ activation time constant `τ_n(V)` (ms).
    pub tau_n_ms: f64,
}

/// Evaluate the Hodgkin–Huxley gating kinetics at membrane potential `v_mv` (mV).
///
/// Each gate `x ∈ {m, h, n}` relaxes toward its steady state
/// `x∞ = α_x/(α_x+β_x)` with a time constant `τ_x = 1/(α_x+β_x)` (ms), from the
/// classic squid-axon `α`/`β` rate functions. These voltage-dependent curves are
/// the heart of the action potential: as the membrane **depolarises**, Na⁺
/// activation `m∞` and K⁺ activation `n∞` rise toward 1 while Na⁺ inactivation
/// `h∞` falls toward 0, and the fast-Na⁺ / slow-K⁺ separation `τ_m ≪ τ_n` is what
/// shapes the spike. Surfaces the gating curves that the integrator
/// ([`HhCompartment`]) otherwise only evaluates internally.
pub fn hh_gating_kinetics(v_mv: f64) -> GatingKinetics {
    let (am, bm) = (alpha_m(v_mv), beta_m(v_mv));
    let (ah, bh) = (alpha_h(v_mv), beta_h(v_mv));
    let (an, bn) = (alpha_n(v_mv), beta_n(v_mv));
    GatingKinetics {
        m_inf: am / (am + bm),
        h_inf: ah / (ah + bh),
        n_inf: an / (an + bn),
        tau_m_ms: 1.0 / (am + bm),
        tau_h_ms: 1.0 / (ah + bh),
        tau_n_ms: 1.0 / (an + bn),
    }
}

/// The **Boltzmann steady-state (in)activation** `x∞(V) = 1/(1 + e^(−(V−V½)/k))`
/// of a voltage-gated channel — the universal two-parameter sigmoid fit to a
/// measured voltage-clamp activation (or inactivation) curve, at membrane
/// potential `voltage_mv` `V` (mV), half-(in)activation voltage
/// `half_activation_mv` `V½` (mV) and slope factor `slope_mv` `k` (mV).
///
/// Where [`hh_gating_kinetics`] gives the *specific* Hodgkin–Huxley gates from
/// the squid-axon `α`/`β` rate equations, this is the *general* empirical form
/// every channel's steady state is reported as (Kᵥ, Naᵥ, Caᵥ, HCN, …). The slope
/// `k` sets both the steepness (smaller `|k|` ⇒ sharper switch) and, through its
/// **sign**, the direction: `k > 0` is **activation** — the gate opens with
/// depolarisation (rises `0 → 1` as `V` increases) — while `k < 0` is
/// **inactivation**, which closes with depolarisation. The curve is exactly
/// `0.5` at `V = V½`, point-symmetric about it (`x∞(V½+Δ) + x∞(V½−Δ) = 1`), and
/// lies strictly in `(0, 1)`. Returns `0.5` (the midpoint) for the degenerate
/// step-function limit `k = 0` and for any non-finite input.
pub fn boltzmann_activation(voltage_mv: f64, half_activation_mv: f64, slope_mv: f64) -> f64 {
    if !voltage_mv.is_finite()
        || !half_activation_mv.is_finite()
        || !slope_mv.is_finite()
        || slope_mv == 0.0
    {
        return 0.5;
    }
    1.0 / (1.0 + (-(voltage_mv - half_activation_mv) / slope_mv).exp())
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

// --- multi-compartment cable -------------------------------------------

/// Per-compartment state derivative `(dV, dm, dh, dn)`.
type Deriv = (f64, f64, f64, f64);

/// A uniform unmyelinated Hodgkin–Huxley **cable**: a chain of compartments
/// coupled by intracellular axial current, with sealed (no-flux) ends. An
/// action potential initiated at one end propagates along it.
pub struct HhCable {
    comps: Vec<HhCompartment>,
    /// Axial coupling coefficient `a / (2 Rᵢ Δx²)` in 1/(Ω·cm²).
    g_c: f64,
}

/// The outcome of a cable run: the peak membrane potential and the time it
/// occurred, per compartment.
pub struct CableRun {
    peak_v: Vec<f64>,
    peak_t_ms: Vec<f64>,
}

impl CableRun {
    /// The time (ms) of peak depolarization at compartment `idx`, or `None`
    /// if that compartment never fired (its peak stayed below 0 mV).
    pub fn peak_time_ms(&self, idx: usize) -> Option<f64> {
        if self.peak_v[idx] > 0.0 {
            Some(self.peak_t_ms[idx])
        } else {
            None
        }
    }

    /// The highest peak membrane potential (mV) reached by any compartment.
    pub fn max_peak_mv(&self) -> f64 {
        self.peak_v.iter().cloned().fold(f64::MIN, f64::max)
    }
}

fn cable_derivs(comps: &[HhCompartment], g_c: f64, ext_drive: &[f64]) -> Vec<Deriv> {
    let n = comps.len();
    (0..n)
        .map(|k| {
            let c = comps[k];
            let vm1 = if k > 0 { comps[k - 1].v } else { c.v };
            let vp1 = if k + 1 < n { comps[k + 1].v } else { c.v };
            // Axial current density (µA/cm²). g_c·ΔV is in mA/cm²
            // (mV / (Ω·cm²)); the 1e3 converts it to µA/cm².
            let i_axial = 1.0e3 * g_c * (vm1 - 2.0 * c.v + vp1);
            // `ext_drive[k]` is the per-compartment external current density
            // (µA/cm²) — either a block stimulus or the extracellular
            // activating term from a field.
            HhCompartment::derivs(c.v, c.m, c.h, c.n, ext_drive[k] + i_axial)
        })
        .collect()
}

fn add_scaled(s0: &[HhCompartment], k: &[Deriv], f: f64) -> Vec<HhCompartment> {
    s0.iter()
        .zip(k)
        .map(|(c, d)| HhCompartment {
            v: c.v + f * d.0,
            m: c.m + f * d.1,
            h: c.h + f * d.2,
            n: c.n + f * d.3,
        })
        .collect()
}

impl HhCable {
    /// Build a uniform cable of `n` compartments, each `dx_um` long, of fiber
    /// radius `a_um`, with intracellular resistivity `ri_ohm_cm` (Ω·cm). All
    /// compartments start at rest.
    pub fn uniform(n: usize, dx_um: f64, a_um: f64, ri_ohm_cm: f64) -> Self {
        let dx_cm = dx_um * 1.0e-4;
        let a_cm = a_um * 1.0e-4;
        let g_c = a_cm / (2.0 * ri_ohm_cm * dx_cm * dx_cm);
        Self {
            comps: vec![HhCompartment::at_rest(); n],
            g_c,
        }
    }

    fn step_rk4(&mut self, ext_drive: &[f64], dt: f64) {
        let g_c = self.g_c;
        let k1 = cable_derivs(&self.comps, g_c, ext_drive);
        let s1 = add_scaled(&self.comps, &k1, 0.5 * dt);
        let k2 = cable_derivs(&s1, g_c, ext_drive);
        let s2 = add_scaled(&self.comps, &k2, 0.5 * dt);
        let k3 = cable_derivs(&s2, g_c, ext_drive);
        let s3 = add_scaled(&self.comps, &k3, dt);
        let k4 = cable_derivs(&s3, g_c, ext_drive);
        for k in 0..self.comps.len() {
            let c = &mut self.comps[k];
            c.v += dt / 6.0 * (k1[k].0 + 2.0 * k2[k].0 + 2.0 * k3[k].0 + k4[k].0);
            c.m += dt / 6.0 * (k1[k].1 + 2.0 * k2[k].1 + 2.0 * k3[k].1 + k4[k].1);
            c.h += dt / 6.0 * (k1[k].2 + 2.0 * k2[k].2 + 2.0 * k3[k].2 + k4[k].2);
            c.n += dt / 6.0 * (k1[k].3 + 2.0 * k2[k].3 + 2.0 * k3[k].3 + k4[k].3);
        }
    }

    /// Stimulate compartment 0 with `stim` and integrate for `duration_ms` at
    /// timestep `dt_ms`, recording each compartment's peak depolarization.
    pub fn stimulate_end(&mut self, stim: StimPulse, duration_ms: f64, dt_ms: f64) -> CableRun {
        let n = self.comps.len();
        // Stimulate a small block at the near end (≈10 % of the cable).
        let n_stim = (n / 10).max(1);
        let mut peak_v: Vec<f64> = self.comps.iter().map(|c| c.v).collect();
        let mut peak_t_ms = vec![0.0; n];
        let n_steps = (duration_ms / dt_ms).round() as usize;
        let mut drive = vec![0.0; n];
        let mut t = 0.0;
        for _ in 0..n_steps {
            let i = stim.current_at(t);
            for (k, d) in drive.iter_mut().enumerate() {
                *d = if k < n_stim { i } else { 0.0 };
            }
            self.step_rk4(&drive, dt_ms);
            t += dt_ms;
            for k in 0..n {
                if self.comps[k].v > peak_v[k] {
                    peak_v[k] = self.comps[k].v;
                    peak_t_ms[k] = t;
                }
            }
        }
        CableRun { peak_v, peak_t_ms }
    }

    /// Integrate with a per-compartment extracellular activating `drive`
    /// (µA/cm²) gated on during `[start_ms, start_ms + width_ms)`.
    pub fn stimulate_extracellular(
        &mut self,
        drive: &[f64],
        start_ms: f64,
        width_ms: f64,
        duration_ms: f64,
        dt_ms: f64,
    ) -> CableRun {
        let n = self.comps.len();
        let mut peak_v: Vec<f64> = self.comps.iter().map(|c| c.v).collect();
        let mut peak_t_ms = vec![0.0; n];
        let n_steps = (duration_ms / dt_ms).round() as usize;
        let mut buf = vec![0.0; n];
        let mut t = 0.0;
        for _ in 0..n_steps {
            let on = t >= start_ms && t < start_ms + width_ms;
            for (k, b) in buf.iter_mut().enumerate() {
                *b = if on { drive[k] } else { 0.0 };
            }
            self.step_rk4(&buf, dt_ms);
            t += dt_ms;
            for k in 0..n {
                if self.comps[k].v > peak_v[k] {
                    peak_v[k] = self.comps[k].v;
                    peak_t_ms[k] = t;
                }
            }
        }
        CableRun { peak_v, peak_t_ms }
    }
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
    fn hh_gating_kinetics_match_the_classic_resting_and_limit_values() {
        // At the resting potential the classic HH gates sit near m≈0.05, h≈0.6, n≈0.32.
        let rest = hh_gating_kinetics(-65.0);
        assert!((rest.m_inf - 0.053).abs() < 0.005, "m∞(−65) {}", rest.m_inf);
        assert!((rest.h_inf - 0.596).abs() < 0.01, "h∞(−65) {}", rest.h_inf);
        assert!((0.30..0.34).contains(&rest.n_inf), "n∞(−65) ≈ 0.32, got {}", rest.n_inf);
        // Na activation is fast, K activation slow — the separation that makes a spike.
        assert!(rest.tau_m_ms < rest.tau_n_ms, "τ_m {} ≪ τ_n {}", rest.tau_m_ms, rest.tau_n_ms);
        assert!(rest.tau_m_ms > 0.0 && rest.tau_h_ms > 0.0 && rest.tau_n_ms > 0.0);

        // Strong hyperpolarisation: activation shut (m∞,n∞→0), Na de-inactivated (h∞→1).
        let lo = hh_gating_kinetics(-100.0);
        assert!(lo.m_inf < 0.02 && lo.n_inf < 0.1, "hyperpol activation {} {}", lo.m_inf, lo.n_inf);
        assert!(lo.h_inf > 0.95, "hyperpol h∞ {}", lo.h_inf);

        // Strong depolarisation: activation on (m∞,n∞→1), Na inactivated (h∞→0).
        let hi = hh_gating_kinetics(50.0);
        assert!(hi.m_inf > 0.98 && hi.n_inf > 0.9, "depol activation {} {}", hi.m_inf, hi.n_inf);
        assert!(hi.h_inf < 0.02, "depol h∞ {}", hi.h_inf);

        // Monotonic: m∞ and n∞ rise with depolarisation, h∞ falls; all stay in [0,1].
        let mid = hh_gating_kinetics(-30.0);
        assert!(lo.m_inf < mid.m_inf && mid.m_inf < hi.m_inf, "m∞ rises with V");
        assert!(lo.n_inf < mid.n_inf && mid.n_inf < hi.n_inf, "n∞ rises with V");
        assert!(lo.h_inf > mid.h_inf && mid.h_inf > hi.h_inf, "h∞ falls with V");
        for x in [rest.m_inf, rest.h_inf, rest.n_inf, lo.m_inf, hi.h_inf, mid.n_inf] {
            assert!((0.0..=1.0).contains(&x), "gate {x} must be in [0,1]");
        }
    }

    #[test]
    fn boltzmann_activation_is_the_universal_voltage_clamp_sigmoid() {
        let (v_half, k) = (-30.0, 9.0); // a typical Kᵥ activation: V½ = −30 mV, k = 9 mV
        // Half-activated exactly at V = V½, for any slope (and either sign).
        assert!((boltzmann_activation(v_half, v_half, k) - 0.5).abs() < 1e-15);
        assert!((boltzmann_activation(v_half, v_half, -4.0) - 0.5).abs() < 1e-15);
        // One slope-factor either side of V½ → the textbook 0.7311 / 0.2689.
        let up = boltzmann_activation(v_half + k, v_half, k);
        assert!((up - 1.0 / (1.0 + (-1.0_f64).exp())).abs() < 1e-12, "closed form");
        assert!((up - 0.7311).abs() < 1e-3, "+k → 0.7311, got {up}");
        assert!((boltzmann_activation(v_half - k, v_half, k) - 0.2689).abs() < 1e-3, "−k → 0.2689");
        // Activation (k>0): opens (rises 0→1) with depolarisation, saturating both ways.
        assert!(boltzmann_activation(30.0, v_half, k) > 0.99, "far-depol → open");
        assert!(boltzmann_activation(-90.0, v_half, k) < 0.01, "far-hyperpol → shut");
        assert!(
            boltzmann_activation(0.0, v_half, k) > boltzmann_activation(-60.0, v_half, k),
            "activation is monotone up in V"
        );
        // Inactivation (k<0): the mirror — closes with depolarisation.
        assert!(boltzmann_activation(30.0, v_half, -k) < 0.01, "inactivation shuts when depolarised");
        assert!(
            boltzmann_activation(0.0, v_half, -k) < boltzmann_activation(-60.0, v_half, -k),
            "inactivation is monotone down in V"
        );
        // Point-symmetric about V½: x∞(V½+Δ) + x∞(V½−Δ) = 1 exactly.
        for d in [1.0_f64, 7.5, 20.0, 50.0] {
            let s = boltzmann_activation(v_half + d, v_half, k)
                + boltzmann_activation(v_half - d, v_half, k);
            assert!((s - 1.0).abs() < 1e-12, "symmetric at Δ={d}: {s}");
        }
        // Always a valid open fraction in (0, 1).
        for v in [-120.0_f64, -50.0, -30.0, 0.0, 50.0] {
            let x = boltzmann_activation(v, v_half, k);
            assert!(x > 0.0 && x < 1.0, "x∞ in (0,1) at V={v}: {x}");
        }
        // The degenerate step (k=0) and any non-finite input → the 0.5 midpoint.
        assert!((boltzmann_activation(0.0, v_half, 0.0) - 0.5).abs() < 1e-15, "k=0 → 0.5");
        assert!((boltzmann_activation(f64::NAN, v_half, k) - 0.5).abs() < 1e-15, "NaN → 0.5");
        assert!((boltzmann_activation(0.0, v_half, f64::INFINITY) - 0.5).abs() < 1e-15, "∞ → 0.5");
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

    #[test]
    fn action_potential_propagates_with_finite_velocity() {
        // squid giant axon: radius 238 µm, Rᵢ = 35.4 Ω·cm, 100-µm compartments
        let mut cable = HhCable::uniform(200, 100.0, 238.0, 35.4);
        let r = cable.stimulate_end(
            StimPulse { amp_ua_cm2: 100.0, start_ms: 1.0, width_ms: 0.5 },
            15.0,
            0.01,
        );
        let t50 = r.peak_time_ms(50).expect("comp 50 fires");
        let t150 = r.peak_time_ms(150).expect("comp 150 fires");
        assert!(t150 > t50, "AP travels 50→150 in +x: t50={t50} t150={t150}");
        let dist_mm = (150 - 50) as f64 * 0.1; // Δx = 100 µm = 0.1 mm → 10 mm
        let vel_m_s = dist_mm / (t150 - t50); // mm/ms ≡ m/s
        assert!(
            (1.0..100.0).contains(&vel_m_s),
            "squid-like conduction velocity expected; got {vel_m_s} m/s"
        );
    }
}
