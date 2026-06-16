//! The Morris–Lecar conductance-based spiking-neuron model.
//!
//! A two-variable model derived from voltage-clamp data on barnacle muscle
//! fibre — unlike the polynomial Izhikevich or the abstract cubic
//! FitzHugh–Nagumo, its currents are *biophysical*: an instantaneous Ca²⁺
//! current, a slower K⁺ current, and a leak, with real reversal potentials and
//! conductances. It is the textbook vehicle for Type-I (continuous f–I, SNIC)
//! excitability.
//!
//! ```text
//! C dV/dt = I − g_L(V − V_L) − g_Ca·m∞(V)·(V − V_Ca) − g_K·w·(V − V_K)
//! dw/dt   = φ · (w∞(V) − w) / τ_w(V)
//!
//! m∞(V) = ½[1 + tanh((V − V₁)/V₂)]
//! w∞(V) = ½[1 + tanh((V − V₃)/V₄)]
//! τ_w(V) = 1 / cosh((V − V₃)/(2 V₄))
//! ```
//!
//! `V` is in mV, `I` in µA/cm², conductances in mS/cm², `C` in µF/cm², time in
//! ms. The default parameters are the standard **Type-I (SNIC)** set
//! (Ermentrout & Terman, *Mathematical Foundations of Neuroscience*).
//!
//! Reference: Morris & Lecar (1981), *Biophys. J.* 35. Integration is RK4.
//!
//! # Honest scope
//!
//! Research/educational. A reduced two-current caricature with the gating
//! curves fit to one preparation — excellent for excitability and bifurcation
//! intuition, not a quantitative model of any specific mammalian neuron (use
//! Hodgkin–Huxley, [`crate::cable`], for detailed biophysics).

/// Morris–Lecar parameters (µA/cm², mS/cm², mV, µF/cm², ms units).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MorrisLecarParams {
    /// Membrane capacitance `C` (µF/cm²).
    pub c: f64,
    /// Leak conductance `g_L`.
    pub g_l: f64,
    /// Calcium conductance `g_Ca`.
    pub g_ca: f64,
    /// Potassium conductance `g_K`.
    pub g_k: f64,
    /// Leak reversal `V_L` (mV).
    pub v_l: f64,
    /// Calcium reversal `V_Ca` (mV).
    pub v_ca: f64,
    /// Potassium reversal `V_K` (mV).
    pub v_k: f64,
    /// Ca-gate midpoint `V₁` (mV).
    pub v1: f64,
    /// Ca-gate slope `V₂` (mV).
    pub v2: f64,
    /// K-gate midpoint `V₃` (mV).
    pub v3: f64,
    /// K-gate slope `V₄` (mV).
    pub v4: f64,
    /// Recovery rate `φ` (1/ms).
    pub phi: f64,
}

impl Default for MorrisLecarParams {
    /// The standard Type-I (SNIC) parameter set.
    fn default() -> Self {
        Self {
            c: 20.0,
            g_l: 2.0,
            g_ca: 4.0,
            g_k: 8.0,
            v_l: -60.0,
            v_ca: 120.0,
            v_k: -84.0,
            v1: -1.2,
            v2: 18.0,
            v3: 12.0,
            v4: 17.4,
            phi: 0.066_7,
        }
    }
}

impl MorrisLecarParams {
    /// Steady-state Ca-gate activation `m∞(V)`.
    pub fn m_inf(&self, v: f64) -> f64 {
        0.5 * (1.0 + ((v - self.v1) / self.v2).tanh())
    }

    /// Steady-state K-gate activation `w∞(V)`.
    pub fn w_inf(&self, v: f64) -> f64 {
        0.5 * (1.0 + ((v - self.v3) / self.v4).tanh())
    }

    /// K-gate time-constant factor `τ_w(V)` (ms).
    pub fn tau_w(&self, v: f64) -> f64 {
        1.0 / ((v - self.v3) / (2.0 * self.v4)).cosh()
    }

    /// Net ionic + applied current at `(V, w)` with input `current`
    /// (the `C dV/dt` right-hand side, µA/cm²).
    fn ionic(&self, v: f64, w: f64, current: f64) -> f64 {
        current
            - self.g_l * (v - self.v_l)
            - self.g_ca * self.m_inf(v) * (v - self.v_ca)
            - self.g_k * w * (v - self.v_k)
    }

    /// The resting equilibrium `(V*, w*)` for a constant `current`: the root of
    /// the steady-state balance with `w = w∞(V)`, found by bisection on
    /// `[V_K, V_Ca]`. (For sub-threshold `current` this is the stable rest
    /// state; in the bistable range it returns one equilibrium on that bracket.)
    pub fn fixed_point(&self, current: f64) -> (f64, f64) {
        let f = |v: f64| self.ionic(v, self.w_inf(v), current);
        let mut lo = self.v_k;
        let mut hi = self.v_ca;
        let (mut flo, _fhi) = (f(lo), f(hi));
        for _ in 0..80 {
            let mid = 0.5 * (lo + hi);
            let fmid = f(mid);
            if (fmid > 0.0) == (flo > 0.0) {
                lo = mid;
                flo = fmid;
            } else {
                hi = mid;
            }
        }
        let v = 0.5 * (lo + hi);
        (v, self.w_inf(v))
    }
}

/// One Morris–Lecar cell's state.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MorrisLecar {
    /// Membrane potential (mV).
    pub v: f64,
    /// K-gate recovery variable.
    pub w: f64,
    /// Model parameters.
    pub params: MorrisLecarParams,
}

impl MorrisLecar {
    /// A cell initialised at its `I = 0` resting equilibrium.
    pub fn at_rest(params: MorrisLecarParams) -> Self {
        let (v, w) = params.fixed_point(0.0);
        Self { v, w, params }
    }

    fn deriv(&self, v: f64, w: f64, current: f64) -> (f64, f64) {
        let dv = self.params.ionic(v, w, current) / self.params.c;
        let dw = self.params.phi * (self.params.w_inf(v) - w) / self.params.tau_w(v);
        (dv, dw)
    }

    /// Advance by `dt` ms with constant input `current` (µA/cm²) using RK4.
    pub fn step(&mut self, dt: f64, current: f64) {
        let (k1v, k1w) = self.deriv(self.v, self.w, current);
        let (k2v, k2w) = self.deriv(self.v + 0.5 * dt * k1v, self.w + 0.5 * dt * k1w, current);
        let (k3v, k3w) = self.deriv(self.v + 0.5 * dt * k2v, self.w + 0.5 * dt * k2w, current);
        let (k4v, k4w) = self.deriv(self.v + dt * k3v, self.w + dt * k3w, current);
        self.v += dt / 6.0 * (k1v + 2.0 * k2v + 2.0 * k3v + k4v);
        self.w += dt / 6.0 * (k1w + 2.0 * k2w + 2.0 * k3w + k4w);
    }
}

/// Count spikes (upward crossings of `threshold` mV by `V`) over a simulation
/// from rest under constant `current`. Module-scoped to avoid colliding with the
/// Hodgkin–Huxley `count_spikes`.
pub fn count_spikes(
    params: MorrisLecarParams,
    current: f64,
    dt: f64,
    duration: f64,
    threshold: f64,
) -> usize {
    let mut cell = MorrisLecar::at_rest(params);
    let steps = (duration / dt).round() as usize;
    let mut spikes = 0;
    let mut prev = cell.v;
    for _ in 0..steps {
        cell.step(dt, current);
        if prev < threshold && cell.v >= threshold {
            spikes += 1;
        }
        prev = cell.v;
    }
    spikes
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn steady_state_gates_are_half_at_their_midpoints() {
        let p = MorrisLecarParams::default();
        assert!((p.m_inf(p.v1) - 0.5).abs() < 1e-12);
        assert!((p.w_inf(p.v3) - 0.5).abs() < 1e-12);
        // Saturating limits.
        assert!(p.m_inf(1000.0) > 0.999 && p.m_inf(-1000.0) < 1e-3);
    }

    #[test]
    fn resting_potential_is_a_true_fixed_point_near_minus_60() {
        let p = MorrisLecarParams::default();
        let (v, w) = p.fixed_point(0.0);
        assert!((-62.0..=-58.0).contains(&v), "rest V = {v}");
        // It is a genuine equilibrium: both derivatives vanish.
        assert!(p.ionic(v, w, 0.0).abs() < 1e-6);
        assert!((p.w_inf(v) - w).abs() < 1e-12);
        // And a simulation started there does not drift or spike.
        let mut cell = MorrisLecar::at_rest(p);
        for _ in 0..20_000 {
            cell.step(0.05, 0.0);
        }
        assert!((cell.v - v).abs() < 1e-3, "drifted to {}", cell.v);
    }

    #[test]
    fn quiescent_at_rest_and_fires_under_strong_drive() {
        let p = MorrisLecarParams::default();
        // No input: stable rest, no spikes.
        assert_eq!(count_spikes(p, 0.0, 0.05, 1000.0, 0.0), 0);
        // Well above the SNIC threshold: repetitive spiking.
        let firing = count_spikes(p, 100.0, 0.05, 1000.0, 0.0);
        assert!(firing >= 3, "I=100 should fire repeatedly, got {firing}");
    }

    #[test]
    fn fires_in_a_window_then_depolarisation_blocks() {
        let p = MorrisLecarParams::default();
        // Mid-range drive gives repetitive firing...
        let mid = count_spikes(p, 60.0, 0.05, 1000.0, 0.0);
        assert!(mid >= 3, "I=60 should fire repeatedly, got {mid}");
        // ...but very strong drive locks the membrane depolarised
        // (depolarisation block): the K current can no longer repolarise, so
        // firing ceases — a real Morris–Lecar property, not monotone f-I.
        let blocked = count_spikes(p, 150.0, 0.05, 1000.0, 0.0);
        assert!(
            blocked <= 2,
            "I=150 should depolarisation-block, got {blocked}"
        );
    }
}
