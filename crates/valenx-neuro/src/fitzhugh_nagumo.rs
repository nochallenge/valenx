//! The FitzHugh–Nagumo excitable-membrane model.
//!
//! A two-variable reduction of Hodgkin–Huxley that keeps the essential
//! qualitative dynamics — a resting state, an excitation threshold, and
//! oscillation under sustained drive — in a form simple enough to reason about
//! with phase-plane geometry:
//!
//! ```text
//! dv/dt = v − v³/3 − w + I        (fast "voltage" variable)
//! dw/dt = ε (v + a − b w)          (slow recovery variable)
//! ```
//!
//! The cubic `v`-nullcline `w = v − v³/3 + I` and the straight `w`-nullcline
//! `w = (v + a)/b` cross at the single equilibrium ([`FitzHughNagumoParams::fixed_point`]).
//! For the classic parameters (`a = 0.7, b = 0.8, ε = 0.08`) that equilibrium is
//! a stable rest state at low `I`; as `I` increases past a Hopf bifurcation
//! (`I ≈ 0.33`) it loses stability and the system settles onto a limit cycle —
//! repetitive spiking.
//!
//! Reference: FitzHugh (1961), *Biophys. J.* 1; Nagumo et al. (1962). Integration
//! is fixed-step RK4.
//!
//! # Honest scope
//!
//! Research/educational. FitzHugh–Nagumo is a dimensionless *qualitative* model:
//! `v` and `w` are not millivolts and microsiemens, and the cubic is a caricature
//! of the fast sodium current. It is the right tool for excitability and
//! phase-plane intuition, not for quantitative membrane biophysics (use
//! Hodgkin–Huxley, [`crate::cable`], for that).

/// FitzHugh–Nagumo parameters.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FitzHughNagumoParams {
    /// Recovery offset `a`.
    pub a: f64,
    /// Recovery–voltage coupling `b`.
    pub b: f64,
    /// Time-scale separation `ε` (small = slow recovery).
    pub epsilon: f64,
}

impl Default for FitzHughNagumoParams {
    /// The classic parameters `a = 0.7, b = 0.8, ε = 0.08`.
    fn default() -> Self {
        Self {
            a: 0.7,
            b: 0.8,
            epsilon: 0.08,
        }
    }
}

impl FitzHughNagumoParams {
    /// The equilibrium `(v*, w*)` for a constant input `current` — the unique
    /// intersection of the cubic `v`-nullcline and the linear `w`-nullcline,
    /// found by Newton's method on
    /// `g(v) = I + v − v³/3 − (v + a)/b` (monotone for the classic `b < 1`,
    /// so the iteration converges to the single real root).
    pub fn fixed_point(&self, current: f64) -> (f64, f64) {
        let mut v = 0.0_f64;
        for _ in 0..100 {
            let g = current + v - v * v * v / 3.0 - (v + self.a) / self.b;
            let g_prime = 1.0 - v * v - 1.0 / self.b;
            v -= g / g_prime;
        }
        (v, (v + self.a) / self.b)
    }
}

/// One FitzHugh–Nagumo cell's state.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FitzHughNagumo {
    /// Fast "voltage" variable.
    pub v: f64,
    /// Slow recovery variable.
    pub w: f64,
    /// Model parameters.
    pub params: FitzHughNagumoParams,
}

impl FitzHughNagumo {
    /// A cell initialised at its `I = 0` resting equilibrium.
    pub fn at_rest(params: FitzHughNagumoParams) -> Self {
        let (v, w) = params.fixed_point(0.0);
        Self { v, w, params }
    }

    fn deriv(&self, v: f64, w: f64, current: f64) -> (f64, f64) {
        let dv = v - v * v * v / 3.0 - w + current;
        let dw = self.params.epsilon * (v + self.params.a - self.params.b * w);
        (dv, dw)
    }

    /// Advance by `dt` with constant input `current` using RK4.
    pub fn step(&mut self, dt: f64, current: f64) {
        let (k1v, k1w) = self.deriv(self.v, self.w, current);
        let (k2v, k2w) = self.deriv(self.v + 0.5 * dt * k1v, self.w + 0.5 * dt * k1w, current);
        let (k3v, k3w) = self.deriv(self.v + 0.5 * dt * k2v, self.w + 0.5 * dt * k2w, current);
        let (k4v, k4w) = self.deriv(self.v + dt * k3v, self.w + dt * k3w, current);
        self.v += dt / 6.0 * (k1v + 2.0 * k2v + 2.0 * k3v + k4v);
        self.w += dt / 6.0 * (k1w + 2.0 * k2w + 2.0 * k3w + k4w);
    }
}

/// Count upward threshold crossings of `v` (spikes) over a simulation from rest
/// under constant `current`. A crossing is counted when `v` rises across
/// `threshold` between two steps.
pub fn count_spikes(
    params: FitzHughNagumoParams,
    current: f64,
    dt: f64,
    duration: f64,
    threshold: f64,
) -> usize {
    let mut cell = FitzHughNagumo::at_rest(params);
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
    fn fixed_point_for_classic_params() {
        let p = FitzHughNagumoParams::default();
        let (v, w) = p.fixed_point(0.0);
        // Analytic resting equilibrium ≈ (−1.1994, −0.6243).
        assert!((v - (-1.1994)).abs() < 1e-3, "v* = {v}");
        assert!((w - (-0.6243)).abs() < 1e-3, "w* = {w}");
        // It really is a root of the nullcline-balance g(v).
        let g = v - v * v * v / 3.0 - (v + p.a) / p.b;
        assert!(g.abs() < 1e-9, "g(v*) = {g}");
    }

    #[test]
    fn rest_is_stable_with_no_drive() {
        let p = FitzHughNagumoParams::default();
        let mut cell = FitzHughNagumo::at_rest(p);
        let (v0, w0) = (cell.v, cell.w);
        for _ in 0..10_000 {
            cell.step(0.05, 0.0);
        }
        assert!((cell.v - v0).abs() < 1e-6, "v drifted to {}", cell.v);
        assert!((cell.w - w0).abs() < 1e-6, "w drifted to {}", cell.w);
    }

    #[test]
    fn sustained_current_drives_a_limit_cycle() {
        let p = FitzHughNagumoParams::default();
        // Above the Hopf threshold (~0.33): repetitive spiking.
        let firing = count_spikes(p, 0.5, 0.01, 400.0, 1.0);
        assert!(firing >= 3, "I=0.5 should oscillate, got {firing} spikes");
        // No drive: the rest state is stable, so no spikes.
        let quiescent = count_spikes(p, 0.0, 0.01, 400.0, 1.0);
        assert_eq!(quiescent, 0, "I=0 should be quiescent, got {quiescent}");
    }

    #[test]
    fn excitable_threshold_separates_small_and_large_responses() {
        let p = FitzHughNagumoParams::default();
        // A small displacement from rest decays back; a large one fires a spike
        // (a big v-excursion) first — the signature of an excitable system.
        let relax = |kick: f64| {
            let mut cell = FitzHughNagumo::at_rest(p);
            cell.v += kick;
            let mut max_v = cell.v;
            for _ in 0..4000 {
                cell.step(0.05, 0.0);
                max_v = max_v.max(cell.v);
            }
            max_v
        };
        let small = relax(0.2);
        let large = relax(1.5);
        assert!(small < 0.5, "sub-threshold kick spiked: max v {small}");
        assert!(
            large > 1.5,
            "supra-threshold kick failed to spike: max v {large}"
        );
    }
}
