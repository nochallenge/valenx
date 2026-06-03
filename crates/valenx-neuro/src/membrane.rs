//! Generic membrane models + an implicit cable integrator.
//!
//! Phase 2 generalises the cable so different membrane models — Hodgkin–Huxley,
//! and (next) myelinated mammalian — share one **numerically stable**
//! integrator. The axial diffusion is solved **implicitly** (a tridiagonal
//! solve each step, Thomas algorithm), making the integrator unconditionally
//! stable: unlike the v1 explicit RK4, it does not diverge at fine (100 µm)
//! compartments in the sub-threshold regime.

use crate::cable::{
    alpha_h, alpha_m, alpha_n, beta_h, beta_m, beta_n, C_M, E_K, E_L, E_NA, G_K, G_L, G_NA, V_REST,
};

/// Solve a tridiagonal linear system by the Thomas algorithm: `sub` is the
/// sub-diagonal (`sub[0]` unused), `diag` the diagonal, `sup` the
/// super-diagonal (`sup[n-1]` unused), `rhs` the right-hand side.
fn thomas(sub: &[f64], diag: &[f64], sup: &[f64], rhs: &[f64]) -> Vec<f64> {
    let n = diag.len();
    let mut c = vec![0.0; n];
    let mut d = vec![0.0; n];
    c[0] = sup[0] / diag[0];
    d[0] = rhs[0] / diag[0];
    for i in 1..n {
        let m = diag[i] - sub[i] * c[i - 1];
        c[i] = sup[i] / m;
        d[i] = (rhs[i] - sub[i] * d[i - 1]) / m;
    }
    let mut x = vec![0.0; n];
    x[n - 1] = d[n - 1];
    for i in (0..n - 1).rev() {
        x[i] = d[i] - c[i] * x[i + 1];
    }
    x
}

/// Exponential-Euler update of a gating variable toward its steady state —
/// stable for any `dt` (the gating ODEs are linear in the gate).
fn exp_euler(x: f64, a: f64, b: f64, dt: f64) -> f64 {
    let inf = a / (a + b);
    let tau = 1.0 / (a + b);
    inf + (x - inf) * (-dt / tau).exp()
}

/// A membrane model: supplies its capacitance + ionic current and advances its
/// own gating state, given the membrane potential. Units: potential mV, current
/// density µA/cm², time ms, capacitance µF/cm².
pub trait Membrane: Clone {
    /// Membrane capacitance (µF/cm²).
    fn c_m(&self) -> f64;
    /// Resting potential (mV).
    fn v_rest(&self) -> f64;
    /// Total ionic current density (µA/cm²) at potential `v`, current gates.
    fn ionic_current(&self, v: f64) -> f64;
    /// Advance the gating variables one `dt` (ms) step at potential `v`.
    fn advance_gates(&mut self, v: f64, dt: f64);
}

/// A Hodgkin–Huxley membrane patch (the v1 squid kinetics) as a [`Membrane`].
#[derive(Clone, Copy, Debug)]
pub struct HhMembrane {
    m: f64,
    h: f64,
    n: f64,
}

impl HhMembrane {
    /// At rest — gates at their steady-state values for the resting potential.
    pub fn at_rest() -> Self {
        let v = V_REST;
        Self {
            m: alpha_m(v) / (alpha_m(v) + beta_m(v)),
            h: alpha_h(v) / (alpha_h(v) + beta_h(v)),
            n: alpha_n(v) / (alpha_n(v) + beta_n(v)),
        }
    }
}

impl Membrane for HhMembrane {
    fn c_m(&self) -> f64 {
        C_M
    }
    fn v_rest(&self) -> f64 {
        V_REST
    }
    fn ionic_current(&self, v: f64) -> f64 {
        G_NA * self.m.powi(3) * self.h * (v - E_NA)
            + G_K * self.n.powi(4) * (v - E_K)
            + G_L * (v - E_L)
    }
    fn advance_gates(&mut self, v: f64, dt: f64) {
        self.m = exp_euler(self.m, alpha_m(v), beta_m(v), dt);
        self.h = exp_euler(self.h, alpha_h(v), beta_h(v), dt);
        self.n = exp_euler(self.n, alpha_n(v), beta_n(v), dt);
    }
}

/// A uniform cable of compartments sharing a [`Membrane`] model, integrated
/// with an **implicit** (backward-Euler) axial solve — unconditionally stable.
pub struct ImplicitCable<M: Membrane> {
    v: Vec<f64>,
    mem: Vec<M>,
    /// Axial coupling `g` (µA/cm² per mV) between adjacent compartments.
    g_axial: f64,
}

impl<M: Membrane> ImplicitCable<M> {
    /// Build a uniform cable of `n` compartments, each `dx_um` long, of fiber
    /// radius `a_um`, intracellular resistivity `ri_ohm_cm` (Ω·cm), all from a
    /// `proto` membrane at rest.
    pub fn uniform(n: usize, proto: M, dx_um: f64, a_um: f64, ri_ohm_cm: f64) -> Self {
        let dx_cm = dx_um * 1.0e-4;
        let a_cm = a_um * 1.0e-4;
        // a / (2 Rᵢ Δx²) in 1/(Ω·cm²); ×1e3 → µA/cm² per mV.
        let g_axial = 1.0e3 * a_cm / (2.0 * ri_ohm_cm * dx_cm * dx_cm);
        let vr = proto.v_rest();
        Self {
            v: vec![vr; n],
            mem: vec![proto; n],
            g_axial,
        }
    }

    /// One implicit step of `dt` (ms) under per-compartment external drive
    /// `ext` (µA/cm²): gates advance explicitly, then a tridiagonal solve
    /// updates the membrane potential (backward-Euler in the axial term).
    fn step(&mut self, ext: &[f64], dt: f64) {
        let n = self.v.len();
        for k in 0..n {
            self.mem[k].advance_gates(self.v[k], dt);
        }
        let g = self.g_axial;
        let mut sub = vec![0.0; n];
        let mut diag = vec![0.0; n];
        let mut sup = vec![0.0; n];
        let mut rhs = vec![0.0; n];
        for k in 0..n {
            let cm = self.mem[k].c_m();
            let left = if k > 0 { g } else { 0.0 };
            let right = if k + 1 < n { g } else { 0.0 };
            sub[k] = -left * dt;
            sup[k] = -right * dt;
            diag[k] = cm + (left + right) * dt;
            let i_ion = self.mem[k].ionic_current(self.v[k]);
            rhs[k] = cm * self.v[k] + dt * (ext[k] - i_ion);
        }
        self.v = thomas(&sub, &diag, &sup, &rhs);
    }

    /// Stimulate the first ≈10 % of compartments with `amp` (µA/cm²) during
    /// `[start, start+width)` ms, integrating `duration` ms at timestep `dt`.
    /// Returns the peak membrane potential reached per compartment (mV).
    pub fn stimulate_block(
        &mut self,
        amp: f64,
        start: f64,
        width: f64,
        duration: f64,
        dt: f64,
    ) -> Vec<f64> {
        let n = self.v.len();
        let n_stim = (n / 10).max(1);
        let mut peak = self.v.clone();
        let mut ext = vec![0.0; n];
        let steps = (duration / dt).round() as usize;
        let mut t = 0.0;
        for _ in 0..steps {
            let on = t >= start && t < start + width;
            for (k, e) in ext.iter_mut().enumerate() {
                *e = if on && k < n_stim { amp } else { 0.0 };
            }
            self.step(&ext, dt);
            t += dt;
            for k in 0..n {
                if self.v[k] > peak[k] {
                    peak[k] = self.v[k];
                }
            }
        }
        peak
    }

    /// Current membrane potential per compartment (mV).
    pub fn potentials(&self) -> &[f64] {
        &self.v
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn thomas_solves_a_known_system() {
        // Tridiagonal [[2,1,0],[1,2,1],[0,1,2]] · x = [1,2,3].
        let x = thomas(&[0.0, 1.0, 1.0], &[2.0, 2.0, 2.0], &[1.0, 1.0, 0.0], &[1.0, 2.0, 3.0]);
        assert!((2.0 * x[0] + x[1] - 1.0).abs() < 1e-12);
        assert!((x[0] + 2.0 * x[1] + x[2] - 2.0).abs() < 1e-12);
        assert!((x[1] + 2.0 * x[2] - 3.0).abs() < 1e-12);
    }

    #[test]
    fn implicit_single_compartment_fires_action_potential() {
        let mut cable = ImplicitCable::uniform(1, HhMembrane::at_rest(), 100.0, 238.0, 35.4);
        let peak = cable.stimulate_block(50.0, 1.0, 0.5, 20.0, 0.01);
        assert!(peak[0] > 20.0, "AP should overshoot 0 mV; peak={}", peak[0]);
    }

    #[test]
    fn implicit_cable_propagates_and_stays_bounded_at_100um() {
        // 200 × 100 µm — the case the v1 *explicit* solver diverged on. The
        // implicit solver must stay finite AND propagate an action potential.
        let mut cable = ImplicitCable::uniform(200, HhMembrane::at_rest(), 100.0, 238.0, 35.4);
        let peak = cable.stimulate_block(100.0, 1.0, 0.5, 15.0, 0.01);
        let maxv = peak.iter().cloned().fold(f64::MIN, f64::max);
        assert!(maxv.is_finite() && maxv < 80.0, "must stay bounded (no blow-up); maxv={maxv}");
        assert!(peak[150] > 0.0, "AP must propagate to compartment 150; peak={}", peak[150]);
    }

    #[test]
    fn implicit_cable_subthreshold_is_stable_at_100um() {
        // The exact configuration that blew the explicit solver to +inf: tiny
        // drive at 100 µm. The implicit solver must simply stay near rest.
        let mut cable = ImplicitCable::uniform(200, HhMembrane::at_rest(), 100.0, 238.0, 35.4);
        let peak = cable.stimulate_block(0.5, 1.0, 0.5, 15.0, 0.01);
        let maxv = peak.iter().cloned().fold(f64::MIN, f64::max);
        assert!(maxv < -40.0, "sub-threshold must stay bounded + not fire; maxv={maxv}");
    }
}
