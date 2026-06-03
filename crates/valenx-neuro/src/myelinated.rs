//! Myelinated mammalian nerve fiber — saltatory conduction.
//!
//! A McNeal/CRRSS-class compartmental model: active **nodes of Ranvier**
//! (Hodgkin–Huxley kinetics) separated by passive **myelinated internodes**
//! (low capacitance, low leak). Internode length scales with fiber diameter
//! (L ≈ 100·D), which makes conduction velocity scale ~linearly with diameter
//! — the hallmark of myelinated fibers. Solved with the implicit integrator.
//!
//! Honest scope: this is the node + single-internode simplification, not the
//! full MRG double-cable (periaxonal space, MYSA/FLUT/STIN segments). It
//! captures saltatory conduction and the CV–diameter scaling; node kinetics are
//! HH (squid-derived) standing in for a mammalian-node model. Units inside the
//! solver: V mV, t ms, capacitance nF, current nA, axial conductance µS.

use crate::cable::V_REST;
use crate::membrane::{thomas, HhMembrane};
use std::f64::consts::PI;

#[derive(Clone)]
enum Seg {
    Node(HhMembrane),
    Internode,
}

impl Seg {
    /// Specific membrane capacitance (µF/cm²). Myelin makes the internode
    /// nearly transparent: with many membrane layers in series the effective
    /// capacitance per area is tiny, so the (long) internode does not load the
    /// cable and the action potential jumps node-to-node.
    fn c_m(&self) -> f64 {
        match self {
            Seg::Node(_) => 1.0,
            Seg::Internode => 1.0e-4,
        }
    }
    /// Ionic current density (µA/cm²).
    fn i_ion(&self, v: f64) -> f64 {
        match self {
            Seg::Node(hh) => {
                use crate::membrane::Membrane;
                hh.ionic_current(v)
            }
            Seg::Internode => 1.0e-5 * (v - V_REST), // negligible myelin leak
        }
    }
    fn advance(&mut self, v: f64, dt: f64) {
        if let Seg::Node(hh) = self {
            use crate::membrane::Membrane;
            hh.advance_gates(v, dt);
        }
    }
}

/// A uniform myelinated fiber: `n_nodes` active nodes joined by passive
/// internodes, with geometry scaled by the fiber diameter.
pub struct MyelinatedFiber {
    v: Vec<f64>,
    seg: Vec<Seg>,
    area_cm2: Vec<f64>,
    g_link_us: Vec<f64>,
    node_idx: Vec<usize>,
    internode_len_cm: f64,
}

impl MyelinatedFiber {
    /// Build a fiber of outer diameter `d_um` (µm) with `n_nodes` nodes.
    pub fn new(d_um: f64, n_nodes: usize) -> Self {
        let ri = 70.0_f64; // Ω·cm intracellular resistivity (mammalian)
        let d_axon_cm = 0.7 * d_um * 1.0e-4; // axon diameter under myelin (cm)
        let a_axon_cm = d_axon_cm / 2.0;
        let node_len_cm = 1.0e-4; // 1 µm node gap
        let inter_len_cm = 100.0 * d_um * 1.0e-4; // L ≈ 100·D µm

        let n = 2 * n_nodes - 1;
        let mut seg = Vec::with_capacity(n);
        let mut area = Vec::with_capacity(n);
        let mut node_idx = Vec::new();
        for k in 0..n {
            if k % 2 == 0 {
                seg.push(Seg::Node(HhMembrane::at_rest()));
                node_idx.push(k);
                area.push(PI * d_axon_cm * node_len_cm);
            } else {
                seg.push(Seg::Internode);
                area.push(PI * d_axon_cm * inter_len_cm);
            }
        }
        let cross = PI * a_axon_cm * a_axon_cm; // axoplasm cross-section (cm²)
        let half_len = |k: usize| {
            if k % 2 == 0 {
                node_len_cm / 2.0
            } else {
                inter_len_cm / 2.0
            }
        };
        let mut g_link = Vec::with_capacity(n - 1);
        for k in 0..n - 1 {
            let r = ri * (half_len(k) + half_len(k + 1)) / cross; // Ω
            g_link.push(1.0e6 / r); // µS
        }
        Self {
            v: vec![V_REST; n],
            seg,
            area_cm2: area,
            g_link_us: g_link,
            node_idx,
            internode_len_cm: inter_len_cm,
        }
    }

    fn step(&mut self, stim_na: &[f64], dt: f64) {
        let n = self.v.len();
        for k in 0..n {
            self.seg[k].advance(self.v[k], dt);
        }
        let mut sub = vec![0.0; n];
        let mut diag = vec![0.0; n];
        let mut sup = vec![0.0; n];
        let mut rhs = vec![0.0; n];
        for k in 0..n {
            let cap_nf = self.seg[k].c_m() * self.area_cm2[k] * 1000.0;
            let i_ion_na = self.seg[k].i_ion(self.v[k]) * self.area_cm2[k] * 1000.0;
            let gl = if k > 0 { self.g_link_us[k - 1] } else { 0.0 };
            let gr = if k + 1 < n { self.g_link_us[k] } else { 0.0 };
            sub[k] = -gl * dt;
            sup[k] = -gr * dt;
            diag[k] = cap_nf + (gl + gr) * dt;
            rhs[k] = cap_nf * self.v[k] + dt * (stim_na[k] - i_ion_na);
        }
        self.v = thomas(&sub, &diag, &sup, &rhs);
    }

    /// Stimulate the first node, integrate, and return the conduction velocity
    /// (m/s) between two downstream nodes — `None` if it does not reach both.
    pub fn conduction_velocity(&mut self) -> Option<f64> {
        let n = self.v.len();
        let dt = 0.002;
        let steps = (5.0_f64 / dt).round() as usize;
        let mut peak_v = self.v.clone();
        let mut peak_t = vec![0.0; n];
        let mut stim = vec![0.0; n];
        let mut t = 0.0;
        for _ in 0..steps {
            stim[0] = if t < 0.1 { 5.0 } else { 0.0 }; // nA into node 0
            self.step(&stim, dt);
            t += dt;
            for k in 0..n {
                if self.v[k] > peak_v[k] {
                    peak_v[k] = self.v[k];
                    peak_t[k] = t;
                }
            }
        }
        let nn = self.node_idx.len();
        let a = self.node_idx[nn / 3];
        let b = self.node_idx[2 * nn / 3];
        if peak_v[a] < 0.0 || peak_v[b] < 0.0 || peak_t[b] <= peak_t[a] {
            return None;
        }
        let n_internodes = ((b - a) / 2) as f64;
        let dist_cm = n_internodes * self.internode_len_cm;
        // cm/ms → m/s is ×10.
        Some(dist_cm / (peak_t[b] - peak_t[a]) * 10.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn conduction_velocity_matches_six_times_diameter_rule() {
        let cv10 = MyelinatedFiber::new(10.0, 21)
            .conduction_velocity()
            .expect("10 µm fiber should fire + propagate");
        let cv20 = MyelinatedFiber::new(20.0, 21)
            .conduction_velocity()
            .expect("20 µm fiber should fire + propagate");
        // Empirical mammalian rule (Hursh 1939 / Rushton 1951): CV ≈ 6·D m/s
        // for D in µm. The simplified node+internode model reproduces it well.
        for (d, cv) in [(10.0_f64, cv10), (20.0_f64, cv20)] {
            let expected = 6.0 * d;
            assert!(
                (cv - expected).abs() / expected < 0.25,
                "CV for {d} µm should be ≈ 6·D = {expected:.0} m/s; got {cv:.1}"
            );
        }
        // CV ∝ D (saltatory) — not ∝ √D as for the unmyelinated v1 cable.
        assert!(
            (cv20 / cv10 - 2.0).abs() < 0.3,
            "CV should ~double when diameter doubles: {cv10:.1} -> {cv20:.1} m/s"
        );
    }
}
