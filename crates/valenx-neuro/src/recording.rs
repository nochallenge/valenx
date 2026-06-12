//! Extracellular recording — the forward model of a neural spike.
//!
//! Stimulation (the rest of this crate) drives neurons; **recording** is the
//! inverse, and the core of a read-out neural interface (Neuralink's N1 *reads*
//! spikes). A firing axon's transmembrane currents set up an extracellular
//! potential that an electrode measures. Each compartment's membrane current is
//! a point current source, so by superposition in the (linear) volume conductor:
//!
//! ```text
//! φ_e(r, t) = 1/(4πσ_e) · Σ_k I_m,k(t) / |r − r_k|
//! ```
//!
//! The membrane current of a compartment is computed as the divergence of the
//! axial current, so the currents **sum to zero** (charge conservation) — and
//! the recorded waveform is the biphasic/triphasic **extracellular action
//! potential** (EAP): microvolt-scale and sign-changing as the spike sweeps
//! past the electrode, exactly as seen on a real micro-electrode.

use crate::membrane::{HhMembrane, ImplicitCable};
use nalgebra::Vector3;
use std::f64::consts::PI;

/// A point electrode recording the extracellular field of a uniform
/// unmyelinated axon lying on the x-axis.
pub struct ExtracellularRecorder {
    sigma_e_s_m: f64,
    dx_m: f64,
    a_m: f64,
    ri_ohm_m: f64,
}

impl ExtracellularRecorder {
    /// Recorder for an axon of compartment length `dx_um`, radius `a_um`, and
    /// intracellular resistivity `ri_ohm_cm`, embedded in extracellular
    /// conductivity `sigma_e` (S/m).
    pub fn new(sigma_e: f64, dx_um: f64, a_um: f64, ri_ohm_cm: f64) -> Self {
        Self {
            sigma_e_s_m: sigma_e,
            dx_m: dx_um * 1.0e-6,
            a_m: a_um * 1.0e-6,
            ri_ohm_m: ri_ohm_cm * 0.01,
        }
    }

    /// Total axial conductance between adjacent compartments (siemens).
    fn g_axial_s(&self) -> f64 {
        PI * self.a_m * self.a_m / (self.ri_ohm_m * self.dx_m)
    }

    /// Trigger an action potential at the left end of an `n_comp`-compartment
    /// axon and record the extracellular potential at `electrode` (metres; the
    /// axon runs along x from 0 to (n−1)·dx) as the spike propagates past it.
    pub fn record(&self, n_comp: usize, electrode: Vector3<f64>) -> Recording {
        let dx_um = self.dx_m * 1.0e6;
        let a_um = self.a_m * 1.0e6;
        let ri_ohm_cm = self.ri_ohm_m * 100.0;
        let mut cable =
            ImplicitCable::uniform(n_comp, HhMembrane::at_rest(), dx_um, a_um, ri_ohm_cm);
        let g = self.g_axial_s();
        let dt = 0.01;
        let n_steps = (15.0_f64 / dt) as usize;
        let n_stim = (n_comp / 10).max(1);

        let mut eap_uv = Vec::with_capacity(n_steps);
        let mut max_abs_current_sum = 0.0_f64;
        let mut ext = vec![0.0; n_comp];
        let mut t = 0.0;
        for _ in 0..n_steps {
            let on = (1.0..1.5).contains(&t);
            for (k, e) in ext.iter_mut().enumerate() {
                *e = if on && k < n_stim { 100.0 } else { 0.0 };
            }
            cable.step(&ext, dt);
            t += dt;

            let v = cable.potentials();
            let mut phi = 0.0;
            let mut isum = 0.0;
            for k in 0..n_comp {
                // Outward membrane current = net axial current *in* (amperes);
                // this is the extracellular point-source strength, and Σ_k = 0.
                let mut dv = 0.0;
                if k > 0 {
                    dv += v[k - 1] - v[k];
                }
                if k + 1 < n_comp {
                    dv += v[k + 1] - v[k];
                }
                let i_m = g * dv * 1.0e-3; // mV → V → A
                isum += i_m;
                let x_k = k as f64 * self.dx_m;
                let r = ((electrode.x - x_k).powi(2)
                    + electrode.y * electrode.y
                    + electrode.z * electrode.z)
                    .sqrt();
                phi += i_m / r;
            }
            max_abs_current_sum = max_abs_current_sum.max(isum.abs());
            eap_uv.push(phi / (4.0 * PI * self.sigma_e_s_m) * 1.0e6); // V → µV
        }
        Recording {
            eap_uv,
            dt_ms: dt,
            max_abs_current_sum,
        }
    }
}

/// A recorded extracellular action potential.
pub struct Recording {
    /// EAP samples (µV).
    pub eap_uv: Vec<f64>,
    /// Sample interval (ms).
    pub dt_ms: f64,
    /// Largest `|Σ membrane current|` over the run (A) — ≈ 0 by charge
    /// conservation; a sanity check on the forward model.
    pub max_abs_current_sum: f64,
}

impl Recording {
    /// Most negative EAP sample (µV).
    pub fn min_uv(&self) -> f64 {
        self.eap_uv.iter().cloned().fold(f64::INFINITY, f64::min)
    }
    /// Most positive EAP sample (µV).
    pub fn max_uv(&self) -> f64 {
        self.eap_uv
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max)
    }
    /// Peak-to-peak amplitude (µV).
    pub fn peak_to_peak_uv(&self) -> f64 {
        self.max_uv() - self.min_uv()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn squid_recorder() -> ExtracellularRecorder {
        // σ_e 0.3 S/m; squid-like axon (matches the validated implicit cable).
        ExtracellularRecorder::new(0.3, 100.0, 238.0, 35.4)
    }

    #[test]
    fn membrane_currents_conserve_charge() {
        // Each membrane current is an axial-current divergence, so they must
        // sum to ~0 at every step — a fundamental check on the forward model.
        let rec = squid_recorder().record(200, Vector3::new(10.0e-3, 1.0e-3, 0.0));
        assert!(
            rec.max_abs_current_sum < 1.0e-12,
            "Σ membrane current must vanish (charge conservation); got {:.3e} A",
            rec.max_abs_current_sum
        );
    }

    #[test]
    fn eap_is_biphasic_with_dominant_negative_phase() {
        // As the spike sweeps past, the extracellular potential is biphasic
        // (sign-changing); the active depolarisation is a current sink, so the
        // dominant deflection is negative — the textbook EAP shape.
        let rec = squid_recorder().record(200, Vector3::new(10.0e-3, 1.0e-3, 0.0));
        let (lo, hi) = (rec.min_uv(), rec.max_uv());
        assert!(
            lo < 0.0 && hi > 0.0,
            "EAP must be biphasic: min={lo:.1} max={hi:.1} µV"
        );
        assert!(
            lo.abs() > hi,
            "negative (sink) phase dominates: |{lo:.1}| vs {hi:.1} µV"
        );
    }

    #[test]
    fn eap_amplitude_falls_off_with_distance() {
        let p2p = |off_mm: f64| {
            squid_recorder()
                .record(200, Vector3::new(10.0e-3, off_mm * 1.0e-3, 0.0))
                .peak_to_peak_uv()
        };
        let (near, mid, far) = (p2p(1.0), p2p(2.0), p2p(4.0));
        assert!(
            near > mid && mid > far,
            "EAP shrinks with electrode distance: {near:.0} > {mid:.0} > {far:.0} µV"
        );
    }
}
