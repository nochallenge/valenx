//! **Harmonic (frequency-response) analysis** — the steady-state response of a
//! damped linear structure to sinusoidal forcing across a frequency sweep.
//! This is the "harmonic response (full method)" analysis a production FEA
//! suite provides, and the piece the modal solver explicitly left out.
//!
//! Given the assembled stiffness `K` and mass `M` (see [`crate::assembly`]) and
//! a real force vector `F`, the steady-state amplitude `X(ω)` at angular
//! frequency `ω` solves the complex dynamic-stiffness system
//!
//! ```text
//! (K − ω²M + iωC) · X = F,   with Rayleigh damping  C = αM + βK.
//! ```
//!
//! We avoid complex arithmetic by solving the equivalent real `2n × 2n` system
//! for the in-phase (`x`) and quadrature (`y`) parts:
//!
//! ```text
//! [ A  −B ] [ x ]   [ F ]
//! [ B   A ] [ y ] = [ 0 ],   A = K − ω²M,  B = ωC.
//! ```
//!
//! Per-DOF amplitude is `√(xᵢ² + yᵢ²)`. Validated against the analytic
//! single-degree-of-freedom receptance (static limit `F/k`, resonance at
//! `ωₙ = √(k/m)` with peak amplitude `F/(2ζk)`).
//!
//! Honest scope: a direct (full-method) harmonic solve on dense `K`/`M` with
//! Rayleigh damping — research / preliminary-design grade, a step toward (not
//! an equal of) a production harmonic-response solver (no modal-superposition
//! acceleration, complex/structural damping, or large-sparse path yet).

use std::f64::consts::TAU;

use nalgebra::{DMatrix, DVector};

/// Steady-state frequency-response of a damped linear structure.
#[derive(Debug, Clone, PartialEq)]
pub struct HarmonicResponse {
    /// The swept excitation frequencies (Hz).
    pub frequencies_hz: Vec<f64>,
    /// `amplitude[f]` is the per-DOF response magnitude at `frequencies_hz[f]`.
    pub amplitude: Vec<DVector<f64>>,
}

impl HarmonicResponse {
    /// Response magnitude at `dof` across the whole sweep.
    pub fn dof_amplitude(&self, dof: usize) -> Vec<f64> {
        self.amplitude.iter().map(|a| a[dof]).collect()
    }

    /// `(frequency_hz, amplitude)` of the largest response at `dof` over the
    /// sweep — the dominant resonance peak. `None` for an empty sweep.
    pub fn resonance_peak(&self, dof: usize) -> Option<(f64, f64)> {
        self.frequencies_hz
            .iter()
            .zip(&self.amplitude)
            .map(|(&f, a)| (f, a[dof]))
            .fold(None, |best, (f, amp)| match best {
                Some((_, b)) if b >= amp => best,
                _ => Some((f, amp)),
            })
    }
}

/// Direct harmonic frequency-response: solve `(K − ω²M + iωC)·X = F` at each
/// frequency in `frequencies_hz`, with Rayleigh damping `C = αM + βK`, and
/// return the per-DOF amplitude at each frequency.
pub fn solve_harmonic(
    stiffness: &DMatrix<f64>,
    mass: &DMatrix<f64>,
    rayleigh_alpha: f64,
    rayleigh_beta: f64,
    force: &DVector<f64>,
    frequencies_hz: &[f64],
) -> HarmonicResponse {
    let n = force.len();
    let damping = mass * rayleigh_alpha + stiffness * rayleigh_beta;
    let mut amplitude = Vec::with_capacity(frequencies_hz.len());
    for &f in frequencies_hz {
        let w = TAU * f;
        let a = stiffness.clone() - mass * (w * w); // A = K − ω²M
        let b = &damping * w; // B = ωC
        let neg_b = -&b;
        let mut big = DMatrix::zeros(2 * n, 2 * n);
        big.view_mut((0, 0), (n, n)).copy_from(&a);
        big.view_mut((0, n), (n, n)).copy_from(&neg_b);
        big.view_mut((n, 0), (n, n)).copy_from(&b);
        big.view_mut((n, n), (n, n)).copy_from(&a);
        let mut rhs = DVector::zeros(2 * n);
        rhs.rows_mut(0, n).copy_from(force);
        let sol = big
            .lu()
            .solve(&rhs)
            .unwrap_or_else(|| DVector::zeros(2 * n));
        let mut amp = DVector::zeros(n);
        for i in 0..n {
            amp[i] = sol[i].hypot(sol[n + i]);
        }
        amplitude.push(amp);
    }
    HarmonicResponse {
        frequencies_hz: frequencies_hz.to_vec(),
        amplitude,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn linspace(start: f64, end: f64, n: usize) -> Vec<f64> {
        (0..n)
            .map(|i| start + (end - start) * i as f64 / (n - 1) as f64)
            .collect()
    }

    #[test]
    fn sdof_receptance_matches_analytic() {
        // m=1, k=400 → ωₙ=20 rad/s, fₙ=3.183 Hz. Rayleigh β=0.005, α=0 ⇒
        // ζ = βk/(2mωₙ) = 0.05.
        let m = DMatrix::from_element(1, 1, 1.0);
        let k = DMatrix::from_element(1, 1, 400.0);
        let force = DVector::from_element(1, 10.0);
        let (alpha, beta) = (0.0, 0.005);
        let freqs = linspace(0.1, 6.0, 500);
        let r = solve_harmonic(&k, &m, alpha, beta, &force, &freqs);

        // Static limit X(0⁺) ≈ F/k.
        let amp_low = r.dof_amplitude(0)[0];
        assert!((amp_low - 10.0 / 400.0).abs() < 2e-3, "static {amp_low}");

        // Resonance: peak near fₙ with amplitude ≈ F/(2ζk).
        let (fp, ap) = r.resonance_peak(0).expect("nonempty sweep");
        let fn_hz = 400.0_f64.sqrt() / TAU;
        assert!((fp - fn_hz).abs() / fn_hz < 0.05, "peak f {fp} vs {fn_hz}");
        let amp_res = 10.0 / (2.0 * 0.05 * 400.0); // = 0.25
        assert!(
            (ap - amp_res).abs() / amp_res < 0.10,
            "peak amp {ap} vs {amp_res}"
        );
    }

    #[test]
    fn two_dof_shows_two_resonances() {
        // Fixed–m–m–fixed spring chain: K = [[2k,−k],[−k,2k]], M = I.
        // Modes at ω=√(k/m) and √(3k/m) → two peaks in the sweep.
        let k = 400.0;
        let stiff = DMatrix::from_row_slice(2, 2, &[2.0 * k, -k, -k, 2.0 * k]);
        let mass = DMatrix::identity(2, 2);
        let force = DVector::from_vec(vec![10.0, 0.0]);
        let freqs = linspace(0.5, 8.0, 1200);
        let r = solve_harmonic(&stiff, &mass, 0.0, 0.0015, &force, &freqs);
        let amps = r.dof_amplitude(0);
        let peaks = (1..amps.len() - 1)
            .filter(|&i| amps[i] > amps[i - 1] && amps[i] > amps[i + 1])
            .count();
        assert!(peaks >= 2, "expected two resonances, found {peaks}");
    }

    #[test]
    fn is_deterministic() {
        let m = DMatrix::identity(2, 2);
        let k = DMatrix::from_row_slice(2, 2, &[800.0, -400.0, -400.0, 800.0]);
        let force = DVector::from_vec(vec![5.0, 0.0]);
        let freqs = linspace(0.5, 8.0, 200);
        let a = solve_harmonic(&k, &m, 0.0, 0.002, &force, &freqs);
        let b = solve_harmonic(&k, &m, 0.0, 0.002, &force, &freqs);
        assert_eq!(a, b);
    }
}
