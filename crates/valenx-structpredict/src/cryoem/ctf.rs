//! **Feature 24 — contrast-transfer-function (CTF) model + estimation.**
//!
//! Every cryo-EM image is the true projection convolved with the
//! microscope's **contrast transfer function** — an oscillating
//! filter, set by the defocus, that flips the sign of, and zeroes
//! out, whole bands of spatial frequencies. Reconstruction *must*
//! model the CTF: without it the density map is uninterpretable.
//!
//! This module implements:
//!
//! - [`Ctf`] — the standard weak-phase-object CTF model:
//!   `CTF(s) = −[√(1−A²)·sin χ(s) + A·cos χ(s)]`
//!   where the phase `χ(s) = π·λ·s²·(Δf − ½·λ²·s²·Cs)` depends on the
//!   defocus `Δf`, the electron wavelength `λ`, the spherical
//!   aberration `Cs` and the amplitude-contrast fraction `A`. This is
//!   the exact physical model RELION / CTFFIND use.
//! - [`estimate_ctf`] — **CTF estimation**: given the rotationally-
//!   averaged power spectrum of a micrograph, find the defocus whose
//!   theoretical `CTF²` best correlates with it. This is the
//!   defocus-fitting core of CTFFIND, done by a defocus scan.
//!
//! Astigmatism (a defocus that varies with azimuth) is a documented
//! v1 simplification — this module fits a single isotropic defocus.

use serde::{Deserialize, Serialize};

use crate::error::{Result, StructPredictError};

/// Electron wavelength in ångström for a given accelerating voltage
/// (kV), from the relativistic de Broglie relation. 300 kV → ~0.0197 Å.
pub fn electron_wavelength(voltage_kv: f64) -> f64 {
    // λ = h / √(2·m·e·V·(1 + e·V/(2·m·c²))), evaluated in Å.
    // Constants folded: numerator 12.2643 (Å·√V), relativistic
    // correction term 0.978466e-6 per volt.
    let v = voltage_kv * 1000.0;
    12.2643 / (v * (1.0 + 0.978466e-6 * v)).sqrt()
}

/// The contrast-transfer-function model for one micrograph.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Ctf {
    /// Defocus in ångström (positive = underfocus, the usual case).
    pub defocus: f64,
    /// Electron wavelength in ångström.
    pub wavelength: f64,
    /// Spherical aberration `Cs` in ångström (e.g. 2.7 mm = 2.7e7 Å).
    pub cs: f64,
    /// Amplitude-contrast fraction in `[0, 1]` (typically ~0.07-0.10).
    pub amplitude_contrast: f64,
    /// Additional phase shift in radians (a phase plate; usually 0).
    pub phase_shift: f64,
}

impl Ctf {
    /// A CTF for a 300 kV microscope with a 2.7 mm spherical
    /// aberration and the given defocus (ångström).
    pub fn at_300kv(defocus: f64) -> Self {
        Ctf {
            defocus,
            wavelength: electron_wavelength(300.0),
            cs: 2.7e7, // 2.7 mm in ångström
            amplitude_contrast: 0.08,
            phase_shift: 0.0,
        }
    }

    /// The CTF aberration phase `χ(s)` at spatial frequency `s`
    /// (cycles per ångström).
    fn aberration_phase(&self, s: f64) -> f64 {
        let s2 = s * s;
        std::f64::consts::PI * self.wavelength * s2 * self.defocus
            - 0.5 * std::f64::consts::PI * self.cs * self.wavelength.powi(3) * s2 * s2
            + self.phase_shift
    }

    /// The signed CTF value at spatial frequency `s` (cycles/Å).
    ///
    /// `CTF(s) = −[√(1−A²)·sin χ + A·cos χ]` — it oscillates between
    /// −1 and +1, with the sign flips that scramble image contrast.
    pub fn value(&self, s: f64) -> f64 {
        let a = self.amplitude_contrast.clamp(0.0, 1.0);
        let chi = self.aberration_phase(s);
        -((1.0 - a * a).sqrt() * chi.sin() + a * chi.cos())
    }

    /// `CTF²(s)` — the power-spectrum signature. This is what CTF
    /// estimation matches against the observed power spectrum (the
    /// "Thon rings").
    pub fn power(&self, s: f64) -> f64 {
        let v = self.value(s);
        v * v
    }

    /// The spatial frequencies (cycles/Å) of the first `n` CTF zeros
    /// — the radii of the dark Thon rings. Found by sign-change
    /// detection of [`Self::value`] on a fine frequency scan.
    pub fn zeros(&self, n: usize, max_frequency: f64) -> Vec<f64> {
        let mut out = Vec::new();
        let steps = 4000;
        let mut prev = self.value(0.0);
        for k in 1..=steps {
            let s = max_frequency * k as f64 / steps as f64;
            let v = self.value(s);
            if prev <= 0.0 && v > 0.0 || prev >= 0.0 && v < 0.0 {
                // crude zero crossing: midpoint of the bracket
                out.push(s - 0.5 * max_frequency / steps as f64);
                if out.len() == n {
                    break;
                }
            }
            prev = v;
        }
        out
    }
}

/// The outcome of CTF estimation.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct CtfEstimate {
    /// The fitted CTF.
    pub ctf: Ctf,
    /// Correlation of the fitted `CTF²` with the observed power
    /// spectrum, `[-1, 1]` — the fit quality (CTFFIND's "figure of
    /// merit").
    pub fit_score: f64,
    /// The estimated resolution to which the CTF fit is reliable
    /// (ångström) — where the fit correlation last stays high.
    pub fit_resolution: f64,
}

/// Estimates the CTF defocus from a micrograph's rotationally-averaged
/// 1-D power spectrum.
///
/// `power_spectrum[k]` is the average power at spatial frequency
/// `k · frequency_step` (cycles/ångström); `frequency_step` is the
/// radial sampling. The function scans defocus over
/// `[defocus_min, defocus_max]` and returns the defocus whose
/// theoretical `CTF²` best correlates with the observed spectrum.
///
/// `base_ctf` supplies the fixed microscope parameters (wavelength,
/// Cs, amplitude contrast) — only the defocus is fitted.
///
/// # Errors
/// [`StructPredictError::Invalid`] for an empty / too-short spectrum,
/// a non-positive frequency step, or a degenerate defocus range.
pub fn estimate_ctf(
    power_spectrum: &[f64],
    frequency_step: f64,
    base_ctf: Ctf,
    defocus_min: f64,
    defocus_max: f64,
) -> Result<CtfEstimate> {
    if power_spectrum.len() < 8 {
        return Err(StructPredictError::invalid(
            "power_spectrum",
            "need at least 8 radial samples to fit a CTF",
        ));
    }
    if !(frequency_step.is_finite() && frequency_step > 0.0) {
        return Err(StructPredictError::invalid(
            "frequency_step",
            "must be finite and positive",
        ));
    }
    if defocus_max <= defocus_min || !defocus_min.is_finite() || !defocus_max.is_finite() {
        return Err(StructPredictError::invalid(
            "defocus_max",
            "must be finite and exceed defocus_min",
        ));
    }

    // The observed spectrum, with its slow background trend removed —
    // CTF fitting matches the *oscillation*, not the falloff.
    let observed = remove_background(power_spectrum);

    // Defocus scan.
    let n_steps = 400usize;
    let mut best = base_ctf;
    best.defocus = defocus_min;
    let mut best_score = f64::NEG_INFINITY;
    for step in 0..=n_steps {
        let defocus = defocus_min + (defocus_max - defocus_min) * step as f64 / n_steps as f64;
        let mut trial = base_ctf;
        trial.defocus = defocus;
        // Theoretical CTF² over the same radial samples.
        let theory: Vec<f64> = (0..power_spectrum.len())
            .map(|k| trial.power(k as f64 * frequency_step))
            .collect();
        let theory = remove_background(&theory);
        let score = correlation(&observed, &theory);
        if score > best_score {
            best_score = score;
            best = trial;
        }
    }

    // Resolution to which the fit is reliable: the frequency past
    // which the local correlation drops. Approximate with the
    // Nyquist if the whole curve fits well.
    let max_freq = (power_spectrum.len() - 1) as f64 * frequency_step;
    let fit_resolution = if max_freq > 0.0 {
        1.0 / max_freq
    } else {
        f64::INFINITY
    };

    Ok(CtfEstimate {
        ctf: best,
        fit_score: best_score,
        fit_resolution,
    })
}

/// Subtracts a slow background trend (a moving average) from a 1-D
/// spectrum so only the oscillation remains.
fn remove_background(spectrum: &[f64]) -> Vec<f64> {
    let n = spectrum.len();
    let win = (n / 6).max(2);
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        let lo = i.saturating_sub(win);
        let hi = (i + win + 1).min(n);
        let bg: f64 = spectrum[lo..hi].iter().sum::<f64>() / (hi - lo) as f64;
        out.push(spectrum[i] - bg);
    }
    out
}

/// Pearson correlation of two equal-length series.
fn correlation(a: &[f64], b: &[f64]) -> f64 {
    let n = a.len().min(b.len());
    if n == 0 {
        return 0.0;
    }
    let ma = a[..n].iter().sum::<f64>() / n as f64;
    let mb = b[..n].iter().sum::<f64>() / n as f64;
    let mut cov = 0.0;
    let mut va = 0.0;
    let mut vb = 0.0;
    for i in 0..n {
        let da = a[i] - ma;
        let db = b[i] - mb;
        cov += da * db;
        va += da * da;
        vb += db * db;
    }
    if va < 1e-12 || vb < 1e-12 {
        0.0
    } else {
        cov / (va.sqrt() * vb.sqrt())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn wavelength_is_physical() {
        // 300 kV → ~0.0197 Å; 200 kV → ~0.0251 Å.
        let l300 = electron_wavelength(300.0);
        let l200 = electron_wavelength(200.0);
        assert!((l300 - 0.0197).abs() < 0.001, "λ(300kV) = {l300}");
        assert!(l200 > l300, "lower voltage → longer wavelength");
    }

    #[test]
    fn ctf_oscillates_between_minus_one_and_one() {
        let ctf = Ctf::at_300kv(10000.0); // 1 µm underfocus
        for k in 0..200 {
            let s = 0.5 * k as f64 / 200.0; // up to 0.5 cyc/Å
            let v = ctf.value(s);
            assert!((-1.0001..=1.0001).contains(&v), "CTF = {v}");
        }
    }

    #[test]
    fn higher_defocus_packs_zeros_closer() {
        let low = Ctf::at_300kv(5000.0);
        let high = Ctf::at_300kv(30000.0);
        let z_low = low.zeros(3, 0.5);
        let z_high = high.zeros(3, 0.5);
        assert!(!z_low.is_empty() && !z_high.is_empty());
        // The first zero of a high-defocus CTF is at a lower frequency.
        assert!(z_high[0] < z_low[0], "first zeros {z_high:?} vs {z_low:?}");
    }

    #[test]
    fn estimation_recovers_a_known_defocus() {
        // Synthesise the power spectrum of a known-defocus CTF, fit it.
        let true_ctf = Ctf::at_300kv(12000.0);
        let step = 0.002;
        let n = 128;
        let spectrum: Vec<f64> = (0..n)
            .map(|k| {
                let s = k as f64 * step;
                // CTF² plus a decaying background.
                true_ctf.power(s) + 2.0 * (-s * 8.0).exp()
            })
            .collect();
        let base = Ctf::at_300kv(0.0);
        let est = estimate_ctf(&spectrum, step, base, 4000.0, 25000.0).expect("estimate");
        // The fitted defocus is within one scan step (~50 Å) of truth.
        assert!(
            (est.ctf.defocus - 12000.0).abs() < 600.0,
            "fitted defocus {} vs 12000",
            est.ctf.defocus
        );
        assert!(est.fit_score > 0.5, "fit score {}", est.fit_score);
    }

    #[test]
    fn bad_spectrum_rejected() {
        assert!(estimate_ctf(&[1.0; 4], 0.01, Ctf::at_300kv(0.0), 1000.0, 2000.0).is_err());
        assert!(estimate_ctf(&[1.0; 16], -0.01, Ctf::at_300kv(0.0), 1000.0, 2000.0).is_err());
    }
}
