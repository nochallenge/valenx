//! The naive (direct-summation) discrete Fourier transform and its
//! inverse, plus the spectrum helpers (bin frequencies, Nyquist,
//! magnitude).
//!
//! ## Conventions
//!
//! For an `N`-point signal `x[0..N]` the forward transform is
//!
//! ```text
//! X[k] = sum_{n=0}^{N-1} x[n] * exp(-2 pi i k n / N),   k = 0..N-1
//! ```
//!
//! and the inverse transform is
//!
//! ```text
//! x[n] = (1 / N) * sum_{k=0}^{N-1} X[k] * exp(+2 pi i k n / N).
//! ```
//!
//! The `1/N` normalisation lives entirely on the inverse side, so
//! `idft(dft(x)) == x` (up to floating-point round-off). This is the
//! same convention used by NumPy's `numpy.fft`.

use crate::complex::Complex;
use crate::error::{validate_sample_rate, FftError};
use std::f64::consts::PI;

/// Forward discrete Fourier transform of a real-valued signal.
///
/// Convenience wrapper around [`dft`] that lifts each real sample into
/// the complex plane (`im = 0`) first.
///
/// # Errors
///
/// Returns [`FftError::EmptySignal`] if `signal` is empty.
pub fn dft_real(signal: &[f64]) -> Result<Vec<Complex>, FftError> {
    let complex: Vec<Complex> = signal.iter().map(|&re| Complex::real(re)).collect();
    dft(&complex)
}

/// Forward discrete Fourier transform.
///
/// Computes `X[k] = sum_n x[n] exp(-2 pi i k n / N)` directly, in
/// `O(N^2)` time. This is the textbook definition, not a fast (FFT)
/// algorithm; it is exact (no radix restrictions) and works for any
/// length `N >= 1`.
///
/// # Errors
///
/// Returns [`FftError::EmptySignal`] if `signal` is empty.
pub fn dft(signal: &[Complex]) -> Result<Vec<Complex>, FftError> {
    let n = signal.len();
    if n == 0 {
        return Err(FftError::EmptySignal);
    }
    let n_f = n as f64;
    let mut out = Vec::with_capacity(n);
    for k in 0..n {
        let mut acc = Complex::ZERO;
        for (sample, n_idx) in signal.iter().zip(0..n) {
            // angle = -2 pi k n / N
            let angle = -2.0 * PI * (k as f64) * (n_idx as f64) / n_f;
            acc += *sample * Complex::expi(angle);
        }
        out.push(acc);
    }
    Ok(out)
}

/// Inverse discrete Fourier transform.
///
/// Computes `x[n] = (1/N) sum_k X[k] exp(+2 pi i k n / N)` directly, in
/// `O(N^2)` time. Applied to the output of [`dft`] it reconstructs the
/// original signal up to floating-point round-off.
///
/// # Errors
///
/// Returns [`FftError::EmptySignal`] if `spectrum` is empty.
pub fn idft(spectrum: &[Complex]) -> Result<Vec<Complex>, FftError> {
    let n = spectrum.len();
    if n == 0 {
        return Err(FftError::EmptySignal);
    }
    let n_f = n as f64;
    let mut out = Vec::with_capacity(n);
    for n_idx in 0..n {
        let mut acc = Complex::ZERO;
        for (coeff, k) in spectrum.iter().zip(0..n) {
            // angle = +2 pi k n / N
            let angle = 2.0 * PI * (k as f64) * (n_idx as f64) / n_f;
            acc += *coeff * Complex::expi(angle);
        }
        out.push(acc * (1.0 / n_f));
    }
    Ok(out)
}

/// Frequency, in hertz, of spectrum bin `k` for an `N`-point transform
/// sampled at `fs` hertz.
///
/// `f_k = k * fs / N`. Note that bins above `N/2` correspond to the
/// negative-frequency half of the spectrum; for a real input the
/// usual physical reading is `f_k = (k - N) * fs / N` there. This
/// helper returns the raw, un-aliased `k * fs / N`.
///
/// # Errors
///
/// Returns [`FftError::EmptySignal`] if `n == 0`,
/// [`FftError::InvalidSampleRate`] if `fs` is not finite and positive,
/// or [`FftError::BinOutOfRange`] if `k >= n`.
pub fn bin_frequency(k: usize, n: usize, fs: f64) -> Result<f64, FftError> {
    if n == 0 {
        return Err(FftError::EmptySignal);
    }
    let fs = validate_sample_rate(fs)?;
    if k >= n {
        return Err(FftError::BinOutOfRange { index: k, len: n });
    }
    Ok((k as f64) * fs / (n as f64))
}

/// All `N` raw bin frequencies `[0, fs/N, 2 fs/N, ..., (N-1) fs/N]`.
///
/// # Errors
///
/// Returns [`FftError::EmptySignal`] if `n == 0` or
/// [`FftError::InvalidSampleRate`] if `fs` is not finite and positive.
pub fn bin_frequencies(n: usize, fs: f64) -> Result<Vec<f64>, FftError> {
    if n == 0 {
        return Err(FftError::EmptySignal);
    }
    let fs = validate_sample_rate(fs)?;
    let n_f = n as f64;
    Ok((0..n).map(|k| (k as f64) * fs / n_f).collect())
}

/// The Nyquist frequency, `fs / 2` — the highest frequency that can be
/// represented without aliasing at sample rate `fs`.
///
/// # Errors
///
/// Returns [`FftError::InvalidSampleRate`] if `fs` is not finite and
/// positive.
pub fn nyquist_frequency(fs: f64) -> Result<f64, FftError> {
    let fs = validate_sample_rate(fs)?;
    Ok(fs / 2.0)
}

/// Magnitude spectrum `|X[k]|` of a transform, one modulus per bin.
///
/// Never fails for a non-empty input, but returns an empty vector for
/// an empty spectrum rather than erroring (the magnitude of "no bins"
/// is "no magnitudes").
pub fn magnitude_spectrum(spectrum: &[Complex]) -> Vec<f64> {
    spectrum.iter().map(|c| c.norm()).collect()
}

/// Total signal energy summed in the time domain, `sum_n |x[n]|^2`.
///
/// Pairs with [`spectral_energy`] to check Parseval's theorem.
pub fn time_energy(signal: &[Complex]) -> f64 {
    signal.iter().map(|c| c.norm_sqr()).sum()
}

/// Total spectral energy `(1/N) sum_k |X[k]|^2`.
///
/// Under the forward / inverse convention used here, Parseval's
/// theorem states `sum_n |x[n]|^2 == (1/N) sum_k |X[k]|^2`, so this
/// equals [`time_energy`] of the originating signal (up to round-off).
///
/// Returns `0.0` for an empty spectrum.
pub fn spectral_energy(spectrum: &[Complex]) -> f64 {
    let n = spectrum.len();
    if n == 0 {
        return 0.0;
    }
    let sum_sq: f64 = spectrum.iter().map(|c| c.norm_sqr()).sum();
    sum_sq / (n as f64)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Generous tolerance for accumulated `O(N^2)` round-off.
    const EPS: f64 = 1e-9;

    fn assert_close(a: f64, b: f64, eps: f64) {
        assert!(
            (a - b).abs() < eps,
            "expected {a} ~= {b} (|diff| = {} >= {eps})",
            (a - b).abs()
        );
    }

    fn assert_complex_close(a: Complex, b: Complex, eps: f64) {
        assert_close(a.re, b.re, eps);
        assert_close(a.im, b.im, eps);
    }

    // ---- input validation -------------------------------------------------

    #[test]
    fn empty_signal_errors() {
        assert_eq!(dft(&[]).unwrap_err(), FftError::EmptySignal);
        assert_eq!(idft(&[]).unwrap_err(), FftError::EmptySignal);
        assert_eq!(dft_real(&[]).unwrap_err(), FftError::EmptySignal);
        assert_eq!(bin_frequencies(0, 8.0).unwrap_err(), FftError::EmptySignal);
    }

    #[test]
    fn bad_sample_rate_errors() {
        assert!(matches!(
            nyquist_frequency(0.0),
            Err(FftError::InvalidSampleRate(_))
        ));
        assert!(matches!(
            nyquist_frequency(-4.0),
            Err(FftError::InvalidSampleRate(_))
        ));
        assert!(matches!(
            nyquist_frequency(f64::NAN),
            Err(FftError::InvalidSampleRate(_))
        ));
        assert!(matches!(
            bin_frequency(0, 4, f64::INFINITY),
            Err(FftError::InvalidSampleRate(_))
        ));
    }

    #[test]
    fn bin_out_of_range_errors() {
        // k == n is invalid (valid range is 0..n).
        assert_eq!(
            bin_frequency(4, 4, 8.0).unwrap_err(),
            FftError::BinOutOfRange { index: 4, len: 4 }
        );
        assert_eq!(bin_frequency(0, 4, 8.0).unwrap(), 0.0);
    }

    // ---- ground-truth: DFT of a delta is flat -----------------------------

    /// The DFT of a unit impulse at n = 0, `x = [1, 0, 0, ...]`, is the
    /// constant spectrum `X[k] = 1` for every `k` (a perfectly flat
    /// magnitude response).
    #[test]
    fn delta_has_flat_spectrum() {
        let n = 16;
        let mut x = vec![Complex::ZERO; n];
        x[0] = Complex::real(1.0);

        let spectrum = dft(&x).unwrap();
        assert_eq!(spectrum.len(), n);
        for bin in &spectrum {
            assert_complex_close(*bin, Complex::new(1.0, 0.0), EPS);
        }

        // Magnitude spectrum is identically 1.
        for m in magnitude_spectrum(&spectrum) {
            assert_close(m, 1.0, EPS);
        }
    }

    /// A delta delayed by one sample, `x = [0, 1, 0, ...]`, has unit
    /// magnitude in every bin (`|X[k]| = 1`) — the delay only rotates
    /// the phase, it does not change the flat magnitude.
    #[test]
    fn delayed_delta_has_flat_magnitude() {
        let n = 8;
        let mut x = vec![Complex::ZERO; n];
        x[1] = Complex::real(1.0);

        let spectrum = dft(&x).unwrap();
        for m in magnitude_spectrum(&spectrum) {
            assert_close(m, 1.0, EPS);
        }
        // Bin 0 (the DC / sum term) is exactly 1 + 0i.
        assert_complex_close(spectrum[0], Complex::new(1.0, 0.0), EPS);
    }

    // ---- ground-truth: forward then inverse round-trips -------------------

    #[test]
    fn dft_then_idft_round_trips() {
        // An arbitrary, non-trivial complex signal.
        let x: Vec<Complex> = (0..12)
            .map(|n| {
                let t = n as f64;
                Complex::new(0.3 * t - 1.1, (0.7 * t).sin())
            })
            .collect();

        let spectrum = dft(&x).unwrap();
        let recovered = idft(&spectrum).unwrap();

        assert_eq!(recovered.len(), x.len());
        for (orig, got) in x.iter().zip(recovered.iter()) {
            assert_complex_close(*orig, *got, EPS);
        }
    }

    #[test]
    fn idft_then_dft_round_trips() {
        // Round-trip the other way: start in the frequency domain.
        let spectrum: Vec<Complex> = (0..10)
            .map(|k| Complex::new((k as f64).cos(), 0.2 * k as f64))
            .collect();

        let time = idft(&spectrum).unwrap();
        let back = dft(&time).unwrap();

        for (orig, got) in spectrum.iter().zip(back.iter()) {
            assert_complex_close(*orig, *got, EPS);
        }
    }

    // ---- ground-truth: a sinusoid peaks at its bin ------------------------

    /// A real cosine at exactly bin `k0`,
    /// `x[n] = cos(2 pi k0 n / N)`, has all its energy in bins `k0`
    /// and `N - k0`, each with magnitude `N/2`. Every other bin is
    /// (numerically) zero. This is the canonical "spectral line"
    /// ground truth.
    #[test]
    fn cosine_peaks_at_its_bin() {
        let n = 32;
        let k0 = 5usize;
        let x: Vec<Complex> = (0..n)
            .map(|nn| {
                let phase = 2.0 * PI * (k0 as f64) * (nn as f64) / (n as f64);
                Complex::real(phase.cos())
            })
            .collect();

        let mag = magnitude_spectrum(&dft(&x).unwrap());

        // The two conjugate-symmetric peaks.
        let expected_peak = (n as f64) / 2.0;
        assert_close(mag[k0], expected_peak, 1e-7);
        assert_close(mag[n - k0], expected_peak, 1e-7);

        // Every other bin is essentially empty.
        for (k, m) in mag.iter().enumerate() {
            if k != k0 && k != n - k0 {
                assert_close(*m, 0.0, 1e-7);
            }
        }

        // The argmax bin is k0 within the first half [0, N/2).
        let (arg_max, _) = mag[..n / 2]
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
            .unwrap();
        assert_eq!(arg_max, k0);
    }

    /// A pure complex exponential `e^{+2 pi i k0 n / N}` puts ALL of
    /// its energy into the single bin `k0` (magnitude `N`), with no
    /// conjugate mirror — the cleanest single-line spectrum.
    #[test]
    fn complex_exponential_is_a_single_line() {
        let n = 24;
        let k0 = 7usize;
        let x: Vec<Complex> = (0..n)
            .map(|nn| {
                let angle = 2.0 * PI * (k0 as f64) * (nn as f64) / (n as f64);
                Complex::expi(angle)
            })
            .collect();

        let mag = magnitude_spectrum(&dft(&x).unwrap());
        for (k, m) in mag.iter().enumerate() {
            if k == k0 {
                assert_close(*m, n as f64, 1e-7);
            } else {
                assert_close(*m, 0.0, 1e-7);
            }
        }
    }

    // ---- ground-truth: bin freq = k fs / N --------------------------------

    #[test]
    fn bin_frequency_is_k_fs_over_n() {
        let n = 8;
        let fs = 16.0; // Hz
                       // f_k = k * fs / N = k * 2 Hz.
        for k in 0..n {
            let expected = (k as f64) * 2.0;
            assert_close(bin_frequency(k, n, fs).unwrap(), expected, EPS);
        }

        // The vector form agrees bin-for-bin.
        let freqs = bin_frequencies(n, fs).unwrap();
        assert_eq!(freqs.len(), n);
        for (k, f) in freqs.iter().enumerate() {
            assert_close(*f, (k as f64) * 2.0, EPS);
        }
    }

    #[test]
    fn dc_bin_is_zero_hz() {
        assert_close(bin_frequency(0, 100, 44_100.0).unwrap(), 0.0, EPS);
    }

    // ---- ground-truth: Nyquist = fs / 2 -----------------------------------

    #[test]
    fn nyquist_is_half_sample_rate() {
        assert_close(nyquist_frequency(44_100.0).unwrap(), 22_050.0, EPS);
        assert_close(nyquist_frequency(2.0).unwrap(), 1.0, EPS);

        // For even N, the Nyquist frequency is exactly the bin at k = N/2.
        let n = 8;
        let fs = 8.0;
        let nyq = nyquist_frequency(fs).unwrap();
        assert_close(bin_frequency(n / 2, n, fs).unwrap(), nyq, EPS);
    }

    // ---- ground-truth: Parseval energy ------------------------------------

    /// Parseval's theorem: `sum_n |x[n]|^2 == (1/N) sum_k |X[k]|^2`.
    #[test]
    fn parseval_energy_is_conserved() {
        let x: Vec<Complex> = (0..20)
            .map(|n| {
                let t = n as f64;
                Complex::new((0.4 * t).cos() + 0.5, (0.9 * t).sin() - 0.2)
            })
            .collect();

        let spectrum = dft(&x).unwrap();

        let e_time = time_energy(&x);
        let e_freq = spectral_energy(&spectrum);

        assert!(e_time > 0.0, "test signal must carry energy");
        assert_close(e_time, e_freq, 1e-7);
    }

    /// Parseval holds for the unit impulse too: time energy 1, and
    /// `(1/N) sum_k 1^2 = (1/N) * N = 1`.
    #[test]
    fn parseval_holds_for_delta() {
        let n = 10;
        let mut x = vec![Complex::ZERO; n];
        x[0] = Complex::real(1.0);

        let spectrum = dft(&x).unwrap();
        assert_close(time_energy(&x), 1.0, EPS);
        assert_close(spectral_energy(&spectrum), 1.0, EPS);
    }

    // ---- DC term equals the sum of samples --------------------------------

    #[test]
    fn dc_bin_equals_sum_of_samples() {
        let x = [1.0, 2.0, 3.0, 4.0, 5.0];
        let spectrum = dft_real(&x).unwrap();
        // X[0] = sum x[n] = 15, purely real.
        assert_complex_close(spectrum[0], Complex::new(15.0, 0.0), EPS);
    }

    #[test]
    fn single_sample_is_its_own_transform() {
        // N = 1: X[0] = x[0], and idft recovers it.
        let x = [Complex::new(3.0, -2.0)];
        let spectrum = dft(&x).unwrap();
        assert_complex_close(spectrum[0], x[0], EPS);
        let back = idft(&spectrum).unwrap();
        assert_complex_close(back[0], x[0], EPS);
    }
}
