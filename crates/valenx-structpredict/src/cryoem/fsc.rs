//! **Feature 29 — resolution estimation (Fourier shell correlation).**
//!
//! How good is a cryo-EM map? The universal answer is the **Fourier
//! shell correlation (FSC)** and the **gold-standard criterion**.
//!
//! - The dataset's particles are split into two independent halves;
//!   each half is reconstructed into its own **half-map**.
//! - The two half-maps are compared in Fourier space: for each
//!   spherical shell of spatial frequency, the **correlation** of the
//!   two maps' Fourier coefficients in that shell is computed. This
//!   correlation-vs-frequency curve is the FSC.
//! - At low frequency the half-maps agree (FSC ≈ 1); at high
//!   frequency they decorrelate into noise (FSC → 0). The frequency
//!   at which the FSC drops through **0.143** is the reported
//!   resolution — the gold-standard threshold (Rosenthal & Henderson
//!   2003), valid precisely *because* the half-maps are independent.
//!
//! This module computes the FSC by a **discrete Fourier transform**
//! of the two volumes, shell-binning the coefficients, and finds the
//! 0.143-crossing. The DFT is the textbook definition (an `O(N²)`
//! transform per axis — correct, not FFT-fast); for the box sizes a
//! classical-v1 reconstruction produces this is entirely adequate.
//! FSC and the gold-standard criterion themselves are *exact*.

use serde::{Deserialize, Serialize};

use crate::cryoem::mrc::Volume3d;
use crate::error::{Result, StructPredictError};

/// A Fourier-shell-correlation curve.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct FscCurve {
    /// Shell spatial frequencies in cycles per ångström, one per
    /// shell (shell `0` is the lowest frequency).
    pub frequencies: Vec<f64>,
    /// The FSC value in each shell, `[-1, 1]`.
    pub fsc: Vec<f64>,
    /// The voxel size the curve was computed with (ångström).
    pub voxel_size: f64,
}

impl FscCurve {
    /// The resolution (ångström) at which the FSC drops through a
    /// threshold — the standard way to read a resolution off an FSC
    /// curve. Returns `None` if the curve never drops that low (the
    /// reconstruction is resolution-limited by the box, not the data).
    pub fn resolution_at(&self, threshold: f64) -> Option<f64> {
        for i in 1..self.fsc.len() {
            if self.fsc[i - 1] >= threshold && self.fsc[i] < threshold {
                // Linear interpolation of the crossing frequency.
                let f0 = self.frequencies[i - 1];
                let f1 = self.frequencies[i];
                let v0 = self.fsc[i - 1];
                let v1 = self.fsc[i];
                let t = if (v0 - v1).abs() > 1e-12 {
                    (v0 - threshold) / (v0 - v1)
                } else {
                    0.0
                };
                let freq = f0 + t * (f1 - f0);
                if freq > 1e-12 {
                    return Some(1.0 / freq);
                }
            }
        }
        None
    }
}

/// A complex number for the DFT (a tiny self-contained type — no
/// external complex crate needed).
#[derive(Copy, Clone)]
struct Complex {
    re: f64,
    im: f64,
}

impl Complex {
    fn zero() -> Self {
        Complex { re: 0.0, im: 0.0 }
    }
}

/// Computes the discrete Fourier transform of a real 3-D volume.
///
/// Returns the complex coefficients in the same `(z, y, x)` layout.
/// This is the separable `O(N⁴)` DFT (a 1-D DFT along each of the
/// three axes) — the textbook definition.
fn dft3(volume: &Volume3d) -> Vec<Complex> {
    let (nx, ny, nz) = (volume.nx, volume.ny, volume.nz);
    let mut data: Vec<Complex> = volume
        .data
        .iter()
        .map(|&v| Complex {
            re: v as f64,
            im: 0.0,
        })
        .collect();

    // 1-D DFT helper: transform a strided line of length `n`.
    let dft_line = |data: &mut Vec<Complex>, start: usize, stride: usize, n: usize| {
        let line: Vec<Complex> = (0..n).map(|i| data[start + i * stride]).collect();
        for k in 0..n {
            let mut acc = Complex::zero();
            for (j, c) in line.iter().enumerate() {
                let angle = -2.0 * std::f64::consts::PI * (k as f64) * (j as f64) / n as f64;
                let (s, co) = angle.sin_cos();
                acc.re += c.re * co - c.im * s;
                acc.im += c.re * s + c.im * co;
            }
            data[start + k * stride] = acc;
        }
    };

    // Transform along x (stride 1).
    for z in 0..nz {
        for y in 0..ny {
            let start = (z * ny + y) * nx;
            dft_line(&mut data, start, 1, nx);
        }
    }
    // Along y (stride nx).
    for z in 0..nz {
        for x in 0..nx {
            let start = z * ny * nx + x;
            dft_line(&mut data, start, nx, ny);
        }
    }
    // Along z (stride nx*ny).
    for y in 0..ny {
        for x in 0..nx {
            let start = y * nx + x;
            dft_line(&mut data, start, nx * ny, nz);
        }
    }
    data
}

/// Computes the Fourier shell correlation between two equally-sized
/// 3-D maps (typically the two gold-standard half-maps).
///
/// Both volumes are DFT-transformed; for each integer-radius shell in
/// Fourier space the correlation
/// `Σ F₁·conj(F₂) / √(Σ|F₁|²·Σ|F₂|²)` is computed.
///
/// # Errors
/// [`StructPredictError::Invalid`] if the volumes differ in shape or
/// are smaller than `4³`.
pub fn fourier_shell_correlation(map_a: &Volume3d, map_b: &Volume3d) -> Result<FscCurve> {
    if map_a.nx != map_b.nx || map_a.ny != map_b.ny || map_a.nz != map_b.nz {
        return Err(StructPredictError::invalid(
            "map_b",
            "the two maps must have identical dimensions",
        ));
    }
    if map_a.nx < 4 || map_a.ny < 4 || map_a.nz < 4 {
        return Err(StructPredictError::invalid(
            "map_a",
            "FSC needs maps of at least 4³ voxels",
        ));
    }
    let (nx, ny, nz) = (map_a.nx, map_a.ny, map_a.nz);
    let fa = dft3(map_a);
    let fb = dft3(map_b);

    // The maximum shell radius — half the smallest box dimension.
    let max_r = (nx.min(ny).min(nz) / 2).max(1);
    let mut cross = vec![0.0f64; max_r + 1]; // Re(Σ Fa·conj(Fb))
    let mut power_a = vec![0.0f64; max_r + 1];
    let mut power_b = vec![0.0f64; max_r + 1];

    let cx = nx as i64 / 2;
    let cy = ny as i64 / 2;
    let cz = nz as i64 / 2;
    for z in 0..nz {
        // Wrapped frequency index (centre the DC term).
        let fz = wrap_freq(z as i64, nz as i64, cz);
        for y in 0..ny {
            let fy = wrap_freq(y as i64, ny as i64, cy);
            for x in 0..nx {
                let fx = wrap_freq(x as i64, nx as i64, cx);
                let r = ((fx * fx + fy * fy + fz * fz) as f64).sqrt().round() as usize;
                if r > max_r {
                    continue;
                }
                let idx = (z * ny + y) * nx + x;
                let a = fa[idx];
                let b = fb[idx];
                // Re(a · conj(b)).
                cross[r] += a.re * b.re + a.im * b.im;
                power_a[r] += a.re * a.re + a.im * a.im;
                power_b[r] += b.re * b.re + b.im * b.im;
            }
        }
    }

    let mut frequencies = Vec::with_capacity(max_r + 1);
    let mut fsc = Vec::with_capacity(max_r + 1);
    let voxel_size = if map_a.voxel_size > 0.0 {
        map_a.voxel_size
    } else {
        1.0
    };
    // The Nyquist frequency corresponds to shell `box/2`.
    let box_dim = nx.min(ny).min(nz) as f64;
    for r in 0..=max_r {
        // Spatial frequency of shell r: r / (box · voxel_size).
        let freq = r as f64 / (box_dim * voxel_size);
        frequencies.push(freq);
        let denom = (power_a[r] * power_b[r]).sqrt();
        fsc.push(if denom > 1e-12 {
            (cross[r] / denom).clamp(-1.0, 1.0)
        } else {
            // An empty shell (e.g. r=0 with zero-mean maps): perfectly
            // correlated by convention.
            1.0
        });
    }

    Ok(FscCurve {
        frequencies,
        fsc,
        voxel_size,
    })
}

/// Maps an array index to a signed, centred frequency index.
fn wrap_freq(i: i64, n: i64, _center: i64) -> i64 {
    if i <= n / 2 {
        i
    } else {
        i - n
    }
}

/// Estimates a map's resolution by the **gold-standard criterion**:
/// the resolution at which the FSC of the two independent half-maps
/// drops through 0.143.
///
/// Returns the resolution in ångström. If the FSC never drops to
/// 0.143 the maps agree to the box's Nyquist limit and that limiting
/// resolution is returned.
///
/// # Errors
/// [`StructPredictError::Invalid`] for mismatched / too-small maps.
pub fn gold_standard_resolution(half_map_a: &Volume3d, half_map_b: &Volume3d) -> Result<f64> {
    let curve = fourier_shell_correlation(half_map_a, half_map_b)?;
    if let Some(res) = curve.resolution_at(0.143) {
        Ok(res)
    } else {
        // Never drops below 0.143 → resolution-limited by the box.
        // The best resolvable frequency is the last (Nyquist) shell.
        let nyq = curve.frequencies.last().copied().unwrap_or(0.0);
        if nyq > 1e-12 {
            Ok(1.0 / nyq)
        } else {
            Err(StructPredictError::invalid(
                "half_map_a",
                "degenerate FSC curve",
            ))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A volume with a smooth low-frequency density bump.
    fn smooth_blob(n: usize) -> Volume3d {
        let mut v = Volume3d::zeros_cube(n);
        v.voxel_size = 1.0;
        let c = (n as f64 - 1.0) / 2.0;
        for z in 0..n {
            for y in 0..n {
                for x in 0..n {
                    let dx = x as f64 - c;
                    let dy = y as f64 - c;
                    let dz = z as f64 - c;
                    let r2 = dx * dx + dy * dy + dz * dz;
                    v.data[(z * n + y) * n + x] = (-r2 / (n as f64)).exp() as f32;
                }
            }
        }
        v
    }

    #[test]
    fn identical_maps_have_perfect_fsc() {
        let v = smooth_blob(8);
        let curve = fourier_shell_correlation(&v, &v).expect("fsc");
        // A map vs itself correlates perfectly in every shell.
        for &f in &curve.fsc {
            assert!(f > 0.999, "FSC shell {f} should be ~1");
        }
    }

    #[test]
    fn noise_decorrelates_at_high_frequency() {
        // Same low-frequency signal, independent high-frequency noise.
        let base = smooth_blob(10);
        let mut a = base.clone();
        let mut b = base.clone();
        let mut sa: u64 = 1;
        let mut sb: u64 = 999;
        for i in 0..a.data.len() {
            sa = sa.wrapping_mul(6364136223846793005).wrapping_add(1);
            sb = sb.wrapping_mul(6364136223846793005).wrapping_add(1);
            // High-frequency-ish noise: alternating sign per voxel.
            let sign = if i % 2 == 0 { 1.0 } else { -1.0 };
            a.data[i] += sign * ((sa >> 40) as f32 / (1u64 << 24) as f32 - 0.5);
            b.data[i] += sign * ((sb >> 40) as f32 / (1u64 << 24) as f32 - 0.5);
        }
        let curve = fourier_shell_correlation(&a, &b).expect("fsc");
        // The lowest-frequency shell (shared signal) correlates better
        // than the highest (independent noise).
        let low = curve.fsc[1];
        let high = *curve.fsc.last().unwrap();
        assert!(
            low > high,
            "low-freq FSC {low} should beat high-freq {high}"
        );
    }

    #[test]
    fn gold_standard_returns_a_resolution() {
        let base = smooth_blob(10);
        let mut a = base.clone();
        let mut b = base.clone();
        let mut sa: u64 = 7;
        let mut sb: u64 = 77;
        for i in 0..a.data.len() {
            sa = sa.wrapping_mul(6364136223846793005).wrapping_add(1);
            sb = sb.wrapping_mul(6364136223846793005).wrapping_add(1);
            let sign = if i % 2 == 0 { 1.0 } else { -1.0 };
            a.data[i] += sign * 0.8 * ((sa >> 40) as f32 / (1u64 << 24) as f32 - 0.5);
            b.data[i] += sign * 0.8 * ((sb >> 40) as f32 / (1u64 << 24) as f32 - 0.5);
        }
        let res = gold_standard_resolution(&a, &b).expect("resolution");
        // A finite, positive resolution in ångström.
        assert!(res.is_finite() && res > 0.0, "resolution {res}");
    }

    #[test]
    fn mismatched_maps_rejected() {
        let a = smooth_blob(8);
        let b = Volume3d::zeros_cube(10);
        assert!(fourier_shell_correlation(&a, &b).is_err());
    }

    #[test]
    fn too_small_maps_rejected() {
        let tiny = Volume3d::zeros_cube(2);
        assert!(fourier_shell_correlation(&tiny, &tiny).is_err());
    }
}
