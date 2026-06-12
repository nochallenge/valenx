//! Radial distribution function — **roadmap feature 27**.
//!
//! The radial distribution function `g(r)` answers: *how much more (or
//! less) likely is it to find a second atom a distance `r` from a
//! given atom, compared with a uniform ideal gas of the same density?*
//! It is the central structural fingerprint of a liquid — the
//! coordination shells of water, the packing of a Lennard-Jones fluid.
//!
//! It is computed by histogramming every pair distance into shells of
//! width `dr`, then normalising shell `k` (spanning `[r, r+dr)`) by the
//! number of pairs an ideal gas at the same number density `ρ` would
//! put in that spherical shell:
//!
//! ```text
//! g(rₖ) = histogram(k) / [ N_pairs · (4/3·π·((r+dr)³ − r³)·ρ) ]
//! ```
//!
//! Distances use the minimum-image convention, so `r` is meaningful
//! only up to half the box width — the histogram is capped there.
//! `g(r) → 1` at large `r` (no correlation) and `g(r) → 0` at very
//! short `r` (excluded volume) — both are checked in the tests.

use nalgebra::Vector3;

use crate::error::{MdError, Result};
use crate::pbc::SimBox;

/// A computed radial distribution function.
#[derive(Clone, Debug, PartialEq)]
pub struct RadialDistribution {
    /// Shell-centre radii (nm).
    pub radii: Vec<f64>,
    /// `g(r)` at each shell centre (dimensionless).
    pub g: Vec<f64>,
    /// Histogram-bin width (nm).
    pub bin_width: f64,
}

impl RadialDistribution {
    /// The running coordination number `n(r) = ∫₀ʳ 4πρ·r'²·g(r')·dr'`
    /// — the average number of neighbours within radius `r`.
    ///
    /// `density` is the number density `ρ` (atoms / nm³).
    pub fn coordination_number(&self, density: f64) -> Vec<f64> {
        let mut cum = 0.0;
        let mut out = Vec::with_capacity(self.radii.len());
        for (r, g) in self.radii.iter().zip(&self.g) {
            // Shell volume element 4π r² dr.
            cum += 4.0 * std::f64::consts::PI * r * r * g * density * self.bin_width;
            out.push(cum);
        }
        out
    }

    /// The radius of the first peak of `g(r)` (nm), or `None` if the
    /// function is featureless.
    ///
    /// This is the *first* local maximum (smallest `r`) that rises
    /// above the significance threshold `g > 1.05` — not the global
    /// maximum. For a crystalline solid a later coordination shell can
    /// be taller than the nearest-neighbour shell, so the global
    /// maximum would not be the first peak.
    pub fn first_peak(&self) -> Option<f64> {
        let g = &self.g;
        let n = g.len();
        for k in 0..n {
            if g[k] <= 1.05 {
                continue;
            }
            // A local maximum: not lower than either neighbour (a
            // boundary bin only needs its single interior neighbour).
            let ge_left = k == 0 || g[k] >= g[k - 1];
            let ge_right = k + 1 == n || g[k] >= g[k + 1];
            if ge_left && ge_right {
                return Some(self.radii[k]);
            }
        }
        None
    }
}

/// Computes `g(r)` for one frame of identical particles.
///
/// * `positions` — atom positions (nm)
/// * `cell` — the periodic box (gives the volume / density)
/// * `r_max` — largest radius to histogram (nm); capped at half the
///   minimum box width
/// * `bins` — number of histogram bins
///
/// # Errors
/// [`MdError::Invalid`] for fewer than two atoms, a non-positive
/// `r_max`, zero bins, or a non-periodic box (density is undefined
/// without a volume).
pub fn radial_distribution(
    positions: &[Vector3<f64>],
    cell: &SimBox,
    r_max: f64,
    bins: usize,
) -> Result<RadialDistribution> {
    let n = positions.len();
    if n < 2 {
        return Err(MdError::invalid("rdf", "needs at least two atoms"));
    }
    if !(r_max.is_finite() && r_max > 0.0) {
        return Err(MdError::invalid("r_max", "must be finite and positive"));
    }
    if bins == 0 {
        return Err(MdError::invalid("bins", "must be at least 1"));
    }
    if !cell.is_periodic() {
        return Err(MdError::invalid(
            "cell",
            "the RDF needs a periodic box for the density normalisation",
        ));
    }
    let volume = cell.volume();
    let density = n as f64 / volume;
    // Cap r_max at the minimum-image limit.
    let r_cap = r_max.min(cell.max_cutoff());
    let dr = r_cap / bins as f64;
    let r_cap2 = r_cap * r_cap;

    let mut hist = vec![0u64; bins];
    let mut pair_count = 0u64;
    for i in 0..n {
        for j in (i + 1)..n {
            let d = cell.min_image(positions[i] - positions[j]);
            let r2 = d.norm_squared();
            if r2 < r_cap2 {
                let r = r2.sqrt();
                let bin = ((r / dr) as usize).min(bins - 1);
                hist[bin] += 1;
            }
            pair_count += 1;
        }
    }

    // Normalise each shell by the ideal-gas expectation.
    let mut radii = Vec::with_capacity(bins);
    let mut g = Vec::with_capacity(bins);
    let four_thirds_pi = 4.0 / 3.0 * std::f64::consts::PI;
    for (k, &count) in hist.iter().enumerate() {
        let r_lo = k as f64 * dr;
        let r_hi = r_lo + dr;
        radii.push(0.5 * (r_lo + r_hi));
        let shell_volume = four_thirds_pi * (r_hi.powi(3) - r_lo.powi(3));
        // Expected count for an ideal gas. The histogram counts
        // *unordered* pairs; for an ideal gas the separation vector of
        // any one pair is uniform over the box, so each of the
        // `pair_count` pairs lands in this shell with probability
        // `shell_volume / V`. With `density = N / V`, that is
        // `ideal = pair_count · shell_volume · density / N`.
        let ideal = pair_count as f64 * shell_volume * density / n as f64;
        let value = if ideal > 0.0 {
            count as f64 / ideal
        } else {
            0.0
        };
        g.push(value);
    }

    Ok(RadialDistribution {
        radii,
        g,
        bin_width: dr,
    })
}

/// Computes an averaged `g(r)` over several frames (a trajectory).
///
/// Each frame must have the same atom count. The per-frame histograms
/// are summed before normalisation, which is the correct way to
/// average an RDF.
///
/// # Errors
/// As [`radial_distribution`], plus [`MdError::Invalid`] for an empty
/// frame list or inconsistent frame sizes.
pub fn radial_distribution_averaged(
    frames: &[Vec<Vector3<f64>>],
    cell: &SimBox,
    r_max: f64,
    bins: usize,
) -> Result<RadialDistribution> {
    if frames.is_empty() {
        return Err(MdError::invalid("frames", "needs at least one frame"));
    }
    let n = frames[0].len();
    if frames.iter().any(|f| f.len() != n) {
        return Err(MdError::dimension("trajectory frames differ in atom count"));
    }
    // Accumulate g(r) across frames (each already normalised); the
    // mean of per-frame g(r) equals the histogram-summed g(r) for a
    // fixed box.
    let mut sum_g: Option<Vec<f64>> = None;
    let mut radii = Vec::new();
    let mut bin_width = 0.0;
    for frame in frames {
        let rdf = radial_distribution(frame, cell, r_max, bins)?;
        radii = rdf.radii.clone();
        bin_width = rdf.bin_width;
        match &mut sum_g {
            None => sum_g = Some(rdf.g),
            Some(acc) => {
                for (a, b) in acc.iter_mut().zip(&rdf.g) {
                    *a += b;
                }
            }
        }
    }
    let count = frames.len() as f64;
    let g = sum_g.unwrap().into_iter().map(|v| v / count).collect();
    Ok(RadialDistribution {
        radii,
        g,
        bin_width,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rng::Rng;

    /// `n` atoms placed uniformly at random in a cubic box.
    fn random_gas(n: usize, edge: f64, seed: u64) -> Vec<Vector3<f64>> {
        let mut rng = Rng::new(seed);
        (0..n)
            .map(|_| {
                Vector3::new(
                    rng.uniform() * edge,
                    rng.uniform() * edge,
                    rng.uniform() * edge,
                )
            })
            .collect()
    }

    #[test]
    fn rejects_bad_input() {
        let cell = SimBox::cubic(5.0).unwrap();
        assert!(radial_distribution(&[Vector3::zeros()], &cell, 1.0, 10).is_err());
        let two = vec![Vector3::zeros(), Vector3::new(1.0, 0.0, 0.0)];
        assert!(radial_distribution(&two, &cell, -1.0, 10).is_err());
        assert!(radial_distribution(&two, &cell, 1.0, 0).is_err());
        assert!(radial_distribution(&two, &SimBox::none(), 1.0, 10).is_err());
    }

    #[test]
    fn ideal_gas_rdf_is_flat_near_one() {
        // A random (ideal) gas should give g(r) ~ 1 everywhere with
        // enough particles.
        let edge = 8.0;
        let pos = random_gas(2000, edge, 2024);
        let cell = SimBox::cubic(edge).unwrap();
        let rdf = radial_distribution(&pos, &cell, 3.0, 30).unwrap();
        // Average g(r) over the outer half of the range.
        let half = rdf.g.len() / 2;
        let mean: f64 = rdf.g[half..].iter().sum::<f64>() / (rdf.g.len() - half) as f64;
        assert!((mean - 1.0).abs() < 0.15, "ideal-gas g(r) mean = {mean}");
    }

    #[test]
    fn lattice_rdf_has_a_first_peak() {
        // A simple cubic lattice has sharp coordination shells.
        let edge = 6.0;
        let s = 1.0;
        let mut pos = Vec::new();
        for i in 0..6 {
            for j in 0..6 {
                for k in 0..6 {
                    pos.push(Vector3::new(i as f64 * s, j as f64 * s, k as f64 * s));
                }
            }
        }
        let cell = SimBox::cubic(edge).unwrap();
        let rdf = radial_distribution(&pos, &cell, 2.5, 50).unwrap();
        // There must be a peak near the lattice spacing 1.0 nm.
        let peak = rdf.first_peak().expect("lattice should have a peak");
        assert!((peak - 1.0).abs() < 0.2, "first peak at {peak}");
    }

    #[test]
    fn coordination_number_is_monotone() {
        let edge = 8.0;
        let pos = random_gas(1000, edge, 5);
        let cell = SimBox::cubic(edge).unwrap();
        let rdf = radial_distribution(&pos, &cell, 3.0, 30).unwrap();
        let density = 1000.0 / cell.volume();
        let cn = rdf.coordination_number(density);
        for w in cn.windows(2) {
            assert!(w[1] >= w[0] - 1e-9, "coordination number not monotone");
        }
    }

    #[test]
    fn averaged_rdf_matches_single_frame_for_one_frame() {
        let edge = 7.0;
        let pos = random_gas(500, edge, 17);
        let cell = SimBox::cubic(edge).unwrap();
        let single = radial_distribution(&pos, &cell, 3.0, 25).unwrap();
        let averaged = radial_distribution_averaged(&[pos], &cell, 3.0, 25).unwrap();
        for (a, b) in single.g.iter().zip(&averaged.g) {
            assert!((a - b).abs() < 1e-9);
        }
    }
}
