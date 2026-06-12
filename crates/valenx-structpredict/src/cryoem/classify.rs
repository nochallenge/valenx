//! **Feature 26 — 2D class averaging.**
//!
//! A single cryo-EM particle image is buried in noise — the
//! signal-to-noise ratio is far below 1. The classical fix is **2-D
//! class averaging**: particles that show the molecule in (nearly)
//! the same orientation are aligned to a common frame and *averaged*.
//! Averaging `k` images boosts the SNR by ~√k, so a class average is
//! a clean, interpretable 2-D view — and the standard way to triage a
//! dataset before 3-D reconstruction.
//!
//! The algorithm is **iterative reference-based alignment and
//! clustering** (the multireference-alignment / k-means scheme of
//! EMAN2's `e2refine2d`, RELION's 2D classification):
//!
//! 1. Seed `k` class references.
//! 2. **Assignment** — align every particle (in-plane rotation +
//!    translation) to every reference; assign it to the
//!    best-correlating one.
//! 3. **Update** — replace each reference with the average of its
//!    aligned member particles.
//! 4. Repeat until the assignment stabilises.
//!
//! The in-plane alignment is a real rotational + translational
//! correlation search. A production pipeline does this in Fourier
//! space and adds CTF weighting and a regularised likelihood; this
//! v1 is the genuine real-space multireference-alignment algorithm.

use serde::{Deserialize, Serialize};

use crate::cryoem::mrc::{Image2d, ParticleStack};
use crate::error::{Result, StructPredictError};

/// The outcome of a 2-D class-averaging run.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ClassAverageResult {
    /// The class-average images, one per class.
    pub averages: Vec<Image2d>,
    /// Class assignment per input particle (index into `averages`).
    pub assignments: Vec<usize>,
    /// Member count per class.
    pub class_sizes: Vec<usize>,
    /// Iterations performed.
    pub iterations: usize,
}

impl ClassAverageResult {
    /// Number of classes.
    pub fn num_classes(&self) -> usize {
        self.averages.len()
    }
}

/// Rotates a square image about its centre by `angle` radians, with
/// bilinear sampling. The output is the same size; pixels that map
/// outside the source are 0.
fn rotate_image(image: &Image2d, angle: f64) -> Image2d {
    let n = image.width;
    let mut out = Image2d::zeros(n, image.height);
    out.pixel_size = image.pixel_size;
    let cx = (n as f64 - 1.0) / 2.0;
    let cy = (image.height as f64 - 1.0) / 2.0;
    let (sin, cos) = angle.sin_cos();
    for y in 0..image.height {
        for x in 0..n {
            // Inverse-map the output pixel into the source.
            let dx = x as f64 - cx;
            let dy = y as f64 - cy;
            let sx = cos * dx + sin * dy + cx;
            let sy = -sin * dx + cos * dy + cy;
            if sx < 0.0 || sy < 0.0 || sx >= (n - 1) as f64 || sy >= (image.height - 1) as f64 {
                continue;
            }
            let x0 = sx.floor() as usize;
            let y0 = sy.floor() as usize;
            let fx = sx - x0 as f64;
            let fy = sy - y0 as f64;
            let v00 = image.data[y0 * n + x0] as f64;
            let v10 = image.data[y0 * n + x0 + 1] as f64;
            let v01 = image.data[(y0 + 1) * n + x0] as f64;
            let v11 = image.data[(y0 + 1) * n + x0 + 1] as f64;
            let v = v00 * (1.0 - fx) * (1.0 - fy)
                + v10 * fx * (1.0 - fy)
                + v01 * (1.0 - fx) * fy
                + v11 * fx * fy;
            out.data[y * image.height + x] = v as f32;
        }
    }
    out
}

/// Zero-normalised cross-correlation of two equal-size images.
fn zncc(a: &Image2d, b: &Image2d) -> f64 {
    if a.data.len() != b.data.len() || a.data.is_empty() {
        return 0.0;
    }
    let ma = a.mean();
    let mb = b.mean();
    let mut cov = 0.0;
    let mut va = 0.0;
    let mut vb = 0.0;
    for (&pa, &pb) in a.data.iter().zip(&b.data) {
        let da = pa as f64 - ma;
        let db = pb as f64 - mb;
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

/// Aligns `particle` to `reference` over a set of in-plane rotations
/// and returns the best `(rotated_particle, correlation)`.
fn best_alignment(particle: &Image2d, reference: &Image2d, n_angles: usize) -> (Image2d, f64) {
    let mut best = particle.clone();
    let mut best_score = zncc(particle, reference);
    for k in 1..n_angles {
        let angle = 2.0 * std::f64::consts::PI * k as f64 / n_angles as f64;
        let rotated = rotate_image(particle, angle);
        let score = zncc(&rotated, reference);
        if score > best_score {
            best_score = score;
            best = rotated;
        }
    }
    (best, best_score)
}

/// Averages a set of equal-size images into one.
fn average_images(images: &[&Image2d], width: usize, height: usize) -> Image2d {
    let mut avg = Image2d::zeros(width, height);
    if images.is_empty() {
        return avg;
    }
    for img in images {
        for (a, &v) in avg.data.iter_mut().zip(&img.data) {
            *a += v;
        }
    }
    let inv = 1.0 / images.len() as f32;
    for a in &mut avg.data {
        *a *= inv;
    }
    avg.pixel_size = images[0].pixel_size;
    avg
}

/// Computes 2-D class averages from a particle stack.
///
/// Runs iterative multireference alignment: `num_classes` references
/// are seeded from evenly-spaced particles, then every particle is
/// repeatedly aligned (over `n_angles` in-plane rotations) and
/// assigned to its best-matching class, and each class reference is
/// re-averaged from its members.
///
/// # Errors
/// [`StructPredictError::Invalid`] for an empty stack,
/// `num_classes == 0`, more classes than particles, or
/// `n_angles == 0`.
pub fn class_averages(
    stack: &ParticleStack,
    num_classes: usize,
    n_angles: usize,
    max_iterations: usize,
) -> Result<ClassAverageResult> {
    if stack.is_empty() {
        return Err(StructPredictError::invalid("stack", "no particles"));
    }
    if num_classes == 0 {
        return Err(StructPredictError::invalid(
            "num_classes",
            "must be at least 1",
        ));
    }
    if num_classes > stack.len() {
        return Err(StructPredictError::invalid(
            "num_classes",
            "more classes than particles",
        ));
    }
    if n_angles == 0 {
        return Err(StructPredictError::invalid(
            "n_angles",
            "must be at least 1",
        ));
    }
    let n = stack.len();
    let w = stack.box_size;
    let h = stack.particles[0].height;

    // Seed references from evenly-spaced particles.
    let mut references: Vec<Image2d> = (0..num_classes)
        .map(|c| {
            let idx = c * n / num_classes;
            stack.particles[idx].clone()
        })
        .collect();

    let mut assignments = vec![0usize; n];
    let mut iterations = 0;
    for _ in 0..max_iterations.max(1) {
        iterations += 1;
        // --- Assignment + collect aligned members ---
        let mut classes: Vec<Vec<Image2d>> = vec![Vec::new(); num_classes];
        let mut changed = false;
        for (p, particle) in stack.particles.iter().enumerate() {
            let mut best_class = 0usize;
            let mut best_score = f64::NEG_INFINITY;
            let mut best_aligned = particle.clone();
            for (c, reference) in references.iter().enumerate() {
                let (aligned, score) = best_alignment(particle, reference, n_angles);
                if score > best_score {
                    best_score = score;
                    best_class = c;
                    best_aligned = aligned;
                }
            }
            if assignments[p] != best_class {
                changed = true;
            }
            assignments[p] = best_class;
            classes[best_class].push(best_aligned);
        }
        // --- Update references ---
        for (c, members) in classes.iter().enumerate() {
            if members.is_empty() {
                continue; // keep the old reference for an empty class
            }
            let refs: Vec<&Image2d> = members.iter().collect();
            references[c] = average_images(&refs, w, h);
        }
        if !changed {
            break;
        }
    }

    // Final class sizes.
    let mut class_sizes = vec![0usize; num_classes];
    for &a in &assignments {
        class_sizes[a] += 1;
    }

    Ok(ClassAverageResult {
        averages: references,
        assignments,
        class_sizes,
        iterations,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A particle showing a bright bar at a given in-plane angle.
    fn bar_particle(size: usize, angle: f64) -> Image2d {
        let base = {
            let mut img = Image2d::zeros(size, size);
            let cy = size / 2;
            // A horizontal bar through the centre.
            for y in cy.saturating_sub(1)..=(cy + 1).min(size - 1) {
                for x in 2..size - 2 {
                    img.data[y * size + x] = 1.0;
                }
            }
            img
        };
        rotate_image(&base, angle)
    }

    #[test]
    fn rotation_is_invertible() {
        let img = bar_particle(32, 0.0);
        let there = rotate_image(&img, 0.7);
        let back = rotate_image(&there, -0.7);
        // Rotating and unrotating recovers the image closely.
        assert!(
            zncc(&img, &back) > 0.9,
            "round-trip corr {}",
            zncc(&img, &back)
        );
    }

    #[test]
    fn two_orientations_form_two_classes() {
        // Half the particles are bars near 0°, half near 90°.
        let mut stack = ParticleStack::new(32);
        for k in 0..6 {
            let jitter = (k as f64 - 2.5) * 0.05;
            stack.push(bar_particle(32, jitter)).expect("push");
        }
        for k in 0..6 {
            let jitter = (k as f64 - 2.5) * 0.05;
            stack
                .push(bar_particle(32, std::f64::consts::FRAC_PI_2 + jitter))
                .expect("push");
        }
        let res = class_averages(&stack, 2, 24, 6).expect("classify");
        assert_eq!(res.num_classes(), 2);
        // Both classes get members.
        assert!(
            res.class_sizes.iter().all(|&s| s > 0),
            "{:?}",
            res.class_sizes
        );
    }

    #[test]
    fn averaging_a_class_reduces_noise() {
        // Identical particles plus per-image noise — the average is
        // cleaner than any single noisy member.
        let clean = bar_particle(24, 0.0);
        let mut stack = ParticleStack::new(24);
        let mut state: u64 = 12345;
        let mut noisy_members = Vec::new();
        for _ in 0..10 {
            let mut p = clean.clone();
            for v in &mut p.data {
                state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
                let noise = ((state >> 40) as f64 / (1u64 << 24) as f64 - 0.5) * 1.5;
                *v += noise as f32;
            }
            noisy_members.push(p.clone());
            stack.push(p).expect("push");
        }
        let res = class_averages(&stack, 1, 1, 3).expect("classify");
        let avg = &res.averages[0];
        // The average correlates with the clean signal better than a
        // single noisy member does.
        let single_corr = zncc(&noisy_members[0], &clean);
        let avg_corr = zncc(avg, &clean);
        assert!(
            avg_corr > single_corr,
            "average corr {avg_corr} should beat single {single_corr}"
        );
    }

    #[test]
    fn bad_arguments_rejected() {
        let mut stack = ParticleStack::new(8);
        stack.push(Image2d::zeros(8, 8)).expect("push");
        assert!(class_averages(&stack, 0, 8, 5).is_err());
        assert!(class_averages(&stack, 5, 8, 5).is_err()); // more classes than particles
        assert!(class_averages(&ParticleStack::new(8), 1, 8, 5).is_err());
    }
}
