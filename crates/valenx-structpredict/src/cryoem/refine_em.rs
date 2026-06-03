//! **Feature 28 — projection-matching iterative refinement.**
//!
//! Reconstruction needs each particle's orientation — but the
//! orientations are unknown. The classical solution is **projection
//! matching**, the iterative refinement engine of every classical
//! cryo-EM package (EMAN, SPIDER, FREALIGN):
//!
//! 1. Start from an initial 3-D map (even a crude one).
//! 2. **Reproject** the current map over a grid of orientations.
//! 3. For every particle, find the reprojection it correlates with
//!    best — that orientation is the particle's new assigned angle.
//! 4. **Reconstruct** a new map from the particles at their new
//!    orientations.
//! 5. Repeat. Each cycle sharpens both the orientations and the map;
//!    the process converges (the orientation assignments stop
//!    changing).
//!
//! This module runs that loop, reusing
//! [`crate::cryoem::reconstruct`]'s forward projection and
//! back-projection. It is the genuine projection-matching algorithm.
//! A production refinement adds CTF correction, a finer angular grid
//! with local refinement, and the regularised-likelihood weighting of
//! RELION; this v1 is the textbook projection-matching core.

use serde::{Deserialize, Serialize};

use crate::cryoem::mrc::{ParticleStack, Volume3d};
use crate::cryoem::reconstruct::{project_volume, reconstruct_3d, Projection};
use crate::error::{Result, StructPredictError};

/// One particle's refined orientation.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ParticleOrientation {
    /// `rot` Euler angle, radians.
    pub rot: f64,
    /// `tilt` Euler angle, radians.
    pub tilt: f64,
    /// `psi` Euler angle, radians.
    pub psi: f64,
    /// The correlation score of the assigned orientation.
    pub score: f64,
}

/// The outcome of a projection-matching refinement run.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct RefinementResult {
    /// The refined 3-D density map.
    pub volume: Volume3d,
    /// The final assigned orientation per particle.
    pub orientations: Vec<ParticleOrientation>,
    /// Mean correlation of the particles to their assigned
    /// reprojections in the final round.
    pub mean_score: f64,
    /// Refinement cycles performed.
    pub iterations: usize,
}

/// A grid of candidate orientations sampled over the sphere.
fn orientation_grid(n_tilt: usize, n_rot: usize) -> Vec<(f64, f64)> {
    let mut grid = Vec::new();
    for ti in 0..n_tilt.max(1) {
        // tilt 0..π
        let tilt = std::f64::consts::PI * ti as f64 / n_tilt.max(1) as f64;
        // The number of rot samples can be scaled by sin(tilt) for a
        // more uniform sphere sampling; a fixed count is simpler and
        // adequate for a v1.
        for ri in 0..n_rot.max(1) {
            let rot = 2.0 * std::f64::consts::PI * ri as f64 / n_rot.max(1) as f64;
            grid.push((rot, tilt));
        }
    }
    grid
}

/// Zero-normalised cross-correlation of two equal-size volumes'
/// projection images.
fn image_zncc(a: &crate::cryoem::mrc::Image2d, b: &crate::cryoem::mrc::Image2d) -> f64 {
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

/// Runs projection-matching iterative refinement.
///
/// `stack` is the particle stack; `initial` is the starting 3-D map
/// (its box must equal the stack's box size). Each cycle reprojects
/// `initial` over an `n_tilt × n_rot` orientation grid, assigns every
/// particle to its best-matching reprojection, and reconstructs a new
/// map; this repeats for up to `max_iterations` cycles or until the
/// orientation assignments stabilise.
///
/// # Errors
/// [`StructPredictError::Invalid`] for an empty stack, a box-size
/// mismatch, a degenerate grid, or `max_iterations == 0`.
pub fn projection_matching(
    stack: &ParticleStack,
    initial: &Volume3d,
    n_tilt: usize,
    n_rot: usize,
    max_iterations: usize,
) -> Result<RefinementResult> {
    if stack.is_empty() {
        return Err(StructPredictError::invalid("stack", "no particles"));
    }
    if max_iterations == 0 {
        return Err(StructPredictError::invalid(
            "max_iterations",
            "must be at least 1",
        ));
    }
    let box_size = stack.box_size;
    if initial.nx != box_size || initial.ny != box_size || initial.nz != box_size {
        return Err(StructPredictError::invalid(
            "initial",
            format!(
                "initial map is {}×{}×{}, not the stack's {box_size}³",
                initial.nx, initial.ny, initial.nz
            ),
        ));
    }
    if n_tilt == 0 || n_rot == 0 {
        return Err(StructPredictError::invalid(
            "grid",
            "n_tilt and n_rot must both be at least 1",
        ));
    }

    let grid = orientation_grid(n_tilt, n_rot);
    let mut current = initial.clone();
    let mut orientations: Vec<ParticleOrientation> = vec![
        ParticleOrientation {
            rot: 0.0,
            tilt: 0.0,
            psi: 0.0,
            score: f64::NEG_INFINITY,
        };
        stack.len()
    ];
    let mut iterations = 0;
    let mut mean_score = 0.0;

    for _ in 0..max_iterations {
        iterations += 1;
        // Reproject the current map over the grid.
        let mut reprojections = Vec::with_capacity(grid.len());
        for &(rot, tilt) in &grid {
            let img = project_volume(&current, rot, tilt, 0.0)?;
            reprojections.push((rot, tilt, img));
        }
        // Assign each particle to its best reprojection.
        let mut changed = false;
        let mut score_sum = 0.0;
        for (p, particle) in stack.particles.iter().enumerate() {
            let mut best = orientations[p];
            best.score = f64::NEG_INFINITY;
            for (rot, tilt, reproj) in &reprojections {
                let score = image_zncc(particle, reproj);
                if score > best.score {
                    best = ParticleOrientation {
                        rot: *rot,
                        tilt: *tilt,
                        psi: 0.0,
                        score,
                    };
                }
            }
            if (best.rot - orientations[p].rot).abs() > 1e-9
                || (best.tilt - orientations[p].tilt).abs() > 1e-9
            {
                changed = true;
            }
            orientations[p] = best;
            score_sum += best.score;
        }
        mean_score = score_sum / stack.len() as f64;

        // Reconstruct a new map from the freshly-assigned orientations.
        let projections: Vec<Projection> = stack
            .particles
            .iter()
            .zip(&orientations)
            .map(|(img, o)| Projection::new(img.clone(), o.rot, o.tilt, o.psi))
            .collect();
        let recon = reconstruct_3d(&projections, box_size)?;
        current = recon.volume;

        if !changed {
            break; // orientations converged
        }
    }

    Ok(RefinementResult {
        volume: current,
        orientations,
        mean_score,
        iterations,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cryoem::mrc::Image2d;

    /// An asymmetric volume — a brick offset from the centre, so that
    /// different views are genuinely distinguishable.
    fn asymmetric_volume(n: usize) -> Volume3d {
        let mut v = Volume3d::zeros_cube(n);
        let c = n / 2;
        for z in 0..n {
            for y in 0..n {
                for x in 0..n {
                    let inside = x > c && x < c + n / 3
                        && y > c.saturating_sub(n / 4)
                        && y < c + 1
                        && z > c.saturating_sub(1)
                        && z < c + 2;
                    if inside {
                        v.data[(z * n + y) * n + x] = 1.0;
                    }
                }
            }
        }
        v
    }

    #[test]
    fn refinement_runs_and_returns_a_map() {
        let n = 12;
        let truth = asymmetric_volume(n);
        // Build a particle stack: reprojections of the true volume at
        // known angles (the particles).
        let mut stack = ParticleStack::new(n);
        for ti in 0..3 {
            for ri in 0..4 {
                let tilt = std::f64::consts::PI * ti as f64 / 3.0;
                let rot = 2.0 * std::f64::consts::PI * ri as f64 / 4.0;
                let img = project_volume(&truth, rot, tilt, 0.0).expect("project");
                stack.push(img).expect("push");
            }
        }
        // Start from a featureless blob.
        let mut initial = Volume3d::zeros_cube(n);
        let c = n / 2;
        for z in 0..n {
            for y in 0..n {
                for x in 0..n {
                    let dx = x as i64 - c as i64;
                    let dy = y as i64 - c as i64;
                    let dz = z as i64 - c as i64;
                    if dx * dx + dy * dy + dz * dz < 9 {
                        initial.data[(z * n + y) * n + x] = 0.5;
                    }
                }
            }
        }
        let res = projection_matching(&stack, &initial, 3, 4, 3).expect("refine");
        assert_eq!(res.orientations.len(), stack.len());
        assert_eq!(res.volume.nx, n);
        // Every particle got a finite assigned score.
        assert!(res.orientations.iter().all(|o| o.score.is_finite()));
    }

    #[test]
    fn box_size_mismatch_rejected() {
        let mut stack = ParticleStack::new(8);
        stack.push(Image2d::zeros(8, 8)).expect("push");
        let wrong = Volume3d::zeros_cube(12);
        assert!(projection_matching(&stack, &wrong, 2, 2, 2).is_err());
    }

    #[test]
    fn empty_stack_rejected() {
        let empty = ParticleStack::new(8);
        let v = Volume3d::zeros_cube(8);
        assert!(projection_matching(&empty, &v, 2, 2, 2).is_err());
    }

    #[test]
    fn degenerate_grid_rejected() {
        let mut stack = ParticleStack::new(8);
        stack.push(Image2d::zeros(8, 8)).expect("push");
        let v = Volume3d::zeros_cube(8);
        assert!(projection_matching(&stack, &v, 0, 2, 2).is_err());
    }
}
