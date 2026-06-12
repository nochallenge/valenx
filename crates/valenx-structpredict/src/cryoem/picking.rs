//! **Feature 25 — particle picking (template / blob-based).**
//!
//! Before any reconstruction, the individual particles must be
//! *located* in each large micrograph. The classical pickers — the
//! ones that predate the neural pickers (Topaz, crYOLO, which are
//! adapter-only here) — are:
//!
//! - **Template matching** — cross-correlate the micrograph with a
//!   reference template (a projection, or a previous class average);
//!   the correlation peaks are particles. This is the
//!   `signature`/`gautomatch`-class picker.
//! - **Blob picking** — when no template exists, correlate against a
//!   simple disk / Gaussian "blob" of roughly the particle's size;
//!   the peaks are blob-like dense regions. EMAN2's
//!   `e2boxer --gauss`, RELION's LoG (Laplacian-of-Gaussian) picker.
//!
//! Both reduce to: compute a correlation map, then extract its
//! significant local maxima with a minimum-separation constraint so
//! the same particle is not picked twice. This module implements
//! that — a real classical picker, with the honest caveat that on
//! low-contrast data a neural picker does substantially better.

use serde::{Deserialize, Serialize};

use crate::cryoem::mrc::Image2d;
use crate::error::{Result, StructPredictError};

/// One picked particle: its centre in the micrograph and the
/// correlation score of the pick.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ParticlePick {
    /// Particle-centre x coordinate (pixels).
    pub x: usize,
    /// Particle-centre y coordinate (pixels).
    pub y: usize,
    /// Normalised correlation score of the pick (higher = stronger).
    pub score: f64,
}

/// Builds a circular-disk blob template of the given diameter, for
/// template-free blob picking.
///
/// The blob is a soft-edged disk: 1.0 in the centre, falling smoothly
/// to 0 at the radius. Used as the matching template when no real
/// reference is available.
pub fn blob_template(diameter: usize) -> Image2d {
    let size = diameter.max(3);
    let mut img = Image2d::zeros(size, size);
    let r = diameter as f64 / 2.0;
    let cx = (size as f64 - 1.0) / 2.0;
    let cy = (size as f64 - 1.0) / 2.0;
    for y in 0..size {
        for x in 0..size {
            let d = ((x as f64 - cx).powi(2) + (y as f64 - cy).powi(2)).sqrt();
            // A soft cosine edge over the outer 20 % of the radius.
            let v = if d <= r * 0.8 {
                1.0
            } else if d <= r {
                let t = (d - r * 0.8) / (r * 0.2);
                0.5 * (1.0 + (std::f64::consts::PI * t).cos())
            } else {
                0.0
            };
            img.data[y * size + x] = v as f32;
        }
    }
    img
}

/// Computes the normalised cross-correlation map of a micrograph with
/// a template.
///
/// At every position the template is overlaid on the micrograph and
/// the **zero-normalised cross-correlation** (ZNCC) is computed — the
/// correlation is invariant to the local micrograph brightness and
/// contrast, which is essential because cryo-EM ice thickness varies
/// across a micrograph. The returned map is the same size as the
/// micrograph (positions where the template would overrun the edge
/// are 0).
fn correlation_map(micrograph: &Image2d, template: &Image2d) -> Vec<f64> {
    let (mw, mh) = (micrograph.width, micrograph.height);
    let (tw, th) = (template.width, template.height);
    let mut map = vec![0.0f64; mw * mh];
    if tw > mw || th > mh {
        return map;
    }
    // Zero-mean the template once.
    let t_mean = template.mean();
    let t_dev: Vec<f64> = template.data.iter().map(|&v| v as f64 - t_mean).collect();
    let t_norm: f64 = t_dev.iter().map(|d| d * d).sum::<f64>().sqrt();
    if t_norm < 1e-12 {
        return map;
    }
    let half_w = tw / 2;
    let half_h = th / 2;
    for cy in half_h..(mh - (th - half_h - 1)) {
        for cx in half_w..(mw - (tw - half_w - 1)) {
            // Window mean.
            let mut wsum = 0.0;
            for ty in 0..th {
                for tx in 0..tw {
                    let mx = cx + tx - half_w;
                    let my = cy + ty - half_h;
                    wsum += micrograph.data[my * mw + mx] as f64;
                }
            }
            let w_mean = wsum / (tw * th) as f64;
            // Covariance and window norm.
            let mut cov = 0.0;
            let mut w_norm = 0.0;
            for ty in 0..th {
                for tx in 0..tw {
                    let mx = cx + tx - half_w;
                    let my = cy + ty - half_h;
                    let wd = micrograph.data[my * mw + mx] as f64 - w_mean;
                    cov += wd * t_dev[ty * tw + tx];
                    w_norm += wd * wd;
                }
            }
            let denom = w_norm.sqrt() * t_norm;
            map[cy * mw + cx] = if denom > 1e-12 { cov / denom } else { 0.0 };
        }
    }
    map
}

/// Extracts the significant local maxima of a correlation map as
/// particle picks, with a minimum-separation constraint.
fn extract_peaks(
    map: &[f64],
    width: usize,
    height: usize,
    threshold: f64,
    min_separation: usize,
) -> Vec<ParticlePick> {
    // Collect all candidate maxima above the threshold.
    let mut candidates: Vec<ParticlePick> = Vec::new();
    for y in 1..height.saturating_sub(1) {
        for x in 1..width.saturating_sub(1) {
            let v = map[y * width + x];
            if v < threshold {
                continue;
            }
            // 8-neighbour local-maximum test.
            let mut is_max = true;
            'nb: for dy in -1i64..=1 {
                for dx in -1i64..=1 {
                    if dx == 0 && dy == 0 {
                        continue;
                    }
                    let nx = (x as i64 + dx) as usize;
                    let ny = (y as i64 + dy) as usize;
                    if map[ny * width + nx] > v {
                        is_max = false;
                        break 'nb;
                    }
                }
            }
            if is_max {
                candidates.push(ParticlePick { x, y, score: v });
            }
        }
    }
    // Greedy non-maximum suppression: take the strongest, drop any
    // candidate within `min_separation` of an already-accepted pick.
    candidates.sort_by(|a, b| {
        b.score
            .partial_cmp(&a.score)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    let mut picks: Vec<ParticlePick> = Vec::new();
    let sep2 = (min_separation * min_separation) as i64;
    for cand in candidates {
        let too_close = picks.iter().any(|p| {
            let dx = p.x as i64 - cand.x as i64;
            let dy = p.y as i64 - cand.y as i64;
            dx * dx + dy * dy < sep2
        });
        if !too_close {
            picks.push(cand);
        }
    }
    picks
}

/// Picks particles from a micrograph by template matching.
///
/// Correlates `micrograph` with `template`, then extracts the
/// correlation peaks above `score_threshold` separated by at least
/// `min_separation` pixels.
///
/// # Errors
/// [`StructPredictError::Invalid`] for an empty micrograph / template
/// or a template larger than the micrograph.
pub fn pick_particles(
    micrograph: &Image2d,
    template: &Image2d,
    score_threshold: f64,
    min_separation: usize,
) -> Result<Vec<ParticlePick>> {
    if micrograph.is_empty() {
        return Err(StructPredictError::invalid("micrograph", "empty image"));
    }
    if template.is_empty() {
        return Err(StructPredictError::invalid("template", "empty template"));
    }
    if template.width > micrograph.width || template.height > micrograph.height {
        return Err(StructPredictError::invalid(
            "template",
            "template is larger than the micrograph",
        ));
    }
    let map = correlation_map(micrograph, template);
    Ok(extract_peaks(
        &map,
        micrograph.width,
        micrograph.height,
        score_threshold,
        min_separation.max(1),
    ))
}

/// Picks particles from a micrograph by template-free blob matching.
///
/// Builds a soft-disk blob of `particle_diameter` pixels and runs
/// [`pick_particles`] with it.
///
/// # Errors
/// [`StructPredictError::Invalid`] as for [`pick_particles`].
pub fn pick_blobs(
    micrograph: &Image2d,
    particle_diameter: usize,
    score_threshold: f64,
    min_separation: usize,
) -> Result<Vec<ParticlePick>> {
    let template = blob_template(particle_diameter);
    pick_particles(micrograph, &template, score_threshold, min_separation)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A micrograph with bright disk-shaped "particles" planted at
    /// known centres on a noisy background.
    fn synthetic_micrograph(size: usize, particles: &[(usize, usize)], radius: f64) -> Image2d {
        let mut img = Image2d::zeros(size, size);
        // A faint deterministic background ripple.
        for y in 0..size {
            for x in 0..size {
                img.data[y * size + x] = 0.05 * ((x as f32 * 0.3).sin() + (y as f32 * 0.21).cos());
            }
        }
        for &(px, py) in particles {
            for y in 0..size {
                for x in 0..size {
                    let d =
                        ((x as f64 - px as f64).powi(2) + (y as f64 - py as f64).powi(2)).sqrt();
                    if d <= radius {
                        img.data[y * size + x] += 1.0;
                    }
                }
            }
        }
        img
    }

    #[test]
    fn blob_template_is_a_soft_disk() {
        let t = blob_template(10);
        // Centre is bright, corners are zero.
        assert!(t.at(5, 5).unwrap() > 0.9);
        assert_eq!(t.at(0, 0).unwrap(), 0.0);
    }

    #[test]
    fn picks_planted_particles() {
        let centres = [(15usize, 15usize), (40, 40), (15, 45)];
        let micro = synthetic_micrograph(64, &centres, 5.0);
        let picks = pick_blobs(&micro, 11, 0.3, 12).expect("pick");
        assert!(
            picks.len() >= centres.len(),
            "found {} picks for {} particles",
            picks.len(),
            centres.len()
        );
        // Every planted centre is near some pick.
        for &(cx, cy) in &centres {
            let near = picks.iter().any(|p| {
                let dx = p.x as i64 - cx as i64;
                let dy = p.y as i64 - cy as i64;
                dx * dx + dy * dy <= 9
            });
            assert!(near, "particle at ({cx},{cy}) was picked");
        }
    }

    #[test]
    fn minimum_separation_prevents_double_picks() {
        // One particle — the minimum-separation NMS must not pick it
        // many times. The score threshold is 0.85: the picker uses
        // *normalised* cross-correlation, which is contrast-invariant,
        // so the deterministic background ripple in the synthetic
        // micrograph also produces shape-correlated peaks (NCC up to
        // ~0.83). A realistic NCC threshold sits above that band — the
        // genuine particle scores ~0.91 — so 0.85 isolates the real
        // particle and the test then exercises exactly the NMS.
        let micro = synthetic_micrograph(48, &[(24, 24)], 5.0);
        let picks = pick_blobs(&micro, 11, 0.85, 15).expect("pick");
        assert_eq!(picks.len(), 1, "single particle picked once");
    }

    #[test]
    fn empty_micrograph_rejected() {
        let empty = Image2d::zeros(0, 0);
        let t = blob_template(8);
        assert!(pick_particles(&empty, &t, 0.5, 4).is_err());
    }

    #[test]
    fn oversized_template_rejected() {
        let micro = Image2d::zeros(10, 10);
        let big = blob_template(40);
        assert!(pick_particles(&micro, &big, 0.5, 4).is_err());
    }
}
