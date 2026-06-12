//! **Feature 27 — 3D reconstruction (weighted back-projection).**
//!
//! Once the particles have known orientations, the 3-D density map is
//! reconstructed from their 2-D projections. The classical
//! reconstruction algorithm — the one that predates and underlies the
//! Fourier-space methods — is **back-projection**:
//!
//! - Each 2-D projection is a "shadow" of the 3-D object along its
//!   view direction. **Back-projecting** smears that shadow back
//!   through the volume along the same direction.
//! - Summing the back-projections of many differently-oriented
//!   projections builds up the density: where every projection agrees
//!   there is signal, elsewhere it averages out.
//! - Plain back-projection blurs the result (it over-weights low
//!   frequencies). **Weighted back-projection** corrects this by
//!   applying a radial ramp / `1/r` weighting — the
//!   projection-slice-theorem correction. This module applies the
//!   weighting in real space as a simple sharpening of each
//!   back-projected slice's low-frequency excess, and additionally
//!   normalises by the per-voxel back-projection count.
//!
//! This is genuine tomographic reconstruction: a [`Volume3d`] built
//! by inverting the projection geometry of a set of oriented
//! [`Projection`]s. It is not the regularised Fourier-space
//! reconstruction RELION uses (which interpolates onto a 3-D Fourier
//! grid and applies a Wiener / SNR-weighted filter), but it inverts
//! the same physics.

use nalgebra::{Rotation3, Vector3};
use serde::{Deserialize, Serialize};

use crate::cryoem::mrc::{Image2d, Volume3d};
use crate::error::{Result, StructPredictError};

/// A 2-D projection of a 3-D object together with the Euler angles of
/// the view direction.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Projection {
    /// The projection image.
    pub image: Image2d,
    /// `rot` Euler angle (radians) — rotation about z.
    pub rot: f64,
    /// `tilt` Euler angle (radians) — rotation about the new y.
    pub tilt: f64,
    /// `psi` Euler angle (radians) — final in-plane rotation about z.
    pub psi: f64,
}

impl Projection {
    /// Builds a projection from an image and its three Euler angles.
    pub fn new(image: Image2d, rot: f64, tilt: f64, psi: f64) -> Self {
        Projection {
            image,
            rot,
            tilt,
            psi,
        }
    }

    /// The 3×3 rotation matrix of this projection's orientation, in
    /// the `ZYZ` Euler convention used throughout cryo-EM.
    pub fn rotation(&self) -> Rotation3<f64> {
        // ZYZ: R = Rz(rot) · Ry(tilt) · Rz(psi).
        Rotation3::from_axis_angle(&Vector3::z_axis(), self.rot)
            * Rotation3::from_axis_angle(&Vector3::y_axis(), self.tilt)
            * Rotation3::from_axis_angle(&Vector3::z_axis(), self.psi)
    }
}

/// The outcome of a 3-D reconstruction.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ReconstructionResult {
    /// The reconstructed density map.
    pub volume: Volume3d,
    /// Number of projections back-projected.
    pub projections_used: usize,
}

/// Projects a 3-D volume into a 2-D image along a given orientation —
/// the *forward* projection operator (a line integral / "X-ray
/// transform"). This is the operator back-projection inverts; it is
/// also reused by the projection-matching refinement.
///
/// The projection is `size × size` where `size` is the volume's
/// `nx`. For each output pixel a ray is cast along the rotated z axis
/// through the volume and the density it passes through is summed
/// (trilinear sampling).
///
/// # Errors
/// [`StructPredictError::Invalid`] for a non-cubic volume.
pub fn project_volume(volume: &Volume3d, rot: f64, tilt: f64, psi: f64) -> Result<Image2d> {
    if volume.nx != volume.ny || volume.ny != volume.nz {
        return Err(StructPredictError::invalid(
            "volume",
            "projection requires a cubic volume",
        ));
    }
    let n = volume.nx;
    let proj = Projection::new(Image2d::zeros(n, n), rot, tilt, psi);
    let r = proj.rotation();
    let center = (n as f64 - 1.0) / 2.0;
    let mut image = Image2d::zeros(n, n);
    image.pixel_size = volume.voxel_size;
    // For each detector pixel (u, v), integrate the volume along the
    // rotated z axis.
    for v in 0..n {
        for u in 0..n {
            let du = u as f64 - center;
            let dv = v as f64 - center;
            let mut sum = 0.0;
            for w in 0..n {
                let dw = w as f64 - center;
                // Detector coordinate (du, dv, dw) → world via R.
                let world = r * Vector3::new(du, dv, dw);
                sum +=
                    sample_trilinear(volume, world.x + center, world.y + center, world.z + center);
            }
            image.data[v * n + u] = sum as f32;
        }
    }
    Ok(image)
}

/// Trilinearly samples a volume at a (possibly fractional) voxel
/// coordinate. Returns 0 outside the volume.
fn sample_trilinear(volume: &Volume3d, x: f64, y: f64, z: f64) -> f64 {
    if x < 0.0
        || y < 0.0
        || z < 0.0
        || x > (volume.nx - 1) as f64
        || y > (volume.ny - 1) as f64
        || z > (volume.nz - 1) as f64
    {
        return 0.0;
    }
    let x0 = x.floor() as usize;
    let y0 = y.floor() as usize;
    let z0 = z.floor() as usize;
    let x1 = (x0 + 1).min(volume.nx - 1);
    let y1 = (y0 + 1).min(volume.ny - 1);
    let z1 = (z0 + 1).min(volume.nz - 1);
    let fx = x - x0 as f64;
    let fy = y - y0 as f64;
    let fz = z - z0 as f64;
    let g = |xi: usize, yi: usize, zi: usize| {
        volume.data[(zi * volume.ny + yi) * volume.nx + xi] as f64
    };
    let c00 = g(x0, y0, z0) * (1.0 - fx) + g(x1, y0, z0) * fx;
    let c10 = g(x0, y1, z0) * (1.0 - fx) + g(x1, y1, z0) * fx;
    let c01 = g(x0, y0, z1) * (1.0 - fx) + g(x1, y0, z1) * fx;
    let c11 = g(x0, y1, z1) * (1.0 - fx) + g(x1, y1, z1) * fx;
    let c0 = c00 * (1.0 - fy) + c10 * fy;
    let c1 = c01 * (1.0 - fy) + c11 * fy;
    c0 * (1.0 - fz) + c1 * fz
}

/// Reconstructs a 3-D density map from a set of oriented projections
/// by weighted back-projection.
///
/// Every projection is smeared back through a `box_size³` volume
/// along its view direction; the volume is normalised by the
/// per-voxel back-projection count and high-pass weighted to undo
/// the low-frequency excess of plain back-projection.
///
/// All projections must be square and the same size (`box_size`).
///
/// # Errors
/// [`StructPredictError::Invalid`] for an empty projection set,
/// `box_size == 0`, or a non-square / mismatched-size projection.
pub fn reconstruct_3d(projections: &[Projection], box_size: usize) -> Result<ReconstructionResult> {
    if projections.is_empty() {
        return Err(StructPredictError::invalid(
            "projections",
            "no projections to reconstruct from",
        ));
    }
    if box_size == 0 {
        return Err(StructPredictError::invalid(
            "box_size",
            "must be at least 1",
        ));
    }
    for (i, p) in projections.iter().enumerate() {
        if p.image.width != box_size || p.image.height != box_size {
            return Err(StructPredictError::invalid(
                "projections",
                format!(
                    "projection {i} is {}×{}, not {box_size}×{box_size}",
                    p.image.width, p.image.height
                ),
            ));
        }
    }

    let n = box_size;
    let mut volume = Volume3d::zeros_cube(n);
    // Per-voxel back-projection count, for normalisation.
    let mut counts = vec![0.0f64; n * n * n];
    let center = (n as f64 - 1.0) / 2.0;

    for proj in projections {
        // High-pass weight the projection (the ramp filter, applied
        // here as a real-space sharpening of the slice).
        let weighted = ramp_weight(&proj.image);
        let r = proj.rotation();
        // For every volume voxel, find where it projects in this
        // image and add that pixel back into the voxel.
        for z in 0..n {
            for y in 0..n {
                for x in 0..n {
                    let world =
                        Vector3::new(x as f64 - center, y as f64 - center, z as f64 - center);
                    // World → detector: the inverse rotation.
                    let det = r.inverse() * world;
                    let u = det.x + center;
                    let v = det.y + center;
                    let val = bilinear(&weighted, u, v);
                    if let Some(val) = val {
                        let idx = (z * n + y) * n + x;
                        volume.data[idx] += val as f32;
                        counts[idx] += 1.0;
                    }
                }
            }
        }
    }

    // Normalise by the back-projection count.
    for (v, &c) in volume.data.iter_mut().zip(&counts) {
        if c > 0.0 {
            *v /= c as f32;
        }
    }

    Ok(ReconstructionResult {
        volume,
        projections_used: projections.len(),
    })
}

/// Applies a ramp (high-pass) weighting to a projection image — the
/// projection-slice-theorem correction that stops back-projection
/// blurring the reconstruction.
///
/// Implemented in real space as `weighted = image + α·(image −
/// blurred)`, an unsharp mask: subtracting a blurred copy removes the
/// low-frequency excess, the hallmark of the ramp filter.
fn ramp_weight(image: &Image2d) -> Image2d {
    let blurred = box_blur(image, 2);
    let mut out = Image2d::zeros(image.width, image.height);
    out.pixel_size = image.pixel_size;
    const ALPHA: f32 = 0.7;
    for i in 0..image.data.len() {
        out.data[i] = image.data[i] + ALPHA * (image.data[i] - blurred.data[i]);
    }
    out
}

/// A simple box blur with the given radius (used by the ramp filter).
fn box_blur(image: &Image2d, radius: usize) -> Image2d {
    let (w, h) = (image.width, image.height);
    let mut out = Image2d::zeros(w, h);
    out.pixel_size = image.pixel_size;
    let r = radius as i64;
    for y in 0..h {
        for x in 0..w {
            let mut sum = 0.0f64;
            let mut count = 0.0f64;
            for dy in -r..=r {
                for dx in -r..=r {
                    let nx = x as i64 + dx;
                    let ny = y as i64 + dy;
                    if nx >= 0 && ny >= 0 && (nx as usize) < w && (ny as usize) < h {
                        sum += image.data[ny as usize * w + nx as usize] as f64;
                        count += 1.0;
                    }
                }
            }
            out.data[y * w + x] = (sum / count) as f32;
        }
    }
    out
}

/// Bilinearly samples an image at a fractional pixel; `None` outside.
fn bilinear(image: &Image2d, x: f64, y: f64) -> Option<f64> {
    if x < 0.0 || y < 0.0 || x > (image.width - 1) as f64 || y > (image.height - 1) as f64 {
        return None;
    }
    let x0 = x.floor() as usize;
    let y0 = y.floor() as usize;
    let x1 = (x0 + 1).min(image.width - 1);
    let y1 = (y0 + 1).min(image.height - 1);
    let fx = x - x0 as f64;
    let fy = y - y0 as f64;
    let g = |xi: usize, yi: usize| image.data[yi * image.width + xi] as f64;
    let c0 = g(x0, y0) * (1.0 - fx) + g(x1, y0) * fx;
    let c1 = g(x0, y1) * (1.0 - fx) + g(x1, y1) * fx;
    Some(c0 * (1.0 - fy) + c1 * fy)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A volume with a bright cube in the centre.
    fn centered_blob(n: usize) -> Volume3d {
        let mut v = Volume3d::zeros_cube(n);
        let c = n / 2;
        let r = (n / 6).max(1) as i64;
        for z in 0..n {
            for y in 0..n {
                for x in 0..n {
                    let dx = x as i64 - c as i64;
                    let dy = y as i64 - c as i64;
                    let dz = z as i64 - c as i64;
                    if dx.abs() <= r && dy.abs() <= r && dz.abs() <= r {
                        v.data[(z * n + y) * n + x] = 1.0;
                    }
                }
            }
        }
        v
    }

    #[test]
    fn projection_of_a_blob_is_brightest_at_the_centre() {
        let vol = centered_blob(20);
        let proj = project_volume(&vol, 0.0, 0.0, 0.0).expect("project");
        let c = proj.width / 2;
        let centre = proj.at(c, c).unwrap();
        let corner = proj.at(0, 0).unwrap();
        assert!(centre > corner, "centre {centre} > corner {corner}");
    }

    #[test]
    fn reconstruction_recovers_a_central_blob() {
        // Project a known blob from many directions, reconstruct,
        // check the density peak is at the centre.
        let n = 16;
        let vol = centered_blob(n);
        let mut projections = Vec::new();
        // A spread of tilt/rot angles.
        for ti in 0..6 {
            for ri in 0..6 {
                let tilt = std::f64::consts::PI * ti as f64 / 6.0;
                let rot = 2.0 * std::f64::consts::PI * ri as f64 / 6.0;
                let img = project_volume(&vol, rot, tilt, 0.0).expect("project");
                projections.push(Projection::new(img, rot, tilt, 0.0));
            }
        }
        let recon = reconstruct_3d(&projections, n).expect("reconstruct");
        assert_eq!(recon.projections_used, 36);
        // The reconstructed density at the centre exceeds the corner.
        let c = n / 2;
        let centre = recon.volume.at(c, c, c).unwrap();
        let corner = recon.volume.at(1, 1, 1).unwrap();
        assert!(
            centre > corner,
            "reconstructed centre {centre} > corner {corner}"
        );
    }

    #[test]
    fn empty_projection_set_rejected() {
        assert!(reconstruct_3d(&[], 16).is_err());
    }

    #[test]
    fn mismatched_projection_size_rejected() {
        let p = Projection::new(Image2d::zeros(8, 8), 0.0, 0.0, 0.0);
        assert!(reconstruct_3d(&[p], 16).is_err());
    }

    #[test]
    fn non_cubic_volume_cannot_be_projected() {
        let v = Volume3d::zeros(8, 8, 4);
        assert!(project_volume(&v, 0.0, 0.0, 0.0).is_err());
    }
}
