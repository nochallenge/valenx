//! Irradiance-volume global illumination — a 3-D grid of light probes
//! for **real-time** indirect lighting.
//!
//! ## What this is
//!
//! A real-time renderer cannot path-trace every pixel, but it still
//! wants the indirect light a path tracer captures — the soft fill, the
//! colour bleed off a red wall. The classic answer (Tatarchuk's
//! "irradiance volumes", and the technique behind every modern engine's
//! baked GI) is to **precompute** the indirect light at a sparse 3-D
//! grid of sample points — *light probes* — and at run time look it up
//! and interpolate.
//!
//! This module ships that:
//!
//! - [`IrradianceVolume`] — a regular 3-D grid of probes spanning an
//!   axis-aligned world box.
//! - **Baking** ([`IrradianceVolume::bake`]) — for each probe, the
//!   incident radiance is gathered by **hemisphere / sphere sampling**
//!   the scene and the result is projected onto **spherical harmonics**
//!   (`L1`, 4 coefficients, or `L2`, 9). The gather is driven by a
//!   caller-supplied radiance closure, so the path tracer (or any
//!   bounded ray-gather) provides the scene radiance and this crate
//!   stays dependency-light.
//! - **Lookup** ([`IrradianceVolume::sample_irradiance`]) — given a
//!   world position and a surface normal, the 8 surrounding probes are
//!   **trilinearly blended** and their SH evaluated in the normal
//!   direction, yielding the irradiance a Lambertian surface there
//!   receives from the indirect environment.
//!
//! ## Spherical harmonics
//!
//! Each probe stores the incident *radiance* field `L(ω)` as a
//! low-order SH expansion — 4 (`L1`) or 9 (`L2`) RGB coefficients.
//! SH is the natural basis for this: the irradiance a Lambertian
//! surface receives is a **cosine-weighted** integral of the radiance,
//! and that convolution is, in the SH basis, just a per-band scalar
//! multiply (Ramamoorthi & Hanrahan 2001 — the `A_l` coefficients
//! `π, 2π/3, π/4`). So an `L2` probe — 9 numbers per channel — captures
//! the entire diffuse-lighting environment essentially exactly, and a
//! lookup is a 9-term dot product.
//!
//! ## Honest scope — the baking is real and fully verifiable
//!
//! The **baking is CPU and exact**: the SH projection, the
//! cosine-convolution lookup, the trilinear blend are all the genuine
//! maths and are unit-tested (uniform scene → uniform irradiance, the
//! SH encode/decode round-trips a constant, a probe by a coloured wall
//! picks up the bleed). What this is *not*: it is a **static** bake
//! (no run-time probe relighting), it has **no visibility term** —
//! a probe blends into a surface even through a thin wall, the classic
//! "light leak" of un-occluded probe GI (a real engine adds a
//! per-probe depth/chebyshev occlusion test — a documented follow-up),
//! and it bakes **one bounce of the supplied radiance** (feed the
//! gather a path-traced radiance for multi-bounce). Each limitation is
//! additive and does not affect the correctness of the SH maths here.

use crate::error::RenderError;

/// The spherical-harmonic order a probe stores its radiance at.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ShOrder {
    /// **L1** — 4 coefficients per colour channel (the constant band
    /// plus the 3 linear bands). Captures a directional + ambient
    /// term; cheap, and enough for soft fill light.
    L1,
    /// **L2** — 9 coefficients per colour channel (adds the 5
    /// quadratic bands). The industry-standard order for baked diffuse
    /// GI — it reconstructs a Lambertian lighting environment to within
    /// a few percent.
    L2,
}

impl ShOrder {
    /// The number of SH coefficients (per colour channel) this order
    /// uses — 4 for [`ShOrder::L1`], 9 for [`ShOrder::L2`].
    #[inline]
    pub fn coeff_count(self) -> usize {
        match self {
            ShOrder::L1 => 4,
            ShOrder::L2 => 9,
        }
    }
}

/// One light probe — the indirect radiance at a grid point, stored as
/// a spherical-harmonic expansion.
///
/// `coeffs` holds `order.coeff_count()` RGB triples: `coeffs[i]` is the
/// `i`-th SH basis function's coefficient. Build one with
/// [`IrradianceVolume::bake`]; a probe is rarely constructed by hand.
#[derive(Clone, Debug, PartialEq)]
pub struct LightProbe {
    /// The SH order — fixes how many `coeffs` are meaningful.
    pub order: ShOrder,
    /// The SH coefficients, one RGB triple per basis function.
    pub coeffs: Vec<[f32; 3]>,
}

impl LightProbe {
    /// A black (all-zero) probe of the given order.
    pub fn zero(order: ShOrder) -> LightProbe {
        LightProbe {
            order,
            coeffs: vec![[0.0; 3]; order.coeff_count()],
        }
    }

    /// Evaluate the stored **radiance** in direction `dir` — the SH
    /// series reconstructed: `Σ coeffs[i]·Y_i(dir)`.
    ///
    /// This is the raw radiance field, *not* the cosine-convolved
    /// irradiance — use [`Self::irradiance`] for the value a
    /// Lambertian surface receives.
    pub fn radiance(&self, dir: [f32; 3]) -> [f32; 3] {
        let basis = sh_basis(dir, self.order);
        let mut out = [0.0f32; 3];
        for (c, &b) in self.coeffs.iter().zip(basis.iter()) {
            out[0] += c[0] * b;
            out[1] += c[1] * b;
            out[2] += c[2] * b;
        }
        out
    }

    /// Evaluate the **irradiance** received by a Lambertian surface
    /// with outward normal `normal` — the cosine-weighted hemisphere
    /// integral of the stored radiance, divided by π so it is directly
    /// the outgoing-radiance scale a diffuse BRDF multiplies its albedo
    /// by.
    ///
    /// In the SH basis the cosine convolution is a per-band scalar
    /// multiply by the Ramamoorthi-Hanrahan `A_l` factors
    /// (`A_0 = π`, `A_1 = 2π/3`, `A_2 = π/4`); dividing the result by
    /// π folds the `1/π` of the Lambert BRDF in, so a uniform white
    /// radiance probe returns irradiance ≈ the radiance itself.
    pub fn irradiance(&self, normal: [f32; 3]) -> [f32; 3] {
        let basis = sh_basis(normal, self.order);
        // Per-band cosine-convolution weights A_l, already divided by π
        // so the output is the diffuse-radiance scale (irradiance/π).
        let inv_pi = std::f32::consts::FRAC_1_PI;
        let a0 = std::f32::consts::PI * inv_pi; // band 0 → 1
        let a1 = (2.0 / 3.0) * std::f32::consts::PI * inv_pi; // band 1 → 2/3
        let a2 = (std::f32::consts::PI / 4.0) * inv_pi; // band 2 → 1/4
        let band_weight = |i: usize| -> f32 {
            match i {
                0 => a0,
                1..=3 => a1,
                _ => a2,
            }
        };
        let mut out = [0.0f32; 3];
        for (i, (c, &b)) in self.coeffs.iter().zip(basis.iter()).enumerate() {
            let w = band_weight(i) * b;
            out[0] += c[0] * w;
            out[1] += c[1] * w;
            out[2] += c[2] * w;
        }
        // Irradiance is non-negative; a low-order SH reconstruction can
        // ring slightly negative — clamp.
        [out[0].max(0.0), out[1].max(0.0), out[2].max(0.0)]
    }
}

/// A 3-D grid of [`LightProbe`]s — the irradiance volume.
///
/// The probes sit on a regular `dims.0 × dims.1 × dims.2` lattice
/// spanning the axis-aligned box `[min, max]`: the corner probes are at
/// `min` and `max`, the rest evenly spaced between. A
/// [`Self::sample_irradiance`] lookup trilinearly blends the 8 probes
/// of the cell containing the query point.
#[derive(Clone, Debug, PartialEq)]
pub struct IrradianceVolume {
    /// Minimum corner of the volume's world-space box.
    pub min: [f32; 3],
    /// Maximum corner of the volume's world-space box.
    pub max: [f32; 3],
    /// Probe-grid resolution along each axis (each ≥ 2).
    pub dims: (usize, usize, usize),
    /// The SH order every probe uses.
    pub order: ShOrder,
    /// `dims.0 · dims.1 · dims.2` probes, indexed `x + y·dx + z·dx·dy`.
    pub probes: Vec<LightProbe>,
}

impl IrradianceVolume {
    /// Create an **unbaked** volume — every probe black — spanning the
    /// box `[min, max]` at the given grid resolution and SH order.
    ///
    /// Call [`Self::bake`] to fill the probes from a scene.
    ///
    /// # Errors
    ///
    /// [`RenderError::BadParameter`] if any grid dimension is below 2
    /// (trilinear blending needs at least a 2×2×2 lattice) or the box
    /// is inverted (`max ≤ min` on any axis).
    pub fn new(
        min: [f32; 3],
        max: [f32; 3],
        dims: (usize, usize, usize),
        order: ShOrder,
    ) -> Result<IrradianceVolume, RenderError> {
        if dims.0 < 2 || dims.1 < 2 || dims.2 < 2 {
            return Err(RenderError::BadParameter {
                name: "dims",
                reason: "irradiance-volume grid must be at least 2×2×2".into(),
            });
        }
        for k in 0..3 {
            if max[k] <= min[k] {
                return Err(RenderError::BadParameter {
                    name: "bounds",
                    reason: "irradiance-volume max must exceed min on every axis".into(),
                });
            }
        }
        let count = dims.0 * dims.1 * dims.2;
        Ok(IrradianceVolume {
            min,
            max,
            dims,
            order,
            probes: vec![LightProbe::zero(order); count],
        })
    }

    /// The world-space position of probe `(ix, iy, iz)`.
    pub fn probe_position(&self, ix: usize, iy: usize, iz: usize) -> [f32; 3] {
        let frac = |i: usize, n: usize| -> f32 {
            if n <= 1 {
                0.0
            } else {
                i as f32 / (n - 1) as f32
            }
        };
        let fx = frac(ix, self.dims.0);
        let fy = frac(iy, self.dims.1);
        let fz = frac(iz, self.dims.2);
        [
            self.min[0] + (self.max[0] - self.min[0]) * fx,
            self.min[1] + (self.max[1] - self.min[1]) * fy,
            self.min[2] + (self.max[2] - self.min[2]) * fz,
        ]
    }

    /// Linear index of probe `(ix, iy, iz)` into [`Self::probes`].
    #[inline]
    fn probe_index(&self, ix: usize, iy: usize, iz: usize) -> usize {
        ix + iy * self.dims.0 + iz * self.dims.0 * self.dims.1
    }

    /// **Bake** the volume — gather every probe's incident radiance
    /// from the scene and project it onto spherical harmonics.
    ///
    /// `radiance` is the scene-radiance closure: `radiance(origin,
    /// direction)` returns the linear-RGB radiance arriving at
    /// `origin` from `direction` — exactly a path-tracer ray-gather, or
    /// any bounded scene query. Decoupling the bake from a concrete
    /// renderer this way keeps the crate dependency-light: the path
    /// tracer (which already depends on this crate) supplies the
    /// closure.
    ///
    /// `samples` directions are drawn over the **full sphere** for each
    /// probe (a probe is omnidirectional — it must capture light from
    /// every direction, not just one hemisphere); each sample is
    /// accumulated into the SH coefficients with the Monte-Carlo
    /// projection `coeff_i += L·Y_i(ω) / (samples·pdf)`, `pdf = 1/4π`
    /// for the uniform-sphere sampler. 64–256 samples is typical.
    ///
    /// The directions are generated by a deterministic Fibonacci-sphere
    /// sequence, so a bake is exactly reproducible.
    pub fn bake<F>(&mut self, samples: usize, mut radiance: F)
    where
        F: FnMut([f32; 3], [f32; 3]) -> [f32; 3],
    {
        let samples = samples.max(1);
        let coeff_count = self.order.coeff_count();
        // Monte-Carlo SH projection: each sample contributes
        // L·Y_i(ω)/(N·pdf); uniform-sphere pdf = 1/(4π).
        let norm = 4.0 * std::f32::consts::PI / samples as f32;

        for iz in 0..self.dims.2 {
            for iy in 0..self.dims.1 {
                for ix in 0..self.dims.0 {
                    let origin = self.probe_position(ix, iy, iz);
                    let mut coeffs = vec![[0.0f32; 3]; coeff_count];
                    for s in 0..samples {
                        let dir = fibonacci_sphere(s, samples);
                        let l = radiance(origin, dir);
                        let basis = sh_basis(dir, self.order);
                        for (c, &b) in coeffs.iter_mut().zip(basis.iter()) {
                            c[0] += l[0] * b * norm;
                            c[1] += l[1] * b * norm;
                            c[2] += l[2] * b * norm;
                        }
                    }
                    let idx = self.probe_index(ix, iy, iz);
                    self.probes[idx] = LightProbe {
                        order: self.order,
                        coeffs,
                    };
                }
            }
        }
    }

    /// Look up the indirect **irradiance** at world `position` for a
    /// surface with outward normal `normal`.
    ///
    /// The cell containing `position` is found, its 8 corner probes
    /// **trilinearly blended** by the fractional position within the
    /// cell, and the blended probe's SH evaluated as cosine-convolved
    /// irradiance in `normal`. The result is the diffuse-radiance scale
    /// a Lambertian BRDF multiplies its albedo by (already `/π`).
    ///
    /// A `position` outside the volume box is clamped to the box, so
    /// the query degrades gracefully to the nearest boundary probes
    /// rather than returning black.
    pub fn sample_irradiance(&self, position: [f32; 3], normal: [f32; 3]) -> [f32; 3] {
        // Continuous grid coordinates of the query point, clamped to
        // the lattice so an out-of-box query uses the boundary probes.
        let grid = |k: usize| -> f32 {
            let span = self.max[k] - self.min[k];
            let t = if span.abs() < 1e-12 {
                0.0
            } else {
                (position[k] - self.min[k]) / span
            };
            let n = match k {
                0 => self.dims.0,
                1 => self.dims.1,
                _ => self.dims.2,
            };
            (t.clamp(0.0, 1.0)) * (n - 1) as f32
        };
        let gx = grid(0);
        let gy = grid(1);
        let gz = grid(2);

        // The cell's lower corner and the fractional position in it.
        let x0 = (gx.floor() as usize).min(self.dims.0 - 2);
        let y0 = (gy.floor() as usize).min(self.dims.1 - 2);
        let z0 = (gz.floor() as usize).min(self.dims.2 - 2);
        let fx = gx - x0 as f32;
        let fy = gy - y0 as f32;
        let fz = gz - z0 as f32;

        // Evaluate the irradiance at all 8 corner probes, then
        // trilinearly blend the *results* (equivalent to blending the
        // coefficients first — SH evaluation is linear).
        let mut acc = [0.0f32; 3];
        for dz in 0..2 {
            for dy in 0..2 {
                for dx in 0..2 {
                    let probe = &self.probes
                        [self.probe_index(x0 + dx, y0 + dy, z0 + dz)];
                    let irr = probe.irradiance(normal);
                    // Trilinear weight of this corner.
                    let wx = if dx == 0 { 1.0 - fx } else { fx };
                    let wy = if dy == 0 { 1.0 - fy } else { fy };
                    let wz = if dz == 0 { 1.0 - fz } else { fz };
                    let w = wx * wy * wz;
                    acc[0] += irr[0] * w;
                    acc[1] += irr[1] * w;
                    acc[2] += irr[2] * w;
                }
            }
        }
        acc
    }
}

/// Evaluate the real spherical-harmonic basis functions in direction
/// `dir` (need not be unit length) up to the requested order.
///
/// Returns 4 values for [`ShOrder::L1`] or 9 for [`ShOrder::L2`], in
/// the canonical order: the constant band, then the 3 linear bands
/// `(y, z, x)`, then the 5 quadratic bands. These are the standard
/// real-SH polynomials with their normalisation constants folded in.
pub fn sh_basis(dir: [f32; 3], order: ShOrder) -> Vec<f32> {
    // Normalise; a degenerate direction collapses to the constant band.
    let len = (dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]).sqrt();
    let (x, y, z) = if len < 1e-12 {
        (0.0, 0.0, 0.0)
    } else {
        (dir[0] / len, dir[1] / len, dir[2] / len)
    };
    // Band 0 (l = 0): the constant 0.5·√(1/π).
    let mut b = vec![0.282_094_8_f32];
    // Band 1 (l = 1): linear, 0.488603·{y, z, x}.
    b.push(0.488_602_5 * y);
    b.push(0.488_602_5 * z);
    b.push(0.488_602_5 * x);
    if order == ShOrder::L2 {
        // Band 2 (l = 2): the 5 quadratics.
        b.push(1.092_548_4 * x * y);
        b.push(1.092_548_4 * y * z);
        b.push(0.315_391_57 * (3.0 * z * z - 1.0));
        b.push(1.092_548_4 * x * z);
        b.push(0.546_274_2 * (x * x - y * y));
    }
    b
}

/// The `i`-th of `n` points of the **Fibonacci sphere** — a
/// deterministic, near-uniform spiral over the unit sphere.
///
/// Used as the bake's sample directions: it gives a low-discrepancy
/// spread (far smoother SH-projection convergence than pseudo-random)
/// and is fully reproducible, so a bake is deterministic.
pub fn fibonacci_sphere(i: usize, n: usize) -> [f32; 3] {
    let n = n.max(1);
    // The golden-angle increment.
    let golden = std::f32::consts::PI * (3.0 - 5.0f32.sqrt());
    // z descends linearly from +1 to −1 across the n points (sampling
    // at stratum centres); the radius of the z-slice is √(1 − z²).
    let z = 1.0 - 2.0 * (i as f32 + 0.5) / n as f32;
    let r = (1.0 - z * z).max(0.0).sqrt();
    let theta = golden * i as f32;
    [r * theta.cos(), r * theta.sin(), z]
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The SH basis is **orthonormal**: `∫ Y_i·Y_j dω = δ_ij`. A
    /// uniform-sphere Monte-Carlo estimate of those integrals must
    /// recover the identity matrix — this confirms the basis
    /// polynomials and their normalisation constants.
    #[test]
    fn sh_basis_is_orthonormal() {
        let n = 20_000;
        let order = ShOrder::L2;
        let m = order.coeff_count();
        let mut gram = vec![0.0f64; m * m];
        let pdf = 1.0 / (4.0 * std::f64::consts::PI);
        for s in 0..n {
            let dir = fibonacci_sphere(s, n);
            let b = sh_basis(dir, order);
            for i in 0..m {
                for j in 0..m {
                    gram[i * m + j] += (b[i] * b[j]) as f64 / (n as f64 * pdf);
                }
            }
        }
        for i in 0..m {
            for j in 0..m {
                let expected = if i == j { 1.0 } else { 0.0 };
                let got = gram[i * m + j];
                assert!(
                    (got - expected).abs() < 0.05,
                    "Gram[{i}][{j}] = {got}, expected {expected}"
                );
            }
        }
    }

    /// SH encode → decode **round-trips a constant**: a constant
    /// radiance field, projected onto SH and reconstructed, gives back
    /// (near) the same constant in every direction. The headline SH
    /// correctness check.
    #[test]
    fn sh_round_trips_a_constant_radiance() {
        let constant = [0.7f32, 0.4, 0.9];
        // Project a constant onto SH by Monte-Carlo (the bake maths).
        let n = 4096;
        let order = ShOrder::L2;
        let norm = 4.0 * std::f32::consts::PI / n as f32;
        let mut coeffs = vec![[0.0f32; 3]; order.coeff_count()];
        for s in 0..n {
            let dir = fibonacci_sphere(s, n);
            let b = sh_basis(dir, order);
            for (c, &bi) in coeffs.iter_mut().zip(b.iter()) {
                c[0] += constant[0] * bi * norm;
                c[1] += constant[1] * bi * norm;
                c[2] += constant[2] * bi * norm;
            }
        }
        let probe = LightProbe { order, coeffs };
        // The reconstructed radiance must equal the constant in every
        // sampled direction.
        for &dir in &[
            [1.0f32, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [-0.577, 0.577, -0.577],
        ] {
            let r = probe.radiance(dir);
            for k in 0..3 {
                assert!(
                    (r[k] - constant[k]).abs() < 0.03,
                    "channel {k}: reconstructed {} vs constant {}",
                    r[k],
                    constant[k]
                );
            }
        }
    }

    /// A probe baked in a **uniformly-lit** scene reports **uniform
    /// irradiance** — the same value for every surface normal. A
    /// constant radiance environment has no directional variation, so
    /// the cosine-convolved irradiance is direction-independent.
    #[test]
    fn uniform_scene_gives_uniform_irradiance() {
        let mut vol = IrradianceVolume::new(
            [-1.0, -1.0, -1.0],
            [1.0, 1.0, 1.0],
            (2, 2, 2),
            ShOrder::L2,
        )
        .unwrap();
        // A scene that radiates a constant 0.5 from every direction.
        vol.bake(256, |_origin, _dir| [0.5, 0.5, 0.5]);
        // Every probe, every normal → the same irradiance ≈ 0.5
        // (irradiance() folds in the 1/π, so a uniform radiance L
        // returns ≈ L).
        for probe in &vol.probes {
            for &n in &[
                [0.0f32, 1.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 0.0, -1.0],
                [0.577, 0.577, 0.577],
            ] {
                let irr = probe.irradiance(n);
                for channel in irr {
                    assert!(
                        (channel - 0.5).abs() < 0.03,
                        "uniform scene: irradiance {channel} should be ≈ 0.5"
                    );
                }
            }
        }
    }

    /// A probe near a **coloured wall** picks up the bleed: bake a
    /// scene where one hemisphere of directions is bright red, and the
    /// probe's irradiance for a normal facing that hemisphere must be
    /// red-dominated.
    #[test]
    fn probe_picks_up_colour_bleed_from_a_wall() {
        let mut vol = IrradianceVolume::new(
            [-1.0, -1.0, -1.0],
            [1.0, 1.0, 1.0],
            (2, 2, 2),
            ShOrder::L2,
        )
        .unwrap();
        // A "red wall" on the +X side: rays heading +X see red, all
        // others see black — a colour-bleed scene.
        vol.bake(512, |_origin, dir| {
            if dir[0] > 0.0 {
                [1.0, 0.0, 0.0]
            } else {
                [0.0, 0.0, 0.0]
            }
        });
        let probe = &vol.probes[0];
        // A surface whose normal faces the red wall (+X) collects the
        // red bleed.
        let facing_wall = probe.irradiance([1.0, 0.0, 0.0]);
        assert!(
            facing_wall[0] > facing_wall[2] + 0.05,
            "a surface facing the red wall should pick up red: {facing_wall:?}"
        );
        assert!(facing_wall[0] > 0.05, "the red bleed should be visible");
        // A surface facing *away* from the wall collects far less red.
        let facing_away = probe.irradiance([-1.0, 0.0, 0.0]);
        assert!(
            facing_away[0] < facing_wall[0],
            "facing away should collect less red bleed: {facing_away:?} vs {facing_wall:?}"
        );
    }

    /// Trilinear blending: a query exactly on a probe returns that
    /// probe's value; a query at the cell centre returns the mean of
    /// the 8 corners.
    #[test]
    fn sample_irradiance_trilinearly_blends_probes() {
        let mut vol = IrradianceVolume::new(
            [0.0, 0.0, 0.0],
            [2.0, 2.0, 2.0],
            (2, 2, 2),
            ShOrder::L1,
        )
        .unwrap();
        // Hand-set the 8 probes: corner (0,0,0) bright, the rest dark.
        // A constant-band-only probe (band 0 coefficient) so
        // irradiance is direction-independent and easy to predict.
        let bright = {
            let mut p = LightProbe::zero(ShOrder::L1);
            // A band-0 coefficient that decodes to radiance 1: the
            // basis constant is 0.2820948, so coeff = 1/0.2820948.
            let c = 1.0 / 0.282_094_8;
            p.coeffs[0] = [c, c, c];
            p
        };
        let dark = LightProbe::zero(ShOrder::L1);
        for iz in 0..2 {
            for iy in 0..2 {
                for ix in 0..2 {
                    let idx = vol.probe_index(ix, iy, iz);
                    vol.probes[idx] = if ix == 0 && iy == 0 && iz == 0 {
                        bright.clone()
                    } else {
                        dark.clone()
                    };
                }
            }
        }
        let n = [0.0, 1.0, 0.0];
        // Exactly on the bright corner → its full irradiance (≈ 1).
        let at_corner = vol.sample_irradiance([0.0, 0.0, 0.0], n);
        assert!(
            at_corner[0] > 0.9,
            "query on the bright probe should read it back: {at_corner:?}"
        );
        // At the cell centre → the mean of 8 corners = 1/8 of the
        // bright probe.
        let at_centre = vol.sample_irradiance([1.0, 1.0, 1.0], n);
        assert!(
            (at_centre[0] - at_corner[0] / 8.0).abs() < 0.05,
            "cell-centre query should be 1/8 of the corner: {} vs {}",
            at_centre[0],
            at_corner[0] / 8.0
        );
        // Halfway along the x-edge from the bright corner → half.
        let at_half = vol.sample_irradiance([1.0, 0.0, 0.0], n);
        assert!(
            (at_half[0] - at_corner[0] / 2.0).abs() < 0.05,
            "edge-midpoint query should be 1/2 of the corner: {} vs {}",
            at_half[0],
            at_corner[0] / 2.0
        );
    }

    /// A query outside the volume box clamps to the boundary rather
    /// than returning black or panicking.
    #[test]
    fn out_of_box_query_clamps_to_the_boundary() {
        let mut vol = IrradianceVolume::new(
            [0.0, 0.0, 0.0],
            [1.0, 1.0, 1.0],
            (2, 2, 2),
            ShOrder::L1,
        )
        .unwrap();
        vol.bake(64, |_o, _d| [0.3, 0.3, 0.3]);
        // A point far outside the box.
        let outside = vol.sample_irradiance([100.0, 100.0, 100.0], [0.0, 1.0, 0.0]);
        assert!(
            (outside[0] - 0.3).abs() < 0.05,
            "out-of-box query should clamp to the boundary probe: {outside:?}"
        );
    }

    /// `new` rejects a too-small grid and an inverted box.
    #[test]
    fn new_rejects_bad_parameters() {
        // A 1×2×2 grid cannot trilinearly blend.
        assert!(IrradianceVolume::new(
            [0.0; 3],
            [1.0; 3],
            (1, 2, 2),
            ShOrder::L1
        )
        .is_err());
        // An inverted box.
        assert!(IrradianceVolume::new(
            [0.0; 3],
            [-1.0, 1.0, 1.0],
            (2, 2, 2),
            ShOrder::L1
        )
        .is_err());
        // A valid call.
        assert!(IrradianceVolume::new(
            [0.0; 3],
            [1.0; 3],
            (3, 3, 3),
            ShOrder::L2
        )
        .is_ok());
    }

    /// The Fibonacci-sphere points are unit length and spread over the
    /// whole sphere (their average is near the origin).
    #[test]
    fn fibonacci_sphere_covers_the_sphere() {
        let n = 1000;
        let mut mean = [0.0f64; 3];
        for i in 0..n {
            let p = fibonacci_sphere(i, n);
            let len = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
            assert!((len - 1.0).abs() < 1e-4, "point {i} not unit: {len}");
            mean[0] += p[0] as f64;
            mean[1] += p[1] as f64;
            mean[2] += p[2] as f64;
        }
        // A uniform sphere sampling averages to ~origin.
        for (k, &axis_sum) in mean.iter().enumerate() {
            assert!(
                (axis_sum / n as f64).abs() < 0.05,
                "Fibonacci sphere axis {k} mean {} should be ~0",
                axis_sum / n as f64
            );
        }
    }

    /// `probe_position` places the corner probes exactly on the box
    /// corners and spaces the interior evenly.
    #[test]
    fn probe_positions_span_the_box() {
        let vol = IrradianceVolume::new(
            [-2.0, 0.0, 1.0],
            [2.0, 4.0, 5.0],
            (3, 3, 3),
            ShOrder::L1,
        )
        .unwrap();
        // The (0,0,0) probe sits at min.
        let p000 = vol.probe_position(0, 0, 0);
        assert_eq!(p000, [-2.0, 0.0, 1.0]);
        // The (2,2,2) probe sits at max.
        let p222 = vol.probe_position(2, 2, 2);
        assert_eq!(p222, [2.0, 4.0, 5.0]);
        // The middle probe sits at the box centre.
        let p111 = vol.probe_position(1, 1, 1);
        assert_eq!(p111, [0.0, 2.0, 3.0]);
    }

    /// A bake is **deterministic** — the Fibonacci-sphere sample set is
    /// fixed, so two bakes of the same scene produce identical probes.
    #[test]
    fn bake_is_deterministic() {
        let make = || {
            let mut v = IrradianceVolume::new(
                [0.0; 3],
                [1.0; 3],
                (2, 2, 2),
                ShOrder::L2,
            )
            .unwrap();
            v.bake(128, |o, d| {
                // A position- and direction-dependent radiance.
                [o[0] + d[1].abs(), 0.5, d[2].abs()]
            });
            v
        };
        let a = make();
        let b = make();
        assert_eq!(a.probes, b.probes, "a bake must be reproducible");
    }
}
