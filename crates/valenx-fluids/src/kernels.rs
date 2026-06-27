//! SPH smoothing kernels (3-D), after Müller, Charypar & Gross,
//! *"Particle-Based Fluid Simulation for Interactive Applications"*,
//! SCA 2003.
//!
//! Smoothed-particle hydrodynamics approximates a field `A` at a point `r` as a
//! sum over neighbour particles `j`, each weighted by a *smoothing kernel*
//! `W(r − r_j, h)` of compact support radius `h`:
//!
//! ```text
//! A(r) ≈ Σ_j  m_j · (A_j / ρ_j) · W(r − r_j, h)
//! ```
//!
//! Müller 2003 uses three different kernels, each chosen for the term it serves:
//!
//! * **poly6** for the **density** sum — smooth and cheap (depends only on
//!   `r²`, so no square root), but its gradient vanishes at the centre, which
//!   makes it a poor choice for a *pressure* force (particles could clump with
//!   no repulsion at very small separation).
//! * **spiky** for the **pressure gradient** — its gradient grows toward the
//!   centre, giving the strong short-range repulsion that resists compression.
//! * **viscosity** for the **viscous Laplacian** — constructed so its Laplacian
//!   is positive everywhere in the support, which keeps viscosity dissipative
//!   (it can only remove relative kinetic energy, never inject it).
//!
//! All three are **compactly supported**: they are exactly `0` for `r > h`, so a
//! particle interacts only with neighbours inside its support ball — the
//! property that makes the spatial hash in [`crate::grid`] effective. Each is
//! normalised in 3-D so that `∫ W dV = 1` over the support; the poly6
//! normalisation is one of the crate's pinned benchmarks (see the tests).
//!
//! # Guarding
//!
//! Every kernel guards its argument: the smoothing length `h` is validated to be
//! finite and strictly positive by the caller ([`SmoothingKernels::new`]), and
//! each scalar kernel returns `0` for `r < 0` or `r > h`. The gradient kernels
//! additionally return the zero vector at `r = 0`, where the radial direction
//! `r̂` is undefined — there is no divide by `r` without first checking `r > 0`.

use nalgebra::Vector3;

use crate::error::FluidError;

/// The three Müller-2003 3-D kernels, precomputed for a fixed smoothing
/// length `h`.
///
/// Constructing a [`SmoothingKernels`] validates `h` once and caches the
/// (otherwise repeatedly recomputed) normalisation constants, so the per-pair
/// kernel evaluations in the hot loop are a handful of multiplies.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SmoothingKernels {
    h: f64,
    h2: f64,
    // Normalisation constants (3-D).
    poly6: f64,      // 315 / (64 π h⁹)
    spiky_grad: f64, // −45 / (π h⁶)   (the scalar in front of r̂, sign included)
    visc_lap: f64,   // 45 / (π h⁶)
}

impl SmoothingKernels {
    /// Build the kernel set for a smoothing length `h`.
    ///
    /// # Errors
    /// [`FluidError::InvalidConfig`] if `h` is not finite and strictly
    /// positive.
    pub fn new(h: f64) -> Result<Self, FluidError> {
        if !(h.is_finite() && h > 0.0) {
            return Err(FluidError::InvalidConfig(format!(
                "smoothing length h must be finite and > 0, got {h}"
            )));
        }
        let pi = std::f64::consts::PI;
        let h3 = h * h * h;
        let h6 = h3 * h3;
        let h9 = h6 * h3;
        Ok(Self {
            h,
            h2: h * h,
            poly6: 315.0 / (64.0 * pi * h9),
            spiky_grad: -45.0 / (pi * h6),
            visc_lap: 45.0 / (pi * h6),
        })
    }

    /// The smoothing length `h` (the compact-support radius).
    #[must_use]
    pub fn h(&self) -> f64 {
        self.h
    }

    /// The poly6 density kernel `W_poly6(r, h)`.
    ///
    /// `W = (315 / 64π h⁹) · (h² − r²)³` for `0 ≤ r ≤ h`, else `0`. Depends only
    /// on `r²`, so callers can pass `r²` directly via [`Self::poly6_r2`] and
    /// avoid the square root.
    #[must_use]
    pub fn poly6(&self, r: f64) -> f64 {
        self.poly6_r2(r * r)
    }

    /// The poly6 kernel evaluated from the *squared* distance `r2 = r²`.
    ///
    /// Returns `0` for `r2 < 0` (nonsensical) or `r2 > h²` (outside support).
    #[must_use]
    pub fn poly6_r2(&self, r2: f64) -> f64 {
        if r2 < 0.0 || r2 > self.h2 {
            return 0.0;
        }
        let d = self.h2 - r2;
        self.poly6 * d * d * d
    }

    /// The gradient of the **spiky** pressure kernel, `∇W_spiky`.
    ///
    /// The spiky kernel is `W_spiky = (15 / π h⁶) · (h − r)³`; its gradient is
    ///
    /// ```text
    /// ∇W_spiky = −(45 / π h⁶) · (h − r)² · r̂
    /// ```
    ///
    /// which points back along `−r̂` (toward the centre particle), giving the
    /// short-range repulsion of the pressure force. `rij` is the displacement
    /// `r_i − r_j`; the returned vector is the gradient with respect to `r_i`.
    /// Returns the zero vector for `r = 0` (the radial direction is undefined
    /// there — **no divide by `r` without `r > 0`**) and for `r > h`.
    #[must_use]
    pub fn spiky_gradient(&self, rij: Vector3<f64>) -> Vector3<f64> {
        let r = rij.norm();
        if r <= 0.0 || r > self.h {
            return Vector3::zeros();
        }
        let d = self.h - r;
        // spiky_grad already carries the −45/(π h⁶) sign.
        let coeff = self.spiky_grad * d * d / r;
        rij * coeff
    }

    /// The **Laplacian** of the viscosity kernel, `∇²W_visc`.
    ///
    /// `∇²W_visc = (45 / π h⁶) · (h − r)` for `0 ≤ r ≤ h`, else `0`. Positive
    /// throughout the support, so the resulting viscous force is strictly
    /// dissipative.
    #[must_use]
    pub fn viscosity_laplacian(&self, r: f64) -> f64 {
        if r < 0.0 || r > self.h {
            return 0.0;
        }
        self.visc_lap * (self.h - r)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rejects_bad_h() {
        assert!(SmoothingKernels::new(0.0).is_err());
        assert!(SmoothingKernels::new(-1.0).is_err());
        assert!(SmoothingKernels::new(f64::NAN).is_err());
        assert!(SmoothingKernels::new(f64::INFINITY).is_err());
        assert!(SmoothingKernels::new(0.1).is_ok());
    }

    #[test]
    fn poly6_is_zero_at_and_beyond_support() {
        let k = SmoothingKernels::new(0.5).unwrap();
        // Exactly at r = h the factor (h² − r²) is 0.
        assert_eq!(k.poly6(0.5), 0.0);
        // Strictly beyond support.
        assert_eq!(k.poly6(0.5 + 1e-9), 0.0);
        assert_eq!(k.poly6(1.0), 0.0);
        // Inside support it is strictly positive and maximal at the centre.
        assert!(k.poly6(0.0) > 0.0);
        assert!(k.poly6(0.1) > k.poly6(0.4));
    }

    // BENCHMARK-PIN (2): the poly6 kernel integrates to ~1 over its support
    // (3-D numerical quadrature in spherical shells), and is exactly 0 beyond h.
    #[test]
    fn poly6_integrates_to_one_over_support() {
        let h = 0.7_f64;
        let k = SmoothingKernels::new(h).unwrap();
        // ∫ W dV = ∫_0^h W(r) · 4π r² dr  (W is radially symmetric).
        let n = 200_000;
        let dr = h / n as f64;
        let mut integral = 0.0;
        for i in 0..n {
            let r = (i as f64 + 0.5) * dr; // midpoint rule
            integral += k.poly6(r) * 4.0 * std::f64::consts::PI * r * r * dr;
        }
        assert!(
            (integral - 1.0).abs() < 1e-4,
            "∫ W_poly6 dV = {integral}, expected ≈ 1"
        );

        // And it is exactly zero past the support, so extending the integral
        // beyond h adds nothing.
        assert_eq!(k.poly6(h * 1.5), 0.0);
    }

    #[test]
    fn spiky_gradient_is_repulsive_and_guarded() {
        let k = SmoothingKernels::new(0.5).unwrap();
        // Zero at the centre (direction undefined) — no divide by zero.
        assert_eq!(k.spiky_gradient(Vector3::zeros()), Vector3::zeros());
        // Zero beyond support.
        assert_eq!(
            k.spiky_gradient(Vector3::new(1.0, 0.0, 0.0)),
            Vector3::zeros()
        );
        // For rij = r_i − r_j pointing +x, the gradient points −x (toward j):
        // the pressure force −∇W then pushes i away from j.
        let g = k.spiky_gradient(Vector3::new(0.1, 0.0, 0.0));
        assert!(g.x < 0.0, "spiky gradient should point back along −r̂");
        assert!(g.y.abs() < 1e-15 && g.z.abs() < 1e-15);
    }

    #[test]
    fn viscosity_laplacian_is_nonnegative_in_support() {
        let k = SmoothingKernels::new(0.5).unwrap();
        assert!(k.viscosity_laplacian(0.0) > 0.0);
        assert!(k.viscosity_laplacian(0.25) > 0.0);
        assert_eq!(k.viscosity_laplacian(0.5), 0.0);
        assert_eq!(k.viscosity_laplacian(0.6), 0.0);
        // Decreases linearly toward the support boundary.
        assert!(k.viscosity_laplacian(0.1) > k.viscosity_laplacian(0.4));
    }
}
