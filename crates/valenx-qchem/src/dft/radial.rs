//! Radial quadrature for the atom-centred molecular grid.
//!
//! The radial half of a molecular DFT grid integrates a function of the
//! distance `r` from a nucleus — `∫₀^∞ f(r) r² dr` — with most points
//! packed near the nucleus where the density is large and varies fast.
//!
//! ## The mapping
//!
//! This module uses the **Treutler-Ahlrichs M4** mapping (Treutler &
//! Ahlrichs, *J. Chem. Phys.* 102, 346 (1995)). A Gauss-Chebyshev
//! second-kind abscissa `x ∈ (−1, 1)` is mapped to a radius by
//!
//! ```text
//! r(x) = (ξ / ln2) · (1 + x)^α · ln(2 / (1 − x))
//! ```
//!
//! with `α = 0.6` and an element-dependent size parameter `ξ`. The
//! Jacobian `dr/dx` folds into the quadrature weight together with the
//! `r²` volume factor, so the returned [`RadialPoint`]s already carry
//! the full `w · r²` weight — a sum `Σ_k w_k f(r_k)` *is* the integral
//! `∫₀^∞ f(r) r² dr`.
//!
//! The Gauss-Chebyshev base rule has the closed-form abscissae
//! `x_i = cos(iπ/(n+1))` and weights `π/(n+1) · sin²(iπ/(n+1))`, so the
//! whole quadrature is exact to machine precision with no tabulated
//! constants.

/// A single radial quadrature point: a radius and its full weight.
///
/// The weight already includes the `r²` spherical volume element and
/// the `dr/dx` Jacobian of the Treutler-Ahlrichs mapping, so
/// `∫₀^∞ f(r) r² dr ≈ Σ_k weight_k · f(r_k)`.
#[derive(Copy, Clone, Debug)]
pub struct RadialPoint {
    /// Radius from the nucleus, in bohr.
    pub radius: f64,
    /// Quadrature weight including `r²` and the mapping Jacobian.
    pub weight: f64,
}

/// The Treutler-Ahlrichs M4 mapping exponent `α`.
const TA_ALPHA: f64 = 0.6;

/// Build an `n`-point Treutler-Ahlrichs radial quadrature with size
/// parameter `xi` (bohr).
///
/// `xi` scales the radial extent: a larger value pushes points further
/// out, appropriate for a more diffuse atom. The returned points run
/// from the innermost radius outward.
pub fn treutler_ahlrichs(n: usize, xi: f64) -> Vec<RadialPoint> {
    let mut out = Vec::with_capacity(n);
    let np1 = (n + 1) as f64;
    let ln2 = std::f64::consts::LN_2;
    for i in 1..=n {
        // Gauss-Chebyshev (2nd kind) abscissa and weight on (-1, 1).
        let theta = i as f64 * std::f64::consts::PI / np1;
        let x = theta.cos();
        let w_cheb = std::f64::consts::PI / np1 * theta.sin().powi(2);
        // The Gauss-Chebyshev rule integrates g(x)·√(1−x²); divide the
        // weight by that factor to get a plain ∫ g(x) dx weight.
        let sqrt_factor = (1.0 - x * x).sqrt();
        let w_plain = if sqrt_factor > 1.0e-14 {
            w_cheb / sqrt_factor
        } else {
            0.0
        };

        // Treutler-Ahlrichs M4 radial map r(x) and its derivative.
        let one_plus = 1.0 + x;
        let one_minus = 1.0 - x;
        let ln_term = (2.0 / one_minus).ln();
        let pref = xi / ln2;
        let r = pref * one_plus.powf(TA_ALPHA) * ln_term;
        // dr/dx = pref [ α (1+x)^{α−1} ln_term + (1+x)^α / (1−x) ].
        let drdx = pref
            * (TA_ALPHA * one_plus.powf(TA_ALPHA - 1.0) * ln_term
                + one_plus.powf(TA_ALPHA) / one_minus);

        // Full weight: w_plain · (dr/dx) · r²  (the r² volume element).
        let weight = w_plain * drdx * r * r;
        out.push(RadialPoint { radius: r, weight });
    }
    out
}

/// A Bragg-Slater-style atomic size parameter `ξ` (bohr) for the
/// Treutler-Ahlrichs radial map, indexed by atomic number `Z = 1..18`.
///
/// These are the Treutler-Ahlrichs `ξ` values (their Table 1) for the
/// first-row atoms, extended with reasonable Bragg-Slater radii for the
/// rest of the supported range; `Z` outside `1..=18` falls back to the
/// carbon value.
pub fn atomic_radial_scale(z: u8) -> f64 {
    // Treutler-Ahlrichs ξ for H..Ne, then Bragg-Slater radii Na..Ar.
    const XI: [f64; 18] = [
        0.8, 0.9, // H, He
        1.8, 1.4, 1.3, 1.1, 0.9, 0.9, 0.9, 0.9, // Li..Ne
        1.4, 1.3, 1.3, 1.2, 1.1, 1.0, 1.0, 1.0, // Na..Ar
    ];
    if z >= 1 && (z as usize) <= XI.len() {
        XI[(z - 1) as usize]
    } else {
        1.1 // carbon-like default
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `∫₀^∞ e^{-r} r² dr = 2`. The radial quadrature integrates a
    /// nucleus-centred exponential to its analytic value.
    #[test]
    fn integrates_exponential_r_squared() {
        let grid = treutler_ahlrichs(60, 1.0);
        let integral: f64 = grid
            .iter()
            .map(|p| p.weight * (-p.radius).exp())
            .sum();
        assert!((integral - 2.0).abs() < 1.0e-4, "∫ = {integral}");
    }

    /// `∫₀^∞ e^{-2r} r² dr = 1/4` — the radial density of a hydrogen 1s
    /// orbital `|ψ|² = e^{-2r}/π` integrates (times 4π) to 1 electron.
    #[test]
    fn integrates_hydrogen_1s_radial_density() {
        let grid = treutler_ahlrichs(80, 1.0);
        let radial: f64 = grid
            .iter()
            .map(|p| p.weight * (-2.0 * p.radius).exp())
            .sum();
        // ∫ e^{-2r} r² dr = 2!/2³ = 1/4.
        assert!((radial - 0.25).abs() < 1.0e-5, "∫ = {radial}");
        // Full normalisation: 4π · (1/π) · ∫ = 4 · (1/4) = 1.
        let norm = 4.0 * std::f64::consts::PI / std::f64::consts::PI * radial;
        assert!((norm - 1.0).abs() < 1.0e-4, "norm = {norm}");
    }

    /// `∫₀^∞ e^{-r²} r² dr = √π / 4` — a Gaussian radial integral.
    #[test]
    fn integrates_gaussian_r_squared() {
        let grid = treutler_ahlrichs(70, 1.0);
        let integral: f64 = grid
            .iter()
            .map(|p| p.weight * (-p.radius * p.radius).exp())
            .sum();
        let exact = std::f64::consts::PI.sqrt() / 4.0;
        assert!((integral - exact).abs() < 1.0e-5, "∫ = {integral}");
    }

    #[test]
    fn all_radii_positive_and_weights_finite() {
        let grid = treutler_ahlrichs(40, 1.2);
        assert_eq!(grid.len(), 40);
        for p in &grid {
            assert!(p.radius > 0.0, "radius {}", p.radius);
            assert!(p.weight.is_finite(), "weight {}", p.weight);
            assert!(p.weight >= 0.0, "weight {} negative", p.weight);
        }
    }

    #[test]
    fn larger_xi_pushes_points_outward() {
        let small = treutler_ahlrichs(30, 0.7);
        let large = treutler_ahlrichs(30, 1.8);
        // The outermost point is further out with the larger scale.
        let r_small = small.iter().map(|p| p.radius).fold(0.0, f64::max);
        let r_large = large.iter().map(|p| p.radius).fold(0.0, f64::max);
        assert!(r_large > r_small);
    }

    #[test]
    fn atomic_scale_is_defined_for_supported_elements() {
        for z in 1..=18u8 {
            assert!(atomic_radial_scale(z) > 0.0);
        }
        // Out-of-range falls back to a positive default.
        assert!(atomic_radial_scale(99) > 0.0);
    }
}
