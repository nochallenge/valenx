//! Built-in black-hole spacetimes.
//!
//! All of the rotating/charged family share a single implementation: the
//! **Kerr–Newman** metric in Boyer–Lindquist coordinates is the master, and
//! Schwarzschild / Kerr / Reissner–Nordström are named constructors obtained
//! by switching off spin and/or charge. One code path, four validated special
//! cases. [`Minkowski`] (flat space) is provided separately as a zero-curvature
//! baseline for testing.
//!
//! Geometrized units (`G = c = 1`) throughout: mass `M`, spin `a = J/M` and
//! charge `Q` all have dimensions of length.
//!
//! Boyer–Lindquist metric, signature `(−,+,+,+)`, with
//! `Σ = r² + a²cos²θ` and `Δ = r² − 2Mr + a² + Q²`:
//!
//! ```text
//! g_tt  = −(Δ − a² sin²θ) / Σ
//! g_tφ  = −a sin²θ (r² + a² − Δ) / Σ
//! g_rr  =  Σ / Δ
//! g_θθ  =  Σ
//! g_φφ  =  ((r² + a²)² − Δ a² sin²θ) sin²θ / Σ
//! ```

use crate::autodiff::Scalar;
use crate::metric::{CoordSystem, Spacetime};

/// The Kerr–Newman family: a black hole of mass `M`, spin `a = J/M` and
/// electric charge `Q`, in Boyer–Lindquist coordinates `(t, r, θ, φ)`.
///
/// Construct the special cases with [`schwarzschild`], [`kerr`] or
/// [`reissner_nordstrom`], or build one directly for the fully general case.
#[derive(Clone, Copy, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct KerrNewman {
    /// Mass `M` (geometrized units; length).
    pub mass: f64,
    /// Spin parameter `a = J/M` (length).
    pub spin: f64,
    /// Electric charge `Q` (geometrized units; length).
    pub charge: f64,
}

impl KerrNewman {
    /// `a² + Q² ≤ M²` — the condition for a genuine black hole (an event
    /// horizon exists). When this fails the solution describes a *naked
    /// singularity* instead. Equality is the extremal limit: a degenerate
    /// horizon (`r+ = r−`) with zero surface gravity and Hawking temperature.
    pub fn is_subextremal(&self) -> bool {
        self.spin * self.spin + self.charge * self.charge <= self.mass * self.mass
    }

    /// `M² − a² − Q²`, the discriminant whose square root sets the horizon
    /// locations `r± = M ± √(M² − a² − Q²)`.
    pub fn horizon_discriminant(&self) -> f64 {
        self.mass * self.mass - self.spin * self.spin - self.charge * self.charge
    }
}

impl Spacetime for KerrNewman {
    fn metric<T: Scalar>(&self, x: [T; 4]) -> [[T; 4]; 4] {
        let m = T::from_f64(self.mass);
        let a = T::from_f64(self.spin);
        let q = T::from_f64(self.charge);
        let two = T::from_f64(2.0);
        let zero = T::from_f64(0.0);

        let r = x[1];
        let theta = x[2];
        let sin_t = theta.sin();
        let sin2 = sin_t.sq();
        let cos_t = theta.cos();
        let cos2 = cos_t.sq();

        let r2 = r.sq();
        let a2 = a.sq();
        let sigma = r2 + a2 * cos2;
        let delta = r2 - two * m * r + a2 + q.sq();
        let inv_sigma = sigma.recip();
        let ra2 = r2 + a2;

        let g_tt = -((delta - a2 * sin2) * inv_sigma);
        let g_tphi = -(a * sin2 * (ra2 - delta) * inv_sigma);
        let g_rr = sigma * delta.recip();
        let g_thth = sigma;
        let g_phph = (ra2.sq() - delta * a2 * sin2) * sin2 * inv_sigma;

        [
            [g_tt, zero, zero, g_tphi],
            [zero, g_rr, zero, zero],
            [zero, zero, g_thth, zero],
            [g_tphi, zero, zero, g_phph],
        ]
    }

    fn coords(&self) -> CoordSystem {
        CoordSystem::BoyerLindquist
    }

    fn mass(&self) -> f64 {
        self.mass
    }
}

/// The **Schwarzschild** black hole: non-rotating, uncharged. The canonical
/// solution, and the one with the most closed-form ground truth.
pub fn schwarzschild(mass: f64) -> KerrNewman {
    KerrNewman {
        mass,
        spin: 0.0,
        charge: 0.0,
    }
}

/// The **Kerr** black hole: rotating (`spin = a = J/M`), uncharged.
pub fn kerr(mass: f64, spin: f64) -> KerrNewman {
    KerrNewman {
        mass,
        spin,
        charge: 0.0,
    }
}

/// The **Reissner–Nordström** black hole: charged (`charge = Q`),
/// non-rotating.
pub fn reissner_nordstrom(mass: f64, charge: f64) -> KerrNewman {
    KerrNewman {
        mass,
        spin: 0.0,
        charge,
    }
}

/// Flat **Minkowski** spacetime in Cartesian coordinates, signature
/// `(−,+,+,+)`. Its curvature is identically zero — a baseline that any correct
/// curvature engine must reproduce exactly.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, serde::Serialize)]
pub struct Minkowski;

impl Spacetime for Minkowski {
    fn metric<T: Scalar>(&self, _x: [T; 4]) -> [[T; 4]; 4] {
        let z = T::from_f64(0.0);
        let p = T::from_f64(1.0);
        let n = T::from_f64(-1.0);
        [[n, z, z, z], [z, p, z, z], [z, z, p, z], [z, z, z, p]]
    }

    fn coords(&self) -> CoordSystem {
        CoordSystem::Cartesian
    }

    fn mass(&self) -> f64 {
        0.0
    }
}
