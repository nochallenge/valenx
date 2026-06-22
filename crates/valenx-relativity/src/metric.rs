//! The [`Spacetime`] trait: anything that defines a metric tensor field.

use crate::autodiff::Scalar;

/// The coordinate chart a metric is written in. Recorded so downstream tools
/// (observables, geodesics) know how to interpret the coordinate components.
#[derive(Clone, Copy, Debug, PartialEq, Eq, serde::Serialize)]
pub enum CoordSystem {
    /// Cartesian `(t, x, y, z)` — used by flat (Minkowski) space.
    Cartesian,
    /// Spherical `(t, r, θ, φ)`.
    Spherical,
    /// Boyer–Lindquist `(t, r, θ, φ)` — the standard rotating black-hole chart.
    BoyerLindquist,
}

/// A spacetime geometry: yields the symmetric, lower-index metric tensor
/// `g_μν` at a coordinate point.
///
/// The method is generic over the scalar type `T`, so the very same definition
/// is evaluated with [`f64`] for plain values and with the automatic-
/// differentiation types ([`crate::Dual`], [`crate::HyperDual`]) to obtain the
/// metric's derivatives exactly. This genericity is what makes the curvature
/// engine work for *any* metric without hand-coding its derivatives.
pub trait Spacetime {
    /// The metric tensor `g_μν` at coordinate point `x`, as a symmetric
    /// `4×4` array.
    fn metric<T: Scalar>(&self, x: [T; 4]) -> [[T; 4]; 4];

    /// The coordinate system in which [`Spacetime::metric`] is expressed.
    fn coords(&self) -> CoordSystem;

    /// The gravitating (ADM) mass `M`, in geometrized units (`G = c = 1`).
    fn mass(&self) -> f64;
}
