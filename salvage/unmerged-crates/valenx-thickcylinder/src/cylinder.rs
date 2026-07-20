//! Validated cylinder geometry and the Lame integration constants.
//!
//! A [`Cylinder`] bundles a strictly-positive inner radius `a` and outer
//! radius `b` with `b > a`. The boundary pressures (internal `p_i`,
//! external `p_o`) close the two-point boundary-value problem and fix the
//! Lame constants `A` and `B`:
//!
//! ```text
//! A = (p_i a^2 - p_o b^2) / (b^2 - a^2)
//! B = a^2 b^2 (p_i - p_o) / (b^2 - a^2)
//! ```
//!
//! These are stored once and reused by the [`crate::stress`] evaluators so
//! the algebra is centralised and tested in exactly one place.

use crate::error::{non_negative, positive, ThickCylinderError};

/// A validated thick-walled cylinder cross-section.
///
/// Units are unconstrained but must be *consistent*: if radii are in
/// millimetres and pressures in megapascals, the returned stresses are in
/// megapascals. The type guarantees `0 < a < b`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Cylinder {
    inner: f64,
    outer: f64,
}

impl Cylinder {
    /// Build a cylinder from inner radius `a` and outer radius `b`.
    ///
    /// # Errors
    ///
    /// Returns [`ThickCylinderError::NonPositive`] if either radius is not
    /// strictly positive, [`ThickCylinderError::NotFinite`] for NaN/Inf,
    /// and [`ThickCylinderError::BadRadii`] unless `b > a`.
    pub fn new(inner: f64, outer: f64) -> Result<Self, ThickCylinderError> {
        let inner = positive("inner_radius", inner)?;
        let outer = positive("outer_radius", outer)?;
        if outer > inner {
            Ok(Self { inner, outer })
        } else {
            Err(ThickCylinderError::BadRadii { inner, outer })
        }
    }

    /// Inner radius `a`.
    #[inline]
    pub fn inner(&self) -> f64 {
        self.inner
    }

    /// Outer radius `b`.
    #[inline]
    pub fn outer(&self) -> f64 {
        self.outer
    }

    /// Wall thickness `t = b - a` (always strictly positive).
    #[inline]
    pub fn thickness(&self) -> f64 {
        self.outer - self.inner
    }

    /// Mean radius `r_m = (a + b) / 2`, used by the thin-wall comparison.
    #[inline]
    pub fn mean_radius(&self) -> f64 {
        0.5 * (self.inner + self.outer)
    }

    /// Diameter ratio `K = b / a`, the natural slenderness parameter of
    /// the Lame solution (`K -> 1` is the thin-wall limit).
    #[inline]
    pub fn diameter_ratio(&self) -> f64 {
        self.outer / self.inner
    }

    /// Check that an evaluation radius lies within the closed wall `[a, b]`.
    ///
    /// # Errors
    ///
    /// [`ThickCylinderError::NotFinite`] for NaN/Inf and
    /// [`ThickCylinderError::OutsideWall`] when `r < a` or `r > b`.
    pub fn check_radius(&self, r: f64) -> Result<f64, ThickCylinderError> {
        let r = crate::error::finite("r", r)?;
        if r < self.inner || r > self.outer {
            Err(ThickCylinderError::OutsideWall {
                inner: self.inner,
                outer: self.outer,
                r,
            })
        } else {
            Ok(r)
        }
    }
}

/// The pair of Lame integration constants `A` and `B` for a given cylinder
/// and boundary pressures.
///
/// With these in hand the stress field is simply
/// `sigma_r(r) = A - B/r^2` and `sigma_theta(r) = A + B/r^2`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LameConstants {
    /// The constant axisymmetric term `A` (mean of hoop and radial stress).
    pub a: f64,
    /// The `1/r^2` amplitude `B`.
    pub b: f64,
}

impl LameConstants {
    /// Solve for `A` and `B` from a [`Cylinder`] under internal pressure
    /// `p_i` and external pressure `p_o`.
    ///
    /// Pressures are taken as non-negative compressive loads (the usual
    /// convention for a contained gauge/fluid pressure). With the sign
    /// convention "tension positive" this yields
    /// `A = (p_i a^2 - p_o b^2)/(b^2 - a^2)` and
    /// `B = a^2 b^2 (p_i - p_o)/(b^2 - a^2)`.
    ///
    /// # Errors
    ///
    /// [`ThickCylinderError::NonPositive`] if a pressure is negative and
    /// [`ThickCylinderError::NotFinite`] for NaN/Inf.
    pub fn solve(
        cyl: &Cylinder,
        p_internal: f64,
        p_external: f64,
    ) -> Result<Self, ThickCylinderError> {
        let p_i = non_negative("p_internal", p_internal)?;
        let p_o = non_negative("p_external", p_external)?;
        let a2 = cyl.inner * cyl.inner;
        let b2 = cyl.outer * cyl.outer;
        let denom = b2 - a2; // strictly positive because b > a > 0
        let a = (p_i * a2 - p_o * b2) / denom;
        let b = a2 * b2 * (p_i - p_o) / denom;
        Ok(Self { a, b })
    }
}
