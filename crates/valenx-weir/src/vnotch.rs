//! Sharp-crested **triangular (V-notch)** weir discharge.
//!
//! A V-notch weir is a triangular opening of full vertex angle `θ` cut
//! into a thin plate; water spills through the notch with upstream head
//! `H` measured from the vertex to the free surface. Because the flow
//! width grows linearly with depth, integrating `√(2 g z)` over the
//! triangular opening yields a **5/2 power** law:
//!
//! ```text
//!   Q = Cd · (8/15) · √(2 g) · tan(θ/2) · H^(5/2)
//! ```
//!
//! The steeper `H^(5/2)` response — versus `H^(3/2)` for a rectangular
//! weir — means a small change in head produces a large, easily-read
//! change in head at **low flows**, which is exactly why a V-notch is
//! the preferred gauging weir for small discharges.

use crate::error::{require_positive, WeirError};
use crate::G_STANDARD;
use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

/// The dimensionless `8/15` prefactor in the V-notch weir equation.
///
/// It is the result of integrating the triangular opening:
/// the local width at depth `z` is `2 (H − z) tan(θ/2)`, and
/// `∫₀ᴴ (H − z) √z dz = (4/15) H^(5/2)`, which with the factor of two
/// from the symmetric notch gives the `8/15` coefficient.
pub const VNOTCH_COEFFICIENT: f64 = 8.0 / 15.0;

/// A sharp-crested triangular (V-notch) weir, validated on construction.
///
/// The struct stores only the three independent quantities that define
/// the discharge: the full vertex angle `θ`, the discharge coefficient
/// `Cd`, and the gravitational acceleration `g`. The derived `tan(θ/2)`
/// is recomputed on demand (see
/// [`half_angle_tangent`](VNotchWeir::half_angle_tangent)) so the type
/// has a single source of truth and serialises losslessly. The head `H`
/// is supplied per evaluation to [`discharge`](VNotchWeir::discharge).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct VNotchWeir {
    /// Full vertex (notch) angle `θ`, in radians. In `(0, π)`.
    vertex_angle_rad: f64,
    /// Dimensionless discharge coefficient `Cd`. Strictly positive;
    /// physically `≈ 0.58`–`0.62` for a fully-contracted V-notch.
    discharge_coefficient: f64,
    /// Gravitational acceleration `g`, in m·s⁻². Strictly positive.
    gravity: f64,
}

impl VNotchWeir {
    /// Construct a V-notch weir at standard gravity
    /// (`g = `[`G_STANDARD`]).
    ///
    /// `vertex_angle_rad` is the **full** notch angle `θ` (e.g. a
    /// 90°-notch is `θ = π/2`); the formula uses `tan(θ/2)`.
    ///
    /// # Errors
    ///
    /// Returns [`WeirError::NotchAngleOutOfRange`] if `vertex_angle_rad`
    /// is not in the open interval `(0, π)`, or a [`WeirError`] if
    /// `discharge_coefficient` is not a finite, strictly-positive
    /// number.
    pub fn new(vertex_angle_rad: f64, discharge_coefficient: f64) -> Result<Self, WeirError> {
        Self::with_gravity(vertex_angle_rad, discharge_coefficient, G_STANDARD)
    }

    /// Construct a 90° V-notch weir at standard gravity — the most
    /// common gauging notch (`θ = π/2`, so `tan(θ/2) = 1`).
    ///
    /// # Errors
    ///
    /// Returns a [`WeirError`] if `discharge_coefficient` is not a
    /// finite, strictly-positive number.
    pub fn ninety_degree(discharge_coefficient: f64) -> Result<Self, WeirError> {
        Self::new(PI / 2.0, discharge_coefficient)
    }

    /// Construct a V-notch weir with an explicit gravitational
    /// acceleration `gravity` (m·s⁻²).
    ///
    /// # Errors
    ///
    /// Returns [`WeirError::NotchAngleOutOfRange`] if `vertex_angle_rad`
    /// is not in `(0, π)`, or a [`WeirError`] if `discharge_coefficient`
    /// or `gravity` is not a finite, strictly-positive number.
    pub fn with_gravity(
        vertex_angle_rad: f64,
        discharge_coefficient: f64,
        gravity: f64,
    ) -> Result<Self, WeirError> {
        if !vertex_angle_rad.is_finite() || vertex_angle_rad <= 0.0 || vertex_angle_rad >= PI {
            return Err(WeirError::NotchAngleOutOfRange {
                radians: vertex_angle_rad,
            });
        }
        Ok(Self {
            vertex_angle_rad,
            discharge_coefficient: require_positive(
                "discharge_coefficient",
                discharge_coefficient,
            )?,
            gravity: require_positive("gravity", gravity)?,
        })
    }

    /// Full vertex (notch) angle `θ`, in radians.
    pub fn vertex_angle_rad(&self) -> f64 {
        self.vertex_angle_rad
    }

    /// Tangent of the half-angle, `tan(θ/2)`.
    ///
    /// Derived from [`vertex_angle_rad`](VNotchWeir::vertex_angle_rad);
    /// for a 90° notch this is exactly `1`.
    pub fn half_angle_tangent(&self) -> f64 {
        (self.vertex_angle_rad / 2.0).tan()
    }

    /// Dimensionless discharge coefficient `Cd`.
    pub fn discharge_coefficient(&self) -> f64 {
        self.discharge_coefficient
    }

    /// Gravitational acceleration `g`, in m·s⁻².
    pub fn gravity(&self) -> f64 {
        self.gravity
    }

    /// Volumetric discharge `Q` (m³·s⁻¹) at upstream head
    /// `head_m` (metres):
    ///
    /// ```text
    ///   Q = Cd · (8/15) · √(2 g) · tan(θ/2) · H^(5/2)
    /// ```
    ///
    /// # Errors
    ///
    /// Returns a [`WeirError`] if `head_m` is not a finite,
    /// strictly-positive number.
    pub fn discharge(&self, head_m: f64) -> Result<f64, WeirError> {
        let h = require_positive("head", head_m)?;
        Ok(self.discharge_coefficient
            * VNOTCH_COEFFICIENT
            * (2.0 * self.gravity).sqrt()
            * self.half_angle_tangent()
            * h.powf(2.5))
    }
}
