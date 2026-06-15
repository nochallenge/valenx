//! Prismatic channel cross-section geometry.
//!
//! ## Model
//!
//! A [`Channel`] is a prismatic open channel whose cross-section is one
//! of the textbook shapes below. For a given flow depth `y` (measured
//! from the channel invert to the free surface) each shape supplies the
//! three geometric quantities that the flow equations consume:
//!
//! - flow area `A` — wetted cross-sectional area,
//! - wetted perimeter `P` — length of the solid boundary in contact
//!   with the water (the free surface is excluded),
//! - top width `T` — width of the free surface.
//!
//! The hydraulic radius then follows as `R = A / P` and the hydraulic
//! (mean) depth as `D = A / T`.
//!
//! ### Rectangular
//!
//! Bottom width `b`, vertical walls. At depth `y`:
//!
//! - `A = b y`
//! - `P = b + 2 y`
//! - `T = b`
//!
//! ### Trapezoidal
//!
//! Bottom width `b` with symmetric side slopes `z` (horizontal run per
//! unit vertical rise, so a wall at run:rise = `z`:1). The rectangular
//! channel is the special case `z = 0`. At depth `y`:
//!
//! - `A = (b + z y) y`
//! - `P = b + 2 y sqrt(1 + z^2)`
//! - `T = b + 2 z y`

use serde::{Deserialize, Serialize};

use crate::error::OpenChannelError;

/// A prismatic open-channel cross-section.
///
/// Build instances with the validated constructors
/// ([`Channel::rectangular`], [`Channel::trapezoidal`]); they reject
/// non-positive or non-finite dimensions so every downstream geometric
/// quantity is well defined.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum Channel {
    /// Rectangular channel with a flat bottom and vertical walls.
    Rectangular {
        /// Bottom width `b` in metres (`> 0`).
        bottom_width_m: f64,
    },
    /// Trapezoidal channel with a flat bottom and symmetric sloping
    /// walls.
    Trapezoidal {
        /// Bottom width `b` in metres (`> 0`).
        bottom_width_m: f64,
        /// Side slope `z` = horizontal run per unit vertical rise
        /// (`>= 0`; `0` degenerates to rectangular).
        side_slope: f64,
    },
}

impl Channel {
    /// Construct a rectangular channel of bottom width `bottom_width_m`.
    ///
    /// # Errors
    ///
    /// Returns [`OpenChannelError`] if `bottom_width_m` is not a finite
    /// strictly-positive number.
    pub fn rectangular(bottom_width_m: f64) -> Result<Self, OpenChannelError> {
        let bottom_width_m = OpenChannelError::non_positive("bottom_width_m", bottom_width_m)?;
        Ok(Self::Rectangular { bottom_width_m })
    }

    /// Construct a trapezoidal channel of bottom width `bottom_width_m`
    /// and symmetric side slope `side_slope` (horizontal:vertical).
    ///
    /// # Errors
    ///
    /// Returns [`OpenChannelError`] if `bottom_width_m` is not finite and
    /// positive, or if `side_slope` is negative / non-finite.
    pub fn trapezoidal(bottom_width_m: f64, side_slope: f64) -> Result<Self, OpenChannelError> {
        let bottom_width_m = OpenChannelError::non_positive("bottom_width_m", bottom_width_m)?;
        let side_slope = OpenChannelError::negative("side_slope", side_slope)?;
        Ok(Self::Trapezoidal {
            bottom_width_m,
            side_slope,
        })
    }

    /// Bottom width `b` of the channel in metres.
    pub fn bottom_width_m(&self) -> f64 {
        match *self {
            Self::Rectangular { bottom_width_m } => bottom_width_m,
            Self::Trapezoidal { bottom_width_m, .. } => bottom_width_m,
        }
    }

    /// Side slope `z` (horizontal run per unit rise). Always `0` for a
    /// rectangular channel.
    pub fn side_slope(&self) -> f64 {
        match *self {
            Self::Rectangular { .. } => 0.0,
            Self::Trapezoidal { side_slope, .. } => side_slope,
        }
    }

    /// Validate a flow depth `y` (must be finite and `> 0`).
    fn check_depth(depth_m: f64) -> Result<f64, OpenChannelError> {
        OpenChannelError::non_positive("depth_m", depth_m)
    }

    /// Flow (wetted) area `A` in m² at flow depth `depth_m`.
    ///
    /// Rectangular: `A = b y`. Trapezoidal: `A = (b + z y) y`.
    ///
    /// # Errors
    ///
    /// Returns [`OpenChannelError`] if `depth_m` is not finite and
    /// positive.
    pub fn area_m2(&self, depth_m: f64) -> Result<f64, OpenChannelError> {
        let y = Self::check_depth(depth_m)?;
        let b = self.bottom_width_m();
        let z = self.side_slope();
        Ok((b + z * y) * y)
    }

    /// Wetted perimeter `P` in m at flow depth `depth_m` (the free
    /// surface is excluded).
    ///
    /// Rectangular: `P = b + 2 y`. Trapezoidal:
    /// `P = b + 2 y sqrt(1 + z^2)`.
    ///
    /// # Errors
    ///
    /// Returns [`OpenChannelError`] if `depth_m` is not finite and
    /// positive.
    pub fn wetted_perimeter_m(&self, depth_m: f64) -> Result<f64, OpenChannelError> {
        let y = Self::check_depth(depth_m)?;
        let b = self.bottom_width_m();
        let z = self.side_slope();
        Ok(b + 2.0 * y * (1.0 + z * z).sqrt())
    }

    /// Free-surface top width `T` in m at flow depth `depth_m`.
    ///
    /// Rectangular: `T = b`. Trapezoidal: `T = b + 2 z y`.
    ///
    /// # Errors
    ///
    /// Returns [`OpenChannelError`] if `depth_m` is not finite and
    /// positive.
    pub fn top_width_m(&self, depth_m: f64) -> Result<f64, OpenChannelError> {
        let y = Self::check_depth(depth_m)?;
        let b = self.bottom_width_m();
        let z = self.side_slope();
        Ok(b + 2.0 * z * y)
    }

    /// Hydraulic radius `R = A / P` in m at flow depth `depth_m`.
    ///
    /// # Errors
    ///
    /// Returns [`OpenChannelError`] if `depth_m` is not finite and
    /// positive.
    pub fn hydraulic_radius_m(&self, depth_m: f64) -> Result<f64, OpenChannelError> {
        let a = self.area_m2(depth_m)?;
        let p = self.wetted_perimeter_m(depth_m)?;
        Ok(a / p)
    }

    /// Hydraulic (mean) depth `D = A / T` in m at flow depth `depth_m`.
    ///
    /// This is the length scale used in the Froude number for a
    /// non-rectangular channel (it equals `y` for a rectangular one).
    ///
    /// # Errors
    ///
    /// Returns [`OpenChannelError`] if `depth_m` is not finite and
    /// positive.
    pub fn hydraulic_depth_m(&self, depth_m: f64) -> Result<f64, OpenChannelError> {
        let a = self.area_m2(depth_m)?;
        let t = self.top_width_m(depth_m)?;
        Ok(a / t)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-12;

    #[test]
    fn rectangular_area_perimeter_topwidth() {
        // b = 3, y = 2  ->  A = 6, P = 3 + 4 = 7, T = 3.
        let ch = Channel::rectangular(3.0).unwrap();
        assert!((ch.area_m2(2.0).unwrap() - 6.0).abs() < EPS);
        assert!((ch.wetted_perimeter_m(2.0).unwrap() - 7.0).abs() < EPS);
        assert!((ch.top_width_m(2.0).unwrap() - 3.0).abs() < EPS);
    }

    #[test]
    fn hydraulic_radius_is_area_over_perimeter() {
        // b = 3, y = 2  ->  R = A/P = 6/7.
        let ch = Channel::rectangular(3.0).unwrap();
        let r = ch.hydraulic_radius_m(2.0).unwrap();
        assert!((r - 6.0 / 7.0).abs() < EPS);
        // Cross-check against the manual A/P composition.
        let a = ch.area_m2(2.0).unwrap();
        let p = ch.wetted_perimeter_m(2.0).unwrap();
        assert!((r - a / p).abs() < EPS);
    }

    #[test]
    fn rectangular_hydraulic_depth_equals_flow_depth() {
        // For a rectangle, D = A/T = (b y)/b = y.
        let ch = Channel::rectangular(4.0).unwrap();
        assert!((ch.hydraulic_depth_m(1.7).unwrap() - 1.7).abs() < EPS);
    }

    #[test]
    fn trapezoidal_geometry_matches_closed_form() {
        // b = 2, z = 1.5, y = 1.
        // A = (2 + 1.5*1)*1 = 3.5
        // P = 2 + 2*1*sqrt(1 + 2.25) = 2 + 2*sqrt(3.25)
        // T = 2 + 2*1.5*1 = 5
        let ch = Channel::trapezoidal(2.0, 1.5).unwrap();
        assert!((ch.area_m2(1.0).unwrap() - 3.5).abs() < EPS);
        let p_expected = 2.0 + 2.0 * 3.25_f64.sqrt();
        assert!((ch.wetted_perimeter_m(1.0).unwrap() - p_expected).abs() < EPS);
        assert!((ch.top_width_m(1.0).unwrap() - 5.0).abs() < EPS);
    }

    #[test]
    fn trapezoid_with_zero_slope_reduces_to_rectangle() {
        let rect = Channel::rectangular(2.5).unwrap();
        let trap = Channel::trapezoidal(2.5, 0.0).unwrap();
        for &y in &[0.5_f64, 1.0, 2.3] {
            assert!((rect.area_m2(y).unwrap() - trap.area_m2(y).unwrap()).abs() < EPS);
            assert!(
                (rect.wetted_perimeter_m(y).unwrap() - trap.wetted_perimeter_m(y).unwrap()).abs()
                    < EPS
            );
            assert!((rect.top_width_m(y).unwrap() - trap.top_width_m(y).unwrap()).abs() < EPS);
        }
    }

    #[test]
    fn constructors_reject_bad_geometry() {
        assert!(Channel::rectangular(0.0).is_err());
        assert!(Channel::rectangular(-1.0).is_err());
        assert!(Channel::trapezoidal(1.0, -0.1).is_err());
        assert!(Channel::trapezoidal(f64::NAN, 1.0).is_err());
    }

    #[test]
    fn non_positive_depth_is_rejected() {
        let ch = Channel::rectangular(2.0).unwrap();
        assert!(ch.area_m2(0.0).is_err());
        assert!(ch.area_m2(-1.0).is_err());
        assert!(ch.hydraulic_radius_m(f64::INFINITY).is_err());
    }
}
