//! Cross-section geometry: second moment of area and extreme-fibre
//! distance.
//!
//! The bending behaviour of an Euler-Bernoulli beam depends on its
//! cross-section only through two scalars:
//!
//! - the **second moment of area** `I` (a.k.a. area moment of inertia)
//!   about the neutral (centroidal, bending) axis, and
//! - the **extreme-fibre distance** `c`, the distance from the neutral
//!   axis to the outermost fibre, which sets the peak bending stress
//!   `sigma = M*c/I`.
//!
//! Two textbook closed forms are provided:
//!
//! - **Solid rectangle** (width `b`, height `h`, bending about the
//!   horizontal centroidal axis): `I = b*h^3 / 12`, `c = h / 2`.
//! - **Solid circle** (diameter `d`): `I = pi*d^4 / 64`, `c = d / 2`.
//!
//! All lengths are in consistent units (e.g. millimetres); `I` then
//! comes out in length^4.

use crate::error::BeamError;
use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

/// A beam cross-section, reduced to the two scalars that drive
/// Euler-Bernoulli bending: the second moment of area `I` and the
/// extreme-fibre distance `c`.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum Section {
    /// Solid rectangle, bending about the horizontal centroidal axis.
    Rectangular {
        /// Width `b` (parallel to the neutral axis), length units.
        width: f64,
        /// Height `h` (in the bending direction), length units.
        height: f64,
    },
    /// Solid circle of the given diameter.
    Circular {
        /// Diameter `d`, length units.
        diameter: f64,
    },
}

impl Section {
    /// Construct a validated solid-rectangle section.
    ///
    /// # Errors
    /// Returns [`BeamError::BadParameter`] if `width` or `height` is not
    /// finite and strictly positive.
    pub fn rectangular(width: f64, height: f64) -> Result<Self, BeamError> {
        BeamError::require_positive("width", width)?;
        BeamError::require_positive("height", height)?;
        Ok(Section::Rectangular { width, height })
    }

    /// Construct a validated solid-circle section.
    ///
    /// # Errors
    /// Returns [`BeamError::BadParameter`] if `diameter` is not finite
    /// and strictly positive.
    pub fn circular(diameter: f64) -> Result<Self, BeamError> {
        BeamError::require_positive("diameter", diameter)?;
        Ok(Section::Circular { diameter })
    }

    /// Second moment of area `I` about the neutral (bending) axis, in
    /// length^4.
    ///
    /// - Rectangle: `I = b*h^3 / 12`.
    /// - Circle: `I = pi*d^4 / 64`.
    pub fn second_moment_area(&self) -> f64 {
        match *self {
            Section::Rectangular { width, height } => width * height.powi(3) / 12.0,
            Section::Circular { diameter } => PI * diameter.powi(4) / 64.0,
        }
    }

    /// Distance `c` from the neutral axis to the extreme fibre, in
    /// length units.
    ///
    /// - Rectangle: `c = h / 2`.
    /// - Circle: `c = d / 2`.
    pub fn extreme_fibre(&self) -> f64 {
        match *self {
            Section::Rectangular { height, .. } => height / 2.0,
            Section::Circular { diameter } => diameter / 2.0,
        }
    }

    /// Cross-sectional area, in length^2.
    ///
    /// Not used by the bending formulas themselves, but handy for
    /// self-weight and reporting.
    ///
    /// - Rectangle: `A = b*h`.
    /// - Circle: `A = pi*d^2 / 4`.
    pub fn area(&self) -> f64 {
        match *self {
            Section::Rectangular { width, height } => width * height,
            Section::Circular { diameter } => PI * diameter * diameter / 4.0,
        }
    }

    /// Elastic section modulus `S = I / c`, in length^3.
    ///
    /// Lets you write the peak stress as `sigma = M / S`.
    ///
    /// # Errors
    /// Returns [`BeamError::DegenerateSection`] if `c` vanishes (only
    /// reachable through a hand-constructed enum variant, since the
    /// validated constructors forbid zero dimensions).
    pub fn section_modulus(&self) -> Result<f64, BeamError> {
        let c = self.extreme_fibre();
        if c.is_finite() && c > 0.0 {
            Ok(self.second_moment_area() / c)
        } else {
            Err(BeamError::DegenerateSection {
                reason: "extreme-fibre distance is zero",
            })
        }
    }
}
