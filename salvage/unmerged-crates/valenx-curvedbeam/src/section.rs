//! Cross-section geometry: area, centroidal radius, neutral-axis radius,
//! eccentricity, and the second moment of area.
//!
//! The whole of Winkler-Bach curved-beam theory rests on one section-level
//! quantity: the radius of the neutral axis,
//!
//! ```text
//! R_n = A / integral(dA / r)
//! ```
//!
//! where the integral runs over the cross-section and `r` is the radial
//! distance of an area element from the centre of curvature. The neutral
//! axis sits *inside* the centroid (closer to the centre of curvature), and
//! the gap
//!
//! ```text
//! e = R - R_n   (>= 0)
//! ```
//!
//! is the eccentricity that appears in the stress formula. This module
//! provides closed-form `integral(dA/r)` evaluations for the two classic
//! textbook sections.

use crate::error::{require_positive, CurvedBeamError};

/// A curved-beam cross-section with closed-form section properties.
#[derive(Copy, Clone, Debug, PartialEq)]
pub enum Section {
    /// Rectangle of radial height `h` and width `b`, spanning inner radius
    /// `r_i` to outer radius `r_o = r_i + h`.
    Rectangle {
        /// Width (perpendicular to the plane of curvature).
        b: f64,
        /// Inner-fibre radius (closest to the centre of curvature).
        r_i: f64,
        /// Outer-fibre radius.
        r_o: f64,
    },
    /// Solid circle of radius `c`, centroid at radius `r_bar` from the
    /// centre of curvature (so inner fibre at `r_bar - c`, outer at
    /// `r_bar + c`).
    Circle {
        /// Section radius.
        c: f64,
        /// Centroidal radius.
        r_bar: f64,
    },
}

/// Derived section properties shared by every stress evaluation.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct SectionProps {
    /// Cross-sectional area `A`.
    pub area: f64,
    /// Centroidal radius `R` (centre of curvature to centroid).
    pub r_centroid: f64,
    /// Neutral-axis radius `R_n`.
    pub r_neutral: f64,
    /// Eccentricity `e = R - R_n` (>= 0).
    pub eccentricity: f64,
    /// Inner-fibre radius `r_i`.
    pub r_inner: f64,
    /// Outer-fibre radius `r_o`.
    pub r_outer: f64,
    /// Second moment of area about the centroidal axis, `I` (used only by
    /// the straight-beam `My/I` cross-check).
    pub inertia: f64,
}

impl Section {
    /// Build a rectangular section from width, inner radius, and outer
    /// radius.
    ///
    /// # Errors
    ///
    /// Returns [`CurvedBeamError::NotPositive`] if `b` or `r_i` is not
    /// strictly positive (or non-finite), and [`CurvedBeamError::BadRadii`]
    /// if `r_i >= r_o`.
    pub fn rectangle(b: f64, r_i: f64, r_o: f64) -> Result<Self, CurvedBeamError> {
        let b = require_positive("b", b)?;
        let r_i = require_positive("r_i", r_i)?;
        let r_o = require_positive("r_o", r_o)?;
        if r_i >= r_o {
            return Err(CurvedBeamError::BadRadii { r_i, r_o });
        }
        Ok(Section::Rectangle { b, r_i, r_o })
    }

    /// Build a solid-circular section from section radius `c` and
    /// centroidal radius `r_bar`.
    ///
    /// # Errors
    ///
    /// Returns [`CurvedBeamError::NotPositive`] if `c` or `r_bar` is not
    /// strictly positive (or non-finite), and
    /// [`CurvedBeamError::SectionTooLarge`] if `r_bar <= c` (the bore would
    /// reach the centre of curvature).
    pub fn circle(c: f64, r_bar: f64) -> Result<Self, CurvedBeamError> {
        let c = require_positive("c", c)?;
        let r_bar = require_positive("r_bar", r_bar)?;
        if r_bar <= c {
            return Err(CurvedBeamError::SectionTooLarge { r_bar, c });
        }
        Ok(Section::Circle { c, r_bar })
    }

    /// Inner-fibre radius `r_i`.
    pub fn r_inner(&self) -> f64 {
        match *self {
            Section::Rectangle { r_i, .. } => r_i,
            Section::Circle { c, r_bar } => r_bar - c,
        }
    }

    /// Outer-fibre radius `r_o`.
    pub fn r_outer(&self) -> f64 {
        match *self {
            Section::Rectangle { r_o, .. } => r_o,
            Section::Circle { c, r_bar } => r_bar + c,
        }
    }

    /// Cross-sectional area `A`.
    pub fn area(&self) -> f64 {
        match *self {
            Section::Rectangle { b, r_i, r_o } => b * (r_o - r_i),
            Section::Circle { c, .. } => std::f64::consts::PI * c * c,
        }
    }

    /// Centroidal radius `R`.
    pub fn r_centroid(&self) -> f64 {
        match *self {
            Section::Rectangle { r_i, r_o, .. } => 0.5 * (r_i + r_o),
            Section::Circle { r_bar, .. } => r_bar,
        }
    }

    /// Neutral-axis radius `R_n = A / integral(dA/r)`.
    ///
    /// Rectangle: `R_n = h / ln(r_o / r_i)` with `h = r_o - r_i`.
    ///
    /// Solid circle: `R_n = (R + sqrt(R^2 - c^2)) / 2`, the closed form of
    /// `c^2 / (2 (R - sqrt(R^2 - c^2)))`.
    pub fn r_neutral(&self) -> f64 {
        match *self {
            Section::Rectangle { r_i, r_o, .. } => (r_o - r_i) / (r_o / r_i).ln(),
            Section::Circle { c, r_bar } => 0.5 * (r_bar + (r_bar * r_bar - c * c).sqrt()),
        }
    }

    /// Second moment of area about the centroidal axis, `I`.
    ///
    /// Rectangle: `b h^3 / 12`. Solid circle: `pi c^4 / 4`. Used only by
    /// the straight-beam limiting cross-check, never by the curved-beam
    /// stress formula itself.
    pub fn inertia(&self) -> f64 {
        match *self {
            Section::Rectangle { b, r_i, r_o } => {
                let h = r_o - r_i;
                b * h * h * h / 12.0
            }
            Section::Circle { c, .. } => std::f64::consts::PI * c * c * c * c / 4.0,
        }
    }

    /// Bundle every derived property in a single struct.
    pub fn props(&self) -> SectionProps {
        let r_centroid = self.r_centroid();
        let r_neutral = self.r_neutral();
        SectionProps {
            area: self.area(),
            r_centroid,
            r_neutral,
            eccentricity: r_centroid - r_neutral,
            r_inner: self.r_inner(),
            r_outer: self.r_outer(),
            inertia: self.inertia(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    #[test]
    fn rectangle_area_and_centroid() {
        // b = 2, r_i = 4, r_o = 6 -> h = 2, A = 4, R = 5.
        let s = Section::rectangle(2.0, 4.0, 6.0).unwrap();
        assert!((s.area() - 4.0).abs() < EPS);
        assert!((s.r_centroid() - 5.0).abs() < EPS);
        assert!((s.inertia() - (2.0 * 8.0 / 12.0)).abs() < EPS);
    }

    #[test]
    fn rectangle_neutral_axis_hand_value() {
        // R_n = h / ln(r_o/r_i) = 2 / ln(6/4) = 2 / 0.4054651081... .
        let s = Section::rectangle(2.0, 4.0, 6.0).unwrap();
        let expected = 2.0 / (6.0_f64 / 4.0).ln();
        assert!((s.r_neutral() - expected).abs() < EPS);
        // Ground truth: the neutral axis lies inside the centroid.
        assert!(s.r_neutral() < s.r_centroid());
        // Independent decimal pin: 2 / ln(1.5) = 4.9326069247... .
        assert!((s.r_neutral() - 4.932606924752863).abs() < 1e-9);
    }

    #[test]
    fn circle_neutral_axis_two_equivalent_forms_agree() {
        // c = 1, R = 5. Both R_n forms must coincide.
        let s = Section::circle(1.0, 5.0).unwrap();
        let closed = 0.5 * (5.0 + (25.0_f64 - 1.0).sqrt());
        let from_integral = 1.0 / (2.0 * (5.0 - (25.0_f64 - 1.0).sqrt()));
        assert!((closed - from_integral).abs() < 1e-9);
        assert!((s.r_neutral() - closed).abs() < EPS);
        assert!(s.r_neutral() < s.r_centroid());
    }

    #[test]
    fn circle_area_and_inertia() {
        let s = Section::circle(2.0, 10.0).unwrap();
        assert!((s.area() - std::f64::consts::PI * 4.0).abs() < EPS);
        assert!((s.inertia() - std::f64::consts::PI * 16.0 / 4.0).abs() < EPS);
        assert!((s.r_inner() - 8.0).abs() < EPS);
        assert!((s.r_outer() - 12.0).abs() < EPS);
    }

    #[test]
    fn eccentricity_is_positive_and_small_relative_to_radius() {
        let s = Section::rectangle(1.0, 9.0, 11.0).unwrap();
        let p = s.props();
        assert!(p.eccentricity > 0.0);
        // For a mildly curved beam e << R.
        assert!(p.eccentricity < 0.05 * p.r_centroid);
    }

    #[test]
    fn rectangle_rejects_bad_radii() {
        let err = Section::rectangle(1.0, 6.0, 4.0).unwrap_err();
        assert_eq!(err.code(), "curvedbeam.bad_radii");
    }

    #[test]
    fn rectangle_rejects_nonpositive_width() {
        let err = Section::rectangle(0.0, 4.0, 6.0).unwrap_err();
        assert_eq!(err.code(), "curvedbeam.not_positive");
    }

    #[test]
    fn circle_rejects_section_too_large() {
        // r_bar = c -> bore reaches centre of curvature.
        let err = Section::circle(5.0, 5.0).unwrap_err();
        assert_eq!(err.code(), "curvedbeam.section_too_large");
    }

    #[test]
    fn rejects_non_finite() {
        let err = Section::rectangle(f64::NAN, 4.0, 6.0).unwrap_err();
        assert_eq!(err.code(), "curvedbeam.not_finite");
    }
}
