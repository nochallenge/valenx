//! Closed-form Euler-Bernoulli analyses for cantilever and
//! simply-supported beams.
//!
//! Each analysis takes a [`Beam`] (length, Young's modulus `E`, and a
//! [`Section`]) plus a [`Load`], and returns a [`BeamResult`] holding the
//! maximum deflection, maximum bending moment and the resulting peak
//! bending stress. Every formula is the standard textbook closed form
//! for a prismatic, homogeneous, linear-elastic beam in small-deflection
//! bending.
//!
//! # Sign and unit conventions
//!
//! - All quantities use one consistent unit system. With `E` in N/mm^2
//!   (= MPa), lengths in mm, point loads `P` in N and distributed loads
//!   `w` in N/mm, deflections come out in mm, moments in N·mm and
//!   stresses in N/mm^2 (MPa).
//! - [`BeamResult::max_deflection`] is the **magnitude** of the largest
//!   transverse deflection (always `>= 0`); the physical deflection is
//!   downward for the gravity-direction loads modelled here.
//! - [`BeamResult::max_moment`] and [`BeamResult::max_stress`] are
//!   likewise reported as magnitudes.
//!
//! # The six cases
//!
//! | Support            | Load                | Max deflection        | Max moment   |
//! |--------------------|---------------------|-----------------------|--------------|
//! | Cantilever         | tip point `P`       | `P L^3 / (3 E I)`     | `P L`        |
//! | Cantilever         | UDL `w`             | `w L^4 / (8 E I)`     | `w L^2 / 2`  |
//! | Simply-supported   | centre point `P`    | `P L^3 / (48 E I)`    | `P L / 4`    |
//! | Simply-supported   | UDL `w`             | `5 w L^4 / (384 E I)` | `w L^2 / 8`  |
//! | Fixed-fixed        | centre point `P`    | `P L^3 / (192 E I)`   | `P L / 8`    |
//! | Fixed-fixed        | UDL `w`             | `w L^4 / (384 E I)`   | `w L^2 / 12` |
//!
//! The fixed-fixed (clamped-clamped) maximum moment is the end moment;
//! the peak bending stress is `sigma = M c / I` in every case, with `c`
//! the section's extreme-fibre distance.

use crate::error::BeamError;
use crate::section::Section;
use serde::{Deserialize, Serialize};

/// Which way the two ends of the beam are held.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Support {
    /// Built in (fixed) at one end, free at the other.
    Cantilever,
    /// Pin / roller supported at both ends (statically determinate).
    SimplySupported,
    /// Built in (fixed / clamped) at *both* ends — statically
    /// indeterminate, and far stiffer than the simply-supported case
    /// (the end moments of fixity carry part of the load).
    FixedFixed,
}

/// The applied load.
///
/// `Point` is applied at the *characteristic* location for the support
/// type (the free tip of a cantilever, or mid-span of a
/// simply-supported or fixed-fixed beam); `Udl` is a uniformly
/// distributed load `w` (force per unit length) along the whole span.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum Load {
    /// Concentrated point load `P` (force).
    Point {
        /// Magnitude `P`, force units.
        force: f64,
    },
    /// Uniformly distributed load `w` (force per unit length).
    Udl {
        /// Intensity `w`, force-per-length units.
        intensity: f64,
    },
}

/// The *kind* of load — concentrated or distributed — without a
/// magnitude.
///
/// Mirrors the two [`Load`] variants but carries no value: it selects
/// point-vs-distributed for the load-capacity inverse
/// [`Beam::allowable_load`], which solves for the magnitude itself.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum LoadKind {
    /// A concentrated point load at the characteristic location
    /// (mirrors [`Load::Point`]).
    Point,
    /// A uniformly distributed load along the whole span (mirrors
    /// [`Load::Udl`]).
    Udl,
}

/// A prismatic, homogeneous, linear-elastic beam.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Beam {
    /// Span `L` between supports (or root-to-tip for a cantilever),
    /// length units.
    pub length: f64,
    /// Young's modulus `E`, force-per-area units.
    pub youngs_modulus: f64,
    /// Cross-section providing `I` and `c`.
    pub section: Section,
}

/// The headline outputs of a beam analysis: peak deflection, peak
/// bending moment and the peak bending stress they imply.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct BeamResult {
    /// Magnitude of the maximum transverse deflection, length units.
    pub max_deflection: f64,
    /// Magnitude of the maximum bending moment, force·length units.
    pub max_moment: f64,
    /// Magnitude of the peak bending stress `sigma = M c / I`,
    /// force-per-area units.
    pub max_stress: f64,
}

impl Beam {
    /// Construct a validated beam.
    ///
    /// # Errors
    /// Returns [`BeamError::BadParameter`] if `length` or
    /// `youngs_modulus` is not finite and strictly positive. The
    /// `section` is assumed already validated (built through
    /// [`Section::rectangular`] / [`Section::circular`]).
    pub fn new(length: f64, youngs_modulus: f64, section: Section) -> Result<Self, BeamError> {
        BeamError::require_positive("length", length)?;
        BeamError::require_positive("youngs_modulus", youngs_modulus)?;
        Ok(Beam {
            length,
            youngs_modulus,
            section,
        })
    }

    /// `E * I`, the flexural rigidity (force·length^2).
    pub fn flexural_rigidity(&self) -> f64 {
        self.youngs_modulus * self.section.second_moment_area()
    }

    /// Analyse this beam for the given `support` and `load`, returning
    /// the peak deflection, moment and stress.
    ///
    /// # Errors
    /// Returns [`BeamError::BadParameter`] if the load magnitude is not
    /// finite and strictly positive, or [`BeamError::DegenerateSection`]
    /// if the section has a vanishing second moment of area (so stress
    /// is undefined).
    pub fn analyze(&self, support: Support, load: Load) -> Result<BeamResult, BeamError> {
        let max_deflection = self.max_deflection(support, load)?;
        let max_moment = self.max_moment(support, load)?;
        let max_stress = self.stress_from_moment(max_moment)?;
        Ok(BeamResult {
            max_deflection,
            max_moment,
            max_stress,
        })
    }

    /// Magnitude of the maximum transverse deflection for `support`
    /// under `load`.
    ///
    /// Closed forms (see the [module docs](crate::beam)):
    /// cantilever point `P L^3 / (3 E I)`, cantilever UDL
    /// `w L^4 / (8 E I)`, simply-supported centre point
    /// `P L^3 / (48 E I)`, simply-supported UDL `5 w L^4 / (384 E I)`,
    /// fixed-fixed centre point `P L^3 / (192 E I)`, fixed-fixed UDL
    /// `w L^4 / (384 E I)`.
    ///
    /// # Errors
    /// Returns [`BeamError::BadParameter`] if the load magnitude is not
    /// finite and strictly positive.
    pub fn max_deflection(&self, support: Support, load: Load) -> Result<f64, BeamError> {
        let ei = self.flexural_rigidity();
        let l = self.length;
        let deflection = match (support, load) {
            (Support::Cantilever, Load::Point { force }) => {
                let p = BeamError::require_positive("force", force)?;
                p * l.powi(3) / (3.0 * ei)
            }
            (Support::Cantilever, Load::Udl { intensity }) => {
                let w = BeamError::require_positive("intensity", intensity)?;
                w * l.powi(4) / (8.0 * ei)
            }
            (Support::SimplySupported, Load::Point { force }) => {
                let p = BeamError::require_positive("force", force)?;
                p * l.powi(3) / (48.0 * ei)
            }
            (Support::SimplySupported, Load::Udl { intensity }) => {
                let w = BeamError::require_positive("intensity", intensity)?;
                5.0 * w * l.powi(4) / (384.0 * ei)
            }
            (Support::FixedFixed, Load::Point { force }) => {
                let p = BeamError::require_positive("force", force)?;
                p * l.powi(3) / (192.0 * ei)
            }
            (Support::FixedFixed, Load::Udl { intensity }) => {
                let w = BeamError::require_positive("intensity", intensity)?;
                w * l.powi(4) / (384.0 * ei)
            }
        };
        Ok(deflection)
    }

    /// Magnitude of the maximum bending moment for `support` under
    /// `load`.
    ///
    /// Closed forms: cantilever point `P L` (at the root), cantilever
    /// UDL `w L^2 / 2` (at the root), simply-supported centre point
    /// `P L / 4` (at mid-span), simply-supported UDL `w L^2 / 8` (at
    /// mid-span), fixed-fixed centre point `P L / 8` (at the ends),
    /// fixed-fixed UDL `w L^2 / 12` (at the ends).
    ///
    /// # Errors
    /// Returns [`BeamError::BadParameter`] if the load magnitude is not
    /// finite and strictly positive.
    pub fn max_moment(&self, support: Support, load: Load) -> Result<f64, BeamError> {
        let l = self.length;
        let moment = match (support, load) {
            (Support::Cantilever, Load::Point { force }) => {
                let p = BeamError::require_positive("force", force)?;
                p * l
            }
            (Support::Cantilever, Load::Udl { intensity }) => {
                let w = BeamError::require_positive("intensity", intensity)?;
                w * l * l / 2.0
            }
            (Support::SimplySupported, Load::Point { force }) => {
                let p = BeamError::require_positive("force", force)?;
                p * l / 4.0
            }
            (Support::SimplySupported, Load::Udl { intensity }) => {
                let w = BeamError::require_positive("intensity", intensity)?;
                w * l * l / 8.0
            }
            (Support::FixedFixed, Load::Point { force }) => {
                // Peak |M| is at the built-in ends, P L / 8 (the centre
                // moment is also P L / 8 but of opposite sign).
                let p = BeamError::require_positive("force", force)?;
                p * l / 8.0
            }
            (Support::FixedFixed, Load::Udl { intensity }) => {
                // Peak |M| is at the built-in ends, w L^2 / 12 (the centre
                // moment is w L^2 / 24).
                let w = BeamError::require_positive("intensity", intensity)?;
                w * l * l / 12.0
            }
        };
        Ok(moment)
    }

    /// Peak bending stress produced by a bending moment of magnitude
    /// `moment`, via the flexure formula `sigma = M c / I`.
    ///
    /// # Errors
    /// Returns [`BeamError::BadParameter`] if `moment` is negative or
    /// non-finite, or [`BeamError::DegenerateSection`] if the section's
    /// second moment of area vanishes.
    pub fn stress_from_moment(&self, moment: f64) -> Result<f64, BeamError> {
        if !moment.is_finite() || moment < 0.0 {
            return Err(BeamError::bad_parameter("moment", moment));
        }
        let i = self.section.second_moment_area();
        if !(i.is_finite() && i > 0.0) {
            return Err(BeamError::DegenerateSection {
                reason: "second moment of area is zero",
            });
        }
        let c = self.section.extreme_fibre();
        Ok(moment * c / i)
    }

    /// The largest load this beam can carry before its peak bending
    /// stress reaches `allowable_stress` — the load-capacity
    /// ("safe load") inverse of [`analyze`](Beam::analyze).
    ///
    /// The most a beam can be bent at the stress limit is the allowable
    /// moment `M_allow = allowable_stress * S`, with `S = I / c` the
    /// elastic [section modulus](Section::section_modulus). Inverting
    /// the moment-from-load relation of [`max_moment`](Beam::max_moment)
    /// then gives the allowable load magnitude:
    ///
    /// | Support          | Point `P_allow`  | Distributed `w_allow` |
    /// |------------------|------------------|-----------------------|
    /// | Cantilever       | `M_allow / L`    | `2 M_allow / L^2`     |
    /// | Simply-supported | `4 M_allow / L`  | `8 M_allow / L^2`     |
    /// | Fixed-fixed      | `8 M_allow / L`  | `12 M_allow / L^2`    |
    ///
    /// The answer is returned as a [`Load`] of the requested `kind`, so
    /// feeding it straight back into [`analyze`](Beam::analyze)
    /// reproduces `allowable_stress` as the peak bending stress.
    ///
    /// # Errors
    /// Returns [`BeamError::BadParameter`] if `allowable_stress` is not
    /// finite and strictly positive, or [`BeamError::DegenerateSection`]
    /// if the section's extreme-fibre distance vanishes (so the section
    /// modulus is undefined).
    pub fn allowable_load(
        &self,
        support: Support,
        kind: LoadKind,
        allowable_stress: f64,
    ) -> Result<Load, BeamError> {
        let sigma = BeamError::require_positive("allowable_stress", allowable_stress)?;
        let m_allow = sigma * self.section.section_modulus()?;
        let l = self.length;
        let load = match (support, kind) {
            (Support::Cantilever, LoadKind::Point) => Load::Point { force: m_allow / l },
            (Support::Cantilever, LoadKind::Udl) => Load::Udl {
                intensity: 2.0 * m_allow / (l * l),
            },
            (Support::SimplySupported, LoadKind::Point) => Load::Point {
                force: 4.0 * m_allow / l,
            },
            (Support::SimplySupported, LoadKind::Udl) => Load::Udl {
                intensity: 8.0 * m_allow / (l * l),
            },
            (Support::FixedFixed, LoadKind::Point) => Load::Point {
                force: 8.0 * m_allow / l,
            },
            (Support::FixedFixed, LoadKind::Udl) => Load::Udl {
                intensity: 12.0 * m_allow / (l * l),
            },
        };
        Ok(load)
    }

    /// The largest load this beam can carry before its peak deflection
    /// reaches `allowable_deflection` — the serviceability (stiffness)
    /// counterpart of [`allowable_load`](Beam::allowable_load), inverting
    /// [`max_deflection`](Beam::max_deflection) instead of the stress.
    ///
    /// Each deflection closed form `delta = k * (load) * L^n / (E I)` is
    /// solved for the load:
    ///
    /// | Support          | Point `P_allow`             | Distributed `w_allow`        |
    /// |------------------|-----------------------------|------------------------------|
    /// | Cantilever       | `3 E I d / L^3`             | `8 E I d / L^4`              |
    /// | Simply-supported | `48 E I d / L^3`            | `384 E I d / (5 L^4)`        |
    /// | Fixed-fixed      | `192 E I d / L^3`           | `384 E I d / L^4`            |
    ///
    /// (`d = allowable_deflection`, `E I` the
    /// [flexural rigidity](Beam::flexural_rigidity)). The answer is a
    /// [`Load`] of the requested `kind`, so feeding it back into
    /// [`max_deflection`](Beam::max_deflection) reproduces
    /// `allowable_deflection`. A real design takes the smaller of this and
    /// [`allowable_load`](Beam::allowable_load) (stress-governed) as the
    /// governing safe load.
    ///
    /// # Errors
    /// Returns [`BeamError::BadParameter`] if `allowable_deflection` is not
    /// finite and strictly positive.
    pub fn deflection_limited_load(
        &self,
        support: Support,
        kind: LoadKind,
        allowable_deflection: f64,
    ) -> Result<Load, BeamError> {
        let d = BeamError::require_positive("allowable_deflection", allowable_deflection)?;
        let ei = self.flexural_rigidity();
        let l = self.length;
        let load = match (support, kind) {
            (Support::Cantilever, LoadKind::Point) => Load::Point {
                force: 3.0 * ei * d / l.powi(3),
            },
            (Support::Cantilever, LoadKind::Udl) => Load::Udl {
                intensity: 8.0 * ei * d / l.powi(4),
            },
            (Support::SimplySupported, LoadKind::Point) => Load::Point {
                force: 48.0 * ei * d / l.powi(3),
            },
            (Support::SimplySupported, LoadKind::Udl) => Load::Udl {
                intensity: 384.0 * ei * d / (5.0 * l.powi(4)),
            },
            (Support::FixedFixed, LoadKind::Point) => Load::Point {
                force: 192.0 * ei * d / l.powi(3),
            },
            (Support::FixedFixed, LoadKind::Udl) => Load::Udl {
                intensity: 384.0 * ei * d / l.powi(4),
            },
        };
        Ok(load)
    }
}
