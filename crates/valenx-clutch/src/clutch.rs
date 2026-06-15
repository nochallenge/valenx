//! Disc / plate friction-clutch torque-capacity models.
//!
//! ## The two classic theories
//!
//! A plate clutch transmits torque by pressing two (or more) annular
//! friction faces together with an axial clamping force `F`. How that
//! force distributes over the face — and therefore how much torque the
//! clutch can carry before it slips — depends on which idealisation you
//! adopt:
//!
//! ### Uniform wear
//!
//! Worn-in clutches wear fastest where the rubbing speed is highest, at
//! the outer edge. Archard-style wear is proportional to `p * r` (local
//! pressure times radius times sliding velocity), so a fully worn face
//! settles to the condition `p * r = constant`: the pressure is highest
//! at the inner edge and falls as `1/r`. Integrating the friction
//! moment over the annulus gives the **mean-radius** result
//!
//! ```text
//! T = mu * F * N * (ro + ri) / 2
//! ```
//!
//! i.e. the whole clamp force acts at the arithmetic mean radius
//! `(ro + ri) / 2`. This is the conservative, design-grade number and
//! is what a clutch delivers after break-in.
//!
//! ### Uniform pressure
//!
//! A brand-new, perfectly flat, rigid face carries a *uniform* contact
//! pressure `p = constant`. Integrating the friction moment over the
//! annulus then gives
//!
//! ```text
//! T = (2/3) * mu * F * N * (ro^3 - ri^3) / (ro^2 - ri^2)
//! ```
//!
//! the clamp force effectively acting at the larger **area-weighted**
//! (centroidal) mean radius `(2/3)(ro^3 - ri^3)/(ro^2 - ri^2)`. Because
//! that lever arm is larger than the arithmetic mean, the
//! uniform-pressure theory always predicts a *higher* torque than
//! uniform wear for the same geometry, friction, and clamp force.
//!
//! In both expressions `N` is the number of friction *surfaces* in
//! contact: a single-plate clutch grips on both faces of the disc, so
//! `N = 2`; a multi-plate clutch with `m` driven plates has `N = 2 m`.
//!
//! ## Transmissible power
//!
//! Once the clutch is fully engaged and the two sides spin together at
//! angular velocity `omega` (rad/s), the power it can transmit is simply
//!
//! ```text
//! P = T * omega
//! ```
//!
//! ## Honest scope
//!
//! These are the **textbook rigid-body, steady-state** capacity
//! formulas (Shigley, Juvinall, Norton). They give the *limiting*
//! torque a dry friction clutch can carry without slipping; they do
//! **not** model engagement transients, slip heating / thermal fade,
//! the friction-coefficient drop with temperature or speed, plate
//! flatness and pressure-distribution departures, lubrication of a wet
//! clutch, hysteresis, or fatigue. Use them for first-order sizing and
//! teaching — not as a substitute for a validated clutch test rig or a
//! production design code.

use nalgebra::DVector;
use serde::{Deserialize, Serialize};

use crate::error::ClutchError;

/// The friction face's annular geometry: an inner and outer radius.
///
/// The contact region is the flat ring `ri <= r <= ro`. Radii are held
/// internally in **metres** (SI) so the torque formulas are
/// unit-consistent, but the ergonomic constructor [`Self::new`] takes
/// the more common engineering unit of **millimetres**.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ClutchGeometry {
    /// Inner radius of the friction annulus, in metres.
    inner_m: f64,
    /// Outer radius of the friction annulus, in metres.
    outer_m: f64,
}

impl ClutchGeometry {
    /// Build a friction-face geometry from inner and outer radii given
    /// in **millimetres**.
    ///
    /// # Errors
    ///
    /// Returns [`ClutchError::InvalidParameter`] if either radius is not
    /// finite or is non-positive, and [`ClutchError::InvertedRadii`] if
    /// the inner radius is not strictly less than the outer radius.
    ///
    /// ```
    /// use valenx_clutch::clutch::ClutchGeometry;
    /// let g = ClutchGeometry::new(50.0, 100.0).unwrap();
    /// assert!((g.inner_radius_m() - 0.050).abs() < 1e-12);
    /// assert!((g.outer_radius_m() - 0.100).abs() < 1e-12);
    /// ```
    pub fn new(inner_mm: f64, outer_mm: f64) -> Result<Self, ClutchError> {
        check_positive_finite("inner_mm", inner_mm)?;
        check_positive_finite("outer_mm", outer_mm)?;
        if inner_mm >= outer_mm {
            return Err(ClutchError::InvertedRadii { inner_mm, outer_mm });
        }
        Ok(Self {
            inner_m: inner_mm * 1.0e-3,
            outer_m: outer_mm * 1.0e-3,
        })
    }

    /// Build a geometry directly from radii in **metres**.
    ///
    /// # Errors
    ///
    /// Same validation as [`Self::new`], reported with the millimetre
    /// equivalents so the `Display` text stays in engineering units.
    pub fn from_metres(inner_m: f64, outer_m: f64) -> Result<Self, ClutchError> {
        check_positive_finite("inner_m", inner_m)?;
        check_positive_finite("outer_m", outer_m)?;
        if inner_m >= outer_m {
            return Err(ClutchError::InvertedRadii {
                inner_mm: inner_m * 1.0e3,
                outer_mm: outer_m * 1.0e3,
            });
        }
        Ok(Self { inner_m, outer_m })
    }

    /// Inner radius `ri` of the friction annulus, in metres.
    pub fn inner_radius_m(&self) -> f64 {
        self.inner_m
    }

    /// Outer radius `ro` of the friction annulus, in metres.
    pub fn outer_radius_m(&self) -> f64 {
        self.outer_m
    }

    /// Effective lever-arm radius (metres) under the **uniform-wear**
    /// assumption: the arithmetic mean `(ro + ri) / 2`.
    ///
    /// Multiplying this by `mu * F * N` yields the uniform-wear torque.
    pub fn mean_radius_uniform_wear_m(&self) -> f64 {
        0.5 * (self.outer_m + self.inner_m)
    }

    /// Effective lever-arm radius (metres) under the **uniform-pressure**
    /// assumption: the area-weighted (centroidal) mean
    /// `(2/3)(ro^3 - ri^3)/(ro^2 - ri^2)`.
    ///
    /// Multiplying this by `mu * F * N` yields the uniform-pressure
    /// torque. It is always strictly larger than the uniform-wear mean
    /// radius for a non-degenerate annulus.
    pub fn mean_radius_uniform_pressure_m(&self) -> f64 {
        let ro = self.outer_m;
        let ri = self.inner_m;
        (2.0 / 3.0) * (ro * ro * ro - ri * ri * ri) / (ro * ro - ri * ri)
    }
}

/// Which idealised contact-pressure distribution to assume when sizing
/// the clutch.
///
/// See the [module documentation](crate::clutch) for the derivation of
/// each. As a rule of thumb: design to [`PressureModel::UniformWear`]
/// (the worn-in, conservative case) and treat
/// [`PressureModel::UniformPressure`] (the new-plate, optimistic case)
/// as the upper bound.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum PressureModel {
    /// Worn-in face, `p * r = const`. Torque acts at the arithmetic mean
    /// radius. Conservative — the standard design assumption.
    UniformWear,
    /// New, rigid, flat face, `p = const`. Torque acts at the larger
    /// area-weighted mean radius. Optimistic upper bound.
    UniformPressure,
}

/// A fully-specified dry friction (plate / disc) clutch.
///
/// Bundles the friction-face [`ClutchGeometry`], the coefficient of
/// friction `mu`, and the number of friction surfaces in contact `N`.
/// Construct with [`Self::new`] (validated), then ask it for its torque
/// capacity under either [`PressureModel`] and the power it can transmit
/// at a given speed.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct FrictionClutch {
    /// Friction-face annular geometry.
    geometry: ClutchGeometry,
    /// Coefficient of (kinetic) friction between the rubbing faces.
    mu: f64,
    /// Number of friction surfaces in contact (single plate => 2).
    surfaces: u32,
}

impl FrictionClutch {
    /// Assemble and validate a friction clutch.
    ///
    /// `geometry` is the friction face, `mu` the coefficient of friction
    /// (dimensionless, must be finite and positive), and `surfaces` the
    /// number of friction faces in contact `N` (must be at least one —
    /// a single-plate clutch has `N = 2`).
    ///
    /// # Errors
    ///
    /// Returns [`ClutchError::InvalidParameter`] if `mu` is not finite or
    /// not positive, and [`ClutchError::InvalidSurfaceCount`] if
    /// `surfaces` is zero.
    ///
    /// ```
    /// use valenx_clutch::clutch::{ClutchGeometry, FrictionClutch};
    /// let g = ClutchGeometry::new(50.0, 100.0).unwrap();
    /// let c = FrictionClutch::new(g, 0.3, 2).unwrap();
    /// assert_eq!(c.surfaces(), 2);
    /// ```
    pub fn new(geometry: ClutchGeometry, mu: f64, surfaces: u32) -> Result<Self, ClutchError> {
        check_positive_finite("mu", mu)?;
        if surfaces == 0 {
            return Err(ClutchError::InvalidSurfaceCount(0.0));
        }
        Ok(Self {
            geometry,
            mu,
            surfaces,
        })
    }

    /// The friction-face geometry.
    pub fn geometry(&self) -> ClutchGeometry {
        self.geometry
    }

    /// Coefficient of friction `mu`.
    pub fn mu(&self) -> f64 {
        self.mu
    }

    /// Number of friction surfaces in contact `N`.
    pub fn surfaces(&self) -> u32 {
        self.surfaces
    }

    /// Torque capacity (newton-metres) under the **uniform-wear**
    /// assumption:
    ///
    /// ```text
    /// T = mu * F * N * (ro + ri) / 2
    /// ```
    ///
    /// where `F` is the axial clamp force in newtons.
    ///
    /// # Errors
    ///
    /// Returns [`ClutchError::InvalidParameter`] if `clamp_force_n` is
    /// not finite or is negative.
    pub fn torque_uniform_wear(&self, clamp_force_n: f64) -> Result<f64, ClutchError> {
        check_nonneg_finite("clamp_force_n", clamp_force_n)?;
        Ok(self.mu
            * clamp_force_n
            * self.surfaces as f64
            * self.geometry.mean_radius_uniform_wear_m())
    }

    /// Torque capacity (newton-metres) under the **uniform-pressure**
    /// assumption:
    ///
    /// ```text
    /// T = (2/3) * mu * F * N * (ro^3 - ri^3) / (ro^2 - ri^2)
    /// ```
    ///
    /// where `F` is the axial clamp force in newtons.
    ///
    /// # Errors
    ///
    /// Returns [`ClutchError::InvalidParameter`] if `clamp_force_n` is
    /// not finite or is negative.
    pub fn torque_uniform_pressure(&self, clamp_force_n: f64) -> Result<f64, ClutchError> {
        check_nonneg_finite("clamp_force_n", clamp_force_n)?;
        Ok(self.mu
            * clamp_force_n
            * self.surfaces as f64
            * self.geometry.mean_radius_uniform_pressure_m())
    }

    /// Torque capacity (newton-metres) under the chosen [`PressureModel`].
    ///
    /// A thin dispatcher over [`Self::torque_uniform_wear`] and
    /// [`Self::torque_uniform_pressure`].
    ///
    /// # Errors
    ///
    /// Propagates the validation error of the underlying model.
    pub fn torque(&self, model: PressureModel, clamp_force_n: f64) -> Result<f64, ClutchError> {
        match model {
            PressureModel::UniformWear => self.torque_uniform_wear(clamp_force_n),
            PressureModel::UniformPressure => self.torque_uniform_pressure(clamp_force_n),
        }
    }

    /// Transmissible power (watts) once fully engaged and spinning at
    /// `omega` radians per second under the chosen [`PressureModel`]:
    ///
    /// ```text
    /// P = T * omega
    /// ```
    ///
    /// # Errors
    ///
    /// Returns [`ClutchError::InvalidParameter`] if `omega_rad_per_s` is
    /// not finite or is negative, and propagates any torque-model
    /// validation error.
    pub fn power(
        &self,
        model: PressureModel,
        clamp_force_n: f64,
        omega_rad_per_s: f64,
    ) -> Result<f64, ClutchError> {
        check_nonneg_finite("omega_rad_per_s", omega_rad_per_s)?;
        let torque = self.torque(model, clamp_force_n)?;
        Ok(torque * omega_rad_per_s)
    }

    /// Sample the torque capacity at each of the supplied clamp forces
    /// (newtons) under the chosen [`PressureModel`], returning the
    /// torques (newton-metres) as an `nalgebra` [`DVector`].
    ///
    /// Convenient for building a torque-vs-clamp-force curve to plot or
    /// to feed downstream linear-algebra code. Because the underlying
    /// model is exactly linear in the clamp force, the returned vector is
    /// a straight line through the origin with slope
    /// `mu * N * r_eff`.
    ///
    /// # Errors
    ///
    /// Returns [`ClutchError::InvalidParameter`] for the first clamp
    /// force in `clamp_forces_n` that is not finite or is negative.
    ///
    /// ```
    /// use valenx_clutch::clutch::{ClutchGeometry, FrictionClutch, PressureModel};
    /// let g = ClutchGeometry::new(50.0, 100.0).unwrap();
    /// let c = FrictionClutch::new(g, 0.3, 2).unwrap();
    /// let curve = c
    ///     .torque_curve(PressureModel::UniformWear, &[0.0, 1000.0, 2000.0])
    ///     .unwrap();
    /// // Linear in clamp force: the 2 kN sample is twice the 1 kN sample.
    /// assert!((curve[2] - 2.0 * curve[1]).abs() < 1e-9);
    /// assert!(curve[0].abs() < 1e-12);
    /// ```
    pub fn torque_curve(
        &self,
        model: PressureModel,
        clamp_forces_n: &[f64],
    ) -> Result<DVector<f64>, ClutchError> {
        let mut out = DVector::zeros(clamp_forces_n.len());
        for (slot, &force) in out.iter_mut().zip(clamp_forces_n.iter()) {
            *slot = self.torque(model, force)?;
        }
        Ok(out)
    }
}

/// Convert a rotational speed in revolutions per minute to angular
/// velocity in radians per second: `omega = 2 * pi * rpm / 60`.
///
/// A small convenience so callers can work in the catalogue unit (rpm)
/// and still feed [`FrictionClutch::power`], which expects rad/s.
///
/// # Errors
///
/// Returns [`ClutchError::InvalidParameter`] if `rpm` is not finite or
/// is negative.
///
/// ```
/// use valenx_clutch::clutch::rpm_to_rad_per_s;
/// let w = rpm_to_rad_per_s(60.0).unwrap();
/// assert!((w - std::f64::consts::TAU).abs() < 1e-12);
/// ```
pub fn rpm_to_rad_per_s(rpm: f64) -> Result<f64, ClutchError> {
    check_nonneg_finite("rpm", rpm)?;
    Ok(rpm * std::f64::consts::TAU / 60.0)
}

/// Validate that a parameter is finite and strictly positive.
fn check_positive_finite(name: &'static str, value: f64) -> Result<(), ClutchError> {
    if !value.is_finite() {
        return Err(ClutchError::InvalidParameter {
            name,
            value,
            reason: "must be a finite number",
        });
    }
    if value <= 0.0 {
        return Err(ClutchError::InvalidParameter {
            name,
            value,
            reason: "must be strictly positive",
        });
    }
    Ok(())
}

/// Validate that a parameter is finite and non-negative.
fn check_nonneg_finite(name: &'static str, value: f64) -> Result<(), ClutchError> {
    if !value.is_finite() {
        return Err(ClutchError::InvalidParameter {
            name,
            value,
            reason: "must be a finite number",
        });
    }
    if value < 0.0 {
        return Err(ClutchError::InvalidParameter {
            name,
            value,
            reason: "must be non-negative",
        });
    }
    Ok(())
}
