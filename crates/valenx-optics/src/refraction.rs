//! Refraction at a planar interface: Snell's law, the critical angle,
//! and total internal reflection.
//!
//! ## What
//!
//! A ray crosses the boundary between two homogeneous, isotropic media
//! of refractive indices `n1` (incident side) and `n2` (transmitted
//! side). All angles are measured from the surface normal.
//!
//! ## Model
//!
//! Snell's law of refraction:
//!
//! ```text
//! n1 * sin(t1) = n2 * sin(t2)
//! ```
//!
//! Solving for the refraction angle gives
//! `sin(t2) = (n1 / n2) * sin(t1)`. When `n1 > n2` (a ray leaving a
//! denser medium for a rarer one) the right-hand side can exceed unity;
//! the incidence angle at which `sin(t2) = 1` exactly is the
//! **critical angle**
//!
//! ```text
//! t_c = asin(n2 / n1)   (defined only for n1 > n2)
//! ```
//!
//! At or beyond `t_c` no real refraction angle exists and the ray is
//! **totally internally reflected** — the interface behaves as a perfect
//! mirror.
//!
//! ## Honest scope
//!
//! First-order ray geometry only. This module returns the refracted ray
//! *direction* (an angle); it does not compute how much light is
//! transmitted versus reflected (the Fresnel equations), nor does it
//! model polarisation, absorption, dispersion (`n` is taken as a single
//! number, not a function of wavelength), or any wave / diffraction
//! effect. Research / educational grade — not a production optical-design
//! tool.

use crate::error::{validate_angle_deg, validate_index, OpticsError};

/// A planar interface between two media, oriented by which side the ray
/// starts on.
///
/// `n1` is the refractive index of the medium the incident ray travels
/// *through* before the surface; `n2` is the index of the medium on the
/// far side. Constructing an [`Interface`] validates both indices.
#[derive(Clone, Copy, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct Interface {
    /// Refractive index of the incident medium (ray starts here).
    pub n1: f64,
    /// Refractive index of the transmitting medium (ray ends here).
    pub n2: f64,
}

impl Interface {
    /// Construct an interface from the two refractive indices.
    ///
    /// # Errors
    ///
    /// Returns [`OpticsError::InvalidParameter`] if either index is not
    /// finite and strictly positive.
    ///
    /// # Examples
    ///
    /// ```
    /// use valenx_optics::refraction::Interface;
    /// let air_to_glass = Interface::new(1.0, 1.5).unwrap();
    /// assert_eq!(air_to_glass.n1, 1.0);
    /// ```
    pub fn new(n1: f64, n2: f64) -> Result<Self, OpticsError> {
        let n1 = validate_index("n1", n1)?;
        let n2 = validate_index("n2", n2)?;
        Ok(Self { n1, n2 })
    }

    /// The interface with the two media swapped (ray travelling the
    /// other way through the same boundary).
    ///
    /// Because both indices were validated at construction, the reversed
    /// interface is always valid, so this is infallible.
    ///
    /// # Examples
    ///
    /// ```
    /// use valenx_optics::refraction::Interface;
    /// let f = Interface::new(1.0, 1.5).unwrap();
    /// let r = f.reversed();
    /// assert_eq!((r.n1, r.n2), (1.5, 1.0));
    /// ```
    pub fn reversed(&self) -> Self {
        Self {
            n1: self.n2,
            n2: self.n1,
        }
    }

    /// Whether this interface can produce total internal reflection.
    ///
    /// TIR is possible only when the ray leaves a denser medium for a
    /// rarer one, i.e. `n1 > n2`.
    ///
    /// # Examples
    ///
    /// ```
    /// use valenx_optics::refraction::Interface;
    /// assert!(Interface::new(1.5, 1.0).unwrap().tir_possible());
    /// assert!(!Interface::new(1.0, 1.5).unwrap().tir_possible());
    /// ```
    pub fn tir_possible(&self) -> bool {
        self.n1 > self.n2
    }

    /// The critical angle of this interface, in degrees, or [`None`] if
    /// total internal reflection cannot occur (`n1 <= n2`).
    ///
    /// When defined, the critical angle is `asin(n2 / n1)`.
    ///
    /// # Examples
    ///
    /// ```
    /// use valenx_optics::refraction::Interface;
    /// // Glass (1.5) to air (1.0): critical angle ~= 41.81 degrees.
    /// let c = Interface::new(1.5, 1.0).unwrap().critical_angle_deg().unwrap();
    /// assert!((c - 41.810_314_9).abs() < 1e-6);
    /// // Air to glass: no critical angle.
    /// assert!(Interface::new(1.0, 1.5).unwrap().critical_angle_deg().is_none());
    /// ```
    pub fn critical_angle_deg(&self) -> Option<f64> {
        critical_angle_deg(self.n1, self.n2)
    }

    /// Refract an incidence angle (degrees) through this interface.
    ///
    /// See [`refract_angle`] for the full contract.
    ///
    /// # Errors
    ///
    /// Propagates the errors of [`refract_angle`]: an out-of-range
    /// incidence angle, or [`OpticsError::TotalInternalReflection`] when
    /// the ray is totally internally reflected.
    pub fn refract(&self, incidence_deg: f64) -> Result<f64, OpticsError> {
        refract_angle(self.n1, self.n2, incidence_deg)
    }
}

/// The critical angle for an interface, in degrees, or [`None`] when no
/// such angle exists.
///
/// The critical angle is defined only for a ray passing from a denser
/// medium into a rarer one (`n1 > n2`); in that case it is
/// `asin(n2 / n1)`. For `n1 <= n2` (entering a denser-or-equal medium)
/// total internal reflection is impossible and the result is [`None`].
///
/// This is a free function (no [`Interface`] needed) for convenience;
/// it does **not** validate the indices, so callers passing raw user
/// input should validate first (or go through [`Interface`]).
///
/// # Examples
///
/// ```
/// use valenx_optics::refraction::critical_angle_deg;
/// // Water (1.33) to air (1.0).
/// let c = critical_angle_deg(1.33, 1.0).unwrap();
/// assert!((c - 48.753_466_6).abs() < 1e-6);
/// // Equal indices: there is no interface to reflect off.
/// assert!(critical_angle_deg(1.5, 1.5).is_none());
/// ```
pub fn critical_angle_deg(n1: f64, n2: f64) -> Option<f64> {
    if n1 > n2 {
        // n2 / n1 is in (0, 1) here, so asin is well defined.
        Some((n2 / n1).asin().to_degrees())
    } else {
        None
    }
}

/// Apply Snell's law to a single incidence angle.
///
/// Given the two refractive indices and an incidence angle `t1` measured
/// from the normal (in degrees), returns the refraction angle `t2` (also
/// in degrees), where `n1 * sin(t1) = n2 * sin(t2)`.
///
/// The map is exactly invertible: refracting `t1` from `n1` into `n2`
/// and then refracting the result back from `n2` into `n1` recovers
/// `t1` (when no total internal reflection intervenes). This round-trip
/// invariance is verified in the module tests.
///
/// # Errors
///
/// - [`OpticsError::InvalidParameter`] if either index is non-finite /
///   non-positive, or `incidence_deg` is outside `[0, 90]`.
/// - [`OpticsError::TotalInternalReflection`] if `n1 > n2` and the
///   incidence angle is at or beyond the critical angle, so Snell's law
///   has no real solution.
///
/// # Examples
///
/// ```
/// use valenx_optics::refraction::refract_angle;
/// // Light entering glass bends toward the normal.
/// let t2 = refract_angle(1.0, 1.5, 30.0).unwrap();
/// assert!(t2 < 30.0);
/// assert!((t2 - 19.471_220_6).abs() < 1e-6);
/// ```
pub fn refract_angle(n1: f64, n2: f64, incidence_deg: f64) -> Result<f64, OpticsError> {
    let n1 = validate_index("n1", n1)?;
    let n2 = validate_index("n2", n2)?;
    let incidence_deg = validate_angle_deg("incidence_deg", incidence_deg)?;

    let sin_t2 = (n1 / n2) * incidence_deg.to_radians().sin();

    // sin_t2 > 1 means the angle would be imaginary -> TIR. Use a tiny
    // tolerance so that the exactly-critical angle (sin_t2 == 1 in exact
    // arithmetic but 1 + epsilon after rounding) is still reported as
    // TIR rather than panicking through asin's NaN.
    if sin_t2 > 1.0 + 1e-12 {
        let critical_deg = critical_angle_deg(n1, n2).unwrap_or(f64::NAN);
        return Err(OpticsError::TotalInternalReflection {
            n1,
            n2,
            incidence_deg,
            critical_deg,
        });
    }

    // Clamp to guard against asin(NaN) from the 1 + epsilon rounding band.
    let sin_t2 = sin_t2.min(1.0);
    Ok(sin_t2.asin().to_degrees())
}

/// The classification of a ray hitting an interface: does it refract, and
/// if so, at what angle?
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum RayOutcome {
    /// The ray refracts; carries the refraction angle in degrees.
    Refracted {
        /// Refraction angle from the normal, in degrees.
        angle_deg: f64,
    },
    /// The ray is totally internally reflected; carries the critical
    /// angle in degrees that it met or exceeded.
    TotallyReflected {
        /// The critical angle for this interface, in degrees.
        critical_deg: f64,
    },
}

/// Classify a ray at an interface into refraction vs. total internal
/// reflection, instead of returning an error for the TIR case.
///
/// This is a convenience wrapper over [`refract_angle`] for callers that
/// treat TIR as an ordinary (non-error) outcome — e.g. a fibre-optics or
/// prism demonstration where reflection is the point.
///
/// # Errors
///
/// Returns [`OpticsError::InvalidParameter`] only for genuinely invalid
/// inputs (bad index or out-of-range angle). The TIR case is reported
/// through [`RayOutcome::TotallyReflected`], not as an error.
///
/// # Examples
///
/// ```
/// use valenx_optics::refraction::{classify_ray, RayOutcome};
/// // Glass to air at 60 deg is past the ~41.8 deg critical angle.
/// match classify_ray(1.5, 1.0, 60.0).unwrap() {
///     RayOutcome::TotallyReflected { critical_deg } => {
///         assert!((critical_deg - 41.810_314_9).abs() < 1e-6);
///     }
///     other => panic!("expected TIR, got {other:?}"),
/// }
/// ```
pub fn classify_ray(n1: f64, n2: f64, incidence_deg: f64) -> Result<RayOutcome, OpticsError> {
    match refract_angle(n1, n2, incidence_deg) {
        Ok(angle_deg) => Ok(RayOutcome::Refracted { angle_deg }),
        Err(OpticsError::TotalInternalReflection { critical_deg, .. }) => {
            Ok(RayOutcome::TotallyReflected { critical_deg })
        }
        Err(other) => Err(other),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::error::ErrorCategory;

    const EPS: f64 = 1e-9;

    #[test]
    fn snell_bends_toward_normal_entering_denser_medium() {
        // n1 < n2: refracted ray is closer to the normal than incident.
        let t2 = refract_angle(1.0, 1.5, 30.0).unwrap();
        assert!(t2 < 30.0);
        // asin(sin(30 deg)/1.5) = asin(0.5/1.5) = asin(1/3).
        let expected = (1.0_f64 / 3.0).asin().to_degrees();
        assert!(
            (t2 - expected).abs() < EPS,
            "t2 = {t2}, expected {expected}"
        );
    }

    #[test]
    fn snell_bends_away_from_normal_entering_rarer_medium() {
        // n1 > n2 but below critical: ray bends away from the normal.
        let t2 = refract_angle(1.5, 1.0, 20.0).unwrap();
        assert!(t2 > 20.0);
        let expected = (1.5 * 20.0_f64.to_radians().sin()).asin().to_degrees();
        assert!(
            (t2 - expected).abs() < EPS,
            "t2 = {t2}, expected {expected}"
        );
    }

    #[test]
    fn normal_incidence_passes_straight_through() {
        // 0 deg in -> 0 deg out, for any index pair.
        let t2 = refract_angle(1.0, 2.4, 0.0).unwrap();
        assert!(t2.abs() < EPS, "t2 = {t2}");
    }

    #[test]
    fn snell_is_invertible_round_trip() {
        // Refract into the medium, then back: recover the original angle.
        for &t1 in &[5.0_f64, 12.5, 30.0, 42.0, 75.0] {
            let into = refract_angle(1.0, 1.5, t1).unwrap();
            let back = refract_angle(1.5, 1.0, into).unwrap();
            assert!(
                (back - t1).abs() < 1e-7,
                "round trip failed: t1 = {t1}, into = {into}, back = {back}"
            );
        }
    }

    #[test]
    fn critical_angle_only_when_n1_greater_than_n2() {
        // n1 > n2: defined.
        assert!(critical_angle_deg(1.5, 1.0).is_some());
        // n1 < n2: undefined.
        assert!(critical_angle_deg(1.0, 1.5).is_none());
        // n1 == n2: undefined (no real interface to reflect from).
        assert!(critical_angle_deg(1.4, 1.4).is_none());
    }

    #[test]
    fn critical_angle_equals_asin_n2_over_n1() {
        // Glass (1.5) -> air (1.0): well-known ~41.81 deg.
        let c = critical_angle_deg(1.5, 1.0).unwrap();
        let expected = (1.0_f64 / 1.5).asin().to_degrees();
        assert!((c - expected).abs() < EPS, "c = {c}, expected {expected}");
        assert!((c - 41.810_314_895_778_6).abs() < 1e-9, "c = {c}");
    }

    #[test]
    fn at_critical_angle_refraction_is_ninety_degrees() {
        // As the incidence approaches the critical angle from below, the
        // refracted ray approaches grazing (90 deg). asin has a vertical
        // tangent at 1, so the approach is slow; verify both that each
        // computed angle matches the exact Snell value and that smaller
        // gaps below critical push the refraction angle monotonically
        // closer to 90 deg.
        let (n1, n2) = (1.5_f64, 1.0_f64);
        let c = critical_angle_deg(n1, n2).unwrap();
        let mut prev = 0.0_f64;
        // Each successive gap is smaller, but note that asin's vertical
        // tangent at 1 means even a 1e-6 deg gap leaves t2 ~0.01 deg shy
        // of 90: the approach is real but slow, and f64 angle resolution
        // near ~41.8 deg caps how close we can numerically push it.
        for &gap in &[1e-2_f64, 1e-4, 1e-6] {
            let inc = c - gap;
            let t2 = refract_angle(n1, n2, inc).unwrap();
            // Exact closed-form: t2 = asin((n1/n2) sin(inc)).
            let expected = ((n1 / n2) * inc.to_radians().sin()).asin().to_degrees();
            assert!(
                (t2 - expected).abs() < EPS,
                "gap = {gap}: t2 = {t2}, expected {expected}"
            );
            assert!(t2 > prev, "gap = {gap}: t2 = {t2} not closer to grazing");
            assert!(t2 < 90.0, "gap = {gap}: t2 = {t2} must stay below 90");
            prev = t2;
        }
        // The closest gap is within a hundredth of a degree of grazing.
        assert!(prev > 89.98, "closest approach t2 = {prev}");
    }

    #[test]
    fn beyond_critical_angle_is_total_internal_reflection() {
        let err = refract_angle(1.5, 1.0, 60.0).unwrap_err();
        match err {
            OpticsError::TotalInternalReflection {
                n1,
                n2,
                incidence_deg,
                critical_deg,
            } => {
                assert!((n1 - 1.5).abs() < EPS);
                assert!((n2 - 1.0).abs() < EPS);
                assert!((incidence_deg - 60.0).abs() < EPS);
                assert!((critical_deg - 41.810_314_895_778_6).abs() < 1e-9);
            }
            other => panic!("expected TIR, got {other:?}"),
        }
        assert_eq!(err.code(), "optics.total_internal_reflection");
        assert_eq!(err.category(), ErrorCategory::Physical);
    }

    #[test]
    fn classify_ray_reports_refraction_and_tir() {
        match classify_ray(1.5, 1.0, 20.0).unwrap() {
            RayOutcome::Refracted { angle_deg } => assert!(angle_deg > 20.0),
            other => panic!("expected refraction, got {other:?}"),
        }
        match classify_ray(1.5, 1.0, 60.0).unwrap() {
            RayOutcome::TotallyReflected { critical_deg } => {
                assert!((critical_deg - 41.810_314_895_778_6).abs() < 1e-9);
            }
            other => panic!("expected TIR, got {other:?}"),
        }
    }

    #[test]
    fn interface_refract_matches_free_function() {
        let iface = Interface::new(1.0, 1.33).unwrap();
        let a = iface.refract(35.0).unwrap();
        let b = refract_angle(1.0, 1.33, 35.0).unwrap();
        assert!((a - b).abs() < EPS, "a = {a}, b = {b}");
    }

    #[test]
    fn interface_reversed_swaps_indices() {
        let iface = Interface::new(1.0, 1.5).unwrap();
        let rev = iface.reversed();
        assert!((rev.n1 - 1.5).abs() < EPS);
        assert!((rev.n2 - 1.0).abs() < EPS);
        assert!(rev.tir_possible());
        assert!(!iface.tir_possible());
    }

    #[test]
    fn invalid_index_is_rejected() {
        assert!(Interface::new(0.0, 1.5).is_err());
        assert!(Interface::new(1.0, -1.0).is_err());
        assert!(refract_angle(f64::NAN, 1.0, 10.0).is_err());
        let err = refract_angle(1.0, 0.0, 10.0).unwrap_err();
        assert_eq!(err.code(), "optics.invalid_parameter");
        assert_eq!(err.category(), ErrorCategory::Input);
    }

    #[test]
    fn out_of_range_angle_is_rejected() {
        assert!(refract_angle(1.0, 1.5, -1.0).is_err());
        assert!(refract_angle(1.0, 1.5, 90.1).is_err());
        // Exactly the boundary values are accepted.
        assert!(refract_angle(1.0, 1.5, 0.0).is_ok());
        assert!(refract_angle(1.0, 1.5, 90.0).is_ok());
    }
}
