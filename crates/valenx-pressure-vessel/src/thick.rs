//! Thick-wall cylinder — the exact Lame elasticity solution.
//!
//! ## Model
//!
//! For a long cylinder of inner radius `a` and outer radius `b` carrying
//! internal pressure `p_i` and external pressure `p_o`, the axisymmetric
//! plane-strain elasticity equations integrate to the **Lame solution**
//! (Lame, 1852; Timoshenko & Goodier, *Theory of Elasticity*, art. 28):
//! the radial and circumferential (hoop) stresses at a wall radius `r`
//! (with `a <= r <= b`) are
//!
//! ```text
//! sigma_r(r) = A - B / r^2
//! sigma_t(r) = A + B / r^2
//! ```
//!
//! where the two Lame constants follow from the surface pressure boundary
//! conditions `sigma_r(a) = -p_i`, `sigma_r(b) = -p_o`:
//!
//! ```text
//! A = (p_i a^2 - p_o b^2) / (b^2 - a^2)
//! B = (p_i - p_o) a^2 b^2 / (b^2 - a^2)
//! ```
//!
//! Two consequences hold for *any* wall and drive the validation tests:
//!
//! - The sum `sigma_r + sigma_t = 2 A` is **constant through the wall**.
//! - With internal pressure only, the hoop stress is largest at the bore
//!   (`r = a`) and falls monotonically outward, while the radial stress
//!   runs from `-p_i` at the bore to `0` at the outer surface. So the
//!   bore is the critical location.
//!
//! For a **closed-ended** cylinder the end caps add a uniform axial
//! stress `sigma_z = (p_i a^2 - p_o b^2) / (b^2 - a^2) = A`, the mean of
//! the radial and hoop stresses (the third principal stress).
//!
//! As the wall becomes thin (`t = b - a ≪ a`) the bore hoop stress tends
//! to the membrane value `p_i r / t`, so [`crate::cylinder`] is exactly
//! the thin-wall limit of this solution — checked in the tests.
//!
//! ## Honest scope
//!
//! Linear-elastic, homogeneous, isotropic, long (plane-strain) cylinder.
//! No yielding / autofrettage / plastic redistribution, no thermal,
//! residual, shrink-fit or end-effect stresses, no code allowables.
//! Open- vs closed-ended only changes the axial term; the radial and hoop
//! fields above are independent of end condition.

use crate::error::{Result, VesselError};
use serde::{Deserialize, Serialize};

/// A thick-wall cylinder defined by its bore radius, outer radius and the
/// internal / external pressures acting on it.
///
/// Construct with [`ThickCylinder::new`] (internal pressure only, the
/// common case) or [`ThickCylinder::with_external`] (both pressures).
/// Stresses are returned in the same unit as the supplied pressures.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct ThickCylinder {
    r_inner: f64,
    r_outer: f64,
    p_inner: f64,
    p_outer: f64,
}

/// The three principal stresses at one wall radius of a thick cylinder.
///
/// `radial` and `hoop` come from the Lame solution; `axial` is the
/// closed-ended end-cap stress (the through-wall mean) and is `None` for
/// an open-ended cylinder.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct StressState {
    /// Radius at which the state was evaluated.
    pub radius: f64,
    /// Radial stress `sigma_r` (compressive / negative on a
    /// pressurised wall).
    pub radial: f64,
    /// Circumferential (hoop) stress `sigma_t`.
    pub hoop: f64,
}

impl ThickCylinder {
    /// Build a thick cylinder under internal pressure only (`p_o = 0`).
    ///
    /// # Errors
    ///
    /// Returns [`VesselError::BadParameter`] if `r_inner`, `r_outer` or
    /// `pressure` is non-finite, zero or negative, and
    /// [`VesselError::Geometry`] if `r_outer <= r_inner`.
    pub fn new(r_inner: f64, r_outer: f64, pressure: f64) -> Result<Self> {
        Self::with_external(r_inner, r_outer, pressure, 0.0)
    }

    /// Build a thick cylinder with both internal and external pressure.
    ///
    /// The internal pressure must be strictly positive; the external
    /// pressure may be zero or positive (a vessel cannot be sucked on by
    /// a negative absolute pressure). Use this for shrink-fit / external-
    /// jacket problems.
    ///
    /// # Errors
    ///
    /// Returns [`VesselError::BadParameter`] if a radius or the internal
    /// pressure is non-finite / non-positive, or the external pressure is
    /// non-finite / negative; and [`VesselError::Geometry`] if
    /// `r_outer <= r_inner`.
    pub fn with_external(
        r_inner: f64,
        r_outer: f64,
        p_inner: f64,
        p_outer: f64,
    ) -> Result<Self> {
        let r_inner = VesselError::require_positive("r_inner", r_inner)?;
        let r_outer = VesselError::require_positive("r_outer", r_outer)?;
        let p_inner = VesselError::require_positive("p_inner", p_inner)?;
        // External pressure may be exactly zero: validate finiteness and
        // non-negativity without the strictly-positive gate.
        if !p_outer.is_finite() || p_outer < 0.0 {
            return Err(VesselError::BadParameter {
                name: "p_outer",
                reason: "must be a finite, non-negative number",
                value: p_outer,
            });
        }
        if r_outer <= r_inner {
            return Err(VesselError::geometry(format!(
                "r_outer ({r_outer}) must be strictly greater than r_inner ({r_inner})"
            )));
        }
        Ok(Self {
            r_inner,
            r_outer,
            p_inner,
            p_outer,
        })
    }

    /// Bore (inner) radius `a`.
    pub fn r_inner(&self) -> f64 {
        self.r_inner
    }

    /// Outer radius `b`.
    pub fn r_outer(&self) -> f64 {
        self.r_outer
    }

    /// Internal pressure `p_i`.
    pub fn p_inner(&self) -> f64 {
        self.p_inner
    }

    /// External pressure `p_o` (zero unless built with
    /// [`Self::with_external`]).
    pub fn p_outer(&self) -> f64 {
        self.p_outer
    }

    /// Wall thickness `t = b - a`.
    pub fn thickness(&self) -> f64 {
        self.r_outer - self.r_inner
    }

    /// The Lame constant `A = (p_i a^2 - p_o b^2) / (b^2 - a^2)`. Also the
    /// (constant) through-wall mean of the radial and hoop stresses and
    /// the closed-ended axial stress.
    fn lame_a(&self) -> f64 {
        let a2 = self.r_inner * self.r_inner;
        let b2 = self.r_outer * self.r_outer;
        (self.p_inner * a2 - self.p_outer * b2) / (b2 - a2)
    }

    /// The Lame constant `B = (p_i - p_o) a^2 b^2 / (b^2 - a^2)`.
    fn lame_b(&self) -> f64 {
        let a2 = self.r_inner * self.r_inner;
        let b2 = self.r_outer * self.r_outer;
        (self.p_inner - self.p_outer) * a2 * b2 / (b2 - a2)
    }

    /// Radial stress `sigma_r(r) = A - B / r^2`.
    ///
    /// # Errors
    ///
    /// Returns [`VesselError::Geometry`] if `r` lies outside the closed
    /// wall interval `[r_inner, r_outer]`.
    pub fn radial_stress(&self, r: f64) -> Result<f64> {
        self.check_radius(r)?;
        Ok(self.lame_a() - self.lame_b() / (r * r))
    }

    /// Circumferential (hoop) stress `sigma_t(r) = A + B / r^2`.
    ///
    /// # Errors
    ///
    /// Returns [`VesselError::Geometry`] if `r` lies outside
    /// `[r_inner, r_outer]`.
    pub fn hoop_stress(&self, r: f64) -> Result<f64> {
        self.check_radius(r)?;
        Ok(self.lame_a() + self.lame_b() / (r * r))
    }

    /// Both Lame stresses at radius `r`, bundled into a [`StressState`].
    ///
    /// # Errors
    ///
    /// Returns [`VesselError::Geometry`] if `r` lies outside
    /// `[r_inner, r_outer]`.
    pub fn stress_at(&self, r: f64) -> Result<StressState> {
        self.check_radius(r)?;
        let b_over_r2 = self.lame_b() / (r * r);
        let a = self.lame_a();
        Ok(StressState {
            radius: r,
            radial: a - b_over_r2,
            hoop: a + b_over_r2,
        })
    }

    /// Maximum in-plane shear stress at radius `r`:
    /// `(sigma_t - sigma_r) / 2 = B / r^2`.
    ///
    /// Largest at the bore (`r = a`), where the hoop stress also peaks.
    ///
    /// # Errors
    ///
    /// Returns [`VesselError::Geometry`] if `r` lies outside
    /// `[r_inner, r_outer]`.
    pub fn max_shear_at(&self, r: f64) -> Result<f64> {
        self.check_radius(r)?;
        Ok(self.lame_b() / (r * r))
    }

    /// Peak hoop stress over the whole wall, which for a vessel under net
    /// internal pressure occurs at the bore:
    /// `sigma_t(a) = (p_i (b^2 + a^2) - 2 p_o b^2) / (b^2 - a^2)`.
    ///
    /// With external pressure zero this is the familiar
    /// `p_i (b^2 + a^2) / (b^2 - a^2)`.
    pub fn max_hoop_stress(&self) -> f64 {
        // Evaluate at the bore; the radius is in range by construction.
        self.lame_a() + self.lame_b() / (self.r_inner * self.r_inner)
    }

    /// Closed-ended axial stress `sigma_z = A`, uniform through the wall.
    /// This is the third principal stress when the cylinder has end caps;
    /// for an open-ended cylinder the axial stress is zero instead.
    pub fn axial_stress_closed(&self) -> f64 {
        self.lame_a()
    }

    /// Validate that `r` lies within `[r_inner, r_outer]`.
    fn check_radius(&self, r: f64) -> Result<()> {
        if !r.is_finite() {
            return Err(VesselError::BadParameter {
                name: "r",
                reason: "must be a finite number",
                value: r,
            });
        }
        if r < self.r_inner || r > self.r_outer {
            return Err(VesselError::geometry(format!(
                "evaluation radius {r} outside wall [{a}, {b}]",
                a = self.r_inner,
                b = self.r_outer
            )));
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cylinder::ThinCylinder;

    /// Loose absolute tolerance for f64 stress comparisons.
    const EPS: f64 = 1e-9;

    #[test]
    fn rejects_bad_geometry_and_inputs() {
        // Outer not greater than inner.
        let err = ThickCylinder::new(2.0, 2.0, 10.0).unwrap_err();
        assert_eq!(err.category(), crate::ErrorCategory::Geometry);
        assert!(ThickCylinder::new(3.0, 2.0, 10.0).is_err());
        // Non-positive radius / pressure.
        assert!(ThickCylinder::new(0.0, 2.0, 10.0).is_err());
        assert!(ThickCylinder::new(1.0, 2.0, 0.0).is_err());
        // Negative external pressure is rejected; zero is fine.
        assert!(ThickCylinder::with_external(1.0, 2.0, 10.0, -1.0).is_err());
        assert!(ThickCylinder::with_external(1.0, 2.0, 10.0, 0.0).is_ok());
    }

    #[test]
    fn radius_out_of_range_errors() {
        let c = ThickCylinder::new(1.0, 2.0, 10.0).unwrap();
        assert!(c.radial_stress(0.99).is_err());
        assert!(c.radial_stress(2.01).is_err());
        assert!(c.hoop_stress(f64::NAN).is_err());
        // Endpoints inclusive.
        assert!(c.hoop_stress(1.0).is_ok());
        assert!(c.hoop_stress(2.0).is_ok());
    }

    #[test]
    fn boundary_conditions_on_radial_stress() {
        // Internal pressure only: sigma_r(a) = -p_i, sigma_r(b) = 0.
        let p = 50.0;
        let c = ThickCylinder::new(1.0, 2.0, p).unwrap();
        assert!((c.radial_stress(1.0).unwrap() + p).abs() < 1e-7);
        assert!(c.radial_stress(2.0).unwrap().abs() < 1e-7);
    }

    #[test]
    fn lame_known_values_internal_pressure() {
        // a=1, b=2, p_i=100, p_o=0. b^2-a^2 = 3.
        // sigma_t(a) = p_i (b^2 + a^2)/(b^2 - a^2) = 100 * 5 / 3.
        // sigma_t(b) = 2 p_i a^2 /(b^2 - a^2)      = 200 / 3.
        let c = ThickCylinder::new(1.0, 2.0, 100.0).unwrap();
        assert!((c.hoop_stress(1.0).unwrap() - 500.0 / 3.0).abs() < 1e-6);
        assert!((c.hoop_stress(2.0).unwrap() - 200.0 / 3.0).abs() < 1e-6);
        // Radial at the bore is -p_i.
        assert!((c.radial_stress(1.0).unwrap() + 100.0).abs() < 1e-6);
    }

    #[test]
    fn hoop_is_maximum_at_inner_radius() {
        // The defining validation: hoop stress peaks at the bore and
        // decreases monotonically outward.
        let c = ThickCylinder::new(1.0, 3.0, 80.0).unwrap();
        let at_inner = c.hoop_stress(1.0).unwrap();
        assert!((at_inner - c.max_hoop_stress()).abs() < EPS);
        let mut prev = at_inner;
        for k in 1..=20 {
            let r = 1.0 + (3.0 - 1.0) * (k as f64) / 20.0;
            let h = c.hoop_stress(r).unwrap();
            assert!(h <= prev + EPS, "hoop not monotone decreasing at r={r}");
            assert!(h <= at_inner + EPS, "found hoop above the bore value");
            prev = h;
        }
    }

    #[test]
    fn stress_sum_is_constant_through_wall() {
        // sigma_r + sigma_t = 2A everywhere -> independent of r.
        let c = ThickCylinder::new(1.0, 2.5, 60.0).unwrap();
        let reference =
            c.radial_stress(1.0).unwrap() + c.hoop_stress(1.0).unwrap();
        for k in 0..=10 {
            let r = 1.0 + 1.5 * (k as f64) / 10.0;
            let s = c.stress_at(r).unwrap();
            assert!(
                ((s.radial + s.hoop) - reference).abs() < 1e-7,
                "sum drifted at r={r}"
            );
        }
        // ...and the closed-ended axial stress is exactly half that sum.
        assert!(
            (c.axial_stress_closed() - reference / 2.0).abs() < 1e-7,
            "axial != mean of radial and hoop"
        );
    }

    #[test]
    fn max_shear_peaks_at_bore_and_equals_half_stress_difference() {
        let c = ThickCylinder::new(1.0, 2.0, 100.0).unwrap();
        let s = c.stress_at(1.0).unwrap();
        let half_diff = (s.hoop - s.radial) / 2.0;
        assert!((c.max_shear_at(1.0).unwrap() - half_diff).abs() < 1e-7);
        // Decreasing outward: bore value is the largest.
        assert!(c.max_shear_at(1.0).unwrap() > c.max_shear_at(2.0).unwrap());
    }

    #[test]
    fn stress_scales_linearly_with_pressure() {
        let a = ThickCylinder::new(1.0, 2.0, 30.0).unwrap();
        let b = ThickCylinder::new(1.0, 2.0, 60.0).unwrap();
        assert!(
            (b.hoop_stress(1.0).unwrap() - 2.0 * a.hoop_stress(1.0).unwrap())
                .abs()
                < 1e-7
        );
        assert!((b.max_hoop_stress() - 2.0 * a.max_hoop_stress()).abs() < 1e-6);
    }

    #[test]
    fn external_pressure_only_puts_wall_in_compression() {
        // p_i tiny, p_o large: hoop stress should be compressive
        // (negative) across the wall.
        let c = ThickCylinder::with_external(1.0, 2.0, 1e-6, 50.0).unwrap();
        assert!(c.hoop_stress(1.0).unwrap() < 0.0);
        assert!(c.hoop_stress(2.0).unwrap() < 0.0);
        // Radial BC: sigma_r(b) = -p_o.
        assert!((c.radial_stress(2.0).unwrap() + 50.0).abs() < 1e-6);
    }

    #[test]
    fn thick_wall_tends_to_thin_wall_for_small_t() {
        // As t = b - a << a, the bore hoop stress -> membrane p r / t.
        // Use a thin wall: a = 100, t = 0.5 (r/t = 200).
        let a = 100.0;
        let t = 0.5;
        let p = 5.0;
        let thick = ThickCylinder::new(a, a + t, p).unwrap();
        let thin = ThinCylinder::new(p, a, t).unwrap();
        let rel_err =
            (thick.max_hoop_stress() - thin.hoop_stress()).abs() / thin.hoop_stress();
        // Membrane uses the inner radius; the exact Lame bore value
        // differs by O(t/a) ~ 0.25%. Require agreement to within 1%.
        assert!(
            rel_err < 1e-2,
            "thin-wall limit off by {rel_err} (thick={}, thin={})",
            thick.max_hoop_stress(),
            thin.hoop_stress()
        );
    }
}
