//! Pneumatic cylinder force from gauge pressure and piston area.
//!
//! ## Model
//!
//! The theoretical (stall) thrust of a pneumatic cylinder is the gauge
//! supply pressure acting over the effective piston area:
//!
//! ```text
//! F = p_gauge * A
//! ```
//!
//! The piston of a single-acting or extending double-acting cylinder
//! presents the full bore area `A_bore = pi/4 * d_bore^2`. On the retract
//! stroke of a double-acting cylinder the rod occupies part of the piston
//! face, so the air acts over the annular area
//! `A_annulus = pi/4 * (d_bore^2 - d_rod^2)`, giving a smaller retract
//! force for the same pressure.
//!
//! Gauge pressure is used deliberately: the atmosphere already pushes on
//! the opposite (or rod-side, vented) face, so only the pressure *above*
//! atmospheric does net work. SI units throughout — metres for diameters
//! and stroke, pascals for pressure, square metres for area, newtons for
//! force.
//!
//! ## Honest scope
//!
//! This is the ideal *static* thrust. A real cylinder delivers less: seal
//! and rod friction (typically a few percent of theoretical force),
//! breakaway / running stiction, back-pressure on the exhaust side, and
//! dynamic effects during motion all reduce the usable force. None of
//! those are modelled here. Treat the result as an upper bound for
//! first-pass sizing, not a guaranteed delivered force.

use crate::error::{PneumaticsError, Result};

/// Mathematical constant pi, to f64 precision.
const PI: f64 = std::f64::consts::PI;

/// The two strokes of a cylinder, which present different effective
/// areas on a double-acting cylinder because the rod occupies part of the
/// piston face on retract.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash, serde::Serialize, serde::Deserialize)]
pub enum Stroke {
    /// Extend stroke: air acts over the full bore (piston) area.
    Extend,
    /// Retract stroke: air acts over the annular area (bore minus rod).
    Retract,
}

/// A pneumatic cylinder defined by its bore and rod diameters.
///
/// Construct with [`Cylinder::single_acting`] (no rod term; both strokes
/// use the full bore area) or [`Cylinder::double_acting`] (the retract
/// stroke subtracts the rod area). All diameters are in metres.
#[derive(Copy, Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct Cylinder {
    /// Bore (piston) diameter, metres. Strictly positive.
    bore: f64,
    /// Rod diameter, metres. Zero for a single-acting cylinder; otherwise
    /// strictly less than `bore`.
    rod: f64,
}

impl Cylinder {
    /// A single-acting cylinder of the given `bore` diameter (metres).
    ///
    /// The rod area is treated as zero, so both [`Stroke`]s use the full
    /// bore area. (A real single-acting cylinder only powers one stroke
    /// and is spring-returned, but its powered-stroke area is the full
    /// bore, which is what this models.)
    ///
    /// # Errors
    ///
    /// [`PneumaticsError::NonPositive`] if `bore` is not finite and `> 0`.
    ///
    /// # Examples
    ///
    /// ```
    /// use valenx_pneumatics::cylinder::Cylinder;
    /// let c = Cylinder::single_acting(0.05).unwrap();
    /// assert!((c.bore_area() - std::f64::consts::PI * 0.025 * 0.025).abs() < 1e-12);
    /// ```
    pub fn single_acting(bore: f64) -> Result<Self> {
        let bore = PneumaticsError::positive("bore", bore)?;
        Ok(Cylinder { bore, rod: 0.0 })
    }

    /// A double-acting cylinder of the given `bore` and `rod` diameters
    /// (metres). The rod must be strictly thinner than the bore.
    ///
    /// # Errors
    ///
    /// - [`PneumaticsError::NonPositive`] if `bore` or `rod` is not finite
    ///   and `> 0`.
    /// - [`PneumaticsError::Geometry`] if `rod >= bore` (which would leave
    ///   a zero or negative annular area).
    ///
    /// # Examples
    ///
    /// ```
    /// use valenx_pneumatics::cylinder::Cylinder;
    /// let c = Cylinder::double_acting(0.05, 0.02).unwrap();
    /// assert!(c.annulus_area() < c.bore_area());
    /// // A rod as wide as the bore is rejected.
    /// assert!(Cylinder::double_acting(0.05, 0.05).is_err());
    /// ```
    pub fn double_acting(bore: f64, rod: f64) -> Result<Self> {
        let bore = PneumaticsError::positive("bore", bore)?;
        let rod = PneumaticsError::positive("rod", rod)?;
        if rod >= bore {
            return Err(PneumaticsError::Geometry(
                "rod diameter must be strictly less than bore diameter",
            ));
        }
        Ok(Cylinder { bore, rod })
    }

    /// The bore (piston) diameter in metres.
    pub fn bore(&self) -> f64 {
        self.bore
    }

    /// The rod diameter in metres (zero for a single-acting cylinder).
    pub fn rod(&self) -> f64 {
        self.rod
    }

    /// The full bore (piston) area, `pi/4 * d_bore^2`, in square metres.
    ///
    /// # Examples
    ///
    /// ```
    /// use valenx_pneumatics::cylinder::Cylinder;
    /// let c = Cylinder::single_acting(0.1).unwrap();
    /// // d = 0.1 m -> A = pi/4 * 0.01 = 0.007853981...
    /// assert!((c.bore_area() - 0.007_853_981_633_974_483).abs() < 1e-15);
    /// ```
    pub fn bore_area(&self) -> f64 {
        0.25 * PI * self.bore * self.bore
    }

    /// The annular (rod-side) area, `pi/4 * (d_bore^2 - d_rod^2)`, in
    /// square metres. Equals [`bore_area`](Self::bore_area) when the rod
    /// diameter is zero (single-acting).
    pub fn annulus_area(&self) -> f64 {
        0.25 * PI * (self.bore * self.bore - self.rod * self.rod)
    }

    /// The effective piston area for the given [`Stroke`], in square
    /// metres: the bore area on [`Stroke::Extend`] and the annular area on
    /// [`Stroke::Retract`].
    pub fn effective_area(&self, stroke: Stroke) -> f64 {
        match stroke {
            Stroke::Extend => self.bore_area(),
            Stroke::Retract => self.annulus_area(),
        }
    }

    /// Theoretical thrust `F = p_gauge * A` for the given stroke, in
    /// newtons, where `gauge_pressure` is the supply pressure *above*
    /// atmospheric in pascals.
    ///
    /// # Errors
    ///
    /// [`PneumaticsError::Negative`] if `gauge_pressure` is negative or
    /// non-finite. (Zero gauge pressure is allowed and yields zero force.)
    ///
    /// # Examples
    ///
    /// ```
    /// use valenx_pneumatics::cylinder::{Cylinder, Stroke};
    /// // 50 mm bore at 6 bar gauge (600 000 Pa) on extend.
    /// let c = Cylinder::single_acting(0.05).unwrap();
    /// let f = c.force(600_000.0, Stroke::Extend).unwrap();
    /// // A = pi/4 * 0.05^2 = 1.9635e-3 m^2; F = 6e5 * A = 1178.1 N.
    /// assert!((f - 1_178.097_245_096_172).abs() < 1e-6);
    /// ```
    pub fn force(&self, gauge_pressure: f64, stroke: Stroke) -> Result<f64> {
        let p = PneumaticsError::non_negative("gauge_pressure", gauge_pressure)?;
        Ok(p * self.effective_area(stroke))
    }

    /// Theoretical extend thrust, `p_gauge * bore_area`. Convenience for
    /// [`force`](Self::force) with [`Stroke::Extend`].
    ///
    /// # Errors
    ///
    /// [`PneumaticsError::Negative`] if `gauge_pressure` is negative or
    /// non-finite.
    pub fn extend_force(&self, gauge_pressure: f64) -> Result<f64> {
        self.force(gauge_pressure, Stroke::Extend)
    }

    /// Theoretical retract thrust, `p_gauge * annulus_area`. Convenience
    /// for [`force`](Self::force) with [`Stroke::Retract`].
    ///
    /// # Errors
    ///
    /// [`PneumaticsError::Negative`] if `gauge_pressure` is negative or
    /// non-finite.
    pub fn retract_force(&self, gauge_pressure: f64) -> Result<f64> {
        self.force(gauge_pressure, Stroke::Retract)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Absolute-difference tolerance for floating comparisons.
    const EPS: f64 = 1e-9;

    #[test]
    fn bore_area_matches_circle_formula() {
        // d = 40 mm -> A = pi/4 * 0.04^2 = 1.256637061e-3 m^2.
        let c = Cylinder::single_acting(0.04).unwrap();
        let expected = 0.25 * PI * 0.04 * 0.04;
        assert!((c.bore_area() - expected).abs() < EPS);
        assert!((c.bore_area() - 1.256_637_061_435_917_2e-3).abs() < 1e-15);
    }

    #[test]
    fn force_equals_gauge_pressure_times_area() {
        // The defining relation: F = p_gauge * A, checked independently.
        let c = Cylinder::single_acting(0.063).unwrap();
        let p = 700_000.0; // 7 bar gauge
        let f = c.force(p, Stroke::Extend).unwrap();
        let expected = p * (0.25 * PI * 0.063 * 0.063);
        assert!((f - expected).abs() < 1e-6);
    }

    #[test]
    fn force_scales_linearly_with_pressure() {
        // Doubling pressure doubles force (linear in p_gauge).
        let c = Cylinder::single_acting(0.05).unwrap();
        let f1 = c.force(300_000.0, Stroke::Extend).unwrap();
        let f2 = c.force(600_000.0, Stroke::Extend).unwrap();
        assert!((f2 - 2.0 * f1).abs() < 1e-6);
    }

    #[test]
    fn force_scales_with_area_quadratically_in_bore() {
        // Doubling the bore quadruples the area and hence the force.
        let small = Cylinder::single_acting(0.025).unwrap();
        let big = Cylinder::single_acting(0.050).unwrap();
        let p = 600_000.0;
        let fs = small.force(p, Stroke::Extend).unwrap();
        let fb = big.force(p, Stroke::Extend).unwrap();
        assert!((fb - 4.0 * fs).abs() < 1e-6);
    }

    #[test]
    fn retract_force_is_less_than_extend_for_double_acting() {
        // Rod steals area on retract, so retract force < extend force.
        let c = Cylinder::double_acting(0.05, 0.02).unwrap();
        let p = 600_000.0;
        let fe = c.extend_force(p).unwrap();
        let fr = c.retract_force(p).unwrap();
        assert!(fr < fe);

        // Quantitatively: ratio = (bore^2 - rod^2) / bore^2.
        let ratio = (0.05 * 0.05 - 0.02 * 0.02) / (0.05 * 0.05);
        assert!((fr / fe - ratio).abs() < EPS);
    }

    #[test]
    fn single_acting_strokes_have_equal_area() {
        // No rod term -> extend and retract areas coincide.
        let c = Cylinder::single_acting(0.032).unwrap();
        assert!((c.effective_area(Stroke::Extend) - c.effective_area(Stroke::Retract)).abs() < EPS);
    }

    #[test]
    fn zero_gauge_pressure_gives_zero_force() {
        let c = Cylinder::single_acting(0.05).unwrap();
        assert!((c.force(0.0, Stroke::Extend).unwrap() - 0.0).abs() < EPS);
    }

    #[test]
    fn rejects_nonpositive_bore() {
        assert!(Cylinder::single_acting(0.0).is_err());
        assert!(Cylinder::single_acting(-0.05).is_err());
        assert!(Cylinder::single_acting(f64::NAN).is_err());
    }

    #[test]
    fn rejects_rod_at_or_above_bore() {
        assert!(Cylinder::double_acting(0.05, 0.05).is_err());
        assert!(Cylinder::double_acting(0.05, 0.06).is_err());
    }

    #[test]
    fn rejects_negative_gauge_pressure() {
        let c = Cylinder::single_acting(0.05).unwrap();
        let err = c.force(-1.0, Stroke::Extend).unwrap_err();
        assert_eq!(err.code(), "pneumatics.negative");
    }

    #[test]
    fn known_textbook_value_80mm_at_6bar() {
        // Festo-style table check: an 80 mm bore at 6 bar (600 kPa) gauge
        // gives a theoretical extend force of ~3016 N.
        // A = pi/4 * 0.08^2 = 5.0265e-3 m^2; F = 6e5 * A = 3015.93 N.
        let c = Cylinder::single_acting(0.08).unwrap();
        let f = c.extend_force(600_000.0).unwrap();
        assert!((f - 3_015.928_947_446_201).abs() < 1e-6, "got {f}");
    }
}
