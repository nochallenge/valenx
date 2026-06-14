//! The helical-compression-spring geometry + material model.
//!
//! [`HelicalSpring`] bundles the four quantities that fully determine a
//! round-wire helical compression spring's linear-elastic behaviour:
//!
//! - wire diameter `d`,
//! - mean coil diameter `D` (centre-to-centre of the wire),
//! - number of *active* coils `N`,
//! - material shear modulus `G`.
//!
//! From these it derives the spring index, the Wahl curvature-stress
//! correction factor, the spring rate, and (given an applied force)
//! the deflection and the corrected torsional shear stress.
//!
//! ## Units
//!
//! The crate is unit-agnostic in the sense that the formulas are
//! dimensionally consistent: pick any coherent set. The doc examples
//! use millimetres for lengths, newtons for force, and megapascals
//! (`N/mm^2`) for the modulus, which yields a rate in `N/mm` and a
//! stress in `MPa`. See [`crate`] for the honest-scope note.

use serde::{Deserialize, Serialize};

use crate::error::{require_positive, SpringError};

/// A round-wire helical compression spring, defined by its geometry
/// and material shear modulus.
///
/// Construct with [`HelicalSpring::new`], which validates every field.
/// All accessor/derived methods are then total (cannot fail) because
/// the invariants — finite, positive, and `D > d` — are guaranteed at
/// construction time.
///
/// # Sign / unit conventions
///
/// Use a single coherent unit system. With mm / N / MPa: rate is
/// `N/mm`, deflection is `mm`, stress is `MPa`.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct HelicalSpring {
    /// Wire diameter `d` (length units, e.g. mm). Strictly positive.
    wire_diameter: f64,
    /// Mean coil diameter `D`, centre-to-centre of the wire cross
    /// section (length units, e.g. mm). Strictly positive and `> d`.
    mean_coil_diameter: f64,
    /// Number of *active* coils `N` (dimensionless). Strictly positive;
    /// may be fractional. Only the coils free to deflect count — seated
    /// / ground end coils are excluded.
    active_coils: f64,
    /// Material shear modulus `G` (pressure units, e.g. MPa). Strictly
    /// positive. For music wire / spring steel `G ~ 79_300 MPa`.
    shear_modulus: f64,
}

impl HelicalSpring {
    /// Build a validated [`HelicalSpring`].
    ///
    /// # Parameters
    ///
    /// - `wire_diameter` — `d`, the wire's diameter.
    /// - `mean_coil_diameter` — `D`, the coil's mean diameter.
    /// - `active_coils` — `N`, the number of active (deflecting) coils.
    /// - `shear_modulus` — `G`, the material's modulus of rigidity.
    ///
    /// # Errors
    ///
    /// Returns [`SpringError::NonPositive`] if any argument is not a
    /// finite, strictly-positive number, and
    /// [`SpringError::Degenerate`] if `mean_coil_diameter` is not
    /// strictly greater than `wire_diameter` (which would make the
    /// spring index `C <= 1` — physically un-windable).
    pub fn new(
        wire_diameter: f64,
        mean_coil_diameter: f64,
        active_coils: f64,
        shear_modulus: f64,
    ) -> Result<Self, SpringError> {
        let wire_diameter =
            require_positive(wire_diameter, "wire_diameter", "wire diameter must be > 0")?;
        let mean_coil_diameter = require_positive(
            mean_coil_diameter,
            "mean_coil_diameter",
            "mean coil diameter must be > 0",
        )?;
        let active_coils = require_positive(
            active_coils,
            "active_coils",
            "active coil count must be > 0",
        )?;
        let shear_modulus =
            require_positive(shear_modulus, "shear_modulus", "shear modulus must be > 0")?;

        if mean_coil_diameter <= wire_diameter {
            return Err(SpringError::Degenerate(format!(
                "mean coil diameter D = {mean_coil_diameter} must exceed wire diameter d = {wire_diameter} (spring index C = D/d must be > 1)"
            )));
        }

        Ok(Self {
            wire_diameter,
            mean_coil_diameter,
            active_coils,
            shear_modulus,
        })
    }

    /// Wire diameter `d`.
    pub fn wire_diameter(&self) -> f64 {
        self.wire_diameter
    }

    /// Mean coil diameter `D`.
    pub fn mean_coil_diameter(&self) -> f64 {
        self.mean_coil_diameter
    }

    /// Number of active coils `N`.
    pub fn active_coils(&self) -> f64 {
        self.active_coils
    }

    /// Material shear modulus `G`.
    pub fn shear_modulus(&self) -> f64 {
        self.shear_modulus
    }

    /// Spring index `C = D / d` (dimensionless).
    ///
    /// The ratio of mean coil diameter to wire diameter. By
    /// construction `C > 1`; well-behaved manufacturable springs
    /// typically fall in the range `4 <= C <= 12`. Values below ~4 are
    /// hard to coil and develop high stress; above ~12 the spring tends
    /// to tangle and buckle. See [`HelicalSpring::index_is_typical`].
    pub fn spring_index(&self) -> f64 {
        self.mean_coil_diameter / self.wire_diameter
    }

    /// Whether the spring index lies in the conventional
    /// manufacturable band `4 <= C <= 12`.
    ///
    /// This is advisory only — an out-of-band index is not an error,
    /// it is merely outside the range most spring-design texts
    /// recommend.
    pub fn index_is_typical(&self) -> bool {
        let c = self.spring_index();
        (4.0..=12.0).contains(&c)
    }

    /// Wahl curvature-correction factor `K_w` (dimensionless).
    ///
    /// Accounts for the stress concentration on the inner face of the
    /// coil due to wire curvature and the direct (transverse) shear:
    ///
    /// ```text
    ///   K_w = (4C - 1)/(4C - 4)  +  0.615 / C
    /// ```
    ///
    /// `K_w` is always `> 1`, is largest for tightly-wound (small-`C`)
    /// springs, and approaches `1` as `C -> infinity` (a nearly
    /// straight wire has negligible curvature correction).
    pub fn wahl_factor(&self) -> f64 {
        let c = self.spring_index();
        (4.0 * c - 1.0) / (4.0 * c - 4.0) + 0.615 / c
    }

    /// Spring rate `k = G·d^4 / (8·D^3·N)` (force per unit length).
    ///
    /// The linear stiffness relating axial force to axial deflection,
    /// `F = k·delta`. Note the strong dependences: rate goes as the
    /// **fourth power** of wire diameter, the **inverse cube** of coil
    /// diameter, and the **inverse** of the active-coil count.
    ///
    /// With mm / N / MPa inputs the result is in `N/mm`.
    pub fn rate(&self) -> f64 {
        let d = self.wire_diameter;
        let big_d = self.mean_coil_diameter;
        let n = self.active_coils;
        self.shear_modulus * d.powi(4) / (8.0 * big_d.powi(3) * n)
    }

    /// Axial deflection `delta = F / k` under applied force `force`.
    ///
    /// # Errors
    ///
    /// Returns [`SpringError::NonPositive`] if `force` is not finite and
    /// `> 0`.
    pub fn deflection(&self, force: f64) -> Result<f64, SpringError> {
        let force = require_positive(force, "force", "applied force must be > 0")?;
        Ok(force / self.rate())
    }

    /// Force `F = k·delta` required to produce axial deflection
    /// `deflection`.
    ///
    /// The inverse of [`HelicalSpring::deflection`]; the two round-trip
    /// exactly (to floating-point precision) for any valid input.
    ///
    /// # Errors
    ///
    /// Returns [`SpringError::NonPositive`] if `deflection` is not
    /// finite and `> 0`.
    pub fn force_for_deflection(&self, deflection: f64) -> Result<f64, SpringError> {
        let deflection = require_positive(deflection, "deflection", "deflection must be > 0")?;
        Ok(self.rate() * deflection)
    }

    /// Wahl-corrected torsional shear stress
    /// `tau = K_w · 8·F·D / (pi·d^3)` under applied force `force`.
    ///
    /// This is the maximum shear stress at the inner fibre of the coil,
    /// including the Wahl curvature/direct-shear correction
    /// ([`HelicalSpring::wahl_factor`]). It scales linearly with the
    /// applied `force` and with the slenderness group `D / d^3`.
    ///
    /// With mm / N inputs the result is in `MPa` (`N/mm^2`).
    ///
    /// # Errors
    ///
    /// Returns [`SpringError::NonPositive`] if `force` is not finite and
    /// `> 0`.
    pub fn shear_stress(&self, force: f64) -> Result<f64, SpringError> {
        let force = require_positive(force, "force", "applied force must be > 0")?;
        let d = self.wire_diameter;
        let big_d = self.mean_coil_diameter;
        let uncorrected = 8.0 * force * big_d / (std::f64::consts::PI * d.powi(3));
        Ok(self.wahl_factor() * uncorrected)
    }

    /// Uncorrected torsional shear stress `tau_0 = 8·F·D / (pi·d^3)`,
    /// i.e. [`HelicalSpring::shear_stress`] with the Wahl factor
    /// omitted.
    ///
    /// Exposed for comparison and teaching; for any real stress check
    /// prefer the Wahl-corrected [`HelicalSpring::shear_stress`], which
    /// is always larger.
    ///
    /// # Errors
    ///
    /// Returns [`SpringError::NonPositive`] if `force` is not finite and
    /// `> 0`.
    pub fn shear_stress_uncorrected(&self, force: f64) -> Result<f64, SpringError> {
        let force = require_positive(force, "force", "applied force must be > 0")?;
        let d = self.wire_diameter;
        let big_d = self.mean_coil_diameter;
        Ok(8.0 * force * big_d / (std::f64::consts::PI * d.powi(3)))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Absolute tolerance for floating-point comparisons. The formulas
    /// are elementary so results match analytic hand-calcs to many
    /// digits; this bound is comfortably loose.
    const EPS: f64 = 1e-9;

    /// Spring steel modulus of rigidity, MPa.
    const G_STEEL: f64 = 79_300.0;

    fn sample() -> HelicalSpring {
        // d = 2 mm, D = 16 mm  -> C = 8, N = 10 active coils.
        HelicalSpring::new(2.0, 16.0, 10.0, G_STEEL).unwrap()
    }

    #[test]
    fn rejects_non_positive_inputs() {
        assert!(HelicalSpring::new(0.0, 16.0, 10.0, G_STEEL).is_err());
        assert!(HelicalSpring::new(2.0, -1.0, 10.0, G_STEEL).is_err());
        assert!(HelicalSpring::new(2.0, 16.0, 0.0, G_STEEL).is_err());
        assert!(HelicalSpring::new(2.0, 16.0, 10.0, f64::NAN).is_err());
        assert!(HelicalSpring::new(f64::INFINITY, 16.0, 10.0, G_STEEL).is_err());
    }

    #[test]
    fn rejects_wire_not_smaller_than_coil() {
        // D == d  -> index 1, un-windable.
        let err = HelicalSpring::new(5.0, 5.0, 4.0, G_STEEL).unwrap_err();
        assert_eq!(err.code(), "spring.degenerate");
        // D < d.
        assert!(HelicalSpring::new(6.0, 5.0, 4.0, G_STEEL).is_err());
    }

    #[test]
    fn spring_index_is_d_over_d() {
        let s = sample();
        assert!((s.spring_index() - 8.0).abs() < EPS);
        assert!(s.index_is_typical());
        // Out-of-band indices flagged but still constructible.
        let tight = HelicalSpring::new(2.0, 6.0, 5.0, G_STEEL).unwrap(); // C = 3
        assert!(!tight.index_is_typical());
        let loose = HelicalSpring::new(1.0, 15.0, 5.0, G_STEEL).unwrap(); // C = 15
        assert!(!loose.index_is_typical());
    }

    #[test]
    fn rate_matches_closed_form() {
        let s = sample();
        // k = G d^4 / (8 D^3 N)
        //   = 79300 * 2^4 / (8 * 16^3 * 10)
        //   = 79300 * 16 / (8 * 4096 * 10)
        //   = 1_268_800 / 327_680
        //   = 3.872070312500 N/mm
        let expected = 79_300.0 * 2.0_f64.powi(4) / (8.0 * 16.0_f64.powi(3) * 10.0);
        assert!((s.rate() - expected).abs() < EPS);
        assert!((s.rate() - 3.8720703125).abs() < 1e-7);
    }

    #[test]
    fn rate_scales_as_d_to_the_fourth() {
        // Doubling the wire diameter multiplies the rate by 2^4 = 16.
        let base = HelicalSpring::new(2.0, 16.0, 10.0, G_STEEL).unwrap();
        let thick = HelicalSpring::new(4.0, 16.0, 10.0, G_STEEL).unwrap();
        let ratio = thick.rate() / base.rate();
        assert!((ratio - 16.0).abs() < 1e-6, "ratio was {ratio}");
    }

    #[test]
    fn rate_scales_as_inverse_d_cubed() {
        // Doubling the mean coil diameter divides the rate by 2^3 = 8.
        let small = HelicalSpring::new(2.0, 16.0, 10.0, G_STEEL).unwrap();
        let big = HelicalSpring::new(2.0, 32.0, 10.0, G_STEEL).unwrap();
        let ratio = small.rate() / big.rate();
        assert!((ratio - 8.0).abs() < 1e-6, "ratio was {ratio}");
    }

    #[test]
    fn rate_scales_as_inverse_coils() {
        // Doubling the active coils halves the rate.
        let few = HelicalSpring::new(2.0, 16.0, 10.0, G_STEEL).unwrap();
        let many = HelicalSpring::new(2.0, 16.0, 20.0, G_STEEL).unwrap();
        let ratio = few.rate() / many.rate();
        assert!((ratio - 2.0).abs() < 1e-9, "ratio was {ratio}");
    }

    #[test]
    fn deflection_is_force_over_rate() {
        let s = sample();
        let f = 40.0;
        let expected = f / s.rate();
        assert!((s.deflection(f).unwrap() - expected).abs() < EPS);
        // Numerically: 40 / 3.8720703125 = 10.330390920... mm
        assert!((s.deflection(f).unwrap() - 10.330390920555).abs() < 1e-6);
    }

    #[test]
    fn force_and_deflection_round_trip() {
        let s = sample();
        let f = 37.5;
        let delta = s.deflection(f).unwrap();
        let back = s.force_for_deflection(delta).unwrap();
        assert!((back - f).abs() < 1e-9, "round-trip force was {back}");
    }

    #[test]
    fn deflection_rejects_bad_force() {
        let s = sample();
        assert!(s.deflection(0.0).is_err());
        assert!(s.deflection(-5.0).is_err());
        assert!(s.deflection(f64::NAN).is_err());
    }

    #[test]
    fn wahl_factor_exceeds_one_and_decays_with_index() {
        // K_w > 1 for every realistic index, and monotonically
        // approaches 1 as C grows.
        let c4 = HelicalSpring::new(2.0, 8.0, 5.0, G_STEEL).unwrap(); // C = 4
        let c8 = HelicalSpring::new(2.0, 16.0, 5.0, G_STEEL).unwrap(); // C = 8
        let c20 = HelicalSpring::new(1.0, 20.0, 5.0, G_STEEL).unwrap(); // C = 20
        assert!(c4.wahl_factor() > 1.0);
        assert!(c8.wahl_factor() > 1.0);
        assert!(c20.wahl_factor() > 1.0);
        assert!(c4.wahl_factor() > c8.wahl_factor());
        assert!(c8.wahl_factor() > c20.wahl_factor());
        // Trends towards 1 at very large index.
        let c200 = HelicalSpring::new(0.1, 20.0, 5.0, G_STEEL).unwrap(); // C = 200
        assert!(c200.wahl_factor() < 1.02);
        assert!(c200.wahl_factor() > 1.0);
    }

    #[test]
    fn wahl_factor_matches_closed_form_at_c8() {
        // K_w(8) = (31/28) + 0.615/8 = 1.1071428571 + 0.076875
        //        = 1.1840178571
        let s = sample();
        let expected = (4.0 * 8.0 - 1.0) / (4.0 * 8.0 - 4.0) + 0.615 / 8.0;
        assert!((s.wahl_factor() - expected).abs() < EPS);
        assert!((s.wahl_factor() - 1.1840178571428).abs() < 1e-9);
    }

    #[test]
    fn shear_stress_matches_closed_form() {
        let s = sample();
        let f = 40.0;
        // tau_0 = 8 F D / (pi d^3) = 8*40*16 / (pi * 8) = 5120 / (8 pi)
        //       = 640 / pi = 203.7183... MPa
        let tau0_expected = 8.0 * f * 16.0 / (std::f64::consts::PI * 2.0_f64.powi(3));
        assert!((s.shear_stress_uncorrected(f).unwrap() - tau0_expected).abs() < 1e-6);
        assert!(
            (s.shear_stress_uncorrected(f).unwrap() - 640.0 / std::f64::consts::PI).abs() < 1e-6
        );
        // Corrected = K_w * tau_0.
        let expected = s.wahl_factor() * tau0_expected;
        assert!((s.shear_stress(f).unwrap() - expected).abs() < 1e-6);
        // And corrected is strictly larger.
        assert!(s.shear_stress(f).unwrap() > s.shear_stress_uncorrected(f).unwrap());
    }

    #[test]
    fn shear_stress_scales_linearly_with_force() {
        let s = sample();
        let a = s.shear_stress(10.0).unwrap();
        let b = s.shear_stress(30.0).unwrap();
        assert!((b / a - 3.0).abs() < 1e-9, "ratio was {}", b / a);
    }

    #[test]
    fn shear_stress_scales_as_d_over_d_cubed() {
        // tau_0 group is D / d^3. Holding force fixed, halving d (with
        // D scaled to keep things constructible) multiplies the
        // uncorrected stress group D/d^3 accordingly. Here: keep D
        // fixed, halve d -> 1/d^3 grows by 8x.
        let f = 25.0;
        let base = HelicalSpring::new(2.0, 16.0, 5.0, G_STEEL).unwrap();
        let thin = HelicalSpring::new(1.0, 16.0, 5.0, G_STEEL).unwrap();
        let ratio =
            thin.shear_stress_uncorrected(f).unwrap() / base.shear_stress_uncorrected(f).unwrap();
        assert!((ratio - 8.0).abs() < 1e-6, "ratio was {ratio}");

        // And doubling D (d fixed) doubles the uncorrected stress.
        let wide = HelicalSpring::new(2.0, 32.0, 5.0, G_STEEL).unwrap();
        let ratio_d =
            wide.shear_stress_uncorrected(f).unwrap() / base.shear_stress_uncorrected(f).unwrap();
        assert!((ratio_d - 2.0).abs() < 1e-9, "ratio was {ratio_d}");
    }

    #[test]
    fn shear_stress_rejects_bad_force() {
        let s = sample();
        assert!(s.shear_stress(0.0).is_err());
        assert!(s.shear_stress(-1.0).is_err());
        assert!(s.shear_stress_uncorrected(f64::NEG_INFINITY).is_err());
    }

    #[test]
    fn serde_round_trip() {
        let s = sample();
        let json = serde_json::to_string(&s).unwrap();
        let back: HelicalSpring = serde_json::from_str(&json).unwrap();
        assert_eq!(s, back);
    }
}
