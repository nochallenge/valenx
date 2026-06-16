//! Differential-pressure flow-meter geometry and the incompressible
//! flow / pressure-drop relations.
//!
//! # Model
//!
//! A differential-pressure (dP) meter places a constriction of bore
//! diameter `d` in a pipe of internal diameter `D` and infers the
//! volumetric flow rate `Q` from the static-pressure difference `dP`
//! measured across that constriction. Applying Bernoulli's equation and
//! continuity to an incompressible, inviscid stream, then folding all
//! the real-fluid losses into a single empirical discharge coefficient
//! `Cd`, gives the textbook working equation
//!
//! ```text
//! Q = Cd * A * sqrt( 2 * dP / (rho * (1 - beta^4)) )
//! ```
//!
//! where
//!
//! ```text
//! A    = pi * d^2 / 4      throat (bore) cross-sectional area  [m^2]
//! beta = d / D             diameter ratio                       [-]
//! rho                      fluid density                        [kg/m^3]
//! dP                       upstream-minus-throat pressure       [Pa]
//! ```
//!
//! The grouping `E = 1 / sqrt(1 - beta^4)` is the *velocity-of-approach
//! factor*; it accounts for the kinetic energy the fluid already carries
//! in the approach pipe and is always greater than one for `0 < beta < 1`.
//!
//! Solving the same equation for `dP` gives the exact inverse used for
//! sizing the differential-pressure transmitter:
//!
//! ```text
//! dP = (rho * (1 - beta^4) / 2) * ( Q / (Cd * A) )^2
//! ```
//!
//! Mass flow follows from the density: `mdot = rho * Q`.
//!
//! Only part of the measured differential is lost for good; the jet
//! re-expands downstream and recovers the rest. The permanently
//! unrecovered fraction of a square-edged orifice plate is the ISO 5167-2
//! ratio
//!
//! ```text
//! dOmega / dP = ( sqrt(1 - beta^4 (1 - Cd^2)) - Cd beta^2 )
//!             / ( sqrt(1 - beta^4 (1 - Cd^2)) + Cd beta^2 ).
//! ```
//!
//! # Honest scope
//!
//! This is a research / educational implementation of the *closed-form*
//! incompressible relation only. It treats `Cd` as a caller-supplied
//! constant and assumes an incompressible fluid (no gas-expansion
//! `epsilon` factor). It does NOT implement the Reynolds-number-dependent
//! Reader-Harris / Gallagher `Cd` correlation, the thermal-expansion
//! corrections, or the installation / tapping requirements of ISO 5167
//! or ASME MFC-3M, and it is not a substitute for an accredited flow
//! calibration. The permanent-loss relation is the square-edged
//! orifice-plate form and overestimates the loss of nozzles and Venturis,
//! which recover pressure in their diffusers. Use it to learn how a dP
//! meter scales, not to bill custody-transfer gas.

use crate::error::{require_non_negative, require_positive, OrificeError};
use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

/// Family of differential-pressure flow meter, each with a representative
/// discharge coefficient.
///
/// The associated [`MeterKind::typical_cd`] values are round textbook
/// figures for fully-developed turbulent flow; they are convenience
/// defaults only. Real meters use a Reynolds-number-dependent `Cd` from
/// the relevant standard, which this crate does not model — see the
/// crate-level honest-scope note.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash, Serialize, Deserialize)]
pub enum MeterKind {
    /// A thin square-edged orifice plate. Highest permanent pressure
    /// loss, lowest discharge coefficient (sharp-edged contraction with a
    /// vena contracta). Typical `Cd` is approximately `0.61`.
    OrificePlate,
    /// A long-radius / ISA flow nozzle. A smoothly converging inlet
    /// recovers more of the ideal flow than an orifice. Typical `Cd` is
    /// approximately `0.97`.
    FlowNozzle,
    /// A classical Venturi tube (converging cone, throat, diverging
    /// cone). The gentle contraction and diffuser give the highest
    /// discharge coefficient and the lowest permanent loss. Typical `Cd`
    /// is approximately `0.98`.
    VenturiTube,
}

impl MeterKind {
    /// A representative discharge coefficient for this meter family in
    /// fully-developed turbulent flow.
    ///
    /// These are nominal textbook values for convenience and for the
    /// ordering test `Venturi > nozzle > orifice`; they are not a
    /// calibrated `Cd` for any specific installation.
    #[must_use]
    pub fn typical_cd(self) -> f64 {
        match self {
            MeterKind::OrificePlate => 0.61,
            MeterKind::FlowNozzle => 0.97,
            MeterKind::VenturiTube => 0.98,
        }
    }

    /// A stable kebab-cased identifier for the variant.
    #[must_use]
    pub fn code(self) -> &'static str {
        match self {
            MeterKind::OrificePlate => "orifice-plate",
            MeterKind::FlowNozzle => "flow-nozzle",
            MeterKind::VenturiTube => "venturi-tube",
        }
    }
}

/// The geometry of a differential-pressure meter: a throat (bore) of
/// diameter `d` inside an upstream pipe of internal diameter `D`.
///
/// Construct with [`MeterGeometry::new`], which validates that both
/// diameters are finite and strictly positive and that the throat is
/// strictly smaller than the pipe (so that `beta = d / D` lies in the
/// open interval `(0, 1)`). All diameters are in metres.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct MeterGeometry {
    /// Throat / bore diameter `d` in metres.
    throat_diameter: f64,
    /// Upstream pipe internal diameter `D` in metres.
    pipe_diameter: f64,
}

impl MeterGeometry {
    /// Build a meter geometry from the throat diameter `d` and the pipe
    /// diameter `D`, both in metres.
    ///
    /// # Errors
    ///
    /// Returns [`OrificeError::NotFinite`] or [`OrificeError::NonPositive`]
    /// if either diameter is not a finite, strictly-positive number, and
    /// [`OrificeError::ThroatNotSmaller`] if `d >= D`.
    pub fn new(throat_diameter: f64, pipe_diameter: f64) -> Result<Self, OrificeError> {
        let throat_diameter = require_positive("throat_diameter", throat_diameter)?;
        let pipe_diameter = require_positive("pipe_diameter", pipe_diameter)?;
        if throat_diameter >= pipe_diameter {
            return Err(OrificeError::ThroatNotSmaller {
                d: throat_diameter,
                pipe: pipe_diameter,
            });
        }
        Ok(Self {
            throat_diameter,
            pipe_diameter,
        })
    }

    /// The throat / bore diameter `d` in metres.
    #[must_use]
    pub fn throat_diameter(self) -> f64 {
        self.throat_diameter
    }

    /// The upstream pipe internal diameter `D` in metres.
    #[must_use]
    pub fn pipe_diameter(self) -> f64 {
        self.pipe_diameter
    }

    /// The diameter ratio `beta = d / D` (dimensionless, in `(0, 1)`).
    #[must_use]
    pub fn beta(self) -> f64 {
        self.throat_diameter / self.pipe_diameter
    }

    /// The throat cross-sectional area `A = pi * d^2 / 4` in square
    /// metres.
    #[must_use]
    pub fn throat_area(self) -> f64 {
        0.25 * PI * self.throat_diameter * self.throat_diameter
    }

    /// The velocity-of-approach factor `E = 1 / sqrt(1 - beta^4)`
    /// (dimensionless).
    ///
    /// Because `0 < beta < 1`, the term `1 - beta^4` lies in `(0, 1)`
    /// and `E` is always strictly greater than one. As `beta -> 0` (a
    /// throat much smaller than the pipe) `E -> 1`; as `beta -> 1` the
    /// factor grows without bound.
    #[must_use]
    pub fn velocity_of_approach_factor(self) -> f64 {
        let beta = self.beta();
        let beta4 = beta * beta * beta * beta;
        1.0 / (1.0 - beta4).sqrt()
    }
}

/// A complete differential-pressure flow meter: a [`MeterGeometry`] plus
/// a discharge coefficient `Cd`.
///
/// Construct with [`Meter::new`] (explicit `Cd`) or
/// [`Meter::with_typical_cd`] (the nominal `Cd` for a [`MeterKind`]).
/// Once built, [`Meter::flow_rate`] and [`Meter::pressure_drop`] are the
/// forward and inverse working relations; they form an exact round trip.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Meter {
    geometry: MeterGeometry,
    discharge_coefficient: f64,
}

impl Meter {
    /// Build a meter from an explicit geometry and discharge coefficient.
    ///
    /// # Errors
    ///
    /// Returns [`OrificeError::NotFinite`] if `Cd` is not finite, or
    /// [`OrificeError::DischargeCoefficientOutOfRange`] if `Cd` is not in
    /// the half-open range `(0, 1]`.
    pub fn new(geometry: MeterGeometry, discharge_coefficient: f64) -> Result<Self, OrificeError> {
        if !discharge_coefficient.is_finite() {
            return Err(OrificeError::NotFinite {
                name: "discharge_coefficient",
            });
        }
        if discharge_coefficient <= 0.0 || discharge_coefficient > 1.0 {
            return Err(OrificeError::DischargeCoefficientOutOfRange {
                value: discharge_coefficient,
            });
        }
        Ok(Self {
            geometry,
            discharge_coefficient,
        })
    }

    /// Build a meter using the nominal discharge coefficient for the
    /// given [`MeterKind`] (see [`MeterKind::typical_cd`]).
    ///
    /// # Errors
    ///
    /// Propagates any error from [`Meter::new`]; in practice the only way
    /// this fails is never, since every `typical_cd` is a valid finite
    /// value in `(0, 1]`, but the signature stays fallible for symmetry.
    pub fn with_typical_cd(geometry: MeterGeometry, kind: MeterKind) -> Result<Self, OrificeError> {
        Meter::new(geometry, kind.typical_cd())
    }

    /// The underlying [`MeterGeometry`].
    #[must_use]
    pub fn geometry(self) -> MeterGeometry {
        self.geometry
    }

    /// The discharge coefficient `Cd` (dimensionless, in `(0, 1]`).
    #[must_use]
    pub fn discharge_coefficient(self) -> f64 {
        self.discharge_coefficient
    }

    /// The flow coefficient `K = Cd * E = Cd / sqrt(1 - beta^4)`
    /// (dimensionless).
    ///
    /// This bundles the discharge coefficient with the velocity-of-
    /// approach factor, so the working equation reads simply
    /// `Q = K * A * sqrt(2 * dP / rho)`. Because `E > 1`, `K` is always
    /// larger than `Cd` alone.
    #[must_use]
    pub fn flow_coefficient(self) -> f64 {
        self.discharge_coefficient * self.geometry.velocity_of_approach_factor()
    }

    /// Volumetric flow rate `Q` (cubic metres per second) for a given
    /// fluid density `rho` (kg/m^3) and pressure drop `dP` (Pa).
    ///
    /// Evaluates the forward working equation
    ///
    /// ```text
    /// Q = Cd * A * sqrt( 2 * dP / (rho * (1 - beta^4)) ).
    /// ```
    ///
    /// A pressure drop of exactly zero yields a flow of exactly zero.
    ///
    /// # Errors
    ///
    /// Returns [`OrificeError::NonPositive`] if `rho <= 0` (or not
    /// finite) and [`OrificeError::Negative`] if `dP < 0` (or not finite).
    pub fn flow_rate(self, density: f64, pressure_drop: f64) -> Result<f64, OrificeError> {
        let density = require_positive("density", density)?;
        let pressure_drop = require_non_negative("pressure_drop", pressure_drop)?;
        let beta = self.geometry.beta();
        let beta4 = beta * beta * beta * beta;
        let area = self.geometry.throat_area();
        let q = self.discharge_coefficient
            * area
            * (2.0 * pressure_drop / (density * (1.0 - beta4))).sqrt();
        Ok(q)
    }

    /// Mass flow rate `mdot` (kilograms per second): `mdot = rho * Q`.
    ///
    /// # Errors
    ///
    /// Same as [`Meter::flow_rate`].
    pub fn mass_flow_rate(self, density: f64, pressure_drop: f64) -> Result<f64, OrificeError> {
        let q = self.flow_rate(density, pressure_drop)?;
        Ok(density * q)
    }

    /// Pressure drop `dP` (Pa) required to drive a given volumetric flow
    /// rate `Q` (m^3/s) through the meter at fluid density `rho`
    /// (kg/m^3).
    ///
    /// Evaluates the exact inverse of [`Meter::flow_rate`]:
    ///
    /// ```text
    /// dP = (rho * (1 - beta^4) / 2) * ( Q / (Cd * A) )^2.
    /// ```
    ///
    /// A flow of exactly zero yields a pressure drop of exactly zero.
    ///
    /// # Errors
    ///
    /// Returns [`OrificeError::NonPositive`] if `rho <= 0` (or not
    /// finite) and [`OrificeError::Negative`] if `Q < 0` (or not finite).
    pub fn pressure_drop(self, density: f64, flow_rate: f64) -> Result<f64, OrificeError> {
        let density = require_positive("density", density)?;
        let flow_rate = require_non_negative("flow_rate", flow_rate)?;
        let beta = self.geometry.beta();
        let beta4 = beta * beta * beta * beta;
        let area = self.geometry.throat_area();
        let ratio = flow_rate / (self.discharge_coefficient * area);
        let dp = 0.5 * density * (1.0 - beta4) * ratio * ratio;
        Ok(dp)
    }

    /// The permanent (unrecovered) pressure-loss ratio `dOmega / dP` of a
    /// square-edged orifice plate — the ISO 5167-2 relation
    ///
    /// ```text
    /// dOmega / dP = ( sqrt(1 - beta^4 (1 - Cd^2)) - Cd beta^2 )
    ///             / ( sqrt(1 - beta^4 (1 - Cd^2)) + Cd beta^2 )
    /// ```
    ///
    /// (dimensionless, in `(0, 1)`). The differential `dP` is what the
    /// transmitter measures across the plate; only this fraction of it is
    /// *permanently* lost downstream — the pumping penalty the meter
    /// imposes — the remainder being recovered as the jet re-expands.
    /// The ratio falls monotonically from one as `beta -> 0` (a tiny bore
    /// nearly blocks the pipe) toward zero as `beta -> 1` (almost no
    /// constriction). In the lossless limit `Cd = 1` it collapses to the
    /// clean form `(1 - beta^2) / (1 + beta^2)`.
    ///
    /// This is the square-edged **orifice-plate** relation. A flow nozzle
    /// or Venturi recovers pressure in its smooth diffuser and so loses
    /// substantially less than this predicts; applying it to those meter
    /// kinds overestimates the permanent loss.
    #[must_use]
    pub fn permanent_pressure_loss_ratio(self) -> f64 {
        let beta2 = self.geometry.beta() * self.geometry.beta();
        let beta4 = beta2 * beta2;
        let cd = self.discharge_coefficient;
        let root = (1.0 - beta4 * (1.0 - cd * cd)).sqrt();
        let cd_beta2 = cd * beta2;
        (root - cd_beta2) / (root + cd_beta2)
    }

    /// The absolute permanent (unrecovered) pressure loss `dOmega` (Pa)
    /// for a volumetric flow rate `Q` (m^3/s) at fluid density `rho`
    /// (kg/m^3): `dOmega = (dOmega / dP) * dP`, combining
    /// [`permanent_pressure_loss_ratio`](Meter::permanent_pressure_loss_ratio)
    /// with the differential [`pressure_drop`](Meter::pressure_drop) that
    /// flow produces. Because the ratio is below one, the permanent loss
    /// is always strictly less than the measured differential.
    ///
    /// As with the ratio, this is the square-edged orifice-plate relation
    /// — see
    /// [`permanent_pressure_loss_ratio`](Meter::permanent_pressure_loss_ratio).
    ///
    /// # Errors
    ///
    /// Returns [`OrificeError::NonPositive`] if `rho <= 0` (or not finite)
    /// and [`OrificeError::Negative`] if `Q < 0` (or not finite).
    pub fn permanent_pressure_loss(
        self,
        density: f64,
        flow_rate: f64,
    ) -> Result<f64, OrificeError> {
        let dp = self.pressure_drop(density, flow_rate)?;
        Ok(self.permanent_pressure_loss_ratio() * dp)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::error::ErrorCategory;

    /// Tight tolerance for comparisons against an independently computed
    /// closed-form value.
    const EPS: f64 = 1e-12;

    fn geom(d: f64, big_d: f64) -> MeterGeometry {
        MeterGeometry::new(d, big_d).expect("valid geometry")
    }

    fn meter(d: f64, big_d: f64, cd: f64) -> Meter {
        Meter::new(geom(d, big_d), cd).expect("valid meter")
    }

    // --- geometry basics ---------------------------------------------------

    #[test]
    fn beta_and_area_match_definitions() {
        let g = geom(0.05, 0.10);
        assert!((g.beta() - 0.5).abs() < EPS, "beta = d/D");
        // A = pi d^2 / 4 computed independently.
        let expected_area = std::f64::consts::PI * 0.05 * 0.05 / 4.0;
        assert!(
            (g.throat_area() - expected_area).abs() < EPS,
            "area {a} vs {expected_area}",
            a = g.throat_area()
        );
        assert!((g.throat_diameter() - 0.05).abs() < EPS);
        assert!((g.pipe_diameter() - 0.10).abs() < EPS);
    }

    #[test]
    fn velocity_of_approach_factor_ground_truth() {
        // beta = 0.5 -> E = 1/sqrt(1 - 0.0625) = 1.0327955589886444.
        let g = geom(0.05, 0.10);
        let e = g.velocity_of_approach_factor();
        assert!(
            (e - 1.032_795_558_988_644_4).abs() < EPS,
            "E = {e}, want 1.0327955589886444"
        );
    }

    #[test]
    fn approach_factor_strictly_exceeds_one_for_all_betas() {
        // For 0 < beta < 1 the factor 1/sqrt(1 - beta^4) is always > 1.
        for &(d, big_d) in &[
            (0.001, 1.0),  // beta ~ 0
            (0.05, 0.10),  // beta = 0.5
            (0.075, 0.10), // beta = 0.75
            (0.099, 0.10), // beta ~ 1
        ] {
            let g = geom(d, big_d);
            let e = g.velocity_of_approach_factor();
            assert!(
                e > 1.0,
                "E must exceed 1, got {e} for beta {b}",
                b = g.beta()
            );
        }
    }

    #[test]
    fn approach_factor_grows_monotonically_with_beta() {
        let small = geom(0.02, 0.10).velocity_of_approach_factor(); // beta 0.2
        let mid = geom(0.05, 0.10).velocity_of_approach_factor(); // beta 0.5
        let large = geom(0.08, 0.10).velocity_of_approach_factor(); // beta 0.8
        assert!(small < mid, "{small} < {mid}");
        assert!(mid < large, "{mid} < {large}");
    }

    #[test]
    fn approach_factor_tends_to_one_as_beta_tends_to_zero() {
        // A throat tiny relative to the pipe -> E ~ 1.
        let e = geom(1.0e-4, 1.0).velocity_of_approach_factor();
        assert!((e - 1.0).abs() < 1e-12, "E -> 1 as beta -> 0, got {e}");
    }

    // --- forward flow: ground truth ----------------------------------------

    #[test]
    fn flow_rate_matches_independent_hand_calc() {
        // d=0.05, D=0.10, Cd=0.61, rho=1000, dP=50000 ->
        // Q = 0.012370124961719516 m^3/s (computed offline).
        let m = meter(0.05, 0.10, 0.61);
        let q = m.flow_rate(1000.0, 50_000.0).expect("flow");
        assert!(
            (q - 0.012_370_124_961_719_516).abs() < 1e-15,
            "Q = {q}, want 0.012370124961719516"
        );
    }

    #[test]
    fn mass_flow_is_density_times_volume_flow() {
        let m = meter(0.05, 0.10, 0.61);
        let q = m.flow_rate(1000.0, 50_000.0).unwrap();
        let mdot = m.mass_flow_rate(1000.0, 50_000.0).unwrap();
        assert!((mdot - 1000.0 * q).abs() < EPS, "mdot = rho * Q");
        // Independent absolute check.
        assert!(
            (mdot - 12.370_124_961_719_517).abs() < 1e-12,
            "mdot = {mdot}"
        );
    }

    #[test]
    fn flow_coefficient_is_cd_times_approach_factor() {
        let m = meter(0.05, 0.10, 0.61);
        let k = m.flow_coefficient();
        let expected = 0.61 * m.geometry().velocity_of_approach_factor();
        assert!((k - expected).abs() < EPS, "K = Cd * E");
        // K = Cd * E always exceeds Cd alone because E > 1.
        assert!(k > m.discharge_coefficient(), "{k} > Cd");
    }

    // --- scaling laws ------------------------------------------------------

    #[test]
    fn flow_scales_as_sqrt_of_pressure_drop() {
        // Quadrupling dP must exactly double Q (Q ~ sqrt(dP)).
        let m = meter(0.04, 0.10, 0.65);
        let q1 = m.flow_rate(998.0, 25_000.0).unwrap();
        let q4 = m.flow_rate(998.0, 100_000.0).unwrap();
        assert!((q4 - 2.0 * q1).abs() < 1e-12, "4x dP -> 2x Q: {q4} vs {q1}");

        // Nine-fold dP -> three-fold Q.
        let q9 = m.flow_rate(998.0, 225_000.0).unwrap();
        assert!((q9 - 3.0 * q1).abs() < 1e-12, "9x dP -> 3x Q");
    }

    #[test]
    fn flow_is_proportional_to_throat_area() {
        // Same beta, density, Cd, dP; double the throat diameter -> 4x area
        // -> 4x flow. Hold beta constant by doubling D too.
        let small = meter(0.02, 0.04, 0.62); // A = pi*0.02^2/4
        let big = meter(0.04, 0.08, 0.62); // A = pi*0.04^2/4 = 4x
        let qs = small.flow_rate(1000.0, 30_000.0).unwrap();
        let qb = big.flow_rate(1000.0, 30_000.0).unwrap();
        let area_ratio = big.geometry().throat_area() / small.geometry().throat_area();
        assert!((area_ratio - 4.0).abs() < 1e-12, "area ratio = 4");
        assert!(
            (qb - 4.0 * qs).abs() < 1e-12,
            "Q proportional to A: {qb} vs {qs}"
        );
    }

    #[test]
    fn flow_is_proportional_to_discharge_coefficient() {
        // Doubling Cd (0.45 -> 0.90) at fixed geometry doubles Q.
        let lo = meter(0.05, 0.10, 0.45);
        let hi = meter(0.05, 0.10, 0.90);
        let q_lo = lo.flow_rate(1000.0, 50_000.0).unwrap();
        let q_hi = hi.flow_rate(1000.0, 50_000.0).unwrap();
        assert!((q_hi - 2.0 * q_lo).abs() < 1e-12, "Q proportional to Cd");
    }

    #[test]
    fn beta_factor_increases_flow_above_the_no_approach_case() {
        // Q with the (1 - beta^4) term must exceed the hypothetical flow
        // computed as if beta were zero (factor E = 1). The ratio equals E.
        let m = meter(0.08, 0.10, 0.7); // beta = 0.8, sizeable approach effect
        let q = m.flow_rate(1000.0, 40_000.0).unwrap();
        let area = m.geometry().throat_area();
        let q_no_approach = 0.7 * area * (2.0 * 40_000.0 / 1000.0_f64).sqrt();
        assert!(
            q > q_no_approach,
            "approach term raises Q: {q} > {q_no_approach}"
        );
        let ratio = q / q_no_approach;
        let e = m.geometry().velocity_of_approach_factor();
        assert!(
            (ratio - e).abs() < 1e-12,
            "flow ratio equals E: {ratio} vs {e}"
        );
        assert!(ratio > 1.0, "ratio 1/sqrt(1 - beta^4) > 1");
    }

    // --- meter-kind ordering ----------------------------------------------

    #[test]
    fn typical_cd_ordering_venturi_gt_nozzle_gt_orifice() {
        let venturi = MeterKind::VenturiTube.typical_cd();
        let nozzle = MeterKind::FlowNozzle.typical_cd();
        let orifice = MeterKind::OrificePlate.typical_cd();
        assert!(venturi > nozzle, "venturi {venturi} > nozzle {nozzle}");
        assert!(nozzle > orifice, "nozzle {nozzle} > orifice {orifice}");
        // Spot-check the nominal magnitudes.
        assert!((venturi - 0.98).abs() < EPS);
        assert!((orifice - 0.61).abs() < EPS);
    }

    #[test]
    fn venturi_passes_more_flow_than_orifice_at_same_geometry() {
        let g = geom(0.05, 0.10);
        let venturi = Meter::with_typical_cd(g, MeterKind::VenturiTube).unwrap();
        let nozzle = Meter::with_typical_cd(g, MeterKind::FlowNozzle).unwrap();
        let orifice = Meter::with_typical_cd(g, MeterKind::OrificePlate).unwrap();
        let qv = venturi.flow_rate(1000.0, 50_000.0).unwrap();
        let qn = nozzle.flow_rate(1000.0, 50_000.0).unwrap();
        let qo = orifice.flow_rate(1000.0, 50_000.0).unwrap();
        assert!(qv > qn, "venturi flow {qv} > nozzle {qn}");
        assert!(qn > qo, "nozzle flow {qn} > orifice {qo}");
    }

    // --- inverse and round trip -------------------------------------------

    #[test]
    fn pressure_drop_matches_closed_form() {
        // Recompute dP independently from Q.
        let m = meter(0.05, 0.10, 0.61);
        let beta4 = 0.5_f64.powi(4);
        let area = m.geometry().throat_area();
        let q = 0.012_370_124_961_719_516;
        let ratio = q / (0.61 * area);
        let expected_dp = 0.5 * 1000.0 * (1.0 - beta4) * ratio * ratio;
        let dp = m.pressure_drop(1000.0, q).unwrap();
        assert!((dp - expected_dp).abs() < 1e-9, "dP {dp} vs {expected_dp}");
        // And it is the 50 kPa that produced this Q.
        assert!((dp - 50_000.0).abs() < 1e-6, "recovers 50 kPa, got {dp}");
    }

    #[test]
    fn flow_then_pressure_round_trips() {
        let m = meter(0.06, 0.12, 0.72);
        for &dp in &[1.0, 100.0, 12_345.0, 250_000.0] {
            let q = m.flow_rate(950.0, dp).unwrap();
            let back = m.pressure_drop(950.0, q).unwrap();
            assert!(
                (back - dp).abs() < 1e-6 * dp.max(1.0),
                "dP round trip: {back} vs {dp}"
            );
        }
    }

    #[test]
    fn pressure_then_flow_round_trips() {
        let m = meter(0.03, 0.05, 0.6);
        for &q in &[1.0e-4, 1.0e-3, 0.01, 0.05] {
            let dp = m.pressure_drop(1000.0, q).unwrap();
            let back = m.flow_rate(1000.0, dp).unwrap();
            assert!(
                (back - q).abs() < 1e-9 * q.max(1.0),
                "Q round trip: {back} vs {q}"
            );
        }
    }

    #[test]
    fn zero_pressure_drop_gives_zero_flow_and_back() {
        let m = meter(0.05, 0.10, 0.61);
        let q = m.flow_rate(1000.0, 0.0).unwrap();
        assert!(q.abs() < EPS, "zero dP -> zero Q, got {q}");
        let dp = m.pressure_drop(1000.0, 0.0).unwrap();
        assert!(dp.abs() < EPS, "zero Q -> zero dP, got {dp}");
    }

    // --- construction validation ------------------------------------------

    #[test]
    fn geometry_rejects_throat_not_smaller_than_pipe() {
        let err = MeterGeometry::new(0.10, 0.10).unwrap_err();
        assert!(matches!(err, OrificeError::ThroatNotSmaller { .. }));
        assert_eq!(err.code(), "orifice.throat_not_smaller");
        assert_eq!(err.category(), ErrorCategory::Domain);

        let err = MeterGeometry::new(0.20, 0.10).unwrap_err();
        assert!(matches!(err, OrificeError::ThroatNotSmaller { .. }));
    }

    #[test]
    fn geometry_rejects_non_positive_diameters() {
        assert!(matches!(
            MeterGeometry::new(0.0, 0.1).unwrap_err(),
            OrificeError::NonPositive {
                name: "throat_diameter",
                ..
            }
        ));
        assert!(matches!(
            MeterGeometry::new(0.05, -0.1).unwrap_err(),
            OrificeError::NonPositive {
                name: "pipe_diameter",
                ..
            }
        ));
    }

    #[test]
    fn geometry_rejects_non_finite_diameters() {
        assert!(matches!(
            MeterGeometry::new(f64::NAN, 0.1).unwrap_err(),
            OrificeError::NotFinite {
                name: "throat_diameter"
            }
        ));
        assert!(matches!(
            MeterGeometry::new(0.05, f64::INFINITY).unwrap_err(),
            OrificeError::NotFinite {
                name: "pipe_diameter"
            }
        ));
    }

    #[test]
    fn meter_rejects_discharge_coefficient_out_of_range() {
        let g = geom(0.05, 0.10);
        assert!(matches!(
            Meter::new(g, 0.0).unwrap_err(),
            OrificeError::DischargeCoefficientOutOfRange { .. }
        ));
        assert!(matches!(
            Meter::new(g, 1.5).unwrap_err(),
            OrificeError::DischargeCoefficientOutOfRange { .. }
        ));
        assert!(matches!(
            Meter::new(g, f64::NAN).unwrap_err(),
            OrificeError::NotFinite {
                name: "discharge_coefficient"
            }
        ));
        // Cd == 1.0 is allowed (the upper bound is inclusive).
        assert!(Meter::new(g, 1.0).is_ok());
    }

    #[test]
    fn flow_and_pressure_reject_bad_fluid_inputs() {
        let m = meter(0.05, 0.10, 0.61);
        assert!(matches!(
            m.flow_rate(0.0, 1.0).unwrap_err(),
            OrificeError::NonPositive {
                name: "density",
                ..
            }
        ));
        assert!(matches!(
            m.flow_rate(1000.0, -1.0).unwrap_err(),
            OrificeError::Negative {
                name: "pressure_drop",
                ..
            }
        ));
        assert!(matches!(
            m.pressure_drop(-1.0, 1.0).unwrap_err(),
            OrificeError::NonPositive {
                name: "density",
                ..
            }
        ));
        assert!(matches!(
            m.pressure_drop(1000.0, -1.0).unwrap_err(),
            OrificeError::Negative {
                name: "flow_rate",
                ..
            }
        ));
    }

    #[test]
    fn meter_kind_codes_are_stable() {
        assert_eq!(MeterKind::OrificePlate.code(), "orifice-plate");
        assert_eq!(MeterKind::FlowNozzle.code(), "flow-nozzle");
        assert_eq!(MeterKind::VenturiTube.code(), "venturi-tube");
    }

    #[test]
    fn serde_round_trip_for_meter() {
        let m = meter(0.05, 0.10, 0.61);
        let json = serde_json::to_string(&m).expect("serialize");
        let back: Meter = serde_json::from_str(&json).expect("deserialize");
        assert_eq!(m, back);
    }

    // --- permanent pressure loss ------------------------------------------

    #[test]
    fn permanent_loss_ratio_matches_iso_closed_form() {
        // beta = 0.5, Cd = 0.61. Recompute the ISO 5167-2 ratio
        // independently and confirm the method agrees, then sanity-check
        // the magnitude (a beta=0.5 orifice loses roughly 70-75% of dP).
        let m = meter(0.05, 0.10, 0.61);
        let beta2 = 0.25_f64;
        let beta4 = beta2 * beta2;
        let cd = 0.61_f64;
        let root = (1.0 - beta4 * (1.0 - cd * cd)).sqrt();
        let expected = (root - cd * beta2) / (root + cd * beta2);
        let r = m.permanent_pressure_loss_ratio();
        assert!((r - expected).abs() < EPS, "ratio {r} vs {expected}");
        assert!(r > 0.70 && r < 0.75, "beta=0.5 orifice loss ~0.73, got {r}");
    }

    #[test]
    fn permanent_loss_ratio_cd_one_is_clean_identity() {
        // GOLD identity: at Cd = 1 the ratio collapses to the exact
        // closed form (1 - beta^2) / (1 + beta^2).
        for &(d, big_d) in &[(0.02, 0.10), (0.05, 0.10), (0.08, 0.10)] {
            let m = Meter::new(geom(d, big_d), 1.0).unwrap();
            let beta2 = m.geometry().beta() * m.geometry().beta();
            let expected = (1.0 - beta2) / (1.0 + beta2);
            let r = m.permanent_pressure_loss_ratio();
            assert!(
                (r - expected).abs() < EPS,
                "Cd=1 identity: {r} vs {expected}"
            );
        }
        // beta = 0.5, Cd = 1 -> (1 - 0.25)/(1 + 0.25) = 0.6 exactly.
        let m = Meter::new(geom(0.05, 0.10), 1.0).unwrap();
        assert!((m.permanent_pressure_loss_ratio() - 0.6).abs() < EPS);
    }

    #[test]
    fn permanent_loss_ratio_is_bounded_and_decreasing_in_beta() {
        // The ratio lies in (0, 1) and falls monotonically as beta rises
        // (less constriction -> less permanent loss).
        let cd = 0.61;
        let small = meter(0.02, 0.10, cd).permanent_pressure_loss_ratio(); // beta 0.2
        let mid = meter(0.05, 0.10, cd).permanent_pressure_loss_ratio(); // beta 0.5
        let large = meter(0.09, 0.10, cd).permanent_pressure_loss_ratio(); // beta 0.9
        for r in [small, mid, large] {
            assert!(r > 0.0 && r < 1.0, "ratio in (0,1), got {r}");
        }
        assert!(small > mid, "loss falls with beta: {small} > {mid}");
        assert!(mid > large, "loss falls with beta: {mid} > {large}");
    }

    #[test]
    fn permanent_loss_tends_to_zero_as_beta_tends_to_one() {
        // A throat almost as large as the pipe permanently loses almost
        // nothing.
        let r = meter(0.0999, 0.10, 0.61).permanent_pressure_loss_ratio();
        assert!(r < 0.05, "near-open meter loses little, got {r}");
    }

    #[test]
    fn permanent_loss_is_a_fraction_of_the_differential() {
        // dOmega = ratio * dP, always strictly below the measured dP.
        let m = meter(0.05, 0.10, 0.61);
        let ratio = m.permanent_pressure_loss_ratio();
        for &q in &[1.0e-3, 0.01, 0.05] {
            let dp = m.pressure_drop(1000.0, q).unwrap();
            let loss = m.permanent_pressure_loss(1000.0, q).unwrap();
            assert!(
                (loss - ratio * dp).abs() < 1e-9 * dp.max(1.0),
                "dOmega = ratio * dP: {loss} vs {}",
                ratio * dp
            );
            assert!(loss < dp, "permanent loss {loss} < differential {dp}");
            assert!(loss > 0.0, "positive loss for positive flow");
        }
    }

    #[test]
    fn permanent_loss_zero_flow_is_zero() {
        let m = meter(0.05, 0.10, 0.61);
        assert!(m.permanent_pressure_loss(1000.0, 0.0).unwrap().abs() < EPS);
    }

    #[test]
    fn permanent_loss_rejects_bad_fluid_inputs() {
        let m = meter(0.05, 0.10, 0.61);
        assert!(matches!(
            m.permanent_pressure_loss(0.0, 1.0).unwrap_err(),
            OrificeError::NonPositive {
                name: "density",
                ..
            }
        ));
        assert!(matches!(
            m.permanent_pressure_loss(1000.0, -1.0).unwrap_err(),
            OrificeError::Negative {
                name: "flow_rate",
                ..
            }
        ));
    }
}
