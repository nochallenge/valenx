//! Singly-reinforced rectangular RC-beam flexure model.
//!
//! This module implements the textbook ultimate-strength flexure
//! equations for a *singly-reinforced rectangular* concrete section
//! using the Whitney equivalent rectangular stress block. All formulas
//! assume the tension steel yields at nominal capacity (the
//! *under-reinforced*, tension-controlled regime), which is the design
//! intent of code-conforming flexural members.
//!
//! ## Unit system
//!
//! The model is unit-agnostic in the sense that it does no internal
//! conversion: it works in any *consistent* set of units. The natural
//! SI choice used throughout the docs and tests is lengths (`b`, `d`,
//! `a`) in millimetres (mm), stresses (`fc'`, `fy`) in megapascals
//! (MPa = N/mm^2), and areas (`As`) in square millimetres (mm^2),
//! which makes [`BeamSection::nominal_moment`] come out in newton-
//! millimetres (N·mm). Divide by `1.0e6` for kN·m. US customary units
//! (in, ksi, in^2 -> kip·in) work equally well as long as every input
//! shares the same system.
//!
//! ## Equations
//!
//! For tension steel area `As` at yield stress `fy`, the tensile force
//! `T = As*fy` is balanced by the concrete compression block of width
//! `b` and uniform stress `0.85*fc'`. Equating `C = T` gives the
//! stress-block depth
//!
//! ```text
//! a = As * fy / (0.85 * fc' * b)
//! ```
//!
//! The nominal moment is the tensile force times the internal lever
//! arm `d - a/2`:
//!
//! ```text
//! Mn = As * fy * (d - a/2)
//! ```
//!
//! The reinforcement (steel) ratio is
//!
//! ```text
//! rho = As / (b * d)
//! ```
//!
//! and the design strength applies the strength-reduction factor
//! `phi` (ACI-318 uses `phi = 0.90` for tension-controlled flexure):
//!
//! ```text
//! phi_Mn = phi * Mn
//! ```

use serde::{Deserialize, Serialize};

use crate::error::RcBeamError;

/// ACI equivalent-stress-block intensity coefficient: the uniform
/// Whitney compressive stress is `0.85 * fc'`.
pub const STRESS_BLOCK_INTENSITY: f64 = 0.85;

/// Default ACI-318 strength-reduction factor for a tension-controlled
/// (under-reinforced) flexural member.
pub const PHI_TENSION_CONTROLLED: f64 = 0.90;

/// A singly-reinforced rectangular reinforced-concrete beam section.
///
/// The fields are the four geometric / material quantities plus the
/// tension-steel area needed to evaluate flexural capacity. Construct
/// with [`BeamSection::new`], which validates that every value is a
/// finite positive number.
///
/// In the natural SI unit system (see the [module docs](crate::beam)):
/// `b`, `d` in mm, `fc_prime`, `fy` in MPa, `area_steel` in mm^2.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct BeamSection {
    /// Section width `b` (compression-face width).
    pub width: f64,
    /// Effective depth `d` — distance from the extreme compression
    /// fibre to the centroid of the tension reinforcement.
    pub effective_depth: f64,
    /// Specified concrete compressive strength `fc'`.
    pub fc_prime: f64,
    /// Specified steel yield strength `fy`.
    pub fy: f64,
    /// Total tension-reinforcement area `As`.
    pub area_steel: f64,
}

impl BeamSection {
    /// Build a validated [`BeamSection`].
    ///
    /// # Arguments
    ///
    /// All five quantities must be finite and strictly positive, in a
    /// single consistent unit system: `width` is the section width `b`,
    /// `effective_depth` the effective depth `d`, `fc_prime` the
    /// concrete strength `fc'`, `fy` the steel yield strength, and
    /// `area_steel` the tension-steel area `As`.
    ///
    /// # Errors
    ///
    /// Returns [`RcBeamError::NonFinite`] if any argument is `NaN` or
    /// infinite, or [`RcBeamError::NonPositive`] if any argument is
    /// `<= 0`.
    pub fn new(
        width: f64,
        effective_depth: f64,
        fc_prime: f64,
        fy: f64,
        area_steel: f64,
    ) -> Result<Self, RcBeamError> {
        let width = RcBeamError::require_positive("b", width)?;
        let effective_depth = RcBeamError::require_positive("d", effective_depth)?;
        let fc_prime = RcBeamError::require_positive("fc", fc_prime)?;
        let fy = RcBeamError::require_positive("fy", fy)?;
        let area_steel = RcBeamError::require_positive("As", area_steel)?;
        Ok(Self {
            width,
            effective_depth,
            fc_prime,
            fy,
            area_steel,
        })
    }

    /// Depth of the Whitney equivalent rectangular stress block,
    /// `a = As*fy / (0.85*fc'*b)`.
    ///
    /// Derived from horizontal force equilibrium `C = T`, with concrete
    /// compression `C = 0.85*fc'*b*a` and steel tension `T = As*fy`.
    /// Because all of `As`, `fy`, `fc'`, `b` are positive (enforced at
    /// construction), `a` is always a finite positive number, so this
    /// method is infallible.
    pub fn stress_block_depth(&self) -> f64 {
        self.area_steel * self.fy / (STRESS_BLOCK_INTENSITY * self.fc_prime * self.width)
    }

    /// Internal lever arm `jd = d - a/2` between the resultant
    /// compressive and tensile forces.
    ///
    /// # Errors
    ///
    /// Returns [`RcBeamError::StressBlockExceedsDepth`] when the stress
    /// block reaches or passes the effective depth (`a >= d`), which
    /// makes the lever arm non-positive and the section degenerate
    /// (grossly over-reinforced — far outside the model's valid,
    /// tension-controlled regime).
    pub fn lever_arm(&self) -> Result<f64, RcBeamError> {
        let a = self.stress_block_depth();
        if a >= self.effective_depth {
            return Err(RcBeamError::StressBlockExceedsDepth {
                a,
                d: self.effective_depth,
            });
        }
        Ok(self.effective_depth - a / 2.0)
    }

    /// Nominal flexural strength `Mn = As*fy*(d - a/2)`.
    ///
    /// In SI inputs (mm, MPa, mm^2) the result is in N·mm.
    ///
    /// # Errors
    ///
    /// Propagates [`RcBeamError::StressBlockExceedsDepth`] from
    /// [`lever_arm`](Self::lever_arm) for a degenerate section.
    pub fn nominal_moment(&self) -> Result<f64, RcBeamError> {
        let jd = self.lever_arm()?;
        Ok(self.area_steel * self.fy * jd)
    }

    /// Design flexural strength `phi*Mn` using the supplied
    /// strength-reduction factor `phi`.
    ///
    /// For a tension-controlled section ACI-318 uses `phi = 0.90`;
    /// [`design_moment_default`](Self::design_moment_default) wires that
    /// in for you.
    ///
    /// # Errors
    ///
    /// Returns [`RcBeamError::PhiOutOfRange`] when `phi` is not in
    /// `(0, 1]` (including `NaN`/`±∞`), and propagates
    /// [`RcBeamError::StressBlockExceedsDepth`] from
    /// [`nominal_moment`](Self::nominal_moment).
    pub fn design_moment(&self, phi: f64) -> Result<f64, RcBeamError> {
        if !phi.is_finite() || phi <= 0.0 || phi > 1.0 {
            return Err(RcBeamError::PhiOutOfRange { value: phi });
        }
        Ok(phi * self.nominal_moment()?)
    }

    /// Design flexural strength using the default tension-controlled
    /// factor [`PHI_TENSION_CONTROLLED`] (`0.90`).
    ///
    /// # Errors
    ///
    /// Propagates [`RcBeamError::StressBlockExceedsDepth`] from
    /// [`nominal_moment`](Self::nominal_moment).
    pub fn design_moment_default(&self) -> Result<f64, RcBeamError> {
        self.design_moment(PHI_TENSION_CONTROLLED)
    }

    /// Reinforcement (steel) ratio `rho = As / (b*d)`.
    ///
    /// Dimensionless. Always finite and positive because `As`, `b`, `d`
    /// are positive, so this method is infallible.
    pub fn reinforcement_ratio(&self) -> f64 {
        self.area_steel / (self.width * self.effective_depth)
    }

    /// Balanced reinforcement ratio `rho_b` for the section, the steel
    /// ratio at which concrete crushing and steel yielding occur
    /// simultaneously:
    ///
    /// ```text
    /// rho_b = 0.85 * beta1 * (fc' / fy) * (eps_cu*Es / (eps_cu*Es + fy))
    /// ```
    ///
    /// using the ultimate concrete strain `eps_cu = 0.003` and the
    /// steel elastic modulus `Es`. A section with `rho < rho_b` is
    /// *under-reinforced* (ductile, steel yields first); see
    /// [`is_under_reinforced`](Self::is_under_reinforced).
    ///
    /// # Arguments
    ///
    /// `beta1` is the stress-block depth factor `beta1 = a/c` (ACI-318:
    /// `0.85` for `fc' <= 28 MPa`, decreasing for higher strengths),
    /// which must be finite and in `(0, 1]`. `steel_modulus` is the
    /// steel elastic modulus `Es` (MPa; ~200_000 for mild steel), which
    /// must be finite and positive.
    ///
    /// # Errors
    ///
    /// Returns [`RcBeamError::PhiOutOfRange`] (reusing the `(0, 1]`
    /// range guard) if `beta1` is out of range, and
    /// [`RcBeamError::NonPositive`] / [`RcBeamError::NonFinite`] if
    /// `steel_modulus` is not finite and positive.
    pub fn balanced_ratio(&self, beta1: f64, steel_modulus: f64) -> Result<f64, RcBeamError> {
        if !beta1.is_finite() || beta1 <= 0.0 || beta1 > 1.0 {
            return Err(RcBeamError::PhiOutOfRange { value: beta1 });
        }
        let es = RcBeamError::require_positive("Es", steel_modulus)?;
        let eps_cu = 0.003_f64;
        let strain_term = eps_cu * es / (eps_cu * es + self.fy);
        Ok(STRESS_BLOCK_INTENSITY * beta1 * (self.fc_prime / self.fy) * strain_term)
    }

    /// Whether the section is *under-reinforced* (ductile): the actual
    /// steel ratio is strictly below the balanced ratio,
    /// `rho < rho_b`.
    ///
    /// Code-conforming flexural members are designed under-reinforced
    /// so that the steel yields — giving warning by deflection and
    /// cracking — well before the concrete crushes.
    ///
    /// # Arguments
    ///
    /// See [`balanced_ratio`](Self::balanced_ratio) for `beta1` and
    /// `steel_modulus`.
    ///
    /// # Errors
    ///
    /// Propagates the validation errors of
    /// [`balanced_ratio`](Self::balanced_ratio).
    pub fn is_under_reinforced(&self, beta1: f64, steel_modulus: f64) -> Result<bool, RcBeamError> {
        Ok(self.reinforcement_ratio() < self.balanced_ratio(beta1, steel_modulus)?)
    }

    /// The tension-steel area `As` required to reach a target **nominal**
    /// moment `Mn` — the design inverse of [`nominal_moment`](Self::nominal_moment).
    ///
    /// Substituting the stress-block depth `a = As*fy/(0.85*fc'*b)` into
    /// `Mn = As*fy*(d - a/2)` makes `Mn` a quadratic in `As`; the
    /// physically meaningful (under-reinforced) root is
    ///
    /// ```text
    /// As = (0.85*fc'*b / fy) * (d - sqrt(d^2 - 2*Mn/(0.85*fc'*b)))
    /// ```
    ///
    /// the standard ACI singly-reinforced sizing formula. `target_moment`
    /// is the nominal moment in the same consistent unit system as the
    /// other arguments (SI: `b`, `d` in mm, `fc'`, `fy` in MPa, `Mn` in
    /// N·mm, giving `As` in mm^2). To size from a *design* moment `Mu`,
    /// pass `Mu / phi`.
    ///
    /// # Errors
    ///
    /// Returns [`RcBeamError::NonFinite`] / [`RcBeamError::NonPositive`]
    /// for any non-finite or non-positive argument, and
    /// [`RcBeamError::MomentExceedsCapacity`] when `target_moment` is
    /// larger than the section can carry singly-reinforced (the
    /// discriminant goes negative at `Mn > 0.85*fc'*b*d^2/2`, the `a = d`
    /// limit).
    pub fn required_steel_area(
        width: f64,
        effective_depth: f64,
        fc_prime: f64,
        fy: f64,
        target_moment: f64,
    ) -> Result<f64, RcBeamError> {
        let b = RcBeamError::require_positive("b", width)?;
        let d = RcBeamError::require_positive("d", effective_depth)?;
        let fc = RcBeamError::require_positive("fc", fc_prime)?;
        let fy = RcBeamError::require_positive("fy", fy)?;
        let mn = RcBeamError::require_positive("target_moment", target_moment)?;
        let cap = STRESS_BLOCK_INTENSITY * fc * b; // 0.85 * fc' * b
        let discriminant = d * d - 2.0 * mn / cap;
        if discriminant < 0.0 {
            return Err(RcBeamError::MomentExceedsCapacity {
                target_moment: mn,
                max_moment: cap * d * d / 2.0,
            });
        }
        Ok((cap / fy) * (d - discriminant.sqrt()))
    }

    /// Build the [`BeamSection`] whose tension steel is sized for a target
    /// **nominal** moment `Mn`, via [`required_steel_area`](Self::required_steel_area).
    ///
    /// The returned section reproduces the target: its
    /// [`nominal_moment`](Self::nominal_moment) equals `target_moment`
    /// (to floating-point tolerance), closing the design loop.
    ///
    /// # Errors
    ///
    /// Propagates the errors of
    /// [`required_steel_area`](Self::required_steel_area) and
    /// [`new`](Self::new).
    pub fn for_nominal_moment(
        width: f64,
        effective_depth: f64,
        fc_prime: f64,
        fy: f64,
        target_moment: f64,
    ) -> Result<Self, RcBeamError> {
        let area_steel =
            Self::required_steel_area(width, effective_depth, fc_prime, fy, target_moment)?;
        Self::new(width, effective_depth, fc_prime, fy, area_steel)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Absolute tolerance for floating-point capacity comparisons
    /// (N·mm scale ~1e8, so 1e-3 is ~1e-11 relative).
    const EPS: f64 = 1.0e-3;
    /// Tighter tolerance for dimensionless ratios.
    const EPS_RATIO: f64 = 1.0e-9;

    /// Canonical worked example used across the suite (SI units):
    /// b = 300 mm, d = 500 mm, fc' = 30 MPa, fy = 420 MPa,
    /// As = 1500 mm^2.
    fn worked() -> BeamSection {
        BeamSection::new(300.0, 500.0, 30.0, 420.0, 1500.0).unwrap()
    }

    #[test]
    fn constructor_rejects_bad_inputs() {
        assert!(BeamSection::new(0.0, 500.0, 30.0, 420.0, 1500.0).is_err());
        assert!(BeamSection::new(300.0, -1.0, 30.0, 420.0, 1500.0).is_err());
        assert!(BeamSection::new(300.0, 500.0, 0.0, 420.0, 1500.0).is_err());
        assert!(BeamSection::new(300.0, 500.0, 30.0, f64::NAN, 1500.0).is_err());
        assert!(BeamSection::new(300.0, 500.0, 30.0, 420.0, 1500.0).is_ok());
    }

    /// Ground truth for the stress-block depth.
    /// a = As*fy / (0.85*fc'*b)
    ///   = 1500*420 / (0.85*30*300)
    ///   = 630_000 / 7_650 = 82.35294117... mm
    #[test]
    fn stress_block_depth_matches_closed_form() {
        let beam = worked();
        let a = beam.stress_block_depth();
        let expected = 630_000.0 / 7_650.0;
        assert!((a - expected).abs() < EPS, "a = {a}, expected {expected}");
        // Independent recomputation straight from the definition.
        let manual = 1500.0 * 420.0 / (0.85 * 30.0 * 300.0);
        assert!((a - manual).abs() < EPS, "a = {a}, manual {manual}");
    }

    /// Ground truth for the lever arm jd = d - a/2.
    /// a = 82.352941..., a/2 = 41.176470..., jd = 458.823529... mm
    #[test]
    fn lever_arm_matches_closed_form() {
        let beam = worked();
        let a = beam.stress_block_depth();
        let jd = beam.lever_arm().unwrap();
        let expected = 500.0 - a / 2.0;
        assert!(
            (jd - expected).abs() < EPS,
            "jd = {jd}, expected {expected}"
        );
        assert!((jd - 458.823_529_411_764_7).abs() < EPS, "jd = {jd}");
    }

    /// Ground truth for nominal moment Mn = As*fy*(d - a/2).
    /// Mn = 1500*420*458.823529... = 289_058_823.5... N·mm
    #[test]
    fn nominal_moment_matches_closed_form() {
        let beam = worked();
        let mn = beam.nominal_moment().unwrap();
        let a = beam.stress_block_depth();
        let expected = 1500.0 * 420.0 * (500.0 - a / 2.0);
        assert!(
            (mn - expected).abs() < EPS,
            "Mn = {mn}, expected {expected}"
        );
        // ~289.06 kN·m.
        assert!(
            (mn / 1.0e6 - 289.058_823_5).abs() < 1.0e-3,
            "Mn = {mn} N·mm"
        );
    }

    /// rho = As / (b*d) = 1500 / (300*500) = 1500/150000 = 0.01.
    #[test]
    fn reinforcement_ratio_matches_closed_form() {
        let beam = worked();
        let rho = beam.reinforcement_ratio();
        assert!((rho - 0.01).abs() < EPS_RATIO, "rho = {rho}");
        let manual = 1500.0 / (300.0 * 500.0);
        assert!((rho - manual).abs() < EPS_RATIO, "rho = {rho}");
    }

    /// phi*Mn with the default phi = 0.90.
    /// phi*Mn = 0.90 * 289_058_823.5... = 260_152_941.1... N·mm.
    #[test]
    fn design_moment_applies_phi() {
        let beam = worked();
        let mn = beam.nominal_moment().unwrap();
        let phi_mn = beam.design_moment_default().unwrap();
        assert!(
            (phi_mn - 0.90 * mn).abs() < EPS,
            "phi*Mn = {phi_mn}, 0.9*Mn = {}",
            0.90 * mn
        );
        // Design strength is strictly below the nominal strength.
        assert!(phi_mn < mn, "phi*Mn {phi_mn} should be < Mn {mn}");
        // Custom phi scales linearly.
        let half = beam.design_moment(0.5).unwrap();
        assert!((half - 0.5 * mn).abs() < EPS, "0.5*Mn check: {half}");
    }

    #[test]
    fn design_moment_rejects_bad_phi() {
        let beam = worked();
        assert!(beam.design_moment(0.0).is_err());
        assert!(beam.design_moment(-0.1).is_err());
        assert!(beam.design_moment(1.5).is_err());
        assert!(beam.design_moment(f64::NAN).is_err());
        // phi = 1.0 is the boundary and is accepted.
        let full = beam.design_moment(1.0).unwrap();
        let mn = beam.nominal_moment().unwrap();
        assert!((full - mn).abs() < EPS, "phi=1 should give Mn");
    }

    /// Monotonicity: more tension steel -> larger nominal moment.
    /// (For under-reinforced sections, increasing As raises both the
    /// tensile force and the stress block; the moment still rises until
    /// the section becomes over-reinforced.)
    #[test]
    fn more_steel_gives_higher_moment() {
        let base = BeamSection::new(300.0, 500.0, 30.0, 420.0, 1000.0).unwrap();
        let more = BeamSection::new(300.0, 500.0, 30.0, 420.0, 2000.0).unwrap();
        let mn_base = base.nominal_moment().unwrap();
        let mn_more = more.nominal_moment().unwrap();
        assert!(
            mn_more > mn_base,
            "Mn(As=2000)={mn_more} should exceed Mn(As=1000)={mn_base}"
        );
        // And a fine sweep is monotone increasing.
        let mut prev = 0.0_f64;
        for as_steel in (500..=2500).step_by(250) {
            let s = BeamSection::new(300.0, 500.0, 30.0, 420.0, as_steel as f64).unwrap();
            let mn = s.nominal_moment().unwrap();
            assert!(mn > prev, "non-monotone at As = {as_steel}: {mn} <= {prev}");
            prev = mn;
        }
    }

    /// Balanced ratio ground truth.
    /// rho_b = 0.85*beta1*(fc'/fy)*(eps_cu*Es/(eps_cu*Es + fy))
    /// with beta1 = 0.85, Es = 200_000, eps_cu = 0.003, fc' = 30,
    /// fy = 420:
    ///   eps_cu*Es = 600
    ///   strain_term = 600 / (600 + 420) = 600/1020 = 0.588235...
    ///   rho_b = 0.85*0.85*(30/420)*0.588235...
    ///         = 0.7225 * 0.0714285... * 0.588235...
    ///         = 0.030364... -> ~0.03036
    #[test]
    fn balanced_ratio_matches_closed_form() {
        let beam = worked();
        let rho_b = beam.balanced_ratio(0.85, 200_000.0).unwrap();
        let strain_term = 600.0 / 1020.0;
        let expected = 0.85 * 0.85 * (30.0 / 420.0) * strain_term;
        assert!(
            (rho_b - expected).abs() < EPS_RATIO,
            "rho_b = {rho_b}, expected {expected}"
        );
        assert!((rho_b - 0.030_364_0).abs() < 1.0e-5, "rho_b = {rho_b}");
    }

    #[test]
    fn balanced_ratio_rejects_bad_inputs() {
        let beam = worked();
        assert!(beam.balanced_ratio(0.0, 200_000.0).is_err());
        assert!(beam.balanced_ratio(1.5, 200_000.0).is_err());
        assert!(beam.balanced_ratio(0.85, 0.0).is_err());
        assert!(beam.balanced_ratio(0.85, f64::INFINITY).is_err());
    }

    /// Under-reinforced check: the worked section (rho = 0.01) sits
    /// below the balanced ratio (~0.0304), so it is under-reinforced.
    #[test]
    fn under_reinforced_check() {
        let beam = worked();
        assert!((beam.reinforcement_ratio() - 0.01).abs() < EPS_RATIO);
        assert!(
            beam.is_under_reinforced(0.85, 200_000.0).unwrap(),
            "rho=0.01 < rho_b~0.0304 must be under-reinforced"
        );

        // Pile in steel until rho exceeds the balanced ratio: now the
        // section is NOT under-reinforced.
        let heavy = BeamSection::new(300.0, 500.0, 30.0, 420.0, 5000.0).unwrap();
        // rho = 5000/150000 = 0.0333... > 0.0304.
        assert!((heavy.reinforcement_ratio() - 0.033_333_333).abs() < 1.0e-6);
        assert!(
            !heavy.is_under_reinforced(0.85, 200_000.0).unwrap(),
            "rho=0.0333 > rho_b~0.0304 must NOT be under-reinforced"
        );
    }

    // ----- Required steel area (design inverse) -----------------------

    /// Sizing the steel for the worked section's own nominal moment
    /// recovers its tension-steel area (As = 1500 mm^2): a forward ->
    /// inverse round trip.
    #[test]
    fn required_steel_area_recovers_worked_section() {
        let beam = worked();
        let mn = beam.nominal_moment().unwrap();
        let as_req = BeamSection::required_steel_area(300.0, 500.0, 30.0, 420.0, mn).unwrap();
        assert!((as_req - 1500.0).abs() < 1e-6, "As = {as_req}");
    }

    /// Independent As -> Mn -> As round trip on a different section.
    #[test]
    fn required_steel_area_round_trips_independent_section() {
        let beam = BeamSection::new(250.0, 450.0, 25.0, 400.0, 2000.0).unwrap();
        let mn = beam.nominal_moment().unwrap();
        let as_req = BeamSection::required_steel_area(250.0, 450.0, 25.0, 400.0, mn).unwrap();
        assert!((as_req - 2000.0).abs() < 1e-6, "As = {as_req}");
    }

    /// Closed-form ground truth plus the design-loop closure: the steel
    /// area matches the ACI formula and a section built with it
    /// reproduces the target moment.
    #[test]
    fn required_steel_area_hand_value_and_closure() {
        // b=250, d=450, fc'=25, fy=400, target Mn = 150e6 N·mm.
        let as_req = BeamSection::required_steel_area(250.0, 450.0, 25.0, 400.0, 150.0e6).unwrap();
        let cap = 0.85 * 25.0 * 250.0; // 5312.5
        let disc = 450.0_f64.powi(2) - 2.0 * 150.0e6 / cap;
        let expected = (cap / 400.0) * (450.0 - disc.sqrt());
        assert!(
            (as_req - expected).abs() < 1e-9,
            "As = {as_req}, expected {expected}"
        );
        // Ballpark: ~901 mm^2.
        assert!((as_req - 901.4).abs() < 0.5, "As = {as_req}");
        // Closure: a section sized this way carries exactly the target.
        let sized = BeamSection::for_nominal_moment(250.0, 450.0, 25.0, 400.0, 150.0e6).unwrap();
        assert!((sized.area_steel - as_req).abs() < 1e-9);
        assert!((sized.nominal_moment().unwrap() - 150.0e6).abs() < EPS);
    }

    /// More required moment needs more tension steel.
    #[test]
    fn required_steel_area_grows_with_target() {
        let lo = BeamSection::required_steel_area(300.0, 500.0, 30.0, 420.0, 100.0e6).unwrap();
        let hi = BeamSection::required_steel_area(300.0, 500.0, 30.0, 420.0, 300.0e6).unwrap();
        assert!(hi > lo, "more moment should need more steel: {hi} vs {lo}");
    }

    /// A target beyond the singly-reinforced limit is rejected, while one
    /// just below it succeeds.
    #[test]
    fn required_steel_area_rejects_overlarge_target() {
        // Max singly-reinforced Mn = 0.85*fc'*b*d^2/2.
        let max = 0.85 * 30.0 * 300.0 * 500.0_f64.powi(2) / 2.0;
        let err =
            BeamSection::required_steel_area(300.0, 500.0, 30.0, 420.0, max * 1.01).unwrap_err();
        assert_eq!(err.code(), "rcbeam.moment-exceeds-capacity");
        assert!(BeamSection::required_steel_area(300.0, 500.0, 30.0, 420.0, max * 0.99).is_ok());
    }

    #[test]
    fn required_steel_area_rejects_bad_inputs() {
        assert!(BeamSection::required_steel_area(0.0, 500.0, 30.0, 420.0, 1.0e8).is_err());
        assert!(BeamSection::required_steel_area(300.0, -1.0, 30.0, 420.0, 1.0e8).is_err());
        assert!(BeamSection::required_steel_area(300.0, 500.0, f64::NAN, 420.0, 1.0e8).is_err());
        assert!(BeamSection::required_steel_area(300.0, 500.0, 30.0, 0.0, 1.0e8).is_err());
        assert!(BeamSection::required_steel_area(300.0, 500.0, 30.0, 420.0, 0.0).is_err());
        assert!(
            BeamSection::required_steel_area(300.0, 500.0, 30.0, 420.0, f64::INFINITY).is_err()
        );
    }

    /// A grossly over-reinforced section drives the stress block past
    /// the effective depth and the capacity calls error out cleanly
    /// rather than returning a bogus negative lever arm.
    #[test]
    fn degenerate_section_errors() {
        // As huge relative to a shallow, weak, narrow section.
        let beam = BeamSection::new(100.0, 100.0, 20.0, 500.0, 8000.0).unwrap();
        // a = 8000*500 / (0.85*20*100) = 4_000_000 / 1700 = 2352.9 mm
        //   >> d = 100 mm.
        let a = beam.stress_block_depth();
        assert!(a > beam.effective_depth, "expected a {a} > d");
        let err = beam.nominal_moment().unwrap_err();
        assert_eq!(err.code(), "rcbeam.stress-block-exceeds-depth");
        assert!(beam.lever_arm().is_err());
        assert!(beam.design_moment_default().is_err());
    }

    /// Round-trip a section through serde JSON.
    #[test]
    fn section_serde_roundtrip() {
        let beam = worked();
        let json = serde_json::to_string(&beam).unwrap();
        let back: BeamSection = serde_json::from_str(&json).unwrap();
        assert_eq!(beam, back);
        let mn_a = beam.nominal_moment().unwrap();
        let mn_b = back.nominal_moment().unwrap();
        assert!((mn_a - mn_b).abs() < EPS);
    }
}
