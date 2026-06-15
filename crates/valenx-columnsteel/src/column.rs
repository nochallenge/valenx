//! Euler-Johnson steel-column buckling.
//!
//! This module implements the classic two-regime column-buckling model
//! used throughout the AISC Allowable Stress Design (ASD) tradition and
//! every undergraduate machine-design / steel-structures text:
//!
//! 1. Compute the **slenderness ratio** `λ = K L / r`.
//! 2. Compute the **column-slenderness transition**
//!    `Cc = sqrt(2 π² E / Fy)`, the slenderness at which Euler buckling
//!    and yielding coincide.
//! 3. For **long** columns (`λ ≥ Cc`) the elastic **Euler** critical
//!    stress governs: `Fcr = π² E / λ²`.
//! 4. For **short / intermediate** columns (`λ < Cc`) the inelastic
//!    **Johnson parabola** governs: `Fcr = Fy [1 − λ² / (2 Cc²)]`.
//! 5. The two curves are *tangent* at `λ = Cc`, where both equal
//!    `Fy / 2`.
//!
//! All stresses are returned in the same units as the supplied modulus
//! and yield stress (the model is dimensionless in the consistent-unit
//! sense — pass psi and get psi, pass MPa and get MPa). Lengths and the
//! radius of gyration only ever appear as the ratio `K L / r`, so any
//! consistent length unit works.

use crate::error::{ColumnError, Result};
use serde::{Deserialize, Serialize};

/// Pi to `f64` precision, named locally so the formulas read like the
/// textbook expressions.
const PI: f64 = std::f64::consts::PI;

/// Which branch of the Euler-Johnson curve governs a column.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash, Serialize, Deserialize)]
pub enum Regime {
    /// Short / intermediate column (`λ < Cc`): the inelastic Johnson
    /// parabola governs.
    Johnson,
    /// Long column (`λ ≥ Cc`): the elastic Euler hyperbola governs.
    Euler,
}

impl Regime {
    /// Stable lowercase identifier (`"johnson"` / `"euler"`).
    pub fn as_str(self) -> &'static str {
        match self {
            Regime::Johnson => "johnson",
            Regime::Euler => "euler",
        }
    }
}

/// Material and geometry of one prismatic steel column, expressed in a
/// single consistent unit system.
///
/// The struct stores only validated, strictly-positive material
/// properties plus the effective slenderness ratio; constructing it
/// through [`Column::new`] is the single validation gate for the whole
/// module.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Column {
    /// Young's modulus `E` (e.g. 29 000 000 psi or 200 000 MPa for
    /// structural steel). Strictly positive.
    youngs_modulus: f64,
    /// Yield stress `Fy` (e.g. 36 000 psi or 250 MPa). Strictly positive.
    yield_stress: f64,
    /// Effective slenderness ratio `λ = K L / r`. Non-negative.
    slenderness: f64,
}

impl Column {
    /// Build a column directly from a known slenderness ratio.
    ///
    /// # Arguments
    ///
    /// `youngs_modulus` and `yield_stress` must be strictly positive and
    /// finite; `slenderness` (`= K L / r`) must be finite and
    /// non-negative.
    ///
    /// # Errors
    ///
    /// Returns [`ColumnError::Invalid`] naming the offending parameter
    /// when any argument is out of range.
    pub fn new(youngs_modulus: f64, yield_stress: f64, slenderness: f64) -> Result<Self> {
        Ok(Column {
            youngs_modulus: ColumnError::require_positive(youngs_modulus, "youngs_modulus")?,
            yield_stress: ColumnError::require_positive(yield_stress, "yield_stress")?,
            slenderness: ColumnError::require_non_negative(slenderness, "slenderness")?,
        })
    }

    /// Build a column from physical geometry, computing the slenderness
    /// ratio `λ = K L / r` from the effective-length factor `K`, the
    /// unbraced length `L`, and the least radius of gyration `r`.
    ///
    /// # Arguments
    ///
    /// `youngs_modulus`, `yield_stress`, `effective_length_factor` (`K`)
    /// and `radius_of_gyration` (`r`) must be strictly positive and
    /// finite; `unbraced_length` (`L`) must be finite and non-negative.
    /// `L` and `r` may use any length unit so long as they share it.
    ///
    /// # Errors
    ///
    /// Returns [`ColumnError::Invalid`] naming the offending parameter.
    pub fn from_geometry(
        youngs_modulus: f64,
        yield_stress: f64,
        effective_length_factor: f64,
        unbraced_length: f64,
        radius_of_gyration: f64,
    ) -> Result<Self> {
        let k = ColumnError::require_positive(effective_length_factor, "effective_length_factor")?;
        let l = ColumnError::require_non_negative(unbraced_length, "unbraced_length")?;
        let r = ColumnError::require_positive(radius_of_gyration, "radius_of_gyration")?;
        Column::new(youngs_modulus, yield_stress, k * l / r)
    }

    /// Young's modulus `E`.
    pub fn youngs_modulus(&self) -> f64 {
        self.youngs_modulus
    }

    /// Yield stress `Fy`.
    pub fn yield_stress(&self) -> f64 {
        self.yield_stress
    }

    /// Slenderness ratio `λ = K L / r`.
    pub fn slenderness(&self) -> f64 {
        self.slenderness
    }

    /// The column-slenderness transition `Cc = sqrt(2 π² E / Fy)`.
    ///
    /// This is the slenderness at which the elastic Euler stress equals
    /// `Fy / 2`; it separates the inelastic (Johnson) and elastic
    /// (Euler) regimes. Always strictly positive for a valid column.
    pub fn cc(&self) -> f64 {
        (2.0 * PI * PI * self.youngs_modulus / self.yield_stress).sqrt()
    }

    /// Whether this column is in the [`Regime::Euler`] (long) or
    /// [`Regime::Johnson`] (short / intermediate) branch.
    ///
    /// The boundary `λ = Cc` is assigned to [`Regime::Euler`]; at that
    /// exact point both formulas return the same value (`Fy / 2`), so
    /// the choice is immaterial to the critical stress.
    pub fn regime(&self) -> Regime {
        if self.slenderness >= self.cc() {
            Regime::Euler
        } else {
            Regime::Johnson
        }
    }

    /// The **Euler** elastic critical buckling stress
    /// `Fcr = π² E / λ²`, valid for long columns (`λ ≥ Cc`).
    ///
    /// # Errors
    ///
    /// Returns [`ColumnError::Invalid`] if the slenderness is zero (the
    /// Euler stress diverges as `λ → 0`), or [`ColumnError::OutOfRange`]
    /// if `λ < Cc` — the column is then in the inelastic Johnson regime,
    /// where this elastic formula over-predicts the real strength.
    pub fn euler_stress(&self) -> Result<f64> {
        if self.slenderness == 0.0 {
            return Err(ColumnError::invalid(
                "slenderness",
                "Euler stress is undefined (infinite) at zero slenderness",
            ));
        }
        if self.slenderness < self.cc() {
            return Err(ColumnError::out_of_range(
                self.slenderness,
                Regime::Euler.as_str(),
                "slenderness is below Cc; the Johnson parabola governs here",
            ));
        }
        Ok(self.euler_stress_unchecked())
    }

    /// The Euler stress `π² E / λ²` without the regime / zero guard.
    ///
    /// Used internally and by callers who want the bare elastic curve
    /// (for plotting the full hyperbola, say). Returns `+∞` at zero
    /// slenderness. Document the domain yourself when you call this.
    pub fn euler_stress_unchecked(&self) -> f64 {
        PI * PI * self.youngs_modulus / (self.slenderness * self.slenderness)
    }

    /// The **Johnson** inelastic parabolic critical buckling stress
    /// `Fcr = Fy [1 − λ² / (2 Cc²)]`, valid for short / intermediate
    /// columns (`λ ≤ Cc`).
    ///
    /// # Errors
    ///
    /// Returns [`ColumnError::OutOfRange`] if `λ > Cc`, where the
    /// parabola would drop below `Fy / 2` and the elastic Euler curve
    /// is the correct (higher-fidelity, lower) bound.
    pub fn johnson_stress(&self) -> Result<f64> {
        if self.slenderness > self.cc() {
            return Err(ColumnError::out_of_range(
                self.slenderness,
                Regime::Johnson.as_str(),
                "slenderness exceeds Cc; the Euler hyperbola governs here",
            ));
        }
        Ok(self.johnson_stress_unchecked())
    }

    /// The Johnson parabola `Fy [1 − λ² / (2 Cc²)]` without the regime
    /// guard.
    ///
    /// Used internally and for plotting the full parabola. Beyond `Cc`
    /// it dips below `Fy / 2` and eventually crosses zero, so document
    /// the domain yourself when you call this directly.
    pub fn johnson_stress_unchecked(&self) -> f64 {
        let cc = self.cc();
        let ratio = self.slenderness / cc;
        self.yield_stress * (1.0 - 0.5 * ratio * ratio)
    }

    /// The governing **critical buckling stress** `Fcr`, automatically
    /// selecting the Euler branch for long columns and the Johnson
    /// branch for short / intermediate columns.
    ///
    /// The result is continuous across `λ = Cc` (both branches meet at
    /// `Fy / 2`), never exceeds the yield stress, and decreases
    /// monotonically as slenderness increases. At zero slenderness it
    /// equals `Fy` (the Johnson parabola's intercept) — a squash load,
    /// not an elastic buckling load.
    ///
    /// This method is total over every validly-constructed [`Column`]:
    /// it never returns an error.
    pub fn critical_stress(&self) -> f64 {
        match self.regime() {
            Regime::Euler => self.euler_stress_unchecked(),
            Regime::Johnson => self.johnson_stress_unchecked(),
        }
    }

    /// The AISC-ASD **variable factor of safety**.
    ///
    /// In the inelastic (Johnson) regime the factor of safety ramps from
    /// `5/3 ≈ 1.667` at `λ = 0` up to `23/12 ≈ 1.917` at `λ = Cc`:
    ///
    /// `FS = 5/3 + (3/8)(λ/Cc) − (1/8)(λ/Cc)³`.
    ///
    /// In the elastic (Euler) regime it is the constant `23/12`.
    /// This is the historical AISC ASD safety factor (the `Fa`
    /// denominators of the 9th-edition ASD specification).
    pub fn factor_of_safety_aisc(&self) -> f64 {
        match self.regime() {
            Regime::Euler => 23.0 / 12.0,
            Regime::Johnson => {
                let x = self.slenderness / self.cc();
                5.0 / 3.0 + (3.0 / 8.0) * x - (1.0 / 8.0) * x * x * x
            }
        }
    }

    /// The **allowable compressive stress** `Fa = Fcr / FS` using the
    /// AISC-ASD variable factor of safety from
    /// [`factor_of_safety_aisc`](Self::factor_of_safety_aisc).
    ///
    /// This is the standard ASD working stress: a member is adequate
    /// when its actual axial stress `P/A` does not exceed `Fa`. Always
    /// strictly less than the critical stress (the factor of safety is
    /// always `> 1`).
    pub fn allowable_stress(&self) -> f64 {
        self.critical_stress() / self.factor_of_safety_aisc()
    }

    /// The **allowable compressive stress** with a caller-supplied
    /// constant factor of safety, `Fa = Fcr / fs`.
    ///
    /// Useful when a code other than AISC ASD (or a teaching example)
    /// prescribes a single safety factor instead of the variable ramp.
    ///
    /// # Errors
    ///
    /// Returns [`ColumnError::Invalid`] if `factor_of_safety` is not
    /// strictly positive and finite.
    pub fn allowable_stress_with(&self, factor_of_safety: f64) -> Result<f64> {
        let fs = ColumnError::require_positive(factor_of_safety, "factor_of_safety")?;
        Ok(self.critical_stress() / fs)
    }

    /// The **allowable axial load** `P_allow = Fa · A` for a gross
    /// cross-sectional area `A`, using the AISC-ASD variable factor of
    /// safety.
    ///
    /// # Errors
    ///
    /// Returns [`ColumnError::Invalid`] if `area` is not strictly
    /// positive and finite.
    pub fn allowable_load(&self, area: f64) -> Result<f64> {
        let a = ColumnError::require_positive(area, "area")?;
        Ok(self.allowable_stress() * a)
    }
}

/// Compute the slenderness ratio `λ = K L / r` directly.
///
/// A free function for callers who only need the ratio and not a full
/// [`Column`].
///
/// # Errors
///
/// Returns [`ColumnError::Invalid`] if `effective_length_factor` or
/// `radius_of_gyration` is not strictly positive, or if
/// `unbraced_length` is negative / non-finite.
pub fn slenderness_ratio(
    effective_length_factor: f64,
    unbraced_length: f64,
    radius_of_gyration: f64,
) -> Result<f64> {
    let k = ColumnError::require_positive(effective_length_factor, "effective_length_factor")?;
    let l = ColumnError::require_non_negative(unbraced_length, "unbraced_length")?;
    let r = ColumnError::require_positive(radius_of_gyration, "radius_of_gyration")?;
    Ok(k * l / r)
}

/// Compute the column-slenderness transition `Cc = sqrt(2 π² E / Fy)`
/// directly from a modulus and yield stress.
///
/// # Errors
///
/// Returns [`ColumnError::Invalid`] if either argument is not strictly
/// positive and finite.
pub fn transition_slenderness(youngs_modulus: f64, yield_stress: f64) -> Result<f64> {
    let e = ColumnError::require_positive(youngs_modulus, "youngs_modulus")?;
    let fy = ColumnError::require_positive(yield_stress, "yield_stress")?;
    Ok((2.0 * PI * PI * e / fy).sqrt())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Absolute tolerance for stress comparisons in psi-scale tests.
    const EPS_PSI: f64 = 1e-3;
    /// Relative tolerance for ratio / dimensionless comparisons.
    const EPS_REL: f64 = 1e-9;

    /// A36 structural steel in US customary units (psi).
    const E_STEEL: f64 = 29_000_000.0;
    const FY_A36: f64 = 36_000.0;

    fn approx(a: f64, b: f64, eps: f64) -> bool {
        (a - b).abs() < eps
    }

    // --- Transition slenderness Cc -----------------------------------

    /// `Cc` for A36 steel is the textbook value ~126.1. AISC tabulates
    /// `Cc = 126.1` for `Fy = 36 ksi`, `E = 29 000 ksi`.
    #[test]
    fn cc_matches_textbook_a36() {
        let cc = transition_slenderness(E_STEEL, FY_A36).unwrap();
        // sqrt(2 * pi^2 * 29e6 / 36e3) = 126.0989...
        assert!(approx(cc, 126.0989, 1e-3), "Cc = {cc}");
    }

    /// `Cc` for Grade-50 steel (`Fy = 50 ksi`) is the tabulated ~107.0.
    #[test]
    fn cc_matches_textbook_grade50() {
        let cc = transition_slenderness(E_STEEL, 50_000.0).unwrap();
        assert!(approx(cc, 107.0, 0.1), "Cc = {cc}");
    }

    /// The `Column::cc` accessor agrees with the free function.
    #[test]
    fn cc_accessor_matches_free_function() {
        let col = Column::new(E_STEEL, FY_A36, 80.0).unwrap();
        let free = transition_slenderness(E_STEEL, FY_A36).unwrap();
        assert!(approx(col.cc(), free, EPS_REL));
    }

    // --- Slenderness ratio -------------------------------------------

    /// `λ = K L / r` is computed exactly. K=1, L=120 in, r=1.5 in -> 80.
    #[test]
    fn slenderness_ratio_is_exact() {
        let lambda = slenderness_ratio(1.0, 120.0, 1.5).unwrap();
        assert!(approx(lambda, 80.0, EPS_REL), "lambda = {lambda}");
    }

    /// The effective-length factor scales the ratio: pinned-fixed K=0.7.
    #[test]
    fn slenderness_ratio_uses_k_factor() {
        let lambda = slenderness_ratio(0.7, 200.0, 2.0).unwrap();
        // 0.7 * 200 / 2 = 70
        assert!(approx(lambda, 70.0, EPS_REL), "lambda = {lambda}");
    }

    /// `from_geometry` and `new(.., slenderness)` agree.
    #[test]
    fn from_geometry_matches_explicit_slenderness() {
        let geo = Column::from_geometry(E_STEEL, FY_A36, 1.0, 120.0, 1.5).unwrap();
        let exp = Column::new(E_STEEL, FY_A36, 80.0).unwrap();
        assert!(approx(geo.slenderness(), exp.slenderness(), EPS_REL));
        assert!(approx(
            geo.critical_stress(),
            exp.critical_stress(),
            EPS_PSI
        ));
    }

    // --- Regime selection: long uses Euler, short uses Johnson -------

    /// A slender column (`λ = 150 > Cc ≈ 126`) is in the Euler regime
    /// and its critical stress equals the bare Euler curve.
    #[test]
    fn long_column_uses_euler() {
        let col = Column::new(E_STEEL, FY_A36, 150.0).unwrap();
        assert_eq!(col.regime(), Regime::Euler);
        let euler = col.euler_stress().unwrap();
        // pi^2 * 29e6 / 150^2 = 12720.8234 psi
        assert!(approx(euler, 12720.8234, 1e-2), "Fe = {euler}");
        assert!(approx(col.critical_stress(), euler, EPS_PSI));
    }

    /// A stocky column (`λ = 50 < Cc`) is in the Johnson regime and its
    /// critical stress equals the Johnson parabola.
    #[test]
    fn short_column_uses_johnson() {
        let col = Column::new(E_STEEL, FY_A36, 50.0).unwrap();
        assert_eq!(col.regime(), Regime::Johnson);
        let johnson = col.johnson_stress().unwrap();
        // Cc = 126.0993; ratio = 50/126.0993 = 0.396512
        // Fcr = 36000 * (1 - 0.5 * 0.396512^2) = 33169.99
        assert!(approx(johnson, 33169.99, 0.1), "Fcr = {johnson}");
        assert!(approx(col.critical_stress(), johnson, EPS_PSI));
    }

    /// Asking for the Euler stress inside the Johnson regime is an error.
    #[test]
    fn euler_rejected_below_cc() {
        let col = Column::new(E_STEEL, FY_A36, 50.0).unwrap();
        let err = col.euler_stress().unwrap_err();
        assert_eq!(err.code(), "columnsteel.out_of_range");
    }

    /// Asking for the Johnson stress beyond the transition is an error.
    #[test]
    fn johnson_rejected_above_cc() {
        let col = Column::new(E_STEEL, FY_A36, 150.0).unwrap();
        let err = col.johnson_stress().unwrap_err();
        assert_eq!(err.code(), "columnsteel.out_of_range");
    }

    // --- Continuity at the transition Cc -----------------------------

    /// At `λ = Cc` the Euler and Johnson curves both equal `Fy / 2`.
    /// This is the defining tangency of the Euler-Johnson construction.
    #[test]
    fn curves_meet_at_cc_equal_half_fy() {
        let cc = transition_slenderness(E_STEEL, FY_A36).unwrap();
        let col = Column::new(E_STEEL, FY_A36, cc).unwrap();

        let euler = col.euler_stress_unchecked();
        let johnson = col.johnson_stress_unchecked();
        let half_fy = FY_A36 / 2.0;

        assert!(approx(euler, half_fy, EPS_PSI), "Euler at Cc = {euler}");
        assert!(
            approx(johnson, half_fy, EPS_PSI),
            "Johnson at Cc = {johnson}"
        );
        assert!(approx(euler, johnson, EPS_PSI), "branches disagree at Cc");
    }

    /// The governing `critical_stress` is continuous across `Cc`: the
    /// value an epsilon below and an epsilon above the transition differ
    /// only by an amount that vanishes with the step. Tests that the
    /// regime switch introduces no jump.
    #[test]
    fn critical_stress_continuous_across_cc() {
        let cc = transition_slenderness(E_STEEL, FY_A36).unwrap();
        let below = Column::new(E_STEEL, FY_A36, cc - 1e-4).unwrap();
        let above = Column::new(E_STEEL, FY_A36, cc + 1e-4).unwrap();
        assert_eq!(below.regime(), Regime::Johnson);
        assert_eq!(above.regime(), Regime::Euler);
        // Both within a few psi of Fy/2 and of each other.
        assert!(approx(
            below.critical_stress(),
            above.critical_stress(),
            1.0
        ));
        assert!(approx(below.critical_stress(), FY_A36 / 2.0, 1.0));
    }

    /// The first derivative also matches at `Cc` (tangency), which is
    /// what makes the Johnson parabola the standard inelastic fit.
    /// Compare central finite differences of the two unchecked curves
    /// straddling the transition.
    #[test]
    fn slopes_match_at_cc_tangency() {
        let cc = transition_slenderness(E_STEEL, FY_A36).unwrap();
        let h = 1e-3;
        let e_plus = Column::new(E_STEEL, FY_A36, cc + h)
            .unwrap()
            .euler_stress_unchecked();
        let e_minus = Column::new(E_STEEL, FY_A36, cc - h)
            .unwrap()
            .euler_stress_unchecked();
        let j_plus = Column::new(E_STEEL, FY_A36, cc + h)
            .unwrap()
            .johnson_stress_unchecked();
        let j_minus = Column::new(E_STEEL, FY_A36, cc - h)
            .unwrap()
            .johnson_stress_unchecked();
        let euler_slope = (e_plus - e_minus) / (2.0 * h);
        let johnson_slope = (j_plus - j_minus) / (2.0 * h);
        // Slopes equal to within finite-difference truncation.
        assert!(
            approx(euler_slope, johnson_slope, 1e-2),
            "Euler slope {euler_slope} vs Johnson slope {johnson_slope}"
        );
    }

    // --- Fcr <= Fy everywhere ----------------------------------------

    /// The critical stress never exceeds the yield stress across the
    /// whole slenderness range, including the squash load at `λ = 0`.
    #[test]
    fn critical_stress_never_exceeds_yield() {
        for n in 0..=300 {
            let lambda = n as f64; // 0 .. 300
            let col = Column::new(E_STEEL, FY_A36, lambda).unwrap();
            let fcr = col.critical_stress();
            assert!(
                fcr <= FY_A36 + EPS_PSI,
                "Fcr {fcr} exceeds Fy {FY_A36} at lambda {lambda}"
            );
        }
    }

    /// At zero slenderness the column squashes at the yield stress: the
    /// Johnson parabola's intercept is exactly `Fy`.
    #[test]
    fn zero_slenderness_yields_fy() {
        let col = Column::new(E_STEEL, FY_A36, 0.0).unwrap();
        assert_eq!(col.regime(), Regime::Johnson);
        assert!(approx(col.critical_stress(), FY_A36, EPS_PSI));
    }

    // --- Monotonicity: higher slenderness -> lower capacity ----------

    /// Critical stress is strictly decreasing in slenderness over the
    /// physical range — a more slender column always buckles at a lower
    /// stress.
    #[test]
    fn capacity_decreases_with_slenderness() {
        let mut prev = f64::INFINITY;
        for n in 0..=250 {
            let lambda = n as f64 * 0.8; // 0 .. 200 in 0.8 steps
            let col = Column::new(E_STEEL, FY_A36, lambda).unwrap();
            let fcr = col.critical_stress();
            assert!(
                fcr < prev + EPS_PSI,
                "capacity rose at lambda {lambda}: {fcr} >= {prev}"
            );
            prev = fcr;
        }
    }

    /// The allowable stress is likewise strictly decreasing in
    /// slenderness.
    #[test]
    fn allowable_decreases_with_slenderness() {
        let a = Column::new(E_STEEL, FY_A36, 40.0)
            .unwrap()
            .allowable_stress();
        let b = Column::new(E_STEEL, FY_A36, 90.0)
            .unwrap()
            .allowable_stress();
        let c = Column::new(E_STEEL, FY_A36, 160.0)
            .unwrap()
            .allowable_stress();
        assert!(a > b, "Fa(40)={a} should exceed Fa(90)={b}");
        assert!(b > c, "Fa(90)={b} should exceed Fa(160)={c}");
    }

    // --- Allowable stress / factor of safety -------------------------

    /// The allowable stress is the critical stress divided by the factor
    /// of safety, and the factor of safety is always greater than one,
    /// so `Fa < Fcr`.
    #[test]
    fn allowable_is_critical_over_safety_factor() {
        let col = Column::new(E_STEEL, FY_A36, 90.0).unwrap();
        let fs = col.factor_of_safety_aisc();
        assert!(fs > 1.0, "factor of safety {fs} must exceed 1");
        let expected = col.critical_stress() / fs;
        assert!(approx(col.allowable_stress(), expected, EPS_PSI));
        assert!(col.allowable_stress() < col.critical_stress());
    }

    /// In the Euler regime the AISC factor of safety is the constant
    /// `23/12 ≈ 1.9167`.
    #[test]
    fn euler_safety_factor_is_23_over_12() {
        let col = Column::new(E_STEEL, FY_A36, 150.0).unwrap();
        assert_eq!(col.regime(), Regime::Euler);
        assert!(approx(col.factor_of_safety_aisc(), 23.0 / 12.0, EPS_REL));
    }

    /// At `λ = 0` the Johnson-branch factor of safety is exactly
    /// `5/3 ≈ 1.6667`.
    #[test]
    fn johnson_safety_factor_at_zero_is_5_over_3() {
        let col = Column::new(E_STEEL, FY_A36, 0.0).unwrap();
        assert!(approx(col.factor_of_safety_aisc(), 5.0 / 3.0, EPS_REL));
    }

    /// At `λ = Cc` the Johnson-branch factor of safety ramps up to
    /// exactly `23/12`, matching the Euler-branch constant — the safety
    /// factor is itself continuous across the transition.
    #[test]
    fn safety_factor_continuous_at_cc() {
        let cc = transition_slenderness(E_STEEL, FY_A36).unwrap();
        // Just inside the Johnson regime.
        let col = Column::new(E_STEEL, FY_A36, cc - 1e-6).unwrap();
        assert!(approx(col.factor_of_safety_aisc(), 23.0 / 12.0, 1e-6));
    }

    /// The AISC allowable stress at `λ = 0` is `Fy / (5/3) = 0.6 Fy`,
    /// the familiar 60%-of-yield short-column working stress.
    #[test]
    fn allowable_at_zero_is_60_percent_yield() {
        let col = Column::new(E_STEEL, FY_A36, 0.0).unwrap();
        assert!(approx(col.allowable_stress(), 0.6 * FY_A36, EPS_PSI));
    }

    /// The constant-FS allowable matches a hand computation and rejects
    /// a non-positive safety factor.
    #[test]
    fn allowable_with_constant_factor() {
        let col = Column::new(E_STEEL, FY_A36, 150.0).unwrap();
        let fa = col.allowable_stress_with(1.92).unwrap();
        assert!(approx(fa, col.critical_stress() / 1.92, EPS_PSI));
        assert!(col.allowable_stress_with(0.0).is_err());
        assert!(col.allowable_stress_with(-1.0).is_err());
    }

    /// The allowable axial load is `Fa · A` and rejects a non-positive
    /// area.
    #[test]
    fn allowable_load_scales_with_area() {
        let col = Column::new(E_STEEL, FY_A36, 80.0).unwrap();
        let area = 10.0; // in^2
        let p = col.allowable_load(area).unwrap();
        assert!(approx(p, col.allowable_stress() * area, EPS_PSI));
        assert!(col.allowable_load(0.0).is_err());
    }

    // --- Worked end-to-end example -----------------------------------

    /// A W-shape column: `K = 1`, `L = 144 in`, `r = 2.10 in`, A36.
    /// Cross-checks slenderness, regime, critical stress, and allowable
    /// stress against a hand calculation.
    #[test]
    fn worked_example_w_shape() {
        let col = Column::from_geometry(E_STEEL, FY_A36, 1.0, 144.0, 2.10).unwrap();
        // lambda = 144 / 2.10 = 68.5714
        assert!(approx(col.slenderness(), 68.5714, 1e-3));
        // Cc = 126.0989 -> short column, Johnson.
        assert_eq!(col.regime(), Regime::Johnson);
        // ratio = 68.5714/126.0993 = 0.543780
        // Fcr = 36000 * (1 - 0.5 * 0.543780^2) = 30677.28
        assert!(approx(col.critical_stress(), 30677.28, 0.1));
        // FS = 5/3 + 3/8*0.543780 - 1/8*0.543780^3 = 1.850487
        assert!(approx(col.factor_of_safety_aisc(), 1.850487, 1e-4));
        // Fa = 30677.28 / 1.850487 = 16577.95
        assert!(approx(col.allowable_stress(), 16577.95, 0.1));
    }

    // --- Unit independence -------------------------------------------

    /// The model is unit-consistent: running A992 steel in SI (MPa)
    /// gives `Cc` independent of the stress unit because `Cc` depends
    /// only on the dimensionless ratio `E / Fy`.
    #[test]
    fn cc_depends_only_on_modulus_ratio() {
        // US: E = 29e6 psi, Fy = 50e3 psi.
        let cc_us = transition_slenderness(29_000_000.0, 50_000.0).unwrap();
        // SI: E = 200_000 MPa, Fy = 344.74 MPa (same E/Fy ratio).
        let cc_si = transition_slenderness(200_000.0, 200_000.0 / (29e6 / 50e3)).unwrap();
        assert!(approx(cc_us, cc_si, 1e-6), "Cc_us {cc_us} != Cc_si {cc_si}");
    }

    // --- Validation guards -------------------------------------------

    #[test]
    fn constructors_reject_nonphysical_inputs() {
        assert!(Column::new(-1.0, FY_A36, 50.0).is_err());
        assert!(Column::new(E_STEEL, 0.0, 50.0).is_err());
        assert!(Column::new(E_STEEL, FY_A36, -1.0).is_err());
        assert!(Column::from_geometry(E_STEEL, FY_A36, 0.0, 100.0, 2.0).is_err());
        assert!(Column::from_geometry(E_STEEL, FY_A36, 1.0, 100.0, 0.0).is_err());
        assert!(Column::from_geometry(E_STEEL, FY_A36, 1.0, -1.0, 2.0).is_err());
        assert!(slenderness_ratio(1.0, 100.0, 0.0).is_err());
        assert!(transition_slenderness(0.0, FY_A36).is_err());
    }

    #[test]
    fn euler_undefined_at_zero_slenderness() {
        let col = Column::new(E_STEEL, FY_A36, 0.0).unwrap();
        assert!(col.euler_stress().is_err());
        // The unchecked variant returns +infinity there.
        assert!(col.euler_stress_unchecked().is_infinite());
    }

    /// The regime discriminator stringifies to its stable identifiers.
    #[test]
    fn regime_string_identifiers() {
        assert_eq!(Regime::Euler.as_str(), "euler");
        assert_eq!(Regime::Johnson.as_str(), "johnson");
    }
}
