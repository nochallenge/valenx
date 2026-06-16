//! Basal-metabolic-rate predictive equations.
//!
//! Three population-fit regressions. Two are closed-form linear
//! functions of body mass, height, age and biological sex:
//!
//! - **Mifflin-St Jeor (1990)** — the modern default; validates against
//!   indirect calorimetry better than older equations for the general
//!   population.
//! - **Harris-Benedict** in the **Roza & Shizgal (1984)** revision —
//!   the classic equation, retained for comparison and because many
//!   nutrition references still quote it.
//!
//! The third regresses on body composition instead:
//!
//! - **Katch-McArdle** — `370 + 21.6·LBM`, a function of lean body mass
//!   alone, so it carries no sex term and needs an accurate body-fat
//!   measurement rather than height and age.
//!
//! The first two take SI anthropometry (kilograms, centimetres, years);
//! Katch-McArdle takes lean mass in kilograms. All three return basal
//! metabolic rate in **kilocalories per day**.
//!
//! # Worked example
//!
//! For a 30-year-old, 80 kg, 180 cm male, Mifflin-St Jeor gives
//! `10·80 + 6.25·180 − 5·30 + 5 = 1780` kcal/day:
//!
//! ```
//! use valenx_bmr::{Sex, mifflin_st_jeor};
//!
//! let bmr = mifflin_st_jeor(Sex::Male, 80.0, 180.0, 30.0).unwrap();
//! assert!((bmr - 1780.0).abs() < 1e-9);
//! ```

use crate::error::{BmrError, Result};
use serde::{Deserialize, Serialize};

/// Biological sex, as required by the BMR regressions.
///
/// The equations were fit on sex-stratified cohorts and carry a
/// sex-specific intercept (and, for Harris-Benedict, sex-specific
/// slopes), so a value here is mandatory rather than optional. This is
/// the binary biological-sex variable the source studies used; it is
/// not a statement about gender identity.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash, Serialize, Deserialize)]
pub enum Sex {
    /// Male coefficients / intercept.
    Male,
    /// Female coefficients / intercept.
    Female,
}

/// Which predictive equation to evaluate.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash, Serialize, Deserialize)]
pub enum BmrEquation {
    /// Mifflin-St Jeor (1990). The modern default.
    MifflinStJeor,
    /// Harris-Benedict, Roza & Shizgal (1984) revision.
    HarrisBenedict,
}

impl BmrEquation {
    /// Evaluate this equation for the given anthropometry.
    ///
    /// `mass_kg` and `height_cm` must be strictly positive; `age_years`
    /// must be strictly positive. Returns kilocalories per day.
    ///
    /// # Errors
    ///
    /// Propagates [`BmrError::OutOfRange`] / [`BmrError::NotFinite`]
    /// from the underlying equation if any input is non-positive or
    /// non-finite.
    pub fn evaluate(self, sex: Sex, mass_kg: f64, height_cm: f64, age_years: f64) -> Result<f64> {
        match self {
            BmrEquation::MifflinStJeor => mifflin_st_jeor(sex, mass_kg, height_cm, age_years),
            BmrEquation::HarrisBenedict => harris_benedict(sex, mass_kg, height_cm, age_years),
        }
    }
}

/// Validate the three anthropometric inputs shared by every equation.
///
/// Returns the validated `(mass, height, age)` triple unchanged so the
/// caller can shadow its bindings.
fn validate_anthropometry(mass_kg: f64, height_cm: f64, age_years: f64) -> Result<(f64, f64, f64)> {
    let mass_kg = BmrError::require_positive("mass_kg", mass_kg)?;
    let height_cm = BmrError::require_positive("height_cm", height_cm)?;
    let age_years = BmrError::require_positive("age_years", age_years)?;
    Ok((mass_kg, height_cm, age_years))
}

/// Mifflin-St Jeor basal metabolic rate, in kcal/day.
///
/// ```text
/// BMR = 10·mass(kg) + 6.25·height(cm) − 5·age(yr) + s
/// ```
///
/// where the sex term `s` is `+5` for [`Sex::Male`] and `−161` for
/// [`Sex::Female`]. Same slopes for both sexes; only the intercept
/// differs.
///
/// # Errors
///
/// [`BmrError::OutOfRange`] / [`BmrError::NotFinite`] if `mass_kg`,
/// `height_cm` or `age_years` is non-positive or non-finite.
pub fn mifflin_st_jeor(sex: Sex, mass_kg: f64, height_cm: f64, age_years: f64) -> Result<f64> {
    let (mass_kg, height_cm, age_years) = validate_anthropometry(mass_kg, height_cm, age_years)?;
    let sex_term = match sex {
        Sex::Male => 5.0,
        Sex::Female => -161.0,
    };
    Ok(10.0 * mass_kg + 6.25 * height_cm - 5.0 * age_years + sex_term)
}

/// Harris-Benedict basal metabolic rate (Roza & Shizgal 1984
/// revision), in kcal/day.
///
/// Unlike Mifflin-St Jeor, this equation uses sex-specific intercepts
/// **and** sex-specific slopes:
///
/// ```text
/// Male:   BMR = 88.362 + 13.397·mass + 4.799·height − 5.677·age
/// Female: BMR = 447.593 +  9.247·mass + 3.098·height − 4.330·age
/// ```
///
/// (mass in kg, height in cm, age in years).
///
/// # Errors
///
/// [`BmrError::OutOfRange`] / [`BmrError::NotFinite`] if `mass_kg`,
/// `height_cm` or `age_years` is non-positive or non-finite.
pub fn harris_benedict(sex: Sex, mass_kg: f64, height_cm: f64, age_years: f64) -> Result<f64> {
    let (mass_kg, height_cm, age_years) = validate_anthropometry(mass_kg, height_cm, age_years)?;
    let bmr = match sex {
        Sex::Male => 88.362 + 13.397 * mass_kg + 4.799 * height_cm - 5.677 * age_years,
        Sex::Female => 447.593 + 9.247 * mass_kg + 3.098 * height_cm - 4.330 * age_years,
    };
    Ok(bmr)
}

/// Validate a body-fat fraction: finite and in the half-open `[0, 1)`.
///
/// A fraction of `1.0` (or more) would leave zero or negative lean
/// mass, so the upper bound is strict; `0.0` (all-lean) is accepted as
/// the limiting case.
///
/// # Errors
///
/// [`BmrError::NotFinite`] if `value` is `NaN` or infinite;
/// [`BmrError::OutOfRange`] if `value` is negative or `>= 1.0`.
fn validate_body_fat_fraction(value: f64) -> Result<f64> {
    if !value.is_finite() {
        return Err(BmrError::NotFinite {
            name: "body_fat_fraction",
            value,
        });
    }
    if !(0.0..1.0).contains(&value) {
        return Err(BmrError::OutOfRange {
            name: "body_fat_fraction",
            value,
            reason: "must be in [0, 1) so lean mass stays positive",
        });
    }
    Ok(value)
}

/// Lean body mass (kg) from total body mass and a body-fat fraction.
///
/// ```text
/// LBM = mass · (1 − body_fat_fraction)
/// ```
///
/// `mass_kg` must be finite and strictly positive; `body_fat_fraction`
/// must be finite and in the half-open `[0, 1)`. The result is the
/// fat-free mass that drives the [`katch_mcardle`] equation.
///
/// ```
/// use valenx_bmr::lean_body_mass;
///
/// // 80 kg at 20% body fat -> 64 kg of lean tissue.
/// let lbm = lean_body_mass(80.0, 0.20).unwrap();
/// assert!((lbm - 64.0).abs() < 1e-9);
/// ```
///
/// # Errors
///
/// [`BmrError::OutOfRange`] / [`BmrError::NotFinite`] if `mass_kg` is
/// non-positive or non-finite, or `body_fat_fraction` is outside
/// `[0, 1)` or non-finite.
pub fn lean_body_mass(mass_kg: f64, body_fat_fraction: f64) -> Result<f64> {
    let mass_kg = BmrError::require_positive("mass_kg", mass_kg)?;
    let body_fat_fraction = validate_body_fat_fraction(body_fat_fraction)?;
    Ok(mass_kg * (1.0 - body_fat_fraction))
}

/// Katch-McArdle basal metabolic rate, in kcal/day.
///
/// Unlike [`mifflin_st_jeor`] and [`harris_benedict`], which regress on
/// mass, height, age and sex, this equation depends only on **lean body
/// mass** — the metabolically active tissue — and so carries no sex
/// term:
///
/// ```text
/// BMR = 370 + 21.6 · lean_body_mass(kg)
/// ```
///
/// When an accurate body-composition measurement is available this is
/// often the most individual-specific of the three, because it removes
/// the confounding effect of fat mass (which is far less metabolically
/// active). It is deliberately **not** a [`BmrEquation`] variant: that
/// enum is parameterised by anthropometry, whereas Katch-McArdle is
/// parameterised by body composition.
///
/// ```
/// use valenx_bmr::katch_mcardle;
///
/// // 60 kg of lean mass: 370 + 21.6·60 = 1666 kcal/day.
/// let bmr = katch_mcardle(60.0).unwrap();
/// assert!((bmr - 1666.0).abs() < 1e-9);
/// ```
///
/// # Errors
///
/// [`BmrError::OutOfRange`] / [`BmrError::NotFinite`] if
/// `lean_body_mass_kg` is non-positive or non-finite.
pub fn katch_mcardle(lean_body_mass_kg: f64) -> Result<f64> {
    let lbm = BmrError::require_positive("lean_body_mass_kg", lean_body_mass_kg)?;
    Ok(370.0 + 21.6 * lbm)
}

/// Katch-McArdle BMR from total mass and a body-fat fraction, in
/// kcal/day.
///
/// Convenience wrapper that computes [`lean_body_mass`] from `mass_kg`
/// and `body_fat_fraction` and feeds it to [`katch_mcardle`] in a single
/// call.
///
/// # Errors
///
/// Propagates the validation errors of [`lean_body_mass`] — a
/// non-positive / non-finite `mass_kg`, or a `body_fat_fraction`
/// outside `[0, 1)` or non-finite.
pub fn katch_mcardle_from_body_fat(mass_kg: f64, body_fat_fraction: f64) -> Result<f64> {
    let lbm = lean_body_mass(mass_kg, body_fat_fraction)?;
    katch_mcardle(lbm)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Tolerance for exact closed-form checks (the formulas are exact
    /// real arithmetic; only float rounding separates us from truth).
    const EPS: f64 = 1e-9;

    #[test]
    fn mifflin_known_male_example() {
        // Canonical worked example: 30 y, 80 kg, 180 cm male.
        // 10*80 + 6.25*180 - 5*30 + 5 = 800 + 1125 - 150 + 5 = 1780.
        let bmr = mifflin_st_jeor(Sex::Male, 80.0, 180.0, 30.0).unwrap();
        assert!((bmr - 1780.0).abs() < EPS, "got {bmr}");
        // Well within a loose ±5 kcal tolerance too.
        assert!((bmr - 1780.0).abs() < 5.0, "got {bmr}");
    }

    #[test]
    fn mifflin_known_female_example() {
        // 30 y, 80 kg, 180 cm female:
        // 10*80 + 6.25*180 - 5*30 - 161 = 2075 - 161 = 1614.
        let bmr = mifflin_st_jeor(Sex::Female, 80.0, 180.0, 30.0).unwrap();
        assert!((bmr - 1614.0).abs() < EPS, "got {bmr}");
    }

    #[test]
    fn female_bmr_below_male_for_same_body() {
        // Same anthropometry: the female intercept is 166 kcal lower
        // (5 − (−161)), so female BMR must be strictly below male.
        let male = mifflin_st_jeor(Sex::Male, 70.0, 175.0, 40.0).unwrap();
        let female = mifflin_st_jeor(Sex::Female, 70.0, 175.0, 40.0).unwrap();
        assert!(female < male, "female {female} should be < male {male}");
        assert!((male - female - 166.0).abs() < EPS, "delta wrong");
    }

    #[test]
    fn harris_known_male_value() {
        // 88.362 + 13.397*70 + 4.799*175 - 5.677*40
        //   = 88.362 + 937.79 + 839.825 - 227.08 = 1638.897.
        let bmr = harris_benedict(Sex::Male, 70.0, 175.0, 40.0).unwrap();
        let expected = 88.362 + 13.397 * 70.0 + 4.799 * 175.0 - 5.677 * 40.0;
        assert!((bmr - expected).abs() < EPS, "got {bmr}");
        assert!((bmr - 1638.897).abs() < 1e-3, "got {bmr}");
    }

    #[test]
    fn harris_known_female_value() {
        let bmr = harris_benedict(Sex::Female, 60.0, 165.0, 35.0).unwrap();
        let expected = 447.593 + 9.247 * 60.0 + 3.098 * 165.0 - 4.330 * 35.0;
        assert!((bmr - expected).abs() < EPS, "got {bmr}");
    }

    #[test]
    fn harris_female_below_male_same_body() {
        let male = harris_benedict(Sex::Male, 75.0, 178.0, 28.0).unwrap();
        let female = harris_benedict(Sex::Female, 75.0, 178.0, 28.0).unwrap();
        assert!(female < male, "female {female} should be < male {male}");
    }

    #[test]
    fn two_equations_are_in_the_same_ballpark() {
        // They are different fits but should agree to within ~10% for a
        // typical adult — a sanity check that neither is grossly wrong.
        let m = mifflin_st_jeor(Sex::Male, 80.0, 180.0, 30.0).unwrap();
        let h = harris_benedict(Sex::Male, 80.0, 180.0, 30.0).unwrap();
        let rel = (m - h).abs() / m;
        assert!(rel < 0.10, "Mifflin {m} vs Harris {h} differ by {rel}");
    }

    #[test]
    fn equation_enum_dispatch_matches_free_fn() {
        let via_enum = BmrEquation::MifflinStJeor
            .evaluate(Sex::Male, 80.0, 180.0, 30.0)
            .unwrap();
        let via_fn = mifflin_st_jeor(Sex::Male, 80.0, 180.0, 30.0).unwrap();
        assert!((via_enum - via_fn).abs() < EPS);

        let via_enum = BmrEquation::HarrisBenedict
            .evaluate(Sex::Female, 60.0, 165.0, 35.0)
            .unwrap();
        let via_fn = harris_benedict(Sex::Female, 60.0, 165.0, 35.0).unwrap();
        assert!((via_enum - via_fn).abs() < EPS);
    }

    #[test]
    fn bmr_increases_with_mass_and_decreases_with_age() {
        let base = mifflin_st_jeor(Sex::Male, 70.0, 175.0, 30.0).unwrap();
        let heavier = mifflin_st_jeor(Sex::Male, 80.0, 175.0, 30.0).unwrap();
        let older = mifflin_st_jeor(Sex::Male, 70.0, 175.0, 50.0).unwrap();
        assert!(heavier > base, "more mass should raise BMR");
        assert!(older < base, "older should lower BMR");
    }

    #[test]
    fn rejects_nonpositive_and_nonfinite_inputs() {
        assert!(mifflin_st_jeor(Sex::Male, 0.0, 180.0, 30.0).is_err());
        assert!(mifflin_st_jeor(Sex::Male, 80.0, -1.0, 30.0).is_err());
        assert!(mifflin_st_jeor(Sex::Male, 80.0, 180.0, 0.0).is_err());
        assert!(harris_benedict(Sex::Female, f64::NAN, 165.0, 35.0).is_err());
        assert!(harris_benedict(Sex::Female, 60.0, f64::INFINITY, 35.0).is_err());
    }

    #[test]
    fn katch_mcardle_known_value() {
        // Canonical closed form: 60 kg lean -> 370 + 21.6*60
        //   = 370 + 1296 = 1666 kcal/day.
        let bmr = katch_mcardle(60.0).unwrap();
        assert!((bmr - 1666.0).abs() < EPS, "got {bmr}");
    }

    #[test]
    fn lean_body_mass_from_fraction() {
        // 80 kg at 20% body fat -> 64 kg lean.
        let lbm = lean_body_mass(80.0, 0.20).unwrap();
        assert!((lbm - 64.0).abs() < EPS, "got {lbm}");
        // Zero body fat (limiting case): the whole mass is lean.
        let all_lean = lean_body_mass(80.0, 0.0).unwrap();
        assert!((all_lean - 80.0).abs() < EPS, "got {all_lean}");
    }

    #[test]
    fn katch_from_body_fat_matches_two_step() {
        let one = katch_mcardle_from_body_fat(80.0, 0.20).unwrap();
        let two = katch_mcardle(lean_body_mass(80.0, 0.20).unwrap()).unwrap();
        assert!((one - two).abs() < EPS, "one {one} vs two {two}");
        // 64 kg lean -> 370 + 21.6*64 = 1752.4.
        assert!((one - 1752.4).abs() < 1e-6, "got {one}");
    }

    #[test]
    fn katch_agrees_with_mifflin_for_typical_male() {
        // Independent cross-check: for an 80 kg, 180 cm, 30 y male at a
        // typical ~18% body fat, the lean-mass equation and the
        // anthropometric one should land within ~10% of each other.
        let mifflin = mifflin_st_jeor(Sex::Male, 80.0, 180.0, 30.0).unwrap();
        let katch = katch_mcardle_from_body_fat(80.0, 0.18).unwrap();
        let rel = (mifflin - katch).abs() / mifflin;
        assert!(
            rel < 0.10,
            "Mifflin {mifflin} vs Katch {katch} differ by {rel}"
        );
    }

    #[test]
    fn katch_mcardle_increases_with_lean_mass() {
        let lean = katch_mcardle(55.0).unwrap();
        let more_muscle = katch_mcardle(70.0).unwrap();
        assert!(more_muscle > lean, "more lean mass should raise BMR");
    }

    #[test]
    fn katch_mcardle_rejects_bad_inputs() {
        assert!(katch_mcardle(0.0).is_err());
        assert!(katch_mcardle(-5.0).is_err());
        assert!(katch_mcardle(f64::NAN).is_err());
        // body_fat_fraction must be in the half-open [0, 1).
        assert!(lean_body_mass(80.0, 1.0).is_err());
        assert!(lean_body_mass(80.0, -0.1).is_err());
        assert!(lean_body_mass(80.0, f64::INFINITY).is_err());
        assert!(lean_body_mass(0.0, 0.2).is_err());
    }

    #[test]
    fn sex_serde_roundtrip() {
        let json = serde_json::to_string(&Sex::Female).unwrap();
        let back: Sex = serde_json::from_str(&json).unwrap();
        assert_eq!(back, Sex::Female);
    }
}
