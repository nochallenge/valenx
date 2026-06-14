//! Energy balance: total daily energy expenditure (TDEE) and a
//! first-order body-mass projection.
//!
//! Once a basal metabolic rate is in hand (see [`crate::bmr`]), the
//! everyday quantities are:
//!
//! - **TDEE** — basal expenditure scaled by a *physical-activity
//!   level* (PAL) multiplier: `TDEE = BMR · activity_factor`. Activity
//!   can only add to basal metabolism, so the multiplier is `>= 1.0`.
//! - **Energy balance** — the daily difference `intake − TDEE`. A
//!   surplus stores energy (mass gain); a deficit draws on stores
//!   (mass loss).
//! - **Mass projection** — converting a *sustained* daily balance over
//!   some number of days into a body-mass change using the
//!   [`KCAL_PER_KG`] energy-density rule of thumb.
//!
//! All energies are kilocalories; masses are kilograms.
//!
//! # Honest scope
//!
//! The mass projection is the classic **linear** energy-balance model
//! (the "7700 kcal per kilogram" / Wishnofsky rule). It deliberately
//! ignores adaptive thermogenesis (metabolism falling as you eat less),
//! the shifting fat-vs-lean composition of the tissue gained or lost,
//! and water-weight swings. Real long-run weight change is sub-linear;
//! treat the projection as a textbook first approximation, never a
//! clinical prescription.

use crate::error::{BmrError, Result};
use serde::{Deserialize, Serialize};

/// Energy density used to convert an energy imbalance into a body-mass
/// change, in kilocalories per kilogram.
///
/// `7700 kcal/kg` is the long-standing rule-of-thumb figure (the
/// Wishnofsky "3500 kcal per pound" rule, ≈ 7716 kcal/kg, rounded). It
/// approximates the energy stored in a kilogram of mixed adipose
/// tissue.
pub const KCAL_PER_KG: f64 = 7700.0;

/// Standard physical-activity levels and their TDEE multipliers.
///
/// These are the widely tabulated factors applied to a Mifflin-St Jeor
/// or Harris-Benedict BMR. Each names a lifestyle band; use
/// [`ActivityLevel::factor`] to get the multiplier or
/// [`ActivityLevel::from_factor`] to bucket an arbitrary value.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash, Serialize, Deserialize)]
pub enum ActivityLevel {
    /// Little or no exercise; desk job. Factor `1.2`.
    Sedentary,
    /// Light exercise 1–3 days per week. Factor `1.375`.
    LightlyActive,
    /// Moderate exercise 3–5 days per week. Factor `1.55`.
    ModeratelyActive,
    /// Hard exercise 6–7 days per week. Factor `1.725`.
    VeryActive,
    /// Very hard exercise plus a physical job or two-a-day training.
    /// Factor `1.9`.
    ExtraActive,
}

impl ActivityLevel {
    /// The TDEE multiplier for this activity level.
    pub fn factor(self) -> f64 {
        match self {
            ActivityLevel::Sedentary => 1.2,
            ActivityLevel::LightlyActive => 1.375,
            ActivityLevel::ModeratelyActive => 1.55,
            ActivityLevel::VeryActive => 1.725,
            ActivityLevel::ExtraActive => 1.9,
        }
    }

    /// Bucket an arbitrary multiplier into the nearest named level,
    /// choosing the band whose canonical factor is closest.
    ///
    /// Useful for labelling a user-entered factor. Ties round toward
    /// the lower (less active) band, which the iteration order makes
    /// natural.
    pub fn from_factor(factor: f64) -> ActivityLevel {
        const LEVELS: [ActivityLevel; 5] = [
            ActivityLevel::Sedentary,
            ActivityLevel::LightlyActive,
            ActivityLevel::ModeratelyActive,
            ActivityLevel::VeryActive,
            ActivityLevel::ExtraActive,
        ];
        let mut best = ActivityLevel::Sedentary;
        let mut best_dist = f64::INFINITY;
        for level in LEVELS {
            let dist = (level.factor() - factor).abs();
            if dist < best_dist {
                best_dist = dist;
                best = level;
            }
        }
        best
    }
}

/// Total daily energy expenditure from a BMR and an explicit activity
/// multiplier.
///
/// `tdee = bmr_kcal · activity_factor`. Because the factor is `>= 1.0`,
/// the result is always at least the BMR (and strictly greater whenever
/// `activity_factor > 1.0`).
///
/// # Errors
///
/// [`BmrError::OutOfRange`] / [`BmrError::NotFinite`] if `bmr_kcal` is
/// non-positive / non-finite, or if `activity_factor` is `< 1.0` /
/// non-finite.
pub fn tdee(bmr_kcal: f64, activity_factor: f64) -> Result<f64> {
    let bmr_kcal = BmrError::require_positive("bmr_kcal", bmr_kcal)?;
    let factor = BmrError::require_activity_factor(activity_factor)?;
    Ok(bmr_kcal * factor)
}

/// Total daily energy expenditure using a named [`ActivityLevel`].
///
/// Convenience wrapper over [`tdee`] that pulls the multiplier from the
/// level. The level's factors are all `>= 1.2`, so this never fails on
/// the activity term.
///
/// # Errors
///
/// [`BmrError::OutOfRange`] / [`BmrError::NotFinite`] if `bmr_kcal` is
/// non-positive or non-finite.
pub fn tdee_for_level(bmr_kcal: f64, level: ActivityLevel) -> Result<f64> {
    tdee(bmr_kcal, level.factor())
}

/// Daily energy balance: `intake_kcal − tdee_kcal`.
///
/// Positive is a surplus (energy stored), negative is a deficit (energy
/// drawn from stores). Both inputs must be finite; intake may not be
/// negative, but the balance itself can be.
///
/// # Errors
///
/// [`BmrError::OutOfRange`] if `intake_kcal` is negative;
/// [`BmrError::NotFinite`] if either argument is non-finite.
pub fn daily_energy_balance(intake_kcal: f64, tdee_kcal: f64) -> Result<f64> {
    if intake_kcal < 0.0 {
        return Err(BmrError::OutOfRange {
            name: "intake_kcal",
            value: intake_kcal,
            reason: "intake cannot be negative",
        });
    }
    let intake = BmrError::require_finite("intake_kcal", intake_kcal)?;
    let expenditure = BmrError::require_finite("tdee_kcal", tdee_kcal)?;
    Ok(intake - expenditure)
}

/// Project the body-mass change from a *sustained* daily energy balance
/// held over `days` days.
///
/// ```text
/// Δmass(kg) = balance_kcal_per_day · days / KCAL_PER_KG
/// ```
///
/// A positive `daily_balance_kcal` (surplus) yields mass **gain**
/// (positive result); a negative balance (deficit) yields mass **loss**
/// (negative result). The balance may be any finite number; `days` must
/// be non-negative and finite.
///
/// # Errors
///
/// [`BmrError::OutOfRange`] if `days` is negative;
/// [`BmrError::NotFinite`] if either argument is non-finite.
pub fn mass_change_kg(daily_balance_kcal: f64, days: f64) -> Result<f64> {
    let balance = BmrError::require_finite("daily_balance_kcal", daily_balance_kcal)?;
    if days < 0.0 {
        return Err(BmrError::OutOfRange {
            name: "days",
            value: days,
            reason: "duration cannot be negative",
        });
    }
    let days = BmrError::require_finite("days", days)?;
    Ok(balance * days / KCAL_PER_KG)
}

/// A complete, validated energy-balance snapshot for one person.
///
/// Produced by [`EnergyBalance::new`]; bundles the BMR, the chosen
/// activity factor, the derived TDEE, intake and the resulting daily
/// balance so a caller (or a UI) gets every figure consistently from a
/// single validated source.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct EnergyBalance {
    /// Basal metabolic rate, kcal/day.
    pub bmr_kcal: f64,
    /// Physical-activity multiplier applied to the BMR (`>= 1.0`).
    pub activity_factor: f64,
    /// Total daily energy expenditure, kcal/day (`bmr · factor`).
    pub tdee_kcal: f64,
    /// Daily energy intake, kcal/day.
    pub intake_kcal: f64,
    /// Daily energy balance, kcal/day (`intake − tdee`).
    pub daily_balance_kcal: f64,
}

impl EnergyBalance {
    /// Build a snapshot from a BMR, an activity factor and a daily
    /// intake, validating every input and deriving TDEE and balance.
    ///
    /// # Errors
    ///
    /// Propagates validation errors from [`tdee`] and
    /// [`daily_energy_balance`].
    pub fn new(bmr_kcal: f64, activity_factor: f64, intake_kcal: f64) -> Result<Self> {
        let tdee_kcal = tdee(bmr_kcal, activity_factor)?;
        let daily_balance_kcal = daily_energy_balance(intake_kcal, tdee_kcal)?;
        Ok(EnergyBalance {
            bmr_kcal,
            activity_factor,
            tdee_kcal,
            intake_kcal,
            daily_balance_kcal,
        })
    }

    /// `true` if intake exceeds expenditure (mass-gain trajectory).
    pub fn is_surplus(&self) -> bool {
        self.daily_balance_kcal > 0.0
    }

    /// `true` if expenditure exceeds intake (mass-loss trajectory).
    pub fn is_deficit(&self) -> bool {
        self.daily_balance_kcal < 0.0
    }

    /// Projected body-mass change (kg) if this daily balance is held
    /// for `days` days. See [`mass_change_kg`].
    ///
    /// # Errors
    ///
    /// [`BmrError::OutOfRange`] if `days` is negative.
    pub fn projected_mass_change_kg(&self, days: f64) -> Result<f64> {
        mass_change_kg(self.daily_balance_kcal, days)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    #[test]
    fn tdee_exceeds_bmr_when_active() {
        // Any activity factor > 1 must lift TDEE above BMR.
        let bmr = 1780.0;
        let out = tdee(bmr, 1.55).unwrap();
        assert!(out > bmr, "TDEE {out} should exceed BMR {bmr}");
        assert!((out - 1780.0 * 1.55).abs() < EPS, "got {out}");
    }

    #[test]
    fn tdee_equals_bmr_at_factor_one() {
        let bmr = 1500.0;
        let out = tdee(bmr, 1.0).unwrap();
        assert!((out - bmr).abs() < EPS, "got {out}");
    }

    #[test]
    fn tdee_rejects_factor_below_one() {
        assert!(tdee(1500.0, 0.9).is_err());
        assert!(tdee(1500.0, f64::NAN).is_err());
        assert!(tdee(-1.0, 1.2).is_err());
    }

    #[test]
    fn activity_level_factors_are_ordered() {
        // The canonical ladder must be strictly increasing.
        let levels = [
            ActivityLevel::Sedentary,
            ActivityLevel::LightlyActive,
            ActivityLevel::ModeratelyActive,
            ActivityLevel::VeryActive,
            ActivityLevel::ExtraActive,
        ];
        for pair in levels.windows(2) {
            assert!(
                pair[0].factor() < pair[1].factor(),
                "{:?} factor {} not < {:?} factor {}",
                pair[0],
                pair[0].factor(),
                pair[1],
                pair[1].factor(),
            );
        }
        // Spot-check the canonical values.
        assert!((ActivityLevel::Sedentary.factor() - 1.2).abs() < EPS);
        assert!((ActivityLevel::ModeratelyActive.factor() - 1.55).abs() < EPS);
        assert!((ActivityLevel::ExtraActive.factor() - 1.9).abs() < EPS);
    }

    #[test]
    fn tdee_for_level_matches_explicit_factor() {
        let a = tdee_for_level(1780.0, ActivityLevel::ModeratelyActive).unwrap();
        let b = tdee(1780.0, 1.55).unwrap();
        assert!((a - b).abs() < EPS, "got {a} vs {b}");
    }

    #[test]
    fn from_factor_buckets_to_nearest() {
        assert_eq!(ActivityLevel::from_factor(1.21), ActivityLevel::Sedentary);
        assert_eq!(
            ActivityLevel::from_factor(1.54),
            ActivityLevel::ModeratelyActive
        );
        assert_eq!(ActivityLevel::from_factor(2.0), ActivityLevel::ExtraActive);
        // Round-trips: each canonical factor maps back to its own level.
        for level in [
            ActivityLevel::Sedentary,
            ActivityLevel::LightlyActive,
            ActivityLevel::ModeratelyActive,
            ActivityLevel::VeryActive,
            ActivityLevel::ExtraActive,
        ] {
            assert_eq!(ActivityLevel::from_factor(level.factor()), level);
        }
    }

    #[test]
    fn surplus_gives_gain_deficit_gives_loss() {
        // +500 kcal/day for 30 days -> +500*30/7700 ≈ +1.948 kg.
        let gain = mass_change_kg(500.0, 30.0).unwrap();
        let expected_gain = 500.0 * 30.0 / KCAL_PER_KG;
        assert!(gain > 0.0, "surplus should gain, got {gain}");
        assert!((gain - expected_gain).abs() < EPS, "got {gain}");
        assert!((gain - 1.9480519).abs() < 1e-6, "got {gain}");

        // −500 kcal/day for 30 days -> symmetric loss.
        let loss = mass_change_kg(-500.0, 30.0).unwrap();
        assert!(loss < 0.0, "deficit should lose, got {loss}");
        assert!((loss + expected_gain).abs() < EPS, "got {loss}");
    }

    #[test]
    fn one_kg_needs_about_7700_kcal() {
        // A net 7700 kcal balance over a single "day" is exactly 1 kg.
        let one = mass_change_kg(KCAL_PER_KG, 1.0).unwrap();
        assert!((one - 1.0).abs() < EPS, "got {one}");
    }

    #[test]
    fn zero_balance_or_zero_days_is_no_change() {
        assert!((mass_change_kg(0.0, 30.0).unwrap()).abs() < EPS);
        assert!((mass_change_kg(500.0, 0.0).unwrap()).abs() < EPS);
    }

    #[test]
    fn mass_change_rejects_negative_days_and_nonfinite() {
        assert!(mass_change_kg(500.0, -1.0).is_err());
        assert!(mass_change_kg(f64::NAN, 30.0).is_err());
        assert!(mass_change_kg(500.0, f64::INFINITY).is_err());
    }

    #[test]
    fn daily_balance_sign_and_validation() {
        let surplus = daily_energy_balance(2500.0, 2000.0).unwrap();
        assert!((surplus - 500.0).abs() < EPS, "got {surplus}");
        let deficit = daily_energy_balance(1800.0, 2200.0).unwrap();
        assert!((deficit + 400.0).abs() < EPS, "got {deficit}");
        // Negative intake is rejected; non-finite is rejected.
        assert!(daily_energy_balance(-1.0, 2000.0).is_err());
        assert!(daily_energy_balance(2000.0, f64::NAN).is_err());
    }

    #[test]
    fn energy_balance_snapshot_is_consistent() {
        let eb = EnergyBalance::new(1780.0, 1.55, 3000.0).unwrap();
        assert!((eb.tdee_kcal - 1780.0 * 1.55).abs() < EPS);
        assert!((eb.daily_balance_kcal - (3000.0 - eb.tdee_kcal)).abs() < EPS);
        assert!(eb.is_surplus());
        assert!(!eb.is_deficit());

        // 30 days of this surplus -> matching projected gain.
        let dm = eb.projected_mass_change_kg(30.0).unwrap();
        let expect = eb.daily_balance_kcal * 30.0 / KCAL_PER_KG;
        assert!((dm - expect).abs() < EPS, "got {dm}");
        assert!(dm > 0.0, "surplus should project a gain, got {dm}");
    }

    #[test]
    fn energy_balance_deficit_branch() {
        let eb = EnergyBalance::new(1780.0, 1.2, 1500.0).unwrap();
        // TDEE = 2136 > 1500 intake -> deficit.
        assert!(eb.is_deficit());
        assert!(!eb.is_surplus());
        let dm = eb.projected_mass_change_kg(60.0).unwrap();
        assert!(dm < 0.0, "deficit should project a loss, got {dm}");
    }

    #[test]
    fn energy_balance_serde_roundtrip() {
        let eb = EnergyBalance::new(1780.0, 1.55, 3000.0).unwrap();
        let json = serde_json::to_string(&eb).unwrap();
        let back: EnergyBalance = serde_json::from_str(&json).unwrap();
        assert_eq!(eb, back);

        let json = serde_json::to_string(&ActivityLevel::VeryActive).unwrap();
        let back: ActivityLevel = serde_json::from_str(&json).unwrap();
        assert_eq!(back, ActivityLevel::VeryActive);
    }
}
