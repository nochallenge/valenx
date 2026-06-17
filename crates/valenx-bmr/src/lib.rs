//! # valenx-bmr
//!
//! Basal metabolic rate and everyday energy-balance arithmetic: predict
//! resting calorie burn from body size and age, scale it to a daily
//! expenditure by activity level, and project the body-mass change a
//! sustained calorie surplus or deficit would produce.
//!
//! ## What
//!
//! Given a person's biological [`Sex`], body mass (kg), height (cm) and
//! age (years), this crate evaluates two standard **basal metabolic
//! rate** (BMR) regressions and turns the result into the figures a
//! calorie planner cares about:
//!
//! - [`mifflin_st_jeor`] / [`harris_benedict`] — BMR in kcal/day, the
//!   resting energy a body burns at complete rest. Dispatchable through
//!   [`BmrEquation`].
//! - [`katch_mcardle`] — a third BMR regression that depends only on
//!   **lean body mass** (no sex term); pair it with [`lean_body_mass`]
//!   (or the [`katch_mcardle_from_body_fat`] shortcut) to derive that
//!   lean mass from total mass and a body-fat fraction.
//! - [`tdee`] / [`tdee_for_level`] — **total daily energy
//!   expenditure**, the BMR scaled by a physical-activity multiplier
//!   ([`ActivityLevel`]) from sedentary (`1.2`) to extra-active
//!   (`1.9`).
//! - [`daily_energy_balance`] and [`mass_change_kg`] — the daily
//!   `intake − TDEE` balance and the body-mass change a sustained
//!   balance implies, bundled together by [`EnergyBalance`]; the inverse
//!   [`days_to_mass_change`] answers the "time to goal" question (and
//!   reports `+infinity` for an unreachable target).
//!
//! ```
//! use valenx_bmr::{ActivityLevel, EnergyBalance, Sex, mifflin_st_jeor, tdee_for_level};
//!
//! // 30-year-old, 80 kg, 180 cm male.
//! let bmr = mifflin_st_jeor(Sex::Male, 80.0, 180.0, 30.0).unwrap();
//! assert!((bmr - 1780.0).abs() < 1e-9); // textbook value
//!
//! let tdee = tdee_for_level(bmr, ActivityLevel::ModeratelyActive).unwrap();
//! assert!(tdee > bmr); // activity always adds to basal burn
//!
//! // Eating 3000 kcal/day against that expenditure -> a surplus.
//! let plan = EnergyBalance::new(bmr, ActivityLevel::ModeratelyActive.factor(), 3000.0).unwrap();
//! assert!(plan.is_surplus());
//! assert!(plan.projected_mass_change_kg(30.0).unwrap() > 0.0); // gains mass
//! ```
//!
//! ## Model
//!
//! - **BMR** uses two population-fit linear regressions of mass,
//!   height, age and sex:
//!   - **Mifflin-St Jeor (1990)** — the modern default
//!     (`10·mass + 6.25·height − 5·age + sex_term`, with the sex term
//!     `+5` male / `−161` female).
//!   - **Harris-Benedict**, in the **Roza & Shizgal (1984)** revision —
//!     the classic equation with sex-specific intercepts and slopes.
//!   - **Katch-McArdle** — `370 + 21.6·LBM`, regressing on lean body
//!     mass alone, so it needs an accurate body-fat measurement instead
//!     of height/age/sex and carries no sex term.
//! - **TDEE** multiplies BMR by a physical-activity-level factor on the
//!   standard sedentary → extra-active ladder. The factor is `>= 1.0`,
//!   so TDEE is always at least the BMR and strictly greater whenever
//!   the factor exceeds `1.0`.
//! - **Mass change** uses the linear energy-balance ("7700 kcal per
//!   kilogram", [`KCAL_PER_KG`]) rule:
//!   `Δmass = (intake − TDEE)·days / 7700`. A surplus projects a gain;
//!   a deficit projects a loss.
//!
//! Inputs are validated through [`BmrError`]'s constructors: masses,
//! heights and ages must be finite and strictly positive, activity
//! factors finite and `>= 1.0`, intake non-negative.
//!
//! ## Honest scope
//!
//! Research/educational grade. These are **textbook closed-form
//! predictive equations** and a **first-order linear energy-balance
//! approximation** — exactly the formulas a physiology or nutrition
//! course teaches — not a clinical, medical, or production dietetics
//! tool. The BMR equations are regressions fit to population cohorts
//! and carry a standard error of roughly ±10% for any individual
//! (Katch-McArdle is moreover only as accurate as the body-fat
//! measurement it is fed); the
//! mass projection ignores adaptive thermogenesis, the fat-vs-lean
//! composition of tissue change, and water-weight swings, so real
//! long-run weight change is sub-linear. Use the numbers to build
//! intuition and ballpark a plan, never to prescribe a diet or make a
//! medical decision.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod bmr;
pub mod energy;
pub mod error;

pub use bmr::{
    harris_benedict, katch_mcardle, katch_mcardle_from_body_fat, lean_body_mass, mifflin_st_jeor,
    BmrEquation, Sex,
};
pub use energy::{
    daily_energy_balance, days_to_mass_change, mass_change_kg, tdee, tdee_for_level, ActivityLevel,
    EnergyBalance, KCAL_PER_KG,
};
pub use error::{BmrError, ErrorCategory, Result};
