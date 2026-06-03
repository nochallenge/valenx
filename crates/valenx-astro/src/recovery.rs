//! Booster-recovery `Δv` budget — the propulsive cost of flying a
//! first stage back to a soft landing.
//!
//! A propulsively-recovered booster spends three retro-propulsive
//! impulses after stage separation:
//!
//! 1. **Boostback** — cancel downrange velocity and head back toward the
//!    landing site (zero for a downrange droneship landing);
//! 2. **Reentry burn** — slow down before the densest air to survive the
//!    aerothermal / aerodynamic loads;
//! 3. **Landing burn** — the [`crate::landing`] hoverslam to touchdown.
//!
//! The total recovery `Δv` is simply their sum, and the propellant it
//! costs follows from the **Tsiolkovsky** rocket equation via
//! [`crate::mass`]: a mass ratio `MR = exp(Δv_total/(Isp·g₀))` and a burned
//! propellant fraction `1 − 1/MR` of the mass at the start of the
//! recovery sequence. This is exactly the reserve that a recoverable
//! first stage must hold back from the ascent — the "recovery tax" on
//! payload.
//!
//! Every relation is an exact algebraic / Tsiolkovsky result, pinned
//! directly by the unit tests.

use serde::{Deserialize, Serialize};

use crate::error::{AstroError, Result};
use crate::mass;

/// The `Δv` components of a booster-recovery sequence (m/s).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct RecoveryBudget {
    /// Boostback-burn `Δv` (m/s). Zero for a downrange landing.
    pub boostback: f64,
    /// Reentry-burn `Δv` (m/s).
    pub reentry_burn: f64,
    /// Landing-burn `Δv` (m/s).
    pub landing_burn: f64,
}

impl RecoveryBudget {
    /// Total recovery `Δv` (m/s): the sum of the three burns.
    ///
    /// # Errors
    ///
    /// Returns [`AstroError::InvalidParameter`] if any component is
    /// non-finite or negative.
    pub fn total_delta_v(&self) -> Result<f64> {
        for value in [self.boostback, self.reentry_burn, self.landing_burn] {
            if !value.is_finite() || value < 0.0 {
                return Err(AstroError::InvalidParameter(
                    "recovery delta-v components must be finite and >= 0",
                ));
            }
        }
        Ok(self.boostback + self.reentry_burn + self.landing_burn)
    }

    /// Propellant mass (kg) the recovery sequence consumes, given the
    /// vehicle mass (kg) at the **start** of the sequence (after stage
    /// separation) and the engine `isp` (s).
    ///
    /// From Tsiolkovsky, the mass ratio over the whole sequence is
    /// `MR = exp(Δv_total/(Isp·g₀))`, so the burned propellant is
    /// `m_start·(1 − 1/MR)`.
    ///
    /// # Errors
    ///
    /// Returns [`AstroError::InvalidParameter`] for a non-physical budget
    /// (see [`RecoveryBudget::total_delta_v`]), a non-finite or
    /// non-positive `mass_at_start`, or a bad `isp` (see
    /// [`crate::mass::mass_ratio`]).
    pub fn propellant_used(&self, mass_at_start: f64, isp: f64) -> Result<f64> {
        if !mass_at_start.is_finite() || mass_at_start <= 0.0 {
            return Err(AstroError::InvalidParameter(
                "mass_at_start must be finite and > 0",
            ));
        }
        let dv = self.total_delta_v()?;
        let mr = mass::mass_ratio(dv, isp)?;
        Ok(mass_at_start * (1.0 - 1.0 / mr))
    }

    /// Fraction of the start-of-sequence mass burned as recovery
    /// propellant: `1 − 1/MR`.
    ///
    /// # Errors
    ///
    /// As [`RecoveryBudget::propellant_used`] (for the budget and `isp`).
    pub fn propellant_fraction(&self, isp: f64) -> Result<f64> {
        let dv = self.total_delta_v()?;
        let mr = mass::mass_ratio(dv, isp)?;
        Ok(1.0 - 1.0 / mr)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::G0;

    #[test]
    fn total_is_the_sum_of_burns() {
        // A Falcon-9-class RTLS budget: boostback ~2000, reentry ~1500,
        // landing ~700 -> 4200 m/s total.
        let b = RecoveryBudget {
            boostback: 2_000.0,
            reentry_burn: 1_500.0,
            landing_burn: 700.0,
        };
        assert!((b.total_delta_v().expect("ok") - 4_200.0).abs() < 1e-9);
        // A downrange (droneship) recovery skips boostback.
        let dr = RecoveryBudget {
            boostback: 0.0,
            ..b
        };
        assert!((dr.total_delta_v().expect("ok") - 2_200.0).abs() < 1e-9);
    }

    #[test]
    fn propellant_matches_tsiolkovsky() {
        // Δv_total = 4200 m/s, Isp = 282 s (Merlin SL), m_start = 30 000 kg.
        // MR = exp(4200/(282·g₀)); burned = m·(1 − 1/MR).
        let b = RecoveryBudget {
            boostback: 2_000.0,
            reentry_burn: 1_500.0,
            landing_burn: 700.0,
        };
        let isp = 282.0;
        let m0 = 30_000.0;
        let mr = (4_200.0 / (isp * G0)).exp();
        let expected = m0 * (1.0 - 1.0 / mr);
        let used = b.propellant_used(m0, isp).expect("ok");
        assert!((used - expected).abs() / expected < 1e-12, "used = {used}");
        // Fraction is consistent with the mass.
        let frac = b.propellant_fraction(isp).expect("ok");
        assert!((frac - used / m0).abs() < 1e-12);
        assert!((0.0..1.0).contains(&frac));
    }

    #[test]
    fn zero_budget_costs_no_propellant() {
        let b = RecoveryBudget {
            boostback: 0.0,
            reentry_burn: 0.0,
            landing_burn: 0.0,
        };
        assert!((b.total_delta_v().expect("ok")).abs() < 1e-12);
        // MR = exp(0) = 1 -> burned fraction 0.
        assert!(b.propellant_used(20_000.0, 300.0).expect("ok").abs() < 1e-9);
        assert!(b.propellant_fraction(300.0).expect("ok").abs() < 1e-12);
    }

    #[test]
    fn more_delta_v_costs_more_propellant() {
        let small = RecoveryBudget {
            boostback: 500.0,
            reentry_burn: 500.0,
            landing_burn: 500.0,
        };
        let big = RecoveryBudget {
            boostback: 2_500.0,
            reentry_burn: 1_500.0,
            landing_burn: 800.0,
        };
        let isp = 300.0;
        let m0 = 25_000.0;
        assert!(
            big.propellant_used(m0, isp).expect("ok")
                > small.propellant_used(m0, isp).expect("ok")
        );
    }

    #[test]
    fn rejects_non_physical_inputs() {
        let good = RecoveryBudget {
            boostback: 1_000.0,
            reentry_burn: 1_000.0,
            landing_burn: 500.0,
        };
        assert!(good.total_delta_v().is_ok());
        let bad = RecoveryBudget {
            boostback: -1.0,
            ..good
        };
        assert!(bad.total_delta_v().is_err());
        assert!(RecoveryBudget {
            reentry_burn: f64::NAN,
            ..good
        }
        .total_delta_v()
        .is_err());
        // Bad mass / isp.
        assert!(good.propellant_used(0.0, 300.0).is_err());
        assert!(good.propellant_used(20_000.0, 0.0).is_err());
    }
}
