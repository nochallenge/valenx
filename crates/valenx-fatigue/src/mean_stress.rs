//! Mean-stress corrections: Goodman, Soderberg, and Gerber.
//!
//! ## Model
//!
//! A stress-life curve is measured under **fully-reversed** loading
//! (mean stress `sm = 0`). Real duty cycles carry a non-zero mean
//! stress, which is damaging: a tensile mean stress reduces the
//! alternating stress amplitude the part can tolerate. The classic
//! constant-life criteria draw a straight failure line on the Haigh
//! diagram (alternating stress `sa` vertical, mean stress `sm`
//! horizontal) between the endurance limit `Se` on the `sa` axis and a
//! static strength on the `sm` axis:
//!
//! ```text
//! Goodman:    sa/Se + sm/Su      = 1/n      (straight line to Su)
//! Soderberg:  sa/Se + sm/Sy      = 1/n      (straight line to Sy)
//! Gerber:     n*sa/Se + (n*sm/Su)^2 = 1     (parabola to Su)
//! ```
//!
//! - `Se` is the (corrected) endurance limit / target fatigue strength.
//! - `Su` is the ultimate tensile strength (Goodman / Gerber intercept).
//! - `Sy` is the yield strength (Soderberg intercept — more
//!   conservative, never permits yielding).
//! - `n` is the design factor of safety (`n = 1` is the failure line).
//! - **Gerber** shares Goodman's `Su` intercept but bows above the
//!   straight line, so for the same mean stress it permits a higher
//!   alternating stress — the least conservative of the three and the
//!   best fit to ductile-material test data.
//!
//! Reading the line two ways:
//!
//! - **Allowable alternating stress** at a given mean stress and design
//!   factor: `sa = (Se / n) * (1 - n*sm/S0)`, where `S0` is `Su`
//!   (Goodman) or `Sy` (Soderberg). At `sm = 0` and `n = 1` this is
//!   exactly `Se`; at `sm = S0` it is `0`.
//! - **Factor of safety** for an applied operating point `(sa, sm)`:
//!   `1/n = sa/Se + sm/S0`, i.e. `n = 1 / (sa/Se + sm/S0)`.
//!
//! ## Honest scope
//!
//! These are the textbook constant-life criteria — the Goodman and
//! Soderberg straight lines plus the Gerber parabola. The crate does not
//! implement the Morrow or Smith-Watson-Topper corrections, and it does
//! not check the separate first-cycle yield line `sa + sm <= Sy`.
//! Research/educational grade, not a production design tool.

use crate::error::{FatigueError, Result};
use serde::{Deserialize, Serialize};

/// Which constant-life criterion to apply.
///
/// The two differ only in the static-strength intercept on the
/// mean-stress axis: Goodman uses the ultimate strength `Su`, Soderberg
/// the (lower) yield strength `Sy`, making Soderberg the more
/// conservative choice.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum MeanStressCriterion {
    /// Goodman line — intercepts the mean-stress axis at the ultimate
    /// tensile strength `Su`.
    Goodman,
    /// Soderberg line — intercepts the mean-stress axis at the yield
    /// strength `Sy` (more conservative; precludes yielding).
    Soderberg,
    /// Gerber **parabola** — `sa/Se + (sm/Su)² = 1/n`, intercepting the
    /// mean-stress axis at the ultimate strength `Su` like Goodman but
    /// bowing *above* the Goodman line, so it permits a higher alternating
    /// stress at the same mean. It is the best fit to experimental data for
    /// ductile materials and the least conservative of the three.
    Gerber,
}

/// The material strengths a mean-stress correction needs.
///
/// All three are stresses in the caller's own consistent unit (MPa,
/// ksi, …); the criteria are dimensionless in that unit. Build one with
/// [`Material::new`], which validates the ordering `Se < Sy <= Su`
/// loosely (it only enforces strict positivity and `Sy <= Su`, since
/// `Se` can in principle approach `Sy`).
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct Material {
    /// Corrected endurance limit / target fatigue strength `Se`.
    pub endurance_limit: f64,
    /// Yield strength `Sy` (Soderberg intercept).
    pub yield_strength: f64,
    /// Ultimate tensile strength `Su` (Goodman intercept).
    pub ultimate_strength: f64,
}

impl Material {
    /// Build a validated material.
    ///
    /// # Errors
    ///
    /// Returns [`FatigueError::Invalid`] if any strength is not strictly
    /// positive / finite, and [`FatigueError::Domain`] if the yield
    /// strength exceeds the ultimate strength (a physically impossible
    /// ordering).
    pub fn new(endurance_limit: f64, yield_strength: f64, ultimate_strength: f64) -> Result<Self> {
        for (name, v) in [
            ("endurance_limit", endurance_limit),
            ("yield_strength", yield_strength),
            ("ultimate_strength", ultimate_strength),
        ] {
            if !v.is_finite() || v <= 0.0 {
                return Err(FatigueError::invalid(
                    name,
                    format!("strength must be finite and > 0, got {v}"),
                ));
            }
        }
        if yield_strength > ultimate_strength {
            return Err(FatigueError::domain(format!(
                "yield strength {yield_strength} must not exceed ultimate \
                 strength {ultimate_strength}"
            )));
        }
        Ok(Material {
            endurance_limit,
            yield_strength,
            ultimate_strength,
        })
    }

    /// The static-strength intercept on the mean-stress axis for a given
    /// criterion: `Su` for Goodman and Gerber, `Sy` for Soderberg.
    fn static_intercept(&self, criterion: MeanStressCriterion) -> f64 {
        match criterion {
            MeanStressCriterion::Goodman | MeanStressCriterion::Gerber => self.ultimate_strength,
            MeanStressCriterion::Soderberg => self.yield_strength,
        }
    }

    /// The allowable alternating stress amplitude at a given mean stress
    /// and design factor of safety.
    ///
    /// For the straight-line criteria solves `sa/Se + sm/S0 = 1/n`, giving
    /// `sa = (Se / n) * (1 - n * sm / S0)`; for [`Gerber`] it solves the
    /// parabola `n*sa/Se + (n*sm/Su)^2 = 1`, giving
    /// `sa = (Se / n) * (1 - (n * sm / Su)^2)`. With `n = 1` and `sm = 0`
    /// every criterion returns exactly the endurance limit `Se`; with
    /// `sm = S0` (the static intercept) it returns `0`.
    ///
    /// [`Gerber`]: MeanStressCriterion::Gerber
    ///
    /// # Errors
    ///
    /// Returns [`FatigueError::Invalid`] for a non-positive / non-finite
    /// design factor, and [`FatigueError::Domain`] if `mean_stress` is
    /// negative or if the operating mean stress already meets or exceeds
    /// the static line at this design factor (so no positive alternating
    /// stress is admissible).
    pub fn allowable_alternating(
        &self,
        criterion: MeanStressCriterion,
        mean_stress: f64,
        design_factor: f64,
    ) -> Result<f64> {
        if !design_factor.is_finite() || design_factor <= 0.0 {
            return Err(FatigueError::invalid(
                "design_factor",
                format!("must be finite and > 0, got {design_factor}"),
            ));
        }
        if !mean_stress.is_finite() || mean_stress < 0.0 {
            return Err(FatigueError::invalid(
                "mean_stress",
                format!("must be finite and >= 0, got {mean_stress}"),
            ));
        }
        let s0 = self.static_intercept(criterion);
        let sa = match criterion {
            // Gerber is parabolic in the mean term: n*sa/Se + (n*sm/Su)^2 = 1.
            MeanStressCriterion::Gerber => {
                let r = design_factor * mean_stress / s0;
                (self.endurance_limit / design_factor) * (1.0 - r * r)
            }
            _ => (self.endurance_limit / design_factor) * (1.0 - design_factor * mean_stress / s0),
        };
        if sa < 0.0 {
            return Err(FatigueError::domain(format!(
                "mean stress {mean_stress} meets or exceeds the static line \
                 (intercept {s0}) at design factor {design_factor}; no \
                 positive alternating stress is admissible"
            )));
        }
        Ok(sa)
    }

    /// The factor of safety for an applied operating point `(sa, sm)`.
    ///
    /// For the straight-line criteria evaluates `n = 1 / (sa/Se + sm/S0)`;
    /// for [`Gerber`](MeanStressCriterion::Gerber) it takes the positive
    /// root of the quadratic `(sm/Su)^2 n^2 + (sa/Se) n - 1 = 0`. A result
    /// `n >= 1` means the point lies on or inside the constant-life line;
    /// `n < 1` means it lies outside (predicted failure).
    ///
    /// # Errors
    ///
    /// Returns [`FatigueError::Invalid`] if `alternating_stress` is not
    /// strictly positive / finite or if `mean_stress` is negative /
    /// non-finite.
    pub fn factor_of_safety(
        &self,
        criterion: MeanStressCriterion,
        alternating_stress: f64,
        mean_stress: f64,
    ) -> Result<f64> {
        if !alternating_stress.is_finite() || alternating_stress <= 0.0 {
            return Err(FatigueError::invalid(
                "alternating_stress",
                format!("must be finite and > 0, got {alternating_stress}"),
            ));
        }
        if !mean_stress.is_finite() || mean_stress < 0.0 {
            return Err(FatigueError::invalid(
                "mean_stress",
                format!("must be finite and >= 0, got {mean_stress}"),
            ));
        }
        let s0 = self.static_intercept(criterion);
        let n = match criterion {
            // Gerber: solve n*sa/Se + (n*sm/Su)^2 = 1, a quadratic
            // B n^2 + A n - 1 = 0 with A = sa/Se, B = (sm/Su)^2; take the
            // positive root. At zero mean it degenerates to n = Se/sa.
            MeanStressCriterion::Gerber => {
                if mean_stress > 0.0 {
                    let a = alternating_stress / self.endurance_limit;
                    let b = (mean_stress / s0).powi(2);
                    (-a + (a * a + 4.0 * b).sqrt()) / (2.0 * b)
                } else {
                    self.endurance_limit / alternating_stress
                }
            }
            _ => 1.0 / (alternating_stress / self.endurance_limit + mean_stress / s0),
        };
        Ok(n)
    }

    /// The **equivalent completely-reversed stress** `sigma_ar` for an
    /// operating point `(sa, sm)` under the chosen criterion: the
    /// fully-reversed (zero-mean) alternating stress that is equally
    /// damaging. This is the value to feed into a fully-reversed S-N curve
    /// ([`SnCurve::cycles_to_failure`](crate::SnCurve::cycles_to_failure))
    /// to predict finite life under a non-zero mean stress — the bridge
    /// between the mean-stress correction and the stress-life curve.
    ///
    /// Straight-line (Goodman / Soderberg): `sigma_ar = sa / (1 - sm/S0)`;
    /// Gerber parabola: `sigma_ar = sa / (1 - (sm/Su)^2)`, where `S0` is
    /// `Su` (Goodman / Gerber) or `Sy` (Soderberg). At zero mean stress it
    /// returns `sa` unchanged; a tensile mean raises it (mean stress is
    /// damaging); and it equals the endurance limit `Se` exactly when the
    /// point lies on the `n = 1` constant-life line.
    ///
    /// # Errors
    ///
    /// Returns [`FatigueError::Invalid`] if `alternating_stress` is not
    /// strictly positive / finite or `mean_stress` is negative /
    /// non-finite, and [`FatigueError::Domain`] if the mean stress meets or
    /// exceeds the static intercept (the knockdown factor is `<= 0`, so no
    /// finite reversed equivalent exists).
    pub fn equivalent_reversed_stress(
        &self,
        criterion: MeanStressCriterion,
        alternating_stress: f64,
        mean_stress: f64,
    ) -> Result<f64> {
        if !alternating_stress.is_finite() || alternating_stress <= 0.0 {
            return Err(FatigueError::invalid(
                "alternating_stress",
                format!("must be finite and > 0, got {alternating_stress}"),
            ));
        }
        if !mean_stress.is_finite() || mean_stress < 0.0 {
            return Err(FatigueError::invalid(
                "mean_stress",
                format!("must be finite and >= 0, got {mean_stress}"),
            ));
        }
        let s0 = self.static_intercept(criterion);
        let knockdown = match criterion {
            MeanStressCriterion::Gerber => {
                let r = mean_stress / self.ultimate_strength;
                1.0 - r * r
            }
            _ => 1.0 - mean_stress / s0,
        };
        if knockdown <= 0.0 {
            return Err(FatigueError::domain(format!(
                "mean stress {mean_stress} meets or exceeds the static \
                 intercept {s0}; no finite completely-reversed equivalent \
                 stress exists"
            )));
        }
        Ok(alternating_stress / knockdown)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(x: f64, y: f64) {
        let tol = 1e-9 * x.abs().max(y.abs()).max(1.0);
        assert!((x - y).abs() < tol, "expected {x} ~= {y}");
    }

    fn steel() -> Material {
        // Se = 200, Sy = 350, Su = 500 (consistent units).
        Material::new(200.0, 350.0, 500.0).unwrap()
    }

    #[test]
    fn constructor_validates() {
        assert!(Material::new(0.0, 350.0, 500.0).is_err());
        assert!(Material::new(200.0, -1.0, 500.0).is_err());
        assert!(Material::new(200.0, 350.0, f64::NAN).is_err());
        // Sy > Su is impossible.
        assert!(Material::new(200.0, 600.0, 500.0).is_err());
        assert!(Material::new(200.0, 350.0, 500.0).is_ok());
    }

    /// Goodman anchor point 1: sm = 0 at n = 1 gives sa = Se.
    #[test]
    fn goodman_zero_mean_gives_endurance_limit() {
        let m = steel();
        let sa = m
            .allowable_alternating(MeanStressCriterion::Goodman, 0.0, 1.0)
            .unwrap();
        close(sa, m.endurance_limit);
    }

    /// Goodman anchor point 2: sm = Su gives sa = 0.
    #[test]
    fn goodman_mean_at_ultimate_gives_zero() {
        let m = steel();
        let sa = m
            .allowable_alternating(MeanStressCriterion::Goodman, m.ultimate_strength, 1.0)
            .unwrap();
        close(sa, 0.0);
    }

    /// Soderberg anchor point 2: sm = Sy gives sa = 0.
    #[test]
    fn soderberg_mean_at_yield_gives_zero() {
        let m = steel();
        let sa = m
            .allowable_alternating(MeanStressCriterion::Soderberg, m.yield_strength, 1.0)
            .unwrap();
        close(sa, 0.0);
    }

    /// Soderberg at zero mean also recovers the endurance limit.
    #[test]
    fn soderberg_zero_mean_gives_endurance_limit() {
        let m = steel();
        let sa = m
            .allowable_alternating(MeanStressCriterion::Soderberg, 0.0, 1.0)
            .unwrap();
        close(sa, m.endurance_limit);
    }

    /// Higher mean stress lowers the allowable alternating stress.
    #[test]
    fn higher_mean_lowers_allowable() {
        let m = steel();
        let lo = m
            .allowable_alternating(MeanStressCriterion::Goodman, 100.0, 1.0)
            .unwrap();
        let hi = m
            .allowable_alternating(MeanStressCriterion::Goodman, 300.0, 1.0)
            .unwrap();
        assert!(hi < lo, "sa at higher mean ({hi}) should be below ({lo})");
    }

    /// Soderberg is more conservative than Goodman at the same mean
    /// stress (lower allowable alternating stress, since Sy < Su).
    #[test]
    fn soderberg_more_conservative_than_goodman() {
        let m = steel();
        let g = m
            .allowable_alternating(MeanStressCriterion::Goodman, 150.0, 1.0)
            .unwrap();
        let s = m
            .allowable_alternating(MeanStressCriterion::Soderberg, 150.0, 1.0)
            .unwrap();
        assert!(s < g, "Soderberg {s} should be below Goodman {g}");
    }

    /// A known mid-line Goodman value: sm = Su/2 at n = 1 gives sa = Se/2.
    #[test]
    fn goodman_midpoint_value() {
        let m = steel();
        let sa = m
            .allowable_alternating(MeanStressCriterion::Goodman, m.ultimate_strength / 2.0, 1.0)
            .unwrap();
        close(sa, m.endurance_limit / 2.0);
    }

    /// A design factor of n scales the allowable alternating stress at
    /// zero mean to Se/n.
    #[test]
    fn design_factor_scales_allowable_at_zero_mean() {
        let m = steel();
        let sa = m
            .allowable_alternating(MeanStressCriterion::Goodman, 0.0, 2.0)
            .unwrap();
        close(sa, m.endurance_limit / 2.0);
    }

    /// Factor of safety inverts the allowable relation: an operating
    /// point exactly on the failure line reports n = 1.
    #[test]
    fn factor_of_safety_on_line_is_one() {
        let m = steel();
        // Pick sm = 200, the Goodman-allowable sa at n = 1 there:
        let sm = 200.0;
        let sa = m
            .allowable_alternating(MeanStressCriterion::Goodman, sm, 1.0)
            .unwrap();
        let n = m
            .factor_of_safety(MeanStressCriterion::Goodman, sa, sm)
            .unwrap();
        close(n, 1.0);
    }

    /// Pure alternating loading at the endurance limit gives n = 1.
    #[test]
    fn factor_of_safety_pure_alternating() {
        let m = steel();
        let n = m
            .factor_of_safety(MeanStressCriterion::Goodman, m.endurance_limit, 0.0)
            .unwrap();
        close(n, 1.0);
    }

    /// A point inside the line has n > 1; outside has n < 1.
    #[test]
    fn factor_of_safety_inside_and_outside() {
        let m = steel();
        let inside = m
            .factor_of_safety(MeanStressCriterion::Goodman, 50.0, 50.0)
            .unwrap();
        let outside = m
            .factor_of_safety(MeanStressCriterion::Goodman, 250.0, 300.0)
            .unwrap();
        assert!(inside > 1.0, "inside point n = {inside} should exceed 1");
        assert!(
            outside < 1.0,
            "outside point n = {outside} should be below 1"
        );
    }

    #[test]
    fn rejects_bad_inputs() {
        let m = steel();
        assert!(m
            .allowable_alternating(MeanStressCriterion::Goodman, 0.0, 0.0)
            .is_err());
        assert!(m
            .allowable_alternating(MeanStressCriterion::Goodman, -10.0, 1.0)
            .is_err());
        // Mean stress beyond the static intercept -> no admissible sa.
        assert!(m
            .allowable_alternating(MeanStressCriterion::Goodman, 600.0, 1.0)
            .is_err());
        assert!(m
            .factor_of_safety(MeanStressCriterion::Goodman, 0.0, 100.0)
            .is_err());
        assert!(m
            .factor_of_safety(MeanStressCriterion::Goodman, 100.0, -1.0)
            .is_err());
    }

    // --- Gerber parabola -------------------------------------------------

    /// Gerber shares the two anchor points: sa = Se at sm = 0 and sa = 0
    /// at sm = Su; and at sm = Su/2 the parabola gives sa = Se*(1 - 1/4).
    #[test]
    fn gerber_anchor_points_and_midpoint() {
        let m = steel();
        close(
            m.allowable_alternating(MeanStressCriterion::Gerber, 0.0, 1.0)
                .unwrap(),
            m.endurance_limit,
        );
        close(
            m.allowable_alternating(MeanStressCriterion::Gerber, m.ultimate_strength, 1.0)
                .unwrap(),
            0.0,
        );
        // sm = Su/2 = 250 -> sa = 200*(1 - 0.25) = 150.
        close(
            m.allowable_alternating(MeanStressCriterion::Gerber, m.ultimate_strength / 2.0, 1.0)
                .unwrap(),
            0.75 * m.endurance_limit,
        );
    }

    /// The Gerber parabola lies above the Goodman line: for the same mean
    /// stress it permits a higher (less conservative) alternating stress.
    #[test]
    fn gerber_is_less_conservative_than_goodman() {
        let m = steel();
        for &sm in &[50.0, 150.0, 300.0, 450.0] {
            let g = m
                .allowable_alternating(MeanStressCriterion::Goodman, sm, 1.0)
                .unwrap();
            let gerber = m
                .allowable_alternating(MeanStressCriterion::Gerber, sm, 1.0)
                .unwrap();
            assert!(
                gerber > g,
                "Gerber {gerber} should exceed Goodman {g} at sm={sm}"
            );
            assert!(gerber <= m.endurance_limit);
        }
    }

    /// Factor of safety inverts the Gerber allowable: a point produced at
    /// n = 1 reports n = 1, and a point produced at n = 2 reports n = 2.
    #[test]
    fn gerber_factor_of_safety_round_trips() {
        let m = steel();
        for &(sm, nd) in &[(250.0, 1.0), (150.0, 1.0), (100.0, 2.0), (300.0, 1.5)] {
            let sa = m
                .allowable_alternating(MeanStressCriterion::Gerber, sm, nd)
                .unwrap();
            let n = m
                .factor_of_safety(MeanStressCriterion::Gerber, sa, sm)
                .unwrap();
            close(n, nd);
        }
    }

    /// At zero mean the Gerber factor of safety reduces to Se/sa.
    #[test]
    fn gerber_pure_alternating_factor_of_safety() {
        let m = steel();
        close(
            m.factor_of_safety(MeanStressCriterion::Gerber, m.endurance_limit, 0.0)
                .unwrap(),
            1.0,
        );
        close(
            m.factor_of_safety(MeanStressCriterion::Gerber, m.endurance_limit / 2.0, 0.0)
                .unwrap(),
            2.0,
        );
    }

    // --- equivalent completely-reversed stress ---------------------------

    const ALL_CRITERIA: [MeanStressCriterion; 3] = [
        MeanStressCriterion::Goodman,
        MeanStressCriterion::Soderberg,
        MeanStressCriterion::Gerber,
    ];

    /// Zero mean stress leaves the alternating stress unchanged for every
    /// criterion: sigma_ar = sa.
    #[test]
    fn equivalent_reversed_at_zero_mean_is_alternating() {
        let m = steel();
        for c in ALL_CRITERIA {
            close(m.equivalent_reversed_stress(c, 120.0, 0.0).unwrap(), 120.0);
        }
    }

    /// A point on the n = 1 constant-life line maps to sigma_ar = Se, for
    /// every criterion — the defining link to allowable_alternating.
    #[test]
    fn equivalent_reversed_on_line_equals_endurance_limit() {
        let m = steel();
        for c in ALL_CRITERIA {
            for &sm in &[0.0_f64, 100.0, 200.0] {
                let sa = m.allowable_alternating(c, sm, 1.0).unwrap();
                close(
                    m.equivalent_reversed_stress(c, sa, sm).unwrap(),
                    m.endurance_limit,
                );
            }
        }
    }

    /// sigma_ar = Se * sa / allowable_alternating(sm, n=1): the equivalent
    /// reversed stress scales linearly with sa relative to the on-line
    /// allowable. Holds for all three criteria.
    #[test]
    fn equivalent_reversed_matches_allowable_ratio() {
        let m = steel();
        let sm = 120.0;
        let sa = 90.0;
        for c in ALL_CRITERIA {
            let allow = m.allowable_alternating(c, sm, 1.0).unwrap();
            let expected = m.endurance_limit * sa / allow;
            close(m.equivalent_reversed_stress(c, sa, sm).unwrap(), expected);
        }
    }

    /// A tensile mean stress raises the equivalent reversed stress above
    /// the raw alternating stress (mean stress is damaging), and Soderberg
    /// (Sy intercept) is more conservative than Goodman (Su intercept).
    #[test]
    fn equivalent_reversed_increases_with_mean_and_soderberg_is_higher() {
        let m = steel();
        let sa = 100.0;
        let sm = 150.0;
        let g = m
            .equivalent_reversed_stress(MeanStressCriterion::Goodman, sa, sm)
            .unwrap();
        let s = m
            .equivalent_reversed_stress(MeanStressCriterion::Soderberg, sa, sm)
            .unwrap();
        assert!(g > sa, "Goodman sigma_ar {g} should exceed sa {sa}");
        assert!(s > g, "Soderberg {s} should exceed Goodman {g}");
    }

    /// Mean stress at/over the static intercept has no finite reversed
    /// equivalent; bad sa / sm are rejected.
    #[test]
    fn equivalent_reversed_rejects_bad_inputs() {
        let m = steel();
        // sm = Su (Goodman intercept) -> knockdown 0 -> Domain.
        assert!(m
            .equivalent_reversed_stress(MeanStressCriterion::Goodman, 100.0, m.ultimate_strength)
            .is_err());
        assert!(m
            .equivalent_reversed_stress(MeanStressCriterion::Goodman, 0.0, 100.0)
            .is_err());
        assert!(m
            .equivalent_reversed_stress(MeanStressCriterion::Goodman, 100.0, -1.0)
            .is_err());
    }

    /// End-to-end with the S-N curve: a non-zero mean stress shortens the
    /// predicted life. Convert (sa, sm) to sigma_ar, then read the Basquin
    /// life — it is shorter than the pure-alternating life at the same sa.
    #[test]
    fn equivalent_reversed_shortens_predicted_life() {
        use crate::sn::SnCurve;
        let m = steel();
        // Basquin fit through (1e3, 0.9*Su=450) and (1e6, Se=200).
        let curve = SnCurve::from_two_points(1.0e3, 450.0, 1.0e6, 200.0).unwrap();
        let sa = 180.0;
        let sm = 150.0;
        let ar = m
            .equivalent_reversed_stress(MeanStressCriterion::Goodman, sa, sm)
            .unwrap();
        assert!(ar > sa);
        let life_mean = curve.cycles_to_failure_unbounded(ar).unwrap();
        let life_pure = curve.cycles_to_failure_unbounded(sa).unwrap();
        assert!(
            life_mean < life_pure,
            "mean-stress life {life_mean} should be below pure-alternating {life_pure}"
        );
    }
}
