//! Simply-supported beam reactions (pin + roller).
//!
//! The canonical determinate problem of introductory statics: a
//! straight horizontal beam resting on two supports —
//!
//! - a **pin** at the left support position `a`, which can react both a
//!   vertical and a horizontal force (2 unknowns, `R_ax` and `R_ay`);
//! - a **roller** at the right support position `b > a`, which reacts
//!   only a vertical force perpendicular to the beam (1 unknown,
//!   `R_by`).
//!
//! That is **3 reaction unknowns**, matched by the **3 planar
//! equilibrium equations** `sum Fx = sum Fy = sum M = 0`, so the beam is
//! statically determinate and the reactions follow in closed form.
//!
//! ## Solution
//!
//! For vertical point loads `P_i` (taken positive **downward**) applied
//! at positions `x_i`, and horizontal point loads `H_i`:
//!
//! - **Horizontal balance** — only the pin resists horizontal load:
//!   ```text
//!   R_ax = -sum H_i.
//!   ```
//! - **Moment about the pin** `a` — the roller's lever arm is
//!   `(b - a)`; each downward load `P_i` makes a clockwise moment with
//!   arm `(x_i - a)`:
//!   ```text
//!   R_by = sum P_i (x_i - a) / (b - a).
//!   ```
//! - **Vertical balance** — recovers the pin's vertical reaction:
//!   ```text
//!   R_ay = sum P_i - R_by.
//!   ```
//!
//! This reproduces the textbook **lever rule**: a load close to a
//! support throws most of its weight onto that support, and a load at
//! mid-span splits evenly (`P/2` to each).
//!
//! Reactions are reported as **upward** vertical forces (the usual sign
//! convention for a downward-loaded beam) and a signed horizontal pin
//! force.

use crate::equilibrium::ForceSystem;
use crate::error::{Result, StaticsError};
use crate::force::AppliedForce;

/// A vertical and/or horizontal point load applied to the beam.
///
/// `position` is the distance along the beam (same axis as the support
/// positions). `down` is the **downward** vertical magnitude (positive
/// pulls the beam down, as gravity loads usually do; a negative value
/// is an uplift). `right` is an optional horizontal component (positive
/// toward `+x`); most beam problems leave it `0`.
#[derive(Debug, Clone, Copy, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct PointLoad {
    /// Position along the beam axis.
    pub position: f64,
    /// Downward vertical magnitude (positive = downward).
    pub down: f64,
    /// Horizontal magnitude (positive = toward `+x`).
    pub right: f64,
}

impl PointLoad {
    /// A purely vertical downward load of `magnitude` at `position`.
    #[must_use]
    pub fn vertical(position: f64, magnitude: f64) -> Self {
        Self {
            position,
            down: magnitude,
            right: 0.0,
        }
    }

    /// A general load with both a downward and a horizontal component.
    #[must_use]
    pub fn new(position: f64, down: f64, right: f64) -> Self {
        Self {
            position,
            down,
            right,
        }
    }
}

/// The solved support reactions of a [`SimpleBeam`].
///
/// All three components follow the same conventions as the inputs:
/// `pin_vertical` and `roller_vertical` are **upward** forces, and
/// `pin_horizontal` is positive toward `+x`.
#[derive(Debug, Clone, Copy, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct Reactions {
    /// Upward vertical reaction at the pin (left support), `R_ay`.
    pub pin_vertical: f64,
    /// Horizontal reaction at the pin, `R_ax` (positive toward `+x`).
    pub pin_horizontal: f64,
    /// Upward vertical reaction at the roller (right support), `R_by`.
    pub roller_vertical: f64,
}

impl Reactions {
    /// Total upward vertical reaction `R_ay + R_by`. In equilibrium
    /// this equals the total downward load.
    #[must_use]
    pub fn total_vertical(&self) -> f64 {
        self.pin_vertical + self.roller_vertical
    }
}

/// The single equivalent force of a beam's vertical loads.
///
/// A set of parallel vertical point loads reduces to one **resultant**:
/// a downward force equal to their sum, acting through their centroid (the
/// load's "centre of gravity"). Returned by
/// [`SimpleBeam::vertical_load_resultant`].
#[derive(Debug, Clone, Copy, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct LoadResultant {
    /// Total downward magnitude `sum P_i` (positive = downward).
    pub down: f64,
    /// Line of action `x_bar = sum P_i x_i / sum P_i`, the load centroid.
    pub position: f64,
}

/// A simply-supported beam: a pin at `support_a` and a roller at
/// `support_b`, carrying any number of [`PointLoad`]s.
///
/// Construct with [`SimpleBeam::new`] (which validates the geometry),
/// add loads, then call [`SimpleBeam::reactions`].
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct SimpleBeam {
    /// Pin (left) support position along the beam axis.
    pub support_a: f64,
    /// Roller (right) support position; must satisfy `support_b > support_a`.
    pub support_b: f64,
    /// Applied point loads.
    pub loads: Vec<PointLoad>,
}

impl SimpleBeam {
    /// Create a beam with a pin at `support_a` and a roller at
    /// `support_b`.
    ///
    /// # Errors
    ///
    /// - [`StaticsError::Invalid`] if either support position is
    ///   non-finite.
    /// - [`StaticsError::Degenerate`] if the supports coincide
    ///   (`support_b <= support_a`), which leaves zero lever arm
    ///   between them and no unique moment solution.
    pub fn new(support_a: f64, support_b: f64) -> Result<Self> {
        if !support_a.is_finite() {
            return Err(StaticsError::invalid("support_a", "must be finite"));
        }
        if !support_b.is_finite() {
            return Err(StaticsError::invalid("support_b", "must be finite"));
        }
        if support_b <= support_a {
            return Err(StaticsError::degenerate(format!(
                "roller at {support_b} must be strictly right of pin at {support_a}"
            )));
        }
        Ok(Self {
            support_a,
            support_b,
            loads: Vec::new(),
        })
    }

    /// The distance between the two supports, `support_b - support_a`.
    /// Always strictly positive for a validly constructed beam.
    #[must_use]
    pub fn span(&self) -> f64 {
        self.support_b - self.support_a
    }

    /// Add a load and return `&mut self` for chaining.
    ///
    /// # Errors
    ///
    /// [`StaticsError::Invalid`] if any load field is non-finite. The
    /// load *may* sit outside `[support_a, support_b]` (an overhang),
    /// which is physically meaningful, so off-span positions are
    /// allowed.
    pub fn add_load(&mut self, load: PointLoad) -> Result<&mut Self> {
        if !load.position.is_finite() {
            return Err(StaticsError::invalid("load.position", "must be finite"));
        }
        if !load.down.is_finite() {
            return Err(StaticsError::invalid("load.down", "must be finite"));
        }
        if !load.right.is_finite() {
            return Err(StaticsError::invalid("load.right", "must be finite"));
        }
        self.loads.push(load);
        Ok(self)
    }

    /// Convenience: add a purely vertical downward load.
    ///
    /// # Errors
    ///
    /// Same as [`add_load`](Self::add_load).
    pub fn add_vertical(&mut self, position: f64, magnitude: f64) -> Result<&mut Self> {
        self.add_load(PointLoad::vertical(position, magnitude))
    }

    /// Total downward vertical load `sum P_i`.
    #[must_use]
    pub fn total_down(&self) -> f64 {
        self.loads.iter().map(|l| l.down).sum()
    }

    /// Total horizontal load `sum H_i` (positive toward `+x`).
    #[must_use]
    pub fn total_right(&self) -> f64 {
        self.loads.iter().map(|l| l.right).sum()
    }

    /// Solve the support reactions in closed form.
    ///
    /// Takes moments about the pin to get the roller reaction, then
    /// uses vertical balance for the pin reaction and horizontal
    /// balance for the pin's horizontal force. See the [module
    /// docs](self) for the equations.
    ///
    /// The span has already been validated as strictly positive by
    /// [`SimpleBeam::new`], so this method does not fail.
    #[must_use]
    pub fn reactions(&self) -> Reactions {
        let span = self.span();

        // Moment about the pin (a): sum of clockwise load moments is
        // balanced by the roller's upward reaction over the span.
        //   R_by * span = sum P_i * (x_i - a)
        let moment_about_pin: f64 = self
            .loads
            .iter()
            .map(|l| l.down * (l.position - self.support_a))
            .sum();
        let roller_vertical = moment_about_pin / span;

        // Vertical balance: R_ay + R_by = sum P_i.
        let pin_vertical = self.total_down() - roller_vertical;

        // Horizontal balance: only the pin resists horizontal load.
        // The reaction opposes the applied horizontal load.
        let pin_horizontal = -self.total_right();

        Reactions {
            pin_vertical,
            pin_horizontal,
            roller_vertical,
        }
    }

    /// Build the full [`ForceSystem`] of applied loads **plus** the
    /// solved reactions, expressed as [`AppliedForce`]s in the `xy`
    /// plane (`x` along the beam, `y` up). Useful for independently
    /// checking that the solved beam satisfies `sum Fx = sum Fy = sum M
    /// = 0`.
    ///
    /// Downward loads become forces with `Fy = -down`; upward reactions
    /// become forces with `Fy = +reaction`.
    #[must_use]
    pub fn equilibrium_system(&self) -> ForceSystem {
        let r = self.reactions();
        let mut sys = ForceSystem::new();

        for l in &self.loads {
            sys.push(AppliedForce::new(l.position, 0.0, l.right, -l.down));
        }
        // Pin reaction at (support_a, 0).
        sys.push(AppliedForce::new(
            self.support_a,
            0.0,
            r.pin_horizontal,
            r.pin_vertical,
        ));
        // Roller reaction at (support_b, 0), vertical only.
        sys.push(AppliedForce::new(
            self.support_b,
            0.0,
            0.0,
            r.roller_vertical,
        ));
        sys
    }

    /// The single equivalent (resultant) of all the beam's vertical loads:
    /// a downward force `sum P_i` acting through the load centroid
    /// `x_bar = sum P_i x_i / sum P_i`.
    ///
    /// Replacing the whole load set with this one force leaves the support
    /// [`Reactions`] unchanged, because both the total vertical load and
    /// its moment about any point are preserved — the defining property of
    /// a force-system resultant.
    ///
    /// Returns `None` when the net vertical load is zero (`sum P_i = 0`):
    /// the loads then form a couple or cancel, and there is no single line
    /// of action.
    #[must_use]
    pub fn vertical_load_resultant(&self) -> Option<LoadResultant> {
        let total = self.total_down();
        if total == 0.0 {
            return None;
        }
        let moment: f64 = self.loads.iter().map(|l| l.down * l.position).sum();
        Some(LoadResultant {
            down: total,
            position: moment / total,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    #[test]
    fn central_load_splits_evenly() {
        // VALIDATE: a central point load P gives equal reactions P/2.
        let mut beam = SimpleBeam::new(0.0, 10.0).unwrap();
        beam.add_vertical(5.0, 100.0).unwrap();
        let r = beam.reactions();
        assert!(
            (r.pin_vertical - 50.0).abs() < EPS,
            "pin {}",
            r.pin_vertical
        );
        assert!(
            (r.roller_vertical - 50.0).abs() < EPS,
            "roller {}",
            r.roller_vertical
        );
    }

    #[test]
    fn reactions_balance_total_load() {
        // VALIDATE: sum Fy = 0 -> R_ay + R_by equals total downward load,
        // for an arbitrary off-centre load.
        let mut beam = SimpleBeam::new(0.0, 8.0).unwrap();
        beam.add_vertical(3.0, 60.0).unwrap();
        let r = beam.reactions();
        assert!(
            (r.total_vertical() - 60.0).abs() < EPS,
            "total {}",
            r.total_vertical()
        );
    }

    #[test]
    fn load_nearer_a_support_loads_it_more() {
        // VALIDATE (lever rule): a load close to the right support puts
        // more of its weight on the right (roller) reaction.
        let mut beam = SimpleBeam::new(0.0, 10.0).unwrap();
        beam.add_vertical(8.0, 100.0).unwrap(); // 8 m from pin, 2 m from roller
        let r = beam.reactions();
        // Roller carries the larger share; pin the smaller.
        assert!(
            r.roller_vertical > r.pin_vertical,
            "roller {} should exceed pin {}",
            r.roller_vertical,
            r.pin_vertical
        );
        // Exact lever rule: R_by = P * a_from_pin / span = 100*8/10 = 80;
        // R_ay = 100 * 2/10 = 20.
        assert!(
            (r.roller_vertical - 80.0).abs() < EPS,
            "roller {}",
            r.roller_vertical
        );
        assert!(
            (r.pin_vertical - 20.0).abs() < EPS,
            "pin {}",
            r.pin_vertical
        );
    }

    #[test]
    fn lever_rule_symmetry() {
        // A load at distance d from the pin and the mirror load at
        // distance d from the roller swap the two reaction values.
        let mut left = SimpleBeam::new(0.0, 12.0).unwrap();
        left.add_vertical(3.0, 90.0).unwrap();
        let rl = left.reactions();

        let mut right = SimpleBeam::new(0.0, 12.0).unwrap();
        right.add_vertical(9.0, 90.0).unwrap(); // mirror position (12 - 3)
        let rr = right.reactions();

        assert!((rl.pin_vertical - rr.roller_vertical).abs() < EPS);
        assert!((rl.roller_vertical - rr.pin_vertical).abs() < EPS);
    }

    #[test]
    fn solved_beam_satisfies_all_three_equilibrium_equations() {
        // VALIDATE: at the solution, sum Fx = sum Fy = sum M = 0.
        let mut beam = SimpleBeam::new(1.0, 7.0).unwrap();
        beam.add_load(PointLoad::new(2.0, 50.0, 12.0)).unwrap();
        beam.add_load(PointLoad::new(5.5, 30.0, -4.0)).unwrap();
        beam.add_vertical(6.0, 20.0).unwrap();

        let sys = beam.equilibrium_system();
        assert!(
            sys.is_in_equilibrium(1e-6),
            "sum Fx {}, sum Fy {}, sum M {}",
            sys.sum_fx(),
            sys.sum_fy(),
            sys.sum_moment_origin()
        );
    }

    #[test]
    fn pin_carries_all_horizontal_load() {
        // The roller takes no horizontal force; the pin reaction
        // exactly opposes the applied horizontal load.
        let mut beam = SimpleBeam::new(0.0, 4.0).unwrap();
        beam.add_load(PointLoad::new(2.0, 0.0, 25.0)).unwrap();
        let r = beam.reactions();
        assert!(
            (r.pin_horizontal + 25.0).abs() < EPS,
            "pin_h {}",
            r.pin_horizontal
        );
        // No vertical load -> both vertical reactions vanish.
        assert!(r.pin_vertical.abs() < EPS);
        assert!(r.roller_vertical.abs() < EPS);
    }

    #[test]
    fn two_loads_superpose() {
        // Reactions are linear: two loads give the sum of what each
        // would give alone. Loads of 40 at 2.5 and 60 at 7.5 on a
        // 10 m beam.
        let mut beam = SimpleBeam::new(0.0, 10.0).unwrap();
        beam.add_vertical(2.5, 40.0).unwrap();
        beam.add_vertical(7.5, 60.0).unwrap();
        let r = beam.reactions();
        // R_by = (40*2.5 + 60*7.5)/10 = (100 + 450)/10 = 55.
        // R_ay = 100 - 55 = 45.
        assert!(
            (r.roller_vertical - 55.0).abs() < EPS,
            "roller {}",
            r.roller_vertical
        );
        assert!(
            (r.pin_vertical - 45.0).abs() < EPS,
            "pin {}",
            r.pin_vertical
        );
        assert!((r.total_vertical() - 100.0).abs() < EPS);
    }

    #[test]
    fn coincident_supports_are_degenerate() {
        let err = SimpleBeam::new(3.0, 3.0).unwrap_err();
        assert_eq!(err.code(), "statics.degenerate");
    }

    #[test]
    fn reversed_supports_are_degenerate() {
        let err = SimpleBeam::new(5.0, 1.0).unwrap_err();
        assert_eq!(err.code(), "statics.degenerate");
    }

    #[test]
    fn non_finite_support_is_invalid() {
        let err = SimpleBeam::new(f64::INFINITY, 2.0).unwrap_err();
        assert_eq!(err.code(), "statics.invalid");
    }

    #[test]
    fn non_finite_load_is_rejected() {
        let mut beam = SimpleBeam::new(0.0, 5.0).unwrap();
        let err = beam.add_vertical(f64::NAN, 10.0).unwrap_err();
        assert_eq!(err.code(), "statics.invalid");
    }

    #[test]
    fn span_is_distance_between_supports() {
        let beam = SimpleBeam::new(2.0, 9.5).unwrap();
        assert!((beam.span() - 7.5).abs() < EPS);
    }

    // -- vertical load resultant -------------------------------------

    #[test]
    fn single_load_resultant_is_that_load() {
        let mut beam = SimpleBeam::new(0.0, 10.0).unwrap();
        beam.add_vertical(3.7, 42.0).unwrap();
        let res = beam.vertical_load_resultant().unwrap();
        assert!((res.down - 42.0).abs() < EPS, "down {}", res.down);
        assert!((res.position - 3.7).abs() < EPS, "x_bar {}", res.position);
    }

    #[test]
    fn resultant_centroid_matches_formula() {
        // x_bar = sum(P_i x_i) / sum(P_i): (40*2 + 60*7) / 100 = 500/100 = 5.
        let mut beam = SimpleBeam::new(0.0, 10.0).unwrap();
        beam.add_vertical(2.0, 40.0).unwrap();
        beam.add_vertical(7.0, 60.0).unwrap();
        let res = beam.vertical_load_resultant().unwrap();
        assert!((res.down - 100.0).abs() < EPS, "down {}", res.down);
        assert!((res.position - 5.0).abs() < EPS, "x_bar {}", res.position);
    }

    #[test]
    fn two_equal_loads_resultant_at_midpoint() {
        let mut beam = SimpleBeam::new(0.0, 12.0).unwrap();
        beam.add_vertical(3.0, 25.0).unwrap();
        beam.add_vertical(9.0, 25.0).unwrap();
        let res = beam.vertical_load_resultant().unwrap();
        assert!((res.position - 6.0).abs() < EPS, "x_bar {}", res.position);
    }

    #[test]
    fn resultant_reproduces_the_same_reactions() {
        // The defining property: a beam carrying just the resultant load
        // has identical support reactions to the full load set.
        let mut beam = SimpleBeam::new(1.0, 9.0).unwrap();
        beam.add_vertical(2.0, 30.0).unwrap();
        beam.add_vertical(5.0, 50.0).unwrap();
        beam.add_vertical(8.0, 20.0).unwrap();
        let full = beam.reactions();
        let res = beam.vertical_load_resultant().unwrap();

        let mut equiv = SimpleBeam::new(1.0, 9.0).unwrap();
        equiv.add_vertical(res.position, res.down).unwrap();
        let single = equiv.reactions();

        assert!(
            (single.pin_vertical - full.pin_vertical).abs() < EPS,
            "pin {} vs {}",
            single.pin_vertical,
            full.pin_vertical
        );
        assert!(
            (single.roller_vertical - full.roller_vertical).abs() < EPS,
            "roller {} vs {}",
            single.roller_vertical,
            full.roller_vertical
        );
    }

    #[test]
    fn zero_net_vertical_load_has_no_resultant() {
        // No loads at all.
        let empty = SimpleBeam::new(0.0, 5.0).unwrap();
        assert!(empty.vertical_load_resultant().is_none());
        // A balanced up/down pair (net zero) is a couple — no resultant.
        let mut couple = SimpleBeam::new(0.0, 5.0).unwrap();
        couple.add_vertical(1.0, 40.0).unwrap();
        couple.add_vertical(4.0, -40.0).unwrap();
        assert!(couple.vertical_load_resultant().is_none());
    }
}
