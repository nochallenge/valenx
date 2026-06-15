//! Ideal rigid levers: classification, mechanical advantage, and the
//! static moment-balance law.
//!
//! A lever is a rigid beam that pivots about a fixed point (the
//! *fulcrum*). An **effort** force is applied at distance `effort_arm`
//! from the fulcrum and balances a **load** force at distance
//! `load_arm`. Under the ideal assumptions of this crate (massless
//! beam, frictionless fulcrum, point forces perpendicular to their
//! arms) the system is in static equilibrium when the moments about the
//! fulcrum cancel:
//!
//! `effort * effort_arm = load * load_arm`
//!
//! The dimensionless **mechanical advantage** is the ratio of the arms,
//! equal in magnitude to the load-to-effort force ratio at balance:
//!
//! `MA = effort_arm / load_arm = load / effort`

use serde::{Deserialize, Serialize};

use crate::error::{validate_arm, validate_force, LeverError};

/// The three classes of lever, distinguished by the relative ordering
/// of the fulcrum, the effort, and the load along the beam.
///
/// The class is a *consequence* of the geometry; this crate also infers
/// it from the arm ratio via [`Lever::class`], so the enum is primarily
/// a labelling aid for callers describing a known arrangement.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub enum LeverClass {
    /// Fulcrum between effort and load (e.g. a seesaw, a pair of
    /// scissors, a crowbar pivoting at its bend). Mechanical advantage
    /// may be greater than, equal to, or less than one depending on
    /// where the fulcrum sits.
    First,
    /// Load between fulcrum and effort (e.g. a wheelbarrow, a nutcracker,
    /// a bottle opener). The effort arm always exceeds the load arm, so
    /// the mechanical advantage is strictly greater than one — these
    /// levers always multiply force.
    Second,
    /// Effort between fulcrum and load (e.g. tweezers, a fishing rod, the
    /// human forearm). The load arm always exceeds the effort arm, so the
    /// mechanical advantage is strictly less than one — these levers trade
    /// force for speed and range of motion.
    Third,
}

impl LeverClass {
    /// A short, stable, lowercase identifier (`"first"`, `"second"`,
    /// `"third"`).
    pub fn label(self) -> &'static str {
        match self {
            LeverClass::First => "first",
            LeverClass::Second => "second",
            LeverClass::Third => "third",
        }
    }
}

/// An ideal rigid lever defined by its two arm lengths.
///
/// Both arms are measured from the fulcrum to the line of action of the
/// corresponding force and are stored as finite, strictly positive
/// lengths in caller-chosen but consistent units (the ratio is
/// dimensionless, so any unit works as long as both arms share it).
/// Construct one with [`Lever::new`], which validates the inputs.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Lever {
    /// Distance from the fulcrum to the applied effort force (> 0).
    pub effort_arm: f64,
    /// Distance from the fulcrum to the resisting load force (> 0).
    pub load_arm: f64,
}

/// The complete static state of a balanced lever.
///
/// Returned by [`Lever::balance_load`] / [`Lever::balance_effort`]; the
/// invariant `effort * effort_arm == load * load_arm` holds to floating
/// point precision.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Balance {
    /// Effort force applied at [`Lever::effort_arm`].
    pub effort: f64,
    /// Load force resisted at [`Lever::load_arm`].
    pub load: f64,
    /// Moment shared by both sides (`effort * effort_arm`), i.e. the
    /// torque about the fulcrum.
    pub moment: f64,
}

impl Lever {
    /// Construct a lever from its effort and load arm lengths.
    ///
    /// # Errors
    ///
    /// Returns [`LeverError::NonPositiveArm`] if either arm is not a
    /// finite, strictly positive number.
    pub fn new(effort_arm: f64, load_arm: f64) -> Result<Self, LeverError> {
        let effort_arm = validate_arm("effort_arm", effort_arm)?;
        let load_arm = validate_arm("load_arm", load_arm)?;
        Ok(Self {
            effort_arm,
            load_arm,
        })
    }

    /// The dimensionless ideal **mechanical advantage**,
    /// `MA = effort_arm / load_arm`.
    ///
    /// `MA > 1` means the lever multiplies the applied effort (a small
    /// effort balances a large load); `MA < 1` means it divides it
    /// (trading force for distance/speed); `MA == 1` is a balance with no
    /// force change.
    pub fn mechanical_advantage(&self) -> f64 {
        self.effort_arm / self.load_arm
    }

    /// Infer the lever [`LeverClass`] implied by the arm ratio.
    ///
    /// `effort_arm > load_arm` (MA > 1) is reported as
    /// [`LeverClass::Second`]; `effort_arm < load_arm` (MA < 1) as
    /// [`LeverClass::Third`]; and an exactly equal pair (MA == 1) as
    /// [`LeverClass::First`] — the only class that admits balanced arms.
    /// Real first-class levers may of course have unequal arms; pass
    /// [`LeverClass::First`] explicitly when you know the arrangement.
    pub fn class(&self) -> LeverClass {
        if self.effort_arm > self.load_arm {
            LeverClass::Second
        } else if self.effort_arm < self.load_arm {
            LeverClass::Third
        } else {
            LeverClass::First
        }
    }

    /// Solve `effort * effort_arm = load * load_arm` for the **load**
    /// that the given `effort` balances.
    ///
    /// `load = effort * MA`. With `MA > 1` the balanced load exceeds the
    /// effort (force multiplied); with `MA < 1` it is smaller.
    ///
    /// # Errors
    ///
    /// Returns [`LeverError::NonFiniteForce`] if `effort` is not finite.
    pub fn balance_load(&self, effort: f64) -> Result<Balance, LeverError> {
        let effort = validate_force("effort", effort)?;
        let load = effort * self.mechanical_advantage();
        Ok(Balance {
            effort,
            load,
            moment: effort * self.effort_arm,
        })
    }

    /// Solve `effort * effort_arm = load * load_arm` for the **effort**
    /// required to balance the given `load`.
    ///
    /// `effort = load / MA`. With `MA > 1` the required effort is less
    /// than the load (the point of a force-multiplying lever); with
    /// `MA < 1` it is greater.
    ///
    /// # Errors
    ///
    /// Returns [`LeverError::NonFiniteForce`] if `load` is not finite.
    pub fn balance_effort(&self, load: f64) -> Result<Balance, LeverError> {
        let load = validate_force("load", load)?;
        let effort = load / self.mechanical_advantage();
        Ok(Balance {
            effort,
            load,
            moment: load * self.load_arm,
        })
    }

    /// Net moment (torque) about the fulcrum for an arbitrary
    /// effort/load pair, taking effort as the positive sense:
    ///
    /// `effort * effort_arm - load * load_arm`
    ///
    /// Zero means the lever is balanced; a positive value means the
    /// effort moment dominates (the load side rises), a negative value
    /// the reverse.
    ///
    /// # Errors
    ///
    /// Returns [`LeverError::NonFiniteForce`] if either force is not
    /// finite.
    pub fn net_moment(&self, effort: f64, load: f64) -> Result<f64, LeverError> {
        let effort = validate_force("effort", effort)?;
        let load = validate_force("load", load)?;
        Ok(effort * self.effort_arm - load * self.load_arm)
    }

    /// Whether an effort/load pair leaves the lever in static balance,
    /// within an absolute moment tolerance `tol`.
    ///
    /// Convenience wrapper over [`Lever::net_moment`]: returns `true`
    /// when its magnitude is `<= tol`.
    ///
    /// # Errors
    ///
    /// Returns [`LeverError::NonFiniteForce`] if either force is not
    /// finite.
    pub fn is_balanced(&self, effort: f64, load: f64, tol: f64) -> Result<bool, LeverError> {
        Ok(self.net_moment(effort, load)?.abs() <= tol)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    #[test]
    fn ma_is_effort_arm_over_load_arm() {
        // ground truth: 30 / 12 = 2.5
        let lever = Lever::new(30.0, 12.0).unwrap();
        assert!((lever.mechanical_advantage() - 2.5).abs() < EPS);
    }

    #[test]
    fn ma_equals_unity_for_equal_arms() {
        let lever = Lever::new(7.0, 7.0).unwrap();
        assert!((lever.mechanical_advantage() - 1.0).abs() < EPS);
        assert_eq!(lever.class(), LeverClass::First);
    }

    #[test]
    fn class_two_has_ma_greater_than_one() {
        // wheelbarrow-like: long effort arm, short load arm
        let lever = Lever::new(100.0, 25.0).unwrap();
        assert_eq!(lever.class(), LeverClass::Second);
        assert!(lever.mechanical_advantage() > 1.0);
    }

    #[test]
    fn class_three_has_ma_less_than_one() {
        // tweezers-like: short effort arm, long load arm
        let lever = Lever::new(20.0, 80.0).unwrap();
        assert_eq!(lever.class(), LeverClass::Third);
        assert!(lever.mechanical_advantage() < 1.0);
    }

    #[test]
    fn longer_effort_arm_gives_more_mechanical_advantage() {
        // Hold the load arm fixed; growing the effort arm must strictly
        // increase MA (monotonicity of the ratio in its numerator).
        let load_arm = 10.0;
        let short = Lever::new(15.0, load_arm).unwrap();
        let long = Lever::new(45.0, load_arm).unwrap();
        assert!(long.mechanical_advantage() > short.mechanical_advantage());
        // closed form: 45/10 = 4.5 vs 15/10 = 1.5
        assert!((short.mechanical_advantage() - 1.5).abs() < EPS);
        assert!((long.mechanical_advantage() - 4.5).abs() < EPS);
    }

    #[test]
    fn balance_equation_holds_for_solved_load() {
        // effort 40 N at arm 2.0 -> moment 80; load arm 0.5 -> load 160 N
        let lever = Lever::new(2.0, 0.5).unwrap();
        let b = lever.balance_load(40.0).unwrap();
        assert!((b.load - 160.0).abs() < EPS);
        // effort * effort_arm == load * load_arm  (the balance law)
        assert!((b.effort * lever.effort_arm - b.load * lever.load_arm).abs() < EPS);
        assert!((b.moment - 80.0).abs() < EPS);
    }

    #[test]
    fn balance_equation_holds_for_solved_effort() {
        // load 200 N at arm 0.5 -> moment 100; effort arm 2.0 -> effort 50 N
        let lever = Lever::new(2.0, 0.5).unwrap();
        let b = lever.balance_effort(200.0).unwrap();
        assert!((b.effort - 50.0).abs() < EPS);
        assert!((b.effort * lever.effort_arm - b.load * lever.load_arm).abs() < EPS);
        assert!((b.moment - 100.0).abs() < EPS);
    }

    #[test]
    fn solve_load_then_effort_round_trips() {
        let lever = Lever::new(3.3, 1.1).unwrap();
        let forward = lever.balance_load(17.0).unwrap();
        let back = lever.balance_effort(forward.load).unwrap();
        assert!((back.effort - 17.0).abs() < EPS);
    }

    #[test]
    fn output_force_multiplied_by_ma() {
        // Mechanical advantage of 4 quadruples the input effort into load.
        let lever = Lever::new(8.0, 2.0).unwrap();
        let b = lever.balance_load(10.0).unwrap();
        assert!((b.load - 40.0).abs() < EPS);
        assert!((b.load / b.effort - lever.mechanical_advantage()).abs() < EPS);
    }

    #[test]
    fn class_three_requires_more_effort_than_load() {
        // MA < 1 means the balancing effort exceeds the load it holds.
        let lever = Lever::new(5.0, 20.0).unwrap();
        let b = lever.balance_effort(10.0).unwrap();
        assert!(b.effort > b.load);
        assert!((b.effort - 40.0).abs() < EPS); // 10 / 0.25
    }

    #[test]
    fn net_moment_is_zero_at_balance() {
        let lever = Lever::new(2.0, 0.5).unwrap();
        let b = lever.balance_load(40.0).unwrap();
        let net = lever.net_moment(b.effort, b.load).unwrap();
        assert!(net.abs() < EPS);
        assert!(lever.is_balanced(b.effort, b.load, EPS).unwrap());
    }

    #[test]
    fn net_moment_sign_tracks_dominant_side() {
        let lever = Lever::new(2.0, 1.0).unwrap();
        // Effort moment (10*2=20) beats load moment (5*1=5) -> positive.
        assert!(lever.net_moment(10.0, 5.0).unwrap() > 0.0);
        // Load moment (50*1=50) beats effort moment (10*2=20) -> negative.
        assert!(lever.net_moment(10.0, 50.0).unwrap() < 0.0);
    }

    #[test]
    fn zero_effort_balances_only_zero_load() {
        let lever = Lever::new(4.0, 2.0).unwrap();
        let b = lever.balance_load(0.0).unwrap();
        assert!((b.load - 0.0).abs() < EPS);
        assert!((b.moment - 0.0).abs() < EPS);
    }

    #[test]
    fn new_rejects_non_positive_arms() {
        assert!(matches!(
            Lever::new(0.0, 1.0),
            Err(LeverError::NonPositiveArm {
                name: "effort_arm",
                ..
            })
        ));
        assert!(matches!(
            Lever::new(1.0, -2.0),
            Err(LeverError::NonPositiveArm {
                name: "load_arm",
                ..
            })
        ));
    }

    #[test]
    fn balance_rejects_non_finite_force() {
        let lever = Lever::new(2.0, 1.0).unwrap();
        assert!(matches!(
            lever.balance_load(f64::NAN),
            Err(LeverError::NonFiniteForce { name: "effort", .. })
        ));
        assert!(matches!(
            lever.balance_effort(f64::INFINITY),
            Err(LeverError::NonFiniteForce { name: "load", .. })
        ));
    }

    #[test]
    fn class_labels_are_stable() {
        assert_eq!(LeverClass::First.label(), "first");
        assert_eq!(LeverClass::Second.label(), "second");
        assert_eq!(LeverClass::Third.label(), "third");
    }
}
