//! Pulley arrangements and their ideal mechanical advantage.
//!
//! ## Model
//!
//! The *ideal mechanical advantage* (IMA) of a rope-and-pulley machine is
//! the number of rope segments that directly support the movable load:
//!
//! ```text
//! IMA = n_supporting_ropes
//! ```
//!
//! This single rule specialises to the three textbook arrangements:
//!
//! - A **fixed** pulley (axle attached to the support) only redirects the
//!   rope. One segment runs to the load, so `IMA = 1`. You pull down with
//!   the same force as the weight, but in a more convenient direction.
//! - A **movable** pulley (axle attached to the load) is held up by two
//!   rope segments — the fixed end and the hauling part — so `IMA = 2`.
//! - A **block and tackle** (a compound of `f` fixed and `m` movable
//!   sheaves reeved by one continuous rope) is supported by `n` segments,
//!   so `IMA = n`. For a system rove to advantage with the hauling part
//!   leaving a movable block, `n = 2 m + 1`; rove to disadvantage it is
//!   `n = 2 m`.
//!
//! All of these are *ideal* (friction-free, massless rope and sheaves)
//! values. The friction-aware *actual* mechanical advantage lives in
//! [`crate::mechanics`], where it is the ideal value scaled by an
//! efficiency `eta` in `(0, 1]`.

use crate::error::{PulleyError, Result};
use serde::{Deserialize, Serialize};

/// The qualitative kind of a single-rope pulley arrangement.
///
/// This is a descriptive tag carried alongside a [`PulleySystem`]; the
/// numeric mechanical advantage always comes from the supporting-rope
/// count, never from this tag alone.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash, Serialize, Deserialize)]
pub enum PulleyKind {
    /// A single pulley whose axle is fixed to the support. Redirects the
    /// rope only; ideal mechanical advantage `1`.
    Fixed,
    /// A single pulley whose axle rides on the load. Supported by two rope
    /// segments; ideal mechanical advantage `2`.
    Movable,
    /// A compound block-and-tackle reeving one continuous rope over
    /// several fixed and movable sheaves. Ideal mechanical advantage
    /// equals the number of rope segments supporting the movable block.
    BlockAndTackle,
}

impl PulleyKind {
    /// A short, stable, lower-case identifier for logging / serialization.
    pub fn id(self) -> &'static str {
        match self {
            PulleyKind::Fixed => "fixed",
            PulleyKind::Movable => "movable",
            PulleyKind::BlockAndTackle => "block_and_tackle",
        }
    }
}

/// A rope-and-pulley machine, reduced to the one quantity that fixes its
/// ideal behaviour: the number of rope segments supporting the load.
///
/// Construct one with [`PulleySystem::fixed`], [`PulleySystem::movable`],
/// or [`PulleySystem::block_and_tackle`]; each validates its inputs and
/// records the matching [`PulleyKind`].
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash, Serialize, Deserialize)]
pub struct PulleySystem {
    kind: PulleyKind,
    supporting_ropes: u32,
}

impl PulleySystem {
    /// A single fixed pulley: ideal mechanical advantage `1`, one
    /// supporting rope segment. Redirects the applied force without
    /// multiplying it.
    pub fn fixed() -> Self {
        PulleySystem {
            kind: PulleyKind::Fixed,
            supporting_ropes: 1,
        }
    }

    /// A single movable pulley: ideal mechanical advantage `2`, two
    /// supporting rope segments.
    pub fn movable() -> Self {
        PulleySystem {
            kind: PulleyKind::Movable,
            supporting_ropes: 2,
        }
    }

    /// A block-and-tackle (compound) system supported by
    /// `supporting_ropes` rope segments, giving ideal mechanical advantage
    /// equal to that count.
    ///
    /// # Errors
    ///
    /// Returns [`PulleyError::Invalid`] if `supporting_ropes` is `0` — a
    /// machine with no segment supporting the load cannot lift it and
    /// would give an undefined (divide-by-zero) mechanical advantage.
    pub fn block_and_tackle(supporting_ropes: u32) -> Result<Self> {
        if supporting_ropes == 0 {
            return Err(PulleyError::invalid(
                "supporting_ropes",
                "a block and tackle needs at least one supporting rope segment",
            ));
        }
        Ok(PulleySystem {
            kind: PulleyKind::BlockAndTackle,
            supporting_ropes,
        })
    }

    /// Build a block-and-tackle from the count of movable sheaves and
    /// whether the hauling part leaves a movable block (*rove to
    /// advantage*) or a fixed block (*rove to disadvantage*).
    ///
    /// With `movable_sheaves = m`:
    ///
    /// - rove to advantage gives `n = 2 m + 1` supporting segments;
    /// - rove to disadvantage gives `n = 2 m` supporting segments.
    ///
    /// # Errors
    ///
    /// Returns [`PulleyError::Invalid`] if the resulting supporting-rope
    /// count would be `0` (i.e. `movable_sheaves == 0` while rove to
    /// disadvantage), which describes no working machine.
    pub fn reeved(movable_sheaves: u32, rove_to_advantage: bool) -> Result<Self> {
        let n = if rove_to_advantage {
            2 * movable_sheaves + 1
        } else {
            2 * movable_sheaves
        };
        Self::block_and_tackle(n)
    }

    /// The [`PulleyKind`] tag for this system.
    pub fn kind(self) -> PulleyKind {
        self.kind
    }

    /// The number of rope segments supporting the load.
    pub fn supporting_ropes(self) -> u32 {
        self.supporting_ropes
    }

    /// The ideal (friction-free) mechanical advantage `IMA = n`, returned
    /// as an exact integer-valued [`f64`]. Equal to the supporting-rope
    /// count and to the velocity ratio.
    pub fn ideal_mechanical_advantage(self) -> f64 {
        f64::from(self.supporting_ropes)
    }

    /// The velocity ratio `VR = distance moved by effort / distance moved
    /// by load`. For an ideal rope machine the inextensible rope ties the
    /// two distances together so that `VR == IMA == n`.
    ///
    /// Pulling in `n` rope segments each by the load's travel `d` requires
    /// hauling `n d` of rope, so the effort end moves `n` times as far as
    /// the load.
    pub fn velocity_ratio(self) -> f64 {
        self.ideal_mechanical_advantage()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-12;

    /// Fixed pulley: MA == 1, one supporting rope, direction-change only.
    #[test]
    fn fixed_has_unit_advantage() {
        let p = PulleySystem::fixed();
        assert_eq!(p.kind(), PulleyKind::Fixed);
        assert_eq!(p.supporting_ropes(), 1);
        assert!((p.ideal_mechanical_advantage() - 1.0).abs() < EPS);
    }

    /// Movable pulley: MA == 2, two supporting ropes.
    #[test]
    fn movable_has_advantage_two() {
        let p = PulleySystem::movable();
        assert_eq!(p.kind(), PulleyKind::Movable);
        assert_eq!(p.supporting_ropes(), 2);
        assert!((p.ideal_mechanical_advantage() - 2.0).abs() < EPS);
    }

    /// Block-and-tackle: MA == n, the supporting-rope count, for a sweep
    /// of values.
    #[test]
    fn block_and_tackle_advantage_equals_n() {
        for n in 1..=12u32 {
            let p = PulleySystem::block_and_tackle(n).unwrap();
            assert_eq!(p.kind(), PulleyKind::BlockAndTackle);
            assert_eq!(p.supporting_ropes(), n);
            let expected = f64::from(n);
            assert!(
                (p.ideal_mechanical_advantage() - expected).abs() < EPS,
                "n = {n}"
            );
        }
    }

    /// Velocity ratio equals mechanical advantage for the ideal machine.
    #[test]
    fn velocity_ratio_equals_mechanical_advantage() {
        for n in 1..=10u32 {
            let p = PulleySystem::block_and_tackle(n).unwrap();
            assert!(
                (p.velocity_ratio() - p.ideal_mechanical_advantage()).abs() < EPS,
                "n = {n}"
            );
        }
    }

    /// A movable pulley is the `n = 2` special case of block-and-tackle.
    #[test]
    fn movable_matches_two_rope_tackle() {
        let movable = PulleySystem::movable();
        let tackle = PulleySystem::block_and_tackle(2).unwrap();
        assert!(
            (movable.ideal_mechanical_advantage() - tackle.ideal_mechanical_advantage()).abs()
                < EPS
        );
    }

    /// A fixed pulley is the `n = 1` special case of block-and-tackle.
    #[test]
    fn fixed_matches_one_rope_tackle() {
        let fixed = PulleySystem::fixed();
        let tackle = PulleySystem::block_and_tackle(1).unwrap();
        assert!(
            (fixed.ideal_mechanical_advantage() - tackle.ideal_mechanical_advantage()).abs() < EPS
        );
    }

    /// Reeving: rove to advantage with `m` movable sheaves gives
    /// `n = 2 m + 1`; rove to disadvantage gives `n = 2 m`.
    #[test]
    fn reeving_counts_segments() {
        // Rove to advantage.
        for m in 0..=6u32 {
            let p = PulleySystem::reeved(m, true).unwrap();
            assert_eq!(p.supporting_ropes(), 2 * m + 1);
        }
        // Rove to disadvantage (m >= 1, else zero segments).
        for m in 1..=6u32 {
            let p = PulleySystem::reeved(m, false).unwrap();
            assert_eq!(p.supporting_ropes(), 2 * m);
        }
    }

    /// One movable sheave rove to advantage is the classic `n = 3`
    /// "luff tackle"; two movable sheaves rove to advantage is the
    /// `n = 5` "double luff".
    #[test]
    fn reeving_named_tackles() {
        assert_eq!(PulleySystem::reeved(1, true).unwrap().supporting_ropes(), 3);
        assert_eq!(PulleySystem::reeved(2, true).unwrap().supporting_ropes(), 5);
    }

    /// A zero-rope block-and-tackle is rejected.
    #[test]
    fn zero_ropes_rejected() {
        let err = PulleySystem::block_and_tackle(0).unwrap_err();
        assert_eq!(err.code(), "pulley.invalid");
        // Rove to disadvantage with no movable sheaves => zero segments.
        let err = PulleySystem::reeved(0, false).unwrap_err();
        assert_eq!(err.code(), "pulley.invalid");
    }

    /// `PulleyKind::id` is stable and distinct.
    #[test]
    fn kind_ids_are_distinct() {
        assert_eq!(PulleyKind::Fixed.id(), "fixed");
        assert_eq!(PulleyKind::Movable.id(), "movable");
        assert_eq!(PulleyKind::BlockAndTackle.id(), "block_and_tackle");
    }
}
