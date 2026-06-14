//! Palmgren-Miner linear cumulative-damage rule.
//!
//! ## Model
//!
//! When a part sees a spectrum of load blocks — `n_i` applied cycles at
//! stress level `i`, where that level alone would fail the part at `N_i`
//! cycles — the **Palmgren-Miner** rule sums the fractional damage of
//! each block:
//!
//! ```text
//! D = sum_i ( n_i / N_i )
//! ```
//!
//! Each ratio `n_i / N_i` is the fraction of the part's life consumed by
//! that block. Failure is predicted when the accumulated damage reaches
//! unity:
//!
//! ```text
//! D = 1  =>  failure
//! ```
//!
//! `D < 1` means remaining life; the remaining-life fraction is `1 - D`.
//! The rule is linear and load-order-independent: it ignores the
//! sequence in which the blocks are applied.
//!
//! A companion calculation finds the **repeats to failure** of a fixed
//! duty block: if one pass of the duty cycle accumulates damage
//! `D_block`, then `1 / D_block` passes reach `D = 1`.
//!
//! ## Honest scope
//!
//! This is the textbook linear damage rule. It deliberately omits
//! load-sequence effects, the empirical sub-/super-unity failure sums
//! often seen in test data, and nonlinear damage models (Marco-Starkey,
//! double-linear, Corten-Dolan). Research/educational grade, not a
//! production design tool.

use crate::error::{FatigueError, Result};
use serde::{Deserialize, Serialize};

/// One block of a load spectrum: `applied_cycles` cycles applied at a
/// level whose stand-alone life is `life_cycles`.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct DamageBlock {
    /// Number of cycles actually applied at this level, `n_i >= 0`.
    pub applied_cycles: f64,
    /// Cycles-to-failure at this level on its own, `N_i > 0`.
    pub life_cycles: f64,
}

impl DamageBlock {
    /// Build a validated block.
    ///
    /// # Errors
    ///
    /// Returns [`FatigueError::Invalid`] if `applied_cycles` is negative
    /// / non-finite, or if `life_cycles` is not strictly positive /
    /// finite.
    pub fn new(applied_cycles: f64, life_cycles: f64) -> Result<Self> {
        if !applied_cycles.is_finite() || applied_cycles < 0.0 {
            return Err(FatigueError::invalid(
                "applied_cycles",
                format!("must be finite and >= 0, got {applied_cycles}"),
            ));
        }
        if !life_cycles.is_finite() || life_cycles <= 0.0 {
            return Err(FatigueError::invalid(
                "life_cycles",
                format!("must be finite and > 0, got {life_cycles}"),
            ));
        }
        Ok(DamageBlock {
            applied_cycles,
            life_cycles,
        })
    }

    /// The fractional damage of this block alone: `n_i / N_i`.
    pub fn damage(&self) -> f64 {
        self.applied_cycles / self.life_cycles
    }
}

/// The accumulated Palmgren-Miner damage over a list of blocks:
/// `D = sum_i n_i / N_i`.
///
/// An empty list returns `0.0` (no damage). The blocks are pre-validated
/// [`DamageBlock`]s, so this cannot fail.
pub fn miner_damage(blocks: &[DamageBlock]) -> f64 {
    blocks.iter().map(DamageBlock::damage).sum()
}

/// Convenience wrapper that builds and sums blocks from two parallel
/// slices of applied cycles `n_i` and stand-alone lives `N_i`.
///
/// # Errors
///
/// Returns [`FatigueError::Dimension`] if the slices differ in length,
/// and propagates any [`FatigueError::Invalid`] from a bad value.
pub fn miner_damage_from_slices(applied: &[f64], lives: &[f64]) -> Result<f64> {
    if applied.len() != lives.len() {
        return Err(FatigueError::dimension(
            applied.len(),
            lives.len(),
            "damage blocks",
        ));
    }
    let mut total = 0.0;
    for (&n, &big_n) in applied.iter().zip(lives.iter()) {
        total += DamageBlock::new(n, big_n)?.damage();
    }
    Ok(total)
}

/// `true` when the accumulated damage reaches or exceeds unity (failure
/// predicted). Uses `>= 1.0` so that a damage of exactly `1.0` counts as
/// failure.
pub fn has_failed(blocks: &[DamageBlock]) -> bool {
    miner_damage(blocks) >= 1.0
}

/// The remaining-life fraction `1 - D`, clamped at zero (a fully damaged
/// or over-damaged spectrum reports `0.0`, never negative).
pub fn remaining_life_fraction(blocks: &[DamageBlock]) -> f64 {
    (1.0 - miner_damage(blocks)).max(0.0)
}

/// How many repeats of a fixed duty block reach failure: `1 / D_block`,
/// where `D_block` is the damage accumulated in one pass of `blocks`.
///
/// # Errors
///
/// Returns [`FatigueError::Domain`] if one pass accumulates zero damage
/// (an empty spectrum, or one with no applied cycles) — failure is never
/// reached, so the repeat count is unbounded.
pub fn repeats_to_failure(blocks: &[DamageBlock]) -> Result<f64> {
    let d = miner_damage(blocks);
    if d <= 0.0 {
        return Err(FatigueError::domain(
            "duty block accumulates zero damage; failure is never reached"
                .to_string(),
        ));
    }
    Ok(1.0 / d)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(x: f64, y: f64) {
        let tol = 1e-9 * x.abs().max(y.abs()).max(1.0);
        assert!((x - y).abs() < tol, "expected {x} ~= {y}");
    }

    #[test]
    fn block_constructor_validates() {
        assert!(DamageBlock::new(-1.0, 100.0).is_err());
        assert!(DamageBlock::new(10.0, 0.0).is_err());
        assert!(DamageBlock::new(10.0, -5.0).is_err());
        assert!(DamageBlock::new(f64::NAN, 100.0).is_err());
        assert!(DamageBlock::new(0.0, 100.0).is_ok()); // zero applied is fine
        assert!(DamageBlock::new(10.0, 100.0).is_ok());
    }

    /// A single block applying exactly its full life is D = 1.
    #[test]
    fn full_life_single_block_is_unit_damage() {
        let b = DamageBlock::new(1.0e6, 1.0e6).unwrap();
        close(miner_damage(&[b]), 1.0);
        assert!(has_failed(&[b]));
    }

    /// Two half-life blocks sum to D = 1 (failure).
    #[test]
    fn two_half_blocks_reach_failure() {
        let b1 = DamageBlock::new(5.0e5, 1.0e6).unwrap(); // 0.5
        let b2 = DamageBlock::new(2.5e5, 5.0e5).unwrap(); // 0.5
        close(miner_damage(&[b1, b2]), 1.0);
        assert!(has_failed(&[b1, b2]));
        close(remaining_life_fraction(&[b1, b2]), 0.0);
    }

    /// Known textbook value: 0.3 + 0.2 + 0.1 = 0.6 damage, 0.4 left.
    #[test]
    fn known_partial_damage_sum() {
        let blocks = [
            DamageBlock::new(300.0, 1000.0).unwrap(), // 0.3
            DamageBlock::new(200.0, 1000.0).unwrap(), // 0.2
            DamageBlock::new(100.0, 1000.0).unwrap(), // 0.1
        ];
        close(miner_damage(&blocks), 0.6);
        assert!(!has_failed(&blocks));
        close(remaining_life_fraction(&blocks), 0.4);
    }

    /// Damage just over unity counts as failure and clamps remaining
    /// life to zero (not negative).
    #[test]
    fn over_unity_is_failure_and_clamps() {
        let blocks = [
            DamageBlock::new(800.0, 1000.0).unwrap(), // 0.8
            DamageBlock::new(500.0, 1000.0).unwrap(), // 0.5
        ];
        let d = miner_damage(&blocks);
        assert!(d > 1.0, "expected over-unity, got {d}");
        assert!(has_failed(&blocks));
        close(remaining_life_fraction(&blocks), 0.0);
    }

    #[test]
    fn empty_spectrum_is_zero_damage() {
        close(miner_damage(&[]), 0.0);
        assert!(!has_failed(&[]));
        close(remaining_life_fraction(&[]), 1.0);
    }

    /// The slice helper agrees with the block API and checks lengths.
    #[test]
    fn slice_helper_matches_and_checks_dims() {
        let applied = [300.0, 200.0, 100.0];
        let lives = [1000.0, 1000.0, 1000.0];
        let d = miner_damage_from_slices(&applied, &lives).unwrap();
        close(d, 0.6);

        // Length mismatch is a dimension error.
        let err = miner_damage_from_slices(&[1.0, 2.0], &[1000.0]).unwrap_err();
        assert_eq!(err.code(), "fatigue.dimension");

        // Bad value still propagates.
        assert!(miner_damage_from_slices(&[1.0], &[0.0]).is_err());
    }

    /// Repeats-to-failure inverts the per-pass damage.
    #[test]
    fn repeats_to_failure_inverts_block_damage() {
        // One pass: 0.25 + 0.05 = 0.30 damage -> 1/0.30 passes to fail.
        let blocks = [
            DamageBlock::new(250.0, 1000.0).unwrap(), // 0.25
            DamageBlock::new(50.0, 1000.0).unwrap(),  // 0.05
        ];
        let r = repeats_to_failure(&blocks).unwrap();
        close(r, 1.0 / 0.30);
        // Sanity: r passes times per-pass damage equals exactly 1.0.
        close(r * miner_damage(&blocks), 1.0);
    }

    #[test]
    fn repeats_to_failure_rejects_zero_damage() {
        assert!(repeats_to_failure(&[]).is_err());
        let no_load = [DamageBlock::new(0.0, 1000.0).unwrap()];
        assert!(repeats_to_failure(&no_load).is_err());
    }
}
