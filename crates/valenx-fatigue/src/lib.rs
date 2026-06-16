//! # valenx-fatigue
//!
//! A stress-life (S-N) fatigue calculator: turn a material's fatigue
//! data and a load history into the engineering answers — predicted life
//! in cycles, the alternating stress a non-zero mean stress still
//! permits, the factor of safety against a constant-life line, and the
//! accumulated damage of a variable-amplitude spectrum.
//!
//! ## What this is
//!
//! Three classical building blocks of high-cycle fatigue analysis:
//!
//! - **Basquin S-N curve** ([`sn`]). [`SnCurve`] holds the power-law
//!   `S = a N^b`. [`SnCurve::stress_at_cycles`] reads the stress at a
//!   life, [`SnCurve::cycles_to_failure`] reads the life at a stress, and
//!   an optional horizontal **endurance limit** ([`SnCurve::with_endurance_limit`])
//!   caps the curve so low stresses return [`Life::Infinite`].
//!   [`SnCurve::from_two_points`] fits the curve to two measured points.
//! - **Mean-stress correction** ([`mean_stress`]). A [`Material`] plus a
//!   [`MeanStressCriterion`] (the Goodman or Soderberg straight line, or
//!   the Gerber parabola) gives the allowable
//!   alternating stress at a given mean stress
//!   ([`Material::allowable_alternating`]) and the factor of safety for
//!   an operating point ([`Material::factor_of_safety`]).
//! - **Cumulative damage** ([`damage`]). [`miner_damage`] sums the
//!   Palmgren-Miner fractions `sum n_i/N_i` over a list of
//!   [`DamageBlock`]s; [`has_failed`] tests for `D >= 1`,
//!   [`remaining_life_fraction`] reports `1 - D`, and
//!   [`repeats_to_failure`] inverts a duty block's per-pass damage.
//!
//! ```
//! use valenx_fatigue::{Life, SnCurve};
//!
//! // A steel: fit the line through (1e3, 0.9*Su) and the fatigue
//! // limit (1e6, Se), then cap it at the endurance limit.
//! let su = 1000.0;
//! let se = 0.5 * su;
//! let curve = SnCurve::from_two_points(1.0e3, 0.9 * su, 1.0e6, se)
//!     .unwrap()
//!     .with_endurance_limit(se)
//!     .unwrap();
//!
//! // Life at a 600-unit stress amplitude is finite; below the limit
//! // it is infinite.
//! assert!(matches!(curve.cycles_to_failure(600.0).unwrap(), Life::Finite(_)));
//! assert_eq!(curve.cycles_to_failure(400.0).unwrap(), Life::Infinite);
//! ```
//!
//! ## Model
//!
//! All three relations are the standard textbook closed forms; each
//! module documents its own equations:
//!
//! - **Basquin** ([`sn`]): `S = a N^b` with `a > 0`, `b < 0`; inverted
//!   `N = (S/a)^(1/b)`; an optional flat endurance cutoff.
//! - **Goodman / Soderberg** ([`mean_stress`]): the straight constant-life
//!   line `sa/Se + sm/S0 = 1/n`, with the static intercept `S0` being the
//!   ultimate strength `Su` (Goodman) or the yield strength `Sy`
//!   (Soderberg). At zero mean stress and unit design factor the
//!   allowable alternating stress is exactly `Se`; at `sm = S0` it is `0`;
//!   a higher mean stress always lowers the allowable alternating stress.
//! - **Palmgren-Miner** ([`damage`]): linear, load-order-independent
//!   damage `D = sum n_i/N_i`, with failure at `D = 1`.
//!
//! Stresses are dimensionless in the caller's own consistent unit
//! (MPa, ksi, …) — every formula is homogeneous in stress, so the crate
//! never assumes a unit system. Float results are compared with absolute
//! tolerances in the test suite, never with exact equality.
//!
//! ## Honest scope
//!
//! Research/educational grade. These are the well-established
//! closed-form stress-life relations taught in a first machine-design or
//! mechanics-of-materials course — useful for learning, estimation, and
//! cross-checking, but **not** a clinical, medical, or production
//! engineering design tool. In particular the crate does **not**:
//!
//! - model strain-life / low-cycle fatigue (Coffin-Manson), the Morrow
//!   or Smith-Watson-Topper mean-stress corrections, or the first-cycle
//!   yield line `sa + sm <= Sy` (the Goodman, Soderberg, and Gerber
//!   constant-life criteria *are* provided);
//! - apply Marin surface / size / loading / temperature / reliability
//!   knock-down factors — the caller supplies an already-corrected curve
//!   and endurance limit;
//! - perform rain-flow cycle counting from a raw load-time history (the
//!   caller supplies the per-block cycle counts and lives); or
//! - account for load-sequence effects, crack-growth (fracture
//!   mechanics / Paris law), or notch / multiaxial stress states.
//!
//! Each omission is a documented, well-understood extension; nothing the
//! crate *does* compute is approximate beyond the stated linear models.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod damage;
pub mod error;
pub mod mean_stress;
pub mod sn;

// --- Convenience re-exports of the most-used types --------------------

pub use error::{ErrorCategory, FatigueError, Result};

pub use sn::{Life, SnCurve};

pub use mean_stress::{Material, MeanStressCriterion};

pub use damage::{
    has_failed, miner_damage, miner_damage_from_slices, remaining_life_fraction,
    repeats_to_failure, DamageBlock,
};

#[cfg(test)]
mod tests {
    use super::*;

    fn close(x: f64, y: f64) {
        let tol = 1e-9 * x.abs().max(y.abs()).max(1.0);
        assert!((x - y).abs() < tol, "expected {x} ~= {y}");
    }

    /// End-to-end: build an S-N curve, read the life of two stress
    /// blocks, and accumulate their Miner damage to a failure verdict.
    #[test]
    fn sn_to_miner_end_to_end() {
        let curve = SnCurve::new(1500.0, -0.1).unwrap();

        // Two operating stress amplitudes and their finite lives.
        let s_hi = 800.0;
        let s_lo = 500.0;
        let n_hi = curve.cycles_to_failure(s_hi).unwrap().finite().unwrap();
        let n_lo = curve.cycles_to_failure(s_lo).unwrap().finite().unwrap();
        assert!(n_lo > n_hi, "lower stress should give longer life");

        // Apply each block to exactly half its life -> D = 1 -> failure.
        let blocks = [
            DamageBlock::new(n_hi / 2.0, n_hi).unwrap(),
            DamageBlock::new(n_lo / 2.0, n_lo).unwrap(),
        ];
        close(miner_damage(&blocks), 1.0);
        assert!(has_failed(&blocks));
    }

    /// End-to-end: a Goodman correction feeds an S-N life lookup. The
    /// mean-stress-corrected equivalent alternating stress has a shorter
    /// life than the raw alternating stress.
    #[test]
    fn mean_stress_then_life() {
        let curve = SnCurve::new(1500.0, -0.1).unwrap();
        let mat = Material::new(250.0, 600.0, 900.0).unwrap();

        // An operating point with tensile mean stress.
        let sa = 300.0;
        let sm = 300.0;

        // Goodman equivalent fully-reversed stress:
        // sa_eq = sa / (1 - sm/Su).
        let sa_eq = sa / (1.0 - sm / mat.ultimate_strength);
        assert!(sa_eq > sa, "mean stress should raise the equivalent stress");

        let n_raw = curve.cycles_to_failure(sa).unwrap().finite().unwrap();
        let n_eq = curve.cycles_to_failure(sa_eq).unwrap().finite().unwrap();
        assert!(
            n_eq < n_raw,
            "mean-stress-corrected life {n_eq} should be below raw {n_raw}"
        );

        // The Goodman factor of safety satisfies its defining relation
        // exactly: 1/n = sa/Se + sm/Su, so n*(sa/Se + sm/Su) == 1. (This
        // is a distinct scalar from the equivalent-stress ratio Se/sa_eq,
        // which only coincides when Se == Su.)
        let n_fos = mat
            .factor_of_safety(MeanStressCriterion::Goodman, sa, sm)
            .unwrap();
        let demand = sa / mat.endurance_limit + sm / mat.ultimate_strength;
        close(n_fos * demand, 1.0);
        // And scaling the operating point by n_fos lands exactly on the
        // Goodman failure line (allowable alternating at design factor 1).
        let allow_on_line = mat
            .allowable_alternating(MeanStressCriterion::Goodman, sm * n_fos, 1.0)
            .unwrap();
        close(allow_on_line, sa * n_fos);
    }

    /// The public re-exports round-trip through serde_json.
    #[test]
    fn public_types_serialize() {
        let curve = SnCurve::new(1200.0, -0.085)
            .unwrap()
            .with_endurance_limit(200.0)
            .unwrap();
        let json = serde_json::to_string(&curve).unwrap();
        let back: SnCurve = serde_json::from_str(&json).unwrap();
        assert_eq!(curve, back);

        let mat = Material::new(200.0, 350.0, 500.0).unwrap();
        let json = serde_json::to_string(&mat).unwrap();
        let back: Material = serde_json::from_str(&json).unwrap();
        assert_eq!(mat, back);

        let block = DamageBlock::new(100.0, 1000.0).unwrap();
        let json = serde_json::to_string(&block).unwrap();
        let back: DamageBlock = serde_json::from_str(&json).unwrap();
        assert_eq!(block, back);
    }
}
