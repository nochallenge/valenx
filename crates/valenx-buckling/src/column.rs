//! Euler column buckling: end conditions, critical load, critical
//! stress, and slenderness.
//!
//! A slender, axially loaded, perfectly straight elastic column fails
//! not by crushing but by *buckling* — it suddenly bows sideways once
//! the compressive load reaches the **Euler critical load**
//!
//! ```text
//! P_cr = pi^2 E I / (K L)^2
//! ```
//!
//! where `E` is Young's modulus, `I` the smallest second moment of area
//! of the cross-section, `L` the unsupported length, and `K` the
//! *effective-length factor* set by how the two ends are restrained.
//! The product `K L` is the **effective length** — the length of the
//! equivalent pinned-pinned column with the same critical load.
//!
//! From `P_cr` two further engineering quantities follow directly:
//!
//! - the **critical stress** `sigma_cr = P_cr / A`, the average axial
//!   stress in the cross-section at the buckling load, and
//! - the **slenderness ratio** `K L / r` with `r = sqrt(I / A)` the
//!   radius of gyration — the dimensionless number that decides whether
//!   a column is "long" (Euler-governed) or "short" (crushing-governed).

use serde::{Deserialize, Serialize};

use crate::error::BucklingError;

/// How the two ends of the column are restrained.
///
/// Each variant carries the textbook **theoretical effective-length
/// factor** `K` (AISC / Hibbeler). The effective length is `K L`, and
/// the critical load scales as `1 / K^2`, so smaller `K` means a stiffer
/// column that carries more load before buckling.
///
/// | End condition | `K`  | Relative `P_cr` |
/// |---------------|------|-----------------|
/// | Fixed-free    | 2.0  | 0.25x pinned    |
/// | Fixed-pinned  | 0.7  | ~2.04x pinned   |
/// | Pinned-pinned | 1.0  | 1.0x (baseline) |
/// | Fixed-fixed   | 0.5  | 4.0x pinned     |
#[derive(Copy, Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum EndCondition {
    /// Both ends pinned (hinged). The reference case, `K = 1`.
    PinnedPinned,
    /// One end fully fixed, the other completely free (a cantilever
    /// column / flagpole). The weakest case, `K = 2`.
    FixedFree,
    /// Both ends fully fixed (built-in / clamped). The strongest of the
    /// four, `K = 0.5`.
    FixedFixed,
    /// One end fixed, the other pinned. `K = 0.7`.
    FixedPinned,
}

impl EndCondition {
    /// Theoretical effective-length factor `K` for this end condition.
    ///
    /// ```
    /// use valenx_buckling::EndCondition;
    ///
    /// assert_eq!(EndCondition::PinnedPinned.factor_k(), 1.0);
    /// assert_eq!(EndCondition::FixedFree.factor_k(), 2.0);
    /// assert_eq!(EndCondition::FixedFixed.factor_k(), 0.5);
    /// assert_eq!(EndCondition::FixedPinned.factor_k(), 0.7);
    /// ```
    pub fn factor_k(self) -> f64 {
        match self {
            EndCondition::PinnedPinned => 1.0,
            EndCondition::FixedFree => 2.0,
            EndCondition::FixedFixed => 0.5,
            EndCondition::FixedPinned => 0.7,
        }
    }

    /// Short human-readable label for UI / reports.
    pub fn label(self) -> &'static str {
        match self {
            EndCondition::PinnedPinned => "pinned-pinned",
            EndCondition::FixedFree => "fixed-free",
            EndCondition::FixedFixed => "fixed-fixed",
            EndCondition::FixedPinned => "fixed-pinned",
        }
    }
}

/// A slender elastic column ready for Euler buckling analysis.
///
/// Units are the caller's responsibility but must be **consistent**. A
/// convenient SI set is: `youngs_modulus` in pascals (`N/m^2`),
/// `second_moment_area` in `m^4`, `length` in metres, and `area` in
/// `m^2`; then [`Column::critical_load`] is in newtons and
/// [`Column::critical_stress`] in pascals. The slenderness ratio is
/// dimensionless in any consistent system.
///
/// Every field is validated on construction via [`Column::new`], so an
/// existing `Column` is always physically well-formed (all quantities
/// strictly positive and finite).
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Column {
    /// Young's modulus `E` (elastic modulus) of the column material.
    pub youngs_modulus: f64,
    /// Smallest second moment of area `I` of the cross-section about the
    /// axis it buckles around (the weak axis).
    pub second_moment_area: f64,
    /// Unsupported length `L` of the column.
    pub length: f64,
    /// Cross-sectional area `A`.
    pub area: f64,
    /// How the two ends are restrained.
    pub end_condition: EndCondition,
}

impl Column {
    /// Construct a validated column.
    ///
    /// `youngs_modulus`, `second_moment_area`, `length`, and `area` must
    /// each be strictly positive and finite; otherwise a
    /// [`BucklingError::NonPositive`] naming the offending quantity is
    /// returned.
    ///
    /// ```
    /// use valenx_buckling::{Column, EndCondition};
    ///
    /// // A 3 m pinned-pinned steel strut.
    /// let col = Column::new(
    ///     200.0e9,   // E  [Pa]
    ///     4.909e-6,  // I  [m^4]
    ///     3.0,       // L  [m]
    ///     3.142e-3,  // A  [m^2]
    ///     EndCondition::PinnedPinned,
    /// )
    /// .unwrap();
    /// assert_eq!(col.length, 3.0);
    ///
    /// // A non-positive length is rejected.
    /// assert!(Column::new(200.0e9, 4.9e-6, 0.0, 3.1e-3, EndCondition::PinnedPinned).is_err());
    /// ```
    pub fn new(
        youngs_modulus: f64,
        second_moment_area: f64,
        length: f64,
        area: f64,
        end_condition: EndCondition,
    ) -> Result<Self, BucklingError> {
        Ok(Column {
            youngs_modulus: BucklingError::require_positive("E", youngs_modulus)?,
            second_moment_area: BucklingError::require_positive("I", second_moment_area)?,
            length: BucklingError::require_positive("L", length)?,
            area: BucklingError::require_positive("A", area)?,
            end_condition,
        })
    }

    /// Effective-length factor `K` of this column's end condition.
    pub fn factor_k(&self) -> f64 {
        self.end_condition.factor_k()
    }

    /// Effective length `K L` — the equivalent pinned-pinned span.
    ///
    /// ```
    /// use valenx_buckling::{Column, EndCondition};
    ///
    /// let col = Column::new(1.0, 1.0, 2.0, 1.0, EndCondition::FixedFixed).unwrap();
    /// // K = 0.5, so effective length is 0.5 * 2.0 = 1.0.
    /// assert!((col.effective_length() - 1.0).abs() < 1e-12);
    /// ```
    pub fn effective_length(&self) -> f64 {
        self.factor_k() * self.length
    }

    /// **Euler critical buckling load** `P_cr = pi^2 E I / (K L)^2`.
    ///
    /// The compressive axial load at which a perfectly straight elastic
    /// column buckles. It scales as `1 / L^2` (doubling the length
    /// quarters the load) and linearly with `I` and `E`.
    ///
    /// ```
    /// use valenx_buckling::{Column, EndCondition};
    ///
    /// // Unit column: P_cr = pi^2 for a pinned-pinned bar with
    /// // E = I = L = 1.
    /// let col = Column::new(1.0, 1.0, 1.0, 1.0, EndCondition::PinnedPinned).unwrap();
    /// let expected = std::f64::consts::PI.powi(2);
    /// assert!((col.critical_load() - expected).abs() < 1e-12);
    /// ```
    pub fn critical_load(&self) -> f64 {
        let kl = self.effective_length();
        std::f64::consts::PI.powi(2) * self.youngs_modulus * self.second_moment_area / (kl * kl)
    }

    /// **Critical (buckling) stress** `sigma_cr = P_cr / A`.
    ///
    /// The average axial stress over the cross-section at the Euler
    /// critical load. For a real column this is only meaningful while it
    /// stays below the material yield strength — above that the column
    /// would yield (crush / Johnson-transition) before reaching the
    /// elastic Euler load, which this textbook model does not capture
    /// (see the crate-level honest-scope note).
    ///
    /// ```
    /// use valenx_buckling::{Column, EndCondition};
    ///
    /// let col = Column::new(1.0, 1.0, 1.0, 2.0, EndCondition::PinnedPinned).unwrap();
    /// // P_cr = pi^2, A = 2  =>  sigma_cr = pi^2 / 2.
    /// let expected = std::f64::consts::PI.powi(2) / 2.0;
    /// assert!((col.critical_stress() - expected).abs() < 1e-12);
    /// ```
    pub fn critical_stress(&self) -> f64 {
        self.critical_load() / self.area
    }

    /// Radius of gyration `r = sqrt(I / A)` of the cross-section.
    ///
    /// ```
    /// use valenx_buckling::{Column, EndCondition};
    ///
    /// // I = 4, A = 1  =>  r = 2.
    /// let col = Column::new(1.0, 4.0, 1.0, 1.0, EndCondition::PinnedPinned).unwrap();
    /// assert!((col.radius_of_gyration() - 2.0).abs() < 1e-12);
    /// ```
    pub fn radius_of_gyration(&self) -> f64 {
        (self.second_moment_area / self.area).sqrt()
    }

    /// **Slenderness ratio** `K L / r`.
    ///
    /// The dimensionless ratio of effective length to radius of
    /// gyration. Larger means more slender (more buckling-prone); the
    /// Euler critical stress falls as `1 / (slenderness)^2`.
    ///
    /// ```
    /// use valenx_buckling::{Column, EndCondition};
    ///
    /// // K L = 1 (pinned, L = 1), r = 1 (I = A = 1)  =>  slenderness 1.
    /// let col = Column::new(1.0, 1.0, 1.0, 1.0, EndCondition::PinnedPinned).unwrap();
    /// assert!((col.slenderness_ratio() - 1.0).abs() < 1e-12);
    /// ```
    pub fn slenderness_ratio(&self) -> f64 {
        self.effective_length() / self.radius_of_gyration()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    const EPS: f64 = 1e-9;

    /// All four textbook effective-length factors.
    #[test]
    fn end_condition_factors_match_textbook() {
        assert!((EndCondition::PinnedPinned.factor_k() - 1.0).abs() < EPS);
        assert!((EndCondition::FixedFree.factor_k() - 2.0).abs() < EPS);
        assert!((EndCondition::FixedFixed.factor_k() - 0.5).abs() < EPS);
        assert!((EndCondition::FixedPinned.factor_k() - 0.7).abs() < EPS);
    }

    /// Unit pinned-pinned column: P_cr = pi^2 exactly.
    #[test]
    fn unit_column_critical_load_is_pi_squared() {
        let col = Column::new(1.0, 1.0, 1.0, 1.0, EndCondition::PinnedPinned).unwrap();
        assert!((col.critical_load() - PI.powi(2)).abs() < EPS);
    }

    /// P_cr scales as 1/L^2: doubling the length quarters the load.
    #[test]
    fn critical_load_scales_inverse_length_squared() {
        let short = Column::new(200.0e9, 5.0e-6, 2.0, 3.0e-3, EndCondition::PinnedPinned).unwrap();
        let long = Column::new(200.0e9, 5.0e-6, 4.0, 3.0e-3, EndCondition::PinnedPinned).unwrap();
        // Length x2 -> load x(1/4).
        let ratio = long.critical_load() / short.critical_load();
        assert!((ratio - 0.25).abs() < EPS, "expected 0.25, got {ratio}");
    }

    /// A longer column buckles at a lower load (monotonic in L).
    #[test]
    fn longer_column_buckles_at_lower_load() {
        let a = Column::new(70.0e9, 2.0e-6, 1.5, 2.0e-3, EndCondition::FixedPinned).unwrap();
        let b = Column::new(70.0e9, 2.0e-6, 3.0, 2.0e-3, EndCondition::FixedPinned).unwrap();
        assert!(
            b.critical_load() < a.critical_load(),
            "longer column should carry less: {} vs {}",
            b.critical_load(),
            a.critical_load()
        );
    }

    /// P_cr scales linearly with the second moment of area I.
    #[test]
    fn critical_load_scales_linearly_with_i() {
        let base = Column::new(200.0e9, 4.0e-6, 2.5, 3.0e-3, EndCondition::PinnedPinned).unwrap();
        let triple =
            Column::new(200.0e9, 12.0e-6, 2.5, 3.0e-3, EndCondition::PinnedPinned).unwrap();
        let ratio = triple.critical_load() / base.critical_load();
        assert!((ratio - 3.0).abs() < EPS, "expected 3.0, got {ratio}");
    }

    /// P_cr scales linearly with Young's modulus E.
    #[test]
    fn critical_load_scales_linearly_with_e() {
        let steel = Column::new(200.0e9, 4.0e-6, 2.5, 3.0e-3, EndCondition::PinnedPinned).unwrap();
        let alu = Column::new(100.0e9, 4.0e-6, 2.5, 3.0e-3, EndCondition::PinnedPinned).unwrap();
        // Half the modulus -> half the load.
        let ratio = alu.critical_load() / steel.critical_load();
        assert!((ratio - 0.5).abs() < EPS, "expected 0.5, got {ratio}");
    }

    /// Fixed-fixed carries 4x the load of an otherwise identical
    /// pinned-pinned column (K = 0.5 vs 1.0, load ~ 1/K^2).
    #[test]
    fn fixed_fixed_carries_four_times_pinned() {
        let pinned = Column::new(200.0e9, 4.0e-6, 3.0, 3.0e-3, EndCondition::PinnedPinned).unwrap();
        let fixed = Column::new(200.0e9, 4.0e-6, 3.0, 3.0e-3, EndCondition::FixedFixed).unwrap();
        let ratio = fixed.critical_load() / pinned.critical_load();
        assert!((ratio - 4.0).abs() < EPS, "expected 4.0, got {ratio}");
    }

    /// Fixed-free is the weakest: 1/4 of the pinned-pinned load
    /// (K = 2 vs 1.0).
    #[test]
    fn fixed_free_is_quarter_of_pinned() {
        let pinned = Column::new(200.0e9, 4.0e-6, 3.0, 3.0e-3, EndCondition::PinnedPinned).unwrap();
        let free = Column::new(200.0e9, 4.0e-6, 3.0, 3.0e-3, EndCondition::FixedFree).unwrap();
        let ratio = free.critical_load() / pinned.critical_load();
        assert!((ratio - 0.25).abs() < EPS, "expected 0.25, got {ratio}");
    }

    /// Ordering of all four end conditions by carrying capacity:
    /// fixed-free < pinned-pinned < fixed-pinned < fixed-fixed.
    #[test]
    fn end_condition_strength_ordering() {
        let mk = |ec| {
            Column::new(200.0e9, 4.0e-6, 3.0, 3.0e-3, ec)
                .unwrap()
                .critical_load()
        };
        let free = mk(EndCondition::FixedFree);
        let pinned = mk(EndCondition::PinnedPinned);
        let fixed_pinned = mk(EndCondition::FixedPinned);
        let fixed_fixed = mk(EndCondition::FixedFixed);
        assert!(free < pinned);
        assert!(pinned < fixed_pinned);
        assert!(fixed_pinned < fixed_fixed);
    }

    /// Critical stress is exactly P_cr / A.
    #[test]
    fn critical_stress_is_load_over_area() {
        let col = Column::new(200.0e9, 4.0e-6, 3.0, 3.0e-3, EndCondition::FixedPinned).unwrap();
        let expected = col.critical_load() / col.area;
        assert!((col.critical_stress() - expected).abs() < EPS * expected.abs());
    }

    /// Effective length equals K * L.
    #[test]
    fn effective_length_is_k_times_l() {
        let col = Column::new(1.0, 1.0, 4.0, 1.0, EndCondition::FixedPinned).unwrap();
        // K = 0.7 -> 2.8.
        assert!((col.effective_length() - 2.8).abs() < EPS);
    }

    /// Radius of gyration r = sqrt(I/A) and slenderness K L / r.
    #[test]
    fn radius_of_gyration_and_slenderness() {
        // I = 9, A = 1 -> r = 3.  Pinned, L = 6 -> K L = 6 -> slenderness 2.
        let col = Column::new(1.0, 9.0, 6.0, 1.0, EndCondition::PinnedPinned).unwrap();
        assert!((col.radius_of_gyration() - 3.0).abs() < EPS);
        assert!((col.slenderness_ratio() - 2.0).abs() < EPS);
    }

    /// Euler critical stress can be written sigma_cr = pi^2 E /
    /// (K L / r)^2; verify it matches P_cr/A independently.
    #[test]
    fn critical_stress_equals_pi2_e_over_slenderness_squared() {
        let col = Column::new(200.0e9, 4.0e-6, 3.0, 3.0e-3, EndCondition::FixedFixed).unwrap();
        let s = col.slenderness_ratio();
        let from_slenderness = PI.powi(2) * col.youngs_modulus / (s * s);
        let from_load = col.critical_stress();
        // Relative agreement.
        assert!(
            (from_slenderness - from_load).abs() < EPS * from_load.abs(),
            "{from_slenderness} vs {from_load}"
        );
    }

    /// A worked textbook case: pinned-pinned A36-like steel strut.
    /// E = 200 GPa, I = 4.909e-6 m^4 (50 mm dia round), L = 3 m.
    /// P_cr = pi^2 * 200e9 * 4.909e-6 / 3^2 ~= 1.0768 MN.
    #[test]
    fn worked_steel_strut_case() {
        let col =
            Column::new(200.0e9, 4.909e-6, 3.0, 3.142e-3, EndCondition::PinnedPinned).unwrap();
        let expected = PI.powi(2) * 200.0e9 * 4.909e-6 / 9.0;
        // ~1.0768e6 N; check to within 1 N.
        assert!((col.critical_load() - expected).abs() < 1.0);
        assert!((col.critical_load() - 1.0768e6).abs() < 2.0e3);
    }

    /// Constructor rejects every non-positive / non-finite quantity and
    /// names the right field.
    #[test]
    fn constructor_rejects_bad_inputs() {
        let ec = EndCondition::PinnedPinned;
        let cases: [(&str, [f64; 4]); 5] = [
            ("E", [0.0, 1.0, 1.0, 1.0]),
            ("I", [1.0, -1.0, 1.0, 1.0]),
            ("L", [1.0, 1.0, 0.0, 1.0]),
            ("A", [1.0, 1.0, 1.0, f64::NAN]),
            ("I", [1.0, f64::INFINITY, 1.0, 1.0]),
        ];
        for (name, [e, i, l, a]) in cases {
            let err = Column::new(e, i, l, a, ec).unwrap_err();
            match err {
                BucklingError::NonPositive { name: got, .. } => {
                    assert_eq!(got, name, "wrong field flagged for case {name}");
                }
            }
        }
    }

    /// Valid inputs construct successfully and round-trip through JSON.
    #[test]
    fn serde_round_trip() {
        let col = Column::new(200.0e9, 4.0e-6, 2.5, 3.0e-3, EndCondition::FixedFixed).unwrap();
        let json = serde_json::to_string(&col).unwrap();
        let back: Column = serde_json::from_str(&json).unwrap();
        assert_eq!(col, back);
        assert!((col.critical_load() - back.critical_load()).abs() < EPS);
    }
}
