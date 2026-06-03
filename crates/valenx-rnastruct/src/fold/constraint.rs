//! Folding constraints: hard structural constraints and soft
//! SHAPE-reactivity pseudo-energies.
//!
//! Both the [`crate::fold::zuker`] minimum-free-energy folder and the
//! [`crate::ensemble`] partition function consult a [`FoldConstraints`]
//! value before allowing a pair `(i, j)` or leaving a position
//! unpaired. Two kinds of constraint are supported:
//!
//! - **Hard constraints** — force a position paired / unpaired, or
//!   force / forbid a specific pair. A configuration violating any
//!   hard constraint is given energy `+∞` so the DP never selects it.
//! - **Soft constraints** — a per-position pseudo-energy added to the
//!   score whenever that base is *unpaired*. This is exactly the
//!   mechanism SHAPE-directed folding uses (the Deigan model): a
//!   chemical reactivity is converted to a free-energy bonus/penalty
//!   that nudges reactive (likely unpaired) bases out of helices.
//!
//! [`FoldConstraints::from_shape`] builds the soft term from a vector
//! of SHAPE reactivities with the Deigan *et al.* (2009) log
//! transform.

use crate::error::{Result, RnaStructError};

/// The Deigan-model slope `m` (kcal/mol) — the default from
/// Deigan *et al.* 2009 and ViennaRNA's `--shapeMethod=D`.
pub const DEIGAN_SLOPE: f64 = 1.8;

/// The Deigan-model intercept `b` (kcal/mol) — the default from
/// Deigan *et al.* 2009.
pub const DEIGAN_INTERCEPT: f64 = -0.6;

/// A per-position hard constraint on a single base.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum BaseConstraint {
    /// No constraint — the folder is free to pair or not pair this
    /// base.
    Free,
    /// This base must be paired (to *some* partner).
    Paired,
    /// This base must remain unpaired.
    Unpaired,
}

/// The full constraint set handed to a folder.
///
/// `n` is the sequence length the constraints were built for; the
/// folder checks that against the sequence it is given.
#[derive(Clone, Debug)]
pub struct FoldConstraints {
    n: usize,
    /// Per-base hard constraint.
    base: Vec<BaseConstraint>,
    /// Forced pairs — `(i, j)` that *must* appear.
    forced: Vec<(usize, usize)>,
    /// Forbidden pairs — `(i, j)` that must *not* appear.
    forbidden: Vec<(usize, usize)>,
    /// Per-base soft pseudo-energy added when the base is unpaired.
    soft_unpaired: Vec<f64>,
}

impl FoldConstraints {
    /// An unconstrained set for a length-`n` sequence.
    pub fn none(n: usize) -> Self {
        FoldConstraints {
            n,
            base: vec![BaseConstraint::Free; n],
            forced: Vec::new(),
            forbidden: Vec::new(),
            soft_unpaired: vec![0.0; n],
        }
    }

    /// The length these constraints were built for.
    pub fn len(&self) -> usize {
        self.n
    }

    /// `true` if built for a zero-length sequence.
    pub fn is_empty(&self) -> bool {
        self.n == 0
    }

    /// `true` if no hard or soft constraint is set (a pure
    /// unconstrained fold).
    pub fn is_unconstrained(&self) -> bool {
        self.forced.is_empty()
            && self.forbidden.is_empty()
            && self.base.iter().all(|c| *c == BaseConstraint::Free)
            && self.soft_unpaired.iter().all(|&e| e == 0.0)
    }

    /// Forces base `i` to be paired.
    ///
    /// # Errors
    /// [`RnaStructError::Invalid`] if `i >= n`.
    pub fn force_paired(&mut self, i: usize) -> Result<()> {
        self.set_base(i, BaseConstraint::Paired)
    }

    /// Forces base `i` to remain unpaired.
    ///
    /// # Errors
    /// [`RnaStructError::Invalid`] if `i >= n`.
    pub fn force_unpaired(&mut self, i: usize) -> Result<()> {
        self.set_base(i, BaseConstraint::Unpaired)
    }

    fn set_base(&mut self, i: usize, c: BaseConstraint) -> Result<()> {
        if i >= self.n {
            return Err(RnaStructError::invalid(
                "position",
                format!("base index {i} out of range for length {}", self.n),
            ));
        }
        self.base[i] = c;
        Ok(())
    }

    /// Forces the pair `(i, j)` to appear in the folded structure.
    ///
    /// # Errors
    /// [`RnaStructError::Invalid`] if an index is out of range or
    /// `i == j`.
    pub fn force_pair(&mut self, i: usize, j: usize) -> Result<()> {
        let (a, b) = self.check_pair(i, j)?;
        self.forced.push((a, b));
        self.base[a] = BaseConstraint::Paired;
        self.base[b] = BaseConstraint::Paired;
        Ok(())
    }

    /// Forbids the pair `(i, j)`.
    ///
    /// # Errors
    /// [`RnaStructError::Invalid`] if an index is out of range or
    /// `i == j`.
    pub fn forbid_pair(&mut self, i: usize, j: usize) -> Result<()> {
        let (a, b) = self.check_pair(i, j)?;
        self.forbidden.push((a, b));
        Ok(())
    }

    fn check_pair(&self, i: usize, j: usize) -> Result<(usize, usize)> {
        if i >= self.n || j >= self.n {
            return Err(RnaStructError::invalid(
                "pair",
                format!("pair ({i}, {j}) out of range for length {}", self.n),
            ));
        }
        if i == j {
            return Err(RnaStructError::invalid(
                "pair",
                "cannot constrain a base to pair with itself",
            ));
        }
        Ok(if i < j { (i, j) } else { (j, i) })
    }

    /// `true` if a pair `(i, j)` (with `i < j`) is allowed by the hard
    /// constraints — neither base is forced-unpaired and the pair is
    /// not forbidden.
    pub fn pair_allowed(&self, i: usize, j: usize) -> bool {
        if i >= self.n || j >= self.n {
            return false;
        }
        if self.base[i] == BaseConstraint::Unpaired
            || self.base[j] == BaseConstraint::Unpaired
        {
            return false;
        }
        let key = if i < j { (i, j) } else { (j, i) };
        !self.forbidden.contains(&key)
    }

    /// `true` if base `i` is allowed to be unpaired by the hard
    /// constraints (not forced-paired).
    pub fn unpaired_allowed(&self, i: usize) -> bool {
        i < self.n && self.base[i] != BaseConstraint::Paired
    }

    /// The list of forced pairs.
    pub fn forced_pairs(&self) -> &[(usize, usize)] {
        &self.forced
    }

    /// The soft pseudo-energy charged when base `i` is unpaired
    /// (0.0 outside range or when no soft term is set).
    pub fn soft_unpaired(&self, i: usize) -> f64 {
        self.soft_unpaired.get(i).copied().unwrap_or(0.0)
    }

    /// Sets the soft unpaired pseudo-energy for base `i` directly.
    ///
    /// # Errors
    /// [`RnaStructError::Invalid`] if `i >= n`.
    pub fn set_soft_unpaired(&mut self, i: usize, energy: f64) -> Result<()> {
        if i >= self.n {
            return Err(RnaStructError::invalid(
                "position",
                "soft-constraint index out of range",
            ));
        }
        self.soft_unpaired[i] = energy;
        Ok(())
    }

    /// Builds a constraint set from a vector of SHAPE reactivities
    /// using the Deigan *et al.* (2009) log model.
    ///
    /// For each position with reactivity `r >= 0`, the per-nucleotide
    /// pseudo-energy is `m · ln(r + 1) + b`. The Deigan model adds
    /// this term to each *paired* nucleotide; equivalently — and this
    /// is the form a pair-free soft constraint takes — we add the
    /// negation to each *unpaired* nucleotide, so a highly reactive
    /// base (large `r`) is rewarded for being unpaired. Negative
    /// reactivities (missing data) contribute zero.
    ///
    /// # Errors
    /// [`RnaStructError::Invalid`] if `reactivities.len() != n`.
    pub fn from_shape(n: usize, reactivities: &[f64]) -> Result<Self> {
        Self::from_shape_with(n, reactivities, DEIGAN_SLOPE, DEIGAN_INTERCEPT)
    }

    /// [`FoldConstraints::from_shape`] with explicit Deigan slope `m`
    /// and intercept `b`.
    ///
    /// # Errors
    /// [`RnaStructError::Invalid`] if `reactivities.len() != n`.
    pub fn from_shape_with(
        n: usize,
        reactivities: &[f64],
        slope: f64,
        intercept: f64,
    ) -> Result<Self> {
        if reactivities.len() != n {
            return Err(RnaStructError::invalid(
                "reactivities",
                format!(
                    "{} SHAPE values for a length-{n} sequence",
                    reactivities.len()
                ),
            ));
        }
        let mut c = Self::none(n);
        for (i, &r) in reactivities.iter().enumerate() {
            if r < 0.0 {
                continue; // missing data
            }
            // Deigan paired-nucleotide term:
            let paired_term = slope * (r + 1.0).ln() + intercept;
            // As a soft *unpaired* bonus this is the negation: a
            // reactive base gains energy (negative => stabilising)
            // for being unpaired.
            c.soft_unpaired[i] = -paired_term;
        }
        Ok(c)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn unconstrained_allows_everything() {
        let c = FoldConstraints::none(10);
        assert!(c.is_unconstrained());
        assert!(c.pair_allowed(0, 9));
        assert!(c.unpaired_allowed(5));
    }

    #[test]
    fn force_unpaired_blocks_pairs() {
        let mut c = FoldConstraints::none(10);
        c.force_unpaired(3).unwrap();
        // a force-unpaired base is still "allowed to be unpaired"
        assert!(c.unpaired_allowed(3));
        // ...but it can no longer take part in any pair
        assert!(!c.pair_allowed(3, 8));
        assert!(c.pair_allowed(2, 8));
    }

    #[test]
    fn force_paired_blocks_unpaired() {
        let mut c = FoldConstraints::none(10);
        c.force_paired(4).unwrap();
        assert!(!c.unpaired_allowed(4));
        assert!(c.pair_allowed(4, 9));
    }

    #[test]
    fn forbid_pair_works() {
        let mut c = FoldConstraints::none(10);
        c.forbid_pair(1, 8).unwrap();
        assert!(!c.pair_allowed(1, 8));
        assert!(!c.pair_allowed(8, 1)); // order-insensitive
        assert!(c.pair_allowed(1, 7));
    }

    #[test]
    fn force_pair_records_and_marks() {
        let mut c = FoldConstraints::none(10);
        c.force_pair(0, 9).unwrap();
        assert_eq!(c.forced_pairs(), &[(0, 9)]);
        assert!(!c.unpaired_allowed(0));
        assert!(!c.unpaired_allowed(9));
    }

    #[test]
    fn out_of_range_rejected() {
        let mut c = FoldConstraints::none(5);
        assert!(c.force_paired(9).is_err());
        assert!(c.forbid_pair(0, 99).is_err());
        assert!(c.force_pair(3, 3).is_err());
    }

    #[test]
    fn shape_reactive_bases_get_unpaired_bonus() {
        // high reactivity at position 2, zero elsewhere
        let react = vec![0.0, 0.0, 3.0, 0.0];
        let c = FoldConstraints::from_shape(4, &react).unwrap();
        // a reactive base gets a stabilising (negative) unpaired bonus
        assert!(c.soft_unpaired(2) < 0.0);
        // an unreactive base: m·ln(1)+b = b = -0.6 -> soft = +0.6
        assert!((c.soft_unpaired(0) - 0.6).abs() < 1e-9);
    }

    #[test]
    fn shape_negative_reactivity_is_ignored() {
        let react = vec![-1.0, -999.0];
        let c = FoldConstraints::from_shape(2, &react).unwrap();
        assert_eq!(c.soft_unpaired(0), 0.0);
        assert_eq!(c.soft_unpaired(1), 0.0);
    }

    #[test]
    fn shape_length_must_match() {
        assert!(FoldConstraints::from_shape(5, &[1.0, 2.0]).is_err());
    }
}
