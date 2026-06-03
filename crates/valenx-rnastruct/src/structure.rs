//! The [`Structure`] model — an RNA secondary structure as a set of
//! base pairs over a fixed-length sequence.
//!
//! An RNA secondary structure is, formally, a set of pairs `(i, j)`
//! with `i < j` over the `0..n` positions of a length-`n` sequence,
//! such that every position appears in at most one pair. A structure
//! is *nested* (pseudoknot-free) if no two pairs `(i, j)` and `(k, l)`
//! "cross" — i.e. `i < k < j < l`. A structure with crossing pairs
//! contains *pseudoknots*.
//!
//! [`Structure`] stores the partner array (`partner[i]` is the index
//! paired with `i`, or `None`) which is the form every algorithm in
//! this crate consumes. Construction always validates: indices in
//! range, `i != j`, no position paired twice.
//!
//! Dot-bracket I/O lives in [`Structure::from_dot_bracket`] /
//! [`Structure::to_dot_bracket`]; ct-file and bpseq I/O live in
//! [`crate::io`].

use crate::error::{Result, RnaStructError};
use serde::{Deserialize, Serialize};

/// A single base pair `(i, j)` with `i < j` (0-based positions).
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash, Serialize, Deserialize)]
pub struct BasePair {
    /// The 5′ (lower-index) partner.
    pub i: usize,
    /// The 3′ (higher-index) partner.
    pub j: usize,
}

impl BasePair {
    /// Builds a pair, ordering the two indices so `i < j`.
    ///
    /// # Errors
    /// [`RnaStructError::Structure`] if `a == b`.
    pub fn new(a: usize, b: usize) -> Result<Self> {
        if a == b {
            return Err(RnaStructError::structure(format!(
                "a base pair cannot pair position {a} with itself"
            )));
        }
        let (i, j) = if a < b { (a, b) } else { (b, a) };
        Ok(BasePair { i, j })
    }

    /// The span `j - i` — the number of positions enclosed plus one.
    pub fn span(&self) -> usize {
        self.j - self.i
    }

    /// `true` if this pair and `other` cross (form a pseudoknot):
    /// `i < k < j < l` (or the mirror). Pairs that nest or are
    /// disjoint do not cross.
    pub fn crosses(&self, other: &BasePair) -> bool {
        let (a, b) = (self, other);
        (a.i < b.i && b.i < a.j && a.j < b.j) || (b.i < a.i && a.i < b.j && b.j < a.j)
    }
}

/// An RNA secondary structure: a partner array over `n` positions.
///
/// `partner[i] == Some(j)` means position `i` pairs with `j` (and then
/// `partner[j] == Some(i)`). `partner[i] == None` means `i` is
/// unpaired. The invariant — symmetry and at-most-one-pair-per-base —
/// is upheld by every constructor and mutator.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Structure {
    /// `partner[i]` = the index paired with `i`, or `None`.
    partner: Vec<Option<usize>>,
}

impl Structure {
    /// An all-unpaired structure of length `n`.
    pub fn empty(n: usize) -> Self {
        Structure {
            partner: vec![None; n],
        }
    }

    /// Builds a structure from a list of base pairs over `n` positions.
    ///
    /// # Errors
    /// [`RnaStructError::Structure`] if any index is `>= n` or any
    /// position is named in more than one pair.
    pub fn from_pairs(n: usize, pairs: &[BasePair]) -> Result<Self> {
        let mut partner = vec![None; n];
        for bp in pairs {
            if bp.i >= n || bp.j >= n {
                return Err(RnaStructError::structure(format!(
                    "base pair ({}, {}) out of range for length {n}",
                    bp.i, bp.j
                )));
            }
            if partner[bp.i].is_some() {
                return Err(RnaStructError::structure(format!(
                    "position {} is paired more than once",
                    bp.i
                )));
            }
            if partner[bp.j].is_some() {
                return Err(RnaStructError::structure(format!(
                    "position {} is paired more than once",
                    bp.j
                )));
            }
            partner[bp.i] = Some(bp.j);
            partner[bp.j] = Some(bp.i);
        }
        Ok(Structure { partner })
    }

    /// Builds a structure directly from a partner array, validating
    /// symmetry. Used by the folding DP traceback.
    ///
    /// # Errors
    /// [`RnaStructError::Structure`] if the array is not symmetric or
    /// any index is out of range.
    pub fn from_partner(partner: Vec<Option<usize>>) -> Result<Self> {
        let n = partner.len();
        for (i, &p) in partner.iter().enumerate() {
            if let Some(j) = p {
                if j >= n {
                    return Err(RnaStructError::structure(format!(
                        "partner index {j} out of range for length {n}"
                    )));
                }
                if j == i {
                    return Err(RnaStructError::structure(format!(
                        "position {i} is paired with itself"
                    )));
                }
                if partner[j] != Some(i) {
                    return Err(RnaStructError::structure(format!(
                        "partner array is not symmetric at position {i}"
                    )));
                }
            }
        }
        Ok(Structure { partner })
    }

    /// Number of positions (the sequence length this structure spans).
    pub fn len(&self) -> usize {
        self.partner.len()
    }

    /// `true` if the structure spans zero positions.
    pub fn is_empty(&self) -> bool {
        self.partner.is_empty()
    }

    /// The partner of position `i`, or `None` if unpaired / out of
    /// range.
    pub fn partner(&self, i: usize) -> Option<usize> {
        self.partner.get(i).copied().flatten()
    }

    /// The raw partner array.
    pub fn partner_array(&self) -> &[Option<usize>] {
        &self.partner
    }

    /// `true` if position `i` is paired.
    pub fn is_paired(&self, i: usize) -> bool {
        self.partner(i).is_some()
    }

    /// All base pairs, sorted ascending by `i`.
    pub fn pairs(&self) -> Vec<BasePair> {
        let mut out = Vec::new();
        for (i, &p) in self.partner.iter().enumerate() {
            if let Some(j) = p {
                if i < j {
                    out.push(BasePair { i, j });
                }
            }
        }
        out
    }

    /// The number of base pairs.
    pub fn n_pairs(&self) -> usize {
        self.partner.iter().filter(|p| p.is_some()).count() / 2
    }

    /// Adds a pair `(i, j)`. Both positions must currently be unpaired.
    ///
    /// # Errors
    /// [`RnaStructError::Structure`] if an index is out of range, the
    /// two indices are equal, or either is already paired.
    pub fn add_pair(&mut self, i: usize, j: usize) -> Result<()> {
        let n = self.partner.len();
        if i >= n || j >= n {
            return Err(RnaStructError::structure(format!(
                "pair ({i}, {j}) out of range for length {n}"
            )));
        }
        if i == j {
            return Err(RnaStructError::structure("cannot pair a position with itself"));
        }
        if self.partner[i].is_some() || self.partner[j].is_some() {
            return Err(RnaStructError::structure(format!(
                "position {i} or {j} is already paired"
            )));
        }
        self.partner[i] = Some(j);
        self.partner[j] = Some(i);
        Ok(())
    }

    /// Removes the pair touching position `i` (if any). Returns the
    /// former partner.
    pub fn remove_pair(&mut self, i: usize) -> Option<usize> {
        let j = self.partner.get(i).copied().flatten()?;
        self.partner[i] = None;
        self.partner[j] = None;
        Some(j)
    }

    /// `true` if the structure is nested (pseudoknot-free): no two
    /// pairs cross.
    pub fn is_nested(&self) -> bool {
        let pairs = self.pairs();
        for a in 0..pairs.len() {
            for b in (a + 1)..pairs.len() {
                if pairs[a].crosses(&pairs[b]) {
                    return false;
                }
            }
        }
        true
    }

    /// `true` if the structure contains at least one crossing pair
    /// (a pseudoknot). The negation of [`Structure::is_nested`].
    pub fn has_pseudoknot(&self) -> bool {
        !self.is_nested()
    }

    /// Parses a dot-bracket string into a structure.
    ///
    /// Recognises the standard nesting brackets `()`, plus three extra
    /// pseudoknot bracket families `[]`, `{}` and `<>` (the
    /// dot-bracket convention used by Rfam / pKiss / ViennaRNA). Each
    /// family is balanced independently, so two different families may
    /// legally cross. `.` marks an unpaired position; spaces are
    /// ignored.
    ///
    /// # Errors
    /// [`RnaStructError::Parse`] on an unbalanced bracket or an
    /// unrecognised character.
    pub fn from_dot_bracket(db: &str) -> Result<Self> {
        // Four bracket families: (open, close).
        const FAMILIES: [(char, char); 4] =
            [('(', ')'), ('[', ']'), ('{', '}'), ('<', '>')];

        let chars: Vec<char> = db.chars().filter(|c| !c.is_whitespace()).collect();
        let n = chars.len();
        let mut partner: Vec<Option<usize>> = vec![None; n];
        // One stack per bracket family.
        let mut stacks: [Vec<usize>; 4] = Default::default();

        for (pos, &c) in chars.iter().enumerate() {
            if c == '.' {
                continue;
            }
            let mut handled = false;
            for (fam, (open, close)) in FAMILIES.iter().enumerate() {
                if c == *open {
                    stacks[fam].push(pos);
                    handled = true;
                    break;
                }
                if c == *close {
                    match stacks[fam].pop() {
                        Some(o) => {
                            partner[o] = Some(pos);
                            partner[pos] = Some(o);
                        }
                        None => {
                            return Err(RnaStructError::parse(
                                "dot-bracket",
                                format!("unmatched `{c}` at position {pos}"),
                            ));
                        }
                    }
                    handled = true;
                    break;
                }
            }
            if !handled {
                return Err(RnaStructError::parse(
                    "dot-bracket",
                    format!("unrecognised character `{c}` at position {pos}"),
                ));
            }
        }
        for (fam, stack) in stacks.iter().enumerate() {
            if !stack.is_empty() {
                return Err(RnaStructError::parse(
                    "dot-bracket",
                    format!(
                        "{} unmatched `{}` bracket(s)",
                        stack.len(),
                        FAMILIES[fam].0
                    ),
                ));
            }
        }
        Ok(Structure { partner })
    }

    /// Renders the structure as a dot-bracket string.
    ///
    /// A pseudoknot-free structure uses only `()`. When pairs cross,
    /// later (by position) crossing pairs are assigned to the `[]`,
    /// `{}` then `<>` families in turn — a greedy colouring that
    /// reproduces the ViennaRNA convention for the common one- and
    /// two-page pseudoknots. If more than four mutually crossing
    /// "pages" are needed the extra pairs fall back to `.` rather than
    /// emit an illegal string; this is rare for real RNA.
    pub fn to_dot_bracket(&self) -> String {
        const OPENS: [char; 4] = ['(', '[', '{', '<'];
        const CLOSES: [char; 4] = [')', ']', '}', '>'];
        let n = self.partner.len();
        let mut out = vec!['.'; n];
        let pairs = self.pairs();

        // Assign each pair to the lowest-numbered page on which it
        // does not cross an already-placed pair (classic interval
        // graph greedy colouring).
        let mut page: Vec<usize> = Vec::with_capacity(pairs.len());
        for (idx, p) in pairs.iter().enumerate() {
            let mut chosen = usize::MAX;
            'page: for cand in 0..4 {
                for (prev_idx, prev) in pairs.iter().enumerate().take(idx) {
                    if page[prev_idx] == cand && p.crosses(prev) {
                        continue 'page;
                    }
                }
                chosen = cand;
                break;
            }
            page.push(chosen);
        }
        for (p, &pg) in pairs.iter().zip(page.iter()) {
            if pg < 4 {
                out[p.i] = OPENS[pg];
                out[p.j] = CLOSES[pg];
            }
            // pg == usize::MAX: more than 4 pages — leave as `.`.
        }
        out.into_iter().collect()
    }

    /// Validates the structure against a sequence of length `seq_len`:
    /// the lengths must match and (optionally) every pair must be a
    /// canonical / wobble pair.
    ///
    /// `require_canonical` enforces that each pair is A-U, G-C or G-U
    /// (in either order) — the only pairs the energy model scores.
    ///
    /// # Errors
    /// [`RnaStructError::Structure`] on a length mismatch or, when
    /// `require_canonical`, a non-canonical pair.
    pub fn validate_against(
        &self,
        seq: &[u8],
        require_canonical: bool,
    ) -> Result<()> {
        if seq.len() != self.partner.len() {
            return Err(RnaStructError::structure(format!(
                "structure length {} does not match sequence length {}",
                self.partner.len(),
                seq.len()
            )));
        }
        if require_canonical {
            for bp in self.pairs() {
                if !crate::fold::energy::can_pair(seq[bp.i], seq[bp.j]) {
                    return Err(RnaStructError::structure(format!(
                        "non-canonical pair {}-{} at positions ({}, {})",
                        seq[bp.i] as char, seq[bp.j] as char, bp.i, bp.j
                    )));
                }
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn base_pair_ordering_and_crossing() {
        let p = BasePair::new(5, 2).unwrap();
        assert_eq!((p.i, p.j), (2, 5));
        assert_eq!(p.span(), 3);
        assert!(BasePair::new(3, 3).is_err());

        let a = BasePair { i: 0, j: 5 };
        let b = BasePair { i: 3, j: 8 }; // crosses a
        let c = BasePair { i: 1, j: 4 }; // nested in a
        assert!(a.crosses(&b));
        assert!(!a.crosses(&c));
    }

    #[test]
    fn from_pairs_validates() {
        let s = Structure::from_pairs(6, &[BasePair { i: 0, j: 5 }, BasePair { i: 1, j: 4 }])
            .unwrap();
        assert_eq!(s.n_pairs(), 2);
        assert_eq!(s.partner(0), Some(5));
        assert!(s.is_paired(1));
        assert!(!s.is_paired(2));

        // double-pairing is rejected
        assert!(Structure::from_pairs(
            6,
            &[BasePair { i: 0, j: 5 }, BasePair { i: 0, j: 3 }]
        )
        .is_err());
        // out of range
        assert!(Structure::from_pairs(4, &[BasePair { i: 0, j: 9 }]).is_err());
    }

    #[test]
    fn dot_bracket_roundtrip_nested() {
        let db = "((..))";
        let s = Structure::from_dot_bracket(db).unwrap();
        assert_eq!(s.n_pairs(), 2);
        assert_eq!(s.partner(0), Some(5));
        assert_eq!(s.partner(1), Some(4));
        assert!(s.is_nested());
        assert_eq!(s.to_dot_bracket(), db);
    }

    #[test]
    fn dot_bracket_pseudoknot() {
        // classic H-type pseudoknot: stem 1 () crosses stem 2 []
        let db = "((..[[..))..]]";
        let s = Structure::from_dot_bracket(db).unwrap();
        assert!(s.has_pseudoknot());
        // round-trips: greedy colouring re-derives the two pages
        let again = Structure::from_dot_bracket(&s.to_dot_bracket()).unwrap();
        assert_eq!(s.pairs(), again.pairs());
    }

    #[test]
    fn dot_bracket_rejects_unbalanced() {
        assert!(Structure::from_dot_bracket("((.)").is_err());
        assert!(Structure::from_dot_bracket("(.))").is_err());
        assert!(Structure::from_dot_bracket("(.x)").is_err());
    }

    #[test]
    fn add_remove_pair() {
        let mut s = Structure::empty(8);
        s.add_pair(0, 7).unwrap();
        assert!(s.is_paired(0));
        assert!(s.add_pair(0, 3).is_err()); // 0 already paired
        assert_eq!(s.remove_pair(7), Some(0));
        assert!(!s.is_paired(0));
    }

    #[test]
    fn validate_against_sequence() {
        let s = Structure::from_dot_bracket("(((...)))").unwrap();
        let seq = b"GGGAAACCC";
        assert!(s.validate_against(seq, true).is_ok());
        // length mismatch
        assert!(s.validate_against(b"GGG", true).is_err());
        // non-canonical pair A-A
        let bad = Structure::from_dot_bracket("(.)").unwrap();
        assert!(bad.validate_against(b"AAA", true).is_err());
    }
}
