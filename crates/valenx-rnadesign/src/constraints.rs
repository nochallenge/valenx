//! Feature 17 — the constrained-design layer.
//!
//! Real RNA design is never unconstrained. A synthesised construct must
//! avoid restriction sites and poly-runs, must fall in a
//! synthesis-friendly GC band, and often has positions that are
//! **locked** — fixed nucleotides the designer may not touch (an
//! engineered aptamer core, a primer-binding footprint, a ribosome
//! contact). This module turns the existing
//! [`DesignConstraints`] into an active
//! constraint layer the principled designers — the ensemble-defect
//! inverse folder ([`crate::inverse`]) and, where applicable, the
//! LinearDesign optimiser ([`crate::lineardesign`]) — honour.
//!
//! A [`DesignConstraintSet`] exposes three things:
//!
//! - **Locked positions** — a per-position `Option<base>`. A locked
//!   position is *never* mutated by a designer and is held at its
//!   fixed base throughout. Locked positions are written into a
//!   [`DesignConstraints`] as `required_subsequences` entries of the
//!   form `"@<index>=<base>"` (e.g. `"@0=G"`), so no schema change is
//!   needed.
//! - **A hard predicate** [`satisfies`](DesignConstraintSet::satisfies)
//!   — `true` when a sequence breaks no constraint (GC in range, no
//!   forbidden motif, no over-long homopolymer, every locked position
//!   correct).
//! - **A soft penalty** [`soft_penalty`](DesignConstraintSet::soft_penalty)
//!   — a non-negative number, `0` when every constraint holds, rising
//!   with the size of each violation. A designer adds this to its
//!   objective so its local search is *steered* toward a legal
//!   sequence rather than only rejecting illegal ones.
//!
//! ## v1 scope
//!
//! Forbidden motifs and forbidden subsequences are matched literally
//! (transcribed to RNA); restriction-enzyme sites are matched via the
//! optimiser's existing scan ([`crate::optimize::forbidden_motif_scan`])
//! when the caller wires it — this module's `soft_penalty` covers the
//! GC band, homopolymer cap, literal forbidden motifs and locked
//! positions, the constraints a sequence-level mutation search can act
//! on directly. The GC and homopolymer penalties are transparent
//! computed terms, not trained models.

use crate::goal::DesignConstraints;
use serde::{Deserialize, Serialize};

/// A constraint set lifted into an active design-time layer.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct DesignConstraintSet {
    /// The underlying declarative constraints.
    constraints: DesignConstraints,
    /// `locked[i] == Some(base)` fixes position `i` to `base` (ASCII
    /// `A C G U`); `None` leaves it free. Length is the design length.
    locked: Vec<Option<u8>>,
}

impl DesignConstraintSet {
    /// Builds a constraint set for a design of `length` nt from a
    /// [`DesignConstraints`]. Locked positions are parsed from the
    /// constraints' `required_subsequences` (`"@<index>=<base>"` form).
    pub fn new(constraints: &DesignConstraints, length: usize) -> Self {
        let mut locked: Vec<Option<u8>> = vec![None; length];
        for (idx, base) in parse_locked(constraints, length) {
            locked[idx] = Some(base);
        }
        DesignConstraintSet {
            constraints: constraints.clone(),
            locked,
        }
    }

    /// The design length the set was built for.
    pub fn len(&self) -> usize {
        self.locked.len()
    }

    /// `true` when the design length is zero.
    pub fn is_empty(&self) -> bool {
        self.locked.is_empty()
    }

    /// The underlying declarative constraints.
    pub fn constraints(&self) -> &DesignConstraints {
        &self.constraints
    }

    /// The locked base at position `i` (`None` if free / out of range).
    pub fn locked_at(&self, i: usize) -> Option<u8> {
        self.locked.get(i).copied().flatten()
    }

    /// `true` if position `i` is free (not locked).
    pub fn is_free(&self, i: usize) -> bool {
        self.locked_at(i).is_none()
    }

    /// The number of locked positions.
    pub fn locked_count(&self) -> usize {
        self.locked.iter().filter(|l| l.is_some()).count()
    }

    /// Applies every locked position to `seq` in place, overriding
    /// whatever base is there — used to seed / repair a candidate.
    ///
    /// `seq` must be at least as long as the constraint set.
    pub fn apply_locks(&self, seq: &mut [u8]) {
        for (i, lk) in self.locked.iter().enumerate() {
            if let (Some(b), Some(slot)) = (lk, seq.get_mut(i)) {
                *slot = *b;
            }
        }
    }

    /// `true` when `seq` breaks **no** constraint: every locked position
    /// holds its base, the GC fraction is in the configured band, no
    /// homopolymer run exceeds the cap, and no forbidden motif /
    /// subsequence occurs.
    pub fn satisfies(&self, seq: &[u8]) -> bool {
        self.locks_ok(seq) && self.gc_ok(seq) && self.homopolymer_ok(seq) && self.motifs_ok(seq)
    }

    /// A non-negative soft penalty for `seq` — `0` when every constraint
    /// holds, rising with each violation. Designed to be *added* to a
    /// designer's objective so its local search is steered toward a
    /// legal sequence. The terms:
    ///
    /// - **locked**: `+1` per locked position holding the wrong base;
    /// - **GC**: `+ (gc excess) · length` (a fractional excess scaled to
    ///   nucleotide units);
    /// - **homopolymer**: `+ (run − cap)` per over-long run;
    /// - **motif**: `+1` per forbidden-motif / subsequence occurrence.
    pub fn soft_penalty(&self, seq: &[u8]) -> f64 {
        let mut penalty = 0.0;

        // Locked positions.
        for (i, lk) in self.locked.iter().enumerate() {
            if let (Some(b), Some(&actual)) = (lk, seq.get(i)) {
                if !actual.eq_ignore_ascii_case(b) {
                    penalty += 1.0;
                }
            }
        }

        // GC band.
        let gc = gc_fraction(seq);
        if gc < self.constraints.gc_min {
            penalty += (self.constraints.gc_min - gc) * seq.len() as f64;
        } else if gc > self.constraints.gc_max {
            penalty += (gc - self.constraints.gc_max) * seq.len() as f64;
        }

        // Homopolymer cap.
        for run in homopolymer_runs(seq) {
            if run > self.constraints.max_homopolymer {
                penalty += (run - self.constraints.max_homopolymer) as f64;
            }
        }

        // Forbidden motifs / subsequences.
        for motif in self
            .constraints
            .forbidden_motifs
            .iter()
            .chain(self.constraints.forbidden_subsequences.iter())
        {
            let needle = to_rna(motif.as_bytes());
            if needle.is_empty() {
                continue;
            }
            penalty += count_occurrences(seq, &needle) as f64;
        }

        penalty
    }

    /// `true` if every locked position holds its base.
    pub fn locks_ok(&self, seq: &[u8]) -> bool {
        self.locked.iter().enumerate().all(|(i, lk)| match lk {
            None => true,
            Some(b) => seq
                .get(i)
                .map(|a| a.eq_ignore_ascii_case(b))
                .unwrap_or(false),
        })
    }

    /// `true` if `seq`'s GC fraction is in the configured band.
    pub fn gc_ok(&self, seq: &[u8]) -> bool {
        let gc = gc_fraction(seq);
        gc >= self.constraints.gc_min - 1e-9 && gc <= self.constraints.gc_max + 1e-9
    }

    /// `true` if no homopolymer run exceeds the configured cap.
    pub fn homopolymer_ok(&self, seq: &[u8]) -> bool {
        homopolymer_runs(seq)
            .into_iter()
            .all(|r| r <= self.constraints.max_homopolymer)
    }

    /// `true` if no forbidden motif / subsequence occurs in `seq`.
    pub fn motifs_ok(&self, seq: &[u8]) -> bool {
        for motif in self
            .constraints
            .forbidden_motifs
            .iter()
            .chain(self.constraints.forbidden_subsequences.iter())
        {
            let needle = to_rna(motif.as_bytes());
            if !needle.is_empty() && count_occurrences(seq, &needle) > 0 {
                return false;
            }
        }
        true
    }
}

/// Parses locked `(position, base)` pairs from a constraint set's
/// `required_subsequences` entries of the form `"@<index>=<base>"`.
///
/// Entries that do not parse, or whose index is out of range for
/// `length`, are ignored.
pub fn parse_locked(constraints: &DesignConstraints, length: usize) -> Vec<(usize, u8)> {
    let mut out = Vec::new();
    for entry in &constraints.required_subsequences {
        if !entry.starts_with('@') {
            continue;
        }
        let rest = &entry[1..];
        let Some(eq) = rest.find('=') else { continue };
        let (idx_str, base_str) = rest.split_at(eq);
        let base_str = &base_str[1..];
        let Ok(idx) = idx_str.trim().parse::<usize>() else {
            continue;
        };
        if idx >= length {
            continue;
        }
        let b = match base_str.trim().to_ascii_uppercase().as_bytes().first() {
            Some(&b @ (b'A' | b'C' | b'G' | b'U')) => b,
            Some(&b'T') => b'U',
            _ => continue,
        };
        out.push((idx, b));
    }
    out
}

/// Builds a `required_subsequences` entry that locks `position` to
/// `base` — the inverse of [`parse_locked`]. A convenience for callers
/// assembling a [`DesignConstraints`] with locked positions.
pub fn lock_entry(position: usize, base: u8) -> String {
    format!("@{}={}", position, base.to_ascii_uppercase() as char)
}

/// GC fraction of an RNA / DNA sequence in `[0, 1]`.
fn gc_fraction(seq: &[u8]) -> f64 {
    if seq.is_empty() {
        return 0.0;
    }
    let gc = seq
        .iter()
        .filter(|&&b| matches!(b.to_ascii_uppercase(), b'G' | b'C'))
        .count();
    gc as f64 / seq.len() as f64
}

/// The lengths of every maximal homopolymer run in `seq`.
fn homopolymer_runs(seq: &[u8]) -> Vec<usize> {
    let mut runs = Vec::new();
    let mut i = 0;
    while i < seq.len() {
        let b = seq[i].to_ascii_uppercase();
        let start = i;
        while i < seq.len() && seq[i].to_ascii_uppercase() == b {
            i += 1;
        }
        runs.push(i - start);
    }
    runs
}

/// Transcribes raw bytes to uppercase RNA (`T`→`U`).
fn to_rna(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .map(|&b| match b.to_ascii_uppercase() {
            b'T' => b'U',
            other => other,
        })
        .collect()
}

/// Counts (overlapping) occurrences of `needle` in `haystack`
/// (case-insensitive).
fn count_occurrences(haystack: &[u8], needle: &[u8]) -> usize {
    if needle.is_empty() || needle.len() > haystack.len() {
        return 0;
    }
    haystack
        .windows(needle.len())
        .filter(|w| w.eq_ignore_ascii_case(needle))
        .count()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A `DesignConstraints` carrying the given `required_subsequences`
    /// entries, everything else default.
    fn with_required(entries: Vec<&str>) -> DesignConstraints {
        DesignConstraints {
            required_subsequences: entries.into_iter().map(String::from).collect(),
            ..DesignConstraints::default()
        }
    }

    #[test]
    fn parses_locked_positions() {
        let c = with_required(vec![
            "@0=G", "@5=c", "@2=T", // T folds to U
            "junk", "@99=A", // out of range
        ]);
        let set = DesignConstraintSet::new(&c, 16);
        assert_eq!(set.locked_at(0), Some(b'G'));
        assert_eq!(set.locked_at(5), Some(b'C'));
        assert_eq!(set.locked_at(2), Some(b'U'));
        assert_eq!(set.locked_at(99), None);
        assert_eq!(set.locked_count(), 3);
    }

    #[test]
    fn lock_entry_round_trips() {
        let entry = lock_entry(7, b'g');
        assert_eq!(entry, "@7=G");
        let c = DesignConstraints {
            required_subsequences: vec![entry],
            ..DesignConstraints::default()
        };
        let set = DesignConstraintSet::new(&c, 10);
        assert_eq!(set.locked_at(7), Some(b'G'));
    }

    #[test]
    fn apply_locks_overrides_bases() {
        let c = with_required(vec!["@0=G", "@3=C"]);
        let set = DesignConstraintSet::new(&c, 6);
        let mut seq = b"AAAAAA".to_vec();
        set.apply_locks(&mut seq);
        assert_eq!(seq, b"GAACAA");
    }

    #[test]
    fn satisfies_a_clean_sequence() {
        let c = DesignConstraints::default().with_gc_range(0.3, 0.7);
        let set = DesignConstraintSet::new(&c, 8);
        // 50 % GC, no runs, no motifs.
        assert!(set.satisfies(b"ACGUACGU"));
    }

    #[test]
    fn detects_gc_violation() {
        let c = DesignConstraints::default().with_gc_range(0.3, 0.6);
        let set = DesignConstraintSet::new(&c, 8);
        assert!(!set.gc_ok(b"GGGGGGGG"));
        assert!(!set.satisfies(b"GGGGGGGG"));
        assert!(set.soft_penalty(b"GGGGGGGG") > 0.0);
    }

    #[test]
    fn detects_forbidden_motif() {
        let c = DesignConstraints::default()
            .with_gc_range(0.0, 1.0)
            .forbid_motif("GAAUUC");
        let set = DesignConstraintSet::new(&c, 12);
        assert!(!set.motifs_ok(b"AAAGAAUUCAAA"));
        assert!(set.soft_penalty(b"AAAGAAUUCAAA") >= 1.0);
        assert!(set.motifs_ok(b"AAAAAAAAAAAA"));
    }

    #[test]
    fn detects_homopolymer_violation() {
        let mut c = DesignConstraints::default().with_gc_range(0.0, 1.0);
        c.max_homopolymer = 4;
        let set = DesignConstraintSet::new(&c, 12);
        assert!(!set.homopolymer_ok(b"ACGAAAAAAACG"));
        assert!(set.homopolymer_ok(b"ACGAAAACGACG"));
    }

    #[test]
    fn detects_locked_violation() {
        let c = DesignConstraints {
            required_subsequences: vec!["@0=G".to_string()],
            ..DesignConstraints::default().with_gc_range(0.0, 1.0)
        };
        let set = DesignConstraintSet::new(&c, 4);
        assert!(set.locks_ok(b"GAAA"));
        assert!(!set.locks_ok(b"AAAA"));
        assert!((set.soft_penalty(b"AAAA") - 1.0).abs() < 1e-9);
        assert!(set.soft_penalty(b"GAAA") < 1e-9);
    }

    #[test]
    fn soft_penalty_is_zero_for_a_clean_sequence() {
        let c = DesignConstraints::default().with_gc_range(0.3, 0.7);
        let set = DesignConstraintSet::new(&c, 8);
        assert!(set.soft_penalty(b"ACGUACGU") < 1e-9);
    }

    #[test]
    fn empty_set() {
        let set = DesignConstraintSet::new(&DesignConstraints::default(), 0);
        assert!(set.is_empty());
        assert_eq!(set.len(), 0);
    }
}
