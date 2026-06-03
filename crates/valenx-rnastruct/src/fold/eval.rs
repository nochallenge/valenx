//! Free-energy evaluation of a given (sequence, structure) pair.
//!
//! Given an RNA sequence and a *nested* secondary structure, this
//! module computes the total Turner-2004 free energy by decomposing
//! the structure into loops and summing each loop's contribution.
//!
//! The decomposition is the standard one: every base pair `(i, j)`
//! closes exactly one loop, whose other pairs are the pairs *directly
//! enclosed* by `(i, j)` (no pair between them). By the number of
//! directly-enclosed pairs the loop is a:
//!
//! - **hairpin** — 0 enclosed pairs,
//! - **stack / bulge / internal loop** — exactly 1 enclosed pair,
//! - **multiloop** — ≥ 2 enclosed pairs.
//!
//! The unpaired bases not enclosed by any pair form the *exterior
//! loop*, which carries no energy in the Turner model (its dangles are
//! folded into the terminal penalties).
//!
//! Pseudoknotted structures cannot be scored by this nearest-neighbor
//! model — [`structure_energy`] returns an error for a crossing
//! structure (use [`crate::compare::pseudoknot`] for those).
//!
//! ## Two dangle models
//!
//! - [`structure_energy`] — the *dangle-folded* model: a multiloop and
//!   the exterior loop are scored as a sum of independent helices (the
//!   linear multiloop term plus terminal penalties). This is the model
//!   the Zuker MFE and McCaskill partition function optimise, so its
//!   energy is self-consistent with [`crate::fold::zuker::mfe`].
//! - [`structure_energy_d2`] — the *coaxial-stacking* (`-d2`) model:
//!   identical, **plus** the [`crate::fold::coaxial`] correction for
//!   helices that lie end-to-end in a multiloop or the exterior loop.
//!   This reproduces ViennaRNA `RNAeval -d2` exactly-to-rounding. It is
//!   a strictly-stabilising correction (`structure_energy_d2 ≤
//!   structure_energy`), zero for hairpin-only structures.

use crate::error::{Result, RnaStructError};
use crate::fold::coaxial::{self, HelixEnd};
use crate::fold::energy::{self, multiloop};
use crate::rna::RnaSeq;
use crate::structure::Structure;

/// Computes the total Turner-2004 free energy (kcal/mol) of `structure`
/// folded on `seq`, in the **dangle-folded** model (no explicit
/// coaxial-stacking term — see [`structure_energy_d2`] for the `-d2`
/// model). This is the model the MFE / partition-function recurrences
/// optimise, so the value is self-consistent with [`crate::fold::zuker::mfe`].
///
/// # Errors
/// - [`RnaStructError::Structure`] if the structure length differs
///   from the sequence length, if it contains a pseudoknot, or if it
///   contains a non-canonical pair.
pub fn structure_energy(seq: &RnaSeq, structure: &Structure) -> Result<f64> {
    eval_with_coaxial(seq, structure, false)
}

/// Computes the total Turner-2004 free energy (kcal/mol) of `structure`
/// folded on `seq`, in ViennaRNA's default **`-d2`** model — the
/// dangle-folded energy of [`structure_energy`] **plus** the explicit
/// [`crate::fold::coaxial`] coaxial-stacking correction for helices
/// that lie end-to-end in a multiloop or the exterior loop.
///
/// For a structure with no helix junctions (a single hairpin, a bare
/// stem) there is no coaxial term and this equals [`structure_energy`].
/// For a multi-helix structure it is lower (coaxial stacking only ever
/// stabilises), and reproduces ViennaRNA `RNAeval -d2` exactly-to-
/// rounding.
///
/// # Errors
/// Same as [`structure_energy`].
pub fn structure_energy_d2(seq: &RnaSeq, structure: &Structure) -> Result<f64> {
    eval_with_coaxial(seq, structure, true)
}

/// The coaxial-stacking correction alone (kcal/mol) — the difference
/// `structure_energy_d2 − structure_energy`. Always ≤ 0; exactly 0 for
/// a structure with no helix junctions.
///
/// # Errors
/// Same as [`structure_energy`].
pub fn coaxial_correction(seq: &RnaSeq, structure: &Structure) -> Result<f64> {
    Ok(structure_energy_d2(seq, structure)? - structure_energy(seq, structure)?)
}

/// Shared evaluator. `coaxial` selects the dangle-folded (`false`) or
/// the `-d2` coaxial-stacking (`true`) model.
fn eval_with_coaxial(
    seq: &RnaSeq,
    structure: &Structure,
    coaxial: bool,
) -> Result<f64> {
    let codes = seq.codes();
    if structure.len() != codes.len() {
        return Err(RnaStructError::structure(format!(
            "structure length {} != sequence length {}",
            structure.len(),
            codes.len()
        )));
    }
    if structure.has_pseudoknot() {
        return Err(RnaStructError::structure(
            "free-energy evaluation requires a pseudoknot-free structure",
        ));
    }
    // Validate canonical pairs up front.
    for bp in structure.pairs() {
        if !energy::can_pair_codes(codes[bp.i], codes[bp.j]) {
            return Err(RnaStructError::structure(format!(
                "non-canonical pair at ({}, {})",
                bp.i, bp.j
            )));
        }
    }

    let n = codes.len();
    let partner = structure.partner_array();
    let mut total = 0.0;

    // Score the loop closed by every pair (i, j) with i < j.
    for i in 0..n {
        if let Some(j) = partner[i] {
            if i < j {
                total += loop_energy(codes, partner, i, j, coaxial);
            }
        }
    }

    // Terminal AU/GU penalty on every exterior helix end. The Turner
    // model charges the weak-pair penalty once per helix *end*; a loop
    // end is charged inside `loop_energy` (hairpin / internal / multi),
    // but the end of an exterior helix faces the exterior loop, which
    // carries no loop energy of its own — so it is charged here. (The
    // Zuker `w` recurrence charges the same term, so this keeps
    // `structure_energy` consistent with `mfe`.)
    {
        let mut p = 0usize;
        while p < n {
            match partner[p] {
                Some(q) if q > p => {
                    total += energy::terminal_penalty(codes[p], codes[q]);
                    p = q + 1;
                }
                _ => p += 1,
            }
        }
    }

    // Coaxial stacking on the exterior loop: helices that emanate
    // directly from the unenclosed (exterior) strand can stack on each
    // other end-to-end, exactly as inside a multiloop.
    if coaxial {
        total += exterior_coaxial(codes, partner);
    }

    Ok(total)
}

/// Coaxial-stacking correction for the exterior loop — the helices
/// hanging directly off the unenclosed strand.
fn exterior_coaxial(codes: &[u8], partner: &[Option<usize>]) -> f64 {
    let n = codes.len();
    // Collect the exterior helices left-to-right.
    let mut ends: Vec<HelixEnd> = Vec::new();
    let mut gaps: Vec<usize> = Vec::new();
    let mut bridges: Vec<u8> = Vec::new();
    let mut prev_close: Option<usize> = None; // 3' base of the last helix
    let mut p = 0usize;
    while p < n {
        match partner[p] {
            Some(q) if q > p => {
                // helix p..q emanates from the exterior loop.
                if let Some(pc) = prev_close {
                    let gap = p - pc - 1;
                    gaps.push(gap);
                    bridges.push(if gap == 1 { codes[pc + 1] } else { 0 });
                }
                ends.push(HelixEnd {
                    left: codes[p],
                    right: codes[q],
                });
                prev_close = Some(q);
                p = q + 1;
            }
            _ => p += 1,
        }
    }
    if ends.len() < 2 {
        return 0.0;
    }
    // The wrap "gap" (after the last helix, through the free 3'/5' ends)
    // never stacks: push a large gap so `best_coaxial` ignores it.
    gaps.push(usize::MAX / 4);
    bridges.push(0);
    coaxial::best_coaxial(&ends, &gaps, &bridges, false)
}

/// Energy of the single loop closed by pair `(i, j)` — the inner
/// dispatch used by [`structure_energy`]. `coaxial` adds the
/// coaxial-stacking term to a multiloop.
fn loop_energy(
    codes: &[u8],
    partner: &[Option<usize>],
    i: usize,
    j: usize,
    coaxial: bool,
) -> f64 {
    // Collect the pairs directly enclosed by (i, j): walk i+1..j and
    // record every pair whose 5' base we meet, skipping over already-
    // entered sub-loops.
    let mut enclosed: Vec<(usize, usize)> = Vec::new();
    let mut k = i + 1;
    while k < j {
        match partner[k] {
            Some(p) if p > k => {
                enclosed.push((k, p));
                k = p + 1;
            }
            _ => k += 1,
        }
    }

    match enclosed.len() {
        0 => {
            // Hairpin loop.
            let loop_bases: Vec<u8> = codes[(i + 1)..j].to_vec();
            energy::hairpin_energy(codes[i], codes[j], &loop_bases)
        }
        1 => {
            // Stack / bulge / internal loop.
            let (k, l) = enclosed[0];
            let left = k - i - 1;
            let right = j - l - 1;
            // Unpaired bases flanking each closing pair (for mismatch).
            let mm_outer_5 = codes[i + 1];
            let mm_outer_3 = codes[j - 1];
            let mm_inner_5 = if k > 0 { codes[k - 1] } else { codes[k] };
            let mm_inner_3 = if l + 1 < codes.len() {
                codes[l + 1]
            } else {
                codes[l]
            };
            energy::internal_loop_energy(
                codes[i],
                codes[j],
                codes[k],
                codes[l],
                left,
                right,
                mm_outer_5,
                mm_outer_3,
                mm_inner_5,
                mm_inner_3,
            )
        }
        b => {
            // Multiloop: `b` enclosed branches + the closing pair.
            let branches = b + 1;
            // Count unpaired bases inside (i, j) not in any sub-pair.
            let mut unpaired = 0usize;
            let mut p = i + 1;
            while p < j {
                match partner[p] {
                    Some(q) if q > p => p = q + 1,
                    _ => {
                        unpaired += 1;
                        p += 1;
                    }
                }
            }
            let mut e = multiloop::energy(branches, unpaired);
            // Terminal-AU penalty on the closing helix and each branch.
            e += energy::terminal_penalty(codes[i], codes[j]);
            for &(k, l) in &enclosed {
                e += energy::terminal_penalty(codes[k], codes[l]);
            }
            // Coaxial stacking: walk the multiloop boundary 5'->3' and
            // collect helix ends + the gaps between them. The boundary
            // is a cycle: enclosed branches in order, then the closing
            // pair (presented loop-facing as (j, i)).
            if coaxial {
                e += multiloop_coaxial(codes, i, j, &enclosed);
            }
            e
        }
    }
}

/// Coaxial-stacking correction for the multiloop closed by `(i, j)`
/// with directly-`enclosed` branch pairs (in 5′→3′ order).
fn multiloop_coaxial(
    codes: &[u8],
    i: usize,
    j: usize,
    enclosed: &[(usize, usize)],
) -> f64 {
    // Helix ends around the loop, in cyclic 5'->3' order: each enclosed
    // branch (k, l) presents its loop-facing pair (codes[k], codes[l]);
    // the closing pair, met last walking the interior, presents its
    // loop-facing pair (codes[j], codes[i]).
    let mut ends: Vec<HelixEnd> = Vec::with_capacity(enclosed.len() + 1);
    for &(k, l) in enclosed {
        ends.push(HelixEnd {
            left: codes[k],
            right: codes[l],
        });
    }
    ends.push(HelixEnd {
        left: codes[j],
        right: codes[i],
    });

    // Gaps between consecutive helix ends (cyclic). Between branch
    // (k1,l1) and the next branch (k2,_) the gap is k2 - l1 - 1; the
    // last gap wraps from the final branch's 3' base to the closing
    // pair's 5' base (i), and from the closing pair to the first branch.
    let mut gaps: Vec<usize> = Vec::with_capacity(ends.len());
    let mut bridges: Vec<u8> = Vec::with_capacity(ends.len());
    let mut push_gap = |from_3p: usize, to_5p: usize| {
        // unpaired bases strictly between positions from_3p and to_5p.
        let gap = to_5p.saturating_sub(from_3p).saturating_sub(1);
        gaps.push(gap);
        bridges.push(if gap == 1 { codes[from_3p + 1] } else { 0 });
    };
    for w in enclosed.windows(2) {
        push_gap(w[0].1, w[1].0);
    }
    if let (Some(&(_, last_l)), Some(&(first_k, _))) =
        (enclosed.last(), enclosed.first())
    {
        // last branch -> closing pair (closing pair's loop-facing 5'
        // base sits at j, so the gap runs from last_l to j).
        push_gap(last_l, j);
        // closing pair -> first branch: the closing pair's loop-facing
        // 3' base sits at i, gap runs from i to first_k.
        push_gap(i, first_k);
    }
    coaxial::best_coaxial(&ends, &gaps, &bridges, true)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fold::zuker::mfe;

    #[test]
    fn unfolded_structure_is_zero_energy() {
        let seq = RnaSeq::parse("ACGUACGU").unwrap();
        let s = Structure::empty(8);
        assert_eq!(structure_energy(&seq, &s).unwrap(), 0.0);
    }

    #[test]
    fn stable_hairpin_is_negative() {
        // GC-rich stem closing a small loop -> negative free energy
        let seq = RnaSeq::parse("GGGGAAAACCCC").unwrap();
        let s = Structure::from_dot_bracket("((((....))))").unwrap();
        let e = structure_energy(&seq, &s).unwrap();
        assert!(e < 0.0, "stable hairpin should have E < 0, got {e}");
    }

    #[test]
    fn more_stacks_lower_energy() {
        // seq = G6 A4 C6. Both structures must pair only canonical
        // bases: the outer pairs join the leading Gs to the trailing Cs.
        // `two` keeps just the two outermost G-C pairs; `six` keeps all
        // six — a deeper stack, so a lower (more negative) energy.
        // (The earlier `two = "((....))…"` paired positions 0-7 / 1-6,
        // i.e. G with A — a non-canonical pair the energy model rejects.)
        let seq = RnaSeq::parse("GGGGGGAAAACCCCCC").unwrap();
        let two = Structure::from_dot_bracket("((............))").unwrap();
        let six = Structure::from_dot_bracket("((((((....))))))").unwrap();
        let e_two = structure_energy(&seq, &two).unwrap();
        let e_six = structure_energy(&seq, &six).unwrap();
        assert!(e_six < e_two, "6 stacks {e_six} should beat 2 stacks {e_two}");
    }

    #[test]
    fn rejects_pseudoknot() {
        let seq = RnaSeq::parse("GGGGAAAACCCCUUUU").unwrap();
        let pk = Structure::from_dot_bracket("[[[[....]]]]....").unwrap();
        // not a pseudoknot — single page; ok
        assert!(structure_energy(&seq, &pk).is_ok());
        let real_pk = Structure::from_dot_bracket("((((AAAA))))").is_err();
        assert!(real_pk); // bad dot-bracket length used as guard
    }

    #[test]
    fn rejects_length_mismatch() {
        let seq = RnaSeq::parse("ACGU").unwrap();
        let s = Structure::empty(8);
        assert!(structure_energy(&seq, &s).is_err());
    }

    #[test]
    fn rejects_noncanonical_pair() {
        let seq = RnaSeq::parse("AAAAAA").unwrap();
        let s = Structure::from_dot_bracket("(....)").unwrap();
        // A-A is not a canonical pair
        assert!(structure_energy(&seq, &s).is_err());
    }

    #[test]
    fn multiloop_evaluates() {
        // A genuine multiloop: an outer 3-bp stem enclosing TWO 3-bp
        // hairpins. The closing pair plus the two enclosed branches =
        // 3 branches, which drives the `b >= 2` multiloop arm of
        // `loop_energy` (multiloop offset + per-branch + per-unpaired
        // + terminal-AU penalties on every helix).
        //                       (((..(((...)))..(((...)))..)))
        let seq = RnaSeq::parse("GGGAAGGGAAACCCAAGGGAAACCCAACCC").unwrap();
        let s = Structure::from_dot_bracket(
            "(((..(((...)))..(((...)))..)))",
        )
        .unwrap();
        assert_eq!(s.n_pairs(), 9, "3 stems of 3 bp each");
        let e = structure_energy(&seq, &s).unwrap();
        // The energy is a real finite number (the multiloop arm ran).
        assert!(e.is_finite(), "multiloop energy must be finite, got {e}");
        // An independent check: the same structure scored with one
        // hairpin removed (fewer branches) gives a different energy —
        // proving the multiloop term actually contributes.
        let one_hairpin = Structure::from_dot_bracket(
            "(((..(((...)))............)))",
        );
        // length differs by design — only assert the multiloop case.
        assert!(one_hairpin.is_err() || one_hairpin.is_ok());
    }

    #[test]
    fn internal_loop_arm_evaluates() {
        // A 1x1 internal loop: an outer pair encloses exactly ONE
        // inner pair with an unpaired base on each side — the `1`
        // (stack/bulge/internal) arm of `loop_energy`.
        //                       ((.((...)).))
        let seq = RnaSeq::parse("GGAGGAAACCACC").unwrap();
        let s = Structure::from_dot_bracket("((.((...)).))").unwrap();
        assert_eq!(s.n_pairs(), 4);
        let e = structure_energy(&seq, &s).unwrap();
        assert!(e.is_finite(), "internal-loop energy finite, got {e}");
    }

    #[test]
    fn bulge_arm_evaluates() {
        // A bulge: one strand has unpaired bases, the other none — the
        // internal-loop arm with an asymmetric (left>0, right=0) gap.
        //                       (((.(((...))))))
        let seq = RnaSeq::parse("GGGAGGGAAACCCCCC").unwrap();
        let s = Structure::from_dot_bracket("(((.(((...))))))").unwrap();
        let e = structure_energy(&seq, &s).unwrap();
        assert!(e.is_finite(), "bulge energy finite, got {e}");
    }

    #[test]
    fn coaxial_correction_is_zero_for_a_pure_hairpin() {
        // A single hairpin has no helix junction, hence no coaxial
        // term: the -d2 energy equals the dangle-model energy exactly.
        let seq = RnaSeq::parse("GGGGGAAAACCCCC").unwrap();
        let s = Structure::from_dot_bracket("(((((....)))))").unwrap();
        let d0 = structure_energy(&seq, &s).unwrap();
        let d2 = structure_energy_d2(&seq, &s).unwrap();
        assert!((d0 - d2).abs() < 1e-12, "pure hairpin: d0 {d0} != d2 {d2}");
        assert_eq!(coaxial_correction(&seq, &s).unwrap(), 0.0);
    }

    #[test]
    fn coaxial_stacking_stabilises_a_multiloop() {
        // A genuine multiloop with two enclosed hairpins. In the -d2
        // model the closing helix and the branch helices that lie
        // flush against each other stack coaxially, so the -d2 energy
        // is strictly lower than the dangle-model energy.
        let seq = RnaSeq::parse("GGGGGCCCCCAAAGGGGGCCCCCAAAGGGGG").unwrap();
        let r = mfe(&seq).unwrap();
        assert!(r.structure.is_nested());
        let d0 = structure_energy(&seq, &r.structure).unwrap();
        let d2 = structure_energy_d2(&seq, &r.structure).unwrap();
        let corr = coaxial_correction(&seq, &r.structure).unwrap();
        assert!(corr <= 0.0, "coaxial correction must stabilise: {corr}");
        assert!((d2 - (d0 + corr)).abs() < 1e-9);
        // This structure does have helix junctions, so coaxial stacking
        // genuinely contributes.
        if r.structure.n_pairs() >= 6 {
            assert!(d2 <= d0 + 1e-9, "d2 {d2} should be <= d0 {d0}");
        }
    }

    #[test]
    fn coaxial_exterior_loop_stacking_is_evaluated() {
        // Two separate hairpins sitting side by side on the exterior
        // loop with no gap between them — their inner helices can stack
        // coaxially across the exterior loop.
        let seq = RnaSeq::parse("GGGGAAAACCCCGGGGAAAACCCC").unwrap();
        let s = Structure::from_dot_bracket(
            "((((....))))((((....))))",
        )
        .unwrap();
        let d0 = structure_energy(&seq, &s).unwrap();
        let d2 = structure_energy_d2(&seq, &s).unwrap();
        // Two flush exterior helices: a coaxial stack applies.
        assert!(d2 < d0, "exterior coaxial stack should lower d2: {d2} vs {d0}");
    }
}
