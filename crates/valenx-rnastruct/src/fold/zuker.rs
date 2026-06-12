//! Zuker minimum-free-energy folding with the Turner-2004 model.
//!
//! This is the standard Zuker-Stiegler (1981) secondary-structure DP:
//! it finds the *pseudoknot-free* structure of minimum total free
//! energy under the nearest-neighbor model in [`crate::fold::energy`].
//!
//! ## Recurrences
//!
//! Four `O(n²)` matrices, filled by increasing sub-sequence length:
//!
//! - `v[i][j]` — minimum energy of `i..=j` **given that `i` pairs
//!   `j`**. The closing pair encloses a hairpin, an internal/bulge
//!   loop (one inner pair `(k,l)`), or a multiloop.
//! - `w[i][j]` — minimum energy of the sub-sequence `i..=j` with no
//!   constraint (the "exterior-loop within `i..=j`" matrix).
//! - `wm[i][j]` — minimum energy of `i..=j` as **part of a
//!   multiloop**, containing **>= 1** branch.
//! - `wm2[i][j]` — like `wm` but containing **>= 2** branches; used to
//!   split a multiloop interior.
//!
//! The internal-loop search is capped at [`MAX_LOOP`] unpaired bases
//! per side, the standard ViennaRNA optimisation that makes the
//! internal-loop case `O(n^2 * MAX_LOOP^2)` rather than `O(n^4)`. The
//! overall complexity is `O(n^3)`.
//!
//! Hard and soft constraints ([`crate::fold::constraint`]) are honored
//! throughout: a forbidden pair is never opened in `v`, a
//! forced-unpaired base contributes its soft pseudo-energy in `w` /
//! `wm`, and `mfe_constrained` finally rejects a structure missing a
//! forced pair.

use crate::error::{Result, RnaStructError};
use crate::fold::constraint::FoldConstraints;
use crate::fold::energy::{self, multiloop};
use crate::fold::nussinov::MIN_HAIRPIN;
use crate::rna::RnaSeq;
use crate::structure::Structure;

/// Maximum number of unpaired bases on either side of an internal
/// loop the DP will consider. Loops larger than this are extremely
/// rare in real structures; capping keeps folding `O(n^3)`.
pub const MAX_LOOP: usize = 30;

/// A "this sub-problem is infeasible" sentinel energy.
const INF: f64 = energy::FORBIDDEN;

/// Result of a minimum-free-energy fold.
#[derive(Clone, Debug, PartialEq)]
pub struct MfeResult {
    /// The minimum-free-energy structure.
    pub structure: Structure,
    /// Its total free energy in kcal/mol.
    pub energy: f64,
}

/// The filled DP matrices — kept so the ensemble suboptimal / sampling
/// code can reuse a Zuker fill. Square `n x n` flat buffers.
pub(crate) struct ZukerTables {
    pub n: usize,
    pub v: Vec<f64>,
    pub w: Vec<f64>,
    pub wm: Vec<f64>,
    pub wm2: Vec<f64>,
}

impl ZukerTables {
    #[inline]
    pub(crate) fn idx(&self, i: usize, j: usize) -> usize {
        i * self.n + j
    }
}

/// Folds `seq` to its minimum-free-energy structure.
pub fn mfe(seq: &RnaSeq) -> Result<MfeResult> {
    mfe_constrained(seq, &FoldConstraints::none(seq.len()))
}

/// Folds `seq` and reports its energy in ViennaRNA's **`-d2`** model —
/// the dangle-folded Zuker MFE *structure*, with its energy re-scored
/// to include the explicit [`crate::fold::coaxial`] coaxial-stacking
/// correction.
///
/// The returned [`MfeResult::structure`] is the optimum of the
/// dangle-folded model (the same one [`mfe`] returns); the
/// [`MfeResult::energy`] is [`crate::fold::eval::structure_energy_d2`]
/// of it, so it equals ViennaRNA `RNAeval -d2` exactly-to-rounding and
/// is self-consistent with `structure_energy_d2`.
///
/// Coaxial stacking only ever stabilises, so the reported energy is
/// ≤ the [`mfe`] energy. For the small minority of sequences where
/// coaxial stacking changes *which* structure is optimal, folding the
/// coaxial term into the recurrence itself (rather than re-scoring the
/// dangle-model optimum) is a documented follow-up — the v1 coaxial
/// model is exact for energy *evaluation* of any given structure.
///
/// # Errors
/// As [`mfe`].
pub fn mfe_d2(seq: &RnaSeq) -> Result<MfeResult> {
    let base = mfe(seq)?;
    let energy = crate::fold::eval::structure_energy_d2(seq, &base.structure)?;
    Ok(MfeResult {
        structure: base.structure,
        energy,
    })
}

/// Folds `seq` to its minimum-free-energy structure subject to
/// `constraints` (hard structural constraints and/or soft SHAPE-style
/// pseudo-energies).
///
/// # Errors
/// - [`RnaStructError::Invalid`] if the constraint length does not
///   match the sequence length.
/// - [`RnaStructError::Structure`] if a forced pair could not be
///   honored (the hard constraints are contradictory / infeasible).
pub fn mfe_constrained(seq: &RnaSeq, constraints: &FoldConstraints) -> Result<MfeResult> {
    let codes = seq.codes();
    let n = codes.len();
    if constraints.len() != n {
        return Err(RnaStructError::invalid(
            "constraints",
            format!(
                "constraint length {} != sequence length {n}",
                constraints.len()
            ),
        ));
    }
    if n == 0 {
        return Ok(MfeResult {
            structure: Structure::empty(0),
            energy: 0.0,
        });
    }

    let tables = fill(codes, constraints);
    let structure = traceback(codes, constraints, &tables)?;

    // Honor forced pairs: if any is missing the constraints were
    // infeasible.
    for &(i, j) in constraints.forced_pairs() {
        if structure.partner(i) != Some(j) {
            return Err(RnaStructError::structure(format!(
                "forced pair ({i}, {j}) could not be satisfied — \
                 the hard constraints are infeasible"
            )));
        }
    }

    let mut energy = tables.w[tables.idx(0, n - 1)];
    if energy >= INF / 2.0 {
        // No feasible nested structure under the constraints.
        return Err(RnaStructError::structure(
            "no feasible structure under the given constraints",
        ));
    }
    // Numerical floor cleanup.
    if energy.abs() < 1e-9 {
        energy = 0.0;
    }
    Ok(MfeResult { structure, energy })
}

/// Fills the four Zuker DP matrices. Pure function of the encoded
/// sequence and constraints.
pub(crate) fn fill(codes: &[u8], cons: &FoldConstraints) -> ZukerTables {
    let n = codes.len();
    let mut t = ZukerTables {
        n,
        v: vec![INF; n * n],
        w: vec![0.0; n * n],
        wm: vec![INF; n * n],
        wm2: vec![INF; n * n],
    };

    // Diagonal: a single unpaired base; w is its soft pseudo-energy.
    for i in 0..n {
        let idx = t.idx(i, i);
        t.w[idx] = if cons.unpaired_allowed(i) {
            cons.soft_unpaired(i)
        } else {
            INF
        };
    }

    // Fill by increasing span.
    for span in 1..n {
        for i in 0..(n - span) {
            let j = i + span;

            // ---- v[i][j]: i pairs j --------------------------------
            let mut vij = INF;
            if span > MIN_HAIRPIN
                && energy::can_pair_codes(codes[i], codes[j])
                && cons.pair_allowed(i, j)
            {
                // Hairpin.
                let loop_bases = &codes[(i + 1)..j];
                vij = vij.min(energy::hairpin_energy(codes[i], codes[j], loop_bases));

                // Internal / bulge / stack: inner pair (k, l).
                let k_max = (i + 1 + MAX_LOOP).min(j.saturating_sub(MIN_HAIRPIN + 1));
                for k in (i + 1)..=k_max {
                    let l_min = (k + MIN_HAIRPIN + 1).max(j.saturating_sub(MAX_LOOP + 1));
                    for l in l_min..j {
                        if l <= k {
                            continue;
                        }
                        let inner = t.v[t.idx(k, l)];
                        if inner >= INF / 2.0 {
                            continue;
                        }
                        let left = k - i - 1;
                        let right = j - l - 1;
                        if left + right != 0 && (left > MAX_LOOP || right > MAX_LOOP) {
                            continue;
                        }
                        let il = energy::internal_loop_energy(
                            codes[i],
                            codes[j],
                            codes[k],
                            codes[l],
                            left,
                            right,
                            codes[i + 1],
                            codes[j - 1],
                            codes[k - 1],
                            codes[l + 1],
                        );
                        vij = vij.min(il + inner);
                    }
                }

                // Multiloop closed by (i, j): interior is wm2 over
                // i+1..j-1 (>= 2 branches), plus closure cost.
                if j >= i + 2 {
                    let interior = t.wm2[t.idx(i + 1, j - 1)];
                    if interior < INF / 2.0 {
                        let closure = multiloop::OFFSET
                            + multiloop::PER_BRANCH
                            + energy::terminal_penalty(codes[i], codes[j]);
                        vij = vij.min(closure + interior);
                    }
                }
            }
            let vidx = t.idx(i, j);
            t.v[vidx] = vij;

            // ---- wm[i][j]: part of a multiloop, >= 1 branch --------
            let mut wmij = INF;
            // (i, j) itself is a branch helix.
            if vij < INF / 2.0 {
                wmij = wmij.min(vij + multiloop::PER_BRANCH);
            }
            // i unpaired, rest is wm.
            if cons.unpaired_allowed(i) {
                let rest = t.wm[t.idx(i + 1, j)];
                if rest < INF / 2.0 {
                    wmij = wmij.min(rest + multiloop::PER_UNPAIRED + cons.soft_unpaired(i));
                }
            }
            // j unpaired, rest is wm.
            if cons.unpaired_allowed(j) {
                let rest = t.wm[t.idx(i, j - 1)];
                if rest < INF / 2.0 {
                    wmij = wmij.min(rest + multiloop::PER_UNPAIRED + cons.soft_unpaired(j));
                }
            }
            // split into two non-empty multiloop pieces.
            for k in i..j {
                let a = t.wm[t.idx(i, k)];
                let b = t.wm[t.idx(k + 1, j)];
                if a < INF / 2.0 && b < INF / 2.0 {
                    wmij = wmij.min(a + b);
                }
            }
            let wm_idx = t.idx(i, j);
            t.wm[wm_idx] = wmij;

            // ---- wm2[i][j]: part of a multiloop, >= 2 branches -----
            let mut wm2ij = INF;
            // i unpaired, rest still >= 2 branches.
            if cons.unpaired_allowed(i) {
                let rest = t.wm2[t.idx(i + 1, j)];
                if rest < INF / 2.0 {
                    wm2ij = wm2ij.min(rest + multiloop::PER_UNPAIRED + cons.soft_unpaired(i));
                }
            }
            // j unpaired.
            if cons.unpaired_allowed(j) {
                let rest = t.wm2[t.idx(i, j - 1)];
                if rest < INF / 2.0 {
                    wm2ij = wm2ij.min(rest + multiloop::PER_UNPAIRED + cons.soft_unpaired(j));
                }
            }
            // split: >= 1 branch (wm) on the left, >= 1 branch (wm) on
            // the right -> together >= 2.
            for k in i..j {
                let a = t.wm[t.idx(i, k)];
                let b = t.wm[t.idx(k + 1, j)];
                if a < INF / 2.0 && b < INF / 2.0 {
                    wm2ij = wm2ij.min(a + b);
                }
            }
            let wm2_idx = t.idx(i, j);
            t.wm2[wm2_idx] = wm2ij;

            // ---- w[i][j]: exterior of i..=j (>= 0 branches) --------
            // all-unpaired baseline (only if every base may be
            // unpaired).
            let mut all_unpaired = 0.0;
            let mut feasible_unpaired = true;
            for p in i..=j {
                if cons.unpaired_allowed(p) {
                    all_unpaired += cons.soft_unpaired(p);
                } else {
                    feasible_unpaired = false;
                    break;
                }
            }
            let mut wij = if feasible_unpaired { all_unpaired } else { INF };
            // i unpaired.
            if cons.unpaired_allowed(i) {
                let rest = t.w[t.idx(i + 1, j)];
                if rest < INF / 2.0 {
                    wij = wij.min(rest + cons.soft_unpaired(i));
                }
            }
            // (i, j) is a closing pair of the whole window.
            if vij < INF / 2.0 {
                wij = wij.min(vij + energy::terminal_penalty(codes[i], codes[j]));
            }
            // split the exterior loop into two.
            for k in i..j {
                let a = t.w[t.idx(i, k)];
                let b = t.w[t.idx(k + 1, j)];
                if a < INF / 2.0 && b < INF / 2.0 {
                    wij = wij.min(a + b);
                }
            }
            let w_idx = t.idx(i, j);
            t.w[w_idx] = wij;
        }
    }
    t
}

/// Reconstructs the MFE structure from the filled tables by an
/// explicit-stack traceback over the same recurrences.
pub(crate) fn traceback(
    codes: &[u8],
    cons: &FoldConstraints,
    t: &ZukerTables,
) -> Result<Structure> {
    let n = t.n;
    let mut partner: Vec<Option<usize>> = vec![None; n];
    if n == 0 {
        return Structure::from_partner(partner);
    }

    /// What kind of sub-problem a traceback frame represents.
    #[derive(Copy, Clone)]
    enum Job {
        W(usize, usize),
        V(usize, usize),
        Wm(usize, usize),
        Wm2(usize, usize),
    }

    let feq = |a: f64, b: f64| (a - b).abs() < 1e-7;
    let mut stack: Vec<Job> = vec![Job::W(0, n - 1)];

    while let Some(job) = stack.pop() {
        match job {
            Job::W(i, j) => {
                if i >= j {
                    continue; // single base or empty: unpaired
                }
                let target = t.w[t.idx(i, j)];
                // i unpaired
                if cons.unpaired_allowed(i) {
                    let rest = t.w[t.idx(i + 1, j)];
                    if rest < INF / 2.0 && feq(target, rest + cons.soft_unpaired(i)) {
                        stack.push(Job::W(i + 1, j));
                        continue;
                    }
                }
                // (i, j) closes
                let vij = t.v[t.idx(i, j)];
                if vij < INF / 2.0
                    && feq(target, vij + energy::terminal_penalty(codes[i], codes[j]))
                {
                    stack.push(Job::V(i, j));
                    continue;
                }
                // split
                let mut split = false;
                for k in i..j {
                    let a = t.w[t.idx(i, k)];
                    let b = t.w[t.idx(k + 1, j)];
                    if a < INF / 2.0 && b < INF / 2.0 && feq(target, a + b) {
                        stack.push(Job::W(i, k));
                        stack.push(Job::W(k + 1, j));
                        split = true;
                        break;
                    }
                }
                if split {
                    continue;
                }
                // else: all-unpaired window — nothing to do.
            }

            Job::V(i, j) => {
                // (i, j) is a pair.
                partner[i] = Some(j);
                partner[j] = Some(i);
                let target = t.v[t.idx(i, j)];

                // hairpin?
                let loop_bases = &codes[(i + 1)..j];
                if feq(
                    target,
                    energy::hairpin_energy(codes[i], codes[j], loop_bases),
                ) {
                    continue;
                }
                // internal / bulge / stack?
                let mut done = false;
                let k_max = (i + 1 + MAX_LOOP).min(j.saturating_sub(MIN_HAIRPIN + 1));
                'inner: for k in (i + 1)..=k_max {
                    let l_min = (k + MIN_HAIRPIN + 1).max(j.saturating_sub(MAX_LOOP + 1));
                    for l in l_min..j {
                        if l <= k {
                            continue;
                        }
                        let inner = t.v[t.idx(k, l)];
                        if inner >= INF / 2.0 {
                            continue;
                        }
                        let left = k - i - 1;
                        let right = j - l - 1;
                        if left + right != 0 && (left > MAX_LOOP || right > MAX_LOOP) {
                            continue;
                        }
                        let il = energy::internal_loop_energy(
                            codes[i],
                            codes[j],
                            codes[k],
                            codes[l],
                            left,
                            right,
                            codes[i + 1],
                            codes[j - 1],
                            codes[k - 1],
                            codes[l + 1],
                        );
                        if feq(target, il + inner) {
                            stack.push(Job::V(k, l));
                            done = true;
                            break 'inner;
                        }
                    }
                }
                if done {
                    continue;
                }
                // multiloop.
                if j >= i + 2 {
                    let interior = t.wm2[t.idx(i + 1, j - 1)];
                    if interior < INF / 2.0 {
                        let closure = multiloop::OFFSET
                            + multiloop::PER_BRANCH
                            + energy::terminal_penalty(codes[i], codes[j]);
                        if feq(target, closure + interior) {
                            stack.push(Job::Wm2(i + 1, j - 1));
                            continue;
                        }
                    }
                }
            }

            Job::Wm(i, j) => {
                if i > j {
                    continue;
                }
                let target = t.wm[t.idx(i, j)];
                let vij = t.v[t.idx(i, j)];
                if vij < INF / 2.0 && feq(target, vij + multiloop::PER_BRANCH) {
                    stack.push(Job::V(i, j));
                    continue;
                }
                if cons.unpaired_allowed(i) && i < j {
                    let rest = t.wm[t.idx(i + 1, j)];
                    if rest < INF / 2.0
                        && feq(
                            target,
                            rest + multiloop::PER_UNPAIRED + cons.soft_unpaired(i),
                        )
                    {
                        stack.push(Job::Wm(i + 1, j));
                        continue;
                    }
                }
                if cons.unpaired_allowed(j) && i < j {
                    let rest = t.wm[t.idx(i, j - 1)];
                    if rest < INF / 2.0
                        && feq(
                            target,
                            rest + multiloop::PER_UNPAIRED + cons.soft_unpaired(j),
                        )
                    {
                        stack.push(Job::Wm(i, j - 1));
                        continue;
                    }
                }
                for k in i..j {
                    let a = t.wm[t.idx(i, k)];
                    let b = t.wm[t.idx(k + 1, j)];
                    if a < INF / 2.0 && b < INF / 2.0 && feq(target, a + b) {
                        stack.push(Job::Wm(i, k));
                        stack.push(Job::Wm(k + 1, j));
                        break;
                    }
                }
            }

            Job::Wm2(i, j) => {
                if i >= j {
                    continue;
                }
                let target = t.wm2[t.idx(i, j)];
                if cons.unpaired_allowed(i) {
                    let rest = t.wm2[t.idx(i + 1, j)];
                    if rest < INF / 2.0
                        && feq(
                            target,
                            rest + multiloop::PER_UNPAIRED + cons.soft_unpaired(i),
                        )
                    {
                        stack.push(Job::Wm2(i + 1, j));
                        continue;
                    }
                }
                if cons.unpaired_allowed(j) {
                    let rest = t.wm2[t.idx(i, j - 1)];
                    if rest < INF / 2.0
                        && feq(
                            target,
                            rest + multiloop::PER_UNPAIRED + cons.soft_unpaired(j),
                        )
                    {
                        stack.push(Job::Wm2(i, j - 1));
                        continue;
                    }
                }
                for k in i..j {
                    let a = t.wm[t.idx(i, k)];
                    let b = t.wm[t.idx(k + 1, j)];
                    if a < INF / 2.0 && b < INF / 2.0 && feq(target, a + b) {
                        stack.push(Job::Wm(i, k));
                        stack.push(Job::Wm(k + 1, j));
                        break;
                    }
                }
            }
        }
    }

    Structure::from_partner(partner)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fold::eval::structure_energy;

    #[test]
    fn empty_and_tiny() {
        let r = mfe(&RnaSeq::parse("A").unwrap()).unwrap();
        assert_eq!(r.structure.n_pairs(), 0);
        assert_eq!(r.energy, 0.0);
    }

    #[test]
    fn folds_a_clean_hairpin() {
        // A strong GC stem closing a 4-nt loop should fold up.
        let seq = RnaSeq::parse("GGGGGAAAACCCCC").unwrap();
        let r = mfe(&seq).unwrap();
        assert!(r.structure.n_pairs() >= 3, "expected a stem, got {r:?}");
        assert!(
            r.energy < 0.0,
            "stable hairpin should have E<0: {}",
            r.energy
        );
        assert!(r.structure.is_nested());
    }

    #[test]
    fn mfe_energy_matches_independent_eval() {
        // The energy the DP reports must equal a fresh evaluation of
        // the structure it returns.
        let seq = RnaSeq::parse("GGGAAACCCAAAGGGAAACCC").unwrap();
        let r = mfe(&seq).unwrap();
        let recomputed = structure_energy(&seq, &r.structure).unwrap();
        assert!(
            (recomputed - r.energy).abs() < 1e-4,
            "DP energy {} != eval {}",
            r.energy,
            recomputed
        );
    }

    #[test]
    fn unpairable_sequence_stays_open() {
        let seq = RnaSeq::parse("AAAAAAAAAAAA").unwrap();
        let r = mfe(&seq).unwrap();
        assert_eq!(r.structure.n_pairs(), 0);
        assert_eq!(r.energy, 0.0);
    }

    #[test]
    fn forbid_pair_changes_fold() {
        let seq = RnaSeq::parse("GGGGGAAAACCCCC").unwrap();
        let free = mfe(&seq).unwrap();
        assert_eq!(free.structure.partner(0), Some(13));
        // forbid the outer pair
        let mut cons = FoldConstraints::none(seq.len());
        cons.forbid_pair(0, 13).unwrap();
        let constrained = mfe_constrained(&seq, &cons).unwrap();
        assert_ne!(constrained.structure.partner(0), Some(13));
    }

    #[test]
    fn force_unpaired_is_respected() {
        let seq = RnaSeq::parse("GGGGGAAAACCCCC").unwrap();
        let mut cons = FoldConstraints::none(seq.len());
        cons.force_unpaired(2).unwrap();
        let r = mfe_constrained(&seq, &cons).unwrap();
        assert!(!r.structure.is_paired(2), "position 2 was forced unpaired");
    }

    #[test]
    fn force_pair_appears_in_result() {
        let seq = RnaSeq::parse("GGGGGAAAAACCCCC").unwrap();
        let mut cons = FoldConstraints::none(seq.len());
        cons.force_pair(1, 13).unwrap();
        let r = mfe_constrained(&seq, &cons).unwrap();
        assert_eq!(r.structure.partner(1), Some(13));
    }

    #[test]
    fn folds_a_multiloop() {
        // two hairpins under one closing helix
        let seq = RnaSeq::parse("GGGGGCCCCCAAAGGGGGCCCCCAAAGGGGG").unwrap();
        let r = mfe(&seq).unwrap();
        assert!(r.structure.is_nested());
        // recomputed energy agrees
        let e = structure_energy(&seq, &r.structure).unwrap();
        assert!((e - r.energy).abs() < 1e-4);
    }
}
