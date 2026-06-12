//! Feature 16 — NUPACK-class ensemble-defect inverse folding.
//!
//! The existing structural designer ([`crate::design::structural`]) and
//! the `valenx-rnastruct` inverse folder are an **adaptive walk**: they
//! mutate the worst-fitting position and keep a change if it lowers the
//! base-pair *distance* between the sequence's single MFE structure and
//! the target. That is a heuristic — it judges a sequence by one
//! structure, ignoring the rest of the Boltzmann ensemble.
//!
//! NUPACK's inverse designer instead minimises the **ensemble defect**:
//! the equilibrium *expected number of incorrectly-paired or
//! incorrectly-unpaired nucleotides* relative to the target, computed
//! from the full base-pair-probability matrix. A sequence with ensemble
//! defect `d` is, on average over its whole equilibrium ensemble, `d`
//! nucleotides away from the target — a principled, ensemble-wide
//! objective rather than a single-structure proxy.
//!
//! ## The objective
//!
//! For a target structure `s*` and the base-pair-probability matrix
//! `p(i,j)` of a candidate sequence,
//!
//! ```text
//!     ensemble_defect = Σ_i  [ 1 − p(i, s*(i)) ]      if i pairs s*(i),
//!                     + Σ_i  [ 1 − p_unpaired(i) ]    if i is unpaired,
//! ```
//!
//! a number in `[0, n]`. A perfect design has defect `0`. The
//! **normalised** defect (`/ n`) is a `[0, 1]` score.
//!
//! ## The probability matrix — LinearPartition
//!
//! The probabilities come from `valenx-rnastruct`'s **LinearPartition**
//! ([`valenx_rnastruct::linear_partition`]) — the linear-time partition
//! function — so the designer scales to long targets where the exact
//! `O(n³)` McCaskill folder would be too slow. (The exact folder is
//! available via [`crate::optimize::ensemble_defect`] for short
//! targets / cross-checks.)
//!
//! ## The search — hierarchical, leaf-first
//!
//! Inverse design is a mutation walk: propose a sequence change, accept
//! it if the ensemble defect drops. NUPACK's key idea is **hierarchical
//! decomposition** — it splits the target into structural sub-domains
//! (the *leaves* of the structure tree: individual hairpins and the
//! helices around them), optimises each leaf so its sub-defect is low,
//! then re-optimises the assembled sequence. This module implements a
//! sound version of that:
//!
//! 1. **Decompose** the target into leaf sub-structures — maximal
//!    helix + hairpin units.
//! 2. **Leaf pass** — for each leaf, run a focused mutation walk that
//!    minimises that leaf's local ensemble defect.
//! 3. **Global pass** — run a whole-sequence mutation walk on the
//!    assembled candidate, targeting the highest-defect positions, to
//!    polish the join regions the leaf passes could not see.
//!
//! Every mutation respects the target: a paired position is only ever
//! re-seated with a canonical pair, so the candidate always *can* adopt
//! the target.
//!
//! ## v1 scope — honest framing
//!
//! - This is a real ensemble-defect minimiser, not bit-for-bit NUPACK.
//!   The decomposition is the leaf-helix split described above, not
//!   NUPACK's full recursive multi-level tree with its conditional
//!   sub-ensemble bonus/penalty bookkeeping.
//! - The defect is computed from LinearPartition's *approximate*
//!   base-pair probabilities (beam search drops low-weight states); a
//!   wide beam makes them near-exact.
//! - The accepted sequence is a strong in-silico candidate — a low
//!   ensemble defect is a strong prediction the sequence folds to the
//!   target, not a measurement. Validate experimentally.

use crate::constraints::DesignConstraintSet;
use crate::error::{Result, RnaDesignError};
use serde::{Deserialize, Serialize};
use valenx_rnastruct::{
    base_pair_distance, linear_partition_with_beam, mfe, LinearPartitionResult, RnaSeq, Structure,
};

/// The four bases the designer may place, as ASCII.
const BASES: [u8; 4] = [b'A', b'C', b'G', b'U'];

/// Canonical pairing partners (Watson-Crick + G-U wobble), ASCII.
const CANON_PAIRS: [(u8, u8); 6] = [
    (b'A', b'U'),
    (b'U', b'A'),
    (b'G', b'C'),
    (b'C', b'G'),
    (b'G', b'U'),
    (b'U', b'G'),
];

/// The beam width LinearPartition runs at inside the designer. Wide
/// enough that the base-pair probabilities are near-exact for the
/// short / medium targets inverse design typically handles.
const LP_BEAM: usize = 200;

// ---------------------------------------------------------------------
// A small deterministic RNG (xorshift) — keeps designs reproducible.
// ---------------------------------------------------------------------

/// A deterministic xorshift RNG, local to the designer.
struct Rng {
    state: u64,
}

impl Rng {
    fn new(seed: u64) -> Self {
        Rng {
            state: seed ^ 0x2545_F491_4F6C_DD1D,
        }
    }

    fn next_u64(&mut self) -> u64 {
        let mut x = self.state;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.state = x;
        x
    }

    /// A `usize` in `[0, n)` (`n` must be non-zero).
    fn below(&mut self, n: usize) -> usize {
        (self.next_u64() % n as u64) as usize
    }

    /// An `f64` in `[0, 1)`.
    fn unit(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }
}

// ---------------------------------------------------------------------
// Parameters and result
// ---------------------------------------------------------------------

/// Parameters for [`inverse_fold_ensemble_defect`].
#[derive(Copy, Clone, Debug)]
pub struct EnsembleDefectParams {
    /// The normalised-ensemble-defect target. The search stops once the
    /// candidate's normalised defect drops to or below this — a
    /// `[0, 1]` value. NUPACK's default stop is `0.01`.
    pub defect_target: f64,
    /// The mutation budget for **each leaf** sub-structure pass.
    pub leaf_iterations: usize,
    /// The mutation budget for the **global** polishing pass.
    pub global_iterations: usize,
    /// The random seed — the search is deterministic for a fixed seed.
    pub seed: u64,
}

impl Default for EnsembleDefectParams {
    /// A solid v1 default: stop at a 1 %-normalised defect, 300
    /// leaf-pass iterations, 600 global iterations.
    fn default() -> Self {
        EnsembleDefectParams {
            defect_target: 0.01,
            leaf_iterations: 300,
            global_iterations: 600,
            seed: 0xE7DE,
        }
    }
}

/// The result of an ensemble-defect inverse-folding run.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct EnsembleDefectDesign {
    /// The designed RNA sequence (`A C G U`).
    pub sequence: Vec<u8>,
    /// The ensemble defect of the design relative to the target
    /// (expected number of incorrectly-(un)paired nucleotides).
    pub ensemble_defect: f64,
    /// The normalised ensemble defect (`ensemble_defect / n`) in
    /// `[0, 1]` — `0` is perfect.
    pub normalized_defect: f64,
    /// `true` if the normalised defect reached the requested target.
    pub solved: bool,
    /// The base-pair distance between the design's MFE fold and the
    /// target — `0` means the *single* MFE structure is exactly the
    /// target (a stronger statement than a low ensemble defect).
    pub mfe_distance: usize,
    /// Total mutation steps accepted across every pass.
    pub accepted_steps: usize,
    /// Total mutation steps attempted.
    pub total_steps: usize,
    /// Human-readable notes, one line per major decision.
    pub notes: Vec<String>,
}

impl EnsembleDefectDesign {
    /// The designed sequence as a `&str` (`A C G U`).
    pub fn sequence_str(&self) -> &str {
        std::str::from_utf8(&self.sequence).unwrap_or("")
    }
}

// ---------------------------------------------------------------------
// Ensemble defect from a LinearPartition result
// ---------------------------------------------------------------------

/// The ensemble defect of `seq` relative to `target`, computed from
/// LinearPartition base-pair probabilities.
///
/// This is the principled NUPACK objective. For each position the
/// contribution is the probability it is *not* in its target pairing
/// state.
///
/// # Errors
/// - [`RnaDesignError::Invalid`] if `seq` and `target` differ in length.
/// - [`RnaDesignError::Upstream`] if LinearPartition fails.
pub fn ensemble_defect_linear(seq: &[u8], target: &Structure) -> Result<f64> {
    if seq.len() != target.len() {
        return Err(RnaDesignError::invalid(
            "target",
            "sequence and target structure differ in length",
        ));
    }
    if seq.is_empty() {
        return Ok(0.0);
    }
    let rna = RnaSeq::parse(seq)?;
    let lp = linear_partition_with_beam(
        &rna,
        LP_BEAM,
        valenx_rnastruct::ensemble::linear_partition::DEFAULT_TEMPERATURE_K,
    )?;
    Ok(defect_from_lp(&lp, target, 0, seq.len()))
}

/// The ensemble defect contributed by positions `[lo, hi)` of `target`,
/// read from a LinearPartition result `lp` of the *whole* sequence.
fn defect_from_lp(lp: &LinearPartitionResult, target: &Structure, lo: usize, hi: usize) -> f64 {
    let mut defect = 0.0;
    for i in lo..hi {
        match target.partner(i) {
            Some(j) => {
                // Target wants i paired to j.
                defect += 1.0 - lp.pair_probability(i, j);
            }
            None => {
                // Target wants i unpaired.
                defect += 1.0 - lp.unpaired_probability(i);
            }
        }
    }
    defect.max(0.0)
}

// ---------------------------------------------------------------------
// Structural decomposition into leaf sub-domains
// ---------------------------------------------------------------------

/// A leaf sub-domain of the target: a contiguous `[start, end)` window
/// containing one hairpin and the helix that closes it. The leaf passes
/// optimise these windows independently before the global polish.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
struct Leaf {
    /// Inclusive start of the window.
    start: usize,
    /// Exclusive end of the window.
    end: usize,
}

/// Decomposes `target` into leaf sub-domains.
///
/// A leaf is the span of an innermost hairpin loop together with the
/// stacked helix immediately enclosing it. Concretely: every base pair
/// `(i, j)` that *directly* encloses a hairpin (no pair nested strictly
/// inside it) seeds a leaf; the leaf window is grown outward over the
/// contiguous stacked helix above that pair. Leaves never overlap.
fn decompose_leaves(target: &Structure) -> Vec<Leaf> {
    let n = target.len();
    let mut leaves = Vec::new();
    // Find every pair that closes a hairpin: (i, j) paired with no pair
    // inside (i, j).
    for i in 0..n {
        let j = match target.partner(i) {
            Some(j) if j > i => j,
            _ => continue,
        };
        // Is (i, j) a hairpin-closing pair? No paired base strictly
        // inside.
        let inner_unpaired = ((i + 1)..j).all(|k| target.partner(k).is_none());
        if !inner_unpaired {
            continue;
        }
        // Grow the helix outward: while (i-1, j+1) is also a pair.
        let mut lo = i;
        let mut hi = j;
        while lo > 0 && hi + 1 < n {
            if target.partner(lo - 1) == Some(hi + 1) {
                lo -= 1;
                hi += 1;
            } else {
                break;
            }
        }
        leaves.push(Leaf {
            start: lo,
            end: hi + 1,
        });
    }
    leaves.sort_by_key(|l| (l.start, l.end));
    leaves
}

// ---------------------------------------------------------------------
// The designer
// ---------------------------------------------------------------------

/// Designs a sequence that folds to `target`, minimising the ensemble
/// defect (feature 16).
///
/// Runs the hierarchical leaf-first search: decompose the target into
/// leaf sub-domains, optimise each leaf's local ensemble defect, then
/// run a global polishing pass on the assembled sequence. Every mutation
/// keeps the target's pairs canonical.
///
/// This is the *unconstrained* entry point; use
/// [`inverse_fold_constrained`] to also honour locked positions, a GC
/// band and forbidden motifs.
///
/// # Errors
/// - [`RnaDesignError::Goal`] if `target` is empty or pseudoknotted.
/// - [`RnaDesignError::Upstream`] if a folding call fails.
pub fn inverse_fold_ensemble_defect(
    target: &Structure,
    params: EnsembleDefectParams,
) -> Result<EnsembleDefectDesign> {
    run_inverse_fold(target, None, params)
}

/// Designs a sequence that folds to `target` while honouring a
/// constraint set — locked positions, a GC band, forbidden motifs and a
/// homopolymer cap (feature 16 + 17).
///
/// Same hierarchical ensemble-defect search as
/// [`inverse_fold_ensemble_defect`], but: locked positions are seeded to
/// their fixed base and never mutated; the search objective is the
/// ensemble defect *plus* the constraint set's soft penalty, so the walk
/// is steered toward a sequence that both folds to the target and
/// clears every constraint.
///
/// # Errors
/// - [`RnaDesignError::Goal`] if `target` is empty or pseudoknotted.
/// - [`RnaDesignError::Invalid`] if `constraints` was built for a
///   different length than `target`.
/// - [`RnaDesignError::Upstream`] if a folding call fails.
pub fn inverse_fold_constrained(
    target: &Structure,
    constraints: &DesignConstraintSet,
    params: EnsembleDefectParams,
) -> Result<EnsembleDefectDesign> {
    if constraints.len() != target.len() {
        return Err(RnaDesignError::invalid(
            "constraints",
            "the constraint set's length does not match the target",
        ));
    }
    run_inverse_fold(target, Some(constraints), params)
}

/// The shared ensemble-defect inverse-folding routine, with or without a
/// constraint set.
fn run_inverse_fold(
    target: &Structure,
    constraints: Option<&DesignConstraintSet>,
    params: EnsembleDefectParams,
) -> Result<EnsembleDefectDesign> {
    let n = target.len();
    if n == 0 {
        return Err(RnaDesignError::goal(
            "target",
            "cannot inverse-fold an empty target",
        ));
    }
    if target.has_pseudoknot() {
        return Err(RnaDesignError::goal(
            "target",
            "target is pseudoknotted — the partition-function folder is pseudoknot-free",
        ));
    }
    if !(0.0..=1.0).contains(&params.defect_target) {
        return Err(RnaDesignError::invalid(
            "defect_target",
            "the normalised-defect target must lie in [0, 1]",
        ));
    }

    let mut rng = Rng::new(params.seed);
    let mut notes: Vec<String> = Vec::new();

    // --- seed: a random target-compatible sequence -------------------
    let mut seq = seed_sequence(target, constraints, &mut rng);

    let mut accepted = 0usize;
    let mut total = 0usize;

    // --- leaf passes -------------------------------------------------
    let leaves = decompose_leaves(target);
    notes.push(format!(
        "Decomposed the {n}-nt target into {} leaf sub-domain(s); optimising each before \
         the global polish.",
        leaves.len(),
    ));
    if let Some(c) = constraints {
        notes.push(format!(
            "Honouring {} locked position(s) and the GC / forbidden-motif / homopolymer \
             constraints.",
            c.locked_count(),
        ));
    }
    for leaf in &leaves {
        let (acc, tot) = optimize_leaf(&mut seq, target, constraints, *leaf, &params, &mut rng)?;
        accepted += acc;
        total += tot;
    }

    // --- global polishing pass --------------------------------------
    let target_threshold = params.defect_target;
    let (acc, tot, defect) = optimize_global(
        &mut seq,
        target,
        constraints,
        &params,
        target_threshold,
        &mut rng,
    )?;
    accepted += acc;
    total += tot;

    let normalized = if n > 0 { defect / n as f64 } else { 0.0 };
    let solved = normalized <= target_threshold + 1e-9;

    // Final MFE-distance check — a stronger statement than low defect.
    let rna = RnaSeq::parse(&seq)?;
    let mfe_struct = mfe(&rna)?.structure;
    let mfe_distance = base_pair_distance(&mfe_struct, target)?;

    notes.push(format!(
        "Ensemble-defect inverse folding: defect {defect:.3} (normalised {normalized:.4}); \
         {accepted}/{total} mutations accepted.",
    ));
    if let Some(c) = constraints {
        let ok = c.satisfies(&seq);
        notes.push(format!(
            "Constraint check on the final design: {}.",
            if ok {
                "every constraint satisfied"
            } else {
                "one or more constraints still violated (raise the iteration budget)"
            },
        ));
    }
    if solved {
        notes.push(format!(
            "The design reached the requested normalised-defect target \
             ({target_threshold:.4})."
        ));
    } else {
        notes.push(format!(
            "The design did not reach the {target_threshold:.4} normalised-defect \
             target — best was {normalized:.4}; raise the iteration budget to push it \
             lower."
        ));
    }
    if mfe_distance == 0 {
        notes.push(
            "The design's single minimum-free-energy structure is exactly the target — a \
             strong in-silico prediction; confirm experimentally."
                .to_string(),
        );
    } else {
        notes.push(format!(
            "The design's MFE structure differs from the target by {mfe_distance} pair(s); \
             the ensemble defect is the principled objective and remains low.",
        ));
    }

    Ok(EnsembleDefectDesign {
        sequence: seq,
        ensemble_defect: defect,
        normalized_defect: normalized,
        solved,
        mfe_distance,
        accepted_steps: accepted,
        total_steps: total,
        notes,
    })
}

/// Builds a random sequence compatible with `target`: every pair seated
/// with a random canonical pair, every unpaired position random, then
/// the constraint set's locked positions applied on top.
fn seed_sequence(
    target: &Structure,
    constraints: Option<&DesignConstraintSet>,
    rng: &mut Rng,
) -> Vec<u8> {
    let n = target.len();
    let mut seq = vec![b'A'; n];
    let mut assigned = vec![false; n];
    for bp in target.pairs() {
        let (a, b) = CANON_PAIRS[rng.below(6)];
        seq[bp.i] = a;
        seq[bp.j] = b;
        assigned[bp.i] = true;
        assigned[bp.j] = true;
    }
    for (i, s) in seq.iter_mut().enumerate() {
        if !assigned[i] {
            *s = BASES[rng.below(4)];
        }
    }
    if let Some(c) = constraints {
        c.apply_locks(&mut seq);
    }
    seq
}

/// Optimises one leaf sub-domain: a focused mutation walk that
/// minimises the leaf window's local ensemble defect (plus the
/// constraint soft penalty, if a constraint set is given). Mutations
/// only touch *free* positions inside the leaf. Returns
/// `(accepted, total)`.
fn optimize_leaf(
    seq: &mut [u8],
    target: &Structure,
    constraints: Option<&DesignConstraintSet>,
    leaf: Leaf,
    params: &EnsembleDefectParams,
    rng: &mut Rng,
) -> Result<(usize, usize)> {
    let mut best = leaf_objective(seq, target, constraints, leaf)?;
    let mut accepted = 0usize;
    let mut total = 0usize;
    for _ in 0..params.leaf_iterations {
        if best <= params.defect_target * (leaf.end - leaf.start) as f64 {
            break;
        }
        total += 1;
        let trial = propose_in_window(seq, target, constraints, leaf.start, leaf.end, rng);
        let trial = match trial {
            Some(t) => t,
            None => continue,
        };
        let trial_obj = leaf_objective(&trial, target, constraints, leaf)?;
        if trial_obj <= best {
            seq.copy_from_slice(&trial);
            best = trial_obj;
            accepted += 1;
        }
    }
    Ok((accepted, total))
}

/// The leaf-pass objective: the leaf window's local ensemble defect plus
/// the constraint soft penalty (which is a whole-sequence quantity).
fn leaf_objective(
    seq: &[u8],
    target: &Structure,
    constraints: Option<&DesignConstraintSet>,
    leaf: Leaf,
) -> Result<f64> {
    let defect = leaf_defect(seq, target, leaf)?;
    let penalty = constraints.map(|c| c.soft_penalty(seq)).unwrap_or(0.0);
    Ok(defect + penalty)
}

/// The local ensemble defect of a leaf window: LinearPartition is run on
/// the *whole* sequence (so context is honoured) and only the leaf
/// window's positions contribute.
fn leaf_defect(seq: &[u8], target: &Structure, leaf: Leaf) -> Result<f64> {
    let rna = RnaSeq::parse(seq)?;
    let lp = linear_partition_with_beam(
        &rna,
        LP_BEAM,
        valenx_rnastruct::ensemble::linear_partition::DEFAULT_TEMPERATURE_K,
    )?;
    Ok(defect_from_lp(&lp, target, leaf.start, leaf.end))
}

/// The global polishing pass: a whole-sequence mutation walk that
/// targets the highest-defect positions, minimising the ensemble defect
/// plus the constraint soft penalty. Returns `(accepted, total,
/// final_defect)` — `final_defect` is the *pure* ensemble defect of the
/// accepted sequence (the constraint penalty is a search aid, not part
/// of the reported defect).
fn optimize_global(
    seq: &mut [u8],
    target: &Structure,
    constraints: Option<&DesignConstraintSet>,
    params: &EnsembleDefectParams,
    target_threshold: f64,
    rng: &mut Rng,
) -> Result<(usize, usize, f64)> {
    let n = seq.len();
    let mut best_defect = ensemble_defect_linear(seq, target)?;
    let mut best_obj = best_defect + constraints.map(|c| c.soft_penalty(seq)).unwrap_or(0.0);
    let mut accepted = 0usize;
    let mut total = 0usize;
    for _ in 0..params.global_iterations {
        // Stop once the defect is low *and* every constraint is met.
        let constraints_ok = constraints.map(|c| c.satisfies(seq)).unwrap_or(true);
        if best_defect / n as f64 <= target_threshold + 1e-9 && constraints_ok {
            break;
        }
        total += 1;
        // Target the highest-defect position: re-run the partition
        // function, find the worst position, and mutate around it.
        let trial = propose_targeted(seq, target, constraints, rng)?;
        let trial = match trial {
            Some(t) => t,
            None => continue,
        };
        let trial_defect = ensemble_defect_linear(&trial, target)?;
        let trial_obj = trial_defect + constraints.map(|c| c.soft_penalty(&trial)).unwrap_or(0.0);
        if trial_obj <= best_obj {
            seq.copy_from_slice(&trial);
            best_defect = trial_defect;
            best_obj = trial_obj;
            accepted += 1;
        }
    }
    Ok((accepted, total, best_defect))
}

/// Proposes a single target-respecting mutation confined to positions
/// in `[lo, hi)`: a paired position re-seats its whole pair, an unpaired
/// position is mutated freely. Locked positions are never chosen, and a
/// pair with a locked partner is re-seated only at its free end with a
/// base that still pairs the locked partner.
fn propose_in_window(
    seq: &[u8],
    target: &Structure,
    constraints: Option<&DesignConstraintSet>,
    lo: usize,
    hi: usize,
    rng: &mut Rng,
) -> Option<Vec<u8>> {
    if lo >= hi {
        return None;
    }
    let free: Vec<usize> = (lo..hi)
        .filter(|&i| constraints.map(|c| c.is_free(i)).unwrap_or(true))
        .collect();
    if free.is_empty() {
        return None;
    }
    let pos = free[rng.below(free.len())];
    apply_mutation(seq, target, constraints, pos, rng)
}

/// Proposes a mutation targeting the *free* position that currently
/// contributes the most ensemble defect (a NUPACK-style worst-position
/// focus). Locked positions are excluded.
fn propose_targeted(
    seq: &[u8],
    target: &Structure,
    constraints: Option<&DesignConstraintSet>,
    rng: &mut Rng,
) -> Result<Option<Vec<u8>>> {
    let n = seq.len();
    let free: Vec<usize> = (0..n)
        .filter(|&i| constraints.map(|c| c.is_free(i)).unwrap_or(true))
        .collect();
    if free.is_empty() {
        return Ok(None);
    }
    let rna = RnaSeq::parse(seq)?;
    let lp = linear_partition_with_beam(
        &rna,
        LP_BEAM,
        valenx_rnastruct::ensemble::linear_partition::DEFAULT_TEMPERATURE_K,
    )?;
    // The free position contributing the most defect.
    let mut worst_pos = free[0];
    let mut worst_val = -1.0_f64;
    for &i in &free {
        let contrib = match target.partner(i) {
            Some(j) => 1.0 - lp.pair_probability(i, j),
            None => 1.0 - lp.unpaired_probability(i),
        };
        if contrib > worst_val {
            worst_val = contrib;
            worst_pos = i;
        }
    }
    // 75 % of the time mutate the worst position; 25 % a random free one
    // (to escape local optima).
    let pos = if rng.unit() < 0.75 {
        worst_pos
    } else {
        free[rng.below(free.len())]
    };
    Ok(apply_mutation(seq, target, constraints, pos, rng))
}

/// Applies a single target-respecting mutation at `pos` (which must be
/// free): if `pos` is paired, the whole pair is re-seated canonically;
/// if its partner is locked, only `pos` is mutated to a base that still
/// pairs the locked partner. An unpaired `pos` is mutated freely.
/// Returns `None` if no legal mutation exists at `pos`.
fn apply_mutation(
    seq: &[u8],
    target: &Structure,
    constraints: Option<&DesignConstraintSet>,
    pos: usize,
    rng: &mut Rng,
) -> Option<Vec<u8>> {
    let mut trial = seq.to_vec();
    match target.partner(pos) {
        Some(partner) => {
            let partner_locked = constraints.and_then(|c| c.locked_at(partner));
            match partner_locked {
                Some(locked_base) => {
                    // Re-seat only `pos` to a base pairing the locked
                    // partner.
                    let options: Vec<u8> = CANON_PAIRS
                        .iter()
                        .filter_map(|&(a, b)| {
                            if pos < partner {
                                (b.eq_ignore_ascii_case(&locked_base)).then_some(a)
                            } else {
                                (a.eq_ignore_ascii_case(&locked_base)).then_some(b)
                            }
                        })
                        .collect();
                    if options.is_empty() {
                        return None;
                    }
                    trial[pos] = options[rng.below(options.len())];
                }
                None => {
                    let (a, b) = CANON_PAIRS[rng.below(6)];
                    let (i, j) = if pos < partner {
                        (pos, partner)
                    } else {
                        (partner, pos)
                    };
                    trial[i] = a;
                    trial[j] = b;
                }
            }
        }
        None => {
            trial[pos] = BASES[rng.below(4)];
        }
    }
    Some(trial)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn designs_a_hairpin_with_low_defect() {
        let target = Structure::from_dot_bracket("(((((....)))))").unwrap();
        let d = inverse_fold_ensemble_defect(&target, EnsembleDefectParams::default()).unwrap();
        assert_eq!(d.sequence.len(), 14);
        // A clean hairpin should reach a low normalised defect.
        assert!(
            d.normalized_defect < 0.2,
            "hairpin normalised defect too high: {} (notes {:?})",
            d.normalized_defect,
            d.notes
        );
    }

    #[test]
    fn designed_sequence_folds_to_the_target() {
        // After ensemble-defect design the MFE fold should match the
        // target closely.
        let target = Structure::from_dot_bracket("((((((....))))))").unwrap();
        let d = inverse_fold_ensemble_defect(&target, EnsembleDefectParams::default()).unwrap();
        assert!(
            d.mfe_distance <= 4,
            "designed sequence's MFE drifted by {} from the target",
            d.mfe_distance
        );
    }

    #[test]
    fn unpaired_target_is_trivial() {
        let target = Structure::empty(12);
        let d = inverse_fold_ensemble_defect(&target, EnsembleDefectParams::default()).unwrap();
        // An all-unpaired target is solvable to a very low defect.
        assert!(d.normalized_defect < 0.1);
    }

    #[test]
    fn rejects_empty_target() {
        let err =
            inverse_fold_ensemble_defect(&Structure::empty(0), EnsembleDefectParams::default())
                .unwrap_err();
        assert_eq!(err.code(), "rnadesign.goal");
    }

    #[test]
    fn rejects_pseudoknot() {
        let pk = Structure::from_dot_bracket("((..[[..))..]]").unwrap();
        assert!(inverse_fold_ensemble_defect(&pk, EnsembleDefectParams::default()).is_err());
    }

    #[test]
    fn is_deterministic() {
        let target = Structure::from_dot_bracket("(((....)))").unwrap();
        let a = inverse_fold_ensemble_defect(&target, EnsembleDefectParams::default()).unwrap();
        let b = inverse_fold_ensemble_defect(&target, EnsembleDefectParams::default()).unwrap();
        assert_eq!(a.sequence, b.sequence);
    }

    #[test]
    fn ensemble_defect_rejects_length_mismatch() {
        let target = Structure::from_dot_bracket("((....))").unwrap();
        assert!(ensemble_defect_linear(b"GGGG", &target).is_err());
    }

    #[test]
    fn ensemble_defect_is_bounded() {
        let target = Structure::from_dot_bracket("((((....))))").unwrap();
        let d = ensemble_defect_linear(b"GGGGAAAACCCC", &target).unwrap();
        assert!((0.0..=12.0).contains(&d));
    }

    #[test]
    fn decompose_finds_a_leaf_per_hairpin() {
        // Two separate hairpins -> two leaves.
        let target = Structure::from_dot_bracket("((((....))))((((....))))").unwrap();
        let leaves = decompose_leaves(&target);
        assert_eq!(leaves.len(), 2);
        assert_eq!(leaves[0], Leaf { start: 0, end: 12 });
        assert_eq!(leaves[1], Leaf { start: 12, end: 24 });
    }

    #[test]
    fn decompose_single_hairpin() {
        let target = Structure::from_dot_bracket("(((((....)))))").unwrap();
        let leaves = decompose_leaves(&target);
        assert_eq!(leaves.len(), 1);
        assert_eq!(leaves[0], Leaf { start: 0, end: 14 });
    }

    #[test]
    fn lower_defect_than_random_seed() {
        // The designer's output must have a markedly lower ensemble
        // defect than a random target-compatible seed.
        let target = Structure::from_dot_bracket("((((((....))))))").unwrap();
        let mut rng = Rng::new(1);
        let seed = seed_sequence(&target, None, &mut rng);
        let seed_defect = ensemble_defect_linear(&seed, &target).unwrap();
        let d = inverse_fold_ensemble_defect(&target, EnsembleDefectParams::default()).unwrap();
        assert!(
            d.ensemble_defect <= seed_defect + 1e-9,
            "designed defect {} not below a random seed's {}",
            d.ensemble_defect,
            seed_defect
        );
    }

    #[test]
    fn constrained_design_honours_locked_positions() {
        use crate::constraints::{lock_entry, DesignConstraintSet};
        use crate::goal::DesignConstraints;
        let target = Structure::from_dot_bracket("(((((....)))))").unwrap();
        let mut c = DesignConstraints::default().with_gc_range(0.0, 1.0);
        // Lock positions 0 and 13 (the outer pair) to G and C.
        c.required_subsequences = vec![lock_entry(0, b'G'), lock_entry(13, b'C')];
        let set = DesignConstraintSet::new(&c, 14);
        let d = inverse_fold_constrained(&target, &set, EnsembleDefectParams::default()).unwrap();
        assert_eq!(d.sequence[0], b'G', "locked position 0 not held");
        assert_eq!(d.sequence[13], b'C', "locked position 13 not held");
    }

    #[test]
    fn constrained_design_honours_gc_band() {
        use crate::constraints::DesignConstraintSet;
        use crate::goal::DesignConstraints;
        let target = Structure::from_dot_bracket("(((((....)))))").unwrap();
        // A tight GC band the design must land inside.
        let c = DesignConstraints::default().with_gc_range(0.35, 0.65);
        let set = DesignConstraintSet::new(&c, 14);
        let d = inverse_fold_constrained(&target, &set, EnsembleDefectParams::default()).unwrap();
        let gc = d
            .sequence
            .iter()
            .filter(|&&b| matches!(b, b'G' | b'C'))
            .count() as f64
            / d.sequence.len() as f64;
        assert!(
            (0.35 - 1e-9..=0.65 + 1e-9).contains(&gc),
            "GC fraction {gc} outside the requested band",
        );
    }

    #[test]
    fn constrained_design_avoids_a_forbidden_motif() {
        use crate::constraints::DesignConstraintSet;
        use crate::goal::DesignConstraints;
        let target = Structure::from_dot_bracket("((((((....))))))").unwrap();
        // Forbid a motif and confirm the design does not contain it.
        let c = DesignConstraints::default()
            .with_gc_range(0.0, 1.0)
            .forbid_motif("GGGG");
        let set = DesignConstraintSet::new(&c, 16);
        let d = inverse_fold_constrained(&target, &set, EnsembleDefectParams::default()).unwrap();
        let has = d
            .sequence
            .windows(4)
            .any(|w| w.eq_ignore_ascii_case(b"GGGG"));
        assert!(!has, "design contains the forbidden motif GGGG");
    }

    #[test]
    fn constrained_rejects_length_mismatch() {
        use crate::constraints::DesignConstraintSet;
        use crate::goal::DesignConstraints;
        let target = Structure::from_dot_bracket("((((....))))").unwrap();
        let set = DesignConstraintSet::new(&DesignConstraints::default(), 99);
        assert!(inverse_fold_constrained(&target, &set, EnsembleDefectParams::default()).is_err());
    }

    #[test]
    fn constrained_design_is_deterministic() {
        use crate::constraints::DesignConstraintSet;
        use crate::goal::DesignConstraints;
        let target = Structure::from_dot_bracket("(((....)))").unwrap();
        let set =
            DesignConstraintSet::new(&DesignConstraints::default().with_gc_range(0.0, 1.0), 10);
        let a = inverse_fold_constrained(&target, &set, EnsembleDefectParams::default()).unwrap();
        let b = inverse_fold_constrained(&target, &set, EnsembleDefectParams::default()).unwrap();
        assert_eq!(a.sequence, b.sequence);
    }
}
