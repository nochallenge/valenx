//! Feature 18 — multi-state (multi-target) sequence design (v1).
//!
//! Single-target inverse folding ([`crate::inverse`]) finds a sequence
//! that adopts *one* structure. Multi-state design — NUPACK's
//! *multi-tube* / multi-complex problem in spirit — finds *one*
//! sequence that adopts **two or more** target structures, each under
//! its own condition. The classic application is a conformational
//! switch: a sequence that folds to structure A in one context and to
//! structure B in another.
//!
//! ## The objective
//!
//! Each target *state* contributes its own **ensemble defect** — the
//! expected number of incorrectly-(un)paired nucleotides relative to
//! that state's structure (see [`crate::inverse`]). Multi-state design
//! minimises the **combined defect**
//!
//! ```text
//!     Σ_state  w_state · ensemble_defect(sequence, state)
//! ```
//!
//! a weighted sum over the states. A sequence with a low combined
//! defect adopts *every* target well; a sequence that folds perfectly
//! to one target but not the others scores poorly. The minimiser is a
//! mutation walk that re-seats positions and keeps a change if it
//! lowers the combined defect.
//!
//! ## Per-state constraints
//!
//! Each [`StateSpec`] carries its own target structure *and* its own
//! [`DesignConstraints`] — the "different constraint set" each state is
//! designed under. A position locked in one state's constraints (a
//! fixed nucleotide) is honoured globally; the GC / motif constraints
//! of every state are folded into the objective as soft penalties so
//! the search is steered toward a sequence legal under all of them.
//!
//! ## v1 scope — honest framing
//!
//! - This is an honest v1: a combined-ensemble-defect mutation walk. It
//!   is **not** a full NUPACK multi-tube designer — it does not model
//!   complex concentrations, strand stoichiometry, or per-tube
//!   partition functions across multiple strands. Every "state" here is
//!   the *same single strand* folding to a different target structure.
//! - The two (or more) targets must be the same length — they are
//!   alternative conformations of one molecule.
//! - A low combined ensemble defect means each target is, by the
//!   energy model, a well-populated conformation; it is a strong
//!   in-silico prediction of multi-stability, not a guarantee the
//!   physical RNA switches. Validate experimentally.
//! - The states' constraints are honoured: locked positions are never
//!   mutated; GC-range / forbidden-motif violations are penalised.

use crate::error::{Result, RnaDesignError};
use crate::goal::DesignConstraints;
use crate::inverse::ensemble_defect_linear;
use serde::{Deserialize, Serialize};
use valenx_rnastruct::{base_pair_distance, mfe, RnaSeq, Structure};

/// The four bases, ASCII.
const BASES: [u8; 4] = [b'A', b'C', b'G', b'U'];

/// Canonical pairing partners, ASCII.
const CANON_PAIRS: [(u8, u8); 6] = [
    (b'A', b'U'),
    (b'U', b'A'),
    (b'G', b'C'),
    (b'C', b'G'),
    (b'G', b'U'),
    (b'U', b'G'),
];

// ---------------------------------------------------------------------
// A small deterministic RNG.
// ---------------------------------------------------------------------

/// A deterministic xorshift RNG, local to the designer.
struct Rng {
    state: u64,
}

impl Rng {
    fn new(seed: u64) -> Self {
        Rng {
            state: seed ^ 0x9E37_79B9_7F4A_7C15,
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

    fn below(&mut self, n: usize) -> usize {
        (self.next_u64() % n as u64) as usize
    }

    fn unit(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }
}

// ---------------------------------------------------------------------
// State specification
// ---------------------------------------------------------------------

/// One target *state* of a multi-state design: a structure the single
/// designed sequence must adopt, with the constraints it must respect
/// in that state and a weight on its contribution to the combined
/// objective.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct StateSpec {
    /// A short human-readable label (`"ligand-free"`, `"bound"`, …).
    pub label: String,
    /// The target secondary structure for this state.
    pub target: Structure,
    /// The design constraints this state imposes — GC range, forbidden
    /// motifs, and locked positions (encoded in `required_subsequences`
    /// and parsed by [`locked_positions`]).
    pub constraints: DesignConstraints,
    /// The weight on this state's ensemble-defect contribution
    /// (`> 0`). Equal weights treat the states equally.
    pub weight: f64,
}

impl StateSpec {
    /// A state with a label, target and unit weight, default
    /// constraints.
    pub fn new(label: impl Into<String>, target: Structure) -> Self {
        StateSpec {
            label: label.into(),
            target,
            constraints: DesignConstraints::default(),
            weight: 1.0,
        }
    }

    /// Sets this state's constraints.
    pub fn with_constraints(mut self, constraints: DesignConstraints) -> Self {
        self.constraints = constraints;
        self
    }

    /// Sets this state's objective weight.
    pub fn with_weight(mut self, weight: f64) -> Self {
        self.weight = weight;
        self
    }
}

// ---------------------------------------------------------------------
// Parameters and result
// ---------------------------------------------------------------------

/// Parameters for [`design_multistate`].
#[derive(Copy, Clone, Debug)]
pub struct MultiStateParams {
    /// The mutation budget for the combined-defect search.
    pub iterations: usize,
    /// The combined-normalised-defect target — the search stops once the
    /// combined defect, divided by `(n · n_states)`, drops to or below
    /// this `[0, 1]` value.
    pub defect_target: f64,
    /// The random seed — the search is deterministic for a fixed seed.
    pub seed: u64,
}

impl Default for MultiStateParams {
    /// 1200 mutation steps, stop at a 5 %-combined-normalised defect.
    fn default() -> Self {
        MultiStateParams {
            iterations: 1200,
            defect_target: 0.05,
            seed: 0x5247,
        }
    }
}

/// The ensemble defect of one state in a finished multi-state design.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct StateDefect {
    /// The state's label.
    pub label: String,
    /// The ensemble defect of the design relative to this state.
    pub ensemble_defect: f64,
    /// The normalised ensemble defect (`/ n`) in `[0, 1]`.
    pub normalized_defect: f64,
    /// The base-pair distance of the design's MFE fold to this state's
    /// target.
    pub mfe_distance: usize,
}

/// The result of a multi-state design (feature 18).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct MultiStateDesign {
    /// The designed RNA sequence (`A C G U`).
    pub sequence: Vec<u8>,
    /// The combined weighted ensemble defect across every state.
    pub combined_defect: f64,
    /// Per-state ensemble-defect breakdown.
    pub states: Vec<StateDefect>,
    /// `true` if the combined normalised defect reached the target.
    pub solved: bool,
    /// Mutation steps accepted.
    pub accepted_steps: usize,
    /// Mutation steps attempted.
    pub total_steps: usize,
    /// Human-readable notes.
    pub notes: Vec<String>,
}

impl MultiStateDesign {
    /// The designed sequence as a `&str`.
    pub fn sequence_str(&self) -> &str {
        std::str::from_utf8(&self.sequence).unwrap_or("")
    }

    /// `true` when *every* state is adopted well — each state's
    /// normalised defect is below `threshold`.
    pub fn all_states_good(&self, threshold: f64) -> bool {
        self.states
            .iter()
            .all(|s| s.normalized_defect <= threshold + 1e-9)
    }
}

// ---------------------------------------------------------------------
// The designer
// ---------------------------------------------------------------------

/// Designs one sequence that adopts every target state (feature 18).
///
/// Minimises the weighted combined ensemble defect across the states by
/// a mutation walk. Each state may carry its own constraints; locked
/// positions are never mutated and GC / motif violations are penalised.
///
/// # Errors
/// - [`RnaDesignError::Goal`] if fewer than two states are given, any
///   target is empty or pseudoknotted, or the targets differ in length.
/// - [`RnaDesignError::Invalid`] if `params.iterations == 0` or a weight
///   is non-positive.
/// - [`RnaDesignError::Upstream`] if a folding call fails.
pub fn design_multistate(
    states: &[StateSpec],
    params: MultiStateParams,
) -> Result<MultiStateDesign> {
    if states.len() < 2 {
        return Err(RnaDesignError::goal(
            "states",
            "multi-state design needs at least two target states",
        ));
    }
    if params.iterations == 0 {
        return Err(RnaDesignError::invalid(
            "iterations",
            "need at least one mutation step",
        ));
    }
    let n = states[0].target.len();
    for s in states {
        if s.target.is_empty() {
            return Err(RnaDesignError::goal("states", "a state's target is empty"));
        }
        if s.target.has_pseudoknot() {
            return Err(RnaDesignError::goal(
                "states",
                format!(
                    "state `{}` is pseudoknotted — the folder is pseudoknot-free",
                    s.label
                ),
            ));
        }
        if s.target.len() != n {
            return Err(RnaDesignError::goal(
                "states",
                "every state's target must be the same length (one molecule, many folds)",
            ));
        }
        if !(s.weight.is_finite() && s.weight > 0.0) {
            return Err(RnaDesignError::invalid(
                "weight",
                "every state weight must be finite and positive",
            ));
        }
    }

    let mut rng = Rng::new(params.seed);
    let mut notes: Vec<String> = Vec::new();

    // Locked positions: a base fixed across the whole design. Collected
    // from every state's constraints (a position locked in any state is
    // locked globally) — see `locked_positions`.
    let locked = collect_locked(states, n);

    // Seed: a sequence compatible with the *first* state's pairs and
    // honouring the locked positions.
    let mut seq = seed_sequence(&states[0].target, &locked, &mut rng);

    let total_weight: f64 = states.iter().map(|s| s.weight).sum();
    let mut best = combined_defect(&seq, states)?;
    let mut accepted = 0usize;
    let mut total = 0usize;

    for _ in 0..params.iterations {
        let combined_norm = best / (n as f64 * total_weight);
        if combined_norm <= params.defect_target + 1e-9 {
            break;
        }
        total += 1;
        let trial = propose_mutation(&seq, states, &locked, &mut rng);
        let trial = match trial {
            Some(t) => t,
            None => continue,
        };
        let trial_defect = combined_defect(&trial, states)?;
        if trial_defect <= best {
            seq = trial;
            best = trial_defect;
            accepted += 1;
        }
    }

    // Per-state breakdown.
    let mut state_defects = Vec::with_capacity(states.len());
    let rna = RnaSeq::parse(&seq)?;
    let mfe_struct = mfe(&rna)?.structure;
    for s in states {
        let d = ensemble_defect_linear(&seq, &s.target)?;
        let mfe_dist = base_pair_distance(&mfe_struct, &s.target)?;
        state_defects.push(StateDefect {
            label: s.label.clone(),
            ensemble_defect: d,
            normalized_defect: d / n as f64,
            mfe_distance: mfe_dist,
        });
    }

    let combined_norm = best / (n as f64 * total_weight);
    let solved = combined_norm <= params.defect_target + 1e-9;

    notes.push(format!(
        "Multi-state design over {} target state(s) of {n} nt each.",
        states.len(),
    ));
    for sd in &state_defects {
        notes.push(format!(
            "State `{}`: ensemble defect {:.3} (normalised {:.4}), MFE distance {}.",
            sd.label, sd.ensemble_defect, sd.normalized_defect, sd.mfe_distance,
        ));
    }
    notes.push(format!(
        "Combined weighted defect {best:.3} (normalised {combined_norm:.4}); \
         {accepted}/{total} mutations accepted.",
    ));
    notes.push(
        "Multi-state design is an honest v1: a combined-ensemble-defect mutation walk over a \
         single strand. It does not model strand concentrations or multi-complex tubes. A low \
         combined defect is a strong in-silico prediction of multi-stability — validate \
         experimentally."
            .to_string(),
    );

    Ok(MultiStateDesign {
        sequence: seq,
        combined_defect: best,
        states: state_defects,
        solved,
        accepted_steps: accepted,
        total_steps: total,
        notes,
    })
}

/// Collects the globally-locked positions from every state's
/// constraints. A position is locked when a state pins an explicit base
/// at it via [`locked_positions`].
fn collect_locked(states: &[StateSpec], n: usize) -> Vec<Option<u8>> {
    let mut locked: Vec<Option<u8>> = vec![None; n];
    for s in states {
        for (pos, base) in locked_positions(&s.constraints, n) {
            if pos < n {
                locked[pos] = Some(base);
            }
        }
    }
    locked
}

/// Extracts locked (position, base) pairs from a constraint set.
///
/// A locked position is expressed as a one-character entry in
/// `required_subsequences` of the form `"@<index>=<base>"` — e.g.
/// `"@0=G"` locks position 0 to `G`. This keeps locked positions inside
/// the existing [`DesignConstraints`] type without a schema change.
/// Entries that do not parse are ignored.
pub fn locked_positions(constraints: &DesignConstraints, n: usize) -> Vec<(usize, u8)> {
    let mut out = Vec::new();
    for entry in &constraints.required_subsequences {
        let bytes = entry.as_bytes();
        if bytes.first() != Some(&b'@') {
            continue;
        }
        // @<index>=<base>
        let rest = &entry[1..];
        let Some(eq) = rest.find('=') else { continue };
        let (idx_str, base_str) = rest.split_at(eq);
        let base_str = &base_str[1..];
        let Ok(idx) = idx_str.parse::<usize>() else {
            continue;
        };
        if idx >= n {
            continue;
        }
        let base = base_str.trim().to_ascii_uppercase();
        let b = match base.as_bytes().first() {
            Some(&b @ (b'A' | b'C' | b'G' | b'U')) => b,
            Some(&b'T') => b'U',
            _ => continue,
        };
        out.push((idx, b));
    }
    out
}

/// Builds a seed sequence compatible with `first_target`'s pairs and the
/// locked positions.
fn seed_sequence(
    first_target: &Structure,
    locked: &[Option<u8>],
    rng: &mut Rng,
) -> Vec<u8> {
    let n = first_target.len();
    let mut seq = vec![b'A'; n];
    let mut assigned = vec![false; n];
    for bp in first_target.pairs() {
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
    // Locked positions override everything.
    for (i, lk) in locked.iter().enumerate() {
        if let Some(b) = lk {
            seq[i] = *b;
        }
    }
    seq
}

/// The combined weighted ensemble defect of `seq` over every state,
/// plus a soft penalty for each state's GC / forbidden-motif violations.
fn combined_defect(seq: &[u8], states: &[StateSpec]) -> Result<f64> {
    let mut total = 0.0;
    for s in states {
        let d = ensemble_defect_linear(seq, &s.target)?;
        total += s.weight * d;
        // Soft constraint penalty: GC out of range and forbidden motifs.
        total += s.weight * constraint_penalty(seq, &s.constraints);
    }
    Ok(total)
}

/// A soft penalty (in "defect units") for a sequence's violations of a
/// constraint set: GC excess and each forbidden-motif occurrence.
fn constraint_penalty(seq: &[u8], constraints: &DesignConstraints) -> f64 {
    let mut penalty = 0.0;
    // GC range.
    let gc = gc_fraction(seq);
    if gc < constraints.gc_min {
        penalty += (constraints.gc_min - gc) * seq.len() as f64;
    } else if gc > constraints.gc_max {
        penalty += (gc - constraints.gc_max) * seq.len() as f64;
    }
    // Forbidden motifs / subsequences.
    for motif in constraints
        .forbidden_motifs
        .iter()
        .chain(constraints.forbidden_subsequences.iter())
    {
        let needle: Vec<u8> = motif
            .bytes()
            .map(|b| match b.to_ascii_uppercase() {
                b'T' => b'U',
                other => other,
            })
            .collect();
        if needle.is_empty() {
            continue;
        }
        penalty += count_occurrences(seq, &needle) as f64;
    }
    penalty
}

/// GC fraction of an RNA sequence.
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

/// Counts non-overlapping... actually overlapping occurrences of
/// `needle` in `haystack` (case-insensitive).
fn count_occurrences(haystack: &[u8], needle: &[u8]) -> usize {
    if needle.is_empty() || needle.len() > haystack.len() {
        return 0;
    }
    haystack
        .windows(needle.len())
        .filter(|w| w.eq_ignore_ascii_case(needle))
        .count()
}

/// Proposes a mutation that respects *every* state: a position paired in
/// any state is re-seated as a canonical pair consistent with the state
/// it is paired in; locked positions are never touched.
fn propose_mutation(
    seq: &[u8],
    states: &[StateSpec],
    locked: &[Option<u8>],
    rng: &mut Rng,
) -> Option<Vec<u8>> {
    let n = seq.len();
    if n == 0 {
        return None;
    }
    // Pick a position that is not locked.
    let free: Vec<usize> = (0..n).filter(|&i| locked[i].is_none()).collect();
    if free.is_empty() {
        return None;
    }
    let pos = free[rng.below(free.len())];

    // Find a state in which `pos` is paired — re-seat that pair so the
    // chosen state stays satisfiable. If `pos` is paired in several
    // states to *different* partners, pick one state at random; the
    // search will balance the others.
    let paired_states: Vec<(usize, usize)> = states
        .iter()
        .enumerate()
        .filter_map(|(si, s)| s.target.partner(pos).map(|p| (si, p)))
        .collect();

    let mut trial = seq.to_vec();
    if paired_states.is_empty() {
        // `pos` is unpaired in every state — mutate it freely.
        trial[pos] = BASES[rng.below(4)];
    } else {
        let (_si, partner) = paired_states[rng.below(paired_states.len())];
        // Re-seat the pair (pos, partner) with a canonical pair, unless
        // the partner is locked — then only mutate `pos` to a base that
        // pairs the locked partner.
        if let Some(locked_partner) = locked[partner] {
            // Choose a base for `pos` that pairs the locked partner.
            let options: Vec<u8> = CANON_PAIRS
                .iter()
                .filter_map(|&(a, b)| {
                    if pos < partner {
                        (b == locked_partner).then_some(a)
                    } else {
                        (a == locked_partner).then_some(b)
                    }
                })
                .collect();
            if options.is_empty() {
                return None;
            }
            trial[pos] = options[rng.below(options.len())];
        } else {
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
    // 20 % of the time also nudge a second random free position — helps
    // escape local optima where one swap alone cannot improve.
    if rng.unit() < 0.2 && free.len() > 1 {
        let pos2 = free[rng.below(free.len())];
        if states.iter().all(|s| s.target.partner(pos2).is_none()) {
            trial[pos2] = BASES[rng.below(4)];
        }
    }
    Some(trial)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn two_states() -> Vec<StateSpec> {
        let a = Structure::from_dot_bracket("((((....))))....").unwrap();
        let b = Structure::from_dot_bracket("....((((....))))").unwrap();
        vec![
            StateSpec::new("state-A", a),
            StateSpec::new("state-B", b),
        ]
    }

    #[test]
    fn designs_a_two_state_sequence() {
        let d = design_multistate(&two_states(), MultiStateParams::default()).unwrap();
        assert_eq!(d.sequence.len(), 16);
        assert_eq!(d.states.len(), 2);
        assert!(!d.notes.is_empty());
    }

    #[test]
    fn combined_defect_is_low_after_design() {
        let d = design_multistate(&two_states(), MultiStateParams::default()).unwrap();
        // Each state should be adopted reasonably — the combined defect
        // should be far below a random seed's.
        let mut rng = Rng::new(7);
        let locked = vec![None; 16];
        let seed = seed_sequence(&two_states()[0].target, &locked, &mut rng);
        let states = two_states();
        let seed_defect = combined_defect(&seed, &states).unwrap();
        assert!(
            d.combined_defect <= seed_defect + 1e-9,
            "designed combined defect {} not below seed {}",
            d.combined_defect,
            seed_defect
        );
    }

    #[test]
    fn rejects_single_state() {
        let one = vec![StateSpec::new(
            "only",
            Structure::from_dot_bracket("((((....))))").unwrap(),
        )];
        assert!(design_multistate(&one, MultiStateParams::default()).is_err());
    }

    #[test]
    fn rejects_length_mismatch() {
        let states = vec![
            StateSpec::new("a", Structure::from_dot_bracket("((((....))))").unwrap()),
            StateSpec::new(
                "b",
                Structure::from_dot_bracket("((((....))))....").unwrap(),
            ),
        ];
        assert!(design_multistate(&states, MultiStateParams::default()).is_err());
    }

    #[test]
    fn rejects_pseudoknot() {
        let states = vec![
            StateSpec::new(
                "a",
                Structure::from_dot_bracket("((..[[..))..]]").unwrap(),
            ),
            StateSpec::new("b", Structure::empty(14)),
        ];
        assert!(design_multistate(&states, MultiStateParams::default()).is_err());
    }

    #[test]
    fn rejects_zero_iterations() {
        let params = MultiStateParams {
            iterations: 0,
            ..MultiStateParams::default()
        };
        assert!(design_multistate(&two_states(), params).is_err());
    }

    #[test]
    fn is_deterministic() {
        let a = design_multistate(&two_states(), MultiStateParams::default()).unwrap();
        let b = design_multistate(&two_states(), MultiStateParams::default()).unwrap();
        assert_eq!(a.sequence, b.sequence);
    }

    #[test]
    fn locked_positions_parse() {
        let c = DesignConstraints {
            required_subsequences: vec![
                "@0=G".to_string(),
                "@3=C".to_string(),
                "garbage".to_string(),
                "@99=A".to_string(), // out of range
            ],
            ..DesignConstraints::default()
        };
        let locked = locked_positions(&c, 16);
        assert_eq!(locked, vec![(0, b'G'), (3, b'C')]);
    }

    #[test]
    fn locked_positions_are_honoured() {
        // Lock position 0 to G in state A's constraints.
        let ca = DesignConstraints {
            required_subsequences: vec!["@0=G".to_string()],
            ..DesignConstraints::default()
        };
        let states = vec![
            StateSpec::new("a", Structure::from_dot_bracket("....((((....))))").unwrap())
                .with_constraints(ca),
            StateSpec::new("b", Structure::from_dot_bracket("((((....))))....").unwrap()),
        ];
        let d = design_multistate(&states, MultiStateParams::default()).unwrap();
        assert_eq!(d.sequence[0], b'G', "locked position 0 was not held at G");
    }

    #[test]
    fn weights_are_validated() {
        let mut states = two_states();
        states[0].weight = 0.0;
        assert!(design_multistate(&states, MultiStateParams::default()).is_err());
    }

    #[test]
    fn all_states_good_predicate() {
        let d = design_multistate(&two_states(), MultiStateParams::default()).unwrap();
        // The predicate is consistent with the per-state defects.
        let thr = 0.5;
        let expected = d.states.iter().all(|s| s.normalized_defect <= thr + 1e-9);
        assert_eq!(d.all_states_good(thr), expected);
    }
}
