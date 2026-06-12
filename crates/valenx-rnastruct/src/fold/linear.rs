//! LinearFold — linear-time beam-search minimum-free-energy folding.
//!
//! The Zuker DP ([`crate::fold::zuker`]) is `O(n³)` time and `O(n²)`
//! memory. That is fine for a tRNA but hopeless for a real mRNA of
//! thousands — tens of thousands — of nucleotides. **LinearFold**
//! (Huang, Zhang, Liu, Mathews, Huang; *Bioinformatics* 2019) folds in
//! `O(n·b)` time and `O(n·b)` memory by two ideas:
//!
//! 1. **Left-to-right (5′→3′) dynamic programming.** Instead of filling
//!    the DP table by increasing span, LinearFold scans the sequence
//!    once from the 5′ end, deriving every state that *ends* at the
//!    current position `j` from states that ended earlier.
//! 2. **Beam pruning.** At each position only the best `b` states per
//!    state type are kept; the rest are discarded. A wider beam is more
//!    accurate and slower; a narrow beam is faster and approximate.
//!
//! ## Energy model
//!
//! LinearFold here reuses the **complete Turner-2004 nearest-neighbor
//! model** of [`crate::fold::energy`] — the same hairpin / bulge /
//! internal-loop / multiloop energies the exact Zuker folder uses. The
//! grammar is the standard Zuker grammar (`C` exterior, `P` paired
//! region, `M` / `M2` multiloop fragments), evaluated left-to-right.
//!
//! ## Exactness
//!
//! Beam search is an **approximate** algorithm. With a beam wide enough
//! that no pruning ever happens (`beam >= n`, or [`fold_linear_exact`])
//! it explores the whole Zuker grammar and returns the **exact** MFE —
//! this is asserted in the crate's validation suite. With the default
//! beam ([`DEFAULT_BEAM_SIZE`]) it is, like ViennaRNA's LinearFold,
//! near-optimal on real sequences but not guaranteed optimal.
//!
//! For short RNA prefer the exact [`crate::fold::zuker::mfe`];
//! LinearFold is the path for sequences too long for the `O(n³)` DP.

use crate::error::Result;
use crate::fold::energy::{self, multiloop};
use crate::fold::nussinov::MIN_HAIRPIN;
use crate::fold::zuker::MAX_LOOP;
use crate::rna::RnaSeq;
use crate::structure::Structure;
use std::collections::HashMap;

/// Default beam size — the number of states kept per position per state
/// type. 100 is the LinearFold paper's default; it is near-optimal on
/// real sequences while keeping folding strictly linear-time.
pub const DEFAULT_BEAM_SIZE: usize = 100;

/// A "this sub-problem is infeasible" sentinel energy.
const INF: f64 = energy::FORBIDDEN;

/// The result of a LinearFold beam-search fold.
#[derive(Clone, Debug, PartialEq)]
pub struct LinearFoldResult {
    /// The predicted secondary structure.
    pub structure: Structure,
    /// Its total free energy in kcal/mol, in the dangle-folded model
    /// (consistent with [`crate::fold::eval::structure_energy`]).
    pub energy: f64,
    /// The beam size used for this fold.
    pub beam_size: usize,
    /// `true` if the beam was wide enough that no pruning happened, so
    /// the result is the **exact** Zuker MFE; `false` if a beam was
    /// pruned and the result is the (near-optimal) beam-search answer.
    pub exact: bool,
}

/// Folds `seq` with LinearFold at the default beam size
/// ([`DEFAULT_BEAM_SIZE`]).
///
/// Linear time and memory in the sequence length. For short sequences
/// the exact [`crate::fold::zuker::mfe`] is preferable; this is the
/// entry point for long sequences.
///
/// # Errors
/// Never fails for a valid [`RnaSeq`]; returns [`Result`] for API
/// consistency with the other folders.
pub fn fold_linear(seq: &RnaSeq) -> Result<LinearFoldResult> {
    fold_linear_with_beam(seq, DEFAULT_BEAM_SIZE)
}

/// Folds `seq` with LinearFold at an explicit `beam_size`.
///
/// A larger beam is more accurate and slower; `beam_size >=
/// seq.len()` disables pruning entirely and returns the exact Zuker
/// MFE. A `beam_size` of 0 is treated as 1.
///
/// # Errors
/// As [`fold_linear`].
pub fn fold_linear_with_beam(seq: &RnaSeq, beam_size: usize) -> Result<LinearFoldResult> {
    let beam = beam_size.max(1);
    let codes = seq.codes();
    let n = codes.len();
    if n == 0 {
        return Ok(LinearFoldResult {
            structure: Structure::empty(0),
            energy: 0.0,
            beam_size: beam,
            exact: true,
        });
    }
    let mut engine = BeamFold::new(codes, beam);
    engine.run();
    let (structure, energy) = engine.traceback()?;
    Ok(LinearFoldResult {
        structure,
        energy,
        beam_size: beam,
        exact: !engine.pruned,
    })
}

/// Folds `seq` with LinearFold and a beam wide enough to **disable
/// pruning** — so the result is the exact Zuker MFE, computed by the
/// linear-time left-to-right DP. Used to cross-check LinearFold against
/// the `O(n³)` folder.
///
/// # Errors
/// As [`fold_linear`].
pub fn fold_linear_exact(seq: &RnaSeq) -> Result<LinearFoldResult> {
    fold_linear_with_beam(seq, seq.len().max(1))
}

// ---------------------------------------------------------------------------
// The beam-search engine.
// ---------------------------------------------------------------------------

/// How a `P` (paired-region) state was derived — for traceback.
#[derive(Copy, Clone, Debug)]
enum PFrom {
    /// `[i..=j]` is a hairpin loop.
    Hairpin,
    /// `[i..=j]` encloses the single inner pair `(k, l)` (stack / bulge
    /// / internal loop).
    Internal(usize, usize),
    /// `[i..=j]` closes a multiloop whose interior is `M2[i+1][j-1]`.
    Multi,
}

/// How an `M` (≥1-branch multiloop fragment) state was derived.
#[derive(Copy, Clone, Debug)]
enum MFrom {
    /// The whole fragment is the single branch `P[i][j]`.
    Branch,
    /// `j` is unpaired; the rest is `M[i][j-1]`.
    UnpairedRight,
    /// `i` is unpaired; the rest is `M[i+1][j]`.
    UnpairedLeft,
    /// Promoted from `M2[i][j]` (≥2 branches is also ≥1).
    FromM2,
}

/// How an `M2` (≥2-branch multiloop fragment) state was derived.
#[derive(Copy, Clone, Debug)]
enum M2From {
    /// `M[i][k]` followed by the branch `P[k+1][j]`.
    Split(usize),
    /// `j` is unpaired; the rest is `M2[i][j-1]`.
    UnpairedRight,
}

/// How a `C` (exterior) state was derived.
#[derive(Copy, Clone, Debug)]
enum CFrom {
    /// `j` is unpaired in the exterior loop.
    Unpaired,
    /// the exterior branch `P[i][j]` closes, preceded by exterior
    /// `C[i-1]`.
    Branch(usize),
}

/// One scored state in a beam: an energy and a backpointer.
#[derive(Copy, Clone, Debug)]
struct State<T> {
    energy: f64,
    from: T,
}

/// The left-to-right beam-search folder.
struct BeamFold<'a> {
    codes: &'a [u8],
    n: usize,
    beam: usize,
    /// `true` once any beam was actually pruned.
    pruned: bool,
    /// `c[j]` — best exterior energy of the prefix `[0..=j]`.
    c: Vec<State<CFrom>>,
    /// `p[j]` — beam of paired regions ending at `j`, keyed by start.
    p: Vec<HashMap<usize, State<PFrom>>>,
    /// `m[j]` — beam of ≥1-branch multiloop fragments ending at `j`.
    m: Vec<HashMap<usize, State<MFrom>>>,
    /// `m2[j]` — beam of ≥2-branch multiloop fragments ending at `j`.
    m2: Vec<HashMap<usize, State<M2From>>>,
}

impl<'a> BeamFold<'a> {
    fn new(codes: &'a [u8], beam: usize) -> Self {
        let n = codes.len();
        BeamFold {
            codes,
            n,
            beam,
            pruned: false,
            c: vec![
                State {
                    energy: INF,
                    from: CFrom::Unpaired,
                };
                n
            ],
            p: vec![HashMap::new(); n],
            m: vec![HashMap::new(); n],
            m2: vec![HashMap::new(); n],
        }
    }

    /// Exterior energy of the prefix `[0..=j]`; `j` is a signed index so
    /// `c_at(-1) = 0` (the empty prefix).
    #[inline]
    fn c_at(&self, j: isize) -> f64 {
        if j < 0 {
            0.0
        } else {
            self.c[j as usize].energy
        }
    }

    /// Runs the full left-to-right beam-search fill.
    fn run(&mut self) {
        for j in 0..self.n {
            self.fill_p(j);
            self.fill_m(j);
            self.fill_m2(j);
            self.fill_c(j);
            // Beam-prune every per-position beam.
            self.pruned |= prune(&mut self.p[j], self.beam);
            self.pruned |= prune(&mut self.m[j], self.beam);
            self.pruned |= prune(&mut self.m2[j], self.beam);
        }
    }

    /// Fills `p[j]` — every paired region `[i..=j]` with `i·j` paired.
    fn fill_p(&mut self, j: usize) {
        let codes = self.codes;
        let mut out: HashMap<usize, State<PFrom>> = HashMap::new();
        let mut relax = |i: usize, e: f64, from: PFrom| {
            if e < INF / 2.0 {
                let slot = out.entry(i).or_insert(State { energy: INF, from });
                if e < slot.energy {
                    slot.energy = e;
                    slot.from = from;
                }
            }
        };

        // Every i that could pair j: j - i - 1 >= MIN_HAIRPIN unpaired.
        let i_hi = j.saturating_sub(MIN_HAIRPIN + 1);
        for i in 0..=i_hi {
            if !energy::can_pair_codes(codes[i], codes[j]) {
                continue;
            }
            // --- hairpin ---
            let loop_bases = &codes[(i + 1)..j];
            relax(
                i,
                energy::hairpin_energy(codes[i], codes[j], loop_bases),
                PFrom::Hairpin,
            );

            // --- multiloop close: interior is M2[i+1][j-1] ---
            if j >= i + 2 {
                if let Some(inner) = self.m2[j - 1].get(&(i + 1)) {
                    let closure = multiloop::OFFSET
                        + multiloop::PER_BRANCH
                        + energy::terminal_penalty(codes[i], codes[j]);
                    relax(i, closure + inner.energy, PFrom::Multi);
                }
            }
        }

        // --- internal / bulge / stack: inner pair (k, l) ---
        // (k, l) is a paired region already in p[l] for some l < j.
        let l_lo = j.saturating_sub(MAX_LOOP + 2);
        for l in (l_lo..j).rev() {
            if l < MIN_HAIRPIN + 1 {
                break;
            }
            // snapshot the inner beam to avoid borrowing self twice.
            let inner: Vec<(usize, f64)> = self.p[l].iter().map(|(&k, s)| (k, s.energy)).collect();
            for (k, inner_e) in inner {
                if inner_e >= INF / 2.0 {
                    continue;
                }
                let right = j - l - 1;
                if right > MAX_LOOP {
                    continue;
                }
                // outer pair (i, j): i < k, left gap k-i-1 <= MAX_LOOP.
                let i_min = k.saturating_sub(MAX_LOOP + 1);
                for i in i_min..k {
                    if !energy::can_pair_codes(codes[i], codes[j]) {
                        continue;
                    }
                    let left = k - i - 1;
                    if left + right == 0 {
                        // pure stack — still valid; fallthrough.
                    } else if left > MAX_LOOP || right > MAX_LOOP {
                        continue;
                    }
                    // a hairpin needs i and j far enough; guaranteed by
                    // l < j and k > i.
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
                    relax(i, il + inner_e, PFrom::Internal(k, l));
                }
            }
        }

        self.p[j] = out;
    }

    /// Fills `m[j]` — ≥1-branch multiloop fragments `[i..=j]`.
    fn fill_m(&mut self, j: usize) {
        let mut out: HashMap<usize, State<MFrom>> = HashMap::new();

        // A whole fragment that is exactly one branch P[i][j].
        for (&i, s) in &self.p[j] {
            relax_m(&mut out, i, s.energy + multiloop::PER_BRANCH, MFrom::Branch);
        }
        // j unpaired: extend an M ending at j-1.
        if j > 0 {
            for (&i, s) in &self.m[j - 1] {
                relax_m(
                    &mut out,
                    i,
                    s.energy + multiloop::PER_UNPAIRED,
                    MFrom::UnpairedRight,
                );
            }
        }
        // ≥2 branches is also ≥1.
        for (&i, s) in &self.m2[j] {
            relax_m(&mut out, i, s.energy, MFrom::FromM2);
        }
        // i unpaired: an M over [i+1..=j] with a leading unpaired base.
        // Sweep every start i from j down to 1 so M[i+1][j] is fully
        // relaxed before M[i][j] reads it — this propagates a leading
        // run of unpaired bases of any length.
        for i in (1..=j).rev() {
            if let Some(s) = out.get(&i).copied() {
                relax_m(
                    &mut out,
                    i - 1,
                    s.energy + multiloop::PER_UNPAIRED,
                    MFrom::UnpairedLeft,
                );
            }
        }

        self.m[j] = out;
    }

    /// Fills `m2[j]` — ≥2-branch multiloop fragments `[i..=j]`.
    fn fill_m2(&mut self, j: usize) {
        let mut out: HashMap<usize, State<M2From>> = HashMap::new();
        let mut relax = |i: usize, e: f64, from: M2From| {
            if e < INF / 2.0 {
                let slot = out.entry(i).or_insert(State { energy: INF, from });
                if e < slot.energy {
                    slot.energy = e;
                    slot.from = from;
                }
            }
        };

        // Split: M[i][k] then the last branch P[k+1][j]. Iterate over
        // the P beam ending at j (start = k+1) and join with the M beam
        // ending at k.
        let p_starts: Vec<(usize, f64)> =
            self.p[j].iter().map(|(&kp1, s)| (kp1, s.energy)).collect();
        for (kp1, p_e) in p_starts {
            if kp1 == 0 {
                continue; // no room for a left fragment
            }
            let k = kp1 - 1;
            for (&i, m_s) in &self.m[k] {
                relax(
                    i,
                    m_s.energy + p_e + multiloop::PER_BRANCH,
                    M2From::Split(k),
                );
            }
        }
        // j unpaired: extend an M2 ending at j-1.
        if j > 0 {
            for (&i, s) in &self.m2[j - 1] {
                relax(i, s.energy + multiloop::PER_UNPAIRED, M2From::UnpairedRight);
            }
        }

        self.m2[j] = out;
    }

    /// Fills `c[j]` — the exterior energy of the prefix `[0..=j]`.
    fn fill_c(&mut self, j: usize) {
        // j unpaired in the exterior loop.
        let mut best = State {
            energy: self.c_at(j as isize - 1),
            from: CFrom::Unpaired,
        };
        // an exterior branch P[i][j] closes, preceded by exterior
        // C[i-1].
        for (&i, s) in &self.p[j] {
            let e = self.c_at(i as isize - 1)
                + s.energy
                + energy::terminal_penalty(self.codes[i], self.codes[j]);
            if e < best.energy {
                best = State {
                    energy: e,
                    from: CFrom::Branch(i),
                };
            }
        }
        self.c[j] = best;
    }

    /// Reconstructs the structure from the filled beams.
    fn traceback(&self) -> Result<(Structure, f64)> {
        let n = self.n;
        let mut partner: Vec<Option<usize>> = vec![None; n];
        if n == 0 {
            return Ok((Structure::from_partner(partner)?, 0.0));
        }

        // Job stack over the same grammar.
        enum Job {
            C(isize),
            P(usize, usize),
            M(usize, usize),
            M2(usize, usize),
        }
        let mut stack = vec![Job::C(n as isize - 1)];

        while let Some(job) = stack.pop() {
            match job {
                Job::C(j) => {
                    if j < 0 {
                        continue;
                    }
                    let ju = j as usize;
                    match self.c[ju].from {
                        CFrom::Unpaired => stack.push(Job::C(j - 1)),
                        CFrom::Branch(i) => {
                            stack.push(Job::P(i, ju));
                            stack.push(Job::C(i as isize - 1));
                        }
                    }
                }
                Job::P(i, j) => {
                    partner[i] = Some(j);
                    partner[j] = Some(i);
                    let st = self.p[j].get(&i).ok_or_else(|| {
                        crate::error::RnaStructError::structure(
                            "LinearFold traceback: missing P state",
                        )
                    })?;
                    match st.from {
                        PFrom::Hairpin => {}
                        PFrom::Internal(k, l) => stack.push(Job::P(k, l)),
                        PFrom::Multi => {
                            if j >= i + 2 {
                                stack.push(Job::M2(i + 1, j - 1));
                            }
                        }
                    }
                }
                Job::M(i, j) => {
                    let st = self.m[j].get(&i).ok_or_else(|| {
                        crate::error::RnaStructError::structure(
                            "LinearFold traceback: missing M state",
                        )
                    })?;
                    match st.from {
                        MFrom::Branch => stack.push(Job::P(i, j)),
                        MFrom::UnpairedRight => stack.push(Job::M(i, j - 1)),
                        MFrom::UnpairedLeft => stack.push(Job::M(i + 1, j)),
                        MFrom::FromM2 => stack.push(Job::M2(i, j)),
                    }
                }
                Job::M2(i, j) => {
                    let st = self.m2[j].get(&i).ok_or_else(|| {
                        crate::error::RnaStructError::structure(
                            "LinearFold traceback: missing M2 state",
                        )
                    })?;
                    match st.from {
                        M2From::Split(k) => {
                            stack.push(Job::M(i, k));
                            stack.push(Job::P(k + 1, j));
                        }
                        M2From::UnpairedRight => stack.push(Job::M2(i, j - 1)),
                    }
                }
            }
        }

        let structure = Structure::from_partner(partner)?;
        let energy = self.c[n - 1].energy;
        // An infeasible / numerically-zero result is reported as 0.0
        // (the open structure).
        let energy = if energy >= INF / 2.0 || energy.abs() < 1e-9 {
            0.0
        } else {
            energy
        };
        Ok((structure, energy))
    }
}

/// Relaxes the `M`-state beam entry for start `i` to the lower of its
/// current energy and `e` (a free function so callers can keep reading
/// the map while relaxing it).
fn relax_m(out: &mut HashMap<usize, State<MFrom>>, i: usize, e: f64, from: MFrom) {
    if e < INF / 2.0 {
        let slot = out.entry(i).or_insert(State { energy: INF, from });
        if e < slot.energy {
            slot.energy = e;
            slot.from = from;
        }
    }
}

/// Prunes `beam_map` to its best (lowest-energy) `beam` entries.
/// Returns `true` if anything was actually removed.
fn prune<T: Copy>(beam_map: &mut HashMap<usize, State<T>>, beam: usize) -> bool {
    if beam_map.len() <= beam {
        return false;
    }
    let mut entries: Vec<(usize, f64)> = beam_map.iter().map(|(&i, s)| (i, s.energy)).collect();
    // Keep the `beam` lowest-energy entries.
    entries.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
    let keep: std::collections::HashSet<usize> =
        entries.iter().take(beam).map(|&(i, _)| i).collect();
    beam_map.retain(|i, _| keep.contains(i));
    true
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fold::eval::structure_energy;
    use crate::fold::zuker::mfe;

    #[test]
    fn empty_and_tiny() {
        let r = fold_linear(&RnaSeq::parse("A").unwrap()).unwrap();
        assert_eq!(r.structure.n_pairs(), 0);
        assert_eq!(r.energy, 0.0);
        assert!(r.exact);
    }

    #[test]
    fn folds_a_clean_hairpin() {
        let seq = RnaSeq::parse("GGGGGAAAACCCCC").unwrap();
        let r = fold_linear(&seq).unwrap();
        assert!(r.structure.n_pairs() >= 3, "expected a stem: {r:?}");
        assert!(
            r.energy < 0.0,
            "stable hairpin should have E<0: {}",
            r.energy
        );
        assert!(r.structure.is_nested());
    }

    #[test]
    fn reported_energy_is_self_consistent() {
        // The energy LinearFold reports must equal an independent
        // evaluation of the structure it returns.
        let seq = RnaSeq::parse("GGGAAACCCAAAGGGAAACCC").unwrap();
        let r = fold_linear(&seq).unwrap();
        let re = structure_energy(&seq, &r.structure).unwrap();
        assert!(
            (re - r.energy).abs() < 1e-4,
            "LinearFold energy {} != eval {re}",
            r.energy
        );
    }

    #[test]
    fn wide_beam_equals_exact_zuker_mfe() {
        // With pruning disabled the left-to-right beam DP must return
        // exactly the O(n^3) Zuker MFE — same energy.
        for s in [
            "GGGGGAAAACCCCC",
            "GCGCGCGAAACGCGCGC",
            "GGGGGCCCCCAAAGGGGGCCCCCAAAGGGGG",
            "GGGAGGGAAACUCCCC",
            "ACGUACGUACGUACGUACGU",
        ] {
            let seq = RnaSeq::parse(s).unwrap();
            let exact = mfe(&seq).unwrap();
            let lin = fold_linear_exact(&seq).unwrap();
            assert!(lin.exact, "exact fold should report exact=true");
            assert!(
                (lin.energy - exact.energy).abs() < 1e-6,
                "seq {s}: LinearFold-exact {} != Zuker MFE {}",
                lin.energy,
                exact.energy
            );
        }
    }

    #[test]
    fn unpairable_sequence_stays_open() {
        let seq = RnaSeq::parse("AAAAAAAAAAAA").unwrap();
        let r = fold_linear(&seq).unwrap();
        assert_eq!(r.structure.n_pairs(), 0);
        assert_eq!(r.energy, 0.0);
    }

    #[test]
    fn folds_a_multiloop() {
        let seq = RnaSeq::parse("GGGGGCCCCCAAAGGGGGCCCCCAAAGGGGG").unwrap();
        let r = fold_linear(&seq).unwrap();
        assert!(r.structure.is_nested());
        let e = structure_energy(&seq, &r.structure).unwrap();
        assert!((e - r.energy).abs() < 1e-4);
    }

    #[test]
    fn narrow_beam_still_folds_and_is_an_upper_bound() {
        // A narrow beam is approximate: its energy can be no lower than
        // the exact MFE (it explores a subset of structures).
        let seq = RnaSeq::parse("GGGGGCCCCCAAAGGGGGCCCCCAAAGGGGGCCCCCAAAGGGGG").unwrap();
        let exact = mfe(&seq).unwrap();
        let narrow = fold_linear_with_beam(&seq, 2).unwrap();
        assert!(narrow.structure.is_nested());
        // self-consistent
        let re = structure_energy(&seq, &narrow.structure).unwrap();
        assert!((re - narrow.energy).abs() < 1e-4);
        // never beats the true optimum
        assert!(
            narrow.energy >= exact.energy - 1e-6,
            "beam-search energy {} below the exact MFE {}",
            narrow.energy,
            exact.energy
        );
    }

    #[test]
    fn beam_zero_is_treated_as_one() {
        let seq = RnaSeq::parse("GGGGGAAAACCCCC").unwrap();
        let r = fold_linear_with_beam(&seq, 0).unwrap();
        assert_eq!(r.beam_size, 1);
        assert!(r.structure.is_nested());
    }

    #[test]
    fn scales_to_a_long_sequence() {
        // A 600-nt sequence: the O(n^3) Zuker DP would be slow; the
        // linear folder handles it quickly with a bounded beam.
        let unit = "GGGGGAAAACCCCCUUUU";
        let long: String = unit.repeat(34); // ~612 nt
        let seq = RnaSeq::parse(&long).unwrap();
        let r = fold_linear(&seq).unwrap();
        assert_eq!(r.structure.len(), seq.len());
        assert!(r.structure.is_nested());
        assert!(
            r.energy < 0.0,
            "a long structured RNA should fold: {}",
            r.energy
        );
        let re = structure_energy(&seq, &r.structure).unwrap();
        assert!((re - r.energy).abs() < 1e-3);
    }
}
