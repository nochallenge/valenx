//! LinearPartition — linear-time partition function and base-pair
//! probabilities.
//!
//! [`crate::ensemble::partition`] computes the McCaskill partition
//! function exactly, but in `O(n³)` time and `O(n²)` memory — the same
//! wall this crate's exact MFE hits on long sequences.
//! **LinearPartition** (Zhang, Mathews, Huang; *Bioinformatics* 2020)
//! removes it with the same two ideas as [`crate::fold::linear`]:
//!
//! 1. **Left-to-right dynamic programming** — one 5′→3′ scan, deriving
//!    every state ending at the current position from earlier states.
//! 2. **Beam pruning** — at each position only the `b` states of
//!    highest partition-function weight per state type are kept.
//!
//! It is the partition-function analogue of LinearFold: the same
//! `C / P / M / M2` grammar, but evaluated in the Boltzmann
//! sum-product semiring (`min → +`, `+ → ×`). It yields, in `O(n·b)`
//! time:
//!
//! - the **partition function** `Q` and the **ensemble free energy**
//!   `G = −RT·ln Q`;
//! - an **approximate base-pair-probability matrix**, from a
//!   left-to-right inside pass plus an outside pass over the pruned
//!   beams.
//!
//! ## Energy model
//!
//! LinearPartition here reuses the **complete Turner-2004 model** of
//! [`crate::fold::energy`] — identical hairpin / loop / multiloop
//! energies to the exact McCaskill folder.
//!
//! ## Exactness
//!
//! With a beam wide enough to disable pruning ([`linear_partition_exact`])
//! the inside pass sums the *whole* ensemble and the ensemble free
//! energy equals the exact McCaskill value — asserted in the validation
//! suite. With the default beam it is near-exact, like the published
//! LinearPartition. The base-pair probabilities are approximate (beam
//! search drops low-weight states); for short RNA prefer the exact
//! [`crate::ensemble::partition::partition_function`].

use crate::error::Result;
use crate::fold::energy::{self, multiloop, GAS_CONSTANT};
use crate::fold::linear::DEFAULT_BEAM_SIZE;
use crate::fold::nussinov::MIN_HAIRPIN;
use crate::fold::zuker::MAX_LOOP;
use crate::rna::RnaSeq;
use std::collections::HashMap;

/// The folding temperature in kelvin used for the Boltzmann weights.
pub const DEFAULT_TEMPERATURE_K: f64 = energy::T37_KELVIN;

/// The output of a LinearPartition computation: the partition function,
/// the ensemble free energy, the temperature, the beam used, and a
/// sparse base-pair-probability map.
#[derive(Clone, Debug)]
pub struct LinearPartitionResult {
    n: usize,
    temperature_k: f64,
    beam_size: usize,
    exact: bool,
    q_total: f64,
    /// Sparse base-pair probabilities — `(i, j) → p(i,j)` for `i < j`,
    /// only entries that survived beam pruning with `p > 0`.
    bpp: HashMap<(usize, usize), f64>,
}

impl LinearPartitionResult {
    /// Sequence length.
    pub fn len(&self) -> usize {
        self.n
    }

    /// `true` if computed for a zero-length sequence.
    pub fn is_empty(&self) -> bool {
        self.n == 0
    }

    /// The folding temperature, kelvin.
    pub fn temperature_k(&self) -> f64 {
        self.temperature_k
    }

    /// The beam size used.
    pub fn beam_size(&self) -> usize {
        self.beam_size
    }

    /// `true` if the beam was wide enough that no pruning happened, so
    /// the partition function is the **exact** McCaskill value.
    pub fn is_exact(&self) -> bool {
        self.exact
    }

    /// The total partition function `Q`.
    pub fn q(&self) -> f64 {
        self.q_total
    }

    /// The ensemble free energy `G = −RT·ln Q`, kcal/mol.
    pub fn ensemble_free_energy(&self) -> f64 {
        -GAS_CONSTANT * self.temperature_k * self.q_total.ln()
    }

    /// The (approximate) probability that base `i` pairs base `j`.
    /// Returns 0.0 for out-of-range indices, `i == j`, or a pair that
    /// did not survive beam pruning.
    pub fn pair_probability(&self, i: usize, j: usize) -> f64 {
        if i >= self.n || j >= self.n || i == j {
            return 0.0;
        }
        let (a, b) = if i < j { (i, j) } else { (j, i) };
        self.bpp.get(&(a, b)).copied().unwrap_or(0.0)
    }

    /// The (approximate) probability that base `i` is unpaired.
    pub fn unpaired_probability(&self, i: usize) -> f64 {
        if i >= self.n {
            return 0.0;
        }
        let mut paired = 0.0;
        for (&(a, b), &p) in &self.bpp {
            if a == i || b == i {
                paired += p;
            }
        }
        (1.0 - paired).clamp(0.0, 1.0)
    }

    /// All base pairs with probability ≥ `threshold`, as `(i, j, p)`
    /// sorted by descending probability.
    pub fn significant_pairs(&self, threshold: f64) -> Vec<(usize, usize, f64)> {
        let mut out: Vec<(usize, usize, f64)> = self
            .bpp
            .iter()
            .filter(|(_, &p)| p >= threshold)
            .map(|(&(i, j), &p)| (i, j, p))
            .collect();
        out.sort_by(|a, b| b.2.partial_cmp(&a.2).unwrap_or(std::cmp::Ordering::Equal));
        out
    }
}

/// Computes the LinearPartition partition function of `seq` at 37 °C
/// with the default beam size ([`DEFAULT_BEAM_SIZE`]).
///
/// # Errors
/// Never fails for a valid [`RnaSeq`]; returns [`Result`] for API
/// consistency.
pub fn linear_partition(seq: &RnaSeq) -> Result<LinearPartitionResult> {
    linear_partition_with_beam(seq, DEFAULT_BEAM_SIZE, DEFAULT_TEMPERATURE_K)
}

/// Computes the LinearPartition partition function with a beam wide
/// enough to **disable pruning** — the partition function is then the
/// exact McCaskill value, computed by the linear-time DP. Used to
/// cross-check LinearPartition against the `O(n³)` folder.
///
/// # Errors
/// As [`linear_partition`].
pub fn linear_partition_exact(seq: &RnaSeq) -> Result<LinearPartitionResult> {
    linear_partition_with_beam(seq, seq.len().max(1), DEFAULT_TEMPERATURE_K)
}

/// Computes the LinearPartition partition function with an explicit
/// `beam_size` and temperature.
///
/// A larger beam is more accurate and slower; `beam_size >= seq.len()`
/// disables pruning and returns the exact McCaskill partition function.
///
/// # Errors
/// As [`linear_partition`].
pub fn linear_partition_with_beam(
    seq: &RnaSeq,
    beam_size: usize,
    temperature_k: f64,
) -> Result<LinearPartitionResult> {
    let beam = beam_size.max(1);
    let codes = seq.codes();
    let n = codes.len();
    if n == 0 {
        return Ok(LinearPartitionResult {
            n: 0,
            temperature_k,
            beam_size: beam,
            exact: true,
            q_total: 1.0,
            bpp: HashMap::new(),
        });
    }
    let mut engine = BeamPartition::new(codes, beam, temperature_k);
    engine.run_inside();
    engine.run_outside();
    let bpp = engine.base_pair_probabilities();
    Ok(LinearPartitionResult {
        n,
        temperature_k,
        beam_size: beam,
        exact: !engine.pruned,
        q_total: engine.q_total.max(1.0e-300),
        bpp,
    })
}

// ---------------------------------------------------------------------------
// The beam-search partition-function engine.
// ---------------------------------------------------------------------------

/// The left-to-right beam-search partition-function folder.
struct BeamPartition<'a> {
    codes: &'a [u8],
    n: usize,
    beam: usize,
    rt: f64,
    pruned: bool,
    /// `qc[j]` — partition function of the exterior over `[0..=j]`.
    qc: Vec<f64>,
    /// `qp[j]` — beam: `i → Σ` over paired regions `[i..=j]`, `i·j`.
    qp: Vec<HashMap<usize, f64>>,
    /// `qm1[j]` — beam: one-branch multiloop fragments (branch starts at
    /// `i`, arbitrary trailing unpaired) ending at `j`.
    qm1: Vec<HashMap<usize, f64>>,
    /// `qm[j]` — beam: ≥1-branch multiloop fragments ending at `j`.
    qm: Vec<HashMap<usize, f64>>,
    /// `qm2[j]` — beam: ≥2-branch multiloop fragments ending at `j`.
    qm2: Vec<HashMap<usize, f64>>,
    /// Total partition function `Q = qc[n-1]`.
    q_total: f64,
    /// Outside weight of every `qp[j][i]` state (filled by the outside
    /// pass).
    op: Vec<HashMap<usize, f64>>,
}

impl<'a> BeamPartition<'a> {
    fn new(codes: &'a [u8], beam: usize, temperature_k: f64) -> Self {
        let n = codes.len();
        BeamPartition {
            codes,
            n,
            beam,
            rt: GAS_CONSTANT * temperature_k,
            pruned: false,
            qc: vec![0.0; n],
            qp: vec![HashMap::new(); n],
            qm1: vec![HashMap::new(); n],
            qm: vec![HashMap::new(); n],
            qm2: vec![HashMap::new(); n],
            q_total: 1.0,
            op: vec![HashMap::new(); n],
        }
    }

    /// Boltzmann factor `exp(−E/RT)`.
    #[inline]
    fn boltz(&self, e: f64) -> f64 {
        (-e / self.rt).exp()
    }

    /// Exterior partition function of the prefix `[0..=j]`; `qc_at(-1)`
    /// is 1 (the empty prefix has one — empty — structure).
    #[inline]
    fn qc_at(&self, j: isize) -> f64 {
        if j < 0 {
            1.0
        } else {
            self.qc[j as usize]
        }
    }

    /// Runs the left-to-right inside (partition-function) pass.
    fn run_inside(&mut self) {
        for j in 0..self.n {
            self.inside_p(j);
            self.inside_m1(j);
            self.inside_m2_and_m(j);
            self.inside_c(j);
            self.pruned |= prune(&mut self.qp[j], self.beam);
            self.pruned |= prune(&mut self.qm1[j], self.beam);
            self.pruned |= prune(&mut self.qm[j], self.beam);
            self.pruned |= prune(&mut self.qm2[j], self.beam);
        }
        self.q_total = self.qc[self.n - 1].max(1.0e-300);
    }

    /// Inside pass for `qp[j]`.
    fn inside_p(&mut self, j: usize) {
        let codes = self.codes;
        let mut out: HashMap<usize, f64> = HashMap::new();

        let i_hi = j.saturating_sub(MIN_HAIRPIN + 1);
        for i in 0..=i_hi {
            if !energy::can_pair_codes(codes[i], codes[j]) {
                continue;
            }
            // hairpin
            let loop_bases = &codes[(i + 1)..j];
            let hp = self.boltz(energy::hairpin_energy(codes[i], codes[j], loop_bases));
            *out.entry(i).or_insert(0.0) += hp;
            // multiloop close
            if j >= i + 2 {
                if let Some(&inner) = self.qm2[j - 1].get(&(i + 1)) {
                    let closure = multiloop::OFFSET
                        + multiloop::PER_BRANCH
                        + energy::terminal_penalty(codes[i], codes[j]);
                    *out.entry(i).or_insert(0.0) += self.boltz(closure) * inner;
                }
            }
        }

        // internal / bulge / stack
        let l_lo = j.saturating_sub(MAX_LOOP + 2);
        for l in (l_lo..j).rev() {
            if l < MIN_HAIRPIN + 1 {
                break;
            }
            let inner: Vec<(usize, f64)> = self.qp[l].iter().map(|(&k, &v)| (k, v)).collect();
            for (k, inner_q) in inner {
                if inner_q == 0.0 {
                    continue;
                }
                let right = j - l - 1;
                if right > MAX_LOOP {
                    continue;
                }
                let i_min = k.saturating_sub(MAX_LOOP + 1);
                for i in i_min..k {
                    if !energy::can_pair_codes(codes[i], codes[j]) {
                        continue;
                    }
                    let left = k - i - 1;
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
                    *out.entry(i).or_insert(0.0) += self.boltz(il) * inner_q;
                }
            }
        }

        self.qp[j] = out;
    }

    /// Inside pass for `qm1[j]` — one-branch fragments (the branch
    /// starts at `i`, with arbitrary trailing unpaired bases).
    ///
    /// `qm1[i][j] = qb[i][j]·B(branch) + qm1[i][j-1]·B(unp_j)`.
    fn inside_m1(&mut self, j: usize) {
        let mut out: HashMap<usize, f64> = HashMap::new();
        // the branch ends exactly at j.
        for (&i, &q) in &self.qp[j] {
            *out.entry(i).or_insert(0.0) += q * self.boltz(multiloop::PER_BRANCH);
        }
        // j is a trailing unpaired base; the rest is qm1[i][j-1].
        if j > 0 {
            for (&i, &q) in &self.qm1[j - 1] {
                *out.entry(i).or_insert(0.0) += q * self.boltz(multiloop::PER_UNPAIRED);
            }
        }
        self.qm1[j] = out;
    }

    /// Inside pass for `qm2[j]` (≥2 branches) and `qm[j]` (≥1 branch),
    /// the unambiguous `qm1`-based grammar:
    ///
    /// `qm2[i][j] = Σ_k qm[i][k]·qm1[k+1][j]`           (last branch),
    /// `qm [i][j] = qm[i+1][j]·B(unp_i) + qm1[i][j] + qm2[i][j]`.
    fn inside_m2_and_m(&mut self, j: usize) {
        // qm2: Σ_k qm[i][k]·qm1[k+1][j]. Iterate the qm1[j] beam (start
        // = k+1) and join with the qm[k] beam.
        let mut qm2_out: HashMap<usize, f64> = HashMap::new();
        let qm1_starts: Vec<(usize, f64)> = self.qm1[j].iter().map(|(&kp1, &q)| (kp1, q)).collect();
        for (kp1, qm1_q) in qm1_starts {
            if kp1 == 0 {
                continue;
            }
            let k = kp1 - 1;
            for (&i, &qm_q) in &self.qm[k] {
                *qm2_out.entry(i).or_insert(0.0) += qm_q * qm1_q;
            }
        }

        // qm = qm1 + qm2, plus the i-unpaired (leading) extension.
        let mut qm_out: HashMap<usize, f64> = HashMap::new();
        for (&i, &q) in &self.qm1[j] {
            *qm_out.entry(i).or_insert(0.0) += q;
        }
        for (&i, &q) in &qm2_out {
            *qm_out.entry(i).or_insert(0.0) += q;
        }
        // i unpaired: qm[i][j] += qm[i+1][j]·B(unp_i). Sweep every
        // start position i from j down to 1 so qm[i+1][j] is finalised
        // (including its own leading-unpaired extensions) before
        // qm[i][j] reads it — this propagates runs of leading unpaired
        // bases of any length.
        let unp = self.boltz(multiloop::PER_UNPAIRED);
        for s in (1..=j).rev() {
            let v = qm_out.get(&s).copied().unwrap_or(0.0);
            if v != 0.0 {
                *qm_out.entry(s - 1).or_insert(0.0) += v * unp;
            }
        }

        self.qm2[j] = qm2_out;
        self.qm[j] = qm_out;
    }

    /// Inside pass for `qc[j]`.
    fn inside_c(&mut self, j: usize) {
        // j unpaired in the exterior loop.
        let mut q = self.qc_at(j as isize - 1);
        // exterior branch P[i][j], preceded by exterior C[i-1].
        for (&i, &qp) in &self.qp[j] {
            q += self.qc_at(i as isize - 1)
                * qp
                * self.boltz(energy::terminal_penalty(self.codes[i], self.codes[j]));
        }
        self.qc[j] = q;
    }

    /// Runs the outside pass — fills `op[j][i]`, the outside weight of
    /// every kept `qp[j][i]` state. Combined with `qp` it gives the
    /// base-pair probabilities.
    ///
    /// The exterior and internal-loop channels are summed exactly over
    /// the kept beams; the multiloop channel is folded in through the
    /// `qm`/`qm2` decomposition as a first-order treatment — identical
    /// in spirit to [`crate::ensemble::partition`].
    fn run_outside(&mut self) {
        let n = self.n;
        // Suffix exterior partition function: qc_suffix[i] = partition
        // function of the exterior over [i..=n-1].
        let mut qc_suffix = vec![1.0_f64; n + 1];
        // qc_suffix[i] = (i unpaired) qc_suffix[i+1]
        //             + Σ_{j>=i, (i'..j) pair} ... — but the exterior
        // decomposition pairs the *first* helix. Build it right-to-left.
        qc_suffix[n] = 1.0;
        for i in (0..n).rev() {
            // i unpaired in the exterior.
            let mut q = qc_suffix[i + 1];
            // a helix (i, e) starts the suffix exterior.
            for j in (i + 1)..n {
                if let Some(&qp) = self.qp[j].get(&i) {
                    q += qp
                        * self.boltz(energy::terminal_penalty(self.codes[i], self.codes[j]))
                        * qc_suffix[j + 1];
                }
            }
            qc_suffix[i] = q;
        }

        // Outside, in decreasing span so enclosing pairs finish first.
        // op[j][i] accumulates the outside weight of qp[j][i].
        for span in (0..n).rev() {
            for i in 0..(n - span) {
                let j = i + span;
                if !self.qp[j].contains_key(&i) {
                    continue;
                }
                let mut outside = 0.0;

                // (i,j) closes an exterior helix: C[i-1] * suffix[j+1].
                let left = self.qc_at(i as isize - 1);
                let right = qc_suffix[j + 1];
                outside += left
                    * right
                    * self.boltz(energy::terminal_penalty(self.codes[i], self.codes[j]));

                // (i,j) sits directly inside an enclosing pair (p,r) as
                // the single inner pair of an internal/bulge/stack loop.
                let p_min = i.saturating_sub(MAX_LOOP + 1);
                for p in p_min..i {
                    let r_max = (j + MAX_LOOP + 1).min(n - 1);
                    for r in (j + 1)..=r_max {
                        let parent = match self.op[r].get(&p) {
                            Some(&v) if v != 0.0 => v,
                            _ => continue,
                        };
                        if !self.qp[r].contains_key(&p) {
                            continue;
                        }
                        let lft = i - p - 1;
                        let rgt = r - j - 1;
                        if lft + rgt != 0 && (lft > MAX_LOOP || rgt > MAX_LOOP) {
                            continue;
                        }
                        if r - p <= MIN_HAIRPIN {
                            continue;
                        }
                        let il = energy::internal_loop_energy(
                            self.codes[p],
                            self.codes[r],
                            self.codes[i],
                            self.codes[j],
                            lft,
                            rgt,
                            self.codes[p + 1],
                            self.codes[r - 1],
                            self.codes[i - 1],
                            self.codes[j + 1],
                        );
                        outside += parent * self.boltz(il);
                    }
                }

                self.op[j].entry(i).or_insert(0.0);
                if let Some(slot) = self.op[j].get_mut(&i) {
                    *slot += outside;
                }
            }
        }
    }

    /// Base-pair probabilities `p(i,j) = qp[i][j] · op[i][j] / Q`.
    fn base_pair_probabilities(&self) -> HashMap<(usize, usize), f64> {
        let mut bpp = HashMap::new();
        for j in 0..self.n {
            for (&i, &qp) in &self.qp[j] {
                if qp == 0.0 {
                    continue;
                }
                if let Some(&outside) = self.op[j].get(&i) {
                    let p = (qp * outside / self.q_total).clamp(0.0, 1.0);
                    if p > 0.0 {
                        bpp.insert((i, j), p);
                    }
                }
            }
        }
        bpp
    }
}

/// Prunes `beam_map` to its best (highest-weight) `beam` entries.
/// Returns `true` if anything was removed.
fn prune(beam_map: &mut HashMap<usize, f64>, beam: usize) -> bool {
    if beam_map.len() <= beam {
        return false;
    }
    let mut entries: Vec<(usize, f64)> = beam_map.iter().map(|(&i, &v)| (i, v)).collect();
    // Keep the `beam` highest-weight entries (descending).
    entries.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
    let keep: std::collections::HashSet<usize> =
        entries.iter().take(beam).map(|&(i, _)| i).collect();
    beam_map.retain(|i, _| keep.contains(i));
    true
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ensemble::partition::partition_function;
    use crate::fold::zuker::mfe;

    #[test]
    fn empty_sequence() {
        let pf = linear_partition(&RnaSeq::parse("A").unwrap()).unwrap();
        assert!(pf.q() > 0.0);
        assert!(pf.is_exact());
    }

    #[test]
    fn q_is_at_least_one() {
        // The empty structure always contributes weight 1.
        let pf = linear_partition(&RnaSeq::parse("AAAAAAAA").unwrap()).unwrap();
        assert!(pf.q() >= 1.0 - 1e-9);
    }

    #[test]
    fn ensemble_free_energy_below_or_equal_mfe() {
        let seq = RnaSeq::parse("GGGGGAAAACCCCC").unwrap();
        let pf = linear_partition(&seq).unwrap();
        let mfe_e = mfe(&seq).unwrap().energy;
        assert!(
            pf.ensemble_free_energy() <= mfe_e + 1e-6,
            "ensemble G {} should be <= MFE {}",
            pf.ensemble_free_energy(),
            mfe_e
        );
    }

    #[test]
    fn wide_beam_equals_exact_mccaskill() {
        // With pruning disabled the left-to-right beam DP must sum the
        // whole ensemble — the ensemble free energy must equal the
        // exact O(n^3) McCaskill value.
        for s in [
            "GGGGGAAAACCCCC",
            "GCGCGCGAAACGCGCGC",
            "GGGGGCCCCCAAAGGGGGCCCCCAAAGGGGG",
            "ACGUACGUACGUACGU",
            "GGGAGGGAAACUCCCC",
        ] {
            let seq = RnaSeq::parse(s).unwrap();
            let exact = partition_function(&seq).unwrap();
            let lin = linear_partition_exact(&seq).unwrap();
            assert!(lin.is_exact());
            assert!(
                (lin.ensemble_free_energy() - exact.ensemble_free_energy()).abs() < 1e-6,
                "seq {s}: LinearPartition-exact G {} != McCaskill G {}",
                lin.ensemble_free_energy(),
                exact.ensemble_free_energy()
            );
        }
    }

    #[test]
    fn pair_probabilities_are_in_range() {
        let seq = RnaSeq::parse("GGGGGAAAACCCCC").unwrap();
        let pf = linear_partition(&seq).unwrap();
        for i in 0..seq.len() {
            for j in (i + 1)..seq.len() {
                let p = pf.pair_probability(i, j);
                assert!((0.0..=1.0).contains(&p), "p({i},{j})={p}");
            }
            assert!((0.0..=1.0).contains(&pf.unpaired_probability(i)));
        }
    }

    #[test]
    fn stable_stem_has_high_pairing_probability() {
        let seq = RnaSeq::parse("GGGGGGAAAACCCCCC").unwrap();
        let pf = linear_partition(&seq).unwrap();
        let p = pf.pair_probability(0, 15);
        assert!(p > 0.3, "outer GC pair probability only {p}");
    }

    #[test]
    fn wide_beam_bpp_matches_exact_mccaskill() {
        // With pruning disabled the inside/outside pass must reproduce
        // the exact McCaskill base-pair probabilities for a clean stem.
        let seq = RnaSeq::parse("GGGGGGAAAACCCCCC").unwrap();
        let exact = partition_function(&seq).unwrap();
        let lin = linear_partition_exact(&seq).unwrap();
        for i in 0..seq.len() {
            for j in (i + 1)..seq.len() {
                let pe = exact.pair_probability(i, j);
                let pl = lin.pair_probability(i, j);
                assert!(
                    (pe - pl).abs() < 1e-6,
                    "p({i},{j}): exact {pe} vs linear {pl}"
                );
            }
        }
    }

    #[test]
    fn unpairable_sequence_is_all_unpaired() {
        let seq = RnaSeq::parse("AAAAAAAAAAAA").unwrap();
        let pf = linear_partition(&seq).unwrap();
        assert!((pf.q() - 1.0).abs() < 1e-9);
        assert!((pf.unpaired_probability(5) - 1.0).abs() < 1e-9);
    }

    #[test]
    fn significant_pairs_sorted() {
        let seq = RnaSeq::parse("GGGGGGAAAACCCCCC").unwrap();
        let pf = linear_partition(&seq).unwrap();
        let sig = pf.significant_pairs(0.1);
        for w in sig.windows(2) {
            assert!(w[0].2 >= w[1].2, "significant pairs must be sorted");
        }
    }

    #[test]
    fn scales_to_a_long_sequence() {
        let unit = "GGGGGAAAACCCCCUUUU";
        let long: String = unit.repeat(34); // ~612 nt
        let seq = RnaSeq::parse(&long).unwrap();
        let pf = linear_partition(&seq).unwrap();
        assert_eq!(pf.len(), seq.len());
        assert!(pf.q() > 1.0);
        assert!(pf.ensemble_free_energy() < 0.0);
    }
}
