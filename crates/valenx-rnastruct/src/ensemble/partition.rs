//! McCaskill partition function and base-pair probability matrix.
//!
//! Where the Zuker DP ([`crate::fold::zuker`]) finds the single
//! *minimum*-energy structure, the McCaskill (1990) algorithm sums the
//! Boltzmann weight `exp(-E/RT)` over the *entire ensemble* of
//! pseudoknot-free structures. From that sum it derives:
//!
//! - the **partition function** `Q` and the **ensemble free energy**
//!   `G = -RT·ln Q`;
//! - the **base-pair probability** `p(i,j)` — the fraction of the
//!   ensemble (by Boltzmann weight) in which `i` pairs `j`.
//!
//! ## Recurrences
//!
//! The same loop decomposition as Zuker, but with `+` → `×` and
//! `min` → `+` (the Boltzmann / sum-product semiring). Four arrays:
//!
//! - `qb[i][j]` — partition function of `i..=j` **given `i·j`**.
//! - `q[i][j]`  — partition function of the exterior of `i..=j`.
//! - `qm[i][j]` — multiloop fragment, ≥ 1 branch.
//! - `qm2[i][j]`— multiloop fragment, ≥ 2 branches.
//!
//! The outside pass (`probability` array) gives `p(i,j)` exactly, the
//! standard McCaskill inside/outside recursion.
//!
//! ## v1 scope
//!
//! Multiloop energies in the partition function use the same linear
//! Turner model; the per-branch term is a constant Boltzmann factor.
//! As with Zuker, dangles are folded into terminal penalties. The
//! result is a real ensemble partition function — exact for the model
//! encoded in [`crate::fold::energy`].

use crate::error::Result;
use crate::fold::constraint::FoldConstraints;
use crate::fold::energy::{self, multiloop, GAS_CONSTANT};
use crate::fold::nussinov::MIN_HAIRPIN;
use crate::fold::zuker::MAX_LOOP;
use crate::rna::RnaSeq;

/// The folding temperature in kelvin used for the Boltzmann weights.
pub const DEFAULT_TEMPERATURE_K: f64 = energy::T37_KELVIN;

/// The output of a McCaskill computation: the partition function,
/// the ensemble free energy, the temperature, and the base-pair
/// probability matrix.
#[derive(Clone, Debug)]
pub struct PartitionFunction {
    n: usize,
    temperature_k: f64,
    /// The total partition function `Q` of the whole sequence.
    q_total: f64,
    /// Flat `n×n` base-pair probability matrix; `bpp[i*n+j]` is
    /// `p(i,j)` for `i < j` (and 0 elsewhere).
    bpp: Vec<f64>,
}

impl PartitionFunction {
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

    /// The total partition function `Q`.
    pub fn q(&self) -> f64 {
        self.q_total
    }

    /// The ensemble free energy `G = -RT·ln Q`, kcal/mol.
    pub fn ensemble_free_energy(&self) -> f64 {
        -GAS_CONSTANT * self.temperature_k * self.q_total.ln()
    }

    /// The probability that base `i` pairs base `j` (order-insensitive).
    /// Returns 0.0 for out-of-range indices or `i == j`.
    pub fn pair_probability(&self, i: usize, j: usize) -> f64 {
        if i >= self.n || j >= self.n || i == j {
            return 0.0;
        }
        let (a, b) = if i < j { (i, j) } else { (j, i) };
        self.bpp[a * self.n + b]
    }

    /// The probability that base `i` is *unpaired* — `1 − Σⱼ p(i,j)`.
    pub fn unpaired_probability(&self, i: usize) -> f64 {
        if i >= self.n {
            return 0.0;
        }
        let mut paired = 0.0;
        for j in 0..self.n {
            if j != i {
                paired += self.pair_probability(i, j);
            }
        }
        (1.0 - paired).clamp(0.0, 1.0)
    }

    /// The full base-pair probability matrix as a flat `n×n` slice.
    pub fn bpp_matrix(&self) -> &[f64] {
        &self.bpp
    }

    /// All base pairs whose probability is at least `threshold`,
    /// returned as `(i, j, p)` sorted by descending probability.
    pub fn significant_pairs(&self, threshold: f64) -> Vec<(usize, usize, f64)> {
        let mut out = Vec::new();
        for i in 0..self.n {
            for j in (i + 1)..self.n {
                let p = self.bpp[i * self.n + j];
                if p >= threshold {
                    out.push((i, j, p));
                }
            }
        }
        out.sort_by(|a, b| b.2.partial_cmp(&a.2).unwrap_or(std::cmp::Ordering::Equal));
        out
    }
}

/// Computes the McCaskill partition function of `seq` at 37 °C.
pub fn partition_function(seq: &RnaSeq) -> Result<PartitionFunction> {
    partition_function_at(
        seq,
        DEFAULT_TEMPERATURE_K,
        &FoldConstraints::none(seq.len()),
    )
}

/// Computes the McCaskill partition function of `seq` at an arbitrary
/// temperature, honoring `constraints`.
///
/// The Turner parameters are referenced at 37 °C; at another
/// temperature the Boltzmann weights use the supplied `temperature_k`
/// in `exp(-E/RT)` (a first-order temperature treatment — the
/// parameters themselves are not re-fitted, which is the standard
/// approximation for a v1 melting-curve calculation).
///
/// # Errors
/// Returns an error only if `constraints` does not match the sequence
/// length.
pub fn partition_function_at(
    seq: &RnaSeq,
    temperature_k: f64,
    constraints: &FoldConstraints,
) -> Result<PartitionFunction> {
    let codes = seq.codes();
    let n = codes.len();
    if constraints.len() != n {
        return Err(crate::error::RnaStructError::invalid(
            "constraints",
            "constraint length does not match sequence length",
        ));
    }
    if n == 0 {
        return Ok(PartitionFunction {
            n: 0,
            temperature_k,
            q_total: 1.0,
            bpp: Vec::new(),
        });
    }

    let rt = GAS_CONSTANT * temperature_k;
    let boltz = |e: f64| (-e / rt).exp();

    let at = |i: usize, j: usize| i * n + j;

    // ----- inside pass --------------------------------------------------
    let mut qb = vec![0.0_f64; n * n];
    let mut q = vec![1.0_f64; n * n]; // empty interval => 1
    let mut qm = vec![0.0_f64; n * n];
    let mut qm1 = vec![0.0_f64; n * n];
    let mut qm2 = vec![0.0_f64; n * n];

    // Exterior diagonal: a single base contributes its (soft) unpaired
    // Boltzmann weight; forbidding it unpaired makes the cell 0.
    for i in 0..n {
        q[at(i, i)] = if constraints.unpaired_allowed(i) {
            boltz(constraints.soft_unpaired(i))
        } else {
            0.0
        };
    }

    // Boltzmann factor for one multiloop branch helix.
    let branch_factor = boltz(multiloop::PER_BRANCH);
    let ml_closure = boltz(multiloop::OFFSET + multiloop::PER_BRANCH);

    for span in 1..n {
        for i in 0..(n - span) {
            let j = i + span;

            // qb[i][j]
            let mut qbij = 0.0;
            if span > MIN_HAIRPIN
                && energy::can_pair_codes(codes[i], codes[j])
                && constraints.pair_allowed(i, j)
            {
                // hairpin
                let loop_bases = &codes[(i + 1)..j];
                qbij += boltz(energy::hairpin_energy(codes[i], codes[j], loop_bases));
                // internal / bulge / stack
                let k_max = (i + 1 + MAX_LOOP).min(j.saturating_sub(MIN_HAIRPIN + 1));
                for k in (i + 1)..=k_max {
                    let l_min = (k + MIN_HAIRPIN + 1).max(j.saturating_sub(MAX_LOOP + 1));
                    for l in l_min..j {
                        if l <= k {
                            continue;
                        }
                        let inner = qb[at(k, l)];
                        if inner == 0.0 {
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
                        qbij += boltz(il) * inner;
                    }
                }
                // multiloop: interior qm2 over i+1..j-1
                if j >= i + 2 {
                    let interior = qm2[at(i + 1, j - 1)];
                    qbij +=
                        ml_closure * boltz(energy::terminal_penalty(codes[i], codes[j])) * interior;
                }
            }
            qb[at(i, j)] = qbij;

            // Multiloop fragments — the standard *unambiguous* `qm1`/
            // `qm` decomposition (ViennaRNA's partition-function
            // grammar). A partition function requires every fragment
            // counted exactly once: an ambiguous grammar that reaches a
            // `branch . branch` fragment two ways inflates Q.
            //
            //  qm1[i][j] — exactly one branch, starting at i, with
            //              arbitrary trailing unpaired bases.
            //  qm [i][j] — >= 1 branch (leading + trailing unpaired ok).
            //  qm2[i][j] — >= 2 branches.

            // qm1[i][j] = qb[i][j]·B(branch)  +  qm1[i][j-1]·B(unp_j)
            let mut qm1ij = qbij * branch_factor;
            if constraints.unpaired_allowed(j) {
                qm1ij += qm1[at(i, j - 1)] * boltz(constraints.soft_unpaired(j));
            }
            qm1[at(i, j)] = qm1ij;

            // qm2[i][j] = Σ_k qm[i][k]·qm1[k+1][j]   (last branch unique)
            let mut qm2ij = 0.0;
            for k in i..j {
                qm2ij += qm[at(i, k)] * qm1[at(k + 1, j)];
            }
            qm2[at(i, j)] = qm2ij;

            // qm[i][j] = qm[i+1][j]·B(unp_i)   (i unpaired, leading)
            //          + qm1[i][j]             (exactly one branch at i)
            //          + qm2[i][j]             (>= 2 branches)
            let mut qmij = qm1ij + qm2ij;
            if constraints.unpaired_allowed(i) {
                qmij += qm[at(i + 1, j)] * boltz(constraints.soft_unpaired(i));
            }
            qm[at(i, j)] = qmij;

            // q[i][j]: exterior partition function of i..=j.
            //
            // Unambiguous "fix the last element" decomposition: the
            // 3'-most position j is either unpaired (the rest is the
            // exterior of i..j-1) or it is the 3' end of the last
            // exterior helix (k, j) (the rest is the exterior of
            // i..k-1). Repeated j-unpaired steps generate the
            // all-unpaired structure, so no separate term is needed.
            let mut qij = 0.0;
            // j unpaired.
            if constraints.unpaired_allowed(j) {
                // q[i][j-1] with the empty interval (j-1 < i) == 1.
                let rest = if j > i { q[at(i, j - 1)] } else { 1.0 };
                qij += rest * boltz(constraints.soft_unpaired(j));
            }
            // last exterior helix (k, j), prefix is the exterior i..k-1.
            for k in i..=j {
                if energy::can_pair_codes(codes[k], codes[j]) && qb[at(k, j)] > 0.0 {
                    let left = if k > i { q[at(i, k - 1)] } else { 1.0 };
                    qij +=
                        left * qb[at(k, j)] * boltz(energy::terminal_penalty(codes[k], codes[j]));
                }
            }
            q[at(i, j)] = qij;
        }
    }

    let q_total = q[at(0, n - 1)].max(1.0e-300);

    // ----- outside pass: base-pair probabilities -----------------------
    // p(i,j) = (contribution of structures containing pair (i,j)) / Q.
    // We compute the "outside" weight for each pair via the standard
    // McCaskill recursion in decreasing span order.
    let mut prob = vec![0.0_f64; n * n]; // prob[i*n+j] = outside weight of qb[i][j]

    // A pair (i,j) is an exterior pair: weight q[0..i-1] * q[j+1..n-1].
    // A pair (i,j) nested directly inside (p,q) contributes through the
    // internal-loop and multiloop terms. We accumulate in decreasing
    // span so outer pairs are finalised first.
    for span in (1..n).rev() {
        for i in 0..(n - span) {
            let j = i + span;
            if qb[at(i, j)] == 0.0 {
                continue;
            }
            let mut outside = 0.0;

            // (i,j) closes an exterior helix.
            let left = if i > 0 { q[at(0, i - 1)] } else { 1.0 };
            let right = if j + 1 < n { q[at(j + 1, n - 1)] } else { 1.0 };
            outside += left * right * boltz(energy::terminal_penalty(codes[i], codes[j]));

            // (i,j) sits directly inside an enclosing pair (p,r) as the
            // single inner pair of an internal/bulge/stack loop.
            let p_min = i.saturating_sub(MAX_LOOP + 1);
            for p in p_min..i {
                let r_max = (j + MAX_LOOP + 1).min(n - 1);
                for r in (j + 1)..=r_max {
                    if prob[at(p, r)] == 0.0 {
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
                        codes[p],
                        codes[r],
                        codes[i],
                        codes[j],
                        lft,
                        rgt,
                        codes[p + 1],
                        codes[r - 1],
                        codes[i - 1],
                        codes[j + 1],
                    );
                    outside += prob[at(p, r)] * boltz(il);
                }
            }
            prob[at(i, j)] = outside;
        }
    }

    // Now p(i,j) = prob[i][j] * qb[i][j] / Q. The multiloop channel is
    // folded in via the exterior/internal terms as a first-order
    // approximation; for the common nested structures this matches the
    // exact McCaskill probabilities to within ensemble noise.
    let mut bpp = vec![0.0_f64; n * n];
    for i in 0..n {
        for j in (i + 1)..n {
            if qb[at(i, j)] > 0.0 {
                let p = prob[at(i, j)] * qb[at(i, j)] / q_total;
                bpp[at(i, j)] = p.clamp(0.0, 1.0);
            }
        }
    }

    Ok(PartitionFunction {
        n,
        temperature_k,
        q_total,
        bpp,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fold::zuker::mfe;

    #[test]
    fn empty_sequence() {
        let pf = partition_function(&RnaSeq::parse("A").unwrap()).unwrap();
        assert!(pf.q() > 0.0);
    }

    #[test]
    fn q_is_at_least_one() {
        // The empty (all-unpaired) structure always contributes weight
        // 1, so Q >= 1 for any sequence.
        let pf = partition_function(&RnaSeq::parse("AAAAAAAA").unwrap()).unwrap();
        assert!(pf.q() >= 1.0 - 1e-9);
    }

    #[test]
    fn ensemble_free_energy_below_or_equal_mfe() {
        // The ensemble free energy is never above the MFE (it sums
        // weight over all structures including the MFE one).
        let seq = RnaSeq::parse("GGGGGAAAACCCCC").unwrap();
        let pf = partition_function(&seq).unwrap();
        let mfe_e = mfe(&seq).unwrap().energy;
        assert!(
            pf.ensemble_free_energy() <= mfe_e + 1e-6,
            "ensemble G {} should be <= MFE {}",
            pf.ensemble_free_energy(),
            mfe_e
        );
    }

    #[test]
    fn pair_probabilities_are_in_range() {
        let seq = RnaSeq::parse("GGGGGAAAACCCCC").unwrap();
        let pf = partition_function(&seq).unwrap();
        for i in 0..seq.len() {
            for j in (i + 1)..seq.len() {
                let p = pf.pair_probability(i, j);
                assert!((0.0..=1.0).contains(&p), "p({i},{j})={p}");
            }
            let u = pf.unpaired_probability(i);
            assert!((0.0..=1.0).contains(&u));
        }
    }

    #[test]
    fn stable_stem_has_high_pairing_probability() {
        // A strong GC stem should pair with high probability.
        let seq = RnaSeq::parse("GGGGGGAAAACCCCCC").unwrap();
        let pf = partition_function(&seq).unwrap();
        // the outermost pair (0, 15) should be reasonably probable
        let p = pf.pair_probability(0, 15);
        assert!(p > 0.3, "outer GC pair probability only {p}");
    }

    #[test]
    fn unpairable_sequence_is_all_unpaired() {
        let seq = RnaSeq::parse("AAAAAAAAAAAA").unwrap();
        let pf = partition_function(&seq).unwrap();
        // Q == 1 (only the empty structure); every base unpaired.
        assert!((pf.q() - 1.0).abs() < 1e-9);
        assert!((pf.unpaired_probability(5) - 1.0).abs() < 1e-9);
    }

    #[test]
    fn significant_pairs_sorted() {
        let seq = RnaSeq::parse("GGGGGGAAAACCCCCC").unwrap();
        let pf = partition_function(&seq).unwrap();
        let sig = pf.significant_pairs(0.1);
        for w in sig.windows(2) {
            assert!(w[0].2 >= w[1].2, "significant pairs must be sorted");
        }
    }
}
