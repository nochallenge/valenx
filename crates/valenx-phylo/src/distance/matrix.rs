//! Pairwise distance matrices and their estimation from an alignment.
//!
//! A [`DistMatrix`] is a labelled, symmetric, zero-diagonal matrix of
//! evolutionary distances. [`distance_matrix`] estimates one from a
//! [`valenx_align::Msa`] under a chosen [`DistanceModel`]:
//!
//! - **p-distance** — the raw proportion of differing sites. An
//!   observed quantity, not an evolutionary distance; it saturates.
//! - **Jukes-Cantor (JC69)** — corrects p for unseen multiple hits
//!   under equal base frequencies and one substitution rate:
//!   `d = -3/4 · ln(1 - 4p/3)`.
//! - **Kimura 2-parameter (K80)** — separates transitions `P` from
//!   transversions `Q`:
//!   `d = ½·ln(1/(1-2P-Q)) + ¼·ln(1/(1-2Q))`.
//! - **Tamura-Nei (TN93)** — distinguishes the two transition classes
//!   (purine A↔G vs pyrimidine C↔T) and uses empirical base
//!   frequencies; the most general of the four.
//!
//! Only aligned columns where *both* sequences have an unambiguous
//! A/C/G/T are counted for a given pair; gaps and ambiguity codes are
//! skipped pairwise (the standard "pairwise deletion" policy).

use crate::error::{PhyloError, Result};
use valenx_align::Msa;

/// Nucleotide-substitution model used to correct an observed
/// p-distance into an evolutionary distance.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DistanceModel {
    /// Uncorrected proportion of differing sites.
    PDistance,
    /// Jukes-Cantor 1969 single-rate correction.
    JukesCantor,
    /// Kimura 1980 two-parameter (transition / transversion).
    Kimura2P,
    /// Tamura-Nei 1993 (two transition classes + base frequencies).
    TamuraNei,
}

/// A symmetric, zero-diagonal matrix of pairwise distances with one
/// label per row/column.
#[derive(Debug, Clone, PartialEq)]
pub struct DistMatrix {
    labels: Vec<String>,
    /// Row-major `n × n` storage.
    data: Vec<f64>,
}

impl DistMatrix {
    /// Builds a matrix from labels and a flat row-major `n × n` buffer.
    ///
    /// # Errors
    /// [`PhyloError::Dimension`] if `data.len() != labels.len()²`.
    pub fn new(labels: Vec<String>, data: Vec<f64>) -> Result<Self> {
        let n = labels.len();
        if data.len() != n * n {
            return Err(PhyloError::dimension(n * n, data.len(), "distance matrix"));
        }
        Ok(DistMatrix { labels, data })
    }

    /// A zero matrix of size `labels.len()`.
    pub fn zeros(labels: Vec<String>) -> Self {
        let n = labels.len();
        DistMatrix {
            data: vec![0.0; n * n],
            labels,
        }
    }

    /// Number of taxa (matrix order).
    pub fn len(&self) -> usize {
        self.labels.len()
    }

    /// `true` if the matrix has no taxa.
    pub fn is_empty(&self) -> bool {
        self.labels.is_empty()
    }

    /// The taxon labels, in row order.
    pub fn labels(&self) -> &[String] {
        &self.labels
    }

    /// Distance between taxa `i` and `j`.
    ///
    /// # Panics
    /// If either index is out of range.
    pub fn get(&self, i: usize, j: usize) -> f64 {
        self.data[i * self.labels.len() + j]
    }

    /// Sets `d(i,j)` and `d(j,i)` together, keeping the matrix
    /// symmetric.
    ///
    /// # Panics
    /// If either index is out of range.
    pub fn set(&mut self, i: usize, j: usize, value: f64) {
        let n = self.labels.len();
        self.data[i * n + j] = value;
        self.data[j * n + i] = value;
    }

    /// `true` if the matrix is symmetric within `tol` and its diagonal
    /// is zero.
    pub fn is_symmetric(&self, tol: f64) -> bool {
        let n = self.labels.len();
        for i in 0..n {
            if self.get(i, i).abs() > tol {
                return false;
            }
            for j in (i + 1)..n {
                if (self.get(i, j) - self.get(j, i)).abs() > tol {
                    return false;
                }
            }
        }
        true
    }
}

/// Maps a nucleotide byte to a 0..3 index (A,C,G,T), or `None` for a
/// gap / ambiguity / non-nucleotide byte.
fn base_index(b: u8) -> Option<usize> {
    match b.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' | b'U' => Some(3),
        _ => None,
    }
}

/// `true` for the purines A (0) and G (2).
fn is_purine(idx: usize) -> bool {
    idx == 0 || idx == 2
}

/// Estimates a pairwise [`DistMatrix`] from a multiple alignment.
///
/// `labels` supplies one name per alignment row (in row order); its
/// length must equal `msa.depth()`.
///
/// # Errors
/// - [`PhyloError::Dimension`] if `labels.len() != msa.depth()`.
/// - [`PhyloError::Invalid`] if the alignment has fewer than two rows.
///
/// Distances that the model cannot compute (saturated p-distance,
/// negative log argument) are clamped to a large finite value rather
/// than `+∞`, so downstream clustering stays numerically well-behaved.
pub fn distance_matrix(msa: &Msa, labels: &[String], model: DistanceModel) -> Result<DistMatrix> {
    let n = msa.depth();
    if n < 2 {
        return Err(PhyloError::invalid(
            "msa",
            "need at least two sequences for a distance matrix",
        ));
    }
    if labels.len() != n {
        return Err(PhyloError::dimension(n, labels.len(), "alignment rows"));
    }
    let mut out = DistMatrix::zeros(labels.to_vec());
    for i in 0..n {
        for j in (i + 1)..n {
            let d = pairwise_distance(&msa.rows[i], &msa.rows[j], model);
            out.set(i, j, d);
        }
    }
    Ok(out)
}

/// Saturated-distance cap — substitutes for `+∞` when a correction's
/// logarithm argument goes non-positive.
const DIST_CAP: f64 = 10.0;

/// Computes the model distance between two aligned rows.
fn pairwise_distance(a: &[u8], b: &[u8], model: DistanceModel) -> f64 {
    // Tally the 4×4 substitution counts over usable columns.
    let mut counts = [[0u64; 4]; 4];
    let len = a.len().min(b.len());
    for k in 0..len {
        if let (Some(x), Some(y)) = (base_index(a[k]), base_index(b[k])) {
            counts[x][y] += 1;
        }
    }
    let mut total = 0u64;
    let mut transitions = 0u64; // A<->G, C<->T
    let mut transversions = 0u64;
    let mut ag = 0u64; // purine transitions
    let mut ct = 0u64; // pyrimidine transitions
    let mut freq = [0.0f64; 4];
    for x in 0..4 {
        for y in 0..4 {
            let c = counts[x][y];
            total += c;
            freq[x] += c as f64;
            freq[y] += c as f64;
            if x == y {
                continue;
            }
            if is_purine(x) == is_purine(y) {
                transitions += c;
                if is_purine(x) {
                    ag += c;
                } else {
                    ct += c;
                }
            } else {
                transversions += c;
            }
        }
    }
    if total == 0 {
        return DIST_CAP;
    }
    let total_f = total as f64;
    for f in &mut freq {
        *f /= 2.0 * total_f;
    }
    let diffs = transitions + transversions;
    let p = diffs as f64 / total_f;

    match model {
        DistanceModel::PDistance => p,
        DistanceModel::JukesCantor => {
            let arg = 1.0 - 4.0 / 3.0 * p;
            if arg <= 0.0 {
                DIST_CAP
            } else {
                (-0.75 * arg.ln()).min(DIST_CAP)
            }
        }
        DistanceModel::Kimura2P => {
            let pp = transitions as f64 / total_f; // transition fraction
            let q = transversions as f64 / total_f; // transversion fraction
            let a1 = 1.0 - 2.0 * pp - q;
            let a2 = 1.0 - 2.0 * q;
            if a1 <= 0.0 || a2 <= 0.0 {
                DIST_CAP
            } else {
                (-0.5 * a1.ln() - 0.25 * a2.ln()).min(DIST_CAP)
            }
        }
        DistanceModel::TamuraNei => tamura_nei(
            &freq,
            ag as f64 / total_f,
            ct as f64 / total_f,
            transversions as f64 / total_f,
        ),
    }
}

/// Tamura-Nei 1993 distance from empirical base frequencies and the
/// three substitution-class proportions.
///
/// `p1` = purine-transition fraction, `p2` = pyrimidine-transition
/// fraction, `q` = transversion fraction. Uses the standard TN93
/// closed form; clamped to [`DIST_CAP`] when any logarithm argument is
/// non-positive.
fn tamura_nei(freq: &[f64; 4], p1: f64, p2: f64, q: f64) -> f64 {
    let (ga, gc, gg, gt) = (freq[0], freq[1], freq[2], freq[3]);
    let gr = ga + gg; // purine frequency
    let gy = gc + gt; // pyrimidine frequency
    if gr <= 0.0 || gy <= 0.0 {
        return DIST_CAP;
    }
    let k1 = 2.0 * ga * gg / gr;
    let k2 = 2.0 * gc * gt / gy;
    let k3 = 2.0 * (gr * gy - ga * gg * gy / gr - gc * gt * gr / gy);

    let w1 = 1.0 - p1 / k1 - q / (2.0 * gr);
    let w2 = 1.0 - p2 / k2 - q / (2.0 * gy);
    let w3 = 1.0 - q / (2.0 * gr * gy);
    if w1 <= 0.0 || w2 <= 0.0 || w3 <= 0.0 {
        return DIST_CAP;
    }
    let d = -k1 * w1.ln() - k2 * w2.ln() - k3 * w3.ln();
    if d.is_finite() {
        d.clamp(0.0, DIST_CAP)
    } else {
        DIST_CAP
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Three short DNA sequences with known pairwise differences.
    fn sample_msa() -> Msa {
        Msa::new(vec![
            b"ACGTACGTAC".to_vec(),
            b"ACGTACGTAC".to_vec(), // identical to row 0
            b"AGGTACGAAC".to_vec(), // 2 differences vs row 0
        ])
        .unwrap()
    }

    #[test]
    fn p_distance_counts_mismatches() {
        let msa = sample_msa();
        let labels: Vec<String> = ["x", "y", "z"].iter().map(|s| s.to_string()).collect();
        let dm = distance_matrix(&msa, &labels, DistanceModel::PDistance).unwrap();
        // Identical rows => distance 0.
        assert!((dm.get(0, 1)).abs() < 1e-12);
        // 2 / 10 sites differ.
        assert!((dm.get(0, 2) - 0.2).abs() < 1e-12);
        assert!(dm.is_symmetric(1e-9));
    }

    #[test]
    fn jc_distance_exceeds_p_distance() {
        let msa = sample_msa();
        let labels: Vec<String> = ["x", "y", "z"].iter().map(|s| s.to_string()).collect();
        let p = distance_matrix(&msa, &labels, DistanceModel::PDistance).unwrap();
        let jc = distance_matrix(&msa, &labels, DistanceModel::JukesCantor).unwrap();
        // JC correction inflates a non-zero distance.
        assert!(jc.get(0, 2) > p.get(0, 2));
        assert!((jc.get(0, 1)).abs() < 1e-12);
    }

    #[test]
    fn k2p_and_tn93_are_finite_and_nonnegative() {
        let msa = sample_msa();
        let labels: Vec<String> = ["x", "y", "z"].iter().map(|s| s.to_string()).collect();
        for model in [DistanceModel::Kimura2P, DistanceModel::TamuraNei] {
            let dm = distance_matrix(&msa, &labels, model).unwrap();
            for i in 0..3 {
                for j in 0..3 {
                    let d = dm.get(i, j);
                    assert!(d.is_finite() && d >= 0.0, "{model:?} d({i},{j}) = {d}");
                }
            }
        }
    }

    #[test]
    fn saturated_distance_is_capped() {
        // Every site differs => p-distance 1, JC log argument negative.
        let msa = Msa::new(vec![b"AAAA".to_vec(), b"TTTT".to_vec()]).unwrap();
        let labels: Vec<String> = ["a", "b"].iter().map(|s| s.to_string()).collect();
        let jc = distance_matrix(&msa, &labels, DistanceModel::JukesCantor).unwrap();
        assert!(jc.get(0, 1) <= DIST_CAP && jc.get(0, 1).is_finite());
    }

    #[test]
    fn dimension_mismatch_is_rejected() {
        let msa = sample_msa();
        let too_few: Vec<String> = vec!["only-one".to_string()];
        assert!(distance_matrix(&msa, &too_few, DistanceModel::PDistance).is_err());
    }

    // ---- Deeper reference-value tests: exact closed-form distances ---
    //
    // JC69 and K80 have exact analytic forms. Each test below builds a
    // two-sequence alignment with a *precisely known* number of
    // transitions and transversions, computes the published closed-form
    // distance by hand in the comment, and asserts the estimator
    // reproduces it to 1e-9 — a genuine reference-value check, not a
    // property/finiteness check.

    /// Build a 2-row alignment with `a` and `b` as the rows.
    fn pair(a: &[u8], b: &[u8]) -> Msa {
        Msa::new(vec![a.to_vec(), b.to_vec()]).unwrap()
    }

    /// The 2-taxon distance under `model`.
    fn dist(a: &[u8], b: &[u8], model: DistanceModel) -> f64 {
        let labels = vec!["a".to_string(), "b".to_string()];
        distance_matrix(&pair(a, b), &labels, model)
            .unwrap()
            .get(0, 1)
    }

    #[test]
    fn jukes_cantor_matches_the_closed_form_exactly() {
        // 20-site alignment: 16 identical + 4 substitutions (any kind).
        //   p = 4/20 = 0.20
        //   JC69:  d = -3/4 · ln(1 - 4p/3)
        //            = -0.75 · ln(1 - 0.8/3)
        //            = -0.75 · ln(0.733333...)
        //            = 0.232616...  kcal-free, just the distance.
        let a = b"AAAAAAAAAAAAAAAAAAAA";
        // 4 changes at positions 0..4: A->G (transition) is fine, the
        // JC model does not distinguish substitution classes.
        let b = b"GGGGAAAAAAAAAAAAAAAA";
        let p: f64 = 4.0 / 20.0;
        let expected = -0.75 * (1.0 - 4.0 / 3.0 * p).ln();
        let got = dist(a, b, DistanceModel::JukesCantor);
        assert!(
            (got - expected).abs() < 1e-9,
            "JC69 distance {got} != closed form {expected}"
        );
        // sanity: it really is bigger than the raw p-distance.
        assert!(got > p);
    }

    #[test]
    fn kimura_2p_matches_the_closed_form_exactly() {
        // 20-site alignment with a precisely known split:
        //   2 transition sites   (A<->G)  -> P = 2/20 = 0.10
        //   2 transversion sites (A<->C)  -> Q = 2/20 = 0.10
        //   16 identical sites
        // K80:  d = -1/2 · ln(1 - 2P - Q) - 1/4 · ln(1 - 2Q)
        //         = -0.5 · ln(1 - 0.2 - 0.1) - 0.25 · ln(1 - 0.2)
        //         = -0.5 · ln(0.70) - 0.25 · ln(0.80)
        let a = b"AAAAAAAAAAAAAAAAAAAA";
        // positions 0,1: A->G transitions; positions 2,3: A->C transversions
        let b = b"GGCCAAAAAAAAAAAAAAAA";
        let p_ts: f64 = 2.0 / 20.0;
        let q_tv: f64 = 2.0 / 20.0;
        let expected = -0.5 * (1.0 - 2.0 * p_ts - q_tv).ln() - 0.25 * (1.0 - 2.0 * q_tv).ln();
        let got = dist(a, b, DistanceModel::Kimura2P);
        assert!(
            (got - expected).abs() < 1e-9,
            "K80 distance {got} != closed form {expected}"
        );
    }

    #[test]
    fn kimura_2p_distinguishes_transitions_from_transversions() {
        // The defining property of K80: it scores an all-transition
        // alignment DIFFERENTLY from an all-transversion alignment of
        // the same raw divergence, whereas JC69 cannot tell them apart.
        // Both K80 values are checked against the closed form.
        let a = b"AAAAAAAAAAAAAAAAAAAA";
        let all_ts = b"GGGGAAAAAAAAAAAAAAAA"; // 4 transitions  -> P=0.2, Q=0
        let all_tv = b"CCCCAAAAAAAAAAAAAAAA"; // 4 transversions -> P=0,   Q=0.2
        let d_ts = dist(a, all_ts, DistanceModel::Kimura2P);
        let d_tv = dist(a, all_tv, DistanceModel::Kimura2P);

        // Closed form, P=0.2 / Q=0:  d = -1/2·ln(0.6) - 1/4·ln(1.0).
        let ts_expected = -0.5_f64 * (1.0 - 0.4_f64).ln();
        // Closed form, P=0 / Q=0.2:  d = -1/2·ln(0.8) - 1/4·ln(0.6).
        let tv_expected = -0.5_f64 * (1.0 - 0.2_f64).ln() - 0.25_f64 * (1.0 - 0.4_f64).ln();
        assert!(
            (d_ts - ts_expected).abs() < 1e-9,
            "K80 all-transition {d_ts} != closed form {ts_expected}"
        );
        assert!(
            (d_tv - tv_expected).abs() < 1e-9,
            "K80 all-transversion {d_tv} != closed form {tv_expected}"
        );
        // The two are genuinely distinct numbers — K80 sees the class.
        assert!(
            (d_ts - d_tv).abs() > 1e-3,
            "K80 must score transitions ({d_ts}) and transversions ({d_tv}) differently"
        );

        // JC69 cannot tell them apart — both give the identical distance.
        let jc_ts = dist(a, all_ts, DistanceModel::JukesCantor);
        let jc_tv = dist(a, all_tv, DistanceModel::JukesCantor);
        assert!((jc_ts - jc_tv).abs() < 1e-12, "JC69 must ignore the class");
    }

    #[test]
    fn identical_sequences_have_zero_distance_under_every_model() {
        let seq = b"ACGTACGTACGTACGT";
        for model in [
            DistanceModel::PDistance,
            DistanceModel::JukesCantor,
            DistanceModel::Kimura2P,
            DistanceModel::TamuraNei,
        ] {
            let d = dist(seq, seq, model);
            assert!(d.abs() < 1e-12, "{model:?} self-distance {d} should be 0");
        }
    }

    #[test]
    fn jukes_cantor_small_divergence_approaches_p_distance() {
        // For very small p, the JC correction d -> p (the leading term
        // of -3/4·ln(1-4p/3) is exactly p). A single change in 100
        // sites: p = 0.01, and d should be within ~1% of p.
        let mut a = vec![b'A'; 100];
        let mut b = vec![b'A'; 100];
        b[0] = b'G';
        let _ = &mut a;
        let p: f64 = 0.01;
        let d = dist(&a, &b, DistanceModel::JukesCantor);
        assert!(
            (d - p).abs() / p < 0.02,
            "at p={p} the JC distance {d} should be within 2% of p"
        );
        assert!(d >= p, "the correction never shrinks the distance");
    }
}
