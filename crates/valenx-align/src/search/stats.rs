//! Karlin-Altschul alignment statistics — E-value and bit-score.
//!
//! The score of a *local* ungapped alignment between random sequences
//! follows an extreme-value (Gumbel) distribution. Karlin & Altschul
//! (1990) showed the expected number of distinct high-scoring segment
//! pairs (HSPs) with score ≥ `S` in a comparison of an `m`-residue
//! query against an `n`-residue database is
//!
//! ```text
//! E = K · m · n · exp(-λ · S)
//! ```
//!
//! and the *bit score* — a score normalised so it is comparable across
//! scoring systems — is
//!
//! ```text
//! S' = (λ · S − ln K) / ln 2
//! ```
//!
//! [`KarlinAltschul`] bundles the two statistical parameters `λ` and
//! `K` and provides the conversions. The `p`-value (probability of at
//! least one such HSP) follows from the Poisson approximation
//! `p = 1 − exp(−E)`.
//!
//! ## v1 scope
//!
//! `λ` and `K` are normally derived by solving the Karlin-Altschul
//! equations for a given scoring matrix and residue composition. This
//! crate ships the published NCBI BLAST values for the common
//! ungapped scoring systems and the standard *gapped* presets for
//! BLOSUM62; an arbitrary scoring system can also have its parameters
//! supplied directly via [`KarlinAltschul::new`].

/// The two Karlin-Altschul statistical parameters plus the natural
/// scale of the scoring system.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct KarlinAltschul {
    /// The scale parameter λ (`lambda`) — natural-log units per score
    /// unit. Always positive.
    pub lambda: f64,
    /// The `K` parameter — a search-space correction factor in `(0,1]`.
    pub k: f64,
    /// Relative entropy `H` of the scoring system in nats (used for
    /// effective-length corrections). Optional context; not required
    /// for the E-value formula.
    pub h: f64,
}

impl KarlinAltschul {
    /// Builds a parameter set from explicit `lambda`, `k` and `h`.
    pub fn new(lambda: f64, k: f64, h: f64) -> Self {
        KarlinAltschul { lambda, k, h }
    }

    /// Published NCBI parameters for **ungapped** BLOSUM62
    /// (`λ = 0.318`, `K = 0.134`, `H = 0.40`).
    // `0.318` is the published Karlin-Altschul λ for ungapped
    // BLOSUM62; clippy mistakes its proximity to 1/π — it is not.
    #[allow(clippy::approx_constant)]
    pub fn blosum62_ungapped() -> Self {
        KarlinAltschul::new(0.318, 0.134, 0.40)
    }

    /// Published NCBI parameters for **gapped** BLOSUM62 with the
    /// default `11/1` gap costs (`λ = 0.267`, `K = 0.041`, `H = 0.14`).
    pub fn blosum62_gapped() -> Self {
        KarlinAltschul::new(0.267, 0.041, 0.14)
    }

    /// Published NCBI parameters for an ungapped `+1 / −3` DNA scoring
    /// system (`λ ≈ 1.374`, `K ≈ 0.711`, `H ≈ 1.31`) — the classic
    /// `megablast` defaults.
    pub fn dna_ungapped() -> Self {
        KarlinAltschul::new(1.374, 0.711, 1.31)
    }

    /// The **bit score** of a raw alignment score `s`.
    ///
    /// `S' = (λ·s − ln K) / ln 2`. Bit scores are comparable across
    /// scoring systems and database sizes.
    pub fn bit_score(&self, s: i32) -> f64 {
        (self.lambda * s as f64 - self.k.ln()) / std::f64::consts::LN_2
    }

    /// The **E-value** of a raw score `s` for an `m`-residue query
    /// against an `n`-residue database (effective lengths).
    ///
    /// `E = K · m · n · exp(−λ·s)`.
    pub fn e_value(&self, s: i32, m: usize, n: usize) -> f64 {
        self.k * m as f64 * n as f64 * (-self.lambda * s as f64).exp()
    }

    /// The E-value computed directly from a bit score and a search
    /// space size — `E = m·n · 2^(−S')`. Equivalent to
    /// [`e_value`](Self::e_value) but parameter-free given the bit
    /// score.
    pub fn e_value_from_bits(bit_score: f64, m: usize, n: usize) -> f64 {
        (m as f64) * (n as f64) * 2f64.powf(-bit_score)
    }

    /// The `p`-value — probability of at least one HSP scoring ≥ `s` —
    /// from the Poisson approximation `p = 1 − e^{−E}`. Always in
    /// `[0, 1]`.
    pub fn p_value(&self, s: i32, m: usize, n: usize) -> f64 {
        let e = self.e_value(s, m, n);
        1.0 - (-e).exp()
    }

    /// The minimum raw score whose E-value is at most `target_e` for
    /// the given search space — the score cutoff a search would apply.
    ///
    /// Inverts the E-value formula:
    /// `s = ceil( (ln(K·m·n) − ln E) / λ )`.
    pub fn score_threshold(&self, target_e: f64, m: usize, n: usize) -> i32 {
        if target_e <= 0.0 {
            return i32::MAX;
        }
        let s = ((self.k * m as f64 * n as f64).ln() - target_e.ln()) / self.lambda;
        s.ceil() as i32
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bit_score_monotone() {
        let ka = KarlinAltschul::blosum62_ungapped();
        let lo = ka.bit_score(20);
        let hi = ka.bit_score(60);
        assert!(hi > lo, "higher raw score => higher bit score");
    }

    #[test]
    fn e_value_decreases_with_score() {
        let ka = KarlinAltschul::blosum62_gapped();
        let e_low = ka.e_value(30, 300, 1_000_000);
        let e_high = ka.e_value(120, 300, 1_000_000);
        assert!(e_high < e_low, "stronger hit => smaller E-value");
        assert!(e_high >= 0.0);
    }

    #[test]
    fn e_value_scales_with_search_space() {
        let ka = KarlinAltschul::blosum62_gapped();
        let small = ka.e_value(50, 100, 1000);
        let big = ka.e_value(50, 100, 1_000_000);
        // 1000x bigger DB => ~1000x bigger E-value.
        assert!((big / small - 1000.0).abs() < 1.0);
    }

    #[test]
    fn p_value_in_unit_interval() {
        let ka = KarlinAltschul::blosum62_ungapped();
        for s in [10, 30, 60, 100] {
            let p = ka.p_value(s, 250, 500_000);
            assert!((0.0..=1.0).contains(&p), "p={p} out of range for s={s}");
        }
    }

    #[test]
    fn bits_and_evalue_consistent() {
        let ka = KarlinAltschul::blosum62_gapped();
        let (s, m, n) = (80, 300, 2_000_000);
        let e_direct = ka.e_value(s, m, n);
        let e_via_bits = KarlinAltschul::e_value_from_bits(ka.bit_score(s), m, n);
        // The two routes must agree closely.
        assert!((e_direct - e_via_bits).abs() / e_direct < 1e-6);
    }

    #[test]
    fn score_threshold_inverts_evalue() {
        let ka = KarlinAltschul::blosum62_gapped();
        let (m, n) = (300, 1_000_000);
        let thr = ka.score_threshold(0.01, m, n);
        // A score at the threshold has E-value <= the target.
        assert!(ka.e_value(thr, m, n) <= 0.01 + 1e-9);
        // One below the threshold exceeds it.
        assert!(ka.e_value(thr - 1, m, n) > 0.01);
    }

    #[test]
    fn dna_params_load() {
        let ka = KarlinAltschul::dna_ungapped();
        assert!(ka.lambda > 1.0);
        assert!(ka.k > 0.0 && ka.k <= 1.0);
    }
}
