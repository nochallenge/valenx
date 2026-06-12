//! A simple Bayesian diploid genotype-likelihood model.
//!
//! Given the read bases observed at one site, a genotype caller must
//! decide between the three diploid genotypes of a biallelic locus —
//! homozygous-reference (`0/0`), heterozygous (`0/1`) and
//! homozygous-alternate (`1/1`). This module implements the standard
//! site-independent model used by samtools / GATK in its simplest
//! form: each read base is an independent Bernoulli observation whose
//! error probability comes from its Phred quality, and the genotype
//! likelihood is the product over reads.
//!
//! The math (for one read base `b`, a genotype with reference-allele
//! fraction `f ∈ {0, 0.5, 1}`, and per-base error `e`):
//!
//! ```text
//! P(b = ref | g) = f·(1 − e)  +  (1 − f)·(e / 3)
//! P(b = alt | g) = (1 − f)·(1 − e)  +  f·(e / 3)
//! ```
//!
//! The `e / 3` term spreads the error mass over the three non-true
//! bases. Likelihoods are accumulated in log-space for stability and
//! combined with a prior to yield posteriors and a phred-scaled `GQ`.

/// The three diploid genotypes of a biallelic site.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash)]
pub enum Genotype {
    /// `0/0` — homozygous reference.
    HomRef,
    /// `0/1` — heterozygous.
    Het,
    /// `1/1` — homozygous alternate.
    HomAlt,
}

impl Genotype {
    /// The reference-allele fraction of this genotype (`1.0`, `0.5`,
    /// `0.0`).
    pub fn ref_fraction(self) -> f64 {
        match self {
            Genotype::HomRef => 1.0,
            Genotype::Het => 0.5,
            Genotype::HomAlt => 0.0,
        }
    }

    /// The VCF `GT` string for this genotype.
    pub fn gt_string(self) -> &'static str {
        match self {
            Genotype::HomRef => "0/0",
            Genotype::Het => "0/1",
            Genotype::HomAlt => "1/1",
        }
    }

    /// All three genotypes, in `0/0`, `0/1`, `1/1` order.
    pub fn all() -> [Genotype; 3] {
        [Genotype::HomRef, Genotype::Het, Genotype::HomAlt]
    }
}

/// One read observation at a site: which allele it supports and its
/// Phred base quality.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct AlleleObs {
    /// `true` when the read base equals the reference allele.
    pub is_ref: bool,
    /// Phred quality of the base (capped internally at 60).
    pub quality: u8,
}

/// Converts a Phred quality to an error probability, clamped to a sane
/// `[1e-6, 0.75]` range so a single perfect or junk base cannot
/// dominate or zero out a likelihood.
fn error_prob(q: u8) -> f64 {
    let q = q.min(60) as f64;
    10f64.powf(-q / 10.0).clamp(1e-6, 0.75)
}

/// The result of genotyping one site.
#[derive(Clone, Debug, PartialEq)]
pub struct GenotypeCall {
    /// The most-probable genotype.
    pub best: Genotype,
    /// Posterior log10 probability of each genotype, in
    /// [`Genotype::all`] order.
    pub log10_posteriors: [f64; 3],
    /// Phred-scaled genotype quality: `−10·log10(1 − P(best))`,
    /// capped at 99.
    pub gq: u8,
    /// Phred-scaled likelihoods normalised so the best is `0` (the VCF
    /// `PL` field), in [`Genotype::all`] order.
    pub pl: [i32; 3],
}

/// Genotypes a site from its allele observations.
///
/// `priors` is the prior probability of each genotype in
/// [`Genotype::all`] order; pass [`default_priors`] for the common
/// site-frequency-neutral choice. An empty observation list yields a
/// `HomRef` call with zero confidence.
pub fn genotype_site(obs: &[AlleleObs], priors: [f64; 3]) -> GenotypeCall {
    // log10 likelihood of each genotype.
    let mut loglik = [0.0f64; 3];
    for (gi, &g) in Genotype::all().iter().enumerate() {
        let f = g.ref_fraction();
        let mut acc = 0.0f64;
        for o in obs {
            let e = error_prob(o.quality);
            let p_ref = f * (1.0 - e) + (1.0 - f) * (e / 3.0);
            let p_alt = (1.0 - f) * (1.0 - e) + f * (e / 3.0);
            let p = if o.is_ref { p_ref } else { p_alt };
            acc += p.max(1e-300).log10();
        }
        loglik[gi] = acc;
    }

    // Posterior ∝ likelihood · prior; normalise in log10 space.
    let mut logpost = [0.0f64; 3];
    for i in 0..3 {
        logpost[i] = loglik[i] + priors[i].max(1e-300).log10();
    }
    let max_lp = logpost.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let denom: f64 = logpost.iter().map(|&lp| 10f64.powf(lp - max_lp)).sum();
    let log_denom = max_lp + denom.log10();
    for lp in &mut logpost {
        *lp -= log_denom;
    }

    // Best genotype.
    let mut best_idx = 0usize;
    for i in 1..3 {
        if logpost[i] > logpost[best_idx] {
            best_idx = i;
        }
    }
    let best = Genotype::all()[best_idx];

    // GQ = -10 log10(1 - P(best)).
    let p_best = 10f64.powf(logpost[best_idx]).clamp(0.0, 1.0);
    let gq = if p_best >= 1.0 {
        99
    } else {
        let q = -10.0 * (1.0 - p_best).max(1e-10).log10();
        q.round().clamp(0.0, 99.0) as u8
    };

    // PL = phred-scaled likelihoods, normalised so best = 0.
    let max_ll = loglik.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let mut pl = [0i32; 3];
    for i in 0..3 {
        // loglik is log10; phred = -10 * log10(L). Normalise.
        let phred = -10.0 * (loglik[i] - max_ll);
        pl[i] = phred.round().clamp(0.0, 255.0) as i32;
    }

    GenotypeCall {
        best,
        log10_posteriors: logpost,
        gq,
        pl,
    }
}

/// The site-frequency-neutral default prior — a mild bias toward the
/// reference (`0.5 / 0.3 / 0.2`), a reasonable starting point when no
/// population allele frequency is known.
pub fn default_priors() -> [f64; 3] {
    [0.5, 0.3, 0.2]
}

/// A flat (uninformative) prior — every genotype equally likely.
pub fn flat_priors() -> [f64; 3] {
    [1.0 / 3.0; 3]
}

#[cfg(test)]
mod tests {
    use super::*;

    fn obs(n_ref: usize, n_alt: usize, q: u8) -> Vec<AlleleObs> {
        let mut v = Vec::new();
        for _ in 0..n_ref {
            v.push(AlleleObs {
                is_ref: true,
                quality: q,
            });
        }
        for _ in 0..n_alt {
            v.push(AlleleObs {
                is_ref: false,
                quality: q,
            });
        }
        v
    }

    #[test]
    fn all_ref_calls_homref() {
        let call = genotype_site(&obs(20, 0, 35), default_priors());
        assert_eq!(call.best, Genotype::HomRef);
        assert!(call.gq > 30);
    }

    #[test]
    fn all_alt_calls_homalt() {
        let call = genotype_site(&obs(0, 20, 35), default_priors());
        assert_eq!(call.best, Genotype::HomAlt);
        assert!(call.gq > 30);
    }

    #[test]
    fn balanced_calls_het() {
        let call = genotype_site(&obs(15, 15, 35), default_priors());
        assert_eq!(call.best, Genotype::Het);
    }

    #[test]
    fn pl_of_best_is_zero() {
        let call = genotype_site(&obs(20, 0, 35), default_priors());
        let best_idx = Genotype::all()
            .iter()
            .position(|&g| g == call.best)
            .unwrap();
        assert_eq!(call.pl[best_idx], 0);
    }

    #[test]
    fn posteriors_sum_to_one() {
        let call = genotype_site(&obs(10, 10, 30), flat_priors());
        let sum: f64 = call.log10_posteriors.iter().map(|&lp| 10f64.powf(lp)).sum();
        assert!((sum - 1.0).abs() < 1e-6, "sum = {sum}");
    }

    #[test]
    fn single_low_quality_base_is_uncertain() {
        // One alt base at quality 2 (high error) -> stays HomRef but
        // not max confidence.
        let mut o = obs(8, 0, 35);
        o.push(AlleleObs {
            is_ref: false,
            quality: 2,
        });
        let call = genotype_site(&o, default_priors());
        assert_eq!(call.best, Genotype::HomRef);
    }

    #[test]
    fn empty_observations_are_safe() {
        let call = genotype_site(&[], default_priors());
        // With no data the prior wins (HomRef has the highest prior).
        assert_eq!(call.best, Genotype::HomRef);
    }

    #[test]
    fn higher_depth_raises_confidence() {
        let low = genotype_site(&obs(5, 0, 30), default_priors());
        let high = genotype_site(&obs(50, 0, 30), default_priors());
        assert!(high.gq >= low.gq);
    }
}
