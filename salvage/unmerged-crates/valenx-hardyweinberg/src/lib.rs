//! # valenx-hardyweinberg
//!
//! ## What
//!
//! Closed-form population-genetics calculators for a single di-allelic
//! locus under the Hardy-Weinberg model:
//!
//! 1. estimate allele frequencies `p` and `q` from genotype counts
//!    ([`AlleleFreq::from_genotype_counts`]),
//! 2. predict the equilibrium genotype frequencies `p^2`, `2pq`, `q^2`
//!    ([`GenotypeFreq`]), and
//! 3. test observed genotype counts against those expectations with a
//!    Pearson chi-square goodness-of-fit test
//!    ([`hwe_chi_square_test`]).
//!
//! Every public function validates its inputs and returns
//! [`Result`](error::Result), rejecting non-finite values and
//! out-of-domain arguments (negative counts, frequencies outside
//! `[0, 1]`, empty samples).
//!
//! ## Model
//!
//! For a locus with alleles `A` and `a`, a diploid sample of `N`
//! individuals carries `2N` allele copies. Gene counting gives the
//! maximum-likelihood allele frequencies
//!
//! ```text
//! p = (2 * n_AA + n_Aa) / (2N)
//! q = (2 * n_aa + n_Aa) / (2N)
//! ```
//!
//! and by construction `p + q = 1`. Under the Hardy-Weinberg
//! assumptions (random mating, no selection / mutation / migration,
//! infinite population) the genotype frequencies after one generation
//! are the binomial expansion of `(p + q)^2`:
//!
//! ```text
//! freq(AA) = p^2,  freq(Aa) = 2pq,  freq(aa) = q^2,
//! ```
//!
//! which sum to `1`. The chi-square goodness-of-fit statistic compares
//! observed counts `O` with the expected counts `E = N * (p^2, 2pq, q^2)`:
//!
//! ```text
//! chi2 = sum_i (O_i - E_i)^2 / E_i,   df = (3 - 1) - 1 = 1
//! ```
//!
//! where one degree of freedom is spent estimating `p`. The upper-tail
//! p-value is the regularised upper incomplete gamma `Q(df/2, chi2/2)`,
//! which for `df = 1` equals `erfc(sqrt(chi2 / 2))`.
//!
//! ## Honest scope
//!
//! This is research/educational grade — standard textbook closed-form
//! and numerical models, validated against analytic ground truth
//! (`p = q = 0.5` gives `0.25 / 0.50 / 0.25`; hand-worked chi-square
//! values; the `df = 1` critical value `3.841459` at `alpha = 0.05`).
//! It is NOT a clinical, medical, diagnostic, or production-certified
//! engineering tool: there are no component tolerances, safety factors,
//! fatigue or thermal limits, code-compliance checks, multi-allelic or
//! multi-locus models, exact (Fisher) HWE tests, or corrections for
//! population structure / inbreeding.
//!
//! ## Example
//!
//! ```rust
//! use valenx_hardyweinberg::{AlleleFreq, GenotypeFreq, hwe_chi_square_test};
//!
//! // 35 AA, 20 Aa, 45 aa individuals.
//! let freq = AlleleFreq::from_genotype_counts(35.0, 20.0, 45.0).unwrap();
//! assert!((freq.p - 0.45).abs() < 1e-12);
//! assert!((freq.q - 0.55).abs() < 1e-12);
//!
//! // Equilibrium genotype frequencies sum to 1.
//! let g = GenotypeFreq::from_alleles(&freq).unwrap();
//! assert!((g.total() - 1.0).abs() < 1e-12);
//!
//! // Goodness-of-fit test for Hardy-Weinberg equilibrium.
//! let test = hwe_chi_square_test(30.0, 30.0, 40.0).unwrap();
//! assert_eq!(test.degrees_of_freedom, 1);
//! assert!(test.rejects_hwe(0.05)); // this sample deviates from HWE
//! ```

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod allele;
pub mod chisquare;
pub mod error;
pub mod genotype;

pub use allele::{major_allele_frequency, minor_allele_frequency, AlleleFreq};
pub use chisquare::{chi_square_upper_tail, hwe_chi_square_test, HweTest};
pub use error::{ErrorCategory, HweError};
pub use genotype::GenotypeFreq;
