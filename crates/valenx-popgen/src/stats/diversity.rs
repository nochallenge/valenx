//! Diversity and neutrality-test statistics.
//!
//! These are the classic single-population summaries, all computable
//! from a [`GenotypeMatrix`] (or, equivalently, its
//! [`crate::stats::Sfs`]):
//!
//! - **Nucleotide diversity pi** ([`nucleotide_diversity`]) — the mean
//!   number of pairwise differences per site.
//! - **Watterson's theta** ([`wattersons_theta`]) — `S / a_n`, where
//!   `S` is the segregating-site count and `a_n` the harmonic number.
//! - **Tajima's D** ([`tajimas_d`]) — the normalised difference
//!   `pi - theta_W`; near 0 under neutrality, negative under an
//!   excess of rare variants (a sweep or expansion), positive under an
//!   excess of intermediate variants (balancing selection, structure).
//! - **Fu & Li's D** ([`fu_li_d`]) — contrasts singleton variation
//!   with total variation; sensitive to recent deleterious mutations.
//! - **Fay & Wu's H** ([`fay_wu_h`]) — contrasts high-frequency
//!   derived variation with intermediate; strongly negative right
//!   after a selective sweep.
//!
//! Fu & Li's D and Fay & Wu's H need the *ancestral* allele known;
//! they assume the matrix's `1` allele is derived.

use crate::error::{PopgenError, Result};
use crate::infer::GenotypeMatrix;

/// `a_n = sum_{i=1}^{n-1} 1/i` — the first harmonic number, the
/// expected coalescent tree length in units of `2N`.
fn a1(n: usize) -> f64 {
    (1..n).map(|i| 1.0 / i as f64).sum()
}

/// `a_{2,n} = sum_{i=1}^{n-1} 1/i^2`.
fn a2(n: usize) -> f64 {
    (1..n).map(|i| 1.0 / (i * i) as f64).sum()
}

/// Mean number of pairwise differences across the sample, summed over
/// all sites (the un-normalised `pi`, often written `k`).
///
/// For a site with `d` derived and `n - d` ancestral copies, the number
/// of differing pairs is `d * (n - d)`, and there are `C(n,2)` pairs,
/// so the site contributes `d(n-d) / C(n,2)`.
///
/// # Errors
/// [`PopgenError::Invalid`] if the matrix has fewer than two samples.
pub fn pairwise_differences(matrix: &GenotypeMatrix) -> Result<f64> {
    let n = matrix.n_samples();
    if n < 2 {
        return Err(PopgenError::invalid(
            "matrix",
            "need at least two samples for diversity",
        ));
    }
    let pairs = (n * (n - 1) / 2) as f64;
    let mut total = 0.0;
    for col in 0..matrix.n_sites() {
        let d = matrix.derived_count(col)? as f64;
        total += d * (n as f64 - d) / pairs;
    }
    Ok(total)
}

/// Nucleotide diversity `pi`: the mean pairwise difference **per
/// site**, dividing [`pairwise_differences`] by the number of sites.
///
/// If `n_sites` is 0 the result is 0.
///
/// # Errors
/// [`PopgenError::Invalid`] if the matrix has fewer than two samples.
pub fn nucleotide_diversity(matrix: &GenotypeMatrix) -> Result<f64> {
    let k = pairwise_differences(matrix)?;
    let sites = matrix.n_sites();
    Ok(if sites == 0 { 0.0 } else { k / sites as f64 })
}

/// Watterson's estimator of `theta`: `S / a_n`, where `S` is the number
/// of segregating sites and `a_n` the harmonic number.
///
/// This is the *total* (segment-wide) estimate; divide by the segment
/// length for a per-site value.
///
/// # Errors
/// [`PopgenError::Invalid`] if the matrix has fewer than two samples.
pub fn wattersons_theta(matrix: &GenotypeMatrix) -> Result<f64> {
    let n = matrix.n_samples();
    if n < 2 {
        return Err(PopgenError::invalid(
            "matrix",
            "need at least two samples for theta",
        ));
    }
    let s = matrix.segregating_sites() as f64;
    Ok(s / a1(n))
}

/// Tajima's D: the difference between `pi` and Watterson's `theta`,
/// normalised by its standard deviation (Tajima 1989).
///
/// `D ~ 0` under neutral equilibrium, `< 0` with an excess of rare
/// alleles, `> 0` with an excess of intermediate-frequency alleles.
/// Returns `0.0` when there is no variance to normalise by (e.g. no
/// segregating sites).
///
/// # Errors
/// [`PopgenError::Invalid`] if the matrix has fewer than four samples
/// (the variance terms are undefined below `n = 4`).
pub fn tajimas_d(matrix: &GenotypeMatrix) -> Result<f64> {
    let n = matrix.n_samples();
    if n < 4 {
        return Err(PopgenError::invalid(
            "matrix",
            "Tajima's D needs at least four samples",
        ));
    }
    let s = matrix.segregating_sites() as f64;
    if s == 0.0 {
        return Ok(0.0);
    }
    let nn = n as f64;
    let a_1 = a1(n);
    let a_2 = a2(n);
    // Tajima 1989 variance coefficients.
    let b1 = (nn + 1.0) / (3.0 * (nn - 1.0));
    let b2 = 2.0 * (nn * nn + nn + 3.0) / (9.0 * nn * (nn - 1.0));
    let c1 = b1 - 1.0 / a_1;
    let c2 = b2 - (nn + 2.0) / (a_1 * nn) + a_2 / (a_1 * a_1);
    let e1 = c1 / a_1;
    let e2 = c2 / (a_1 * a_1 + a_2);
    let variance = e1 * s + e2 * s * (s - 1.0);
    if variance <= 0.0 {
        return Ok(0.0);
    }
    let pi = pairwise_differences(matrix)?;
    let theta_w = s / a_1;
    Ok((pi - theta_w) / variance.sqrt())
}

/// Fu & Li's D (Fu & Li 1993): contrasts the number of singletons
/// `eta_e` with the total number of segregating sites, normalised.
///
/// Strongly negative when recent (typically deleterious) mutations
/// inflate the singleton class.
///
/// # Errors
/// [`PopgenError::Invalid`] if the matrix has fewer than four samples.
pub fn fu_li_d(matrix: &GenotypeMatrix) -> Result<f64> {
    let n = matrix.n_samples();
    if n < 4 {
        return Err(PopgenError::invalid(
            "matrix",
            "Fu & Li's D needs at least four samples",
        ));
    }
    let s = matrix.segregating_sites() as f64;
    if s == 0.0 {
        return Ok(0.0);
    }
    let nn = n as f64;
    let a_1 = a1(n);
    let a_2 = a2(n);
    // Count derived singletons (the external-branch mutations).
    let mut eta_e = 0.0;
    for col in 0..matrix.n_sites() {
        if matrix.derived_count(col)? == 1 {
            eta_e += 1.0;
        }
    }
    // Fu & Li 1993 coefficients.
    let c = if (n - 1) == 0 {
        0.0
    } else {
        2.0 * (nn * a_1 - 2.0 * (nn - 1.0)) / ((nn - 1.0) * (nn - 2.0))
    };
    let v = 1.0 + (a_1 * a_1 / (a_2 + a_1 * a_1)) * (c - (nn + 1.0) / (nn - 1.0));
    let u = a_1 - 1.0 - v;
    let variance = u * s + v * s * s;
    if variance <= 0.0 {
        return Ok(0.0);
    }
    Ok((s - a_1 * eta_e) / variance.sqrt())
}

/// Fay & Wu's H (Fay & Wu 2000): the difference between `pi` and
/// `theta_H`, an estimator that weights high-frequency derived variants
/// heavily.
///
/// `theta_H = sum_i (2 * xi_i * i^2) / (n(n-1))`. A sweep leaves an
/// excess of high-frequency derived alleles, driving H sharply
/// negative.
///
/// This returns the *un-normalised* `pi - theta_H` (the classic Fay &
/// Wu statistic; the normalised version requires an outgroup-error
/// model out of v1 scope).
///
/// # Errors
/// [`PopgenError::Invalid`] if the matrix has fewer than two samples.
pub fn fay_wu_h(matrix: &GenotypeMatrix) -> Result<f64> {
    let n = matrix.n_samples();
    if n < 2 {
        return Err(PopgenError::invalid(
            "matrix",
            "Fay & Wu's H needs at least two samples",
        ));
    }
    let nn = n as f64;
    let mut theta_h = 0.0;
    for col in 0..matrix.n_sites() {
        let i = matrix.derived_count(col)? as f64;
        if i > 0.0 && i < nn {
            theta_h += 2.0 * i * i;
        }
    }
    theta_h /= nn * (nn - 1.0);
    let pi = pairwise_differences(matrix)?;
    Ok(pi - theta_h)
}

/// Expected heterozygosity (gene diversity) `He` — the mean over sites of `2·p·(1−p)`, where
/// `p` is the derived-allele frequency at each site. Measures the diversity expected under
/// random mating from allele frequencies alone; distinct from nucleotide diversity π (the mean
/// pairwise difference per site). Returns `0.0` when there are no sites.
///
/// # Errors
/// Propagates [`PopgenError`] if a site's allele frequency cannot be computed.
pub fn expected_heterozygosity(matrix: &GenotypeMatrix) -> Result<f64> {
    let sites = matrix.n_sites();
    if sites == 0 {
        return Ok(0.0);
    }
    let mut sum = 0.0;
    for col in 0..sites {
        let p = matrix.frequency(col)?;
        sum += 2.0 * p * (1.0 - p);
    }
    Ok(sum / sites as f64)
}

/// Mean minor-allele frequency (MAF) — the mean over sites of `min(p, 1−p)`, where `p` is the
/// derived-allele frequency at each site. Measures the mean burden of the rarer allele; distinct
/// from expected heterozygosity (which uses the product 2·p·(1−p)) and from π. Returns `0.0`
/// when there are no sites.
///
/// # Errors
/// Propagates [`PopgenError`] if a site's allele frequency cannot be computed.
pub fn minor_allele_frequency(matrix: &GenotypeMatrix) -> Result<f64> {
    let sites = matrix.n_sites();
    if sites == 0 {
        return Ok(0.0);
    }
    let mut sum = 0.0;
    for col in 0..sites {
        let p = matrix.frequency(col)?;
        sum += p.min(1.0 - p);
    }
    Ok(sum / sites as f64)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn matrix(rows: Vec<Vec<u8>>) -> GenotypeMatrix {
        let cols = rows[0].len();
        let pos: Vec<f64> = (0..cols).map(|i| i as f64).collect();
        GenotypeMatrix::from_rows(rows, pos).unwrap()
    }

    #[test]
    fn pairwise_differences_on_a_known_site() {
        // 4 samples, one site with 2 derived: differing pairs =
        // 2*2 = 4, over C(4,2) = 6 -> 0.6667.
        let m = matrix(vec![vec![1], vec![1], vec![0], vec![0]]);
        let k = pairwise_differences(&m).unwrap();
        assert!((k - 4.0 / 6.0).abs() < 1e-12);
    }

    #[test]
    fn wattersons_theta_is_s_over_harmonic() {
        // 4 samples, 2 segregating sites. a_4 = 1 + 1/2 + 1/3.
        let m = matrix(vec![vec![1, 1], vec![1, 0], vec![0, 0], vec![0, 0]]);
        let theta = wattersons_theta(&m).unwrap();
        let a = 1.0 + 0.5 + 1.0 / 3.0;
        assert!((theta - 2.0 / a).abs() < 1e-12);
    }

    #[test]
    fn nucleotide_diversity_is_per_site() {
        let m = matrix(vec![vec![1, 1], vec![1, 0], vec![0, 0], vec![0, 0]]);
        let pi = nucleotide_diversity(&m).unwrap();
        let k = pairwise_differences(&m).unwrap();
        assert!((pi - k / 2.0).abs() < 1e-12);
    }

    #[test]
    fn tajimas_d_is_near_zero_for_a_balanced_spectrum() {
        // A spectrum with one singleton and one doubleton in n=4 is
        // roughly neutral-shaped; D should be small in magnitude.
        let m = matrix(vec![vec![1, 1], vec![0, 1], vec![0, 0], vec![0, 0]]);
        let d = tajimas_d(&m).unwrap();
        assert!(d.abs() < 3.0, "D = {d}");
    }

    #[test]
    fn tajimas_d_is_negative_with_excess_singletons() {
        // Many singletons, no intermediate variants -> pi < theta_W
        // -> D negative.
        let mut rows = vec![vec![0u8; 8]; 8];
        for (i, row) in rows.iter_mut().enumerate().take(8) {
            row[i] = 1; // each site a private singleton
        }
        let m = matrix(rows);
        let d = tajimas_d(&m).unwrap();
        assert!(d < 0.0, "expected negative D, got {d}");
    }

    #[test]
    fn tajimas_d_zero_when_no_segregating_sites() {
        let m = matrix(vec![vec![0], vec![0], vec![0], vec![0]]);
        assert_eq!(tajimas_d(&m).unwrap(), 0.0);
    }

    #[test]
    fn fu_li_d_runs_and_is_finite() {
        let mut rows = vec![vec![0u8; 6]; 6];
        for (i, row) in rows.iter_mut().enumerate().take(4) {
            row[i] = 1;
        }
        rows[0][4] = 1;
        rows[1][4] = 1; // a doubleton
        let m = matrix(rows);
        let d = fu_li_d(&m).unwrap();
        assert!(d.is_finite());
    }

    #[test]
    fn fay_wu_h_is_negative_with_high_frequency_derived() {
        // 6 samples, sites with 5 derived copies (high frequency).
        let m = matrix(vec![
            vec![1, 1],
            vec![1, 1],
            vec![1, 1],
            vec![1, 1],
            vec![1, 1],
            vec![0, 0],
        ]);
        let h = fay_wu_h(&m).unwrap();
        assert!(h < 0.0, "expected negative H, got {h}");
    }

    #[test]
    fn rejects_too_small_samples() {
        let m = matrix(vec![vec![1], vec![0]]);
        assert!(tajimas_d(&m).is_err());
        assert!(fu_li_d(&m).is_err());
        // pi and theta only need n >= 2.
        assert!(nucleotide_diversity(&m).is_ok());
    }

    #[test]
    fn expected_heterozygosity_from_allele_frequencies() {
        // 2 samples, 1 site, p = 0.5 → He = 2·0.5·0.5 = 0.5.
        let m = matrix(vec![vec![1], vec![0]]);
        assert!((expected_heterozygosity(&m).unwrap() - 0.5).abs() < 1e-12);
        // monomorphic site (all ancestral) → He = 0.
        let m0 = matrix(vec![vec![0], vec![0], vec![0]]);
        assert_eq!(expected_heterozygosity(&m0).unwrap(), 0.0);
        // 4 samples, 2 sites: site0 p=0.5 → 0.5, site1 p=0.75 → 0.375; mean = 0.4375.
        let m2 = matrix(vec![vec![1, 1], vec![1, 1], vec![0, 1], vec![0, 0]]);
        assert!((expected_heterozygosity(&m2).unwrap() - 0.4375).abs() < 1e-12);
    }

    #[test]
    fn minor_allele_frequency_is_mean_of_min() {
        // 2 samples, 1 site, p = 0.5 → min(0.5, 0.5) = 0.5.
        let m = matrix(vec![vec![1], vec![0]]);
        assert!((minor_allele_frequency(&m).unwrap() - 0.5).abs() < 1e-12);
        // 4 samples, 2 sites both p = 0.75 → min(0.75, 0.25) = 0.25.
        let m2 = matrix(vec![vec![1, 1], vec![1, 1], vec![1, 1], vec![0, 0]]);
        assert!((minor_allele_frequency(&m2).unwrap() - 0.25).abs() < 1e-12);
        // Distinct from He: for p = 0.75, He = 0.375 but MAF = 0.25.
        assert!(
            (minor_allele_frequency(&m2).unwrap() - expected_heterozygosity(&m2).unwrap()).abs()
                > 0.1
        );
    }

    #[test]
    fn tajimas_d_matches_the_hand_computed_value() {
        // GROUND TRUTH (Tajima 1989). A 4-sample matrix with 3 segregating
        // sites of derived counts 1, 2, 1 — every quantity hand-computed from
        // first principles. The existing Tajima tests only check the SIGN
        // (near-zero / negative / zero); this pins the exact statistic, so it
        // validates the variance-coefficient pipeline (c1, c2, e1, e2):
        //   pairs C(4,2)=6;  pi = (1·3 + 2·2 + 1·3)/6 = 10/6 = 5/3
        //   S = 3;  a1 = 1+1/2+1/3 = 11/6;  a2 = 1+1/4+1/9 = 49/36
        //   b1 = (n+1)/(3(n−1)) = 5/9;  b2 = 2(n²+n+3)/(9n(n−1)) = 23/54
        //   c1 = b1 − 1/a1 = 0.0101010;  c2 = b2 − (n+2)/(a1·n) + a2/a1² = 0.0127028
        //   e1 = c1/a1 = 0.0055096;  e2 = c2/(a1²+a2) = 0.0026899
        //   Var = e1·S + e2·S(S−1) = 0.0165289 + 0.0161395 = 0.0326684
        //   theta_w = S/a1 = 18/11;  D = (pi − theta_w)/√Var = 0.030303/0.180744 = 0.16765
        let m = matrix(vec![
            vec![1, 1, 0],
            vec![0, 1, 0],
            vec![0, 0, 1],
            vec![0, 0, 0],
        ]);
        // Independent component checks (separate code paths), exact:
        assert_eq!(m.segregating_sites(), 3);
        assert!((pairwise_differences(&m).unwrap() - 5.0 / 3.0).abs() < 1e-12);
        assert!((wattersons_theta(&m).unwrap() - 18.0 / 11.0).abs() < 1e-12);
        // The full statistic, pinned to the first-principles hand value.
        let d = tajimas_d(&m).unwrap();
        assert!((d - 0.16765).abs() < 1e-4, "Tajima's D {d} != 0.16765 (hand)");
    }
}
