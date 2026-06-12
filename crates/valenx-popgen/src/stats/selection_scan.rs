//! Haplotype-based selection scans: EHH, iHH and iHS.
//!
//! A recent selective sweep leaves a long stretch of *homozygous*
//! haplotype around the swept allele — selection carries one
//! haplotype to high frequency faster than recombination can break it
//! up. EHH-class statistics detect that signature.
//!
//! - **EHH** ([`ehh`]) — Extended Haplotype Homozygosity: the
//!   probability that two haplotypes carrying a given core allele are
//!   identical from the core out to a test site (Sabeti et al. 2002).
//!   It decays from 1 at the core toward 0 with distance; the decay is
//!   *slow* for a swept allele.
//! - **iHH** ([`integrated_ehh`]) — the integral of EHH over genomic
//!   distance, the area under the EHH-decay curve.
//! - **iHS** ([`ihs`]) — the standardised log-ratio of the iHH of the
//!   ancestral versus derived core allele (Voight et al. 2006). A large
//!   `|iHS|` flags a site whose one allele sits on unusually long
//!   haplotypes.
//!
//! ## v1 simplifications
//!
//! Haplotypes are the rows of a [`GenotypeMatrix`]; "identical" means
//! identical allelic states over the spanned columns. iHS here is the
//! *un-genome-standardised* log(iHH_A / iHH_D) — genome-wide
//! frequency-binned standardisation (the final Voight step) needs a
//! whole-genome scan and is left to the caller.

use crate::error::{PopgenError, Result};
use crate::infer::GenotypeMatrix;

/// EHH for one core allele, evaluated at every site of the matrix.
///
/// Returns `(position, ehh)` pairs. `core` is the column index of the
/// core site; `core_allele` (`0` or `1`) selects which haplotype set to
/// follow. EHH at the core itself is `1.0` by definition; at a far site
/// it is the fraction of carrier-haplotype pairs still identical from
/// the core to there.
///
/// # Errors
/// [`PopgenError::Invalid`] if `core` is out of range, `core_allele` is
/// not `0`/`1`, or fewer than two haplotypes carry the core allele.
pub fn ehh(matrix: &GenotypeMatrix, core: usize, core_allele: u8) -> Result<Vec<(f64, f64)>> {
    if core >= matrix.n_sites() {
        return Err(PopgenError::invalid("core", "site index out of range"));
    }
    if core_allele > 1 {
        return Err(PopgenError::invalid("core_allele", "must be 0 or 1"));
    }
    // Haplotype rows carrying the core allele.
    let carriers: Vec<usize> = (0..matrix.n_samples())
        .filter(|&r| matrix.get(r, core) == core_allele)
        .collect();
    if carriers.len() < 2 {
        return Err(PopgenError::invalid(
            "core_allele",
            "fewer than two haplotypes carry the core allele",
        ));
    }
    let mut out = Vec::with_capacity(matrix.n_sites());
    // For each test site compute EHH over the span [core, test].
    for (test, &pos) in matrix.positions().iter().enumerate() {
        let (lo, hi) = if core <= test {
            (core, test)
        } else {
            (test, core)
        };
        out.push((pos, homozygosity(matrix, &carriers, lo, hi)));
    }
    Ok(out)
}

/// The haplotype homozygosity of the carrier set over columns
/// `[lo, hi]`: the probability two random carriers are identical there.
///
/// Computed by partitioning the carriers into distinct-haplotype groups
/// and summing `C(group, 2) / C(total, 2)`.
fn homozygosity(matrix: &GenotypeMatrix, carriers: &[usize], lo: usize, hi: usize) -> f64 {
    use std::collections::HashMap;
    let mut groups: HashMap<Vec<u8>, usize> = HashMap::new();
    for &r in carriers {
        let hap: Vec<u8> = (lo..=hi).map(|c| matrix.get(r, c)).collect();
        *groups.entry(hap).or_insert(0) += 1;
    }
    let total = carriers.len();
    if total < 2 {
        return 1.0;
    }
    let pairs = (total * (total - 1) / 2) as f64;
    let same: f64 = groups
        .values()
        .map(|&g| (g * (g.saturating_sub(1)) / 2) as f64)
        .sum();
    same / pairs
}

/// Integrated EHH (iHH): the area under the EHH-vs-position curve for a
/// core allele, by the trapezoidal rule.
///
/// # Errors
/// See [`ehh`].
pub fn integrated_ehh(matrix: &GenotypeMatrix, core: usize, core_allele: u8) -> Result<f64> {
    let curve = ehh(matrix, core, core_allele)?;
    Ok(trapezoid(&curve))
}

/// The integrated haplotype score (iHS) at a core site.
///
/// `iHS = ln(iHH_ancestral / iHH_derived)`. A large positive value
/// means the *ancestral* allele sits on long haplotypes; a large
/// negative value means the *derived* allele does (the classic sweep
/// signal). Returns `0.0` if either allele has fewer than two carriers.
///
/// # Errors
/// [`PopgenError::Invalid`] if `core` is out of range.
pub fn ihs(matrix: &GenotypeMatrix, core: usize) -> Result<f64> {
    if core >= matrix.n_sites() {
        return Err(PopgenError::invalid("core", "site index out of range"));
    }
    let ihh_anc = integrated_ehh(matrix, core, 0);
    let ihh_der = integrated_ehh(matrix, core, 1);
    match (ihh_anc, ihh_der) {
        (Ok(a), Ok(d)) if a > 1e-12 && d > 1e-12 => Ok((a / d).ln()),
        _ => Ok(0.0),
    }
}

/// Trapezoidal integral of a `(x, y)` curve, sorted by `x`.
fn trapezoid(points: &[(f64, f64)]) -> f64 {
    let mut pts = points.to_vec();
    pts.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    let mut area = 0.0;
    for w in pts.windows(2) {
        let (x0, y0) = w[0];
        let (x1, y1) = w[1];
        area += (x1 - x0) * (y0 + y1) / 2.0;
    }
    area
}

#[cfg(test)]
mod tests {
    use super::*;

    fn matrix(rows: Vec<Vec<u8>>, positions: Vec<f64>) -> GenotypeMatrix {
        GenotypeMatrix::from_rows(rows, positions).unwrap()
    }

    #[test]
    fn ehh_is_one_at_the_core() {
        let m = matrix(
            vec![vec![1, 1, 0], vec![1, 0, 1], vec![0, 1, 1], vec![0, 0, 0]],
            vec![0.0, 100.0, 200.0],
        );
        let curve = ehh(&m, 0, 1).unwrap();
        // EHH at the core column is 1.
        let at_core = curve.iter().find(|&&(p, _)| p == 0.0).unwrap().1;
        assert!((at_core - 1.0).abs() < 1e-12);
    }

    #[test]
    fn ehh_decays_away_from_the_core_when_haplotypes_differ() {
        // Carriers of the derived core allele have divergent flanks.
        let m = matrix(
            vec![vec![1, 0, 0], vec![1, 1, 1], vec![0, 0, 0], vec![0, 1, 1]],
            vec![0.0, 50.0, 100.0],
        );
        let curve = ehh(&m, 0, 1).unwrap();
        let at_core = curve.iter().find(|&&(p, _)| p == 0.0).unwrap().1;
        let far = curve.iter().find(|&&(p, _)| p == 100.0).unwrap().1;
        assert!(far <= at_core, "EHH rose away from the core");
    }

    #[test]
    fn swept_allele_keeps_high_ehh() {
        // The derived core carriers share an identical long haplotype
        // -> EHH stays 1 across the whole span.
        let m = matrix(
            vec![
                vec![1, 1, 1, 1],
                vec![1, 1, 1, 1],
                vec![1, 1, 1, 1],
                vec![0, 0, 1, 0],
            ],
            vec![0.0, 100.0, 200.0, 300.0],
        );
        let curve = ehh(&m, 0, 1).unwrap();
        for &(_, e) in &curve {
            assert!((e - 1.0).abs() < 1e-12, "EHH dropped: {e}");
        }
    }

    #[test]
    fn integrated_ehh_is_larger_for_a_longer_haplotype() {
        // Allele 1: identical flanks (long haplotype). Allele 0:
        // divergent flanks (short).
        let m = matrix(
            vec![vec![1, 1, 1], vec![1, 1, 1], vec![0, 1, 0], vec![0, 0, 1]],
            vec![0.0, 100.0, 200.0],
        );
        let ihh_der = integrated_ehh(&m, 0, 1).unwrap();
        let ihh_anc = integrated_ehh(&m, 0, 0).unwrap();
        assert!(
            ihh_der > ihh_anc,
            "derived iHH {ihh_der} not above ancestral {ihh_anc}"
        );
    }

    #[test]
    fn ihs_is_negative_when_derived_haplotype_is_long() {
        let m = matrix(
            vec![vec![1, 1, 1], vec![1, 1, 1], vec![0, 1, 0], vec![0, 0, 1]],
            vec![0.0, 100.0, 200.0],
        );
        let score = ihs(&m, 0).unwrap();
        // ln(iHH_anc / iHH_der) < 0 when derived haplotypes are longer.
        assert!(score < 0.0, "iHS = {score}");
    }

    #[test]
    fn ehh_rejects_bad_input() {
        let m = matrix(vec![vec![1], vec![0]], vec![0.0]);
        assert!(ehh(&m, 9, 1).is_err());
        assert!(ehh(&m, 0, 3).is_err());
        // Only one carrier of allele 1 -> error.
        assert!(ehh(&m, 0, 1).is_err());
    }
}
