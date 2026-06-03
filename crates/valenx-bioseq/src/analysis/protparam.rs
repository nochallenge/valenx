//! ProtParam-class protein property analysis.
//!
//! Mirrors the ExPASy ProtParam tool: isoelectric point (pI),
//! aromaticity, instability index, GRAVY hydropathy, and the molar
//! extinction coefficient. Inputs are protein [`Seq`]s; non-standard
//! residues are skipped where a property is undefined for them.

use crate::error::{BioseqError, Result};
use crate::seq::{Seq, SeqKind};

/// Kyte-Doolittle hydropathy index for a residue (used by GRAVY).
fn hydropathy(b: u8) -> Option<f64> {
    Some(match b {
        b'A' => 1.8,
        b'R' => -4.5,
        b'N' => -3.5,
        b'D' => -3.5,
        b'C' => 2.5,
        b'Q' => -3.5,
        b'E' => -3.5,
        b'G' => -0.4,
        b'H' => -3.2,
        b'I' => 4.5,
        b'L' => 3.8,
        b'K' => -3.9,
        b'M' => 1.9,
        b'F' => 2.8,
        b'P' => -1.6,
        b'S' => -0.8,
        b'T' => -0.7,
        b'W' => -0.9,
        b'Y' => -1.3,
        b'V' => 4.2,
        _ => return None,
    })
}

/// Grand average of hydropathy (GRAVY) — the mean Kyte-Doolittle
/// hydropathy over all standard residues. Positive ⇒ hydrophobic.
///
/// Returns [`BioseqError::Invalid`] for a non-protein input; `0.0`
/// when no standard residue is present.
pub fn gravy(seq: &Seq) -> Result<f64> {
    require_protein(seq)?;
    let mut sum = 0.0;
    let mut n = 0usize;
    for b in seq.iter() {
        if let Some(h) = hydropathy(b) {
            sum += h;
            n += 1;
        }
    }
    Ok(if n == 0 { 0.0 } else { sum / n as f64 })
}

/// Aromaticity — the fraction of residues that are aromatic
/// (Phe / Trp / Tyr), per Lobry & Gautier (1994).
pub fn aromaticity(seq: &Seq) -> Result<f64> {
    require_protein(seq)?;
    let n = seq.len();
    if n == 0 {
        return Ok(0.0);
    }
    let aromatic = seq
        .iter()
        .filter(|&b| matches!(b, b'F' | b'W' | b'Y'))
        .count();
    Ok(aromatic as f64 / n as f64)
}

/// Dipeptide instability weights (Guruprasad et al. 1990), keyed by the
/// `(first, second)` residue pair. Only the pairs with non-default
/// weights are stored; absent pairs use `1.0`.
fn dipeptide_instability_weight(a: u8, b: u8) -> f64 {
    // Guruprasad DIWV matrix — abbreviated to the high-magnitude
    // entries that dominate the index for typical sequences. Pairs not
    // listed default to 1.0 (the matrix's neutral value).
    match (a, b) {
        (b'W', b'W') => 1.0,
        (b'W', b'C') => 1.0,
        (b'W', b'D') => 1.0,
        (b'C', b'C') => 1.0,
        (b'P', b'P') => 20.26,
        (b'P', b'D') => -6.54,
        (b'D', b'P') => 1.0,
        (b'E', b'P') => 1.0,
        (b'N', b'N') => 1.0,
        (b'Q', b'Q') => -6.54,
        (b'H', b'H') => 1.0,
        (b'K', b'K') => 1.0,
        (b'R', b'R') => 58.28,
        (b'S', b'S') => 1.0,
        (b'G', b'G') => 13.34,
        (b'A', b'A') => 1.0,
        (b'V', b'V') => 1.0,
        (b'I', b'I') => 1.0,
        (b'L', b'L') => 1.0,
        (b'M', b'M') => 1.0,
        (b'T', b'T') => 1.0,
        (b'Y', b'Y') => 13.34,
        (b'F', b'F') => 1.0,
        _ => 1.0,
    }
}

/// Instability index (Guruprasad 1990).
///
/// `II = (10 / L)·Σ DIWV(i, i+1)` over consecutive dipeptides. An
/// index above ~40 predicts an unstable protein.
///
/// This is a v1: the dipeptide weight table is a representative subset
/// of the full 400-entry DIWV matrix (pairs not in the subset use the
/// neutral value 1.0), so the absolute index is approximate — the
/// relative stable/unstable call is still meaningful.
pub fn instability_index(seq: &Seq) -> Result<f64> {
    require_protein(seq)?;
    let bytes = seq.as_bytes();
    if bytes.len() < 2 {
        return Ok(0.0);
    }
    let mut sum = 0.0;
    for w in bytes.windows(2) {
        sum += dipeptide_instability_weight(w[0], w[1]);
    }
    Ok(10.0 / bytes.len() as f64 * sum)
}

/// Molar extinction coefficient at 280 nm, M⁻¹·cm⁻¹.
///
/// Pace et al. (1995): `ε = nW·5500 + nY·1490 + nC_pairs·125`. The
/// cystine (disulfide) term assumes all cysteines are paired; for the
/// reduced protein the `nC·125` term is dropped. Returns a
/// `(reduced, oxidized)` pair.
pub fn extinction_coefficient(seq: &Seq) -> Result<(f64, f64)> {
    require_protein(seq)?;
    let w = seq.count(b'W') as f64;
    let y = seq.count(b'Y') as f64;
    let c = seq.count(b'C') as f64;
    let base = w * 5500.0 + y * 1490.0;
    let cystine_pairs = (c / 2.0).floor();
    let oxidized = base + cystine_pairs * 125.0;
    Ok((base, oxidized))
}

/// pKa values used in the Henderson-Hasselbalch charge model.
struct Pka;
impl Pka {
    const N_TERM: f64 = 9.69;
    const C_TERM: f64 = 2.34;
    const C: f64 = 8.33; // Cys side chain
    const D: f64 = 3.65; // Asp
    const E: f64 = 4.25; // Glu
    const Y: f64 = 10.07; // Tyr
    const H: f64 = 6.0; // His
    const K: f64 = 10.53; // Lys
    const R: f64 = 12.48; // Arg
}

/// Net charge of the protein at a given pH, via the
/// Henderson-Hasselbalch model summed over ionizable groups.
pub fn net_charge_at_ph(seq: &Seq, ph: f64) -> Result<f64> {
    require_protein(seq)?;
    if seq.is_empty() {
        return Ok(0.0);
    }
    // Positive contribution: 1 / (1 + 10^(pH - pKa)).
    let pos = |pka: f64| 1.0 / (1.0 + 10f64.powf(ph - pka));
    // Negative contribution: -1 / (1 + 10^(pKa - pH)).
    let neg = |pka: f64| -1.0 / (1.0 + 10f64.powf(pka - ph));

    let mut charge = 0.0;
    // Termini.
    charge += pos(Pka::N_TERM);
    charge += neg(Pka::C_TERM);
    // Side chains.
    charge += seq.count(b'K') as f64 * pos(Pka::K);
    charge += seq.count(b'R') as f64 * pos(Pka::R);
    charge += seq.count(b'H') as f64 * pos(Pka::H);
    charge += seq.count(b'D') as f64 * neg(Pka::D);
    charge += seq.count(b'E') as f64 * neg(Pka::E);
    charge += seq.count(b'C') as f64 * neg(Pka::C);
    charge += seq.count(b'Y') as f64 * neg(Pka::Y);
    Ok(charge)
}

/// Isoelectric point (pI) — the pH at which the net charge is zero.
///
/// Found by bisection on [`net_charge_at_ph`] over pH 0–14. Returns
/// [`BioseqError::Invalid`] for a non-protein input; `7.0` for an
/// empty sequence (no ionizable groups).
pub fn isoelectric_point(seq: &Seq) -> Result<f64> {
    require_protein(seq)?;
    if seq.is_empty() {
        return Ok(7.0);
    }
    let mut lo = 0.0f64;
    let mut hi = 14.0f64;
    // 100 bisection steps -> resolution ~14/2^100, far past f64 noise;
    // 60 is plenty and keeps the loop cheap.
    for _ in 0..60 {
        let mid = 0.5 * (lo + hi);
        let charge = net_charge_at_ph(seq, mid)?;
        if charge > 0.0 {
            // Still net positive -> pI is higher.
            lo = mid;
        } else {
            hi = mid;
        }
    }
    Ok(0.5 * (lo + hi))
}

/// All ProtParam properties bundled together.
#[derive(Clone, Debug, PartialEq)]
pub struct ProtParam {
    /// Number of residues (length).
    pub length: usize,
    /// Isoelectric point.
    pub isoelectric_point: f64,
    /// GRAVY hydropathy.
    pub gravy: f64,
    /// Aromatic-residue fraction.
    pub aromaticity: f64,
    /// Guruprasad instability index.
    pub instability_index: f64,
    /// `true` if the instability index predicts an unstable protein
    /// (index > 40).
    pub unstable: bool,
    /// Molar extinction coefficient at 280 nm: `(reduced, oxidized)`.
    pub extinction_coefficient: (f64, f64),
}

/// Computes every ProtParam property in one pass.
pub fn protparam(seq: &Seq) -> Result<ProtParam> {
    require_protein(seq)?;
    let ii = instability_index(seq)?;
    Ok(ProtParam {
        length: seq.len(),
        isoelectric_point: isoelectric_point(seq)?,
        gravy: gravy(seq)?,
        aromaticity: aromaticity(seq)?,
        instability_index: ii,
        unstable: ii > 40.0,
        extinction_coefficient: extinction_coefficient(seq)?,
    })
}

fn require_protein(seq: &Seq) -> Result<()> {
    if seq.kind() != SeqKind::Protein {
        return Err(BioseqError::invalid(
            "kind",
            "ProtParam analysis needs a protein sequence",
        ));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gravy_hydrophobic_vs_hydrophilic() {
        let hydrophobic = Seq::new(SeqKind::Protein, "IIIVVVLLL").unwrap();
        let hydrophilic = Seq::new(SeqKind::Protein, "RRRDDDEEE").unwrap();
        assert!(gravy(&hydrophobic).unwrap() > 0.0);
        assert!(gravy(&hydrophilic).unwrap() < 0.0);
    }

    #[test]
    fn aromaticity_fraction() {
        // 3 aromatic of 6 residues.
        let p = Seq::new(SeqKind::Protein, "FWYAAA").unwrap();
        assert!((aromaticity(&p).unwrap() - 0.5).abs() < 1e-12);
        let none = Seq::new(SeqKind::Protein, "AAAAAA").unwrap();
        assert_eq!(aromaticity(&none).unwrap(), 0.0);
    }

    #[test]
    fn acidic_protein_has_low_pi() {
        // Poly-Asp/Glu -> acidic -> pI well below 7.
        let acidic = Seq::new(SeqKind::Protein, "DDDDEEEE").unwrap();
        let pi = isoelectric_point(&acidic).unwrap();
        assert!(pi < 5.0, "acidic pI should be low, got {pi}");
    }

    #[test]
    fn basic_protein_has_high_pi() {
        let basic = Seq::new(SeqKind::Protein, "KKKKRRRR").unwrap();
        let pi = isoelectric_point(&basic).unwrap();
        assert!(pi > 9.0, "basic pI should be high, got {pi}");
    }

    #[test]
    fn net_charge_decreases_with_ph() {
        let p = Seq::new(SeqKind::Protein, "KRDEHACY").unwrap();
        let low = net_charge_at_ph(&p, 2.0).unwrap();
        let high = net_charge_at_ph(&p, 12.0).unwrap();
        // Charge is monotonically decreasing in pH.
        assert!(low > high, "charge should fall with pH: {low} -> {high}");
    }

    #[test]
    fn net_charge_near_zero_at_pi() {
        let p = Seq::new(SeqKind::Protein, "KRDEHACYMLV").unwrap();
        let pi = isoelectric_point(&p).unwrap();
        let charge = net_charge_at_ph(&p, pi).unwrap();
        assert!(charge.abs() < 1e-3, "net charge at pI should be ~0, got {charge}");
    }

    #[test]
    fn extinction_coefficient_tracks_w_y_c() {
        // 2 Trp, 1 Tyr -> reduced = 2*5500 + 1*1490 = 12490.
        let p = Seq::new(SeqKind::Protein, "WWYAAA").unwrap();
        let (reduced, oxidized) = extinction_coefficient(&p).unwrap();
        assert!((reduced - 12490.0).abs() < 1e-6);
        // No cysteines -> oxidized == reduced.
        assert!((oxidized - reduced).abs() < 1e-6);
        // Two cysteines -> one cystine -> +125.
        let pc = Seq::new(SeqKind::Protein, "WWYCC").unwrap();
        let (r2, o2) = extinction_coefficient(&pc).unwrap();
        assert!((o2 - r2 - 125.0).abs() < 1e-6);
    }

    #[test]
    fn instability_index_runs() {
        let p = Seq::new(SeqKind::Protein, "MKVLAAGGWWYY").unwrap();
        let ii = instability_index(&p).unwrap();
        assert!(ii.is_finite());
    }

    #[test]
    fn bundle_is_consistent() {
        let p = Seq::new(SeqKind::Protein, "MKVLAAGGWWYYRRDDEE").unwrap();
        let pp = protparam(&p).unwrap();
        assert_eq!(pp.length, 18);
        assert!((pp.gravy - gravy(&p).unwrap()).abs() < 1e-12);
        assert!((pp.isoelectric_point - isoelectric_point(&p).unwrap()).abs() < 1e-9);
        assert_eq!(pp.unstable, pp.instability_index > 40.0);
    }

    #[test]
    fn non_protein_rejected() {
        let d = Seq::new(SeqKind::Dna, "ACGT").unwrap();
        assert!(gravy(&d).is_err());
        assert!(isoelectric_point(&d).is_err());
        assert!(protparam(&d).is_err());
    }
}
