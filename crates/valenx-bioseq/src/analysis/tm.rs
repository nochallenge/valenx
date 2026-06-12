//! Oligonucleotide melting temperature (Tm).
//!
//! Two models are provided:
//!
//! - [`tm_wallace`] — the Wallace ("2+4") rule, `Tm = 2·(A+T) +
//!   4·(G+C)`. A quick estimate, valid only for short oligos (<14 nt).
//! - [`tm_nearest_neighbor`] — the SantaLucia (1998) unified
//!   nearest-neighbor thermodynamic model with monovalent + Mg²⁺ +
//!   dNTP salt corrections, terminal-AT/GC initiation, and a symmetry
//!   correction for self-complementary duplexes. The accurate choice
//!   for primer design.
//!
//! The nearest-neighbor model is implemented by
//! [`crate::analysis::thermo`]; this module exposes it under the
//! historical `tm_*` names.

use crate::analysis::thermo::{self, SaltParams, StrandParams};
use crate::error::{BioseqError, Result};
use crate::seq::{Seq, SeqKind};

/// Wallace-rule melting temperature in °C.
///
/// `Tm = 2·(nA + nT) + 4·(nG + nC)`. Ambiguity codes are ignored.
/// Returns [`BioseqError::Invalid`] for a non-DNA sequence.
pub fn tm_wallace(seq: &Seq) -> Result<f64> {
    if seq.kind() != SeqKind::Dna {
        return Err(BioseqError::invalid("kind", "Tm needs a DNA sequence"));
    }
    let a = seq.count(b'A');
    let t = seq.count(b'T');
    let g = seq.count(b'G');
    let c = seq.count(b'C');
    Ok(2.0 * (a + t) as f64 + 4.0 * (g + c) as f64)
}

/// Parameters for [`tm_nearest_neighbor`].
#[derive(Copy, Clone, Debug)]
pub struct NnParams {
    /// Total strand concentration, mol/L. Typical PCR primer: 0.25 µM
    /// = `0.25e-6`.
    pub strand_conc: f64,
    /// Monovalent cation (Na⁺ / K⁺) concentration, mol/L. Typical
    /// PCR buffer: 0.05 M.
    pub na_conc: f64,
    /// Divalent Mg²⁺ concentration, mol/L. Typical PCR: 1.5 mM.
    pub mg_conc: f64,
    /// Total dNTP concentration, mol/L. Typical PCR: 0.2 mM. dNTPs
    /// chelate Mg²⁺, so [Mg²⁺]_free is taken as max(0, mg - dntp).
    pub dntp_conc: f64,
}

impl Default for NnParams {
    fn default() -> Self {
        NnParams {
            strand_conc: 0.25e-6,
            na_conc: 0.05,
            mg_conc: 1.5e-3,
            dntp_conc: 0.2e-3,
        }
    }
}

impl From<NnParams> for (SaltParams, StrandParams) {
    fn from(p: NnParams) -> Self {
        (
            SaltParams {
                na_conc: p.na_conc,
                mg_conc: p.mg_conc,
                dntp_conc: p.dntp_conc,
            },
            StrandParams {
                strand_conc: p.strand_conc,
            },
        )
    }
}

/// SantaLucia nearest-neighbor melting temperature in °C.
///
/// Implemented by [`crate::analysis::thermo::duplex_tm`]: sums the 10
/// nearest-neighbor ΔH°/ΔS° increments, adds the initiation terms
/// (terminal A/T vs. G/C), applies the self-complementary symmetry
/// correction when applicable, and folds the monovalent + Mg²⁺ + dNTP
/// salt environment into the entropy via the von Ahsen 1999
/// effective-monovalent formula. Returns [`BioseqError::Invalid`] for
/// a non-DNA sequence, a sequence shorter than 2 nt, or any base that
/// is not a canonical `ACGT` (the model is undefined for ambiguity
/// codes).
pub fn tm_nearest_neighbor(seq: &Seq, params: NnParams) -> Result<f64> {
    if seq.kind() != SeqKind::Dna {
        return Err(BioseqError::invalid("kind", "Tm needs a DNA sequence"));
    }
    let bytes = seq.as_bytes();
    if bytes.len() < 2 {
        return Err(BioseqError::invalid(
            "sequence",
            "nearest-neighbor Tm needs at least 2 nt",
        ));
    }
    for &b in bytes {
        if !matches!(b, b'A' | b'C' | b'G' | b'T') {
            return Err(BioseqError::invalid(
                "sequence",
                format!(
                    "non-canonical base `{}` — NN model needs pure ACGT",
                    b as char
                ),
            ));
        }
    }
    let (salt, strand) = params.into();
    thermo::duplex_tm(bytes, salt, strand)
        .ok_or_else(|| BioseqError::invalid("params", "degenerate thermodynamic denominator"))
}

/// Nearest-neighbor Tm with the default PCR-primer parameters.
pub fn tm_nearest_neighbor_default(seq: &Seq) -> Result<f64> {
    tm_nearest_neighbor(seq, NnParams::default())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn wallace_rule_arithmetic() {
        // 4 G/C + 4 A/T -> 4*4 + 4*2 = 24.
        let s = Seq::new(SeqKind::Dna, "GGCCAATT").unwrap();
        assert!((tm_wallace(&s).unwrap() - 24.0).abs() < 1e-12);
    }

    #[test]
    fn wallace_all_gc() {
        let s = Seq::new(SeqKind::Dna, "GGGG").unwrap();
        assert!((tm_wallace(&s).unwrap() - 16.0).abs() < 1e-12);
    }

    #[test]
    fn wallace_rejects_non_dna() {
        let r = Seq::new(SeqKind::Rna, "ACGU").unwrap();
        assert!(tm_wallace(&r).is_err());
    }

    #[test]
    fn nn_tm_is_reasonable_for_a_typical_primer() {
        // A standard ~20-mer primer melts in the 50-65 °C range.
        let s = Seq::new(SeqKind::Dna, "ACGTACGTACGTACGTACGT").unwrap();
        let tm = tm_nearest_neighbor_default(&s).unwrap();
        assert!(tm > 30.0 && tm < 80.0, "Tm out of plausible range: {tm}");
    }

    #[test]
    fn nn_tm_gc_rich_higher_than_at_rich() {
        let gc = Seq::new(SeqKind::Dna, "GCGCGCGCGCGCGCGC").unwrap();
        let at = Seq::new(SeqKind::Dna, "ATATATATATATATAT").unwrap();
        let tm_gc = tm_nearest_neighbor_default(&gc).unwrap();
        let tm_at = tm_nearest_neighbor_default(&at).unwrap();
        assert!(
            tm_gc > tm_at,
            "GC-rich Tm {tm_gc} should exceed AT-rich {tm_at}"
        );
    }

    #[test]
    fn nn_tm_longer_oligo_higher() {
        let short = Seq::new(SeqKind::Dna, "ACGTACGTAC").unwrap();
        let long = Seq::new(SeqKind::Dna, "ACGTACGTACGTACGTACGTACGT").unwrap();
        assert!(
            tm_nearest_neighbor_default(&long).unwrap()
                > tm_nearest_neighbor_default(&short).unwrap()
        );
    }

    #[test]
    fn nn_tm_salt_dependence() {
        let s = Seq::new(SeqKind::Dna, "ACGTACGTACGTACGT").unwrap();
        let low_salt = tm_nearest_neighbor(
            &s,
            NnParams {
                na_conc: 0.01,
                mg_conc: 0.0,
                dntp_conc: 0.0,
                ..Default::default()
            },
        )
        .unwrap();
        let high_salt = tm_nearest_neighbor(
            &s,
            NnParams {
                na_conc: 0.2,
                mg_conc: 0.0,
                dntp_conc: 0.0,
                ..Default::default()
            },
        )
        .unwrap();
        // Higher monovalent salt stabilizes the duplex -> higher Tm.
        assert!(high_salt > low_salt, "low {low_salt} high {high_salt}");
    }

    #[test]
    fn nn_tm_magnesium_raises_tm() {
        let s = Seq::new(SeqKind::Dna, "ACGTACGTACGTACGT").unwrap();
        let no_mg = tm_nearest_neighbor(
            &s,
            NnParams {
                na_conc: 0.05,
                mg_conc: 0.0,
                dntp_conc: 0.0,
                ..Default::default()
            },
        )
        .unwrap();
        let with_mg = tm_nearest_neighbor(
            &s,
            NnParams {
                na_conc: 0.05,
                mg_conc: 3e-3,
                dntp_conc: 0.2e-3,
                ..Default::default()
            },
        )
        .unwrap();
        assert!(
            with_mg > no_mg,
            "Mg should raise Tm; no={no_mg} with={with_mg}"
        );
    }

    #[test]
    fn nn_tm_rejects_short_and_ambiguous() {
        let one = Seq::new(SeqKind::Dna, "A").unwrap();
        assert!(tm_nearest_neighbor_default(&one).is_err());
        let amb = Seq::new(SeqKind::Dna, "ACGTN").unwrap();
        assert!(tm_nearest_neighbor_default(&amb).is_err());
    }
}
