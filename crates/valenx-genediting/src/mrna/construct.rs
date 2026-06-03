//! Feature 16 — the mRNA construct model.
//!
//! A therapeutic mRNA is an assembly of five parts, 5′→3′:
//!
//! 1. a **5′ cap** — an `m7G` cap analog that recruits the ribosome
//!    and protects the 5′ end;
//! 2. a **5′UTR** — an untranslated leader carrying the Kozak context;
//! 3. the **CDS** — the open reading frame, `ATG`…stop;
//! 4. a **3′UTR** — an untranslated trailer carrying stability
//!    determinants;
//! 5. a **poly-A tail** — a run of `A`s that sets mRNA half-life.
//!
//! This module gives a [`MrnaConstruct`] data model and a
//! [`MrnaConstructBuilder`] that assembles and **validates** one — the
//! CDS must start with a start codon, end with a stop codon and have a
//! length divisible by three; the cap and tail must be plausible.
//!
//! The construct is stored as RNA (`A C G U`); the builder transcribes
//! DNA input automatically.

use crate::error::{GeneditingError, Result};
use crate::sequtil::{is_acgu, transcribe, upper};
use serde::{Deserialize, Serialize};
use valenx_bioseq::ops::translate::GeneticCode;

/// The 5′ cap chemistry of an mRNA construct.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum CapType {
    /// Cap 0 — `m7GpppN`; an `m7G` cap with no 2′-O-methylation.
    Cap0,
    /// Cap 1 — `m7GpppNm`; 2′-O-methylated first transcribed
    /// nucleotide. The standard for therapeutic mRNA (evades RIG-I /
    /// IFIT innate sensing).
    Cap1,
    /// A co-transcriptional trinucleotide cap analog (CleanCap-class).
    CleanCap,
    /// An anti-reverse cap analog (ARCA) used in older T7 IVT.
    Arca,
}

impl CapType {
    /// Human-readable name.
    pub fn name(self) -> &'static str {
        match self {
            CapType::Cap0 => "Cap 0 (m7GpppN)",
            CapType::Cap1 => "Cap 1 (m7GpppNm)",
            CapType::CleanCap => "CleanCap (co-transcriptional cap 1)",
            CapType::Arca => "ARCA (anti-reverse cap analog)",
        }
    }

    /// `true` when the cap evades innate immune RNA sensors well — the
    /// 2′-O-methylated cap-1 chemistries.
    pub fn is_innate_immune_silent(self) -> bool {
        matches!(self, CapType::Cap1 | CapType::CleanCap)
    }
}

/// A therapeutic-mRNA construct: the five parts plus their chemistry.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct MrnaConstruct {
    /// The 5′ cap chemistry.
    pub cap: CapType,
    /// The 5′UTR sequence (RNA, `A C G U`).
    pub utr5: Vec<u8>,
    /// The coding sequence (RNA, `A C G U`; `AUG`…stop).
    pub cds: Vec<u8>,
    /// The 3′UTR sequence (RNA, `A C G U`).
    pub utr3: Vec<u8>,
    /// The poly-A tail length in nucleotides.
    pub poly_a_len: usize,
}

impl MrnaConstruct {
    /// The full transcript body (5′UTR + CDS + 3′UTR + poly-A), RNA.
    /// The cap is a chemistry, not a base sequence, so it is not part
    /// of this string.
    pub fn transcript(&self) -> Vec<u8> {
        let mut t = Vec::with_capacity(
            self.utr5.len() + self.cds.len() + self.utr3.len() + self.poly_a_len,
        );
        t.extend_from_slice(&self.utr5);
        t.extend_from_slice(&self.cds);
        t.extend_from_slice(&self.utr3);
        t.extend(std::iter::repeat_n(b'A', self.poly_a_len));
        t
    }

    /// Total transcript length in nucleotides (excluding the cap).
    pub fn len(&self) -> usize {
        self.utr5.len() + self.cds.len() + self.utr3.len() + self.poly_a_len
    }

    /// `true` when the transcript body is empty.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// The number of codons in the CDS (its length divided by three).
    pub fn codon_count(&self) -> usize {
        self.cds.len() / 3
    }

    /// The 0-based index where the CDS starts within
    /// [`transcript`](Self::transcript) — i.e. the 5′UTR length.
    pub fn cds_start(&self) -> usize {
        self.utr5.len()
    }
}

/// A builder that assembles and validates an [`MrnaConstruct`].
///
/// Accepts DNA or RNA for each part (DNA is transcribed); the CDS is
/// validated on [`build`](Self::build).
#[derive(Clone, Debug, Default)]
pub struct MrnaConstructBuilder {
    cap: Option<CapType>,
    utr5: Vec<u8>,
    cds: Vec<u8>,
    utr3: Vec<u8>,
    poly_a_len: usize,
}

impl MrnaConstructBuilder {
    /// A fresh builder (cap defaults to [`CapType::Cap1`], the
    /// therapeutic standard, if never set).
    pub fn new() -> Self {
        MrnaConstructBuilder::default()
    }

    /// Sets the 5′ cap chemistry.
    pub fn cap(mut self, cap: CapType) -> Self {
        self.cap = Some(cap);
        self
    }

    /// Sets the 5′UTR (DNA or RNA; transcribed to RNA).
    pub fn utr5(mut self, seq: impl AsRef<[u8]>) -> Self {
        self.utr5 = transcribe(seq.as_ref());
        self
    }

    /// Sets the CDS (DNA or RNA; transcribed to RNA).
    pub fn cds(mut self, seq: impl AsRef<[u8]>) -> Self {
        self.cds = transcribe(seq.as_ref());
        self
    }

    /// Sets the 3′UTR (DNA or RNA; transcribed to RNA).
    pub fn utr3(mut self, seq: impl AsRef<[u8]>) -> Self {
        self.utr3 = transcribe(seq.as_ref());
        self
    }

    /// Sets the poly-A tail length in nucleotides.
    pub fn poly_a(mut self, len: usize) -> Self {
        self.poly_a_len = len;
        self
    }

    /// Assembles and validates the construct.
    ///
    /// # Errors
    /// - [`GeneditingError::InvalidTarget`] if the CDS is empty, has a
    ///   length not divisible by three, does not begin with a start
    ///   codon, does not end with a stop codon, or any part contains a
    ///   non-ACGU base.
    /// - [`GeneditingError::Invalid`] for an absurd poly-A length
    ///   (`> 500` nt — beyond any physiological tail).
    pub fn build(self) -> Result<MrnaConstruct> {
        let cds = upper(&self.cds);
        if cds.is_empty() {
            return Err(GeneditingError::invalid_target("cds", "CDS is empty"));
        }
        if !is_acgu(&cds) {
            return Err(GeneditingError::invalid_target(
                "cds",
                "CDS must be A/C/G/U after transcription",
            ));
        }
        if cds.len() % 3 != 0 {
            return Err(GeneditingError::invalid_target(
                "cds",
                "CDS length must be a multiple of 3",
            ));
        }
        let code = GeneticCode::standard();
        // Validate the start codon — translate works on DNA codons, so
        // reverse-transcribe the first codon.
        let first = dna_codon(&cds[0..3]);
        if !code.is_start_codon(&first) {
            return Err(GeneditingError::invalid_target(
                "cds",
                "CDS does not begin with a start codon (AUG)",
            ));
        }
        let last = dna_codon(&cds[cds.len() - 3..]);
        if !code.is_stop_codon(&last) {
            return Err(GeneditingError::invalid_target(
                "cds",
                "CDS does not end with a stop codon",
            ));
        }
        for (name, part) in [("utr5", &self.utr5), ("utr3", &self.utr3)] {
            if !part.is_empty() && !is_acgu(part) {
                return Err(GeneditingError::invalid_target(
                    "region",
                    format!("{name} must be A/C/G/U after transcription"),
                ));
            }
        }
        if self.poly_a_len > 500 {
            return Err(GeneditingError::invalid(
                "poly_a_len",
                "poly-A tail length is beyond any physiological value (> 500 nt)",
            ));
        }
        Ok(MrnaConstruct {
            cap: self.cap.unwrap_or(CapType::Cap1),
            utr5: upper(&self.utr5),
            cds,
            utr3: upper(&self.utr3),
            poly_a_len: self.poly_a_len,
        })
    }
}

/// Reverse-transcribes a 3-base RNA codon to a DNA codon array so the
/// `valenx-bioseq` genetic-code tables (which key on DNA) can read it.
fn dna_codon(rna: &[u8]) -> [u8; 3] {
    let mut out = [b'N'; 3];
    for (i, &b) in rna.iter().take(3).enumerate() {
        out[i] = match b.to_ascii_uppercase() {
            b'U' => b'T',
            other => other,
        };
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn builds_a_valid_construct() {
        // ATG + one sense codon + stop.
        let c = MrnaConstructBuilder::new()
            .cap(CapType::Cap1)
            .utr5(b"GGGACC")
            .cds(b"ATGGCCTAA")
            .utr3(b"AAUAAA")
            .poly_a(120)
            .build()
            .unwrap();
        assert_eq!(c.cap, CapType::Cap1);
        assert_eq!(c.codon_count(), 3);
        assert_eq!(c.poly_a_len, 120);
        // Stored as RNA — the DNA CDS was transcribed.
        assert!(!c.cds.contains(&b'T'));
    }

    #[test]
    fn transcript_concatenates_parts() {
        let c = MrnaConstructBuilder::new()
            .utr5(b"GG")
            .cds(b"ATGTAA")
            .utr3(b"CC")
            .poly_a(5)
            .build()
            .unwrap();
        let t = c.transcript();
        assert_eq!(t.len(), 2 + 6 + 2 + 5);
        assert!(t.ends_with(b"AAAAA"));
        assert_eq!(c.cds_start(), 2);
    }

    #[test]
    fn rejects_cds_without_start_codon() {
        let err = MrnaConstructBuilder::new()
            .cds(b"GCCGCCTAA")
            .build()
            .unwrap_err();
        assert_eq!(err.code(), "genediting.invalid_target");
    }

    #[test]
    fn rejects_cds_without_stop_codon() {
        let err = MrnaConstructBuilder::new()
            .cds(b"ATGGCCGCC")
            .build()
            .unwrap_err();
        assert_eq!(err.code(), "genediting.invalid_target");
    }

    #[test]
    fn rejects_cds_not_multiple_of_three() {
        let err = MrnaConstructBuilder::new()
            .cds(b"ATGGCCTA")
            .build()
            .unwrap_err();
        assert_eq!(err.code(), "genediting.invalid_target");
    }

    #[test]
    fn rejects_empty_cds() {
        assert!(MrnaConstructBuilder::new().build().is_err());
    }

    #[test]
    fn rejects_absurd_poly_a() {
        let err = MrnaConstructBuilder::new()
            .cds(b"ATGTAA")
            .poly_a(9999)
            .build()
            .unwrap_err();
        assert_eq!(err.code(), "genediting.invalid");
    }

    #[test]
    fn cap_default_is_cap1() {
        let c = MrnaConstructBuilder::new().cds(b"ATGTAA").build().unwrap();
        assert_eq!(c.cap, CapType::Cap1);
    }

    #[test]
    fn cap1_is_innate_immune_silent() {
        assert!(CapType::Cap1.is_innate_immune_silent());
        assert!(CapType::CleanCap.is_innate_immune_silent());
        assert!(!CapType::Cap0.is_innate_immune_silent());
    }

    #[test]
    fn accepts_rna_input_directly() {
        let c = MrnaConstructBuilder::new()
            .cds(b"AUGGCCUAA")
            .build()
            .unwrap();
        assert_eq!(c.codon_count(), 3);
    }
}
