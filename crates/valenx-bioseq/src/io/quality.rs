//! Phred / Solexa quality-score codecs.
//!
//! FASTQ stores per-base quality as printable ASCII. Three historical
//! encodings exist, differing in the ASCII offset and (for Solexa) the
//! score definition:
//!
//! | Variant | Offset | Score | Range typically |
//! |---|---|---|---|
//! | [`QualityEncoding::Phred33`] | 33 | Phred | Sanger / Illumina ≥1.8 |
//! | [`QualityEncoding::Phred64`] | 64 | Phred | Illumina 1.3–1.7 |
//! | [`QualityEncoding::Solexa64`] | 64 | Solexa (log-odds) | Solexa / Illumina <1.3 |
//!
//! Phred quality `Q` relates to base-call error probability `p` by
//! `Q = -10·log10(p)`. Solexa quality uses the odds ratio:
//! `Q_sol = -10·log10(p / (1-p))`.

use crate::error::{BioseqError, Result};

/// A FASTQ quality-string encoding.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub enum QualityEncoding {
    /// Phred+33 — Sanger and Illumina 1.8+. The modern default.
    Phred33,
    /// Phred+64 — Illumina 1.3 through 1.7.
    Phred64,
    /// Solexa+64 — Solexa / early Illumina (pre-1.3). Log-odds score.
    Solexa64,
}

impl QualityEncoding {
    /// The ASCII offset subtracted on decode / added on encode.
    pub fn offset(self) -> u8 {
        match self {
            QualityEncoding::Phred33 => 33,
            QualityEncoding::Phred64 | QualityEncoding::Solexa64 => 64,
        }
    }

    /// `true` for the Solexa log-odds score (vs. the Phred score).
    pub fn is_solexa(self) -> bool {
        matches!(self, QualityEncoding::Solexa64)
    }
}

/// Converts a Phred quality score to a base-call error probability
/// `p = 10^(-Q/10)`.
pub fn phred_to_error_prob(q: f64) -> f64 {
    10f64.powf(-q / 10.0)
}

/// Converts a base-call error probability to a Phred quality score
/// `Q = -10·log10(p)`. `p` is clamped away from 0 so the result is
/// finite.
pub fn error_prob_to_phred(p: f64) -> f64 {
    let p = p.clamp(1e-300, 1.0);
    -10.0 * p.log10()
}

/// Converts a Solexa quality score to a base-call error probability.
///
/// The Solexa score is defined as `Q = -10·log10(p / (1-p))`, so
/// `p / (1-p) = 10^(-Q/10)` and therefore
/// `p = 10^(-Q/10) / (1 + 10^(-Q/10))`. (The exponent is `-Q/10`, not
/// `+Q/10`: a *higher* Solexa score must map to a *lower* `p`.)
pub fn solexa_to_error_prob(q: f64) -> f64 {
    let t = 10f64.powf(-q / 10.0);
    t / (1.0 + t)
}

/// Converts an error probability to a Solexa quality score.
/// `Q = -10·log10(p / (1-p))`.
pub fn error_prob_to_solexa(p: f64) -> f64 {
    let p = p.clamp(1e-300, 1.0 - 1e-12);
    -10.0 * (p / (1.0 - p)).log10()
}

/// Converts a Phred score to the equivalent Solexa score (same error
/// probability).
pub fn phred_to_solexa(q: f64) -> f64 {
    error_prob_to_solexa(phred_to_error_prob(q))
}

/// Converts a Solexa score to the equivalent Phred score.
pub fn solexa_to_phred(q: f64) -> f64 {
    error_prob_to_phred(solexa_to_error_prob(q))
}

/// Decodes a FASTQ quality string into per-base **Phred** scores.
///
/// Solexa-encoded input is converted to the Phred scale so callers
/// always work in one unit. Returns [`BioseqError::Parse`] on a
/// character below the encoding's offset.
pub fn decode(qual: &[u8], encoding: QualityEncoding) -> Result<Vec<u8>> {
    let off = encoding.offset();
    let mut out = Vec::with_capacity(qual.len());
    for &c in qual {
        if c < off {
            return Err(BioseqError::parse(
                "fastq",
                format!(
                    "quality char {:?} (0x{c:02x}) below offset {off} for {encoding:?}",
                    c as char
                ),
            ));
        }
        let raw = (c - off) as i32;
        let phred = if encoding.is_solexa() {
            // Solexa scores can be negative; round to the Phred scale.
            solexa_to_phred(raw as f64).round().max(0.0) as i32
        } else {
            raw
        };
        out.push(phred.clamp(0, 93) as u8);
    }
    Ok(out)
}

/// Encodes per-base **Phred** scores into a FASTQ quality string.
///
/// For a Solexa encoding the Phred scores are first converted to the
/// Solexa scale. Scores are clamped to the printable ASCII range of
/// the chosen encoding.
pub fn encode(phred: &[u8], encoding: QualityEncoding) -> Vec<u8> {
    let off = encoding.offset() as i32;
    phred
        .iter()
        .map(|&q| {
            let score = if encoding.is_solexa() {
                phred_to_solexa(q as f64).round() as i32
            } else {
                q as i32
            };
            // Printable ASCII is 33..=126.
            let c = (score + off).clamp(33, 126);
            c as u8
        })
        .collect()
}

/// Mean Phred quality of a decoded score vector. Returns `0.0` for an
/// empty slice.
pub fn mean_quality(phred: &[u8]) -> f64 {
    if phred.is_empty() {
        return 0.0;
    }
    let sum: u32 = phred.iter().map(|&q| q as u32).sum();
    sum as f64 / phred.len() as f64
}

/// Heuristically guesses the encoding of a quality string by its
/// character range. Returns `None` when the range is ambiguous (a
/// common situation — many strings are valid under both Phred+33 and
/// Phred+64).
pub fn guess_encoding(qual: &[u8]) -> Option<QualityEncoding> {
    let min = *qual.iter().min()?;
    let max = *qual.iter().max()?;
    if min < 59 {
        // Below ';' — only Phred+33 reaches here.
        Some(QualityEncoding::Phred33)
    } else if max > 104 {
        // Above 'h' — only Phred+64-family reaches here. Solexa vs.
        // Phred64 is indistinguishable from range alone; assume Phred64.
        Some(QualityEncoding::Phred64)
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn phred_error_prob_roundtrip() {
        // Q20 -> p = 0.01.
        assert!((phred_to_error_prob(20.0) - 0.01).abs() < 1e-12);
        assert!((error_prob_to_phred(0.01) - 20.0).abs() < 1e-9);
        // Q30 -> p = 0.001.
        assert!((phred_to_error_prob(30.0) - 0.001).abs() < 1e-12);
    }

    #[test]
    fn phred33_decode() {
        // '!' = 33 -> Q0 ; 'I' = 73 -> Q40.
        let q = decode(b"!I", QualityEncoding::Phred33).unwrap();
        assert_eq!(q, vec![0, 40]);
    }

    #[test]
    fn phred64_decode() {
        // '@' = 64 -> Q0 ; 'h' = 104 -> Q40.
        let q = decode(b"@h", QualityEncoding::Phred64).unwrap();
        assert_eq!(q, vec![0, 40]);
    }

    #[test]
    fn encode_decode_roundtrip_phred33() {
        let scores = vec![0u8, 10, 20, 30, 40];
        let enc = encode(&scores, QualityEncoding::Phred33);
        let dec = decode(&enc, QualityEncoding::Phred33).unwrap();
        assert_eq!(dec, scores);
    }

    #[test]
    fn below_offset_is_error() {
        // ' ' = 32, below the Phred+33 offset.
        assert!(decode(b" ", QualityEncoding::Phred33).is_err());
    }

    #[test]
    fn solexa_conversion_high_quality_matches_phred() {
        // At high quality Solexa and Phred scores converge.
        let p = phred_to_solexa(40.0);
        assert!((p - 40.0).abs() < 0.001, "got {p}");
        // At low quality they diverge: Solexa(Q10) < 10.
        assert!(phred_to_solexa(10.0) < 10.0);
    }

    #[test]
    fn solexa_roundtrip_through_error_prob() {
        let q = 25.0;
        let back = solexa_to_phred(phred_to_solexa(q));
        assert!((back - q).abs() < 1e-6, "got {back}");
    }

    #[test]
    fn mean_quality_computation() {
        assert_eq!(mean_quality(&[10, 20, 30]), 20.0);
        assert_eq!(mean_quality(&[]), 0.0);
    }

    #[test]
    fn guess_encoding_extremes() {
        assert_eq!(guess_encoding(b"!!!!"), Some(QualityEncoding::Phred33));
        assert_eq!(
            guess_encoding(b"iiii"),
            Some(QualityEncoding::Phred64)
        );
        // Ambiguous mid-range -> None.
        assert_eq!(guess_encoding(b"BBBB"), None);
    }
}
