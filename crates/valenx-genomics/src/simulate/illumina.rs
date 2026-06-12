//! Illumina-style short-read simulator (ART-class).
//!
//! ART simulates Illumina reads by drawing fragments from a reference,
//! then corrupting each base with a **position-specific** substitution
//! error model and assigning a per-base quality from a position
//! profile. This module implements that model:
//!
//! - a per-cycle error-rate profile (errors rise toward the 3′ end —
//!   the universal Illumina pattern);
//! - a per-cycle quality profile, with the simulated quality string
//!   *consistent* with the error decision (a corrupted base gets a low
//!   quality, a correct one a high quality);
//! - substitution-only errors (Illumina's dominant mode), the error
//!   base drawn uniformly from the three non-true bases.
//!
//! Everything is seeded through [`crate::util::rng::Rng`] for
//! reproducibility.
//!
//! ## v1 scope
//!
//! Substitution errors only — Illumina indels are rare and ART itself
//! defaults them off. The error / quality profiles are smooth analytic
//! curves (a documented simplification of ART's empirically-trained
//! per-instrument profiles); they capture the right *shape*. GC bias
//! in fragment selection is not modelled.

use crate::error::{GenomicsError, Result};
use crate::util::rng::Rng;
use valenx_bioseq::alphabet::SeqKind;
use valenx_bioseq::io::fastq::FastqRecord;
use valenx_bioseq::record::SeqRecord;
use valenx_bioseq::seq::Seq;

/// Configuration for the Illumina simulator.
#[derive(Clone, Debug, PartialEq)]
pub struct IlluminaProfile {
    /// Read length in bases.
    pub read_length: usize,
    /// Mean per-base substitution-error rate at the 5′ end.
    pub error_5prime: f64,
    /// Mean per-base substitution-error rate at the 3′ end (typically
    /// higher — the Illumina quality droop).
    pub error_3prime: f64,
    /// Phred quality assigned to a *correct* base near the 5′ end.
    pub quality_high: u8,
    /// Phred quality assigned to a *correct* base near the 3′ end.
    pub quality_low: u8,
}

impl IlluminaProfile {
    /// A HiSeq-2500-like 150 bp profile.
    pub fn hiseq_150() -> Self {
        IlluminaProfile {
            read_length: 150,
            error_5prime: 0.001,
            error_3prime: 0.01,
            quality_high: 40,
            quality_low: 28,
        }
    }

    /// A MiSeq-like 250 bp profile.
    pub fn miseq_250() -> Self {
        IlluminaProfile {
            read_length: 250,
            error_5prime: 0.0015,
            error_3prime: 0.02,
            quality_high: 38,
            quality_low: 25,
        }
    }

    /// Per-cycle error rate at 0-based cycle `i` — a linear ramp from
    /// [`error_5prime`](Self::error_5prime) to
    /// [`error_3prime`](Self::error_3prime).
    pub fn error_rate(&self, i: usize) -> f64 {
        if self.read_length <= 1 {
            return self.error_5prime;
        }
        let t = i as f64 / (self.read_length - 1) as f64;
        self.error_5prime + t * (self.error_3prime - self.error_5prime)
    }

    /// Phred quality for a *correct* base at 0-based cycle `i` — a
    /// linear droop from [`quality_high`](Self::quality_high) to
    /// [`quality_low`](Self::quality_low).
    pub fn correct_quality(&self, i: usize) -> u8 {
        if self.read_length <= 1 {
            return self.quality_high;
        }
        let t = i as f64 / (self.read_length - 1) as f64;
        let q = self.quality_high as f64 - t * (self.quality_high as f64 - self.quality_low as f64);
        q.round().clamp(2.0, 60.0) as u8
    }

    fn validate(&self) -> Result<()> {
        if self.read_length == 0 {
            return Err(GenomicsError::invalid("read_length", "must be positive"));
        }
        for (name, v) in [
            ("error_5prime", self.error_5prime),
            ("error_3prime", self.error_3prime),
        ] {
            if !(0.0..=1.0).contains(&v) {
                return Err(GenomicsError::invalid(name, "must be in [0, 1]"));
            }
        }
        Ok(())
    }
}

/// The complement of a DNA base (uppercased; non-ACGT → `N`).
fn complement(b: u8) -> u8 {
    match b.to_ascii_uppercase() {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        _ => b'N',
    }
}

/// Reverse-complements a base slice.
fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| complement(b)).collect()
}

/// Picks an error base uniformly from the three bases other than
/// `true_base`.
fn error_base(true_base: u8, rng: &mut Rng) -> u8 {
    const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];
    let tb = true_base.to_ascii_uppercase();
    let choices: Vec<u8> = BASES.iter().copied().filter(|&b| b != tb).collect();
    if choices.is_empty() {
        return tb;
    }
    choices[rng.below(choices.len())]
}

/// Applies the position-specific error model to one true read,
/// returning the corrupted bases and a consistent quality vector.
///
/// `true_read` must already be in read orientation. For each cycle,
/// the base is substituted with probability `profile.error_rate(i)`; a
/// substituted base is given a quality reflecting that error rate, a
/// correct base gets `profile.correct_quality(i)`.
pub fn apply_error_model(
    true_read: &[u8],
    profile: &IlluminaProfile,
    rng: &mut Rng,
) -> (Vec<u8>, Vec<u8>) {
    let mut bases = Vec::with_capacity(true_read.len());
    let mut quals = Vec::with_capacity(true_read.len());
    for (i, &tb) in true_read.iter().enumerate() {
        let e = profile.error_rate(i);
        if rng.chance(e) {
            bases.push(error_base(tb, rng));
            // Quality of an erroneous base: derived from the local
            // error rate, floored low.
            let q = (-10.0 * e.max(1e-6).log10()).round().clamp(2.0, 20.0) as u8;
            quals.push(q);
        } else {
            bases.push(tb.to_ascii_uppercase());
            quals.push(profile.correct_quality(i));
        }
    }
    (bases, quals)
}

/// Simulates `n_reads` single-end Illumina reads from a reference
/// sequence.
///
/// Each read starts at a uniformly-random position that admits a full
/// read length, takes a random strand, and is corrupted by
/// [`apply_error_model`]. The reference must be at least
/// `profile.read_length` bases long.
pub fn simulate_reads(
    reference: &[u8],
    profile: &IlluminaProfile,
    n_reads: usize,
    seed: u64,
) -> Result<Vec<FastqRecord>> {
    profile.validate()?;
    if reference.len() < profile.read_length {
        return Err(GenomicsError::invalid(
            "reference",
            format!(
                "reference length {} < read length {}",
                reference.len(),
                profile.read_length
            ),
        ));
    }
    let upper: Vec<u8> = reference.iter().map(|b| b.to_ascii_uppercase()).collect();
    let max_start = upper.len() - profile.read_length;
    let mut rng = Rng::new(seed);
    let mut out = Vec::with_capacity(n_reads);

    for i in 0..n_reads {
        let start = rng.below(max_start + 1);
        let frag = &upper[start..start + profile.read_length];
        let reverse = rng.chance(0.5);
        let true_read: Vec<u8> = if reverse {
            revcomp(frag)
        } else {
            frag.to_vec()
        };
        let (bases, quals) = apply_error_model(&true_read, profile, &mut rng);
        let seq = Seq::new(SeqKind::Dna, &bases)
            .map_err(|e| GenomicsError::invalid("simulated_read", e.to_string()))?;
        let id = format!("sim_read_{i}");
        let desc = format!(
            "pos={} strand={}",
            start + 1,
            if reverse { '-' } else { '+' }
        );
        let mut rec = SeqRecord::new(id, seq);
        rec.description = desc;
        out.push(FastqRecord {
            record: rec,
            quality: quals,
        });
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ref_seq(len: usize) -> Vec<u8> {
        // A pseudo-random but deterministic reference.
        let bases = [b'A', b'C', b'G', b'T'];
        (0..len).map(|i| bases[(i * 7 + 3) % 4]).collect()
    }

    #[test]
    fn error_rate_ramps_up() {
        let p = IlluminaProfile::hiseq_150();
        assert!(p.error_rate(0) < p.error_rate(149));
        assert!((p.error_rate(0) - p.error_5prime).abs() < 1e-9);
        assert!((p.error_rate(149) - p.error_3prime).abs() < 1e-9);
    }

    #[test]
    fn quality_droops() {
        let p = IlluminaProfile::hiseq_150();
        assert!(p.correct_quality(0) >= p.correct_quality(149));
    }

    #[test]
    fn simulate_yields_requested_count() {
        let reference = ref_seq(5000);
        let reads = simulate_reads(&reference, &IlluminaProfile::hiseq_150(), 100, 42).unwrap();
        assert_eq!(reads.len(), 100);
        for r in &reads {
            assert_eq!(r.len(), 150);
            assert_eq!(r.quality.len(), 150);
        }
    }

    #[test]
    fn simulation_is_deterministic() {
        let reference = ref_seq(5000);
        let p = IlluminaProfile::hiseq_150();
        let a = simulate_reads(&reference, &p, 50, 7).unwrap();
        let b = simulate_reads(&reference, &p, 50, 7).unwrap();
        assert_eq!(a, b);
    }

    #[test]
    fn rejects_short_reference() {
        let reference = ref_seq(50);
        assert!(simulate_reads(&reference, &IlluminaProfile::hiseq_150(), 10, 1).is_err());
    }

    #[test]
    fn error_model_keeps_length() {
        let mut rng = Rng::new(1);
        let true_read = ref_seq(150);
        let (bases, quals) = apply_error_model(&true_read, &IlluminaProfile::hiseq_150(), &mut rng);
        assert_eq!(bases.len(), 150);
        assert_eq!(quals.len(), 150);
    }

    #[test]
    fn high_error_rate_corrupts_some_bases() {
        // A degenerate profile with 100% error: every base substituted.
        let p = IlluminaProfile {
            read_length: 100,
            error_5prime: 1.0,
            error_3prime: 1.0,
            quality_high: 40,
            quality_low: 30,
        };
        let true_read = vec![b'A'; 100];
        let mut rng = Rng::new(3);
        let (bases, _) = apply_error_model(&true_read, &p, &mut rng);
        // No base may remain 'A' under a forced substitution.
        assert!(bases.iter().all(|&b| b != b'A'));
    }

    #[test]
    fn zero_error_keeps_reference_bases() {
        let p = IlluminaProfile {
            read_length: 100,
            error_5prime: 0.0,
            error_3prime: 0.0,
            quality_high: 40,
            quality_low: 30,
        };
        let true_read = ref_seq(100);
        let mut rng = Rng::new(3);
        let (bases, _) = apply_error_model(&true_read, &p, &mut rng);
        assert_eq!(bases, true_read);
    }

    #[test]
    fn rejects_bad_profile() {
        let mut p = IlluminaProfile::hiseq_150();
        p.read_length = 0;
        assert!(simulate_reads(&ref_seq(1000), &p, 1, 1).is_err());
        let mut p = IlluminaProfile::hiseq_150();
        p.error_3prime = 2.0;
        assert!(simulate_reads(&ref_seq(1000), &p, 1, 1).is_err());
    }
}
