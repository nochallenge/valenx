//! Long-read simulator (PacBio / Nanopore-style — pbsim / Badread-class).
//!
//! Third-generation reads differ from Illumina in two ways the
//! simulator must capture: reads are **long** with a broad length
//! distribution, and the error mode is **indel-heavy** — insertions
//! and deletions dominate over substitutions, the opposite of
//! Illumina.
//!
//! This module simulates that:
//!
//! - read lengths drawn from a (clamped) normal length distribution;
//! - a per-base error process with three outcomes — substitution,
//!   insertion, deletion — whose relative weights match the chosen
//!   chemistry ([`LongReadTech::PacBioClr`],
//!   [`LongReadTech::PacBioHifi`], [`LongReadTech::Nanopore`]);
//! - a per-base quality that reflects the chemistry's accuracy.
//!
//! All randomness is seeded through [`crate::util::rng::Rng`].
//!
//! ## v1 scope
//!
//! The error process is per-base i.i.d. — it does not model Nanopore's
//! homopolymer-length systematic bias or PacBio's strand-correlated
//! errors. Read lengths are normal, not the empirical
//! log-normal/gamma mixtures of pbsim3 / Badread. It is a real
//! indel-heavy long-read generator, not those tools' trained models.

use crate::error::{GenomicsError, Result};
use crate::util::rng::Rng;
use valenx_bioseq::alphabet::SeqKind;
use valenx_bioseq::io::fastq::FastqRecord;
use valenx_bioseq::record::SeqRecord;
use valenx_bioseq::seq::Seq;

/// Long-read sequencing chemistry.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum LongReadTech {
    /// PacBio continuous-long-read (CLR) — ~10-15 % error, indel-heavy.
    PacBioClr,
    /// PacBio HiFi / CCS — ~0.5-1 % error after consensus.
    PacBioHifi,
    /// Oxford Nanopore — ~5-10 % error, indel-heavy with a deletion
    /// skew.
    Nanopore,
}

impl LongReadTech {
    /// `(substitution, insertion, deletion)` per-base error
    /// probabilities for this chemistry.
    pub fn error_rates(self) -> (f64, f64, f64) {
        match self {
            // CLR: ~13 % total, insertions slightly dominant.
            LongReadTech::PacBioClr => (0.01, 0.07, 0.05),
            // HiFi: ~0.6 % total.
            LongReadTech::PacBioHifi => (0.003, 0.0015, 0.0015),
            // Nanopore: ~8 % total, deletion-skewed.
            LongReadTech::Nanopore => (0.02, 0.025, 0.035),
        }
    }

    /// Mean Phred quality for a base from this chemistry.
    pub fn base_quality(self) -> u8 {
        match self {
            LongReadTech::PacBioClr => 9,
            LongReadTech::PacBioHifi => 30,
            LongReadTech::Nanopore => 12,
        }
    }
}

/// Configuration for the long-read simulator.
#[derive(Clone, Debug, PartialEq)]
pub struct LongReadProfile {
    /// The sequencing chemistry.
    pub tech: LongReadTech,
    /// Mean read length.
    pub mean_length: usize,
    /// Standard deviation of the read length.
    pub length_sd: f64,
    /// Reads shorter than this are re-drawn (a hard floor).
    pub min_length: usize,
}

impl LongReadProfile {
    /// A PacBio-CLR-like profile: ~10 kb reads.
    pub fn pacbio_clr() -> Self {
        LongReadProfile {
            tech: LongReadTech::PacBioClr,
            mean_length: 10_000,
            length_sd: 4_000.0,
            min_length: 500,
        }
    }

    /// A PacBio-HiFi-like profile: ~15 kb reads.
    pub fn pacbio_hifi() -> Self {
        LongReadProfile {
            tech: LongReadTech::PacBioHifi,
            mean_length: 15_000,
            length_sd: 3_000.0,
            min_length: 1_000,
        }
    }

    /// A Nanopore-like profile: long, very broad length spread.
    pub fn nanopore() -> Self {
        LongReadProfile {
            tech: LongReadTech::Nanopore,
            mean_length: 8_000,
            length_sd: 6_000.0,
            min_length: 200,
        }
    }

    /// Draws a read length from the (clamped) normal length
    /// distribution.
    pub fn draw_length(&self, rng: &mut Rng) -> usize {
        let mut len = rng
            .next_normal(self.mean_length as f64, self.length_sd)
            .round() as i64;
        if len < self.min_length as i64 {
            len = self.min_length as i64;
        }
        len.max(1) as usize
    }
}

const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];

fn random_base(rng: &mut Rng) -> u8 {
    BASES[rng.below(4)]
}

fn other_base(b: u8, rng: &mut Rng) -> u8 {
    let bu = b.to_ascii_uppercase();
    let choices: Vec<u8> = BASES.iter().copied().filter(|&x| x != bu).collect();
    choices[rng.below(choices.len())]
}

fn complement(b: u8) -> u8 {
    match b.to_ascii_uppercase() {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        _ => b'N',
    }
}

/// Applies the indel-heavy long-read error process to a true fragment.
///
/// Walks the fragment base by base. At each base, with the chemistry's
/// `(sub, ins, del)` probabilities one of: a substitution (emit a
/// wrong base), an insertion (emit a random base, do **not** advance
/// the reference), a deletion (advance the reference, emit nothing), or
/// — most often — a clean copy. Insertions and substitutions get a low
/// quality; clean bases get the chemistry quality.
pub fn apply_long_read_errors(
    fragment: &[u8],
    tech: LongReadTech,
    rng: &mut Rng,
) -> (Vec<u8>, Vec<u8>) {
    let (p_sub, p_ins, p_del) = tech.error_rates();
    let clean_q = tech.base_quality();
    let err_q = (clean_q / 3).max(2);

    let mut bases = Vec::with_capacity(fragment.len());
    let mut quals = Vec::with_capacity(fragment.len());

    let mut i = 0usize;
    while i < fragment.len() {
        let roll = rng.next_f64();
        if roll < p_ins {
            // Insertion: emit a random base, stay on the same ref base.
            bases.push(random_base(rng));
            quals.push(err_q);
            // No `i` advance — the same reference base is processed
            // again next iteration.
        } else if roll < p_ins + p_del {
            // Deletion: skip this reference base, emit nothing.
            i += 1;
        } else if roll < p_ins + p_del + p_sub {
            // Substitution.
            bases.push(other_base(fragment[i], rng));
            quals.push(err_q);
            i += 1;
        } else {
            // Clean copy.
            bases.push(fragment[i].to_ascii_uppercase());
            quals.push(clean_q);
            i += 1;
        }
    }
    (bases, quals)
}

/// Simulates `n_reads` long reads from a reference sequence.
///
/// Each read draws a length from the profile's length distribution,
/// picks a random start admitting that length (a length longer than
/// the reference is capped to the reference length), picks a strand,
/// and is corrupted by [`apply_long_read_errors`].
pub fn simulate_long_reads(
    reference: &[u8],
    profile: &LongReadProfile,
    n_reads: usize,
    seed: u64,
) -> Result<Vec<FastqRecord>> {
    if reference.is_empty() {
        return Err(GenomicsError::invalid("reference", "reference is empty"));
    }
    if profile.mean_length == 0 {
        return Err(GenomicsError::invalid("mean_length", "must be positive"));
    }
    let upper: Vec<u8> = reference.iter().map(|b| b.to_ascii_uppercase()).collect();
    let mut rng = Rng::new(seed);
    let mut out = Vec::with_capacity(n_reads);

    for i in 0..n_reads {
        let want = profile.draw_length(&mut rng).min(upper.len());
        let max_start = upper.len() - want;
        let start = rng.below(max_start + 1);
        let frag = &upper[start..start + want];
        let reverse = rng.chance(0.5);
        let true_frag: Vec<u8> = if reverse {
            frag.iter().rev().map(|&b| complement(b)).collect()
        } else {
            frag.to_vec()
        };
        let (bases, quals) = apply_long_read_errors(&true_frag, profile.tech, &mut rng);
        // An all-deletion fragment could empty the read; guard it.
        if bases.is_empty() {
            continue;
        }
        let seq = Seq::new(SeqKind::Dna, &bases)
            .map_err(|e| GenomicsError::invalid("simulated_read", e.to_string()))?;
        let id = format!("sim_lread_{i}");
        let mut rec = SeqRecord::new(id, seq);
        rec.description = format!(
            "pos={} len={} strand={}",
            start + 1,
            want,
            if reverse { '-' } else { '+' }
        );
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
        let bases = [b'A', b'C', b'G', b'T'];
        (0..len).map(|i| bases[(i * 5 + 1) % 4]).collect()
    }

    #[test]
    fn error_rates_distinguish_chemistries() {
        let hifi: f64 = {
            let (s, i, d) = LongReadTech::PacBioHifi.error_rates();
            s + i + d
        };
        let clr: f64 = {
            let (s, i, d) = LongReadTech::PacBioClr.error_rates();
            s + i + d
        };
        assert!(hifi < clr, "HiFi must be more accurate than CLR");
    }

    #[test]
    fn simulate_yields_reads() {
        let reference = ref_seq(100_000);
        let reads =
            simulate_long_reads(&reference, &LongReadProfile::pacbio_clr(), 20, 42).unwrap();
        // Some reads may be dropped if they end up empty, but CLR at
        // 10kb should always survive.
        assert_eq!(reads.len(), 20);
        for r in &reads {
            assert_eq!(r.len(), r.quality.len());
            assert!(r.len() > 100);
        }
    }

    #[test]
    fn simulation_is_deterministic() {
        let reference = ref_seq(100_000);
        let p = LongReadProfile::nanopore();
        let a = simulate_long_reads(&reference, &p, 10, 7).unwrap();
        let b = simulate_long_reads(&reference, &p, 10, 7).unwrap();
        assert_eq!(a, b);
    }

    #[test]
    fn indel_errors_change_read_length() {
        // CLR error model is indel-heavy; a long fragment's read length
        // should differ from the fragment length for at least one read.
        let reference = ref_seq(50_000);
        let reads =
            simulate_long_reads(&reference, &LongReadProfile::pacbio_clr(), 30, 11).unwrap();
        let any_length_changed = reads.iter().any(|r| {
            // Parse "len=" from the description.
            r.record
                .description
                .split_whitespace()
                .find_map(|t| t.strip_prefix("len="))
                .and_then(|s| s.parse::<usize>().ok())
                .map(|frag_len| frag_len != r.len())
                .unwrap_or(false)
        });
        assert!(any_length_changed, "indel model should perturb length");
    }

    #[test]
    fn hifi_is_near_lossless() {
        // HiFi error rate is tiny; a short fragment usually copies
        // through unchanged.
        let frag = ref_seq(200);
        let mut rng = Rng::new(1);
        let (bases, _) = apply_long_read_errors(&frag, LongReadTech::PacBioHifi, &mut rng);
        // Length should be within a few bases of the input.
        let diff = (bases.len() as i64 - frag.len() as i64).abs();
        assert!(diff < 20, "HiFi perturbed length by {diff}");
    }

    #[test]
    fn length_distribution_respects_minimum() {
        let p = LongReadProfile {
            tech: LongReadTech::Nanopore,
            mean_length: 1000,
            length_sd: 5000.0, // huge spread can push below the floor
            min_length: 300,
        };
        let mut rng = Rng::new(5);
        for _ in 0..200 {
            assert!(p.draw_length(&mut rng) >= 300);
        }
    }

    #[test]
    fn rejects_empty_reference() {
        assert!(simulate_long_reads(&[], &LongReadProfile::pacbio_clr(), 1, 1).is_err());
    }
}
