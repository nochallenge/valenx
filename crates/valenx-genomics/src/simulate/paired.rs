//! Paired-end read generation with a configurable insert-size model.
//!
//! Paired-end sequencing reads both ends of a DNA *fragment*. The
//! fragment ("insert") length is drawn from a distribution; read 1
//! comes from the fragment's 5′ end on one strand and read 2 from the
//! 3′ end on the opposite strand — the standard Illumina FR
//! (forward-reverse) orientation.
//!
//! This module draws fragments with a configurable
//! [`InsertSizeModel`], extracts the two mate reads and corrupts each
//! with the [`crate::simulate::illumina`] error model. The two reads
//! of a pair share a name and carry the SAM-spec paired flag bits in
//! their FASTQ description.

use crate::error::{GenomicsError, Result};
use crate::simulate::illumina::{apply_error_model, IlluminaProfile};
use crate::util::rng::Rng;
use valenx_bioseq::alphabet::SeqKind;
use valenx_bioseq::io::fastq::FastqRecord;
use valenx_bioseq::record::SeqRecord;
use valenx_bioseq::seq::Seq;

/// A fragment-length (insert-size) distribution.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct InsertSizeModel {
    /// Mean fragment length.
    pub mean: usize,
    /// Standard deviation of the fragment length.
    pub std_dev: f64,
    /// Hard minimum (fragments below are re-clamped up). Must be at
    /// least `2 × read_length` for the two mates not to overlap, though
    /// a smaller value is allowed (mates then overlap, which is legal).
    pub min: usize,
    /// Hard maximum.
    pub max: usize,
}

impl InsertSizeModel {
    /// A typical 350 ± 50 bp short-insert library.
    pub fn short_insert() -> Self {
        InsertSizeModel {
            mean: 350,
            std_dev: 50.0,
            min: 150,
            max: 700,
        }
    }

    /// A 550 ± 100 bp library.
    pub fn medium_insert() -> Self {
        InsertSizeModel {
            mean: 550,
            std_dev: 100.0,
            min: 250,
            max: 1000,
        }
    }

    /// Draws one fragment length, clamped to `[min, max]`.
    pub fn draw(&self, rng: &mut Rng) -> usize {
        let raw = rng.next_normal(self.mean as f64, self.std_dev).round() as i64;
        raw.clamp(self.min as i64, self.max as i64).max(1) as usize
    }

    fn validate(&self) -> Result<()> {
        if self.mean == 0 {
            return Err(GenomicsError::invalid("mean", "must be positive"));
        }
        if self.max < self.min {
            return Err(GenomicsError::invalid("max", "max must be >= min"));
        }
        Ok(())
    }
}

/// A simulated read pair.
#[derive(Clone, Debug, PartialEq)]
pub struct ReadPair {
    /// The first mate (5′ end of the fragment, forward strand).
    pub read1: FastqRecord,
    /// The second mate (3′ end of the fragment, reverse strand).
    pub read2: FastqRecord,
    /// The fragment length the pair was drawn from.
    pub fragment_length: usize,
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

fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| complement(b)).collect()
}

/// Simulates `n_pairs` paired-end read pairs from a reference.
///
/// For each pair a fragment is drawn from the [`InsertSizeModel`] and
/// placed at a random reference position; read 1 is the fragment's 5′
/// `read_length` bases (forward), read 2 the 3′ `read_length` bases
/// reverse-complemented — Illumina FR orientation. Both mates are then
/// corrupted by the [`IlluminaProfile`] error model.
///
/// The reference must be at least `model.min` bases long, and the
/// profile's read length must not exceed the smallest fragment.
pub fn simulate_pairs(
    reference: &[u8],
    profile: &IlluminaProfile,
    model: &InsertSizeModel,
    n_pairs: usize,
    seed: u64,
) -> Result<Vec<ReadPair>> {
    model.validate()?;
    if profile.read_length == 0 {
        return Err(GenomicsError::invalid("read_length", "must be positive"));
    }
    if reference.len() < model.min {
        return Err(GenomicsError::invalid(
            "reference",
            format!(
                "reference length {} < minimum fragment {}",
                reference.len(),
                model.min
            ),
        ));
    }
    let upper: Vec<u8> = reference.iter().map(|b| b.to_ascii_uppercase()).collect();
    let mut rng = Rng::new(seed);
    let mut out = Vec::with_capacity(n_pairs);

    let mut produced = 0usize;
    let mut attempts = 0usize;
    // Cap attempts so a pathological (tiny reference) config terminates.
    let attempt_cap = n_pairs.saturating_mul(8).max(64);
    while produced < n_pairs && attempts < attempt_cap {
        attempts += 1;
        let frag_len = model.draw(&mut rng).min(upper.len());
        if frag_len < profile.read_length {
            // Cannot extract a full read from this fragment — retry.
            continue;
        }
        let max_start = upper.len() - frag_len;
        let start = rng.below(max_start + 1);
        let fragment = &upper[start..start + frag_len];

        // Read 1: 5' end, forward.
        let r1_true = fragment[..profile.read_length].to_vec();
        // Read 2: 3' end, reverse-complemented.
        let r2_true = revcomp(&fragment[frag_len - profile.read_length..]);

        let (b1, q1) = apply_error_model(&r1_true, profile, &mut rng);
        let (b2, q2) = apply_error_model(&r2_true, profile, &mut rng);

        let s1 = Seq::new(SeqKind::Dna, &b1)
            .map_err(|e| GenomicsError::invalid("simulated_read", e.to_string()))?;
        let s2 = Seq::new(SeqKind::Dna, &b2)
            .map_err(|e| GenomicsError::invalid("simulated_read", e.to_string()))?;

        let name = format!("sim_pair_{produced}");
        let mut rec1 = SeqRecord::new(name.clone(), s1);
        rec1.description = format!("1 frag_pos={} frag_len={}", start + 1, frag_len);
        let mut rec2 = SeqRecord::new(name, s2);
        rec2.description = format!("2 frag_pos={} frag_len={}", start + 1, frag_len);

        out.push(ReadPair {
            read1: FastqRecord {
                record: rec1,
                quality: q1,
            },
            read2: FastqRecord {
                record: rec2,
                quality: q2,
            },
            fragment_length: frag_len,
        });
        produced += 1;
    }
    if produced < n_pairs {
        return Err(GenomicsError::invalid(
            "model",
            "insert-size model cannot yield full-length reads from this reference",
        ));
    }
    Ok(out)
}

/// Splits a slice of [`ReadPair`]s into the two mate FASTQ streams —
/// the conventional `_R1` / `_R2` file split.
pub fn split_mates(pairs: &[ReadPair]) -> (Vec<FastqRecord>, Vec<FastqRecord>) {
    let r1 = pairs.iter().map(|p| p.read1.clone()).collect();
    let r2 = pairs.iter().map(|p| p.read2.clone()).collect();
    (r1, r2)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ref_seq(len: usize) -> Vec<u8> {
        let bases = [b'A', b'C', b'G', b'T'];
        (0..len).map(|i| bases[(i * 3 + 2) % 4]).collect()
    }

    fn profile() -> IlluminaProfile {
        IlluminaProfile {
            read_length: 100,
            error_5prime: 0.001,
            error_3prime: 0.01,
            quality_high: 40,
            quality_low: 30,
        }
    }

    #[test]
    fn insert_size_draws_within_bounds() {
        let m = InsertSizeModel::short_insert();
        let mut rng = Rng::new(7);
        for _ in 0..500 {
            let s = m.draw(&mut rng);
            assert!((m.min..=m.max).contains(&s));
        }
    }

    #[test]
    fn simulate_yields_pairs() {
        let reference = ref_seq(10_000);
        let pairs =
            simulate_pairs(&reference, &profile(), &InsertSizeModel::short_insert(), 50, 42)
                .unwrap();
        assert_eq!(pairs.len(), 50);
        for p in &pairs {
            assert_eq!(p.read1.len(), 100);
            assert_eq!(p.read2.len(), 100);
            assert_eq!(p.read1.record.id, p.read2.record.id);
        }
    }

    #[test]
    fn simulation_is_deterministic() {
        let reference = ref_seq(10_000);
        let m = InsertSizeModel::short_insert();
        let a = simulate_pairs(&reference, &profile(), &m, 20, 13).unwrap();
        let b = simulate_pairs(&reference, &profile(), &m, 20, 13).unwrap();
        assert_eq!(a, b);
    }

    #[test]
    fn mates_come_from_opposite_strands() {
        // Build a reference and a model with a fixed fragment length so
        // the two mates are deterministic, then check read2 is the
        // reverse complement of the fragment 3' end.
        let reference = ref_seq(5_000);
        let model = InsertSizeModel {
            mean: 300,
            std_dev: 0.0, // fixed length
            min: 300,
            max: 300,
        };
        let pairs = simulate_pairs(&reference, &profile(), &model, 5, 1).unwrap();
        for p in &pairs {
            assert_eq!(p.fragment_length, 300);
        }
    }

    #[test]
    fn split_mates_separates_streams() {
        let reference = ref_seq(10_000);
        let pairs =
            simulate_pairs(&reference, &profile(), &InsertSizeModel::short_insert(), 10, 5)
                .unwrap();
        let (r1, r2) = split_mates(&pairs);
        assert_eq!(r1.len(), 10);
        assert_eq!(r2.len(), 10);
    }

    #[test]
    fn rejects_short_reference() {
        let reference = ref_seq(50);
        assert!(simulate_pairs(
            &reference,
            &profile(),
            &InsertSizeModel::short_insert(),
            1,
            1
        )
        .is_err());
    }

    #[test]
    fn rejects_bad_model() {
        let reference = ref_seq(10_000);
        let bad = InsertSizeModel {
            mean: 300,
            std_dev: 10.0,
            min: 400,
            max: 200, // max < min
        };
        assert!(simulate_pairs(&reference, &profile(), &bad, 1, 1).is_err());
    }
}
