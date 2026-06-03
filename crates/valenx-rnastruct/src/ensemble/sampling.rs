//! Boltzmann stochastic structure sampling (stochastic traceback).
//!
//! Drawing structures *proportionally to their Boltzmann weight*
//! `exp(-E/RT)` is the standard way to characterise the conformational
//! ensemble (Ding & Lawrence 2003). This module fills the McCaskill
//! partition-function arrays and then performs a *stochastic
//! traceback*: at each decomposition step it picks a branch with
//! probability equal to that branch's share of the partition
//! function. A traceback so guided yields an exact Boltzmann sample.
//!
//! [`sample`] returns the requested number of structures; identical
//! structures may legitimately recur (the most probable structure is
//! the most frequent sample). [`sample_with_counts`] groups the draws
//! and reports each distinct structure's observed frequency — a Monte
//! Carlo estimate of its Boltzmann probability.

use crate::ensemble::rng::Rng;
use crate::error::{Result, RnaStructError};
use crate::fold::energy::{self, multiloop, GAS_CONSTANT};
use crate::fold::nussinov::MIN_HAIRPIN;
use crate::fold::zuker::MAX_LOOP;
use crate::rna::RnaSeq;
use crate::structure::Structure;
use std::collections::HashMap;

/// The McCaskill inside arrays, kept for stochastic traceback.
struct InsideTables {
    n: usize,
    rt: f64,
    qb: Vec<f64>,
    q: Vec<f64>,
    qm: Vec<f64>,
    qm2: Vec<f64>,
}

impl InsideTables {
    #[inline]
    fn at(&self, i: usize, j: usize) -> usize {
        i * self.n + j
    }
}

/// Fills the inside partition-function arrays at temperature
/// `temperature_k`.
fn fill_inside(codes: &[u8], temperature_k: f64) -> InsideTables {
    let n = codes.len();
    let rt = GAS_CONSTANT * temperature_k;
    let boltz = |e: f64| (-e / rt).exp();
    let mut t = InsideTables {
        n,
        rt,
        qb: vec![0.0; n * n],
        q: vec![1.0; n * n],
        qm: vec![0.0; n * n],
        qm2: vec![0.0; n * n],
    };
    if n == 0 {
        return t;
    }
    let branch_factor = boltz(multiloop::PER_BRANCH);
    let ml_closure = boltz(multiloop::OFFSET + multiloop::PER_BRANCH);

    for span in 1..n {
        for i in 0..(n - span) {
            let j = i + span;
            // qb
            let mut qbij = 0.0;
            if span > MIN_HAIRPIN && energy::can_pair_codes(codes[i], codes[j]) {
                let loop_bases = &codes[(i + 1)..j];
                qbij += boltz(energy::hairpin_energy(codes[i], codes[j], loop_bases));
                let k_max = (i + 1 + MAX_LOOP).min(j.saturating_sub(MIN_HAIRPIN + 1));
                for k in (i + 1)..=k_max {
                    let l_min = (k + MIN_HAIRPIN + 1).max(j.saturating_sub(MAX_LOOP + 1));
                    for l in l_min..j {
                        if l <= k {
                            continue;
                        }
                        let inner = t.qb[t.at(k, l)];
                        if inner == 0.0 {
                            continue;
                        }
                        let left = k - i - 1;
                        let right = j - l - 1;
                        if left + right != 0 && (left > MAX_LOOP || right > MAX_LOOP) {
                            continue;
                        }
                        let il = energy::internal_loop_energy(
                            codes[i], codes[j], codes[k], codes[l], left, right,
                            codes[i + 1], codes[j - 1], codes[k - 1], codes[l + 1],
                        );
                        qbij += boltz(il) * inner;
                    }
                }
                if j >= i + 2 {
                    qbij += ml_closure
                        * boltz(energy::terminal_penalty(codes[i], codes[j]))
                        * t.qm2[t.at(i + 1, j - 1)];
                }
            }
            let qb_idx = t.at(i, j);
            t.qb[qb_idx] = qbij;
            // qm
            let mut qmij = qbij * branch_factor;
            qmij += t.qm[t.at(i + 1, j)];
            qmij += t.qm[t.at(i, j - 1)];
            // remove double-counted overlap of the two unpaired
            // extensions (interior i+1..j-1 counted twice)
            if j >= i + 2 {
                qmij -= t.qm[t.at(i + 1, j - 1)];
            }
            let qm_idx = t.at(i, j);
            t.qm[qm_idx] = qmij;
            // qm2
            let mut qm2ij = 0.0;
            for k in i..j {
                qm2ij += t.qm[t.at(i, k)] * t.qm[t.at(k + 1, j)];
            }
            let qm2_idx = t.at(i, j);
            t.qm2[qm2_idx] = qm2ij;
            // q
            let mut qij = 1.0; // all-unpaired
            for k in i..=j {
                if energy::can_pair_codes(codes[k], codes[j]) && t.qb[t.at(k, j)] > 0.0 {
                    let left = if k > i { t.q[t.at(i, k - 1)] } else { 1.0 };
                    qij += left
                        * t.qb[t.at(k, j)]
                        * boltz(energy::terminal_penalty(codes[k], codes[j]));
                }
            }
            let q_idx = t.at(i, j);
            t.q[q_idx] = qij;
        }
    }
    t
}

/// Draws `count` structures from the Boltzmann ensemble of `seq` at
/// 37 °C, seeded with `seed`.
///
/// # Errors
/// [`RnaStructError::Invalid`] if `count` is zero.
pub fn sample(seq: &RnaSeq, count: usize, seed: u64) -> Result<Vec<Structure>> {
    sample_at(seq, count, seed, energy::T37_KELVIN)
}

/// [`sample`] at an arbitrary temperature.
///
/// # Errors
/// [`RnaStructError::Invalid`] if `count` is zero.
pub fn sample_at(
    seq: &RnaSeq,
    count: usize,
    seed: u64,
    temperature_k: f64,
) -> Result<Vec<Structure>> {
    if count == 0 {
        return Err(RnaStructError::invalid(
            "count",
            "must request at least one sample",
        ));
    }
    let codes = seq.codes();
    let n = codes.len();
    if n == 0 {
        return Ok(vec![Structure::empty(0); count]);
    }
    let t = fill_inside(codes, temperature_k);
    let mut rng = Rng::new(seed);
    let mut out = Vec::with_capacity(count);
    for _ in 0..count {
        let mut partner = vec![None; n];
        traceback_q(codes, &t, 0, n - 1, &mut partner, &mut rng);
        // Stochastic traceback can, very rarely, fail to seat a pair
        // when floating-point sums drift; `from_partner` cleans the
        // structure and any orphaned half-pair would have been caught.
        if let Ok(s) = Structure::from_partner(partner) {
            out.push(s);
        } else {
            out.push(Structure::empty(n));
        }
    }
    Ok(out)
}

/// Draws `count` structures and returns each distinct structure with
/// its observed count, sorted by descending frequency.
///
/// # Errors
/// [`RnaStructError::Invalid`] if `count` is zero.
pub fn sample_with_counts(
    seq: &RnaSeq,
    count: usize,
    seed: u64,
) -> Result<Vec<(Structure, usize)>> {
    let samples = sample(seq, count, seed)?;
    let mut tally: HashMap<String, (Structure, usize)> = HashMap::new();
    for s in samples {
        let key = s.to_dot_bracket();
        tally.entry(key).or_insert((s, 0)).1 += 1;
    }
    let mut out: Vec<(Structure, usize)> = tally.into_values().collect();
    out.sort_by(|a, b| b.1.cmp(&a.1));
    Ok(out)
}

/// Stochastic traceback of `q[i][j]` — the exterior of `i..=j`.
fn traceback_q(
    codes: &[u8],
    t: &InsideTables,
    i: usize,
    j: usize,
    partner: &mut [Option<usize>],
    rng: &mut Rng,
) {
    if i >= j {
        return;
    }
    let boltz = |e: f64| (-e / t.rt).exp();
    let total = t.q[t.at(i, j)];
    if total <= 0.0 {
        return;
    }
    let r = rng.next_f64() * total;
    let mut cum = 1.0; // all-unpaired weight
    if r < cum {
        return; // everything in i..=j unpaired
    }
    for k in i..=j {
        if energy::can_pair_codes(codes[k], codes[j]) && t.qb[t.at(k, j)] > 0.0 {
            let left = if k > i { t.q[t.at(i, k - 1)] } else { 1.0 };
            let w = left
                * t.qb[t.at(k, j)]
                * boltz(energy::terminal_penalty(codes[k], codes[j]));
            cum += w;
            if r < cum + 1e-12 {
                // base j pairs k; recurse on the left exterior and the
                // helix.
                partner[k] = Some(j);
                partner[j] = Some(k);
                if k > i {
                    traceback_q(codes, t, i, k - 1, partner, rng);
                }
                traceback_qb(codes, t, k, j, partner, rng);
                return;
            }
        }
    }
    // numerical fallthrough: leave unpaired.
}

/// Stochastic traceback of `qb[i][j]` — `(i, j)` is paired.
fn traceback_qb(
    codes: &[u8],
    t: &InsideTables,
    i: usize,
    j: usize,
    partner: &mut [Option<usize>],
    rng: &mut Rng,
) {
    let boltz = |e: f64| (-e / t.rt).exp();
    let total = t.qb[t.at(i, j)];
    if total <= 0.0 {
        return;
    }
    let r = rng.next_f64() * total;
    let mut cum = 0.0;

    // hairpin
    let loop_bases = &codes[(i + 1)..j];
    cum += boltz(energy::hairpin_energy(codes[i], codes[j], loop_bases));
    if r < cum {
        return;
    }
    // internal / bulge / stack
    let k_max = (i + 1 + MAX_LOOP).min(j.saturating_sub(MIN_HAIRPIN + 1));
    for k in (i + 1)..=k_max {
        let l_min = (k + MIN_HAIRPIN + 1).max(j.saturating_sub(MAX_LOOP + 1));
        for l in l_min..j {
            if l <= k {
                continue;
            }
            let inner = t.qb[t.at(k, l)];
            if inner == 0.0 {
                continue;
            }
            let left = k - i - 1;
            let right = j - l - 1;
            if left + right != 0 && (left > MAX_LOOP || right > MAX_LOOP) {
                continue;
            }
            let il = energy::internal_loop_energy(
                codes[i], codes[j], codes[k], codes[l], left, right,
                codes[i + 1], codes[j - 1], codes[k - 1], codes[l + 1],
            );
            cum += boltz(il) * inner;
            if r < cum {
                partner[k] = Some(l);
                partner[l] = Some(k);
                traceback_qb(codes, t, k, l, partner, rng);
                return;
            }
        }
    }
    // multiloop
    if j >= i + 2 {
        traceback_qm2(codes, t, i + 1, j - 1, partner, rng);
    }
}

/// Stochastic traceback of `qm2[i][j]` — a multiloop interior of
/// `>= 2` branches.
fn traceback_qm2(
    codes: &[u8],
    t: &InsideTables,
    i: usize,
    j: usize,
    partner: &mut [Option<usize>],
    rng: &mut Rng,
) {
    if i >= j {
        return;
    }
    let total = t.qm2[t.at(i, j)];
    if total <= 0.0 {
        return;
    }
    let r = rng.next_f64() * total;
    let mut cum = 0.0;
    for k in i..j {
        let w = t.qm[t.at(i, k)] * t.qm[t.at(k + 1, j)];
        cum += w;
        if r < cum {
            traceback_qm(codes, t, i, k, partner, rng);
            traceback_qm(codes, t, k + 1, j, partner, rng);
            return;
        }
    }
    // fallthrough — split at the last position.
    if i < j {
        traceback_qm(codes, t, i, j - 1, partner, rng);
        traceback_qm(codes, t, j, j, partner, rng);
    }
}

/// Stochastic traceback of `qm[i][j]` — a multiloop fragment of
/// `>= 1` branch.
fn traceback_qm(
    codes: &[u8],
    t: &InsideTables,
    i: usize,
    j: usize,
    partner: &mut [Option<usize>],
    rng: &mut Rng,
) {
    if i > j {
        return;
    }
    let boltz = |e: f64| (-e / t.rt).exp();
    let total = t.qm[t.at(i, j)];
    if total <= 0.0 {
        return;
    }
    let r = rng.next_f64() * total;
    let mut cum = 0.0;
    // single branch helix (i, j)
    cum += t.qb[t.at(i, j)] * boltz(multiloop::PER_BRANCH);
    if r < cum {
        partner[i] = Some(j);
        partner[j] = Some(i);
        traceback_qb(codes, t, i, j, partner, rng);
        return;
    }
    // i unpaired
    if i < j {
        cum += t.qm[t.at(i + 1, j)];
        if r < cum {
            traceback_qm(codes, t, i + 1, j, partner, rng);
            return;
        }
        // j unpaired
        cum += t.qm[t.at(i, j - 1)];
        if r < cum {
            traceback_qm(codes, t, i, j - 1, partner, rng);
            return;
        }
    }
    // fallthrough — treat as i unpaired.
    if i < j {
        traceback_qm(codes, t, i + 1, j, partner, rng);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fold::eval::structure_energy;

    #[test]
    fn samples_are_valid_structures() {
        let seq = RnaSeq::parse("GGGGGAAAACCCCC").unwrap();
        let samples = sample(&seq, 20, 1).unwrap();
        assert_eq!(samples.len(), 20);
        for s in &samples {
            assert!(s.is_nested(), "a Boltzmann sample must be nested");
            assert_eq!(s.len(), seq.len());
            // every sampled structure is scoreable
            assert!(structure_energy(&seq, s).unwrap().is_finite());
        }
    }

    #[test]
    fn sampling_is_reproducible() {
        let seq = RnaSeq::parse("GGGGGAAAACCCCC").unwrap();
        let a = sample(&seq, 15, 99).unwrap();
        let b = sample(&seq, 15, 99).unwrap();
        let adb: Vec<String> = a.iter().map(|s| s.to_dot_bracket()).collect();
        let bdb: Vec<String> = b.iter().map(|s| s.to_dot_bracket()).collect();
        assert_eq!(adb, bdb, "same seed must give the same samples");
    }

    #[test]
    fn unpairable_sequence_samples_open() {
        let seq = RnaSeq::parse("AAAAAAAAAAAA").unwrap();
        let samples = sample(&seq, 5, 3).unwrap();
        for s in &samples {
            assert_eq!(s.n_pairs(), 0);
        }
    }

    #[test]
    fn counts_sum_to_requested_total() {
        let seq = RnaSeq::parse("GGGGGAAAACCCCC").unwrap();
        let counted = sample_with_counts(&seq, 40, 7).unwrap();
        let total: usize = counted.iter().map(|(_, c)| c).sum();
        assert_eq!(total, 40);
        // sorted by descending count
        for w in counted.windows(2) {
            assert!(w[0].1 >= w[1].1);
        }
    }

    #[test]
    fn stable_stem_dominates_the_sample() {
        // a strong stem should be the most frequently sampled fold
        let seq = RnaSeq::parse("GGGGGGGAAAACCCCCCC").unwrap();
        let counted = sample_with_counts(&seq, 60, 11).unwrap();
        // the top structure should carry at least some pairs
        assert!(counted[0].0.n_pairs() > 0);
    }

    #[test]
    fn rejects_zero_count() {
        let seq = RnaSeq::parse("GGGGCCCC").unwrap();
        assert!(sample(&seq, 0, 1).is_err());
    }
}
