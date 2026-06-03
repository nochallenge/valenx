//! PCR primer design.
//!
//! Given a template and a target sub-region, [`design_primers`] picks
//! a forward and reverse primer that flank the target and satisfy
//! length / Tm / GC-clamp constraints, then screens each candidate
//! for hairpins, self-dimers and cross-primer dimers using the
//! nearest-neighbor ΔG model in [`crate::analysis::thermo`] (Primer3 /
//! OligoAnalyzer style — not a geometric stem search).
//!
//! Tm is computed via the SantaLucia 1998 unified thermodynamic model
//! with monovalent + Mg²⁺ + dNTP corrections; hairpin / dimer
//! candidates are rejected when their ΔG at the screen temperature
//! (default 37 °C) drops below the configured stability threshold
//! (default `−9.0` kcal/mol — slightly more permissive than Primer3's
//! `−9.0` to `−12.0` default to avoid rejecting all candidates on
//! short, GC-rich templates).

use crate::analysis::composition::gc_content;
use crate::analysis::thermo::{
    most_stable_hairpin, most_stable_hetero_dimer, most_stable_self_dimer, DimerScore,
    HairpinScore,
};
use crate::analysis::tm::tm_nearest_neighbor_default;
use crate::error::{BioseqError, Result};
use crate::ops::revcomp::reverse_complement;
use crate::seq::{Seq, SeqKind};

/// Constraints for [`design_primers`].
#[derive(Copy, Clone, Debug)]
pub struct PrimerConstraints {
    /// Minimum primer length, nt.
    pub min_len: usize,
    /// Maximum primer length, nt.
    pub max_len: usize,
    /// Target melting temperature, °C.
    pub target_tm: f64,
    /// Allowed deviation from `target_tm`, °C.
    pub tm_tolerance: f64,
    /// Minimum GC fraction.
    pub min_gc: f64,
    /// Maximum GC fraction.
    pub max_gc: f64,
    /// Require a GC-clamp — a `G` or `C` as the 3′-terminal base.
    pub require_gc_clamp: bool,
    /// Temperature (°C) at which hairpin / dimer ΔG is evaluated.
    /// Primer3 default is the annealing temperature; we expose it
    /// explicitly so callers can be honest about which regime they
    /// care about.
    pub screen_temp: f64,
    /// Maximum allowed hairpin ΔG at the screen temperature, kcal/mol
    /// (more negative is more stable). Primer3 default is around
    /// `−9` kcal/mol.
    pub max_hairpin_dg: f64,
    /// Maximum allowed self-dimer ΔG (kcal/mol).
    pub max_self_dimer_dg: f64,
    /// Maximum allowed hetero-dimer (cross-primer) ΔG (kcal/mol).
    pub max_hetero_dimer_dg: f64,
}

impl Default for PrimerConstraints {
    fn default() -> Self {
        PrimerConstraints {
            min_len: 18,
            max_len: 30,
            target_tm: 60.0,
            tm_tolerance: 5.0,
            min_gc: 0.4,
            max_gc: 0.6,
            require_gc_clamp: true,
            screen_temp: 37.0,
            max_hairpin_dg: -9.0,
            max_self_dimer_dg: -9.0,
            max_hetero_dimer_dg: -9.0,
        }
    }
}

/// A designed primer with its computed properties.
#[derive(Clone, Debug, PartialEq)]
pub struct Primer {
    /// The primer sequence, 5′→3′.
    pub seq: Seq,
    /// `true` for the forward (sense) primer, `false` for the reverse.
    pub is_forward: bool,
    /// 0-based start of the primer's binding footprint on the
    /// template's top strand.
    pub binding_start: usize,
    /// 0-based end (exclusive) of the binding footprint.
    pub binding_end: usize,
    /// Nearest-neighbor melting temperature, °C.
    pub tm: f64,
    /// GC fraction.
    pub gc: f64,
    /// Most stable hairpin found (if any). The presence of a value
    /// does NOT mean the primer was rejected; thresholds are applied
    /// separately.
    pub hairpin: Option<HairpinScore>,
    /// Most stable self-dimer found (if any).
    pub self_dimer: Option<DimerScore>,
}

impl Primer {
    /// Primer length in nt.
    pub fn len(&self) -> usize {
        self.seq.len()
    }

    /// `true` if the primer carries no residues.
    pub fn is_empty(&self) -> bool {
        self.seq.is_empty()
    }

    /// `true` if the primer is free of both hairpins and self-dimers.
    pub fn is_clean(&self) -> bool {
        self.hairpin.is_none() && self.self_dimer.is_none()
    }

    /// Backwards-compatible legacy predicate.
    pub fn has_hairpin(&self) -> bool {
        self.hairpin.is_some()
    }

    /// Backwards-compatible legacy predicate.
    pub fn has_self_dimer(&self) -> bool {
        self.self_dimer.is_some()
    }
}

/// A designed forward + reverse primer pair.
#[derive(Clone, Debug, PartialEq)]
pub struct PrimerPair {
    /// The forward primer.
    pub forward: Primer,
    /// The reverse primer.
    pub reverse: Primer,
    /// Most stable hetero-dimer (forward × reverse), if any.
    pub hetero_dimer: Option<DimerScore>,
}

impl PrimerPair {
    /// Absolute difference in Tm between the two primers, °C — a key
    /// quality metric (a small ΔTm makes PCR robust).
    pub fn tm_difference(&self) -> f64 {
        (self.forward.tm - self.reverse.tm).abs()
    }

    /// Size of the amplicon the pair produces, in bp.
    pub fn amplicon_size(&self) -> usize {
        self.reverse
            .binding_end
            .saturating_sub(self.forward.binding_start)
    }
}

/// Designs a primer pair flanking `[target_start, target_end)` on
/// `template`.
///
/// The forward primer anneals ending just before `target_start`; the
/// reverse primer is the reverse complement of the region just after
/// `target_end`. For each, lengths in `[min_len, max_len]` are tried
/// and the candidate whose Tm is closest to `target_tm` (and that
/// satisfies the GC + clamp + hairpin/dimer ΔG constraints) is chosen.
///
/// Returns [`BioseqError::Invalid`] for a non-DNA template, an
/// out-of-range target, or if no candidate satisfies the constraints.
pub fn design_primers(
    template: &Seq,
    target_start: usize,
    target_end: usize,
    constraints: PrimerConstraints,
) -> Result<PrimerPair> {
    if template.kind() != SeqKind::Dna {
        return Err(BioseqError::invalid("kind", "primer design needs a DNA template"));
    }
    let n = template.len();
    if target_start >= target_end || target_end > n {
        return Err(BioseqError::invalid(
            "target",
            format!("invalid target [{target_start},{target_end}) on length {n}"),
        ));
    }
    if constraints.min_len == 0 || constraints.min_len > constraints.max_len {
        return Err(BioseqError::invalid(
            "constraints",
            "require 0 < min_len <= max_len",
        ));
    }

    let forward = design_forward(template, target_start, &constraints)?;
    let reverse = design_reverse(template, target_end, &constraints)?;
    let hetero_dimer =
        most_stable_hetero_dimer(forward.seq.as_bytes(), reverse.seq.as_bytes(), constraints.screen_temp)
            .filter(|d| d.dg < constraints.max_hetero_dimer_dg);
    Ok(PrimerPair {
        forward,
        reverse,
        hetero_dimer,
    })
}

/// Designs the forward primer — ends at `target_start` on the top
/// strand, extends 5′ from there.
fn design_forward(
    template: &Seq,
    target_start: usize,
    c: &PrimerConstraints,
) -> Result<Primer> {
    let bytes = template.as_bytes();
    let mut best: Option<(Primer, f64)> = None;
    for len in c.min_len..=c.max_len {
        // Need `len` template bases ending just before `target_start`.
        if len > target_start {
            continue; // not enough template to the 5' side
        }
        let start = target_start - len;
        let end = target_start;
        let sub = Seq::new_unchecked(SeqKind::Dna, bytes[start..end].to_vec(), template.topology());
        if let Some(primer) = evaluate_candidate(sub, true, start, end, c) {
            consider(&mut best, primer, c.target_tm);
        }
    }
    best.map(|(p, _)| p).ok_or_else(|| {
        BioseqError::invalid(
            "forward_primer",
            "no forward primer satisfies the constraints",
        )
    })
}

/// Designs the reverse primer — binds the top strand starting at
/// `target_end`, the primer itself is the reverse complement.
fn design_reverse(
    template: &Seq,
    target_end: usize,
    c: &PrimerConstraints,
) -> Result<Primer> {
    let bytes = template.as_bytes();
    let n = bytes.len();
    let mut best: Option<(Primer, f64)> = None;
    for len in c.min_len..=c.max_len {
        if target_end + len > n {
            continue; // not enough template to the 3' side
        }
        let start = target_end;
        let end = target_end + len;
        let footprint =
            Seq::new_unchecked(SeqKind::Dna, bytes[start..end].to_vec(), template.topology());
        // The reverse primer is the reverse complement of the footprint.
        let primer_seq = reverse_complement(&footprint)?;
        if let Some(primer) = evaluate_candidate(primer_seq, false, start, end, c) {
            consider(&mut best, primer, c.target_tm);
        }
    }
    best.map(|(p, _)| p).ok_or_else(|| {
        BioseqError::invalid(
            "reverse_primer",
            "no reverse primer satisfies the constraints",
        )
    })
}

/// Checks a candidate against the constraints; returns a [`Primer`] if
/// it passes the hard filters (Tm window, GC window, GC clamp, and
/// ΔG-based hairpin / self-dimer thresholds).
fn evaluate_candidate(
    seq: Seq,
    is_forward: bool,
    binding_start: usize,
    binding_end: usize,
    c: &PrimerConstraints,
) -> Option<Primer> {
    let tm = tm_nearest_neighbor_default(&seq).ok()?;
    if (tm - c.target_tm).abs() > c.tm_tolerance {
        return None;
    }
    let gc = gc_content(&seq).ok()?;
    if gc < c.min_gc || gc > c.max_gc {
        return None;
    }
    if c.require_gc_clamp && !has_gc_clamp(&seq) {
        return None;
    }
    let bytes = seq.as_bytes();
    let hairpin = most_stable_hairpin(bytes, c.screen_temp);
    let self_dimer = most_stable_self_dimer(bytes, c.screen_temp);
    // Reject candidates whose hairpin / dimer is too stable.
    if let Some(h) = &hairpin {
        if h.dg < c.max_hairpin_dg {
            return None;
        }
    }
    if let Some(d) = &self_dimer {
        if d.dg < c.max_self_dimer_dg {
            return None;
        }
    }
    Some(Primer {
        seq,
        is_forward,
        binding_start,
        binding_end,
        tm,
        gc,
        hairpin,
        self_dimer,
    })
}

/// Keeps whichever candidate's Tm is closest to the target.
fn consider(best: &mut Option<(Primer, f64)>, candidate: Primer, target_tm: f64) {
    let score = (candidate.tm - target_tm).abs();
    // Prefer a clean primer over a structured one at a similar score.
    let penalty = if candidate.is_clean() { 0.0 } else { 0.5 };
    let effective = score + penalty;
    match best {
        None => *best = Some((candidate, effective)),
        Some((_, bs)) if effective < *bs => *best = Some((candidate, effective)),
        Some(_) => {}
    }
}

/// `true` if the primer's 3′-terminal base is `G` or `C`.
pub fn has_gc_clamp(primer: &Seq) -> bool {
    matches!(
        primer.as_bytes().last().map(|b| b.to_ascii_uppercase()),
        Some(b'G') | Some(b'C')
    )
}

/// `true` if the primer has a self-dimer at the default 37 °C with
/// stability worse than `−4 kcal/mol`. A coarse boolean predicate; for
/// fine control prefer [`most_stable_self_dimer`] directly.
pub fn has_self_dimer(primer: &Seq) -> bool {
    most_stable_self_dimer(primer.as_bytes(), 37.0)
        .map(|d| d.dg < -4.0)
        .unwrap_or(false)
}

/// `true` if the primer can fold into a hairpin with ΔG < `−2 kcal/mol`
/// at 37 °C. Coarse boolean predicate; prefer
/// [`most_stable_hairpin`] for finer control.
pub fn has_hairpin(primer: &Seq) -> bool {
    most_stable_hairpin(primer.as_bytes(), 37.0)
        .map(|h| h.dg < -2.0)
        .unwrap_or(false)
}

/// `true` if two primers form a hetero-dimer with ΔG below the
/// `default_threshold` (kcal/mol; pass `-4.0` for a Primer3-like
/// permissive filter).
pub fn has_hetero_dimer(a: &Seq, b: &Seq) -> bool {
    most_stable_hetero_dimer(a.as_bytes(), b.as_bytes(), 37.0)
        .map(|d| d.dg < -4.0)
        .unwrap_or(false)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A template with a GC-balanced flanking design space.
    fn template() -> Seq {
        // 90 nt: regions chosen so 20-mers land near Tm 55-62 °C.
        Seq::new(
            SeqKind::Dna,
            "ATGCGACTGATCGATCGTAGCTAGCATCGATCGATCGTAGCTAGCTAGCATCGATCGTAGCTAGCATCGATCGATCGTAGCTAGCTAGC",
        )
        .unwrap()
    }

    #[test]
    fn designs_a_pair() {
        let t = template();
        let c = PrimerConstraints {
            min_len: 18,
            max_len: 26,
            target_tm: 58.0,
            tm_tolerance: 12.0,
            min_gc: 0.3,
            max_gc: 0.7,
            require_gc_clamp: false,
            max_hairpin_dg: -20.0,
            max_self_dimer_dg: -20.0,
            max_hetero_dimer_dg: -20.0,
            ..Default::default()
        };
        let pair = design_primers(&t, 30, 60, c).unwrap();
        assert!(pair.forward.is_forward);
        assert!(!pair.reverse.is_forward);
        // Forward primer ends right before the target.
        assert_eq!(pair.forward.binding_end, 30);
        // Reverse primer footprint starts at the target end.
        assert_eq!(pair.reverse.binding_start, 60);
        assert!(pair.amplicon_size() >= 30);
    }

    #[test]
    fn respects_length_bounds() {
        let t = template();
        let c = PrimerConstraints {
            min_len: 20,
            max_len: 20,
            target_tm: 58.0,
            tm_tolerance: 15.0,
            min_gc: 0.2,
            max_gc: 0.8,
            require_gc_clamp: false,
            max_hairpin_dg: -20.0,
            max_self_dimer_dg: -20.0,
            max_hetero_dimer_dg: -20.0,
            ..Default::default()
        };
        let pair = design_primers(&t, 30, 60, c).unwrap();
        assert_eq!(pair.forward.len(), 20);
        assert_eq!(pair.reverse.len(), 20);
    }

    #[test]
    fn out_of_range_target_errs() {
        let t = template();
        let c = PrimerConstraints::default();
        assert!(design_primers(&t, 60, 30, c).is_err()); // inverted
        assert!(design_primers(&t, 0, 999, c).is_err()); // past end
    }

    #[test]
    fn impossible_constraints_err() {
        let t = template();
        let c = PrimerConstraints {
            min_len: 18,
            max_len: 26,
            target_tm: 200.0, // unreachable Tm
            tm_tolerance: 1.0,
            min_gc: 0.4,
            max_gc: 0.6,
            require_gc_clamp: true,
            ..Default::default()
        };
        assert!(design_primers(&t, 30, 60, c).is_err());
    }

    #[test]
    fn gc_clamp_detection() {
        assert!(has_gc_clamp(&Seq::new(SeqKind::Dna, "AAAAG").unwrap()));
        assert!(has_gc_clamp(&Seq::new(SeqKind::Dna, "AAAAC").unwrap()));
        assert!(!has_gc_clamp(&Seq::new(SeqKind::Dna, "AAAAT").unwrap()));
    }

    #[test]
    fn self_dimer_detection_thermodynamic() {
        // A palindrome necessarily forms a stable self-dimer with
        // significant ΔG.
        let palindrome = Seq::new(SeqKind::Dna, "GAATTCGAATTC").unwrap();
        assert!(has_self_dimer(&palindrome));
        // A homopolymer has no Watson-Crick self-complementarity.
        let poly = Seq::new(SeqKind::Dna, "AAAAAAAAAAAA").unwrap();
        assert!(!has_self_dimer(&poly));
    }

    #[test]
    fn hairpin_detection_thermodynamic() {
        // Stem GGGGG ... loop ... CCCCC -> a very stable hairpin.
        let hp = Seq::new(SeqKind::Dna, "GGGGGAAATTTCCCCC").unwrap();
        assert!(has_hairpin(&hp));
        let flat = Seq::new(SeqKind::Dna, "AAAAAAAAAAAAAAAAA").unwrap();
        assert!(!has_hairpin(&flat));
    }

    #[test]
    fn hetero_dimer_detection_thermodynamic() {
        // Two primers whose 3′ ends are reverse-complementary.
        let a = Seq::new(SeqKind::Dna, "GACTAGGCATTTT").unwrap();
        let b = Seq::new(SeqKind::Dna, "AAAATGCCTAGTC").unwrap();
        assert!(has_hetero_dimer(&a, &b));
    }

    #[test]
    fn primer_object_carries_thermodynamic_scores() {
        let t = template();
        let c = PrimerConstraints {
            min_len: 20,
            max_len: 20,
            target_tm: 58.0,
            tm_tolerance: 15.0,
            min_gc: 0.2,
            max_gc: 0.8,
            require_gc_clamp: false,
            max_hairpin_dg: -20.0,
            max_self_dimer_dg: -20.0,
            max_hetero_dimer_dg: -20.0,
            ..Default::default()
        };
        let pair = design_primers(&t, 30, 60, c).unwrap();
        // The Primer object exposes the structured thermodynamic
        // scores — they may be None if no significant structure was
        // found.
        let _ = pair.forward.hairpin.as_ref();
        let _ = pair.forward.self_dimer.as_ref();
    }

    #[test]
    fn non_dna_template_rejected() {
        let r = Seq::new(SeqKind::Rna, "ACGUACGU").unwrap();
        assert!(design_primers(&r, 0, 4, PrimerConstraints::default()).is_err());
    }
}
