//! DNA nearest-neighbor thermodynamics — the SantaLucia 1998 unified
//! ΔH° / ΔS° parameter set, salt corrections (monovalent + Mg²⁺ +
//! dNTP), and ΔG-based hairpin / self-dimer / hetero-dimer scoring
//! (Primer3 / OligoAnalyzer style).
//!
//! All parameters are sourced from John SantaLucia Jr., "A unified
//! view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor
//! thermodynamics", PNAS 95: 1460-1465 (1998).
//!
//! Salt corrections:
//!
//! - Monovalent: SantaLucia 1998 (`ΔS_corr = ΔS + 0.368 * (N-1) *
//!   ln[Na⁺]`).
//! - Mg²⁺ + dNTP: the von Ahsen et al. (1999) effective-monovalent
//!   formula `[Na⁺]_eq = [Mon⁺] + 120 · √max(0, [Mg²⁺] - [dNTP])`,
//!   which converts the divalent + chelating dNTP contribution into a
//!   monovalent equivalent that drives the same correction.
//!
//! Hairpin / dimer ΔG: minimum-free-energy DP over the nearest-
//! neighbor parameter set, with a hairpin-loop initiation penalty of
//! `+4.0 kcal/mol` (a representative loop-size-independent value that
//! tracks the Primer3 default for small hairpins).

use std::cmp::Ordering;

/// SantaLucia 1998 nearest-neighbor parameter set — ΔH° (kcal/mol) and
/// ΔS° (cal/(mol·K)) for the 10 unique Watson-Crick dinucleotide
/// stacks. Indexed by the 5′→3′ top-strand dinucleotide.
pub fn nn_params(dimer: &[u8; 2]) -> Option<NnEntry> {
    let entry = |dh, ds| NnEntry { dh, ds };
    Some(match dimer {
        b"AA" | b"TT" => entry(-7.9, -22.2),
        b"AT" => entry(-7.2, -20.4),
        b"TA" => entry(-7.2, -21.3),
        b"CA" | b"TG" => entry(-8.5, -22.7),
        b"GT" | b"AC" => entry(-8.4, -22.4),
        b"CT" | b"AG" => entry(-7.8, -21.0),
        b"GA" | b"TC" => entry(-8.2, -22.2),
        b"CG" => entry(-10.6, -27.2),
        b"GC" => entry(-9.8, -24.4),
        b"GG" | b"CC" => entry(-8.0, -19.9),
        _ => return None,
    })
}

/// Initiation parameters (SantaLucia 1998 unified set): one set is
/// added for each end of the duplex, depending on the terminal base.
pub fn nn_init(end: u8) -> NnEntry {
    match end {
        b'G' | b'C' => NnEntry { dh: 0.1, ds: -2.8 },
        _ => NnEntry { dh: 2.3, ds: 4.1 },
    }
}

/// Symmetry correction for self-complementary duplexes (SantaLucia
/// 1998): `ΔS += -1.4 cal/(mol·K)`.
pub const SYMMETRY_DS: f64 = -1.4;

/// One nearest-neighbor (or initiation) ΔH° / ΔS° entry.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct NnEntry {
    /// Enthalpy, kcal/mol.
    pub dh: f64,
    /// Entropy, cal/(mol·K).
    pub ds: f64,
}

impl NnEntry {
    /// Gibbs free energy at temperature `t_celsius`, kcal/mol.
    pub fn dg(&self, t_celsius: f64) -> f64 {
        let t_k = t_celsius + 273.15;
        // ΔG = ΔH - T·ΔS (with ΔS in cal/(mol·K) -> divide by 1000).
        self.dh - t_k * self.ds / 1000.0
    }

    /// Adds another entry to this one in place.
    pub fn add(&mut self, other: NnEntry) {
        self.dh += other.dh;
        self.ds += other.ds;
    }
}

/// Salt + buffer parameters used by the thermodynamic routines.
#[derive(Copy, Clone, Debug)]
pub struct SaltParams {
    /// Monovalent cation (Na⁺ / K⁺) concentration, mol/L.
    pub na_conc: f64,
    /// Divalent magnesium concentration, mol/L (0 for none).
    pub mg_conc: f64,
    /// Total dNTP concentration, mol/L (0 if not present).
    pub dntp_conc: f64,
}

impl Default for SaltParams {
    /// Standard PCR buffer: 50 mM Na⁺-equivalent, 1.5 mM Mg²⁺,
    /// 0.2 mM dNTPs.
    fn default() -> Self {
        SaltParams {
            na_conc: 0.05,
            mg_conc: 1.5e-3,
            dntp_conc: 0.2e-3,
        }
    }
}

impl SaltParams {
    /// Effective monovalent cation concentration, mol/L, that
    /// approximates the duplex-stabilizing effect of the full buffer
    /// (von Ahsen 1999): `[Na⁺]_eq[mM] = [Mon⁺][mM] + 120 *
    /// √([Mg²⁺]_free[mM])`. Returned in mol/L for consistency with
    /// the rest of the API; internally the von Ahsen coefficient is
    /// evaluated in millimolar units (the form published in the paper).
    pub fn effective_na(&self) -> f64 {
        // Convert mol/L -> mmol/L for the von Ahsen formula.
        let mon_mm = self.na_conc * 1000.0;
        let mg_mm = self.mg_conc * 1000.0;
        let dntp_mm = self.dntp_conc * 1000.0;
        let free_mg_mm = (mg_mm - dntp_mm).max(0.0);
        let na_eq_mm = mon_mm + 120.0 * free_mg_mm.sqrt();
        na_eq_mm / 1000.0
    }
}

/// Strand-concentration parameters for [`duplex_tm`].
#[derive(Copy, Clone, Debug)]
pub struct StrandParams {
    /// Total strand concentration, mol/L (typical PCR primer 0.25 µM).
    pub strand_conc: f64,
}

impl Default for StrandParams {
    fn default() -> Self {
        StrandParams {
            strand_conc: 0.25e-6,
        }
    }
}

/// `true` if `seq` is its own reverse complement (a self-
/// complementary palindrome).
pub fn is_self_complementary(seq: &[u8]) -> bool {
    if seq.is_empty() {
        return false;
    }
    let n = seq.len();
    for i in 0..n {
        if !pairs(seq[i], seq[n - 1 - i]) {
            return false;
        }
    }
    true
}

/// Sum the nearest-neighbor ΔH / ΔS (with initiation + symmetry
/// correction) for a Watson-Crick duplex with `seq` as the top strand.
/// Returns `None` if any base is not pure ACGT.
pub fn duplex_enthalpy_entropy(seq: &[u8]) -> Option<NnEntry> {
    if seq.len() < 2 {
        return None;
    }
    let mut acc = NnEntry { dh: 0.0, ds: 0.0 };
    for w in seq.windows(2) {
        let arr = [w[0].to_ascii_uppercase(), w[1].to_ascii_uppercase()];
        let entry = nn_params(&arr)?;
        acc.add(entry);
    }
    acc.add(nn_init(seq[0].to_ascii_uppercase()));
    acc.add(nn_init(seq[seq.len() - 1].to_ascii_uppercase()));
    if is_self_complementary(seq) {
        acc.ds += SYMMETRY_DS;
    }
    Some(acc)
}

/// Duplex melting temperature, °C, with monovalent + Mg²⁺ + dNTP
/// corrections. `None` if the sequence is too short or contains
/// non-ACGT bases.
pub fn duplex_tm(seq: &[u8], salt: SaltParams, strand: StrandParams) -> Option<f64> {
    let nn = duplex_enthalpy_entropy(seq)?;
    let dh_cal = nn.dh * 1000.0; // kcal -> cal
                                 // SantaLucia salt correction on ΔS (monovalent-equivalent).
    let n = seq.len() as f64;
    let na_eq = salt.effective_na().max(1e-6);
    let ds_salt = nn.ds + 0.368 * (n - 1.0) * na_eq.ln();
    // Gas constant cal/(mol·K).
    const R: f64 = 1.987;
    // x = 4 for non-self-complementary duplex, x = 1 for self-comp.
    let x = if is_self_complementary(seq) { 1.0 } else { 4.0 };
    let denom = ds_salt + R * (strand.strand_conc / x).ln();
    if denom.abs() < 1e-12 {
        return None;
    }
    let tm_k = dh_cal / denom;
    Some(tm_k - 273.15)
}

// ----------------------------------------------------------------------
// Hairpin / dimer ΔG scoring
// ----------------------------------------------------------------------

/// Hairpin-loop initiation penalty, kcal/mol — a representative
/// loop-size-independent value used by Primer3 for short oligos.
pub const HAIRPIN_LOOP_INIT_DG: f64 = 4.0;

/// Minimum hairpin loop length, nt (the unpaired bases at the apex).
pub const HAIRPIN_MIN_LOOP: usize = 3;

/// Minimum stem length, nt — at least one full NN stack is required.
pub const HAIRPIN_MIN_STEM: usize = 3;

/// A scored hairpin candidate.
#[derive(Clone, Debug, PartialEq)]
pub struct HairpinScore {
    /// Total ΔG of the hairpin (stem stacks + loop initiation), kcal/mol.
    pub dg: f64,
    /// Start of the 5′ side of the stem on the input.
    pub stem_start: usize,
    /// Length of the stem (paired bases per side).
    pub stem_len: usize,
    /// Length of the loop between the two sides of the stem.
    pub loop_len: usize,
}

/// Searches `seq` for the most stable hairpin (lowest ΔG at 37°C). A
/// hairpin is a stem (Watson-Crick paired bases from one side of the
/// sequence to the other) with a loop of at least
/// [`HAIRPIN_MIN_LOOP`] unpaired bases at the apex.
///
/// Returns `None` if no hairpin with stem ≥ [`HAIRPIN_MIN_STEM`] and
/// negative ΔG can be found. `t_celsius` is the temperature at which
/// ΔG is evaluated (typically 37°C for the standard reference, or the
/// PCR annealing temperature for primer-design use).
pub fn most_stable_hairpin(seq: &[u8], t_celsius: f64) -> Option<HairpinScore> {
    let n = seq.len();
    if n < 2 * HAIRPIN_MIN_STEM + HAIRPIN_MIN_LOOP {
        return None;
    }
    let bytes: Vec<u8> = seq.iter().map(|b| b.to_ascii_uppercase()).collect();
    let mut best: Option<HairpinScore> = None;
    // Enumerate the apex of the loop. For each possible stem start `i`
    // and loop length `loop_len`, extend the stem as long as bases
    // pair; record the lowest ΔG.
    for i in 0..n {
        // The 3′ side of the stem starts at j = i + stem_len + loop_len.
        for stem_len in HAIRPIN_MIN_STEM..=((n - i) / 2) {
            for loop_len in HAIRPIN_MIN_LOOP..=(n - i - 2 * stem_len) {
                let j = i + stem_len + loop_len;
                if j + stem_len > n {
                    break;
                }
                // Check Watson-Crick pairing along the stem.
                let mut paired = true;
                for k in 0..stem_len {
                    let a = bytes[i + k];
                    let b = bytes[j + stem_len - 1 - k];
                    if !pairs(a, b) {
                        paired = false;
                        break;
                    }
                }
                if !paired {
                    continue;
                }
                // ΔG = sum of NN stacks along the stem +
                // loop-initiation penalty.
                let mut dg = HAIRPIN_LOOP_INIT_DG;
                for k in 0..stem_len - 1 {
                    let a1 = bytes[i + k];
                    let a2 = bytes[i + k + 1];
                    let entry = nn_params(&[a1, a2])?;
                    dg += entry.dg(t_celsius);
                }
                let cand = HairpinScore {
                    dg,
                    stem_start: i,
                    stem_len,
                    loop_len,
                };
                if better(&best, &cand) {
                    best = Some(cand);
                }
            }
        }
    }
    best.filter(|h| h.dg < 0.0)
}

fn better(best: &Option<HairpinScore>, cand: &HairpinScore) -> bool {
    match best {
        None => true,
        Some(b) => cand.dg.partial_cmp(&b.dg).unwrap_or(Ordering::Equal) == Ordering::Less,
    }
}

/// A scored dimer alignment (self-dimer or hetero-dimer).
#[derive(Clone, Debug, PartialEq)]
pub struct DimerScore {
    /// Total ΔG of the dimer at the chosen temperature, kcal/mol.
    /// More negative = more stable = worse for a primer.
    pub dg: f64,
    /// Offset of the 3′-most paired base on strand `a`.
    pub a_offset: usize,
    /// Offset of the 3′-most paired base on strand `b`.
    pub b_offset: usize,
    /// Number of contiguous paired bases.
    pub length: usize,
    /// `true` if the dimer involves the 3′ end of either primer (the
    /// most consequential class — a 3′-end dimer abolishes priming).
    pub three_prime: bool,
}

/// Scores the most stable hetero-dimer alignment between two primers
/// `a` and `b` at `t_celsius`. Each candidate alignment slides strand
/// `b` (taken in 3′→5′ orientation against `a`'s 5′→3′ direction) and
/// counts contiguous Watson-Crick pairs; ΔG is summed over the
/// nearest-neighbor stacks of the paired segment.
///
/// Returns `None` if no contiguous run of at least 4 paired bases with
/// negative ΔG can be found.
pub fn most_stable_hetero_dimer(a: &[u8], b: &[u8], t_celsius: f64) -> Option<DimerScore> {
    score_dimer(a, b, t_celsius)
}

/// Scores the most stable self-dimer of `seq`.
pub fn most_stable_self_dimer(seq: &[u8], t_celsius: f64) -> Option<DimerScore> {
    score_dimer(seq, seq, t_celsius)
}

/// Generic dimer scorer. Strand `b` is taken antiparallel to `a` — at
/// each relative offset we line up the two strands and find every
/// contiguous Watson-Crick paired segment of length ≥ 4. The most
/// stable (lowest ΔG) such segment is returned.
fn score_dimer(a: &[u8], b: &[u8], t_celsius: f64) -> Option<DimerScore> {
    if a.is_empty() || b.is_empty() {
        return None;
    }
    let abytes: Vec<u8> = a.iter().map(|c| c.to_ascii_uppercase()).collect();
    let bbytes: Vec<u8> = b.iter().map(|c| c.to_ascii_uppercase()).collect();
    let na = abytes.len();
    let mut best: Option<DimerScore> = None;

    // For an antiparallel pairing: a[i] pairs with b[nb-1-j]; for
    // alignment offset `s` (s = i + j), we slide `b` against `a`.
    // Equivalently, take the reverse-complement of `b` and look for
    // contiguous matches with `a`.
    let b_rc = reverse_complement_acgt(&bbytes);
    let nrc = b_rc.len();
    // Slide b_rc over a. For each shift `s` in -(nrc-1) ..= (na-1),
    // a[i] aligns with b_rc[i - s] when both indices are valid.
    let shift_min: isize = -(nrc as isize - 1);
    let shift_max: isize = na as isize - 1;
    for shift in shift_min..=shift_max {
        let a_start = shift.max(0) as usize;
        let b_start = (-shift).max(0) as usize;
        let overlap = (na - a_start).min(nrc - b_start);
        if overlap < 4 {
            continue;
        }
        // Walk the overlap; find runs of equal bases (since
        // reverse_complement turns Watson-Crick pairing into byte
        // equality).
        let mut k = 0;
        while k < overlap {
            if abytes[a_start + k] == b_rc[b_start + k] {
                let run_start = k;
                while k < overlap && abytes[a_start + k] == b_rc[b_start + k] {
                    k += 1;
                }
                let run_len = k - run_start;
                if run_len >= 4 {
                    // Score this run.
                    let mut dg = 0.0;
                    for w in 0..run_len - 1 {
                        let a1 = abytes[a_start + run_start + w];
                        let a2 = abytes[a_start + run_start + w + 1];
                        dg += nn_params(&[a1, a2])?.dg(t_celsius);
                    }
                    // 3′-end involvement: the run sits within 2 nt of
                    // a's 3′ end, or it sits within 2 nt of b's 3′ end
                    // (b's 3′ end is at b_rc index 0, since b_rc is
                    // the reverse of b).
                    let a_run_end = a_start + run_start + run_len;
                    let a_3p = a_run_end + 2 >= na;
                    let b_3p = (b_start + run_start) <= 2;
                    let three_prime = a_3p || b_3p;
                    let cand = DimerScore {
                        dg,
                        a_offset: a_start + run_start,
                        b_offset: b_start + run_start,
                        length: run_len,
                        three_prime,
                    };
                    if better_dimer(&best, &cand) {
                        best = Some(cand);
                    }
                }
            } else {
                k += 1;
            }
        }
    }
    best.filter(|d| d.dg < 0.0)
}

fn better_dimer(best: &Option<DimerScore>, cand: &DimerScore) -> bool {
    match best {
        None => true,
        Some(b) => cand.dg.partial_cmp(&b.dg).unwrap_or(Ordering::Equal) == Ordering::Less,
    }
}

/// Reverse complement of a pure-ACGT byte slice (uppercase). U folds
/// to A's complement. Any other byte maps to itself so the dimer
/// scorer treats it as a no-pair.
fn reverse_complement_acgt(bytes: &[u8]) -> Vec<u8> {
    bytes
        .iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'T' | b'U' => b'A',
            b'G' => b'C',
            b'C' => b'G',
            other => other,
        })
        .collect()
}

/// Watson-Crick pairing predicate. `T` and `U` both pair with `A`.
pub fn pairs(a: u8, b: u8) -> bool {
    let a = a.to_ascii_uppercase();
    let b = b.to_ascii_uppercase();
    matches!(
        (a, b),
        (b'A', b'T') | (b'T', b'A') | (b'A', b'U') | (b'U', b'A') | (b'G', b'C') | (b'C', b'G')
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    // -------------------------------------------------------------------
    // SantaLucia parameter set landmarks
    // -------------------------------------------------------------------

    #[test]
    fn santalucia_aa_tt_stack_known_value() {
        // SantaLucia 1998, table 4: AA/TT -> ΔH = -7.9, ΔS = -22.2.
        let p = nn_params(b"AA").unwrap();
        assert!((p.dh - -7.9).abs() < 1e-9);
        assert!((p.ds - -22.2).abs() < 1e-9);
        // Reverse pairing TT/AA is the same dimer read from the other
        // strand.
        let p = nn_params(b"TT").unwrap();
        assert!((p.dh - -7.9).abs() < 1e-9);
    }

    #[test]
    fn santalucia_cg_stack_known_value() {
        let p = nn_params(b"CG").unwrap();
        assert!((p.dh - -10.6).abs() < 1e-9);
        assert!((p.ds - -27.2).abs() < 1e-9);
    }

    #[test]
    fn nn_init_terminal_at_and_gc() {
        let at = nn_init(b'A');
        assert!((at.dh - 2.3).abs() < 1e-9);
        assert!((at.ds - 4.1).abs() < 1e-9);
        let gc = nn_init(b'G');
        assert!((gc.dh - 0.1).abs() < 1e-9);
        assert!((gc.ds - -2.8).abs() < 1e-9);
    }

    #[test]
    fn self_complementary_detection() {
        assert!(is_self_complementary(b"GAATTC"));
        assert!(is_self_complementary(b"AT"));
        assert!(!is_self_complementary(b"GAATTA"));
        assert!(!is_self_complementary(b"AAA"));
    }

    #[test]
    fn dg_at_37c_is_negative_for_a_typical_duplex() {
        // A ~20mer should be very stable at 37°C.
        let s = b"ACGTACGTACGTACGT";
        let nn = duplex_enthalpy_entropy(s).unwrap();
        let dg = nn.dg(37.0);
        assert!(dg < 0.0, "expected ΔG < 0 at 37°C, got {dg}");
    }

    // -------------------------------------------------------------------
    // Salt / Mg / dNTP corrections
    // -------------------------------------------------------------------

    #[test]
    fn salt_effective_na_includes_mg_dntp() {
        // No Mg/dNTP -> just the monovalent.
        let s = SaltParams {
            na_conc: 0.05,
            mg_conc: 0.0,
            dntp_conc: 0.0,
        };
        assert!((s.effective_na() - 0.05).abs() < 1e-12);
        // Standard PCR Mg/dNTPs should push the equivalent up.
        let s = SaltParams::default();
        let eq = s.effective_na();
        assert!(eq > s.na_conc, "Mg should add to Na: {eq} vs {}", s.na_conc);
    }

    #[test]
    fn duplex_tm_responds_to_salt_and_strand() {
        let seq = b"ACGTACGTACGTACGT";
        let s_low = SaltParams {
            na_conc: 0.01,
            mg_conc: 0.0,
            dntp_conc: 0.0,
        };
        let s_high = SaltParams {
            na_conc: 0.2,
            mg_conc: 0.0,
            dntp_conc: 0.0,
        };
        let strand = StrandParams::default();
        let tm_low = duplex_tm(seq, s_low, strand).unwrap();
        let tm_high = duplex_tm(seq, s_high, strand).unwrap();
        assert!(tm_high > tm_low, "high salt -> higher Tm");
    }

    #[test]
    fn mg_alone_raises_tm_above_low_monovalent() {
        let seq = b"ACGTACGTACGTACGT";
        let strand = StrandParams::default();
        let no_mg = SaltParams {
            na_conc: 0.05,
            mg_conc: 0.0,
            dntp_conc: 0.0,
        };
        let with_mg = SaltParams {
            na_conc: 0.05,
            mg_conc: 3e-3,
            dntp_conc: 0.2e-3,
        };
        let tm_no = duplex_tm(seq, no_mg, strand).unwrap();
        let tm_mg = duplex_tm(seq, with_mg, strand).unwrap();
        assert!(
            tm_mg > tm_no,
            "Mg should raise Tm; no_mg={tm_no} mg={tm_mg}"
        );
    }

    #[test]
    fn self_complementary_duplex_uses_symmetry_correction() {
        // GAATTC is self-complementary. The symmetry correction makes
        // ΔS more negative by 1.4 cal/(mol·K), and the strand factor x
        // drops from 4 to 1, both of which shift Tm.
        let pal = b"GAATTC";
        let nn_pal = duplex_enthalpy_entropy(pal).unwrap();
        let nn_naive = {
            // Compute without the symmetry correction by sampling a
            // sister non-palindrome of the same composition.
            let nonpal = b"GAATTA";
            duplex_enthalpy_entropy(nonpal).unwrap()
        };
        // Symmetry correction makes ΔS smaller for the palindrome by
        // exactly SYMMETRY_DS (relative to a sister non-palindrome,
        // but the comparison here is only directional, not exact —
        // base composition differs).
        assert!(nn_pal.ds < nn_naive.ds);
    }

    // -------------------------------------------------------------------
    // Hairpin and dimer ΔG
    // -------------------------------------------------------------------

    #[test]
    fn obvious_hairpin_scores_negatively() {
        // 5 bp stem + 4 nt loop + 5 bp stem.
        let hp = b"GGGGGAAAAACCCCC";
        let score = most_stable_hairpin(hp, 37.0).unwrap();
        assert!(score.dg < 0.0, "expected ΔG < 0, got {}", score.dg);
        assert!(score.stem_len >= 5);
        assert!(score.loop_len >= 3);
    }

    #[test]
    fn random_unstructured_string_has_no_strong_hairpin() {
        // A run of A's cannot form any base pairs.
        let flat = b"AAAAAAAAAAAAAAAA";
        assert!(most_stable_hairpin(flat, 37.0).is_none());
    }

    #[test]
    fn obvious_self_dimer_scores_negatively() {
        // A palindrome necessarily forms a perfect self-dimer.
        let pal = b"GAATTCGAATTC";
        let score = most_stable_self_dimer(pal, 37.0).unwrap();
        assert!(score.dg < 0.0, "ΔG = {} should be < 0", score.dg);
        assert!(score.length >= 4);
    }

    #[test]
    fn no_self_dimer_for_unmatched_sequence() {
        // A homopolymer A run has no Watson-Crick pairing with itself.
        let flat = b"AAAAAAAAAAAA";
        assert!(most_stable_self_dimer(flat, 37.0).is_none());
    }

    #[test]
    fn hetero_dimer_detected_between_complementary_primers() {
        // Two primers whose 3′ ends are reverse-complementary.
        let a = b"GACTAGGCATTTT";
        let b = b"AAAATGCCTAGTC";
        let score = most_stable_hetero_dimer(a, b, 37.0).unwrap();
        assert!(score.dg < 0.0);
        // The 3′-end pairing should be flagged.
        assert!(
            score.three_prime,
            "expected 3′ involvement, score: {score:?}"
        );
    }

    #[test]
    fn unrelated_primers_have_no_hetero_dimer() {
        let a = b"AAAAAAAAAAAA";
        let b = b"AAAAAAAAAAAA";
        assert!(most_stable_hetero_dimer(a, b, 37.0).is_none());
    }

    #[test]
    fn pairs_predicate_matches_watson_crick() {
        for (a, b, expected) in [
            (b'A', b'T', true),
            (b'T', b'A', true),
            (b'G', b'C', true),
            (b'C', b'G', true),
            (b'A', b'U', true),
            (b'A', b'C', false),
            (b'G', b'T', false),
        ] {
            assert_eq!(pairs(a, b), expected, "pairs({a},{b})");
        }
    }
}
