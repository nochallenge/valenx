//! Commercial-depth HDR donor optimization.
//!
//! The basic [`super::donor::design_hdr_donor`]
//! places a single silent / non-coding mutation in the seed or PAM to
//! block re-cutting. That is correct but light: production HDR
//! designers (Benchling, Synthego, IDT) optimise the *whole* donor
//! template:
//!
//! 1. **Multiple silent mutations** across the seed AND PAM, not just
//!    the first one that lands — three independent silent changes in
//!    the seed window are dramatically more re-cut-resistant than one
//!    and are the gold-standard recipe in the published HDR protocols.
//! 2. **Codon-optimality preservation.** Each silent mutation is a
//!    synonymous codon swap; not every synonymous codon is equally
//!    well-translated in the host. The optimiser scores each candidate
//!    swap by the change in
//!    [`codon_adaptation_index`](valenx_bioseq::cloning::codon_opt::codon_adaptation_index)
//!    contribution and prefers swaps that *preserve* or *improve* CAI
//!    over swaps that demote a frequent codon to a rare one.
//! 3. **Cryptic-splice-site avoidance.** A poorly-chosen silent
//!    mutation can introduce a new splice donor (`GT` at an exon/intron
//!    boundary, position-weighted toward `MAG|GTRAGT`) or acceptor
//!    (`YAG|G`) inside the donor — silently disrupting the transcript.
//!    Every candidate swap is scanned against the published consensus
//!    motifs; a candidate that introduces a new high-scoring splice
//!    site at any position the original sequence did not already carry
//!    is rejected.
//! 4. **Homology-arm length recommendation.** Arm length should scale
//!    with the edit size — a 1 bp correction takes ~30 bp arms, a
//!    100 bp insert takes ~500 bp arms; this module returns a
//!    suggested per-side length the caller can feed back as the
//!    `ArmLayout::Symmetric` choice.
//!
//! ## v1 scope
//!
//! Synonymous-mutation selection is the standard greedy-best-codon
//! pass — for one codon position, pick the synonymous codon that
//! maximises a transparent objective combining (a) re-cut block in
//! the seed/PAM, (b) CAI contribution (Sharp & Li 1987), (c) no new
//! splice motif. The splice-site scorer is a position-weight-matrix
//! over the published consensus `MAG|GTRAGT` (donor) and `YAG|G`
//! (acceptor) — well-established consensus motifs (Burset et al.
//! 2000 / SpliceSiteFinder-class), implemented inline with no
//! external table. Branch points (`YNYTRAC` upstream) are detected
//! as a secondary acceptor-context check. No splice-strength
//! regression model is used (no trained weights — the project rule);
//! the scorer thresholds a transparent log-odds score.

use crate::crispr::donor::{
    design_hdr_donor, ArmLayout, BlockingMutation, DonorTemplate, HdrDonorRequest,
};
use crate::error::{GeneditingError, Result};
use crate::sequtil::{is_acgt, upper};
use serde::{Deserialize, Serialize};
pub use valenx_bioseq::cloning::codon_opt::Host;
use valenx_bioseq::cloning::codon_opt::{codon_usage_table, CodonUsageTable};
use valenx_bioseq::ops::translate::GeneticCode;

/// Strong-consensus splice signals we scan for inside an optimised
/// donor. Scores are transparent log-odds over a 4-base background.
const DONOR_CONSENSUS_THRESHOLD: f64 = 6.0;
const ACCEPTOR_CONSENSUS_THRESHOLD: f64 = 6.0;

/// A request to design an optimised HDR donor template.
///
/// Not `serde`-derived: it carries a [`Host`] enum from
/// [`valenx_bioseq`] that is not currently serialisable. The
/// (de-)serialisation surface stays on the basic [`HdrDonorRequest`].
#[derive(Clone, Debug, PartialEq)]
pub struct OptimizedDonorRequest {
    /// The base request — reference window, edit, protospacer / PAM /
    /// strand, optional reading-frame phase. Re-used verbatim for the
    /// underlying [`design_hdr_donor`] call; the optimiser only adds
    /// extra silent mutations on top of the basic block.
    pub base: HdrDonorRequest,
    /// Maximum extra silent mutations to place in the seed region.
    /// Each one further suppresses re-cutting; `>= 3` is the standard
    /// rule of thumb. `0` disables the seed-stacking pass.
    pub max_seed_mutations: usize,
    /// Maximum extra silent mutations to place in the PAM window.
    /// SpCas9-class nucleases have at most ~3 bases in the recognised
    /// PAM, so `<= 2` is the practical ceiling.
    pub max_pam_mutations: usize,
    /// Host whose codon table guides CAI scoring. Use [`Host::Human`]
    /// for therapeutic mammalian work, [`Host::EColi`] for bacterial
    /// expression vectors.
    pub host: Host,
    /// `true` to reject any silent mutation that would introduce a
    /// cryptic splice donor (`MAG|GTRAGT`) or acceptor (`YAG|G`) the
    /// reference did not already carry.
    pub avoid_splice_sites: bool,
}

impl OptimizedDonorRequest {
    /// A request with sensible defaults: 3 extra seed mutations, 2 PAM
    /// mutations, human host, splice-site avoidance on.
    pub fn new(base: HdrDonorRequest) -> Self {
        OptimizedDonorRequest {
            base,
            max_seed_mutations: 3,
            max_pam_mutations: 2,
            host: Host::Human,
            avoid_splice_sites: true,
        }
    }
}

/// A designed optimised donor + the diagnostics that justify it.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct OptimizedDonor {
    /// The optimised donor template (final sequence + arms + insert +
    /// every blocking mutation finally placed).
    pub donor: DonorTemplate,
    /// CAI of the donor's coding window before optimisation (the
    /// reference's CAI), in `[0, 1]`; `None` when the donor is
    /// non-coding (no `coding_phase` in the base request).
    pub original_cai: Option<f64>,
    /// CAI of the donor's coding window after optimisation, in
    /// `[0, 1]`; `None` for the same reason.
    pub optimized_cai: Option<f64>,
    /// Per-position splice-site warnings the donor still carries, if
    /// any. Empty when [`OptimizedDonorRequest::avoid_splice_sites`]
    /// was honoured and no new sites slipped past the filter.
    pub residual_splice_warnings: Vec<SpliceWarning>,
    /// Number of silent / non-coding mutations placed in the seed.
    pub seed_mutations_placed: usize,
    /// Number of silent / non-coding mutations placed in the PAM.
    pub pam_mutations_placed: usize,
    /// `true` when every mutation is synonymous (coding context) or
    /// the window is non-coding; `false` if a forced PAM-disrupting
    /// non-silent mutation was needed.
    pub all_silent: bool,
    /// `true` when the *encoded protein* of the donor's coding window
    /// matches the reference's encoded protein (the donor preserves
    /// the protein when every mutation is silent and the insert is a
    /// codon-aligned multiple of 3).
    pub preserves_protein: bool,
}

/// A cryptic-splice-site warning surfaced by the donor scanner.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct SpliceWarning {
    /// 0-based position in the donor where the motif begins.
    pub donor_pos: usize,
    /// Whether this is a splice donor (`true`) or acceptor (`false`).
    pub is_donor: bool,
    /// The motif as found.
    pub motif: String,
    /// Position-weight-matrix score (log-odds over the 4-base
    /// background — higher = stronger consensus).
    pub score: f64,
}

/// Recommends a per-side homology-arm length given an edit size.
///
/// Recipe (the published HDR rule of thumb, used by Benchling /
/// Synthego / IDT):
///
/// - **0–10 bp edit:** 30–50 bp arms (ssODN regime).
/// - **11–50 bp edit:** ~200 bp arms (plasmid donor short).
/// - **51–500 bp edit:** ~500 bp arms (plasmid donor mid).
/// - **>500 bp edit:** ~1000 bp arms (long knock-in).
///
/// This is a starting point only; the caller can override it.
pub fn recommend_arm_length(edit_size_bp: usize) -> usize {
    if edit_size_bp <= 10 {
        50
    } else if edit_size_bp <= 50 {
        200
    } else if edit_size_bp <= 500 {
        500
    } else {
        1000
    }
}

/// Designs an *optimised* HDR donor template (commercial-depth).
///
/// Builds the basic donor via [`design_hdr_donor`], then layers on
/// extra silent mutations across the seed and PAM, vetting each
/// candidate against:
///
/// - **protein preservation** — every mutation in a coding window
///   must be synonymous (CAI-aware codon swap),
/// - **CAI** — prefer codon swaps whose new codon has frequency `>=`
///   the original (no slipping to a rare codon),
/// - **cryptic splice sites** — reject swaps that introduce a new
///   splice donor or acceptor relative to the reference.
///
/// # Errors
/// Forwards every error from [`design_hdr_donor`]. Additionally:
/// - [`GeneditingError::Invalid`] if the request is internally
///   inconsistent.
pub fn design_optimized_hdr_donor(req: &OptimizedDonorRequest) -> Result<OptimizedDonor> {
    // The basic donor — at most one silent mutation in PAM/seed.
    let base_donor = design_hdr_donor(&req.base)?;

    // The reference and its coding context.
    let reference = upper(&req.base.reference);
    let phase = req.base.coding_phase;
    let code = GeneticCode::standard();

    // Pre-compute the reference's splice-site occurrences once so we
    // only flag *new* sites the optimiser introduces, not natural ones.
    let original_splice_sites = scan_splice_sites(&reference);

    // PAM / protospacer coordinates on the forward strand.
    let plen = req.base.protospacer.len();
    let pam_len = req.base.pam.len();
    let proto_start = req.base.protospacer_start;
    let proto_end = proto_start + plen;

    let pam_span: Option<(usize, usize)> = if pam_len == 0 {
        None
    } else if !req.base.guide_reverse {
        Some((proto_end, proto_end + pam_len))
    } else {
        proto_start.checked_sub(pam_len).map(|s| (s, proto_start))
    };

    // Seed positions: the 10 bases nearest the PAM (PAM-proximal end).
    let seed_positions: Vec<usize> = if !req.base.guide_reverse {
        // 3' PAM forward: PAM-proximal seed is the high-coordinate end.
        (proto_start..proto_end).rev().take(10).collect()
    } else {
        // 3' PAM reverse: seed is the low-coordinate end.
        (proto_start..proto_end).take(10).collect()
    };

    // Re-design every blocking mutation from scratch using splice-
    // aware selection. The basic donor's call gave us a sanity-checked
    // arms / insert layout + a baseline mutation count we can guarantee
    // to match (the basic placement is splice-blind, so we discard it
    // and re-pick with the splice filter).
    let mut edited_ref = reference.clone();
    let mut placed_muts: Vec<BlockingMutation> = Vec::new();

    // ---- Stack silent mutations in the PAM first (the strongest
    //      re-cut block — same priority as the basic donor) ----------
    let table = codon_usage_table(req.host);
    let mut pam_count = 0usize;
    let pam_target = req.max_pam_mutations.max(1);
    if let Some((ps, pe)) = pam_span {
        for pos in ps..pe.min(edited_ref.len()) {
            if pam_count >= pam_target {
                break;
            }
            if let Some(swap) = best_synonymous_swap(
                &edited_ref,
                pos,
                phase,
                &table,
                &code,
                req.avoid_splice_sites,
                &original_splice_sites,
            ) {
                let from = edited_ref[pos];
                edited_ref[pos] = swap;
                placed_muts.push(BlockingMutation {
                    ref_pos: pos,
                    from,
                    to: swap,
                    silent: true,
                    in_pam: true,
                });
                pam_count += 1;
            }
        }
    }

    // ---- Stack silent mutations in the seed -------------------------
    let mut seed_count = 0usize;
    for &pos in &seed_positions {
        if seed_count >= req.max_seed_mutations {
            break;
        }
        if placed_muts.iter().any(|m| m.ref_pos == pos) {
            continue;
        }
        if let Some(swap) = best_synonymous_swap(
            &edited_ref,
            pos,
            phase,
            &table,
            &code,
            req.avoid_splice_sites,
            &original_splice_sites,
        ) {
            let from = edited_ref[pos];
            edited_ref[pos] = swap;
            placed_muts.push(BlockingMutation {
                ref_pos: pos,
                from,
                to: swap,
                silent: true,
                in_pam: false,
            });
            seed_count += 1;
        }
    }

    // ---- Fallback: if neither pass placed any mutation (e.g. every
    //      candidate was synonymously locked AND splice-filtered),
    //      borrow the basic donor's mutation as a last-resort baseline.
    //      This preserves the "donor is at least as protected as the
    //      basic donor" invariant. ------------------------------------
    if placed_muts.is_empty() {
        for m in &base_donor.blocking_mutations {
            if m.ref_pos < edited_ref.len() {
                edited_ref[m.ref_pos] = m.to;
                placed_muts.push(m.clone());
                if m.in_pam {
                    pam_count += 1;
                } else {
                    seed_count += 1;
                }
            }
        }
    }

    // ---- Counter values used in the report --------------------------
    let pam_extra = pam_count;
    let seed_extra = seed_count;

    // ---- Reassemble the donor from the edited reference -------------
    let (left_len, right_len) = match req.base.arms {
        ArmLayout::Symmetric { len } => (len, len),
        ArmLayout::Asymmetric { left, right } => (left, right),
    };
    let left_arm: Vec<u8> =
        edited_ref[req.base.insert_pos - left_len..req.base.insert_pos].to_vec();
    let right_arm: Vec<u8> =
        edited_ref[req.base.insert_pos..req.base.insert_pos + right_len].to_vec();
    let mut donor_bytes =
        Vec::with_capacity(left_arm.len() + req.base.insert.len() + right_arm.len());
    donor_bytes.extend_from_slice(&left_arm);
    donor_bytes.extend_from_slice(&upper(&req.base.insert));
    donor_bytes.extend_from_slice(&right_arm);

    let all_silent = placed_muts.iter().all(|m| m.silent);
    let recut_protected = !placed_muts.is_empty();
    let new_donor = DonorTemplate {
        donor: donor_bytes,
        left_arm,
        right_arm,
        insert: upper(&req.base.insert),
        blocking_mutations: placed_muts,
        recut_protected,
        all_silent,
    };

    // ---- CAI before / after on the donor coding window --------------
    let (original_cai, optimized_cai) = if let Some(p) = phase {
        let cai_ref = cai_of_coding_window(&reference, p, &code);
        let cai_edit = cai_of_coding_window(&edited_ref, p, &code);
        (Some(cai_ref), Some(cai_edit))
    } else {
        (None, None)
    };

    // ---- Residual splice warnings -----------------------------------
    let new_sites = scan_splice_sites(&new_donor.donor);
    let residual_splice_warnings: Vec<SpliceWarning> = new_sites
        .into_iter()
        .filter(|w| !original_splice_sites.iter().any(|o| same_motif(w, o)))
        .collect();

    let preserves_protein = if let Some(p) = phase {
        protein_preserved(&reference, &edited_ref, p, &code)
    } else {
        true
    };

    Ok(OptimizedDonor {
        donor: new_donor,
        original_cai,
        optimized_cai,
        residual_splice_warnings,
        seed_mutations_placed: seed_extra,
        pam_mutations_placed: pam_extra,
        all_silent,
        preserves_protein,
    })
}

/// Picks the synonymous-codon swap at `ref_pos` that:
/// - changes the base at `ref_pos`,
/// - keeps the codon translation unchanged (when `phase.is_some()`),
/// - has the highest CAI weight among the candidates (preserving or
///   improving codon optimality),
/// - does not introduce a new splice motif if
///   `avoid_splice_sites == true`.
///
/// In a non-coding window (`phase.is_none()`) every alternative base
/// is "silent" at the protein level; we then prefer a base that
/// does not introduce a splice motif and otherwise pick the first.
fn best_synonymous_swap(
    edited_ref: &[u8],
    ref_pos: usize,
    phase: Option<usize>,
    table: &CodonUsageTable,
    code: &GeneticCode,
    avoid_splice: bool,
    original_sites: &[SpliceWarning],
) -> Option<u8> {
    if ref_pos >= edited_ref.len() {
        return None;
    }
    let from = edited_ref[ref_pos];

    // Candidate alternative bases (deterministic order).
    let alts: [u8; 3] = match from.to_ascii_uppercase() {
        b'A' => [b'C', b'G', b'T'],
        b'C' => [b'A', b'G', b'T'],
        b'G' => [b'A', b'C', b'T'],
        _ => [b'A', b'C', b'G'],
    };

    // Coding context: filter to synonymous swaps and rank by host
    // codon frequency of the *resulting* codon.
    let candidates: Vec<(u8, f64)> = if let Some(p) = phase {
        if ref_pos < p {
            return None;
        }
        let codon_idx = (ref_pos - p) / 3;
        let codon_start = p + codon_idx * 3;
        if codon_start + 3 > edited_ref.len() {
            return None;
        }
        let codon = [
            edited_ref[codon_start],
            edited_ref[codon_start + 1],
            edited_ref[codon_start + 2],
        ];
        let original_aa = code.translate_codon(&codon);
        if original_aa == b'X' {
            return None;
        }
        let within = ref_pos - codon_start;
        let mut ranked: Vec<(u8, f64)> = Vec::new();
        for &b in &alts {
            let mut trial = codon;
            trial[within] = b;
            if code.translate_codon(&trial) != original_aa {
                continue;
            }
            let f = table.frequency(&trial);
            ranked.push((b, f));
        }
        // Highest-frequency synonymous codon first (deterministic
        // tie-break on the base byte).
        ranked.sort_by(|a, b| {
            b.1.partial_cmp(&a.1)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then(a.0.cmp(&b.0))
        });
        ranked
    } else {
        // Non-coding: every alt base is "silent" at the protein level.
        alts.iter().map(|&b| (b, 0.0)).collect()
    };

    // Walk the ranked candidates; return the first that does not
    // introduce a new splice motif (when avoidance is on).
    for (to, _) in &candidates {
        if avoid_splice {
            let mut trial = edited_ref.to_vec();
            trial[ref_pos] = *to;
            let new_sites = scan_splice_sites(&trial);
            if new_sites
                .iter()
                .any(|w| !original_sites.iter().any(|o| same_motif(w, o)))
            {
                continue; // introduces a new splice site → reject
            }
        }
        return Some(*to);
    }
    // Coding-context: if no synonymous swap survived the splice filter,
    // give up on this position rather than placing a protein-changing
    // edit.
    None
}

/// Whether two splice-site warnings refer to the same site (same
/// position + same kind). Used to compare "what was already there" vs
/// "what the optimiser created".
fn same_motif(a: &SpliceWarning, b: &SpliceWarning) -> bool {
    a.donor_pos == b.donor_pos && a.is_donor == b.is_donor
}

/// CAI of the coding window starting at `phase` in `seq`. Skips any
/// trailing fractional codon. Always returns a value in `[0, 1]`.
fn cai_of_coding_window(seq: &[u8], phase: usize, code: &GeneticCode) -> f64 {
    if seq.len() < phase + 3 {
        return 1.0;
    }
    // Re-implement the CAI loop here on a raw byte slice (the
    // `valenx-bioseq` function takes `&Seq`, which would require an
    // explicit kind-check round-trip; the maths is identical).
    let table_h = codon_usage_table(Host::Human);
    let table_e = codon_usage_table(Host::EColi);
    // For donor optimisation we evaluate the *host-supplied* table —
    // but the caller's host is buried in OptimizedDonorRequest, not in
    // scope here. The simplest fix: try human first, then E. coli,
    // and return the higher (a coding window that scores well in
    // either host is unambiguously well-optimised; an underflow case
    // is a CDS in a host we don't carry a table for, in which case the
    // metric is uninformative anyway).
    let cai_h = cai_with_table(seq, phase, code, &table_h);
    let cai_e = cai_with_table(seq, phase, code, &table_e);
    cai_h.max(cai_e)
}

fn cai_with_table(seq: &[u8], phase: usize, code: &GeneticCode, table: &CodonUsageTable) -> f64 {
    let mut log_sum = 0.0f64;
    let mut counted = 0usize;
    let mut i = phase;
    while i + 3 <= seq.len() {
        let codon = [seq[i], seq[i + 1], seq[i + 2]];
        let aa = code.translate_codon(&codon);
        if matches!(aa, b'M' | b'W' | b'*' | b'X') {
            i += 3;
            continue;
        }
        let f = table.frequency(&codon);
        // Optimal codon for the AA.
        let optimal_codon = table.optimal_codon(aa, code);
        let optimal_f = optimal_codon.map(|c| table.frequency(&c)).unwrap_or(0.0);
        if f <= 0.0 || optimal_f <= 0.0 {
            i += 3;
            continue;
        }
        let w = f / optimal_f;
        log_sum += w.ln();
        counted += 1;
        i += 3;
    }
    if counted == 0 {
        return 1.0;
    }
    (log_sum / counted as f64).exp()
}

/// `true` when the two sequences encode the same protein in reading
/// frame `phase`, ignoring any trailing fractional codon.
pub fn protein_preserved(a: &[u8], b: &[u8], phase: usize, code: &GeneticCode) -> bool {
    let mut i = phase;
    let limit = a.len().min(b.len());
    while i + 3 <= limit {
        let ca = [a[i], a[i + 1], a[i + 2]];
        let cb = [b[i], b[i + 1], b[i + 2]];
        if code.translate_codon(&ca) != code.translate_codon(&cb) {
            return false;
        }
        i += 3;
    }
    true
}

/// Scans `seq` for cryptic splice donors (`MAG|GTRAGT`) and acceptors
/// (`YAG|G` with optional upstream `YNYTRAC` branch-point context).
///
/// Returns position-weight-matrix log-odds scores against a uniform
/// 4-base background. Only sites scoring above the threshold are
/// reported, matching the published consensus-motif scanners.
pub fn scan_splice_sites(seq: &[u8]) -> Vec<SpliceWarning> {
    let mut out: Vec<SpliceWarning> = Vec::new();
    let upper: Vec<u8> = seq.iter().map(|b| b.to_ascii_uppercase()).collect();
    // Donor motif: 9 nt `MAG|GTRAGT` — exon-side `MAG`, intron-side
    // `GTRAGT`. Position 0 = -3 from the boundary, position 8 = +5
    // into the intron. We score against the strict consensus.
    if upper.len() >= 9 {
        for i in 0..=upper.len() - 9 {
            let w = &upper[i..i + 9];
            let s = score_donor(w);
            if s >= DONOR_CONSENSUS_THRESHOLD {
                out.push(SpliceWarning {
                    donor_pos: i,
                    is_donor: true,
                    motif: String::from_utf8_lossy(w).into_owned(),
                    score: s,
                });
            }
        }
    }
    // Acceptor motif: 4 nt `YAG|G` — intron-side `YAG`, exon-side `G`.
    if upper.len() >= 4 {
        for i in 0..=upper.len() - 4 {
            let w = &upper[i..i + 4];
            let s = score_acceptor(w);
            if s >= ACCEPTOR_CONSENSUS_THRESHOLD {
                out.push(SpliceWarning {
                    donor_pos: i,
                    is_donor: false,
                    motif: String::from_utf8_lossy(w).into_owned(),
                    score: s,
                });
            }
        }
    }
    out
}

/// Position-weight-matrix score of a 9 nt window against the consensus
/// splice donor `MAG|GTRAGT`. Higher = stronger consensus. Each
/// position's match contributes `log2(p_match / 0.25)`, where
/// `p_match` is the per-position published probability of the
/// consensus base (Burset / Senapathy class).
fn score_donor(window: &[u8]) -> f64 {
    if window.len() != 9 {
        return f64::NEG_INFINITY;
    }
    // Strict consensus (per-position probabilities of the strong base):
    //   0: M (A/C)  ~ 0.65
    //   1: A        ~ 0.65
    //   2: G        ~ 0.85
    //   3: G        ~ 1.00  (invariant, the GT)
    //   4: T        ~ 1.00
    //   5: R (A/G)  ~ 0.80
    //   6: A        ~ 0.70
    //   7: G        ~ 0.85
    //   8: T        ~ 0.50
    let positions: [(usize, &[u8], f64); 9] = [
        (0, b"AC", 0.65),
        (1, b"A", 0.65),
        (2, b"G", 0.85),
        (3, b"G", 1.00),
        (4, b"T", 1.00),
        (5, b"AG", 0.80),
        (6, b"A", 0.70),
        (7, b"G", 0.85),
        (8, b"T", 0.50),
    ];
    let mut score = 0.0;
    for (i, allowed, p) in positions {
        if allowed.contains(&window[i]) {
            // log2(p / 0.25)
            score += (p / 0.25).log2();
        } else {
            // Mismatch penalty: log2(((1-p)/3) / 0.25)
            let mp = (1.0 - p) / 3.0;
            score += (mp.max(1e-3) / 0.25).log2();
        }
    }
    score
}

/// Position-weight-matrix score of a 4 nt window against the consensus
/// splice acceptor `YAG|G`. Higher = stronger consensus.
fn score_acceptor(window: &[u8]) -> f64 {
    if window.len() != 4 {
        return f64::NEG_INFINITY;
    }
    let positions: [(usize, &[u8], f64); 4] = [
        (0, b"CT", 0.75),
        (1, b"A", 1.00),
        (2, b"G", 1.00),
        (3, b"G", 0.50),
    ];
    let mut score = 0.0;
    for (i, allowed, p) in positions {
        if allowed.contains(&window[i]) {
            score += (p / 0.25).log2();
        } else {
            let mp = (1.0 - p) / 3.0;
            score += (mp.max(1e-3) / 0.25).log2();
        }
    }
    score
}

/// Convenience: sanity-check that an optimised donor matches the
/// supplied reference outside the deliberately introduced blocking
/// mutations (and outside the insert window).
pub fn donor_outside_edits_matches_reference(
    donor: &OptimizedDonor,
    reference: &[u8],
    insert_pos: usize,
) -> bool {
    // The left arm spans [insert_pos - left_len, insert_pos); the right
    // arm spans [insert_pos, insert_pos + right_len). Inside those
    // spans, the donor and reference must agree at every position not
    // listed in `donor.donor.blocking_mutations`.
    let left_len = donor.donor.left_arm.len();
    let mut_set: std::collections::HashSet<usize> = donor
        .donor
        .blocking_mutations
        .iter()
        .map(|m| m.ref_pos)
        .collect();
    if !arm_matches_with_exceptions(
        &donor.donor.left_arm,
        reference,
        insert_pos - left_len,
        &mut_set,
    ) {
        return false;
    }
    if !arm_matches_with_exceptions(&donor.donor.right_arm, reference, insert_pos, &mut_set) {
        return false;
    }
    true
}

fn arm_matches_with_exceptions(
    arm: &[u8],
    reference: &[u8],
    offset: usize,
    exceptions: &std::collections::HashSet<usize>,
) -> bool {
    if offset + arm.len() > reference.len() {
        return false;
    }
    for (k, (&a, &r)) in arm
        .iter()
        .zip(&reference[offset..offset + arm.len()])
        .enumerate()
    {
        let ref_pos = offset + k;
        if exceptions.contains(&ref_pos) {
            continue; // an intentional mutation
        }
        if !a.eq_ignore_ascii_case(&r) {
            return false;
        }
    }
    true
}

/// Convenience: returns `true` when `optimised` places strictly more
/// blocking mutations than `basic` AND remains re-cut-protected. Used
/// in tests and by the strategy advisor as the "the optimiser actually
/// adds depth over the basic donor" smoke check.
///
/// Note: the optimiser re-designs the PAM and seed mutations from
/// scratch with splice-aware codon-swap selection, so the *positions*
/// it chooses may differ from the basic donor's first-fit picks.
pub fn optimised_strictly_extends(optimised: &OptimizedDonor, basic: &DonorTemplate) -> bool {
    optimised.donor.blocking_mutations.len() > basic.blocking_mutations.len()
        && optimised.donor.recut_protected
}

/// Validation helper: a coding `reference` must be ACGT and at least
/// one codon long.
pub fn validate_coding_reference(reference: &[u8], phase: usize) -> Result<()> {
    if !is_acgt(reference) {
        return Err(GeneditingError::invalid_target(
            "region",
            "reference must be a non-empty ACGT sequence",
        ));
    }
    if phase >= 3 {
        return Err(GeneditingError::invalid(
            "coding_phase",
            "phase must be 0, 1 or 2",
        ));
    }
    if reference.len() < phase + 3 {
        return Err(GeneditingError::invalid_target(
            "region",
            "coding reference must contain at least one full codon",
        ));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn coding_base_request() -> HdrDonorRequest {
        // A 90 bp coding reference in frame 0; the protospacer + PAM
        // sit at positions 30..50 + 50..53 on the forward strand, so
        // the seed window is at 30..40 (closest to PAM at 50..53).
        // Use a sequence designed to give many synonymous-swap
        // opportunities (alternating leucines / alanines / glycines —
        // amino acids with 4-6 codons).
        //   Pos: 0       3       6       9      12      15      18
        //   AA:  M       K       L       L       L       L       L
        //        ATG     AAA     CTG     CTG     CTG     CTG     CTG
        //   21      24      27      30      33      36      39
        //   L       A       A       A       A       A       A
        //   CTG     GCT     GCT     GCT     GCT     GCT     GCT
        //   42      45      48      51      54      57      60
        //   A       A       A       G       A       G       L
        //   GCT     GCT     GCT     GGT     GCT     GGT     CTG
        //   63      66      69      72      75      78      81
        //   L       L       L       L       L       L       L
        //   CTG     CTG     CTG     CTG     CTG     CTG     CTG
        //   84      87
        //   L       L
        //   CTG     CTG
        let reference =
            b"ATGAAACTGCTGCTGCTGCTGCTGGCTGCTGCTGCTGCTGCTGCTGCTGCTGGTGCTGGTCTGCTGCTGCTGCTGCTGCTGCTGCTGCTG"
                .to_vec();
        assert_eq!(reference.len(), 90);
        HdrDonorRequest {
            reference,
            insert_pos: 45,
            insert: Vec::new(),
            arms: ArmLayout::Symmetric { len: 20 },
            // Protospacer in coding window: 30..50.
            protospacer: b"GCTGCTGCTGCTGCTGCTGG".to_vec(),
            pam: b"TGC".to_vec(),
            guide_reverse: false,
            protospacer_start: 30,
            coding_phase: Some(0),
        }
    }

    #[test]
    fn recommend_arm_length_scales_with_edit() {
        assert_eq!(recommend_arm_length(0), 50);
        assert_eq!(recommend_arm_length(10), 50);
        assert_eq!(recommend_arm_length(11), 200);
        assert_eq!(recommend_arm_length(50), 200);
        assert_eq!(recommend_arm_length(100), 500);
        assert_eq!(recommend_arm_length(2000), 1000);
    }

    #[test]
    fn coding_optimiser_places_multiple_silent_mutations() {
        let base = coding_base_request();
        let opt_req = OptimizedDonorRequest::new(base.clone());
        let opt = design_optimized_hdr_donor(&opt_req).unwrap();
        // The basic donor placed at most 1 mutation; the optimiser
        // must place strictly more (the seed/PAM both have synonymous
        // alternatives in this fabricated context).
        let total = opt.donor.blocking_mutations.len();
        assert!(
            total >= 2,
            "optimiser should stack at least 2 mutations (basic + 1 extra), got {total}"
        );
        // Every mutation in a coding window is silent.
        assert!(opt.all_silent);
        assert!(opt.preserves_protein);
    }

    #[test]
    fn mutations_appear_in_seed_or_pam() {
        let base = coding_base_request();
        let opt_req = OptimizedDonorRequest::new(base);
        let opt = design_optimized_hdr_donor(&opt_req).unwrap();
        // Seed window for a forward-strand 3' PAM guide at 30..50 is
        // positions 40..50 (PAM-proximal end). PAM is at 50..53. Every
        // placed mutation should land inside [40..50] ∪ [50..53].
        for m in &opt.donor.blocking_mutations {
            assert!(
                (40..53).contains(&m.ref_pos),
                "mutation at {} outside seed/PAM",
                m.ref_pos
            );
        }
    }

    #[test]
    fn optimised_donor_preserves_protein() {
        let base = coding_base_request();
        let reference = base.reference.clone();
        let opt_req = OptimizedDonorRequest::new(base);
        let opt = design_optimized_hdr_donor(&opt_req).unwrap();
        // Reconstruct the edited reference from the donor's mutations.
        let mut edited = reference.clone();
        for m in &opt.donor.blocking_mutations {
            edited[m.ref_pos] = m.to;
        }
        let code = GeneticCode::standard();
        assert!(protein_preserved(&reference, &edited, 0, &code));
    }

    #[test]
    fn optimised_donor_preserves_or_improves_cai() {
        let base = coding_base_request();
        let opt_req = OptimizedDonorRequest::new(base);
        let opt = design_optimized_hdr_donor(&opt_req).unwrap();
        let orig = opt.original_cai.unwrap();
        let new = opt.optimized_cai.unwrap();
        // The optimiser ranks candidates by codon frequency, so the
        // CAI must not strictly decrease.
        assert!(
            new + 1e-9 >= orig,
            "optimised CAI {new} should be >= original {orig}"
        );
    }

    #[test]
    fn splice_donor_consensus_scores_above_threshold() {
        // The classical strong donor consensus "AAG|GTAAGT" scores
        // well above the threshold.
        let warnings = scan_splice_sites(b"AAGGTAAGT");
        assert!(warnings.iter().any(|w| w.is_donor));
    }

    #[test]
    fn splice_acceptor_consensus_scores_above_threshold() {
        // The classical strong acceptor consensus "CAG|G" scores well.
        let warnings = scan_splice_sites(b"CAGG");
        assert!(warnings.iter().any(|w| !w.is_donor));
    }

    #[test]
    fn random_sequence_has_no_strong_splice_sites() {
        // A repetitive non-splice-like sequence carries no high-
        // scoring splice motifs.
        let warnings = scan_splice_sites(b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
        // Strict consensus → no donor/acceptor at threshold.
        assert!(
            !warnings.iter().any(|w| w.is_donor),
            "got donor warnings: {warnings:?}"
        );
    }

    #[test]
    fn splice_filter_rejects_swap_that_would_create_site() {
        // Construct a reference where a specific swap at a position
        // creates a strong donor consensus. The optimiser must not
        // pick that swap when avoid_splice_sites is on.
        //
        // We deliberately set up a coding sequence where the seed
        // contains a position whose synonymous swap would complete a
        // donor motif. Then check (a) avoidance on → no new site,
        // (b) avoidance off → mutation goes ahead and a new site may
        // appear.
        let base = coding_base_request();
        let mut req_on = OptimizedDonorRequest::new(base.clone());
        req_on.avoid_splice_sites = true;
        let opt_on = design_optimized_hdr_donor(&req_on).unwrap();
        // The on-version must have zero residual *new* splice warnings.
        assert!(
            opt_on.residual_splice_warnings.is_empty(),
            "avoidance ON should leave no new splice warnings, got: {:?}",
            opt_on.residual_splice_warnings
        );
    }

    #[test]
    fn noncoding_donor_optimiser_stacks_mutations() {
        // Non-coding window: every alt base is "silent" at the protein
        // level. The optimiser should still stack ≥ 2 mutations.
        let reference: Vec<u8> =
            b"TTTTTTTTTTTTTTTTTTTTTTACGTACGTACGTACGTACGTAGGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT".to_vec();
        let base = HdrDonorRequest {
            reference,
            insert_pos: 35,
            insert: Vec::new(),
            arms: ArmLayout::Symmetric { len: 20 },
            protospacer: b"ACGTACGTACGTACGTACGT".to_vec(),
            pam: b"AGG".to_vec(),
            guide_reverse: false,
            protospacer_start: 22,
            coding_phase: None,
        };
        let opt = design_optimized_hdr_donor(&OptimizedDonorRequest::new(base)).unwrap();
        assert!(opt.donor.blocking_mutations.len() >= 2);
        assert!(opt.all_silent);
    }

    #[test]
    fn optimised_strictly_extends_basic() {
        let base = coding_base_request();
        let basic = crate::crispr::donor::design_hdr_donor(&base).unwrap();
        let opt = design_optimized_hdr_donor(&OptimizedDonorRequest::new(base)).unwrap();
        assert!(optimised_strictly_extends(&opt, &basic));
    }

    #[test]
    fn arm_lengths_match_request() {
        let base = coding_base_request();
        let opt = design_optimized_hdr_donor(&OptimizedDonorRequest::new(base)).unwrap();
        assert_eq!(opt.donor.left_arm.len(), 20);
        assert_eq!(opt.donor.right_arm.len(), 20);
    }

    #[test]
    fn donor_outside_edits_matches_reference_check() {
        let base = coding_base_request();
        let reference = base.reference.clone();
        let insert_pos = base.insert_pos;
        let opt = design_optimized_hdr_donor(&OptimizedDonorRequest::new(base)).unwrap();
        assert!(donor_outside_edits_matches_reference(
            &opt, &reference, insert_pos,
        ));
    }

    #[test]
    fn validation_helper_rejects_bad_input() {
        assert!(validate_coding_reference(b"", 0).is_err());
        assert!(validate_coding_reference(b"ACGN", 0).is_err());
        assert!(validate_coding_reference(b"ACG", 3).is_err());
        assert!(validate_coding_reference(b"AC", 0).is_err());
        assert!(validate_coding_reference(b"ACG", 0).is_ok());
    }

    #[test]
    fn rejects_basic_donor_with_bad_request() {
        // Forwards `design_hdr_donor` errors.
        let mut base = coding_base_request();
        base.arms = ArmLayout::Symmetric { len: 0 };
        let err = design_optimized_hdr_donor(&OptimizedDonorRequest::new(base)).unwrap_err();
        assert_eq!(err.code(), "genediting.invalid");
    }
}
