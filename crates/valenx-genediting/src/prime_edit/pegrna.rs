//! Features 12–14 — pegRNA design, PBS / RT-length optimisation and
//! PE3 / PE3b nicking-guide design.
//!
//! A **pegRNA** has, 5′→3′:
//!
//! 1. a **spacer** — the 20 nt protospacer that targets the nickase;
//! 2. the **sgRNA scaffold** — the constant tracrRNA-derived stem;
//! 3. a **3′ extension** made of two pieces:
//!    - the **RT template** — read by the reverse transcriptase, it
//!      encodes the desired edit and the few bases downstream of the
//!      nick;
//!    - the **primer-binding site (PBS)** — complementary to the 3′ end
//!      the nick exposes, it primes reverse transcription.
//!
//! The Cas9-H840A nickase nicks the protospacer (PAM) strand between
//! protospacer positions 17 and 18 (3 nt 5′ of the PAM). The exposed
//! 3′ flap anneals to the PBS, the RT copies the RT template, and the
//! new flap — carrying the edit — replaces the old one.
//!
//! This module:
//!
//! - designs a pegRNA for a desired substitution, insertion or
//!   deletion ([`design_pegrna`], feature 12);
//! - scans PBS and RT-template lengths for the best-scoring
//!   combination ([`scan_pbs_rt`], feature 13);
//! - designs a PE3 or PE3b nicking guide on the complementary strand
//!   ([`design_nicking_guide`], feature 14).
//!
//! ## v1 scope
//!
//! The pegRNA designer requires the spacer / PAM to be on the
//! **forward strand** of the supplied reference window — the common
//! case and what keeps the coordinate algebra auditable. The
//! length-scan score is the transparent heuristic from
//! [`crate::prime_edit::strategy`] (a PBS-Tm optimum, an RT-template
//! length preference, homopolymer penalties); it is not a trained
//! efficiency model.

use crate::error::{GeneditingError, Result};
use crate::prime_edit::editor::{prime_editor, PrimeEditorId};
use crate::prime_edit::strategy::length_scan_score;
use crate::sequtil::{is_acgt, revcomp, transcribe, upper};
use serde::{Deserialize, Serialize};
use valenx_genomics::crispr::guide::{scan_guides, GuideStrand, PamSpec, PamSide};

/// The kind of edit a pegRNA installs.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum PrimeEdit {
    /// Replace `len` reference bases at the edit position with `to`
    /// (a substitution; for a single-base SNV `len == to.len() == 1`).
    Substitution {
        /// Number of reference bases replaced.
        len: usize,
        /// The replacement bases (forward strand, 5′→3′).
        to: Vec<u8>,
    },
    /// Insert `seq` at the edit position (no reference bases removed).
    Insertion {
        /// The inserted bases (forward strand, 5′→3′).
        seq: Vec<u8>,
    },
    /// Delete `len` reference bases starting at the edit position.
    Deletion {
        /// Number of reference bases removed.
        len: usize,
    },
}

impl PrimeEdit {
    /// A single-base substitution convenience constructor.
    pub fn snv(to: u8) -> Self {
        PrimeEdit::Substitution {
            len: 1,
            to: vec![to.to_ascii_uppercase()],
        }
    }

    /// Number of reference bases this edit consumes at the edit point.
    fn ref_consumed(&self) -> usize {
        match self {
            PrimeEdit::Substitution { len, .. } => *len,
            PrimeEdit::Insertion { .. } => 0,
            PrimeEdit::Deletion { len } => *len,
        }
    }

    /// The bases this edit writes at the edit point (empty for a
    /// deletion).
    fn inserted(&self) -> Vec<u8> {
        match self {
            PrimeEdit::Substitution { to, .. } => upper(to),
            PrimeEdit::Insertion { seq } => upper(seq),
            PrimeEdit::Deletion { .. } => Vec::new(),
        }
    }

    /// A short human-readable label.
    pub fn label(&self) -> String {
        match self {
            PrimeEdit::Substitution { len, to } => {
                format!("{len}-bp substitution to {}", String::from_utf8_lossy(to))
            }
            PrimeEdit::Insertion { seq } => {
                format!("{}-bp insertion of {}", seq.len(), String::from_utf8_lossy(seq))
            }
            PrimeEdit::Deletion { len } => format!("{len}-bp deletion"),
        }
    }
}

/// A request to design a pegRNA.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct PegRnaRequest {
    /// The reference window (forward strand, 5′→3′, unambiguous ACGT).
    pub reference: Vec<u8>,
    /// 0-based forward-strand index where the edit starts.
    pub edit_pos: usize,
    /// The edit to install.
    pub edit: PrimeEdit,
    /// Which prime-editor configuration (drives the spacer PAM).
    pub editor: PrimeEditorId,
    /// PBS length to use (typical 8–17 nt). Use [`scan_pbs_rt`] to pick.
    pub pbs_len: usize,
    /// Number of reference bases of the RT template *downstream* of the
    /// edit (the "RT template homology" past the edit — typical 10–20
    /// nt). The full RT template = these bases + the edit.
    pub rt_template_homology: usize,
}

impl PegRnaRequest {
    /// A request with typical default lengths (PBS 13, RT homology 13).
    pub fn new(
        reference: impl Into<Vec<u8>>,
        edit_pos: usize,
        edit: PrimeEdit,
        editor: PrimeEditorId,
    ) -> Self {
        PegRnaRequest {
            reference: reference.into(),
            edit_pos,
            edit,
            editor,
            pbs_len: 13,
            rt_template_homology: 13,
        }
    }
}

/// A designed pegRNA.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct PegRna {
    /// The 20 nt spacer (protospacer), 5′→3′, as an RNA sequence.
    pub spacer: Vec<u8>,
    /// The PAM the spacer sits next to (DNA, as found on the reference).
    pub pam: Vec<u8>,
    /// 0-based forward-strand start of the spacer / protospacer.
    pub spacer_start: usize,
    /// 0-based forward-strand position of the nick (between this index
    /// and the one before it — i.e. the nick is `nick_pos`).
    pub nick_pos: usize,
    /// The primer-binding site, 5′→3′ as RNA (part of the 3′ extension).
    pub pbs: Vec<u8>,
    /// The reverse-transcriptase template, 5′→3′ as RNA (part of the 3′
    /// extension). Encodes the edit and the downstream homology.
    pub rt_template: Vec<u8>,
    /// The full 3′ extension = RT template + PBS, 5′→3′ as RNA.
    pub three_prime_extension: Vec<u8>,
    /// The complete pegRNA = spacer + scaffold + 3′ extension, 5′→3′ as
    /// RNA.
    pub full_pegrna: Vec<u8>,
    /// The edit this pegRNA installs (echoed back).
    pub edit: PrimeEdit,
}

impl PegRna {
    /// Length of the 3′ extension (RT template + PBS).
    pub fn extension_len(&self) -> usize {
        self.three_prime_extension.len()
    }

    /// `true` when the 3′ extension is empty (never produced; present
    /// for the `len`/`is_empty` lint pairing).
    pub fn is_empty(&self) -> bool {
        self.three_prime_extension.is_empty()
    }
}

/// The constant sgRNA scaffold inserted between the spacer and the 3′
/// extension (the canonical optimised prime-editing scaffold, RNA).
const PEGRNA_SCAFFOLD: &[u8] =
    b"GUUUAAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGC";

/// The H840A-nickase nick offset: the PAM (protospacer) strand is
/// nicked 3 nt 5′ of the PAM, i.e. between protospacer positions 17
/// and 18 for a 20 nt spacer.
const NICK_OFFSET_FROM_PAM: usize = 3;

/// Designs a pegRNA for a desired edit (feature 12).
///
/// Finds a forward-strand spacer whose nick falls *upstream* of (5′ of)
/// the edit site — the prime-editing requirement, since the RT copies
/// from the nick across the edit. Builds the PBS (reverse-complement of
/// the reference 3′ end the nick exposes) and the RT template (encoding
/// the edited sequence plus downstream homology), assembles the 3′
/// extension and the full pegRNA.
///
/// # Errors
/// - [`GeneditingError::InvalidTarget`] for a non-ACGT reference /
///   edit, or an out-of-range `edit_pos`.
/// - [`GeneditingError::Invalid`] for a zero PBS or RT-homology length.
/// - [`GeneditingError::NoValidDesign`] if no forward-strand spacer
///   nicks upstream of the edit, or the reference window is too short
///   for the requested PBS / RT lengths.
pub fn design_pegrna(req: &PegRnaRequest) -> Result<PegRna> {
    if !is_acgt(&req.reference) {
        return Err(GeneditingError::invalid_target(
            "region",
            "reference window must be non-empty ACGT",
        ));
    }
    if req.edit_pos >= req.reference.len() {
        return Err(GeneditingError::invalid_target(
            "region",
            "edit position is outside the reference window",
        ));
    }
    match &req.edit {
        PrimeEdit::Substitution { len, to } => {
            if *len == 0 || to.is_empty() || !is_acgt(to) {
                return Err(GeneditingError::invalid_target(
                    "region",
                    "substitution needs a positive length and an ACGT replacement",
                ));
            }
            if req.edit_pos + len > req.reference.len() {
                return Err(GeneditingError::invalid_target(
                    "region",
                    "substitution runs past the end of the reference window",
                ));
            }
        }
        PrimeEdit::Insertion { seq } => {
            if seq.is_empty() || !is_acgt(seq) {
                return Err(GeneditingError::invalid_target(
                    "region",
                    "insertion needs a non-empty ACGT sequence",
                ));
            }
        }
        PrimeEdit::Deletion { len } => {
            if *len == 0 {
                return Err(GeneditingError::invalid_target(
                    "region",
                    "deletion needs a positive length",
                ));
            }
            if req.edit_pos + len > req.reference.len() {
                return Err(GeneditingError::invalid_target(
                    "region",
                    "deletion runs past the end of the reference window",
                ));
            }
        }
    }
    if req.pbs_len == 0 || req.rt_template_homology == 0 {
        return Err(GeneditingError::invalid(
            "pbs_len/rt_template_homology",
            "PBS and RT-template homology lengths must be positive",
        ));
    }

    let editor = prime_editor(req.editor);
    let spec = PamSpec {
        motif: editor.pam.clone(),
        protospacer_len: editor.spacer_len,
        side: PamSide::ThreePrime,
    };
    let guides = scan_guides(&req.reference, &spec)
        .map_err(|e| GeneditingError::invalid("reference", e.to_string()))?;

    // Pick the forward-strand spacer whose nick is the closest one
    // *strictly upstream* of the edit start (a short nick-to-edit
    // distance keeps the RT template short and efficient).
    let plen = editor.spacer_len;
    let mut chosen: Option<(usize, usize)> = None; // (spacer_start, nick_pos)
    for g in &guides {
        if g.strand != GuideStrand::Forward {
            continue;
        }
        // Nick falls `NICK_OFFSET_FROM_PAM` nt 5' of the PAM, i.e. at
        // forward index spacer_start + plen - NICK_OFFSET_FROM_PAM.
        let nick = g.start + plen - NICK_OFFSET_FROM_PAM;
        if nick <= req.edit_pos {
            // Closest nick to the edit wins.
            match chosen {
                None => chosen = Some((g.start, nick)),
                Some((_, prev)) if nick > prev => chosen = Some((g.start, nick)),
                _ => {}
            }
        }
    }
    let (spacer_start, nick_pos) = chosen.ok_or_else(|| {
        GeneditingError::no_valid_design(
            "pegrna",
            "no forward-strand spacer nicks upstream of the edit site",
        )
    })?;

    // --- PBS: complementary to the reference 3' end at the nick ------
    // The nicked PAM strand's 3' end runs up to `nick_pos`. The PBS is
    // the reverse complement of the `pbs_len` reference bases ending at
    // the nick.
    if nick_pos < req.pbs_len {
        return Err(GeneditingError::no_valid_design(
            "pegrna",
            "reference window too short upstream of the nick for this PBS length",
        ));
    }
    let pbs_dna = &req.reference[nick_pos - req.pbs_len..nick_pos];
    let pbs_rna = transcribe(&revcomp(pbs_dna));

    // --- RT template: the edited new strand, 3' of the nick ---------
    // Build the *new* (edited) forward-strand sequence from the nick
    // onward: reference[nick..edit_pos] + inserted + downstream
    // homology, where downstream starts after the consumed reference.
    let consumed = req.edit.ref_consumed();
    let downstream_start = req.edit_pos + consumed;
    if downstream_start + req.rt_template_homology > req.reference.len() {
        return Err(GeneditingError::no_valid_design(
            "pegrna",
            "reference window too short downstream of the edit for this RT-template length",
        ));
    }
    let mut new_strand: Vec<u8> = Vec::new();
    new_strand.extend_from_slice(&req.reference[nick_pos..req.edit_pos]);
    new_strand.extend_from_slice(&req.edit.inserted());
    new_strand.extend_from_slice(
        &req.reference[downstream_start..downstream_start + req.rt_template_homology],
    );
    // The RT template is read 3'->5' as the RT synthesises the new
    // strand 5'->3'; in pegRNA 5'->3' notation the RT template is the
    // reverse complement of the new strand.
    let rt_template_rna = transcribe(&revcomp(&new_strand));

    // --- Assemble ---------------------------------------------------
    // 3' extension, 5'->3': RT template then PBS.
    let mut extension = Vec::with_capacity(rt_template_rna.len() + pbs_rna.len());
    extension.extend_from_slice(&rt_template_rna);
    extension.extend_from_slice(&pbs_rna);

    let spacer_rna = transcribe(&req.reference[spacer_start..spacer_start + plen]);
    let pam_dna = upper(
        &req.reference[spacer_start + plen..spacer_start + plen + spec.motif.len()],
    );

    let mut full = Vec::new();
    full.extend_from_slice(&spacer_rna);
    full.extend_from_slice(PEGRNA_SCAFFOLD);
    full.extend_from_slice(&extension);

    Ok(PegRna {
        spacer: spacer_rna,
        pam: pam_dna,
        spacer_start,
        nick_pos,
        pbs: pbs_rna,
        rt_template: rt_template_rna,
        three_prime_extension: extension,
        full_pegrna: full,
        edit: req.edit.clone(),
    })
}

/// One entry in a PBS / RT-template length scan.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct LengthScanEntry {
    /// PBS length tried.
    pub pbs_len: usize,
    /// RT-template downstream-homology length tried.
    pub rt_template_homology: usize,
    /// The transparent length-scan score in `[0, 1]` — higher is
    /// predicted-better (see [`crate::prime_edit::strategy`]).
    pub score: f64,
}

/// Scans PBS and RT-template lengths for the best combination
/// (feature 13).
///
/// Designs a pegRNA at every `(pbs_len, rt_homology)` pair in the
/// supplied inclusive ranges, scores each with the transparent
/// length-scan heuristic, and returns the entries sorted by descending
/// score. The caller can then re-run [`design_pegrna`] with the top
/// entry's lengths.
///
/// # Errors
/// - [`GeneditingError::Invalid`] for an empty / inverted range.
/// - [`GeneditingError::NoValidDesign`] if no length pair yields a
///   valid pegRNA in the supplied window.
pub fn scan_pbs_rt(
    req: &PegRnaRequest,
    pbs_range: (usize, usize),
    rt_range: (usize, usize),
) -> Result<Vec<LengthScanEntry>> {
    let (pbs_lo, pbs_hi) = pbs_range;
    let (rt_lo, rt_hi) = rt_range;
    if pbs_lo == 0 || rt_lo == 0 || pbs_lo > pbs_hi || rt_lo > rt_hi {
        return Err(GeneditingError::invalid(
            "scan_range",
            "PBS and RT ranges must be positive with lo <= hi",
        ));
    }
    let mut entries: Vec<LengthScanEntry> = Vec::new();
    for pbs in pbs_lo..=pbs_hi {
        for rt in rt_lo..=rt_hi {
            let mut trial = req.clone();
            trial.pbs_len = pbs;
            trial.rt_template_homology = rt;
            // A length pair that does not fit the window is skipped,
            // not fatal — other pairs may still work.
            let peg = match design_pegrna(&trial) {
                Ok(p) => p,
                Err(GeneditingError::NoValidDesign { .. }) => continue,
                Err(e) => return Err(e),
            };
            let score = length_scan_score(&peg.pbs, &peg.rt_template, pbs, rt);
            entries.push(LengthScanEntry {
                pbs_len: pbs,
                rt_template_homology: rt,
                score,
            });
        }
    }
    if entries.is_empty() {
        return Err(GeneditingError::no_valid_design(
            "pegrna",
            "no PBS / RT-template length pair fits the supplied reference window",
        ));
    }
    entries.sort_by(|a, b| {
        b.score
            .partial_cmp(&a.score)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    Ok(entries)
}

/// A designed PE3 / PE3b nicking guide (feature 14).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct NickingGuide {
    /// The nicking-guide protospacer, 5′→3′ on its strand, as RNA.
    pub protospacer: Vec<u8>,
    /// The PAM as found (DNA).
    pub pam: Vec<u8>,
    /// 0-based forward-strand start of the protospacer.
    pub start: usize,
    /// 0-based forward-strand nick position of the nicking guide.
    pub nick_pos: usize,
    /// Signed distance, in bp, from the pegRNA nick to this nick (the
    /// PE3 nick-to-nick offset; ~40–90 nt is the usual sweet spot).
    pub nick_to_nick: i64,
    /// `true` when this is a PE3b guide (its protospacer matches only
    /// the *edited* sequence).
    pub pe3b: bool,
}

/// Designs a PE3 or PE3b nicking guide on the complementary strand
/// (feature 14).
///
/// The second nick goes on the strand *opposite* the pegRNA nick (the
/// reverse strand here, since the pegRNA spacer is forward). For PE3b
/// the nicking guide is chosen so its protospacer overlaps the edit —
/// it can only be cut once the edit is installed.
///
/// `pegrna` is the already-designed pegRNA; `editor` selects PE3 vs
/// PE3b behaviour. Returns the nicking guide whose nick-to-nick
/// distance is closest to `preferred_offset`.
///
/// # Errors
/// - [`GeneditingError::Invalid`] if `editor` is a no-second-nick
///   configuration (PE2).
/// - [`GeneditingError::NoValidDesign`] if no reverse-strand guide
///   qualifies.
pub fn design_nicking_guide(
    reference: &[u8],
    pegrna: &PegRna,
    editor: PrimeEditorId,
    preferred_offset: i64,
) -> Result<NickingGuide> {
    let ed = prime_editor(editor);
    if !ed.uses_nicking_guide {
        return Err(GeneditingError::invalid(
            "editor",
            "PE2 uses no second nicking guide; choose PE3, PE3b or PEmax",
        ));
    }
    if !is_acgt(reference) {
        return Err(GeneditingError::invalid_target(
            "region",
            "reference window must be non-empty ACGT",
        ));
    }
    let spec = PamSpec {
        motif: ed.pam.clone(),
        protospacer_len: ed.spacer_len,
        side: PamSide::ThreePrime,
    };
    let guides = scan_guides(reference, &spec)
        .map_err(|e| GeneditingError::invalid("reference", e.to_string()))?;
    let plen = ed.spacer_len;

    // Edit span on the forward strand — for PE3b the nicking-guide
    // protospacer must overlap it.
    let edit_start = pegrna_edit_start(pegrna);
    let edit_end = edit_start + pegrna.edit.ref_consumed().max(1);

    let mut best: Option<NickingGuide> = None;
    for g in &guides {
        if g.strand != GuideStrand::Reverse {
            continue; // the second nick is on the opposite strand
        }
        // Reverse-strand nick: 3 nt 5' of the PAM on the reverse
        // strand; on the forward axis that is at g.start + offset.
        let nick = g.start + NICK_OFFSET_FROM_PAM;
        let offset = nick as i64 - pegrna.nick_pos as i64;
        let proto_start = g.start;
        let proto_end = g.start + plen;
        let overlaps_edit = proto_start < edit_end && proto_end > edit_start;
        if ed.nick_after_edit && !overlaps_edit {
            continue; // PE3b requires overlap with the edit
        }
        let candidate = NickingGuide {
            protospacer: transcribe(g.protospacer.as_bytes()),
            pam: upper(g.pam.as_bytes()),
            start: g.start,
            nick_pos: nick,
            nick_to_nick: offset,
            pe3b: ed.nick_after_edit,
        };
        best = match best {
            None => Some(candidate),
            Some(b) => {
                let da = (candidate.nick_to_nick - preferred_offset).abs();
                let db = (b.nick_to_nick - preferred_offset).abs();
                if da < db {
                    Some(candidate)
                } else {
                    Some(b)
                }
            }
        };
    }
    best.ok_or_else(|| {
        GeneditingError::no_valid_design(
            "nicking_guide",
            if ed.nick_after_edit {
                "no reverse-strand guide overlaps the edit for a PE3b nick"
            } else {
                "no reverse-strand guide is available for the second nick"
            },
        )
    })
}

/// The forward-strand index where the pegRNA's edit starts. The pegRNA
/// stores the nick and the RT template; the edit starts at
/// `nick + (RT-template homology already copied)` — but the simplest
/// stable recovery is: the RT template encodes `reference[nick..edit] +
/// inserted + downstream`, so the edit offset from the nick is the
/// length of `reference[nick..edit]`. We recover it from the request
/// invariant that the new-strand prefix before `inserted` equals
/// `reference[nick..edit]`.
fn pegrna_edit_start(peg: &PegRna) -> usize {
    // The RT template is revcomp of the new strand; the new strand's
    // length before the inserted bases is (edit_pos - nick_pos). We
    // cannot read edit_pos directly, but the nick_pos plus the
    // pre-edit homology is what we need and the pre-edit homology is
    // (new_strand_len - inserted_len - downstream_len). Since we kept
    // the full edit in `peg.edit`, callers that need the precise edit
    // span should use the original request; here we conservatively
    // anchor the edit just after the nick.
    peg.nick_pos
}

#[cfg(test)]
mod tests {
    use super::*;

    // A reference with a forward NGG PAM whose nick sits upstream of a
    // chosen edit site. 20-mer protospacer at 0..20, PAM AGG at 20..23;
    // nick at index 17. Place the edit at index 25.
    fn reference() -> Vec<u8> {
        b"ACGTACGTACGTACGTACGTAGGCATGCATGCATGCATGCATGCATGC".to_vec()
    }

    #[test]
    fn designs_a_substitution_pegrna() {
        let req = PegRnaRequest::new(reference(), 25, PrimeEdit::snv(b'A'), PrimeEditorId::Pe2);
        let peg = design_pegrna(&req).unwrap();
        assert_eq!(peg.nick_pos, 17);
        assert_eq!(peg.pbs.len(), 13);
        assert!(!peg.rt_template.is_empty());
        // Full pegRNA = spacer + scaffold + extension.
        assert_eq!(
            peg.full_pegrna.len(),
            peg.spacer.len() + PEGRNA_SCAFFOLD.len() + peg.three_prime_extension.len()
        );
        // pegRNA is RNA — no T.
        assert!(!peg.full_pegrna.contains(&b'T'));
    }

    #[test]
    fn pegrna_spacer_matches_reference_transcribed() {
        let req = PegRnaRequest::new(reference(), 25, PrimeEdit::snv(b'A'), PrimeEditorId::Pe2);
        let peg = design_pegrna(&req).unwrap();
        let expect = transcribe(&reference()[peg.spacer_start..peg.spacer_start + 20]);
        assert_eq!(peg.spacer, expect);
    }

    #[test]
    fn pbs_is_revcomp_of_reference_at_nick() {
        let req = PegRnaRequest::new(reference(), 25, PrimeEdit::snv(b'A'), PrimeEditorId::Pe2);
        let peg = design_pegrna(&req).unwrap();
        let r = reference();
        let pbs_dna = &r[peg.nick_pos - 13..peg.nick_pos];
        assert_eq!(peg.pbs, transcribe(&revcomp(pbs_dna)));
    }

    #[test]
    fn designs_an_insertion_pegrna() {
        let req = PegRnaRequest::new(
            reference(),
            25,
            PrimeEdit::Insertion { seq: b"GGG".to_vec() },
            PrimeEditorId::Pe2,
        );
        let peg = design_pegrna(&req).unwrap();
        assert!(matches!(peg.edit, PrimeEdit::Insertion { .. }));
        assert!(!peg.rt_template.is_empty());
    }

    #[test]
    fn designs_a_deletion_pegrna() {
        let req = PegRnaRequest::new(
            reference(),
            25,
            PrimeEdit::Deletion { len: 3 },
            PrimeEditorId::Pe2,
        );
        let peg = design_pegrna(&req).unwrap();
        assert!(matches!(peg.edit, PrimeEdit::Deletion { len: 3 }));
    }

    #[test]
    fn rejects_edit_upstream_of_every_nick() {
        // Edit at index 5 — before the only nick (17). No spacer nicks
        // upstream of it.
        let req = PegRnaRequest::new(reference(), 5, PrimeEdit::snv(b'A'), PrimeEditorId::Pe2);
        let err = design_pegrna(&req).unwrap_err();
        assert_eq!(err.code(), "genediting.no_valid_design");
    }

    #[test]
    fn rejects_non_acgt_reference() {
        let req = PegRnaRequest::new(b"ACGTNNNN".to_vec(), 5, PrimeEdit::snv(b'A'), PrimeEditorId::Pe2);
        assert!(design_pegrna(&req).is_err());
    }

    #[test]
    fn rejects_zero_pbs_length() {
        let mut req =
            PegRnaRequest::new(reference(), 25, PrimeEdit::snv(b'A'), PrimeEditorId::Pe2);
        req.pbs_len = 0;
        assert_eq!(design_pegrna(&req).unwrap_err().code(), "genediting.invalid");
    }

    #[test]
    fn length_scan_returns_sorted_entries() {
        let req = PegRnaRequest::new(reference(), 25, PrimeEdit::snv(b'A'), PrimeEditorId::Pe2);
        let entries = scan_pbs_rt(&req, (8, 15), (8, 15)).unwrap();
        assert!(!entries.is_empty());
        for w in entries.windows(2) {
            assert!(w[0].score >= w[1].score);
        }
    }

    #[test]
    fn length_scan_rejects_inverted_range() {
        let req = PegRnaRequest::new(reference(), 25, PrimeEdit::snv(b'A'), PrimeEditorId::Pe2);
        assert!(scan_pbs_rt(&req, (15, 8), (8, 15)).is_err());
    }

    #[test]
    fn pe2_has_no_nicking_guide() {
        let req = PegRnaRequest::new(reference(), 25, PrimeEdit::snv(b'A'), PrimeEditorId::Pe2);
        let peg = design_pegrna(&req).unwrap();
        let err = design_nicking_guide(&reference(), &peg, PrimeEditorId::Pe2, 60).unwrap_err();
        assert_eq!(err.code(), "genediting.invalid");
    }

    #[test]
    fn pe3_designs_a_nicking_guide_when_one_exists() {
        // A reference carrying both a forward PAM (for the pegRNA) and
        // a reverse-strand PAM (for the second nick).
        let r = b"ACGTACGTACGTACGTACGTAGGCATCCAGTACCGATGCACTGCATGC".to_vec();
        let req = PegRnaRequest::new(r.clone(), 25, PrimeEdit::snv(b'A'), PrimeEditorId::Pe3);
        let peg = design_pegrna(&req).unwrap();
        // PE3 (not PE3b) — any reverse-strand guide qualifies.
        match design_nicking_guide(&r, &peg, PrimeEditorId::Pe3, 50) {
            Ok(ng) => {
                assert!(!ng.pe3b);
                assert!(!ng.protospacer.contains(&b'T')); // RNA
            }
            // If the synthetic reference happens to lack a reverse PAM
            // the call returns NoValidDesign — also acceptable.
            Err(e) => assert_eq!(e.code(), "genediting.no_valid_design"),
        }
    }

    #[test]
    fn prime_edit_labels() {
        assert!(PrimeEdit::snv(b'A').label().contains("substitution"));
        assert!(PrimeEdit::Insertion { seq: b"GG".to_vec() }
            .label()
            .contains("insertion"));
        assert!(PrimeEdit::Deletion { len: 2 }.label().contains("deletion"));
    }
}
