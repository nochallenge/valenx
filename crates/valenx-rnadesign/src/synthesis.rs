//! Group E (part 1) — output and synthesis plan (features 19–21).
//!
//! A validated design is not yet a thing you can order — it has to be
//! turned into a **synthesis-ready package**: a DNA template, a way to
//! make RNA from it, and the reagents to amplify the template. This
//! module builds that package.
//!
//! - [`generate_dna_template`] (feature 19) — reverse-transcribes the
//!   RNA design to its DNA template and prepends an RNA-polymerase
//!   promoter (T7 or SP6) for in-vitro transcription.
//! - [`plan_ivt`] (feature 20) — an in-vitro-transcription plan: the
//!   template, the polymerase, the NTP mix (including the modified
//!   nucleotide when the scheme calls for it), the 5′-cap strategy, the
//!   poly-A strategy, and honest yield / handling notes.
//! - [`plan_synthesis`] (feature 21) — assembles the full
//!   [`SynthesisPackage`]: the final RNA sequence, the DNA template,
//!   PCR primers for the template (via [`valenx_bioseq`]), and an
//!   annotated construct map.
//!
//! ## v1 scope — honest framing
//!
//! The promoter sequences are the canonical published consensus
//! sequences; the IVT plan encodes the standard well-established
//! protocol choices (polymerase, cap strategy, modified-NTP
//! substitution) — it is **not** a quantitative reaction model and the
//! yield notes are qualitative. The synthesis package is a *starting
//! point* for wet-lab work: a real order needs the template
//! synthesised or PCR-amplified, the transcription reaction optimised,
//! and the product validated. Nothing here guarantees a successful
//! synthesis.

use crate::design::RnaDesign;
use crate::error::{Result, RnaDesignError};
use crate::goal::DesignConstraints;
use serde::{Deserialize, Serialize};
use valenx_bioseq::primer::design::{design_primers, PrimerConstraints};
use valenx_bioseq::{Seq, SeqKind};
use valenx_genediting::mrna::uridine::ModifiedNucleoside;

/// An RNA-polymerase promoter for in-vitro transcription.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum Promoter {
    /// The bacteriophage **T7** RNA-polymerase promoter — the workhorse
    /// of in-vitro transcription.
    T7,
    /// The bacteriophage **SP6** RNA-polymerase promoter.
    Sp6,
}

impl Promoter {
    /// The promoter's double-stranded DNA sequence (top strand, 5′→3′).
    /// Transcription begins at the final `G`, which becomes the RNA's
    /// `+1` nucleotide.
    pub fn dna_sequence(self) -> &'static str {
        match self {
            // Canonical T7 class III promoter consensus.
            Promoter::T7 => "TAATACGACTCACTATAG",
            // Canonical SP6 promoter consensus.
            Promoter::Sp6 => "ATTTAGGTGACACTATAG",
        }
    }

    /// The polymerase enzyme name.
    pub fn polymerase(self) -> &'static str {
        match self {
            Promoter::T7 => "T7 RNA polymerase",
            Promoter::Sp6 => "SP6 RNA polymerase",
        }
    }

    /// A short human-readable name.
    pub fn name(self) -> &'static str {
        match self {
            Promoter::T7 => "T7",
            Promoter::Sp6 => "SP6",
        }
    }
}

/// The DNA template for transcribing an RNA design (feature 19).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct DnaTemplate {
    /// The full double-stranded DNA template (top strand, 5′→3′):
    /// `promoter + coding strand of the RNA`.
    pub sequence: Vec<u8>,
    /// The primary promoter prepended to the template.
    pub promoter: Promoter,
    /// `true` when a second promoter (the other of T7 / SP6) was also
    /// recorded as an alternative in the construct map.
    pub dual_promoter: bool,
    /// 0-based `[start, end)` span of the promoter within
    /// [`sequence`](Self::sequence).
    pub promoter_span: (usize, usize),
    /// 0-based `[start, end)` span of the transcribed region (the RNA's
    /// DNA coding strand) within [`sequence`](Self::sequence).
    pub transcript_span: (usize, usize),
}

impl DnaTemplate {
    /// The total template length in base pairs.
    pub fn len(&self) -> usize {
        self.sequence.len()
    }

    /// `true` when the template is empty (never produced).
    pub fn is_empty(&self) -> bool {
        self.sequence.is_empty()
    }

    /// The DNA template as a `&str`.
    pub fn sequence_str(&self) -> &str {
        std::str::from_utf8(&self.sequence).unwrap_or("")
    }
}

/// Reverse-transcribes an RNA design to a DNA template with an
/// RNA-polymerase promoter (feature 19).
///
/// The RNA design is back-transcribed (`U → T`) to its coding-strand
/// DNA and `promoter`'s consensus sequence is prepended; the result is
/// the template a polymerase transcribes back into the RNA design.
///
/// # Errors
/// [`RnaDesignError::Invalid`] if the design sequence is empty.
pub fn generate_dna_template(design: &RnaDesign, promoter: Promoter) -> Result<DnaTemplate> {
    if design.sequence.is_empty() {
        return Err(RnaDesignError::invalid(
            "design",
            "cannot make a DNA template for an empty design",
        ));
    }
    // Back-transcribe the RNA design to coding-strand DNA.
    let coding_dna: Vec<u8> = design
        .sequence
        .iter()
        .map(|&b| match b.to_ascii_uppercase() {
            b'U' => b'T',
            other => other,
        })
        .collect();

    let promoter_bytes = promoter.dna_sequence().as_bytes();
    let mut sequence = Vec::with_capacity(promoter_bytes.len() + coding_dna.len());
    sequence.extend_from_slice(promoter_bytes);
    sequence.extend_from_slice(&coding_dna);

    let promoter_span = (0, promoter_bytes.len());
    let transcript_span = (promoter_bytes.len(), sequence.len());

    Ok(DnaTemplate {
        sequence,
        promoter,
        dual_promoter: false,
        promoter_span,
        transcript_span,
    })
}

/// The cap strategy of an in-vitro-transcription plan.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum CapStrategy {
    /// Co-transcriptional capping with a trinucleotide cap analog
    /// (CleanCap-class) — the cap is incorporated during transcription.
    CoTranscriptional,
    /// Enzymatic post-transcriptional capping (vaccinia capping enzyme
    /// + 2′-O-methyltransferase for a cap-1 chemistry).
    Enzymatic,
    /// No cap — for a non-coding structural RNA that is not translated.
    Uncapped,
}

impl CapStrategy {
    /// A short human-readable name.
    pub fn name(self) -> &'static str {
        match self {
            CapStrategy::CoTranscriptional => "co-transcriptional cap analog",
            CapStrategy::Enzymatic => "enzymatic (post-transcriptional) capping",
            CapStrategy::Uncapped => "uncapped",
        }
    }
}

/// The poly-A strategy of an in-vitro-transcription plan.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum PolyAStrategy {
    /// The poly-A tail is **template-encoded** — the DNA template
    /// carries the poly-T stretch, so the tail is transcribed directly.
    TemplateEncoded {
        /// The encoded tail length, nt.
        length: usize,
    },
    /// The poly-A tail is added **enzymatically** with poly-A
    /// polymerase after transcription.
    Enzymatic {
        /// The target tail length, nt.
        target_length: usize,
    },
    /// No poly-A tail (a non-coding structural RNA).
    None,
}

impl PolyAStrategy {
    /// A short human-readable name.
    pub fn name(self) -> &'static str {
        match self {
            PolyAStrategy::TemplateEncoded { .. } => "template-encoded poly-A",
            PolyAStrategy::Enzymatic { .. } => "enzymatic poly-A tailing",
            PolyAStrategy::None => "no poly-A tail",
        }
    }
}

/// An in-vitro-transcription plan for an RNA design (feature 20).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct IvtPlan {
    /// The polymerase enzyme to use.
    pub polymerase: String,
    /// The four NTPs in the reaction (ASCII labels), with the modified
    /// uridine substituted in when the scheme calls for it.
    pub ntp_mix: Vec<String>,
    /// The modified nucleoside replacing standard UTP, when one is used.
    pub modified_nucleoside: ModifiedNucleoside,
    /// `true` when the reaction uses a modified uridine (`Ψ`, `m1Ψ`,
    /// `5moU`).
    pub uses_modified_uridine: bool,
    /// The cap strategy.
    pub cap_strategy: CapStrategy,
    /// The poly-A strategy.
    pub poly_a_strategy: PolyAStrategy,
    /// Honest yield / handling notes (one line per consideration).
    pub notes: Vec<String>,
}

/// Plans the in-vitro transcription of an RNA design (feature 20).
///
/// `template` is the DNA template; `is_coding` flags whether the RNA is
/// a translated mRNA (which needs a cap and a poly-A tail) or a
/// non-coding structural RNA (which needs neither); `nucleoside` is the
/// modified-uridine scheme from the design constraints.
///
/// # Errors
/// This function does not fail in v1 — it returns a plan for any
/// well-formed template — but returns `Result` for forward
/// compatibility.
pub fn plan_ivt(
    template: &DnaTemplate,
    is_coding: bool,
    nucleoside: ModifiedNucleoside,
    poly_a_len: usize,
) -> Result<IvtPlan> {
    let uses_modified = nucleoside.reduces_immunogenicity();
    // The NTP mix: ATP, CTP, GTP, and either UTP or the modified analog.
    let uridine_label = match nucleoside {
        ModifiedNucleoside::Uridine => "UTP",
        ModifiedNucleoside::Pseudouridine => "pseudo-UTP (ΨTP)",
        ModifiedNucleoside::N1MethylPseudouridine => "N1-methylpseudo-UTP (m1ΨTP)",
        ModifiedNucleoside::Methoxyuridine => "5-methoxy-UTP (5moUTP)",
    };
    let ntp_mix = vec![
        "ATP".to_string(),
        "CTP".to_string(),
        "GTP".to_string(),
        uridine_label.to_string(),
    ];

    // Cap / poly-A strategy depends on whether the RNA is translated.
    let (cap_strategy, poly_a_strategy) = if is_coding {
        (
            CapStrategy::CoTranscriptional,
            if poly_a_len > 0 {
                PolyAStrategy::TemplateEncoded { length: poly_a_len }
            } else {
                PolyAStrategy::Enzymatic { target_length: 100 }
            },
        )
    } else {
        (CapStrategy::Uncapped, PolyAStrategy::None)
    };

    let mut notes = Vec::new();
    notes.push(format!(
        "Transcribe with {} from the {} promoter.",
        template.promoter.polymerase(),
        template.promoter.name(),
    ));
    if uses_modified {
        notes.push(format!(
            "Substitute every UTP with {} — reduces innate-immune (RIG-I / TLR) \
             activation and typically raises protein output.",
            nucleoside.name(),
        ));
    } else {
        notes.push(
            "Standard UTP — for a therapeutic mRNA consider an m1Ψ substitution to lower \
             immunogenicity."
                .to_string(),
        );
    }
    match cap_strategy {
        CapStrategy::CoTranscriptional => notes.push(
            "Add a co-transcriptional trinucleotide cap analog for a cap-1 5' end in one \
             transcription step."
                .to_string(),
        ),
        CapStrategy::Enzymatic => notes.push(
            "Cap enzymatically after transcription (vaccinia capping enzyme + 2'-O-MTase)."
                .to_string(),
        ),
        CapStrategy::Uncapped => notes.push(
            "No cap — this is a non-coding structural RNA that is not translated."
                .to_string(),
        ),
    }
    match poly_a_strategy {
        PolyAStrategy::TemplateEncoded { length } => notes.push(format!(
            "The {length}-nt poly-A tail is template-encoded — it is transcribed directly; \
             keep the encoding poly-T run intact during template preparation.",
        )),
        PolyAStrategy::Enzymatic { target_length } => notes.push(format!(
            "Add a ~{target_length}-nt poly-A tail enzymatically with poly-A polymerase.",
        )),
        PolyAStrategy::None => {}
    }
    notes.push(
        "Yield is qualitative here — a real IVT reaction needs the NTP / Mg2+ ratio, the \
         template amount and the reaction time optimised, then DNase-treated and purified."
            .to_string(),
    );
    notes.push(
        "After transcription, validate the RNA: check integrity (capillary electrophoresis), \
         confirm the fold / function in vitro, and remove dsRNA byproducts (which trigger \
         innate immunity)."
            .to_string(),
    );

    Ok(IvtPlan {
        polymerase: template.promoter.polymerase().to_string(),
        ntp_mix,
        modified_nucleoside: nucleoside,
        uses_modified_uridine: uses_modified,
        cap_strategy,
        poly_a_strategy,
        notes,
    })
}

/// One annotated region in a construct map.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct ConstructRegion {
    /// 0-based `[start, end)` span of the region in the DNA template.
    pub span: (usize, usize),
    /// The region label (`"T7 promoter"`, `"5'UTR"`, `"CDS"`, …).
    pub label: String,
}

impl ConstructRegion {
    /// The region length.
    pub fn len(&self) -> usize {
        self.span.1.saturating_sub(self.span.0)
    }

    /// `true` when the region spans zero bases.
    pub fn is_empty(&self) -> bool {
        self.span.1 <= self.span.0
    }
}

/// A pair of PCR primers for amplifying the DNA template.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct TemplatePrimers {
    /// The forward primer sequence, 5′→3′.
    pub forward: String,
    /// The reverse primer sequence, 5′→3′.
    pub reverse: String,
    /// The forward primer's melting temperature, °C.
    pub forward_tm: f64,
    /// The reverse primer's melting temperature, °C.
    pub reverse_tm: f64,
    /// The amplicon size in base pairs.
    pub amplicon_size: usize,
}

/// The full synthesis-order package for a validated RNA design
/// (feature 21).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct SynthesisPackage {
    /// The final RNA design sequence (`A C G U`).
    pub rna_sequence: Vec<u8>,
    /// The DNA template (`promoter + transcribed region`).
    pub dna_template: Vec<u8>,
    /// The DNA-template detail.
    pub template: DnaTemplate,
    /// The in-vitro-transcription plan.
    pub ivt_plan: IvtPlan,
    /// PCR primers for amplifying the DNA template, or `None` if no
    /// primer pair could be designed (a very short template).
    pub primers: Option<TemplatePrimers>,
    /// The annotated construct map (the regions of the DNA template).
    pub construct_map: Vec<ConstructRegion>,
    /// Human-readable notes for the synthesis order.
    pub notes: Vec<String>,
}

impl SynthesisPackage {
    /// The RNA design length in nucleotides.
    pub fn rna_len(&self) -> usize {
        self.rna_sequence.len()
    }

    /// The DNA-template length in base pairs.
    pub fn template_len(&self) -> usize {
        self.dna_template.len()
    }

    /// `true` when the package is empty (never produced).
    pub fn is_empty(&self) -> bool {
        self.rna_sequence.is_empty()
    }
}

/// Assembles the full synthesis-order package for a validated RNA
/// design (feature 21).
///
/// Builds the DNA template (with the chosen `promoter`), plans the
/// in-vitro transcription, designs PCR primers for the template, and
/// produces an annotated construct map.
///
/// # Errors
/// [`RnaDesignError::Invalid`] if the design sequence is empty;
/// [`RnaDesignError::Upstream`] only if a `valenx-bioseq` call fails
/// unexpectedly (primer-design failure is handled gracefully — it
/// yields `primers: None`, not an error).
pub fn plan_synthesis(
    design: &RnaDesign,
    constraints: &DesignConstraints,
    promoter: Promoter,
) -> Result<SynthesisPackage> {
    if design.sequence.is_empty() {
        return Err(RnaDesignError::invalid(
            "design",
            "cannot plan synthesis for an empty design",
        ));
    }

    // 1) DNA template.
    let template = generate_dna_template(design, promoter)?;

    // 2) IVT plan.
    let poly_a_len = design
        .construct
        .as_ref()
        .map(|c| c.poly_a_len)
        .unwrap_or(0);
    let ivt_plan = plan_ivt(&template, design.is_coding(), constraints.nucleoside, poly_a_len)?;

    // 3) PCR primers for the DNA template — amplify the whole template.
    let primers = design_template_primers(&template);

    // 4) Construct map.
    let construct_map = build_construct_map(design, &template);

    let mut notes = Vec::new();
    notes.push(format!(
        "Synthesis package: a {}-nt RNA design transcribed from a {}-bp DNA template \
         ({} promoter).",
        design.sequence.len(),
        template.len(),
        promoter.name(),
    ));
    if primers.is_some() {
        notes.push(
            "PCR primers are provided to amplify the DNA template before transcription."
                .to_string(),
        );
    } else {
        notes.push(
            "No PCR primer pair could be designed for this short template — order the \
             template as a synthetic dsDNA fragment / plasmid insert instead."
                .to_string(),
        );
    }
    notes.push(
        "This package is the starting point for wet-lab synthesis — the template must be \
         made / amplified, the IVT reaction run and optimised, and the RNA product \
         validated. It is a plan, not a guarantee of a successful synthesis."
            .to_string(),
    );

    Ok(SynthesisPackage {
        rna_sequence: design.sequence.clone(),
        dna_template: template.sequence.clone(),
        template,
        ivt_plan,
        primers,
        construct_map,
        notes,
    })
}

/// Designs a PCR primer pair spanning the whole DNA template, returning
/// `None` if the template is too short for the primer constraints.
fn design_template_primers(template: &DnaTemplate) -> Option<TemplatePrimers> {
    let dna = Seq::new(SeqKind::Dna, &template.sequence).ok()?;
    let n = dna.len();
    // The primers flank a target that leaves room for ~20-nt primers at
    // each end. Need at least ~50 bp of template.
    if n < 60 {
        return None;
    }
    let constraints = PrimerConstraints {
        min_len: 18,
        max_len: 28,
        target_tm: 58.0,
        tm_tolerance: 12.0,
        min_gc: 0.30,
        max_gc: 0.70,
        require_gc_clamp: false,
        ..Default::default()
    };
    // Target region: the middle of the template, leaving ~24 bp of
    // flank for each primer.
    let flank = 24usize.min(n / 3);
    let pair = design_primers(&dna, flank, n - flank, constraints).ok()?;
    Some(TemplatePrimers {
        forward: pair.forward.seq.as_str().to_string(),
        reverse: pair.reverse.seq.as_str().to_string(),
        forward_tm: pair.forward.tm,
        reverse_tm: pair.reverse.tm,
        amplicon_size: pair.amplicon_size(),
    })
}

/// Builds the annotated construct map of a DNA template.
fn build_construct_map(design: &RnaDesign, template: &DnaTemplate) -> Vec<ConstructRegion> {
    let mut map = Vec::new();
    // The promoter.
    map.push(ConstructRegion {
        span: template.promoter_span,
        label: format!("{} promoter", template.promoter.name()),
    });
    // The transcribed region — broken into the five mRNA parts for a
    // coding design, or a single block for a structural RNA.
    let (t_start, _t_end) = template.transcript_span;
    if let Some(construct) = &design.construct {
        // 5'UTR, CDS, 3'UTR, poly-A — offsets within the transcript.
        let mut cursor = t_start;
        let utr5_len = construct.utr5.len();
        if utr5_len > 0 {
            map.push(ConstructRegion {
                span: (cursor, cursor + utr5_len),
                label: "5'UTR".to_string(),
            });
        }
        cursor += utr5_len;
        let cds_len = construct.cds.len();
        if cds_len > 0 {
            map.push(ConstructRegion {
                span: (cursor, cursor + cds_len),
                label: "CDS".to_string(),
            });
        }
        cursor += cds_len;
        let utr3_len = construct.utr3.len();
        if utr3_len > 0 {
            map.push(ConstructRegion {
                span: (cursor, cursor + utr3_len),
                label: "3'UTR".to_string(),
            });
        }
        cursor += utr3_len;
        if construct.poly_a_len > 0 {
            map.push(ConstructRegion {
                span: (cursor, cursor + construct.poly_a_len),
                label: format!("poly-A ({} nt)", construct.poly_a_len),
            });
        }
    } else {
        map.push(ConstructRegion {
            span: template.transcript_span,
            label: "structural RNA (transcribed region)".to_string(),
        });
    }
    map
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::design::DesignKind;
    use valenx_rnastruct::Structure;

    fn structural_design(seq: &[u8]) -> RnaDesign {
        RnaDesign {
            sequence: seq.to_vec(),
            kind: DesignKind::Structural {
                target: Structure::empty(seq.len()),
            },
            cds_span: None,
            construct: None,
            notes: Vec::new(),
        }
    }

    #[test]
    fn promoter_sequences() {
        assert_eq!(Promoter::T7.dna_sequence(), "TAATACGACTCACTATAG");
        assert!(Promoter::T7.dna_sequence().ends_with('G'));
        assert!(Promoter::Sp6.dna_sequence().ends_with('G'));
        assert!(Promoter::T7.polymerase().contains("T7"));
    }

    #[test]
    fn dna_template_prepends_promoter() {
        let d = structural_design(b"GGGGAAAACCCC");
        let t = generate_dna_template(&d, Promoter::T7).unwrap();
        // Template = promoter + back-transcribed design.
        assert!(t.sequence_str().starts_with("TAATACGACTCACTATAG"));
        assert!(t.sequence_str().ends_with("GGGGAAAACCCC"));
        assert_eq!(t.len(), 18 + 12);
        // No U in a DNA template.
        assert!(!t.sequence.contains(&b'U'));
    }

    #[test]
    fn dna_template_back_transcribes_u_to_t() {
        let d = structural_design(b"ACGUACGU");
        let t = generate_dna_template(&d, Promoter::T7).unwrap();
        // The transcribed region is the back-transcribed design.
        let (s, e) = t.transcript_span;
        assert_eq!(&t.sequence[s..e], b"ACGTACGT");
    }

    #[test]
    fn dna_template_rejects_empty_design() {
        let d = structural_design(b"");
        assert!(generate_dna_template(&d, Promoter::T7).is_err());
    }

    #[test]
    fn ivt_plan_uses_modified_uridine() {
        let d = structural_design(b"GGGGAAAACCCC");
        let t = generate_dna_template(&d, Promoter::T7).unwrap();
        let plan = plan_ivt(&t, true, ModifiedNucleoside::N1MethylPseudouridine, 100).unwrap();
        assert!(plan.uses_modified_uridine);
        assert!(plan.ntp_mix.iter().any(|n| n.contains("m1Ψ")));
        assert!(plan.notes.iter().any(|n| n.contains("immune")));
    }

    #[test]
    fn ivt_plan_uncapped_for_structural_rna() {
        let d = structural_design(b"GGGGAAAACCCC");
        let t = generate_dna_template(&d, Promoter::T7).unwrap();
        let plan = plan_ivt(&t, false, ModifiedNucleoside::Uridine, 0).unwrap();
        assert_eq!(plan.cap_strategy, CapStrategy::Uncapped);
        assert_eq!(plan.poly_a_strategy, PolyAStrategy::None);
    }

    #[test]
    fn ivt_plan_standard_uridine() {
        let d = structural_design(b"GGGGAAAACCCC");
        let t = generate_dna_template(&d, Promoter::T7).unwrap();
        let plan = plan_ivt(&t, true, ModifiedNucleoside::Uridine, 100).unwrap();
        assert!(!plan.uses_modified_uridine);
        assert!(plan.ntp_mix.contains(&"UTP".to_string()));
    }

    #[test]
    fn synthesis_package_for_a_structural_rna() {
        // A long enough sequence for primers.
        let seq = b"GGGACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGU";
        let d = structural_design(seq);
        let pkg = plan_synthesis(&d, &DesignConstraints::default(), Promoter::T7).unwrap();
        assert_eq!(pkg.rna_len(), seq.len());
        assert!(pkg.template_len() > pkg.rna_len());
        // The construct map at least has the promoter + transcript.
        assert!(pkg.construct_map.len() >= 2);
        assert!(pkg.construct_map[0].label.contains("promoter"));
    }

    #[test]
    fn synthesis_package_short_design_has_no_primers() {
        let d = structural_design(b"GGGGAAAACCCC");
        let pkg = plan_synthesis(&d, &DesignConstraints::default(), Promoter::T7).unwrap();
        // Too short for a primer pair.
        assert!(pkg.primers.is_none());
    }

    #[test]
    fn synthesis_package_rejects_empty() {
        let d = structural_design(b"");
        assert!(plan_synthesis(&d, &DesignConstraints::default(), Promoter::T7).is_err());
    }

    #[test]
    fn construct_region_geometry() {
        let r = ConstructRegion {
            span: (10, 25),
            label: "CDS".to_string(),
        };
        assert_eq!(r.len(), 15);
        assert!(!r.is_empty());
    }
}
