//! Group E (part 2) — export (feature 22).
//!
//! The final step of the workflow is **export**: turning a finished
//! design and its synthesis package into shareable files. Three
//! formats:
//!
//! - [`export_fasta`] — FASTA of the RNA design and its DNA template
//!   (via [`valenx_bioseq`]'s FASTA writer).
//! - [`export_genbank`] — an annotated GenBank flat file of the DNA
//!   template, with the construct map written out as feature
//!   annotations (via [`valenx_bioseq`]'s GenBank writer).
//! - [`export_design_report`] — a human-readable plain-text design
//!   report bundling the goal, the design, the validation verdict and
//!   the synthesis plan, with the honest in-silico-prediction
//!   disclaimer.
//!
//! ## v1 scope
//!
//! Export reuses the `valenx-bioseq` FASTA / GenBank writers, so the
//! file formats are exactly that crate's (a real v1 GenBank writer that
//! covers the universal fields). The text report is a fixed
//! human-readable layout — it is a summary for a human reader, not a
//! machine-parsed format (use the typed [`crate::driver::RnaDesignReport`]
//! and its `serde` derives for structured export).

use crate::design::{DesignKind, RnaDesign};
use crate::error::Result;
use crate::synthesis::SynthesisPackage;
use crate::validate::{ValidationReport, ValidationVerdict};
use valenx_bioseq::io::fasta::{write as fasta_write, FastaWriteOptions};
use valenx_bioseq::io::genbank;
use valenx_bioseq::record::{Location, SeqFeature, SeqRecord};
use valenx_bioseq::{Seq, SeqKind};

/// Exports the RNA design and its DNA template as a multi-record FASTA
/// string (feature 22).
///
/// The FASTA holds two records: the RNA design (`A C G U`) and the DNA
/// template (`A C G T`, with the promoter).
///
/// # Errors
/// [`crate::error::RnaDesignError::Upstream`] if either sequence is not
/// valid for its alphabet (the designers only ever emit valid
/// sequences, so this should not happen).
pub fn export_fasta(
    design: &RnaDesign,
    synthesis: &SynthesisPackage,
    id_prefix: &str,
) -> Result<String> {
    let rna = Seq::new(SeqKind::Rna, &design.sequence)?;
    let dna = Seq::new(SeqKind::Dna, &synthesis.dna_template)?;

    let rna_rec = SeqRecord::new(format!("{id_prefix}_rna"), rna)
        .with_description(format!("designed RNA ({})", design.kind.name()));
    let dna_rec = SeqRecord::new(format!("{id_prefix}_template"), dna).with_description(
        format!("DNA template ({} promoter)", synthesis.template.promoter.name()),
    );

    Ok(fasta_write(
        &[rna_rec, dna_rec],
        FastaWriteOptions::default(),
    ))
}

/// Exports the DNA template as an annotated GenBank flat file
/// (feature 22).
///
/// The construct map is written out as GenBank feature annotations (the
/// promoter, the UTRs, the CDS, the poly-A tail).
///
/// # Errors
/// [`crate::error::RnaDesignError::Upstream`] if the DNA template is not
/// a valid DNA sequence.
pub fn export_genbank(
    design: &RnaDesign,
    synthesis: &SynthesisPackage,
    id: &str,
) -> Result<String> {
    let dna = Seq::new(SeqKind::Dna, &synthesis.dna_template)?;
    let mut record = SeqRecord::new(id, dna).with_description(format!(
        "Valenx synthetic-RNA design DNA template — {} ({} promoter). \
         In-silico design; lab-validate before use.",
        design.kind.name(),
        synthesis.template.promoter.name(),
    ));
    record
        .annotations
        .insert("organism".to_string(), "synthetic construct".to_string());

    // Write the construct map out as features.
    for region in &synthesis.construct_map {
        let (start, end) = region.span;
        let feature_type = construct_feature_type(&region.label);
        let feature = SeqFeature::new(feature_type, Location::single(start, end))
            .with_qualifier("label", region.label.clone());
        record.add_feature(feature);
    }

    Ok(genbank::write(&record))
}

/// Maps a construct-map region label to a GenBank feature-type string.
fn construct_feature_type(label: &str) -> &'static str {
    let l = label.to_ascii_lowercase();
    if l.contains("promoter") {
        "promoter"
    } else if l.contains("cds") {
        "CDS"
    } else if l.contains("utr") {
        "5'UTR" // GenBank has 5'UTR / 3'UTR; a generic UTR maps here
    } else if l.contains("poly-a") {
        "polyA_signal"
    } else {
        "misc_feature"
    }
}

/// Exports a human-readable plain-text design report (feature 22).
///
/// The report bundles the design provenance, the sequence, the
/// validation verdict and the synthesis plan, prefixed with the honest
/// in-silico-prediction disclaimer. It is a summary for a human reader.
pub fn export_design_report(
    design: &RnaDesign,
    validation: &ValidationReport,
    synthesis: &SynthesisPackage,
    title: &str,
) -> String {
    let mut out = String::new();

    // --- header ------------------------------------------------------
    out.push_str("================================================================\n");
    out.push_str(&format!("  Valenx synthetic-RNA design report — {title}\n"));
    out.push_str("================================================================\n\n");

    // --- the honest disclaimer, front and centre ---------------------
    out.push_str("IN-SILICO DESIGN — READ THIS\n");
    out.push_str("----------------------------\n");
    out.push_str(
        "This is a designed candidate, predicted by a nearest-neighbor RNA\n\
         energy model and rule-based heuristics. It is a strong, validated-\n\
         in-silico hypothesis — NOT a guarantee of correct in-vivo behaviour.\n\
         The physical RNA must be synthesised (chemical synthesis or in-vitro\n\
         transcription) and validated in the wet lab: folded, assayed, and\n\
         tested for function. Treat this report as the start of wet-lab work.\n\n",
    );

    // --- design ------------------------------------------------------
    out.push_str("DESIGN\n");
    out.push_str("------\n");
    out.push_str(&format!("  Provenance : {}\n", design.kind.name()));
    out.push_str(&format!("  Length     : {} nt\n", design.len()));
    if let DesignKind::Structural { target } = &design.kind {
        out.push_str(&format!("  Target     : {}\n", target.to_dot_bracket()));
    }
    if let Some((s, e)) = design.cds_span {
        out.push_str(&format!("  CDS span   : {s}..{e} ({} nt)\n", e - s));
    }
    out.push_str("  Sequence (RNA, 5'->3'):\n");
    append_wrapped(&mut out, design.sequence_str(), 60, "    ");
    out.push('\n');
    if !design.notes.is_empty() {
        out.push_str("  Design notes:\n");
        for n in &design.notes {
            append_bullet(&mut out, n);
        }
    }
    out.push('\n');

    // --- validation --------------------------------------------------
    out.push_str("VALIDATION (in silico)\n");
    out.push_str("----------------------\n");
    out.push_str(&format!(
        "  Overall verdict : {}\n",
        verdict_label(validation.verdict)
    ));
    out.push_str(&format!(
        "  Checks          : {} pass, {} warn, {} fail\n",
        validation.count_with(ValidationVerdict::Pass),
        validation.count_with(ValidationVerdict::Warn),
        validation.count_with(ValidationVerdict::Fail),
    ));
    for check in &validation.constraint_checks {
        out.push_str(&format!(
            "    [{}] {} — {}\n",
            check.verdict.name(),
            check.name,
            check.detail,
        ));
    }
    if let Some(fb) = &validation.fold_back {
        out.push_str(&format!(
            "  Fold-back       : {:.0}% of target pairs recovered, MFE {:.2} kcal/mol\n",
            fb.structure_match_percent, fb.mfe_energy,
        ));
    }
    if let Some(en) = &validation.ensemble {
        out.push_str(&format!(
            "  Ensemble        : target probability {:.1}%, ensemble defect {:.2}\n",
            en.target_probability * 100.0,
            en.ensemble_defect,
        ));
    }
    if let Some(r) = &validation.robustness {
        out.push_str(&format!(
            "  Robustness      : mutational robustness {:.0}%, co-transcriptional {}\n",
            r.mutational_robustness * 100.0,
            if r.cotranscriptional_ok { "ok" } else { "flagged" },
        ));
    }
    out.push('\n');

    // --- synthesis plan ----------------------------------------------
    out.push_str("SYNTHESIS PLAN\n");
    out.push_str("--------------\n");
    out.push_str(&format!(
        "  DNA template : {} bp ({} promoter)\n",
        synthesis.template_len(),
        synthesis.template.promoter.name(),
    ));
    out.push_str(&format!(
        "  Polymerase   : {}\n",
        synthesis.ivt_plan.polymerase,
    ));
    out.push_str(&format!(
        "  NTP mix      : {}\n",
        synthesis.ivt_plan.ntp_mix.join(", "),
    ));
    out.push_str(&format!(
        "  Cap strategy : {}\n",
        synthesis.ivt_plan.cap_strategy.name(),
    ));
    out.push_str(&format!(
        "  Poly-A       : {}\n",
        synthesis.ivt_plan.poly_a_strategy.name(),
    ));
    match &synthesis.primers {
        Some(p) => {
            out.push_str(&format!(
                "  PCR primers  : fwd {} (Tm {:.1} C), rev {} (Tm {:.1} C), amplicon {} bp\n",
                p.forward, p.forward_tm, p.reverse, p.reverse_tm, p.amplicon_size,
            ));
        }
        None => {
            out.push_str(
                "  PCR primers  : none designed — order the template as synthetic dsDNA\n",
            );
        }
    }
    out.push_str("  Construct map:\n");
    for region in &synthesis.construct_map {
        out.push_str(&format!(
            "    {:>5}..{:<5}  {}\n",
            region.span.0, region.span.1, region.label,
        ));
    }
    out.push('\n');

    out.push_str("DNA template (5'->3'):\n");
    append_wrapped(
        &mut out,
        std::str::from_utf8(&synthesis.dna_template).unwrap_or(""),
        60,
        "  ",
    );
    out.push('\n');

    out.push_str("----------------------------------------------------------------\n");
    out.push_str(
        "End of report. Reminder: this is an in-silico prediction — synthesise\n\
         and lab-validate the design before relying on it.\n",
    );

    out
}

/// A human-readable label for a validation verdict.
fn verdict_label(verdict: ValidationVerdict) -> &'static str {
    match verdict {
        ValidationVerdict::Pass => "PASS (in-silico checks cleared)",
        ValidationVerdict::Warn => "WARN (passed with concerns)",
        ValidationVerdict::Fail => "FAIL (a hard requirement is violated)",
    }
}

/// Appends `text` to `out`, wrapped to `width` columns, each line
/// prefixed with `indent`.
fn append_wrapped(out: &mut String, text: &str, width: usize, indent: &str) {
    let bytes = text.as_bytes();
    if bytes.is_empty() {
        out.push_str(indent);
        out.push_str("(empty)\n");
        return;
    }
    for chunk in bytes.chunks(width.max(1)) {
        out.push_str(indent);
        out.push_str(std::str::from_utf8(chunk).unwrap_or(""));
        out.push('\n');
    }
}

/// Appends a `note` to `out` as a wrapped bullet point.
fn append_bullet(out: &mut String, note: &str) {
    // Simple word-wrap at ~72 columns for readability.
    const WIDTH: usize = 68;
    let mut line = String::new();
    let mut first = true;
    for word in note.split_whitespace() {
        if !line.is_empty() && line.len() + 1 + word.len() > WIDTH {
            out.push_str(if first { "    - " } else { "      " });
            out.push_str(&line);
            out.push('\n');
            line.clear();
            first = false;
        }
        if !line.is_empty() {
            line.push(' ');
        }
        line.push_str(word);
    }
    if !line.is_empty() {
        out.push_str(if first { "    - " } else { "      " });
        out.push_str(&line);
        out.push('\n');
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::design::DesignKind;
    use crate::goal::DesignConstraints;
    use crate::synthesis::{plan_synthesis, Promoter};
    use crate::validate::validate_design;
    use valenx_rnastruct::Structure;

    fn structural_design(seq: &[u8], db: &str) -> RnaDesign {
        RnaDesign {
            sequence: seq.to_vec(),
            kind: DesignKind::Structural {
                target: Structure::from_dot_bracket(db).unwrap(),
            },
            cds_span: None,
            construct: None,
            notes: vec!["inverse-folded".to_string()],
        }
    }

    #[test]
    fn fasta_export_has_two_records() {
        let d = structural_design(b"GGGGGGAAAACCCCCC", "((((((....))))))");
        let pkg = plan_synthesis(&d, &DesignConstraints::default(), Promoter::T7).unwrap();
        let fasta = export_fasta(&d, &pkg, "design1").unwrap();
        // Two records => two '>' headers.
        assert_eq!(fasta.matches('>').count(), 2);
        assert!(fasta.contains("design1_rna"));
        assert!(fasta.contains("design1_template"));
    }

    #[test]
    fn genbank_export_is_a_flat_file() {
        let d = structural_design(b"GGGGGGAAAACCCCCC", "((((((....))))))");
        let pkg = plan_synthesis(&d, &DesignConstraints::default(), Promoter::T7).unwrap();
        let gb = export_genbank(&d, &pkg, "DESIGN1").unwrap();
        assert!(gb.starts_with("LOCUS"));
        assert!(gb.contains("FEATURES"));
        assert!(gb.contains("ORIGIN") || gb.contains("//"));
        // The promoter feature is annotated.
        assert!(gb.contains("promoter"));
    }

    #[test]
    fn text_report_carries_the_disclaimer() {
        let d = structural_design(b"GGGGGGAAAACCCCCC", "((((((....))))))");
        let report = validate_design(&d, &DesignConstraints::default()).unwrap();
        let pkg = plan_synthesis(&d, &DesignConstraints::default(), Promoter::T7).unwrap();
        let text = export_design_report(&d, &report, &pkg, "Test design");
        // The honest framing is present.
        assert!(text.contains("IN-SILICO DESIGN"));
        assert!(text.contains("lab"));
        assert!(text.contains("NOT a guarantee"));
        // It has the major sections.
        assert!(text.contains("DESIGN"));
        assert!(text.contains("VALIDATION"));
        assert!(text.contains("SYNTHESIS PLAN"));
    }

    #[test]
    fn text_report_shows_the_sequence() {
        let d = structural_design(b"GGGGGGAAAACCCCCC", "((((((....))))))");
        let report = validate_design(&d, &DesignConstraints::default()).unwrap();
        let pkg = plan_synthesis(&d, &DesignConstraints::default(), Promoter::T7).unwrap();
        let text = export_design_report(&d, &report, &pkg, "Test");
        assert!(text.contains("GGGGGGAAAACCCCCC"));
    }

    #[test]
    fn construct_feature_type_mapping() {
        assert_eq!(construct_feature_type("T7 promoter"), "promoter");
        assert_eq!(construct_feature_type("CDS"), "CDS");
        assert_eq!(construct_feature_type("poly-A (100 nt)"), "polyA_signal");
    }

    #[test]
    fn wrap_helper_wraps() {
        let mut s = String::new();
        append_wrapped(&mut s, "ACGUACGUACGU", 4, "  ");
        // 12 chars / 4 = 3 lines.
        assert_eq!(s.lines().count(), 3);
    }
}
