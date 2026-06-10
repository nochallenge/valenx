//! Panel 1 — **Sequence** (`valenx-bioseq`).
//!
//! Load / edit a nucleotide sequence, then run reverse-complement,
//! translation, ORF finding, GC / composition analysis, a restriction
//! digest with a virtual-gel fragment list, Tm-targeted primer design
//! and in-silico PCR — all native `valenx-bioseq` calls.

use eframe::egui;

use valenx_bioseq::analysis::composition::{gc_content, gc_skew, residue_counts, shannon_entropy};
use valenx_bioseq::analysis::weight::{molecular_weight_nucleic, molecular_weight_protein};
use valenx_bioseq::analysis::tm::tm_wallace;
use valenx_bioseq::cloning::digest_map::{digest_to_gel, restriction_map};
use valenx_bioseq::cloning::restriction::{enzyme_by_name, Enzyme};
use valenx_bioseq::ops::orf::{find_orfs, OrfOptions};
use valenx_bioseq::ops::revcomp::reverse_complement;
use valenx_bioseq::ops::translate::{translate, GeneticCode, TranslateOptions};
use valenx_bioseq::primer::design::{design_primers, PrimerConstraints};
use valenx_bioseq::primer::pcr::{simulate_pcr, PcrOptions};
use valenx_bioseq::{Seq, SeqKind};

use super::common;
use crate::ValenxApp;

/// Which sub-tool of the Sequence panel is showing.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Default)]
enum Tool {
    #[default]
    Analyze,
    Translate,
    Orf,
    Restriction,
    Primers,
    Pcr,
}

/// Form + result state for the Sequence panel.
pub struct SequencePanel {
    /// The working sequence (raw user text — cleaned at run time).
    seq_text: String,
    /// Sequence kind selected by the radio row.
    kind: SeqKind,
    /// Active sub-tool.
    tool: Tool,
    /// NCBI genetic-code table id for translation / ORFs.
    code_id: u8,
    /// Minimum protein length for the ORF finder.
    orf_min_len: usize,
    /// Restriction enzyme names (comma-separated).
    enzymes: String,
    /// Primer-design target window (start, end).
    primer_start: usize,
    primer_end: usize,
    primer_target_tm: f64,
    /// PCR primers.
    pcr_fwd: String,
    pcr_rev: String,
    /// Last error (if any).
    error: Option<String>,
    /// Rendered result text.
    result: String,
    /// Undo / redo history over the editable sequence text. Pushed
    /// just before any mutating action (Reverse-complement, Load
    /// FASTA, multiline-edit committed). Wires into the host's
    /// Ctrl+Z / Ctrl+Y via [`SequencePanel::undo_edit`].
    history: crate::undo::History<String>,
}

impl Default for SequencePanel {
    fn default() -> Self {
        // The first 65 nt are the original MCS demo (EcoRI / BamHI /
        // HindIII / KpnI / SacI sites — exercises the restriction and
        // PCR sub-tools). The trailing 36 nt are a palindrome-free,
        // GC-balanced 3' flank that gives the primer-design sub-tool
        // a clean reverse-primer window — without it every reverse
        // candidate landed in the MCS and failed the self-dimer /
        // hairpin ΔG screen (no candidate satisfies the constraints).
        SequencePanel {
            seq_text:
                "ATGGCATTAGCCGGTAAATCGGATCCGAATTCAAGCTTGGTACCGAGCTCGGATCCAAGCTTTAA\
                 GTACGATCAACGTCAATGACGTTAACGCATCAAGAC"
                    .to_string(),
            kind: SeqKind::Dna,
            tool: Tool::Analyze,
            code_id: 1,
            orf_min_len: 20,
            enzymes: "EcoRI, BamHI, HindIII".to_string(),
            // The forward primer anneals to the template *before*
            // `primer_start` and the reverse primer to the template
            // *after* `primer_end`, so the target window must leave a
            // primer-sized flank on both sides. `primer_start = 21`
            // puts the forward primer in the clean 5' coding region
            // (template[0..21], ending in `G` for the GC clamp);
            // `primer_end = 65` puts the reverse primer in the clean
            // appended 3' flank (template[65..] starts with `G`,
            // giving the reverse primer a `C` GC clamp).
            primer_start: 21,
            primer_end: 65,
            primer_target_tm: 60.0,
            pcr_fwd: "ATGGCATTAGCCGGTAAATCG".to_string(),
            pcr_rev: "TTAAAGCTTGGATCCGAGCTC".to_string(),
            error: None,
            result: String::new(),
            history: crate::undo::History::new(),
        }
    }
}

impl SequencePanel {
    /// Undo the most recent sequence-text edit. Returns `true` if a
    /// snapshot was popped.
    pub fn undo_edit(&mut self) -> bool {
        if let Some(prev) = self.history.undo(self.seq_text.clone()) {
            self.seq_text = prev;
            self.error = None;
            true
        } else {
            false
        }
    }

    /// Redo the most recently undone sequence-text edit.
    pub fn redo_edit(&mut self) -> bool {
        if let Some(next) = self.history.redo(self.seq_text.clone()) {
            self.seq_text = next;
            self.error = None;
            true
        } else {
            false
        }
    }

    /// Push the current sequence text onto the undo stack so the
    /// next mutation is reversible. Called by the mutating buttons
    /// (Reverse-complement, Load FASTA) just before they overwrite
    /// the text.
    fn record(&mut self) {
        self.history.record(self.seq_text.clone());
    }

    /// Whether the most recent mutation can be undone.
    pub fn can_undo(&self) -> bool {
        self.history.can_undo()
    }
    /// Whether the most recently undone mutation can be redone.
    pub fn can_redo(&self) -> bool {
        self.history.can_redo()
    }

    /// Build a [`valenx_bioseq::record::SeqRecord`] from the current
    /// sequence text for the 2D DNA viewport. Returns `None` when the
    /// text is empty or fails to parse — the viewport then falls back to
    /// its demo plasmid. The record carries no feature annotations yet
    /// (the panel edits a bare sequence), so the viewport draws a
    /// length-accurate backbone + bp ruler for the user's real sequence.
    pub fn to_seq_record(&self) -> Option<valenx_bioseq::record::SeqRecord> {
        use valenx_bioseq::record::{Location, SeqFeature, SeqRecord, Span};

        let cleaned = common::clean_sequence(&self.seq_text);
        if cleaned.is_empty() {
            return None;
        }
        let seq = Seq::new(self.kind, cleaned.as_bytes()).ok()?;

        // Auto-annotate with detected ORFs (native bioseq analysis) so the
        // 2D viewport draws real feature arrows for the user's sequence
        // rather than a bare backbone. Nucleotide-only, and skipped on very
        // long sequences to keep this (called per repaint) cheap. Borrow
        // `seq` for ORF finding before it moves into the record.
        let mut features = Vec::new();
        if self.kind != SeqKind::Protein && seq.len() <= 200_000 {
            if let Ok(code) = GeneticCode::by_id(self.code_id) {
                if let Ok(orfs) = find_orfs(&seq, &code, OrfOptions::default()) {
                    for (i, orf) in orfs.iter().enumerate() {
                        let span =
                            Span::with_strand(orf.span.start, orf.span.end, orf.strand);
                        features.push(
                            SeqFeature::new("CDS", Location::Single(span))
                                .with_qualifier("label", format!("ORF {}", i + 1)),
                        );
                    }
                }
            }
        }

        let mut record = SeqRecord::new("Sequence", seq);
        record.features = features;
        Some(record)
    }
}

/// Keyboard-shortcut entry: run the active sub-tool. Mirrors what
/// clicking the Run button does, with no UI side-effects.
pub fn run_primary_shortcut(app: &mut crate::ValenxApp) {
    let p = &mut app.genetics.sequence;
    match p.tool {
        Tool::Analyze => run_analyze(p),
        Tool::Translate => run_translate(p),
        Tool::Orf => run_orf(p),
        Tool::Restriction => run_restriction(p),
        Tool::Primers => run_primers(p),
        Tool::Pcr => run_pcr(p),
    }
}

/// Parse the panel's text into a validated [`Seq`].
fn parse_seq(panel: &SequencePanel) -> Result<Seq, String> {
    let cleaned = common::clean_sequence(&panel.seq_text);
    if cleaned.is_empty() {
        return Err("sequence is empty".into());
    }
    Seq::new(panel.kind, cleaned.as_bytes()).map_err(|e| e.to_string())
}

/// Turn a raw `valenx-bioseq` error string into a friendlier user
/// message, when we recognise the failure pattern. Unknown messages
/// pass through unchanged so we never lose information — the goal is
/// to dress up the recognised cases, not hide the unrecognised ones.
fn friendly_error(raw: &str, action: &str) -> String {
    let low = raw.to_ascii_lowercase();
    if low.contains("no primer") || low.contains("no candidate") {
        format!(
            "No primer satisfied the constraints during {action}. \
             Try widening the Target Tm range (±3 °C), or moving the \
             primer window away from restriction sites."
        )
    } else if low.contains("invalid character") || low.contains("not a") && low.contains("base") {
        format!(
            "{raw}\n→ Check the input sequence for ambiguous characters; \
             FASTA headers + whitespace are auto-stripped, but mixed-case \
             IUPAC codes like N / W / S are not."
        )
    } else if low.contains("empty") {
        format!(
            "{action} needs at least one residue. \
             Paste a sequence into the Residues box or click Load FASTA…"
        )
    } else {
        format!("{action}: {raw}")
    }
}

/// Render the Sequence panel.
pub fn draw(app: &mut ValenxApp, ui: &mut egui::Ui) {
    let p = &mut app.genetics.sequence;

    common::section(ui, "Sequence input");
    ui.horizontal(|ui| {
        ui.label("Kind:")
            .on_hover_text("Pick the sequence alphabet — affects validation + which tools are valid.");
        ui.radio_value(&mut p.kind, SeqKind::Dna, "DNA")
            .on_hover_text("A / C / G / T over the four DNA bases.");
        ui.radio_value(&mut p.kind, SeqKind::Rna, "RNA")
            .on_hover_text("A / C / G / U over the four RNA bases.");
        ui.radio_value(&mut p.kind, SeqKind::Protein, "Protein")
            .on_hover_text("Amino-acid sequence in IUPAC single-letter code (ACDEFGHIKLMNPQRSTVWY).");
    });
    common::seq_input(ui, "seq_panel_input", "Residues (FASTA or bare):", &mut p.seq_text, 4);
    ui.horizontal(|ui| {
        let cleaned = common::clean_sequence(&p.seq_text);
        ui.label(format!("length: {} nt/aa", cleaned.len()))
            .on_hover_text("Cleaned residue count (FASTA headers + whitespace + digits stripped).");
        if ui
            .small_button("Load FASTA…")
            .on_hover_text(
                "Read a FASTA file from disk. Multi-record files load the first record.",
            )
            .clicked()
        {
            if let Some(path) = rfd::FileDialog::new()
                .add_filter("FASTA", &["fa", "fasta", "fna", "txt"])
                .pick_file()
            {
                // Round-21 H1: see biostruct loader.
                match valenx_core::io_caps::read_capped_to_string(
                    &path,
                    valenx_core::io_caps::MAX_GENETICS_FILE_BYTES as usize,
                ) {
                    Ok(text) => {
                        // Push the current text onto the undo stack
                        // so Ctrl+Z restores the user's pre-load state.
                        p.record();
                        let recs = common::parse_fasta(&text);
                        if let Some((_, body)) = recs.into_iter().next() {
                            p.seq_text = body;
                        } else {
                            p.seq_text = text;
                        }
                    }
                    Err(e) => p.error = Some(format!("read: {e}")),
                }
            }
        }
        if ui
            .small_button("Reverse-complement")
            .on_hover_text(
                "Replace the sequence with its reverse-complement. \
                 Nucleotide sequences only. Undo with Ctrl+Z.",
            )
            .clicked()
        {
            match parse_seq(p) {
                Ok(s) if s.kind() != SeqKind::Protein => match reverse_complement(&s) {
                    Ok(rc) => {
                        // Snapshot before mutating so Ctrl+Z reverses
                        // the operation.
                        p.record();
                        p.seq_text = rc.as_str().to_string();
                        p.error = None;
                    }
                    Err(e) => p.error = Some(friendly_error(&e.to_string(), "reverse-complement")),
                },
                Ok(_) => {
                    p.error = Some(
                        "Reverse-complement only works on nucleotide sequences. \
                         Switch the Kind radio to DNA or RNA."
                            .into(),
                    )
                }
                Err(e) => p.error = Some(friendly_error(&e, "reverse-complement")),
            }
        }
        // Undo / redo affordance — small buttons. Mirrors Ctrl+Z /
        // Ctrl+Y, but discoverable for first-time users who haven't
        // hit the cheat-sheet yet.
        ui.separator();
        if ui
            .add_enabled(p.can_undo(), egui::Button::new("↶"))
            .on_hover_text("Undo last edit (Ctrl+Z)")
            .clicked()
        {
            p.undo_edit();
        }
        if ui
            .add_enabled(p.can_redo(), egui::Button::new("↷"))
            .on_hover_text("Redo last edit (Ctrl+Y / Ctrl+Shift+Z)")
            .clicked()
        {
            p.redo_edit();
        }
    });
    ui.separator();

    common::section(ui, "Tool");
    ui.horizontal_wrapped(|ui| {
        ui.selectable_value(&mut p.tool, Tool::Analyze, "Composition");
        ui.selectable_value(&mut p.tool, Tool::Translate, "Translate");
        ui.selectable_value(&mut p.tool, Tool::Orf, "ORF finder");
        ui.selectable_value(&mut p.tool, Tool::Restriction, "Restriction");
        ui.selectable_value(&mut p.tool, Tool::Primers, "Primer design");
        ui.selectable_value(&mut p.tool, Tool::Pcr, "in-silico PCR");
    });
    ui.add_space(4.0);

    match p.tool {
        Tool::Analyze => draw_analyze(p, ui),
        Tool::Translate => draw_translate(p, ui),
        Tool::Orf => draw_orf(p, ui),
        Tool::Restriction => draw_restriction(p, ui),
        Tool::Primers => draw_primers(p, ui),
        Tool::Pcr => draw_pcr(p, ui),
    }

    common::error_line(ui, &p.error);
    if !p.result.is_empty() {
        ui.separator();
        common::section(ui, "Result");
        common::mono_output(ui, "seq_panel_result", &p.result, 12);
    }
}

fn draw_analyze(p: &mut SequencePanel, ui: &mut egui::Ui) {
    if common::run_button(ui, "Analyze composition") {
        run_analyze(p);
    }
}

/// Run the composition analysis — extracted from the button closure so
/// it is callable from the headless UI tests.
fn run_analyze(p: &mut SequencePanel) {
    p.error = None;
    match parse_seq(p) {
        Ok(s) => {
            let mut out = String::new();
            out.push_str(&format!("Length      : {} residues\n", s.len()));
            if s.kind() != SeqKind::Protein {
                match gc_content(&s) {
                    Ok(gc) => out.push_str(&format!("GC content  : {:.2} %\n", gc * 100.0)),
                    Err(e) => out.push_str(&format!("GC content  : (n/a — {e})\n")),
                }
                match gc_skew(&s) {
                    Ok(sk) => out.push_str(&format!("GC skew     : {sk:.4}\n")),
                    Err(e) => out.push_str(&format!("GC skew     : (n/a — {e})\n")),
                }
                match molecular_weight_nucleic(&s) {
                    Ok(mw) => out.push_str(&format!("MW (nucleic): {mw:.2} Da\n")),
                    Err(e) => out.push_str(&format!("MW (nucleic): (n/a — {e})\n")),
                }
                match tm_wallace(&s) {
                    Ok(tm) => out.push_str(&format!("Tm (Wallace): {tm:.1} °C\n")),
                    Err(e) => out.push_str(&format!("Tm (Wallace): (n/a — {e})\n")),
                }
            } else {
                match molecular_weight_protein(&s) {
                    Ok(mw) => out.push_str(&format!("MW (protein): {mw:.2} Da\n")),
                    Err(e) => out.push_str(&format!("MW (protein): (n/a — {e})\n")),
                }
            }
            out.push_str(&format!("Shannon H   : {:.4} bits\n", shannon_entropy(&s)));
            out.push_str("\nResidue counts:\n");
            for (res, n) in residue_counts(&s) {
                let pct = 100.0 * n as f64 / s.len().max(1) as f64;
                out.push_str(&format!("  {res}  {n:>6}  ({pct:5.1} %)\n"));
            }
            p.result = out;
        }
        Err(e) => p.error = Some(e),
    }
}

fn draw_translate(p: &mut SequencePanel, ui: &mut egui::Ui) {
    ui.horizontal(|ui| {
        ui.label("NCBI code table:");
        egui::ComboBox::from_id_source("seq_code_table")
            .selected_text(format!("table {}", p.code_id))
            .show_ui(ui, |ui| {
                for id in [1u8, 2, 4, 5, 6, 9, 11, 16] {
                    ui.selectable_value(&mut p.code_id, id, format!("table {id}"));
                }
            });
    });
    if common::run_button(ui, "Translate (frame +1)") {
        run_translate(p);
    }
}

/// Run the frame +1 translation — extracted for the headless UI tests.
fn run_translate(p: &mut SequencePanel) {
    p.error = None;
    match (parse_seq(p), GeneticCode::by_id(p.code_id)) {
        (Ok(s), Ok(code)) if s.kind() != SeqKind::Protein => {
            let opts = TranslateOptions { to_stop: false, start_as_met: true };
            match translate(&s, &code, opts) {
                Ok(prot) => {
                    p.result = format!(
                        "Genetic code : {} (table {})\nProtein ({} aa):\n{}",
                        code.name,
                        p.code_id,
                        prot.len(),
                        prot.as_str(),
                    );
                }
                Err(e) => p.error = Some(e.to_string()),
            }
        }
        (Ok(_), Ok(_)) => {
            p.error = Some("translation needs a nucleotide sequence".into())
        }
        (Err(e), _) => p.error = Some(e),
        (_, Err(e)) => p.error = Some(e.to_string()),
    }
}

fn draw_orf(p: &mut SequencePanel, ui: &mut egui::Ui) {
    ui.horizontal(|ui| {
        ui.label("Min protein length:");
        ui.add(egui::DragValue::new(&mut p.orf_min_len).range(1..=10_000));
        ui.label("code:");
        ui.add(egui::DragValue::new(&mut p.code_id).range(1..=33));
    });
    if common::run_button(ui, "Find ORFs (6 frames)") {
        run_orf(p);
    }
}

/// Run the 6-frame ORF finder — extracted for the headless UI tests.
fn run_orf(p: &mut SequencePanel) {
    p.error = None;
    match (parse_seq(p), GeneticCode::by_id(p.code_id)) {
        (Ok(s), Ok(code)) if s.kind() != SeqKind::Protein => {
            let opts = OrfOptions {
                atg_only: true,
                min_protein_len: p.orf_min_len,
                allow_no_stop: false,
            };
            match find_orfs(&s, &code, opts) {
                Ok(orfs) => {
                    let mut out = format!("{} ORFs (longest first):\n", orfs.len());
                    for (i, orf) in orfs.iter().take(40).enumerate() {
                        out.push_str(&format!(
                            "#{:<3} frame {} {:<7} {:>5}..{:<5} {:>5} aa  {}\n",
                            i + 1,
                            orf.frame,
                            format!("{:?}", orf.strand),
                            orf.span.start,
                            orf.span.end,
                            orf.protein_len(),
                            if orf.has_stop { "stop" } else { "no-stop" },
                        ));
                    }
                    p.result = out;
                }
                Err(e) => p.error = Some(e.to_string()),
            }
        }
        (Ok(_), Ok(_)) => {
            p.error = Some("ORF finding needs a nucleotide sequence".into())
        }
        (Err(e), _) => p.error = Some(e),
        (_, Err(e)) => p.error = Some(e.to_string()),
    }
}

fn draw_restriction(p: &mut SequencePanel, ui: &mut egui::Ui) {
    ui.horizontal(|ui| {
        ui.label("Enzymes:");
        ui.text_edit_singleline(&mut p.enzymes);
    });
    ui.label(
        egui::RichText::new("e.g. EcoRI, BamHI, HindIII, NotI, XhoI")
            .weak()
            .small(),
    );
    if common::run_button(ui, "Digest + virtual gel") {
        run_restriction(p);
    }
}

/// Run the restriction digest + virtual gel — extracted for the
/// headless UI tests.
fn run_restriction(p: &mut SequencePanel) {
    p.error = None;
    match parse_seq(p) {
        Ok(s) if s.kind() == SeqKind::Dna => {
                let names: Vec<String> = p
                    .enzymes
                    .split(',')
                    .map(|n| n.trim().to_string())
                    .filter(|n| !n.is_empty())
                    .collect();
                let mut enzymes: Vec<&Enzyme> = Vec::new();
                let mut unknown: Vec<String> = Vec::new();
                for n in &names {
                    match enzyme_by_name(n) {
                        Some(e) => enzymes.push(e),
                        None => unknown.push(n.clone()),
                    }
                }
                if enzymes.is_empty() {
                    p.error = Some(format!("no known enzymes ({} unrecognised)", unknown.len()));
                    return;
                }
                match (restriction_map(&s, &enzymes), digest_to_gel(&s, &enzymes)) {
                    (Ok(map), Ok(gel)) => {
                        let mut out = String::new();
                        if !unknown.is_empty() {
                            out.push_str(&format!("(unrecognised: {})\n", unknown.join(", ")));
                        }
                        out.push_str(&format!("{} cut sites:\n", map.cut_sites.len()));
                        for c in map.cut_sites.iter().take(40) {
                            out.push_str(&format!(
                                "  {:<10} site@{:<6} cut@{}\n",
                                c.enzyme, c.site_start, c.top_cut_pos,
                            ));
                        }
                        out.push_str(&format!(
                            "\nVirtual gel — {} fragments, {} distinct bands:\n",
                            gel.band_sizes.len(),
                            gel.distinct_bands(),
                        ));
                        let mut bands = gel.band_sizes.clone();
                        bands.sort_unstable_by(|a, b| b.cmp(a));
                        for bp in &bands {
                            let bar = "#".repeat((*bp as f64).log10().max(0.0) as usize * 4 + 1);
                            out.push_str(&format!("  {bp:>7} bp  {bar}\n"));
                        }
                        p.result = out;
                    }
                    (Err(e), _) | (_, Err(e)) => p.error = Some(e.to_string()),
                }
            }
            Ok(_) => p.error = Some("restriction digest needs a DNA sequence".into()),
            Err(e) => p.error = Some(e),
        }
}

fn draw_primers(p: &mut SequencePanel, ui: &mut egui::Ui) {
    ui.horizontal(|ui| {
        ui.label("Target window:");
        ui.add(egui::DragValue::new(&mut p.primer_start).prefix("start "));
        ui.add(egui::DragValue::new(&mut p.primer_end).prefix("end "));
    });
    ui.horizontal(|ui| {
        ui.label("Target Tm (°C):");
        ui.add(egui::DragValue::new(&mut p.primer_target_tm).range(40.0..=80.0).speed(0.5));
    });
    if common::run_button(ui, "Design primer pair") {
        run_primers(p);
    }
}

/// Run the primer-pair design — extracted for the headless UI tests.
fn run_primers(p: &mut SequencePanel) {
    p.error = None;
    match parse_seq(p) {
        Ok(s) if s.kind() == SeqKind::Dna => {
            let constraints = PrimerConstraints {
                target_tm: p.primer_target_tm,
                ..PrimerConstraints::default()
            };
            let end = p.primer_end.min(s.len());
            match design_primers(&s, p.primer_start, end, constraints) {
                Ok(pair) => {
                    p.result = format!(
                        "Forward : {}\n  binding {}..{}  Tm {:.1} °C  GC {:.0} %{}\n\
                         Reverse : {}\n  binding {}..{}  Tm {:.1} °C  GC {:.0} %{}\n\n\
                         Amplicon size : {} bp\nΔTm           : {:.2} °C",
                        pair.forward.seq.as_str(),
                        pair.forward.binding_start,
                        pair.forward.binding_end,
                        pair.forward.tm,
                        pair.forward.gc * 100.0,
                        if pair.forward.is_clean() { "  clean" } else { "  ⚠ dimer/hairpin" },
                        pair.reverse.seq.as_str(),
                        pair.reverse.binding_start,
                        pair.reverse.binding_end,
                        pair.reverse.tm,
                        pair.reverse.gc * 100.0,
                        if pair.reverse.is_clean() { "  clean" } else { "  ⚠ dimer/hairpin" },
                        pair.amplicon_size(),
                        pair.tm_difference(),
                    );
                }
                Err(e) => p.error = Some(e.to_string()),
            }
        }
        Ok(_) => p.error = Some("primer design needs a DNA template".into()),
        Err(e) => p.error = Some(e),
    }
}

fn draw_pcr(p: &mut SequencePanel, ui: &mut egui::Ui) {
    ui.label("Forward primer (5'→3'):");
    ui.add(
        egui::TextEdit::singleline(&mut p.pcr_fwd)
            .font(egui::TextStyle::Monospace)
            .desired_width(f32::INFINITY),
    );
    ui.label("Reverse primer (5'→3'):");
    ui.add(
        egui::TextEdit::singleline(&mut p.pcr_rev)
            .font(egui::TextStyle::Monospace)
            .desired_width(f32::INFINITY),
    );
    if common::run_button(ui, "Simulate PCR") {
        run_pcr(p);
    }
}

/// Run the in-silico PCR — extracted for the headless UI tests.
fn run_pcr(p: &mut SequencePanel) {
    p.error = None;
    let fwd = common::clean_sequence(&p.pcr_fwd);
    let rev = common::clean_sequence(&p.pcr_rev);
    let parsed = (
        parse_seq(p),
        Seq::new(SeqKind::Dna, fwd.as_bytes()),
        Seq::new(SeqKind::Dna, rev.as_bytes()),
    );
    match parsed {
        (_, Err(e), _) | (_, _, Err(e)) => p.error = Some(format!("primer: {e}")),
        (Err(e), _, _) => p.error = Some(e),
        (Ok(tpl), _, _) if tpl.kind() != SeqKind::Dna => {
            p.error = Some("in-silico PCR needs a DNA template".into())
        }
        (Ok(tpl), Ok(f), Ok(r)) => {
            match simulate_pcr(&tpl, &f, &r, PcrOptions::default()) {
                Ok(amplicons) => {
                    let mut out = format!("{} amplicon(s):\n", amplicons.len());
                    for (i, a) in amplicons.iter().enumerate() {
                        out.push_str(&format!(
                            "#{:<3} {:>5}..{:<5}  {} bp{}\n",
                            i + 1,
                            a.start,
                            a.end,
                            a.size,
                            if a.wraps_origin { "  (wraps origin)" } else { "" },
                        ));
                    }
                    if let Some(first) = amplicons.first() {
                        let prod = first.product.as_str();
                        let preview: String = prod.chars().take(120).collect();
                        out.push_str(&format!("\nFirst product (preview):\n{preview}"));
                        if prod.len() > 120 {
                            out.push_str(" …");
                        }
                    }
                    p.result = out;
                }
                Err(e) => p.error = Some(e.to_string()),
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_sequence_parses() {
        let p = SequencePanel::default();
        let s = parse_seq(&p).expect("default seq must parse");
        assert!(s.len() > 30);
        assert_eq!(s.kind(), SeqKind::Dna);
    }

    #[test]
    fn empty_text_is_rejected() {
        let p = SequencePanel {
            seq_text: "  \n>just a header\n".to_string(),
            ..SequencePanel::default()
        };
        assert!(parse_seq(&p).is_err());
    }

    #[test]
    fn translate_default_table_produces_protein() {
        let p = SequencePanel::default();
        let s = parse_seq(&p).unwrap();
        let code = GeneticCode::by_id(1).unwrap();
        let prot = translate(&s, &code, TranslateOptions::default()).unwrap();
        assert_eq!(prot.kind(), SeqKind::Protein);
        assert!(!prot.as_str().is_empty());
    }
}

/// Headless egui UI-logic tests for the Sequence panel.
///
/// These render the panel into a windowless [`egui::Context`] across
/// representative states and exercise every extracted Run action — they
/// never open an OS window and never reach `rfd::FileDialog` (the
/// Load-FASTA button is a layout-only element that is drawn but never
/// clicked).
#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;
    use crate::genetics_workbench::GeneticsPanel;
    use crate::ValenxApp;

    /// Draw the Sequence panel once in a headless context. Panics
    /// inside the closure surface as a failed test.
    fn draw_headless(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            egui::CentralPanel::default().show(ctx, |ui| {
                super::draw(app, ui);
            });
        });
    }

    /// A fresh app with the Sequence panel active.
    fn app_with_panel() -> ValenxApp {
        let mut app = ValenxApp::default();
        app.genetics.active = GeneticsPanel::Sequence;
        app
    }

    #[test]
    fn draws_every_tool_without_panic() {
        // Fresh / valid-populated state: draw the panel with each
        // sub-tool selected. The draw fn must never panic.
        for tool in [
            Tool::Analyze,
            Tool::Translate,
            Tool::Orf,
            Tool::Restriction,
            Tool::Primers,
            Tool::Pcr,
        ] {
            let mut app = app_with_panel();
            app.genetics.sequence.tool = tool;
            draw_headless(&mut app);
        }
    }

    #[test]
    fn draws_post_run_and_error_states_without_panic() {
        // Post-run: a populated result string must render.
        let mut app = app_with_panel();
        app.genetics.sequence.result = "Length : 64 residues\n".to_string();
        draw_headless(&mut app);
        // Error state: a surfaced error must render.
        let mut app = app_with_panel();
        app.genetics.sequence.error = Some("sequence is empty".to_string());
        draw_headless(&mut app);
        // Empty input drawn (the default text cleared).
        let mut app = app_with_panel();
        app.genetics.sequence.seq_text.clear();
        draw_headless(&mut app);
    }

    #[test]
    fn run_analyze_produces_a_sane_result() {
        // Valid DNA input → composition analysis calls the real
        // valenx-bioseq API and produces a correctly-formatted result.
        let mut p = SequencePanel::default();
        run_analyze(&mut p);
        assert!(p.error.is_none(), "analyze errored: {:?}", p.error);
        assert!(p.result.contains("Length"));
        assert!(p.result.contains("GC content"));
        assert!(p.result.contains("Residue counts"));
    }

    #[test]
    fn run_translate_produces_a_protein() {
        let mut p = SequencePanel::default();
        run_translate(&mut p);
        assert!(p.error.is_none(), "translate errored: {:?}", p.error);
        assert!(p.result.contains("Protein"));
        assert!(p.result.contains("Genetic code"));
    }

    #[test]
    fn run_orf_finds_open_reading_frames() {
        let mut p = SequencePanel::default();
        run_orf(&mut p);
        assert!(p.error.is_none(), "ORF finder errored: {:?}", p.error);
        assert!(p.result.contains("ORFs"));
    }

    #[test]
    fn run_restriction_produces_a_virtual_gel() {
        let mut p = SequencePanel::default();
        run_restriction(&mut p);
        assert!(p.error.is_none(), "restriction errored: {:?}", p.error);
        assert!(p.result.contains("cut sites"));
        assert!(p.result.contains("Virtual gel"));
    }

    #[test]
    fn run_primers_designs_a_pair() {
        let mut p = SequencePanel::default();
        run_primers(&mut p);
        assert!(p.error.is_none(), "primer design errored: {:?}", p.error);
        assert!(p.result.contains("Forward"));
        assert!(p.result.contains("Reverse"));
        assert!(p.result.contains("Amplicon size"));
    }


    #[test]
    fn run_pcr_amplifies_the_template() {
        let mut p = SequencePanel::default();
        run_pcr(&mut p);
        assert!(p.error.is_none(), "PCR errored: {:?}", p.error);
        assert!(p.result.contains("amplicon"));
    }

    #[test]
    fn run_actions_surface_errors_on_empty_input() {
        // Empty input → every Run action sets an error rather than
        // panicking or producing a result.
        let blank = || SequencePanel {
            seq_text: String::new(),
            ..SequencePanel::default()
        };
        let mut p = blank();
        run_analyze(&mut p);
        assert!(p.error.is_some(), "analyze should error on empty input");
        let mut p = blank();
        run_translate(&mut p);
        assert!(p.error.is_some(), "translate should error on empty input");
        let mut p = blank();
        run_orf(&mut p);
        assert!(p.error.is_some(), "ORF finder should error on empty input");
        let mut p = blank();
        run_restriction(&mut p);
        assert!(p.error.is_some(), "restriction should error on empty input");
        let mut p = blank();
        run_primers(&mut p);
        assert!(p.error.is_some(), "primer design should error on empty input");
        let mut p = blank();
        run_pcr(&mut p);
        assert!(p.error.is_some(), "PCR should error on empty input");
    }

    #[test]
    fn run_actions_reject_a_protein_for_nucleotide_tools() {
        // A protein sequence is malformed input for translate / ORF /
        // restriction — they must surface an error, not panic.
        let protein = || SequencePanel {
            seq_text: "MKVLAGDRENT".to_string(),
            kind: SeqKind::Protein,
            ..SequencePanel::default()
        };
        let mut p = protein();
        run_translate(&mut p);
        assert!(p.error.is_some());
        let mut p = protein();
        run_orf(&mut p);
        assert!(p.error.is_some());
        let mut p = protein();
        run_restriction(&mut p);
        assert!(p.error.is_some());
    }
}
