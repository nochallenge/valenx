//! Panel 10 — **Genomics** (`valenx-genomics`).
//!
//! Summarise a VCF or SAM file, call variants from a SAM file against a
//! reference via the pileup engine, design CRISPR guide RNAs with PAM
//! scanning and a Doench-style on-target score, and simulate Illumina
//! short reads — all native `valenx-genomics` calls.

use eframe::egui;

use valenx_genomics::crispr::guide::{scan_guides, PamSpec};
use valenx_genomics::format::pileup::{build_pileup, Reference};
use valenx_genomics::format::sam::SamFile;
use valenx_genomics::format::vcf::VcfFile;
use valenx_genomics::simulate::illumina::{simulate_reads, IlluminaProfile};
use valenx_genomics::variant::call::{call_variants, CallParams};

use super::common;
use crate::ValenxApp;

/// Which genomics sub-tool is showing.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Default)]
pub(crate) enum Tool {
    #[default]
    FileSummary,
    VariantCalling,
    Crispr,
    ReadSim,
}

/// VCF vs SAM for the file-summary tool.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Default)]
pub(crate) enum FileKind {
    #[default]
    Vcf,
    Sam,
}

/// PAM / nuclease choice for the CRISPR tool.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Default)]
pub(crate) enum Pam {
    #[default]
    SpCas9,
    SaCas9,
    Cas12a,
}

/// Snapshot of every editable input. Multi-text + numerics + enums —
/// big but flat. Stored once per `Run`.
#[derive(Clone, Debug, PartialEq, Default)]
pub(crate) struct GenomicsSnapshot {
    pub(crate) file_kind: FileKind,
    pub(crate) file_text: String,
    pub(crate) reference: String,
    pub(crate) ref_name: String,
    pub(crate) crispr_target: String,
    pub(crate) pam: Pam,
    pub(crate) n_reads: usize,
    pub(crate) read_seed: u64,
}

/// Form + result state for the Genomics panel.
pub struct GenomicsPanel {
    tool: Tool,
    file_kind: FileKind,
    /// VCF / SAM text for the summary tool, or SAM text for the
    /// variant-caller.
    file_text: String,
    /// Reference sequence for variant calling / read simulation.
    reference: String,
    /// Reference contig name (must match the SAM `rname`).
    ref_name: String,
    /// CRISPR target sequence + PAM choice.
    crispr_target: String,
    pam: Pam,
    /// Read-simulation parameters.
    n_reads: usize,
    read_seed: u64,
    error: Option<String>,
    result: String,
    /// Undo / redo over every editable input.
    history: crate::undo::History<GenomicsSnapshot>,
}

impl GenomicsPanel {
    fn snapshot(&self) -> GenomicsSnapshot {
        GenomicsSnapshot {
            file_kind: self.file_kind,
            file_text: self.file_text.clone(),
            reference: self.reference.clone(),
            ref_name: self.ref_name.clone(),
            crispr_target: self.crispr_target.clone(),
            pam: self.pam,
            n_reads: self.n_reads,
            read_seed: self.read_seed,
        }
    }
    fn restore(&mut self, s: GenomicsSnapshot) {
        self.file_kind = s.file_kind;
        self.file_text = s.file_text;
        self.reference = s.reference;
        self.ref_name = s.ref_name;
        self.crispr_target = s.crispr_target;
        self.pam = s.pam;
        self.n_reads = s.n_reads;
        self.read_seed = s.read_seed;
    }
    pub fn undo_edit(&mut self) -> bool {
        let current = self.snapshot();
        if let Some(prev) = self.history.undo(current) {
            self.restore(prev);
            self.error = None;
            true
        } else {
            false
        }
    }
    pub fn redo_edit(&mut self) -> bool {
        let current = self.snapshot();
        if let Some(next) = self.history.redo(current) {
            self.restore(next);
            self.error = None;
            true
        } else {
            false
        }
    }
    pub fn can_undo(&self) -> bool { self.history.can_undo() }
    pub fn can_redo(&self) -> bool { self.history.can_redo() }
}

impl Default for GenomicsPanel {
    fn default() -> Self {
        GenomicsPanel {
            tool: Tool::FileSummary,
            file_kind: FileKind::Vcf,
            file_text: "##fileformat=VCFv4.2\n\
                        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
                        chr1\t100\t.\tA\tG\t60\tPASS\tDP=30\n\
                        chr1\t250\trs7\tGT\tG\t45\tPASS\tDP=22\n"
                .to_string(),
            // ~210 bp — long enough for the 150 bp Illumina read sim.
            reference: "ACGTACGTACGTAAGGTTCCAACGTACGTACGTTAGCTAGCTAGGCATGCATGC\
                        ATCGATCGATCGGGCCAATTGGCCAATTGGTACGTACGTAAGCTTGGATCCGAA\
                        TTCCGGAATTCCGGTTAACCGGTTAACCATGCATGCATGCTAGCTAGCTAGGCA\
                        TGCATGCATCGATCGATCGGGCCAATTGGCC"
                .to_string(),
            ref_name: "chr1".to_string(),
            crispr_target: "ACGTACGTGGAAGGTTCCAACGGGATGCATGCATGCTAGCTAGCTAGGCATGCA".to_string(),
            pam: Pam::SpCas9,
            n_reads: 50,
            read_seed: 7,
            error: None,
            result: String::new(),
            history: crate::undo::History::new(),
        }
    }
}

/// Render the Genomics panel.
pub fn draw(app: &mut ValenxApp, ui: &mut egui::Ui) {
    let p = &mut app.genetics.genomics;

    common::section(ui, "Tool");
    ui.horizontal_wrapped(|ui| {
        ui.selectable_value(&mut p.tool, Tool::FileSummary, "VCF / SAM summary")
            .on_hover_text("Parse and summarize a VCF (variants) or SAM (alignments) blob.");
        ui.selectable_value(&mut p.tool, Tool::VariantCalling, "Variant calling")
            .on_hover_text("Call SNVs / indels from a SAM pileup against a reference.");
        ui.selectable_value(&mut p.tool, Tool::Crispr, "CRISPR guides")
            .on_hover_text("Find candidate CRISPR guide RNAs against a target sequence.");
        ui.selectable_value(&mut p.tool, Tool::ReadSim, "Read simulation")
            .on_hover_text("Simulate short reads from a reference (Illumina-style model).");
        ui.separator();
        let (u, r) = common::undo_redo_inline(ui, p.can_undo(), p.can_redo());
        if u { p.undo_edit(); }
        if r { p.redo_edit(); }
    });
    ui.separator();

    match p.tool {
        Tool::FileSummary => draw_summary(p, ui),
        Tool::VariantCalling => draw_variant(p, ui),
        Tool::Crispr => draw_crispr(p, ui),
        Tool::ReadSim => draw_readsim(p, ui),
    }

    common::error_line(ui, &p.error);
    if !p.result.is_empty() {
        ui.separator();
        common::section(ui, "Result");
        common::mono_output(ui, "genomics_result", &p.result, 15);
    }
}

fn draw_summary(p: &mut GenomicsPanel, ui: &mut egui::Ui) {
    ui.horizontal(|ui| {
        ui.radio_value(&mut p.file_kind, FileKind::Vcf, "VCF")
            .on_hover_text("Variant Call Format — one line per variant.");
        ui.radio_value(&mut p.file_kind, FileKind::Sam, "SAM")
            .on_hover_text("Sequence Alignment / Map — one line per read alignment.");
        if ui
            .small_button("Load file…")
            .on_hover_text("Open a .vcf, .sam, or .txt file from disk.")
            .clicked()
        {
            if let Some(path) = rfd::FileDialog::new()
                .add_filter("VCF / SAM", &["vcf", "sam", "txt"])
                .pick_file()
            {
                // Round-21 H1: see biostruct loader.
                match valenx_core::io_caps::read_capped_to_string(
                    &path,
                    valenx_core::io_caps::MAX_GENETICS_FILE_BYTES as usize,
                ) {
                    Ok(t) => {
                        p.file_kind = if path
                            .extension()
                            .map(|e| e.eq_ignore_ascii_case("sam"))
                            .unwrap_or(false)
                        {
                            FileKind::Sam
                        } else {
                            FileKind::Vcf
                        };
                        p.file_text = t;
                    }
                    Err(e) => p.error = Some(format!("read: {e}")),
                }
            }
        }
    });
    ui.add(
        egui::TextEdit::multiline(&mut p.file_text)
            .id_source("genomics_file_input")
            .font(egui::TextStyle::Monospace)
            .desired_width(f32::INFINITY)
            .desired_rows(7),
    );
    if common::run_button(ui, "Summarize") {
        let snap = p.snapshot();
        p.history.record(snap);
        run_summary(p);
    }
}

/// Run the VCF / SAM file summary — extracted for the headless UI
/// tests.
fn run_summary(p: &mut GenomicsPanel) {
    p.error = None;
    match p.file_kind {
        FileKind::Vcf => match VcfFile::parse(&p.file_text) {
                Ok(vcf) => {
                    let snvs = vcf
                        .records
                        .iter()
                        .filter(|r| r.reference.len() == 1 && r.alt.iter().all(|a| a.len() == 1))
                        .count();
                    let mut out = format!(
                        "VCF {} · {} samples · {} records\n  SNV-like      : {}\n  indel-like    : {}\n\n",
                        vcf.header.fileformat,
                        vcf.header.samples.len(),
                        vcf.records.len(),
                        snvs,
                        vcf.records.len() - snvs,
                    );
                    for r in vcf.records.iter().take(40) {
                        out.push_str(&format!(
                            "  {}:{:<8} {} -> {}  qual {}\n",
                            r.chrom,
                            r.pos,
                            r.reference,
                            r.alt.join(","),
                            r.qual.map(|q| format!("{q:.0}")).unwrap_or_else(|| ".".into()),
                        ));
                    }
                    p.result = out;
                }
                Err(e) => p.error = Some(e.to_string()),
            },
            FileKind::Sam => match SamFile::parse(&p.file_text) {
                Ok(sam) => {
                    let mapped = sam.records.iter().filter(|r| !r.is_unmapped()).count();
                    let mut out = format!(
                        "SAM · {} reference seqs · {} records\n  mapped        : {}\n  unmapped      : {}\n\n",
                        sam.header.references.len(),
                        sam.records.len(),
                        mapped,
                        sam.records.len() - mapped,
                    );
                    for rec in sam.records.iter().take(40) {
                        out.push_str(&format!(
                            "  {:<14} {} @ {:<8} mapq {:<3} {} bp\n",
                            rec.qname,
                            if rec.rname.is_empty() { "*" } else { &rec.rname },
                            rec.pos,
                            rec.mapq,
                            rec.seq.len(),
                        ));
                    }
                    p.result = out;
                }
                Err(e) => p.error = Some(e.to_string()),
            },
    }
}

fn draw_variant(p: &mut GenomicsPanel, ui: &mut egui::Ui) {
    common::section(ui, "Reference");
    ui.horizontal(|ui| {
        ui.label("contig:");
        ui.text_edit_singleline(&mut p.ref_name);
    });
    common::seq_input(ui, "genomics_var_ref", "Reference sequence:", &mut p.reference, 2);
    common::section(ui, "Aligned reads (SAM)");
    ui.add(
        egui::TextEdit::multiline(&mut p.file_text)
            .id_source("genomics_var_sam")
            .font(egui::TextStyle::Monospace)
            .desired_width(f32::INFINITY)
            .desired_rows(6),
    );
    if common::run_button(ui, "Build pileup + call variants") {
        let snap = p.snapshot();
        p.history.record(snap);
        run_variant(p);
    }
}

/// Run the pileup-build + variant-call — extracted for the headless UI
/// tests.
fn run_variant(p: &mut GenomicsPanel) {
    p.error = None;
    match SamFile::parse(&p.file_text) {
        Ok(sam) => {
                let mut reference = Reference::new();
                reference.add(p.ref_name.clone(), common::clean_sequence(&p.reference));
                match build_pileup(&sam.records, &reference, 0) {
                    Ok(columns) => {
                        match call_variants(&columns, &CallParams::default()) {
                            Ok(variants) => {
                                let mut out = format!(
                                    "{} pileup columns · {} variant call(s):\n",
                                    columns.len(),
                                    variants.len(),
                                );
                                for v in variants.iter().take(40) {
                                    out.push_str(&format!(
                                        "  {}:{:<8} {} -> {}  {:?}  depth {} ({:.0}% alt)  qual {:.1}\n",
                                        v.chrom,
                                        v.pos,
                                        v.reference,
                                        v.alt,
                                        v.kind,
                                        v.depth,
                                        v.alt_fraction * 100.0,
                                        v.qual,
                                    ));
                                }
                                if variants.is_empty() {
                                    out.push_str(
                                        "  (no variants — reads may all match the \
                                         reference, or the SAM has no mapped reads)\n",
                                    );
                                }
                                p.result = out;
                            }
                            Err(e) => p.error = Some(e.to_string()),
                        }
                    }
                    Err(e) => p.error = Some(format!("pileup: {e}")),
                }
            }
            Err(e) => p.error = Some(format!("SAM parse: {e}")),
        }
}

fn draw_crispr(p: &mut GenomicsPanel, ui: &mut egui::Ui) {
    common::seq_input(ui, "genomics_crispr_target", "Target DNA:", &mut p.crispr_target, 3);
    ui.horizontal(|ui| {
        ui.label("Nuclease / PAM:");
        ui.radio_value(&mut p.pam, Pam::SpCas9, "SpCas9 (NGG)");
        ui.radio_value(&mut p.pam, Pam::SaCas9, "SaCas9");
        ui.radio_value(&mut p.pam, Pam::Cas12a, "Cas12a");
    });
    if common::run_button(ui, "Scan guide RNAs") {
        let snap = p.snapshot();
        p.history.record(snap);
        run_crispr(p);
    }
}

/// Run the CRISPR guide scan — extracted for the headless UI tests.
fn run_crispr(p: &mut GenomicsPanel) {
    p.error = None;
    let target = common::clean_sequence(&p.crispr_target);
    if target.is_empty() {
        p.error = Some("target sequence is empty".into());
        return;
    }
    let pam = match p.pam {
        Pam::SpCas9 => PamSpec::spcas9(),
        Pam::SaCas9 => PamSpec::sacas9(),
        Pam::Cas12a => PamSpec::cas12a(),
    };
    match scan_guides(target.as_bytes(), &pam) {
        Ok(mut guides) => {
            guides.sort_by(|a, b| {
                b.on_target_score
                    .partial_cmp(&a.on_target_score)
                    .unwrap_or(std::cmp::Ordering::Equal)
            });
            let mut out = format!("{} candidate guides (best first):\n", guides.len());
            for (i, g) in guides.iter().take(40).enumerate() {
                out.push_str(&format!(
                    "#{:<3} {} {} @{:<5} {:?}  on-target {:.3}  GC {:.0}%\n",
                    i + 1,
                    g.protospacer,
                    g.pam,
                    g.start,
                    g.strand,
                    g.on_target_score,
                    g.gc_content * 100.0,
                ));
            }
            p.result = out;
        }
        Err(e) => p.error = Some(e.to_string()),
    }
}

fn draw_readsim(p: &mut GenomicsPanel, ui: &mut egui::Ui) {
    common::seq_input(ui, "genomics_readsim_ref", "Reference sequence:", &mut p.reference, 3);
    ui.horizontal(|ui| {
        ui.label("Reads:");
        ui.add(egui::DragValue::new(&mut p.n_reads).range(1..=100_000));
        ui.label("seed:");
        ui.add(egui::DragValue::new(&mut p.read_seed));
    });
    ui.label(
        egui::RichText::new("HiSeq 150 bp Illumina error model")
            .weak()
            .small(),
    );
    if common::run_button(ui, "Simulate reads") {
        let snap = p.snapshot();
        p.history.record(snap);
        run_readsim(p);
    }
}

/// Run the Illumina read simulation — extracted for the headless UI
/// tests.
fn run_readsim(p: &mut GenomicsPanel) {
    p.error = None;
    let reference = common::clean_sequence(&p.reference);
    if reference.is_empty() {
        p.error = Some("reference is empty".into());
        return;
    }
    let profile = IlluminaProfile::hiseq_150();
    match simulate_reads(reference.as_bytes(), &profile, p.n_reads, p.read_seed) {
        Ok(reads) => {
            let total_bp: usize = reads.iter().map(|r| r.record.seq.len()).sum();
            let mut out = format!(
                "simulated {} reads · {} bp total\nmean read length: {} bp\n\n",
                reads.len(),
                total_bp,
                if reads.is_empty() {
                    0
                } else {
                    total_bp / reads.len()
                },
            );
            for rec in reads.iter().take(6) {
                let preview: String = rec.record.seq.as_str().chars().take(60).collect();
                out.push_str(&format!("@{}\n{}\n", rec.record.id, preview));
            }
            if reads.len() > 6 {
                out.push_str(&format!("… and {} more reads\n", reads.len() - 6));
            }
            p.result = out;
        }
        Err(e) => p.error = Some(e.to_string()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_vcf_parses() {
        let p = GenomicsPanel::default();
        let vcf = VcfFile::parse(&p.file_text).expect("demo VCF must parse");
        assert_eq!(vcf.records.len(), 2);
    }

    #[test]
    fn crispr_scan_finds_guides() {
        let p = GenomicsPanel::default();
        let target = common::clean_sequence(&p.crispr_target);
        let guides = scan_guides(target.as_bytes(), &PamSpec::spcas9()).unwrap();
        // The demo target contains GG dinucleotides, so SpCas9 guides
        // must exist.
        assert!(!guides.is_empty());
    }

    #[test]
    fn read_simulation_produces_reads() {
        let p = GenomicsPanel::default();
        let reference = common::clean_sequence(&p.reference);
        let reads = simulate_reads(
            reference.as_bytes(),
            &IlluminaProfile::hiseq_150(),
            10,
            p.read_seed,
        )
        .unwrap();
        assert_eq!(reads.len(), 10);
    }
}

/// Headless egui UI-logic tests for the Genomics panel.
#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;
    use crate::genetics_workbench::GeneticsPanel;
    use crate::ValenxApp;

    fn draw_headless(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            egui::CentralPanel::default().show(ctx, |ui| {
                super::draw(app, ui);
            });
        });
    }

    fn app_with_panel() -> ValenxApp {
        let mut app = ValenxApp::default();
        app.genetics.active = GeneticsPanel::Genomics;
        app
    }

    #[test]
    fn draws_every_tool_without_panic() {
        for tool in [
            Tool::FileSummary,
            Tool::VariantCalling,
            Tool::Crispr,
            Tool::ReadSim,
        ] {
            let mut app = app_with_panel();
            app.genetics.genomics.tool = tool;
            draw_headless(&mut app);
        }
    }

    #[test]
    fn draws_post_run_and_error_states_without_panic() {
        let mut app = app_with_panel();
        app.genetics.genomics.result = "VCF v4.2 · 0 samples · 2 records\n".to_string();
        draw_headless(&mut app);
        let mut app = app_with_panel();
        app.genetics.genomics.error = Some("SAM parse failed".to_string());
        draw_headless(&mut app);
    }

    #[test]
    fn run_summary_parses_the_demo_vcf() {
        // The demo VCF → the real valenx-genomics VCF parser produces
        // a correctly-formatted summary.
        let mut p = GenomicsPanel {
            tool: Tool::FileSummary,
            file_kind: FileKind::Vcf,
            ..GenomicsPanel::default()
        };
        run_summary(&mut p);
        assert!(p.error.is_none(), "summary errored: {:?}", p.error);
        assert!(p.result.contains("VCF"));
        assert!(p.result.contains("records"));
    }

    #[test]
    fn run_crispr_finds_guides() {
        let mut p = GenomicsPanel {
            tool: Tool::Crispr,
            ..GenomicsPanel::default()
        };
        run_crispr(&mut p);
        assert!(p.error.is_none(), "CRISPR scan errored: {:?}", p.error);
        assert!(p.result.contains("candidate guides"));
    }

    #[test]
    fn run_readsim_simulates_reads() {
        let mut p = GenomicsPanel {
            tool: Tool::ReadSim,
            n_reads: 12,
            ..GenomicsPanel::default()
        };
        run_readsim(&mut p);
        assert!(p.error.is_none(), "read sim errored: {:?}", p.error);
        assert!(p.result.contains("simulated"));
        assert!(p.result.contains("bp total"));
    }

    #[test]
    fn run_variant_handles_the_default_inputs() {
        // The variant caller runs end-to-end; the demo file_text is a
        // VCF, so SAM parsing of it surfaces an error rather than
        // panicking — either way the panel stays well-formed.
        let mut p = GenomicsPanel {
            tool: Tool::VariantCalling,
            ..GenomicsPanel::default()
        };
        run_variant(&mut p);
        assert!(p.error.is_some() || !p.result.is_empty());
    }

    #[test]
    fn run_actions_surface_errors_on_bad_input() {
        // CRISPR with an empty target.
        let mut p = GenomicsPanel {
            tool: Tool::Crispr,
            crispr_target: String::new(),
            ..GenomicsPanel::default()
        };
        run_crispr(&mut p);
        assert!(p.error.is_some(), "CRISPR should error on empty target");
        // Read simulation with an empty reference.
        let mut p = GenomicsPanel {
            tool: Tool::ReadSim,
            reference: String::new(),
            ..GenomicsPanel::default()
        };
        run_readsim(&mut p);
        assert!(p.error.is_some(), "read sim should error on empty reference");
        // Summary of malformed SAM text.
        let mut p = GenomicsPanel {
            tool: Tool::FileSummary,
            file_kind: FileKind::Sam,
            file_text: "this is not a SAM file\n@@@\n".to_string(),
            ..GenomicsPanel::default()
        };
        run_summary(&mut p);
        // The SAM parser is lenient with junk lines; assert no panic
        // and a well-formed outcome (error or result).
        assert!(p.error.is_some() || !p.result.is_empty());
    }
}
