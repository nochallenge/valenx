//! Per-panel contextual help text. F1 opens an [`egui::Window`]
//! whose body is the entry from [`help_for_panel`] for the
//! currently-active panel.
//!
//! The text is intentionally short (≤ ~10 lines per panel) so the
//! help popup is glance-readable rather than a dump of the user
//! manual. Format is:
//!
//! ```text
//! # PanelName
//!
//! What it does — one paragraph
//!
//! ## Key controls
//! - **Input X** — what it means
//! - **Button Y** — what it does
//!
//! ## Quick start
//! 1. ...
//! 2. ...
//! ```
//!
//! Rendered with [`render_help_window`] which parses the markdown
//! enough to render headings + bullets — no full markdown impl,
//! just `#`, `##`, `-`, blank-line paragraphs.

use eframe::egui;

/// Render the contextual-help window for the given panel name. The
/// `open` parameter follows egui's window-open convention — set to
/// `false` to hide. Returns `true` if any markdown-ish content was
/// rendered (always — the catalogue has a generic fallback).
pub fn render_help_window(ctx: &egui::Context, open: &mut bool, panel_name: &str) -> bool {
    let body = help_for_panel(panel_name);
    let title = format!("Help — {panel_name}");
    let mut rendered = false;
    egui::Window::new(title)
        .open(open)
        .collapsible(false)
        .resizable(true)
        .default_width(480.0)
        .default_height(360.0)
        .show(ctx, |ui| {
            egui::ScrollArea::vertical().show(ui, |ui| {
                for line in body.lines() {
                    let line = line.trim_end();
                    if let Some(rest) = line.strip_prefix("# ") {
                        ui.heading(rest);
                    } else if let Some(rest) = line.strip_prefix("## ") {
                        ui.add_space(4.0);
                        ui.label(egui::RichText::new(rest).strong());
                    } else if let Some(rest) = line.strip_prefix("- ") {
                        // Bullet line — render with a leading dot
                        // and italicise the bolded segment if any.
                        ui.horizontal(|ui| {
                            ui.label("•");
                            ui.label(rest);
                        });
                    } else if line.is_empty() {
                        ui.add_space(2.0);
                    } else {
                        ui.label(line);
                    }
                }
                rendered = true;
                ui.add_space(8.0);
                ui.separator();
                ui.weak("F1 to reopen this help. Ctrl+P to search every action.");
            });
        });
    rendered
}

/// One-line summary of the panel — shown as the tab tooltip when
/// the user hovers a workbench-selector chip.
pub fn short_summary(panel: &str) -> &'static str {
    match panel {
        // Genetics panels.
        "Sequence" => "Edit DNA/RNA/protein, translate, find ORFs, run primer / PCR.",
        "Alignment" => "Pairwise / multiple alignment, k-mer search (valenx-align).",
        "Phylogenetics" => "Build distance trees, render Newick / ASCII (valenx-phylo).",
        "Population Genetics" => "Wright-Fisher / coalescent, summary stats (valenx-popgen).",
        "RNA Structure" => "Fold RNA into secondary structure (valenx-rnastruct).",
        "RNA Designer" => "Guided synthetic-RNA design wizard (valenx-rnadesign).",
        "Molecular Dynamics" => "Classical MD simulation (valenx-md).",
        "Cheminformatics" => "Descriptors + fingerprints + similarity (valenx-cheminf).",
        "Macromolecular Structure" => "Ramachandran / superposition / contacts (valenx-biostruct).",
        "Quantum Chemistry" => "Hartree-Fock SCF on small molecules (valenx-qchem).",
        "Genomics" => "VCF parsing, CRISPR target search, read simulation (valenx-genomics).",
        "Systems Biology" => "ODE reaction networks (valenx-sysbio).",
        "Docking" => "Ligand docking + virtual screening (valenx-dock-screen).",
        "Gene Editing" => "CRISPR / base / prime / mRNA design (valenx-genediting).",
        // Aero workbench sections.
        "Body" => "Pick the body source (demo box, primitive, loaded mesh).",
        "Wind" => "Set freestream speed, AoA, density, viscosity.",
        "Ground & Wheels" => "Toggle moving-ground + rotating-wheel models.",
        "Tunnel & Mesh" => "Tunnel dimensions + cut-cell mesh resolution.",
        "Solver" => "k-ε / k-ω SST, iteration cap, residual target.",
        "Run" => "Steady solve or angle-sweep polar.",
        "Results" => "Drag / lift / moment coefficients, residual plot, AoA polar.",
        "Visualization" => "Push Cp / velocity / pressure / Q-criterion fields to viewport.",
        // Mesh toolbox (a few representative).
        "Part" => "CAD primitive creation + boolean ops.",
        "Draft" => "2D drafting (lines, circles, dimensions).",
        "TechDraw" => "Generate 2D engineering drawings from solids.",
        "Assembly" => "Constraints + parts list.",
        "Surface" => "NURBS + BRep surface tooling.",
        "CAM" => "Toolpath generation + post-processing.",
        "Arch" => "AEC / BIM architectural elements.",
        "Spreadsheet" => "Parameter table that drives the model.",
        "Dock" => "Layout manager for floating panels.",
        "Sketcher" => "2D sketch constraint solver.",
        "Part Design" => "Sketch-driven feature modelling.",
        _ => "(no help text available for this panel yet — F1 again to dismiss)",
    }
}

/// Full help body for the panel. Markdown-ish — `#`, `##`, `-` lines
/// + blank-line paragraphs.
///
/// Brevity matters here: this text fills a popover, not a doc page.
/// Aim for ~10 lines including the heading and the quick-start
/// recipe. Panels without bespoke text fall through to a generic
/// "use Ctrl+P to find actions" hint.
pub fn help_for_panel(panel: &str) -> &'static str {
    match panel {
        "Sequence" => "\
# Sequence — valenx-bioseq

Edit a DNA / RNA / protein sequence and run six sub-tools:
analyze composition, translate, find ORFs (6 frames), digest with
restriction enzymes, design primers, simulate PCR.

## Key controls
- **Residues** — multiline input (FASTA headers are stripped).
- **Kind radio** — DNA / RNA / Protein.
- **Tool selector** — picks which sub-tool the Run button drives.
- **Run button** — executes the active sub-tool (Ctrl+R).

## Quick start
1. Paste a DNA sequence or click Load FASTA…
2. Pick the Tool you want (e.g. ORF finder).
3. Click Run — results appear in the Result box.
4. Ctrl+Z to undo a paste; Ctrl+R to re-run.
",
        "Alignment" => "\
# Alignment — valenx-align

Pairwise alignment (global / local / semi-global), multiple-sequence
alignment, and k-mer seed search.

## Key controls
- **Sequence A / B** — the two input sequences.
- **Mode radio** — Global (Needleman-Wunsch) / Local (Smith-Waterman) / Semi-global.
- **Match / Mismatch / Gap** — substitution matrix tweaks.
- **Run** — executes the selected mode (Ctrl+R).

## Quick start
1. Paste two sequences into Sequence A / Sequence B.
2. Choose Mode — Local for finding a shared motif, Global for whole-sequence.
3. Click Run; the alignment renders in the Result box.
",
        "RNA Structure" => "\
# RNA Structure — valenx-rnastruct

Fold an RNA sequence into a minimum-free-energy secondary structure;
output is dot-bracket notation + the ΔG (kcal/mol).

## Key controls
- **Sequence** — RNA in IUPAC letters.
- **Algorithm** — Nussinov (max-pairs) or Zuker (MFE).
- **Temperature** — folding temperature (°C).
- **Run** — folds the sequence (Ctrl+R).

## Quick start
1. Paste an RNA sequence into the input box.
2. Pick the Zuker algorithm for MFE-based folding.
3. Click Run — the dot-bracket structure appears below.
",
        "RNA Designer" => "\
# RNA Designer — valenx-rnadesign

A guided 6-step wizard for designing synthetic RNAs end-to-end:
target structure → codon optimization → folding check → restriction-site
masking → primer design → assembly file.

## Navigation
- **Next / Back** — step through the wizard (Ctrl+R = next step).
- **Reset** — clear the wizard state (after confirmation).

## Quick start
1. Step 1 — paste the target dot-bracket structure.
2. Step 2 — pick the codon-optimization target organism.
3. Step 3 → 6 — accept the defaults or tune as needed.
",
        "Molecular Dynamics" => "\
# Molecular Dynamics — valenx-md

Classical MD on a small system (Lennard-Jones + harmonic bonds).
Outputs a trajectory the Visualization tab can render.

## Key controls
- **Timestep (ps)** — integration step.
- **Total time (ns)** — total simulated time.
- **Temperature (K)** — thermostat setpoint.
- **Run** — integrates the trajectory (Ctrl+R).
",
        "Cheminformatics" => "\
# Cheminformatics — valenx-cheminf

Compute molecular descriptors (LogP, TPSA, MW), Morgan fingerprints,
and Tanimoto similarity between two molecules.

## Key controls
- **SMILES** — input molecule(s) in SMILES notation.
- **Tool** — descriptors / fingerprint / similarity.
- **Run** — executes the selected tool (Ctrl+R).
",
        "Genomics" => "\
# Genomics — valenx-genomics

Parse VCF variant calls, design CRISPR guides, simulate short reads.

## Key controls
- **VCF text** — paste a VCF blob (header optional).
- **Tool** — summary / CRISPR design / read simulator.
- **Run** — runs the selected tool (Ctrl+R).
",
        "Body" => "\
# Body — Wind Tunnel

Pick the body the solver puts in the tunnel: a built-in demo box,
a CAD primitive, or the currently-loaded mesh.

## Tip
The body's projected frontal area is auto-computed from the source —
it's the reference area for Cd / Cl coefficients.
",
        "Wind" => "\
# Wind conditions — Wind Tunnel

Sets the freestream: speed (m/s), angle of attack (deg), air density
(kg/m³), and dynamic viscosity (Pa·s). Reynolds number is derived
from speed × chord ÷ ν and displayed read-only.

## Quick start
- Road-car defaults: 30 m/s, AoA 0, sea-level air.
- Aircraft defaults: bump speed to 60 m/s; pick a chord-scale body.
",
        "Run" => "\
# Run — Wind Tunnel

Steady solve runs once; Angle sweep loops through a range of AoA and
reports a Cd/Cl polar.

## Controls
- **Mode** — Steady / Sweep.
- **Sweep range** — AoA start / stop / step (deg).
- **Run / Cancel** — Ctrl+R to start, Esc to abort.
",
        "Results" => "\
# Results — Wind Tunnel

Reports the converged Cd / Cl / Cm coefficients, residual history,
and (for sweep runs) the AoA polar plot. Push fields into the
3-D viewport from the Visualization section.
",
        "Visualization" => "\
# Flow visualization — Wind Tunnel

Push a flow field into the 3-D viewport's colour overlay: pressure
coefficient, velocity magnitude, static pressure, or Q-criterion
vortex-marker. Use **Clear overlay** to return to a plain wireframe.
",
        _ => "\
# Panel help

This panel doesn't have a dedicated help entry yet. Try:

- **Hover any control** — every interactive control has a tooltip.
- **Ctrl+P** — open the command palette and fuzzy-search every action.
- **?** — toggle the keyboard-shortcut cheat-sheet.
",
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn short_summary_returns_non_empty_for_known_panels() {
        for p in [
            "Sequence",
            "Alignment",
            "RNA Structure",
            "Body",
            "Wind",
            "Run",
        ] {
            assert!(!short_summary(p).is_empty(), "{p} has empty summary");
        }
    }

    #[test]
    fn short_summary_falls_back_gracefully() {
        let s = short_summary("DefinitelyNotARealPanel");
        assert!(!s.is_empty());
        assert!(s.contains("no help"));
    }

    #[test]
    fn help_for_panel_starts_with_heading() {
        // Every catalogue entry must lead with a `# Heading` so the
        // popup renders a clear title block.
        for p in [
            "Sequence",
            "Alignment",
            "RNA Structure",
            "Molecular Dynamics",
            "Body",
            "Wind",
            "Run",
            "Results",
            "Visualization",
        ] {
            let body = help_for_panel(p);
            assert!(
                body.trim_start().starts_with("# "),
                "{p} help does not start with a # heading: {:?}",
                body.lines().next()
            );
        }
    }

    #[test]
    fn fallback_help_text_offers_a_path_forward() {
        let body = help_for_panel("UnregisteredPanel");
        assert!(body.contains("Ctrl+P"));
    }

    #[test]
    fn render_help_window_runs_without_panic() {
        // Headless: render the popup in a windowless ctx for every
        // catalogue entry — must not panic.
        let ctx = egui::Context::default();
        let mut open = true;
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            for p in [
                "Sequence",
                "Alignment",
                "Wind",
                "Run",
                "UnregisteredPanel",
            ] {
                render_help_window(ctx, &mut open, p);
            }
        });
    }
}
