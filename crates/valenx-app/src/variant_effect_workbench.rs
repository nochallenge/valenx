//! Variant-effect workbench — HGVS variant parsing on `valenx-variant-effect`.
//!
//! A right-side panel that parses an HGVS-style substitution string (protein
//! `p.R273H` / `p.Arg273His`, or coding `c.817C>T`) into a structured
//! [`Variant`] and reports the parsed components — or a clear parse error.
//! The entry point to the variant-effect orchestration crate. Headless-testable.

use eframe::egui;

use valenx_variant_effect::{parse, Variant};

use crate::ValenxApp;

/// Persistent state for the variant-effect workbench.
pub struct VariantEffectWorkbenchState {
    input: String,
    parsed: Option<Result<Variant, String>>,
}

impl Default for VariantEffectWorkbenchState {
    fn default() -> Self {
        Self { input: "p.R273H".to_string(), parsed: None }
    }
}

/// A human-readable description of a parsed variant.
fn describe(v: &Variant) -> String {
    match v {
        Variant::ProteinSub { wt, pos, mt } => {
            format!("protein substitution {}{}{} — residue {pos}", wt.as_char(), pos, mt.as_char())
        }
        Variant::CodingSub { pos, wt, mt } => {
            format!("coding substitution c.{pos}{}>{}", wt.as_char(), mt.as_char())
        }
    }
}

/// Draw the variant-effect workbench (a no-op unless toggled on via
/// View → Variant Effect).
pub fn draw_variant_effect_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_variant_effect_workbench {
        return;
    }
    egui::SidePanel::right("valenx_variant_effect_workbench")
        .resizable(true)
        .default_width(340.0)
        .width_range(300.0..=560.0)
        .show(ctx, |ui| {
            ui.heading("Variant Effect");
            ui.label(
                egui::RichText::new("HGVS variant parser · valenx-variant-effect")
                    .weak()
                    .small(),
            );
            ui.separator();
            let s = &mut app.variant_effect;
            ui.horizontal(|ui| {
                ui.label("variant");
                ui.add(
                    egui::TextEdit::singleline(&mut s.input)
                        .desired_width(160.0)
                        .hint_text("p.R273H or c.817C>T"),
                );
            });
            if ui.button("▶ Parse").clicked() {
                let r = parse(&s.input).map_err(|e| e.to_string());
                s.parsed = Some(r);
            }
            ui.label(
                egui::RichText::new("Examples:  p.R273H · p.Arg273His · c.817C>T")
                    .weak()
                    .small(),
            );
            if let Some(res) = &s.parsed {
                ui.separator();
                match res {
                    Ok(v) => {
                        ui.colored_label(egui::Color32::from_rgb(80, 220, 120), describe(v));
                    }
                    Err(e) => {
                        ui.colored_label(egui::Color32::from_rgb(220, 120, 80), format!("parse error: {e}"));
                    }
                }
            }
        });
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_a_protein_substitution() {
        let v = parse("p.R273H").expect("parse");
        assert!(matches!(v, Variant::ProteinSub { pos: 273, .. }), "got {v:?}");
        assert!(describe(&v).contains("273"));
    }

    #[test]
    fn parses_a_coding_substitution() {
        let v = parse("c.817C>T").expect("parse");
        assert!(matches!(v, Variant::CodingSub { pos: 817, .. }), "got {v:?}");
    }

    #[test]
    fn rejects_garbage() {
        assert!(parse("not a variant").is_err());
    }
}
