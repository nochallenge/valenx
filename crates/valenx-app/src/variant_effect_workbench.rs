//! Variant-effect workbench — HGVS variant parsing on `valenx-variant-effect`.
//!
//! A right-side panel that parses HGVS-style substitution strings (protein
//! `p.R273H` / `p.Arg273His`, or coding `c.817C>T`) into structured
//! [`Variant`]s. Accepts a **batch** — one variant per line — and reports
//! each line's parsed components or a clear parse error, with an OK/error
//! tally. Headless-testable.

use eframe::egui;

use valenx_variant_effect::{parse, Variant};

use crate::ValenxApp;

/// Persistent state for the variant-effect workbench.
pub struct VariantEffectWorkbenchState {
    /// One HGVS variant per line.
    input: String,
    /// Per-line parse results (the trimmed line + its parse), built on "Parse".
    results: Option<Vec<(String, Result<Variant, String>)>>,
}

impl Default for VariantEffectWorkbenchState {
    fn default() -> Self {
        Self {
            input: "p.R273H\np.Arg249Ser\nc.817C>T".to_string(),
            results: None,
        }
    }
}

/// A human-readable description of a parsed variant.
fn describe(v: &Variant) -> String {
    match v {
        Variant::ProteinSub { wt, pos, mt } => {
            format!(
                "protein substitution {}{}{} — residue {pos}",
                wt.as_char(),
                pos,
                mt.as_char()
            )
        }
        Variant::CodingSub { pos, wt, mt } => {
            format!(
                "coding substitution c.{pos}{}>{}",
                wt.as_char(),
                mt.as_char()
            )
        }
    }
}

/// Parse a batch of variants — one per line, blank lines skipped — into a
/// `(line, result)` list. Extracted from the draw closure so it is
/// unit-testable.
fn parse_batch(input: &str) -> Vec<(String, Result<Variant, String>)> {
    input
        .lines()
        .map(str::trim)
        .filter(|l| !l.is_empty())
        .map(|l| (l.to_string(), parse(l).map_err(|e| e.to_string())))
        .collect()
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
            if crate::workbench_ui::header(
                ui,
                "Variant Effect",
                "HGVS variant parser · valenx-variant-effect",
            ) {
                app.show_variant_effect_workbench = false;
            }
            let s = &mut app.variant_effect;
            ui.label(egui::RichText::new("Variants (one per line)").strong());
            ui.add(
                egui::TextEdit::multiline(&mut s.input)
                    .desired_rows(4)
                    .desired_width(f32::INFINITY)
                    .hint_text("p.R273H\np.Arg249Ser\nc.817C>T")
                    .font(egui::TextStyle::Monospace),
            );
            if ui.button("\u{25B6} Parse").clicked() {
                s.results = Some(parse_batch(&s.input));
            }
            ui.label(
                egui::RichText::new("Examples:  p.R273H · p.Arg273His · c.817C>T")
                    .weak()
                    .small(),
            );
            if let Some(results) = &s.results {
                ui.separator();
                let ok = results.iter().filter(|(_, r)| r.is_ok()).count();
                let err = results.len() - ok;
                ui.label(
                    egui::RichText::new(format!(
                        "{} variant(s) · {ok} parsed · {err} error(s)",
                        results.len()
                    ))
                    .strong(),
                );
                egui::ScrollArea::vertical()
                    .auto_shrink([false, false])
                    .max_height(220.0)
                    .show(ui, |ui| {
                        for (line, res) in results {
                            match res {
                                Ok(v) => {
                                    ui.colored_label(
                                        egui::Color32::from_rgb(80, 220, 120),
                                        format!("{line}  \u{2192}  {}", describe(v)),
                                    );
                                }
                                Err(e) => {
                                    ui.colored_label(
                                        egui::Color32::from_rgb(220, 120, 80),
                                        format!("{line}  \u{2192}  parse error: {e}"),
                                    );
                                }
                            }
                        }
                    });
            }
        });
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_a_protein_substitution() {
        let v = parse("p.R273H").expect("parse");
        assert!(
            matches!(v, Variant::ProteinSub { pos: 273, .. }),
            "got {v:?}"
        );
        assert!(describe(&v).contains("273"));
    }

    #[test]
    fn parses_a_coding_substitution() {
        let v = parse("c.817C>T").expect("parse");
        assert!(
            matches!(v, Variant::CodingSub { pos: 817, .. }),
            "got {v:?}"
        );
    }

    #[test]
    fn rejects_garbage() {
        assert!(parse("not a variant").is_err());
    }

    #[test]
    fn parse_batch_splits_lines_and_tallies() {
        // Four lines (one blank skipped → 3 entries): 2 valid, 1 garbage.
        let out = parse_batch("p.R273H\n\nc.817C>T\nnot a variant");
        assert_eq!(out.len(), 3);
        assert!(out[0].1.is_ok());
        assert!(out[1].1.is_ok());
        assert!(out[2].1.is_err());
        // The line text is preserved (trimmed).
        assert_eq!(out[0].0, "p.R273H");
    }

    #[test]
    fn parse_batch_default_is_all_valid() {
        let out = parse_batch(&VariantEffectWorkbenchState::default().input);
        assert_eq!(out.len(), 3);
        assert!(out.iter().all(|(_, r)| r.is_ok()));
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_variant_effect_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_variant_effect_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_batch_results_without_panic() {
        let mut app = ValenxApp::default();
        app.show_variant_effect_workbench = true;
        app.variant_effect.results = Some(parse_batch("p.R273H\nnot a variant"));
        draw_workbench(&mut app);
    }
}
