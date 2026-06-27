//! Variant-effect workbench — HGVS variant parsing on `valenx-variant-effect`.
//!
//! A right-side panel that parses HGVS-style substitution strings (protein
//! `p.R273H` / `p.Arg273His`, or coding `c.817C>T`) into structured
//! [`Variant`]s. Accepts a **batch** — one variant per line — and reports
//! each line's parsed components or a clear parse error, with an OK/error
//! tally. Headless-testable.

use eframe::egui;

use valenx_variant_effect::{parse, Variant};

use crate::agent_commands::AgentValue;
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

impl VariantEffectWorkbenchState {
    /// The user-visible captions of every control the agent bridge can set via
    /// `SetControl` (see [`crate::agent_commands`]). Returned by `ListControls`
    /// so an agent can discover the name space. This panel has a single settable
    /// control — the multi-line variant batch — addressed by the heading the
    /// panel draws above it. The `Parse` button is an action, not a settable
    /// control, so it is not listed.
    pub fn agent_control_names() -> &'static [&'static str] {
        &["Variants (one per line)"]
    }

    /// Set the variant-batch text field by its caption, for the agent
    /// `SetControl` bridge. Fail-loud: an unknown caption or a non-string value
    /// returns `Err(String)` (the bridge turns it into a `warn` feed note) —
    /// never a panic, and no field is written on error. The batch is free text
    /// ([`AgentValue::as_str`]); the agent then triggers `Parse` to evaluate it
    /// (so setting this clears the stale `results` until the next parse).
    pub fn agent_set(&mut self, name: &str, value: &AgentValue) -> Result<(), String> {
        match name {
            "Variants (one per line)" => {
                self.input = value.as_str()?.to_string();
                // The displayed results no longer correspond to the new input;
                // drop them so the panel does not show stale parses.
                self.results = None;
            }
            other => return Err(format!("unknown variant-effect control: {other:?}")),
        }
        Ok(())
    }

    /// Read-only readout for the agent `ReadReadout` bridge and the product
    /// self-test ([`crate::self_test`]): one line per parsed variant
    /// (`<line> → <description>` on success, `<line> → error: …` on failure), or
    /// `None` before the first parse. Mirrors the panel's per-line result list.
    pub fn agent_readout(&self) -> Option<String> {
        let results = self.results.as_ref()?;
        Some(
            results
                .iter()
                .map(|(line, r)| match r {
                    Ok(v) => format!("{line} → {}", describe(v)),
                    Err(e) => format!("{line} → error: {e}"),
                })
                .collect::<Vec<_>>()
                .join("\n"),
        )
    }
}

/// Run the HGVS batch parse (the in-panel **Parse** action). Factored out so the
/// button and the product self-test ([`crate::self_test`]) share one path.
pub(crate) fn run(app: &mut ValenxApp) {
    app.variant_effect.results = Some(parse_batch(&app.variant_effect.input));
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
    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_variant_effect_workbench",
        "Variant Effect",
        |app, ui| {
            ui.label(
                egui::RichText::new("HGVS variant parser · valenx-variant-effect")
                    .weak()
                    .small(),
            );
            ui.separator();
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
        },
    );
    if close {
        app.show_variant_effect_workbench = false;
    }
}

/// Build the **Variant Effect** result card for the Workbench+Agent bridge — a
/// DATA-ONLY [`crate::WorkspaceProduct`] (`mesh: None`) whose `lines` are the
/// genuine HGVS parse results ([`parse_batch`]) for the canonical default batch
/// (`p.R273H` / `p.Arg249Ser` / `c.817C>T`): a parsed/error tally followed by one
/// described row per variant. Registered as the `"variant_effect"` producer in
/// [`crate::products_registry::lookup`]; the tile renders it as a text card, not
/// a 3-D view. The rows mirror the panel's own per-variant readout.
pub(crate) fn variant_effect_product() -> crate::WorkspaceProduct {
    let results = parse_batch(&VariantEffectWorkbenchState::default().input);
    let ok = results.iter().filter(|(_, r)| r.is_ok()).count();
    let err = results.len() - ok;
    let mut lines = vec![format!(
        "{} variant(s) · {ok} parsed · {err} error(s)",
        results.len()
    )];
    for (line, res) in &results {
        match res {
            Ok(v) => lines.push(format!("{line}  \u{2192}  {}", describe(v))),
            Err(e) => lines.push(format!("{line}  \u{2192}  parse error: {e}")),
        }
    }
    crate::WorkspaceProduct {
        title: "Variant Effect".into(),
        lines,
        mesh: None,
        vertex_colors: None,
        camera: Default::default(),
        kind2d: None,
        last_export: None,
        image: None,
        image_texture: None,
        animation: None,
    }
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

    // ---- agent_set / agent_control_names (the SetControl bridge) ----

    #[test]
    #[allow(clippy::field_reassign_with_default)]
    fn agent_set_sets_batch_and_rejects_unknown_and_typemismatch() {
        let mut s = VariantEffectWorkbenchState::default();
        // Seed a stale results table to confirm agent_set clears it.
        s.results = Some(parse_batch("p.R273H"));

        // The settable text control, verified via state.
        s.agent_set(
            "Variants (one per line)",
            &AgentValue::Str("p.G12D\nc.35G>A".into()),
        )
        .expect("set the variant batch");
        assert_eq!(s.input, "p.G12D\nc.35G>A");
        assert!(s.results.is_none(), "stale results must be cleared on set");
        // The new text parses through the same pipeline the panel uses.
        let out = parse_batch(&s.input);
        assert_eq!(out.len(), 2);
        assert!(out.iter().all(|(_, r)| r.is_ok()));

        // Unknown caption -> Err (not a panic).
        assert!(s.agent_set("nope", &AgentValue::Str("x".into())).is_err());
        // Type mismatch: the text caption fed a number -> Err.
        assert!(s
            .agent_set("Variants (one per line)", &AgentValue::Int(7))
            .is_err());

        // The single advertised control name is settable.
        assert_eq!(
            VariantEffectWorkbenchState::agent_control_names(),
            &["Variants (one per line)"]
        );
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
