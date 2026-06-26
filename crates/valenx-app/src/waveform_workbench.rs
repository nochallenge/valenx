//! The right-side **Waveform (VCD viewer)** workbench — an in-house digital
//! waveform-capture panel over the [`valenx_waveform`] crate (a
//! self-contained Value Change Dump parser; Valenx's digital-oscilloscope /
//! logic-analyzer input).
//!
//! The user pastes (or edits) VCD source into a labelled, AI-settable
//! multiline text box, seeded with a small clock+counter sample. **Parse**
//! calls [`valenx_waveform::Waveform::parse`] and lists every signal (name +
//! bit width), each signal's transition count, and the overall time range.
//!
//! It mirrors the other real-time workbenches (`thermo_workbench`): a
//! [`crate::workbench_chrome::workbench_shell`] panel gated on
//! [`crate::ValenxApp::show_waveform_workbench`], toggled from the View menu
//! and openable by the agent bridge under the workbench id `"waveform"`. The
//! bridge can set the controls (`agent_set` / `agent_control_names`), read a
//! status line (`agent_readout`), and fire the parse via the RunCommand id
//! `waveform.parse`.

use eframe::egui;

use valenx_waveform::Waveform;

use crate::ValenxApp;

// ---------------------------------------------------------------------------
// Seed sample
// ---------------------------------------------------------------------------

/// A small, valid VCD sample — a 1-bit clock toggling and a 2-bit counter,
/// dumped at `#0 / #5 / #10` — used as the initial textbox contents (the
/// same trace the crate documents). Parses to two signals over time range
/// `(0, 10)`.
const SAMPLE_VCD: &str = "\
$timescale 1ns $end
$scope module top $end
$var wire 1 ! clk $end
$var wire 2 # cnt [1:0] $end
$upscope $end
$enddefinitions $end
#0
0!
b00 #
#5
1!
b01 #
#10
0!
b10 #
";

// ---------------------------------------------------------------------------
// Workbench state
// ---------------------------------------------------------------------------

/// One parsed signal row for the readout: its name, bit width, and the
/// number of value-change transitions recorded for it.
struct SignalRow {
    /// The signal's declared name.
    name: String,
    /// The signal's bit width (1 for a scalar wire).
    width: u32,
    /// How many value-change transitions the signal has.
    transitions: usize,
}

/// The result of a parse: the per-signal rows + the overall time range.
struct WaveformResult {
    /// One row per parsed signal (name / width / transition count).
    signals: Vec<SignalRow>,
    /// The `(first, last)` change time across all signals, or `None` for an
    /// empty trace.
    time_range: Option<(u64, u64)>,
}

/// Persistent state for the Waveform workbench: the editable VCD source and
/// the latest parse result.
pub struct WaveformWorkbenchState {
    /// The VCD source text the user edits / pastes (seeded with [`SAMPLE_VCD`]).
    pub source: String,

    /// Latest parse result, or `None` before the first parse.
    result: Option<WaveformResult>,
    /// Last error string (shown in the panel), cleared on a good parse.
    error: Option<String>,
}

impl Default for WaveformWorkbenchState {
    fn default() -> Self {
        // Seeded with the clock+counter sample. Parse is NOT run until the
        // user (or the bridge) presses Parse.
        Self {
            source: SAMPLE_VCD.to_string(),
            result: None,
            error: None,
        }
    }
}

impl WaveformWorkbenchState {
    /// Parse the current VCD source and store the result (or an error).
    /// Factored out so the in-panel **Parse** button and the
    /// `waveform.parse` bridge id share one path.
    fn parse_now(&mut self) {
        match self.try_parse() {
            Ok(res) => {
                self.result = Some(res);
                self.error = None;
            }
            Err(e) => {
                self.error = Some(e);
                // Keep any previous result on screen so a failed edit doesn't
                // blank the readout; the error line explains why.
            }
        }
    }

    /// The fallible parse pipeline. Separated so `parse_now` can map a typed
    /// error into the panel's error line.
    fn try_parse(&self) -> Result<WaveformResult, String> {
        let wf = Waveform::parse(&self.source).map_err(|e| e.to_string())?;
        let signals = wf
            .signals()
            .iter()
            .map(|s| SignalRow {
                name: s.name.clone(),
                width: s.width,
                transitions: s.transitions().len(),
            })
            .collect();
        Ok(WaveformResult {
            signals,
            time_range: wf.time_range(),
        })
    }

    /// The user-visible captions of every control the agent bridge can set
    /// via `SetControl` (returned by `ListControls`). The single settable
    /// control is the VCD source text.
    pub fn agent_control_names() -> &'static [&'static str] {
        &["VCD source"]
    }

    /// Set one labelled control by its caption, for the agent `SetControl`
    /// bridge. Fail-loud on an unknown caption / wrong type; no state is
    /// written on error and nothing panics. `VCD source` reads a string (the
    /// raw VCD text to parse).
    pub fn agent_set(
        &mut self,
        name: &str,
        value: &crate::agent_commands::AgentValue,
    ) -> Result<(), String> {
        match name {
            "VCD source" => {
                self.source = value.as_str()?.to_string();
            }
            other => return Err(format!("unknown waveform control: {other:?}")),
        }
        Ok(())
    }

    /// The current readout text for the agent `ReadReadout` bridge: the
    /// parsed signal count, each signal's name/width/transition count, and
    /// the overall time range once a parse exists. `Some` once parsed,
    /// `None` before the first parse (or after a parse error with no prior
    /// result).
    pub fn agent_readout(&self) -> Option<String> {
        if let Some(err) = &self.error {
            if self.result.is_none() {
                return Some(format!("Waveform parse failed: {err}"));
            }
        }
        let r = self.result.as_ref()?;
        let range = match r.time_range {
            Some((lo, hi)) => format!("[{lo}, {hi}]"),
            None => "empty".to_string(),
        };
        let sigs = r
            .signals
            .iter()
            .map(|s| format!("{}(w{}, {}x)", s.name, s.width, s.transitions))
            .collect::<Vec<_>>()
            .join(", ");
        Some(format!(
            "Waveform \u{00B7} {} signal(s) \u{00B7} time {range} \u{00B7} {sigs}",
            r.signals.len(),
        ))
    }
}

// ---------------------------------------------------------------------------
// Bridge run action (parse)
// ---------------------------------------------------------------------------

/// Run the parse (the in-panel **Parse** action). Factored out so the button
/// and the `waveform.parse` bridge id share one path.
pub(crate) fn run(app: &mut ValenxApp) {
    app.waveform.parse_now();
}

// ---------------------------------------------------------------------------
// Entry point
// ---------------------------------------------------------------------------

/// Draw the Waveform workbench. A no-op unless toggled on via
/// View → Waveform (VCD viewer).
pub fn draw_waveform_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_waveform_workbench {
        return;
    }
    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_waveform_workbench",
        "Waveform (VCD viewer) \u{2014} digital logic-analyzer / oscilloscope capture",
        waveform_workbench_body,
    );
    if close {
        app.show_waveform_workbench = false;
    }
}

// ---------------------------------------------------------------------------
// Workbench body
// ---------------------------------------------------------------------------

fn waveform_workbench_body(app: &mut ValenxApp, ui: &mut egui::Ui) {
    egui::ScrollArea::vertical()
        .auto_shrink([false, false])
        .show(ui, |ui| waveform_workbench_body_inner(app, ui));
}

fn waveform_workbench_body_inner(app: &mut ValenxApp, ui: &mut egui::Ui) {
    ui.label(
        egui::RichText::new(
            "In-house digital waveform capture \u{00B7} a self-contained Value Change Dump (VCD) \
             parser [paste or edit a VCD trace, then Parse to list the signals (name + bit width), \
             each signal's transition count, and the overall time range].",
        )
        .weak()
        .small(),
    );
    ui.separator();

    let s = &mut app.waveform;

    // --- VCD source ---------------------------------------------------------
    let lbl = ui.label("VCD source");
    egui::ScrollArea::vertical()
        .id_source("waveform_vcd_scroll")
        .max_height(220.0)
        .show(ui, |ui| {
            ui.add(
                egui::TextEdit::multiline(&mut s.source)
                    .code_editor()
                    .desired_rows(12)
                    .desired_width(f32::INFINITY),
            )
            .labelled_by(lbl.id);
        });

    // --- Parse --------------------------------------------------------------
    ui.add_space(6.0);
    if ui
        .button("\u{25B6} Parse")
        .on_hover_text(
            "Parse the VCD source and list every signal (name + bit width), each signal's \
             transition count, and the overall time range.",
        )
        .clicked()
    {
        s.parse_now();
    }

    ui.add_space(6.0);
    ui.separator();
    draw_result(s, ui);
}

// ---------------------------------------------------------------------------
// Result render
// ---------------------------------------------------------------------------

fn draw_result(s: &WaveformWorkbenchState, ui: &mut egui::Ui) {
    if let Some(err) = &s.error {
        ui.label(
            egui::RichText::new(format!("Parse error: {err}"))
                .color(egui::Color32::from_rgb(220, 110, 90))
                .strong(),
        );
        ui.add_space(4.0);
    }

    let Some(r) = s.result.as_ref() else {
        ui.label(
            egui::RichText::new(
                "No result yet \u{2014} edit / paste a VCD trace, then press Parse.",
            )
            .italics()
            .weak(),
        );
        return;
    };

    let range = match r.time_range {
        Some((lo, hi)) => format!("[{lo}, {hi}]"),
        None => "empty".to_string(),
    };
    ui.label(
        egui::RichText::new(format!(
            "Parsed {} signal(s) \u{00B7} time range {range}",
            r.signals.len()
        ))
        .strong()
        .color(egui::Color32::from_rgb(150, 200, 230)),
    );
    ui.add_space(2.0);
    egui::Grid::new("waveform_signals")
        .num_columns(3)
        .striped(true)
        .show(ui, |ui| {
            ui.label(egui::RichText::new("signal").strong());
            ui.label(egui::RichText::new("width").strong());
            ui.label(egui::RichText::new("transitions").strong());
            ui.end_row();
            for sig in &r.signals {
                ui.label(&sig.name);
                ui.label(format!("{}", sig.width));
                ui.label(format!("{}", sig.transitions));
                ui.end_row();
            }
        });
}

// ---------------------------------------------------------------------------
// Tests (unit)
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_seed_is_the_sample_and_unparsed() {
        let s = WaveformWorkbenchState::default();
        assert!(s.source.contains("$timescale"));
        assert!(s.source.contains("clk"));
        assert!(s.result.is_none());
    }

    /// Ground truth: the seeded clock+counter sample parses to exactly the
    /// expected two signals over time range `(0, 10)`, with the clock
    /// carrying three transitions (the crate's documented trace).
    #[test]
    fn sample_parses_to_expected_signal_count() {
        let mut s = WaveformWorkbenchState::default();
        s.parse_now();
        assert!(s.error.is_none(), "the sample should parse: {:?}", s.error);
        let r = s.result.as_ref().expect("a result after parse");
        assert_eq!(r.signals.len(), 2, "the sample has clk + cnt");
        assert_eq!(r.time_range, Some((0, 10)), "time range");
        // The parser prefixes the declared `$scope module top` and re-attaches
        // any `[..]` bit-select — names are the fully-qualified `top.clk` /
        // `top.cnt[1:0]`; match on the leaf base (strip scope + bit-select),
        // mirroring the crate's own `signal_by_name`.
        let leaf = |row: &&SignalRow, want: &str| {
            let leaf = row.name.rsplit('.').next().unwrap_or(&row.name);
            leaf.split('[').next().map(str::trim) == Some(want)
        };
        let clk = r
            .signals
            .iter()
            .find(|s| leaf(s, "clk"))
            .expect("clk signal");
        assert_eq!(clk.width, 1, "clk is a 1-bit wire");
        assert_eq!(clk.transitions, 3, "clk: 0@0, 1@5, 0@10");
        let cnt = r
            .signals
            .iter()
            .find(|s| leaf(s, "cnt"))
            .expect("cnt signal");
        assert_eq!(cnt.width, 2, "cnt is a 2-bit wire");
    }

    #[test]
    fn parse_produces_a_readout() {
        let mut s = WaveformWorkbenchState::default();
        s.parse_now();
        let readout = s.agent_readout().expect("readout after parse");
        assert!(readout.contains("Waveform"), "readout: {readout}");
        assert!(
            readout.contains("2 signal"),
            "readout names count: {readout}"
        );
        assert!(readout.contains("clk"), "readout names clk: {readout}");
    }

    #[test]
    fn agent_set_replaces_source() {
        use crate::agent_commands::AgentValue;
        let mut s = WaveformWorkbenchState::default();
        let tiny = "$var wire 1 ! a $end\n$enddefinitions $end\n#0\n0!\n".to_string();
        s.agent_set("VCD source", &AgentValue::Str(tiny.clone()))
            .expect("source");
        assert_eq!(s.source, tiny);
        s.parse_now();
        let r = s.result.as_ref().expect("parses");
        assert_eq!(r.signals.len(), 1, "the tiny trace has one signal");
    }

    #[test]
    fn agent_set_rejects_bad_values() {
        use crate::agent_commands::AgentValue;
        let mut s = WaveformWorkbenchState::default();
        // A numeric value is the wrong type for the text control.
        assert!(s.agent_set("VCD source", &AgentValue::Float(1.0)).is_err());
        assert!(s.agent_set("bogus", &AgentValue::Str("x".into())).is_err());
    }

    #[test]
    fn parse_error_surfaces_for_garbage_without_panic() {
        let mut s = WaveformWorkbenchState::default();
        s.source = "this is not a vcd file at all".to_string();
        s.parse_now();
        // Either a clean parse with no signals or a surfaced error — never a
        // panic. The garbage here has no `$enddefinitions`, so it errors.
        assert!(
            s.error.is_some() || s.result.as_ref().is_some_and(|r| r.signals.is_empty()),
            "garbage should not yield a populated result"
        );
    }

    #[test]
    fn readout_is_none_before_parse() {
        let s = WaveformWorkbenchState::default();
        assert!(s.agent_readout().is_none(), "no readout before a parse");
    }

    #[test]
    fn run_bridge_helper_parses_through_app() {
        let mut app = ValenxApp::default();
        run(&mut app);
        assert!(
            app.waveform.result.is_some(),
            "the waveform.parse bridge helper should produce a result"
        );
    }

    #[test]
    fn control_names_are_listed() {
        let names = WaveformWorkbenchState::agent_control_names();
        assert!(names.contains(&"VCD source"), "missing VCD source control");
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;
    use egui::accesskit::{Node, NodeId, Role};

    fn draw_and_collect_nodes(app: &mut ValenxApp) -> Vec<(NodeId, Node)> {
        let ctx = egui::Context::default();
        ctx.enable_accesskit();
        let out = ctx.run(egui::RawInput::default(), |ctx| {
            draw_waveform_workbench(app, ctx);
        });
        out.platform_output
            .accesskit_update
            .expect("accesskit tree is produced when enabled")
            .nodes
    }

    fn has_named_node(nodes: &[(NodeId, Node)], name: &str) -> bool {
        nodes.iter().any(|(_, n)| n.name() == Some(name))
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_waveform_workbench);
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_waveform_workbench(&mut app, ctx);
        });
        // No panic = pass.
    }

    #[test]
    fn workbench_draws_when_shown() {
        let mut app = ValenxApp::default();
        app.show_waveform_workbench = true;
        let nodes = draw_and_collect_nodes(&mut app);
        assert!(!nodes.is_empty(), "a shown workbench produces a11y nodes");
    }

    #[test]
    fn numeric_controls_are_named_and_associated() {
        // The single interactive control is the multiline VCD `TextEdit`
        // (a TextInput / MultilineTextInput); it must be `labelled_by` its
        // caption so an AI / screen reader can find it by caption text.
        let mut app = ValenxApp::default();
        app.show_waveform_workbench = true;
        let nodes = draw_and_collect_nodes(&mut app);

        let by_id: std::collections::HashMap<NodeId, &Node> =
            nodes.iter().map(|(id, n)| (*id, n)).collect();

        let text_inputs: Vec<&Node> = nodes
            .iter()
            .map(|(_, n)| n)
            .filter(|n| matches!(n.role(), Role::TextInput | Role::MultilineTextInput))
            .collect();
        assert!(
            !text_inputs.is_empty(),
            "expected the VCD source text box as a text input"
        );
        assert!(
            text_inputs.iter().all(|n| !n.labelled_by().is_empty()),
            "the VCD TextEdit must be labelled_by a caption (AI-drivable name)"
        );
        assert!(
            text_inputs.iter().all(|n| {
                n.labelled_by()
                    .iter()
                    .any(|id| by_id.get(id).is_some_and(|t| t.name().is_some()))
            }),
            "the VCD TextEdit's labelled_by must point at a named caption node"
        );

        assert!(
            has_named_node(&nodes, "VCD source"),
            "caption 'VCD source' must be a named node in the a11y tree"
        );
    }
}
