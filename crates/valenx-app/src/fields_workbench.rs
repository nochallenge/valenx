//! The right-side **Field Statistics Workbench** panel — descriptive
//! statistics over a pasted number list, via `valenx-fields`.
//!
//! Mirrors the springs / gears / geomatics / piping / collision / sheet
//! metal workbenches: a resizable [`egui::SidePanel`] gated on
//! `crate::ValenxApp::show_fields_workbench`, toggled from the View
//! menu. The form takes a whitespace/comma-separated list of numbers; the
//! "Compute" button parses them into a scalar [`valenx_fields::Field`] and
//! reports the central-tendency, dispersion, and shape statistics
//! (mean / median / variance / std dev / rms / skewness / excess kurtosis
//! / coefficient of variation / min·max) as a monospace readout.

use eframe::egui;

use valenx_fields::integrals::{
    field_coefficient_of_variation, field_excess_kurtosis, field_mean, field_median, field_min_max,
    field_peak_to_peak, field_rms, field_skewness, field_std_dev, field_sum, field_variance,
};
use valenx_fields::units::DIMENSIONLESS;
use valenx_fields::{Field, FieldKind, Location, RegionRef, TimeKey};

use crate::ValenxApp;

/// Persistent form + result state for the Field Statistics Workbench.
pub struct FieldsWorkbenchState {
    /// The raw text the user pastes — numbers separated by whitespace
    /// and/or commas.
    text: String,
    /// Formatted statistics readout (empty until the first compute).
    result: String,
    /// Validation / parse error, if any.
    error: Option<String>,
}

impl Default for FieldsWorkbenchState {
    fn default() -> Self {
        Self {
            text: "1 2 3 4 5".to_string(),
            result: String::new(),
            error: None,
        }
    }
}

impl FieldsWorkbenchState {
    /// The user-visible captions of every control the agent bridge can set via
    /// `SetControl` (see [`crate::agent_commands`]). The only editable control is
    /// the numbers text box (the rest of the panel is computed output).
    pub fn agent_control_names() -> &'static [&'static str] {
        &["numbers"]
    }

    /// Set one labelled control by caption, for the agent `SetControl` bridge.
    /// The single settable field is the `numbers` text box (whitespace/comma-
    /// separated list); it reads `AgentValue::as_str`. Fail-loud on an unknown
    /// caption or a value of the wrong type — never a panic, no field written on
    /// error.
    pub fn agent_set(
        &mut self,
        name: &str,
        value: &crate::agent_commands::AgentValue,
    ) -> Result<(), String> {
        match name {
            "numbers" => self.text = value.as_str()?.to_string(),
            other => return Err(format!("unknown Field Statistics control: {other:?}")),
        }
        Ok(())
    }
}

/// Wrap a `Vec<f64>` in the canonical dimensionless scalar [`Field`] (the
/// same shape the crate's own tests use). `valenx-fields` exposes no
/// dedicated constructor, but every field is public, so this is a plain
/// struct literal of public items.
fn scalar_field(data: Vec<f64>) -> Field {
    Field {
        name: "workbench".to_string(),
        kind: FieldKind::Scalar,
        location: Location::OnNode,
        region: RegionRef("workbench".to_string()),
        units: DIMENSIONLESS,
        time: TimeKey::Steady,
        data,
        range: None,
    }
}

/// Draw the Field Statistics Workbench right-side panel. A no-op when the
/// `show_fields_workbench` toggle is off.
pub fn draw_fields_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_fields_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_fields_workbench",
        "Field Statistics",
        |app, ui| {
            ui.label(
                egui::RichText::new("descriptive statistics · valenx-fields")
                    .weak()
                    .small(),
            );
            ui.separator();

            let s = &mut app.fields;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(
                        egui::RichText::new("Numbers (whitespace or comma separated)").strong(),
                    );
                    ui.add(
                        egui::TextEdit::multiline(&mut s.text)
                            .desired_rows(4)
                            .font(egui::TextStyle::Monospace),
                    );

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("\u{25B6} Compute").strong())
                        .clicked()
                    {
                        run_fields(s);
                    }

                    if let Some(e) = &s.error {
                        ui.add_space(4.0);
                        ui.colored_label(egui::Color32::from_rgb(220, 90, 90), e);
                    }

                    if !s.result.is_empty() {
                        ui.separator();
                        ui.label(egui::RichText::new("Statistics").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        },
    );
    if close {
        app.show_fields_workbench = false;
    }
}

/// Parse the pasted numbers, build a scalar [`Field`], compute the
/// descriptive statistics, and format the readout. Extracted from the
/// draw closure so it is unit-testable.
fn run_fields(s: &mut FieldsWorkbenchState) {
    s.error = None;

    let tokens: Vec<&str> = s
        .text
        .split(|c: char| c.is_whitespace() || c == ',')
        .filter(|t| !t.is_empty())
        .collect();
    if tokens.is_empty() {
        s.error = Some("enter at least one number".into());
        return;
    }
    let mut data = Vec::with_capacity(tokens.len());
    for t in &tokens {
        match t.parse::<f64>() {
            Ok(v) if v.is_finite() => data.push(v),
            _ => {
                s.error = Some(format!("could not parse '{t}' as a finite number"));
                return;
            }
        }
    }

    let field = scalar_field(data);
    let n = field.data.len();
    let (min, max) = field_min_max(&field).unwrap_or((0.0, 0.0));

    s.result = format!(
        "samples        : {}\n\
         sum            : {:.6}\n\
         mean           : {:.6}\n\
         median         : {:.6}\n\
         min / max      : {:.6} / {:.6}\n\
         peak-to-peak   : {:.6}\n\
         variance (pop) : {:.6}\n\
         std dev        : {:.6}\n\
         rms            : {:.6}\n\
         skewness       : {:.6}\n\
         excess kurtosis: {:.6}\n\
         coeff of var   : {:.6}",
        n,
        field_sum(&field),
        field_mean(&field),
        field_median(&field),
        min,
        max,
        field_peak_to_peak(&field),
        field_variance(&field),
        field_std_dev(&field),
        field_rms(&field),
        field_skewness(&field),
        field_excess_kurtosis(&field),
        field_coefficient_of_variation(&field),
    );
}

/// Build the **Field Statistics** result card for the Workbench+Agent bridge —
/// a DATA-ONLY [`crate::WorkspaceProduct`] (`mesh: None`) whose `lines` are the
/// genuine descriptive statistics ([`run_fields`]) over the canonical default
/// dataset (`1 2 3 4 5`). Registered as the `"fields"` producer in
/// [`crate::products_registry::lookup`]; the tile renders it as a text card (the
/// same mesh-less path the `dna` card uses), not a 3-D view.
pub(crate) fn fields_product() -> crate::WorkspaceProduct {
    let mut s = FieldsWorkbenchState::default();
    run_fields(&mut s);
    crate::WorkspaceProduct {
        title: "Field Statistics".into(),
        lines: crate::products_registry::lines_from_readout(&s.result),
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
    fn default_state_is_idle() {
        let s = FieldsWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn compute_default_dataset() {
        let mut s = FieldsWorkbenchState::default();
        run_fields(&mut s);
        assert!(s.error.is_none());
        assert!(!s.result.is_empty());
        assert!(s.result.contains("mean"));
        assert!(s.result.contains("median"));
        assert!(s.result.contains("variance"));
        assert!(s.result.contains("skewness"));
        // [1,2,3,4,5]: mean 3, median 3 → "3.000000"; population variance 2
        // → "2.000000"; 5 samples.
        assert!(s.result.contains("samples        : 5"));
        assert!(s.result.contains("3.000000"));
        assert!(s.result.contains("2.000000"));
        // Recompute via the backend on the same scalar Field to confirm the
        // construction path + the worked values.
        let f = scalar_field(vec![1.0, 2.0, 3.0, 4.0, 5.0]);
        assert!((field_mean(&f) - 3.0).abs() < 1e-12);
        assert!((field_median(&f) - 3.0).abs() < 1e-12);
        assert!((field_variance(&f) - 2.0).abs() < 1e-12);
    }

    #[test]
    fn compute_accepts_comma_and_newline_separators() {
        let mut s = FieldsWorkbenchState {
            text: "10, 20\n30".to_string(),
            ..Default::default()
        };
        run_fields(&mut s);
        assert!(s.error.is_none());
        assert!(s.result.contains("samples        : 3"));
        // mean of 10,20,30 = 20.
        assert!(s.result.contains("20.000000"));
    }

    #[test]
    fn agent_set_sets_numbers_unknown_and_type_mismatch_err() {
        use crate::agent_commands::AgentValue;
        let mut s = FieldsWorkbenchState::default();
        // The numbers text box accepts a string and is then computable.
        s.agent_set("numbers", &AgentValue::Str("10 20 30".into()))
            .unwrap();
        assert_eq!(s.text, "10 20 30");
        run_fields(&mut s);
        assert!(s.error.is_none());
        assert!(s.result.contains("samples        : 3"));
        // Unknown caption -> Err.
        assert!(s
            .agent_set("no such control", &AgentValue::Str("x".into()))
            .is_err());
        // Type mismatch (number into the string field) -> Err, field untouched.
        assert!(s.agent_set("numbers", &AgentValue::Int(7)).is_err());
        assert_eq!(s.text, "10 20 30", "rejected set leaves field untouched");
    }

    #[test]
    fn compute_rejects_empty_and_non_numeric() {
        let mut empty = FieldsWorkbenchState {
            text: "   ".to_string(),
            ..Default::default()
        };
        run_fields(&mut empty);
        assert!(empty.error.is_some());
        assert!(empty.result.is_empty());

        let mut bad = FieldsWorkbenchState {
            text: "1 2 foo 3".to_string(),
            ..Default::default()
        };
        run_fields(&mut bad);
        assert!(bad.error.is_some());
        assert!(bad.result.is_empty());
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    /// Render the whole workbench panel once in a headless egui context.
    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_fields_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_fields_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_fields_workbench = true;
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_a_result_and_error_without_panic() {
        let mut app = ValenxApp::default();
        app.show_fields_workbench = true;
        run_fields(&mut app.fields);
        app.fields.error = Some("invalid input".to_string());
        draw_workbench(&mut app);
    }
}
